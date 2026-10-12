# 理论与数值模型

[English](theory.md) | 中文  
[返回中文 README](../README.zh-CN.md)

> **作者:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **更新日期:** 2026-10-11
> **版本:** 8.5.2

本文档承接原长版 README 中的理论、模型和数值方法说明。仓库首页现在保持精简；实际运行步骤见 [`../USER_GUIDE.md`](../USER_GUIDE.md)，版本历史见 [`../CHANGELOG.md`](../CHANGELOG.md)。

## 1. 物理问题

仿真对象是在 Stokes 流中沉降的细菌胞体及其柔性纤毛 / 鞭毛：

- **胞体:** 刚性胶囊体、近似球体或导入的自定义三维点云。
- **纤毛 / 鞭毛:** 一根或多根不可伸长 Kirchhoff 弹性杆。
- **流体:** 零雷诺数 Stokes 流。
- **耦合:** 胞体与所有纤毛节点之间的水动力相互作用通过全局迁移率求解保留。

默认参数化胞体由两个半球端盖和一个圆柱中段组成。自定义几何模式下，胞体与纤毛必须同时来自用户导入的点云；当混用字段进入校验器时（内联/程序化 `cfg`，或经 `run_simulation.m` 传入的 JSON `cfg`），`validate_simulation_params.m` 会拒绝“自定义 + 参数化”的双向混用。注意：JSON 参数化模式下，`load_sim_config.m` **不会**将顶层 `custom_geometry` 读入 `cfg`，因此 JSON 中残留的该字段会在加载阶段被**静默忽略**而非触发校验——应从 JSON 中删除该块，或改用 `geometry_mode: "custom"`。

导入的纤毛曲线只转换一次为 canonical 等弦长折线：先去除连续重复点并建立 PCHIP 参数曲线，再用嵌套标量根求解获得参数单调、端点保持且每段弦长相同的节点。参考节点、全部材料 directors 和内禀曲率都来自同一折线。初始胞体姿态对全部参考 directors 施加同一个刚体旋转后再恢复 Euler 角，从而保证弯曲无应力杆的客观性。所有纤毛统一使用全局 \(E\) 与 \(r_{\mathrm{pili}}\)，不存在逐纤毛刚度参数。

参数化纤毛由 `pili_type_building.m` 生成沿 \(+x\) 的模板（弯曲平面为 \(xy\)），再由 `pili_building.m` 放置到胞体 14 个附着位点之一。自 v8.3.0 起，每个位点变换都是真旋转 \(R=[\mathbf{c}_1,\mathbf{c}_2,\mathbf{c}_1\times\mathbf{c}_2]\)（\(\mathbf{c}_1\) 为外法向、\(\mathbf{c}_2\) 为弯曲方向）加平移，可选的偏转矩阵 \(R_z(\theta)\)、\(R_x(\gamma)\) 均为右手约定。旧版在位点 1、2、5、7、10、11、14 采用 \(\det=-1\) 的坐标交换并使用左手 \(R_x\)，会镜像手性模板（锥形螺旋）及 `deflection_x` 的扭转方向；不带 `deflection_x` 的平面模板在两版中放置完全一致。

## 2. 主要方法

| 组件 | 方法 | 主要来源 |
| --- | --- | --- |
| 胞体水动力 | 胞体表面点上的正则化 Stokeslet 方法 (RSM) | Cortez et al. (2005) |
| 纤毛水动力 | 正则化 Stokeslet 线段法 (RSS)，线段精确积分 | Cortez (2018) |
| 纤毛扭转-流动耦合 | 正则化 Rotlet 线段 | Cortez (2018) |
| 纤毛弹性 | Euler 角 director 表示的 Kirchhoff 弹性杆 | Walker (2020) |
| 胞体姿态 | 四元数运动学 | 本项目实现 |
| 时间积分 | MATLAB `ode15s` 与事件驱动分段积分 | 本项目实现 |

RSS 的关键优势是对每根直线纤毛段上的正则化 Stokeslet 核函数进行精确积分，使段长与正则化半径解耦。典型纤毛可使用约 `10` 段每微米的离散密度，同时仍保留较完整的水动力耦合。

## 3. 无量纲化

以胞体半长轴作为特征长度：

| 物理量 | 符号 | 定义 |
| --- | --- | --- |
| 长度 | $L_0$ | 参数化：`major_axis`；自定义：`L_0_override` 或体几何建议值 |
| 弯曲刚度 | $EI$ | $E\pi r_{\mathrm{pili}}^4/4$ |
| 时间 | $T_0$ | $8\pi\mu L_0^4/EI$ |
| 速度 | $U_0$ | $L_0/T_0$ |
| 力 | $F_0$ | $EI/L_0^2$ |

关键无量纲参数包括：

- 弹性-重力数：$E_g = F_{\mathrm{sink}}L_0^2/EI$。
- 纤毛细长比：$\varepsilon = r_{\mathrm{pili}}/L_0$。
- 泊松比：$\sigma$，由用户给定。

沉降力（有量纲）为

```math
F_{\mathrm{sink}} = (\rho_{\mathrm{bac}} - \rho_{\mathrm{fluid}})\, V_{\mathrm{body}}\, g,
```

其中 $V_{\mathrm{body}}$ 在参数化模式下为胶囊体积，在自定义模式下为 alpha-shape 体积。重力通过 ODE 右端 $\mathbf{R}_{1:3} = -E_g\,\hat{\mathbf{g}}$ 进入，$\hat{\mathbf{g}} = [0,0,-1]^\mathsf{T}$。

根部与胞体间隙（有量纲）为

```math
\mathrm{dis} = 1.1\,(\varepsilon_{\mathrm{body}} + \varepsilon_{\mathrm{pili}})\,L_0,
```

传入 `pili_building` / `custom_pili_building`，避免纤毛根部正则化核重叠。

其他物理输入使用同一套尺度：

- 背景流速度：$\mathbf{u}_{\mathrm{bg}}^\ast = \mathbf{u}_{\mathrm{bg}}/U_0$。
- 积分后的主动力矩：$\mathbf{M}^\ast = \mathbf{M}_{\mathrm{phys}}L_{\mathrm{pili}}/EI$。
- 内禀曲率：使用与 Euler 角曲率场一致的无量纲弧长尺度。

## 4. 状态向量与几何重构

完整状态向量包含胞体位置、胞体姿态和每根纤毛的杆状态：

```math
\mathbf{Z} =
\left[
\mathbf{X}_{\mathrm{cm}},
\mathbf{q},
\boldsymbol{\theta}_1,
\boldsymbol{\phi}_1,
\boldsymbol{\psi}_1,
\ldots,
\boldsymbol{\theta}_M,
\boldsymbol{\phi}_M,
\boldsymbol{\psi}_M
\right].
```

胞体四元数 $\mathbf{q}$ 表示体坐标到实验室坐标的旋转，其演化为：

```math
\frac{d\mathbf{q}}{dt}
= \frac{1}{2}
\begin{bmatrix}
0 \\
\boldsymbol{\Omega}
\end{bmatrix}
\otimes \mathbf{q}.
```

其中 $\boldsymbol{\Omega}$ 为**实验室坐标系下的胞体角速度**，取自耦合求解的 `z_dot(4:6)`（`quat_kinematics.m`）。

每根纤毛由 Euler 角定义材料标架 `(d1, d2, d3)`，其中 `d3` 是纤毛段切向。[`reconstruct_geometry.m`](../src/geometry/reconstruct_geometry.m) 将状态 `Z` 转换为实验室坐标系下的胞体表面点、纤毛节点、director 标架和水动力求解所需的辅助几何量。

## 5. 水动力模型

### 5.1 胞体 RSM

胞体表面力使用 Cortez 正则化 Stokeslet 核函数：

```math
\mathbf{S}^{\varepsilon}(\mathbf{r})
=
\frac{
\left(\lVert\mathbf{r}\rVert^2 + 2\varepsilon^2\right)\mathbf{I}
+ \mathbf{r}\mathbf{r}^{\mathsf{T}}
}{
\left(\lVert\mathbf{r}\rVert^2 + \varepsilon^2\right)^{3/2}
}.
```

胞体正则化长度默认值为：

```math
\varepsilon_{\mathrm{body}} = h_{\mathrm{mean}}/2.6667.
```

其中 `h_mean` 是胞体表面点平均最近邻距离的无量纲值。该默认值来自 `N_body x eps_body` Stokes 球沉降标定扫描，也可通过 JSON 中的 `body.eps_body` 或 MATLAB 中的 `cfg.eps_body` 覆盖。

### 5.2 纤毛 RSS

对于从 $\mathbf{y}_0$ 到 $\mathbf{y}_1$ 的直线段，RSS 计算线性变化力密度产生的速度：

```math
\mathbf{u}(\mathbf{x})
=
\mathbf{M}_1\mathbf{f}_0
+ \mathbf{M}_2\mathbf{f}_1.
```

矩阵 `M1` 和 `M2` 已包含线段精确积分权重，因此 RSS 线段列不再额外乘梯形长度权重。所需的 `T_{n,q}` 积分量通过 Cortez (2018) 递推公式计算。

纤毛正则化长度通过标定因子 $\beta$（`hydro_radius_factor`）与物理半径相联系：

```math
\varepsilon_{\mathrm{pili}} = \beta\, r_{\mathrm{pili}}/L_0,
\qquad \beta = \sqrt{2}\ \text{（默认）}.
```

与胞体 RSM 正则化不同，该长度不是应趋近于零的纯数值参数，它决定纤毛的水动力半径。中心线上的正则化奇点线不是无滑移圆柱，其平动与轴向转动的等效半径并不相同，因此 8.5.1 及以前使用的 $\beta=1$ 并没有经过标定。弯曲刚度始终使用物理半径 $r_{\mathrm{pili}}$，根部间隙 `dis` 也由它计算。对 Cortez (2018) 7/2 blob，在直丝中点处的检验结果（$r_{\mathrm{pili}}/L_{\mathrm{pili}}=3\times10^{-3}$，100 段）：

- **轴向扭转。** 正则化 rotlet 线段的单位长度阻力为 $\xi_t=2\pi\mu\varepsilon^2$（数值核对，比值 $1.0000$），而无滑移圆柱为 $4\pi\mu r^2$（Walker 等 2020 式 24；Maxian 等 2022 §II A 3）。两者在 $\varepsilon=\sqrt2\,r$ 时一致；$\beta=1$ 时扭转阻力只有圆柱值的一半。
- **平动。** Rotne–Prager–Yamakawa 线的正则化半径 $\hat a=r\,e^{3/2}/4\approx1.12\,r$ 与半径 $r$ 的圆柱的细长体理论一致（Maxian 等 2021–2023）。对本 blob，取 $\varepsilon=\sqrt{e}\,r\approx1.649\,r$ 时垂直与平行平动迁移率都与之相同（解析与数值）。$\beta=\sqrt2$ 时丝中点的迁移率比该参考高 $2.4\,\%$（垂直）和 $2.9\,\%$（平行），$\beta=1$ 时分别高 $7.9\,\%$ 和 $9.4\,\%$。

因此 $\beta=\sqrt2$ 是默认折中：扭转阻力与圆柱完全一致，平动偏差约 $3\,\%$。要平动精确可取 $\beta=\sqrt e$（此时扭转阻力比圆柱值高 $36\,\%$）；取 $\beta=1$ 可复现 8.5.1 及以前的版本。单根菌毛的纯沉降算例中，$\beta=1$ 与 $\beta=\sqrt2$ 的沉降距离相差 $0.13\,\%$；在马达算例（`D4`）中，整个运行期间胞体的反向转动角相差 $1.95$ 倍（$0.0976\to0.190$ rad）。

### 5.3 Rotlet 线段

纤毛扭转诱导的流动通过正则化 Rotlet 线段表示：

```math
\mathbf{u}(\mathbf{x})
=
-L
\left(
c_a[\mathbf{x}_0]_{\times}
+ c_b[\mathbf{v}]_{\times}
\right)
\boldsymbol{\tau}.
```

其中 `tau` 为线段扭矩密度。实现中将 Cortez (2018) 的全旋度 rotlet 线段核保存为 `M_rotlet`；当该分块与 Stokeslet 平动互易耦合时，通过 `0.5*M_rotlet` 因子转换为物理标准 rotlet 约定。

## 6. 耦合线性系统

每次 ODE 右端项计算都会组装：

```math
\mathbf{L}\dot{\mathbf{z}} = \mathbf{R},
```

其中：

```math
\mathbf{L} = -\mathbf{B}\mathbf{A}\mathbf{Q}.
```

自 v7.6.0 起，轴向扭转不再使用局部 RFT 扭转阻力。平动与轴向自旋未知量通过如下大迁移率耦合：

```math
\begin{bmatrix}
\mathbf{u} \\
\omega_\parallel
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{M}_{\mathrm{weighted}} & \mathbf{C}_{\mathrm{forward}} \\
\mathbf{C}_{\mathrm{reverse}} & \mathbf{M}_{\mathrm{twist}}
\end{bmatrix}
\begin{bmatrix}
\mathbf{f} \\
\mathbf{m}
\end{bmatrix},
```

其中 `M_twist` 为非局部线段积分轴向转动迁移率，且

```math
\mathbf{C}_{\mathrm{forward}}
=
\tfrac{1}{2}\mathbf{M}_{\mathrm{rotlet}}\mathbf{D}_{\mathrm{twist}},
\qquad
\mathbf{C}_{\mathrm{reverse}}
=
\mathbf{W}_m^{-1}\mathbf{C}_{\mathrm{forward}}^{\mathsf T}\mathbf{W}_f.
```

其中 \(\mathbf{W}_f\) 将胞体/RSS 求积权重重复到每个笛卡尔力分量，
\(\mathbf{W}_m=\operatorname{diag}(L_s)\) 使用各线段实际长度。于是
\(\mathbf{W}_f\mathbf{C}_{\mathrm{forward}}
=\mathbf{C}_{\mathrm{reverse}}^{\mathsf T}\mathbf{W}_m\)，
这正是 Lorentz 互易在离散虚功内积下的形式；只有两个测度完全相同时才可直接用普通转置。`1/2` 因子将 Cortez 的全旋度 rotlet 线段核转换为物理标准 rotlet 约定。代码通过 Schur 互补施加逆作用：

```math
\mathbf{S}_{\mathrm{twist}}
=
\mathbf{M}_{\mathrm{twist}}
-
\mathbf{C}_{\mathrm{reverse}}
\mathbf{M}_{\mathrm{weighted}}^{-1}
\mathbf{C}_{\mathrm{forward}},
```

因此主线性系统与背景流修正使用同一个耦合算子。

`M_twist` 将扭矩密度映射为角速度，因此它也必须在线测度下自伴：
\(\mathbf{W}_m\mathbf{M}_{\mathrm{twist}}
=\mathbf{M}_{\mathrm{twist}}^{\mathsf T}\mathbf{W}_m\)。
实现使用每一段端点的实际距离强制该关系，不假设不同纤毛具有相同段长。

主要矩阵块为：

- `Q_global`: 将广义速率映射为胞体和纤毛节点速度。
- `M_weighted`: 全局平动水动力迁移率矩阵，包含胞体 RSM 和纤毛 RSS 分块。
- `M_twist`: 非局部轴向扭转迁移率。
- `D_twist`: 将标量轴向扭矩密度嵌入为沿线段方向的 3D rotlet 向量。
- `C_forward`、`C_reverse`: 在线测度下互为加权伴随、并采用物理标准 rotlet 缩放的平动/扭转耦合块。
- `B_global`: 将分布力和线积分力矩映射到全局力平衡、力矩平衡和纤毛力矩平衡方程。
- `R`: 重力、弹性力矩、主动力矩、背景流修正和马达力矩项。

**全局角动量平衡（自 v7.6.0 起）。** 游动体作用于流体的总水动力矩为 Stokeslet 力矩与线积分后的净 rotlet 力矩之和：

```math
\mathbf{T} = \sum_i \mathbf{r}_i \times \mathbf{f}_i + \sum_s L_s \mathbf{m}_s,
```

因为 `m_s` 是第 `s` 段上的扭矩密度。因此 `B_global` 的力矩行同时包含 `sum(r x f)` 与 `sum(L_s m_s)`。这使受力/力矩自由的游动体严格力矩为零，并使马达驱动或自旋纤毛给出正确的胞体反向旋转。

求解使用 MATLAB decomposition / LU 分解，不显式计算矩阵逆。

## 7. 弹性杆力矩

纤毛采用 Kirchhoff 弹性杆模型。Darboux 向量写作：

```math
\mathbf{U}
=
\kappa_1\mathbf{d}_1
+ \kappa_2\mathbf{d}_2
+ \kappa_3\mathbf{d}_3,
```

director 演化关系为：

```math
\mathbf{d}_{\alpha}' = \mathbf{U}\times\mathbf{d}_{\alpha}.
```

[`compute_pili_elastic_moments.m`](../src/physics/compute_pili_elastic_moments.m) 使用的曲率分量为：

```math
\kappa_1
=
\mathbf{d}_3\cdot\mathbf{d}_2'
- \kappa_{1,0},
```

```math
\kappa_2
=
\mathbf{d}_1\cdot\mathbf{d}_3'
- \kappa_{2,0},
```

```math
\kappa_3
=
\frac{
\mathbf{d}_2\cdot\mathbf{d}_1'
- \kappa_{3,0}
}{1+\sigma}.
```

每根纤毛的内禀曲率可配置为 `frozen`、`zero`、`constant` 或 `helix`。主动力矩可关闭、设为行波形式，或由自定义配置提供。

### 离散曲率差分模板（v8.3.0）

Directors 按分段存储（第 \(i\) 段覆盖 \(s\in[(i-1)\Delta s,\ i\Delta s]\)，\(\Delta s=1/N\)），而力矩平衡写在节点上。共享算子 [`compute_pili_curvature.m`](../src/physics/compute_pili_curvature.m) 在节点 \(i\) 处按

```math
\kappa_1(s_i) = \mathbf{d}_3^{(i)}\cdot\frac{\mathbf{d}_{2,i}-\mathbf{d}_{2,i-1}}{\Delta s},\qquad
\kappa_2(s_i) = \mathbf{d}_1^{(i)}\cdot\frac{\mathbf{d}_{3,i}-\mathbf{d}_{3,i-1}}{\Delta s},\qquad
\kappa_3(s_i) = \mathbf{d}_2^{(i)}\cdot\frac{\mathbf{d}_{1,i}-\mathbf{d}_{1,i-1}}{\Delta s}
```

计算三个分量，其中 \(\mathbf{d}_{\alpha,i}\) 为**分段** directors，\(\mathbf{d}_\alpha^{(i)}\) 为由节点两侧分段坐标架经 SVD 极分解平均得到的**节点** director（[`compute_node_directors.m`](../src/physics/compute_node_directors.m)）。第 \(i-1\)、\(i\) 段的中心位于 \(s_i\mp\Delta s/2\)，因此该差商是 \(\mathbf{d}_\alpha'(s_i)\) 的中心估计，且与其投影所用的节点 director 共位，格式为 \(O(\Delta s^2)\)。8.3.0 之前的实现对节点 directors 本身做后向差分，导数中心落在 \(s_i-\Delta s/2\)，使投影退化为 \(O(\Delta s)\)（在解析的常 Darboux 矢量坐标架上数值验证：观测收敛阶 2.00 对 0.98）。节点 1 借用节点 2 的模板并随后被根部边界条件覆盖；节点 \(N+1\) 无力矩、不参与计算。`init_coupled_params.m` 中的内禀曲率 `kb0` 由同一算子生成，因此 `frozen` 纤毛在 \(t=0\) 处为精确的零力矩平衡（残差 \(\sim 10^{-13}\)）。

## 8. 根部耦合与马达力矩

每根纤毛根部支持两类边界条件：

- `clamped`: 根部弯曲角和扭转角均锁定到胞体坐标系。
- `free_torsion`: 根部弯曲角仍受约束，但根部扭转角 `psi_1` 自由演化。

在 `free_torsion` 模式下，根部 `d3` 力矩平衡行是整根鞭毛关于根部的轴向力矩平衡，其右端为胞体经根部传递的边界扭矩（马达扭矩；无摩擦轴承时为零），见下文“根部扭转行”。自 v7.3.0 起，代码不再向 `R_global(4:6)` 显式加入马达反力矩。原因是这些行表示整个“胞体 + 纤毛”耦合系统的角动量平衡，马达在胞体与纤毛之间交换的力矩属于系统内部力矩，在整体求和中相互抵消。胞体角响应仍然存在，但它由全局 Stokes 耦合求解和根部运动学约束隐式给出。

### 马达扭矩–转速关系（v8.2.0）

有效马达扭矩采用细菌鞭毛马达的线性扭矩–转速关系（Berg & Turner 1993；Berg 2003），以转子相对胞体的转速表述：

$$\tau_{\mathrm{eff}} = \tau_{\mathrm{motor}}\left(1 - \frac{\omega_{\mathrm{rel}}}{\omega_{\max}}\right),\qquad
\omega_{\mathrm{rel}} = \operatorname{sign}(\tau_{\mathrm{motor}})\,\Big[\big(\dot\psi_1 + \cos\theta_1\,\dot\phi_1\big) - \boldsymbol{\Omega}_{\mathrm{body}}\cdot\mathbf{d}_3^{(1)}\Big].$$

$\dot\psi_1 + \cos\theta_1\dot\phi_1$ 是 Walker Euler 角运动学中第 1 段的完整轴向角速度（与线性系统扭转率行所用表达式相同）；减去 $\boldsymbol{\Omega}_{\mathrm{body}}\cdot\mathbf{d}_3^{(1)}$ 后该速度与参考系无关，等于马达实际感受到的转子–定子相对转动。符号因子使转速沿驱动方向度量，因此同一公式适用于 CW 与 CCW 马达。该关系是非对称的：被迫反转（$\omega_{\mathrm{rel}}<0$）时扭矩高于失速值，$\omega_{\mathrm{rel}}>\omega_{\max}$ 时马达变为制动。由于 $\dot{z}$ 关于马达扭矩是仿射的，$\omega_{\mathrm{rel}} = a + B\tau$（$a$ 为无马达解，$B$ 为单位扭矩响应，二者来自 `lin_sys` 的同一次分解），自洽扭矩由小型线性系统 $(I + \operatorname{diag}(|\tau_{\mathrm{motor}}|/\omega_{\max})B)\,\tau = \tau_{\mathrm{motor}} - \operatorname{diag}(|\tau_{\mathrm{motor}}|/\omega_{\max})\,a$ 精确求得。

符号约定：`R_global` 存储的是所施加广义力的负值，因此马达扭矩以 $-\tau$ 进入其 $\psi_1$ 行；正的 `torque_pN_nm` 使转子绕向外根部切向量 $\mathbf{d}_3$ 做右手旋转（从鞭毛末端朝胞体看为逆时针）。

### `free_torsion` 的根部扭转行（v8.4.0）

纤毛力矩平衡的第 $i$ 行是杆段 $[s_i,1]$ 关于节点 $i$ 的平衡：

```math
\mathbf{M}(s_i) + \int_{s_i}^{1}\big(\mathbf{x}-\mathbf{x}_i\big)\times\mathbf{f}\,\mathrm{d}s + \int_{s_i}^{1}\mathbf{m}\,\mathrm{d}s = \mathbf{0},
```

其中 $\mathbf{M}(s_i)$ 是 $[0,s_i)$ 作用于 $[s_i,1]$ 的内力矩。对内部节点它是 Kirchhoff 力矩 $EI(\kappa_1\mathbf{d}_1+\kappa_2\mathbf{d}_2)+GJ\kappa_3\mathbf{d}_3$ 加上积分后的主动力矩。在根部节点，$[0,s_1)$ 就是胞体，故 $\mathbf{M}(0)$ 是胞体施加于杆的扭矩：`clamped` 时为未知反力（三行改为运动学约束）；`free_torsion` 时 $\mathbf{M}(0)\cdot\mathbf{d}_3=\tau_{\mathrm{eff}}$（无摩擦轴承时为零）。节点 1 处可用的离散 $\kappa_3$ 实际是在第 1、2 段之间求值的，它已经是节点 2 平衡所用的内力矩。8.4.0 之前它仍被加入根部行；将根部行与节点 2 行相减可见，马达扭矩当时只由第 1 段的水动力扭矩平衡，因而直杆的稳态转速被放大为物理值的 $N$ 倍（验证：$N=10$ 时修复前后分别为 $8.47\times10^4$ 与 $8.84\times10^3$ rad/s；扭转分布自节点 2 到末端线性递减，$GJ\kappa_3(s_2)=\tau_{\mathrm{eff}}(1-1/N)$）。现在 `compute_pili_elastic_moments.m` 与 `coupled_active_moment.m` 的根部行在两种耦合模式下均为零，`free_torsion` 的 $\psi_1$ 行只携带边界扭矩。由此：无马达的自由轴承要求鞭毛所受净轴向水动力扭矩为零（纯内部扭转波无法使其整体旋转）；有马达时直杆稳态转速为 $\omega=\tau_{\mathrm{eff}}/(\xi_t L)$，其中对 Cortez (2018) 7/2 blob 的正则化 rotlet 线段 $\xi_t=2\pi\mu\epsilon_{\mathrm{pili}}^2$；默认 $\epsilon_{\mathrm{pili}}=\sqrt2\,r_{\mathrm{pili}}$（§5.2）时它等于无滑移圆柱值 $4\pi\mu r_{\mathrm{pili}}^2$，而 $\epsilon_{\mathrm{pili}}=r_{\mathrm{pili}}$（8.5.1 及以前）只有其一半；模拟器复现系数 $2\pi\mu\epsilon^2$ 的误差 $<10^{-3}$。

## 9. 背景流

背景流作为给定稳态速度场进入，并通过水动力修正项写入右端：

```math
\mathbf{R}
\leftarrow
\mathbf{R}
- \mathbf{B}\mathbf{A}^{-1}\mathbf{u}_{\mathrm{bg}}.
```

支持的背景流类型包括：

- `off`: 关闭背景流。
- `uniform`: 用户指定方向的均匀流。
- `shear_z`: 用户指定流向、剪切轴和厚度的有限厚度剪切层。不可压缩性要求归一化流向与剪切方向正交；验证器会拒绝 `abs(dot)>1e-12` 的配置。

在未截断的剪切层内部，
\(\mathbf{u}=U\,\hat{\mathbf d}(\hat{\mathbf s}\cdot\mathbf{x})/H\)，
因此 \(\nabla\cdot\mathbf{u}=(U/H)\hat{\mathbf d}\cdot\hat{\mathbf s}\)。
所以正交条件来自质量守恒，而不是单纯的坐标选择偏好。

## 10. 布朗动力学

启用布朗运动时，观测窗口被划分为长度为 `dt_B` 的 Brownian 区间（最后一个区间使用实际剩余时长）。自 v8.5.0 起，body-only 热运动以**在整个区间内恒定的广义热力**作用于胞体来实现（"力实现"方案），而不再在确定性段之间施加位移 kick：

1. 在每个区间起点，通过探测完整确定性胞体-纤毛系统计算有效胞体迁移率（§10.1）。
2. 抽取冻结迁移率下的 Euler–Maruyama 增量 \(\delta\mathbf{U}\)，并将其换算为在冻结系统上恰好产生该增量的广义力 \(\mathbf{f}_{\mathrm{th}}\)（§10.2）。
3. 将 \(\mathbf{f}_{\mathrm{th}}\) 加到右端项的六个胞体行（`ode_coupled_rhs` 第 8a 步），用 `ode15s` 在该区间内积分耦合系统；纤毛通过与重力、马达扭矩相同的水动力与弹性耦合响应胞体的热运动。

### 10.1 有效胞体迁移率

有效迁移率通过 6 个单位广义胞体载荷探测耦合线性系统得到：

```math
\mathbf{M}_{\mathrm{eff}}^{\ast}(:,i)
=
\left(\mathbf{L}^{-1}(-\mathbf{e}_i)\right)_{1:6},
\qquad i=1,\ldots,6.
```

对称正定检查及舍入尺度的谱修复在无量纲迁移率上完成；若出现具有物理量级的负特征值，则视为装配错误而不是直接截断。修复后的矩阵再通过合同变换转换到物理单位：

```math
\mathbf{D}_{\mathrm{eff}}
=
k_B T\,\mathbf{M}_{\mathrm{eff,dim}},
```

```math
\mathbf{L}_{\mathrm{eff}}
\mathbf{L}_{\mathrm{eff}}^{\mathsf{T}}
=
2\mathbf{D}_{\mathrm{eff}},
```

```math
\delta\mathbf{U}
=
\mathbf{L}_{\mathrm{eff}}
\boldsymbol{\eta}
\sqrt{\Delta t_B},
\qquad
\boldsymbol{\eta}\sim\mathcal{N}(\mathbf{0},\mathbf{I}_6).
```

Brownian 抽样来自本次运行独立的随机流，不受几何缓存创建影响；每个区间恰好消耗 6 个正态随机数。

### 10.2 增量的力实现（v8.5.0）

以区间起点冻结的有效迁移率为准，在 \(\Delta t_B\) 内恰好产生 \(\delta\mathbf{U}\) 的广义力为

```math
\mathbf{f}_{\mathrm{th}}
=
\mathbf{M}_{\mathrm{eff}}^{-1}\,
\frac{\delta\mathbf{U}}{\Delta t_B},
```

在耦合系统的无量纲广义单位中求值（第 1–3 行以 \(F_0\) 计，第 4–6 行以 \(F_0 L_0\) 计；`sample_brownian_body_force`）。它以与重力相同的符号约定进入右端项，\(\mathbf{R}(1{:}6) \leftarrow \mathbf{R}(1{:}6) - \mathbf{f}_{\mathrm{th}}\)，并在 `ode15s` 积分该区间期间保持恒定。因此实际的胞体增量是沿真实耦合轨迹的 \(\int \mathbf{M}_{\mathrm{eff}}(t)\,\mathbf{f}_{\mathrm{th}}\,dt\)，而非冻结迁移率下的值。

**根部不再出现折点的原因。** v8.5.0 之前，胞体姿态在每个节点被瞬时移动，随后第 1 段纤毛被重新投影到移动后的胞体上，而第 \(2\ldots N\) 段保持原取向。因此每次 kick 都在第 1、2 段之间注入一个幅值与 \(|\delta\boldsymbol{\Theta}|/\ell\) 同阶的曲率间断；当 \(\ell \ll 1\) 时这会激发高波数弹性模态，其松弛时间尺度为 \(\ell^4\)，迫使 `ode15s` 采用极小步长。改为力实现后，胞体连续运动，纤毛通过耦合求解跟随；根部约束从未被破坏，也就不会激发伪模态。

**隐式 Ito 漂移。** 迁移率依赖于构型的过阻尼 SDE 在 Ito 约定下需要漂移项 \(k_B T\,\nabla\!\cdot\!\mathbf{M}\)。在一个区间内展开耦合 ODE，

```math
\delta\mathbf{U}_{\mathrm{real}}
=
\mathbf{M}_n\mathbf{f}_{\mathrm{th}}\Delta t_B
+ \tfrac12\,(\partial_j\mathbf{M})\,(\mathbf{M}_n\mathbf{f}_{\mathrm{th}})_j\,\mathbf{f}_{\mathrm{th}}\,\Delta t_B^2
+ \mathcal{O}(\Delta t_B^{3/2}),
```

并利用 \(\mathbb{E}[\mathbf{f}_{\mathrm{th}}\mathbf{f}_{\mathrm{th}}^{\mathsf T}]\,\Delta t_B = 2k_BT\,\mathbf{M}_n^{-1}\)，得到

```math
\mathbb{E}\big[\delta\mathbf{U}_{\mathrm{real}}\big]
=
k_B T\,\nabla\!\cdot\!\mathbf{M}\,\Delta t_B
+ \mathcal{O}(\Delta t_B^{3/2}),
\qquad
\operatorname{Cov}\big[\delta\mathbf{U}_{\mathrm{real}}\big]
=
2k_BT\,\mathbf{M}_n\,\Delta t_B
+ \mathcal{O}(\Delta t_B^{3/2}).
```

沿轨迹采样到的迁移率梯度自动生成了 Ito 漂移：该方案是 body-only Ito SDE 的一致弱一阶积分器，v7.5.0–v8.4.0 使用的显式有限差分 \(k_B T\,\nabla\!\cdot\!\mathbf{M}\) 修正（每个节点额外 12 次耦合求解）不再施加——若再加入将导致漂移被重复计入。

对当前无壁面模型，漂移有闭式表达。在无界流体中，有效迁移率仅通过取向依赖于胞体姿态，\(\mathbf{M}_{\mathrm{eff}} = \mathbf{T}\,\mathbf{M}_{\mathrm{shape}}\,\mathbf{T}^{\mathsf T}\)，\(\mathbf{T}=\mathrm{blkdiag}(\mathbf{R},\mathbf{R})\)，所以平动坐标不进入散度。对绕实验室系轴 \(\mathbf{e}_j\) 的无穷小转动，每个分块按 \(\mathbf{A}\to\mathbf{A}+\varepsilon([\mathbf{e}_j]_\times\mathbf{A}-\mathbf{A}[\mathbf{e}_j]_\times)\) 变化，第 \(i\) 行的漂移为 \(k_BT\sum_j\partial M_{i,R_j}/\partial\theta_j\)。**转动行**为零（\(\mathbf{M}\) 对称时）；**平动行一般不为零**：

```math
\big(k_BT\,\nabla\!\cdot\!\mathbf{M}\big)_{T,i}
= k_BT\sum_{j,n}\epsilon_{ijn}\,(\mathbf{M}_{TR})_{nj},
```

即 \(k_BT\) 乘以平动–转动耦合块反对称部分的轴向量。当 \(\mathbf{M}_{TR}\) 对称时为零，特别是跟踪点处有对称中心的物体（无菌毛的胶囊体、对称布置的菌毛）；当 \(\mathbf{X}_{cm}\) 不是细胞的流体力学中心时不为零：此时胞体转动，跟踪点向该中心松弛。其量级约为 \(2D_R\) 乘以偏移量；对默认胶囊体加一根 3 um 菌毛，约为 0.02 um/s（约为沉降速度的 15%）。该矢量固定在体坐标系中，对取向平均为零，不产生长时间净漂移。力实现格式会自动产生这一项，因为它沿区间积分的是实际的 \(\mathbf{M}(\theta(t))\)；对 \(\mathbf{M}_{TR}\neq0\) 的二维刚体，按代码实现做的蒙特卡洛检验在 1.3 个标准误内复现了 \(k_BT\,\mathbf{R}\mathbf{J}\mathbf{b}\)，\(\mathbf{M}_{TR}=0\) 时为零。因此无需额外添加，如前所述，再加显式的 \(k_BT\,\nabla\!\cdot\!\mathbf{M}\) 会重复计入。区间内纤毛在热力作用下变形导致的 \(\mathbf{M}_{\mathrm{eff}}\) 变化属于耦合动力学本身，而非漂移修正。

**Brownian 近似的适用范围。** 所有纤毛的水动力拖曳都会影响探测得到的 6 x 6 胞体迁移率，但纤毛内部 Euler 角自由度不接受随机驱动，只在确定性弹性、水动力、主动力矩、马达力矩以及由胞体传递而来的热力作用下演化。因此，该模型适合在忽略纤毛内部热涨落时研究胞体质心和取向扩散；它不是完整耦合热平衡过程，也不能再现纤毛内部模态的能量均分。

**已知近似 —— 确定性算子与随机算子的一致性。**

确定性求解使用原始水动力迁移率矩阵 `M_weighted`，当 `eps_body != eps_pili` 时该矩阵在结构上非对称：胞体 RSM 块（源正则化长度 `eps = eps_body`）与纤毛 RSS 块（源正则化长度 `eps = eps_pili`）对同一对胞体-纤毛点对施加了不同的正则化长度，因此 Lorentz 互易定理在非对角的胞体↔纤毛耦合块上不严格成立。相反，布朗协方差路径使用对称化后的有效胞体迁移率 `M_eff = 0.5*(M + M')`，以满足涨落耗散定理对对称扩散张量的要求。

因此该算符分裂依赖"非互易分量足够小"这一近似。实现在每个布朗步测量相对不对称度 `rel_asym = ||M_eff - M_eff^T||_F / ||M_eff||_F`，当 `rel_asym > 1e-3` 时发出警告。典型混合正则化配置下 `rel_asym ~ 1e-5..1e-3`；纯 RSM（无纤毛）系统降至 `<1e-10`（仅 LU 舍入误差）。从完整 `M_weighted` 经 Schur 互补提取到 6×6 胞体块的过程会进一步衰减非对称性，故对平衡分布的实际影响以 `rel_asym` 本身为上界。当 `eps_body ≈ eps_pili` 时该近似成为精确。追求最严格统计保真度的用户应监控 `[BROWN]` 不对称度日志，若 `rel_asym > 1e-2`，可考虑在保证胞体表面分辨率 `N_body` 足以将 RSM 泄漏控制在容差内（见 §5.1）的前提下，将 `eps_body` 调向 `eps_pili`（或反向调整）。

### 10.3 Brownian ODE 重启中的滞后 Jacobian（v8.5.1）

本小节只涉及时间积分器，§10.1–10.2 的模型不变。每个 Brownian 区间都会重启 `ode15s`。对隐式 BDF/NDF 步 \(\mathbf{y}_{n+1}-\gamma h\,\mathbf{f}(\mathbf{y}_{n+1})=\boldsymbol{\psi}\)，求解器用修正牛顿迭代，迭代矩阵为 \(\mathbf{W}=\mathbf{I}-\gamma h\,\tilde{\mathbf{J}}\)，其中 \(\tilde{\mathbf{J}}\) 是 \(\partial\mathbf{f}/\partial\mathbf{y}\) 的任意近似。迭代收敛当且仅当 \(\rho\big(\mathbf{I}-\mathbf{W}^{-1}\mathbf{W}_\text{true}\big)<1\)；收敛后无论 \(\tilde{\mathbf{J}}\) 取何值，求解的都是同一个隐式方程，局部误差估计与步长/阶数选择均基于收敛后的解。因此陈旧的 \(\tilde{\mathbf{J}}\) 只改变收敛速度，不改变离散格式。每次重启重新做稠密有限差分 Jacobian 需要 \(n_\text{ode}\) 次右端项计算，是 Brownian 节点耗时的主要部分（4 菌毛约 70–85%），因此复用上一次的精确 Jacobian（`LaggedJacobian.m`）：当 `ode15s` 在同一次调用中第二次索要 \(\mathbf{J}\)（牛顿失败）、Euler 基旋转之后、以及缓存达到 `cfg.jacobian_max_age` 个节点时刷新。在 4 菌毛验证算例（N_body = 200，5 段/µm，刚性 \(|\lambda|_\max\approx9\times10^{7}\)，单位 \(T_0^{-1}\)）中，实际步长（\(h\le1.7\times10^{-3}\)）下，缓存至多 64 个节点的 Jacobian 对应的收缩矩阵谱半径 \(\le0.17\)；同种子轨迹与新鲜 Jacobian 的解相差约 \(10^{-5}\) 相对量。`cfg.jacobian_reuse=false` 与 8.5.0 逐位一致。

## 11. 流场可视化

[`compute_plot_flowfield.m`](../src/visualization/compute_plot_flowfield.m) 在末状态重构几何，在胞体质心附近的 XZ 和 YZ 平面上评估速度场，遮罩胞体内部区域，并导出 2 x 2 面板图：

- 3D XZ 平面视图。
- 3D YZ 平面视图。
- 2D XZ 速度大小和归一化方向箭头。
- 2D YZ 速度大小和归一化方向箭头。

自定义胞体在可用时使用保存的三角化表面和 alpha-shape 内部遮罩。

## 12. Walker 标度一致性

当耦合问题退化为无胞体的单根自由弹性纤毛时，本项目的无量纲结构设计为与 Walker 风格的 Kirchhoff 弹性杆公式保持一致：

- 空间 $Q$ 矩阵作用一致：将体坐标系或杆坐标系中的广义速率映射为纤毛上的平动和角速度。
- 力平衡行使用梯形求积权重，纤毛力矩平衡行使用每根杆自身的长度标度。
- 弹性力矩和主动力矩均存储在每根纤毛的归一化弧长坐标 $s_k\in[0,1]$ 上。

对于扭转阻力，有量纲线段关系正比于局部段长：

```math
\tau_{\mathrm{twist},j}
\sim
c_t\,\Delta s_k\,\dot{\psi}_j,
\qquad
c_t \sim 4\pi\mu r_{\mathrm{pili}}^2.
```

使用项目全局尺度后，无量纲扭转系数为：

```math
\frac{c_t\,\Delta s_k/T_0}{EI/L_0}
\sim
\frac{\varepsilon^2}{2}\,\Delta \hat{s}_k.
```

代码中每根纤毛使用自身的归一化弧长，因此从全局胞体尺度转换到每根纤毛的 Walker 尺度时，会引入长度比 $\ell_k = L_{\mathrm{pili},k}/L_0$。每根纤毛的无量纲细长度为 $\varepsilon_k = \varepsilon/\ell_k = r_{\mathrm{pili}}/L_{\mathrm{pili},k}$（`build_global_resistance.m`）。

```math
\widetilde{\mathbf{A}}_{\mathrm{twist},k}
=
\frac{\ell_k\varepsilon_k^2}{2}
\widehat{\mathbf{A}}_k
=
\frac{\varepsilon^2}{2\ell_k}
\widehat{\mathbf{A}}_k.
```

两种写法代数等价，因为 $\ell_k\varepsilon_k^2 = \varepsilon^2/\ell_k$。

纤毛力矩平衡行中配套使用的 $\ell_k^3$ 因子会把每根杆尺度下的力矩方程转换回全局尺度，从而保证组装后的系统

```math
\mathbf{L}\dot{\mathbf{z}} = \mathbf{R}
```

在量纲上保持一致。

## 13. 验证范围

代码层面已围绕以下模型不变量进行过检查：

- 四元数约定与胞体姿态运动学。
- RSM 核函数尺度与胞体体积一致性。
- RSS 递推项和线段矩阵对称性。
- Rotlet 反对称性和递推行为。
- Walker 风格 Euler 角 `Q` 矩阵与力矩平衡尺度。
- 固支与自由扭转根部边界条件。
- 多纤毛状态布局、几何重构和矩阵维度。
- 完整耦合有效胞体迁移率与 body-only Brownian Cholesky 正定性检查。
- 平动/扭转加权互易关系与对称功率耗散检查。
- 尺度无关线性求解后向残差、选定状态 `rcond`、有限值及多马达真实不动点残差检查。
- MAT、CSV、图像和日志输出单位。

后向残差超过 \(10^{-10}\) 时告警，超过 \(10^{-6}\) 时中止。`rcond<1e-10` 与 `rcond<1e-14` 仅告警并写入诊断；仅当后向残差同时进入误差带（\(>10^{-6}\)）且 `rcond<1e-14` 时才联合中止。实测极值与固定阈值都会写入 `results.diagnostics`。

详细缺陷修复和版本历史属于 [`../CHANGELOG.md`](../CHANGELOG.md)，不再放在理论文档中。

## 14. 参考文献

- Cortez, R. (2005). The method of regularized Stokeslets.
- Cortez, R. (2018). Regularized Stokeslet segments.
- Walker, B. J. (2020). Elastic filament / bacterial appendage formulation used as the rod-model reference.

## 15. 许可证

MIT License。详见 [`../LICENSE`](../LICENSE)。
