# 理论与数值模型

[English](theory.md) | 中文  
[返回中文 README](../README.zh-CN.md)

> **作者:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **更新日期:** 2026-11-24
> **版本:** 10.2.1

本文档说明物理模型和数值方法。实际运行步骤见 [`../USER_GUIDE.md`](../USER_GUIDE.md)，版本历史见 [`../CHANGELOG.md`](../CHANGELOG.md)。

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
\mathrm{dis} = 1.1\,(\varepsilon_{\mathrm{body}}L_0 + r_{\mathrm{pili}}),
```

传入 `pili_building` / `custom_pili_building`，避免纤毛根部正则化核重叠。式中用的是 §5.1 的胞体正则化长度和纤毛的物理半径，因此间隙与 `hydro_radius_factor`（§5.2）无关。

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

其中 \(\mathbf{Q}\) 把广义速率映射为各点速度和各段自旋角速度；\(\mathbf{A}\) 把这些速度映射为作用在流体上的力密度和扭矩密度（即下文大迁移率的逆，由线性求解施加，不显式构造）；\(\mathbf{B}\) 把这些密度汇总进各平衡方程。

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

### 6.1 胞体框架 Block Schur 互补与冻结水动力三色雅可比（v10.0.0）

由于胞体为刚体，其自相互作用 RSM 迁移率块 \(\mathbf{M}_{BB}^{\mathrm{loc}}\)（\(3N_{\mathrm{body}} \times 3N_{\mathrm{body}}\)，对称正定）在胞体局部坐标系中恒定不变。初始化阶段由 `build_body_precomp` 一次性完成 Cholesky 分解 \(\mathbf{M}_{BB}^{\mathrm{loc}} = \mathbf{L}_{BB}\mathbf{L}_{BB}^{\mathsf T}\) 并预求 6 列刚体运动学作用 \(\mathbf{U}_{QB}^{\mathrm{loc}} = (\mathbf{M}_{BB}^{\mathrm{loc}})^{-1}\mathbf{Q}_{B,6}^{\mathrm{loc}}\)。在每次 RHS 评估中，`build_coupled_lin_sys` 将纤毛坐标与切向量转入胞体框架，通过 `build_pili_mobility_blocks` 仅组装随构型变化的分块 \((\mathbf{M}_{PB}, \mathbf{M}_{BP}, \mathbf{M}_{PP}, \mathbf{C}_{\mathrm{twist}}, \mathbf{M}_{\mathrm{twist}})\)，并通过分块 Schur 互补代数消去全部 \(3N_{\mathrm{body}}\) 个胞体未知量：

```math
\mathbf{S}_{PP} = \mathbf{M}_{PP} - \mathbf{M}_{PB}(\mathbf{M}_{BB}^{\mathrm{loc}})^{-1}\mathbf{M}_{BP}.
```

该 reformulation 将每步稠密 LU 分解复杂度从 \(O((3N_{\mathrm{body}} + 3N_{\mathrm{pili}})^3)\) 降至 \(O((3N_{\mathrm{pili}})^3)\)，且与稠密路径严格代数等价（在条件数对应的机器精度内一致，见 `validation/cpc_performance_benchmark.m`）。无纤毛算例（\(M=0\)）则直接旋转预计算的 \(6\times 6\) 阻力矩阵（\(N_{\mathrm{body}}=300\) 时提速约 99×）。

对于刚性 `ode15s` 牛顿迭代，主导的 \(O(N_k^4)\) 刚性来自局部 Kirchhoff 弹性恢复力矩 \(\partial \mathbf{R}/\partial \mathbf{Z}\)，而水动力阻力算子 \(\mathrm{lin\_sys}(\mathbf{Z})\) 仅以 \(O(1)\) 平滑变化。因此 `ode_coupled_jacobian` 仅组装并 LU 分解一次 \(\mathrm{lin\_sys}(\mathbf{Z})\)，对局部模板采用步长 3 图着色（3-colour）求 \(\partial \mathbf{R}/\partial(\theta,\phi,\psi)\)（同时逐列处理根段并包含精确线性化马达扭矩–转速反馈），通过单次多右端回代得到雅可比矩阵，实现 75×–517× 雅可比加速，且不改变 ODE 右端函数和误差控制；如牛顿迭代失败则自动回退有限差分雅可比。布朗分段重启期间，`LaggedJacobian` 在 `jacobian_max_age` 步内复用缓存雅可比，并按需或在换基后通过 `ode_coupled_jacobian` 刷新。

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
- \mathbf{B}\mathbf{A}\,\mathbf{u}_{\mathrm{bg}},
```

其中 \(\mathbf{A}\) 即 §6 的算子：流体载荷取决于相对背景流的速度，\(\mathbf{A}(\mathbf{Q}\dot{\mathbf{z}} - \mathbf{u}_{\mathrm{bg}})\)。

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

**Brownian 近似的适用范围。** 所有纤毛的水动力拖曳都会影响探测得到的 6 x 6 胞体迁移率，但在 body-only 方案中，纤毛内部 Euler 角自由度不接受随机驱动，只在确定性弹性、水动力、主动力矩、马达力矩以及由胞体传递而来的热力作用下演化。因此，body-only 方案适合在忽略纤毛内部热涨落时研究胞体质心和取向扩散；它不是完整耦合热平衡过程，也不能再现纤毛内部模态的能量均分。可选的纤毛慢弯曲模态热驱动见 §10.4（v9.0.0）。

**已知近似 —— 确定性算子与随机算子的一致性。**

确定性求解使用原始水动力迁移率矩阵 `M_weighted`，当 `eps_body != eps_pili` 时该矩阵在结构上非对称：胞体 RSM 块（源正则化长度 `eps = eps_body`）与纤毛 RSS 块（源正则化长度 `eps = eps_pili`）对同一对胞体-纤毛点对施加了不同的正则化长度，因此 Lorentz 互易定理在非对角的胞体↔纤毛耦合块上不严格成立。相反，布朗协方差路径使用对称化后的有效胞体迁移率 `M_eff = 0.5*(M + M')`，以满足涨落耗散定理对对称扩散张量的要求。

因此该算符分裂依赖"非互易分量足够小"这一近似。实现在每个布朗步测量相对不对称度 `rel_asym = ||M_eff - M_eff^T||_F / ||M_eff||_F`，当 `rel_asym > 1e-3` 时发出警告。典型混合正则化配置下 `rel_asym ~ 1e-5..1e-3`；纯 RSM（无纤毛）系统降至 `<1e-10`（仅 LU 舍入误差）。从完整 `M_weighted` 经 Schur 互补提取到 6×6 胞体块的过程会进一步衰减非对称性，故对平衡分布的实际影响以 `rel_asym` 本身为上界。当 `eps_body ≈ eps_pili` 时该近似成为精确。追求最严格统计保真度的用户应监控 `[BROWN]` 不对称度日志，若 `rel_asym > 1e-2`，可考虑在保证胞体表面分辨率 `N_body` 足以将 RSM 泄漏控制在容差内（见 §5.1）的前提下，将 `eps_body` 调向 `eps_pili`（或反向调整）。

### 10.3 Brownian ODE 重启中的滞后 Jacobian（v8.5.1）

本小节只涉及时间积分器，§10.1–10.2 的模型不变。每个 Brownian 区间都会重启 `ode15s`。对隐式 BDF/NDF 步 \(\mathbf{y}_{n+1}-\gamma h\,\mathbf{f}(\mathbf{y}_{n+1})=\boldsymbol{\psi}\)，求解器用修正牛顿迭代，迭代矩阵为 \(\mathbf{W}=\mathbf{I}-\gamma h\,\tilde{\mathbf{J}}\)，其中 \(\tilde{\mathbf{J}}\) 是 \(\partial\mathbf{f}/\partial\mathbf{y}\) 的任意近似。迭代收敛当且仅当 \(\rho\big(\mathbf{I}-\mathbf{W}^{-1}\mathbf{W}_\text{true}\big)<1\)；收敛后无论 \(\tilde{\mathbf{J}}\) 取何值，求解的都是同一个隐式方程，局部误差估计与步长/阶数选择均基于收敛后的解。因此陈旧的 \(\tilde{\mathbf{J}}\) 只改变收敛速度，不改变离散格式。每次重启重新做稠密有限差分 Jacobian 需要 \(n_\text{ode}\) 次右端项计算，是 Brownian 节点耗时的主要部分（4 菌毛约 70–85%），因此复用上一次的精确 Jacobian（`LaggedJacobian.m`）：当 `ode15s` 在同一次调用中第二次索要 \(\mathbf{J}\)（牛顿失败）、Euler 基旋转之后、以及缓存达到 `cfg.jacobian_max_age` 个节点时刷新。在 4 菌毛验证算例（N_body = 200，5 段/µm，刚性 \(|\lambda|_\max\approx9\times10^{7}\)，单位 \(T_0^{-1}\)）中，实际步长（\(h\le1.7\times10^{-3}\)）下，缓存至多 64 个节点的 Jacobian 对应的收缩矩阵谱半径 \(\le0.17\)；同种子轨迹与新鲜 Jacobian 的解相差约 \(10^{-5}\) 相对量。`cfg.jacobian_reuse=false` 与 8.5.0 逐位一致。

### 10.4 纤毛的热涨落（v9.0.0）

设 `brownian.pili_thermal = true`（`cfg.enable_pili_brown`）后，纤毛最慢的弯曲模态也接受热力。默认是每根纤毛的第一、二阶弯曲，每阶两个偏振方向（`brownian.pili_mode_orders = 2`）。扭转和更高阶弯曲仍为确定性。§10.1–10.2 的 body-only 方案不变，仍是默认方案。

**热模态**（`init_pili_brownian_modes`，每次运行在无应力初始状态计算一次）。记 \(\mathbf{x}\) 为所有纤毛内部节点处材料曲率和扭率相对其本征值的偏差（`pili_bending_coords`）。离散 Kirchhoff 杆的弹性能在这组坐标下严格为二次型：

```math
E = \tfrac12\,\mathbf{x}^{\mathsf T}\mathbf{K}\,\mathbf{x},
\qquad
\mathbf{K} = \operatorname{diag}(w_K w_E),
```

其中 \(w_K = (1, 1, 1/(1+\sigma))\)，\(w_E = 1/(N\ell)\)。对每个坐标施加单位广义力探测耦合系统，得到曲率迁移率 \(\dot{\mathbf{x}} = \mathbf{M}_x\mathbf{f}\)；每次探测中胞体和其余节点都不受力，所以 \(\mathbf{M}_x\) 包含纤毛之间以及纤毛与运动胞体之间的水动力耦合。自由松弛为 \(\dot{\mathbf{x}} = -\mathbf{M}_x\mathbf{K}\mathbf{x}\)。取对称部分 \(\mathbf{M}_s\)：

```math
\mathbf{K}^{1/2}\mathbf{M}_s\mathbf{K}^{1/2}\,\mathbf{u}_j = \lambda_j\,\mathbf{u}_j,
\qquad
a_j = \mathbf{u}_j^{\mathsf T}\mathbf{K}^{1/2}\mathbf{x},
\qquad
E = \tfrac12\sum_j a_j^2 .
```

最慢的 \(n = 2\,n_o M\) 个以弯曲为主的模态是热模态（\(n_o\) 为弯曲阶数，\(M\) 为纤毛根数）。它们是整个细胞的模态：多根等长纤毛的同阶弯曲会被水动力耦合混合。平衡时每个模态的能量为 \(k_BT/2\)，即 \(\langle a_j^2\rangle = k_BT/(F_0L_0)\)。模态形状按节点存放在材料坐标系中，随纤毛一起运动，与 Euler 基的选择无关。

本节的数值都对应 `configs/config_example.json` 的细胞：胶囊 \(L_0 = 1.25\) µm、半径 0.4 µm，\(E = 3.279\times10^{8}\) Pa，\(r_{\mathrm{pili}} = 3\) nm，于是 \(EI = 2.09\times10^{-26}\) N m²，\(T_0 = 2.94\) s，293.15 K 时 \(k_BT/(F_0L_0) = 0.2425\)；流体为水（\(\mu = 10^{-3}\) Pa s）。松弛时间与 \(\mu/EI\) 成正比，并大致与纤毛长度的四次方成正比。

| 细胞（3 µm 纤毛，30 段） | 一阶 | 二阶 | 不加噪声的最慢模态 |
| --- | --- | --- | --- |
| 胶囊 + 1 根 | 0.210 s | 13.4 ms | 2.14 ms |
| 胶囊 + 4 根 | 0.075–0.62 s | 9.0–16.9 ms | 2.30 ms |

**联合迁移率**（`compute_effective_coupled_brownian`，每个区间计算一次）。热坐标为 \(\mathbf{y} = (\mathbf{X}_{cm}, \boldsymbol{\theta}_{body}, a_1,\ldots,a_n)\)。在当前状态下用 6 个单位胞体载荷和 \(n\) 个单位模态力探测，得到 \((6+n)\times(6+n)\) 迁移率 \(\dot{\mathbf{y}} = \mathbf{M}_{joint}\,(\mathbf{f}_{body}, \mathbf{g})\)。单位模态力 \(g_j\) 是耦合系统弯曲行上的力矩分布 \(\sqrt{w_K/w_E}\,\mathbf{u}_j\)。矩阵的对称化和正定检查与 §10.1 相同。

每次探测都固定所有不加噪声的纤毛模态，即上述模态分析中的全部扭转模态和全部更快的弯曲模态：对这些模态施加保持力 \(\boldsymbol{\lambda}\)，使它们的幅值不变，\(\mathbf{A}_{hold}^{\mathsf T}\dot{\mathbf{x}} = 0\)。这样做的依据是时间尺度分离。确定性模态的松弛时间从几十微秒（扭转）到约 2 ms（三阶弯曲），热模态的寿命则为 9 ms 到 0.6 s。作用在热坐标上的力如果同时带动了某个快模态，该模态的弹性力矩几乎立刻就把它顶回去；所以在热模态的相关时间内，快模态相当于约束，支配热坐标的是固定快模态时的迁移率（阻力矩阵限制在热坐标上的逆）。让快模态不受力的探测给出的是瞬时迁移率。在直线静止状态，两者在模态块上相同，因为模态此时正是本征向量。纤毛弯曲后两者不同：

- *扭转。* 材料坐标系中的力矩会让弯曲纤毛的各段绕自身轴线自转。自转几乎没有粘性阻力，却使两个曲率分量在材料坐标系中互相转换。胶囊加 1 根纤毛、弯曲幅度为热幅度的五分之一时，不固定扭转得到的模态迁移率是弯曲松弛速率的 14 到 196 倍。
- *快弯曲模态。* 纤毛弯曲后，静止状态的模态形状会与快弯曲模态混合。同一细胞在 293.15 K 的一条热轨迹上，只固定扭转的探测把一阶模态迁移率高估了 16%（中位数），最大达 5 倍；时间积分在 0.2 ms 内的实际响应比该探测低，最多低 50%。固定全部非热模态后，探测与时间积分相差 1–3%。

按偏大的迁移率抽取的噪声，相对于时间积分实际产生的耗散偏弱（涨落–耗散平衡被破坏），热模态达不到能量均分；因此采用完整的固定。于是 \(\mathbf{M}_{joint}\) 的胞体块是“非热模态固定、热模态自由”时的胞体迁移率。它介于 §10.1 的完全柔性迁移率（适用于比快松弛时间更短的时间）和刚性形状迁移率（观测时间超过最慢热松弛时间后达到）之间。

**抽样**（`sample_brownian_coupled_force`）。以区间起点冻结的联合迁移率为准：

```math
\delta\mathbf{Y} = \mathbf{L}\,\boldsymbol{\eta}\sqrt{\Delta t_B},
\qquad
\mathbf{L}\mathbf{L}^{\mathsf T} = 2k_BT\,\mathbf{M}_{joint},
\qquad
\begin{pmatrix}\mathbf{f}_{th}\\ \mathbf{g}\end{pmatrix}
= \mathbf{M}_{joint}^{-1}\,\frac{\delta\mathbf{Y}}{\Delta t_B},
```

每个区间消耗 \(6+n\) 个正态随机数。\(\mathbf{M}_{joint}\) 的交叉块给出胞体与纤毛之间的相关。与 §10.2 相同，力在整个区间内保持恒定，Ito 漂移由沿轨迹采样到的迁移率产生。

**热力的实现**（`pili_thermal_rhs`，`ode_coupled_rhs` 第 8e 步）。模态力不以无约束白噪声力矩的形式施加（那样会把弯曲纤毛拧成螺旋）。实际做法是：每个热力（6 个胞体分量和各模态力）都配有该坐标受约束探测运动中纤毛所受的完整粘性载荷——即胞体坐标系下的中心线力密度 \(\mathbf{f}_c = (\mathbf{A}\mathbf{Q}_{\mathrm{trans}}\,\mathbf{z}_c)_{pili}\)（区间内随胞体共同转动），以及沿各段当前切向 \(\mathbf{d}_3\) 施加的平衡轴向标量力矩密度 \(m_c = (\mathbf{A}\mathbf{Q}_{\mathrm{twist}}\,\mathbf{z}_c)\)，\(c = 1,\ldots,6+n\)。二者合并后在纤毛行上的广义力可在机器精度（\(\sim 10^{-15}\)）下精确复现受约束探测右端项（模态力矩分布与所有非热模态的保持力），即便在大变形 3D 弯曲纤毛上也是如此。保留探测运动自带的平衡轴向力矩密度可抵消弯曲纤毛的粘性搓转阻力而不激发扭转（将弯曲状态下的扭转能再降约两个数量级，并恢复小 \(\Delta t_B\) 下的单调协方差收敛）。该载荷在 6 个胞体行上的合力不重复施加（第 1–6 行置零）；模态部分对整个细胞不产生净力和净力矩。

**不等长多纤毛的模态筛选。** 默认按 `pili_mode_orders` 选取最慢的 \(2\,n_o M\) 个弯曲模态。对包含不等长纤毛的细胞（松弛时间 \(\tau \propto L^4\)），可使用可选的时间尺度门槛 `brownian.pili_tau_cut_s` 选取所有满足 \(\tau_j \ge \tau_{\mathrm{cut}}\) 的慢弯曲模态，并可通过 `brownian.pili_protect_pairs` 在能隙比 \(< 1.15\) 时自动扩充截断边界，避免劈开近简并偏振对。

**保持载荷带来的漂移。** §10.2 中隐式产生 Ito 漂移的前提，是广义力在区间内保持不变。固定在胞体坐标系中的中心线载荷不满足这一点：纤毛运动时，力矩平衡中的力臂和标架随之改变，载荷产生的模态力也跟着变。记 \(g_j(\mathbf{P}_k;\mathbf{y})\) 为区间起点构造的第 \(k\) 个模态的载荷 \(\mathbf{P}_k\) 在构型 \(\mathbf{y}\) 下对模态 \(j\) 的力。把一个区间内的位移展开到二阶，并对热力求平均（\(\langle F_mF_n\rangle = 2k_BT(\mathbf{M}_{joint}^{-1})_{mn}/\Delta t_B\)），得到

```math
\langle\delta\mathbf{Y}\rangle
= k_BT\,\Delta t_B\left(\nabla_{\mathbf{y}}\!\cdot\mathbf{M}_{joint}
  + \mathbf{M}_{joint}\,\mathbf{c}\right),
\qquad
c_j = \sum_{k=1}^{n}\frac{\partial g_j(\mathbf{P}_k;\mathbf{y})}{\partial a_k}.
```

第一项是热坐标的 Ito 漂移。第二项是载荷实现方式带来的假象：一个多余的广义力 \(k_BT\,\mathbf{c}\)。对位于胞体轴线上的单根纤毛，它在直线静止状态为零，随变形线性增大；293.15 K 时其效果相当于二阶弯曲的刚度变化约 \(+4\%\)，一阶约 \(-3\%\)。对不在轴线上的纤毛，它在静止状态就不为零：探测运动中胞体在动，载荷带有沿纤毛方向的分量，纤毛一偏转，这部分载荷的力矩就有一阶变化。胶囊加 4 根纤毛时，静止状态下对称的一阶模态 \(c = -0.31\)；不补偿的话，293.15 K 时该模态的均值会偏移其热幅度（均方根）的 0.15。胞体坐标没有贡献，因为载荷随胞体转动，且其合力不施加。抽样时在该区间的模态力上加 \(-k_BT\,\mathbf{c}\)，把这一项消去。\(\mathbf{c}\) 在区间起点沿各热模态的运动学方向（使用对称化并正定修复后的 \(\mathbf{M}_{joint}\) 计算 \(\mathbf{z}_{probe}\mathbf{M}_{joint}^{-1}\)，非热模态固定，胞体不动）用中心差分计算；在 \(2n\) 个位移后的状态上无需重建水动力迁移率或完整稠密矩阵，而是通过 `compute_effective_coupled_brownian` 内部与 `build_force_balance` 代数等价的 \(O(N_k)\) 无矩阵子程序 `apply_pili_B_force_rows` 直接计算纤毛力矩平衡行。加上这一项后，热模态的平稳分布是幅值 \(a_j\) 上的 \(\exp(-E/k_BT)\)，即能量均分。

**Brownian 区间。** 本方案是建立在时间尺度分离基础上的粗粒化有效模型，因此 \(\Delta t_B\) 需落在双向合理窗口内：

- *上限*（\(\Delta t_B \lesssim \tau_{\mathrm{slow}}/4\)）：力在 \(\Delta t_B\) 内保持恒定，会使松弛速率为 \(\lambda\) 的模态平稳方差乘以 \((2/x)\tanh(x/2)\)，\(x = \lambda\Delta t_B\)（\(x = 0.25\) 时低 0.5%，\(x = 0.5\) 时低 2%）。\(\Delta t_B\) 应不超过最快热模态松弛时间的四分之一，对上表的 3 µm 纤毛即 0.002 s（`run_simulation` 在 \(x > 0.5\) 时发出警告，\(x > 2.0\) 时默认拦截报错）。
- *下限*（\(\Delta t_B \gtrsim \tau_{\mathrm{det}}\)）：若固定 `pili_mode_orders` 而将 \(\Delta t_B\) 缩小至远低于最慢非热化弯曲模态的松弛时间 \(\tau_{\mathrm{det}}\)（上表为 2.14–2.30 ms），单步内非热化弯曲模态来不及绝热松弛（即便平衡轴向力矩密度已保证单步协方差单调收敛且无扭转泄露，仍会残留少量来自 \(\mathbf{M}_{yb}\) 的漂移误差）；因此本方案不能通过固定模态数而盲目令 \(\Delta t_B \to 0\) 来检验收敛性（当 \(\Delta t_B / \tau_{\mathrm{det}} < 0.5\) 时 `run_simulation` 会发出警告）。如需更细时间分辨率，应同步提高 `pili_mode_orders` 或降低 `pili_tau_cut_s`。

**胞体扩散。** 下表是线性化模型的胞体扩散系数（本方案的精确协方差递推，\(\Delta t_B = 0.002\) s）。\(\parallel\) 和 \(\perp\) 相对胞体长轴；\(\boldsymbol{\mu}_{flex}\) 是 §10.1 的完全柔性胞体迁移率，\(\boldsymbol{\mu}_{rigid}\) 是纤毛形状冻结时细胞的迁移率。

| | 1 根：\(D_{t\perp}\) | \(D_{r\perp}\) | 4 根：\(D_{t\parallel}\) | \(D_{t\perp}\) | \(D_{r\parallel}\) | \(D_{r\perp}\) |
|---|---|---|---|---|---|---|
| 时间间隔 \(\Delta t_B\)，\(D / k_BT\mu_{flex}\)：只对胞体加噪声 | 0.97 | 0.87 | 0.91 | 0.95 | 0.64 | 0.88 |
| 时间间隔 \(\Delta t_B\)，\(D / k_BT\mu_{flex}\)：一、二阶 | 0.98 | 0.89 | 0.93 | 0.96 | 0.66 | 0.90 |
| 长时间，\(D / k_BT\mu_{rigid}\)：只对胞体加噪声 | 0.98 | 0.41 | 0.55 | 0.78 | 0.05 | 0.26 |
| 长时间，\(D / k_BT\mu_{rigid}\)：一、二阶 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |

纤毛模态热化后，热坐标构成一个满足涨落–耗散关系的封闭系统：噪声按系统实际响应的迁移率抽取，噪声和胞体运动都不向确定性模态输入能量。这时长时间扩散等于 \(k_BT\,\boldsymbol{\mu}_{rigid}\)，即形状涨落被平均掉之后的值，与热化的阶数无关。body-only 方案中纤毛没有噪声却仍在耗散，所以长时间扩散低于该值；4 根纤毛的细胞绕长轴的自转低到二十分之一。在一个 Brownian 区间的时间间隔上，两种方案相差不超过 3%。两者都低于完全柔性的值，因为比区间更快的模态此时已经松弛：body-only 方案中靠弹性松弛，热涨落方案中则按定义处理（探测时固定）；后者 \(\mathbf{M}_{joint}\) 的胞体块是完全柔性迁移率的 0.90（1 根，\(D_{r\perp}\)）和 0.70（4 根，\(D_{r\parallel}\)）。从 \(k_BT\,\boldsymbol{\mu}_{flex}\) 过渡到该值发生在这些松弛时间之内（最长约 2 ms），本方案不解析这一过程。

**适用范围。** 本方案针对被动、本征为直线、根部固支的纤毛。`root_coupling = free_torsion` 会被拒绝：纤毛绕根部轴线的自转是一个没有弹性回复力矩的坐标，既不加噪声也不被固定，纤毛一旦弯曲，联合迁移率就失去对称性和正定性。本征曲率覆盖和主动力矩可以使用但会给出警告：热模态取自被动的初始状态，下面的能量均分检验不覆盖这些情形。

**验收（v9.0.0）。** 在生产分辨率下（306 个胞体点，每根纤毛 30 段）对 1 根和 4 根纤毛做了算子检验：联合迁移率不对称度 \(2.7\times10^{-4}\) 和 \(1.4\times10^{-3}\)；模态对角元与松弛速率相差不超过 0.2%；右端项对热力的响应与 \(\mathbf{M}_{joint}\mathbf{F}\) 相差不超过 \(3.7\times10^{-3}\)；恒定模态力下的静态变形与该力相差不超过 \(2.6\times10^{-3}\)；由 ODE 测得的松弛速率与本征值相差不超过 0.2%；更换 Euler 基后的差别小于 \(10^{-7}\)。在一条热轨迹的 6 个状态上（1 根纤毛，293.15 K），时间积分的响应是探测迁移率的 0.97–1.00。对静止状态的 4 根纤毛细胞，用正反两个方向的力各积分一次，测出各模态在一个区间内的平均位移：它与 \(k_BT\,\Delta t\,(\nabla\!\cdot\mathbf{M}_{joint} + \mathbf{M}_{joint}\mathbf{c})\) 相差不超过 2%，其中保持载荷这一项单独复现的系数为 0.993。线性协方差递推给出每个热模态的 \(\langle a_j^2\rangle/k_BT\) 在 0.996 到 1.000 之间，body-only 方案为 0–0.51。

通过 `run_simulation` 运行的蒙特卡洛系综（116 个胞体点，\(\Delta t_B = 0.002\) s，无重力）在平稳时间窗内给出的 \(\langle a_j^2\rangle/k_BT\) 如下，格式为“模拟值 / 同一时间窗的线性理论值”：

| 系综 | 运行 | 一阶 | 二阶 |
|---|---|---|---|
| 1 根，293.15 K | 48 × 1.5 s | 1.044 ± 0.061 / 0.998 | 0.989 ± 0.016 / 0.998 |
| 1 根，2.93 K | 24 × 1.2 s | 1.039 ± 0.106 / 0.997 | 0.992 ± 0.027 / 0.998 |
| 4 根，293.15 K | 24 × 1.0 s | 0.938 ± 0.050 / 0.957 | 0.976 ± 0.016 / 0.998 |
| 4 根，2.93 K | 24 × 1.0 s | 0.903 ± 0.063 / 0.957 | 1.012 ± 0.014 / 0.998 |
| 1 根，293.15 K，只对胞体加噪声 | 48 × 1.5 s | 0.403 ± 0.024 / 0.456 | 0.169 ± 0.002 / 0.174 |
| 4 根，293.15 K，只对胞体加噪声 | 24 × 1.0 s | 0.267 ± 0.023 / 0.231 | 0.106 ± 0.002 / 0.104 |

加纤毛噪声的 8 个值与理论值的偏差都不超过 1.4 倍标准误差。（4 根纤毛的运行时长是最慢松弛时间 0.62 s 的 1.6 倍，所以一阶的理论值是 0.957 而不是 1。）2.93 K 时变形只有 293.15 K 时的十分之一，胞体扩散在所有时间间隔上都与线性理论一致（24 个值的 \(\chi^2\) 为 28.9 和 20.3）。293.15 K 时，凡是细胞形状起作用的地方就有偏差：该纤毛的持久长度 \(EI/k_BT = 5.2\) µm，与其长度 3 µm 相当，热变形并不小，而弯曲的纤毛会改变细胞的迁移率。4 根纤毛时，时间间隔为 \(\Delta t_B\) 的绕长轴自转扩散是线性值的 0.937 ± 0.016；\(\mathbf{M}_{joint}\) 的胞体块对平衡形状取平均给出 0.951 ± 0.008，2.93 K 的系综给出 0.990 ± 0.014。详见 `CHANGELOG.md`（9.0.0）。

## 11. 流场可视化

[`compute_plot_flowfield.m`](../src/visualization/compute_plot_flowfield.m) 在末状态重构几何，在过胞体质心的 XZ（Y = Y_cm）、YZ（X = X_cm）和 XY（Z = Z_cm）平面上评估速度场，遮罩胞体内部区域，并导出 2 x 3 面板图：

- 上排：3D XZ、YZ、XY 平面视图。
- 下排：2D XZ、YZ、XY 速度大小和归一化方向箭头。

自定义胞体在可用时使用保存的三角化表面和 alpha-shape 内部遮罩。

自 v8.7.0 起，每个已保存输出点还记录 `results.flow_state`：完整状态、菌毛局部基，以及产生该输出的 Brownian 热力。[`render_flowfield_from_results.m`](../post_processing/render_flowfield_from_results.m) 用这些数据绘制末状态、指定时刻或动画。速度核与仿真结束时的末状态图相同。v8.7.1 将胞体内部留白，坐标范围盖住整个细菌，并用光滑的表面截线作轮廓。v8.7.3 增加 XY 平面，色标放在子图 f 右侧，99 分位色标计入三个平面。`opts.drop_first_frame`（默认 false）会在 `times`、`animation` 和 `all` 重绘中去掉已保存输出 1。v8.7.0 之前的文件没有这套历史数据。

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

## 14. 实测步态主动力矩

[`../chlamydomonas/`](../chlamydomonas/) 中的双鞭毛游动模型不使用预先给定的行波力矩。实测平面步态作为已知的状态和速率代入。刚体速度来自与正向求解器相同的耦合线性系统的力、力矩行（§6）。杆的各行则给出使该运动成为离散模型精确解的内主动力矩

\[
\mathbf{M}_{\mathrm{act}} = \boldsymbol{\kappa}_{\mathrm{dev}} - \mathrm{lin\_sys}\,\dot{\mathbf{z}}.
\]

只保留绕 \(\mathbf{d}_2\) 的弯曲分量。`chlamy_inverse_moment` 把它存成每个节点上的时间傅里叶级数。闭环算例使用 `chlamy_fit_moment` 的二次谐波解析拟合；鞭毛 2 是鞭毛 1 的镜像。拍动平面倾斜把力矩绕每根鞭毛自己的根部切线转动，两根鞭毛的幅值、相位和频率比可以不同；可选的近端 1ω+2ω 面外弯矩包络（`chlamy_make_perp_envelope`）沿 \(\mathbf{d}_3 \times \hat{\mathbf{n}}\) 施加且严格零寄生扭矩。沉降和马达驱动算例不走这条路径。

## 15. 参考文献

- Cortez, R., Fauci, L. & Medovikov, A. (2005). The method of regularized Stokeslets in three dimensions: analysis, validation, and application to helical swimming. *Physics of Fluids* 17, 031504. 胞体 RSM 核函数。
- Cortez, R. (2018). Regularized Stokeslet segments. *Journal of Computational Physics* 375, 783–796. 纤毛 RSS 与 rotlet 线段。
- Walker, B. J., Ishimoto, K. & Gaffney, E. A. (2020). Efficient simulation of filament elastohydrodynamics in three dimensions. *Physical Review Fluids* 5, 123103. Euler 角弹性杆列式、自适应基、圆柱的扭转阻力。
- Maxian, O., Mogilner, A. & Donev, A. (2021). Integral-based spectral method for inextensible slender fibers in Stokes flow. *Physical Review Fluids* 6, 014102.
- Maxian, O., Sprinkle, B., Peskin, C. S. & Donev, A. (2022). Hydrodynamics of a twisting, bending, inextensible fiber in Stokes flow. *Physical Review Fluids* 7, 074101.
- Maxian, O., Sprinkle, B. & Donev, A. (2023). Bending fluctuations in semiflexible, inextensible, slender filaments in Stokes flow: towards a spectral discretization. *Journal of Chemical Physics* 158, 154114. 以上三篇 Maxian 等的论文给出纤毛水动力半径的细长体参考值（§5.2）。
- Berg, H. C. & Turner, L. (1993). Torque generated by the flagellar motor of *Escherichia coli*. *Biophysical Journal* 65, 2201–2216.
- Berg, H. C. (2003). The rotary motor of bacterial flagella. *Annual Review of Biochemistry* 72, 19–54. 马达的线性扭矩–转速关系（§8）。

## 16. 许可证

MIT License。详见 [`../LICENSE`](../LICENSE)。
