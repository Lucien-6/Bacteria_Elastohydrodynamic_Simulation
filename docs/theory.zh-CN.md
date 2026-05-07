# 理论与数值模型

[English](theory.md) | 中文  
[返回中文 README](../README.zh-CN.md)

> **作者:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **更新日期:** 2026-05-06
> **版本:** 7.4.2

本文档承接原长版 README 中的理论、模型和数值方法说明。仓库首页现在保持精简；实际运行步骤见 [`../USER_GUIDE.md`](../USER_GUIDE.md)，版本历史见 [`../CHANGELOG.md`](../CHANGELOG.md)。

## 1. 物理问题

仿真对象是在 Stokes 流中沉降的细菌胞体及其柔性纤毛 / 鞭毛：

- **胞体:** 刚性胶囊体、近似球体或导入的自定义三维点云。
- **纤毛 / 鞭毛:** 一根或多根不可伸长 Kirchhoff 弹性杆。
- **流体:** 零雷诺数 Stokes 流。
- **耦合:** 胞体与所有纤毛节点之间的水动力相互作用通过全局迁移率求解保留。

默认参数化胞体由两个半球端盖和一个圆柱中段组成。自定义几何模式下，胞体与纤毛必须同时来自用户导入的点云；代码会拒绝“自定义胞体 + 参数化纤毛”或反向混用。

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
| 长度 | $L_0$ | 胞体 `major_axis` |
| 弯曲刚度 | $EI$ | $E\pi r_{\mathrm{pili}}^4/4$ |
| 时间 | $T_0$ | $8\pi\mu L_0^4/EI$ |
| 速度 | $U_0$ | $L_0/T_0$ |
| 力 | $F_0$ | $EI/L_0^2$ |

关键无量纲参数包括：

- 弹性-重力数：$E_g = F_{\mathrm{sink}}L_0^2/EI$。
- 纤毛细长比：$\varepsilon = r_{\mathrm{pili}}/L_0$。
- 泊松比：$\sigma$，由用户给定。

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

每根纤毛由 Euler 角定义材料标架 `(d1, d2, d3)`，其中 `d3` 是纤毛段切向。[`reconstruct_geometry.m`](../reconstruct_geometry.m) 将状态 `Z` 转换为实验室坐标系下的胞体表面点、纤毛节点、director 标架和水动力求解所需的辅助几何量。

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

纤毛正则化长度是物理参数：

```math
\varepsilon_{\mathrm{pili}} = r_{\mathrm{pili}}/L_0.
```

它代表纤毛有限半径，而不是应趋近于零的纯数值参数。

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

其中 `tau` 为线段扭矩密度。对典型细菌纤毛而言，该耦合量级很小，但保留它可以使公式结构更完整。

## 6. 耦合线性系统

每次 ODE 右端项计算都会组装：

```math
\mathbf{L}\dot{\mathbf{z}} = \mathbf{R},
```

其中：

```math
\mathbf{L} = -\mathbf{B}\mathbf{A}\mathbf{Q},
```

```math
\mathbf{A}\mathbf{Q}
=
\begin{bmatrix}
\mathbf{M}_{\mathrm{weighted}}^{-1}\mathbf{Q}_{\mathrm{trans,corr}} \\
\widetilde{\mathbf{A}}_{\mathrm{twist}}\mathbf{Q}_{\mathrm{twist}}
\end{bmatrix},
```

```math
\mathbf{Q}_{\mathrm{trans,corr}}
=
\mathbf{Q}_{\mathrm{trans}}
- \mathbf{M}_{\mathrm{rotlet}}
\widetilde{\mathbf{A}}_{\mathrm{twist}}
\mathbf{Q}_{\mathrm{twist}}.
```

主要矩阵块为：

- `Q_global`: 将广义速率映射为胞体和纤毛节点速度。
- `M_weighted`: 全局水动力迁移率矩阵，包含胞体 RSM 和纤毛 RSS 分块。
- `A_twist`: 扭转阻力算子。
- `B_global`: 将分布力和扭矩映射到全局力平衡、力矩平衡和纤毛力矩平衡方程。
- `R`: 重力、弹性力矩、主动力矩、背景流修正和马达力矩项。

求解使用 MATLAB 反斜杠 / LU 分解，不显式计算矩阵逆。

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

[`compute_pili_elastic_moments.m`](../compute_pili_elastic_moments.m) 使用的曲率分量为：

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

## 8. 根部耦合与马达力矩

每根纤毛根部支持两类边界条件：

- `clamped`: 根部弯曲角和扭转角均锁定到胞体坐标系。
- `free_torsion`: 根部弯曲角仍受约束，但根部扭转角 `psi_1` 自由演化。

在 `free_torsion` 模式下，根部 `d3` 力矩平衡行可以接受弹性扭转、主动 `d3` 力矩和恒定马达力矩。自 v7.3.0 起，代码不再向 `R_global(4:6)` 显式加入马达反力矩。原因是这些行表示整个“胞体 + 纤毛”耦合系统的角动量平衡，马达在胞体与纤毛之间交换的力矩属于系统内部力矩，在整体求和中相互抵消。胞体角响应仍然存在，但它由全局 Stokes 耦合求解和根部运动学约束隐式给出。

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
- `shear_z`: 用户指定流向、剪切轴和厚度的有限厚度剪切层。

## 10. 布朗动力学

启用布朗运动时，确定性弹性流体动力学与随机胞体运动采用算符分裂：

1. 使用 `ode15s` 积分确定性耦合系统。
2. 在每个布朗时间节点计算全系统有效胞体迁移率。
3. 施加 6 维平移-旋转相关布朗扰动。

有效迁移率通过 6 个单位广义胞体载荷探测耦合线性系统得到：

```math
\mathbf{M}_{\mathrm{eff}}^{\ast}(:,i)
=
\left(\mathbf{L}^{-1}(-\mathbf{e}_i)\right)_{1:6},
\qquad i=1,\ldots,6.
```

量纲转换后：

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

可选 `fixman_midpoint` 会在随机中点重新计算迁移率，以近似构型依赖 Ito 漂移。代价是每个布朗步需要约两次迁移率评估。

## 11. 流场可视化

[`compute_plot_flowfield.m`](../compute_plot_flowfield.m) 在末状态重构几何，在胞体质心附近的 XZ 和 YZ 平面上评估速度场，遮罩胞体内部区域，并导出 2 x 2 面板图：

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

代码中每根纤毛使用自身的归一化弧长，因此从全局胞体尺度转换到每根纤毛的 Walker 尺度时，会引入长度比 $\ell_k = L_{\mathrm{pili},k}/L_0$：

```math
\widetilde{\mathbf{A}}_{\mathrm{twist},k}
=
\frac{\ell_k\varepsilon_k^2}{2}
\widehat{\mathbf{A}}_k
=
\frac{\varepsilon^2}{2\ell_k}
\widehat{\mathbf{A}}_k.
```

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
- 全系统布朗有效迁移率和 Cholesky 正定性检查。
- MAT、CSV、图像和日志输出单位。

详细缺陷修复和版本历史属于 [`../CHANGELOG.md`](../CHANGELOG.md)，不再放在理论文档中。

## 14. 参考文献

- Cortez, R. (2005). The method of regularized Stokeslets.
- Cortez, R. (2018). Regularized Stokeslet segments.
- Walker, B. J. (2020). Elastic filament / bacterial appendage formulation used as the rod-model reference.

## 15. 许可证

MIT License。详见 [`../LICENSE`](../LICENSE)。
