# 细菌胞体-纤毛耦合弹性流体动力学仿真

> **作者**: Lucien (lucien-6@qq.com)  
> **创建日期**: 2026-02-06  
> **更新日期**: 2026-02-12  
> **版本**: 1.6

---

## 目录

1. [项目概述](#1-项目概述)
2. [快速开始](#2-快速开始)
3. [架构与数据流](#3-架构与数据流)
4. [物理问题与建模策略](#4-物理问题与建模策略)
5. [无量纲化体系](#5-无量纲化体系)
6. [状态向量设计](#6-状态向量设计)（含 6.4 初始姿态参数 `alpha_0`）
7. [几何重构](#7-几何重构)
8. [流体力学：正则化 Stokeslet 方法](#8-流体力学正则化-stokeslet-方法)
9. [运动学矩阵 Q](#9-运动学矩阵-q)
10. [力/力矩平衡矩阵 B](#10-力力矩平衡矩阵-b)
11. [弹性恢复力矩（Kirchhoff 弹性杆）](#11-弹性恢复力矩kirchhoff-弹性杆)
12. [耦合线性系统的组装与求解](#12-耦合线性系统的组装与求解)
13. [固接边界条件](#13-固接边界条件)
14. [自适应基底旋转（欧拉角奇异性处理）](#14-自适应基底旋转欧拉角奇异性处理)
15. [时间积分方案](#15-时间积分方案)
16. [输出系统：日志、数据导出与动画](#16-输出系统日志数据导出与动画)
17. [文件结构](#17-文件结构)
18. [缩放一致性证明：与 Walker 的等价关系](#18-缩放一致性证明与-walker-的等价关系)
19. [多纤毛（Multi-Pili）支持](#19-多纤毛multi-pili支持)
20. [审查与验证汇总](#20-审查与验证汇总)
21. [参考文献](#21-参考文献)
22. [附录 A：关键缩放因子推导](#附录-a关键缩放因子推导)
    - [A.4 扭转阻力的额外 ℓ_k 因子](#a4-为什么扭转阻力需要额外的-ell_k-因子)
    - [A.5 逐行量纲一致性验证](#a5-线性系统-ldotzr-逐行量纲一致性验证)
23. [许可证](#许可证)

---

## 1. 项目概述

本项目模拟一个带有弹性纤毛（pili）的胶囊形细菌在 Stokes 流中因重力沉降的运动。胞体被视为刚性胶囊，每根纤毛建模为不可伸长的 Kirchhoff 弹性杆。所有水动力学相互作用（胞体-胞体、胞体-纤毛、纤毛-纤毛）均通过统一的正则化 Stokeslet 方法（RSM）求解。

**核心特性：**

- **刚性胞体**：胶囊形（半球帽 + 圆柱中段），四元数参数化姿态（无万向节死锁）
- **弹性纤毛**：Kirchhoff 弹性杆 + 欧拉角参数化（Walker 2020）
- **水动力学**：正则化 Stokeslet 方法（RSM），自适应逐点正则化参数
- **时间积分**：MATLAB `ode15s`（隐式刚性求解器）
- **自适应基底**：逐纤毛欧拉角基底旋转以回避奇异性
- **固接边界条件**：双重保障（运动学率级约束 + 几何投影修正）
- **输出系统**：双通道日志、`.mat` 数据持久化、静态图、3D MP4 动画
- **可视化**：参数化光滑曲面胞体（半透明，Gouraud 着色）+ 逐纤毛独立配色折线 + 时序渐变色
- **初始姿态**：用户可配置 `alpha_0`（XZ 平面内胞体长轴与 +X 的夹角），范围 `[-π/2, π/2]`
- **案例管理**：用户自定义案例名称 `case_name`，自动命名输出目录
- **计时系统**：全局/分段计时、ETA 预估、专业进度条、人类可读时间格式化

---

## 2. 快速开始

```matlab
cd Coupled_Body_Pili
main_coupled_simulation   % 以默认参数运行
```

仿真输出自动保存至 `output/<case_name>/` 或 `output/run_YYYYMMDD_HHMMSS/`（取决于 `case_name` 是否设置）：

```
output/<case_name>/                 % 用户自定义案例名称
output/run_YYYYMMDD_HHMMSS/         % 默认（case_name 为空时）
├── simulation_results.mat    % results + params 结构体
├── summary_plot.png          % 4 子图概览（300 DPI）
├── bacteria_3d_motion.mp4    % 3D 运动动画
└── simulation.log            % 完整运行日志
```

**用户可配置参数**（在 `main_coupled_simulation.m` 中设置）：

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `case_name` | `''`（空） | 案例名称（非空时用作输出目录名） |
| `alpha_0` | `pi/2` | 初始姿态角 [rad]：XZ 平面内长轴与 +X 的夹角 |
| `N_body` | 1000 | 胞体表面力作用点数 |
| `seg_per_um` | 150 | 纤毛每微米离散段数 |
| `T_obs` | 180.0 | 观测时间 [s] |
| `animation_skip` | 1 | 动画帧跳跃间隔 |

---

## 3. 架构与数据流

### 3.1 模块调用关系图

```mermaid
graph TD
    A[main_coupled_simulation.m] --> B[init_coupled_params.m]
    A --> C[init_coupled_state.m]
    A --> D[adaptive_basis_check.m]
    A --> CG[enforce_clamped_geometry.m]
    A --> E[ode15s]
    E --> F[ode_coupled_rhs.m]
    F --> G[reconstruct_geometry.m]
    F --> H[build_global_resistance.m]
    F --> I[compute_pili_elastic_moments.m]
    G --> J[unpack_state.m]
    G --> K[pili_directors.m]
    G --> L[pili_spatial_coords.m]
    G --> M[quat_to_rotmat.m]
    D --> N[pili_euler_angles.m]
    CG --> N
    F --> O[clamped_bc_rates.m]
    F --> P[pack_state_rates.m]
    F --> Q[quat_kinematics.m]
    Q --> R[quat_multiply.m]
    F --> S[skew_matrix.m]
    B --> T[capsule_body_building.m]
    B --> U[pili_building.m]
    U --> V[pili_type_building.m]
    U --> W[compute_curvature_binormal.m]
```

### 3.2 分层模块设计

| 层次 | 文件 | 职责 |
|------|------|------|
| 主脚本层 | `main_coupled_simulation.m` | 参数设置、积分循环、后处理、输出 |
| 初始化层 | `init_coupled_params.m`, `init_coupled_state.m` | 特征尺度、几何构建、初始状态 |
| 核心计算层 | `ode_coupled_rhs.m` | Q、A、B 矩阵组装、线性系统求解 |
| 几何工具层 | `reconstruct_geometry.m`, `pili_*.m`, `quat_*.m` | 状态 ↔ 空间坐标转换 |
| 流体力学层 | `build_global_resistance.m` | RSM 迁移率/阻力矩阵 |
| 弹性力学层 | `compute_pili_elastic_moments.m` | Kirchhoff 杆曲率 |
| 边界条件层 | `clamped_bc_rates.m`, `enforce_clamped_geometry.m` | 固接运动学约束 + 几何投影 |
| 数值稳定性层 | `adaptive_basis_check.m` | 欧拉角奇异性回避 |
| 几何生成层 | `capsule_body_building.m`, `pili_building.m`, `pili_type_building.m`, `compute_curvature_binormal.m` | 胞体/纤毛离散化 |
| 进度显示层 | `odetpbar.m` | ODE 段内专业进度条（ode15s OutputFcn） |

---

## 4. 物理问题与建模策略

### 4.1 问题描述

模拟一个带有弹性纤毛的胶囊形细菌在 Stokes 流中因重力沉降的运动：

- **细菌胞体**：刚性胶囊形（两个半球帽 + 圆柱中段），表面离散化为 $N_{\text{body}}$ 个力作用点（默认 1000）
- **弹性纤毛**：$M$ 根不可伸长的 Kirchhoff 弹性杆，各含 $N_k$ 个线段（$N_k + 1$ 个节点）
- **流体**：Stokes 流（零 Reynolds 数），所有水动力学相互作用通过 RSM 处理

### 4.2 方法来源

| 组件 | 来源 | 方法 |
|------|------|------|
| 胞体形状、布点、纤毛附着位 | 父目录 `Stokes Brown Flexible` | `capsule_body_building`, `pili_building` |
| 纤毛弹性杆模型 | `3D-Flexible-Filaments` (Walker 2020) | Kirchhoff 杆 + 欧拉角参数化 |
| 水动力学（平动） | Walker (2020) | 正则化 Stokeslet 方法（RSM） |
| 水动力学（扭转） | Walker (2020) | 阻力力理论（RFT） |
| 胞体姿态 | 本方案新增 | 四元数参数化（避免万向节死锁） |
| 时间积分 | Walker (2020) | MATLAB `ode15s` 隐式求解器 |

### 4.3 核心耦合机制

1. **运动学耦合**：纤毛根部固接于胞体表面，根节点速度 = 胞体刚体运动速度
2. **水动力学耦合**：全局阻力矩阵包含胞体-胞体、胞体-纤毛、纤毛-纤毛的所有相互作用
3. **弹性耦合**：纤毛的弹性恢复力矩通过水动力学反作用力影响胞体运动

---

## 5. 无量纲化体系

### 5.1 特征尺度

所有量均以胞体半长轴 $a$（`major_axis`）为基准进行无量纲化：

| 特征量 | 符号 | 定义 | 物理含义 | 已验证 |
|--------|------|------|----------|--------|
| 长度尺度 | $L_0$ | $a$（半长轴） | 胞体几何特征长度 | ✓ |
| 时间尺度 | $T_0$ | $\displaystyle\frac{8\pi\mu L_0^4}{EI}$ | 弹性弛豫时间 | ✓ |
| 速度尺度 | $U_0$ | $\displaystyle\frac{L_0}{T_0} = \frac{EI}{8\pi\mu L_0^3}$ | 弹性速度 | ✓ |
| 力尺度 | $F_0$ | $\displaystyle\frac{EI}{L_0^2}$ | 弯曲力 | ✓ |
| 力线密度尺度 | $F_{\ell}$ | $\displaystyle\frac{EI}{L_0^3}$ | 弯曲力密度 | ✓ |

其中 $\mu$ 为流体动力粘度，$EI$ 为纤毛弯曲刚度（$EI = E \cdot \pi r_{\text{pili}}^4 / 4$）。

有量纲 Stokeslet 方程中的 $1/(8\pi\mu)$ 前因子被无量纲化**完全吸收**，无量纲形式为 $\hat{\mathbf{u}} = \hat{S}\hat{\mathbf{f}}$。✓

### 5.2 关键无量纲参数

| 参数 | 符号 | 定义 | 含义 | 已验证 |
|------|------|------|------|--------|
| 弹性-重力数 | $Eg$ | $\displaystyle\frac{F_{\text{sink}} L_0^2}{EI}$ | 重力与弯曲刚度之比 | ✓ |
| 全局细长比 | $\varepsilon$ | $\displaystyle\frac{r_{\text{pili}}}{L_0}$ | 纤毛半径/胞体尺度 | ✓ |
| 泊松比 | $\sigma$ | — | 纤毛材料泊松比 | ✓ |
| 纤毛无量纲长度 | $\ell_k$ | $\displaystyle\frac{L_{\text{pili},k}}{L_0}$ | 第 $k$ 根纤毛的无量纲长度 | ✓ |

沉降力：$F_{\text{sink}} = \Delta\rho \cdot V_{\text{capsule}} \cdot g$

胶囊（spherocylinder）精确体积：

$$
V_{\text{capsule}} = \pi b^2 \left(\frac{4}{3}b + 2(a - b)\right)
$$

其中 $b$ 为半短轴（球冠半径），$a$ 为半长轴，$2(a-b)$ 为中间圆柱段长度。

### 5.3 无量纲坐标

$$
\mathbf{x}^* = \frac{\mathbf{x}}{L_0}, \quad t^* = \frac{t}{T_0}, \quad \mathbf{v}^* = \frac{\mathbf{v}}{U_0}, \quad \mathbf{f}^* = \frac{\mathbf{f}}{F_\ell}
$$

以下省略星号 ${}^*$，所有方程均在无量纲坐标下表述。

### 5.4 有量纲输出单位

后处理阶段将无量纲结果映射回物理单位：

| 量 | 单位 | 转换公式 | 已验证 |
|------|------|------|--------|
| 时间 | s（秒） | $t = t^* \cdot T_0$ | ✓ |
| 长度 | μm（微米） | $x = x^* \cdot L_0 \times 10^6$ | ✓ |
| 力 | pN（皮牛） | $F = F^* \cdot F_0 \times 10^{12}$ | ✓ |
| 速度 | μm/s | $v = v^* \cdot U_0 \times 10^6$ | ✓ |
| 角度/角速度 | rad, rad/s | 无量纲化不改变 | ✓ |

---

## 6. 状态向量设计

### 6.1 全局状态向量

$$
\mathbf{Z} = \begin{pmatrix} \mathbf{X}_{\text{cm}} \\ \mathbf{q} \\ \boldsymbol{\theta}_1 \\ \boldsymbol{\phi}_1 \\ \boldsymbol{\psi}_1 \\ \vdots \\ \boldsymbol{\theta}_M \\ \boldsymbol{\phi}_M \\ \boldsymbol{\psi}_M \end{pmatrix} \in \mathbb{R}^{n_{\text{state}}}
$$

- $\mathbf{X}_{\text{cm}} \in \mathbb{R}^3$：胞体质心位置
- $\mathbf{q} = (q_0, q_1, q_2, q_3)^T \in \mathbb{R}^4$：胞体姿态四元数（$|\mathbf{q}| = 1$）
- $\boldsymbol{\theta}_k, \boldsymbol{\phi}_k, \boldsymbol{\psi}_k \in \mathbb{R}^{N_k}$：第 $k$ 根纤毛各段的欧拉角

总自由度：$n_{\text{state}} = 7 + 3\sum_{k=1}^{M} N_k$

### 6.2 广义速率向量

线性系统的核心未知量（速度/角速度空间，而非 $\dot{\mathbf{Z}}$）：

$$
\dot{\mathbf{z}} = \begin{pmatrix} \mathbf{V} \\ \boldsymbol{\Omega} \\ \dot{\boldsymbol{\theta}}_1 \\ \dot{\boldsymbol{\phi}}_1 \\ \dot{\boldsymbol{\psi}}_1 \\ \vdots \end{pmatrix} \in \mathbb{R}^{n_{\text{dof}}}
$$

自由度总数（不含四元数冗余维度）：$n_{\text{dof}} = 6 + 3\sum_{k=1}^{M} N_k$

### 6.3 四元数运动学

$$
\dot{\mathbf{q}} = \frac{1}{2} \begin{pmatrix} 0 \\ \boldsymbol{\Omega} \end{pmatrix} \otimes \mathbf{q}
$$

左乘约定，$\boldsymbol{\Omega}$ 在实验室坐标系中。投影归一化维持 $|\mathbf{q}| = 1$：

$$
\dot{\mathbf{q}} \leftarrow \dot{\mathbf{q}} - (\dot{\mathbf{q}} \cdot \mathbf{q})\,\mathbf{q}
$$

四元数到旋转矩阵的转换（$R$：体坐标系 → 实验室坐标系）：

$$
R(\mathbf{q}) = \begin{pmatrix} 1-2(q_2^2+q_3^2) & 2(q_1 q_2 - q_0 q_3) & 2(q_1 q_3 + q_0 q_2) \\ 2(q_1 q_2 + q_0 q_3) & 1-2(q_1^2+q_3^2) & 2(q_2 q_3 - q_0 q_1) \\ 2(q_1 q_3 - q_0 q_2) & 2(q_2 q_3 + q_0 q_1) & 1-2(q_1^2+q_2^2) \end{pmatrix}
$$

### 6.4 初始姿态参数 `alpha_0`

用户通过 `alpha_0`（单位：rad）控制胞体的初始姿态——XZ 平面内胞体长轴与 +X 轴正方向的夹角，取值范围 $[-\pi/2, \pi/2]$。

胞体长轴在体坐标系中沿 Z 轴。为使其在实验室坐标系中指向 $(\cos\alpha_0, 0, \sin\alpha_0)$，需绕 Y 轴旋转 $\theta_Y = \pi/2 - \alpha_0$，对应初始四元数：

$$
\mathbf{q}_0 = \left(\cos\frac{\theta_Y}{2}, \; 0, \; \sin\frac{\theta_Y}{2}, \; 0\right)^T
$$

| `alpha_0` | 含义 | 初始四元数 |
|-----------|------|------------|
| $\pi/2$（默认） | 竖直向上（+Z） | $(1, 0, 0, 0)^T$ |
| $0$ | 水平（+X） | $(\cos\frac{\pi}{4}, 0, \sin\frac{\pi}{4}, 0)^T$ |
| $-\pi/2$ | 竖直向下（−Z） | $(0, 0, 1, 0)^T$ |

初始化时（`init_coupled_state.m`），纤毛的初始几何坐标通过 $R(\mathbf{q}_0)$ 旋转至匹配胞体姿态的实验室坐标系后，再计算欧拉角。`init_coupled_params.m` 中包含范围校验。

---

## 7. 几何重构

每次 ODE 右端函数调用时，从状态向量 $\mathbf{Z}$ 重构所有空间坐标。

### 7.1 胞体表面点

$$
\mathbf{x}_i^{\text{lab}} = R(\mathbf{q})\,\mathbf{x}_i^{\text{local}} + \mathbf{X}_{\text{cm}}
$$

### 7.2 纤毛节点坐标

**根节点位置**：$\mathbf{x}_{\text{root},k} = \mathbf{X}_{\text{cm}} + R(\mathbf{q})\,\mathbf{r}_k^{\text{local}}$

**欧拉角到指向向量 directors**（自适应基底 $R_k$ 中的 ZYZ 约定）：

$$
\mathbf{d}_3^{(j)} = \begin{pmatrix} \sin\theta_j \cos\phi_j \\ \sin\theta_j \sin\phi_j \\ \cos\theta_j \end{pmatrix}, \quad \mathbf{d}_1^{(j)} = \begin{pmatrix} -\sin\phi_j\cos\psi_j - \cos\theta_j\cos\phi_j\sin\psi_j \\ \cos\phi_j\cos\psi_j - \cos\theta_j\sin\phi_j\sin\psi_j \\ \sin\theta_j\sin\psi_j \end{pmatrix}
$$

$$
\mathbf{d}_2^{(j)} = \begin{pmatrix} \sin\phi_j\sin\psi_j - \cos\theta_j\cos\phi_j\cos\psi_j \\ -\cos\phi_j\sin\psi_j - \cos\theta_j\sin\phi_j\cos\psi_j \\ \sin\theta_j\cos\psi_j \end{pmatrix}
$$

**节点坐标**（实验室坐标系）：$\mathbf{x}_{j+1}^{(k)} = \mathbf{x}_j^{(k)} + \frac{\ell_k}{N_k}\, R_k^T \,\mathbf{d}_3^{(j)}$

---

## 8. 流体力学：正则化 Stokeslet 方法

### 8.1 RSM 核函数（Cortez et al. 2005）

$$
S_{ij}^{\epsilon} = \frac{(|\mathbf{r}|^2 + 2\epsilon_{ij}^2)\,\mathbf{I} + \mathbf{r}\mathbf{r}^T}{(|\mathbf{r}|^2 + \epsilon_{ij}^2)^{3/2}}
$$

其中 $\mathbf{r} = \mathbf{x}_i - \mathbf{x}_j$，$\epsilon_{ij}$ 为自适应逐对正则化参数。

### 8.2 加权迁移率矩阵

$$
M_{\text{weighted}} = S^{\epsilon} \cdot \text{diag}(\mathbf{w}), \quad w_j = \begin{cases} 1 & \text{胞体点（点力）} \\ \Delta s_k = \ell_k / N_k & \text{纤毛 } k \text{ 节点（力线密度）} \end{cases}
$$

### 8.3 平动阻力矩阵

$\mathbf{A}_{\text{trans}} = M_{\text{weighted}}^{-1}$ — 实际计算通过 LU 分解直接求解 `M_weighted \ Q_trans`，避免显式 $O(N^3)$ 矩阵求逆。

### 8.4 扭转阻力矩阵（RFT）

纤毛的扭转使用阻力力理论，对每根纤毛 $k$：

$$
\tilde{A}_k = \frac{\ell_k\,\varepsilon_k^2}{2}\,\hat{A}_k, \quad [\hat{A}_k]_{3(j-1)+1:3j, \, j} = \mathbf{d}_3^{(j)}
$$

其中 $\varepsilon_k = \varepsilon / \ell_k = r_{\text{pili}} / L_{\text{pili},k}$ 为第 $k$ 根纤毛的局部细长比（纤毛物理半径与其自身长度之比）。$\hat{A}_k$ 将标量扭转速率映射为 3D 转矩密度向量，$\mathbf{d}_3^{(j)}$ 使用**实验室坐标系**中的切线方向。

> **额外因子 $\ell_k$ 的来源**：$\mathbf{B}_k$ 的平动力列涉及 $\Delta s \cdot x$（两个长度因子 $\sim \ell_k^2$），而扭转力矩列仅涉及 $\Delta s$（一个长度因子 $\sim \ell_k$）。统一的 $\ell_k^3$ 缩放对平动列正确，但扭转列需要额外补偿一个 $\ell_k$，通过在 $\tilde{A}_k$ 中引入该因子实现。详见 [附录 A.4](#a4-为什么扭转阻力需要额外的-ell_k-因子)。

### 8.5 自适应正则化参数

> **RSM 力点重合原则**：在正则化 Stokeslet 方法中，力的施加点与速度的评估点天然重合（force point ≡ evaluation point），因此几何生成模块仅输出一组表面/节点坐标，无需独立的力点数据。

$$
\epsilon_j = \begin{cases} h_{\text{body}}^\text{min} / 2 & \text{胞体表面点} \\ r_{\text{pili}} / L_0 & \text{纤毛节点} \end{cases}, \quad \epsilon_{ij} = \epsilon_j \text{（力点决定规则）}
$$

纤毛根部间距：$d_{\text{gap}} = 1.1 \times (\epsilon_{\text{body}} \cdot L_0 + r_{\text{pili}})$

---

## 9. 运动学矩阵 Q

$\mathbf{Q}_{\text{global}}$ 将广义速率 $\dot{\mathbf{z}}$ 映射为所有节点的速度（平动 + 扭转）：

$$
\mathbf{v}_{\text{all}} = \mathbf{Q}_{\text{global}}\,\dot{\mathbf{z}}, \quad \mathbf{Q}_{\text{global}} \in \mathbb{R}^{n_{\text{vel}} \times n_{\text{dof}}}
$$

其中 $n_{\text{vel}} = 3 N_{\text{total}} + \sum_k N_k$。

### 9.1 胞体点速度

$$
\mathbf{v}_i = \mathbf{V} - [\mathbf{r}_i]_\times\,\boldsymbol{\Omega}
$$

其中 $\mathbf{r}_i = \mathbf{x}_i - \mathbf{X}_{\text{cm}}$，$[\cdot]_\times$ 为反对称矩阵。

### 9.2 纤毛节点速度

> **关键缩放因子**：欧拉角速率列使用 $\ell_k / N_k$（而非 Walker 原始代码的 $1/N_k$），确保纤毛速度单位为 $L_0/T_0$，与胞体速度一致。

**扭转速率**：$\dot{\psi}_j^{\text{twist}} = \dot{\psi}_j + \cos\theta_j\,\dot{\phi}_j$

**基底 → 实验室坐标变换**：$\mathbf{Q}_k^{\text{lab}} = \text{blkdiag}(\text{kron}(\mathbf{I}_{N_k+1}, R_k^T), \mathbf{I}_{N_k}) \cdot \mathbf{Q}_k^{\text{ext}}$

**根部耦合**：$\dot{\mathbf{x}}_{\text{root}} = \mathbf{V} + \boldsymbol{\Omega} \times \mathbf{r}_{\text{root}}$

---

## 10. 力/力矩平衡矩阵 B

$\mathbf{B}_{\text{global}} \in \mathbb{R}^{n_{\text{dof}} \times n_{\text{force}}}$ 将力/力矩密度映射为平衡方程。

### 10.1 全局力平衡（第 1–3 行）

采用梯形积分：端点权重 $\Delta s_k / 2$，内部权重 $\Delta s_k$。

> 这与 Walker 的 $1/(2N^2)$ 权重不同。Walker 的权重在其系统中可行（力平衡行 RHS=0），但本模型含非零重力 RHS（$Eg\hat{\mathbf{g}}$），必须使用物理正确的梯形权重。

### 10.2 全局力矩平衡（第 4–6 行）

$$
\sum_{\text{all}} w_i (\mathbf{x}_i - \mathbf{X}_{\text{cm}}) \times \mathbf{f}_i = \mathbf{0}
$$

### 10.3 纤毛弯矩平衡（第 7 行起）

每根纤毛 $3N_k$ 个方程。纤毛节点坐标归一化至单位长度后构造 $\mathbf{B}_k$，使其与 Walker 的 $\mathbf{B}_W$ 结构完全一致。

**$\ell_k^3$ 缩放**：$\mathbf{B}_{k}^{\text{lab}} = \ell_k^3 \cdot \mathbf{B}_W \cdot \mathbf{R}_f$ — 严格推导见 [§18](#18-缩放一致性证明与-walker-的等价关系)。

> **平动力列与扭转力矩列的缩放差异**：$\ell_k^3$ 统一应用于 $\mathbf{B}_k$ 的所有列，对平动力列正确（涉及 $\Delta s \cdot x$ 两个长度因子）。扭转力矩列仅涉及 $\Delta s$（一个长度因子），需要额外一个 $\ell_k$ 才能与 $\ell_k^4 = EH_k$ 一致——该补偿因子已移入 $\tilde{A}_k$（见 [§8.4](#84-扭转阻力矩阵rft)）。

---

## 11. 弹性恢复力矩（Kirchhoff 弹性杆）

三分量曲率偏差（相对于固有形状）：

$$
\kappa_1^{(j)} = \mathbf{d}_3^{NC}(j,:) \cdot \frac{d\mathbf{d}_2^{NC}}{ds}\bigg|_j - \kappa_{1,0}^{(j)}
$$

$$
\kappa_2^{(j)} = \mathbf{d}_1^{NC}(j,:) \cdot \frac{d\mathbf{d}_3^{NC}}{ds}\bigg|_j - \kappa_{2,0}^{(j)}
$$

$$
\kappa_3^{(j)} = \frac{1}{1+\sigma}\left(\mathbf{d}_2^{NC}(j,:) \cdot \frac{d\mathbf{d}_1^{NC}}{ds}\bigg|_j - \kappa_{3,0}^{(j)}\right)
$$

$\kappa_3$ 的 $1/(1+\sigma)$ 因子反映了扭转刚度 $GJ = EI/(1+\sigma)$ 的差异。导数中乘以 $N_k$ 因子对应归一化坐标下的有限差分（$\Delta s = 1/N_k$）。

---

## 12. 耦合线性系统的组装与求解

### 12.1 线性系统

$$
\mathbf{L}\,\dot{\mathbf{z}} = \mathbf{R}, \quad \mathbf{L} = -\mathbf{B}_{\text{global}} \cdot \mathbf{A}_{\text{ext}} \cdot \mathbf{Q}_{\text{global}}
$$

实际计算采用分块直接求解策略（不显式构造 $\mathbf{A}_{\text{ext}}$）：

```matlab
AQ_trans = M_weighted \ Q_trans;     % 平动：LU 直接求解
AQ_twist = A_twist_total * Q_twist;  % 扭转：显式矩阵乘
AQ       = [AQ_trans; AQ_twist];
lin_sys  = -B_global * AQ;
```

**性能**：LU 分解 $O(N^3/3)$ + 回代 $O(N^2 \cdot n_{\text{col}})$，相比全矩阵求逆的 $O(N^3)$ 提速约 **2–3 倍**。

### 12.2 右端向量

$$
\mathbf{R} = \begin{pmatrix} -Eg\,\hat{\mathbf{g}} \\ \mathbf{0}_{3\times 1} \\ \kappa_1^{(1)},\, \kappa_2^{(1)},\, \kappa_3^{(1)},\, \ldots \end{pmatrix}
$$

| 行 | 方程 | RHS |
|----|------|-----|
| 1–3 | 力平衡 | $-Eg\,\hat{\mathbf{g}}$（RSM 符号约定） |
| 4–6 | 力矩平衡 | $\mathbf{0}$ |
| 7+ | 纤毛弯矩平衡 | $\kappa_\alpha^{(j)}$（弹性恢复力矩） |

---

## 13. 固接边界条件

双重保障：**率级约束**（线性系统行替换）+ **几何投影**（`enforce_clamped_geometry`）。

### 13.1 首段欧拉角速率

将胞体角速度变换到纤毛基底：$\boldsymbol{\omega}_{\text{basis}} = R_k \boldsymbol{\Omega}_{\text{lab}}$

$$
\dot{\theta}_1 = -\omega_x \sin\phi_1 + \omega_y \cos\phi_1
$$

$$
\dot{\phi}_1 = \omega_z - \frac{\cos\theta_1}{\sin\theta_1}(\omega_x \cos\phi_1 + \omega_y \sin\phi_1)
$$

$$
\dot{\psi}_1 = \frac{\omega_x \cos\phi_1 + \omega_y \sin\phi_1}{\sin\theta_1}
$$

Jacobian：$\mathbf{J}_{\text{euler}} = \mathbf{J}_{\text{basis}} \cdot R_k$，其中

$$
\mathbf{J}_{\text{basis}} = \begin{pmatrix} -\sin\phi_1 & \cos\phi_1 & 0 \\ -\frac{\cos\theta_1\cos\phi_1}{\sin\theta_1} & -\frac{\cos\theta_1\sin\phi_1}{\sin\theta_1} & 1 \\ \frac{\cos\phi_1}{\sin\theta_1} & \frac{\sin\phi_1}{\sin\theta_1} & 0 \end{pmatrix}
$$

在线性系统中替换对应行，设置约束 $\dot{\theta}_1 - \mathbf{J}_{\text{euler}}(1,:) \cdot \boldsymbol{\Omega} = 0$，RHS = 0。

### 13.2 几何投影修正

每个积分区间结束后，`enforce_clamped_geometry` 将首段欧拉角投影至满足 $D_{\text{lab}}^{(1)} = R_{\text{body}} \cdot D_{\text{local}}^{(k)}$ 的精确值，校正数值漂移。

---

## 14. 自适应基底旋转（欧拉角奇异性处理）

### 14.1 检测与旋转

当任意纤毛的 $\theta_j$ 接近 $0$ 或 $\pi$（阈值 $\delta = \pi/20$）时：

1. 在 $10^4$ 个球面采样点中，寻找与所有 $\pm\mathbf{d}_3$ 方向距离最大的最优新 z 轴方向 $\mathbf{v}_{\text{opt}}$
2. 构造旋转矩阵：$R_{\text{local}} = \mathbf{I} + 2\mathbf{K}^2$，其中 $\mathbf{K} = [\mathbf{v}_{\text{opt}}]_\times$（等价于 $R_{\text{local}} = 2\mathbf{v}\mathbf{v}^T - \mathbf{I}$，即绕 $\mathbf{v}$ 旋转 $\pi$ 的 Rodrigues 公式）
3. 更新基底：$R_k^{\text{new}} = R_{\text{local}} \cdot R_k^{\text{old}}$

### 14.2 固有曲率 $\kappa_0$ 的基底不变性

正交变换保持点积不变：$(R\mathbf{a}) \cdot (R\mathbf{b}) = \mathbf{a} \cdot \mathbf{b}$。因此 $\kappa_0$ 在基底变换下**精确不变**，无需重新计算。✓

### 14.3 ODE 事件函数

使用 `ode15s` 的事件机制，在 $\theta$ 接近奇异性（$< \delta/2$）时终止当前积分段，执行基底旋转后重新启动。

---

## 15. 时间积分方案

### 15.1 求解器

MATLAB `ode15s`（隐式多步法），适用于**刚性**弹性流体动力学系统。

### 15.2 主循环流程

```
while T_achieved < T_obs_dimless
    1. 自适应基底检查（在 ode15s 外部）
    2. 设置事件函数（奇异性检测 + 超时）
    3. 调用 ode15s 积分到终止时间或事件触发
    4. 归一化四元数
    5. 执行几何投影修正（enforce_clamped_geometry）
    6. 在输出时间点插值保存结果
end
```

> 自适应基底检查在 `ode15s` **外部**执行，避免积分器内部状态的不连续性。

### 15.3 ODE 右端函数评估流程

```
ode_coupled_rhs(t, Z, params):
  1.  从 Z 解包 → X_cm, q, {θ,φ,ψ}_k
  2.  重构几何 → X_all, directors
  3.  构建加权迁移率矩阵 M_weighted 和扭转阻力 A_twist
  4.  构建运动学矩阵 Q_global
  5.  分块求解 AQ = [M_weighted\Q_trans; A_twist*Q_twist]
  6.  构建力平衡矩阵 B_global
  7.  计算弹性力矩 κ（每根纤毛）
  8.  组装线性系统 L = -B·AQ, R = [-Eg·ĝ; 0; κ]
  9.  施加固接边界条件（替换首段方程）
  10. 求解 L·ż = R
  11. 从 ż 提取 V, Ω → 计算 dq/dt
  12. 打包 dZ = [V; dq/dt; dθ; dφ; dψ; ...]
```

---

## 16. 输出系统：日志、数据导出与动画

### 16.1 输出目录与案例管理

用户可通过 `case_name` 参数自定义案例名称。当 `case_name` 非空时，输出目录以案例名称命名；留空时回退至时间戳命名：

```
output/<case_name>/                 % case_name = 'test_vertical'
output/run_YYYYMMDD_HHMMSS/         % case_name = ''（默认）
├── simulation_results.mat    % results + params 结构体（-v7.3 格式）
├── summary_plot.png          % 4 子图概览（300 DPI，Arial 字体）
├── bacteria_3d_motion.mp4    % 3D 运动动画（MPEG-4，95% 品质）
└── simulation.log            % 完整双通道日志
```

`case_name` 定义在脚本顶部（`CASE IDENTIFICATION` 区域），位于所有参数设置之前，确保日志和输出路径的一致性。日志中同时记录 `Run ID`（案例名称或自动生成的 `run_` 前缀时间戳）和 `Timestamp`（精确执行时刻）。

### 16.2 双通道日志系统

所有日志消息同时输出到控制台和 `simulation.log` 文件，通过 `dual_log` 局部函数实现。结构化标签：

| 标签 | 含义 |
|------|------|
| `[PARAMS]` | 参数加载 |
| `[INIT]` | 初始化过程 |
| `[ODE]` | ODE 积分进度 |
| `[POST]` | 后处理 |
| `[SAVE]` | 数据保存 |
| `[PLOT]` | 绘图生成 |
| `[ANIM]` | 动画渲染 |

### 16.3 数据持久化

`results` 结构体保存的字段：

| 字段 | 说明 | 单位 |
|------|------|------|
| `t` | 物理时间 | s |
| `t_dimless` | 无量纲时间 | — |
| `X_cm`, `X_cm_um` | 质心轨迹 | 无量纲 / μm |
| `q` | 四元数历史 | — |
| `X_pili`, `X_pili_um` | 纤毛节点坐标 | 无量纲 / μm |
| `V_sink` | 沉降速度 | μm/s |
| `F_sink_pN` | 沉降力 | pN |
| `units` | 单位元数据结构体 | — |

### 16.4 3D 动画与可视化

使用 MATLAB `VideoWriter` 生成 MPEG-4 动画（95% 品质）。

**胞体渲染**：由参数化球柱体（spherocylinder）曲面网格组成，使用 `surf` 绘制：

- 网格结构：左半球（pole → equator）+ 圆柱段 + 右半球（equator → pole），三段拼接为连续曲面
- 分辨率：圆周 36 点 × 经向 34 点 = 1,224 顶点
- 面颜色：森林绿 `[0.20, 0.65, 0.33]`，透明度 `FaceAlpha = 0.35`
- Gouraud 光照着色（双光源：`headlight` + `left`），`material('dull')`
- 每帧通过 `R(q) * capsule_flat + X_cm` 向量化变换，无需重新生成网格

**纤毛渲染**：连续折线 `plot3`，`LineWidth = 2.5`，根部（附着点）圆点标记。每根纤毛分配独立颜色（从 10 色可区分色板循环取色），支持多纤毛场景下的视觉区分：

| 索引 | 颜色 | RGB |
|------|------|-----|
| 1 | Medium blue | `[0.20, 0.40, 0.80]` |
| 2 | Burnt orange | `[0.85, 0.33, 0.10]` |
| 3 | Olive green | `[0.47, 0.67, 0.19]` |
| 4 | Crimson | `[0.64, 0.08, 0.18]` |
| 5 | Purple | `[0.49, 0.18, 0.56]` |
| 6 | Goldenrod | `[0.93, 0.69, 0.13]` |
| 7 | Sky blue | `[0.30, 0.75, 0.93]` |
| 8 | Grey | `[0.60, 0.60, 0.60]` |
| 9 | Dark blue | `[0.00, 0.45, 0.74]` |
| 10 | Rose pink | `[0.85, 0.47, 0.66]` |

静态快照图中，每根纤毛的时序渐变色基于其独立基色，从浅色（早期时刻）渐变至深色（晚期时刻）。

**质心轨迹**：浅红色 `[1.00, 0.55, 0.55]` 折线，从起始时刻到当前帧的累积轨迹，清晰展示沉降路径。

**期刊级排版**：Arial 字体，μ 为正体（TeX `{\mu}`），坐标轴标签 12pt，标题 13pt 加粗，刻度朝外。

**其他特性**：

- `animation_skip`：用户自行设置的帧跳跃间隔（每隔 N 个输出点取一帧）
- **等比例立方体包围盒**：遍历全部帧计算全局坐标范围，构建等纵横比的固定坐标轴范围，确保画面完全不跳动
- OpenGL 渲染器确保透明度正确
- 离屏渲染（`Visible='off'`），不干扰交互

### 16.5 计时与 ETA 预估系统

全面的计时系统覆盖仿真全周期，所有计时信息同步输出至控制台和日志文件：

| 阶段 | 计时器 | ETA | 输出示例 |
|------|--------|-----|----------|
| 初始化 | `t_init_elapsed` | — | `[INIT] Initialisation completed in 2.3 s` |
| ODE 段内 | `odetpbar` OutputFcn | ✓ 线性外推 | `[ODE]   [==========..........] 50.0%  Elapsed: 23 s  ETA: 23 s` |
| ODE 全局 | `tic_total`，逐段 `seg_wall` | ✓ 线性外推 | `[ODE]   Progress: 35.2% ǀ Elapsed: 1 min 23.5 s ǀ ETA: 2 min 34.8 s` |
| 动画渲染 | `t_anim_start` | ✓ 线性外推 | `[ANIM] Frame 40/201 (20%) ǀ Elapsed: 8.2 s ǀ ETA: 32.8 s` |
| 全局 | `tic_sim_start` | — | 最终摘要中的 `Total wall time` |

**ODE 段内进度条**（`odetpbar.m`）：自包含的专业进度条，作为 `ode15s` 的 `OutputFcn` 运行。使用 `=` / `.` 字符绘制 30 格宽度进度条，实时显示百分比、已用时间和 ETA。输出节流（≥ 0.25 s 间隔更新）避免过多 I/O。采用 `\b` 退格机制原地刷新，仅输出至控制台（不写入日志文件），段结束后打印最终 100% 行并换行。

**段间进度**：每个积分段完成后，通过双通道日志 `logf` 输出全局进度百分比、累计已用时间和剩余 ETA 估计。

**最终摘要**包含分阶段耗时明细（初始化 / ODE 积分 / 动画渲染 / 总计），以及 ODE 积分段的统计信息（mean / min / max）。

时间格式化使用 `format_duration()` 辅助函数，根据量级自动选择 `s` / `min s` / `h min s` 格式。

---

## 17. 文件结构

```
Coupled_Body_Pili/
│
├── main_coupled_simulation.m       % 主仿真脚本（入口）
│
├── init_coupled_params.m           % 参数初始化（无量纲化、几何构建）
├── init_coupled_state.m            % 初始状态向量 Z0
│
├── ode_coupled_rhs.m               % 核心 ODE 右端函数（Q/A/B 组装 + 线性求解）
├── reconstruct_geometry.m          % 从状态向量重构所有空间坐标
├── build_global_resistance.m       % RSM 迁移率矩阵 + 扭转 RFT
├── compute_pili_elastic_moments.m  % 弹性曲率偏差计算
│
├── enforce_clamped_geometry.m      % 固接几何投影（约束漂移校正）
├── clamped_bc_rates.m              % 固接 BC 首段欧拉角速率
├── adaptive_basis_check.m          % 逐纤毛自适应基底旋转管理
│
├── capsule_body_building.m         % 胶囊胞体表面网格生成（RSM：仅位置点）
├── pili_building.m                 % 纤毛结构初始化（弯曲刚度、曲率等）
├── pili_type_building.m            % 纤毛形态生成（5 种类型，RSM：仅位置点）
├── compute_curvature_binormal.m    % 离散弹性杆曲率二法向计算
│
├── quat_to_rotmat.m                % 四元数 → 旋转矩阵
├── quat_multiply.m                 % Hamilton 四元数乘法
├── quat_kinematics.m               % dq/dt = ½[0;Ω]⊗q
├── skew_matrix.m                   % 3×3 反对称矩阵 [v]×
│
├── pili_directors.m                % 欧拉角 → directors (d1, d2, d3)
├── pili_euler_angles.m             % directors → 欧拉角（反函数）
├── pili_spatial_coords.m           % 根位置 + 欧拉角 → 节点坐标
├── unpack_state.m                  % 从 Z 提取各分量
├── pack_state_rates.m              % 将各分量打包为 dZ
│
├── odetpbar.m                      % ODE 段内专业进度条（ode15s OutputFcn）
├── textprogressbar.m               % 文本进度条工具（备用）
│
├── README.md                       % 本文件
├── coupled_scheme_documentation.md % 理论方案文档
├── review_report.md                % 全面审查报告
│
└── output/                         % 仿真输出目录（自动创建）
    └── <case_name>/ 或 run_YYYYMMDD_HHMMSS/
        ├── simulation_results.mat
        ├── summary_plot.png
        ├── bacteria_3d_motion.mp4
        └── simulation.log
```

**命名规范**：所有 MATLAB 函数文件统一使用 `snake_case`。头注释采用标准格式：`%FUNC_NAME  Brief description`，包含 `Description:`, `Syntax:`, `Inputs:`, `Outputs:`, `Author:`, `Created:`, `Updated:` 等字段。

---

## 18. 缩放一致性证明：与 Walker 的等价关系

本节严格证明耦合系统的方程与 Walker et al. (2020) 在数学上完全等价。

### 18.1 Stokeslet 坐标缩放性质

$$
\hat{S}^{\ell_k \hat{\epsilon}_W}(\ell_k \hat{\mathbf{r}}) = \frac{1}{\ell_k} \hat{S}^{\hat{\epsilon}_W}(\hat{\mathbf{r}})
$$

### 18.2 A 矩阵等价

本项目的加权迁移率（纤毛 $k$ 部分）：

$$
M_{\text{pili}} = \hat{S}_{L_0} \cdot \Delta s_k = \frac{1}{\ell_k} \hat{S}_W \cdot \frac{\ell_k}{N_k} = \frac{\hat{S}_W}{N_k}
$$

$$
A_{\text{pili}} = M_{\text{pili}}^{-1} = N_k \hat{S}_W^{-1} = A_W \quad \checkmark
$$

### 18.3 $EH = \ell_k^4$ 等价性

Walker 的弹性流体动力学数：$EH = 8\pi\mu L_{\text{pili}}^4 / (EI \cdot T_{\text{ref}})$

以本项目的时间尺度 $T_0 = 8\pi\mu L_0^4 / EI$ 代入：

$$
EH_k = \frac{L_{\text{pili},k}^4}{L_0^4} = \ell_k^4
$$

### 18.4 B·A·Q 乘积缩放验证（平动部分）

| 矩阵 | 本项目定义 | 与 Walker 的比值 |
|-------|-----------|-----------------|
| $\mathbf{Q}_{\text{trans}}$ | $(\ell_k/N_k) \cdot Q_{\text{base}}$ | $\times \ell_k$ |
| $\mathbf{A}_{\text{trans}}$ | $A_W$（已证） | $\times 1$ |
| $\mathbf{B}_{k,\text{trans}}$ | $\ell_k^3 \cdot B_{W,\text{trans}}$ | $\times \ell_k^3$ |

乘积：

$$
\mathbf{B}_{k,\text{trans}} \cdot \mathbf{A}_{\text{trans}} \cdot \mathbf{Q}_{\text{trans}} = \ell_k^3 B_W \cdot A_W \cdot \ell_k Q_W = \ell_k^4 \cdot B_W A_W Q_W = EH_k \cdot B_W A_W Q_W \quad \checkmark
$$

### 18.5 B·A·Q 乘积缩放验证（扭转部分）

扭转部分的缩放需要独立验证，因为 $\mathbf{B}$ 矩阵中平动力列与扭转力矩列涉及不同数量的长度因子：

| 矩阵 | 本项目定义 | 缩放 |
|-------|-----------|------|
| $\mathbf{Q}_{\text{twist}}$ | $\hat{\omega}_j = \cos\theta_j\dot{\phi}_j + \dot{\psi}_j$ | $\times 1$ |
| $\tilde{A}_k$ | $\ell_k\varepsilon_k^2/2 \cdot \hat{A}_k$ | $\times \ell_k\varepsilon_k^2/2$ |
| $\mathbf{B}_{k,\text{twist}}$ | $\ell_k^3 \cdot B_{W,\text{twist}}$, 其中 $B_{W,\text{twist}}$ 条目 $\sim d_\alpha/N_k$ | $\times \ell_k^3$ |

乘积（单项 $j$）：

$$
\mathbf{B}_{k,\text{twist}} \cdot \tilde{A}_k \cdot \mathbf{Q}_{\text{twist}} = \ell_k^3 \cdot \frac{d_\alpha}{N_k} \cdot \frac{\ell_k\varepsilon_k^2}{2}\hat{\omega}_j d_3^j = \frac{\varepsilon_k^2 \ell_k^4}{2N_k}(d_\alpha \cdot d_3)\hat{\omega}_j
$$

**有量纲验证**：扭转力矩密度 $\tau_j = 4\pi\mu r_{\text{pili}}^2\omega_j d_3^j$，Walker Eq. (19) 的贡献：

$$
d_\alpha \cdot \Delta s_k\,\tau_j = \frac{\ell_k L_0}{N_k}\cdot\frac{4\pi\mu(\varepsilon_k\ell_k L_0)^2}{T_0}\hat{\omega}_j(d_\alpha\cdot d_3) = \frac{\varepsilon_k^2\ell_k^3 EI}{2N_k L_0}\hat{\omega}_j(d_\alpha\cdot d_3)
$$

除以弯矩平衡方程的尺度 $EI/(\ell_k L_0)$：

$$
\frac{\varepsilon_k^2\ell_k^3 EI/(2N_k L_0)}{EI/(\ell_k L_0)} = \frac{\varepsilon_k^2\ell_k^4}{2N_k} \quad \checkmark
$$

与上述乘积**完全一致**。

**结论**：弯矩平衡方程的平动和扭转两部分均与 Walker 的 $-EH \cdot B \cdot \mathcal{A} \cdot Q \cdot \dot{z} = \kappa$ **完全等价**。✓

---

## 19. 多纤毛（Multi-Pili）支持

### 19.1 使用方法

在 `main_coupled_simulation.m` 中将 `Pili_Matrix` 扩展为多列，每列定义一根纤毛：

```matlab
Pili_Matrix = [1,   1,   3;        % Type: 直线, 直线, 圆弧
               3,   2,   4;        % Length: 3 um, 2 um, 4 um
               0,   0,   pi/6;    % Deflection from Z
               0,   pi/4, 0;      % Deflection from X
               NaN, NaN, pi/3;    % Shape param 1
               NaN, NaN, NaN;     % Shape param 2
               NaN, NaN, NaN;     % Shape param 3
               5,   6,   11];     % Attachment: Pole 1, Pole 2, Middle y+
```

`N_segments` 会根据各纤毛长度自动计算（`round(L_um * seg_per_um)`）。

### 19.2 理论基础

多纤毛支持在理论框架层面是天然的：

| 耦合机制 | 实现方式 | 说明 |
|----------|----------|------|
| 纤毛-纤毛流体耦合 | 全局 RSM 迁移率矩阵 $S^{\epsilon}$ | $N_{\text{total}} = N_{\text{body}} + \sum_k(N_k+1)$ 个点的全对相互作用 |
| 菌体-纤毛流体耦合 | 同上 | body-pili 交叉项自动包含在 $S^{\epsilon}$ 中 |
| 弹性力矩 | 每根独立计算 Kirchhoff 曲率 | $\kappa_k$ 仅依赖第 $k$ 根纤毛自身的 directors |
| 固接边界条件 | 每根纤毛独立约束到菌体 | 不同纤毛的约束行互不干扰 |
| 扭转阻力 | 块对角 $\tilde{A} = \text{blkdiag}(\tilde{A}_1, \ldots, \tilde{A}_M)$ | 局部 RFT，纤毛间无扭转交叉耦合（物理正确） |
| 自适应基底旋转 | 逐纤毛独立检测与旋转 | 每根纤毛维护自己的 $R_k$ |

### 19.3 状态向量布局

对 $M$ 根纤毛、分段数 $N_1, N_2, \ldots, N_M$：

$$
\mathbf{Z} = \underbrace{[\mathbf{X}_{\text{cm}}(3), \; \mathbf{q}(4)]}_{7} \; \underbrace{[\boldsymbol{\theta}_1(N_1), \boldsymbol{\phi}_1(N_1), \boldsymbol{\psi}_1(N_1)]}_{3N_1} \; \cdots \; \underbrace{[\boldsymbol{\theta}_M(N_M), \boldsymbol{\phi}_M(N_M), \boldsymbol{\psi}_M(N_M)]}_{3N_M}
$$

总自由度：$n = 7 + 3\sum_{k=1}^M N_k$

### 19.4 逐模块多纤毛验证结果

| 模块 | 文件 | 多纤毛处理 | 验证状态 |
|------|------|-----------|----------|
| 状态打包/解包 | `unpack_state.m`, `pack_state_rates.m` | `idx` 偏移累加循环 | ✓ |
| 参数初始化 | `init_coupled_params.m` | `M_pili = size(Pili_Matrix, 2)`, cell 数组 | ✓ |
| 初始状态 | `init_coupled_state.m` | 逐纤毛 Euler 角计算 + 重建验证 | ✓ |
| 几何重建 | `reconstruct_geometry.m` | 循环 k=1:M, `X_all` 正确拼接 | ✓ |
| RSM 迁移率 | `build_global_resistance.m` | 全局矩阵含所有交叉耦合 | ✓ |
| Q 矩阵 | `ode_coupled_rhs.m` | `pt_offset`, `col_offset`, `twist_row` 正确累加 | ✓ |
| B 矩阵 | `ode_coupled_rhs.m` | 力/矩/弯矩平衡行正确偏移 | ✓ |
| 弹性力矩 | `compute_pili_elastic_moments.m` | 每根独立计算，正确索引 | ✓ |
| 固接 BC | `ode_coupled_rhs.m` + `enforce_clamped_geometry.m` | 逐纤毛行替换，互不干扰 | ✓ |
| 自适应基底 | `adaptive_basis_check.m` | 独立 per-pili 检测与旋转 | ✓ |
| 奇异性事件 | `singularity_event` (main) | 遍历所有纤毛取最小 margin | ✓ |
| 结果存储 | main | `cell(M_pili, 1)` 数组 | ✓ |
| 可视化 | main | 逐纤毛独立配色 | ✓ |

### 19.5 附着点冲突检测

`init_coupled_params.m` 在构建纤毛几何前验证 `Pili_Matrix(8,:)` 中是否有重复的附着位编号。若检测到冲突，发出 `warning` 警告（不中断执行），提示用户多根纤毛共享同一附着位可能在物理上不合理。

### 19.6 维度一致性

| 矩阵 | 维度 | 多纤毛公式 |
|-------|------|-----------|
| $\mathbf{Z}$ | $n$ | $7 + 3\sum_k N_k$ |
| $\dot{\mathbf{z}}$ | $n_{\text{dof}}$ | $6 + 3\sum_k N_k$ |
| $\mathbf{v}_{\text{all}}$ | $n_{\text{vel}}$ | $3(N_{\text{body}} + \sum_k(N_k+1)) + \sum_k N_k$ |
| $\mathbf{f}_{\text{all}}$ | $n_{\text{force}}$ | $3(N_{\text{body}} + \sum_k(N_k+1)) + 3\sum_k N_k$ |
| $\mathbf{Q}_{\text{global}}$ | $n_{\text{vel}} \times n_{\text{dof}}$ | — |
| $\mathbf{B}_{\text{global}}$ | $n_{\text{dof}} \times n_{\text{force}}$ | — |
| $\mathbf{L}$ | $n_{\text{dof}} \times n_{\text{dof}}$ | 方阵 |

所有维度经过完整追踪验证，M=1 与 M>1 情况一致。

---

## 20. 审查与验证汇总

### 20.1 已验证项目（28 项）

| # | 审查项 | 涉及文件 | 状态 |
|---|--------|----------|------|
| 1 | $1/(8\pi\mu)$ 吸收 | `build_global_resistance.m` | ✓ |
| 2 | 特征尺度自洽性 | `init_coupled_params.m` | ✓ |
| 3 | 胶囊体精确体积 | `init_coupled_params.m` | ✓ |
| 4 | 四元数左乘约定 | `quat_kinematics.m` | ✓ |
| 5 | Q 矩阵缩放、变换、根部耦合 | `ode_coupled_rhs.m` | ✓ |
| 6 | B 矩阵力平衡梯形权重 | `ode_coupled_rhs.m` | ✓ |
| 7 | B 矩阵弯矩平衡 | `ode_coupled_rhs.m` | ✓ |
| 8 | $\ell_k^3$ 缩放因子 | `ode_coupled_rhs.m` | ✓ |
| 9 | A = Walker 的 A（精确等价） | `build_global_resistance.m` | ✓ |
| 10 | 均匀列权重（Walker 一致） | `build_global_resistance.m` | ✓ |
| 11 | 线性系统符号约定 | `ode_coupled_rhs.m` | ✓ |
| 12 | 固接 BC Jacobian 推导 | `ode_coupled_rhs.m` | ✓ |
| 13 | 几何投影修正 | `enforce_clamped_geometry.m` | ✓ |
| 14 | 自适应旋转矩阵 | `adaptive_basis_check.m` | ✓ |
| 15 | $\kappa_0$ 基底不变性 | `adaptive_basis_check.m` | ✓ |
| 16 | 欧拉角正/反函数一致性 | `pili_directors.m`, `pili_euler_angles.m` | ✓ |
| 17 | `M_weighted \ Q_trans` 优化等价性 | `ode_coupled_rhs.m` | ✓ |
| 18 | 胶囊体/纤毛几何构建 | `capsule_body_building.m` 等 | ✓ |
| 19 | 有量纲输出单位映射 | `main_coupled_simulation.m` | ✓ |
| 20 | 输出系统（日志、数据、动画） | `main_coupled_simulation.m` | ✓ |
| 21 | 计时与 ETA 预估系统 | `main_coupled_simulation.m` | ✓ |
| 22 | 参数化胞体曲面可视化（surf + 透明度 + 光照） | `main_coupled_simulation.m` | ✓ |
| 23 | RSM 力点/位置点重合：死代码清理 | `capsule_body_building.m`, `pili_type_building.m` | ✓ |
| 24 | 逐行量纲一致性检查（$L\dot{z}=R$ 全行） | `ode_coupled_rhs.m`, `build_global_resistance.m` | ✓ |
| 25 | 初始姿态 `alpha_0` 参数化（四元数 + 纤毛几何旋转） | `init_coupled_state.m`, `init_coupled_params.m` | ✓ |
| 26 | 案例命名 `case_name` + 输出目录自动管理 | `main_coupled_simulation.m` | ✓ |
| 27 | 多纤毛（Multi-Pili）全模块支持与维度一致性 | 全部文件 | ✓ |
| 28 | 附着点冲突检测 + 逐纤毛独立配色可视化 | `init_coupled_params.m`, `main_coupled_simulation.m` | ✓ |

### 20.2 Bug 修复历史（全部已修复）

| # | 问题 | 修复日期 | 修复内容 |
|---|------|----------|----------|
| 1 | 纤毛端点列权重 $\Delta s/2$ 不一致 | 2026-02-07 | 改为均匀 $\Delta s$（Walker 一致性） |
| 2 | `inv()` 数值稳定性 | 2026-02-07 | 改为 `\` 运算符 |
| 3 | 全矩阵求逆性能瓶颈 | 2026-02-10 | 返回 `M_weighted`，使用 LU 直接求解 |
| 4 | BUG-3：`A_twist` 坐标系不一致 | 2026-02-10 | `d3_basis` → `d3`（实验室坐标系） |
| 5 | 4 个 PascalCase 文件名不规范 | 2026-02-11 | 重命名为 `snake_case`，更新所有引用 |
| 6 | `capsule_body_building` 含 RSM 无用的 `shift`/力点数据 | 2026-02-11 | 移除 `shift` 参数、`SGLoc/CGLoc`、`gxB/gyB/gzB` |
| 7 | `pili_type_building` 含 RSM 无用的力点字段 | 2026-02-11 | 移除 `gxtail/gytail/gztail`，精简输出结构 |
| 8 | 胞体以散点渲染，不够形象 | 2026-02-11 | 改为参数化球柱体曲面 `surf`，半透明 + Gouraud 光照 |
| 9 | $\tilde{A}_k$ 扭转阻力缺少 $\ell_k$ 因子 | 2026-02-12 | `eps_k^2/2` → `ell_k*eps_k^2/2`（量纲一致性修复） |
| 10 | 初始化中 `reconstruct_geometry` 被循环内重复调用 | 2026-02-12 | 移至循环外单次调用（效率优化） |
| 11 | 多纤毛可视化使用统一蓝色，难以区分 | 2026-02-12 | 引入 10 色可区分色板，逐纤毛独立配色 |
| 12 | 多纤毛附着点冲突无检测 | 2026-02-12 | 在 `init_coupled_params.m` 中添加附着点重复检测与 `warning` |

### 20.3 功能增强历史

| # | 功能 | 日期 | 说明 |
|---|------|------|------|
| 1 | 初始姿态参数化 `alpha_0` | 2026-02-12 | 用户可配置 XZ 平面内胞体倾斜角；`init_coupled_state.m` 从角度计算四元数 + 旋转纤毛几何 |
| 2 | 案例命名 `case_name` | 2026-02-12 | 用户自定义案例名称，自动命名输出目录；留空回退至时间戳 |
| 3 | 专业进度条 `odetpbar.m` | 2026-02-12 | 重写为自包含实现：30 格可视进度条 + 百分比 + 已用时间 + ETA；节流更新 |
| 4 | 多纤毛独立配色可视化 | 2026-02-12 | 10 色可区分色板，per-pili 时序渐变；静态图 + 动画双覆盖 |
| 5 | 附着点冲突检测 | 2026-02-12 | 初始化时检测重复附着位，发出 `warning` 提示 |

### 20.4 综合质量评定

| 维度 | 评分 | 说明 |
|------|------|------|
| **理论正确性** | ★★★★★ | 核心理论完全正确；BUG-3 已修复 |
| **实现正确性** | ★★★★★ | 全部源文件逐行审查，与 Walker 精确对比 |
| **无量纲化** | ★★★★★ | 特征尺度自洽，有量纲输出单位明确 |
| **胞体-纤毛耦合** | ★★★★★ | Q 根部耦合、B 力/矩平衡、固接 BC 双重保障 |
| **多纤毛支持** | ★★★★★ | 全模块维度一致性验证、附着点冲突检测、独立配色 |
| **数值稳定性** | ★★★★★ | LU 直接求解、四元数归一化、几何投影、自适应基底 |
| **代码质量** | ★★★★★ | 统一 `snake_case` 命名、完整头注释、规范排版 |
| **输出系统** | ★★★★★ | 双通道日志、`.mat` 数据、PNG 图、MP4 动画、专业进度条、计时/ETA |
| **可视化** | ★★★★★ | 期刊级排版、绿色胞体/逐纤毛独立色/浅红轨迹、时序渐变、等比例固定坐标轴 |
| **用户体验** | ★★★★★ | 可配置初始姿态、案例命名、自动输出管理、清晰进度反馈 |
| **性能** | ★★★★☆ | LU 优化提速 2–3 倍；可进一步使用 Cholesky |

### 20.5 结论

> 全部 MATLAB 源文件经过多轮全面审查（2026-02-06 至 2026-02-12）。理论框架与数值实现**正确且已完整验证**。已识别的 12 个问题全部修复（含 2026-02-12 发现的扭转阻力缩放量纲不一致、初始化效率优化、多纤毛可视化改进、附着点冲突检测）。多纤毛（Multi-Pili）支持经过逐模块全面验证（状态向量、几何重建、RSM 耦合、Q/B 矩阵、固接 BC、自适应基底旋转），维度一致性在 M=1 和 M>1 情况下均正确。线性系统 $L\dot{z}=R$ 每一行均已通过有量纲追溯验证。当前代码无已知错误。

---

## 21. 参考文献

1. **Walker, B. J., Ishimoto, K., & Gaffney, E. A.** (2020). "Filament mechanics in a half-space via regularised Stokeslet segments." *Physical Review E*, 102(4), 043103.

2. **Cortez, R., Fauci, L., & Medovikov, A.** (2005). "The method of regularized Stokeslets in three dimensions: Analysis, validation, and application to helical swimming." *Physics of Fluids*, 17(3), 031504.

3. **Cortez, R.** (2018). "Regularized stokeslet segments." *Journal of Computational Physics*, 375, 783–796.

4. **Perrin, F.** (1934). "Mouvement brownien d'un ellipsoide." *Journal de Physique et le Radium*, 5(10), 497–511.

---

## 附录 A：关键缩放因子推导

### A.1 为什么 Q 的欧拉角列使用 $\ell_k/N_k$ 而非 $1/N_k$

Walker 在单位长度纤毛（$L=1$）下使用 $1/N$ 作为段长。在本项目以 $L_0$ 为基准的坐标系中，纤毛实际长度为 $\ell_k = L_{\text{pili}}/L_0$，段长为 $\ell_k/N_k$。因此速度计算中的正确缩放因子为 $\ell_k/N_k$，比 Walker 多出一个 $\ell_k$ 倍。

### A.2 为什么 B 的弯矩行使用 $\ell_k^3$ 缩放

归一化后 $\mathbf{B}_k = \mathbf{B}_W$，而 $\mathbf{A}_{\text{our}} = \mathbf{A}_W$，$\mathbf{Q}_{\text{our}} = \ell_k \cdot \mathbf{Q}_W$。乘积：

$$
\ell_k^3 \cdot \mathbf{B}_W \cdot \mathbf{A}_W \cdot (\ell_k \cdot \mathbf{Q}_W) = \ell_k^4 \cdot \mathbf{B}_W \cdot \mathbf{A}_W \cdot \mathbf{Q}_W
$$

其中 $\ell_k^4 = EH_k$（Walker 的弹性流体动力学数）。

### A.3 为什么力平衡使用 $\ell_k/N_k$ 梯形权重

本模型含重力（RHS $\neq 0$），因此必须使用物理正确的梯形积分权重 $[\Delta s/2, \Delta s, \ldots, \Delta s, \Delta s/2]$，而非 Walker 的 $1/(2N^2)$ 权重。Walker 的权重在其系统中可行（力平衡行被固接 BC 替换或 RHS=0），但不适用于含重力的耦合系统。

### A.4 为什么扭转阻力需要额外的 $\ell_k$ 因子

**问题**：$\mathbf{B}_k$ 的 $\ell_k^3$ 统一缩放对平动力列正确，但对扭转力矩列缺少一个 $\ell_k$。

**根本原因**：$\mathbf{B}_k$ 在 Walker 的单位长度坐标下构造时，平动力列和扭转力矩列涉及**不同数量的长度因子**：

| 列类型 | Walker Eq. | $\mathbf{B}$ 条目中的长度因子 | 归一化后的有量纲尺度 |
|--------|------------|-------------------------------|---------------------|
| 平动力 | (17)–(18) | $\Delta s \cdot (x_j - x_i) \sim L^2$ | $\ell_k^2 L_0^2$ |
| 扭转力矩 | (19) | $\Delta s \sim L$ | $\ell_k L_0$ |

用 $\ell_k^3$ 统一缩放：

- **平动乘积**：$\ell_k^3 \cdot 1 \cdot \ell_k = \ell_k^4 = EH_k$ ✓（B×A×Q，其中 A 项 ×1，Q 项 ×$\ell_k$）
- **扭转乘积**：$\ell_k^3 \cdot (\varepsilon_k^2/2) \cdot 1 = \varepsilon_k^2\ell_k^3/2$ ✗（应为 $\varepsilon_k^2\ell_k^4/2$）

**有量纲验证**：

有量纲扭转力矩贡献（一项）：

$$
d_\alpha \cdot \Delta s_k \tau_j = \frac{\ell_k L_0}{N_k}\cdot\frac{4\pi\mu(\varepsilon_k\ell_k L_0)^2}{T_0}\hat{\omega}_j(d_\alpha\cdot d_3) = \frac{\varepsilon_k^2\ell_k^3 EI}{2N_k L_0}\hat{\omega}_j(d_\alpha\cdot d_3)
$$

弯矩平衡方程的尺度为 $EI/(\ell_k L_0)$，因此正确的无量纲扭转项为：

$$
\frac{\varepsilon_k^2\ell_k^3 EI/(2N_k L_0)}{EI/(\ell_k L_0)} = \frac{\varepsilon_k^2\ell_k^4}{2N_k}
$$

而原代码给出 $\varepsilon_k^2\ell_k^3/(2N_k)$，**差一个 $\ell_k$**。

**修复**：将缺失的 $\ell_k$ 补入 $\tilde{A}_k$（保持 $\mathbf{B}_k$ 的统一 $\ell_k^3$ 缩放不变）：

$$
\tilde{A}_k = \frac{\varepsilon_k^2}{2}\hat{A}_k \quad\longrightarrow\quad \tilde{A}_k = \frac{\ell_k\,\varepsilon_k^2}{2}\hat{A}_k
$$

验证：$\ell_k^3 \cdot (\ell_k\varepsilon_k^2/2) \cdot 1 / N_k = \varepsilon_k^2\ell_k^4/(2N_k)$ ✓

**影响评估**：由于 $\varepsilon_k \sim 10^{-3}$，扭转贡献 $O(\varepsilon_k^2) \sim 10^{-6}$ 为弯曲贡献的百万分之一。修复前的 $O(\ell_k)$ 误差（$\ell_k \sim 1\text{–}5$）改变了一个已经极小的量，**对仿真结果的实际影响可忽略**，但量纲一致性现已严格保证。

### A.5 线性系统 $L\dot{z}=R$ 逐行量纲一致性验证

| 行 | 方程 | 两端尺度 | 有量纲追溯结果 | 状态 |
|----|------|----------|---------------|------|
| 1–3 | 力平衡 | $F_0 = EI/L_0^2$ | 胞体 $\sum F_0\hat{f}$ + 纤毛 $\sum(\ell_k/N_k)F_\ell\hat{f} \cdot L_0 = \sum(\ell_k/N_k)F_0\hat{f}$；RHS = $Eg \cdot F_0$ | ✓ |
| 4–6 | 力矩平衡 | $F_0 L_0 = EI/L_0$ | 胞体 $F_0 L_0(\hat{r}\times\hat{f})$ + 纤毛 $(\ell_k/N_k)F_0 L_0(\hat{r}\times\hat{f})$；RHS = 0 | ✓ |
| 7+ (平动) | 纤毛弯矩—力贡献 | $EI/(\ell_k L_0)$ | $\ell_k^3 B_W \hat{f}^P = B_W(\ell_k^3\hat{f}^P) = B_W\hat{f}^W = \hat{\kappa}$ | ✓ |
| 7+ (扭转) | 纤毛弯矩—扭矩贡献 | $EI/(\ell_k L_0)$ | $\ell_k^3 \cdot \ell_k\varepsilon_k^2/(2N_k) = \varepsilon_k^2\ell_k^4/(2N_k)$（修复后） | ✓ |

---

## 许可证

MIT License (c) 2026 Lucien
