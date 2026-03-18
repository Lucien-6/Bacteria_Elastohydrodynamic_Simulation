# 细菌胞体-纤毛耦合弹性流体动力学仿真

> **作者**: Lucien (lucien-6@qq.com)  
> **创建日期**: 2026-02-06  
> **更新日期**: 2026-03-18  
> **版本**: 3.3

---

## 目录

1. [项目概述](#1-项目概述)
2. [快速开始](#2-快速开始)
3. [架构与数据流](#3-架构与数据流)
4. [物理问题与建模](#4-物理问题与建模)
5. [无量纲化](#5-无量纲化)
6. [状态向量设计](#6-状态向量设计)
7. [几何重构](#7-几何重构)
8. [流体力学：正则化 Stokeslet 线段法 (RSS)](#8-流体力学正则化-stokeslet-线段法-rss)
9. [运动学矩阵 Q](#9-运动学矩阵-q)
10. [力/力矩平衡矩阵 B](#10-力力矩平衡矩阵-b)
11. [弹性恢复力矩（Kirchhoff 弹性杆）](#11-弹性恢复力矩kirchhoff-弹性杆)
12. [耦合线性系统组装与求解](#12-耦合线性系统组装与求解)
13. [固支边界条件](#13-固支边界条件)
14. [自适应基底旋转（Euler 角奇异性规避）](#14-自适应基底旋转euler-角奇异性规避)
15. [时间积分](#15-时间积分)
16. [胞体布朗动力学（分裂步 Euler-Maruyama）](#16-胞体布朗动力学分裂步-euler-maruyama)
17. [输出系统：日志、数据导出与动画](#17-输出系统日志数据导出与动画)
18. [流场可视化](#18-流场可视化)
19. [文件结构](#19-文件结构)
20. [与 Walker (2020) 的标度一致性](#20-与-walker-2020-的标度一致性)
21. [多纤毛支持](#21-多纤毛支持)
22. [验证总结](#22-验证总结)
23. [参考文献](#23-参考文献)
24. [附录 A：关键标度推导](#附录-a关键标度推导)
25. [许可证](#许可证)

---

## 1. 项目概述

本项目模拟一个具有弹性纤毛（以 Kirchhoff 弹性杆建模）的胶囊形细菌在 Stokes 流中重力沉降的过程。细菌胞体被视为刚性胶囊体；每根纤毛为不可伸长的弹性杆，采用 Euler 角参数化。

**流体力学方法：**

| 组件                                 | 方法                           | 参考文献             |
| ------------------------------------ | ------------------------------ | -------------------- |
| 胞体表面（$N_{\text{body}}$ 个力点） | RSM 点力法                     | Cortez et al. (2005) |
| 纤毛丝（每根 $N_k$ 段）              | RSS 精确线段积分 + Rotlet 线段 | Cortez (2018)        |

RSS（Regularized Stokeslet Segments，正则化 Stokeslet 线段法）对正则化 Stokeslet 核函数沿每段直线段进行精确积分，将正则化参数 $\varepsilon$ 与线段长度 $h$ 解耦。这使得纤毛仅需约 15 段/微米即可保持高精度，与逐点方法相比大幅降低了计算成本。

**核心特性：**

- **刚体胞体**：胶囊几何体（半球端盖 + 圆柱中段），四元数参数化姿态（无万向锁）
- **弹性纤毛**：Kirchhoff 不可伸长弹性杆 + Euler 角参数化（Walker 2020）
- **流体力学**：胞体 RSM 点力 + 纤毛 RSS Stokeslet 线段 + Rotlet 线段耦合（扭转—平移）
- **布朗动力学**：胞体热涨落（分裂步 Euler-Maruyama，RSM 摩擦矩阵 → Schur 补解耦 → Cholesky 扩散张量）
- **时间积分**：MATLAB `ode15s`（隐式刚性求解器）+ 固定 `dt_B` 布朗节点
- **自适应基底**：逐纤毛 Euler 角基底旋转以规避奇异性
- **固支边界条件**：双重施加（运动学速率约束 + 几何投影）
- **输出**：双通道日志、`.mat` 数据、PNG 图像、3D MP4 动画、JSON 配置快照
- **可视化**：参数化光滑胞体表面（半透明，Gouraud 着色）+ 逐纤毛着色折线，含时间梯度
- **流场可视化**：末状态双平面流场图（XZ + YZ，各含 3D 表面视图与 2D 着色图；速度大小 jet 着色 + 归一化方向箭头；胞体内部遮罩）
- **动画布局**：三面板 1×3（3D + XZ 投影 + YZ 投影），动态图像尺寸
- **初始姿态**：用户可配置 `alpha_0`（胞体在 XZ 平面内的倾斜角）
- **外部配置**：支持 JSON 配置文件，运行时自动保存配置快照至输出目录（可追溯复现）
- **参数验证**：启动前自动执行参数合理性预检（错误聚合 + 警告分级）
- **无人值守模式**：`auto_confirm = true` 跳过所有交互对话框，支持批量参数扫描
- **纤毛结构体接口**：`parse_pili_config` 提供命名字段的用户友好纤毛定义方式
- **函数化封装**：`run_simulation()` 将全部仿真逻辑封装为可编程调用的函数
- **进度反馈**：动画渲染和流场计算配备实时文本进度条（含 ETA）
- **输出目录保护**：检测已有仿真数据时弹出覆盖/重命名/中止对话框
- **参数扫描套件**（`batch sweep/`）：6 个独立脚本分别扫描初始倾斜角、胞体长轴、纤毛长度、纤毛附着点位、纤毛数量、纤毛杨氏模量，均采用 headless 无头模式批量调用 `run_simulation()`
- **JSON 批量运行器**：`batch_run_from_json.m` 自动检索用户指定文件夹内全部 JSON 配置文件并依次执行，单工况失败自动跳过继续，运行结束输出对齐汇总表及批次日志

---

## 2. 快速开始

### 2.1 方式一：编辑内联参数（默认）

```matlab
% 1. 在 MATLAB 中打开 main_coupled_simulation.m
% 2. 在 "Inline parameter definitions" 区域调整 cfg 结构体字段：

cfg.case_name  = 'MyCase';
cfg.alpha_0    = pi/4;
cfg.N_body     = 1000;

% 3. 定义纤毛配置（命名结构体，推荐）：
pili_config(1).type      = 'linear';
pili_config(1).length_um = 3;
pili_config(1).site      = 1;
cfg.Pili_Matrix = parse_pili_config(pili_config);

% 4. 运行：
main_coupled_simulation
```

### 2.2 方式二：从 JSON 文件加载配置

```matlab
% 在 main_coupled_simulation.m 中设置：
config_file = 'configs/my_case.json';    % 指向 JSON 配置文件
% 运行即可，所有参数从 JSON 加载（参见 config_example.json 模板）
```

### 2.3 方式三：编程式调用 / 批量扫描

```matlab
% 直接调用 run_simulation() 函数：
cfg = load_sim_config('configs/base.json');
cfg.alpha_0      = pi/4;
cfg.auto_confirm = true;     % 无人值守模式
[results, params, output_dir] = run_simulation(cfg);

% 参数扫描脚本参见 batch sweep/ 文件夹（6 个扫描脚本 + 1 个 JSON 批量运行器）
```

**输出**（保存至 `output/<case_name>/`）：

```
output/MyCase/
├── simulation_results.mat    % 完整结果 + 参数（v7.3 / HDF5）
├── timeseries.csv            % 标量时间序列（10 列，纯文本）
├── config.json               % 配置快照（可直接复用以复现运行）
├── summary_plot.png          % 静态总览图（600 dpi）
├── initial_pose.png          % 初始姿态预览图（600 dpi，预览开启时）
├── bacteria_3d_motion.mp4    % 3D 运动动画（三面板 1×3，动态尺寸，动画开启时）
├── flowfield_panels.png      % 双平面流场 2×2 面板图（600 dpi）
└── simulation.log            % 完整运行日志
```

**用户可配置参数**（位于 `main_coupled_simulation.m` 或 JSON 配置文件）：

| 参数 | 默认值 | 说明 |
| ---- | ------ | ---- |
| `case_name` | `'Test'` | 案例名称（作为输出文件夹名；为空则使用时间戳） |
| `minor_axis` | `0.4e-6` | 胞体椭球短半轴（体宽半径）[m] |
| `major_axis` | `1.25e-6` | 胞体椭球长半轴（体长半径）[m] |
| `alpha_0` | `0` | XZ 平面内初始胞体倾斜角 [rad] |
| `N_body` | `1000` | 胞体表面力点数 |
| `seg_per_um` | `15` | 纤毛每微米段数（RSS 精确积分） |
| `Pili_Matrix` | 单根直线 | 纤毛配置矩阵 [8 × M]（可用 `parse_pili_config` 从结构体生成） |
| `E` | `1.768e7` | 杨氏模量 [Pa] |
| `r_pili` | `3e-9` | 纤毛半径 [m] |
| `sigma_p` | `0.2` | 泊松比 |
| `mu` | `1e-3` | 流体动力粘度 [Pa·s] |
| `rho_bac` | `1180.0` | 细菌密度 [kg/m³] |
| `rho_fluid` | `998.232` | 流体密度 [kg/m³] |
| `T_obs` | `30.0` | 总观测时间 [s] |
| `AbsTol` / `RelTol` | `1e-4` | ODE 求解器容差 |
| `delta` | `pi/20` | Euler 角奇异性阈值 [rad] |
| `tlim` | `120` | ODE 段超时限制 [s] |
| `preview_initial` | `true` | 积分前显示初始姿态预览 |
| `generate_animation` | `true` | 生成 MP4 动画 |
| `animation_fps` | `30` | 动画帧率 |
| `animation_skip` | `1` | 动画帧跳步间隔 |
| `enable_brown` | `false` | 启用胞体布朗运动 |
| `Temper` | `293.15` | 流体绝对温度 [K]（默认 20°C） |
| `dt_B` | `dt/2` | 布朗时间步 [s]（`enable_brown = true` 时生效） |
| `auto_confirm` | `false` | 无人值守模式（跳过所有交互对话框） |
| `flowfield_grid_spacing` | `0.05` | 流场网格间距 [μm] |
| `flowfield_grid_range` | `10.0` | 流场网格范围 [μm] |
| `flowfield_quiver_skip` | `5` | 流场箭头降采样倍率 |
| `flowfield_chunk_size` | `2000` | 流场评估分块大小 |

---

## 3. 架构与数据流

### 3.1 模块层级

```mermaid
graph TD
    ENTRY[main_coupled_simulation.m<br/>Entry Point — cfg 定义] --> PC[parse_pili_config.m]
    ENTRY --> LC[load_sim_config.m]
    ENTRY --> RUN[run_simulation.m<br/>Core Engine]

    RUN --> VAL[validate_simulation_params.m]
    RUN --> SC[save_sim_config.m]
    RUN --> B[init_coupled_params.m]
    RUN --> C[init_coupled_state.m]
    RUN --> D[ode15s + ode_coupled_rhs.m]
    RUN --> E[adaptive_basis_check.m]
    RUN --> F[enforce_clamped_geometry.m]
    RUN --> V[compute_plot_flowfield.m]
    RUN --> PB[render_progress_bar.m]
    RUN --> G[Output: .mat / .png / .mp4 / .log / .json]

    B --> H[capsule_body_building.m]
    B --> I[pili_building.m]
    B --> W[compute_body_brownian_params.m]
    I --> J[pili_type_building.m]
    I --> K[compute_curvature_binormal.m]

    C --> L[reconstruct_geometry.m]
    C --> F

    D --> L
    D --> M[build_global_resistance.m]
    D --> N[compute_pili_elastic_moments.m]
    D --> P[pack_state_rates.m]

    V --> L
    V --> D
    V --> M
    V --> PB

    L --> Q[pili_directors.m]
    L --> R[pili_spatial_coords.m]
    L --> S[unpack_state.m]
    L --> T[quat_to_rotmat.m]

    E --> Q
    E --> U[pili_euler_angles.m]
    E --> S

    BATCH["batch sweep/ 文件夹<br/>batch_sweep_alpha / major_axis / pili_length<br/>pili_site / pili_count / pili_E<br/>batch_run_from_json"] --> RUN
```

### 3.2 功能分层

| 层级       | 文件                                                                                                 | 职责                                          |
| ---------- | ---------------------------------------------------------------------------------------------------- | --------------------------------------------- |
| 入口脚本   | `main_coupled_simulation.m`                                                                          | cfg 结构体定义、配置来源选择、调用 `run_simulation` |
| 仿真引擎   | `run_simulation.m`                                                                                   | 积分循环、布朗积分、后处理、输出（可编程调用）|
| 配置管理   | `load_sim_config.m`, `save_sim_config.m`, `parse_pili_config.m`                                      | JSON 读写、纤毛结构体→矩阵转换               |
| 参数验证   | `validate_simulation_params.m`                                                                       | 启动前参数合理性预检（错误聚合 + 警告分级）   |
| 初始化     | `init_coupled_params.m`, `init_coupled_state.m`                                                      | 特征尺度、几何构建、布朗参数、初始状态        |
| 核心计算   | `ode_coupled_rhs.m`                                                                                  | Q、A、B 矩阵组装，线性系统求解               |
| 几何工具   | `reconstruct_geometry.m`, `pili_*.m`, `quat_*.m`                                                     | 状态 ↔ 空间坐标转换                          |
| 流体力学   | `build_global_resistance.m`                                                                          | RSS 迁移率 + Rotlet 线段迁移率 + 扭转 RFT    |
| 弹性力学   | `compute_pili_elastic_moments.m`                                                                     | Kirchhoff 弹性杆曲率                          |
| 布朗参数   | `compute_body_brownian_params.m`                                                                     | RSM 摩擦矩阵 → Schur 补 → Cholesky 扩散因子  |
| 边界条件   | `enforce_clamped_geometry.m`                                                                         | 固支运动学约束 + 几何投影                     |
| 数值稳定性 | `adaptive_basis_check.m`                                                                             | Euler 角奇异性规避                            |
| 几何生成   | `capsule_body_building.m`, `pili_building.m`, `pili_type_building.m`, `compute_curvature_binormal.m` | 胞体/纤毛离散化                               |
| 流场可视化 | `compute_plot_flowfield.m`                                                                           | 末状态双平面流场评估与绘图（2×2 面板）        |
| 进度显示   | `odetpbar.m`, `render_progress_bar.m`                                                                | ODE 段内进度条 + 通用文本进度条（动画/流场）  |
| 批量扫描   | `batch sweep/batch_sweep_alpha.m`, `batch sweep/batch_sweep_major_axis.m`, `batch sweep/batch_sweep_pili_length.m`, `batch sweep/batch_sweep_pili_site.m`, `batch sweep/batch_sweep_pili_count.m`, `batch sweep/batch_sweep_pili_E.m`, `batch sweep/batch_run_from_json.m` | 6 个独立参数扫描脚本（初始角/胞体长轴/纤毛长度/附着点位/数量/杨氏模量）+ JSON 批量运行器 |
| 工具函数   | `skew_matrix.m`, `pack_state_rates.m`, `unpack_state.m`                                              | 线性代数辅助、状态打包                        |

---

## 4. 物理问题与建模

### 4.1 问题描述

模拟具有弹性纤毛的细菌在 Stokes 流中重力沉降：

- **胞体**：刚性胶囊体（两个半球端盖 + 圆柱中段），离散为 $N_{\text{body}}$ 个表面力点（默认 1000）
- **纤毛**：$M$ 根不可伸长的 Kirchhoff 弹性杆，每根含 $N_k$ 段（$N_k + 1$ 个节点）
- **流体**：Stokes 流（零雷诺数）；所有流体力学相互作用通过 RSS（纤毛）+ RSM（胞体）实现

### 4.2 方法来源

| 组件                       | 来源                           | 方法                                                 |
| -------------------------- | ------------------------------ | ---------------------------------------------------- |
| 胞体形状、表面点、附着位点 | 母项目 `Stokes Brown Flexible` | `capsule_body_building`, `pili_building`             |
| 纤毛弹性杆模型             | Walker (2020)                  | Kirchhoff 弹性杆 + Euler 角参数化                    |
| 流体力学（平移）           | Cortez (2005, 2018)            | 胞体：RSM 点力；纤毛：RSS 精确线段积分 + Rotlet 线段 |
| 流体力学（扭转）           | Walker (2020)                  | 阻力理论 (RFT)                                     |
| 胞体姿态                   | 本项目                         | 四元数参数化（规避万向锁）                           |
| 时间积分                   | Walker (2020)                  | MATLAB `ode15s` 隐式求解器                           |

### 4.3 核心耦合机制

1. **运动学耦合**：纤毛根部固支于胞体表面；根节点速度 = 胞体刚体运动速度
2. **流体力学耦合**：全局迁移率矩阵（胞体 RSM + 纤毛 RSS + Rotlet）捕获所有胞体—胞体、胞体—纤毛和纤毛—纤毛相互作用
3. **弹性耦合**：纤毛弹性恢复力矩通过流体力学反力反馈影响胞体运动

---

## 5. 无量纲化

### 5.1 特征尺度

所有量均以胞体半长轴 $a$（`major_axis`）进行无量纲化：

| 尺度     | 符号  | 定义                                |
| -------- | ----- | ----------------------------------- |
| 长度     | $L_0$ | $a$                                 |
| 弯曲刚度 | $EI$  | $E \cdot \pi r_{\text{pili}}^4 / 4$ |
| 时间     | $T_0$ | $8\pi\mu L_0^4 / EI$                |
| 速度     | $U_0$ | $L_0 / T_0$                         |
| 力       | $F_0$ | $EI / L_0^2$                        |

### 5.2 关键无量纲参数

| 参数        | 符号          | 定义                         |
| ----------- | ------------- | ---------------------------- |
| 弹性-重力数 | $Eg$          | $F_{\text{sink}} L_0^2 / EI$ |
| 细长比      | $\varepsilon$ | $r_{\text{pili}} / L_0$      |
| 泊松比      | $\sigma$      | 用户指定                     |

---

## 6. 状态向量设计

### 6.1 布局

$$
\mathbf{Z} = \underbrace{[\mathbf{X}_{\text{cm}}(3), \; \mathbf{q}(4)]}_{7} \; \underbrace{[\boldsymbol{\theta}_1, \boldsymbol{\phi}_1, \boldsymbol{\psi}_1]}_{3N_1} \; \cdots \; \underbrace{[\boldsymbol{\theta}_M, \boldsymbol{\phi}_M, \boldsymbol{\psi}_M]}_{3N_M}
$$

总自由度数：$n = 7 + 3\sum_{k=1}^M N_k$

### 6.2 速率向量（约化）

$$
\dot{\mathbf{z}} = [\mathbf{V}(3), \; \boldsymbol{\Omega}(3), \; \dot{\boldsymbol{\theta}}_1, \dot{\boldsymbol{\phi}}_1, \dot{\boldsymbol{\psi}}_1, \; \ldots]
$$

自由度数：$n_{\text{dof}} = 6 + 3\sum_k N_k$

### 6.3 初始姿态参数 `alpha_0`

角度 `alpha_0` 定义胞体长轴（体坐标系 X 轴）在 XZ 平面内相对于 +X 方向的初始倾斜角：

- `pi/2` → 竖直（+Z 方向，标准沉降姿态）
- `0` → 水平（+X 方向）
- `-pi/2` → 倒置（-Z 方向）

---

## 7. 几何重构

`reconstruct_geometry.m` 将状态向量 $\mathbf{Z}$ 转换为实验室坐标系中的位置：

1. **胞体**：$\mathbf{X}_{\text{cm}}$ + 四元数 $\mathbf{q}$ → 旋转矩阵 $R$ → 胞体表面点 $\mathbf{X}_i = \mathbf{X}_{\text{cm}} + R \cdot \mathbf{x}_i^{\text{local}}$
2. **纤毛**：由胞体表面获取根部位置 + Euler 角 → `pili_spatial_coords` → 实验室坐标系节点位置
3. **标架**：`pili_directors` → 材料标架 $(\mathbf{d}_1, \mathbf{d}_2, \mathbf{d}_3)$，其中 $\mathbf{d}_3$ 为切线方向

输出：`geom` 结构体，包含 `X_all`、`X_cm`、`X_body`、`X_pili`、标架、角度等。

---

## 8. 流体力学：正则化 Stokeslet 线段法 (RSS)

### 8.1 RSM 核函数 (Cortez 2005) — 胞体表面

$$
S_{ij}^{\varepsilon} = \frac{(|\mathbf{r}|^2 + 2\varepsilon_{ij}^2)\,\mathbf{I} + \mathbf{r}\mathbf{r}^T}{(|\mathbf{r}|^2 + \varepsilon_{ij}^2)^{3/2}}
$$

仅应用于迁移率矩阵中胞体表面对应的列。

### 8.2 RSS 线段核函数 (Cortez 2018)

对于线段 $\mathbf{y}(\alpha) = \mathbf{y}_0 - \alpha\mathbf{v}$（$\mathbf{v} = \mathbf{y}_0 - \mathbf{y}_1$，$L = |\mathbf{v}|$），力密度从 $\mathbf{f}_0$ 到 $\mathbf{f}_1$ 线性变化：

$$
\mathbf{u}(\hat{\mathbf{x}}) = \mathbf{M}_1\,\mathbf{f}_0 + \mathbf{M}_2\,\mathbf{f}_1
$$

$$
\mathbf{M}_2 = L\left[(T_{1,-1} + \varepsilon^2 T_{1,-3})\mathbf{I} + T_{1,-3}\,\mathbf{x}_0\mathbf{x}_0^T + T_{2,-3}(\mathbf{x}_0\mathbf{v}^T + \mathbf{v}\mathbf{x}_0^T) + T_{3,-3}\,\mathbf{v}\mathbf{v}^T\right]
$$

$$
\mathbf{M}_1 = L\left[(T_{0,-1} + \varepsilon^2 T_{0,-3})\mathbf{I} + \cdots\right] - \mathbf{M}_2
$$

$T_{n,q}$ 值通过递推关系计算（Cortez 2018, Eq. 7），起始于：

- $T_{0,-1} = \frac{1}{L}\ln\frac{H_1}{H_0}$，其中 $H = LR + \mathbf{x}_0 \cdot \mathbf{v}$
- $T_{0,-3} = \frac{1}{R_0 H_0} - \frac{1}{R_1 H_1}$

**核心优势**：$\mathbf{M}_1$、$\mathbf{M}_2$ 已包含线段长度积分权重——无需额外的 $\Delta s$ 列权重。

### 8.3 全局迁移率矩阵

$$
M_{\text{weighted}} = \begin{pmatrix} M_{BB} & M_{BP} \\ M_{PB} & M_{PP} \end{pmatrix}
$$

| 分块     | 评估点 | 力源 | 方法     | $\varepsilon$         |
| -------- | ------ | ---- | -------- | --------------------- |
| $M_{BB}$ | 胞体   | 胞体 | RSM 点力 | $h_{\text{mean}}/3$   |
| $M_{PB}$ | 纤毛   | 胞体 | RSM 点力 | $h_{\text{mean}}/3$   |
| $M_{BP}$ | 胞体   | 纤毛 | RSS 线段 | $r_{\text{pili}}/L_0$ |
| $M_{PP}$ | 纤毛   | 纤毛 | RSS 线段 | $r_{\text{pili}}/L_0$ |

**配对规则**：$\varepsilon_{ij} = \varepsilon_j$（由力源点决定）。

**向量化组装**：所有 $T_{n,q}$ 值和 $\mathbf{M}_1$/$\mathbf{M}_2$ 分量以 $[N_{\text{total}} \times N_{\text{seg,all}}]$ 矩阵形式计算，然后散射至全局迁移率矩阵中。

### 8.4 Rotlet 线段迁移率

捕获纤毛上扭矩密度到所有评估点平移速度的流体力学耦合（Cortez 2018, Sec. 2.4）：

$$
\mathbf{u}(\hat{\mathbf{x}}) = -L\left[c_a\,[\mathbf{x}_0]_\times + c_b\,[\mathbf{v}]_\times\right]\boldsymbol{\tau}
$$

其中 $c_a = 2T_{0,-3} + 3\varepsilon^2 T_{0,-5}$，$c_b = 2T_{1,-3} + 3\varepsilon^2 T_{1,-5}$。

$T_{n,-5}$ 通过额外的递推关系计算，需要 $d^2 = |\mathbf{x}_0 \times \mathbf{v}|^2 + L^2\varepsilon^2$。

$M_{\text{rotlet}} \in \mathbb{R}^{3N_{\text{total}} \times 3N_{\text{seg}}}$ 将每段的三维扭矩密度映射到所有点的平移速度。

**线性系统中的集成**：

$$
\mathbf{f} = M_{\text{stokeslet}}^{-1}\left(\mathbf{Q}_{\text{trans}} - M_{\text{rotlet}} \cdot \tilde{A}_{\text{twist}} \cdot \mathbf{Q}_{\text{twist}}\right) \dot{\mathbf{z}}
$$

**物理影响**：对于典型的细菌纤毛（$\varepsilon \sim 10^{-3}$），Rotlet 耦合量级为 $O(\varepsilon^4) \sim 10^{-12}$，实际上可忽略，但确保了理论的完备性。

### 8.5 正则化参数

$$
\varepsilon_j = \begin{cases} h_{\text{body}}^{\text{mean}} / 3 & \text{胞体表面点} \\ r_{\text{pili}} / L_0 & \text{纤毛节点} \end{cases}
$$

**胞体**（RSM，力源在封闭曲面上）：$h_{\text{body}}^{\text{mean}}$ 为胞体表面所有点最近邻距离的均值（平均网格间距 $\Delta s_{\text{avg}}$），$h_{\text{mean}}/3$ 对应比值 $\varepsilon/\Delta s_{\text{avg}} \approx 0.333$，为针对当前胶囊几何（$a/b = 3.125$，$N_{\text{body}} = 1000$）数值调试所得的经验最优值，在阻力矩阵精度与数值稳定性之间取得最佳平衡。$\varepsilon \to 0$ 的极限存在（退化为奇异 BEM），故 $\varepsilon$ 为纯数值参数，应通过与解析解（如 Stokes 球阻力）对比验证。

**纤毛**（RSS，力源在一维曲线上）：$\varepsilon = r_{\text{pili}}/L_0$ 是**物理参数**，代表力作用区域的空间延展半径（纤毛横截面半径）。Cortez et al. (2005) 明确指出，当力作用于曲线时 $\varepsilon \to 0$ 的极限不存在（速度场发散），因此 $\varepsilon$ 不可趋向零，而应取纤毛的物理横截面半径。细长比 $\varepsilon/\ell = r_{\text{pili}}/L_{\text{pili}} \approx 10^{-3}$ 远小于 1，满足细长体假设；段长与正则化参数之比 $\Delta s/\varepsilon = \Delta s \cdot L_0 / r_{\text{pili}} \approx 22$，RSS 精确积分在此范围内精度充分。

纤毛根部与胞体表面的间隙：$d_{\text{gap}} = 1.1 \times (\varepsilon_{\text{body}} \cdot L_0 + r_{\text{pili}})$。

### 8.6 平移阻力

$\mathbf{A}_{\text{trans}} = M_{\text{weighted}}^{-1}$ — 通过 LU 分解计算：`M_weighted \ Q_trans`（避免显式 $O(N^3)$ 矩阵求逆）。

### 8.7 扭转阻力 (RFT)

逐纤毛扭转阻力（阻力理论）：

$$
\tilde{A}_k = \frac{\ell_k\,\varepsilon_k^2}{2}\,\hat{A}_k, \quad [\hat{A}_k]_{3(j-1)+1:3j, \, j} = \mathbf{d}_3^{(j)}
$$

---

## 9. 运动学矩阵 Q

$\mathbf{Q}_{\text{global}} \in \mathbb{R}^{n_{\text{vel}} \times n_{\text{dof}}}$ 将速率向量 $\dot{\mathbf{z}}$ 映射为物理速度：

- **胞体点**：$\mathbf{v}_i = \mathbf{V} + \boldsymbol{\Omega} \times (\mathbf{x}_i - \mathbf{X}_{\text{cm}})$
- **纤毛节点**：由 Walker 的空间 Q 矩阵构建，通过自适应旋转 $R_k$ 从纤毛基底变换至实验室坐标系
- **扭转速率**：$\psi$ 速率通过 $\cos\theta$ 与 $\theta$ 速率耦合（Euler 角运动学）

平移子块：$\mathbf{Q}_{\text{trans}} = \mathbf{Q}[1:3N_{\text{total}}, :]$
扭转子块：$\mathbf{Q}_{\text{twist}} = \mathbf{Q}[3N_{\text{total}}+1:\text{end}, :]$

---

## 10. 力/力矩平衡矩阵 B

$\mathbf{B}_{\text{global}} \in \mathbb{R}^{n_{\text{dof}} \times n_{\text{force}}}$ 将分布力/力矩映射至速率向量方程：

| 行  | 方程                                     | 右端项                  |
| --- | ---------------------------------------- | ----------------------- |
| 1–3 | 力平衡                                   | $-Eg\,\hat{\mathbf{g}}$ |
| 4–6 | 力矩平衡                                 | 0（无外力矩）           |
| 7+  | 纤毛力矩平衡（每根 $N_k$ 段 × 3 个分量） | 弹性曲率 $\kappa$       |

胞体和纤毛的力贡献使用梯形求积权重。纤毛力矩平衡遵循 Walker (2020) 的公式，采用 $\ell_k^3$ 标度。

---

## 11. 弹性恢复力矩（Kirchhoff 弹性杆）

`compute_pili_elastic_moments.m` 计算离散曲率与固有曲率的偏差：

$$
\kappa_\alpha^{(j)} = \mathbf{d}_\alpha^{(j)} \cdot \mathbf{d}_3'^{(j)} - \kappa_{\alpha,0}^{(j)}, \quad \alpha = 1, 2
$$

$$
\kappa_3^{(j)} = \frac{1}{1+\sigma}\left(\mathbf{d}_2^{(j)} \cdot \mathbf{d}_1'^{(j)} - \kappa_{3,0}^{(j)}\right)
$$

---

## 12. 耦合线性系统组装与求解

每个 ODE 时间步的核心方程：

$$
\mathbf{L}\,\dot{\mathbf{z}} = \mathbf{R}
$$

其中：

- $\mathbf{L} = -\mathbf{B} \cdot \mathbf{AQ}$
- $\mathbf{AQ} = \begin{bmatrix} M_{\text{weighted}}^{-1} \cdot \mathbf{Q}_{\text{trans,corr}} \\ \tilde{A}_{\text{twist}} \cdot \mathbf{Q}_{\text{twist}} \end{bmatrix}$
- $\mathbf{Q}_{\text{trans,corr}} = \mathbf{Q}_{\text{trans}} - M_{\text{rotlet}} \cdot \tilde{A}_{\text{twist}} \cdot \mathbf{Q}_{\text{twist}}$
- $\mathbf{R} = [-Eg\,\hat{\mathbf{g}}; \; \mathbf{0}; \; \boldsymbol{\kappa}_1; \ldots; \boldsymbol{\kappa}_M]$

通过 LU 分解求解：`lin_sys \ R_global`。

---

## 13. 固支边界条件

纤毛根部固支于胞体表面。双重施加：

1. **速率层面**（在 `ode_coupled_rhs.m` 步骤 7 中）：将每根纤毛力矩平衡块的前 3 行替换为：

$$
[\dot{\theta}_1,\; \dot{\phi}_1,\; \dot{\psi}_1]^T = J_e \cdot \boldsymbol{\Omega}
$$

其中 $J_e \in \mathbb{R}^{3 \times 3}$ 为自适应基底下 Euler 角速率到角速度的雅可比矩阵，将胞体角速度 $\boldsymbol{\Omega}$ 映射为根节点的三个 Euler 角速率。

2. **位置层面**（在 `enforce_clamped_geometry.m` 中）：每个 ODE 步后，几何投影修正纤毛—胞体连接处的约束漂移。

---

## 14. 自适应基底旋转（Euler 角奇异性规避）

`adaptive_basis_check.m` 监控每根纤毛的 $|\theta_1^{(k)}|$。当 $\theta_1$ 趋近 0 或 $\pi$（在阈值 $\delta$ 范围内）时，更新旋转基底 $R_k$ 以将工作点移离奇异性。所有 Euler 角在新基底下重新计算。

---

## 15. 时间积分

```
算法：
  对每个时间段 [t_n, t_{n+1}]：
      1. 解包 Z → X_cm, q, {θ,φ,ψ}_k
      2. 重构几何 → X_all, 标架
      3. 构建迁移率 M_weighted（胞体 RSM + 纤毛 RSS）、M_rotlet、A_twist
      4. 构建 Q_global（运动学）
      5. 求解 AQ = [M_weighted\Q_trans_corr; A_twist*Q_twist]
      6. 构建 B_global（力/力矩平衡）
      7. 计算弹性力矩 κ（逐纤毛）
      8. 组装线性系统 L = -B·AQ, R = [-Eg·ĝ; 0; κ]
      9. 施加固支边界条件（替换首段行）
     10. 求解 L·ż = R
     11. 提取 V, Ω → 计算 dq/dt = ½[0;Ω]⊗q
     12. 打包 dZ = [V; dq/dt; dθ; dφ; dψ; ...]
```

求解器：`ode15s`，配合事件驱动的分段分裂（超时事件和奇异性事件）。

---

## 16. 胞体布朗动力学（分裂步 Euler-Maruyama）

### 16.1 方法概述

胞体布朗运动通过**算符分裂（Operator Splitting）**与确定性弹性流体动力学解耦：

- **Step A（确定性部分）**：`ode15s` 完整求解耦合胞体-纤毛弹性流体动力学；RHS 函数 `ode_coupled_rhs.m` 无需任何修改
- **Step B（随机部分）**：每隔固定的布朗时间步 `dt_B`，向胞体状态施加 Euler-Maruyama 布朗位移"踢脚"

该策略避免了将随机噪声注入隐式自适应求解器 `ode15s` 所导致的 Jacobian 污染和误差控制失效问题。分裂步引入的额外误差为 $O(\text{dt}_B)$，与 Euler-Maruyama 本身弱一阶精度一致，不降低整体精度。

### 16.2 物理模型

过阻尼 Langevin 方程（Stokes 流极限）：

$$
\boldsymbol{\xi}_T \dot{\mathbf{X}} = \mathbf{F}_{\text{det}} + \sqrt{2k_BT\,\boldsymbol{\xi}_T}\;\boldsymbol{\eta}_T(t)
$$

$$
\boldsymbol{\xi}_R \dot{\boldsymbol{\Theta}} = \mathbf{T}_{\text{det}} + \sqrt{2k_BT\,\boldsymbol{\xi}_R}\;\boldsymbol{\eta}_R(t)
$$

其中 $\boldsymbol{\xi}_T$、$\boldsymbol{\xi}_R$ 分别为平移和旋转摩擦张量，$\boldsymbol{\eta}$ 为高斯白噪声，$k_B$ 为 Boltzmann 常数，$T$ 为绝对温度。

### 16.3 摩擦矩阵与扩散张量

`compute_body_brownian_params.m` 计算胞体 6×6 摩擦矩阵并提取 Cholesky 扩散因子：

1. **RSM 迁移率矩阵** $M_{BB} \in \mathbb{R}^{3N_b \times 3N_b}$：胞体表面点之间的正则化 Stokeslet 核函数
2. **无滑移约束矩阵** $N_{\text{ns}} \in \mathbb{R}^{3N_b \times 6}$：将 6 个刚体自由度映射为 $3N_b$ 个表面点速度
3. **力/力矩提取矩阵** $E \in \mathbb{R}^{6 \times 3N_b}$：将分布力求和为总力和总力矩
4. **摩擦矩阵**：$\Xi = E \cdot M_{BB}^{-1} \cdot N_{\text{ns}} \in \mathbb{R}^{6 \times 6}$
5. **Schur 补解耦**：将 $\Xi$ 的平移—旋转耦合消除，得到解耦的 $\boldsymbol{\xi}_T$ 和 $\boldsymbol{\xi}_R$
6. **扩散张量**：$D_T = k_BT\,\boldsymbol{\xi}_T^{-1}$，$D_R = k_BT\,\boldsymbol{\xi}_R^{-1}$
7. **Cholesky 分解**：$2D_T = L_T L_T^T$，$2D_R = L_R L_R^T$

输出结构体 `params.brown` 包含：`L_trans`、`L_rot`（Cholesky 因子）、`D_trans`、`D_rot`（扩散张量）、`FM`（完整摩擦矩阵）、`DFM_trans`、`DFM_rot`（解耦摩擦张量）。

### 16.4 布朗扰动施加

在每个布朗节点 $t = n \cdot \text{dt}_B$ 处：

**平移扰动**（体坐标系 → 实验室坐标系）：

$$
\Delta\mathbf{X}_{\text{lab}} = R_{\text{body}} \cdot L_T \cdot \boldsymbol{\eta}_T \cdot \sqrt{\text{dt}_B}
$$

**旋转扰动**（四元数合成）：

$$
\Delta\boldsymbol{\Theta}_{\text{lab}} = R_{\text{body}} \cdot L_R \cdot \boldsymbol{\eta}_R \cdot \sqrt{\text{dt}_B}
$$

$$
\Delta q = \left[\cos\frac{|\Delta\boldsymbol{\Theta}|}{2};\; \frac{\sin(|\Delta\boldsymbol{\Theta}|/2)}{|\Delta\boldsymbol{\Theta}|}\,\Delta\boldsymbol{\Theta}_{\text{lab}}\right]
$$

$$
q_{\text{new}} = \Delta q \otimes q_{\text{old}} \quad (\text{归一化})
$$

每次扰动后重新施加固支几何投影（`enforce_clamped_geometry`）以修正约束漂移。

### 16.5 主循环积分策略

当布朗运动启用时，外层 `while` 循环的积分目标从 `T_obs_dimless` 变为 `min(T_brown_next, T_obs_dimless)`，其中 `T_brown_next` 为下一个布朗节点（以 `dt_B / T_0` 为间隔的等距序列）。每当 `ode15s` 到达布朗节点时施加扰动；若被奇点事件提前终止则不施加扰动，继续积分直到节点。

### 16.6 用户配置

| 参数           | 默认值   | 说明                                     |
| -------------- | -------- | ---------------------------------------- |
| `enable_brown` | `false`  | 启用/禁用胞体布朗运动                   |
| `Temper`       | `293.15` | 流体绝对温度 [K]（默认 20°C）           |
| `dt_B`         | `dt/2`   | 布朗时间步 [s]（应远小于最短弛豫时间）   |

**启用布朗运动时**，`init_coupled_params` 会自动调用 `compute_body_brownian_params` 计算扩散参数，并将结果存入 `params.brown` 结构体。主循环在积分过程中自动检测弛豫时间并发出步长合理性警告。

---

## 17. 输出系统：日志、数据导出与动画

### 17.1 输出目录

- 非空 `case_name` → `output/<case_name>/`
- 为空 → `output/run_YYYYMMDD_HHMMSS/`

**覆盖保护**（v3.1）：当输出目录已包含仿真数据文件（`.mat`、`.csv`、`.log`、`.mp4`、`.png`）时，弹出三选一对话框：

| 选项 | 行为 |
| ---- | ---- |
| Overwrite | 清空目录内所有文件，重新运行 |
| Rename | 为当前运行追加时间戳创建新目录（`<name>_YYYYMMDD_HHmmss`），原目录保留不变 |
| Abort | 终止仿真，返回空结果 |

`auto_confirm = true`（无人值守模式）下默认执行 Rename 操作。

### 17.2 输出文件

| 文件                     | 格式            | 内容                                                  |
| ------------------------ | --------------- | ----------------------------------------------------- |
| `simulation_results.mat` | MAT v7.3 (HDF5) | `results` + `params` 完整结构体                                                      |
| `timeseries.csv`         | 纯文本          | 10 列标量时间序列（见下表）                                                          |
| `config.json`            | JSON            | 配置参数快照（可直接作为输入复用以精确复现运行）                                     |
| `summary_plot.png`       | PNG 600 dpi     | 静态总览：CM 轨迹、沉降速度、倾斜角、快照三视图                                      |
| `initial_pose.png`       | PNG 600 dpi     | 初始姿态预览三视图（`preview_initial = true` 时生成）                                |
| `bacteria_3d_motion.mp4` | MPEG-4 质量 100 | 3D 动画，三面板 1×3 布局（3D + XZ + YZ 投影），动态图像尺寸，含光滑胞体表面 + 着色纤毛 + CM 轨迹 |
| `flowfield_panels.png`   | PNG 600 dpi     | 双平面流场 2×2 面板图：XZ + YZ，3D 表面视图 + 2D 着色图；jet 色图 + 归一化箭头       |
| `simulation.log`         | 纯文本          | 带时间戳的完整双通道运行日志                                                         |

### 17.3 results 结构体字段（物理单位）

| 字段         | 单位 | 尺寸           | 说明                                    |
| ------------ | ---- | -------------- | --------------------------------------- |
| `t`          | s    | `[1 × NT]`     | 物理时间序列                            |
| `t_dimless`  | —    | `[1 × NT]`     | 无量纲时间序列                          |
| `X_cm_um`    | μm   | `[3 × NT]`     | 质心位置 (X/Y/Z)                        |
| `V_cm_um`    | μm/s | `[3 × NT]`     | 质心速度三分量                          |
| `V_sink`     | μm/s | `[1 × NT]`     | 沉降速度（Z 分量）                      |
| `tilt_angle` | rad  | `[1 × NT]`     | 胞体长轴偏仰角（与 `alpha_0` 约定一致） |
| `q`          | —    | `[4 × NT]`     | 方向四元数                              |
| `X_pili_um`  | μm   | cell `{M × 1}` | 各纤毛节点坐标 `[Nseg+1 × 3 × NT]`      |
| `F_sink_pN`  | pN   | 标量           | 净沉降力（重力浮力差）                  |
| `units`      | —    | struct         | 单位及特征尺度元数据                    |

**`timeseries.csv` 列说明：**

| 列名                            | 单位 | 说明               |
| ------------------------------- | ---- | ------------------ |
| `t_s`                           | s    | 时间               |
| `X_um`, `Y_um`, `Z_um`          | μm   | 质心坐标           |
| `Vx_um_s`, `Vy_um_s`, `Vz_um_s` | μm/s | 质心速度三分量     |
| `Vsink_um_s`                    | μm/s | 沉降速度（Z 方向） |
| `tilt_rad`                      | rad  | 胞体倾斜角         |
| `tilt_deg`                      | °    | 胞体倾斜角（度）   |

### 17.4 单位约定

| 物理量 | 输出单位       |
| ------ | -------------- |
| 时间   | 秒 [s]         |
| 长度   | 微米 [μm]      |
| 速度   | 微米/秒 [μm/s] |
| 力     | 皮牛 [pN]      |
| 角度   | 弧度 [rad]     |

---

## 18. 流场可视化

`compute_plot_flowfield.m` 在仿真末状态计算并绘制 XZ 平面（Y = Y_cm）和 YZ 平面（X = X_cm）的流场图像，以 2×2 面板布局呈现。

### 18.1 计算流程

1. **状态还原**：调用 `reconstruct_geometry` 重建末状态几何；调用 `ode_coupled_rhs` 获取约束速率 `dZ`
2. **运动学提取**：从 `dZ` 提取胞体平移/旋转速度 (V, Ω) 及纤毛角速率 (dθ, dφ, dψ)
3. **节点速度**：利用刚体公式和运动链计算所有胞体/纤毛节点的平移速度 `v_trans` 与扭转速率 `v_twist`
4. **力系数求解**：重建迁移率矩阵 `M_weighted`、Rotlet 矩阵 `M_rotlet`、扭转阻力 `A_twist`；计算扭矩密度 `tau`；对 `v_trans` 施加 Rotlet 校正后求解 `f_coeffs = M_weighted \ v_trans_corr`
5. **双平面网格评估**：分别在 XZ 和 YZ 平面上，以胞体质心为中心的 10×10 μm 网格（间距 0.05 μm）上，利用 RSM（胞体）、RSS（纤毛）和 Rotlet（纤毛扭矩）核函数向量化评估流体速度，采用共享的局部函数 `evaluate_plane_velocity`，分块策略（`chunk_size = 2000`）控制内存
6. **胞体遮罩**：将各网格点变换至胞体局部坐标系，计算其到胶囊体中轴线的距离，标记位于胞体内部的点不参与绘制

### 18.2 可视化内容

2×2 面板布局：

| 位置 | 视图       | 说明                                                                                         |
| ---- | ---------- | -------------------------------------------------------------------------------------------- |
| 左上 | 3D — XZ 面 | XZ 流场以纹理表面（`texturemap`）铺设于 3D 视图中，半透明胞体 + 纤毛叠加                     |
| 右上 | 3D — YZ 面 | YZ 流场以纹理表面铺设于 3D 视图中，半透明胞体 + 纤毛叠加                                     |
| 左下 | 2D — XZ 面 | 速度大小 `imagesc` + 归一化等长 `quiver` 箭头（白色光晕 + 深色核心）；胞体白色填充 + 灰色描边 |
| 右下 | 2D — YZ 面 | 速度大小 `imagesc` + 归一化等长 `quiver` 箭头（白色光晕 + 深色核心）；胞体白色填充 + 灰色描边 |

**共享特性**：

- `jet` 色图，统一色标（99 百分位数截断）
- 胞体内部遮罩（胶囊体 SDF 距离场 + `contourc` 轮廓提取）
- 纤毛穿越可视化平面处以圆圈标记

### 18.3 配置参数

流场参数已提升至主脚本 `cfg` 结构体，通过 `flow_opts` 传入 `compute_plot_flowfield`。缺省值在函数内部通过 `get_flow_opt` 回退提供：

| cfg 字段 | flow_opts 字段 | 默认值 | 说明 |
| -------- | -------------- | ------ | ---- |
| `flowfield_grid_spacing` | `grid_spacing_um` | 0.05 | 网格间距 [μm] |
| `flowfield_grid_range` | `grid_range_um` | 10 | 网格总范围 [μm]（10×10） |
| `flowfield_quiver_skip` | `quiver_skip` | 5 | 箭头降采样倍率（每 5 个点一箭） |
| `flowfield_chunk_size` | `chunk_size` | 2000 | 核函数评估分块大小 |

流场计算过程现配备 `render_progress_bar` 实时进度条（含 ETA），分别显示 `[FLOW-XZ]` 和 `[FLOW-YZ]` 标签。

---

## 19. 文件结构

```
Coupled_Body_Pili/
│
├── main_coupled_simulation.m           % 入口脚本：cfg 结构体定义 + 调用 run_simulation
├── run_simulation.m                    % 仿真引擎：积分循环、布朗积分、后处理、输出
│
├── load_sim_config.m                   % JSON 配置文件 → cfg 结构体
├── save_sim_config.m                   % cfg 结构体 → JSON 配置快照
├── parse_pili_config.m                 % 纤毛命名结构体 → Pili_Matrix [8×M]
├── validate_simulation_params.m        % 参数合理性预检（错误聚合 + 警告分级）
│
├── init_coupled_params.m               % 参数初始化（无量纲化、几何、布朗参数）
├── init_coupled_state.m                % 初始状态向量 Z0
│
├── ode_coupled_rhs.m                   % 核心 ODE 右端项（Q/A/B 组装 + 线性求解 + Rotlet 耦合）
├── reconstruct_geometry.m              % 状态向量 → 实验室坐标系坐标
├── build_global_resistance.m           % RSS 迁移率 + Rotlet 线段 + 扭转 RFT
├── compute_pili_elastic_moments.m      % Kirchhoff 弹性杆弹性曲率偏差
│
├── enforce_clamped_geometry.m          % 固支几何投影（约束漂移修正）
├── adaptive_basis_check.m              % 逐纤毛自适应 Euler 角基底旋转
│
├── compute_body_brownian_params.m      % 胞体布朗参数：RSM 摩擦矩阵 → Schur 补 → Cholesky 扩散因子
│
├── capsule_body_building.m             % 胶囊体表面网格生成
├── pili_building.m                     % 纤毛结构初始化（刚度、曲率）
├── pili_type_building.m                % 纤毛形态生成（5 种类型）
├── compute_curvature_binormal.m        % 离散弹性杆曲率副法线 (Bergou 2008)
│
├── quat_to_rotmat.m                    % 四元数 → 旋转矩阵
├── quat_multiply.m                     % Hamilton 四元数乘法
├── quat_kinematics.m                   % dq/dt = ½[0;Ω]⊗q
├── skew_matrix.m                       % 3×3 反对称叉积矩阵
│
├── pili_directors.m                    % Euler 角 → 标架 (d1, d2, d3)
├── pili_euler_angles.m                 % 标架 → Euler 角（逆运算）
├── pili_spatial_coords.m               % 根部 + Euler 角 → 节点位置
├── unpack_state.m                      % Z → 分量提取
├── pack_state_rates.m                  % 分量 → dZ 打包
│
├── compute_plot_flowfield.m            % 末状态双平面流场评估与绘图（2×2 面板）
│
├── odetpbar.m                          % ODE 段内进度条（ode15s OutputFcn）
├── render_progress_bar.m               % 通用文本进度条（动画渲染/流场计算）
│
├── batch sweep/                        % 自动化参数扫描与批量运行脚本集
│   ├── batch_sweep_alpha.m             % 扫描初始倾斜角 alpha_0（0°~90°，7 工况）
│   ├── batch_sweep_major_axis.m        % 扫描胞体长半轴（0.75~2.00 μm，6 工况）
│   ├── batch_sweep_pili_length.m       % 扫描纤毛长度（1~6 μm，6 工况）
│   ├── batch_sweep_pili_site.m         % 扫描纤毛附着点位（1~14，14 工况）
│   ├── batch_sweep_pili_count.m        % 扫描纤毛数量（1~6 根，6 工况）
│   ├── batch_sweep_pili_E.m            % 扫描纤毛杨氏模量（1e6~1e9 Pa，5 工况）
│   └── batch_run_from_json.m           % JSON 批量运行器：检索文件夹→依次执行→失败跳过
│
├── configs/                            % JSON 配置文件目录
│   └── OnePili-P1-L3.json             % 单根纤毛、P1 附着点、3 μm 长度的参考工况配置
│
├── config_example.json                 % JSON 配置文件模板
│
├── README.md                           % 本文件
├── Brownian_Implementation_Plan.md     % 布朗动力学实施方案文档
└── Cortez - 2018 - Regularized stokeslet segments.pdf
```

总计：**38 个 MATLAB 源文件 + 1 个 JSON 配置模板 + 1 个 JSON 参考工况配置**（31 个核心模块 + 7 个批量扫描脚本，无冗余代码，无未使用函数）。

---

## 20. 与 Walker (2020) 的标度一致性

无量纲化体系经设计可在简化为无胞体的单根纤毛时精确复现 Walker 的框架。关键等价关系：

- $Q$ 矩阵：结构相同（基底速度 + Euler 角速率）
- $B$ 矩阵：力平衡使用梯形求积权重；力矩平衡使用 $\ell^3$ 标度
- $\tilde{A}_{\text{twist}}$：包含额外的 $\ell_k$ 因子以确保量纲一致性
- 弹性力矩：曲率偏差公式与 Walker 的离散算子一致

---

## 21. 多纤毛支持

### 21.1 配置

`Pili_Matrix` 的每一列定义一根纤毛。系统支持 $M \geq 1$ 根纤毛，各自独立设定：

- 长度、形状类型、附着位点
- Euler 角状态和自适应基底
- 弹性属性（共享材料参数，独立的固有曲率）

### 21.2 耦合

| 机制              | 实现方式                                                                                    |
| ----------------- | ------------------------------------------------------------------------------------------- |
| 纤毛—纤毛流体力学 | 全局迁移率矩阵中的完全交叉耦合（$M_{PP}$ 块）                                               |
| 胞体—纤毛流体力学 | $M_{BP}$ 和 $M_{PB}$ 块中的交叉耦合                                                         |
| 弹性力矩          | 逐纤毛独立计算（Kirchhoff 曲率）                                                            |
| 固支边界条件      | 逐纤毛独立约束（不同附着位点）                                                              |
| 扭转阻力          | 块对角 $\tilde{A} = \text{blkdiag}(\tilde{A}_1, \ldots, \tilde{A}_M)$（局部 RFT，物理正确） |
| 自适应基底        | 逐纤毛独立检测与旋转                                                                        |

### 21.3 附着位点验证

`init_coupled_params.m` 检查 `Pili_Matrix(8,:)` 中是否存在重复附着位点，若检测到则发出警告。

### 21.4 矩阵维度汇总

| 矩阵                         | 维度                                                      |
| ---------------------------- | --------------------------------------------------------- |
| $\mathbf{Z}$                 | $7 + 3\sum_k N_k$                                         |
| $\dot{\mathbf{z}}$           | $6 + 3\sum_k N_k$                                         |
| $M_{\text{weighted}}$        | $3N_{\text{total}} \times 3N_{\text{total}}$              |
| $M_{\text{rotlet}}$          | $3N_{\text{total}} \times 3\sum_k N_k$                    |
| $\mathbf{Q}_{\text{global}}$ | $(3N_{\text{total}} + \sum_k N_k) \times n_{\text{dof}}$  |
| $\mathbf{B}_{\text{global}}$ | $n_{\text{dof}} \times (3N_{\text{total}} + 3\sum_k N_k)$ |

---

## 22. 验证总结

### 22.1 已验证项目 (30+)

| 类别                | 详情                                                                         |
| ------------------- | ---------------------------------------------------------------------------- |
| RSM 胞体核函数      | $1/(8\pi\mu)$ 吸收、特征尺度自洽性、胶囊体积                                 |
| 四元数代数          | 左乘约定、运动学 $\dot{q} = \frac{1}{2}[0;\Omega]\otimes q$、归一化          |
| Q 矩阵              | 标度、坐标变换、根部耦合、交错排序                                           |
| B 矩阵              | 力的梯形权重、力矩平衡、$\ell^3$ 标度                                        |
| RSS $T_{n,q}$ 递推  | 全部 6 项（$T_{0,-1}$ 至 $T_{3,-3}$）已对照 Cortez 2018 Eq. 7 验证           |
| RSS $M_1$/$M_2$ 块  | Eq. 12–13，对称性检查（机器精度），自相互作用安全                            |
| RSS $T_{n,-5}$ 递推 | Rotlet 递推、$d^2$ 安全措施、Sec. 2.4 公式                                   |
| Rotlet 反对称块     | 反对称性已验证（$\|R+R^T\|/\|R\| = 0$），收敛率 2.0                          |
| 向量化组装          | 所有 $[N_{\text{total}} \times N_{\text{seg,all}}]$ 数组与标量核函数输出一致 |
| 扭转阻力            | $\ell_k \varepsilon_k^2 / 2$ 标度，量纲一致                                  |
| 固支边界条件        | 雅可比矩阵推导、速率层面 + 位置层面双重施加                                  |
| 自适应基底          | 逐纤毛检测与旋转、Euler 角重新计算                                           |
| 多纤毛              | 状态向量、几何、Q/B 矩阵、边界条件均已验证（$M \geq 1$）                     |
| 线性系统            | $L\dot{z}=R$ 每行逐一量纲追踪                                                |
| 布朗动力学          | RSM 摩擦矩阵、Schur 补解耦、Cholesky 正定性检查、FDT 一致性                  |

### 22.2 缺陷修复历史（均已解决）

| #   | 问题                                    | 日期       | 修复                                     |
| --- | --------------------------------------- | ---------- | ---------------------------------------- |
| 1   | 纤毛端点列权重不一致                    | 2026-02-07 | 统一 $\Delta s$（Walker 一致性）         |
| 2   | `inv()` 数值稳定性                      | 2026-02-07 | 替换为 `\` 运算符                        |
| 3   | 全矩阵求逆性能                          | 2026-02-10 | 返回 $M_{\text{weighted}}$，使用 LU 求解 |
| 4   | $\tilde{A}_{\text{twist}}$ 坐标系不匹配 | 2026-02-10 | 使用实验室坐标系 `d3`                    |
| 5   | PascalCase 文件命名                     | 2026-02-11 | 重命名为 `snake_case`                    |
| 6   | 几何构建器中未使用的力点数据            | 2026-02-11 | 移除冗余代码                             |
| 7   | 胞体以散点渲染                          | 2026-02-11 | 参数化 `surf` + Gouraud 着色             |
| 8   | 扭转阻力缺少 $\ell_k$ 因子              | 2026-02-12 | `eps_k^2/2` → `ell_k*eps_k^2/2`          |
| 9   | `reconstruct_geometry` 冗余调用         | 2026-02-12 | 移至循环外                               |
| 10  | 多纤毛统一着色                          | 2026-02-12 | 10 色可区分调色板                        |
| 11  | 无附着位点冲突检测                      | 2026-02-12 | 添加重复位点警告                         |
| 12  | `acos` 因浮点域错误产生 NaN             | 2026-02-12 | 钳位至 $[-1, 1]$                         |
| 13  | 四元数初始化假设胞体 Z 为长轴           | 2026-02-12 | 修正为胞体 X 轴                          |
| 14  | 移除纯 RSM 纤毛模式 (v2.0)              | 2026-02-14 | 统一为 RSS + Rotlet；删除 `use_rss` 标志 |
| 15  | 纤毛线宽与质心标记过大 (v2.1)           | 2026-02-15 | 静态图与动画中减半                       |
| 16  | 坐标轴范围无最小下限 (v2.1)             | 2026-02-15 | 增加 `4*major_axis` 最小范围约束         |
| 17  | 动画中 CM 轨迹断续 (v2.1)               | 2026-02-15 | 使用全分辨率轨迹数据替代帧采样稀疏数据   |
| 18  | 重力标识偏离坐标轴框 (v2.1)             | 2026-02-15 | 缩小间距；添加 `drawnow` 确保位置准确    |
| 19  | `TightInset(3)` 对三维坐标轴高估，重力标识与轴框间距仍然过大 (v2.3→v2.4) | 2026-02-19 | 将 `TightInset(3)` 截断至 0.020（图宽 2%）；三维摄像机投影会把 tick/axis 标签投影到 Position 框外，MATLAB 将其全部计入 TightInset 导致虚高 |
| 20  | 动画中重力标识在坐标轴全属性配置前放置，`Position` 在第一帧后收缩导致标识垂直偏移 (v2.4) | 2026-02-19 | 在 `add_gravity_indicator_outside` 调用前预配置 `xlim/ylim/zlim`、`DataAspectRatio`、`view`、`xlabel/ylabel/zlabel` 及 `title`（字体大小一致）；`cla` 不清除 `HandleVisibility='off'` 的装饰对象，预配置状态全程稳定 |
| 21  | 根本原因：`DataAspectRatio=[1 1 1]` 等比例约束使极细长轴框（Z>>X,Y）居中于 Position 矩形，实际轴框右边缘与 Position 右边缘相差数百像素；之前所有基于 Position 右边缘的方案均无法修复此问题 (v2.4) | 2026-02-19 | 彻底重写 `add_gravity_indicator_outside`：通过正射影投影（`right=[cos az, sin az, 0]`）将数据边界框投影到屏幕水平方向，计算等比例缩放后的真实渲染宽度，以居中关系求得轴框真实右边界，完全移除对 `TightInset` 的依赖 |
| 22  | 静态快照中胞体与纤毛颜色偏暗，无法呈现预设明亮渐变色 (v2.4→v2.5) | 2026-02-19 | 移除三个快照视图的 `camlight`/`lighting('gouraud')`/`material('dull')`，改为 `lighting(axi,'none')` 均匀环境光，曲面以 `FaceColor` 原色全亮度渲染，消除单向定向光照造成的暗区 |
| 23  | 3D 动画中 CM 轨迹线遮挡半透明胞体和纤毛 (v2.4→v2.5) | 2026-02-19 | 调整绘制顺序（胞体→纤毛→CM 标记→轨迹），绘制轨迹后立即调用 `uistack(h_traj,'bottom')` 将其移至 `Children` 数组首位，强制 Painters 渲染器最先绘制轨迹，使其位于所有对象之下 |
| 24  | 胞体正则化参数 `eps_body` 由 `h_min/2` 改为 `h_mean/6` (v2.6) | 2026-03-01 | 原策略受球帽—圆柱拼接处局部几何（CDis/2 ≈ 38.6 nm）驱动，使 ε/Δs 系统性偏大（≈ 0.24–0.35）；新策略基于全局平均最近邻间距，ε/Δs ≈ 0.167，与 Cortez et al. (2005) 最优比值 0.154 吻合，正则化误差降低约 40% |
| 25  | 胞体正则化参数 `eps_body` 由 `h_mean/6` 调整为 `h_mean/3` (v2.7) | 2026-03-01 | 经数值测试验证，ε/Δs ≈ 0.333 对当前胶囊几何（$a/b=3.125$，$N_\text{body}=1000$）为经验最优值，阻力矩阵精度与数值稳定性均优于 h_mean/6 |
| 26  | 动画 CM 轨迹在 OpenGL 渲染器下遮挡纤毛 (v2.8) | 2026-03-03 | `uistack('bottom')` 仅对 Painters 渲染器有效；动画中因 `FaceAlpha` 曲面强制启用 OpenGL，渲染顺序由深度缓冲（z-buffer）决定，轨迹线在几何位置靠近相机时必然遮挡纤毛。解决方案：将 `traj_color` 由 `[R G B]` 改为 4 分量 `[R G B 0.30]`（RGBA，alpha = 0.30），即使几何位置在前也不完全遮挡，同时移除无效的 `uistack` 调用 |
| 27  | 动画视频左右两侧大面积留白 (v2.8) | 2026-03-03 | `DataAspectRatio=[1 1 1]` 将等比例数据框居中于 Position 矩形，细菌沉降时 Z 范围远大于 X/Y，投影内容呈竖向细长形，固定 1080×810 横向图像导致两侧大量空白。解决方案：按正射影投影（az=45°, el=25°）将数据边界框投影到屏幕，推算内容宽高比，动态设定图像尺寸（固定高度 1440 px，宽度 = `max(1440 × proj_w/proj_h × fill_ratio, 900)`），视频分辨率随数据自适应，留白消除 |
| 28  | `validate_simulation_params.m` 重复附着位点检测索引错误 (v3.1→v3.2) | 2026-03-17 | `sites(accumarray(ic,1)>1)` 以逻辑向量（长度=唯一位点数）索引 `sites`（长度=总纤毛数），维度不匹配导致崩溃；改为 `u_sites(counts > 1)` 正确过滤唯一位点 |
| 29  | `run_simulation.m` "覆盖" 选项未清除旧文件 (v3.1→v3.2) | 2026-03-17 | 输出目录覆盖对话框的 "Overwrite" 分支仅记录日志但未删除已有文件，旧输出残留；添加 `delete` 循环显式清除所有已有文件 |
| 30  | `load_sim_config.m` `get_json` 死代码 (v3.1→v3.2) | 2026-03-17 | `isnumeric(default) && islogical(default)` 恒为 `false`（`islogical` 为 `true` 时 `isnumeric` 必为 `false`），`logical(v)` 转换永远不执行；移除 `isnumeric(default) &&` |
| 31  | `main_coupled_simulation.m` 案例名拼写错误 (v3.1→v3.2) | 2026-03-17 | `'Brwon-Test-Soft'` → `'Brown-Test-Soft'` |
| 32  | `render_progress_bar.m` 进度条无法原地更新 (v3.1→v3.2) | 2026-03-17 | `prev_len_container` 数值数组按值传递，函数内修改不回传至调用方，`\b` 退格数始终为 0，进度信息逐行堆叠而非原地刷新；改用 `persistent prev_msg_len` 管理状态（与 `odetpbar.m` 同模式），`current==1` 时重置 |

### 22.3 功能历史

| #   | 功能                                                              | 日期       |
| --- | ----------------------------------------------------------------- | ---------- |
| 1   | 初始姿态参数化 `alpha_0`                                          | 2026-02-12 |
| 2   | 案例命名 `case_name`                                              | 2026-02-12 |
| 3   | 专业进度条 `odetpbar.m`                                           | 2026-02-12 |
| 4   | 多纤毛可区分着色                                                  | 2026-02-12 |
| 5   | 附着位点冲突检测                                                  | 2026-02-12 |
| 6   | ODE 事件触发日志                                                  | 2026-02-12 |
| 7   | 初始姿态交互式预览                                                | 2026-02-12 |
| 8   | 纤毛 RSS 流体力学 (Cortez 2018)                                   | 2026-02-14 |
| 9   | 向量化 RSS 组装（消除线段循环）                                   | 2026-02-14 |
| 10  | Rotlet 线段迁移率 $M_{\text{rotlet}}$                             | 2026-02-14 |
| 11  | 纯 RSS 模式——移除 RSM 纤毛回退 (v2.0)                             | 2026-02-14 |
| 12  | 纤毛线宽与质心标记减半 (v2.1)                                     | 2026-02-15 |
| 13  | 坐标轴范围最小下限 `4*major_axis` (v2.1)                          | 2026-02-15 |
| 14  | 动画中全分辨率 CM 轨迹 (v2.1)                                     | 2026-02-15 |
| 15  | 重力标识紧贴坐标轴框 (v2.1)                                       | 2026-02-15 |
| 16  | `results.tilt_angle`：胞体长轴偏仰角时间序列 (v2.2)               | 2026-02-18 |
| 17  | `results.V_cm_um`：质心速度三分量时间序列 (v2.2)                  | 2026-02-18 |
| 18  | 初始姿态预览图自动保存为 `initial_pose.png` (v2.2)                | 2026-02-18 |
| 19  | `timeseries.csv` 导出：10 列标量时间序列，便于外部工具调用 (v2.2) | 2026-02-18 |
| 20  | summary_plot 第 3 子图替换为胞体倾斜角时间演化图 (v2.2)           | 2026-02-18 |
| 21  | 修复重力标识位置偏移：`sgtitle` 后调用以获取最终 Position (v2.3)  | 2026-02-18 |
| 22  | 重力标识定位改用 `TightInset(3)` 精确计算内容右边界 (v2.3)        | 2026-02-18 |
| 23  | 根治重力标识右偏：`TightInset(3)` 截断至 0.020，消除三维轴标签余量虚高 (v2.4) | 2026-02-19 |
| 24  | 动画轴完整预配置（`xlim/ylim/zlim` + `DataAspectRatio` + `view` + 轴标签 + `title`）后再放置重力标识，`Position` 全帧稳定 (v2.4) | 2026-02-19 |
| 25  | 彻底根治：正射影投影法计算 `DataAspectRatio=[1 1 1]` 等比例约束下的轴框真实渲染宽度，完全移除对 `TightInset`/`Position` 右边缘的依赖；对极细长数据范围（细菌沉降 Z>>X,Y）准确有效 (v2.4) | 2026-02-19 |
| 26  | 静态快照均匀光照：三视图改用 `lighting('none')` 环境光，消除 Gouraud 单向光照对非正对摄像机曲面的暗化，胞体/纤毛以真实颜色渐变完整呈现 (v2.5) | 2026-02-19 |
| 27  | 动画 CM 轨迹置底层渲染：绘制完成后调用 `uistack(h_traj,'bottom')`，轨迹线始终位于胞体与纤毛之下，不影响半透明胞体及纤毛的视觉效果 (v2.5) | 2026-02-19 |
| 28  | 胞体正则化参数优化为 `h_mean/6`（v2.6）：正则化参数策略由局部最小间距 `h_min/2` 改为全局平均最近邻间距 `h_mean/6`，ε/Δs ≈ 0.167，符合 Cortez et al. (2005) 对球体算例的最优配置，兼顾正则化精度与数值稳定性 | 2026-03-01 |
| 29  | 胞体正则化参数进一步调整为 `h_mean/3`（v2.7）：经数值测试，ε/Δs ≈ 0.333 对当前胶囊几何为经验最优，精度与稳定性优于 `h_mean/6` | 2026-03-01 |
| 30  | PNG 输出分辨率升级至 600 dpi（v2.8）：`initial_pose.png` 由 200 dpi 升至 600 dpi；`summary_plot.png` 由 300 dpi 升至 600 dpi，满足期刊投稿图像质量要求 | 2026-03-03 |
| 31  | 坐标轴最小范围由 4 倍半长轴（5 μm）改为 4 倍胞体全长（10 μm，= `4 × 2 × major_axis`）（v2.8）：确保三视图在细菌几何尺度下具有充足的空间背景；同时为 `initial_pose.png` 新增最小范围约束（此前无约束），三图规则统一 | 2026-03-03 |
| 32  | 动画视频分辨率与质量提升（v2.8）：图像高度由 900 px 提升至 1440 px（2K 级纵向分辨率），最小宽度由 560 px 提升至 900 px；MPEG-4 压缩质量由 95 升至最高值 100，视频更清晰锐利 | 2026-03-03 |
| 33  | XZ 平面流场可视化（v2.9）：新增 `compute_flowfield_xz.m`，在仿真末状态计算并绘制 XZ 平面流场。利用 RSM/RSS/Rotlet 核函数在 10×10 μm 网格（间距 0.1 μm）上向量化评估流体速度，jet 色图背景着色速度大小 + 归一化等长 quiver 箭头表示流向；胞体内部由球冠柱距离场遮罩；分块策略控制内存；输出 `flowfield_xz.png`（600 dpi） | 2026-03-06 |
| 34  | 胞体布朗动力学（v3.0）：新增分裂步算符分裂（Operator Splitting）框架，`ode15s` 确定性求解器与固定间隔 `dt_B` Euler-Maruyama 布朗扰动解耦。RSM 摩擦矩阵 → Schur 补解耦 → Cholesky 扩散因子；四元数旋转布朗扰动；布朗节点后自动重施加固支几何投影 | 2026-03-14 |
| 35  | `compute_body_brownian_params.m` 新建（v3.0）：基于 RSM 迁移率矩阵构建胶囊体 6×6 刚体摩擦矩阵，Schur 补解耦平移/旋转子块，Cholesky 分解生成 Euler-Maruyama 所需的 $L_T$、$L_R$ 因子 | 2026-03-14 |
| 36  | 双平面流场 2×2 面板可视化（v3.0）：`compute_flowfield_xz.m` 重构为 `compute_plot_flowfield.m`，同时计算 XZ 和 YZ 平面流场；2×2 面板布局（3D 表面视图 + 2D imagesc，各含 XZ 和 YZ）；统一色标（99 百分位截断）；`contourc` 轮廓提取胞体遮罩轮廓；输出 `flowfield_panels.png` | 2026-03-14 |
| 37  | 流场网格间距由 0.1 μm 细化至 0.05 μm（v3.0）：分辨率翻倍，流场细节更清晰 | 2026-03-14 |
| 38  | 动画三面板 1×3 布局（v3.0）：3D 视图（左）+ XZ 投影（中）+ YZ 投影（右），每帧三个面板同步更新；XZ/YZ 面板采用正射影投影（`Projection = orthographic`），手动补充封闭轴框线 | 2026-03-14 |
| 39  | 动画图像尺寸自适应（v3.0）：固定高度 1200 px，宽度按正射影投影数据边界框宽高比动态计算（最小 1800 px），适配 1×3 三面板布局 | 2026-03-14 |
| 40  | 流场参数提升至主脚本（v3.1）：`grid_spacing_um`、`grid_range_um`、`quiver_skip`、`chunk_size` 从 `compute_plot_flowfield` 硬编码提升至 `cfg` 结构体，通过 `flow_opts` 传入，函数内通过 `get_flow_opt` 回退默认值 | 2026-03-17 |
| 41  | 纤毛结构体接口（v3.1）：新增 `parse_pili_config.m`，支持以命名字段（`type`, `length_um`, `site` 等）定义纤毛配置，自动转换为 `Pili_Matrix [8×M]`；保留原始矩阵格式向后兼容 | 2026-03-17 |
| 42  | 无人值守模式（v3.1）：新增 `auto_confirm` 参数（默认 `false`），设为 `true` 时跳过初始姿态预览确认和输出目录覆盖对话框，支持无交互批量运行 | 2026-03-17 |
| 43  | 参数合理性预检（v3.1）：新增 `validate_simulation_params.m`，在仿真启动前对所有用户参数执行合理性检查；错误聚合（一次性报告全部违规）+ 警告分级（次优但可接受的参数值） | 2026-03-17 |
| 44  | 输出目录覆盖保护（v3.1）：检测到输出目录已包含仿真数据文件时，弹出三选一对话框（覆盖/重命名/中止）；`auto_confirm` 模式下默认重命名（追加时间戳） | 2026-03-17 |
| 45  | 通用文本进度条（v3.1→v3.2 修复）：新增 `render_progress_bar.m`，为动画渲染（`[ANIM]`）和流场计算（`[FLOW-XZ]`/`[FLOW-YZ]`）提供带 ETA 的实时进度反馈；0.25s 刷新节流，`\b` 原地更新；v3.2 改用 `persistent` 变量管理状态（解决值传递导致无法原地刷新的问题） | 2026-03-17 |
| 46  | 外部 JSON 配置系统（v3.1）：新增 `load_sim_config.m` / `save_sim_config.m`，支持从 JSON 文件加载参数 + 运行时自动保存 `config.json` 至输出目录；配置快照可直接复用以精确复现运行 | 2026-03-17 |
| 47  | 函数化封装 + 批量扫描（v3.1）：核心仿真逻辑封装为 `run_simulation(cfg)` 函数，支持编程式调用；`main_coupled_simulation.m` 简化为入口脚本；新增 `batch_sweep_example.m` 参数扫描示例 | 2026-03-17 |
| 48  | 参数扫描套件 `batch sweep/`（v3.3）：`batch_sweep_example.m` 重命名为 `batch_sweep_alpha.m` 并与 5 个新脚本统一归入 `batch sweep/` 子文件夹；各脚本分别扫描：初始倾斜角 `alpha_0`（7 工况）、胞体长半轴 `major_axis`（6 工况）、纤毛长度 `length_um`（6 工况）、纤毛附着点位 `site`（14 工况）、纤毛数量（6 工况）、纤毛杨氏模量 `E`（5 工况）；全部采用 headless 无头模式；代码审查修复杨氏模量脚本注释错误及 case_name 小数点问题（`1p77e07` 格式） | 2026-03-18 |
| 49  | JSON 批量运行器 `batch_run_from_json.m`（v3.3）：`uigetdir` 交互选择文件夹（支持硬编码路径旁路）→ 自动发现并字母序排列全部 `*.json` → `load_sim_config` 加载并强制 headless 覆盖 → `try-catch` 捕获异常跳过失败工况继续执行 → 运行结束输出对齐汇总表（状态 / 稳态 V_sink / 耗时 / 错误原因）→ 自动生成 `batch_summary_*.log`；ETA 基于历史均值动态估算；使用 `mean(...,'omitnan')` 消除 Statistics Toolbox 依赖 | 2026-03-18 |
| 50  | 批量脚本代码审查修复（v3.3）：修复 `batch_run_from_json.m` 的死代码（`cfg.generate_animation` 自赋值）与冗余赋值（`blog_fid=-1` 在 `fopen` 已返回 -1 后的重复写入） | 2026-03-18 |

---

## 23. 参考文献

1. **Walker, B.J. et al. (2020)**. "Filament mechanics in a half-space via regularised Stokeslet segments." _Physical Review E_, 102(2), 023101.
2. **Cortez, R. et al. (2005)**. "The method of regularized Stokeslets in three dimensions: analysis, validation, and application to helical swimming." _Physics of Fluids_, 17(3), 031504.
3. **Cortez, R. (2018)**. "Regularized Stokeslet segments." _Journal of Computational Physics_, 375, 783–796.
4. **Bergou, M. et al. (2008)**. "Discrete elastic rods." _ACM Transactions on Graphics_, 27(3), 63.

---

## 附录 A：关键标度推导

### A.1 时间尺度

$$
T_0 = \frac{8\pi\mu L_0^4}{EI} = \frac{8\pi\mu a^4}{E \cdot \pi r^4/4} = \frac{32\mu a^4}{E r^4}
$$

### A.2 弹性-重力数

$$
Eg = \frac{F_{\text{sink}} L_0^2}{EI} = \frac{\Delta\rho \cdot V_{\text{capsule}} \cdot g \cdot a^2}{EI}
$$

### A.3 扭转阻力 $\ell_k$ 因子

纤毛每段的有量纲扭转阻力为 $\tilde{c}_t = 4\pi\mu$。无量纲化后：

$$
\tilde{A}_k = \frac{\ell_k \varepsilon_k^2}{2} \hat{A}_k
$$

额外的 $\ell_k$ 因子确保线性系统中的量纲一致性。

### A.4 线性系统 $L\dot{z}=R$ 逐行量纲检查

所有行均已验证量纲一致：

- 第 1–3 行（力）：两侧均为 $[F_0]$
- 第 4–6 行（力矩）：两侧均为 $[F_0 \cdot L_0]$
- 纤毛力矩行：两侧均为 $[EI/L_0]$

---

## 许可证

MIT License

Copyright (c) 2026 Lucien

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
