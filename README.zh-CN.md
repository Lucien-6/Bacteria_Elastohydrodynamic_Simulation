# 细菌弹性流体动力学仿真

[English](README.md) | 中文

> **作者:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **创建日期:** 2026-02-06  
> **更新日期:** 2026-09-04
> **版本:** 8.5.0
> **运行环境:** 纯 MATLAB 工程；运行仿真不需要 Python 或 pip 依赖。

本项目用于模拟刚性细菌胞体与一根或多根柔性纤毛 / 鞭毛在低雷诺数 Stokes 流中的重力沉降动力学。胞体通过表面正则化 Stokeslet 点云表示；每根纤毛作为不可伸长 Kirchhoff 弹性杆，并通过正则化 Stokeslet 线段法与流体耦合。

## 核心能力

- 支持参数化胶囊胞体或导入自定义胞体点云。
- 支持参数化纤毛或导入自定义纤毛中心线。
- RSM/RSS 混合流体力学模型，并包含非局部线段积分扭转迁移率、互易 Schur 互补耦合和物理标准 rotlet 缩放。
- Kirchhoff 弹性杆模型，支持内禀曲率和主动力矩。
- 可选 body-only 胞体布朗运动：胞体 6 自由度协方差来自完整胞体-纤毛耦合有效迁移率，并采用无量纲 SPD 修复及独立可复现随机流；每个区间的增量以恒定广义热力经耦合 ODE 实现（隐式 Ito 漂移、无根部折点）。纤毛内部弹性模态本身不接受热噪声。
- 可选背景流：关闭、均匀流或有限厚度剪切流。
- 支持自由扭转根部边界条件和恒定马达力矩，用于马达驱动鞭毛。
- 支持 MATLAB 内联配置、JSON 配置和程序化批处理运行。
- 支持面向论文图像的导出：MAT、CSV、日志、PNG、矢量 PDF、可选交互式 HTML 和 MP4 动画。
- 可选**并行**批量驱动（[`batch_sweep/batch_repeat_brownian.m`](batch_sweep/batch_repeat_brownian.m)、[`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m)）：在已安装 **Parallel Computing Toolbox** 时使用 `parfor` 与自动池规模估计；无工具箱或池创建失败时自动回退为串行循环。
- **Linux 超算集群运行支持**（7.10.0 新增）：[`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m) / [`batch_sweep/batch_repeat_brownian_HPC_Linux.m`](batch_sweep/batch_repeat_brownian_HPC_Linux.m) 是专为无图形界面、独占整节点的 Linux 集群加固的维护分支：平台自适应视频导出、Linux 原生 CPU/内存探测、并行池任务存储目录迁移到节点本地、抽样降低单 worker 内存占用、实时子任务进度上报。原有 Windows/交互式脚本不受影响。

- 推荐 MATLAB R2020b 或更新版本。
- PNG 与 PDF 图像导出依赖 `exportgraphics`，MATLAB R2020a 起可用。
- 可选交互式 HTML 图像导出需要 MATLAB R2026a 或更新版本。
- 胶囊胞体表面布点需要 MATLAB Add-On / File Exchange 中 Anton Semechko 的 **S2 Sampling Toolbox**。安装后需确保 `ParticleSampleSphere`、`SpiralSampleSphere`、`IcosahedronMesh`、`QuadCubeMesh` 和 `SubdivideSphericalMesh` 在 MATLAB path 中。
- 运行仿真不需要 Python / Conda / pip。详见 [`requirements.txt`](requirements.txt)。
- **Parallel Computing Toolbox** 为可选依赖：在使用 [`batch_sweep/batch_repeat_brownian.m`](batch_sweep/batch_repeat_brownian.m) 与 [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m) 的默认并行分支时需要（`use_parallel = true`）。

## 快速开始

### 方式 A：内联参数

1. 打开 [`main_coupled_simulation.m`](main_coupled_simulation.m)。
2. 修改 `cfg` 配置块，包括几何、纤毛、材料、仿真时长、布朗运动和输出选项。
3. 运行：

```matlab
main_coupled_simulation
```

### 方式 B：JSON 配置

1. 将 [`configs/config_example.json`](configs/config_example.json) 复制到本地配置路径，例如 `configs/my_case.json`。
2. 在 [`main_coupled_simulation.m`](main_coupled_simulation.m) 中设置 `config_file`：

```matlab
config_file = 'configs/my_case.json';
main_coupled_simulation
```

自定义胞体和纤毛几何可参考 [`configs/config_custom_example.json`](configs/config_custom_example.json)；自由扭转马达鞭毛可参考 [`configs/config_motor_flagella_example.json`](configs/config_motor_flagella_example.json)。

### 方式 C：程序化调用

```matlab
cfg = load_sim_config('configs/config_example.json');
cfg.auto_confirm = true;
[results, params, output_dir] = run_simulation(cfg);
```

批量扫描**驱动脚本**位于 [`batch_sweep/`](batch_sweep/)；运行后 CSV/MAT 等汇总输出写入 `batch sweep/`（带空格，由脚本自动创建）。对文件夹内全部 `*.json` 进行批跑（可选并行、非空 `case_name` 去重、实时 `batch_summary_*.log`）请使用 [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m)。多工况后处理入口为 [`post_processing/postprocess_cases.m`](post_processing/postprocess_cases.m)：三维平移 MSD 拟合（\(\mathrm{MSD}=6D\tau+V^2\tau^2\)，拟合窗 \(0<\tau\le T/2\)）与 XZ 平面转动扩散 MSAD（\(\mathrm{MSAD}=2D_\mathrm{rot}\tau\)），导出 Excel 及 PNG/PDF/（可选 HTML）汇总图。

## 主要输出

结果写入 `output/<case_name>/` 或自动生成的时间戳运行目录。

| 文件 | 用途 |
| --- | --- |
| `simulation_results.mat` | 完整 `results` 与 `params` 结构体。 |
| `timeseries.csv` | 便于外部分析的标量时间序列。 |
| `config.json` | 本次运行配置快照，可用于复现。 |
| `summary_plot.png` / `.pdf` / 可选 `.html` | 轨迹与汇总图。 |
| `flowfield_panels.png` / `.pdf` / 可选 `.html` | XZ / YZ 双平面流场图。 |
| `initial_pose.png` / `.pdf` / 可选 `.html` | 初始姿态预览图，开启预览时生成。 |
| `bacteria_3d_motion.mp4` | 可选三维动画。 |
| `simulation.log` | 控制台与文件双通道运行日志。 |
| `brownian_rng_seed.txt` | 仅在启用布朗运动时写出随机种子。 |

`results.completed` 与 `results.termination_reason` 用于区分完整运行和超时截断运行。任一积分段超时后，程序保存已获得的 MAT/CSV/配置/日志数据并正常返回；批处理与汇总后处理不会把该部分轨迹纳入稳态统计。

## 核心文件

| 路径 | 作用 |
| --- | --- |
| [`main_coupled_simulation.m`](main_coupled_simulation.m) | 主脚本与内联配置入口。 |
| [`run_simulation.m`](run_simulation.m) | 仿真引擎、积分循环、导出和日志。 |
| [`src/config/load_sim_config.m`](src/config/load_sim_config.m) | JSON 配置读取。 |
| [`src/config/save_sim_config.m`](src/config/save_sim_config.m) | JSON 配置快照写出。 |
| [`src/config/validate_simulation_params.m`](src/config/validate_simulation_params.m) | 运行前参数校验。 |
| [`src/config/init_coupled_params.m`](src/config/init_coupled_params.m) | 几何、无量纲化、布朗参数和模型参数初始化。 |
| [`src/config/init_coupled_state.m`](src/config/init_coupled_state.m) | 初始状态向量构建。 |
| [`src/core/ode_coupled_rhs.m`](src/core/ode_coupled_rhs.m) | 耦合弹性流体动力学 ODE 右端项。 |
| [`src/core/build_coupled_lin_sys.m`](src/core/build_coupled_lin_sys.m) | 耦合线性系统组装。 |
| [`src/physics/compute_pili_curvature.m`](src/physics/compute_pili_curvature.m) | 共享的节点共位二阶曲率/扭转算子（用于力矩与 `kb0`）。 |
| [`src/physics/compute_pili_elastic_moments.m`](src/physics/compute_pili_elastic_moments.m) | Kirchhoff 弹性恢复力矩及根部边界规则。 |
| [`src/physics/compute_effective_body_brownian.m`](src/physics/compute_effective_body_brownian.m) | 完整耦合有效胞体迁移率、Brownian 协方差因子与无量纲迁移率。 |
| [`src/physics/sample_brownian_body_force.m`](src/physics/sample_brownian_body_force.m) | 每个 Brownian 区间的广义热力（布朗增量的力实现）。 |
| [`src/physics/background_flow_velocity.m`](src/physics/background_flow_velocity.m) | 共享背景流速度求值器（求解器与流场绘图共用）。 |
| [`src/visualization/compute_plot_flowfield.m`](src/visualization/compute_plot_flowfield.m) | 在实际达到的终止时间进行末状态流场计算与图像导出。 |

## 文档导航

- [`README.md`](README.md)：英文入口页。
- [`USER_GUIDE.md`](USER_GUIDE.md)：英文用户指南。
- [`USER_GUIDE.zh-CN.md`](USER_GUIDE.zh-CN.md)：中文用户指南，包含运行步骤、输出说明、自定义几何和批处理说明。
- [`docs/theory.md`](docs/theory.md)：英文理论文档，包含模型、无量纲化、线性系统、布朗动力学和验证说明。
- [`docs/theory.zh-CN.md`](docs/theory.zh-CN.md)：中文理论文档。
- [`CHANGELOG.md`](CHANGELOG.md)：版本变更与兼容性记录。
- [`docs/index.html`](docs/index.html)：项目网页。

## 版本说明

当前版本为 **8.5.0 (2026-09-04)**，改变 body-only 布朗运动的施加方式：每个 Brownian 区间的热增量以恒定广义力 `f_th = M_eff \ (dU/dt_B)` 经耦合 ODE 作用于胞体（`sample_brownian_body_force.m`、`ode_coupled_rhs` 第 8a 步），取代“瞬时位移 kick + 根部重投影”。纤毛现在连续跟随胞体的热运动，因此每次 kick 在第 1、2 段之间注入的伪曲率折点——以及它强加给 `ode15s` 的刚性高波数松弛——不再出现；Ito 漂移 `kBT div(M)` 由沿轨迹的迁移率梯度隐式生成，显式有限差分漂移（每个节点额外 12 次耦合求解）已移除。即使随机种子相同，Brownian 轨迹也与旧版本不同；确定性算例不受影响。2 型抛物线菌毛模板亦精确锚定于原点。**8.4.0 (2026-09-04)** 修正 `free_torsion` 根部扭转方程：根部 $\psi_1$ 行是整根鞭毛的轴向力矩平衡，现在只携带边界（马达）扭矩；此前内部 Kirchhoff 扭转力矩也被加入该行，导致马达扭矩仅与第 1 段的阻力平衡，稳态转子转速被放大为分段数 $N$ 倍。**8.3.0 (2026-09-04)** 为离散精度与几何一致性更新：Kirchhoff 曲率差分模板改为节点共位的二阶格式（共享算子 `compute_pili_curvature.m`，同时用于 `kb0`）；`pili_building.m` 的 14 个附着位点变换全部改为真旋转（此前位点 1、2、5、7、10、11、14 会镜像手性菌毛），`deflection_x` 改为右手约定；平面菌毛的放置与旧版完全一致。上一版本 **8.2.0 (2026-09-04)** 为全面数理审查后的物理正确性更新：鞭毛马达改用以转子**相对胞体**转速（完整轴向速率 `dpsi_1 + cos(theta_1) dphi_1 - Omega_body . d3`）表述的带符号（非对称）线性扭矩–转速关系，并明确符号约定（正扭矩 = 绕向外根部切向量右手旋转）；流场后处理改为由相对环境流的速度求解力密度，并将背景流叠加到绘制的流场；`validate_simulation_params` 对 2–5 型菌毛缺失/非法的 `shape_params` 给出明确报错；用户指南补充了全部曲率与力矩输入的 `1/L_pili`、`EI/L_pili` 单位说明。

**8.1.0 (2026-07-25)** 为批量后处理运动学功能更新：[`post_processing/postprocess_cases.m`](post_processing/postprocess_cases.m) 在 \(0<\tau\le T/2\) 上拟合三维平移扩散，并新增 XZ 平面转动扩散 \(D_\mathrm{rot}\)（优先体四元数，回退 `tilt_rad`）。流体/弹性/布朗仿真核心相对 8.0.1 **未改**。

**8.0.1 (2026-07-21)** 为在 8.0.0 之上的数值诊断加固补丁。稠密加权迁移率耗散检验改用全谱 `eig`（`eigs` 仅用于稀疏矩阵并带回退）；孤立 `rcond < 1e-14` 以 `NearSingularMatrix` 告警而不再中止；`SingularMatrix` 需同时满足大后向残差与近奇异 `rcond` 的联合门控；健康体检通过 `check_linear_solve` 探针评估求解质量。诊断快照字段为 `rcond_near_singular_threshold`。

**8.0.0 (2026-07-13)** 为面向数学物理正确性和运行契约的重大升级：不等段长平动—扭转加权互易、自定义菌毛等弦长 canonical 几何、完整 director 客观初始化、无量纲 body-only Brownian SPD 修复与余数 kick/独立 RNG、严格不可压缩 `shear_z` 正交性、超时部分结果安全返回，以及经 `simulation.log` / `results.diagnostics` 可审计的数值防线。

**自 8.0.0 起保留的破坏性变更：** `Pili_Matrix` 必须精确为 8 行；自定义或具名配置中的 `bending_stiffness` 会被拒绝，刚度统一由全局 `material.E` 与 `r_pili` 决定；非正交 `shear_z` 配置会被拒绝；超时运行在保存部分轨迹后返回 `results.completed=false` 与 `termination_reason='timeout'`。

**7.10.0 (2026-07-03)** 新增 Linux 超算集群运行能力，并修复两处跨平台批处理问题；该版本引入维护分支 **[`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m)** / **[`batch_sweep/batch_repeat_brownian_HPC_Linux.m`](batch_sweep/batch_repeat_brownian_HPC_Linux.m)**。

**7.9.6 (2026-07-03)** 为控制台输出与异常健壮性修正版：修复了自 v7.9.3 起导致所有正常（非批处理）运行的 `logf` 控制台输出全部消失的 `dual_log` 极性反转 bug，以及导致每次自定义几何运行都在 `[PLOT] Generating summary figure` 步骤必然崩溃的 `trisurf` 坐标轴指定问题，并新增了贯穿整个运行体的 `try/catch`，确保 `simulation.log` 始终记录完整的 `[FATAL]` 分节而不是静默截断。

**7.9.5 (2026-07-03)** 将静默模式扩展到**自定义几何**：`custom_body_building` 与 `custom_pili_building` 接受注入的 `logf` 句柄（`[Custom Body]` / `[Custom Pili]` 诊断信息亦受 `cfg.silent` 约束），且 summary 图在 silent 下使用 `Visible='off'`（与动画图一致）。**`cfg.silent`** 为**仅运行时**批处理开关（不写入 JSON）：为 true 时抑制本 run 内 `logf` 控制台回显、直接 `fprintf`、ODE/动画进度条与 warning；完整文本仍写入 `simulation.log`；error 从不抑制。并行批处理驱动（如 `batch_repeat_brownian`）在 worker 上设 `cfg.silent = true`。

**7.9.4 (2026-07-02)** 封堵了剩余的 `cfg.silent` 控制台输出缺口（直接 `fprintf`、`odetpbar`/`render_progress_bar` 进度条、16 个 `src/` 文件中约 24 处 `warning`）。7.9.3 引入了 `cfg.silent` 本身及增量式批次汇总写入；7.9.2 为数值健壮性版；7.9.1 为审查后 Bug 修复；7.9.0 为参数验证深度重构；7.8.0 修正了参考曲率 `kb0` 离散化与常值电机力矩无量纲化。完整历史详见 [`CHANGELOG.md`](CHANGELOG.md)。

## 许可证

MIT License。详见 [`LICENSE`](LICENSE)。
