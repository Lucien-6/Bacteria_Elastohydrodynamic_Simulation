# BacEHD — 细菌弹性流体动力学仿真

[English](README.md) | 中文

> **程序名称:** BacEHD（详见 [`cpc/PROGRAM_SUMMARY.md`](cpc/PROGRAM_SUMMARY.md) 与 [`CITATION.cff`](CITATION.cff)）  
> **作者:** 陆尚平 (Lucien, [lucien-6@qq.com](mailto:lucien-6@qq.com)), 赵坤 ([kzhao@uestc.edu.cn](mailto:kzhao@uestc.edu.cn))  
> **仓库:** [GitHub](https://github.com/Lucien-6/Bacteria_Elastohydrodynamic_Simulation)  
> **创建日期:** 2026-02-06  
> **更新日期:** 2026-11-25  
> **版本:** 10.2.2  
> **运行环境:** 纯基础 MATLAB R2022a+（核心仿真无需 Statistics Toolbox、MEX 或 Python；可选论文重绘脚本位于 `post_processing/`）。

本项目用于模拟刚性细胞胞体与一根或多根柔性纤毛 / 鞭毛在低雷诺数 Stokes 流中的动力学：重力沉降、马达力矩或主动力矩驱动的游动、背景流中的运动，以及布朗运动。胞体通过表面正则化 Stokeslet 点云表示；每根纤毛作为不可伸长 Kirchhoff 弹性杆，并通过正则化 Stokeslet 线段法与流体耦合。

## 核心能力

- 支持参数化胶囊胞体或导入自定义胞体点云。
- 支持参数化纤毛或导入自定义纤毛中心线。
- RSM/RSS 混合流体力学模型，并包含非局部线段积分扭转迁移率、互易 Schur 互补耦合和物理标准 rotlet 缩放，以及经标定的纤毛水动力半径（`material.hydro_radius_factor`）。
- Kirchhoff 弹性杆模型，支持内禀曲率和主动力矩。
- 可选 body-only 胞体布朗运动：胞体 6 自由度协方差来自完整胞体-纤毛耦合有效迁移率，并采用无量纲 SPD 修复及独立可复现随机流；每个区间的增量以恒定广义热力经耦合 ODE 实现（隐式 Ito 漂移、无根部折点）。未开启下一项时，纤毛内部弹性模态本身不接受热噪声。
- 可选纤毛热涨落（9.0.0+，10.0.0 全面加固，`brownian.pili_thermal`）：纤毛慢弯曲模态与胞体按联合迁移率共同抽取热力，并以受约束探测的中心线力密度及平衡轴向标量力矩密度（`brown_pili_tau`）精确实现，使热模态达到能量均分且扭转泄露可忽略。支持按时间尺度门槛筛选（`brownian.pili_tau_cut_s`）、近简并偏振对保护（`brownian.pili_protect_pairs`）及独立分纤毛随机流（`brownian.pili_rng_seeds`）。扭转和更高阶弯曲模态仍为确定性。适用于被动、本征为直线、根部固支的纤毛。见 [`docs/theory.zh-CN.md`](docs/theory.zh-CN.md) §10.4。
- 胞体框架 Block Schur 互补水动力预计算（`build_body_precomp`、`build_pili_mobility_blocks`）与快速冻结水动力三色雅可比（`ode_coupled_jacobian`，协同 `LaggedJacobian`），显著加速多纤毛与布朗仿真（10.0.0）。
- 自动化 CPC 审稿人验证测试套件（[`cpc/run_cpc_tests.m`](cpc/run_cpc_tests.m)、[`cpc/PROGRAM_SUMMARY.md`](cpc/PROGRAM_SUMMARY.md)）、外部细长体/悬臂梁校验驱动（[`validation/`](validation/)）及算子/并行/GPU 扩展基准（[`validation/cpc_performance_benchmark.m`](validation/cpc_performance_benchmark.m)、[`validation/cpc_parallel_and_flowfield_benchmark.m`](validation/cpc_parallel_and_flowfield_benchmark.m)）。
- 可选背景流：关闭、均匀流或有限厚度剪切流。
- 支持自由扭转根部边界条件和恒定马达力矩，用于马达驱动鞭毛。
- 双鞭毛莱茵衣藻三维游动：在实测双波数 1ω+2ω 平面驱动基础上叠加 17° 手性拍平面倾角与零寄生扭矩近端 1ω+2ω 面外弯矩，复现 1–2 Hz 自旋与左手超螺旋轨迹。详见 [`chlamydomonas/README.zh-CN.md`](chlamydomonas/README.zh-CN.md) 与 [`chlamydomonas/cases/run_3d_helix_2s.m`](chlamydomonas/cases/run_3d_helix_2s.m)。
- 支持 MATLAB 内联配置、JSON 配置和程序化批处理运行。
- 支持面向论文图像的导出：MAT、CSV、日志、PNG、矢量 PDF、可选交互式 HTML 和 MP4 动画。
- 高性能六面板三维/二维流线流场重建与动画（`render_flowfield_from_results`，10.2.0）：直接复用 Block Schur 力系数求解、三平面合批 FP64 GPU/BLAS3 速度核（`use_gpu = 'auto'`）、精确解析胶囊截线、线程池并行流线预计算（`use_parallel = 'auto'`）以及轻量缓存（~1.6–4.0 MiB/帧，视网格域大小而定）。
- 可选**并行**批量驱动（[`batch_sweep/batch_repeat_brownian.m`](batch_sweep/batch_repeat_brownian.m)、[`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m)）：在已安装 **Parallel Computing Toolbox** 时使用 `parfor` 与自动池规模估计；无工具箱或池创建失败时自动回退为串行循环。
- **Linux 超算集群运行支持**：[`run_simulation.m`](run_simulation.m) 现已原生统一支持无界面/Linux 环境（经 `render_motion_animation.m` 平台自适应 MP4/AVI 导出、`.avi` 检测、`cfg.silent` 与 `DataQueue` 实时进度心跳）；[`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m) 作为薄包装转发别名保留，向下兼容 [`batch_sweep/batch_repeat_brownian_HPC_Linux.m`](batch_sweep/batch_repeat_brownian_HPC_Linux.m)（Linux 原生 CPU/内存探测、节点本地 `parpool` 存储与重型输出抽样）。

## 目录架构

```text
<project_root>/
├── main_coupled_simulation.m    # 交互式内联 / JSON 主入口脚本
├── run_simulation.m             # 统一仿真主引擎（Windows 与 Linux HPC）
├── src/                         # 物种无关核心求解器模块
│   ├── config/                  #   JSON 读写、聚合式参数预校验与初始化
│   ├── geometry/                #   胶囊/自定义胞体与等弦长纤毛构建 + S2 网格缓存
│   ├── core/                    #   Block Schur 线性系统、RHS、三色雅可比与健康诊断
│   ├── physics/                 #   Kirchhoff 弹性、主动力矩、两层级耦合布朗动力学
│   ├── math/                    #   四元数、方向架、欧拉角与自适应 SO(3) 换基
│   ├── visualization/           #   四面板运动动画与 GPU/线程池六面板流场重建
│   └── utils/                   #   版本号、状态打包与 ETA/进度估计
├── configs/                     # 示例与典型物种 JSON 配置模板
├── geometry/                    # 示例、大肠杆菌与衣藻点云 CSV 及生成脚本
├── chlamydomonas/               # 可插拔双鞭毛衣藻逆向/三维螺旋游动模块
├── cpc/                         # CPC Program Summary 与 18+2 自动化审稿测试套件
├── validation/                  # 外部 SBT/悬臂梁校验与算子/并行/GPU 性能基准
├── batch_sweep/                 # 多核 parfor 与 Linux HPC 批量系综/扫参驱动
├── post_processing/             # 离线流场/动画重绘、MSD/MSAD 统计与可选 plot_*.py
├── docs/                        # 理论与数值模型文档（中英文）
└── MATLAB Add-On/               # 内置 MIT 许可 S2 Sampling Toolbox（A. Semechko）
```

## 环境要求

- MATLAB R2022a 或更新版本（不兼容 GNU Octave）；本项目在 Windows 11 R2026a 上开发测试，并在 Linux 超算 R2023a 上运行验证。更早版本不支持：每次运行用 `jsonencode(..., 'PrettyPrint', true)`（R2021a）写出 `config.json`，默认流场图使用 `clim`（R2022a）。流场色标采用纯基础 MATLAB 线性分位数插值，不依赖 Statistics and Machine Learning Toolbox。
- PNG 与 PDF 图像由 `exportgraphics` 导出；可选交互式 HTML 图像导出需要 MATLAB R2026a 或更新版本。
- 胶囊胞体表面布点使用 Anton Semechko 的 MIT 许可 **S2 Sampling Toolbox**，已随本仓库直接打包在 `MATLAB Add-On/Suite of functions to perform uniform sampling of a sphere/` 中（`cpc/run_cpc_tests.m` 与验证脚本自动加入路径，无需额外下载）。
- 运行核心 MATLAB 仿真不需要 Python / Conda / pip；仅可选的论文重绘脚本（`post_processing/plot_*.py`）使用 [`requirements.txt`](requirements.txt) 中列出的 Python 绘图包。
- **Parallel Computing Toolbox**（及可选 CUDA GPU）为可选依赖：用于批量并行扫描（`batch_sweep/`）以及多帧流场动画自动加速（`use_gpu = 'auto'`、`use_parallel = 'auto'`），无工具箱或无 GPU 时自动回退为串行 CPU。

## 快速开始

运行前请把项目根目录设为 MATLAB 当前文件夹，或把它加入 MATLAB 路径。

### 审稿复现测试与性能基准（CPC Test Run）

```matlab
addpath('cpc'); addpath('validation');
run_cpc_tests('SkipExternal', true)      % 18 项内置物理与数值核验（~22.5 s）
run_cpc_tests                            % 18 项内置 + 2 项外部文献基准（~1-2 min）
cpc_performance_benchmark                % Block Schur 与三色雅可比加速基准（~17 s）
cpc_parallel_and_flowfield_benchmark     % parfor 强扩展与 GPU 流场加速基准
```

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

自定义胞体和纤毛几何可参考 [`configs/config_custom_example.json`](configs/config_custom_example.json)；自由扭转马达鞭毛可参考 [`configs/config_motor_flagella_example.json`](configs/config_motor_flagella_example.json)。接近真实的单根尾部鞭毛大肠杆菌游动案例（胞体 2.5 µm × 0.88 µm，左手螺旋鞭毛丝、螺距 2.22 µm，马达扭矩 1300 pN·nm）见 [`configs/config_ecoli_single_flagellum.json`](configs/config_ecoli_single_flagellum.json)。双鞭毛莱茵衣藻三维游动可在 MATLAB 中直接打开 [`chlamydomonas/cases/run_3d_helix_2s.m`](chlamydomonas/cases/run_3d_helix_2s.m)（三维螺旋游动，2.0 s）或 [`chlamydomonas/cases/run_planar_sym50.m`](chlamydomonas/cases/run_planar_sym50.m)（平面全对称参考，50 拍）后点 Run，完整说明见 [`chlamydomonas/README.zh-CN.md`](chlamydomonas/README.zh-CN.md)。

本机专用的工况请复制到 `configs/` 下的本地路径（例如 `configs/my_case.json`）。仓库跟踪 `configs/` 下的示例 JSON，该目录不在 git 的忽略范围内，自建的工况文件会显示为未跟踪文件；运行结果写入已被忽略的 `output/`。

### 方式 C：程序化调用

```matlab
cfg = load_sim_config('configs/config_example.json');
cfg.auto_confirm = true;
[results, params, output_dir] = run_simulation(cfg);
```

批量扫描**驱动脚本**位于 [`batch_sweep/`](batch_sweep/)；运行后 CSV/MAT 等汇总输出写入 `batch sweep/`（带空格，由脚本自动创建）。对文件夹内全部 `*.json` 进行批跑（可选并行、非空 `case_name` 去重、实时 `batch_summary_*.log`）请使用 [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m)。多工况后处理入口为 [`post_processing/postprocess_cases.m`](post_processing/postprocess_cases.m)：三维平移 MSD 拟合（\(\mathrm{MSD}=6D\tau+V^2\tau^2\)，拟合窗 \(0<\tau\le T/2\)）与 XZ 平面转动扩散 MSAD（\(\mathrm{MSAD}=2D_\mathrm{rot}\tau\)），导出 Excel 及 PNG/PDF/（可选 HTML）汇总图。若只想根据已保存的 `simulation_results.mat` 重新渲染运动动画（更换抽帧间隔或帧率）而不重跑 ODE，请使用 [`post_processing/render_animation_from_results.m`](post_processing/render_animation_from_results.m)。

## 主要输出

结果写入 `output/<case_name>/` 或自动生成的时间戳运行目录。

| 文件 | 用途 |
| --- | --- |
| `simulation_results.mat` | 完整 `results` 与 `params` 结构体。 |
| `timeseries.csv` | 便于外部分析的标量时间序列。 |
| `config.json` | 本次运行配置快照，可用于复现。 |
| `summary_plot.png` / `.pdf` / 可选 `.html` | 轨迹与汇总图。 |
| `flowfield_panels.png` / `.pdf` / 可选 `.html` | XZ / YZ / XY 流场图（2×3）。 |
| `initial_pose.png` / `.pdf` / 可选 `.html` | 初始姿态预览图，开启预览时生成。 |
| `bacteria_3d_motion.mp4` | 可选三维动画（Linux 上为 `bacteria_3d_motion.avi`）。 |
| `bacteria_3d_motion_rerender.mp4` | 可选；由 `render_animation_from_results` 根据已保存结果重新渲染的动画。 |
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
| [`src/core/ode_coupled_jacobian.m`](src/core/ode_coupled_jacobian.m) | 面向 `ode15s` 的快速冻结水动力三色雅可比（10.0.0）。 |
| [`src/core/LaggedJacobian.m`](src/core/LaggedJacobian.m) | Brownian `ode15s` 重启所用的滞后雅可比句柄，默认由 `ode_coupled_jacobian` 刷新并自动回退有限差分。 |
| [`src/core/build_body_precomp.m`](src/core/build_body_precomp.m) | 胞体局部框架下刚体 RSM 块的 Cholesky 预计算，用于 Block Schur 加速（10.0.0）。 |
| [`src/core/build_pili_mobility_blocks.m`](src/core/build_pili_mobility_blocks.m) | Block Schur 路径所用的胞体框架非 BB 水动力迁移率子块（10.0.0）。 |
| [`src/core/build_coupled_lin_sys.m`](src/core/build_coupled_lin_sys.m) | 耦合线性系统组装（Block Schur 快速路径 + 稠密 LU 路径）。 |
| [`src/core/build_force_balance.m`](src/core/build_force_balance.m) | 耦合系统的力/力矩平衡算子，从组装中拆出（9.0.0）。 |
| [`src/physics/compute_pili_curvature.m`](src/physics/compute_pili_curvature.m) | 共享的节点共位二阶曲率/扭转算子（用于力矩与 `kb0`）。 |
| [`src/physics/compute_pili_elastic_moments.m`](src/physics/compute_pili_elastic_moments.m) | Kirchhoff 弹性恢复力矩及根部边界规则。 |
| [`src/physics/compute_effective_body_brownian.m`](src/physics/compute_effective_body_brownian.m) | 完整耦合有效胞体迁移率、Brownian 协方差因子与无量纲迁移率。 |
| [`src/physics/sample_brownian_body_force.m`](src/physics/sample_brownian_body_force.m) | 每个 Brownian 区间的广义热力（布朗增量的力实现）。 |
| [`src/physics/init_pili_brownian_modes.m`](src/physics/init_pili_brownian_modes.m) | 接受热噪声的纤毛慢弯曲模态（支持按阶数或 `tau_cut_s` 筛选及偏振对保护，9.0.0–10.0.0）。 |
| [`src/physics/compute_effective_coupled_brownian.m`](src/physics/compute_effective_coupled_brownian.m) | 胞体与纤毛热模态联合迁移率、中心线力与平衡轴向力矩模板（`P_load`/`P_tau`）及无矩阵漂移系数（9.0.0–10.0.0）。 |
| [`src/physics/sample_brownian_coupled_force.m`](src/physics/sample_brownian_coupled_force.m) | 每个 Brownian 区间内胞体与纤毛慢模态的热力及平衡轴向力矩抽样（9.0.0–10.0.0）。 |
| [`src/physics/pili_thermal_rhs.m`](src/physics/pili_thermal_rhs.m) | 纤毛热力密度与平衡轴向力矩密度在耦合系统右端项中的贡献（9.0.0–10.0.0）。 |
| [`src/physics/background_flow_velocity.m`](src/physics/background_flow_velocity.m) | 共享背景流速度求值器（求解器与流场绘图共用）。 |
| [`src/visualization/compute_plot_flowfield.m`](src/visualization/compute_plot_flowfield.m) | 在实际达到的终止时间进行末状态流场计算与图像导出。 |
| [`post_processing/render_flowfield_from_results.m`](post_processing/render_flowfield_from_results.m) | 从 v8.7.0 的 `simulation_results.mat` 重算末状态、指定时刻或全时序流场动画。 |
| [`src/visualization/render_motion_animation.m`](src/visualization/render_motion_animation.m) | 共用的四面板运动动画渲染器（3D 与 XZ 总览，加上跟随质心的 3D 与 XZ 窗口，含图例、比例尺、腰部一侧的深绿色短标记、带播放速度的时间戳）：图形对象只创建一次，逐帧仅更新数据；输出分阶段耗时日志。跟随窗口的刻度间隔不细于 0.5 μm。 |
| [`src/visualization/add_gravity_indicator_outside.m`](src/visualization/add_gravity_indicator_outside.m) | 预览图、汇总图和动画图共用的重力箭头标注（可选参数 `head_scale`）。 |
| [`chlamydomonas/chlamy_run_preset.m`](chlamydomonas/chlamy_run_preset.m) | 莱茵衣藻三维/平面标准预设入口（详见 [`chlamydomonas/README.zh-CN.md`](chlamydomonas/README.zh-CN.md)）。 |
| [`chlamydomonas/core/chlamy_run_case.m`](chlamydomonas/core/chlamy_run_case.m) | 通用衣藻算例运行器：实测步态表 + 三维面外弯矩与手性倾角覆写。 |
| [`chlamydomonas/core/chlamy_swim_analysis.m`](chlamydomonas/core/chlamy_swim_analysis.m) | 解缠 SO(3) 自旋、多圈螺旋参数、胞体间隙与三维轨道后处理分析。 |

## 文档导航

- [`README.md`](README.md)：英文入口页。
- [`USER_GUIDE.md`](USER_GUIDE.md)：英文用户指南。
- [`USER_GUIDE.zh-CN.md`](USER_GUIDE.zh-CN.md)：中文用户指南，包含运行步骤、输出说明、自定义几何和批处理说明。
- [`docs/theory.md`](docs/theory.md)：英文理论文档，包含模型、无量纲化、线性系统、布朗动力学和验证说明。
- [`docs/theory.zh-CN.md`](docs/theory.zh-CN.md)：中文理论文档。
- [`CHANGELOG.md`](CHANGELOG.md)：版本变更与兼容性记录。
- [`docs/index.html`](docs/index.html)：项目网页。

## 许可证

MIT License。详见 [`LICENSE`](LICENSE)。
