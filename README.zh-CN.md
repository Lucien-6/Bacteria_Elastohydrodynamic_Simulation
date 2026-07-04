# 细菌弹性流体动力学仿真

[English](README.md) | 中文

> **作者:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **创建日期:** 2026-02-06  
> **更新日期:** 2026-07-03
> **版本:** 7.10.0
> **运行环境:** 纯 MATLAB 工程；运行仿真不需要 Python 或 pip 依赖。

本项目用于模拟刚性细菌胞体与一根或多根柔性纤毛 / 鞭毛在低雷诺数 Stokes 流中的重力沉降动力学。胞体通过表面正则化 Stokeslet 点云表示；每根纤毛作为不可伸长 Kirchhoff 弹性杆，并通过正则化 Stokeslet 线段法与流体耦合。

## 核心能力

- 支持参数化胶囊胞体或导入自定义胞体点云。
- 支持参数化纤毛或导入自定义纤毛中心线。
- RSM/RSS 混合流体力学模型，并包含非局部线段积分扭转迁移率、互易 Schur 互补耦合和物理标准 rotlet 缩放。
- Kirchhoff 弹性杆模型，支持内禀曲率和主动力矩。
- 可选胞体布朗运动，扩散张量来自完整胞体-纤毛耦合有效迁移率，并使用 SPD 投影协方差与显式严格 Ito 热漂移。
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

批量扫描**驱动脚本**位于 [`batch_sweep/`](batch_sweep/)；运行后 CSV/MAT 等汇总输出写入 `batch sweep/`（带空格，由脚本自动创建）。对文件夹内全部 `*.json` 进行批跑（可选并行、非空 `case_name` 去重、实时 `batch_summary_*.log`）请使用 [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m)。多工况后处理入口为 [`post_processing/postprocess_cases.m`](post_processing/postprocess_cases.m)。

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
| [`src/physics/compute_effective_body_brownian.m`](src/physics/compute_effective_body_brownian.m) | 全系统有效布朗迁移率、SPD 扩散张量和 Cholesky 噪声因子计算。 |
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

当前版本为 **7.10.0 (2026-07-03)**。本版新增 Linux 超算集群运行能力，并修复了同一轮审查中发现的两处跨平台批处理 bug。针对将 `batch_repeat_brownian.m` 部署到独占整节点、采用分布式/网络存储、MATLAB R2023a 的 Linux 超算分区所做的专项审查发现：多处 Windows 专属 API（`memory()`、`wmic`）在 Linux 上要么直接报错、要么静默把并行池规模算错；`VideoWriter` 的默认编码 profile 仅 Windows/Mac 支持；`parpool` 的 `JobStorageLocation` 默认路径在家目录为网络存储的集群上是已知会导致池启动竞态的隐患；而且完全看不到正在运行中的子任务进度。本版新增维护分支 **[`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m)** / **[`batch_sweep/batch_repeat_brownian_HPC_Linux.m`](batch_sweep/batch_repeat_brownian_HPC_Linux.m)**，针对 Linux 逐一解决了以上问题（平台自适应视频导出、`lscpu`/`/proc/meminfo` 资源探测、CPU/内存余量独立设置、抽样降低单 worker 内存占用、并行池任务存储目录迁移到节点本地、带节流的实时进度心跳），同时**完全不改动**原有面向 Windows 的脚本。此外，同一轮审查还发现 `capture_partial_outdir`（v7.9.5 引入的失败运行 `outdir`/`log_hint` 字段）在两个平台上其实都从未真正生效过——因为 MATLAB 函数抛出异常时不会把返回值传回调用者——现已在 `batch_repeat_brownian.m` 与 `batch_run_from_json.m` 中改为根据确定性的 `case_name` 反推预期路径来修复。**仿真物理/数值输出无任何变化。**

**7.9.6 (2026-07-03)** 为控制台输出与异常健壮性修正版：修复了自 v7.9.3 起导致所有正常（非批处理）运行的 `logf` 控制台输出全部消失的 `dual_log` 极性反转 bug，以及导致每次自定义几何运行都在 `[PLOT] Generating summary figure` 步骤必然崩溃的 `trisurf` 坐标轴指定问题，并新增了贯穿整个运行体的 `try/catch`，确保 `simulation.log` 始终记录完整的 `[FATAL]` 分节而不是静默截断。

**7.9.5 (2026-07-03)** 将静默模式扩展到**自定义几何**：`custom_body_building` 与 `custom_pili_building` 接受注入的 `logf` 句柄（`[Custom Body]` / `[Custom Pili]` 诊断信息亦受 `cfg.silent` 约束），且 summary 图在 silent 下使用 `Visible='off'`（与动画图一致）。**`cfg.silent`** 为**仅运行时**批处理开关（不写入 JSON）：为 true 时抑制本 run 内 `logf` 控制台回显、直接 `fprintf`、ODE/动画进度条与 warning；完整文本仍写入 `simulation.log`；error 从不抑制。并行批处理驱动（如 `batch_repeat_brownian`）在 worker 上设 `cfg.silent = true`。

**7.9.4 (2026-07-02)** 封堵了剩余的 `cfg.silent` 控制台输出缺口（直接 `fprintf`、`odetpbar`/`render_progress_bar` 进度条、16 个 `src/` 文件中约 24 处 `warning`）。7.9.3 引入了 `cfg.silent` 本身及增量式批次汇总写入；7.9.2 为数值健壮性版；7.9.1 为审查后 Bug 修复；7.9.0 为参数验证深度重构；7.8.0 修正了参考曲率 `kb0` 离散化与常值电机力矩无量纲化。完整历史详见 [`CHANGELOG.md`](CHANGELOG.md)。

## 许可证

MIT License。详见 [`LICENSE`](LICENSE)。
