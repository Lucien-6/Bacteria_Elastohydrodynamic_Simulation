# 细菌弹性流体动力学仿真

[English](README.md) | 中文

> **作者:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **创建日期:** 2026-02-06  
> **更新日期:** 2026-07-03
> **版本:** 7.9.6
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

当前版本为 **7.9.6 (2026-07-03)**。本版为控制台输出与异常健壮性修正版：对 `run_simulation.m` 日志显示链路的一次专项审查发现，`dual_log` 的控制台回显开关自 v7.9.3 起被写反——导致所有正常（非批处理）运行的 `logf` 控制台输出**全部消失**，而批处理/并行 worker 反而全部回显，与 `cfg.silent` 文档约定完全相反；`simulation.log` 文件本身未受影响。审查还发现自定义几何（`geometry_mode='custom'`）运行会在无法跳过的 `[PLOT] Generating summary figure` 步骤必然崩溃，原因是 `trisurf(ax, F, x, y, z, ...)` 不像 `surf`/`plot3` 那样支持把坐标轴句柄作为首个参数。以上两处均已修复。本版还新增了贯穿整个运行体的 `try/catch`，确保 `simulation.log` 始终记录完整堆栈的 `[FATAL]` 分节而不是静默截断；补齐了 `compute_plot_flowfield` 流场进度条上最后一处未受 `cfg.silent` 约束的缺口；把 `init_coupled_params.m`、`custom_body_building.m`、`custom_pili_building.m` 中9处此前仅走原生 `warning()` 的告警恢复为项目标准的 `logf_h` + `warning` 双通道模式；并对布朗运动逐段日志做了节流，避免长时程运行把 `simulation.log` 刷爆。**仿真物理/数值输出无任何变化。**

**7.9.5 (2026-07-03)** 将静默模式扩展到**自定义几何**：`custom_body_building` 与 `custom_pili_building` 接受注入的 `logf` 句柄（`[Custom Body]` / `[Custom Pili]` 诊断信息亦受 `cfg.silent` 约束），且 summary 图在 silent 下使用 `Visible='off'`（与动画图一致）。**`cfg.silent`** 为**仅运行时**批处理开关（不写入 JSON）：为 true 时抑制本 run 内 `logf` 控制台回显、直接 `fprintf`、ODE/动画进度条与 warning；完整文本仍写入 `simulation.log`；error 从不抑制。并行批处理驱动（如 `batch_repeat_brownian`）在 worker 上设 `cfg.silent = true`。

**7.9.4 (2026-07-02)** 是在 v7.9.3 基础上的分层静默完整性更新：v7.9.3 的 `cfg.silent` 字段仅门控了 `dual_log`（约 100 处 `logf` 调用），仍遗漏了 `run_simulation.m` 中 5 处直接 `fprintf`（ notably 每个 run 末尾的 `Log saved:` 与 `[OUTPUT]` 重命名提示）、16 个 `src/` 子函数中约 24 处 `warning` 调用，以及两个用 `\b` 退格符的就地进度条（`odetpbar`、`render_progress_bar`）——并行模式下这些仍会产生数千行交错的控制台输出与退格乱码。v7.9.4 通过分层设计封堵全部四个控制台输出通道：(1) 直接 `fprintf` 调用包裹 `if ~silent`；(2) `silent` 时不再挂载 `odetpbar` 作为 ODE `OutputFcn`；(3) `silent` 时不再调用 `render_progress_bar`；(4) 子函数 `warning` 调用通过 `warning('off','all')` 抑制，并用 `onCleanup` 守卫在任何退出路径（正常返回、error、Ctrl-C）恢复调用方的 warning 设置。**仿真输出无变化；调用方的 warning 设置始终恢复。** v7.9.3 为批处理日志健壮性与并行 stdout 整洁性版；v7.9.2 为数值健壮性与文档准确性版；v7.9.1 为审查后 Bug 修复；v7.9.0 为参数验证深度重构；v7.8.0 修正了参考曲率 `kb0` 离散化与常值电机力矩无量纲化。完整历史详见 [`CHANGELOG.md`](CHANGELOG.md)。

## 许可证

MIT License。详见 [`LICENSE`](LICENSE)。
