# 细菌弹性流体动力学仿真

[English](README.md) | 中文

> **作者:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **创建日期:** 2026-02-06  
> **更新日期:** 2026-11-04
> **版本:** 8.7.2
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

自定义胞体和纤毛几何可参考 [`configs/config_custom_example.json`](configs/config_custom_example.json)；自由扭转马达鞭毛可参考 [`configs/config_motor_flagella_example.json`](configs/config_motor_flagella_example.json)。接近真实的单根尾部鞭毛大肠杆菌游动案例（胞体 2.5 µm × 0.88 µm，左手螺旋鞭毛丝、螺距 2.22 µm，马达扭矩 1300 pN·nm）见 [`configs/config_ecoli_single_flagellum.json`](configs/config_ecoli_single_flagellum.json)。

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
| `flowfield_panels.png` / `.pdf` / 可选 `.html` | XZ / YZ 双平面流场图。 |
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
| [`src/core/LaggedJacobian.m`](src/core/LaggedJacobian.m) | Brownian `ode15s` 重启所用的滞后有限差分 Jacobian 句柄（v8.5.1）。 |
| [`src/core/build_coupled_lin_sys.m`](src/core/build_coupled_lin_sys.m) | 耦合线性系统组装。 |
| [`src/physics/compute_pili_curvature.m`](src/physics/compute_pili_curvature.m) | 共享的节点共位二阶曲率/扭转算子（用于力矩与 `kb0`）。 |
| [`src/physics/compute_pili_elastic_moments.m`](src/physics/compute_pili_elastic_moments.m) | Kirchhoff 弹性恢复力矩及根部边界规则。 |
| [`src/physics/compute_effective_body_brownian.m`](src/physics/compute_effective_body_brownian.m) | 完整耦合有效胞体迁移率、Brownian 协方差因子与无量纲迁移率。 |
| [`src/physics/sample_brownian_body_force.m`](src/physics/sample_brownian_body_force.m) | 每个 Brownian 区间的广义热力（布朗增量的力实现）。 |
| [`src/physics/background_flow_velocity.m`](src/physics/background_flow_velocity.m) | 共享背景流速度求值器（求解器与流场绘图共用）。 |
| [`src/visualization/compute_plot_flowfield.m`](src/visualization/compute_plot_flowfield.m) | 在实际达到的终止时间进行末状态流场计算与图像导出。 |
| [`post_processing/render_flowfield_from_results.m`](post_processing/render_flowfield_from_results.m) | 从 v8.7.0 的 `simulation_results.mat` 重算末状态、指定时刻或全时序流场动画。 |
| [`src/visualization/render_motion_animation.m`](src/visualization/render_motion_animation.m) | 共用的四面板运动动画渲染器（3D 与 XZ 总览，加上跟随质心的 3D 与 XZ 窗口，含图例、比例尺、腰部一侧的深绿色短标记、带播放速度的时间戳）：图形对象只创建一次，逐帧仅更新数据；输出分阶段耗时日志。跟随窗口的刻度间隔不细于 0.5 μm。 |
| [`src/visualization/add_gravity_indicator_outside.m`](src/visualization/add_gravity_indicator_outside.m) | 预览图、汇总图和动画图共用的重力箭头标注（可选参数 `head_scale`）。 |

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
