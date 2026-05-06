# 细菌胞体-纤毛耦合弹性流体动力学仿真

[English](README.md) | 中文

> **作者:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **创建日期:** 2026-02-06  
> **更新日期:** 2026-05-06  
> **版本:** 7.4.0  
> **运行环境:** 纯 MATLAB 工程；运行仿真不需要 Python 或 pip 依赖。

本项目用于模拟刚性细菌胞体与一根或多根柔性纤毛 / 鞭毛在低雷诺数 Stokes 流中的重力沉降动力学。胞体通过表面正则化 Stokeslet 点云表示；每根纤毛作为不可伸长 Kirchhoff 弹性杆，并通过正则化 Stokeslet 线段法与流体耦合。

## 核心能力

- 支持参数化胶囊胞体或导入自定义胞体点云。
- 支持参数化纤毛或导入自定义纤毛中心线。
- RSM/RSS 混合流体力学模型，并包含 Rotlet 线段扭转耦合。
- Kirchhoff 弹性杆模型，支持内禀曲率和主动力矩。
- 可选胞体布朗运动，扩散张量来自完整胞体-纤毛耦合有效迁移率。
- 可选背景流：关闭、均匀流或有限厚度剪切流。
- 支持自由扭转根部边界条件和恒定马达力矩，用于马达驱动鞭毛。
- 支持 MATLAB 内联配置、JSON 配置和程序化批处理运行。
- 支持面向论文图像的导出：MAT、CSV、日志、PNG、矢量 PDF、可选交互式 HTML 和 MP4 动画。
- 可选**并行**布朗重复批量（[`batch sweep/batch_repeat_brownian.m`](batch%20sweep/batch_repeat_brownian.m)）：在已安装 **Parallel Computing Toolbox** 时使用 `parfor` 与自动池规模估计；无工具箱或池创建失败时自动回退为串行循环。

- 推荐 MATLAB R2020b 或更新版本。
- PNG 与 PDF 图像导出依赖 `exportgraphics`，MATLAB R2020a 起可用。
- 可选交互式 HTML 图像导出需要 MATLAB R2026a 或更新版本。
- 胶囊胞体表面布点需要 MATLAB Add-On / File Exchange 中 Anton Semechko 的 **S2 Sampling Toolbox**。安装后需确保 `ParticleSampleSphere`、`SpiralSampleSphere`、`IcosahedronMesh`、`QuadCubeMesh` 和 `SubdivideSphericalMesh` 在 MATLAB path 中。
- 运行仿真不需要 Python / Conda / pip。详见 [`requirements.txt`](requirements.txt)。
- **Parallel Computing Toolbox** 为可选依赖：仅在使用 [`batch sweep/batch_repeat_brownian.m`](batch%20sweep/batch_repeat_brownian.m) 的并行分支时需要（默认 `use_parallel = true`）。

## 快速开始

### 方式 A：内联参数

1. 打开 [`main_coupled_simulation.m`](main_coupled_simulation.m)。
2. 修改 `cfg` 配置块，包括几何、纤毛、材料、仿真时长、布朗运动和输出选项。
3. 运行：

```matlab
main_coupled_simulation
```

### 方式 B：JSON 配置

1. 将 [`config_example.json`](config_example.json) 复制到本地配置路径，例如 `configs/my_case.json`。
2. 在 [`main_coupled_simulation.m`](main_coupled_simulation.m) 中设置 `config_file`：

```matlab
config_file = 'configs/my_case.json';
main_coupled_simulation
```

自定义胞体和纤毛几何可参考 [`config_custom_example.json`](config_custom_example.json)；自由扭转马达鞭毛可参考 [`config_motor_flagella_example.json`](config_motor_flagella_example.json)。

### 方式 C：程序化调用

```matlab
cfg = load_sim_config('config_example.json');
cfg.auto_confirm = true;
[results, params, output_dir] = run_simulation(cfg);
```

批量扫描脚本位于 [`batch sweep/`](batch%20sweep/)，多工况后处理入口为 [`post_processing/postprocess_cases.m`](post_processing/postprocess_cases.m)。

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
| [`load_sim_config.m`](load_sim_config.m) | JSON 配置读取。 |
| [`save_sim_config.m`](save_sim_config.m) | JSON 配置快照写出。 |
| [`validate_simulation_params.m`](validate_simulation_params.m) | 运行前参数校验。 |
| [`init_coupled_params.m`](init_coupled_params.m) | 几何、无量纲化、布朗参数和模型参数初始化。 |
| [`init_coupled_state.m`](init_coupled_state.m) | 初始状态向量构建。 |
| [`ode_coupled_rhs.m`](ode_coupled_rhs.m) | 耦合弹性流体动力学 ODE 右端项。 |
| [`build_coupled_lin_sys.m`](build_coupled_lin_sys.m) | 耦合线性系统组装。 |
| [`compute_effective_body_brownian.m`](compute_effective_body_brownian.m) | 全系统有效布朗迁移率与扩散计算。 |
| [`compute_plot_flowfield.m`](compute_plot_flowfield.m) | 末状态流场计算与图像导出。 |

## 文档导航

- [`README.md`](README.md)：英文入口页。
- [`USER_GUIDE.md`](USER_GUIDE.md)：英文用户指南。
- [`USER_GUIDE.zh-CN.md`](USER_GUIDE.zh-CN.md)：中文用户指南，包含运行步骤、输出说明、自定义几何和批处理说明。
- [`docs/theory.md`](docs/theory.md)：英文理论文档，包含模型、无量纲化、线性系统、布朗动力学和验证说明。
- [`docs/theory.zh-CN.md`](docs/theory.zh-CN.md)：中文理论文档。
- [`CHANGELOG.md`](CHANGELOG.md)：版本变更与兼容性记录。
- [`docs/index.html`](docs/index.html)：项目网页。

## 版本说明

当前版本为 **7.4.0 (2026-05-06)**。在 **`batch_repeat_brownian.m`** 中新增可选 `parfor` 并行调度（独立布朗重复、资源上限约束的 worker 数、无 PCT 或建池失败时回退串行）。单次调用的 **`run_simulation`** 物理路径与 v7.3.x 一致。完整历史见 [`CHANGELOG.md`](CHANGELOG.md)。

## 许可证

MIT License。详见 [`LICENSE`](LICENSE)。
