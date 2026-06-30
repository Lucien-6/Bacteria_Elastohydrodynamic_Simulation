# 细菌弹性流体动力学仿真

[English](README.md) | 中文

> **作者:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **创建日期:** 2026-02-06  
> **更新日期:** 2026-06-30
> **版本:** 7.9.2
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

当前版本为 **7.9.2 (2026-06-30)**。本版是在 v7.9.1 基础上的数值健壮性与文档准确性更新：(1) 在 `reconstruct_geometry.m` 与 `quat_to_rotmat.m` 中加入防御性单位四元数归一化，防止 `ode15s` 内部试探步产生 `|q| != 1` 时通过非正交旋转矩阵悄悄使胞体变形；(2) 更正 Euler 角约定标注——实际为 Walker (2020) 的 `R = Rz(phi)·Ry(-theta)·Rz(psi)`（标准 ZYZ 主动旋转的转置），而非此前文档所写的 ZXZ，仅注释/文档变更，数值输出逐位一致；(3) 在 `theory.md §10` 中显式声明布朗路径的确定性/随机算子分裂近似（确定性求解用原始 `M_weighted`，布朗协方差用对称化 `M_eff`），并增加 `rel_asym > 1e-2` 量化诊断与用户指引；(4) 移除 `M_eff_dim` 的冗余对称化（代数精确——对称矩阵逐元素除以对称外积 `S * S'` 保持对称性）。**对于正常运行的仿真，四元数守卫为空操作、对称化移除为代数恒等，故输出与 v7.9.1 完全一致。** v7.9.1 为审查后 Bug 修复；v7.9.0 为参数验证深度重构；v7.8.0 修正了参考曲率 `kb0` 离散化与常值电机力矩无量纲化。完整历史详见 [`CHANGELOG.md`](CHANGELOG.md)。

## 许可证

MIT License。详见 [`LICENSE`](LICENSE)。
