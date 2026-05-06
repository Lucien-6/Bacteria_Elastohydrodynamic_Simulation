# 用户指南 — 细菌胞体-纤毛耦合仿真

[English](USER_GUIDE.md) | 中文

> **作者:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **版本:** 7.4.1  
> **更新日期:** 2026-05-06

本文档面向实际使用：如何运行仿真、如何选择配置方式、输出文件在哪里、以及如何做批量后处理。模型公式与数值方法见 [`docs/theory.zh-CN.md`](docs/theory.zh-CN.md)，版本历史见 [`CHANGELOG.md`](CHANGELOG.md)。

## 1. 运行前准备

### 运行需求

- 推荐 MATLAB R2020b 或更新版本。
- PNG 和 PDF 图像导出依赖 `exportgraphics`；可选交互式 HTML 导出需要 MATLAB R2026a 或更新版本。
- 胶囊胞体表面布点需要 Anton Semechko 的 **S2 Sampling Toolbox**，可通过 MATLAB Add-On Explorer / File Exchange 安装。
- 运行仿真不需要 Python、Conda 或 pip 环境；[`requirements.txt`](requirements.txt) 已明确说明。
- **Parallel Computing Toolbox** 为可选：仅在使用 [`batch sweep/batch_repeat_brownian.m`](batch%20sweep/batch_repeat_brownian.m) 默认并行分支时需要（`use_parallel = true`）。

### 推荐工作目录

建议从仓库根目录打开 MATLAB，或在 MATLAB 中运行：

```matlab
cd('path/to/Coupled_Body_Pili')
```

这样 JSON 配置、几何文件和输出目录的相对路径会更稳定。

## 2. 选择运行方式

| 目标 | 推荐方式 |
| --- | --- |
| 快速测试或修改默认示例 | 编辑 [`main_coupled_simulation.m`](main_coupled_simulation.m)。 |
| 复现或分享一个工况 | 从 [`config_example.json`](config_example.json) 复制 JSON 配置。 |
| 导入自定义胞体和纤毛几何 | 使用 [`config_custom_example.json`](config_custom_example.json)。 |
| 测试自由扭转马达鞭毛 | 使用 [`config_motor_flagella_example.json`](config_motor_flagella_example.json)。 |
| 无人值守运行多个工况 | 调用 [`run_simulation.m`](run_simulation.m)，并设置 `cfg.auto_confirm = true`。 |

## 3. 使用内联参数运行

1. 打开 [`main_coupled_simulation.m`](main_coupled_simulation.m)。
2. 修改 `cfg` 配置块和 `pili_config` 条目。
3. 运行：

```matlab
main_coupled_simulation
```

当你正在交互式开发新工况，并希望配置直接保留在入口脚本附近时，适合使用这种方式。

## 4. 使用 JSON 配置运行

1. 将 [`config_example.json`](config_example.json) 复制到本地目录，例如 `configs/my_case.json`。
2. 在 [`main_coupled_simulation.m`](main_coupled_simulation.m) 中设置 `config_file`：

```matlab
config_file = 'configs/my_case.json';
main_coupled_simulation
```

仓库默认忽略 `configs/` 文件夹，因此本地路径和探索性参数不会被意外提交。

### JSON 顶层字段

完整 JSON 配置建议包含以下顶层字段：

- `case_name`
- `geometry_mode`
- `body`
- `pili`
- `material`
- `fluid`
- `simulation`
- `brownian`
- `output`
- `flowfield`
- `background_flow`

权威格式以仓库中的示例 JSON 文件和 [`load_sim_config.m`](load_sim_config.m) 的实现为准。

## 5. 程序化运行与批量运行

脚本或参数扫描可直接调用仿真引擎：

```matlab
cfg = load_sim_config('config_example.json');
cfg.auto_confirm = true;
[results, params, output_dir] = run_simulation(cfg);
```

批量扫描脚本位于 [`batch sweep/`](batch%20sweep/)。其中包含胞体分辨率、胞体正则化参数、纤毛长度、纤毛刚度、纤毛段密度、初始姿态角和布朗重复实验等扫描入口。

### 5.1 布朗重复批量（`batch_repeat_brownian.m`）

从仓库根目录打开 [`batch sweep/batch_repeat_brownian.m`](batch%20sweep/batch_repeat_brownian.m)，对**同一套**布朗开启配置进行多次独立重复，输出写入 `output/<case_name_base>/001`、`002`、…（零填充子目录）。

- **并行（默认）：** `use_parallel = true` 在已许可 **Parallel Computing Toolbox** 且能创建 `parpool` 时使用 `parfor`。worker 数量由探测到的物理核心数与可用内存共同约束，默认最多使用约 `max_resource_frac`（默认 `0.75`）的“资源预算”，并配合 `mem_per_worker_gb` 估算内存上限。若无工具箱或建池失败，脚本会**告警并回退**到与旧版等价的串行 `for` 循环。
- **串行：** 设置 `use_parallel = false` 可强制单线程循环，便于与 v7.4 之前的行为逐字节对照。
- **随机种子：** 并行模式下，`brown_rng_seed_mode = 'fixed'` 对每个重复索引 `k` 设置 `cfg.brown_rng_seed = rng_seed_base + k`；`'timestamp'` 则使用时间基底 + `k`，避免多 worker 撞同一种子。**串行** `'timestamp'` 仍像单次交互运行一样不在每轮传入 `cfg.brown_rng_seed`。

工况基目录下会写出带时间戳的 `batch_repeat_summary_*.txt`，记录每轮状态、耗时以及并行诊断字段（`n_pool`、`core_probe`、CPU/内存受限原因等）。

## 6. 自定义几何模式

当胞体和纤毛都来自用户点云文件时，使用自定义几何模式。

1. 复制 [`config_custom_example.json`](config_custom_example.json)。
2. 设置 `"geometry_mode": "custom"`。
3. 在 `"custom_geometry"` 中提供胞体和纤毛文件。
4. 参考 [`geometry/`](geometry/) 中的示例文件格式。

### 自定义几何规则

- 自定义胞体和自定义纤毛必须同时使用。
- 胞体点云支持 `.csv`、`.txt` 和 `.mat` 文件，且应包含三维空间坐标。
- 纤毛中心线必须预先离散为等长段；允许的相对段长离散度不超过 $10^{-3}$。
- 连续重复点或近重复点会在校验前自动移除。
- 对任意初始姿态，优先使用 `initial_orientation`，不要依赖只绕单轴的 `alpha_0`。
- 自定义模式下 `deflection_z` 表示绕对齐后切向的扭转，`deflection_x` 表示倾斜；参数化模式下二者的几何角色相反。

## 7. 马达与自由扭转运行

建议从 [`config_motor_flagella_example.json`](config_motor_flagella_example.json) 开始。

每根纤毛的关键字段包括：

- `root_coupling`: `"clamped"` 或 `"free_torsion"`。
- `motor_torque.type`: `"off"` 或 `"constant"`。
- `motor_torque.torque_pN_nm`: 马达力矩，单位 pN nm。
- `motor_torque.max_omega_Hz`: 可选停转速度限制，单位 Hz。

在 `free_torsion` 模式下，根部扭转角 $\psi_1$ 会动态演化。马达力矩只施加到纤毛根部 $d_3$ 力矩平衡行；不会额外加入 `R_global(4:6)`，因为这些行表示完整“胞体 + 纤毛”耦合系统的角动量平衡。

## 8. 布朗动力学

在 JSON 中设置 `brownian.enable`，或在 MATLAB 中设置 `cfg.enable_brown`。

重要选项包括：

- `brownian.rng_seed`: 固定正整数种子，用于可复现运行。
- `brownian.dt_B_ratio`: 布朗时间步相对于输出时间步的比例。
- `brownian.fixman_midpoint`: 启用 Fixman 中点修正，用于处理构型依赖迁移率。

启用布朗运动后，实际使用的随机种子会保存到 `config.json`、`brownian_rng_seed.txt`、`simulation.log` 和 `simulation_results.mat`。

**迁移率条件数告警（自 v7.4.1 起）：** 每次布朗节拍会在 [`compute_effective_body_brownian.m`](compute_effective_body_brownian.m) 内对已对称化的**无量纲**有效迁移率 `M_eff_dimless` 做 `rcond` 诊断；因此正常运行不应再只因 `M_eff_dim` 各分块的 SI 量纲差异而出现大量虚假的 `IllConditioned`。若 MATLAB 仍以极小的 **`M_eff_dimless` rcond** 抛出该告警，需排查欧拉根是否接近 \(\theta=\pi\)（钳制根部行退化）、体形极端纵横比等与几何/动力学相关的问题，并结合 `simulation.log` 与轨迹核对。

## 9. 主要输出

输出写入 `output/<case_name>/` 或自动生成的时间戳运行目录。

| 文件 | 用途 |
| --- | --- |
| `simulation_results.mat` | 完整 `results` 和 `params` 结构体。 |
| `timeseries.csv` | 便于外部分析的标量时间序列。 |
| `config.json` | 本次运行的配置快照，可用于复现。 |
| `summary_plot.png` / `.pdf` / 可选 `.html` | 主汇总图。 |
| `flowfield_panels.png` / `.pdf` / 可选 `.html` | XZ 和 YZ 双平面流场图。 |
| `initial_pose.png` / `.pdf` / 可选 `.html` | 开启预览时生成的初始姿态图。 |
| `bacteria_3d_motion.mp4` | 可选三维动画。 |
| `simulation.log` | 完整运行日志。 |
| `brownian_rng_seed.txt` | 启用布朗运动时写出的随机种子和复现提示。 |

## 10. 后处理已完成工况

汇总多个已完成工况：

1. 在 MATLAB 中运行 `addpath('post_processing')`。
2. 运行：

```matlab
postprocess_cases()
```

1. 选择包含各工况文件夹的父目录。
2. 选择要纳入统计的工况文件夹。
3. 选择输出 `.xlsx` 文件路径。

后处理脚本会读取 `timeseries.csv` 和 `config.json`，计算运动学指标，写出 Excel 工作簿，并导出 PNG、PDF 和可选 HTML 格式的汇总图。

## 11. 常见检查

| 现象 | 建议检查 |
| --- | --- |
| MATLAB 找不到函数 | 确认仓库根目录已加入 MATLAB path。 |
| 胶囊体生成失败 | 确认 S2 Sampling Toolbox 已安装且在 MATLAB path 中。 |
| JSON 工况运行后没有纤毛 | 检查 `pili` 是否为合法 JSON 数组，且每根纤毛字段完整。 |
| 自定义几何被拒绝 | 检查点云方向、重复点和纤毛等长段要求。 |
| 输出目录已存在 | 交互式选择处理方式，或设置 `auto_confirm = true` 使用无人值守重命名策略。 |
| 布朗重复结果不一致 | 确认种子模式为固定种子还是时间戳种子。 |
| 布朗重复批量未并行执行 | 安装 **Parallel Computing Toolbox**、检查 `parpool` 集群配置，或在 `batch_repeat_brownian.m` 中设置 `use_parallel = false`。 |
| `IllConditioned` 布朗迁移率告警（`M_eff_dimless rcond`） | v7.4.1+ 在无量纲迁移率上做诊断；若仍出现告警，多半是真实的病态几何/欧拉根问题，核对 §8、`simulation.log` 及 \(\theta\approx\pi\) 配置。 |

## 12. 文档导航

- [`README.md`](README.md): 英文项目入口。
- [`README.zh-CN.md`](README.zh-CN.md): 中文项目入口。
- [`USER_GUIDE.md`](USER_GUIDE.md): 英文用户指南。
- [`docs/theory.md`](docs/theory.md): 英文理论与数值模型。
- [`docs/theory.zh-CN.md`](docs/theory.zh-CN.md): 中文理论与数值模型。
- [`CHANGELOG.md`](CHANGELOG.md): 版本说明和兼容性变更。
- [`docs/index.html`](docs/index.html): 项目网页。

## 13. 许可证

MIT License。详见 [`LICENSE`](LICENSE)。
