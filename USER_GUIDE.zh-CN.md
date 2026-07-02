# 用户指南 — 细菌弹性流体动力学仿真

[English](USER_GUIDE.md) | 中文

> **作者:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **版本:** 7.9.4
> **更新日期:** 2026-07-02

本文档面向实际使用：如何运行仿真、如何选择配置方式、输出文件在哪里、以及如何做批量后处理。模型公式与数值方法见 [`docs/theory.zh-CN.md`](docs/theory.zh-CN.md)，版本历史见 [`CHANGELOG.md`](CHANGELOG.md)。

版本提示：v7.9.4 为分层静默完整性版：v7.9.3 的 `cfg.silent` 字段仅门控了 `dual_log`（约 100 处 `logf` 调用），仍遗漏了 `run_simulation.m` 中 5 处直接 `fprintf`（notably 每个 run 末尾的 `Log saved:` 与 `[OUTPUT]` 重命名提示）、16 个 `src/` 子函数中约 24 处 `warning` 调用，以及两个用 `\b` 退格符的就地进度条（`odetpbar`、`render_progress_bar`）。v7.9.4 通过分层设计封堵全部四个控制台输出通道：(1) 直接 `fprintf` 调用包裹 `if ~silent`；(2) `silent` 时不再挂载 `odetpbar` 作为 ODE `OutputFcn`；(3) `silent` 时不再调用 `render_progress_bar`；(4) 子函数 `warning` 调用通过 `warning('off','all')` 抑制，并用 `onCleanup` 守卫在任何退出路径恢复调用方的 warning 设置。**仿真输出无变化；`src/` 子函数未被修改；调用方的 warning 设置始终恢复。** v7.9.3 为批处理日志健壮性与并行 stdout 整洁性版；v7.9.2 为数值健壮性与文档准确性版；v7.9.1 为审查后 Bug 修复（MSD `D_trans` 下界保护、批处理错误栈追踪捕获）；v7.9.0 将所有起飞前参数检查统一收口到 `validate_simulation_params.m`（直接接受完整 `cfg` 结构，P1–P9）；v7.8.0 将内禀（参考）曲率 `kb0` 改用与运行时曲率相同的 SVD 投影节点 director 算子进行离散，使 `frozen` 曲率纤毛在参考构型下成为精确的零力矩平衡态，并将常值电机力矩改用 `EI/L_pili` 无量纲化（与主动力矩约定一致）。

## 1. 运行前准备

### 运行需求

- 推荐 MATLAB R2020b 或更新版本。
- PNG 和 PDF 图像导出依赖 `exportgraphics`；可选交互式 HTML 导出需要 MATLAB R2026a 或更新版本。
- 胶囊胞体表面布点需要 Anton Semechko 的 **S2 Sampling Toolbox**，可通过 MATLAB Add-On Explorer / File Exchange 安装。
- 运行仿真不需要 Python、Conda 或 pip 环境；[`requirements.txt`](requirements.txt) 已明确说明。
- **Parallel Computing Toolbox** 为可选：在使用 [`batch_sweep/batch_repeat_brownian.m`](batch_sweep/batch_repeat_brownian.m) 与 [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m) 默认并行分支时需要（`use_parallel = true`）。

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
| 复现或分享一个工况 | 从 [`configs/config_example.json`](configs/config_example.json) 复制 JSON 配置。 |
| 导入自定义胞体和纤毛几何 | 使用 [`configs/config_custom_example.json`](configs/config_custom_example.json)。 |
| 测试自由扭转马达鞭毛 | 使用 [`configs/config_motor_flagella_example.json`](configs/config_motor_flagella_example.json)。 |
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

1. 将 [`configs/config_example.json`](configs/config_example.json) 复制到本地目录，例如 `configs/my_case.json`。
2. 在 [`main_coupled_simulation.m`](main_coupled_simulation.m) 中设置 `config_file`：

```matlab
config_file = 'configs/my_case.json';
main_coupled_simulation
```

仓库在 `configs/` 下提供示例 JSON；本地工况请复制为 `configs/my_case.json` 等路径使用。`.gitignore` 中仅忽略 `output/`。

### 按几何模式划分的 JSON 字段

**参数化模式**（`geometry_mode: "parametric"`）：使用 `body`、`pili`、`material`、`fluid`、`simulation`、`brownian`、`output`、`flowfield`、`background_flow`。不得设置非空 `custom_geometry`（自 v7.4.4 起）。权威示例：[`configs/config_example.json`](configs/config_example.json)。

**自定义模式**（`geometry_mode: "custom"`）：提供 `custom_geometry`（导入胞体与纤毛点云）。不得同时设置参数化胞体字段（`major_axis`、`N_body` 等）或非空 `Pili_Matrix`。权威示例：[`configs/config_custom_example.json`](configs/config_custom_example.json)。

自 v7.4.4 起，当 `cfg` 上同时存在混用字段时，`validate_simulation_params.m` 会**双向**拒绝 custom/parametric 混用并给出明确错误。

**JSON 说明：** 当 `geometry_mode` 为 `"parametric"` 时，`load_sim_config.m` **不会**将顶层 `custom_geometry` 映射到 `cfg`。参数化 JSON 中残留的 `custom_geometry` 因此会在**加载阶段被静默忽略**，不会进入校验——请从 JSON 中删除该块，或改为 `"geometry_mode": "custom"`。

典型顶层布局包括：

- `case_name`、`geometry_mode`
- 参数化：`body` + `pili`；自定义：`custom_geometry`
- `material`、`fluid`、`simulation`、`brownian`、`output`、`flowfield`、`background_flow`

### 配置枚举（与代码一致）

| 字段 | 合法值 | 代码来源 |
| --- | --- | --- |
| `intrinsic_curvature` | 裸字符串（`frozen`、`zero`、`constant`、`helix`）**或**含 `.type` 的结构体（v7.9.0 P7） | `validate_simulation_params.m` |
| `intrinsic_curvature`（`constant`） | 含 3 元素 `kappa` 的结构体 | `validate_simulation_params.m` |
| `intrinsic_curvature`（`helix`） | 含有限标量 `kappa_d2`、`kappa_d3` 的结构体 | `validate_simulation_params.m` |
| `active_moment` | 裸字符串（`off`、`traveling_wave`、`custom`）**或**含 `.type` 的结构体（v7.9.0 P7） | `validate_simulation_params.m` |
| `motor_torque` | 裸字符串（`off`、`constant`）**或**含 `.type` 的结构体（v7.9.0 P7/P8） | `validate_simulation_params.m` |
| `background_flow.type` | `off`、`uniform`、`shear_z` | `coupled_background_flow.m` |
| `root_coupling` | `clamped`、`free_torsion` | `build_coupled_lin_sys.m` |

### 起飞前参数检查（v7.9.0）

[`validate_simulation_params.m`](src/config/validate_simulation_params.m) 在初始化前运行，并将所有错误聚合成单份报告。v7.9.0 新增或加强的检查：

| 检查项 | 条件 | 严重度 |
| --- | --- | --- |
| `pili_ext` 数量一致性（P1） | `numel(pili_ext)` != 纤毛数（parametric：`size(Pili_Matrix,2)`；custom：`numel(custom_geometry.pili)`） | 错误 |
| `delta` 范围（P2） | `delta <= 0` | 错误 |
| `delta` 范围（P2） | `delta >= pi/4`（奇异事件过早触发） | 警告 |
| `N_body` 上限（P4） | `N_body > 2000`（`O(N^2)` 内存） | 警告 |
| 背景流正交性（P5） | `shear_z` 且 `direction` 与 `shear_axis` 近平行（`|cos(theta)| > 0.99`，零涡量） | 警告 |
| `alpha_0` 范围（P6/A1） | 超出 `[-pi/2, pi/2]`；在 parametric 模式及 custom 模式下 `initial_orientation` 为空时生效 | 错误 |
| `AbsTol` / `RelTol` 措辞（P9） | 超出 `(0, 1]`；推荐 `1e-6 ~ 1e-3` | 警告 |

布朗 JSON 字段由 [`load_sim_config.m`](src/config/load_sim_config.m) 映射：`brownian.enable`、`brownian.temperature` -> `Temper`、`brownian.rng_seed`、`brownian.dt_B_ratio`。

权威格式以仓库中的示例 JSON 文件和 [`src/config/load_sim_config.m`](src/config/load_sim_config.m) 的实现为准。

## 5. 程序化运行与批量运行

脚本或参数扫描可直接调用仿真引擎：

```matlab
cfg = load_sim_config('configs/config_example.json');
cfg.auto_confirm = true;
[results, params, output_dir] = run_simulation(cfg);
```

批量扫描**驱动脚本**位于 [`batch_sweep/`](batch_sweep/)。运行后 CSV/MAT 等汇总输出写入仓库根目录下的 **`batch sweep/`**（带空格，由脚本自动创建；与源码目录 `batch_sweep/` 不同）。

### 5.1 布朗重复批量（`batch_repeat_brownian.m`）

从仓库根目录打开 [`batch_sweep/batch_repeat_brownian.m`](batch_sweep/batch_repeat_brownian.m)，对**同一套**布朗开启配置进行多次独立重复，输出写入 `output/<case_name_base>/001`、`002`、…（零填充子目录）。

- **并行（默认）：** `use_parallel = true` 在已许可 **Parallel Computing Toolbox** 且能创建 `parpool` 时使用 `parfor`。worker 数量由探测到的物理核心数与可用内存共同约束，默认最多使用约 `max_resource_frac`（默认 `0.75`）的“资源预算”，并配合 `mem_per_worker_gb` 估算内存上限。若无工具箱或建池失败，脚本会**告警并回退**到与旧版等价的串行 `for` 循环。
- **串行：** 设置 `use_parallel = false` 可强制单线程循环，便于与 v7.4 之前的行为逐字节对照。
- **随机种子：** 并行模式下，`brown_rng_seed_mode = 'fixed'` 对每个重复索引 `k` 设置 `cfg.brown_rng_seed = rng_seed_base + k`；`'timestamp'` 则使用时间基底 + `k`，避免多 worker 撞同一种子。**串行** `'timestamp'` 仍像单次交互运行一样不在每轮传入 `cfg.brown_rng_seed`。
- **worker 控制台静默（仅并行）：** parfor 分支设置 `cfg_k.silent = true`，因此每个 worker 的 `run_simulation` 只把完整日志写入 `output/<case>/<NNN>/simulation.log`，**不**回显到共享 stdout。这样可避免约 15000 行交错的 `[ODE] Segment N completed` 淹没有用的 `[k/100] OK` 进度行。串行分支不变（保持 verbose），便于交互式调试。

工况基目录下会写出带时间戳的 `batch_repeat_summary_*.txt`，记录每轮状态、耗时以及并行诊断字段（`n_pool`、`core_probe`、CPU/内存受限原因等）。per-run 记录**增量写入**（每次 `fopen(...,'a')` / `fclose` 强制刷盘），因此批次中途崩溃时摘要文件中仍保留头部与全部已完成 run 的记录，仅尾部（`elapsed_total_s`、`success_count`、加速比）丢失。并行模式下 per-run 行按 worker 完成顺序排列（行首 `NNN` 编号即 `k`，便于事后排序）；串行模式按 `k` 顺序写入。

### 5.2 JSON 文件夹批跑（`batch_run_from_json.m`）

从仓库根目录运行 [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m)，选择含 `*.json` 的文件夹（或在脚本顶部设置 `BATCH_JSON_DIR`）。每个入队文件经 `load_sim_config` 加载，强制无界面（`auto_confirm = true`，`preview_initial = false`），再调用 `run_simulation`。布朗 RNG、动画与流场开关均遵循各 JSON。

- **预扫描：** 按文件名字母序；**非空**重复 `case_name` 标记为 `SKIPPED`（排序第一个运行，其余不入队）。空 `case_name` 不去重（由 `run_simulation` 分配唯一时间戳输出目录）。
- **并行（默认）：** `use_parallel = true`，`mem_per_worker_gb = 3.0`，`max_resource_frac = 0.75`，`worker_blas_threads = 1`。已许可 **Parallel Computing Toolbox** 时使用 `parfor`；否则告警并回退串行。
- **实时日志：** 在 JSON 文件夹写入 `batch_summary_<时间戳>.log`；并行路径经 `parallel.pool.DataQueue` 推送，串行路径逐轮写入，均 `fflush` 刷盘。
- **汇总：** 控制台与日志报告 JSON 总数、入队数、成功/失败/跳过、各文件稳态 `V_sink`；并行时附带加速比估计。

```matlab
run('batch_sweep/batch_run_from_json.m');
```

## 6. 自定义几何模式

当胞体和纤毛都来自用户点云文件时，使用自定义几何模式。

1. 复制 [`configs/config_custom_example.json`](configs/config_custom_example.json)。
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

建议从 [`configs/config_motor_flagella_example.json`](configs/config_motor_flagella_example.json) 开始。

每根纤毛的关键字段包括：

- `root_coupling`: `"clamped"` 或 `"free_torsion"`。
- `motor_torque.type`: `"off"` 或 `"constant"`。
- `motor_torque.torque_pN_nm`: 马达力矩，单位 pN nm。
- `motor_torque.max_omega_Hz`: 可选停转速度限制，单位 Hz。

在 `free_torsion` 模式下，根部扭转角 $\psi_1$ 会动态演化。马达力矩只施加到纤毛根部 $d_3$ 力矩平衡行；不会额外加入 `R_global(4:6)`，因为这些行表示完整“胞体 + 纤毛”耦合系统的角动量平衡。

## 8. 布朗动力学

在 JSON 中设置 `brownian.enable`，或在 MATLAB 中设置 `cfg.enable_brown`。

自 v7.7.0 起，Brownian 更新采用严格 body-only Ito 形式：通过有限差分显式计算热漂移 `kBT*div(M_eff)`，并使用 SPD 投影后的扩散协方差。当前不再提供用户可配置的中点漂移选项。

重要选项包括：

- `brownian.temperature`: 正的流体绝对温度，单位 K；启用布朗运动时必需。
- `brownian.rng_seed`: 固定正整数种子，用于可复现运行。
- `brownian.dt_B_ratio`: 布朗时间步相对于输出时间步的比例。

启用布朗运动后，实际使用的随机种子会保存到 `config.json`、`brownian_rng_seed.txt`、`simulation.log` 和 `simulation_results.mat`。

**迁移率条件数告警（自 v7.4.1 起）：** 每次布朗节拍会在 [`src/physics/compute_effective_body_brownian.m`](src/physics/compute_effective_body_brownian.m) 内对已对称化的**无量纲**有效迁移率 `M_eff_dimless` 做 `rcond` 诊断；因此正常运行不应再只因 `M_eff_dim` 各分块的 SI 量纲差异而出现大量虚假的 `IllConditioned`。若 MATLAB 仍以极小的 **`M_eff_dimless` rcond** 抛出该告警，需排查欧拉根是否接近 \(\theta=\pi\)（钳制根部行退化）、体形极端纵横比等与几何/动力学相关的问题，并结合 `simulation.log` 与轨迹核对。

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
| `geometry_mode` 混用报错 | 仅使用 parametric **或** custom 字段；勿将 `custom_geometry` 与 `body`/`pili` 矩阵字段混用（v7.4.4+）。 |
| `helix` 内禀曲率报错 | 须提供结构体：`{"type":"helix","kappa_d2":…,"kappa_d3":…}`，且为有限标量；裸 `"helix"` 无效（v7.4.4+）。 |
| 输出目录已存在 | 交互式选择处理方式，或设置 `auto_confirm = true` 使用无人值守重命名策略。 |
| 布朗重复结果不一致 | 确认种子模式为固定种子还是时间戳种子。 |
| 布朗重复批量未并行执行 | 安装 **Parallel Computing Toolbox**、检查 `parpool` 集群配置，或在 `batch_repeat_brownian.m` 中设置 `use_parallel = false`。 |
| JSON 文件夹批跑未并行执行 | 同上，针对 `batch_run_from_json.m`，或设置 `use_parallel = false`。 |
| JSON 批跑因重复 `case_name` 被跳过 | 非空 `case_name` 仅字母序首个 JSON 会运行；请使用唯一 `case_name` 或留空以自动时间戳目录。 |
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
