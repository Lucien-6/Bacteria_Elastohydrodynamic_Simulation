# 用户指南 — BacEHD（细菌弹性流体动力学仿真）

[English](USER_GUIDE.md) | 中文

> **作者:** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **版本:** 10.2.1
> **更新日期:** 2026-11-24

本文档面向实际使用：如何运行仿真、如何选择配置方式、输出文件在哪里、以及如何做批量后处理。模型公式与数值方法见 [`docs/theory.zh-CN.md`](docs/theory.zh-CN.md)，版本历史见 [`CHANGELOG.md`](CHANGELOG.md)。

## 1. 运行前准备

### 运行需求

- MATLAB R2022a 或更新版本。本项目在 Windows 11 R2026a 上开发和测试，并在 Linux 超算 R2023a 上运行验证。更早版本不支持：每次运行都用 `jsonencode(..., 'PrettyPrint', true)`（R2021a）写出 `config.json`；默认开启的流场图使用 `clim`（R2022a）。流场色标采用纯基础 MATLAB 分位数插值，不依赖 Statistics Toolbox。
- PNG 和 PDF 图像导出依赖 `exportgraphics`；可选交互式 HTML 导出需要 MATLAB R2026a 或更新版本。
- 胶囊胞体表面布点需要 Anton Semechko 的 **S2 Sampling Toolbox**，可通过 MATLAB Add-On Explorer / File Exchange 安装。
- 运行仿真不需要 Python、Conda 或 pip 环境；[`requirements.txt`](requirements.txt) 已明确说明。
- **Parallel Computing Toolbox** 为可选：在使用 [`batch_sweep/batch_repeat_brownian.m`](batch_sweep/batch_repeat_brownian.m) 与 [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m) 默认并行分支时需要（`use_parallel = true`）。
- **在 Linux 超算集群上运行（7.10.0 起）：** 请使用 [`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m) / [`batch_sweep/batch_repeat_brownian_HPC_Linux.m`](batch_sweep/batch_repeat_brownian_HPC_Linux.m) 而非原版脚本——详见 §5.3。

### 推荐工作目录

建议从仓库根目录打开 MATLAB，或在 MATLAB 中运行：

```matlab
cd('path/to/project_root')
```

这样 JSON 配置、几何文件和输出目录的相对路径会更稳定。

## 2. 选择运行方式

| 目标 | 推荐方式 |
| --- | --- |
| 快速测试或修改默认示例 | 编辑 [`main_coupled_simulation.m`](main_coupled_simulation.m)。 |
| 复现或分享一个工况 | 从 [`configs/config_example.json`](configs/config_example.json) 复制 JSON 配置。 |
| 模拟接近真实的大肠杆菌（单根尾部鞭毛、马达驱动） | 使用 [`configs/config_ecoli_single_flagellum.json`](configs/config_ecoli_single_flagellum.json)；几何文件由 [`geometry/generate_ecoli_geometry.m`](geometry/generate_ecoli_geometry.m) 生成。 |
| 导入自定义胞体和纤毛几何 | 使用 [`configs/config_custom_example.json`](configs/config_custom_example.json)。 |
| 测试自由扭转马达鞭毛 | 使用 [`configs/config_motor_flagella_example.json`](configs/config_motor_flagella_example.json)。 |
| 模拟三维双鞭毛莱茵衣藻游动 | 打开 [`chlamydomonas/cases/run_3d_helix_2s.m`](chlamydomonas/cases/run_3d_helix_2s.m) 或 [`chlamydomonas/cases/run_planar_sym50.m`](chlamydomonas/cases/run_planar_sym50.m) 后点 Run（详见 [`chlamydomonas/README.zh-CN.md`](chlamydomonas/README.zh-CN.md)）。 |
| 无人值守运行多个工况 | 调用 [`run_simulation.m`](run_simulation.m)，并设置 `cfg.auto_confirm = true`。 |

### 双鞭毛莱茵衣藻

几何配置为 [`configs/config_chlamy_biflagellate_gz.json`](configs/config_chlamy_biflagellate_gz.json)：5 µm × 4 µm 旋转椭球体与两根 13° 固支鞭毛。JSON 内的行波项仅作占位，运行时由 `chlamydomonas/data/fit_1p00A.mat`（57.75 Hz 双波数 1ω+2ω 实测拟合表）与可选的三维面外弯矩（`chlamydomonas/core/chlamy_make_perp_envelope.m`）替换。完整模型说明与参数见 [`chlamydomonas/README.zh-CN.md`](chlamydomonas/README.zh-CN.md)。

```matlab
run('chlamydomonas/cases/run_3d_helix_2s.m')    % 三维螺旋游动，2.0 s，含动画
run('chlamydomonas/cases/run_planar_sym50.m')   % 纯平面全对称参考，50 拍，含动画
```

`tilt_deg`、`drive_amp`、`drive_phase`、`drive_freq_ratio`、`drive_amp_perp` 和 `perp_sign` 是 `chlamy_run_case` 对每根鞭毛的覆写项。预设脚本会在算例目录内自动保存 `swim_analysis.mat`。

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

仓库在 `configs/` 下提供示例 JSON；本地工况请复制为 `configs/my_case.json` 等路径使用。`configs/` 下的文件不在 git 的忽略范围内，自建的工况文件会显示为未跟踪文件；运行结果写入已被忽略的 `output/`。

### 按几何模式划分的 JSON 字段

**参数化模式**（`geometry_mode: "parametric"`）：使用 `body`、`pili`、`material`、`fluid`、`simulation`、`brownian`、`output`、`flowfield`、`background_flow`。不得设置非空 `custom_geometry`（自 v7.4.4 起）。权威示例：[`configs/config_example.json`](configs/config_example.json)。

**自定义模式**（`geometry_mode: "custom"`）：提供 `custom_geometry`（导入胞体与纤毛点云）。不得同时设置参数化胞体字段（`major_axis`、`N_body` 等）或非空 `Pili_Matrix`。权威示例：[`configs/config_custom_example.json`](configs/config_custom_example.json)。

自 v7.4.4 起，当 `cfg` 上同时存在混用字段时，`validate_simulation_params.m` 会**双向**拒绝 custom/parametric 混用并给出明确错误。

**JSON 说明：** 当 `geometry_mode` 为 `"parametric"` 时，`load_sim_config.m` **不会**将顶层 `custom_geometry` 映射到 `cfg`。参数化 JSON 中残留的 `custom_geometry` 因此会在**加载阶段被静默忽略**，不会进入校验——请从 JSON 中删除该块，或改为 `"geometry_mode": "custom"`。

典型顶层布局包括：

- `case_name`、`geometry_mode`
- 参数化：`body` + `pili`；自定义：`custom_geometry`
- `material`、`fluid`、`simulation`、`brownian`、`output`、`flowfield`、`background_flow`

`material` 段还接受可选键 `hydro_radius_factor`（默认 `sqrt(2)`，对应 `cfg.hydro_radius_factor`）：菌毛水动力正则化长度为 `eps_pili = hydro_radius_factor * r_pili`，而弯曲刚度始终使用物理半径 `r_pili`。`sqrt(2)` 使扭转阻力与无滑移圆柱精确一致、平动迁移率偏差约 3% 以内；`sqrt(e)` 使平动精确；`1` 复现 8.5.1 及以前的版本。超出 [1, 3] 的值会触发校验警告。

### 配置枚举（与代码一致）

| 字段 | 合法值 | 代码来源 |
| --- | --- | --- |
| `intrinsic_curvature` | 裸字符串（`frozen`、`zero`、`constant`、`helix`）**或**含 `.type` 的结构体 | `validate_simulation_params.m` |
| `intrinsic_curvature`（`constant`） | 含 3 元素 `kappa` 的结构体 | `validate_simulation_params.m` |
| `intrinsic_curvature`（`helix`） | 含有限标量 `kappa_d2`、`kappa_d3` 的结构体 | `validate_simulation_params.m` |
| `active_moment` | 裸字符串（`off`、`traveling_wave`、`custom`）**或**含 `.type` 的结构体 | `validate_simulation_params.m` |
| `motor_torque` | 裸字符串（`off`、`constant`）**或**含 `.type` 的结构体 | `validate_simulation_params.m` |
| `background_flow.type` | `off`、`uniform`、`shear_z` | `coupled_background_flow.m` |
| `root_coupling` | `clamped`、`free_torsion` | `build_coupled_lin_sys.m` |

### 菌毛曲率与力矩输入的单位

每根菌毛的弹性杆方程以**菌毛长度** $L_{\mathrm{pili}}$（而非全局体尺度 $L_0$）无量纲化：弧长取 $s \in [0, 1]$，节点间距为 $1/N_{\mathrm{seg}}$。因此所有曲率量都是以 $1/L_{\mathrm{pili}}$ 为单位的无量纲数，所有分布力矩都以 $EI/L_{\mathrm{pili}}$ 为单位。

| 输入 | 含义 | 单位 / 换算 |
| --- | --- | --- |
| `intrinsic_curvature.kappa = [k1, k2, k3]`（`constant`） | 参考 $\kappa_1 = \mathbf{d}_3\cdot\mathbf{d}_2'$（绕 $\mathbf{d}_1$ 弯曲）、$\kappa_2 = \mathbf{d}_1\cdot\mathbf{d}_3'$（绕 $\mathbf{d}_2$ 弯曲）、$\kappa_3 = \mathbf{d}_2\cdot\mathbf{d}_1'$（扭转） | 无量纲；物理值 $= k_i / L_{\mathrm{pili}}$ [1/m] |
| `intrinsic_curvature.kappa_d2`、`kappa_d3`（`helix`） | 均匀的绕 $\mathbf{d}_2$ 参考弯曲与参考扭转 | 无量纲；对半径 $R$、约化螺距 $p$（螺距 $=2\pi p$）的螺旋：`kappa_d2` $= L_{\mathrm{pili}}\,R/(R^2+p^2)$，`kappa_d3` $= L_{\mathrm{pili}}\,p/(R^2+p^2)$ |
| `active_moment.amplitude_Nm` | 内部主动力矩幅值 | 物理量 [N m]；`amplitude_dimless = amplitude_Nm * L_pili / EI` |
| `motor_torque.torque_pN_nm` | 马达失速扭矩 | 物理量 [pN nm]；`torque_dimless = torque / (EI/L_pili)` |
| `motor_torque.max_omega_Hz` | 零扭矩（空载）转速 | 物理量 [Hz]；`max_omega_dimless = 2*pi*max_omega_Hz*T_0` |

由于曲率输入按菌毛长度缩放，同一个 `kappa` 值用于不同 `length_um` 的菌毛时对应不同的物理曲率。`intrinsic_curvature = "frozen"`（默认）直接取初始几何为参考曲率，无需换算。

### 运行前参数检查

[`validate_simulation_params.m`](src/config/validate_simulation_params.m) 在初始化前运行，并将所有错误聚合成单份报告。最常拦下运行的检查：

| 检查项 | 条件 | 严重度 |
| --- | --- | --- |
| `Pili_Matrix` 形状 | 行数不是 8（每根纤毛一列）；不接受第 9 行的逐纤毛刚度，所有纤毛共用同一组 `E` 和 `r_pili` | 错误 |
| `pili_ext` 数量一致性 | `numel(pili_ext)` != 纤毛数（parametric：`size(Pili_Matrix,2)`；custom：`numel(custom_geometry.pili)`） | 错误 |
| `delta` 范围 | `delta <= 0` | 错误 |
| `delta` 范围 | `delta >= pi/4`（奇异事件过早触发） | 警告 |
| `N_body` 上限 | `N_body > 2000`（`O(N^2)` 内存） | 警告 |
| 背景流正交性 | `shear_z` 要求 `direction` 与 `shear_axis` 均为有限非零向量，且 `abs(dot(direction_hat,shear_axis_hat)) <= 1e-12` | 错误 |
| `alpha_0` 范围 | 超出 `[-pi/2, pi/2]`；在 parametric 模式及 custom 模式下 `initial_orientation` 为空时生效 | 错误 |
| `AbsTol` / `RelTol` 范围 | 超出 `(0, 1]`；推荐 `1e-6 ~ 1e-3` | 警告 |

布朗 JSON 字段由 [`load_sim_config.m`](src/config/load_sim_config.m) 映射：`brownian.enable`、`brownian.temperature` -> `Temper`、`brownian.rng_seed`、`brownian.dt_B_ratio`、`brownian.pili_thermal` -> `enable_pili_brown`、`brownian.pili_mode_orders` -> `pili_brown_orders`。

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

- **并行（默认）：** `use_parallel = true` 在已许可 **Parallel Computing Toolbox** 且能创建 `parpool` 时使用 `parfor`。worker 数量由探测到的物理核心数与可用内存共同约束，默认最多使用约 `max_resource_frac`（默认 `0.90`）的“资源预算”，并配合 `mem_per_worker_gb`（默认 `2.5`）估算内存上限。若无工具箱或建池失败，脚本会**告警并回退**到与旧版等价的串行 `for` 循环。
- **串行：** 设置 `use_parallel = false` 可强制单线程循环，便于与 v7.4 之前的行为逐字节对照。
- **随机种子：** 并行模式下，`brown_rng_seed_mode = 'fixed'` 对每个重复索引 `k` 设置 `cfg.brown_rng_seed = rng_seed_base + k`；`'timestamp'` 则使用时间基底 + `k`，避免多 worker 撞同一种子。**串行** `'timestamp'` 仍像单次交互运行一样不在每轮传入 `cfg.brown_rng_seed`。
- **worker 控制台静默（仅并行）：** parfor 分支设置 `cfg_k.silent = true`，因此每个 worker 的 `run_simulation` 只把完整日志写入 `output/<case>/<NNN>/simulation.log`，**不**回显到共享 stdout。这样可避免约 15000 行交错的 `[ODE] Segment N completed` 淹没有用的 `[k/100] OK` 进度行。自 v7.9.5 起，自定义几何 run 在 worker 上亦完全静默。串行分支不变（保持 verbose），便于交互式调试。

工况基目录下会写出带时间戳的 `batch_repeat_summary_<batch_id>.txt`，记录批次元数据、每轮状态、耗时与并行诊断。**头部（v7.9.5+）：** `batch_id`、人读格式 `started`、`config_source_type` / `config_json` 或 `config_source_note`、指向 `batch_config_snapshot_<batch_id>.json` 的路径（完整 inline `cfg` 副本，便于复现）、`project_root`、`script_path`、`brown_rng_policy`、`parallel_requested`、`execution_mode`、`summary_record_order`、可选 `[POOL] Resource probe` 字段（串行回退时亦写入）、以及串行执行但已预计算并行种子时的 `seeds_precomputed_unused`。**per-run 记录** **增量写入**（每轮 `fopen(...,'a')` / `fclose`）：OK 行含 `out:`（相对路径）；FAILED 行含 `Identifier` / `Message` / `Location`，并在可用时给出 `outdir` / `log_hint` 指向部分输出。**尾部：** `finished`、`elapsed_total_s`、`success_count`、并行加速比字段（如适用）、`failed_count` / `failed_indices` / `failed_runs`，以及异常退出后仍为 `PENDING` 的 `pending_count` / `pending_indices`。批次中途崩溃时摘要仍保留头部与全部已完成 run；仅尾部可能缺失。并行模式下 per-run 行按 worker 完成顺序排列（行首 `NNN` 即 `k`）；串行模式按 `k` 顺序写入。

### 5.2 JSON 文件夹批跑（`batch_run_from_json.m`）

从仓库根目录运行 [`batch_sweep/batch_run_from_json.m`](batch_sweep/batch_run_from_json.m)，选择含 `*.json` 的文件夹（或在脚本顶部设置 `BATCH_JSON_DIR`）。每个入队文件经 `load_sim_config` 加载，强制无界面（`auto_confirm = true`，`preview_initial = false`），再调用 `run_simulation`。布朗 RNG、动画与流场开关均遵循各 JSON。

- **预扫描：** 按文件名字母序；**非空**重复 `case_name` 标记为 `SKIPPED`（排序第一个运行，其余不入队）。空 `case_name` 不去重（由 `run_simulation` 分配唯一时间戳输出目录）。
- **并行（默认）：** `use_parallel = true`，`mem_per_worker_gb = 3.0`，`max_resource_frac = 0.75`，`worker_blas_threads = 1`。已许可 **Parallel Computing Toolbox** 时使用 `parfor`；否则告警并回退串行。
- **实时日志：** 在 JSON 文件夹写入 `batch_summary_<时间戳>.log`；并行路径经 `parallel.pool.DataQueue` 推送，串行路径逐轮写入，均 `fflush` 刷盘。
- **汇总：** 控制台与日志报告 JSON 总数、入队数、成功/失败/跳过、各文件稳态 `V_sink`；并行时附带加速比估计。

```matlab
run('batch_sweep/batch_run_from_json.m');
```

### 5.3 Linux 超算集群批跑（`batch_repeat_brownian_HPC_Linux.m`，7.10.0 起）

如需在无图形界面、独占整节点的 Linux 超算集群上运行，请使用 [`batch_sweep/batch_repeat_brownian_HPC_Linux.m`](batch_sweep/batch_repeat_brownian_HPC_Linux.m)。自 v10.0.0 起，所有 Linux/无界面求解器特性（平台自适应 AVI/MP4 导出、`.avi` 检测、`cfg.silent` 与 `DataQueue` 进度心跳）已统一合入 `run_simulation.m`，[`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m) 保留为转发至 `run_simulation.m` 的薄包装别名。用法与 §5.1 完全一致（同样的 `case_name_base` / `n_repeats` / `brown_rng_seed_mode` / 内联或 JSON 配置方式）。

与 `batch_repeat_brownian.m` 的差异：

- **平台自适应视频导出：** `VideoWriter` 在 Windows/Mac 上仍用 `'MPEG-4'`（`.mp4`），但在 Linux 上改用 `'Motion JPEG AVI'`（`.avi`）——MathWorks 官方文档明确 `'MPEG-4'` 仅支持 Windows/Mac，在任何 Linux 系统上都会报错，与 MATLAB 版本、是否有图形界面无关。
- **Linux 原生资源探测：** 当 Windows 专属的 `wmic`/`memory()` 探测不可用时，物理核心数会退回 `lscpu -p=Core,Socket`，可用内存会退回解析 `/proc/meminfo` 的 `MemAvailable` 字段，取代了此前"假设只有 4GB 可用"这个会在真实 Linux 节点上把并行池规模严重算小的兜底逻辑。
- **CPU/内存余量独立设置：** `cpu_resource_frac`（默认 `0.90`）与 `mem_resource_frac`（默认 `0.65`，更保守）取代了原来单一的 `max_resource_frac`。
- **重度输出抽样：** 无论内联配置还是 JSON 怎么设置，`generate_animation` / `compute_flowfield` 默认只在 `heavy_output_sample_idx`（默认第一个、中间一个、最后一个）指定的重复里打开，其余全部强制关闭，以控制单个 worker 的峰值内存。确认内存余量足够安全后，可将 `heavy_output_sample_idx` 设为 `1:n_repeats` 以恢复每次重复都输出完整内容。
- **并行池任务存储目录迁移到节点本地：** 在 Linux 下，本地并行池的 `JobStorageLocation` 会被迁移到 `tempdir()`（自动识别作业调度系统设置的 `$TMPDIR`，通常是节点本地临时存储）下一个按 `SLURM_JOB_ID` / `PBS_JOBID` 命名的专属子目录，避免 MATLAB 默认路径（用户家目录下）在家目录为分布式/网络存储的集群上已知会出现的并行池启动竞态问题。若迁移本身失败会自动退回原有方式，Windows 完全不受影响。
- **实时子任务进度心跳：** 尽管 `cfg.silent = true` 会抑制某次重复运行的全部正常控制台回显，仍会按每次重复至少间隔 30 秒的节流频率打印一行 `[HB] k=NNN running  PP.P%  elapsed=...`，让长时间运行的重复也能看到进度，无需等它跑完。该心跳只打印到控制台，不写入 summary 文件。

**在真正信任集群上自动估算出的并行池规模之前**，`mem_per_worker_gb`（默认 `2.5` GB）**并未针对任何具体配置做过实测标定**。建议先单独跑一次，用 `/usr/bin/time -v matlab -batch ...` 或 Linux 下的 `top`/`ps` 测出真实的内存峰值，再把 `mem_per_worker_gb` 设为该测量值再加 30%~50% 余量。这个值预算过低，是导致 worker 被 Linux OOM Killer 中途杀死的最可能根因——由于 `parfor` 没有内置的超时或掉线 worker 恢复机制，这可能会让整批任务无限期卡死而不是报错退出。上面的重度输出抽样能降低触发概率，但不能完全杜绝。

## 6. 自定义几何模式

当胞体和纤毛都来自用户点云文件时，使用自定义几何模式。

1. 复制 [`configs/config_custom_example.json`](configs/config_custom_example.json)。
2. 设置 `"geometry_mode": "custom"`。
3. 在 `"custom_geometry"` 中提供胞体和纤毛文件。
4. 参考 [`geometry/`](geometry/) 中的示例文件格式。

### 自定义几何规则

- 自定义胞体和自定义纤毛必须同时使用。
- 胞体点云支持 `.csv`、`.txt` 和 `.mat` 文件，且应包含三维空间坐标。
- 导入的纤毛中心线允许原始点距不均匀。程序先删除连续重复点，再将 PCHIP 参数曲线转换为保持两端点的唯一等弦长 canonical 折线。
- 不支持逐纤毛 `bending_stiffness`；所有纤毛统一使用全局 `material.E` 与 `r_pili`，旧配置中的非空字段会被明确拒绝。
- 对任意初始姿态，优先使用 `initial_orientation`，不要依赖只绕单轴的 `alpha_0`。
- 自定义模式下 `deflection_z` 表示绕对齐后切向的扭转，`deflection_x` 表示倾斜；参数化模式下二者的几何角色相反。

## 7. 马达与自由扭转运行

建议从 [`configs/config_motor_flagella_example.json`](configs/config_motor_flagella_example.json) 开始。

每根纤毛的关键字段包括：

- `root_coupling`: `"clamped"` 或 `"free_torsion"`。
- `motor_torque.type`: `"off"` 或 `"constant"`。
- `motor_torque.torque_pN_nm`: 马达失速扭矩，单位 pN nm（带符号，见下）。
- `motor_torque.max_omega_Hz`: 线性扭矩–转速关系的零扭矩转速，单位 Hz。

在 `free_torsion` 模式下，根部扭转角 $\psi_1$ 会动态演化。其方程是**整根鞭毛关于根部**的轴向（$d_3$）力矩平衡：鞭毛所受水动力扭矩由胞体经根部传递的扭矩平衡，即马达扭矩（`motor_torque.type = "off"` 时为无摩擦轴承，扭矩为零）。内部 Kirchhoff 扭转力矩不出现在该行——它由内部节点承担（v8.4.0；此前内部扭转力矩也被加入根部行，使马达扭矩仅与第 1 段的阻力平衡，稳态转子转速被放大为分段数 $N$ 倍）。不会额外向 `R_global(4:6)` 加入反力矩，因为这些行表示完整“胞体 + 纤毛”耦合系统的角动量平衡。

稳态下直鞭毛的转子转速由整根鞭毛的扭转阻力决定：$\omega \approx \tau_{\mathrm{eff}} / (\xi_t L)$，其中 $\xi_t$ 为正则化 Stokeslet 线段模型的单位长度轴向转动阻力（默认 $\epsilon = \sqrt{2}\,r_{\mathrm{pili}}$（由 `material.hydro_radius_factor` 设定）时 $\xi_t = 2\pi\mu\epsilon^2 = 4\pi\mu r^2$，即无滑移圆柱值；8.5.1 及以前该因子为 1，$\xi_t$ 只有其一半。`validate_simulation_params.m` 的运行前失速诊断采用同一系数）。

**扭矩–转速关系（v8.2.0）。** 有效马达扭矩采用细菌鞭毛马达实测的线性扭矩–转速关系：

$$\tau_{\mathrm{eff}} = \tau_{\mathrm{motor}}\left(1 - \frac{\omega_{\mathrm{rel}}}{\omega_{\max}}\right),\qquad
\omega_{\mathrm{rel}} = \operatorname{sign}(\tau_{\mathrm{motor}})\,\big[(\dot\psi_1 + \cos\theta_1\,\dot\phi_1) - \boldsymbol{\Omega}_{\mathrm{body}}\cdot\mathbf{d}_3^{(1)}\big].$$

- $\omega_{\mathrm{rel}}$ 是转子（第 1 段）**相对胞体**的自旋速率，以马达驱动方向为正；它使用第 1 段完整的轴向角速度，而不是单独的 $\dot\psi_1$。
- 该关系是非对称的：$\omega_{\mathrm{rel}} = 0$ 给出失速扭矩，$\omega_{\mathrm{rel}} = \omega_{\max}$ 给出零扭矩，被迫反转（$\omega_{\mathrm{rel}} < 0$）时扭矩高于失速值，转速超过 $\omega_{\max}$ 时马达变为制动器（$\tau_{\mathrm{eff}} < 0$）。旧版的 $|\omega|$ 形式对两个转向一视同仁。
- 由于关系是线性的，自洽扭矩在每次 RHS 调用中通过一个小型线性方程组精确求得（任意数量的马达）。

**符号约定。** 正的 `torque_pN_nm` 使转子绕指向远离胞体的根部切向量 $\mathbf{d}_3$ 做右手旋转，即从鞭毛末端朝胞体看为逆时针（左手螺旋鞭毛的大肠杆菌“直行”转向）；负值为顺时针。v8.2.0 之前的配置采用相反且未记录的转向。

## 8. 布朗动力学

在 JSON 中设置 `brownian.enable`，或在 MATLAB 中设置 `cfg.enable_brown`。

Brownian 动力学默认为 body-only：6×6 有效迁移率 `M_eff` 来自完整胞体-纤毛确定性耦合系统，但热驱动只施加到胞体的三个平移和三个转动自由度；纤毛内部 Euler 角模态本身不接受热噪声。因此该方案不是完整耦合热平衡模型，不适用于研究纤毛内部模态的能量均分。选项 `brownian.pili_thermal`（v9.0.0，见下）为纤毛的慢弯曲模态加上热驱动。

**力实现方案（v8.5.0）。** 每个 Brownian 区间 `dt_B` 开始时探测 `M_eff`，抽取 Euler–Maruyama 增量 `dU = L_eff * eta * sqrt(dt_B)`，并将其换算为恒定的广义热力 `f_th = M_eff \ (dU / dt_B)`，在整个区间内加到 ODE 右端项的胞体行上。胞体由耦合求解器推动，纤毛通过水动力与弹性耦合连续跟随——与它们跟随重力或马达扭矩的方式完全一致。两点实际影响：

- 根部不再出现伪折点。旧版本在每个节点瞬时移动胞体，并只将第 1 段纤毛重新投影到移动后的胞体上，每次 kick 都在第 1、2 段之间注入曲率跳变，激发刚性高波数弹性模态，迫使 `ode15s` 在每个节点之后采用极小步长。改为力实现后根部约束从未被破坏，Brownian 算例不再需要更短的 `dt_B` 或更细的分段来保持光滑。
- 隐式 Ito 漂移。沿轨迹采样到的迁移率梯度在 `dt_B` 一阶精度上复现热漂移 `kBT*div(M_eff)`（见 `docs/theory.md` §10.2），因此不再计算显式有限差分漂移（v7.5.0–v8.4.0 每个节点额外 12 次耦合求解）。`[BROWN]` 日志在前 5 个区间及之后每 100 个区间起点报告抽取的 `|F_th|`、`|tau_th|` 与冻结迁移率下的目标增量。

SPD 修复和 Cholesky 分解在无量纲广义坐标中完成，再转换到 SI 单位。当 `T_obs` 不能被 `dt_B` 整除时，最后一个 Brownian 区间使用实际剩余时间。当前不提供用户可配置的中点漂移选项。固定 `brownian.rng_seed` 的算例仍可复现，但由于随机增量的实现方式不同，轨迹与 v8.4.0 及更早版本不一致。

重要选项包括：

- `brownian.temperature`: 正的流体绝对温度，单位 K；启用布朗运动时必需。
- `brownian.rng_seed`: 固定正整数种子，用于可复现运行。
- `brownian.dt_B_ratio`: 布朗时间步相对于输出时间步的比例。
- `brownian.pili_thermal`（`cfg.enable_pili_brown`，默认 `false`）：纤毛热涨落，见下。
- `brownian.pili_mode_orders`（`cfg.pili_brown_orders`，默认 `2`）：每根纤毛接受热噪声的弯曲阶数。
- `brownian.pili_tau_cut_s`（`cfg.pili_brown_tau_cut_s`，可选，默认不启用）：物理时间尺度门槛（单位秒）——自动保留所有松弛时间 `>= tau_cut_s` 的慢弯曲模态（推荐用于不等长多纤毛，默认开启偏振对保护）。
- `brownian.pili_protect_pairs`（`cfg.pili_brown_protect_pairs`，可选）：当截断处能隙比 `< 1.15` 时自动扩充至下一个清晰能隙，避免劈开近简并偏振对。
- `brownian.pili_rng_seeds`（`cfg.brown_pili_rng_seeds`，可选）：每根纤毛独立的正整数随机种子（与 `brownian.rng_seed` 互异）；在 `batch_repeat_brownian*.m` 批处理中自动为每次重复分配互不冲突的独立种子段。

**纤毛热涨落（v9.0.0，v10.0.0 全面加固升级）。** 设 `brownian.pili_thermal = true` 后，纤毛的慢弯曲模态（默认为每根纤毛第一、二阶弯曲，每阶两个偏振方向）也接受热力。这些模态在初始状态由耦合算子计算一次。每个区间内，胞体与热模态上的热力按联合迁移率一起抽取，并通过受约束探测的中心线力密度及平衡轴向标量力矩密度（`brown_pili_tau`）精确实现，保留胞体–纤毛相关且不向非热化/扭转模态泄露能量，每个热模态的弹性能松弛到 `kT/2`。扭转和更高阶弯曲模态仍为确定性。该选项要求 `brownian.enable = true` 且至少有一根纤毛。

- **在日志给出的双向窗口内选 `dt_B`。** 启动时的 `[BROWN]` 日志列出各热模态及最慢非热化弯曲模态（`tau_det`）的松弛时间。`dt_B` 上限不超过最快热模态松弛时间的四分之一（采用 `configs/config_example.json` 材料参数时，3 µm 纤毛取 0.002 s；`x_fast > 0.5` 时发出警告，`x_fast > 2.0` 默认拦截报错，除非显式设置 `cfg.allow_large_pili_dt_B = true`）；同时不应在固定阶数下把 `dt_B` 缩至远低于 `tau_det`（`dt_B / tau_det < 0.5` 时发出非绝热警告，如需更细步长请提高 `pili_mode_orders` 或降低 `pili_tau_cut_s`）。末段极短残余步长（`< 5% dt_B`）自动并入最后一个区间。
- **输出。** `results.pili_modes` 包含模态幅值 `a`（`energy_kT = a.^2 / (2 kT)`，平衡均值 1/2）、总弯曲能与总扭转能（`elastic_energy_kT`）、各模态松弛时间、阶数、`dom_pilus`、`pilus_order`、`selection_mode`、`pilus_share`，以及运行期最大残差诊断（`max_pattern_resid`、`max_hold_pattern_resid`、`max_twist_row_resid`、`max_hold_force`）。`results.flow_state.brown_mode_force`、`brown_pili_force` 与 `brown_pili_tau` 保存每个输出时刻的模态力、中心线力密度及平衡轴向力矩密度，用于精确重绘流场。
- **适用范围。** 本方案针对本征为直线、根部固支、无主动驱动的纤毛，并只在这类纤毛上验证过，且无重力、无背景流。`root_coupling = free_torsion`（可转动的鞭毛，无论有无马达）会直接报错；这类细胞请用 body-only 方案。本征曲率覆盖和主动力矩可以运行，但会给出警告：热模态取自被动的初始状态；中心线上的力只对本征为直线的纤毛才不引起扭转；这些情形下的能量均分没有验证过。使用前请检查 `elastic_energy_kT(2, :)`（扭转能）。抽取噪声时，比热模态更快的模态按已松弛处理，所以约 2 ms（3 µm 纤毛三阶弯曲的松弛时间）以内的胞体扩散不被解析。从时间间隔 `dt_B` 起，胞体扩散与线性理论一致；长时间下等于形状冻结的细胞的值（数值见 `docs/theory.zh-CN.md` §10.4）。把 `pili_mode_orders` 提到 3 需要约 0.0005 s 的 `dt_B`。

**滞后 Jacobian（v8.5.1）。** Brownian 运行在每个 `dt_B` 节点重启 `ode15s`，而每次重启过去都要做一次稠密有限差分 Jacobian（额外 `n_ode` 次右端项计算，占 4 菌毛算例耗时的 70–85%）。自 v8.5.1 起复用上一次的精确 Jacobian（[`src/core/LaggedJacobian.m`](src/core/LaggedJacobian.m)）：它只进入隐式步的牛顿迭代矩阵，因此解与此前结果的差异保持在积分容差以内，而每个 Brownian 节点的耗时降低约 2–3.7 倍（1 根菌毛约 1.4 倍）。当 `ode15s` 在同一次调用中再次索要 Jacobian（牛顿失败）、发生 Euler 基旋转，或缓存已达 50 个节点时，会重新计算精确 Jacobian；若某个滞后段报错或异常提前停止，则以默认 Jacobian 重算该段。`simulation.log` 的 `[JAC]` 行和 `results.jacobian_reuse` 给出计数（`n_fresh`、`n_lagged`、`n_on_demand`、`n_invalidated`、`n_fallback`）。仅运行时选项（不写入 JSON）：`cfg.jacobian_reuse`（默认 `true`；`false` 与 v8.5.0 逐位一致）和 `cfg.jacobian_max_age`（默认 50 个节点）。确定性算例不受影响。

启用布朗运动后，实际使用的随机种子会保存到 `config.json`、`brownian_rng_seed.txt`、`simulation.log` 和 `simulation_results.mat`。该种子只控制本次运行的独立随机流，因此几何缓存是否命中以及 MATLAB 全局随机流中的其他调用都不会改变 Brownian 轨迹。

**迁移率条件数告警（自 v7.4.1 起）：** 每次布朗节拍会在 [`src/physics/compute_effective_body_brownian.m`](src/physics/compute_effective_body_brownian.m) 内对已对称化的**无量纲**有效迁移率 `M_eff_dimless` 做 `rcond` 诊断；因此正常运行不应再只因 `M_eff_dim` 各分块的 SI 量纲差异而出现大量虚假的 `IllConditioned`。若 MATLAB 仍以极小的 **`M_eff_dimless` rcond** 抛出该告警，需排查欧拉根是否接近 \(\theta=\pi\)（钳制根部行退化）、体形极端纵横比等与几何/动力学相关的问题，并结合 `simulation.log` 与轨迹核对。

## 9. 主要输出

输出写入 `output/<case_name>/` 或自动生成的时间戳运行目录。

| 文件 | 用途 |
| --- | --- |
| `simulation_results.mat` | 完整 `results` 和 `params` 结构体。 |
| `timeseries.csv` | 便于外部分析的标量时间序列。 |
| `config.json` | 本次运行的配置快照，可用于复现。 |
| `summary_plot.png` / `.pdf` / 可选 `.html` | 主汇总图。 |
| `flowfield_panels.png` / `.pdf` / 可选 `.html` | 仿真结束时的 XZ、YZ 和 XY 末状态流场图（2×3）。 |
| `flowfield_panels_rerender_final.png` | 可选；由 `render_flowfield_from_results` 重绘的末状态流场图。 |
| `flowfield_panels_tNNNNNN.png` | 可选；指定时刻的流场图。 |
| `flowfield_animation.mp4` | 可选流场动画（Linux 上为 `flowfield_animation.avi`）。 |
| `initial_pose.png` / `.pdf` / 可选 `.html` | 开启预览时生成的初始姿态图。 |
| `bacteria_3d_motion.mp4` | 可选三维动画（Linux 上为 `bacteria_3d_motion.avi`）。 |
| `bacteria_3d_motion_rerender.mp4` | 可选；由 `render_animation_from_results` 写出（见下文）。 |
| `simulation.log` | 完整运行日志。 |
| `brownian_rng_seed.txt` | 启用布朗运动时写出的随机种子和复现提示。 |

动画由 `src/visualization/render_motion_animation.m` 渲染：所有图形对象只创建一次，逐帧仅更新数据。画面为固定 16:9 画布（2560 x 1440 逻辑像素；编码后的实际尺寸随显示缩放而变，例如 150% 缩放时为 3840 x 2160）上的四个面板：**a** 3D 总览（正交投影、固定范围、重力箭头），**b** XZ 总览（固定范围、2 μm 比例尺），**c** 跟随质心的 3D 窗口，**d** 跟随质心的 XZ 窗口。窗口 c、d 跟随质心并保持恒定放大倍数，且不是立方体或正方形：沿每个实验室坐标轴，窗口范围取所有帧中胞体和全部菌毛相对质心触及的范围再加 8% 余量（至少 0.1 μm），使细菌占据窗口主体，X、Y、Z 三轴比例相同（不拉伸、不压缩）；均带实验室坐标系下的刻度、刻度标签和网格（因此网格会从细菌旁边滑过）；面板 c、d 的刻度间隔不细于 0.5 μm。四个面板的排版自动选取：各面板的宽高比和相对大小由数据决定，程序在多种排版（并排/上下堆叠的嵌套分割，阅读顺序 a-b-c-d、a-c-b-d、a-c-d-b）中选出公共缩放系数最大的一种，并写入 `simulation.log`（例如 `Layout H(a,H(b,V(c,d)))`），因此细长、扁平或紧凑的细菌都不会在画面上留下大片空白。质心轨迹仅在 a、b 中绘制。胞体为均匀浅绿，腰部一侧有一小块深绿色（约占体长 18%，方位角 ±32°），随胞体一起转动，因此在各面板中都能看出绕长轴的滚转。画布上有图例、全局时间戳 *t* 和相对真实时间的播放速度（每视频秒对应的模拟秒数 = 帧间隔 x `animation_fps`，标注“slow motion”或“accelerated”），面板标题 a–d 紧贴各图像左上角。`render_motion_animation` 还接受 `opts.scale_bar_um`（默认 2），并在返回的 `info` 结构中给出 `window_half_um`（窗口半宽 `[hx hy hz]`）、`window_offset_um` 和 `playback_speed`。动画阶段结束时，`simulation.log` 会输出分阶段耗时表（坐标范围预扫描、图形搭建、数据更新、`drawnow`、`getframe`、帧缩放、`writeVideo`），用于判断渲染时间花在哪里。

若要在不重跑 ODE 的情况下根据已完成算例重新生成动画（例如更换抽帧间隔或帧率），先执行 `addpath('post_processing')`，再运行 `render_animation_from_results()`（弹出文件对话框）或 `render_animation_from_results(mat_file, fps, skip)`。未给出 `fps` / `skip` 时，取该算例文件夹 `config.json` 中的 `output.animation_fps` / `output.animation_skip`，否则为 30 和 1。结果保存在 `.mat` 文件旁，文件名为 `bacteria_3d_motion_rerender.mp4`（Linux 上为 `.avi`），不会覆盖原视频。

v8.7.0 及之后的结果可以用 `render_flowfield_from_results()` 或 `render_flowfield_from_results(mat_file, mode, requested_times_s, fps, skip)` 重算流场。`mode` 为 `final`（默认）、`times`、`animation` 或 `all`。时刻单位为秒，并匹配最近的已保存输出点；多个请求落到同一点时只计算一次。未给出 `fps` / `skip` 时，与运动动画一样读取 `config.json`。末状态色标只用该帧，指定时刻只用被请求的帧，动画使用全部已保存输出点；分位数计入 XZ、YZ 和 XY 三个平面。自 v8.7.3 起，图为 2500 x 1400 画布上的六面板：上排是 XZ、YZ、XY 的 3D 视图，下排是对应的 2D 视图，色标在子图 f 右侧。第六个参数 `opts.drop_first_frame = true` 时，`times`、`animation` 和 `all` 去掉已保存输出 1（`final` 不变）。自 v10.2.0 起，流场计算自动启用 Block Schur 力系数直接提取、三平面合批 FP64 GPU 速度核（`opts.use_gpu = 'auto'`，设 `false` 可强制 CPU）、线程池并行力求解与 2D 流线预计算（`opts.use_parallel = 'auto'`，设 `false` 为串行）以及轻量级网格缓存（~1.6 MiB/帧，2 GiB 可存约 1,300 帧）和逐帧编码后即时释放；`skip > 1` 时可在 `opts` 或 `config.json` 的 `flowfield` 段中设 `scale_from_video_frames = true` 以仅用视频帧计算色标分位数。动画在计算流场时显示一条原地刷新的 Compute 进度条，写视频时显示一条 Video 进度条。文件写在 `.mat` 旁边，不覆盖 `flowfield_panels.png`。动画在 Windows 和 macOS 上为 `flowfield_animation.mp4`，在 Linux 上为 `flowfield_animation.avi`。无图形界面时必须传入 MAT 路径，不能依赖文件选择框。

若某个积分段超时，`results.completed=false` 且 `results.termination_reason='timeout'`。程序会保存并正常返回已经获得的部分轨迹，但跳过动画、流场以及批处理稳态统计。

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

后处理脚本读取 `timeseries.csv` 与 `config.json`（可选从 `simulation_results.mat` 读取 `results.q` 以获得取向）。每个选中的 Case 文件夹对应一条轨迹，计算：

| 指标 | 方法 |
| --- | --- |
| 沉降速度 | \(-\texttt{Vsink\_um\_s}\) 的均值 / 末 5% 稳态 / 标准差 |
| 最大水平位移、三维曲折度 | 由 \(X,Y,Z\) 轨迹几何量 |
| 平移 \(D_\mathrm{trans}\)、\(V_\mathrm{msd}\) | 三维 TAMSD；在 \(0<\tau\le T/2\) 拟合 \(\mathrm{MSD}=6D\tau+V^2\tau^2\)（滞后覆盖整段 \(T\)） |
| 转动 \(D_\mathrm{rot}\) | XZ 体轴 MSAD；在 \(0<\tau\le T/2\) 拟合 \(\mathrm{MSAD}=2D_\mathrm{rot}\tau\) |

结果写入单个 Excel，并导出多面板汇总图（PNG、PDF；R2026a+ 可选 HTML）。超时部分轨迹（`results.completed=false`）不纳入扩散类汇总指标。

## 11. 常见检查

| 现象 | 建议检查 |
| --- | --- |
| MATLAB 找不到函数 | 确认仓库根目录已加入 MATLAB path。 |
| 胶囊体生成失败 | 确认 S2 Sampling Toolbox 已安装且在 MATLAB path 中。 |
| JSON 工况运行后没有纤毛 | 检查 `pili` 是否为合法 JSON 数组，且每根纤毛字段完整。 |
| 自定义几何被拒绝 | 检查点云方向、重复/退化区间，以及 canonical 等弦长构造是否可收敛。 |
| `geometry_mode` 混用报错 | 仅使用 parametric **或** custom 字段；勿将 `custom_geometry` 与 `body`/`pili` 矩阵字段混用（v7.4.4+）。 |
| `helix` 内禀曲率报错 | 须提供结构体：`{"type":"helix","kappa_d2":…,"kappa_d3":…}`，且为有限标量；裸 `"helix"` 无效（v7.4.4+）。 |
| 输出目录已存在 | 交互式选择处理方式，或设置 `auto_confirm = true` 使用无人值守重命名策略。 |
| 布朗重复结果不一致 | 确认种子模式为固定种子还是时间戳种子。 |
| 布朗重复批量未并行执行 | 安装 **Parallel Computing Toolbox**、检查 `parpool` 集群配置，或在 `batch_repeat_brownian.m` 中设置 `use_parallel = false`。 |
| JSON 文件夹批跑未并行执行 | 同上，针对 `batch_run_from_json.m`，或设置 `use_parallel = false`。 |
| JSON 批跑因重复 `case_name` 被跳过 | 非空 `case_name` 仅字母序首个 JSON 会运行；请使用唯一 `case_name` 或留空以自动时间戳目录。 |
| `IllConditioned` 布朗迁移率告警（`M_eff_dimless rcond`） | v7.4.1+ 在无量纲迁移率上做诊断；若仍出现告警，多半是真实的病态几何/欧拉根问题，核对 §8、`simulation.log` 及 \(\theta\approx\pi\) 配置。 |
| `NearSingularMatrix` 告警（`rcond < 1e-14`） | v8.0.1+：仅记录/告警，**不会单独中止**。中止需同时出现大后向残差（`> 1e-6`）→ `SingularMatrix`。查看 `results.diagnostics` / `[DIAG]` 日志。 |
| 警告 `numerical_diagnostics:DissipationEigensolveFallback` | 无害。耗散检查默认用稠密 `eig`；只有算子是稀疏矩阵且超过 600 行时才尝试迭代 `eigs`，不收敛时给出该警告并回退到 `eig`。不影响计算。 |
| `run_simulation:SilentPreviewDeadlock` 报错（v7.9.6+） | `cfg.silent = true` 且 `preview_initial = true` 时必须同时 `auto_confirm = true`（无头/并行 worker 无法应答交互式 `questdlg`）。设置 `auto_confirm = true` 或 `preview_initial = false`。 |
| Linux 上 `VideoWriter` 报错 `"The specified profile is not valid."` | 改用 [`run_simulation_HPC_Linux.m`](run_simulation_HPC_Linux.m)（§5.3）而非 `run_simulation.m`——`'MPEG-4'` 仅支持 Windows/Mac，HPC_Linux 分支会在 Linux 上自动切换为 `'Motion JPEG AVI'`。 |
| 批跑过程中 `top` 里 worker 进程数量突然变少、整批卡死不动（Linux） | 大概率是 Linux OOM Killer；`parfor` 没有内置超时或掉线 worker 恢复机制。在 `batch_repeat_brownian_HPC_Linux.m`（§5.3）上，用实测得到的真实内存峰值重新标定 `mem_per_worker_gb`（不要用默认值），和/或缩小 `heavy_output_sample_idx` 覆盖范围。可用 `dmesg \| grep -i "killed process"` 或 `journalctl -k \| grep -i "out of memory"` 确认。 |
| 家目录为分布式存储的集群上 `parpool` 启动很慢或直接失败 | `batch_repeat_brownian_HPC_Linux.m`（§5.3）的 `ensure_parpool` 已经会在 Linux 下把 `JobStorageLocation` 迁移到节点本地的 `tempdir()`；确认控制台是否打印了 `[INFO] parpool JobStorageLocation (node-local): ...`。如果打印的是迁移失败的警告，说明回退到了 MATLAB 默认路径，而这在家目录为分布式/网络存储的集群上是已知的启动慢/竞态风险点。 |

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
