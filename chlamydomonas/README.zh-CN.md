# 莱茵衣藻双鞭毛三维游动模拟模块

[English](README.md) | 中文

> **作者：** Lucien ([lucien-6@qq.com](mailto:lucien-6@qq.com))  
> **创建日期：** 2026-11-24  
> **适用引擎：** `Bacteria Elastohydrodynamic Simulation`（Kirchhoff 弹性细杆 + 正则化 Stokeslet 耦合求解器）

本模块依托项目主体的低雷诺数流固耦合引擎，建立并标定了野生型莱茵衣藻（*Chlamydomonas reinhardtii*，CC-124 / CC-125）的双鞭毛三维自由游动模型。方案在实测平面步态反演驱动的基础上，引入**手性拍平面倾角 + 近端 1ω/2ω 零寄生扭矩面外主动弯矩 + 5% 生理振幅不对称**，无需穿体接触即可复现活细胞典型轴向自旋（`1–2 Hz`）、左手超螺旋轨迹、三维非平面尖端闭环轨道与近端鞭毛瞬时 3D 挠率。

---

## 1. 目录结构与主工程衔接关系

本模块与主工程保持统一的职责分层（`configs/`、`geometry/`、`src/`、`output/` 共享主干，`chlamydomonas/` 内部按数据、核心物理、一键案例与绘图清晰隔离）：

```text
<project_root>/
├── configs/
│   └── config_chlamy_biflagellate_gz.json   # 权威双鞭毛配置（5×4 µm 椭球体，13° 根部附着）
├── geometry/
│   ├── chlamy_body_spheroid.csv             # 长椭球胞体表面点云（400 个 Riesz 采样点）
│   ├── chlamy_flagellum_ccw_gz.csv          # +Y 侧鞭毛实测平均静息形状（L = 11.70 µm）
│   ├── chlamy_flagellum_cw_gz.csv           # -Y 侧鞭毛镜像静息形状
│   └── generate_chlamydomonas_geometry.m    # 实测几何 CSV 生成函数
└── chlamydomonas/
    ├── README.md                            # 英文文档
    ├── README.zh-CN.md                      # 本文档（中文）
    ├── chlamy_run_preset.m                  # 统一预设运行入口（三维基准与平面参考）
    ├── data/
    │   └── fit_1p00A.mat                    # 受版本控制的双波数 1ω+2ω 标定驱动表（57.75 Hz）
    ├── core/
    │   ├── chlamy_build_params.m            # 参数与初值构建器（与 run_simulation 一致）
    │   ├── chlamy_inverse_moment.m          # 逆问题求解器：实测步态 -> 主动弯矩 q(s,t)
    │   ├── chlamy_fit_moment.m              # 双波数 1ω+2ω 解析谐波拟合
    │   ├── chlamy_make_perp_envelope.m      # 三维近端 1ω+2ω 面外弯矩构造器
    │   ├── chlamy_table_moment_fun.m        # 自定义主动矩句柄工厂（值捕获 par，防存档失效）
    │   ├── chlamy_table_moment_eval.m       # 独立主动矩评估器（节点局部轴、零寄生扭矩）
    │   ├── chlamy_run_case.m                # 通用无头算例运行器（支持全参数覆写）
    │   └── chlamy_swim_analysis.m           # 解缠 SO(3) 旋转与多圈螺旋/3D 轨道后处理分析
    ├── cases/
    │   ├── run_3d_helix_2s.m                # 一键运行：三维超螺旋 2.0 s 长时案例（~116 拍）
    │   ├── run_planar_sym50.m               # 一键运行：纯平面全对称 50 拍参考案例
    │   └── run_planar_a75f75.m              # 一键运行：纯平面振幅/频率不对称对照（1.0 s）
    ├── plotting/
    │   ├── plot_measured_gait_1p00A.py      # 绘制 1p00 A 实测实验步态图
    │   └── plot_gait_comparison_sym50.py    # 对比 SYM50 模拟与实测步态
    └── Refers/                              # 参考文献与 Dryad 实验波形原始数据
```

**与主求解器 `src/` 的解耦机制**：主求解器本身不含任何衣藻硬编码。`chlamy_run_case.m` 读取 `configs/config_chlamy_biflagellate_gz.json`，通过 `custom_geometry.pili(k).deflection_z` 施加根部拍平面倾角，并用 `type = 'custom'` 句柄（指向 `chlamy_table_moment_eval.m`）替换 JSON 中的占位行波项。所有运行结果 `simulation_results.mat` 均自带参数值捕获，在新会话中可经 `src/visualization/private/restore_saved_moment_handles.m` 直接加载后处理。

---

## 2. 快速上手

建议在项目根目录下运行（或直接在 MATLAB 中打开 `chlamydomonas/cases/` 下任意脚本点 **Run**，脚本会自动定位项目根目录）。

### 方式一：一键案例脚本

```matlab
% 三维螺旋游动长时案例（2.0 s ~ 116 拍，约 2.2 个完整螺旋周期，含 3D 动画）
% 输出目录：output/Chlamy-3D-Helix-2s/
run('chlamydomonas/cases/run_3d_helix_2s.m')

% 纯平面全对称 50 拍参考案例
run('chlamydomonas/cases/run_planar_sym50.m')
```

### 方式二：统一预设入口（`chlamy_run_preset`）

```matlab
addpath('chlamydomonas');

% 短时高精度 3 拍基准（RelTol = AbsTol = 1e-4）
s_sym = chlamy_run_preset('SYM3D');       % C2 对称三维蛙泳（直线前进 + ~1.17 Hz 自旋）
s_hlx = chlamy_run_preset('HELIX3D');     % 5% 左右振幅差（左手超螺旋）

% 长时三维螺旋 + MP4 动画
s_20b = chlamy_run_preset('HELIX3D_20');  % 20 拍（约 0.35 s）
s_2s  = chlamy_run_preset('HELIX3D_2S');  % 2.0 s（约 116 拍，约 2.2 圈完整螺旋）

% 干跑静态验证（~3 秒检查几何、驱动句柄与初始 ODE 右端项，不跑积分）
diag = chlamy_run_preset('HELIX3D_2S', [], struct('dry_run', true));
```

### 方式三：自定义参数扫参（`chlamy_run_case`）

```matlab
addpath('chlamydomonas/core');
S = load('chlamydomonas/data/fit_1p00A.mat', 'fit');
tab3d = chlamy_make_perp_envelope(S.fit, struct( ...
    'r_perp', 0.25, 'r_perp2', 0.05, 's_w', 0.20, ...
    'delta_perp', pi/2, 'delta_perp2', pi/2));

ov = struct( ...
    'case_name',      'My-Chlamy-3D-Test', ...
    'moment_table',   tab3d, ...
    'tilt_deg',       [17, 17], ...
    'perp_sign',      [1, -1], ...
    'drive_amp',      [1, 0.95], ...
    'T_obs',          3 / S.fit.f_Hz, ...
    'rel_tol',        1e-4, ...
    'abs_tol',        1e-4);

[results, params, out_dir] = chlamy_run_case('configs/config_chlamy_biflagellate_gz.json', ov);
s = chlamy_swim_analysis(results, struct('f_Hz', S.fit.f_Hz, 'n_skip', 1));
```

### 长时间仿真实时进度查看

`run_simulation.m` 会每 10 秒向 `output/<case_name>/progress.txt` 写入一行精确的实时进度（含百分比、当前无量纲时间与 ETA）：

```text
pct=42.50 t=1.4536e+01 T=3.4203e+01 elapsed=14820.1 eta=20050.7 seg=1
```

在 PowerShell 中随时查看（不影响 MATLAB 求解）：

```powershell
Get-Content output\Chlamy-3D-Helix-2s\progress.txt
```

---

## 3. 物理与数学模型（三维摆动方案）

### 3.1 几何与弹性本构

- **胞体**：长椭球体，沿体轴 \(+X\)（前进方向）长半轴 \(c = 5.0\,\mu\text{m}\)，短半轴 \(a = 4.0\,\mu\text{m}\)，表面采用 400 点缓存 Riesz 采样网格。
- **双鞭毛**：两根固支 Kirchhoff 弹性杆，长 \(L_k = 11.70\,\mu\text{m}\)（`seg_per_um = 6` 对应 \(N_{\mathrm{seg}} = 71\) 段），水动力半径 \(r = 0.125\,\mu\text{m}\)，抗弯刚度 \(EI = 8.40 \times 10^{-22}\,\text{N}\cdot\text{m}^2\)（无量纲弯矩单位 \(EI/L_k = 7.18 \times 10^{-17}\,\text{N}\cdot\text{m}\)）。附着点极角与根部切线半角均为 \(13^\circ\)（基体间距 \(1.80\,\mu\text{m}\)）。
- **静息形状**：取水相自由游动细胞实测时间平均曲率轮廓（Qin et al. 2015 / Gutierrez et al. 2025，静态总转角约 \(124^\circ\)），以 `intrinsic_curvature = 'frozen'` 冻结在参考构型中。
- **本构符号**：杆方程残差为 \(R_{\mathrm{rod}} = \kappa_{\mathrm{dev}} - \mathbf{M}_{\mathrm{act}}\)。

### 3.2 双正交分量三维主动弯矩

<!-- markdownlint-disable MD049 -->
在第 \(k\) 根鞭毛（\(k \in \{1, 2\}\)）、平衡节点 \(j \in \{2,\dots,N_{\mathrm{seg}}\}\)（无量纲弧长 \(s = (j-1)/N_{\mathrm{seg}}\)，根节点 \(j=1\) 恒为 0）处，主动内弯矩矢量为：

\[
\mathbf{M}_k(s, t) = q_{\parallel,k}(s, t)\,\hat{\mathbf{n}}_k(t) + q_{\perp,k}(s, t)\,\hat{\mathbf{b}}_k(s, t)
\]

1. **倾斜拍平面法向 \(\hat{\mathbf{n}}_k(t)\)**：将参考法向 \(\hat{\mathbf{z}}\) 绕第 \(k\) 根鞭毛根部外切向 \(\hat{\mathbf{t}}_{0,k}\) 右手旋转 \(\gamma_k\) 得到，固定在根部局部坐标系并随胞体刚体旋转。
   - **重要符号约定（`tilt_deg = [+17, +17]`）**：两鞭毛原本左右镜像弯曲（分别朝 \(+Y\) 与 \(-Y\)），因此**同号**右手旋转 `[+17, +17]` 恰好把 +Y 侧鞭毛抬向 `+Z`（尖端 `+2.89 µm`）、把 -Y 侧鞭毛压向 `-Z`（尖端 `-2.89 µm`），对应 Cortese & Wan（*PRL* 2021）的 \(C_2\) 手性轨道倾角（\(\beta \approx 17^\circ\)）。若误写异号 `[+17, -17]`，两根鞭毛会同时翘向 `+Z`，净自旋抵消为零。
2. **节点局部面外弯曲轴 \(\hat{\mathbf{b}}_k(s, t)\)**：
   \[
   \hat{\mathbf{b}}_k(s, t) = \frac{\mathbf{d}_{3,k}(s, t) \times \hat{\mathbf{n}}_k(t)}{\|\mathbf{d}_{3,k}(s, t) \times \hat{\mathbf{n}}_k(t)\|}
   \]
   由于 \(\hat{\mathbf{b}}_k\) 在每个节点上都严格垂直于当前局部切向 \(\mathbf{d}_{3,k}\)，无论鞭毛弯曲角度多大，面外驱动都只产生**纯面外弯曲，严格零寄生扭转矩**（\(M_{3,k}^\perp \equiv 0\)）。

### 3.3 平面主驱动 \(q_{\parallel,k}(s, t)\)（双波数 1ω+2ω 测量拟合）

来自实测细胞 `1p00 A`（\(f = 57.75\,\text{Hz}\)）精确逆问题反演的紧致解析拟合（`chlamydomonas/data/fit_1p00A.mat`，单位 \(EI/L_k\)）：

\[
q_{\parallel,k}(s, t) = (-1)^{k-1}\Big\{ c_0(s) + \alpha_k\big[A(s)\cos(\omega t + \phi(s)) + B(s)\cos(2\omega t + 2\phi(s) + \delta(s))\big] \Big\}
\]

- \(A(s) = A_0(1 - d_A s) + A_{\mathrm{pk}}\exp\!\big(-((s - s_p)/w_p)^2\big)\)
- \(\phi(s) = \phi_0 - k_1 s - (k_2 - k_1)\tau\ln\!\big(1 + \exp((s - s_b)/\tau)\big)\)（近端波数 \(k_1 = 7.61\)、远端波数 \(k_2 = 3.97\,\text{rad}/L\)）
- \(B(s) = B_0 + B_{\mathrm{pk}}\exp\!\big(-((s - s_p)/w_p)^2\big)\)，\(\delta(s) = d_0 + \dfrac{d_1 - d_0}{1 + \exp(-(s - s_d)/w_d)}\)
- \(c_0(s)\)：8 节点分段线性静态修正矩（不受 \(\alpha_k\) 缩放）。

### 3.4 三维面外驱动 \(q_{\perp,k}(s, t)\)

无静态偏置（\(c_{0,\perp} \equiv 0\)），相位锁相到主驱动的一阶相位 \(\phi(s)\)，采用近端指数包络：

\[
q_{\perp,k}(s, t) = \sigma_k\,\mathrm{e}^{-s/s_w}\Big[ A_\perp \cos(\omega t + \phi(s) + \delta_\perp) + B_\perp \cos(2\omega t + 2\phi(s) + \delta_{\perp 2}) \Big]
\]

其中 \(A_\perp = r_\perp\, q_{\mathrm{rms}}\)，\(B_\perp = r_{\perp 2}\, q_{\mathrm{rms}}\)，主振荡均方根 \(q_{\mathrm{rms}} = 3.421\,EI/L_k\)（\(2.46 \times 10^{-16}\,\text{N}\cdot\text{m}\)）。

**手性单次乘法原则**：表内 `tab.coef_perp = {c1, c1}` 对两根鞭毛存完全相同的空间包络，cis/trans 反向手性 \(\sigma = [+1, -1]\) **只在运行时由 `ov.perp_sign` 乘一次**，确保对称驱动下严格满足 \(\mathbf{M}_2 = R_x(\pi)\mathbf{M}_1\) 的 \(C_2\) 旋转对称。
<!-- markdownlint-enable MD049 -->

### 3.5 三维方案标定参数汇总

| 类别 | 参数 | 取值 | 物理含义 |
| --- | --- | --- | --- |
| **基础尺度** | \(f\)、\(L_k\)、\(EI/L_k\) | `57.75 Hz`、`11.70 µm`、`7.18e-17 N·m` | 拍打频率、鞭毛长度、弯矩尺度 |
| **平面主驱动** | \(A_0, d_A, A_{\mathrm{pk}}, s_p, w_p\) | `2.997, 0.363, 12.578, 0.097, 0.046` | 一阶幅值分布（\(EI/L_k\)） |
| (`fit_1p00A`) | \(\phi_0, k_1, k_2, s_b, \tau\) | `1.410, 7.606, 3.975, 0.435, 0.030` | 双波数相位分布 |
| | \(B_0, B_{\mathrm{pk}}, d_0, d_1, s_d, w_d\) | `1.083, 2.334, 2.583, 5.831, 0.310, 0.048` | 二阶谐波幅值与相移 |
| **几何手性** | `tilt_deg` (\(\gamma_1, \gamma_2\)) | `[+17°, +17°]` | \(C_2\) 手性轨道倾角 |
| **面外驱动** | `r_perp`、`r_perp2` | **`0.25`**、**`0.05`** | 1ω（`6.14e-17 N·m`）与 2ω（`1.23e-17 N·m`）比例 |
| (**1ω+2ω**) | `s_w`、`delta_perp`、`delta_perp2` | **`0.20`**、**`+π/2`**、**`+π/2`** | 近端衰减长度与正交相移 |
| | `perp_sign` (\(\sigma_1, \sigma_2\)) | **`[+1, -1]`** | cis/trans 手性符号（与倾角同向增强自旋） |
| **超螺旋控制** | `drive_amp` (\(\alpha_1, \alpha_2\)) | `[1, 1]` 或 **`[1, 0.95]`** | 5% 生理振幅差 -> 左手超螺旋 |

---

## 4. 标定基准与实验观测对比

全套标定汇总表与四面板总结图存放在 `output/chlamy_oop/`（`oop_summary.csv`、`oop_step1_summary.csv`、`phase2_3d_summary.pdf/.png`）。

| 算例 | 配置要点 | `U_fwd` [µm/s] | `spin` [Hz] | 螺旋半径 `R` [µm] | 螺距 [µm] | 尖端 `z_pp` [µm] | 中心线 3D 非平面度 | `f_body_min` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **活细胞实验** | 野生型水相文献 | ~80–128（50 Hz） | **1–2 Hz** | ~1–5（左手） | ~50–100 | 三维闭环 | 近端手性扭转 | > 1 |
| `OOP0`（基线） | `tilt=[17,17]`，`r_perp=0` | 80.1 | -0.685 | 0.03 | -116.9 | 2.49 | 0.67% | 1.124 |
| `OOP10F`（对照） | `tilt=[0,0]`，`r_perp=0.10` | 81.1 | -0.195 | 0.18 | -414.8 | 0.62 | 0.54% | 1.124 |
| **`SYM3D`**（三维，`10⁻⁴`） | `drive_amp=[1, 1]` | **79.6** | **-1.166** | 0.00 | -68.2 | 3.54 | **1.88%** | **1.124** |
| **`HELIX3D`**（三维，`10⁻⁴`） | `drive_amp=[1, 0.95]` | **77.6** | **-1.092** | **1.49** | **-69.9** | 3.41 | **1.82%** | **1.124** |

关键验证结论：

- **时间积分收敛性**：`HELIX3D` 在 `tol = 1e-3` 与 `1e-4` 之间所有核心量偏差均 `< 0.7%`（`ΔU = 0.19%`、`Δspin = 0.09%`、`ΔR = 0.67%`）。
- **精确线性叠加**：纯几何倾角（`-0.685 Hz`）+ 纯动力学面外（`-0.195 Hz`）= `-0.880 Hz`，与联合仿真（`-0.875 Hz`）误差仅 `0.6%`；反转 `perp_sign` 为 `[-1, +1]` 精确抵消至 `-0.497 Hz`，同号 `[+1, +1]` 自旋零增益而全部转为胞体俯仰晃动。
- **2ω 项作用**：纯 1ω 空间闭曲线在数学上必共面，叠加 5% 二阶项（`r_perp2 = 0.05`）在不损失推进力和安全间隙的前提下，将尖端轨道真三维非平面度（`sv3/sv1`）提升了 `+46%`。

---

## 5. 重新生成几何与标定表（进阶用法）

仓库已随附可直接运行的几何 CSV 和 `chlamydomonas/data/fit_1p00A.mat`。如需从 `chlamydomonas/Refers/` 原始实验数据重新构建：

```matlab
% 1. 重新生成几何 CSV（写入 geometry/）
addpath('geometry');
generate_chlamydomonas_geometry();

% 2. 由实测步态反演主动弯矩 -> 双波数 1ω+2ω 解析拟合
addpath('chlamydomonas/core');
raw_mat = ['chlamydomonas/Refers/doi_10_5061_dryad_sbcc2frj6__v20250313/' ...
           'Quantification_of_Flagellar_Gait_Changes_Gutierrez_Data/Data/' ...
           'Newtonian 1p00 A_Viscosity1.00.mat'];
tab = chlamy_inverse_moment(raw_mat);
[fit, P] = chlamy_fit_moment(tab);
save('chlamydomonas/data/fit_1p00A.mat', 'fit', 'P');
```

---

## 6. 输出文件与后处理说明

每个预设会在 `output/<case_name>/` 下生成标准求解器产物（`simulation_results.mat`、`timeseries.csv`、`config.json`、`simulation.log`、`progress.txt`，可选 `bacteria_3d_motion.mp4`）以及由 `chlamy_swim_analysis.m` 输出的 `swim_analysis.mat`。

- **长时多圈旋转无折叠**：`chlamy_swim_analysis.m` 通过逐帧累加 \(\mathrm{SO}(3)\) 李代数旋转向量 \(\log(R_m^\top R_{m+1})\) 与中点体架位移投影来计算 `Omega_body_rad_s`、`spin_Hz`、`total_spin_turns` 和螺旋参数，对多秒多圈（如 `HELIX3D_2S` 的 2.0 s / 2.2 圈）完全精确，不会发生 `< π` 角度折叠。
- **间隙安全检查**：当前求解器不含接触排斥力，分析新算例时务必确认 `s.f_body_min > 1`（鞭毛不穿胞体）且 `s.d_flag_min_um > 0.25`（双鞭毛不重叠）。
