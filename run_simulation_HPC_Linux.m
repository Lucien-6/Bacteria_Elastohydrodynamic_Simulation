function [results, params, output_dir] = run_simulation_HPC_Linux(cfg_or_path)
%RUN_SIMULATION_HPC_LINUX  HPC/Linux entry alias delegating to run_simulation.
%
%  Description:
%    Backward-compatible entry point retained for HPC/Linux cluster batch
%    drivers (e.g. batch_sweep/batch_repeat_brownian_HPC_Linux.m).  All
%    Linux/headless features formerly maintained in a separate fork —
%    platform-aware AVI/MP4 video encoding (in render_motion_animation.m),
%    .avi existing-output detection, cfg.silent, cfg.brown_pili_rng_seeds,
%    and cfg.progress_queue/cfg.progress_tag DataQueue heartbeats — are
%    unified in run_simulation.m so there is a single source of truth.
%
%  Syntax:
%    [results, params, output_dir] = run_simulation_HPC_Linux(cfg)
%    [results, params, output_dir] = run_simulation_HPC_Linux(json_path)
%
%  Inputs / Outputs:
%    Identical to run_simulation.  See: help run_simulation
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-03
%  Updated: 2026-11-24 — unified with run_simulation.m (thin forwarding wrapper).

[results, params, output_dir] = run_simulation(cfg_or_path);

end
