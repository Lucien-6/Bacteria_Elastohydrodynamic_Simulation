function info = call_saved_flowfield(results, params, output_dir, opts, engine_dir)
%CALL_SAVED_FLOWFIELD  Call the src renderer without shadowing this folder.
%
%  The post-processing entry has the same name as the computational
%  renderer.  Changing to the renderer directory makes that file the one
%  resolved by the shared name.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-28

old_dir = pwd;
cleanup_dir = onCleanup(@() cd(old_dir));
cd(engine_dir);
info = render_flowfield_from_results(results, params, output_dir, opts);
end
