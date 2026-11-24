function info = call_saved_flowfield(results, params, output_dir, opts)
%CALL_SAVED_FLOWFIELD  Call the computational flow-field engine.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-28
%  Updated: 2026-11-18 — call render_flowfield_engine by name.

info = render_flowfield_engine(results, params, output_dir, opts);
end
