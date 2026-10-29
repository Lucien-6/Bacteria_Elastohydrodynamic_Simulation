function info = open_saved_flowfield_ui(varargin)
%OPEN_SAVED_FLOWFIELD_UI  Reach the post-processing entry by folder lookup.
%
%  src/visualization is normally ahead of post_processing on the MATLAB
%  path.  This function changes to post_processing before using the shared
%  function name, so the file browser is the function that runs.
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-10-28

this_dir = fileparts(mfilename('fullpath'));
proj_root = fileparts(fileparts(fileparts(this_dir)));
post_dir = fullfile(proj_root, 'post_processing');
old_dir = pwd;
cleanup_dir = onCleanup(@() cd(old_dir));
cd(post_dir);
info = render_flowfield_from_results(varargin{:});
end
