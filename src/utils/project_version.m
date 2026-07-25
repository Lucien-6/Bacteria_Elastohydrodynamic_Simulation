function v = project_version()
%PROJECT_VERSION  Single source of truth for the project version string.
%
%  Description:
%    Returns the current project version as a char row vector.  Referenced
%    by run_simulation.m's log banner (and any other runtime code that
%    needs to report the version) so that a release bump only requires
%    editing this one file instead of hunting down every hardcoded
%    version string scattered across the codebase.
%
%    Documentation files (README.md, CHANGELOG.md, USER_GUIDE.md, etc.)
%    are still updated manually as part of the release process; this
%    function only covers strings emitted at runtime.
%
%  Syntax:
%    v = project_version()
%
%  Output:
%    v - char row vector, e.g. '8.1.0'
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-07-03
%  Updated: 2026-07-25

v = '8.1.0';

end
