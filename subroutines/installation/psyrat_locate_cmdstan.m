function home = psyrat_locate_cmdstan()
%Single CmdStan discovery point for the toolbox (G43).
%
%home = psyrat_locate_cmdstan()
%
%Background:
% The repo used to carry FIVE divergent CmdStan discovery orders: the bundled
% MatlabStan resolver (+mstan/stan_home.m: env vars > toolbox-managed home
% installs > MATLAB path), psyrat_detect_cmdstan_version (path > env, with
% the env names in reversed order and no home globs), the
% psyrat_checkversionsofdeps version read (path-only), the repair template
% inside psyrat_computevarcomp (which() fallbacks ahead of the home globs),
% and -- found during the G43 fix itself -- psyrat_find_cmdstan_home, the
% helper that injects a stan_home argument at the stan() call (no STAN_HOME
% in its env list, no ~/cmdstan-* glob, lexicographic glob order). A setup
% that only one of them could see produced inconsistent behavior -- most
% visibly, an environment-variable-only CmdStan install that estimation
% could use perfectly well still hit the startup install dialog and a
% spurious version warning.
%
% This function is the one discovery point the toolbox-side callers consult
% (the psyrat_start gate, psyrat_detect_cmdstan_version,
% psyrat_checkversionsofdeps, and psyrat_find_cmdstan_home's first tier
% before its legacy fallbacks). It DELEGATES to mstan.stan_home -- the
% resolver estimation actually uses at fit time -- rather than carrying yet
% another copy of the order, so what the gate and the version checks see is
% by construction what a fit would use. One known limit rides on that
% delegation: which('mstan.stan_home') answers for whichever MatlabStan is
% first on the path, and an UNPATCHED upstream copy (no env tiers, no
% globs, no sort) shadowing the bundled one would silently change the order
% -- the PSYRAT_PATCHLEVEL.txt stamp documents that two copies can be
% pathed at once, and psyrat_matlabstan_bundleinfo is the detector for it.
% The canonical order and the validity criterion live in
% +mstan/stan_home.m:
%   1. Environment variables (PSYRAT_TEST_CMDSTAN_PATH, STAN_HOME,
%      CMDSTAN_HOME, CMDSTAN).
%   2. Toolbox-managed install dirs under the user home
%      (PsyRATDependents/cmdstan-*, then ~/cmdstan-*), newest parseable
%      version first within each location (G43 ruling).
%   3. Whatever is discoverable on the MATLAB path.
% A candidate counts only when it is a VALID CmdStan root (a makefile plus a
% bin/ directory), so a stale environment variable pointing at a non-CmdStan
% directory reads as "not found" here rather than as a present-but-broken
% install.
%
%Required Inputs:
% No inputs are required.
%
%Outputs:
% home - the discovered CmdStan home directory as a char row, or '' when no
%  valid install was found (including when the bundled MatlabStan is not on
%  the MATLAB path -- estimation could not run in that state either, so
%  "not found" is the honest answer).

% Copyright (C) 2016-2026 Peter E. Clayson
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program (gpl.txt). If not, see
%     <http://www.gnu.org/licenses/>.
%

home = '';

%mstan.stan_home is a package function inside the bundled (patched)
%MatlabStan; which() is the reliable existence test for package functions.
if isempty(which('mstan.stan_home'))
    return;
end

try
    home = mstan.stan_home();
catch
    %stan_home is written to return '' rather than error, so a throw means a
    %corrupted file (the psyrat_repair_mstan_stan_home_file class of damage);
    %report "not found" and let the estimation-time repair handle the file.
    home = '';
end

if isstring(home)
    home = char(home);
end
if ~ischar(home)
    home = '';
end
end
