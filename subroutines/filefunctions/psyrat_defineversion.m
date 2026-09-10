function psyratver = psyrat_defineversion
%Set version for the PsyRAT Toolbox
%
%

%This function just sets the version number for loading by various scripts
%
%Input
% No inputs required in the command line
%
%Output
% psyratver - PsyRAT Toolbox version
%

% Copyright (C) 2016-2025 Peter E. Clayson
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

%SINGLE SOURCE OF TRUTH for the toolbox version. Every consumer routes through
%this function -- the startup banner, psyrat_prefs.ver, the run-config sidecar,
%the export-table headers, the project writer and the GitHub release check --
%so this literal is the only place a version is declared.
%
%The '-beta' suffix is deliberate (B4 / W7). The bare '0.1.0' understated the
%toolbox: the ERA-lineage Version History reaches 0.5.2 (Oct 2020), and PsyRAT
%has since added the gamma/chi-square family, the dynamic/conditional designs
%and the person-specific DoD retest designs on top of that. '-beta' states that
%the API and outputs are still moving, which is the honest claim while the
%copula designs (analyses 7/8/10) ship experimental with convergence not
%validated in-repo (S11).
%
%SAFE FOR THE SEMVER COMPARISON, measured rather than assumed:
%sscanf('0.1.0-beta','%d.%d.%d') returns [0 1 0] (numel 3), so
%local_compare_semver in psyrat_checkgithubreleases parses it as 0.1.0 and its
%numel(...) ~= 3 guard does not fire. The same holds for the test helper
%localStepVersion. Keep any future suffix in the '<x.y.z>-<tag>' shape: a suffix
%that perturbs the leading x.y.z would trip 'psyrat:versionparse'.
psyratver = '0.1.0-beta';


end
