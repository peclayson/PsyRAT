function depvers = psyrat_dependentsversions
%Specifies the versions used for the dependents of the PsyRAT Toolbox
%
%

%This function just sets the version numbers for CmdStan, MatlabStan, and
% MatlabProcessManager
%
%Input
% No inputs required in the command line
%
%Output
% depvers - versions for each dependent
%   cmdstan - version of CmdStan
%   matlabstan - version of MatlabStan
%   matlabprocessmanager - version of MatlabProcessManager
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

depvers = struct;
depvers.cmdstan = '2.38.0';
depvers.matlabstan = '2.15.1.0';

%PsyRAT patches its vendored MatlabStan, and MatlabStan itself records no
%version anywhere in its files, so the release string above cannot tell a
%patched bundle from an unpatched one. This is the patch level PsyRAT expects,
%compared against the stamp in bundled_dependents/MatlabStan-<ver>/
%PSYRAT_PATCHLEVEL.txt. Bump BOTH in the same commit; see that file for when.
depvers.matlabstanpatch = '2';
depvers.matlabprocessmanager = '0.5.1';

end

