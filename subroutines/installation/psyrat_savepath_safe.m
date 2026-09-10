function ok = psyrat_savepath_safe()
%Persist the MATLAB search path, tolerating read-only/managed installations.
%
%psyrat_savepath_safe
%
%
%Background:
% The toolbox calls savepath after adding PsyRAT and its dependents to the
% MATLAB path so the additions persist across sessions. On centrally managed
% or admin-locked MATLAB installations (common on shared lab, cluster, or
% campus machines), pathdef.m is not writable, and savepath then fails. The
% stock behavior is to abort with an error even though the install otherwise
% succeeded. This helper downgrades that failure to an informative warning so
% the current session still works and the user is told how to make the path
% persist.
%
%Required Inputs:
% No inputs are required.
%
%Output:
% ok - logical scalar; true if the path was saved, false if it could not be
%  saved (a warning is issued in that case).

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

ok = false;

%savepath returns 0 on success and 1 on failure; it can also throw. Handle
%both so a read-only path never aborts an otherwise successful install.
try
    status = savepath;
    if status == 0
        ok = true;
    else
        warning('psyrat:savepath',...
            strcat('Unable to save the MATLAB search path (savepath reported a\n',...
            'failure). The dependents were added for this session but may not\n',...
            'persist after MATLAB restarts.\n',...
            'How to fix: run MATLAB with write access to its pathdef.m, or add\n',...
            'the dependents folder to your startup.m or userpath so it loads\n',...
            'automatically each session.'));
    end
catch err
    %Pass err.message as a %s argument, not inside the format string, so a
    %percent sign or backslash in the message cannot corrupt the warning.
    warning('psyrat:savepath',...
        strcat('Unable to save the MATLAB search path: %s\n',...
        'The dependents were added for this session but may not persist after\n',...
        'MATLAB restarts.\n',...
        'How to fix: run MATLAB with write access to its pathdef.m, or add the\n',...
        'dependents folder to your startup.m or userpath so it loads\n',...
        'automatically each session.'),...
        err.message);
end

end
