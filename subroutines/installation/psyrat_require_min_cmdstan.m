function ok = psyrat_require_min_cmdstan(minver,mode)
%Preflight check that the installed CmdStan meets PsyRAT's minimum version.
%
%psyrat_require_min_cmdstan(minver,mode)
%
%
%Background:
% The Stan models PsyRAT generates use the array[N] T declaration syntax that
% was introduced in Stan 2.26. CmdStan older than 2.26 cannot compile them and
% fails with a cryptic Stan syntax error partway through compilation. This
% function detects the installed CmdStan version and surfaces a clear,
% actionable message before that happens.
%
% The guard is intentionally fail-open: if the CmdStan version cannot be
% located or parsed, it does nothing (returns true) rather than blocking an
% otherwise valid run. It only acts when a version below the minimum is
% positively detected.
%
%Required Inputs:
% No inputs are required.
%
%Optional Inputs:
% minver - minimum acceptable CmdStan version as 'major.minor.patch'
%  (default '2.26.0').
% mode - 'error' (default) to raise an error when CmdStan is too old, or
%  'warn' to issue a warning instead (used at startup so the GUI still opens).
%
%Output:
% ok - logical scalar; true if the installed CmdStan meets the minimum or the
%  version could not be determined, false if it was found to be too old.

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

if nargin < 1 || isempty(minver)
    minver = '2.26.0';
end
if nargin < 2 || isempty(mode)
    mode = 'error';
end

%default to not blocking: only a positively-detected old version trips it
ok = true;

verstr = psyrat_detect_cmdstan_version();
if isempty(verstr)
    %could not determine the version; do not block (fail-open)
    return;
end

%compare major.minor.patch via the shared gate. psyrat_version_ge is fail-open
%(returns true when either version is unparseable), so an unreadable version
%leaves isbelow false and ok true, matching the prior fail-open behavior.
isbelow = ~psyrat_version_ge(verstr,minver);

if isbelow
    ok = false;
    %Build with sprintf, not strcat: strcat strips trailing whitespace from
    %char arguments, which would delete the spaces before the interpolated
    %version strings (yielding 'CmdStan2.25.0'). The %% / escaped text below is
    %literal here, so passing the finished string to warning/error is safe.
    msg = sprintf(['Detected CmdStan %s, but PsyRAT requires CmdStan %s or ' ...
        'newer.\n\n' ...
        'The generated Stan models use the array[N] T declaration syntax\n' ...
        'introduced in Stan 2.26, which older CmdStan cannot compile. How to\n' ...
        'fix: install CmdStan %s or newer (psyrat_installdependents\n' ...
        'installs the tested version), then restart MATLAB.'],verstr,minver,minver);
    %pass msg as a literal %s argument: it is already fully formatted, so this
    %avoids warning/error re-interpreting any character as a format specifier.
    if strcmpi(mode,'warn')
        warning('PsyRAT:cmdstanVersion','%s',msg);
    else
        error('PsyRAT:cmdstanVersion','%s',msg);
    end
end

end
