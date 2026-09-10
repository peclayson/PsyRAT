function psyrat_printversions
%Print the versions of Matlab, PsyRAT Toolbox, and dependents to command
%window
%
%

%This function just will examine whether the dependents are up to date. It
%currently does not work for MatlabStan because the version info is not
%contained in any of the downloaded files.
%
%Input
% No inputs required in the command line
%
%Output
% The versions are printed in the command window
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

versions = struct;
versions.matlab = [];
versions.PsyRAT = [];
versions.cmdstan = [];
versions.matlabstan = [];
versions.matlabprocessmanager = [];

%store verion of Matlab
versions.matlab = version;

%store version of the PsyRAT Toolbox
versions.PsyRAT = psyrat_defineversion;

%check the version of cmdstan
%CmdStan stopped shipping cmdstan-guide.tex around 2.22, so which() returns
%empty on any modern install and the unguarded fileread that used to be here
%errored - taking the whole report down before it printed anything, including
%the versions that WERE available. psyrat_checkversionsofdeps handles this by
%asking stanc directly (it guards the same which() call); that logic already
%exists as psyrat_detect_cmdstan_version, so call it rather than keep a third
%copy of the parse. The .tex read stays as the fallback for genuinely old
%installs, which is the case stanc detection does not cover.
versions.cmdstan = psyrat_detect_cmdstan_version();

if isempty(versions.cmdstan)
    cs_filepath = which('cmdstan-guide.tex');
    if ~isempty(cs_filepath)
        %load the file and pull the version information
        cs_file = fileread(cs_filepath);
        C = strsplit(cs_file,'\');
        str = C{5};
        c_beg = strfind(str,'{');
        versions.cmdstan = str(c_beg+1:end-2);
    else
        versions.cmdstan = 'unknown';
    end
end


%check the version of MatlabStan

%MatlabStan carries no version of its own, so PsyRAT reads the stamp it puts
%in its vendored bundle. Reported from whichever MatlabStan Matlab would
%actually load, so a stale copy on the path shows as stale here.
msinfo = psyrat_matlabstan_bundleinfo;
if msinfo.found
    versions.matlabstan = sprintf('%s (PsyRAT patch level %s)', ...
        msinfo.version,msinfo.patchlevel);
elseif isempty(msinfo.dir)
    versions.matlabstan = 'not found';
else
    versions.matlabstan = 'unstamped (predates PsyRAT patch levels)';
end

%check the version of MatlabProcessManager
%find the file that has the information about the version
mpm_filepath = which('processManager.m');

%load the file and pull the version information
mpm_file = fileread(mpm_filepath);
C = strsplit(mpm_file,'version = ');
str = C{2};
c_beg = strfind(str,'''');
versions.matlabprocessmanager = str(c_beg(1)+1:c_beg(2)-1);

%print information to command window
fprintf('\nMatlab version: %s', versions.matlab);
fprintf('\nPsyRAT Toolbox version: %s', versions.PsyRAT);
fprintf('\nCmdStan version: %s', versions.cmdstan);
fprintf('\nMatlabStan version: %s', versions.matlabstan);
fprintf('\nMatlab Process Manager version: %s\n\n',...
    versions.matlabprocessmanager);

end

