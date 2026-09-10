function versions = psyrat_checkversionsofdeps
%Check the versions of the dependents for PsyRAT Toolbox
%
%

%This function just will examine whether the dependents are up to date.
%MatlabStan carries no version information in its own files; PsyRAT stamps its
%vendored bundle instead, and psyrat_matlabstan_bundleinfo reads that stamp.
%
%Input
% No inputs required in the command line
%
%Output
% versions - array with information regarding whether versions of the
% dependents are 0-old, 1-up to date, or 2-untested new versions.
% -1 means the version could not be determined; only matlabstan returns it,
% when no MatlabStan is on the path or its stamp cannot be parsed.
% 3 means the in-use copy is PsyRAT's own in-tree bundle with a would-be-old
% stamp (G56d); only matlabstan returns it, and it deliberately fails every
% "== 0" needs-update check because a dependents update cannot replace the
% bundle. See psyrat_matlabstan_bundleinfo for the full scale.
%   cmdstan - information for cmdstan
%   matlabstan - information for matlab stan
%   matlabprocessmanager - information for matlab process manager
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
versions.cmdstan = [];
versions.matlabstan = [];
versions.matlabprocessmanager = [];

%pull the current versions of the dependents used by the toolbox
depvers = psyrat_dependentsversions;

%start off checking the version of cmdstan. Discovery goes through the single
%toolbox resolver (G43): psyrat_locate_cmdstan delegates to mstan.stan_home,
%the resolver estimation actually uses, so an env-var-only or managed-home
%install is found here the same way a fit would find it -- the old path-only
%which() lookups made exactly those setups trip the spurious "Unable to
%determine CmdStan version" warning below. The cmdstan-guide.tex read stays
%as a legacy fallback (pre-~2.22 installs record the version there) for
%installs the resolver cannot see.
cs_home = psyrat_locate_cmdstan;
cs_filepath = which('cmdstan-guide.tex');

if isempty(cs_home) && ~isempty(cs_filepath)
    %load the file and pull the version information
    cs_file = fileread(cs_filepath);
    C = strsplit(cs_file,'\');
    str = C{5};
    c_beg = strfind(str,'{');
    cs_ver = str(c_beg+1:end-2);
    
else
    %cs_home from the resolver above; when it is empty (and no guide.tex was
    %found either) the command below fails and the tested-version fallback
    %with its warning still applies, same as before. Quote the path and
    %handle the Windows .exe the same way psyrat_detect_cmdstan_version
    %does -- the adversarial wave caught the two centralized readers
    %disagreeing here, which let a space-bearing home pass the startup gate
    %and still trip the warning below.
    stanc = fullfile(cs_home,'bin','stanc');
    if ispc && exist([stanc '.exe'],'file') == 2
        stanc = [stanc '.exe'];
    end
    command = ['"' stanc '" --version'];
    cs_ver = '';
    stdout_text = '';

    try
        p = processManager('id','stanc version','command',command,...
            'keepStdout',true,...
            'printStdout',false,...
            'pollInterval',0.005);

        %Give the command enough time to emit version text.
        p.block(1);

        if isprop(p,'stdout') && ~isempty(p.stdout)
            stdout_text = strjoin(p.stdout,' ');
        end
    catch
        stdout_text = '';
    end

    %Fallback for environments where processManager stdout is empty.
    if isempty(strtrim(stdout_text))
        [~, cmdout] = system(command);
        stdout_text = cmdout;
    end

    vers_match = regexp(stdout_text,'\d+\.\d+\.\d+','match','once');
    if ~isempty(vers_match)
        cs_ver = vers_match;
    else
        %If version cannot be parsed, assume tested version to avoid hard
        %failure during startup.
        cs_ver = depvers.cmdstan;
        warning('Unable to determine CmdStan version from stanc output; using expected version for compatibility check.');
    end
end

%format version information for cmdstan that is used by the toolbox
cs_used_parts = sscanf(depvers.cmdstan,'%d.%d.%d')';

%format version information for cmdstan that is found in the matlab path
cs_found_parts = sscanf(cs_ver,'%d.%d.%d')';

%compare cmdstan versions
for ii = 1:3
    if cs_used_parts(ii) > cs_found_parts(ii)
        versions.cmdstan = 0;
        break;
    elseif cs_used_parts(ii) < cs_found_parts(ii)
        versions.cmdstan = 2;
        break;
    elseif cs_used_parts(ii) == cs_found_parts(ii)
        versions.cmdstan = 1;
    end
end


%check the version of MatlabStan

%MatlabStan records no release version in any of its own files, so it cannot be
%checked the way CmdStan and MatlabProcessManager are. PsyRAT therefore stamps
%its vendored bundle with PSYRAT_PATCHLEVEL.txt, which carries both the upstream
%release and a PsyRAT patch level, and the stamp is read from whichever
%MatlabStan directory Matlab would actually load. This replaces the -1 "not
%useful" sentinel that used to be returned here.
%
%A directory with no stamp is reported as OLD rather than unknown: it predates
%stamping and therefore cannot contain PsyRAT's local patches to MatlabStan.
%-1 is still returned when no MatlabStan can be found at all, or when a stamp is
%present but unparseable.
msinfo = psyrat_matlabstan_bundleinfo;
versions.matlabstan = msinfo.status;


%check the version of MatlabProcessManager
%find the file that has the information about the version
mpm_filepath = which('processManager.m');

%load the file and pull the version information
mpm_file = fileread(mpm_filepath);
C = strsplit(mpm_file,'version = ');
str = C{2};
c_beg = strfind(str,'''');
mpm_ver = str(c_beg(1)+1:c_beg(2)-1);

%format version information for cmdstan that is used by the toolbox
mpm_used_parts = sscanf(depvers.matlabprocessmanager,'%d.%d.%d')';

%format version information for cmdstan that is found in the matlab path
mpm_found_parts = sscanf(mpm_ver,'%d.%d.%d')';

%compare cmdstan versions
for ii = 1:3
    if mpm_used_parts(ii) > mpm_found_parts(ii)
        versions.matlabprocessmanager = 0;
        break;
    elseif mpm_used_parts(ii) < mpm_found_parts(ii)
        versions.matlabprocessmanager = 2;
        break;
    elseif mpm_used_parts(ii) == mpm_found_parts(ii)
        versions.matlabprocessmanager = 1;
    end
end


end
