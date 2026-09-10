function psyrat_installdependents(varargin)
%Install the Matlab dependents. User will be prompted when input necessary
%
%psyrat_installdependents
%
%
%Required Inputs:
% No inputs are required.
% Some features only work for certain versions of Mac or Windows. When
%  automatic installation is not possible, the error text states what must
%  be installed manually.
%Optional Inputs:
% depvercheck - structure array with three dependents. The script will
%  update all of the dependents with a value of 0 (indicating that the
%  particular dependnet is out of date).
%
%Output:
% No data are outputted to the Matlab command window. However, software
%  will be installed and a new Matlab path will be saved.

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

if nargin > 0 && ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of 
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)  
        error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',... 
        'Make sure each variable input is paired with a value \n',...
        'See help psyrat_installdependents for more information on optional inputs'));
    end
    
    %check if the dataset is present
    ind = find(strcmp('depvercheck',varargin),1);
    if ~isempty(ind)
        depvercheck = varargin{ind+1}; 
    else 
        depvercheck = [];
    end
elseif nargin ~= 2
    depvercheck = [];
end

%load the versions of the dependents used by the toolbox
depvers = psyrat_dependentsversions;

%CmdStan is always fetched as the self-contained release tarball, which
%bundles the Stan/math submodules required to build it. The GitHub source
%archive (/archive/v*.zip) omits those submodules and cannot be built, so
%it is no longer used.
url_cs_tar = strcat('https://github.com/stan-dev/cmdstan/releases/download/v',...
    depvers.cmdstan,'/cmdstan-',depvers.cmdstan,'.tar.gz');
url_ms = strcat('https://github.com/brian-lau/MatlabStan/archive/v',...
    depvers.matlabstan);
url_mp = strcat('https://github.com/brian-lau/MatlabProcessManager/archive/v',...
    depvers.matlabprocessmanager);

urls = struct;
urls.cmdstan_tar = url_cs_tar;
urls.ms_zip = strcat(url_ms,'.zip');
urls.ms_tar = strcat(url_ms,'.tar.gz');
urls.mpm_zip = strcat(url_mp,'.zip');
urls.mpm_tar = strcat(url_mp,'.tar.gz');

%determine the version of OS that is being used
if ismac
    sys = 1; %because Apple is number 1, obviously
elseif ispc
    sys = 2; %because Windows OS is inferior
elseif isunix && ~ismac
    sys = 3; %Linux (EXPERIMENTAL: auto-install is not verified by the maintainer)
end

%check whether the version of the OS is supported
if sys == 1 %mac
    [~,cmdout] = system('sw_vers');
    parseout = strsplit(cmdout);
    parseOS = strsplit(parseout{6},'.');
    
    if isnan(str2double(parseOS{1})) || ~contains(char(parseOS{1}),'.')
        parseOS = strsplit(parseout{4},'.');
    end
    
    if str2double(parseOS{1}) >= 10
        if str2double(parseOS{2}) >= 9 || str2double(parseOS{1}) >= 11
            
            if isempty(depvercheck)
            
                [~, savepath] = uiputfile('PsyRATDependents',...
                'Where would you like to save the directory for the dependents?');

                %if the user does not select a file, then take the user back to psyrat_start    
                if savepath == 0
                    %relaunch first, dialog last (B20 creation order); modal per B21
                    psyrat_start;
                    errordlg('Location not selected','File Error','modal');
                    return;
                end

                wrkdir = fullfile(savepath,'PsyRATDependents');
                
            elseif ~isempty(depvercheck)

                wrkdir = fileparts(fileparts(which('runCmdStanTests.py')));

                %G56 (partial hardening). which() sees only the MATLAB
                %path, so wrkdir is empty whenever the resolvable CmdStan
                %lives off-path -- including the toolbox's OWN managed
                %install when a failed savepath dropped it from the saved
                %path (psyrat_savepath_safe downgrades that failure to a
                %warning). "No install at all" cannot arrive here: the
                %startup gate diverts to the fresh-install prompt when the
                %resolver finds nothing. An empty wrkdir is never usable
                %downstream: mkdir('') is a silent no-op (MEASURED R2025b,
                %status 1 / MATLAB:MKDIR:EmptyDirectoryName), the status
                %GUI renders, and the first hard failure is cd('') in the
                %download block (MATLAB:cd:NonExistentFolder, MEASURED),
                %with non-download copies resolving fullfile('',...)
                %relative to pwd. Stop BEFORE rmolddeps so a doomed update
                %cannot delete working dependents first, and exit by the
                %file's B20 convention (relaunch first, dialog last, modal
                %per B21) instead of erroring out of a callback that has
                %already closed every PsyRAT window. Deriving wrkdir from
                %the resolver instead of which() is part of the G56 design
                %ruling still open in the ledger.
                if isempty(wrkdir)
                    reshome = psyrat_locate_cmdstan();
                    if isempty(reshome)
                        whereline = 'No CmdStan install could be resolved.';
                    else
                        whereline = sprintf('CmdStan was found at: %s', reshome);
                    end
                    dlgtext = {...
                        'The update flow could not find the dependents directory.'; ...
                        whereline; ...
                        'That install is not on the MATLAB path, so there is no directory to update in place.'; ...
                        'Add that directory to the MATLAB path (addpath + savepath) and rerun the update,'; ...
                        'run psyrat_installdependents from the Command Window to choose a fresh install directory,'; ...
                        'or update an environment-variable (STAN_HOME/CMDSTAN) install manually.'};
                    %relaunch first, dialog last (B20 creation order); modal per B21
                    psyrat_start;
                    errordlg(dlgtext,'Update location not found','modal');
                    return;
                end

                %G56c: an env-var pin outranks any fresh download at fit
                %time; disclose the inert update and let the user cancel
                %BEFORE anything is removed.
                if ~local_envpin_update_gate(depvercheck)
                    return;
                end

                rmolddeps(depvercheck);

            end

            %the third argument carries the version verdict ([] outside the
            %update flow) -- the per-OS installers use it to say plainly
            %when a CmdStan flagged old survived removal and nothing will
            %be downloaded (G56)
            psyrat_macdepsinstall(wrkdir,urls,depvercheck);
            
        else
            error('mac:oldver',... %Error code and associated error
            strcat('WARNING: Automatic installation only works for OSX 10.9 and newer \n'));
        end
    else
        error('mac:oldver',... %Error code and associated error
        strcat('WARNING: Automatic installation only works for Mac OS 10.9 or newer \n'));
    end
elseif sys == 2 %windows
    
    if isempty(depvercheck)
    
        fprintf('Warning: automatic installation has only been tested using Windows 7\n');
        [~, savepath] = uiputfile('PsyRATDependents',...
        'Where would you like to save the directory for the dependents?');

        %if the user does not select a file, then take the user back to psyrat_start    
        if savepath == 0
            %relaunch first, dialog last (B20 creation order); modal per B21
            psyrat_start;
            errordlg('Location not selected','File Error','modal');
            return;
        end

        wrkdir = fullfile(savepath,'PsyRATDependents');
     
    elseif ~isempty(depvercheck)

        wrkdir = fileparts(fileparts(which('runCmdStanTests.py')));

        %G56 (partial hardening). Same guard as the macOS arm (see the
        %comment there): an empty which() means no in-place directory to
        %update; stop before rmolddeps, by the B20 exit convention.
        if isempty(wrkdir)
            reshome = psyrat_locate_cmdstan();
            if isempty(reshome)
                whereline = 'No CmdStan install could be resolved.';
            else
                whereline = sprintf('CmdStan was found at: %s', reshome);
            end
            dlgtext = {...
                'The update flow could not find the dependents directory.'; ...
                whereline; ...
                'That install is not on the MATLAB path, so there is no directory to update in place.'; ...
                'Add that directory to the MATLAB path (addpath + savepath) and rerun the update,'; ...
                'run psyrat_installdependents from the Command Window to choose a fresh install directory,'; ...
                'or update an environment-variable (STAN_HOME/CMDSTAN) install manually.'};
            %relaunch first, dialog last (B20 creation order); modal per B21
            psyrat_start;
            errordlg(dlgtext,'Update location not found','modal');
            return;
        end

        %G56c: an env-var pin outranks any fresh download at fit time;
        %disclose the inert update and let the user cancel BEFORE anything
        %is removed.
        if ~local_envpin_update_gate(depvercheck)
            return;
        end

        rmolddeps(depvercheck);

    end

    %the third argument carries the version verdict ([] outside the
    %update flow; see the macOS call above; G56)
    psyrat_windepsinstall(wrkdir,urls,depvercheck);
elseif sys == 3 %linux (EXPERIMENTAL - unverified by the maintainer)

    warning('os:linuxExperimental',... %not an error: Linux install is experimental
        strcat('EXPERIMENTAL: automatic dependency installation on Linux has not\n',...
        'been verified by the maintainer. It mirrors the macOS path (download the\n',...
        'CmdStan release tarball and build it with make/g++). If it fails, install\n',...
        'CmdStan manually and point the toolbox at it with the STAN_HOME or CMDSTAN\n',...
        'environment variable. See the dependencies support matrix for details.\n'));

    if isempty(depvercheck)

        [~, savepath] = uiputfile('PsyRATDependents',...
        'Where would you like to save the directory for the dependents?');

        %if the user does not select a file, then take the user back to psyrat_start
        if savepath == 0
            %relaunch first, dialog last (B20 creation order); modal per B21
            psyrat_start;
            errordlg('Location not selected','File Error','modal');
            return;
        end

        wrkdir = fullfile(savepath,'PsyRATDependents');

    elseif ~isempty(depvercheck)

        wrkdir = fileparts(fileparts(which('runCmdStanTests.py')));

        %G56 (partial hardening). Same guard as the macOS arm (see the
        %comment there): an empty which() means no in-place directory to
        %update; stop before rmolddeps, by the B20 exit convention.
        if isempty(wrkdir)
            reshome = psyrat_locate_cmdstan();
            if isempty(reshome)
                whereline = 'No CmdStan install could be resolved.';
            else
                whereline = sprintf('CmdStan was found at: %s', reshome);
            end
            dlgtext = {...
                'The update flow could not find the dependents directory.'; ...
                whereline; ...
                'That install is not on the MATLAB path, so there is no directory to update in place.'; ...
                'Add that directory to the MATLAB path (addpath + savepath) and rerun the update,'; ...
                'run psyrat_installdependents from the Command Window to choose a fresh install directory,'; ...
                'or update an environment-variable (STAN_HOME/CMDSTAN) install manually.'};
            %relaunch first, dialog last (B20 creation order); modal per B21
            psyrat_start;
            errordlg(dlgtext,'Update location not found','modal');
            return;
        end

        %G56c: an env-var pin outranks any fresh download at fit time;
        %disclose the inert update and let the user cancel BEFORE anything
        %is removed.
        if ~local_envpin_update_gate(depvercheck)
            return;
        end

        rmolddeps(depvercheck);

    end

    %the third argument carries the version verdict ([] outside the
    %update flow; see the macOS call above; G56)
    psyrat_linuxdepsinstall(wrkdir,urls,depvercheck);
end

end

function psyrat_macdepsinstall(wrkdir,urls,depvercheck)

%get the current directory so it can be reverted to after installation
startdir = cd;

%create a structure array that will store which dependents are not properly
%installed
%0 - not installed, 1 - installed
depcheck = struct;

%first check for the XCode command line tools
if exist('/Library/Developer/CommandLineTools','file') == 7
    depcheck.CLT = 1;
    fprintf('XCode Command line tools is installed\n');
else
    depcheck.CLT = 0;
    fprintf('XCode Command line tools does not appear to be installed\n');
end

%check for cmdstan.
%G51: the exist() tests below see only the MATLAB path, but a CmdStan the
%startup gate accepts can live entirely off-path (an environment-variable or
%home-glob install; see psyrat_locate_cmdstan/G43). Consult the shared
%resolver FIRST so such a user is not told to install a CmdStan that
%estimation can already use; the path-only checks remain as the fallback for
%installs the resolver cannot see. The resolver validates makefile+bin/, so
%built-ness is still checked separately, against the resolved home.
cmdstanhome = psyrat_locate_cmdstan();
if ~isempty(cmdstanhome)

    %CmdStan has been installed (found by the shared resolver)
    depcheck.cmdstan = 1;
    fprintf('CmdStan is installed\n');

    %Record the resolved home as THE CmdStan directory for this run. The
    %downstream consumers (the make-build repair, the MatlabStan
    %stan_home pointing) fall back to which('runCmdStanTests.py') only
    %when cmdstandir is unset, and that which() probe is empty for the
    %off-path installs this arm exists for -- without this line the gate
    %could declare the resolved home unbuilt and then cd('') or build a
    %DIFFERENT path-visible install (adversarial wave, 2026-08-29).
    cmdstandir = cmdstanhome;

    %check whether the resolved CmdStan has been properly built
    if exist(fullfile(cmdstanhome,'bin','stansummary'),'file') ~= 2 || ...
        exist(fullfile(cmdstanhome,'bin','stanc'),'file') ~= 2

        %CmdStan has not been built
        depcheck.cmdbuild = 0;
        fprintf('CmdStan has not been properly built\n');
    else
        %CmdStan has been built
        depcheck.cmdbuild = 1;
        fprintf('CmdStan has been properly built\n');
    end

elseif exist('makefile','file') ~= 2 || ...
    exist('test-all.sh','file') ~= 2 || ...
    exist('runCmdStanTests.py','file') ~= 2

    %CmdStan has not been installed
    depcheck.cmdstan = 0;

    %CmdStan has not been built
    depcheck.cmdbuild = 0;
    fprintf('CmdStan needs to be installed and built\n');
else

    %CmdStan has been installed
    depcheck.cmdstan = 1;
    fprintf('CmdStan is installed\n');

    %check whether CmdStan has been properly built
    if exist('stansummary','file') ~= 2 || ...
        exist('stanc','file') ~= 2

        %CmdStan has not been built
        depcheck.cmdbuild = 0;
        fprintf('CmdStan has not been properly built\n');
    else
        %CmdStan has been built
        depcheck.cmdbuild = 1;
        fprintf('CmdStan has been properly built\n');
    end
end

%G56 (partial hardening). The old-version verdict used to be discarded
%before this point: the gate above re-derives installed-ness from scratch
%AFTER rmolddeps, so a CmdStan the verdict flagged OLD that survived
%removal reads as installed and the download block below is skipped --
%previously with no mention, an all-Installed status GUI (plain text;
%psyrat_guistatus colors nothing), and the version check
%flagging the same install again at next launch. Say so plainly. The
%realistic survivor classes are a deletion psyrat_safe_rmdepdir refused
%(reported on the console above) and a SECOND install the gate resolves
%after the first was removed; the off-path class the first draft named
%cannot reach here any more, because the dispatcher's empty-wrkdir guard
%returns before the per-OS installers run. The verdict conditioning
%matters: rmolddeps targets CmdStan only when depvercheck.cmdstan == 0,
%so a CURRENT CmdStan survives legitimately (an update prompted by the
%other dependents) and gets no notice. Force-reinstall semantics and what
%"update" should mean for env-var-pinned installs are the open ledger
%item G56.
if ~isempty(depvercheck) && depvercheck.cmdstan == 0 && depcheck.cmdstan == 1
    if ~exist('cmdstandir','var') || isempty(cmdstandir)
        %single fileparts: runCmdStanTests.py sits at the CmdStan root, so
        %its folder IS the home. The double-fileparts spelling the wrkdir
        %derivations use names the PARENT dependents directory, which must
        %never appear in a message that invites manual deletion (it holds
        %MatlabStan and MatlabProcessManager too).
        survivor = fileparts(which('runCmdStanTests.py'));
    else
        survivor = cmdstandir;
    end
    fprintf(strcat('NOTE: an update was requested, but an existing CmdStan install\n',...
        '(%s)\n',...
        'survived the removal step, so no new CmdStan will be downloaded.\n',...
        'If a removal refusal was printed above, it explains why the\n',...
        'install was kept. At fit time the toolbox resolves CmdStan\n',...
        'through its normal tiers (environment variables, then managed\n',...
        'home directories, then the MATLAB path), so estimation may use\n',...
        'this or another install. To force a fresh install, remove that\n',...
        'CmdStan directory manually and rerun the update, or install\n',...
        'CmdStan yourself and point STAN_HOME or CMDSTAN at it.\n'), survivor);
end

%check for MatlabProcessManager (MPM)
if exist('processManager.m','file') ~= 2 || ...
    exist('processState.m','file') ~= 2

    %MPM not installed
    depcheck.mpm = 0;
    fprintf('MatlabProcessManager is not installed\n');
else
    %MPM is installed
    depcheck.mpm = 1;
    fprintf('MatlabProcessManager is installed\n');
end

%check for MatlabStan
if exist('StanFit.m','file') ~= 2 || ...
    exist('StanModel.m','file') ~= 2

    %mstan not installed
    depcheck.mstan = 0;
    fprintf('MatlabStan is not installed\n');
else
    %mstan is installed
    depcheck.mstan = 1;
    fprintf('MatlabStan is installed\n');
end

%pull up a gui to display the status of the PsyRAT toolbox dependents
psyrat_guistatus(depcheck);

%check whether the directory to save the files exists, if not create it
if exist(wrkdir,'file') ~=7
    %make the directory where the files will be saved
    mkdir(wrkdir);
end

%if needed, install CLT
if depcheck.CLT == 0
    %user will need to interacte with gui for installation of CLT
    system('xcode-select --install');
    
    fsize = get(0,'DefaultTextFontSize') + 3;
    figwidth = 400; figheight = 200;
    %initialize gui
    f = psyrat_newfigure('unit','pix','Visible','off',...
      'position',[400 400 figwidth figheight],...
      'menub','no',...
      'numbertitle','off',...
      'resize','off');

    movegui(f,'center');

    %Write text
    uicontrol(f,'Style','text','fontsize',fsize+2,...
        'HorizontalAlignment','center',...
        'String','Please indicate when the XCode installation is complete',...
        'Position',[0 100 figwidth figheight/3]);          

    %Create a button that will resume installations
    uicontrol(f,'Style','push','fontsize',fsize,...
        'HorizontalAlignment','center',...
        'String',{'Installation';'Complete'},...
        'Position', [(figwidth/5)/2 25 (figwidth/5)*4 75],...
        'Callback','uiresume(gcbf)'); 
    
    %display gui
    set(f,'Visible','on');
    %tag gui
    %
    %B49: this was the one-character Tag 'f'. findobj with no handle argument
    %searches the whole graphics root and is not restricted to figures, so any
    %HANDLE-VISIBLE object tagged 'f' -- including a uicontrol inside somebody
    %else's figure -- matched the lookup below.
    %
    %AND THE CONSEQUENCE WAS WORSE THAN A STRAY WINDOW CLOSING. close() rejects
    %non-figure handles, so a matched uicontrol was NOT quietly closed: the call
    %threw. Measured in R2025b -- findobj returns the uicontrol, then close()
    %raises MATLAB:close:InvalidFigureHandle and the uicontrol survives. With a
    %mixed match close() aborts before closing anything, so even this figure
    %stayed up. There is no try/catch between psyrat_macdepsinstall's call site
    %and here, so the exception escaped psyrat_installdependents entirely,
    %killing the Command Line Tools install the instant the user clicked
    %"Installation Complete".
    %
    %Namespaced like the rest of the toolbox's guis, which also brings it
    %under psyrat_closepsyratfigs' 'psyrat_' prefix (B48).
    f.Tag = 'psyrat_gui_cltwait';
    
    uiwait(f);
    
    %check if the Command Line Tools wait prompt is still open.
    f = findobj('Tag','psyrat_gui_cltwait');
    if ~isempty(f)
        close(f);
    end
    
    %Re-check the install location rather than trusting the launch-time exit
    %status of `xcode-select --install`: that command returns nonzero when the
    %tools are already present and returns 0 merely when the installer dialog
    %launches (not when installation actually completes). Confirming the
    %directory after the user clicks "Installation Complete" is reliable.
    if exist('/Library/Developer/CommandLineTools','file') == 7
        depcheck.CLT = 1;
        fprintf('XCode Command Line Tools successfully installed\n');
        %pull up a gui to display the status
        psyrat_guistatus(depcheck);
    else
        error('CLT:notinstalled',... %Error code and associated error
        strcat('WARNING: Automatic installation of XCode CLT failed \n\n',...
        'Please install the command line tools manually\n'));
    end
end

%check whether cmdstan needs to be installed
if depcheck.cmdstan == 0

    %url for the cmdstan release tarball (bundles the Stan submodules)
    loc = urls.cmdstan_tar;

    %change the working directory to where the files will be installed
    cd(wrkdir);

    %download the cmdstan tarball
    fileout = websave('cmdstan.tar.gz',loc);

    %unpack the tarball using the system tar. Matlab's untar was found not to
    %unpack the cmdstan tarball correctly on macOS, so system tar stays the
    %primary extractor; untar is only a fallback if system tar is unavailable.
    %Plain (POSIX) tar flags are used so the command works with both BSD tar
    %(the default on macOS) and GNU tar; the GNU-only --no-same-owner flag is
    %intentionally avoided.
    status = system('tar -x -z -f cmdstan.tar.gz');
    if status ~= 0
        try
            untar(fileout,wrkdir);
        catch
            error('CmdStan:notinstalled',... %Error code and associated error
                strcat('ERROR: The system failed to unpack the cmdstan tarball\n\n',...
                'Please verify that a C++ development environment is installed\n',...
                'The recommended C++ development environment is XCode\n'));
        end
    end

    %delete the tarball after it's been unpacked
    delete(fullfile(wrkdir,'cmdstan.tar.gz'));

    %now find the path to the new cmdstan directory
    ls = dir(wrkdir);
    ind = find(strncmp({ls.name}, 'cmdstan', 7)==1);

    if isempty(ind)
        error('CmdStan:notinstalled',... %Error code and associated error
            strcat('ERROR: No cmdstan folder was found after unpacking the tarball\n\n',...
            'Please verify the download succeeded and rerun psyrat_installdependents'));
    end

    %if more than one cmdstan folder exists, keep the most recent one
    if length(ind) > 1
        [~,newind] = max([ls(ind).datenum]);
        ind = ind(newind);
    end

    cmdstandir = fullfile(wrkdir,ls(ind).name);

    %update depcheck and display status
    depcheck.cmdstan = 1;
    fprintf('CmdStan successfully downloaded and unpacked\n');
    psyrat_guistatus(depcheck);
end

%check whether CmdStan needs to be built
if depcheck.cmdbuild == 0
    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir
    if ~exist('cmdstandir','var')
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    if isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    
    %change the cd to where cmdstan is installed so it can be built
    cd(cmdstandir);
    
    %find out how many cores are available for building stan
    ncores = feature('numCores');
    if ncores > 2
        ncores = ncores - 1;
    end

    %execute command to build stan
    buildcmdstan = strcat('make build -j',num2str(ncores));
    status = system(buildcmdstan);

    %update depcheck and display status
    if status == 0
        depcheck.cmdbuild = 1;
        fprintf('CmdStan successfully built\n');
        psyrat_guistatus(depcheck);
    else
        error('CmdStan:notbuilt',... %Error code and associated error
            strcat('ERROR: Automatic build of CmdStan failed \n\n',...
            'Please install CmdStan manually\n'));
    end
    
    newcmdstanbuild = 1;
    
end

%check whether the Matlab Process Manager needs to be installed
if depcheck.mpm == 0
    %use the bundled MatlabProcessManager if available, otherwise download it
    bundledmpm = psyrat_copybundledmpm(wrkdir);

    if bundledmpm == 0
        %the bundled copy could not be found, so fall back to the upstream
        %download. Say so: the download is the unpatched upstream and the
        %patch-level check at the next launch would otherwise be the first
        %notice (audit 2026-09-05).
        warning('psyrat:unbundledDependency',...
            ['The bundled MatlabProcessManager was not found under bundled_dependents, ',...
            'so it is being downloaded from %s instead.'],urls.mpm_tar);
        %url for Matlab Process Manager
        loc = urls.mpm_tar;

        %download tarball
        fileout = websave('mpm.tar.gz',loc);

        %unpack tarball
        untar(fileout,wrkdir);
    else
        fprintf('Using bundled MatlabProcessManager\n');
    end
    
    %change status of depcheck and display gui
    depcheck.mpm = 1;
    fprintf('Matlab Process manager successfully installed\n');
    psyrat_guistatus(depcheck);
end

%check whether MatlabStan is installed
if depcheck.mstan == 0
    %use bundled MatlabStan if available, otherwise download it
    [bundledms,msdir] = psyrat_copybundledmatlabstan(wrkdir);
    
    if bundledms == 0
        %the bundled, PATCHED copy could not be found, so fall back to the
        %upstream download. The upstream copy lacks the patches PsyRAT relies
        %on (see bundled_dependents/MatlabStan-*/PSYRAT_PATCHLEVEL.txt), so say
        %so now rather than at the next launch's patch-level check (audit
        %2026-09-05).
        warning('psyrat:unbundledDependency',...
            ['The bundled, patched MatlabStan was not found under bundled_dependents, ',...
            'so the UNPATCHED upstream copy is being downloaded from %s. Restore the ',...
            'bundled_dependents folder and rerun the installer to get the patched copy.'],...
            urls.ms_tar);
        %url for MatlabStan
        loc = urls.ms_tar;

        %download tarball
        fileout = websave('ms.tar.gz',loc);

        %unpack tarball
        untar(fileout,wrkdir);
        
        %get dir of MatlabStan so stan_home can be edited
        msdir = psyrat_findmatlabstandir(wrkdir);
    else
        fprintf('Using bundled MatlabStan\n');
    end
    
    %delete the odd file that can be created some times
    if exist(fullfile(wrkdir,'pax_global_header'),'file') ~= 0
        delete(fullfile(wrkdir,'pax_global_header'));
    end
    
    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir. Guard ~exist first: on a reinstall where CmdStan is already
    %installed AND built (depcheck.cmdstan/cmdbuild both 1) but MatlabStan is
    %missing, neither assignment path above ran, so cmdstandir is undefined
    %here. The Windows/Linux installers use this same combined guard.
    if ~exist('cmdstandir','var') || isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    cd(cmdstandir);

    stanhomefile = fullfile(msdir,'+mstan','stan_home.m');
    psyrat_set_stanhome_default(stanhomefile,cmdstandir);
    
    
    %open the StanModel.m file and change the version to current cmdstan
    % version
    fid=fopen(fullfile(msdir,'StanModel.m'));
    storefile = {};
    while 1
        tline = fgetl(fid);
        storefile{end+1} = tline;
        if ~ischar(tline), break, end
    end
    fclose(fid);

    ind = strcmp(storefile,"            ver = cellfun(@str2num,regexp(str{3},'\.','split'));");
    
    %Only rewrite the version line when the legacy MatlabStan template that
    %contains it is present. The bundled PsyRAT fork detects the CmdStan
    %version at runtime (StanModel.stan_version) and has no such line, so a
    %blind assignment here errors on an all-false logical index (MATLAB:
    %index:expected_one_output_for_assignment).
    if any(ind)
        jvers = psyrat_dependentsversions;
        storefile{ind} = ['            ver = ''' jvers.cmdstan ''';'];
    end

    %overwrite the existing stan_home.m file with the new information
    fid=fopen(fullfile(msdir,'StanModel.m'),'w');
    for i=1:(length(storefile)-1)
        storefile{i} = strrep(storefile{i},'%','%%');
        storefile{i} = strrep(storefile{i},'\','\\');
        fprintf(fid,[storefile{i} '\n']);
    end
    fclose(fid);
    
    %update depcheck and display gui
    depcheck.mstan = 1;
    fprintf('MatlabStan successfully installed\n');
    psyrat_guistatus(depcheck);
  
%this is necessary in case cmdstan was updated, but a new version  
%matlabstan was not installed. Matlabstan needs to be pointed to cmdstan.
elseif depcheck.mstan == 1 && exist('newcmdstanbuild','var')
    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir
    if isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    cd(cmdstandir);
    
    %get dir of MatlabStan so stan_home can be edited
    msdir = psyrat_findmatlabstandir(wrkdir);
    if isempty(msdir)
        p = which('StanModel.m');
        if ~isempty(p)
            msdir = fileparts(p);
        end
    end
    
    stanhomefile = fullfile(msdir,'+mstan','stan_home.m');
    if isempty(msdir) || exist(stanhomefile,'file') ~= 2
        error('MatlabStan:notfound',... %Error code and associated error
            strcat('ERROR: MatlabStan was detected but stan_home.m could not be found\n\n',...
            'Please verify MatlabStan is installed and included in the Matlab path\n',...
            'then rerun psyrat_installdependents'));
    end

    psyrat_set_stanhome_default(stanhomefile,cmdstandir);
    
    %update depcheck and display gui
    depcheck.mstan = 1;
    fprintf('MatlabStan successfully installed\n');
    psyrat_guistatus(depcheck);
end

%add files to the Matlab path and save it
addpath(genpath(wrkdir));
psyrat_savepath_safe;
fprintf('Directories for dependents added to Matlab path\n');

%close gui
psyrat_instgui = findobj('Tag','psyrat_instgui');

if ~isempty(psyrat_instgui)
    
    close(psyrat_instgui);
    
end

%go back to the starting directory
cd(startdir);

%rerun PsyRAT Toolbox
psyrat_start;

end


function psyrat_windepsinstall(wrkdir,urls,depvercheck)
%install dependents on Windows Operating System

%get the current directory so it can be reverted to after installation
startdir = cd;

%create a structure array that will store which dependents are not properly
%installed
%0 - not installed, 1 - installed
depcheck = struct;

%first check whether the make and g++ compilers exist (CmdStan needs both)
[~,makeout] = system('make -v');
makematch = strncmp(makeout,'GNU Make',8);

%g++ ships with the Rtools toolchain; `g++ --version` exits 0 and prints a
%banner mentioning g++/GCC when it is on the PATH. A missing compiler instead
%returns a nonzero status and a shell "not recognized" message.
[gppstatus,gppout] = system('g++ --version');
gppmatch = (gppstatus == 0) && ~isempty(regexpi(gppout,'g\+\+|GCC','once'));

if makematch == 1 && gppmatch
    depcheck.CLT = 1;
    fprintf('Command line tools installed\n');
else
    error('Rtools:notinstalled',... %Error code and associated error
        strcat('WARNING: Command line tools not installed\n\n',...
        'Automatic installation of command line tools not supported for Windows\n',...
        'CmdStan requires both GNU make and the g++ compiler on the PATH\n',...
        'Rtools contains the necessary command line tools needed for CmdStan\n',...
        'Sorry I was unable to automate this process for Windows users!\n\n'));
end

%check for cmdstan.
%G51: consult the shared resolver first (see the macOS gate above); the
%path-only checks remain as the fallback. Built-ness on Windows keys on the
%.exe binaries under the resolved home's bin directory.
cmdstanhome = psyrat_locate_cmdstan();
if ~isempty(cmdstanhome)

    %CmdStan has been installed (found by the shared resolver)
    depcheck.cmdstan = 1;
    fprintf('CmdStan is installed\n');

    %Record the resolved home as THE CmdStan directory for this run. The
    %downstream consumers (the make-build repair, the MatlabStan
    %stan_home pointing) fall back to which('runCmdStanTests.py') only
    %when cmdstandir is unset, and that which() probe is empty for the
    %off-path installs this arm exists for -- without this line the gate
    %could declare the resolved home unbuilt and then cd('') or build a
    %DIFFERENT path-visible install (adversarial wave, 2026-08-29).
    cmdstandir = cmdstanhome;

    %check whether the resolved CmdStan has been properly built
    if exist(fullfile(cmdstanhome,'bin','stansummary.exe'),'file') ~= 2 || ...
        exist(fullfile(cmdstanhome,'bin','stanc.exe'),'file') ~= 2

        %CmdStan has not been built
        depcheck.cmdbuild = 0;
        fprintf('CmdStan has not been properly built\n');
    else
        %CmdStan has been built
        depcheck.cmdbuild = 1;
        fprintf('CmdStan has been properly built\n');
    end

elseif exist('makefile','file') ~= 2 || ...
    exist('test-all.sh','file') ~= 2 || ...
    exist('runCmdStanTests.py','file') ~= 2

    %CmdStan has not been installed
    depcheck.cmdstan = 0;

    %CmdStan has not been built
    depcheck.cmdbuild = 0;
    fprintf('CmdStan needs to be installed and built\n');
else

    %CmdStan has been installed
    depcheck.cmdstan = 1;
    fprintf('CmdStan is installed\n');

    %check whether CmdStan has been properly built
    if exist('stansummary.exe','file') ~= 2 || ...
        exist('stanc.exe','file') ~= 2

        %CmdStan has not been built
        depcheck.cmdbuild = 0;
        fprintf('CmdStan has not been properly built\n');
    else
        %CmdStan has been built
        depcheck.cmdbuild = 1;
        fprintf('CmdStan has been properly built\n');
    end
end

%G56 (partial hardening). The old-version verdict used to be discarded
%before this point: the gate above re-derives installed-ness from scratch
%AFTER rmolddeps, so a CmdStan the verdict flagged OLD that survived
%removal reads as installed and the download block below is skipped --
%previously with no mention, an all-Installed status GUI (plain text;
%psyrat_guistatus colors nothing), and the version check
%flagging the same install again at next launch. Say so plainly. The
%realistic survivor classes are a deletion psyrat_safe_rmdepdir refused
%(reported on the console above) and a SECOND install the gate resolves
%after the first was removed; the off-path class the first draft named
%cannot reach here any more, because the dispatcher's empty-wrkdir guard
%returns before the per-OS installers run. The verdict conditioning
%matters: rmolddeps targets CmdStan only when depvercheck.cmdstan == 0,
%so a CURRENT CmdStan survives legitimately (an update prompted by the
%other dependents) and gets no notice. Force-reinstall semantics and what
%"update" should mean for env-var-pinned installs are the open ledger
%item G56.
if ~isempty(depvercheck) && depvercheck.cmdstan == 0 && depcheck.cmdstan == 1
    if ~exist('cmdstandir','var') || isempty(cmdstandir)
        %single fileparts: runCmdStanTests.py sits at the CmdStan root, so
        %its folder IS the home. The double-fileparts spelling the wrkdir
        %derivations use names the PARENT dependents directory, which must
        %never appear in a message that invites manual deletion (it holds
        %MatlabStan and MatlabProcessManager too).
        survivor = fileparts(which('runCmdStanTests.py'));
    else
        survivor = cmdstandir;
    end
    fprintf(strcat('NOTE: an update was requested, but an existing CmdStan install\n',...
        '(%s)\n',...
        'survived the removal step, so no new CmdStan will be downloaded.\n',...
        'If a removal refusal was printed above, it explains why the\n',...
        'install was kept. At fit time the toolbox resolves CmdStan\n',...
        'through its normal tiers (environment variables, then managed\n',...
        'home directories, then the MATLAB path), so estimation may use\n',...
        'this or another install. To force a fresh install, remove that\n',...
        'CmdStan directory manually and rerun the update, or install\n',...
        'CmdStan yourself and point STAN_HOME or CMDSTAN at it.\n'), survivor);
end

%check for MatlabProcessManager (MPM)
if exist('processManager.m','file') ~= 2 || ...
    exist('processState.m','file') ~= 2 
    
    %MPM not installed
    depcheck.mpm = 0;
    fprintf('MatlabProcessManager is not installed\n');
else
    %MPM is installed
    depcheck.mpm = 1;
    fprintf('MatlabProcessManager is installed\n');
end

%check for MatlabStan
if exist('StanFit.m','file') ~= 2 || ...
    exist('StanModel.m','file') ~= 2 
    
    %mstan not installed
    depcheck.mstan = 0;
    fprintf('MatlabStan is not installed\n');
else
    %mstan is installed
    depcheck.mstan = 1;
    fprintf('MatlabStan is installed\n');
end

%pull up a gui to display the status of the PsyRAT toolbox dependents
psyrat_guistatus(depcheck);

%check whether the directory to save the files exists, if not create it
if exist(wrkdir,'file') ~=7
    %make the directory where the files will be saved
    mkdir(wrkdir);
end

%check whether cmdstan needs to be installed
if depcheck.cmdstan == 0 

    %use release tarball (includes Stan submodules)
    loc = urls.cmdstan_tar;

    %change the working directory to where the files will be installed
    cd(wrkdir);

    %download the cmdstan tarball
    fileout = websave('cmdstan.tar.gz',loc);

    %unpack the tarball
    %prefer Matlab untar but fall back to system tar in case untar fails
    try
        untar(fileout,wrkdir);
    catch
        status = system('tar -x -z -f cmdstan.tar.gz');
        if status ~= 0
            error('CmdStan:notinstalled',... %Error code and associated error
                strcat('ERROR: The system failed to unpack the cmdstan tarball\n\n',...
                'Please verify that C++ development tools are installed\n'));
        end
    end

    %delete the tarball after it's been unpacked
    delete(fullfile(wrkdir,'cmdstan.tar.gz'));

    %now find the path to the new cmdstan directory
    ls = dir(wrkdir);
    ind = find(strncmp({ls.name}, 'cmdstan', 7)==1);

    if isempty(ind)
        error('CmdStan:notinstalled',... %Error code and associated error
            strcat('ERROR: No cmdstan folder was found after unpacking the tarball\n\n',...
            'Please verify the download succeeded and rerun psyrat_installdependents'));
    end

    %if more than one cmdstan folder exists, keep the most recent one
    if length(ind) > 1
        [~,newind] = max([ls(ind).datenum]);
        ind = ind(newind);
    end
    
    cmdstandir = fullfile(wrkdir,ls(ind).name);
    
    %update depcheck and display status
    depcheck.cmdstan = 1;
    fprintf('CmdStan successfully installed\n');
    psyrat_guistatus(depcheck);
    
end

%check whether CmdStan needs to be built
if depcheck.cmdbuild == 0
    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir
    if ~exist('cmdstandir','var')
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    if isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    
    %change the cd to where cmdstan is installed so it can be built
    cd(cmdstandir);
    
    %find out how many cores are available for building stan
    ncores = feature('numCores');
    if ncores > 2
        ncores = ncores - 1;
    end

    %execute command to build stan
    buildcmdstan = strcat('make build -j',num2str(ncores));
    status = system(buildcmdstan);

    %update depcheck and display status
    if status == 0
        depcheck.cmdbuild = 1;
        fprintf('CmdStan successfully built\n');
        psyrat_guistatus(depcheck);
    else
        error('CmdStan:notbuilt',... %Error code and associated error
            strcat('WARNING: Automatic build of CmdStan failed \n\n',...
            'Please install CmdStan manually\n'));
    end
    
    cmdstannewbuild = 1;
end

%check whether the Matlab Process Manager needs to be installed
if depcheck.mpm == 0
    %use the bundled MatlabProcessManager if available, otherwise download it
    bundledmpm = psyrat_copybundledmpm(wrkdir);

    if bundledmpm == 0
        %url for Matlab Process Manager
        loc = urls.mpm_zip;

        %download zip
        fileout = websave('mpm.zip',loc);

        %unpack zip
        unzip(fileout,wrkdir);
    else
        fprintf('Using bundled MatlabProcessManager\n');
    end
    
    %change status of depcheck and display gui
    depcheck.mpm = 1;
    fprintf('Matlab Process manager successfully installed\n');
    psyrat_guistatus(depcheck);
end

%check whether MatlabStan is installed
if depcheck.mstan == 0
    %use bundled MatlabStan if available, otherwise download it
    [bundledms,msdir] = psyrat_copybundledmatlabstan(wrkdir);
    
    if bundledms == 0
        %url for MatlabStan (zip archive on Windows)
        loc = urls.ms_zip;

        %download zip
        fileout = websave('ms.zip',loc);

        %unpack zip
        unzip(fileout,wrkdir);
        
        %get dir of MatlabStan so stan_home can be edited
        msdir = psyrat_findmatlabstandir(wrkdir);
    else
        fprintf('Using bundled MatlabStan\n');
    end
    
    %delete the odd file that can be created some times
    if exist(fullfile(wrkdir,'pax_global_header'),'file') ~= 0
        delete(fullfile(wrkdir,'pax_global_header'));
    end
    
    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir
    if ~exist('cmdstandir','var') || isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    elseif exist('cmdstandir','var') && isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    cd(cmdstandir);

    stanhomefile = fullfile(msdir,'+mstan','stan_home.m');
    psyrat_set_stanhome_default(stanhomefile,cmdstandir);
    
    %open the StanModel.m file and change the version to current cmdstan
    % version
    fid=fopen(fullfile(msdir,'StanModel.m'));
    storefile = {};
    while 1
        tline = fgetl(fid);
        storefile{end+1} = tline;
        if ~ischar(tline), break, end
    end
    fclose(fid);

    ind = strcmp(storefile,"            ver = cellfun(@str2num,regexp(str{3},'\.','split'));");
    
    %Only rewrite the version line when the legacy MatlabStan template that
    %contains it is present. The bundled PsyRAT fork detects the CmdStan
    %version at runtime (StanModel.stan_version) and has no such line, so a
    %blind assignment here errors on an all-false logical index (MATLAB:
    %index:expected_one_output_for_assignment).
    if any(ind)
        jvers = psyrat_dependentsversions;
        storefile{ind} = ['            ver = ''' jvers.cmdstan ''';'];
    end

    %overwrite the existing stan_home.m file with the new information
    fid=fopen(fullfile(msdir,'StanModel.m'),'w');
    for i=1:(length(storefile)-1)
        storefile{i} = strrep(storefile{i},'%','%%');
        storefile{i} = strrep(storefile{i},'\','\\');
        fprintf(fid,[storefile{i} '\n']);
    end
    fclose(fid);
    
    %update depcheck and display gui
    depcheck.mstan = 1;
    fprintf('MatlabStan successfully installed\n');
    psyrat_guistatus(depcheck);
    
%this is necessary if cmdstan is updated. Matlabstan needs to have the 
%location of cmdstan updated.
elseif depcheck.mstan == 1 && exist('cmdstannewbuild','var')
    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir
    if isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    cd(cmdstandir);
    
    %get dir of MatlabStan so stan_home can be edited
    msdir = psyrat_findmatlabstandir(wrkdir);
    if isempty(msdir)
        p = which('StanModel.m');
        if ~isempty(p)
            msdir = fileparts(p);
        end
    end
    
    stanhomefile = fullfile(msdir,'+mstan','stan_home.m');
    if isempty(msdir) || exist(stanhomefile,'file') ~= 2
        error('MatlabStan:notfound',... %Error code and associated error
            strcat('ERROR: MatlabStan was detected but stan_home.m could not be found\n\n',...
            'Please verify MatlabStan is installed and included in the Matlab path\n',...
            'then rerun psyrat_installdependents'));
    end

    psyrat_set_stanhome_default(stanhomefile,cmdstandir);
    
    %update depcheck and display gui
    depcheck.mstan = 1;
    fprintf('MatlabStan successfully installed\n');
    psyrat_guistatus(depcheck);
end

%add files to the Matlab path and save it
addpath(genpath(wrkdir));
psyrat_savepath_safe;
fprintf('Directories for dependents added to Matlab path\n');

%close gui
psyrat_instgui = findobj('Tag','psyrat_instgui');

if ~isempty(psyrat_instgui)
    
    close(psyrat_instgui);
    
end

%go back to the starting directory
cd(startdir);

%rerun PsyRAT Toolbox
psyrat_start;

end

function psyrat_linuxdepsinstall(wrkdir,urls,depvercheck)
%install dependents on Linux (EXPERIMENTAL - unverified by the maintainer)
%
%This branch mirrors the macOS path: it checks for a C++ toolchain, downloads
%the CmdStan release tarball, builds it with make, and installs
%MatlabProcessManager and the bundled MatlabStan. It has not been verified on
%a Linux machine. If anything fails, install CmdStan manually and point the
%toolbox at it with the STAN_HOME or CMDSTAN environment variable.

%get the current directory so it can be reverted to after installation
startdir = cd;

%create a structure array that will store which dependents are not properly
%installed
%0 - not installed, 1 - installed
depcheck = struct;

%first check whether the make and g++ compilers exist (CmdStan needs both)
[~,makeout] = system('make -v');
makematch = strncmp(makeout,'GNU Make',8);

[gppstatus,gppout] = system('g++ --version');
gppmatch = (gppstatus == 0) && ~isempty(regexpi(gppout,'g\+\+|GCC','once'));

if makematch == 1 && gppmatch
    depcheck.CLT = 1;
    fprintf('Command line tools (make and g++) found\n');
else
    error('Buildtools:notinstalled',... %Error code and associated error
        strcat('WARNING: A C++ build toolchain was not found\n\n',...
        'CmdStan requires both GNU make and the g++ compiler on the PATH\n',...
        'Install them with your package manager, for example:\n',...
        '  Debian/Ubuntu: sudo apt-get install build-essential\n',...
        '  Fedora/RHEL:   sudo dnf groupinstall "Development Tools"\n',...
        'then rerun psyrat_installdependents\n\n'));
end

%check for cmdstan.
%G51: consult the shared resolver first (see the macOS gate above); the
%path-only checks remain as the fallback.
cmdstanhome = psyrat_locate_cmdstan();
if ~isempty(cmdstanhome)

    %CmdStan has been installed (found by the shared resolver)
    depcheck.cmdstan = 1;
    fprintf('CmdStan is installed\n');

    %Record the resolved home as THE CmdStan directory for this run. The
    %downstream consumers (the make-build repair, the MatlabStan
    %stan_home pointing) fall back to which('runCmdStanTests.py') only
    %when cmdstandir is unset, and that which() probe is empty for the
    %off-path installs this arm exists for -- without this line the gate
    %could declare the resolved home unbuilt and then cd('') or build a
    %DIFFERENT path-visible install (adversarial wave, 2026-08-29).
    cmdstandir = cmdstanhome;

    %check whether the resolved CmdStan has been properly built
    if exist(fullfile(cmdstanhome,'bin','stansummary'),'file') ~= 2 || ...
        exist(fullfile(cmdstanhome,'bin','stanc'),'file') ~= 2

        %CmdStan has not been built
        depcheck.cmdbuild = 0;
        fprintf('CmdStan has not been properly built\n');
    else
        %CmdStan has been built
        depcheck.cmdbuild = 1;
        fprintf('CmdStan has been properly built\n');
    end

elseif exist('makefile','file') ~= 2 || ...
    exist('test-all.sh','file') ~= 2 || ...
    exist('runCmdStanTests.py','file') ~= 2

    %CmdStan has not been installed
    depcheck.cmdstan = 0;

    %CmdStan has not been built
    depcheck.cmdbuild = 0;
    fprintf('CmdStan needs to be installed and built\n');
else

    %CmdStan has been installed
    depcheck.cmdstan = 1;
    fprintf('CmdStan is installed\n');

    %check whether CmdStan has been properly built
    if exist('stansummary','file') ~= 2 || ...
        exist('stanc','file') ~= 2

        %CmdStan has not been built
        depcheck.cmdbuild = 0;
        fprintf('CmdStan has not been properly built\n');
    else
        %CmdStan has been built
        depcheck.cmdbuild = 1;
        fprintf('CmdStan has been properly built\n');
    end
end

%G56 (partial hardening). The old-version verdict used to be discarded
%before this point: the gate above re-derives installed-ness from scratch
%AFTER rmolddeps, so a CmdStan the verdict flagged OLD that survived
%removal reads as installed and the download block below is skipped --
%previously with no mention, an all-Installed status GUI (plain text;
%psyrat_guistatus colors nothing), and the version check
%flagging the same install again at next launch. Say so plainly. The
%realistic survivor classes are a deletion psyrat_safe_rmdepdir refused
%(reported on the console above) and a SECOND install the gate resolves
%after the first was removed; the off-path class the first draft named
%cannot reach here any more, because the dispatcher's empty-wrkdir guard
%returns before the per-OS installers run. The verdict conditioning
%matters: rmolddeps targets CmdStan only when depvercheck.cmdstan == 0,
%so a CURRENT CmdStan survives legitimately (an update prompted by the
%other dependents) and gets no notice. Force-reinstall semantics and what
%"update" should mean for env-var-pinned installs are the open ledger
%item G56.
if ~isempty(depvercheck) && depvercheck.cmdstan == 0 && depcheck.cmdstan == 1
    if ~exist('cmdstandir','var') || isempty(cmdstandir)
        %single fileparts: runCmdStanTests.py sits at the CmdStan root, so
        %its folder IS the home. The double-fileparts spelling the wrkdir
        %derivations use names the PARENT dependents directory, which must
        %never appear in a message that invites manual deletion (it holds
        %MatlabStan and MatlabProcessManager too).
        survivor = fileparts(which('runCmdStanTests.py'));
    else
        survivor = cmdstandir;
    end
    fprintf(strcat('NOTE: an update was requested, but an existing CmdStan install\n',...
        '(%s)\n',...
        'survived the removal step, so no new CmdStan will be downloaded.\n',...
        'If a removal refusal was printed above, it explains why the\n',...
        'install was kept. At fit time the toolbox resolves CmdStan\n',...
        'through its normal tiers (environment variables, then managed\n',...
        'home directories, then the MATLAB path), so estimation may use\n',...
        'this or another install. To force a fresh install, remove that\n',...
        'CmdStan directory manually and rerun the update, or install\n',...
        'CmdStan yourself and point STAN_HOME or CMDSTAN at it.\n'), survivor);
end

%check for MatlabProcessManager (MPM)
if exist('processManager.m','file') ~= 2 || ...
    exist('processState.m','file') ~= 2

    %MPM not installed
    depcheck.mpm = 0;
    fprintf('MatlabProcessManager is not installed\n');
else
    %MPM is installed
    depcheck.mpm = 1;
    fprintf('MatlabProcessManager is installed\n');
end

%check for MatlabStan
if exist('StanFit.m','file') ~= 2 || ...
    exist('StanModel.m','file') ~= 2

    %mstan not installed
    depcheck.mstan = 0;
    fprintf('MatlabStan is not installed\n');
else
    %mstan is installed
    depcheck.mstan = 1;
    fprintf('MatlabStan is installed\n');
end

%pull up a gui to display the status of the PsyRAT toolbox dependents
psyrat_guistatus(depcheck);

%check whether the directory to save the files exists, if not create it
if exist(wrkdir,'dir') ~= 7
    %make the directory where the files will be saved
    mkdir(wrkdir);
end

%check whether cmdstan needs to be installed
if depcheck.cmdstan == 0

    %url for the cmdstan release tarball (bundles the Stan submodules)
    loc = urls.cmdstan_tar;

    %change the working directory to where the files will be installed
    cd(wrkdir);

    %download the cmdstan tarball
    fileout = websave('cmdstan.tar.gz',loc);

    %unpack the tarball using the system tar (primary), falling back to
    %Matlab's untar only if system tar is unavailable. POSIX tar flags keep
    %this working with GNU tar.
    status = system('tar -x -z -f cmdstan.tar.gz');
    if status ~= 0
        try
            untar(fileout,wrkdir);
        catch
            error('CmdStan:notinstalled',... %Error code and associated error
                strcat('ERROR: The system failed to unpack the cmdstan tarball\n\n',...
                'Please verify that a C++ development environment is installed\n',...
                'then rerun psyrat_installdependents'));
        end
    end

    %delete the tarball after it's been unpacked
    delete(fullfile(wrkdir,'cmdstan.tar.gz'));

    %now find the path to the new cmdstan directory
    ls = dir(wrkdir);
    ind = find(strncmp({ls.name}, 'cmdstan', 7)==1);

    if isempty(ind)
        error('CmdStan:notinstalled',... %Error code and associated error
            strcat('ERROR: No cmdstan folder was found after unpacking the tarball\n\n',...
            'Please verify the download succeeded and rerun psyrat_installdependents'));
    end

    %if more than one cmdstan folder exists, keep the most recent one
    if length(ind) > 1
        [~,newind] = max([ls(ind).datenum]);
        ind = ind(newind);
    end

    cmdstandir = fullfile(wrkdir,ls(ind).name);

    %update depcheck and display status
    depcheck.cmdstan = 1;
    fprintf('CmdStan successfully downloaded and unpacked\n');
    psyrat_guistatus(depcheck);
end

%check whether CmdStan needs to be built
if depcheck.cmdbuild == 0
    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir
    if ~exist('cmdstandir','var')
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    if isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end

    %change the cd to where cmdstan is installed so it can be built
    cd(cmdstandir);

    %find out how many cores are available for building stan
    ncores = feature('numCores');
    if ncores > 2
        ncores = ncores - 1;
    end

    %execute command to build stan
    buildcmdstan = strcat('make build -j',num2str(ncores));
    status = system(buildcmdstan);

    %update depcheck and display status
    if status == 0
        depcheck.cmdbuild = 1;
        fprintf('CmdStan successfully built\n');
        psyrat_guistatus(depcheck);
    else
        error('CmdStan:notbuilt',... %Error code and associated error
            strcat('ERROR: Automatic build of CmdStan failed \n\n',...
            'Please install CmdStan manually\n',...
            'then rerun psyrat_installdependents'));
    end

    newcmdstanbuild = 1;

end

%check whether the Matlab Process Manager needs to be installed
if depcheck.mpm == 0
    %use the bundled MatlabProcessManager if available, otherwise download it
    bundledmpm = psyrat_copybundledmpm(wrkdir);

    if bundledmpm == 0
        %the bundled copy could not be found, so fall back to the upstream
        %download. Say so: the download is the unpatched upstream and the
        %patch-level check at the next launch would otherwise be the first
        %notice (audit 2026-09-05).
        warning('psyrat:unbundledDependency',...
            ['The bundled MatlabProcessManager was not found under bundled_dependents, ',...
            'so it is being downloaded from %s instead.'],urls.mpm_tar);
        %url for Matlab Process Manager
        loc = urls.mpm_tar;

        %download tarball
        fileout = websave('mpm.tar.gz',loc);

        %unpack tarball
        untar(fileout,wrkdir);
    else
        fprintf('Using bundled MatlabProcessManager\n');
    end

    %change status of depcheck and display gui
    depcheck.mpm = 1;
    fprintf('Matlab Process manager successfully installed\n');
    psyrat_guistatus(depcheck);
end

%check whether MatlabStan is installed
if depcheck.mstan == 0
    %use bundled MatlabStan if available, otherwise download it
    [bundledms,msdir] = psyrat_copybundledmatlabstan(wrkdir);

    if bundledms == 0
        %the bundled, PATCHED copy could not be found, so fall back to the
        %upstream download. The upstream copy lacks the patches PsyRAT relies
        %on (see bundled_dependents/MatlabStan-*/PSYRAT_PATCHLEVEL.txt), so say
        %so now rather than at the next launch's patch-level check (audit
        %2026-09-05).
        warning('psyrat:unbundledDependency',...
            ['The bundled, patched MatlabStan was not found under bundled_dependents, ',...
            'so the UNPATCHED upstream copy is being downloaded from %s. Restore the ',...
            'bundled_dependents folder and rerun the installer to get the patched copy.'],...
            urls.ms_tar);
        %url for MatlabStan
        loc = urls.ms_tar;

        %download tarball
        fileout = websave('ms.tar.gz',loc);

        %unpack tarball
        untar(fileout,wrkdir);

        %get dir of MatlabStan so stan_home can be edited
        msdir = psyrat_findmatlabstandir(wrkdir);
    else
        fprintf('Using bundled MatlabStan\n');
    end

    %delete the odd file that can be created some times
    if exist(fullfile(wrkdir,'pax_global_header'),'file') ~= 0
        delete(fullfile(wrkdir,'pax_global_header'));
    end

    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir
    if ~exist('cmdstandir','var') || isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    cd(cmdstandir);

    stanhomefile = fullfile(msdir,'+mstan','stan_home.m');
    psyrat_set_stanhome_default(stanhomefile,cmdstandir);

    %open the StanModel.m file and change the version to current cmdstan
    % version
    fid=fopen(fullfile(msdir,'StanModel.m'));
    storefile = {};
    while 1
        tline = fgetl(fid);
        storefile{end+1} = tline;
        if ~ischar(tline), break, end
    end
    fclose(fid);

    ind = strcmp(storefile,"            ver = cellfun(@str2num,regexp(str{3},'\.','split'));");

    %Only rewrite the version line when the legacy MatlabStan template that
    %contains it is present. The bundled PsyRAT fork detects the CmdStan
    %version at runtime (StanModel.stan_version) and has no such line, so a
    %blind assignment here errors on an all-false logical index.
    if any(ind)
        jvers = psyrat_dependentsversions;
        storefile{ind} = ['            ver = ''' jvers.cmdstan ''';'];
    end

    %overwrite the existing StanModel.m file with the new information
    fid=fopen(fullfile(msdir,'StanModel.m'),'w');
    for i=1:(length(storefile)-1)
        storefile{i} = strrep(storefile{i},'%','%%');
        storefile{i} = strrep(storefile{i},'\','\\');
        fprintf(fid,[storefile{i} '\n']);
    end
    fclose(fid);

    %update depcheck and display gui
    depcheck.mstan = 1;
    fprintf('MatlabStan successfully installed\n');
    psyrat_guistatus(depcheck);

%this is necessary in case cmdstan was updated, but a new version
%matlabstan was not installed. Matlabstan needs to be pointed to cmdstan.
elseif depcheck.mstan == 1 && exist('newcmdstanbuild','var')
    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir
    if isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    cd(cmdstandir);

    %get dir of MatlabStan so stan_home can be edited
    msdir = psyrat_findmatlabstandir(wrkdir);
    if isempty(msdir)
        p = which('StanModel.m');
        if ~isempty(p)
            msdir = fileparts(p);
        end
    end

    stanhomefile = fullfile(msdir,'+mstan','stan_home.m');
    if isempty(msdir) || exist(stanhomefile,'file') ~= 2
        error('MatlabStan:notfound',... %Error code and associated error
            strcat('ERROR: MatlabStan was detected but stan_home.m could not be found\n\n',...
            'Please verify MatlabStan is installed and included in the Matlab path\n',...
            'then rerun psyrat_installdependents'));
    end

    psyrat_set_stanhome_default(stanhomefile,cmdstandir);

    %update depcheck and display gui
    depcheck.mstan = 1;
    fprintf('MatlabStan successfully installed\n');
    psyrat_guistatus(depcheck);
end

%add files to the Matlab path and save it
addpath(genpath(wrkdir));
psyrat_savepath_safe;
fprintf('Directories for dependents added to Matlab path\n');

%close gui
psyrat_instgui = findobj('Tag','psyrat_instgui');

if ~isempty(psyrat_instgui)

    close(psyrat_instgui);

end

%go back to the starting directory
cd(startdir);

%rerun PsyRAT Toolbox
psyrat_start;

end

function psyrat_set_stanhome_default(stanhomefile,cmdstandir)

if exist(stanhomefile,'file') ~= 2
    error('MatlabStan:notfound',... %Error code and associated error
        strcat('ERROR: stan_home.m could not be found\n\n',...
        'Path: ',stanhomefile));
end

fid = fopen(stanhomefile);
if fid == -1
    error('MatlabStan:fileopen',... %Error code and associated error
        strcat('ERROR: Unable to open stan_home.m for reading\n\n',...
        'Path: ',stanhomefile));
end

storefile = {};
while 1
    tline = fgetl(fid);
    storefile{end+1} = tline;
    if ~ischar(tline), break, end
end
fclose(fid);

escaped = strrep(cmdstandir,'\','\\');
newline = ['d = ''' escaped ''';'];
replace_idx = [];
for i = 1:(length(storefile)-1)
    if ~isempty(regexp(storefile{i},'^\s*d\s*=\s*''[^'']*''\s*;\s*$','once'))
        replace_idx = i;
        break;
    end
end

if isempty(replace_idx)
    %fallback for unexpected templates: insert after function header
    replace_idx = min(4,length(storefile)-1);
end
storefile{replace_idx} = newline;

fid = fopen(stanhomefile,'w');
if fid == -1
    error('MatlabStan:fileopen',... %Error code and associated error
        strcat('ERROR: Unable to open stan_home.m for writing\n\n',...
        'Path: ',stanhomefile));
end
for i=1:(length(storefile)-1)
    storefile{i} = strrep(storefile{i},'%','%%');
    storefile{i} = strrep(storefile{i},'\','\\');
    fprintf(fid,[storefile{i} '\n']);
end
fclose(fid);

end

function [copied,msdir] = psyrat_copybundledmatlabstan(wrkdir)

copied = 0;
msdir = [];

depvers = psyrat_dependentsversions;
toolboxroot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
bundledsrc = fullfile(toolboxroot,'bundled_dependents',...
    strcat('MatlabStan-',depvers.matlabstan));

if exist(bundledsrc,'dir') == 7
    [~,foldername] = fileparts(bundledsrc);
    msdir = fullfile(wrkdir,foldername);
    
    if exist(msdir,'dir') ~= 7
        copyfile(bundledsrc,msdir);
    end
    
    copied = 1;
end

end

function copied = psyrat_copybundledmpm(wrkdir)
%Use the vendored MatlabProcessManager under bundled_dependents/ instead of
%fetching it from GitHub. Mirrors psyrat_copybundledmatlabstan; the process
%manager needs no post-copy editing.

copied = 0;

depvers = psyrat_dependentsversions;
toolboxroot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
bundledsrc = fullfile(toolboxroot,'bundled_dependents',...
    strcat('MatlabProcessManager-',depvers.matlabprocessmanager));

if exist(bundledsrc,'dir') == 7
    [~,foldername] = fileparts(bundledsrc);
    mpmdir = fullfile(wrkdir,foldername);

    if exist(mpmdir,'dir') ~= 7
        copyfile(bundledsrc,mpmdir);
    end

    copied = 1;
end

end

function msdir = psyrat_findmatlabstandir(wrkdir)

msdir = [];
ls = dir(wrkdir);
ind = find(strncmp({ls.name}, 'MatlabStan', 10)==1);

if ~isempty(ind)
    if length(ind)>1
        [~,newind] = max([ls(ind).datenum]);
        ind = ind(newind);
    end
    
    msdir = fullfile(wrkdir,ls(ind).name);
end

end


function psyrat_guistatus(depcheck)

%define parameters for figure position
figwidth = 600;
figheight = 450;
fsize = get(0,'DefaultTextFontSize') + 3;

%define space between rows and first row location
rowspace = 35;
row = figheight - rowspace*2;

%define locations of column 1 and 2
lcol = 30;
rcol = (figwidth/8)*5;

psyrat_instgui = findobj('Tag','psyrat_instgui');

if isempty(psyrat_instgui)

    %create the gui
    psyrat_instgui= psyrat_newfigure('unit','pix','Visible','off',...
      'position',[400 400 figwidth figheight],...
      'menub','no',...
      'name','Installation Summary',...
      'numbertitle','off',...
      'resize','off');

    movegui(psyrat_instgui,'center');
    
end

%Print the name of the gui
uicontrol(psyrat_instgui,'Style','text','fontsize',fsize+4,...
    'HorizontalAlignment','center',...
    'String','Installation Status',...
    'Position',[0 row figwidth 25]);     

row = row - rowspace*2;

%CLI status
uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String','Command Line Tools',...
    'Position', [lcol row figwidth/3 25]);  

if depcheck.CLT == 0
    uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
        'String','Not Installed',...
        'Position', [rcol row figwidth/4 25]);  
elseif depcheck.CLT == 1
    uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
        'String','Installed',...
        'Position', [rcol row figwidth/4 25]);  
end

row = row - rowspace;

%CmdStan status
uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String','CmdStan Files',...
    'Position', [lcol row figwidth/3 25]);  

if depcheck.cmdstan == 0
    uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
        'String','Not Installed',...
        'Position', [rcol row figwidth/4 25]);  
elseif depcheck.cmdstan == 1
    uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
        'String','Installed',...
        'Position', [rcol row figwidth/4 25]);  
end

row = row - rowspace;

%CmdStan build status
uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String','CmdStan Build',...
    'Position', [lcol row figwidth/3 25]);  

if depcheck.cmdbuild == 0
    uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
        'String','Not Built',...
        'Position', [rcol row figwidth/4 25]);  
elseif depcheck.cmdbuild == 1
    uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
        'String','Built',...
        'Position', [rcol row figwidth/4 25]);  
end

row = row - rowspace;

%Matlab Process Manager status
uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String','Matlab Process Manager',...
    'Position', [lcol row figwidth/3 25]);  

if depcheck.mpm == 0
    uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
        'String','Not Installed',...
        'Position', [rcol row figwidth/4 25]);  
elseif depcheck.mpm == 1
    uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
        'String','Installed',...
        'Position', [rcol row figwidth/4 25]);  
end

row = row - rowspace;

%mstan status
uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String','Matlab Stan',...
    'Position', [lcol row figwidth/3 25]);  

if depcheck.mstan == 0
    uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
        'String','Not Installed',...
        'Position', [rcol row figwidth/4 25]);  
elseif depcheck.mstan == 1
    uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
        'String','Installed',...
        'Position', [rcol row figwidth/4 25]);  
end

row = row - rowspace*2;

%mstan status
uicontrol(psyrat_instgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','Installation may take a while...',...
    'Position',[0 row figwidth 25]);

%make sure the gui is displayed if it has not been shown. Tag the gui as
%well
if strcmp(get(psyrat_instgui,'Visible'),'off')
    set(psyrat_instgui,'Visible','on');
    
    %tag gui
    psyrat_instgui.Tag = 'psyrat_instgui';
end

%pause for the gui to be displayed
pause(.02);

end

function rmolddeps(depvercheck)

fprintf('Removing old versions of dependents\nThis may take a while\n\n');

%temporarily turn off warnings that filepaths are being removed
warning('off','MATLAB:rmpath:DirNotFound');
warning('off','MATLAB:RMDIR:RemovedFromPath');

%Each target is a directory resolved from the Matlab path, which is not
%necessarily an installed copy: PsyRAT puts its own vendored bundled_dependents
%tree on the path and the installed copies carry identical folder names, so a
%bare which() can resolve to PsyRAT's own source. psyrat_safe_rmdepdir refuses
%any target inside the toolbox tree and reports rather than throwing, so one
%refusal does not abandon the other dependents.
if depvercheck.cmdstan == 0
    %delete old cmdstan files
    local_report_rmdep(fileparts(which('runCmdStanTests.py')),'CmdStan');
end

if depvercheck.matlabstan == 0
    %delete old matlabstan files
    local_report_rmdep(fileparts(which('mcmc.m')),'MatlabStan');
end

if depvercheck.matlabprocessmanager == 0
    %delete old matlabprocessmanager files
    local_report_rmdep(fileparts(which('processManager.m')), ...
        'MatlabProcessManager');
end

%turn warnings back on
warning('on','MATLAB:rmpath:DirNotFound');
warning('on','MATLAB:RMDIR:RemovedFromPath');

end

function local_report_rmdep(target,label)
%LOCAL_REPORT_RMDEP Remove one dependency directory and say what happened.
%   A refusal is printed rather than raised: it is a normal outcome when the
%   resolved copy is PsyRAT's own bundled one, which needs no reinstalling.

[removed,reason] = psyrat_safe_rmdepdir(target,label);
if removed
    fprintf('Removed old %s directory: %s\n',label,target);
else
    fprintf('%s\n',reason);
end

end

function proceed = local_envpin_update_gate(depvercheck)
%LOCAL_ENVPIN_UPDATE_GATE G56c: disclose an inert CmdStan update up front.
%   At fit time the CmdStan resolver's environment tier outranks every other
%   tier (psyrat_locate_cmdstan: env vars first), so when the resolver's
%   winning home IS an env-var pin that this update cannot replace, the
%   freshly downloaded CmdStan will never be selected while that pin stands.
%   Say so BEFORE anything is removed and let the user proceed or cancel.
%
%   The gate fires only when ALL of the following hold, each mirroring the
%   condition of the mechanism it reports (the G56 surviving-install lesson):
%   1. depvercheck.cmdstan == 0 -- the same verdict rmolddeps reads before
%      touching CmdStan. On a MatlabStan/MPM-only update nothing CmdStan is
%      removed or downloaded, so a CmdStan-pin dialog would be noise.
%   2. A set env var, normalized the way the resolver normalizes candidates
%      (local_envpin_normalize mirrors +mstan/stan_home's
%      local_normalize_cmdstan_home: strtrim; a FILE resolves to its parent,
%      stanc/stanc.exe to its grandparent; a directory named bin with a
%      sibling makefile to its parent), names the resolver's winning home.
%      Comparing the RAW env value instead missed genuine pins spelled as
%      STAN_HOME=<home>/makefile or <home>/bin (adversarial wave).
%   3. The pinned home is NOT the directory rmolddeps is about to delete
%      (fileparts(which('runCmdStanTests.py')), rmolddeps' own target). When
%      it is, the pin dies with the deletion, the env candidate goes invalid,
%      the resolver falls through to the fresh download, and the update IS
%      effective -- a dialog claiming otherwise would be false. This is also
%      the state the dev sweep manufactures by seeding
%      PSYRAT_TEST_CMDSTAN_PATH with the resolver's own answer.
%
%   A cancel is a chosen outcome, not an error: print why, then relaunch
%   psyrat_start so the user is not left windowless (the update callback
%   already closed every PsyRAT window before this runs). The relaunch
%   cannot precede the dialog here -- questdlg must block for the answer
%   that decides whether to relaunch -- so this is deliberately NOT the
%   B20 relaunch-first creation order the file's error exits use.

proceed = true;

%condition 1: the verdict rmolddeps reads (fields are set by
%psyrat_checkversionsofdeps; the update arms always pass a populated struct)
if ~isstruct(depvercheck) || ~isfield(depvercheck,'cmdstan') || ...
        depvercheck.cmdstan ~= 0
    return;
end

home = psyrat_locate_cmdstan();
if isempty(home)
    return; %no resolvable install, so nothing can be pinned
end

%condition 2: the resolver's env tier, in its own precedence order
pinvars = {'PSYRAT_TEST_CMDSTAN_PATH','STAN_HOME','CMDSTAN_HOME','CMDSTAN'};
pinned = '';
for ii = 1:numel(pinvars)
    val = local_envpin_normalize(getenv(pinvars{ii}));
    if ~isempty(val) && local_samedir(val,home)
        pinned = pinvars{ii};
        break;
    end
end

if isempty(pinned)
    return; %the winner is not an env pin, so the update can take effect
end

%condition 3: a pin naming rmolddeps' own removal target dies with the
%removal, so the update is effective and needs no dialog
rmtarget = fileparts(which('runCmdStanTests.py'));
if ~isempty(rmtarget) && local_samedir(rmtarget,home)
    return;
end

fprintf(['\nCmdStan is pinned by the environment variable %s (%s).\n' ...
    'A downloaded CmdStan will not be selected while that pin stands;\n' ...
    'any other dependents this update covers are updated normally.\n\n'], ...
    pinned,home);

choice = questdlg({ ...
    sprintf('CmdStan is pinned by the environment variable %s:',pinned); ...
    home; ...
    ''; ...
    'PsyRAT selects that install ahead of any newly downloaded CmdStan,'; ...
    'so this update will not change which CmdStan estimation uses'; ...
    'while the pin stands.'; ...
    ''; ...
    'Any other dependents this update covers (MatlabStan,'; ...
    'MatlabProcessManager) are updated normally either way.'}, ...
    'CmdStan is pinned by an environment variable', ...
    'Update anyway','Cancel','Cancel');

if ~strcmp(choice,'Update anyway')
    proceed = false;
    fprintf('\nDependents will not be updated (environment-variable pin).\n\n');
    psyrat_start;
end

end

function p = local_envpin_normalize(p)
%LOCAL_ENVPIN_NORMALIZE Env value -> the home the resolver would derive.
%   Mirrors +mstan/stan_home's local_normalize_cmdstan_home, which is the
%   normalization every resolver candidate -- env values included -- passes
%   through before the validity test, so the gate's comparison sees the same
%   home the resolver's winner was derived from. Kept in the same order:
%   strtrim; a FILE resolves to its parent, except stanc/stanc.exe to its
%   grandparent (they live in bin/); a directory named bin with a sibling
%   makefile resolves to its parent; anything nonexistent to ''.

if isstring(p)
    p = char(p);
end
if isempty(p) || ~ischar(p)
    p = '';
    return;
end
p = strtrim(p);
if isempty(p)
    return;
end

if exist(p,'file') == 2
    [pp,nm,ex] = fileparts(p);
    token = lower([nm ex]);
    if strcmp(token,'stanc') || strcmp(token,'stanc.exe')
        p = fileparts(pp);
    else
        p = pp;
    end
elseif exist(p,'dir') == 7
    [pp,nm] = fileparts(p);
    if strcmpi(nm,'bin') && exist(fullfile(pp,'makefile'),'file') == 2
        p = pp;
    end
else
    p = '';
end

end

function tf = local_samedir(a,b)
%LOCAL_SAMEDIR True when two paths name the same directory.
%   Canonicalizes through the filesystem (dir() reports the resolved parent
%   for a directory, the same trick psyrat_safe_rmdepdir uses) and compares
%   with a case sensitivity matching the platform's filesystem, so a path
%   spelled with a trailing separator still matches the resolver's winner.
%   Symlinked spellings are NOT promised to match (whether dir() reports the
%   physical path is unmeasured here); a missed match only skips the
%   disclosure, it removes nothing.

a = local_canondir(a);
b = local_canondir(b);
if isempty(a) || isempty(b)
    tf = false;
elseif ispc || ismac
    tf = strcmpi(a,b);
else
    tf = strcmp(a,b);
end

end

function p = local_canondir(p)
%LOCAL_CANONDIR Absolute directory path without a trailing separator.

if isempty(p)
    p = '';
    return;
end
p = char(p);

d = dir(p);
if ~isempty(d) && isfield(d,'folder') && ~isempty(d(1).folder)
    if strcmp(d(1).name,'.')
        p = d(1).folder;
    end
end

while numel(p) > 1 && endsWith(p,filesep)
    p = p(1:end-1);
end

end
