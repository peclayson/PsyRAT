function psyrat_updatepsyrat
%Update the PsyRAT Toolbox from the PsyRAT.zip asset of the latest stable
% GitHub release (resolved through the releases API). This currently only
% works for some newer versions of Mac or Windows. If using Linux, the PsyRAT
% will need to be updated manually.
%
%psyrat_updatepsyrat
%
%
%Required Inputs:
% No inputs are required.
% Some features only work for certain versions of Mac or Windows.
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

%Resolve the latest STABLE release through the GitHub releases API, the same
%path psyrat_start's startup check uses, and download the PsyRAT.zip asset
%that release/PUBLISHING.md attaches to every release (audit finding G87,
%owner ruling 2026-09-09). Until that ruling this function fetched the
%releases web page, parsed a heading whose markup GitHub has since changed
%and whose title the runbook never wrote, and appended a guessed asset name
%to a download URL no release carried; because the update prompt's Yes
%callback closes every PsyRAT window before calling here, that failure left
%the user with no window at all.
%A check that THROWS (it wraps only its own web request) must not escape
%here: the update prompt's Yes callback has already closed every PsyRAT
%window, so route it into the same refusal below, with the releases index.
try
    info = psyrat_checkgithubreleases('verbose', false);
catch me
    %only the fields the refusal arm reads (statusCode -1 never reaches the
    %tag or the version)
    info = struct('statusCode', -1, 'latestReleaseAssetURL', '', ...
        'latestReleaseURL', 'https://github.com/peclayson/PsyRAT/releases', ...
        'message', ['The release check failed: ' me.message]);
end

%the GitHub page every message below names: the release's own page when a
%stable release was found, the releases index otherwise
urlstr = info.latestReleaseURL;

if info.statusCode ~= 0 || isempty(info.latestReleaseAssetURL)
    %Nothing downloadable: no newer stable release (statusCode 1, 2 or 3),
    %the check itself failed (-1), or the release carries no PsyRAT.zip
    %asset. Relaunch FIRST and create the modal dialog LAST (B20 creation
    %order; modal per B21), exactly as the cancel arms below do.
    if info.statusCode == 0
        why = sprintf(['Release %s does not provide a PsyRAT.zip download, ', ...
            'so the toolbox cannot be updated automatically. Download the ', ...
            'release from the page below and install it by hand.'], ...
            info.latestReleaseTag);
    else
        why = info.message;
        if isempty(why)
            why = 'No downloadable PsyRAT release was found.';
        end
    end
    psyrat_start;
    errordlg({why; ['Download page: ' urlstr]}, 'Update unavailable', 'modal');
    return;
end

%the asset's own download URL and the semantic version of that release
asseturl = info.latestReleaseAssetURL;
ver = info.latestReleaseVersion;

%version number will be appended to directory
psyrat_dirname = strcat('PsyRAT_Toolbox_v',ver);

str_cantdl = 'The new version of the toolbox is available on ';
str_cantdl = [str_cantdl... 
    '<a href="matlab:web(''https://github.com/peclayson/PsyRAT/releases/latest'',''-browser'')">Github</a>'];

%determine the version of OS that is being used
if ismac
    sys = 1; %because Apple is number 1, obviously
elseif ispc
    sys = 2; %because Windows OS is inferior
elseif isunix && ~ismac
    sys = 3; %Linux (EXPERIMENTAL: auto-update is not verified by the maintainer)
end

%check whether the version of the OS is supported
if sys == 1 %mac
    [~,cmdout] = system('sw_vers');
    parseout = strsplit(cmdout);
    parseOS = strsplit(parseout{6},'.');

    %sw_vers token positions shifted when the macOS ProductName changed
    %from "Mac OS X" to "macOS"; if the build-version token was grabbed
    %instead of the product version, fall back to the ProductVersion token.
    %(Mirrors the detection in psyrat_installdependents so the toolbox
    %self-update also works on macOS 11 and newer.)
    if isnan(str2double(parseOS{1})) || ~contains(char(parseOS{1}),'.')
        parseOS = strsplit(parseout{4},'.');
    end

    if str2double(parseOS{1}) >= 10
        if str2double(parseOS{2}) >= 9 || str2double(parseOS{1}) >= 11
            [~, savepath] = uiputfile('PsyRAT',...
            'Where would you like to save the PsyRAT?');
        
            %if the user does not select a file, then take the user back to psyrat_start    
            if savepath == 0
                %relaunch first, dialog last (B20 creation order); modal per B21
                psyrat_start;
                errordlg('Location not selected','File Error','modal');
                return;
            end
            
            wrkdir = fullfile(savepath,psyrat_dirname);
            psyrat_toolboxinstall(wrkdir,asseturl);
            
        else
            error('mac:oldver',... %Error code and associated error
            strcat('WARNING: Automatic installation only works for OSX 10.9 and newer \n\n',...
            str_cantdl));
        end
    else
        error('mac:oldver',... %Error code and associated error
        strcat('WARNING: Automatic installation only works for Mac OS 10.9 or newer \n\n',...
            str_cantdl));
    end
elseif sys == 2 %windows
    fprintf('Warning: automatic installation has only been tested using Windows 7\n');
    [~, savepath] = uiputfile('PsyRAT',...
    'Where would you like to save the directory for the updated PsyRAT?');

    %if the user does not select a file, then take the user back to psyrat_start    
    if savepath == 0
        %relaunch first, dialog last (B20 creation order); modal per B21
        psyrat_start;
        errordlg('Location not selected','File Error','modal');
        return;
    end

    wrkdir = fullfile(savepath,psyrat_dirname);

    psyrat_toolboxinstall(wrkdir,asseturl);
elseif sys == 3 %linux (EXPERIMENTAL - unverified by the maintainer)
    warning('os:linuxExperimental',... %not an error: Linux update is experimental
        strcat('EXPERIMENTAL: automatic PsyRAT toolbox update on Linux has not\n',...
        'been verified by the maintainer. The toolbox itself is pure MATLAB, so\n',...
        'the download/unzip update is expected to work, but proceed with caution.\n'));
    [~, savepath] = uiputfile('PsyRAT',...
    'Where would you like to save the directory for the updated PsyRAT?');

    %if the user does not select a file, then take the user back to psyrat_start
    if savepath == 0
        %relaunch first, dialog last (B20 creation order); modal per B21
        psyrat_start;
        errordlg('Location not selected','File Error','modal');
        return;
    end

    wrkdir = fullfile(savepath,psyrat_dirname);

    psyrat_toolboxinstall(wrkdir,asseturl);
end
        
end

function psyrat_toolboxinstall(wrkdir,asseturl)
%install the new toolbox from the PsyRAT.zip asset at ASSETURL, the direct
%download URL the releases API reported; the archive holds one top-level
%PsyRAT folder, which is what the movefile below expects

%get the current directory so it can be changed after installation
startdir = cd;

%get the directory of the current psyrat_start file
old_psyratdir = which('psyrat_start.m');
old_psyratdir = fileparts(old_psyratdir);

%download zip
fileout = websave('psyrat_toolbox.zip',asseturl);

%get the installation directory
installdir = fileparts(fileout);

%unzip
unzip(fileout,installdir);

%move the files to where they're supposed to go
movefile(fullfile(installdir,'PsyRAT'),wrkdir)

%delete the zip file after it's been unpacked
delete(fileout);

%temporarily turn off warnings that filepaths are being removed
warning('off','MATLAB:rmpath:DirNotFound');

%remove old PsyRAT files from the path
rmpath(genpath(old_psyratdir)); 

%turn warnings back on
warning('on','MATLAB:rmpath:DirNotFound');

%add files to the Matlab path and save it
addpath(genpath(wrkdir));
psyrat_savepath_safe;
fprintf('Directories for PsyRAT Toolbox added and saved to Matlab path\n');

%go back to the starting directory
if strcmp(genpath(old_psyratdir),startdir)
    cd(startdir);
else
    cd(wrkdir);
end

psyrat_start;

end
