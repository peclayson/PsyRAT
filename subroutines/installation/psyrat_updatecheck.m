function psyratinstalledver = psyrat_updatecheck(psyratver)
%Check whether a new release has been posted on Github
%
%psyrat_updatecheck
%
%Lasted Updated 8/5/19
%
%Required Input:
% psyratver - PsyRAT Toolbox version
%
%Outputs:
% psyratinstalledver - contains information regarding whether the toolbox is
%  up to date 0 - old, 1 - current, 2 - beta
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

try
    %pull webpage from github
    urlstr = 'https://www.peterclayson.com/wp-content/uploads/psyrat_currentversion.txt';
    webraw = webread(urlstr,'text','html');

    %extract semantic version even if extra text/html is returned
    ver = local_extractversion(webraw);
    
    %format version of the currently installed toolbox
    used_parts = sscanf(psyratver,'%d.%d.%d')';
    if numel(used_parts) ~= 3
        error('psyrat_updatecheck:versionparse', ...
            'Unable to parse installed version: %s', psyratver);
    end

    %format version of the toolbox found online
    found_parts = sscanf(ver,'%d.%d.%d')';
    if numel(found_parts) ~= 3
        error('psyrat_updatecheck:versionparse', ...
            'Unable to parse available version: %s', ver);
    end

    %compare cmdstan versions
    for ii = 1:3
        if used_parts(ii) < found_parts(ii)
            psyratinstalledver = 0;
            break;
        elseif used_parts(ii) > found_parts(ii)
            psyratinstalledver = 2;
            break;
        elseif used_parts(ii) == found_parts(ii)
            psyratinstalledver = 1;
        end
    end
    
    %let user know what was found
    switch psyratinstalledver
        case 1
            str = 'You are running the most up-to-date version of the toolbox';
            fprintf('\n%s\n',str);
        case 0
            str = 'There is a new version of the toolbox available on ';
            str = [str... 
                '<a href="matlab:web(''https://github.com/peclayson/PsyRAT/releases/latest'',''-browser'')">Github</a>'];
            fprintf('\n%s\n',str);
        case 2
            warning('You are running a non-stable release of the toolbox');
            fprintf(strcat(...
                'You likely cloned github, rather than installed the latest stable release\n',...
                'As a result, complete functionality cannot be guaranteed\n',...
                'The lastest stable release can be downloaded at'));
            str = '<a href="matlab:web(''https://github.com/peclayson/PsyRAT/releases/latest'',''-browser'')">Github</a>';
            fprintf('\n\n%s\n',str);
    end
catch
    fprintf('\nUnable to connect to internet to check for new releases\n');
    psyratinstalledver = [];
end

end

function ver = local_extractversion(webraw)
token = regexp(char(webraw), '\d+\.\d+\.\d+', 'match', 'once');
if isempty(token)
    error('psyrat_updatecheck:versiontoken', ...
        'No semantic version token found in update check response.');
end
ver = strtrim(token);
end
