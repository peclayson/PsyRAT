function [psyrat_prefs,psyrat_data] = psyrat_findprefsdata(varin)
%Somersault through inputs to find psyrat_prefs and psyrat_data
%
%[psyrat_prefs,psyrat_data] = psyrat_findprefsdata(varargin)
%
%
%Input
% varin - varargin from guis
%
%output
% psyrat_prefs - PsyRAT Toolbox structure array containing preferences
% psyrat_data - PsyRAT Toolbox structure array containing data
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

ind = find(strcmp('psyrat_prefs',varin),1);
if ~isempty(ind)
    psyrat_prefs = varin{ind+1}; 
else
    psyrat_prefs = '';
end

%check if psyrat_data has been defined
ind = find(strcmp('psyrat_data',varin),1);
if ~isempty(ind)
    psyrat_data = varin{ind+1}; 
else
    psyrat_data = '';
end

end

