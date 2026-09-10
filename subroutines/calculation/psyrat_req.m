function value = psyrat_req(args,name,msg)
%Fetch a required name-value argument from a varargin cell array.
%
% value = psyrat_req(args,name,msg)
%
% Looks up 'name' (case-insensitive) in the name-value cell array 'args' and
% returns the paired value. If 'name' is absent, raises an error with id
% 'varargin:<lowername>' and the supplied message.
%
% This is the shared implementation of the required-argument getter that was
% previously duplicated as an identical local function in psyrat_cutscore_sing
% and psyrat_rel_sing. See psyrat_opt for the optional-argument counterpart.
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

ind = find(strcmpi(name,args),1);
if ~isempty(ind)
    value = args{ind+1};
else
    error(['varargin:' lower(name)],... %Error code and associated error
        strcat('WARNING: ',msg,'\n'));
end
end
