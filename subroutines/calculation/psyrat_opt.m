function value = psyrat_opt(args,name,default)
%Fetch an optional name-value argument from a varargin cell array.
%
% value = psyrat_opt(args,name,default)
%
% Looks up 'name' (case-insensitive) in the name-value cell array 'args' and
% returns the paired value, or 'default' when 'name' is absent. This is the
% shared form of the find(strcmpi(...))-with-default idiom used across the
% calculation functions. See psyrat_req for the required-argument counterpart.
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
    value = default;
end
end
