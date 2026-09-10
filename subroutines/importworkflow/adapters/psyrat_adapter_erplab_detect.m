function out = psyrat_adapter_erplab_detect(filepath)
%Detect ERPLAB ERP files.

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

[~,~,ext] = fileparts(filepath);
ext = lower(ext);

if strcmp(ext,'.erp')
    out = struct('confidence',0.95,'reason','File extension .erp matches ERPLAB ERP file.');
else
    out = struct('confidence',0,'reason','File extension is not .erp.');
end

end
