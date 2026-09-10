function out = psyrat_adapter_epmat_detect(filepath)
%Detect EP Toolkit exported MATLAB struct files (.mat/.ept).

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

if ~any(strcmp(ext,{'.mat','.ept'}))
    out = struct('confidence',0,'reason','File extension is not .mat or .ept.');
    return;
end

confidence = 0.15;
reason = 'Generic MATLAB/EP Toolkit file; EP structure not yet confirmed.';

try
    vars = whos('-file',filepath);
    names = lower(string({vars.name}));
    if any(names == "epdata")
        confidence = 0.98;
        reason = 'File contains variable EPdata.';
    elseif any(contains(names,'ep'))
        confidence = 0.45;
        reason = 'MAT file variable names suggest EP Toolkit content.';
    end
catch
    confidence = 0.05;
    reason = 'Could not inspect MAT file variables.';
end

out = struct('confidence',confidence,'reason',reason);

end
