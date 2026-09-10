function det = psyrat_import_detect_format(filepath)
%Detect best matching import adapter for a source file.

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

if exist(filepath,'file') ~= 2
    error('psyrat_import:fileMissing','Source file was not found: %s',filepath);
end

adapters = psyrat_import_registry();
allDet = struct('adapter',{},'display_name',{},'confidence',{},'reason',{});

for i = 1:numel(adapters)
    out = adapters(i).detect_fn(filepath);
    if isempty(out)
        continue;
    end

    entry = struct();
    entry.adapter = adapters(i).name;
    entry.display_name = adapters(i).display_name;
    entry.confidence = localBoundConfidence(out.confidence);
    entry.reason = out.reason;
    allDet(end+1) = entry; %#ok<AGROW>
end

if isempty(allDet)
    det = struct('ok',false,'adapter','','display_name','', ...
        'confidence',0,'reason','No adapter recognized this file.', ...
        'alternates',struct([]));
    return;
end

[~,idx] = sort([allDet.confidence],'descend');
allDet = allDet(idx);

best = allDet(1);
det = struct('ok',best.confidence > 0, ...
    'adapter',best.adapter, ...
    'display_name',best.display_name, ...
    'confidence',best.confidence, ...
    'reason',best.reason, ...
    'alternates',allDet(2:end));

end

function c = localBoundConfidence(c)
if isempty(c) || ~isnumeric(c) || ~isfinite(c)
    c = 0;
    return;
end
c = max(0,min(1,c));
end
