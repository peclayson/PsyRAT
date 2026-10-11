function lines = psyrat_bootstrap_lines(rel)
%Per-fit bootstrap count and validity disclosure shared by viewer and exports.
% Copyright (C) 2016-2026 Peter E. Clayson
% SPDX-License-Identifier: GPL-3.0-or-later
lines = {};
if ~isfield(rel,'engine') || ~strcmpi(rel.engine,'fitlme') || ...
        ~isfield(rel,'engine_meta') || ~isfield(rel.engine_meta,'bootstrap_ok')
    return;
end
m = rel.engine_meta;
if isfield(m,'fits'); f = m.fits; else; f = m; end
for k = 1:numel(f)
    label = '';
    if numel(f)>1 && isfield(f,'label'); label = [' [' f(k).label ']']; end
    failed = f(k).bootstrap_reps-f(k).bootstrap_ok;
    line = sprintf('Bootstrap replicates%s: %d of %d succeeded', ...
        label,f(k).bootstrap_ok,f(k).bootstrap_reps);
    if failed>0
        line = sprintf('%s (%d failed refits dropped)',line,failed);
    end
    lines{end+1} = line; %#ok<AGROW>
    if isfield(f,'original') && ~f(k).original.valid
        lines{end+1} = sprintf('Original REML fit%s failed the validity check; its bootstrap does not validate inference.',label); %#ok<AGROW>
    end
end
if isfield(m,'draw_layout')
    lines{end+1} = ['Draw layout: ' m.draw_layout '. Each fit retains its own successful replicate indices; columns do not establish joint cross-fit uncertainty.'];
end
end
