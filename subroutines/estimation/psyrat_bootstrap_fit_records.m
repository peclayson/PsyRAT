function records = psyrat_bootstrap_fit_records(records,fit)
%Collect per-fit evidence before an outer estimation loop replaces its fit.
% Copyright (C) 2016-2026 Peter E. Clayson
% SPDX-License-Identifier: GPL-3.0-or-later
if ~isa(fit,'PsyRATNativeFit') || ~strcmp(fit.engine,'fitlme')
    return;
end
m = fit.meta;
if isfield(m,'fits')
    f = m.fits;
else
    f = m;
end
for k = 1:numel(f)
    r = struct('label','','formula',f(k).formula,'point',f(k).point, ...
        'original',f(k).original,'bootstrap_reps',f(k).bootstrap_reps, ...
        'bootstrap_ok',f(k).bootstrap_ok,'bootstrap_failed',f(k).bootstrap_failed, ...
        'bootstrap_success_mask',f(k).bootstrap_success_mask, ...
        'bootstrap_draw_indices',f(k).bootstrap_draw_indices, ...
        'bootstrap_records',f(k).bootstrap_records, ...
        'bootstrap_fail_ids',{f(k).bootstrap_fail_ids}, ...
        'seed',m.seed,'design',m.design);
    if isempty(records); records = r; else; records(end+1) = r; end %#ok<AGROW>
end
end
