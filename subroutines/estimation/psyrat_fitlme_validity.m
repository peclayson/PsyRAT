function record = psyrat_fitlme_validity(hat,warning_id,warning_message)
%Same admissibility rule for original REML estimates and bootstrap refits.
%Zero variance is an admissible boundary; nonfinite/complex estimates and
%recognized optimizer warnings are not. Records preserve the last warning
%reported by MATLAB (not a claim to capture every warning from the optimizer).
% Copyright (C) 2016-2026 Peter E. Clayson
% SPDX-License-Identifier: GPL-3.0-or-later
record = struct('valid',true,'warning_id',warning_id, ...
    'warning_message',warning_message,'error_id','','error_message','');
markers = {'Singular','IllConditioned','rankDeficient','PerfectFit', ...
    'Converg','Unable','NonConverg'};
if any(contains(warning_id,markers,'IgnoreCase',true))
    record.valid = false;
    record.error_id = 'PsyRAT:fitlmeInvalidWarning';
    record.error_message = 'Optimizer warning invalidated the REML estimates.';
end
names = fieldnames(hat);
for k = 1:numel(names)
    v = hat.(names{k});
    if isnumeric(v) && (~isreal(v) || any(~isfinite(v(:))))
        record.valid = false;
        record.error_id = 'PsyRAT:fitlmeInvalidEstimate';
        record.error_message = 'Nonfinite or complex REML estimate.';
        return;
    end
end
end
