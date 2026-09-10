function flags = psyrat_copula_flag_counts(cop)
%PSYRAT_COPULA_FLAG_COUNTS Stage-2 diagnostic counts for ONE stratum's modular
%cut-copula output (the struct psyrat_gamma_cut_copula_stage stores in
%REL.out.copula{s}).
%
%Shared by psyrat_dynrel_summary (per-stratum summary.strata(s).copula_flags)
%and psyrat_provenance_lines (the aggregated disclosure line) so the two
%disclosure carriers read the flags with ONE set of semantics and cannot
%drift - the 2026-08-11 code review found the two prior inline copies already
%disagreed on missing-field handling at birth.
%
%Output fields
% recorded       - false when cop is not a struct or has no sel (a malformed or
%                  pre-release stored result). The counts are then all zero and
%                  MUST NOT be read as an affirmative "none flagged": nothing
%                  was summarizable.
% flags_complete - true only when BOTH per-draw flag fields were present. A
%                  missing flag field counts as zero flagged (best-effort, so
%                  an older stored result still summarizes) but clears this,
%                  keeping absence of evidence distinguishable from evidence
%                  of absence.
% n_draws        - numel(cop.sel) when recorded, 0 otherwise.
% n_invalid      - sum of cop.invalidity_flag (numerical trouble: excessive
%                  PIT clipping or copula posterior mass at the rho support
%                  bound), 0 when the field is absent.
% n_calibration  - sum of cop.calibration_flag (descriptive PIT mean/SD
%                  departure), 0 when the field is absent.
%
%Neither flag ever blocks computation (the stage's diagnose-not-hide
%contract); counting them is what makes them visible downstream.

% Copyright (C) 2016-2026 Peter E. Clayson
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

flags = struct('recorded',false,'flags_complete',false,...
    'n_draws',0,'n_invalid',0,'n_calibration',0);
if ~isstruct(cop) || ~isfield(cop,'sel')
    return;
end
flags.recorded = true;
flags.n_draws = numel(cop.sel);
hasinv = isfield(cop,'invalidity_flag');
hascal = isfield(cop,'calibration_flag');
flags.flags_complete = hasinv && hascal;
if hasinv
    flags.n_invalid = sum(cop.invalidity_flag(:));
end
if hascal
    flags.n_calibration = sum(cop.calibration_flag(:));
end

end
