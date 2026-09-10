function flags = psyrat_dod_flag_counts(ssrel, ndraws)
%PSYRAT_DOD_FLAG_COUNTS Signed-projection invalidity counts for ONE stratum's
%nonconcurrent difference-of-differences per-participant table (the output of
%PSYRAT_SSREL_DODIFFDYNREL_GAMMA_LS / _TRT_).
%
%   flags = psyrat_dod_flag_counts(ssrel, ndraws)
%
%Shared by psyrat_dynrel_summary (per-stratum summary.strata(s).dod_flags)
%and psyrat_provenance_lines (the aggregated disclosure line) so the two
%disclosure carriers read the flags with ONE set of semantics and cannot
%drift - the same construction as PSYRAT_COPULA_FLAG_COUNTS.
%
%The calculators flag a draw x participant evaluation when a signed composite
%is materially negative or a coefficient leaves [0,1] beyond tolerance
%(record-never-repair; the reference's posterior_invalid_* convention), and
%store the per-participant FRACTION. The counts here are those fractions
%rescaled by the draw count: exact integers up to floating point, restored by
%rounding.
%
%Output fields
% recorded       - false when ssrel is not a table carrying invalid_frac (a
%                  malformed or pre-release stored result). The counts are
%                  then zero and MUST NOT be read as "none flagged".
% flags_complete - true only when BOTH designs' fractions were present
%                  (invalid_frac and invalid_frac_std). A missing column
%                  counts as zero flagged but clears this, keeping absence of
%                  evidence distinguishable from evidence of absence.
% n_evals        - participants x draws when recorded, 0 otherwise.
% n_invalid      - flagged evaluations under the ACTUAL D-study design.
% n_invalid_std  - flagged evaluations under the STANDARDIZED design.
%
%Neither flag ever blocks computation; counting them is what makes them
%visible downstream.
%
%See also PSYRAT_DOD_LABELS, PSYRAT_COPULA_FLAG_COUNTS.

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
    'n_evals',0,'n_invalid',0,'n_invalid_std',0);
if ~istable(ssrel) || ~ismember('invalid_frac',ssrel.Properties.VariableNames) ...
        || ~isscalar(ndraws) || ~isfinite(ndraws) || ndraws < 1
    return;
end
flags.recorded = true;
flags.n_evals = height(ssrel) * ndraws;
flags.n_invalid = round(sum(ssrel.invalid_frac) * ndraws);
hasstd = ismember('invalid_frac_std',ssrel.Properties.VariableNames);
flags.flags_complete = hasstd;
if hasstd
    flags.n_invalid_std = round(sum(ssrel.invalid_frac_std) * ndraws);
end

end
