function labels = psyrat_dod_labels(out)
%PSYRAT_DOD_LABELS The estimand-disclosure labels from a REL.out struct of the
%person-specific dynamic NONCONCURRENT difference-of-differences designs
%(analyses 28/29), or [] when absent (any other design, or a pre-release
%stored result).
%
%The labels are run-constants the engine writes once into REL.out.dod_labels
%(the analog of the reference bundle's settings labels): the inference
%framework, the no-dependence-module statement, the person-specific
%scale parameterization, and the estimand version. "Well-formed" requires
%inference_framework - the field every consumer quotes - so the summary and
%the provenance text can never disagree about whether labels exist (the same
%acceptance rule PSYRAT_COPULA_LABELS adopted after the 2026-08-11 review
%found two inline scans disagreeing).
%
%Shared by psyrat_dynrel_summary (summary.dod_labels) and
%psyrat_provenance_lines (the nonconcurrent-DoD disclosure block), which are
%both PRESENCE-keyed on this function's result - never on analysis strings.
%
%See also PSYRAT_COPULA_LABELS, PSYRAT_DOD_FLAG_COUNTS.

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

labels = [];
if ~isstruct(out) || ~isfield(out,'dod_labels')
    return;
end
cand = out.dod_labels;
if isstruct(cand) && isfield(cand,'inference_framework')
    labels = cand;
end

end
