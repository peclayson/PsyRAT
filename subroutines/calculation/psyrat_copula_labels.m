function labels = psyrat_copula_labels(copcell)
%PSYRAT_COPULA_LABELS The stage-archived provenance labels from a REL.out.copula
%cell array (analyses 19/20 modular cut-copula), or [] when no stratum carries
%them (a pre-release stored result).
%
%The labels are run-constants written identically into every stratum's stage
%output; the first stratum that carries a well-formed labels struct wins.
%"Well-formed" requires inference_framework - the field every consumer quotes -
%so the summary and the provenance text can never disagree about whether labels
%exist (the 2026-08-11 code review found the two prior inline scans used
%different acceptance tests).
%
%Shared by psyrat_dynrel_summary (summary.copula_labels) and
%psyrat_provenance_lines (the disclosure line). The stage
%(psyrat_gamma_cut_copula_stage) remains the single source of truth for the
%strings themselves; this function only locates them.

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
if ~iscell(copcell)
    return;
end
for c = 1:numel(copcell)
    cop = copcell{c};
    if isstruct(cop) && isfield(cop,'labels') && ...
            isfield(cop.labels,'inference_framework')
        labels = cop.labels;
        return;
    end
end

end
