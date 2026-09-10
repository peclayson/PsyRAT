function lbl = psyrat_trt_occlabel(psyrat_data)
%Occasion-estimand label for test-retest (TRT) reliability outputs.
%
%lbl = psyrat_trt_occlabel(psyrat_data)
%
%Distinguishes the two user-selectable TRT occasion estimands (author
%estimand convention, scientific audit sec.14):
% - the default single-occasion score (n'_o = 1): reliability of a score
%   from n'_i trials at one occasion, generalizing over the occasion facet
% - a multi-occasion composite score (n'_o = k): reliability of a score
%   averaged across k occasions
%
%The two coefficients are never paired with each other's SEM, so they are
%labeled distinctly wherever TRT reliability is reported. Older result
%structs produced before this option existed carry no occasion setting and
%fall back to the single-occasion label.
%
%Inputs
% psyrat_data - the toolbox results struct (uses psyrat_data.relsummary)
%
%Outputs
% lbl - a human-readable label for the selected TRT occasion estimand

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
%

%default to the single-occasion estimand for back-compatibility
noccmode = 1;
nocc = 1;

if isfield(psyrat_data,'relsummary')
    rs = psyrat_data.relsummary;
    if isfield(rs,'noccmode') && ~isempty(rs.noccmode)
        noccmode = rs.noccmode;
    end
    if isfield(rs,'nocc') && ~isempty(rs.nocc)
        nocc = rs.nocc;
    end
end

if noccmode == 2
    lbl = sprintf('Multi-occasion composite TRT reliability, n''_o = %d',nocc);
else
    lbl = sprintf('Single-occasion TRT reliability, n''_o = 1');
end

end
