function s = psyrat_dod_plugin_cell(bp,bt,er,obs,metric,cvec)
%Evaluate one plug-in cell of the DoD trial-count heatmap and record why it is (not) drawn
%
%s = psyrat_dod_plugin_cell(bp,bt,er,obs,metric)
%s = psyrat_dod_plugin_cell(bp,bt,er,obs,metric,cvec)
%
%Inputs
% bp - 4 x 4 posterior MEAN between-person covariance block
% bt - 4 x 4 posterior MEAN between-trial covariance block
% er - 1 x 4 posterior mean residual VARIANCES (exp(log-SD).^2 averaged over draws)
% obs - [n1 n2 n3 n4] trial counts the cell is evaluated at
% metric - 1 dependability (absolute error), 2 generalizability (relative error)
%
%Optional Inputs
% cvec - signed contrast weights (default [1 -1 -1 1], i.e. (ERP_1 - ERP_2) - (ERP_3 - ERP_4))
%
%Outputs (fields of s)
% value - the plug-in coefficient exactly as the heatmap has always computed it:
%   0 when any count is <= 0 (the zero-count CONVENTION, not a computed
%   coefficient), NaN when the denominator is <= 0 or not finite, otherwise
%   u / denominator with nothing clipped
% u - universe-score contrast variance c' bp c
% bt_scaled - between-trial contrast variance projected to the counts: own
%   count on the diagonal, pairwise HARMONIC mean of the two counts on the
%   off-diagonals (the retained Rocha et al., 2026, Table 6 convention)
% er_scaled - residual contrast variance, each event divided by its own count
% bt_unscaled - between-trial contrast variance at a single observation, c' bt c
%   (the ICC record; it does not depend on the counts)
% denominator - the raw denominator of the selected metric, retained even when
%   it is <= 0 (value is NaN then, the denominator is not)
% zero_count - true when the 0-by-convention branch was taken; every quantity
%   above is NaN on that branch because nothing was computed
% admissible - the heatmap mask value: true = the cell is drawn
% reason - '' for a computed, admissible cell; otherwise the reasons joined by
%   '; '. A zero-count cell carries the convention marker and stays drawn.
%
%What this is
% The plug-in evaluation of the heatmap grid in psyrat_startview_dodiff, moved
% here VERBATIM from that viewer's former local_dod_point_metric so the plotted
% values are bit-identical to before, plus the admissibility record the grid
% never had (2B / G11, audit 2026-09-12 PBA-07 residual). It is a plug-in: each
% variance component enters at its posterior MEAN, and the ratio of means is
% NOT the posterior mean of the draw-level ratio the adjacent table reports
% (two equally weighted draws with universe variance 1 and 9 against error 1
% give .7 for the table and 5/6 here). No draw-level quantity is estimated from
% component means; the record below describes THIS plug-in evaluation only.
%
%Which checks, and why the same four as the table
% psyrat_dodiffrel records four quantities that must be variances -- u,
% bt_scaled, er_scaled and bt_unscaled -- BEFORE it branches on the metric, so
% its record (and the dialog the viewer raises from it) describes the
% evaluation at those counts, whatever coefficient is then formed. This helper
% applies the same four records to the same quantities at this cell's counts,
% under either metric. A negative between-trial contrast variance therefore
% masks a generalizability cell too, although that term is not in the
% generalizability denominator; the reason text says so. The final ratio lying
% in [0,1] proves nothing about the components: the frozen anchor (PSD block
% [[1,2],[2,4]] at counts [1,100,1,1]) has bt_scaled = -.98 and D = .330033.
%
% The rule is the STRICT computed sign, with no tolerance and no clipping, on
% the arithmetic exactly as the loop below performs it. Exact zero is a valid
% boundary (a zero universe variance is a legitimate coefficient of 0; zero
% error components are a legitimate, extreme, coefficient of 1), so zero is
% never a reason here, matching psyrat_admissibility. At the theoretical
% boundary of the anchor (n2 = 2 n1, where bt_scaled is exactly 0 in exact
% arithmetic) floating evaluation can return a tiny negative such as
% -1.39e-17; that cell is then masked, the double itself is kept in the
% bt_scaled field, and the reason prints it at 15 significant digits (readable,
% never rounded to -0; the raw fields, not the text, are the full-precision
% record). That is the honest report of what was computed. Do not add a
% tolerance to hide it.
%
% A component that is not finite is reported as such; it is not a variance and
% the coefficient built from it is not interpretable. The range rule the grid
% has carried since PBA-07 (denominator <= 0, coefficient outside [0,1]) stays
% as a second layer, so nothing that was masked before is drawn now.

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

if nargin < 6 || isempty(cvec)
    cvec = [1; -1; -1; 1];
end
cvec = cvec(:);

s = struct('value',NaN,'u',NaN,'bt_scaled',NaN,'er_scaled',NaN,...
    'bt_unscaled',NaN,'denominator',NaN,'zero_count',false,...
    'admissible',false,'reason','');

%Zero-count convention (unchanged): the cell shows 0 and nothing is computed.
%The components stay NaN so a 0 here cannot be mistaken for a calculated
%variance of zero, and the marker says which branch produced the 0.
if any(obs <= 0)
    s.value = 0;
    s.zero_count = true;
    s.admissible = true;
    s.reason = 'not computed: zero trial count (0 shown by convention)';
    return;
end

%---- plug-in arithmetic, verbatim from the former local_dod_point_metric ----
u = cvec' * bp * cvec;

bt_scaled = 0;
for i = 1:4
    bt_scaled = bt_scaled + (cvec(i)^2) * (bt(i,i) / obs(i));
end

for i = 1:3
    for j = i+1:4
        nharm = 2 / ((1/obs(i)) + (1/obs(j)));
        bt_scaled = bt_scaled + (2 * cvec(i) * cvec(j) * bt(i,j) / nharm);
    end
end

er_scaled = 0;
for i = 1:4
    er_scaled = er_scaled + (cvec(i)^2) * (er(i) / obs(i));
end

if metric == 1
    denom = u + bt_scaled + er_scaled;
else
    denom = u + er_scaled;
end

if denom <= 0 || ~isfinite(denom)
    val = NaN;
else
    val = u / denom;
end
%---- end of the verbatim arithmetic ----------------------------------------

%the single-observation (ICC) record, as psyrat_dodiffrel forms it (:133)
bt_unscaled = cvec' * bt * cvec;

s.value = val;
s.u = u;
s.bt_scaled = bt_scaled;
s.er_scaled = er_scaled;
s.bt_unscaled = bt_unscaled;
s.denominator = denom;

%Reasons. Labels are the ones psyrat_dodiffrel records under, so the grid and
%the table's note name one quantity the same way. The metric notes make the
%evaluation-level nature of the record explicit where a quantity is outside
%the denominator actually drawn.
if metric == 1
    bt_note = '';
else
    bt_note = ['; not in the generalizability denominator, recorded for the ',...
        'evaluation as psyrat_dodiffrel does'];
end
icc_note = '; single-observation record, not in the heatmap denominator';

reasons = {};
reasons = local_component_reason(reasons,'universe-score variance',u,'');
reasons = local_component_reason(reasons,'between-trial contrast variance',...
    bt_scaled,bt_note);
reasons = local_component_reason(reasons,'residual contrast variance',er_scaled,'');
reasons = local_component_reason(reasons,...
    'between-trial contrast variance at a single observation (ICC)',...
    bt_unscaled,icc_note);

%the retained range layer (PBA-07): kept after the component records so a
%cell that fails both says both
if ~isfinite(denom)
    reasons{end+1} = sprintf('denominator not finite (%s)',local_num(denom));
elseif denom <= 0
    reasons{end+1} = sprintf('denominator <= 0 (%s)',local_num(denom));
end
if isfinite(val) && (val < 0 || val > 1)
    reasons{end+1} = sprintf('coefficient outside [0,1] (%s)',local_num(val));
end

s.admissible = isempty(reasons);
s.reason = strjoin(reasons,'; ');

end

function reasons = local_component_reason(reasons,label,v,note)
%STRICT computed sign: negative is a reason, exact zero is not, nonfinite is
%said out loud. No tolerance.
if ~isfinite(v)
    reasons{end+1} = sprintf('%s not finite (%s)%s',label,local_num(v),note);
elseif v < 0
    reasons{end+1} = sprintf('%s negative (%s)%s',label,local_num(v),note);
end
end

function txt = local_num(v)
%15 significant digits: a readable numeric reason in which a tiny raw negative
%such as -1.39e-17 prints as itself and never as -0, and -.98 prints as -0.98.
%Not a round-trip of every double; the raw numeric fields carry those.
txt = sprintf('%.15g',v);
end
