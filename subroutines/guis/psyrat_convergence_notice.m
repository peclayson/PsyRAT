function msg = psyrat_convergence_notice(rel)
%Human-readable notice flagging estimation-diagnostic problems (failed MCMC
%convergence and/or divergent transitions), for surfacing in the results viewer
%so a GUI-only user cannot read an untrustworthy reliability coefficient without
%an on-screen warning.
%
%msg = psyrat_convergence_notice(rel)
%
%Input:
% rel - the reliability results struct (psyrat_data.rel). Two fields are read,
%  both stamped by the estimation core:
%   rel.out.conv.converged   - 0/1 convergence flag from psyrat_checkconv. 0
%    means at least one parameter failed the R-hat/n_eff criteria.
%   rel.out.conv.ndivergent  - count of HMC divergent transitions on the CmdStan
%    path (psyrat_storeconv, a local function inside psyrat_computevarcomp).
%    NaN means "not assessed": native engines (hmcSampler
%    exposes no divergence API; fitlme is not MCMC) and legacy files.
%
%Output:
% msg - a notice string to show the user (e.g. via warndlg), or '' when there is
%  nothing to flag. Both checks are guarded with isfield so legacy results files
%  saved before these fields existed produce no notice (no error), leaving the
%  default workflow unchanged.
%
%Convergence and divergences are reported independently: divergences are
%report-only and never alter the converged flag, so a clean run with divergences
%and a non-converged run are flagged separately and both, when both apply.

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

msg = '';

%nothing to read on a malformed or legacy struct
if ~isstruct(rel) || ~isfield(rel,'out') || ~isstruct(rel.out) ...
        || ~isfield(rel.out,'conv') || ~isstruct(rel.out.conv)
    return;
end
conv = rel.out.conv;

parts = {};

%non-convergence: converged is stamped 0/1 by psyrat_checkconv for every engine
if isfield(conv,'converged') && ~isempty(conv.converged) && conv.converged == 0
    parts{end+1} = ['The estimation model did not converge for at least one ' ...
        'measurement. The reliability estimates shown here may be unreliable. ' ...
        'Consider rerunning with more warmup and sampling iterations.'];
end

%divergent transitions: a positive, non-NaN count (CmdStan path only)
if isfield(conv,'ndivergent') && ~isempty(conv.ndivergent) ...
        && isnumeric(conv.ndivergent) && ~isnan(conv.ndivergent) ...
        && conv.ndivergent > 0
    parts{end+1} = sprintf(['Sampling produced %d divergent transition(s). ' ...
        'Variance components and reliability estimates may be biased even when ' ...
        'R-hat looks acceptable. Consider rerunning with more warmup iterations ' ...
        'or stronger priors.'], conv.ndivergent);
end

if ~isempty(parts)
    msg = strjoin(parts, sprintf('\n\n'));
end

end
