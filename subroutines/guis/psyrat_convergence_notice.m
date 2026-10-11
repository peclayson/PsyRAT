function [msg,requires_ack] = psyrat_convergence_notice(rel)
%Human-readable estimation diagnostic scope and warning notice.
%
%[msg,requires_ack] = psyrat_convergence_notice(rel)
%
%The results struct supplies saved convergence, divergence, exclusion and
%source-quality records. Missing records are disclosed as unrecorded; they are
%not reconstructed or certified by this presentation helper. A malformed
%struct without a convergence record returns an empty notice.
%
%msg always includes assessment scope and final-readiness limits when a
%convergence record exists. requires_ack is true for failed checks, positive
%divergences, recorded exclusions or recorded degraded diagnostic formats.
%Scope-only and legacy-unrecorded information can be displayed without blocking
%the viewer. The optional second output does not change one-output callers.
%Divergences remain report-only and never alter the convergence flag.

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

msg = ''; requires_ack=false;

%nothing to read on a malformed or legacy struct
if ~isstruct(rel) || ~isfield(rel,'out') || ~isstruct(rel.out) ...
        || ~isfield(rel.out,'conv') || ~isstruct(rel.out.conv)
    return;
end
conv = rel.out.conv;

parts = psyrat_diagnostic_lines(rel);

%non-convergence: converged is stamped 0/1 by psyrat_checkconv for every engine
if isfield(conv,'converged') && ~isempty(conv.converged) && conv.converged == 0
    requires_ack=true;
    if isfield(rel,'engine') && strcmpi(rel.engine,'fitlme')
        parts{end+1} = ['The REML validity or bootstrap success-count check failed. ' ...
            'Inspect per-fit optimizer warnings, invalid estimates and refit counts.'];
    else
    parts{end+1} = ['The estimation diagnostic check failed for at least one ' ...
        'measurement. The reliability estimates shown here may be unreliable. ' ...
        'Consider rerunning with more warmup and sampling iterations.'];
    end
end

%divergent transitions: a positive, non-NaN count (CmdStan path only)
if isfield(conv,'ndivergent') && ~isempty(conv.ndivergent) ...
        && isnumeric(conv.ndivergent) && ~isnan(conv.ndivergent) ...
        && conv.ndivergent > 0
    requires_ack=true;
    parts{end+1} = sprintf(['Sampling produced %d divergent transition(s). ' ...
        'Variance components and reliability estimates may be biased even when ' ...
        'R-hat looks acceptable. Consider rerunning with more warmup iterations ' ...
        'and inspect the model geometry.'], conv.ndivergent);
end

%Saved exclusion names remain visible. Historical files may predate the
%model-backed classifier, so the notice does not retroactively certify them.
if isfield(conv,'excluded') && iscell(conv.excluded) && ~isempty(conv.excluded)
    requires_ack=true;
    ex = cellfun(@char, conv.excluded(:)', 'UniformOutput', false);
    parts{end+1} = sprintf(['%d parameter(s) were recorded as excluded from the ' ...
        'convergence check: %s.'], ...
        numel(ex), strjoin(ex, ', '));
end

if isfield(conv,'diagnostics') && iscell(conv.diagnostics)
    for k=1:numel(conv.diagnostics)
        m=conv.diagnostics{k};
        if isstruct(m)&&isfield(m,'degraded')&&isequal(m.degraded,true)&& ...
                isfield(m,'format')&&~strcmp(m.format,'unrecorded')
            requires_ack=true;
        end
    end
end

if ~isempty(parts)
    msg = strjoin(parts, sprintf('\n\n'));
end

end
