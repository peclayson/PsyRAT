function lines = psyrat_diagnostic_lines(rel)
%Common saved/view/export diagnostic scope and quality disclosure.
% Copyright (C) 2016-2026 Peter E. Clayson
% SPDX-License-Identifier: GPL-3.0-or-later
lines={};
if ~isstruct(rel); return; end
c=struct();
if isfield(rel,'out')&&isfield(rel.out,'conv'); c=rel.out.conv; end
native=isfield(rel,'engine')&&strcmpi(rel.engine,'fitlme');
if native; scope='REML optimizer validity and bootstrap count floor';
else; scope='Early-workflow Rhat/ESS check'; end
verdict='unassessed';
if isfield(c,'converged')
    if isequal(c.converged,1); verdict='passed'; else; verdict='failed'; end
end
lines{end+1,1}=[scope ': ' verdict];
lines{end+1,1}='Final-inference readiness: unassessed (quantity-specific MCSE and uncertainty calibration required)';
if native; return; end
if ~isfield(c,'diagnostics')||isempty(c.diagnostics)
    lines{end+1,1}='Diagnostic source/precision: unrecorded; re-diagnose saved chains'; return;
end
for k=1:numel(c.diagnostics)
    m=c.diagnostics{k};
    text=sprintf('Fit %d diagnostics: %s; precision %g significant digits',k,m.format,m.precision_significant_digits);
    if m.degraded; text=[text '; degraded: ' m.degradation_reason]; end %#ok<AGROW>
    lines{end+1,1}=text;
end
end
