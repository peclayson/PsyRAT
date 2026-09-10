function line = psyrat_iterline(rel)
%PSYRAT_ITERLINE One provenance line reporting chains and the iteration split.
%
% line = psyrat_iterline(rel)
%
%Returns 'Chains: N, Warmup: W, Sampling: S' for exported-table headers. An
%unqualified 'Iterations' total reads, in Stan reporting convention, as
%post-warmup draws and overstates the run by 2x at the defaults, so the
%combined count is never printed without its label. Legacy results that
%carry only .niter fall back to an explicitly labeled total.
%
%Input
% rel - psyrat_data.rel after estimation (.nchains, .nwarmup, .nsampling,
%  legacy .niter).
%
%Output
% line - char row for the header cell/fprintf pipelines.

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

if isfield(rel,'nwarmup') && isfield(rel,'nsampling') && ...
        ~isempty(rel.nwarmup) && ~isempty(rel.nsampling)
    line = sprintf('Chains: %d, Warmup: %d, Sampling: %d', ...
        rel.nchains, rel.nwarmup, rel.nsampling);
else
    line = sprintf('Chains: %d, Iterations (warmup + sampling): %d', ...
        rel.nchains, rel.niter);
end
end
