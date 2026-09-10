function ops = psyrat_dynrel_grid()
%Shared grid-shaping helpers for the dynamic-reliability surface calculators.
%
% ops = psyrat_dynrel_grid()
%
%Returns a struct of function handles so the four psyrat_rel_*dynrel* surface
%builders (psyrat_rel_dynrel, psyrat_rel_dynrel_trt, psyrat_rel_diffdynrel,
%psyrat_rel_diffdynrel_trt) share one definition of these byte-identical grid
%utilities instead of each re-declaring them as local functions:
%
%  ops.empty(nz1,nz2)    - preallocate an ll/pt/ul reliability-surface struct of
%                          NaNs (nz1 x nz2 each), with rows indexing z1 and
%                          columns indexing z2.
%  ops.squeeze(g,ndim)   - collapse a one-dimension surface struct to 1 x nz1 row
%                          vectors (no-op when ndim == 2).
%  ops.squeeze1d(v,ndim) - collapse a one-dimension vector to a 1 x nz1 row vector
%                          (no-op when ndim == 2).
%
%These shape the surface so the figure (psyrat_dynrelplot) and viewer consume a
%1 x nz1 row for one dimension and an nz1 x nz2 grid for two, identically across
%all four variants. The helpers are pure (no state), so routing the four builders
%through them is output-identical to the previous per-file local functions.

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

ops = struct('empty',@local_emptygrid,...
    'squeeze',@local_squeezegrid,...
    'squeeze1d',@local_squeeze1d);
end

function g = local_emptygrid(nz1,nz2)
%Preallocate an ll/pt/ul reliability-surface struct of NaNs (nz1 x nz2 each).
g = struct('ll',nan(nz1,nz2),'pt',nan(nz1,nz2),'ul',nan(nz1,nz2));
end

function g = local_squeezegrid(g,ndim)
%For one dimension the grid is a single row; return 1 x nz1 row vectors.
if ndim == 1
    g.ll = g.ll(:)';
    g.pt = g.pt(:)';
    g.ul = g.ul(:)';
end
end

function v = local_squeeze1d(v,ndim)
%For one dimension collapse the vector to a 1 x nz1 row.
if ndim == 1
    v = v(:)';
end
end
