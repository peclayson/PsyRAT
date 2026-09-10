function [L, lp] = psyrat_corr_chol_lkj(y, K, eta)
%General K x K LKJ correlation-matrix Cholesky transform for the native HMC engine.
%
%[L, lp] = psyrat_corr_chol_lkj(y, K, eta)
%
%Builds a K x K correlation-matrix Cholesky factor L (lower-triangular, unit row
%norms, so R = L*L' is a valid correlation matrix) from the K(K-1)/2 unconstrained
%canonical partial correlations (CPCs) y, and returns the LKJ(eta) log-density
%over y: the LKJ density on R PLUS the change-of-variables Jacobian from the
%strictly-lower-triangular Cholesky entries to y. This is the general-K
%generalization of the 2x2 atanh(rho) form used by the diff family
%(L_id / L_trl / Lrescor); the > 2-D correlation blocks (e.g. the case-13 joint
%4-D person block coupling location and log-residual events) need it.
%
%Inputs:
% y   - K(K-1)/2 x 1 unconstrained CPC vector, strict-lower-triangle row-major:
%       for i = 2..K, j = 1..i-1 take the next y, i.e.
%       y = [cpc(L[2,1]); cpc(L[3,1]); cpc(L[3,2]); cpc(L[4,1]); ...]. May be a
%       plain double vector or a dlarray (AD path); the construction uses only
%       tanh/sqrt/log and per-row horizontal concatenation (no in-place indexed
%       assignment, which dlarray AD does not support).
% K   - correlation-matrix dimension (>= 1; K = 1 returns L = 1, lp = 0).
% eta - LKJ shape parameter (eta = 1 uniform over correlation matrices; eta = 2
%       matches the diff family's lkj_corr_cholesky(2)).
%
%Outputs:
% L   - K x K lower-triangular correlation Cholesky factor (R = L*L').
% lp  - scalar LKJ(eta) log-density in the UNCONSTRAINED y-space (drops the
%       parameter-independent normalizing constant), ready to add to an HMC
%       log-posterior target.
%
%Derivation (matches Stan's lkj_corr_cholesky up to a constant): with w = tanh(y)
%the CPCs and the row-wise construction L[i,j] = w_ij * s_ij,
%s_ij = sqrt(1 - sum_{k<j} L[i,k]^2), L[i,i] = s_ii, the LKJ(eta) target in
%y-space is
%  sum_{i=2..K} (K - i + 2*eta - 2) * log(L[i,i])   % lkj_corr_cholesky density
%  + sum_{i>j} log(s_ij)                             % |det d(L_free)/dw| part
%  + sum_{i>j} log(1 - w_ij^2)                       % tanh Jacobian (dw/dy)
%The first term is the density on the constrained Cholesky factor; the last two
%are log|det d(L_free)/dy|. This reduces to 2*log(1-rho^2) for K = 2, eta = 2 (the
%validated 2x2 form), and its implied constrained-correlation marginals match the
%known LKJ(eta) Beta marginals (verified CmdStan-free by a prior-marginal Monte
%Carlo check; see tests/TestCorrCholLkj).
%
%See also: psyrat_native_hmc, psyrat_hmc_sample

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

if K <= 1
    if isa(y,'dlarray'); L = dlarray(1); else; L = 1; end
    lp = 0;
    return;
end

isdl = isa(y,'dlarray');
if ~isdl
    %Plain-double fast path. The cell-concatenation construction below exists only
    %because dlarray AD forbids in-place indexed assignment; for ordinary doubles
    %(the always-double posterior-draw transforms, which call this thousands of
    %times) a preallocated indexed-assignment build is far cheaper and numerically
    %identical (scalar tanh/sqrt/log, same accumulation order).
    [L,lp] = local_lkj_double(y,K,eta);
    return;
end
%dlarray (AD) path only past this point: build with per-row concatenation, which
%dlarray AD supports (in-place indexed assignment is not).
one = dlarray(1); zero = dlarray(0);
rows = cell(K,1);
rows{1} = [one, repmat(zero,1,K-1)];   % first row = [1 0 ... 0]
lp = 0;
idx = 0;
for i = 2:K
    rowi = cell(1,K);
    sumsq = zero;                      % running sum of L[i,k]^2 for k < j
    for j = 1:(i-1)
        idx = idx + 1;
        w = tanh(y(idx));
        s = sqrt(1 - sumsq);           % s_ij = sqrt(1 - sum_{k<j} L[i,k]^2)
        Lij = w .* s;
        rowi{j} = Lij;
        sumsq = sumsq + Lij.^2;
        lp = lp + log(s) + log(1 - w.^2);          % Jacobian terms
    end
    Lii = sqrt(1 - sumsq);             % diagonal s_ii
    rowi{i} = Lii;
    for j = (i+1):K
        rowi{j} = zero;
    end
    rows{i} = horzcat(rowi{:});
    lp = lp + (K - i + 2*eta - 2) .* log(Lii);      % lkj_corr_cholesky density
end
L = vertcat(rows{:});
end

function [L,lp] = local_lkj_double(y,K,eta)
%Plain-double LKJ Cholesky build (K >= 2). Identical math to the dlarray
%cell-concatenation path in the parent function, but uses a preallocated K x K
%matrix with indexed assignment instead of per-row horzcat/vertcat of cells.
%Returns the same L and the same y-space LKJ(eta) log-density lp.
L = zeros(K,K);
L(1,1) = 1;                            % first row = [1 0 ... 0]
lp = 0;
idx = 0;
for i = 2:K
    sumsq = 0;                         % running sum of L(i,k)^2 for k < j
    for j = 1:(i-1)
        idx = idx + 1;
        w = tanh(y(idx));
        s = sqrt(1 - sumsq);           % s_ij = sqrt(1 - sum_{k<j} L(i,k)^2)
        Lij = w * s;
        L(i,j) = Lij;
        sumsq = sumsq + Lij^2;
        lp = lp + log(s) + log(1 - w^2);            % Jacobian terms
    end
    Lii = sqrt(1 - sumsq);             % diagonal s_ii
    L(i,i) = Lii;
    lp = lp + (K - i + 2*eta - 2) * log(Lii);       % lkj_corr_cholesky density
end
end
