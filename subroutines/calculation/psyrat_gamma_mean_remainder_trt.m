function R = psyrat_gamma_mean_remainder_trt(alpha, s_p, s_o, s_t, s_po, s_pt, s_ot)
%PSYRAT_GAMMA_MEAN_REMAINDER_TRT Population highest-order (person x occasion x
% trial) OBSERVED-score-scale mean-surface remainder of the crossed Gamma(log
% link) test-retest design, R = K_PIO (Q8/PBA-08 ruling, 2026-09-13).
%
%   R = psyrat_gamma_mean_remainder_trt(alpha, s_p, s_o, s_t, s_po, s_pt, s_ot)
%
% For the log-linked mean
%   log mu_{pot} = alpha + u_p + u_o + u_t + u_po + u_pt + u_ot,
% with independent centered normal effects whose log-scale variances are
% v_A = s_A^2, put V = sum of the six variances and B = exp(2*alpha + V) = E[mu]^2.
% Expanding mu = sqrt(B) * prod_A (1 + D_A) with E[D_A] = 0, E[D_A^2] = q_A =
% expm1(v_A) gives mutually orthogonal products, each assigned to the UNION of
% its factor indices. The highest-order component (union = {p,o,t}) is
%
%   R = B * [ e^V - e^(v_p+v_o+v_po) - e^(v_p+v_t+v_pt) - e^(v_o+v_t+v_ot)
%             + e^(v_p) + e^(v_o) + e^(v_t) - 1 ],
%
% which reduces to B * q_p * q_o * q_t when the three log two-way variances are
% zero. It is the same quantity PSYRAT_GAMMA_VARCOMPS_TRT leaves inside
% sigma_res2 by subtraction (sigma_res2 = R + e_cond_var, up to the float-noise
% clip). It is NOT computed here by that closed form: both the six-way
% subtraction and the inclusion-exclusion bracket cancel catastrophically when
% the log variances are small (the additive limit, where R is many orders of
% magnitude below the terms being differenced). Instead the support-union sum
% is evaluated DIRECTLY - every product of q_A over a factor subset whose index
% union is {p,o,t} - so every term is nonnegative, nothing cancels, and a facet
% whose q's are zero drops out exactly.
%
% WHAT IT IS FOR. Under the Q8 ruling every per-participant two-facet gamma
% residual (analyses 25 and 27, both scale parameterizations) is
%   E_p = R + C_p,
% the POPULATION mean-surface remainder plus the participant's OWN conditional
% observation variance, and the typical-person reference is E_typ = R + C_typ.
% This is the two-facet counterpart of the one-facet S17/S20 pooling of the
% person x trial term (pi_mean). R is a population quantity: it is built from
% alpha (at the participant's own z for the dynamic design) and the six
% log-mean SDs, never from a participant's own effects. It is NOT the person-
% conditional remainder R_p, nor the conditional occasion x trial component
% H_p, and it must not be added to a residual that already carries it.
%
% Inputs (scalars or equal-length column vectors of posterior draws; scalars
% broadcast against vectors), all on the LOG expected-score scale:
%   alpha - log expected-score intercept (at the z of interest)
%   s_p   - person SD;  s_o - occasion SD;  s_t - trial SD
%   s_po  - person x occasion SD;  s_pt - person x trial SD;
%   s_ot  - occasion x trial SD
%
% Output
%   R - [ndraws x 1] observed-score-scale remainder, nonnegative by
%       construction (a sum of nonnegative products; no clip is needed).
%
% Zero-facet identities (pinned by tests/TestGammaTwoFacetResidualPooling):
% removing an entire facet (v_o = v_po = v_ot = 0, or the trial analogue) gives
% R = 0 exactly (every product with that facet in its union carries a zero
% factor), while the person x other-facet component survives on its own.
% Setting only a main effect to zero while its interactions stay does NOT give
% zero, so it is not a zero-facet control.
%
% See also PSYRAT_GAMMA_VARCOMPS_TRT, PSYRAT_GAMMA_VARCOMPS_TRT_LS,
% PSYRAT_SSREL_DYNREL_TRT_GAMMA, PSYRAT_SSREL_DYNREL_TRT_GAMMA_LS.

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

% Column-vector form so scalars and draw vectors mix cleanly.
alpha = alpha(:);
s_p   = s_p(:);
s_o   = s_o(:);
s_t   = s_t(:);
s_po  = s_po(:);
s_pt  = s_pt(:);
s_ot  = s_ot(:);

% Broadcast every input to the common draw length up front (same reason as
% PSYRAT_GAMMA_VARCOMPS_TRT_LS): a scalar/vector mix must return one column of
% the draw length, and a length mismatch must fail loudly rather than broadcast
% into a wrong shape.
lengths = [numel(alpha) numel(s_p) numel(s_o) numel(s_t) numel(s_po) ...
    numel(s_pt) numel(s_ot)];
ndraw = max(lengths);
if any(lengths ~= 1 & lengths ~= ndraw)
    error('psyrat_gamma_mean_remainder_trt:draws',... %Error code and associated error
        ['Posterior-draw inputs must be scalars or all the same length. '...
        'How to fix: pass each of alpha, s_p, s_o, s_t, s_po, s_pt and s_ot '...
        'as either a scalar or a column vector of the same number of draws.']);
end
grow  = ones(ndraw,1);
alpha = alpha .* grow;
s_p   = s_p   .* grow;
s_o   = s_o   .* grow;
s_t   = s_t   .* grow;
s_po  = s_po  .* grow;
s_pt  = s_pt  .* grow;
s_ot  = s_ot  .* grow;

% Squared log-scale SDs (variances of each independent effect on the log mean).
vp  = s_p.^2;
vo  = s_o.^2;
vt  = s_t.^2;
vpo = s_po.^2;
vpt = s_pt.^2;
vot = s_ot.^2;

% V = variance of log mu at a single cell (the same sum as the converters).
V = vp + vo + vt + vpo + vpt + vot;

% Common prefactor B = exp(2*alpha + V) = E[mu]^2.
B = exp(2.*alpha + V);

% Normalized second moments of the centered multiplicative factors,
% q_A = E[D_A^2] = exp(v_A) - 1, one column per factor in the order
% {P, O, T, PO, PT, OT}; expm1 keeps small variances accurate.
q = [expm1(vp) expm1(vo) expm1(vt) expm1(vpo) expm1(vpt) expm1(vot)];   % [ndraws x 6]

% Index set of each factor over the three facets (P = 1, O = 2, T = 3).
factor_idx = {1, 2, 3, [1 2], [1 3], [2 3]};

% Support-union enumeration: sum, over every nonempty subset of the six
% factors whose index union is the full {P,O,T}, of the product of that
% subset's q_A. This is the definition of K_S for S = {p,o,t} (the union of
% the factors multiplied, not merely their count, decides which component a
% product belongs to). 63 subsets; only nonnegative terms are ever added.
total = zeros(ndraw,1);
for mask = 1:63
    members = find(bitand(mask, 2.^(0:5)) > 0);   % factors in this subset
    covered = false(1,3);
    for a = members
        covered(factor_idx{a}) = true;
    end
    if all(covered)
        total = total + prod(q(:,members), 2);
    end
end

R = B .* total;

end
