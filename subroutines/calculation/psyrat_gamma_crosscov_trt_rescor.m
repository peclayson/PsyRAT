function cov = psyrat_gamma_crosscov_trt_rescor(alpha_x, s_p_x, s_o_x, s_t_x, ...
    s_po_x, s_pt_x, s_ot_x, alpha_y, s_p_y, s_o_y, s_t_y, ...
    s_po_y, s_pt_y, s_ot_y, cor_p, cor_o, cor_t, cor_po, cor_pt, cor_ot, ...
    cov_e_obs)
%PSYRAT_GAMMA_CROSSCOV_TRT_RESCOR Cross-cell observed-score covariance components
% for the CONCURRENT two-facet (trial x occasion) gamma difference designs,
% including the residual term the non-concurrent helper fixes at zero.
%
%   cov = psyrat_gamma_crosscov_trt_rescor(alpha_x, s_p_x, s_o_x, s_t_x, ...
%             s_po_x, s_pt_x, s_ot_x, alpha_y, s_p_y, s_o_y, s_t_y, ...
%             s_po_y, s_pt_y, s_ot_y, cor_p, cor_o, cor_t, cor_po, cor_pt, ...
%             cor_ot, cov_e_obs)
%
% Inputs are the twenty arguments of PSYRAT_GAMMA_CROSSCOV_TRT, plus:
%   cov_e_obs - the OBSERVATION-level residual cross-cell covariance on the
%               observed score scale, supplied by the caller (draws vector or
%               scalar). For the modular cut-copula designs this is the
%               Gauss-Hermite conversion of the copula correlation, rescaled by
%               the two components' expected means; see
%               PSYRAT_GAMMA_COPULA_RESCOV.
%
% Output struct cov. All eight of PSYRAT_GAMMA_CROSSCOV_TRT's fields are passed
% through unchanged and under the same names (seven cross-covariances plus the
% total):
%   cov_p, cov_o, cov_t     - main-effect cross-cell covariances
%   cov_po, cov_pt, cov_ot  - two-way interaction cross-cell covariances
%   cov_res                 - MEAN-SURFACE three-way (p x o x t) cross term
%   cov_mu_total            - total mean-surface cross-cell covariance
% and this file adds:
%   cov_e_obs               - the observation-level residual cross term. The
%                             supplied value, expanded to the draw length when a
%                             scalar was passed.
%   cov_res_e_obs           - cov_res + cov_e_obs, the pooled residual channel
%
% WHY THIS EXISTS RATHER THAN AN EDIT TO PSYRAT_GAMMA_CROSSCOV_TRT. That helper
% serves the NON-concurrent two-facet designs, where the scaled-chi-square
% observations of two distinct cells are independent and its cov_res is therefore
% the complete residual cross-covariance. Changing it would move shipped numbers
% (analysis 10). This file leaves it untouched, calls it for the seven terms it
% already computes - so every shared quantity is bit-identical to the shipped
% path - and adds only what concurrency introduces. Same construction as the
% one-facet PSYRAT_GAMMA_CROSSCOV_RESCOR.
%
% WHAT IS POOLED HERE, AND WHAT DELIBERATELY IS NOT. Rast & Clayson (in press)
% page 4 specifies BOTH designs, and they differ. For one facet the residual is
% a combined quantity - "the combined person-by-task interaction and
% unsystematic measurement error". For two facets the same page extends the
% coefficient to person x task x occasion, keeps sigma^2_pt and sigma^2_po as
% their own components, and defines the residual separately: "the term
% sigma^2_pto,e is the residual variance, encompassing the three-way interaction
% and other unmodelled noise", with
%
%   G = sigma^2_p / (sigma^2_p + sigma^2_pt/n'_t + sigma^2_po/n'_o
%                              + sigma^2_pto,e/(n'_t*n'_o)).
%
% So the split implemented here is the paper's own two-facet decomposition, not
% an adaptation of its one-facet one. CITE IT BY PAGE: the two-facet display is
% unnumbered, and the only numbered equation on page 4 is (2), which is the
% ONE-facet G - so "Equation (2)" would name the wrong design. Do NOT cite the
% paper's page-7 sentence here either: it concerns the one-facet
% condition-specific residual sigma^2_pt,e(t), says only that it "may include a
% person-by-condition interaction", and is about conditions rather than
% occasions.
%
% Concretely, realizing the residual is not the same operation in the two
% designs, and reading the one-facet sibling across without this paragraph would
% be wrong:
%
%   ONE-FACET. Person x trial and measurement error are confounded by design -
%   there is one observation per person x trial cell - so the toolbox carries
%   them in a single slot (PSYRAT_GAMMA_VARCOMPS_LS's sigma_pi_e2), and the
%   pooling question is only whether the mean-surface half is dropped.
%
%   TWO-FACET. Person x trial is SEPARATELY ESTIMABLE, is reported as its own
%   component, and PSYRAT_DIFFREL_TRT gives it its own divisor n_i - different
%   from the residual's n_i*n_o. Folding it into the residual would change the
%   D-study arithmetic, not just a label. It therefore stays out. The quantity
%   that pools here is the highest-order term: the three-way mean-surface
%   remainder plus the observation dispersion.
%
% The reference implementation makes exactly this split (R/reliability.R:241-246
% selects residual_name = "pio_e" and mean_name = "pio_mean" for the two-facet
% branch where the one-facet branch selects "pi_e"/"pi_mean"; :266-268 then form
% covariance$pio_e = covariance$pio_mean + residual_covariance). cov_res_e_obs is
% that line.
%
% Its diagonal counterpart is PSYRAT_GAMMA_VARCOMPS_TRT_LS's sigma_res2, which
% that function recovers by subtraction and which is therefore always
% pio_mean plus whatever observation variance the scale argument encodes. WHICH
% observation variance depends on the caller, and the two callers differ:
%   - PSYRAT_SSREL_DIFFDYNREL_TRT_RESCOR_GAMMA_LS passes the PARTICIPANT'S own
%     log residual SD with s_sd = 0, giving pio_mean + sigma_person^2 - the
%     reference's per-person quantity.
%   - PSYRAT_REL_DIFFDYNREL_TRT_GAMMA_LS passes the POPULATION log residual SD
%     with s_sd = the person log-sigma SD, giving pio_mean + E[sigma_p^2]
%     marginalized over the person population, which is the typical-person
%     surface's convention.
% Do not describe either as the diagonal counterpart without saying which.
%
% NUMERICS - A MEASURED LIMITATION, NOT A SILENT FIX. cov_res arrives from the
% shipped helper as an eight-term alternating inclusion-exclusion,
%
%   P^-1 * cov_res = ABCDEF - ABD - ACE - BCF + A + B + C - 1
%
% writing P for the shared prefactor and A = exp(cov_p_log) and so on. As the
% log-scale covariances shrink toward zero EVERY one of those eight terms tends
% to +/-1, while their sum is second order in the covariances - at a common
% log-covariance of 1e-3 the terms are all ~1.00 and the sum is 6.0e-6, a
% cancellation of about 1.7e5. That is where the significant digits go. An
% algebraically identical regrouping,
%
%   P^-1 * cov_res = A*B*D*expm1(c+e+f) - A*expm1(c+e) - B*expm1(c+f) + expm1(c)
%
% brings the terms down to first order, so at the same point they are ~2e-3
% against a 6.0e-6 sum: still a cancellation, but ~5e2 rather than ~1.7e5. It is
% BETTER CONDITIONED, not exact - it loses a couple of digits itself. It is NOT
% used here. The eight-term form is what analysis 10 ships with, and replacing it
% would silently move those numbers; the one-facet sibling chose its product form
% only because its mean-surface cross term did not exist in the shipped helper at
% all and had to be written fresh.
%
% The size of the discrepancy was MEASURED, not assumed
% (TestGammaDiffDynrelTrtLsConversion/testResidualCrossCovCancellation holds it
% to a ceiling just above each measurement): at log-scale covariances of 1e-3 the
% two forms agree to 3.8e-12 relative, at 1e-4 to 2.8e-9, and at 1e-5 to only
% 1.3e-7. Realistic ERP designs sit near 1e-3 - facet
% SDs around 0.05-0.10 with cross-event correlations of 0.2-0.5 - so the loss is
% immaterial where this function is actually used, and becomes visible only for
% very weak facet correlations or very small facet SDs. That is why the shipped
% form is kept. If a design ever lands in the degraded regime, it is a finding
% for the maintainer to adjudicate against the papers (Golden Rule 1), not a
% reason to edit a shipped helper from here.
%
% See also PSYRAT_GAMMA_CROSSCOV_TRT, PSYRAT_GAMMA_CROSSCOV_RESCOR,
% PSYRAT_GAMMA_COPULA_RESCOV, PSYRAT_GAMMA_VARCOMPS_TRT_LS.

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

% Every shared term comes from the shipped helper, unchanged, so the concurrent
% and non-concurrent two-facet paths cannot drift apart on them.
base = psyrat_gamma_crosscov_trt(alpha_x, s_p_x, s_o_x, s_t_x, ...
    s_po_x, s_pt_x, s_ot_x, alpha_y, s_p_y, s_o_y, s_t_y, ...
    s_po_y, s_pt_y, s_ot_y, cor_p, cor_o, cor_t, cor_po, cor_pt, cor_ot);

cov_e_obs = cov_e_obs(:);
if isscalar(cov_e_obs)
    cov_e_obs = repmat(cov_e_obs, size(base.cov_res));
end
if numel(cov_e_obs) ~= numel(base.cov_res)
    error('psyrat_gamma_crosscov_trt_rescor:size', ...
        ['The residual cross-covariance must be a scalar or match the number '...
        'of draws (%d); got %d.'], numel(base.cov_res), numel(cov_e_obs));
end

cov = base;
cov.cov_e_obs = cov_e_obs;
% The pooled residual channel: the three-way mean-surface remainder plus the
% observation-level covariance. See the header for why person x trial is absent.
cov.cov_res_e_obs = base.cov_res + cov_e_obs;

end
