function [mu, sigma] = psyrat_gamma_margin_predictors(draws, standata, d, m)
%PSYRAT_GAMMA_MARGIN_PREDICTORS Conditional mean and residual SD of every
% observation, for one posterior draw of one component of the Gamma
% location-scale margins model.
%
%   [mu, sigma] = psyrat_gamma_margin_predictors(draws, standata, d, m)
%
% Inputs
%   draws    - stage-1 posterior draw arrays (see PSYRAT_GAMMA_CUT_COPULA_STAGE
%              for the field contract)
%   standata - the data the fit was run on, same contract
%   d        - index of the draw to reconstruct
%   m        - component index, 1 or 2
%
% Outputs
%   mu, sigma - N-by-1 conditional mean and residual SD, in observation order
%
% THE TWO SUBMODELS ARE NOT SYMMETRIC, and that asymmetry is the model rather
% than an omission:
%
%   log mu    = offset + intercept + dimension slopes + person location effect
%               + EVERY facet effect (trial; plus occasion, person x trial,
%               person x occasion and trial x occasion in the two-facet design)
%   log sigma = offset + intercept + dimension slopes + person SCALE effect
%
% Facets load on the mean only. The residual SD varies across people (that is
% what makes person-specific reliability possible) and with the dimension
% predictors (that is what makes it dynamic), but not across trials or
% occasions.
%
% The person block carries four columns per draw in a fixed order: the two
% components' location effects first, then their two log-residual-SD effects.
% Reading columns 3 and 4 as locations, or transposing the block, produces a
% model that still runs and still returns plausible reliabilities, so this
% mapping is pinned by a test against the reference implementation rather than
% left to inspection.
%
% The offsets are DATA-DERIVED centering constants carried in standata, not
% parameters: the fitted intercepts are deviations from them, so omitting them
% silently rescales every mean and residual SD. This is the trap recorded for
% the earlier external-validation harness in
% documentation/gamma_locationscale_rollout.md.
%
% Separated from PSYRAT_GAMMA_CUT_COPULA_STAGE so it can be pinned directly, and
% because the person-specific reliability stage needs the same reconstruction.
%
% See also PSYRAT_GAMMA_CUT_COPULA_STAGE, PSYRAT_GAMMA_PIT.

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

if ~isscalar(m) || ~any(m == [1 2])
    error('psyrat_gamma_margin_predictors:component', ...
        'The component index must be 1 or 2.');
end

twofacet = isfield(standata, 'o_id') && ~isempty(standata.o_id);

% Person block: columns 1-2 are the location effects, 3-4 the log-residual-SD
% effects, so component m reads column m for the mean and column 2+m for scale.
u_p = reshape(draws.u_p(d, :, :), size(draws.u_p, 2), 4);

% ---- Location submodel -------------------------------------------------
beta_mu = reshape(draws.beta_mu(d, m, :), [], 1);
eta = standata.mu_offset(m) + draws.alpha_mu(d, m) + standata.Z * beta_mu ...
    + u_p(standata.p_id, m) ...
    + localFacet(draws.u_i, d, m, standata.i_id);

if twofacet
    eta = eta ...
        + localFacet(draws.u_o, d, m, standata.o_id) ...
        + localFacet(draws.u_pi, d, m, standata.pi_id) ...
        + localFacet(draws.u_po, d, m, standata.po_id) ...
        + localFacet(draws.u_io, d, m, standata.io_id);
end

% ---- Scale submodel ----------------------------------------------------
beta_sigma = reshape(draws.beta_sigma(d, m, :), [], 1);
log_sigma = standata.log_sigma_offset(m) + draws.alpha_sigma(d, m) ...
    + standata.Z * beta_sigma + u_p(standata.p_id, 2 + m);

mu = exp(eta);
sigma = exp(log_sigma);

if ~all(isfinite(mu)) || any(mu <= 0) || ~all(isfinite(sigma)) || any(sigma <= 0)
    error('psyrat_gamma_margin_predictors:reconstruction', ...
        ['A reconstructed conditional mean or residual SD is nonpositive or '...
        'non-finite for draw %d, component %d. That is a slot-mapping or '...
        'overflow failure rather than a data problem: both are exponentials '...
        'of a linear predictor and cannot be nonpositive unless the predictor '...
        'itself is corrupt.'], d, m);
end

end

% ---------------------------------------------------------------------------
% Local helpers
% ---------------------------------------------------------------------------

function contribution = localFacet(effect, d, m, index)
%LOCALFACET One two-column facet effect, expanded to the observation rows.
levels = size(effect, 2);
slice = reshape(effect(d, :, :), levels, 2);
contribution = slice(index, m);
end
