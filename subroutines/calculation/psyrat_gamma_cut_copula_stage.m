function out = psyrat_gamma_cut_copula_stage(draws, standata, seed, varargin)
%PSYRAT_GAMMA_CUT_COPULA_STAGE Stage 2 of the modular cut-copula workflow: the
% conditional posterior of the Gaussian-copula residual correlation, given
% retained draws from a Stan fit of the Gamma location-scale MARGINS.
%
%   out = psyrat_gamma_cut_copula_stage(draws, standata, seed)
%   out = psyrat_gamma_cut_copula_stage(draws, standata, seed, 'name', value, ...)
%
% WHAT THIS IS, AND WHAT IT IS NOT. The stage-1 Stan program contains the two
% components' exact Gamma marginal likelihoods and every location-scale random
% effect, but NO copula density and no rho_e. This function then estimates the
% copula correlation conditionally on each retained marginal draw. Copula
% information is deliberately never fed back into the marginal parameters, so
% the pair (marginal draw, rho_e draw) is a sample from the CUT (modular)
% posterior
%
%   p_cut(psi, rho | y) = p_margin(psi | y) * p_copula(rho | psi, y),
%
% NOT from the full joint posterior. Marginal uncertainty propagates into rho_e;
% the arrow back from the copula to the margins is severed. This is closely
% related to inference-functions-for-margins methodology. Results must be
% described as modular/two-stage in any methods or limitations section, and the
% provenance labels returned in out.labels must not be stripped from downstream
% output.
%
% WHY IT IS DONE THIS WAY. The full joint likelihood must differentiate two
% incomplete-Gamma CDFs and two inverse-normal transforms for every paired
% observation at every leapfrog step. At the scale this design targets (tens of
% thousands of pairs) that is months of sampling per chain. Moving the copula
% probability transforms outside HMC makes the analysis feasible without
% replacing the person-specific residual SDs with arithmetic sample variances.
% See documentation/modular_cut_copula_workflow.md.
%
% Inputs
%   draws    - struct of posterior draw arrays from the stage-1 margins fit, with
%              the leading dimension indexing draws:
%                alpha_mu     [nd x 2]        component log-mean intercepts
%                alpha_sigma  [nd x 2]        component log-residual-SD intercepts
%                beta_mu      [nd x 2 x K]    dimension slopes on log mean
%                beta_sigma   [nd x 2 x K]    dimension slopes on log residual SD
%                u_p          [nd x P x 4]    person effects; columns 1-2 are the
%                                             two components' location effects,
%                                             3-4 their log-residual-SD effects
%                u_i          [nd x I x 2]    trial effects (location only)
%              and, for the two-facet design only:
%                u_o, u_pi, u_po, u_io        [nd x level x 2] location effects
%   standata - struct describing the data the fit was run on:
%                N, K, P, I                   counts
%                Z            [N x K]         standardized dimension predictors
%                y_1, y_2     [N x 1]         the two components' positive scores
%                p_id, i_id   [N x 1]         person and trial index per row
%                mu_offset, log_sigma_offset  [1 x 2] data-derived centering
%              plus, two-facet only: O, PI, PO, IO and o_id, pi_id, po_id, io_id
%   seed     - seed for THIS invocation; the only random quantity here is the
%              uniform variate used to draw rho_e from its conditional
%              posterior. The estimation call sites derive it per group
%              stratum as run seed + (stratum - 1), so multi-group strata
%              draw decoupled rho_e sequences rather than comonotone ones;
%              the received value is archived in out.settings.seed
%
% Name/value options (defaults match the owner's validated reference run):
%   'max_draws'   500    cap on retained draws carried into the copula stage
%   'grid_points' 20001  grid resolution for the conditional posterior of rho
%   'prior_sd'    0.5    normal(0, prior_sd) prior on rho
%   'clip'        1e-12  probability floor for the Gamma transform
%   'verbose'     0      print progress every 25 draws when nonzero
%
% Output struct out:
%   sel            indices of the input draws that were retained
%   rho_e          one conditional draw of the copula correlation per retained draw
%   rho_mean/sd/median/q025/q975/mode      conditional posterior summaries
%   rho_boundary_mass                      mass pressed against the +/-0.95 bound
%   pit_*                                  normal-score calibration diagnostics
%   invalidity_flag                        numerical failure: excessive clipping
%                                          or material mass at the rho bound
%   calibration_flag                       descriptive PIT mean/SD departure
%   settings, labels                       archived constants and provenance
%
% THE TWO FLAGS ARE DIFFERENT KINDS OF CLAIM and are kept apart deliberately.
% invalidity_flag marks numerical trouble. calibration_flag is a descriptive
% model check: conditional PITs use random effects estimated from the same
% observations, so their plug-in mean and SD are informative but are not a
% posterior-validity proof. Neither blocks computation here; reliability is still
% calculated so that either problem can be diagnosed rather than hidden.
%
% DRAW SELECTION IS DETERMINISTIC - evenly spaced indices, matching the
% reference. It deliberately does NOT follow the seeded randperm convention used
% by psyrat_gamma_extract_diff_copula: even spacing spreads the retained draws
% across chains and iterations without consuming the random stream, and it makes
% the retained set reproducible from the draw count alone.
%
% BASE MATLAB ONLY. See PSYRAT_GAMMA_PIT and PSYRAT_GAMMA_COPULA_RESCOV for the
% tail-accuracy reasoning behind the transforms.
%
% See also PSYRAT_GAMMA_PIT, PSYRAT_COPULA_SUFFSTATS, PSYRAT_COPULA_RHO_GRID.

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

opt = localOptions(varargin);
localValidate(draws, standata);

nd = size(draws.alpha_mu, 1);

% ---- Draw selection ----------------------------------------------------
% Evenly spaced retained draws. Round-half-to-EVEN reproduces R's rounding, so
% the retained set matches the reference for the same draw count rather than
% differing on exact .5 midpoints (MATLAB's round() goes half away from zero).
if nd > opt.max_draws
    sel = unique(localRoundHalfToEven(linspace(1, nd, opt.max_draws)))';
else
    sel = (1:nd)';
end
nsel = numel(sel);

% ---- Random stream -----------------------------------------------------
% The ONLY random quantity in this stage. Everything else - the grid, the
% transforms, the conditional summaries - is deterministic given a draw.
stream = RandStream('mt19937ar', 'Seed', seed);
uniforms = rand(stream, nsel, 1);

% ---- Per-draw storage --------------------------------------------------
z = @() zeros(nsel, 1);
out = struct( ...
    'sel', sel, 'rho_e', z(), 'rho_mean', z(), 'rho_sd', z(), ...
    'rho_median', z(), 'rho_q025', z(), 'rho_q975', z(), 'rho_mode', z(), ...
    'rho_boundary_mass', z(), ...
    'pit_z_1_mean', z(), 'pit_z_1_sd', z(), 'pit_z_2_mean', z(), ...
    'pit_z_2_sd', z(), 'pit_normal_score_correlation', z(), ...
    'pit_1_abs_gt_1_96', z(), 'pit_2_abs_gt_1_96', z(), ...
    'pit_joint_upper_1_96', z(), 'pit_joint_lower_1_96', z(), ...
    'pit_1_lower_clipped', z(), 'pit_1_upper_clipped', z(), ...
    'pit_2_lower_clipped', z(), 'pit_2_upper_clipped', z(), ...
    'pit_any_clip_fraction', z(), ...
    'invalidity_flag', false(nsel, 1), 'calibration_flag', false(nsel, 1));

nobs = standata.N;

for k = 1:nsel
    d = sel(k);

    if opt.verbose && (k == 1 || mod(k, 25) == 0)
        fprintf('  modular cut-copula draw %d/%d\n', k, nsel);
    end

    % ---- Reconstruct every observation's conditional mean and residual SD ----
    % Facets load on the mean; only the person effect and the dimension slopes
    % move the residual SD. See PSYRAT_GAMMA_MARGIN_PREDICTORS.
    [mu_1, sigma_1] = psyrat_gamma_margin_predictors(draws, standata, d, 1);
    [mu_2, sigma_2] = psyrat_gamma_margin_predictors(draws, standata, d, 2);

    pit_1 = psyrat_gamma_pit(standata.y_1, mu_1, sigma_1, opt.clip);
    pit_2 = psyrat_gamma_pit(standata.y_2, mu_2, sigma_2, opt.clip);

    sufficient = psyrat_copula_suffstats(pit_1.z, pit_2.z);
    post = psyrat_copula_rho_grid(sufficient, uniforms(k), opt.grid_points, ...
        opt.prior_sd);

    out.rho_e(k) = post.rho_e;
    out.rho_mean(k) = post.mean;
    out.rho_sd(k) = post.sd;
    out.rho_median(k) = post.median;
    out.rho_q025(k) = post.q025;
    out.rho_q975(k) = post.q975;
    out.rho_mode(k) = post.mode;
    out.rho_boundary_mass(k) = post.boundary_mass;

    % ---- Normal-score calibration diagnostics ----
    % If the fitted margins describe the data, the normal scores are standard
    % normal: mean 0, SD 1, ~5% beyond +/-1.96. The joint tail frequencies are
    % reported because copula dependence shows up there first.
    z1 = pit_1.z;
    z2 = pit_2.z;
    out.pit_z_1_mean(k) = mean(z1);
    out.pit_z_1_sd(k) = std(z1);
    out.pit_z_2_mean(k) = mean(z2);
    out.pit_z_2_sd(k) = std(z2);
    out.pit_normal_score_correlation(k) = localCorrelation(z1, z2);
    out.pit_1_abs_gt_1_96(k) = mean(abs(z1) > 1.96);
    out.pit_2_abs_gt_1_96(k) = mean(abs(z2) > 1.96);
    out.pit_joint_upper_1_96(k) = mean(z1 > 1.96 & z2 > 1.96);
    out.pit_joint_lower_1_96(k) = mean(z1 < -1.96 & z2 < -1.96);

    out.pit_1_lower_clipped(k) = pit_1.lower_clipped;
    out.pit_1_upper_clipped(k) = pit_1.upper_clipped;
    out.pit_2_lower_clipped(k) = pit_2.lower_clipped;
    out.pit_2_upper_clipped(k) = pit_2.upper_clipped;
    clipped = pit_1.lower_clipped + pit_1.upper_clipped + ...
        pit_2.lower_clipped + pit_2.upper_clipped;
    out.pit_any_clip_fraction(k) = clipped / (2 * nobs);

    % ---- Flags (recorded, never blocking) ----
    out.calibration_flag(k) = abs(out.pit_z_1_mean(k)) > 0.10 || ...
        abs(out.pit_z_2_mean(k)) > 0.10 || ...
        out.pit_z_1_sd(k) < 0.90 || out.pit_z_1_sd(k) > 1.10 || ...
        out.pit_z_2_sd(k) < 0.90 || out.pit_z_2_sd(k) > 1.10;
    out.invalidity_flag(k) = out.pit_any_clip_fraction(k) > 1e-4 || ...
        post.boundary_mass > 0.01;
end

% ---- Archived constants and provenance --------------------------------
% Every tunable that could move the estimand is recorded with the result, so a
% stored fit carries the configuration it was produced under.
out.settings = struct( ...
    'max_draws', opt.max_draws, 'grid_points', opt.grid_points, ...
    'prior_sd', opt.prior_sd, 'clip', opt.clip, 'seed', seed, ...
    'n_draws_available', nd, 'n_draws_used', nsel, ...
    'draw_selection', 'evenly_spaced_deterministic', ...
    'rho_support', [-0.95 0.95]);

% These three strings state what kind of posterior produced the numbers. They
% must ride through to every output surface; do not remove them in reports.
out.labels = struct( ...
    'inference_framework', ...
        'modular_cut_two_stage_gamma_margins_gaussian_copula_v1', ...
    'copula_feedback_to_margins', false, ...
    'joint_posterior_status', 'not_a_full_joint_bayesian_posterior');

end

% ---------------------------------------------------------------------------
% Local helpers
% ---------------------------------------------------------------------------

function r = localCorrelation(a, b)
%LOCALCORRELATION Pearson correlation without the Statistics Toolbox.
a = a(:) - mean(a);
b = b(:) - mean(b);
denominator = sqrt(sum(a.^2) * sum(b.^2));
if denominator <= 0
    r = NaN;
else
    r = sum(a .* b) / denominator;
end
end

function r = localRoundHalfToEven(x)
%LOCALROUNDHALFTOEVEN Round to nearest, ties to even (R's rounding).
%   MATLAB's round() breaks ties away from zero, so an index sequence that lands
%   on an exact .5 would select a different draw than the reference does.
r = round(x);
half = abs(x - fix(x)) == 0.5;
r(half) = 2 * round(x(half) / 2);
end

function opt = localOptions(args)
%LOCALOPTIONS Name/value options with the reference-validated defaults.
opt = struct('max_draws', 500, 'grid_points', 20001, 'prior_sd', 0.5, ...
    'clip', 1e-12, 'verbose', 0);

if mod(numel(args), 2) ~= 0
    error('psyrat_gamma_cut_copula_stage:options', ...
        'Options must be supplied as name/value pairs.');
end
for k = 1:2:numel(args)
    name = args{k};
    if ~ischar(name) && ~isstring(name)
        error('psyrat_gamma_cut_copula_stage:options', ...
            'Option names must be character vectors.');
    end
    name = char(name);
    if ~isfield(opt, name)
        error('psyrat_gamma_cut_copula_stage:options', ...
            'Unrecognized option ''%s''. Valid options: %s.', ...
            name, strjoin(fieldnames(opt)', ', '));
    end
    opt.(name) = args{k+1};
end

if ~isscalar(opt.max_draws) || opt.max_draws < 1
    error('psyrat_gamma_cut_copula_stage:options', ...
        'max_draws must be a positive scalar.');
end
%clip is re-validated by psyrat_gamma_pit per call, but rejecting it here means
%the failure names the documented public option instead of surfacing per draw.
if ~isscalar(opt.clip) || ~isfinite(opt.clip) || opt.clip <= 0 || opt.clip >= 0.5
    error('psyrat_gamma_cut_copula_stage:options', ...
        'clip must be a finite scalar strictly between 0 and 0.5.');
end

end

function localValidate(draws, standata)
%LOCALVALIDATE Shape agreement between the draws and the data they came from.

needed_draws = {'alpha_mu', 'alpha_sigma', 'beta_mu', 'beta_sigma', 'u_p', 'u_i'};
missing = needed_draws(~isfield(draws, needed_draws));
if ~isempty(missing)
    error('psyrat_gamma_cut_copula_stage:draws', ...
        ['The draw struct is missing the field(s): %s. How to fix: check the '...
        'extraction that builds it from the stage-1 fit.'], ...
        strjoin(missing, ', '));
end

needed_data = {'N', 'K', 'Z', 'y_1', 'y_2', 'p_id', 'i_id', 'mu_offset', ...
    'log_sigma_offset'};
missing = needed_data(~isfield(standata, needed_data));
if ~isempty(missing)
    error('psyrat_gamma_cut_copula_stage:standata', ...
        'The standata struct is missing the field(s): %s.', ...
        strjoin(missing, ', '));
end

nd = size(draws.alpha_mu, 1);
if size(draws.beta_mu, 1) ~= nd || size(draws.u_p, 1) ~= nd || ...
        size(draws.u_i, 1) ~= nd || size(draws.alpha_sigma, 1) ~= nd || ...
        size(draws.beta_sigma, 1) ~= nd
    error('psyrat_gamma_cut_copula_stage:draws', ...
        'Every draw array must share the same leading (draw) dimension.');
end
if size(draws.u_p, 3) ~= 4
    error('psyrat_gamma_cut_copula_stage:draws', ...
        ['The person block must carry four columns per draw: two location '...
        'effects then two log-residual-SD effects.']);
end
if size(draws.beta_mu, 3) ~= standata.K || size(draws.beta_sigma, 3) ~= standata.K
    error('psyrat_gamma_cut_copula_stage:draws', ...
        'The dimension slopes must have one column per predictor (K = %d).', ...
        standata.K);
end
if size(standata.Z, 1) ~= standata.N || size(standata.Z, 2) ~= standata.K
    error('psyrat_gamma_cut_copula_stage:standata', ...
        'Z must be N-by-K.');
end
if numel(standata.y_1) ~= standata.N || numel(standata.y_2) ~= standata.N ...
        || numel(standata.p_id) ~= standata.N || numel(standata.i_id) ~= standata.N
    error('psyrat_gamma_cut_copula_stage:standata', ...
        ['The two score vectors and the person/trial index vectors must each '...
        'have one entry per observation. Concurrent pairing means row j of '...
        'each component is the SAME physical observation.']);
end

% Two-facet data need their four extra location effects, and each index vector
% needs its matching draw array. Catching this here names the missing piece;
% without it the reconstruction would fail on a bare field-access error.
if isfield(standata, 'o_id') && ~isempty(standata.o_id)
    pairs = {'o_id', 'u_o'; 'pi_id', 'u_pi'; 'po_id', 'u_po'; 'io_id', 'u_io'};
    for k = 1:size(pairs, 1)
        if ~isfield(standata, pairs{k,1}) || ~isfield(draws, pairs{k,2})
            error('psyrat_gamma_cut_copula_stage:twofacet', ...
                ['The two-facet design needs both standata.%s and draws.%s. '...
                'How to fix: supply all four of the occasion, person x trial, '...
                'person x occasion and trial x occasion effects, or omit '...
                'o_id entirely for the one-facet design.'], ...
                pairs{k,1}, pairs{k,2});
        end
        if numel(standata.(pairs{k,1})) ~= standata.N
            error('psyrat_gamma_cut_copula_stage:standata', ...
                'standata.%s must have one entry per observation.', pairs{k,1});
        end
    end
end

end
