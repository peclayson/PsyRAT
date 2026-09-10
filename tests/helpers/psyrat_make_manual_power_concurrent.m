function [T, P] = psyrat_make_manual_power_concurrent(outfile, seed)
%PSYRAT_MAKE_MANUAL_POWER_CONCURRENT Generate the manual's concurrent power dataset.
%
%   T = PSYRAT_MAKE_MANUAL_POWER_CONCURRENT() regenerates the tutorial dataset
%   test_data/manual_power_concurrent.csv and returns it as a table.
%
%   T = PSYRAT_MAKE_MANUAL_POWER_CONCURRENT(OUTFILE) writes to OUTFILE instead.
%
%   [T, P] = PSYRAT_MAKE_MANUAL_POWER_CONCURRENT(OUTFILE, SEED) draws with a
%   different generator seed (used only by the seed-selection screen described
%   under SEED SELECTION) and also returns the design constants in the struct
%   P, so that the production ledger reads the generating parameters from this
%   file rather than from a copy.
%
% PURPOSE. This file is the dataset behind chapter 15's gamma-family material:
% the exported-table provenance header of section 15.4 and the modular
% cut-posterior disclosure of section 15.5 (also chapter 9, section 9.5). Its
% scores are strictly positive and right-skewed, which is what the gamma
% (scaled chi-square) family is for, and its two events are measured on the
% SAME trials, which is what the residual-coupling design is for. It is fully
% simulated: every generating parameter is listed below. The committed CSV and
% this generator are welded together by tests/TestManualDatasets.m. Do not
% edit the CSV by hand; edit the parameters here, rerun, and re-capture.
%
% DESIGN. Single-trial theta-band and alpha-band power (microvolts squared) at
% a frontal-midline site, both scored from the same error trials of a simulated
% flanker task. 40 participants (ids 501 to 540), one group. Each participant
% contributes n = randi([20 40]) error trials, and every trial yields one theta
% score and one alpha score, so the file holds 2n rows per participant and the
% two events have equal counts within every participant; that equality is the
% concurrent design's requirement, because the toolbox pairs the two events by
% their position within each participant. Each participant also carries one
% questionnaire-style worry score (an integer on a 16 to 80 scale), the
% Dimension 1 covariate, standardized once across the 40 participants.
%
% GENERATING MODEL (the location-scale gamma margins that the modular
% cut-copula workflow fits in its first stage, coupled per trial by a Gaussian
% copula, which is the coupling its second stage estimates):
%
%   log mu_e(p, t)  = a_e + b_e * z(p) + u(p, e) + v(t, e)
%   log sigma_e(p)  = c_e + d_e * z(p) + w(p, e)
%   Y_e(p, t) ~ Gamma(shape = (mu / sigma)^2, scale = sigma^2 / mu)
%                                          so that E(Y) = mu, Var(Y) = sigma^2
%
%   * e indexes the band (1 theta, 2 alpha); z(p) is the standardized worry
%     score, (worry - mean) / SD across the 40 participants with the sample SD.
%   * (u(p, 1), u(p, 2), w(p, 1), w(p, 2)) is 4-variate normal, mean zero, SDs
%     [0.30 0.30 0.25 0.25]: the person's log-mean effects for the two bands and
%     the person's log-residual-SD effects for the two bands. Correlations:
%     0.50 between the two mean effects, 0.35 between the two scale effects,
%     0.30 (theta) and 0.25 (alpha) between a band's mean and scale effect.
%     The scale effects are what make per-participant reliability possible.
%   * v(t, e): trial main effects on the log mean, SDs [0.11 0.13], correlation
%     0.20 across the bands, drawn once per trial index and shared by every
%     participant at that index (persons crossed with trials).
%   * a = log([8 5]) and c = log([4 2.5]): at z = 0 a typical participant's
%     theta power has mean 8 and SD 4, alpha power mean 5 and SD 2.5 (a
%     coefficient of variation of 0.5 in both bands). b = [0.15 0.10] and
%     d = [0.25 0.20] are the dimension slopes on the log mean and on the log
%     residual SD.
%   * The copula: for each trial, (n1, n2) is bivariate standard normal with
%     correlation RHO_E = 0.40; u_e = Phi(n_e); Y_e is the gamma quantile of
%     u_e under that trial's shape and scale. The two bands therefore share
%     trial-level noise with residual correlation 0.40 on the normal scale.
%
% TRUE VALUES. The manual prints no numbers from this file; the planted values
% above are recorded so the production ledger's Session C criteria can be
% judged (the copula residual correlation's credible interval must exclude
% zero; the toolbox's own documentation says its estimate attenuates at
% designs of this size, so the point estimate is expected below 0.40).
%
% REPRODUCIBILITY. rng(SEED, 'twister') below fixes everything. The gamma
% quantile is gammaincinv (a base-MATLAB built-in) times the scale, and Phi is
% written with erfc, so no toolbox is needed. The measurement column is
% rounded to 6 decimals before writing so that the committed CSV and a
% regenerated table agree exactly after a round trip; the rounded scores are
% checked to be strictly positive, which the gamma family requires.
%
% SEED SELECTION (disclosed). The committed seed is the first in the fixed
% sequence 20260829, 20260830, ... whose realization met the manual's
% pre-registered Session C criteria (MANUAL_PRODUCTION.md, section 5: both
% producer fits on this file converged with zero divergent transitions, the
% gamma header carries the family and location-scale lines, and the copula
% residual correlation's credible interval excludes zero). The committed seed,
% 20260829, is position 1. This is a pedagogical selection; it is not a
% validation of the estimator.
%
% This generator is standalone base MATLAB: no toolboxes, no test framework.

% ---- fixed design parameters ------------------------------------------------
SEED    = 20260829;               % committed generator seed (distinct from the
                                  % toolbox's estimation seed, 12345, on purpose)
NSUB    = 40;                     % participants (one group)
EVENTS  = {'theta', 'alpha'};     % the two bands, scored from the same trials
NTRIALS = [20 40];                % error trials per participant, randi range

% log-mean intercepts, dimension slopes on the log mean (per band)
A_MU    = log([8.0, 5.0]);
B_MU    = [0.15, 0.10];

% log-residual-SD intercepts, dimension slopes on the log residual SD (per band)
C_SIG   = log([4.0, 2.5]);
D_SIG   = [0.25, 0.20];

% person block [u_theta u_alpha w_theta w_alpha]: SDs and correlations
SD_PERSON = [0.30, 0.30, 0.25, 0.25];
R_PERSON  = [ 1.00, 0.50, 0.30, 0.00 ;
              0.50, 1.00, 0.00, 0.25 ;
              0.30, 0.00, 1.00, 0.35 ;
              0.00, 0.25, 0.35, 1.00 ];

% trial main effects on the log mean: SDs per band, correlation across bands
SD_TRIAL = [0.11, 0.13];
R_TRIAL  = 0.20;

% Gaussian-copula residual correlation between the two bands on one trial
RHO_E   = 0.40;

% the covariate: an integer questionnaire total, one per participant
WORRY_MEAN  = 45;
WORRY_SD    = 12;
WORRY_RANGE = [16 80];

% ---- design constants, returned for the ledger ------------------------------
P = struct('seed', SEED, 'nsub', NSUB, 'events', {EVENTS}, ...
    'ntrials', NTRIALS, 'a_mu', A_MU, 'b_mu', B_MU, 'c_sig', C_SIG, ...
    'd_sig', D_SIG, 'sd_person', SD_PERSON, 'r_person', R_PERSON, ...
    'sd_trial', SD_TRIAL, 'r_trial', R_TRIAL, 'rho_e', RHO_E, ...
    'worry_mean', WORRY_MEAN, 'worry_sd', WORRY_SD, 'worry_range', WORRY_RANGE);

% ---- output location --------------------------------------------------------
if nargin < 1 || isempty(outfile)
    helpersdir = fileparts(mfilename('fullpath'));        % tests/helpers
    repo = fileparts(fileparts(helpersdir));              % repo root
    outfile = fullfile(repo, 'test_data', 'manual_power_concurrent.csv');
end

% ---- simulate ---------------------------------------------------------------
if nargin < 2 || isempty(seed)
    seed = SEED;                              % the committed realization
end
P.seed = seed;
rng(seed, 'twister');

% Covariate first, standardized across the sample the way the toolbox does.
worry = round(WORRY_MEAN + WORRY_SD * randn(NSUB, 1));
worry = min(max(worry, WORRY_RANGE(1)), WORRY_RANGE(2));
z     = (worry - mean(worry)) / std(worry);

% Person block: one 4-vector per participant from the correlated normal.
Sp = (SD_PERSON' * SD_PERSON) .* R_PERSON;
Lp = chol(Sp, 'lower');
U  = (Lp * randn(4, NSUB))';                  % NSUB x 4

% Trial main effects on the log mean, over the largest possible count, shared
% by every participant at the same trial index.
St = (SD_TRIAL' * SD_TRIAL) .* [1, R_TRIAL; R_TRIAL, 1];
V  = randn(NTRIALS(2), 2) * chol(St);         % maxn x 2

% Copula transform: correlated standard normals to uniforms.
Lc = [1, 0; RHO_E, sqrt(1 - RHO_E^2)];

id = []; event = {}; meas = []; wcol = [];
for p = 1:NSUB
    pid = 500 + p;                            % ids 501..540
    n   = randi(NTRIALS);

    % One correlated normal pair per trial, then Phi to (0, 1).
    N = (Lc * randn(2, n))';                  % n x 2, correlation RHO_E
    Uc = 0.5 * erfc(-N / sqrt(2));

    % The two bands are emitted as one block each, in trial order, so that the
    % toolbox's positional pairing within participant matches trial to trial.
    for e = 1:numel(EVENTS)
        logmu  = A_MU(e) + B_MU(e) * z(p) + U(p, e) + V(1:n, e);
        logsig = C_SIG(e) + D_SIG(e) * z(p) + U(p, 2 + e);
        mu     = exp(logmu);                  % n x 1
        sig    = exp(logsig);                 % scalar for this participant
        shape  = (mu / sig).^2;
        scale  = sig^2 ./ mu;
        y      = gammaincinv(Uc(:, e), shape) .* scale;

        id    = [id;    repmat(pid, n, 1)];              %#ok<AGROW>
        event = [event; repmat(EVENTS(e), n, 1)];        %#ok<AGROW>
        meas  = [meas;  y];                              %#ok<AGROW>
        wcol  = [wcol;  repmat(worry(p), n, 1)];         %#ok<AGROW>
    end
end

% Round the measurements so the committed CSV and a regenerated table agree
% exactly after writetable/readtable, and confirm the gamma precondition.
meas = round(meas, 6);
if any(meas <= 0)
    error('psyrat_make_manual_power_concurrent:nonpositive', ...
        'a rounded score is not strictly positive; choose another seed');
end

T = table(id, event, meas, wcol, ...
    'VariableNames', {'id', 'event', 'meas', 'worry'});

writetable(T, outfile);
fprintf('psyrat_make_manual_power_concurrent: wrote %d rows to %s\n', ...
    height(T), outfile);
end
