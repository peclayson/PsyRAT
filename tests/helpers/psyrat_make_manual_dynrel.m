function [T, P] = psyrat_make_manual_dynrel(outfile, seed)
%PSYRAT_MAKE_MANUAL_DYNREL Generate the manual's dynamic-reliability dataset.
%
%   T = PSYRAT_MAKE_MANUAL_DYNREL() regenerates the tutorial dataset
%   test_data/manual_ern_dynrel.csv and returns it as a table.
%
%   T = PSYRAT_MAKE_MANUAL_DYNREL(OUTFILE) writes to OUTFILE instead.
%
%   [T, P] = PSYRAT_MAKE_MANUAL_DYNREL(OUTFILE, SEED) draws with a different
%   generator seed (used only by the seed-selection screen described under
%   SEED SELECTION) and also returns the design constants in the struct P, so
%   that the manual's truth statements and the production ledger read the
%   generating parameters from this file rather than from a copy.
%
% PURPOSE. This file is the dataset behind the manual's dynamic (conditional)
% reliability figures (chapter 9, section 9.5, and chapter 15, section 15.2):
% a reliability surface drawn over a between-person covariate. It is fully
% simulated: every generating parameter is listed below. The committed CSV and
% this generator are welded together by tests/TestManualDatasets.m, which
% regenerates the table and requires it to match the committed file. Do not
% edit the CSV by hand; edit the parameters here, rerun, and re-capture the
% figures.
%
% DESIGN. A simulated flanker study measuring error-related negativity (ERN),
% scored as a mean amplitude in microvolts. One session; 80 participants in
% two independent groups (GroupA, GroupB; between-person), 40 per group; two
% events, inc_err and inc_cor (the two events the single-session tutorial
% selects), with trial counts drawn per participant per event as in the
% single-session dataset (errors randi([12 36]), corrects randi([24 40])).
% Each participant also carries ONE questionnaire-style score, worry, an
% integer on a 16 to 80 scale (the range of a 16-item, 1 to 5 scale), repeated
% down that participant's rows. That column is the Dimension 1 covariate. The
% toolbox standardizes it once across all 80 participants before fitting, so
% z = 0 is the sample mean of worry and z = 1 is one sample SD above it.
%
% GENERATING MODEL (the location-scale model psyrat_build_stan_dynrel fits,
% one fit per group x event stratum, persons crossed with trials):
%
%   meas(p, e, t) = mu(g, e) + b_mean(e) * z(p) + person(p, e) + trial(g, e, t)
%                 + exp(log_s0(g, e) + b_sig * z(p) + delta(p)) * eps(p, e, t)
%
%   * z(p): the participant's standardized worry score, (worry - mean) / SD
%     across the 80 participants with the sample SD, which is exactly the
%     standardization the toolbox applies, so the slopes below are stated on
%     the axis the surface is drawn over.
%   * mu(g, e), the person-effect SDs sigma_p, the trial main-effect SD
%     sigma_i = 1.0, and the residual SDs at z = 0 (s0) are the single-session
%     generator's inc_err and inc_cor cells; GroupB's error trials are noisier.
%     person(p, :) is bivariate normal with r = 0.75 between the two events.
%   * b_mean: -0.5 microvolts per SD of worry on error trials (higher worry,
%     more negative ERN), 0 on correct trials. It moves the mean surface only.
%   * b_sig = +0.25 on the log residual SD in every stratum: one SD of worry
%     multiplies the trial-to-trial SD by exp(0.25) = 1.28. This is the point
%     the figure makes: reliability falls as worry rises, because the same
%     person variance sits behind a growing error variance.
%   * delta(p): one person-level log-residual effect per participant, SD 0.20,
%     shared by both events and independent of person(p, :). It separates the
%     "typical person" surface (delta = 0) from the population average, and it
%     is what the model's sigma_delta_p term estimates.
%
% TRUE VALUES. With one fit per stratum the toolbox's generalizability surface
% at n trials is G(z) = sigma_p^2 / (sigma_p^2 + sigma_e(z)^2 / n) with
% sigma_e(z) = s0 * exp(0.25 z); dependability adds sigma_i^2 / n to the error
% term and sits about .01 lower. At the median drawn trial counts (about 24
% error trials, 32 correct trials):
%
%     stratum            sigma_p   s0    G(z=-2)   G(z=0)   G(z=+2)
%     GroupA inc_err       2.0     5.5     .90       .76      .54
%     GroupB inc_err       2.2     7.5     .85       .67      .43
%     inc_cor (both)       1.7     5.0     .91       .79      .58
%
%   So every panel of the surface falls by roughly .3 to .4 across the
%   observed range of worry, which is the slope the chapter teaches the
%   reader to look for, and the two groups share one predictor axis.
%
% REPRODUCIBILITY. rng(SEED, 'twister') below fixes everything. The
% measurement column is rounded to 6 decimals before writing so that the
% committed CSV and a regenerated table agree exactly after a round trip.
%
% SEED SELECTION (disclosed). The committed seed is the first in the fixed
% sequence 20260829, 20260830, ... whose realization met the manual's
% pre-registered Session C criteria (MANUAL_PRODUCTION.md, section 5: the
% producer fit converged with zero divergent transitions, and in each error
% stratum the generalizability surface at the median trial count falls by at
% least .10 from the bottom of the covariate grid to the top). The committed
% seed, 20260829, is position 1. This is a pedagogical selection so that the
% figure shows the lesson; it is not a validation of the estimator, and the
% manual says so.
%
% This generator is standalone base MATLAB: no toolboxes, no test framework.

% ---- fixed design parameters ------------------------------------------------
SEED    = 20260829;               % committed generator seed (distinct from the
                                  % toolbox's estimation seed, 12345, on purpose;
                                  % see SEED SELECTION above)
NPERGRP = 40;                     % participants per group (80 total)
GROUPS  = {'GroupA', 'GroupB'};
EVENTS  = {'inc_err', 'inc_cor'};
RPERS   = 0.75;                   % person-effect correlation between events

% cell means in microvolts, rows = groups, cols = events (order matches EVENTS)
MU      = [ -5.0, -1.5 ;          % GroupA
            -4.0, -1.5 ];         % GroupB

% person-effect SDs, same layout
SIGMA_P = [  2.0,  1.7 ;
             2.2,  1.7 ];

% residual (trial-noise) SDs at z = 0, same layout
SIGMA_E0 = [ 5.5,  5.0 ;
             7.5,  5.0 ];

% trial main-effect SD (sigma_i), one value for every cell
SIGMA_I = 1.0;

% dimension slopes: on the mean (per event) and on the log residual SD (all
% strata), per SD of the standardized covariate
B_MEAN  = [ -0.5,  0.0 ];
B_SIG   = 0.25;

% person-level log-residual SD (sigma_delta_p), one effect per participant
SIGMA_DELTA = 0.20;

% the covariate: an integer questionnaire total, one per participant
WORRY_MEAN  = 45;
WORRY_SD    = 12;
WORRY_RANGE = [16 80];

% trial-count ranges per event type: errors scarce, corrects plentiful
NTRIALS_ERR = [12 36];
NTRIALS_COR = [24 40];
ISERR       = [true false];       % which EVENTS entries are error events

% ---- design constants, returned for the ledger and the truth table ----------
P = struct('seed', SEED, 'npergrp', NPERGRP, 'groups', {GROUPS}, ...
    'events', {EVENTS}, 'rpers', RPERS, 'mu', MU, 'sigma_p', SIGMA_P, ...
    'sigma_e0', SIGMA_E0, 'sigma_i', SIGMA_I, 'b_mean', B_MEAN, ...
    'b_sig', B_SIG, 'sigma_delta', SIGMA_DELTA, 'worry_mean', WORRY_MEAN, ...
    'worry_sd', WORRY_SD, 'worry_range', WORRY_RANGE, ...
    'ntrials_err', NTRIALS_ERR, 'ntrials_cor', NTRIALS_COR, 'iserr', ISERR);

% ---- output location --------------------------------------------------------
if nargin < 1 || isempty(outfile)
    helpersdir = fileparts(mfilename('fullpath'));        % tests/helpers
    repo = fileparts(fileparts(helpersdir));              % repo root
    outfile = fullfile(repo, 'test_data', 'manual_ern_dynrel.csv');
end

% ---- simulate ---------------------------------------------------------------
if nargin < 2 || isempty(seed)
    seed = SEED;                              % the committed realization
end
P.seed = seed;
rng(seed, 'twister');

% The covariate is drawn for all 80 participants first, then standardized
% across the whole sample the way psyrat_zscore_subjectlevel does (sample SD),
% so that the generating slopes apply to the axis the toolbox draws.
nsub  = NPERGRP * numel(GROUPS);
worry = round(WORRY_MEAN + WORRY_SD * randn(nsub, 1));
worry = min(max(worry, WORRY_RANGE(1)), WORRY_RANGE(2));
z     = (worry - mean(worry)) / std(worry);

id = []; group = {}; event = {}; meas = []; wcol = [];
nextid = 401;                                 % ids 401..480
psub = 0;                                     % running participant index
for g = 1:numel(GROUPS)
    % Person-effect covariance for this group: a compound-symmetric correlation
    % (RPERS off the diagonal) scaled by the per-event SDs.
    D = diag(SIGMA_P(g, :));
    R = RPERS * ones(numel(EVENTS)) + (1 - RPERS) * eye(numel(EVENTS));
    C = chol(D * R * D, 'lower');

    % Trial main effects for this group: one draw per trial index per event,
    % over the largest count the event's range can produce, so every
    % participant at index t shares trial t's effect (persons crossed with
    % trials). Drawn before the participant loop to keep the stream simple.
    ieff = cell(1, numel(EVENTS));
    for e = 1:numel(EVENTS)
        if ISERR(e)
            maxn = NTRIALS_ERR(2);
        else
            maxn = NTRIALS_COR(2);
        end
        ieff{e} = SIGMA_I * randn(maxn, 1);
    end

    for p = 1:NPERGRP
        psub = psub + 1;
        pid = nextid; nextid = nextid + 1;
        peff  = C * randn(numel(EVENTS), 1);  % this person's 2 event effects
        delta = SIGMA_DELTA * randn;          % this person's log-residual effect

        for e = 1:numel(EVENTS)
            if ISERR(e)
                n = randi(NTRIALS_ERR);
            else
                n = randi(NTRIALS_COR);
            end
            sig_e  = exp(log(SIGMA_E0(g, e)) + B_SIG * z(psub) + delta);
            scores = MU(g, e) + B_MEAN(e) * z(psub) + peff(e) + ieff{e}(1:n) + ...
                sig_e * randn(n, 1);

            id    = [id;    repmat(pid, n, 1)];              %#ok<AGROW>
            group = [group; repmat(GROUPS(g), n, 1)];        %#ok<AGROW>
            event = [event; repmat(EVENTS(e), n, 1)];        %#ok<AGROW>
            meas  = [meas;  scores];                         %#ok<AGROW>
            wcol  = [wcol;  repmat(worry(psub), n, 1)];      %#ok<AGROW>
        end
    end
end

% Round the measurements so the committed CSV and a regenerated table agree
% exactly after writetable/readtable.
meas = round(meas, 6);

T = table(id, group, event, meas, wcol, ...
    'VariableNames', {'id', 'group', 'event', 'meas', 'worry'});

writetable(T, outfile);
fprintf('psyrat_make_manual_dynrel: wrote %d rows to %s\n', height(T), outfile);
end
