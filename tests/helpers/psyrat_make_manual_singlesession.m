function [T, P] = psyrat_make_manual_singlesession(outfile, seed)
%PSYRAT_MAKE_MANUAL_SINGLESESSION Generate the manual's single-session dataset.
%
%   T = PSYRAT_MAKE_MANUAL_SINGLESESSION() regenerates the tutorial dataset
%   test_data/manual_ern_singlesession.csv and returns it as a table.
%
%   T = PSYRAT_MAKE_MANUAL_SINGLESESSION(OUTFILE) writes to OUTFILE instead.
%
%   [T, P] = PSYRAT_MAKE_MANUAL_SINGLESESSION(OUTFILE, SEED) draws with a
%   different generator seed (used only by the seed-selection screen described
%   under SEED SELECTION) and also returns the design constants in the struct
%   P, so that the manual's truth tables and the drift checker read the
%   generating parameters from this file rather than from a copy.
%
% PURPOSE. This file is the follow-along dataset for the user manual's
% single-session tutorials (single-score reliability, two-event difference
% scores, and the four-event difference-of-differences example). It is fully
% simulated: every generating parameter is listed below, so the manual can
% print the true values next to the estimates the toolbox recovers. The
% committed CSV and this generator are welded together by
% tests/TestManualDatasets.m, which regenerates the table and requires it to
% match the committed file. Do not edit the CSV by hand; edit the parameters
% here, rerun, and update the manual's expected-value tables.
%
% DESIGN. A simulated flanker study measuring error-related negativity (ERN),
% scored as a mean amplitude in microvolts. One session; 80 participants in
% two independent groups (GroupA, GroupB; between-person), 40 per group, which
% keeps every group above the toolbox's below-20 advisory and gives the
% cell-level person variances enough participants to land near their
% generating values. Four events cross
% flanker congruency with response accuracy:
%
%   inc_err  incongruent-trial errors      inc_cor  incongruent-trial corrects
%   con_err  congruent-trial errors        con_cor  congruent-trial corrects
%
% Trial counts are drawn per participant per event (errors are scarcer than
% corrects, and counts differ across participants on purpose: the toolbox is
% built for unequal trial counts, and the tutorial says so).
%
% GENERATING MODEL (one-facet, persons crossed with trials within each
% group x event cell):
%
%   meas(p, e, t) = mu(g, e) + person(p, e) + trial(g, e, t) + noise(p, e, t)
%
%   * mu(g, e): cell means in microvolts. Error trials are more negative than
%     correct trials (the ERN), incongruent errors most negative, and GroupB's
%     error-trial means are attenuated by +1.0 microvolt:
%         GroupA: inc_err -5.0, inc_cor -1.5, con_err -4.0, con_cor -2.0
%         GroupB: inc_err -4.0, inc_cor -1.5, con_err -3.0, con_cor -2.0
%   * person(p, :): one draw per participant from a 4-variate normal with
%     mean zero, per-event SDs sigma_p (below), and a common correlation of
%     r = 0.75 between every pair of events. The correlation matters: the
%     two-event difference tutorial teaches that difference scores lose the
%     variance the two events share, and 0.75 is a realistic value for
%     correct/error ERN person effects.
%   * trial(g, e, :): a TRIAL MAIN EFFECT, SD sigma_i = 1.0, drawn once per
%     group x event cell over the cell's trial indices and shared by every
%     participant at the same index (persons crossed with trials). It is what
%     separates dependability from generalizability on real output: absolute
%     error carries sigma_i^2 / n and relative error does not, so the two
%     coefficients differ by roughly .005 to .01 at the observed trial counts,
%     enough to show at two decimals in the correct-trial cells (.78 vs .79 at
%     n = 32) though the error-trial cells can still round to the same value
%     (.74 vs .74 at n = 22). Participants with n trials take the first n
%     indices.
%   * noise: independent normal, SD sigma_e per group x event cell. Error
%     trials are noisier, and GroupB's error trials are noisier still, so the
%     tutorial's group comparison has something honest to show.
%
% TRUE VALUES (the manual's "truth tables" are computed from these):
%
%     cell            sigma_p   sigma_i   sigma_e     n per person
%     GroupA err        2.0       1.0       5.5       randi([12 36])
%     GroupA cor        1.7       1.0       5.0       randi([24 40])
%     GroupB err        2.2       1.0       7.5       randi([12 36])
%     GroupB cor        1.7       1.0       5.0       randi([24 40])
%
%   One-facet dependability at n trials is sigma_p^2 / (sigma_p^2 +
%   (sigma_i^2 + sigma_e^2) / n) and generalizability drops the sigma_i^2
%   term (Rocha et al., 2026, Table 2). For example, at n = 22 error trials,
%   GroupA inc_err: dependability 4.00 / (4.00 + 31.25/22) = 0.738 and
%   generalizability 4.00 / (4.00 + 30.25/22) = 0.744; GroupB inc_err:
%   4.84 / (4.84 + 57.25/22) = 0.650. At n = 32 correct trials, either
%   group's cor cells: dependability 2.89 / (2.89 + 26.00/32) = 0.781,
%   generalizability 2.89 / (2.89 + 25.00/32) = 0.787.
%   The .70 walkthrough stays inside the observed error-trial counts: at
%   dependability .70, GroupA's error cells need 31.25 / (4.00 * 0.3/0.7)
%   = 18.2, so 19 trials, and GroupB's need 57.25 / (4.84 * 0.3/0.7) = 27.6,
%   so 28 trials -- both within the drawn range of 12-36, with room for the
%   estimated GroupB cutoff to land above its true value and still fall inside
%   the data. At .80 the same cells need 32 and 48 trials: GroupA's sits at the
%   top of its range and GroupB's is beyond it, which is the case the manual
%   walks into on purpose.
%   The two-event difference dependability follows the same logic with
%   sigma^2_diff = sigma_p1^2 + sigma_p2^2 - 2 r sigma_p1 sigma_p2 in the
%   numerator; with r = 0.75 the difference is far less dependable than
%   either constituent, which is the point of that tutorial.
%
% REPRODUCIBILITY. rng(SEED, 'twister') below fixes everything. The
% measurement column is rounded to 6 decimals before writing so that the
% committed CSV and a regenerated table agree exactly after a round trip.
%
% SEED SELECTION (disclosed). The committed seed is the first in the fixed
% sequence 20260829, 20260830, ... whose realization met the manual's
% pre-registered realization criteria (MANUAL_PRODUCTION.md, section 5:
% every tutorial fit converged with zero divergent transitions; at .70 every
% cell's trial cutoff lies inside the observed counts with GroupA retaining
% at least 12 of 40 and GroupB at least 1; at .80 GroupB's error cells
% extrapolate; the two-event difference cutoff extrapolates in both groups).
% Those criteria were written before the first regenerated fit, and the
% record there names every seed that was fitted. The committed seed,
% 20260836, is position 8: positions 1-7 failed a
% method-of-moments pre-screen of the same criteria, and position 20 was
% fitted under a miscalibrated pre-screen and failed the .80 retention
% criterion. This is a pedagogical selection so that the tutorial walks
% through a realization that shows each lesson; it is not a validation of
% the estimator, and the manual says so.
%
% This generator is standalone base MATLAB: no toolboxes, no test framework.

% ---- fixed design parameters ------------------------------------------------
SEED    = 20260836;               % committed generator seed (distinct from the
                                  % toolbox's estimation seed, 12345, on purpose;
                                  % see SEED SELECTION above)
NPERGRP = 40;                     % participants per group (80 total)
GROUPS  = {'GroupA', 'GroupB'};
EVENTS  = {'inc_err', 'inc_cor', 'con_err', 'con_cor'};
RPERS   = 0.75;                   % person-effect correlation between events

% cell means, rows = groups, cols = events (order matches EVENTS)
MU      = [ -5.0, -1.5, -4.0, -2.0 ;     % GroupA
            -4.0, -1.5, -3.0, -2.0 ];    % GroupB

% person-effect SDs, same layout
SIGMA_P = [  2.0,  1.7,  2.0,  1.7 ;
             2.2,  1.7,  2.2,  1.7 ];

% trial-noise SDs, same layout
SIGMA_E = [  5.5,  5.0,  5.5,  5.0 ;
             7.5,  5.0,  7.5,  5.0 ];

% trial main-effect SD (sigma_i), one value for every cell: large enough to
% separate dependability from generalizability on screen, small enough that
% it never dominates the person variance
SIGMA_I = 1.0;

% trial-count ranges per event type: errors scarce, corrects plentiful
NTRIALS_ERR = [12 36];
NTRIALS_COR = [24 40];
ISERR       = [true false true false];   % which EVENTS entries are error events

% ---- design constants, returned for the manual's truth tables ---------------
P = struct('seed', SEED, 'npergrp', NPERGRP, 'groups', {GROUPS}, ...
    'events', {EVENTS}, 'rpers', RPERS, 'mu', MU, 'sigma_p', SIGMA_P, ...
    'sigma_e', SIGMA_E, 'sigma_i', SIGMA_I, 'ntrials_err', NTRIALS_ERR, ...
    'ntrials_cor', NTRIALS_COR, 'iserr', ISERR);

% ---- output location --------------------------------------------------------
if nargin < 1 || isempty(outfile)
    helpersdir = fileparts(mfilename('fullpath'));        % tests/helpers
    repo = fileparts(fileparts(helpersdir));              % repo root
    outfile = fullfile(repo, 'test_data', 'manual_ern_singlesession.csv');
end

% ---- simulate ---------------------------------------------------------------
if nargin < 2 || isempty(seed)
    seed = SEED;                              % the committed realization
end
P.seed = seed;
rng(seed, 'twister');

id = []; group = {}; event = {}; meas = [];
nextid = 101;                                 % ids 101..180
for g = 1:numel(GROUPS)
    % Build the person-effect covariance for this group once: a compound-
    % symmetric correlation (RPERS off the diagonal) scaled by the per-event
    % SDs. chol() gives the transform from iid normals to correlated effects.
    D = diag(SIGMA_P(g, :));
    R = RPERS * ones(numel(EVENTS)) + (1 - RPERS) * eye(numel(EVENTS));
    C = chol(D * R * D, 'lower');

    % Trial main effects for this group: one draw per trial index per event,
    % over the largest count the event's range can produce, so every
    % participant at index t shares trial t's effect (persons crossed with
    % trials). Drawn before the participant loop to keep the stream layout
    % simple.
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
        pid = nextid; nextid = nextid + 1;
        peff = C * randn(numel(EVENTS), 1);   % this person's 4 event effects

        for e = 1:numel(EVENTS)
            if ISERR(e)
                n = randi(NTRIALS_ERR);
            else
                n = randi(NTRIALS_COR);
            end
            scores = MU(g, e) + peff(e) + ieff{e}(1:n) + ...
                SIGMA_E(g, e) * randn(n, 1);

            id    = [id;    repmat(pid, n, 1)];              %#ok<AGROW>
            group = [group; repmat(GROUPS(g), n, 1)];        %#ok<AGROW>
            event = [event; repmat(EVENTS(e), n, 1)];        %#ok<AGROW>
            meas  = [meas;  scores];                         %#ok<AGROW>
        end
    end
end

% Round the measurements so the committed CSV and a regenerated table agree
% exactly after writetable/readtable (writetable prints the shortest decimal
% that round-trips, so rounding first removes any representation ambiguity).
meas = round(meas, 6);

T = table(id, group, event, meas, ...
    'VariableNames', {'id', 'group', 'event', 'meas'});

writetable(T, outfile);
fprintf('psyrat_make_manual_singlesession: wrote %d rows to %s\n', ...
    height(T), outfile);
end
