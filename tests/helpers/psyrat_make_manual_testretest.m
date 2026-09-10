function [T, P] = psyrat_make_manual_testretest(outfile, seed)
%PSYRAT_MAKE_MANUAL_TESTRETEST Generate the manual's test-retest dataset.
%
%   T = PSYRAT_MAKE_MANUAL_TESTRETEST() regenerates the tutorial dataset
%   test_data/manual_ern_testretest.csv and returns it as a table.
%
%   T = PSYRAT_MAKE_MANUAL_TESTRETEST(OUTFILE) writes to OUTFILE instead.
%
%   [T, P] = PSYRAT_MAKE_MANUAL_TESTRETEST(OUTFILE, SEED) draws with another
%   generator seed (seed-selection screen only) and returns the design
%   constants in P for the manual's truth tables and the drift checker.
%
% PURPOSE. Follow-along data for the manual's test-retest tutorial. Fully
% simulated with the generating parameters listed below, so the manual can
% show the toolbox recovering known variance components. Welded to the
% committed CSV by tests/TestManualDatasets.m: edit parameters here, rerun,
% and update the manual's expected values; never hand-edit the CSV.
%
% DESIGN. A simulated study in which 30 participants (ids 201-230, no group
% column) complete the same flanker task in two sessions two weeks apart
% (time = t1, t2). Two events (err, cor). ERN scored as mean amplitude in
% microvolts. Trial counts are balanced by design here (24 error and 40
% correct trials per participant per occasion): the crossed generating model
% below needs a common trial index, and the unequal-counts lesson already
% belongs to the single-session tutorial.
%
% GENERATING MODEL. Within each event, the fully crossed
% person x occasion x trial decomposition used throughout the toolbox's
% test-retest validation (see localSimTrt in tests/TestNativeEngine.m, the
% template this generator mirrors):
%
%   meas = mu + p(id) + o(occ) + t(trl) + pt(id,trl) + po(id,occ)
%             + ot(occ,trl) + noise
%
% with independent normal effects of SD (microvolts):
%
%     component                 err     cor
%     sigma_p    (person)       3.2     2.8
%     sigma_o    (occasion)     0.6     0.6
%     sigma_t    (trial)        1.0     1.0
%     sigma_pt   (person x trial)   1.0     0.8
%     sigma_po   (person x occ)     1.7     1.4
%     sigma_ot   (occ x trial)      0.5     0.5
%     sigma_e    (residual)     5.5     5.0
%     mu         (mean)        -4.5    -1.8
%
% Why these values (design note). With two occasions the occasion SD is
% barely identified and its posterior is prior-dominated (the toolbox's
% half-Cauchy(0, 5) prior on two occasion effects lands well above the truth
% whatever the truth is: 2.2 to 2.6 microvolts on the committed fit against a
% generating 0.6; 1.8 was the design-time projection). Stability and
% equivalence-and-stability DEPENDABILITY count
% sigma_o^2 as error, so a design with the ERN-typical person SD of 2.0
% produces estimates of those two coefficients far below their true values
% and below the .50 threshold the tutorial walks through. The person SDs
% are set so that both coefficients clear .50 in the estimate even with the
% inflated occasion term, and the person x occasion SDs so that, at the true
% values, neither coefficient reaches .80 at any trial count (the chapter's
% lesson that more trials cannot repair day-to-day instability). sigma_t and
% sigma_ot are nonzero so that a trial and an occasion x trial main effect
% exist to estimate; on the committed fit the occasion x trial SD posterior
% still reaches near zero, the two-level occasion SD spans 0.1 to 7, and the
% sampler reports divergent transitions (909 at the tutorial's settings) that
% no retune of these components removed (MANUAL_PRODUCTION.md section 5).
% Equivalence dependability and generalizability agree to two decimals at
% the observed counts and separate only at the trial cutoffs and in the SEMs.
%
% Person effects are correlated r = 0.75 between the two events (a shared /
% event-specific decomposition preserving each event's marginal SD), matching
% the single-session dataset's convention.
%
% TRUE COEFFICIENTS. Equivalence and stability coefficients follow from these
% components via the toolbox's test-retest formulas (Rocha et al., 2026,
% Table 3, with the toolbox's n'_o = 1 default). The manual's tutorial
% chapter derives and prints those true values next to the recovered
% estimates; this header records only the generating components, which are
% the ground truth everything else is computed from.
%
% REPRODUCIBILITY. rng(SEED, 'twister'); measurements rounded to 6 decimals
% before writing (exact CSV round trip). Standalone base MATLAB.
%
% SEED SELECTION (disclosed). The committed seed, 20260829, is position 1 in
% the fixed sequence 20260829, 20260830, ... . Its realization passed the
% cutoff criterion of the manual's pre-registered realization criteria
% (MANUAL_PRODUCTION.md, section 5, R5: the equivalence .70 cutoffs and the
% stability and equivalence-and-stability .50 cutoffs all lie inside the
% design's trial counts, all 30 participants retained) but NOT the
% zero-divergence criterion (R1: 909 divergent transitions at the tutorial's
% settings). The owner committed it as measured on 2026-09-03; the tutorial
% chapter states the count and explains its cause. Pedagogical selection, not
% estimator validation.

% ---- fixed design parameters ------------------------------------------------
SEED   = 20260829;                 % committed generator seed (see SEED SELECTION)
NSUB   = 30;                       % participants, ids 201..230
NOCC   = 2;                        % occasions t1, t2
OCCLAB = {'t1', 't2'};
EVENTS = {'err', 'cor'};
NTRL   = [24 40];                  % trials per participant per occasion, by event
RPERS  = 0.75;                     % person-effect correlation between events

%                 err    cor
MU      = [      -4.5,  -1.8 ];
SIGMA_P = [       3.2,   2.8 ];
SIGMA_O = [       0.6,   0.6 ];
SIGMA_T = [       1.0,   1.0 ];
SIGMA_PT = [      1.0,   0.8 ];
SIGMA_PO = [      1.7,   1.4 ];
SIGMA_OT = [      0.5,   0.5 ];
SIGMA_E = [       5.5,   5.0 ];

% ---- design constants, returned for the manual's truth tables ---------------
P = struct('seed', SEED, 'nsub', NSUB, 'nocc', NOCC, 'occlab', {OCCLAB}, ...
    'events', {EVENTS}, 'ntrl', NTRL, 'rpers', RPERS, 'mu', MU, ...
    'sigma_p', SIGMA_P, 'sigma_o', SIGMA_O, 'sigma_t', SIGMA_T, ...
    'sigma_pt', SIGMA_PT, 'sigma_po', SIGMA_PO, 'sigma_ot', SIGMA_OT, ...
    'sigma_e', SIGMA_E);

% ---- output location --------------------------------------------------------
if nargin < 1 || isempty(outfile)
    helpersdir = fileparts(mfilename('fullpath'));        % tests/helpers
    repo = fileparts(fileparts(helpersdir));              % repo root
    outfile = fullfile(repo, 'test_data', 'manual_ern_testretest.csv');
end

% ---- simulate ---------------------------------------------------------------
if nargin < 2 || isempty(seed)
    seed = SEED;                   % the committed realization
end
P.seed = seed;
rng(seed, 'twister');

% Correlated person effects: p_event = sqrt(r)*shared + sqrt(1-r)*specific,
% with shared and specific both unit-normal, then scaled by the event's
% sigma_p. Marginal SD per event is exactly sigma_p; the between-event
% correlation is exactly RPERS.
shared = randn(NSUB, 1);
pers = zeros(NSUB, numel(EVENTS));
for e = 1:numel(EVENTS)
    specific = randn(NSUB, 1);
    pers(:, e) = SIGMA_P(e) * (sqrt(RPERS) * shared + sqrt(1 - RPERS) * specific);
end

id = []; time = {}; event = {}; meas = [];
for e = 1:numel(EVENTS)
    n = NTRL(e);

    % Crossed effects for this event stratum (mirrors localSimTrt): every
    % participant sees trial indices 1..n on both occasions, so occasion,
    % trial, and all two-way interactions are identifiable.
    occeff = SIGMA_O(e) * randn(NOCC, 1);
    trleff = SIGMA_T(e) * randn(n, 1);
    pteff  = SIGMA_PT(e) * randn(NSUB, n);
    poeff  = SIGMA_PO(e) * randn(NSUB, NOCC);
    oteff  = SIGMA_OT(e) * randn(NOCC, n);

    for p = 1:NSUB
        pid = 200 + p;
        for o = 1:NOCC
            noise = SIGMA_E(e) * randn(n, 1);
            scores = MU(e) + pers(p, e) + occeff(o) + trleff ...
                + pteff(p, :)' + poeff(p, o) + oteff(o, :)' + noise;

            id    = [id;    repmat(pid, n, 1)];               %#ok<AGROW>
            time  = [time;  repmat(OCCLAB(o), n, 1)];         %#ok<AGROW>
            event = [event; repmat(EVENTS(e), n, 1)];         %#ok<AGROW>
            meas  = [meas;  scores];                          %#ok<AGROW>
        end
    end
end

meas = round(meas, 6);   % exact CSV round trip (see single-session generator)

T = table(id, time, event, meas, ...
    'VariableNames', {'id', 'time', 'event', 'meas'});

writetable(T, outfile);
fprintf('psyrat_make_manual_testretest: wrote %d rows to %s\n', ...
    height(T), outfile);
end
