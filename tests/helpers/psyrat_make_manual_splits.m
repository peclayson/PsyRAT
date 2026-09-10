function [T, P] = psyrat_make_manual_splits(outfile, seed)
%PSYRAT_MAKE_MANUAL_SPLITS Generate the manual's reliability-from-splits dataset.
%
%   T = PSYRAT_MAKE_MANUAL_SPLITS() regenerates the tutorial dataset
%   test_data/manual_ern_splits.csv and returns it as a table.
%
%   T = PSYRAT_MAKE_MANUAL_SPLITS(OUTFILE) writes to OUTFILE instead.
%
%   [T, P] = PSYRAT_MAKE_MANUAL_SPLITS(OUTFILE, SEED) draws with another
%   generator seed (seed-selection screen only) and returns the design
%   constants in P for the manual's truth tables and the drift checker.
%
% PURPOSE. Follow-along data for the manual's splits tutorial (reliability
% when only split means, not single trials, are available). Fully simulated;
% welded to the committed CSV by tests/TestManualDatasets.m.
%
% DESIGN. 36 participants (ids 301-336). Each contributes 16 split means:
% each split is the average of n_i single-trial ERN scores, and n_i is
% recorded in the REQUIRED weight column (a positive integer; the splits
% input contract). The n_i are deliberately unequal, randi([4 40]), so the
% file exercises the Nonparallel path; the tutorial notes that equal n_i
% would instead qualify for the Parallel declaration.
%
% Why 16 splits and a 4-40 weight range (design note). The observed-design
% splits model (case 23) separates the person x split SD from the per-item
% residual only through the weight pattern, and it runs at CmdStan's default
% target acceptance. With 6 splits and weights 8-30 the sampler reported
% hundreds of divergent transitions on clean simulated data whatever the
% split components were set to (2026-09-03 design experiments, recorded in
% MANUAL_PRODUCTION.md section 5); 16 splits with weights 4-40 was the least
% change that sampled with none.
%
% GENERATING MODEL (matches the toolbox's observed-design splits model, case
% 23: persons crossed with splits, and the residual of a split MEAN shrinks
% with the number of items averaged into it):
%
%   meas(p, s) = mu + person(p) + split(s) + ps(p, s) + e(p, s),
%                e ~ Normal(0, sigma_item^2 / n_i)
%
%     mu          -3.5   microvolts (grand mean)
%     sigma_p      2.0   person SD
%     sigma_s      0.5   split main-effect SD (one draw per split index,
%                        shared by every participant at that index)
%     sigma_ps     0.7   person x split interaction SD
%     sigma_item   5.0   single-item (single-trial) SD
%
% The two split-level components are small against the person SD and the
% split-mean noise (about 5.0 / sqrt(22) = 1.07 at the mean weight), but they
% are not zero: the toolbox estimates both, and a generator that sets them to
% exactly zero is no better behaved (the 2026-09-02 six-split realization
% produced 414 divergent transitions that way). The item-level residual
% sigma_item is the ground truth the split-level estimation targets;
% dependability computed from these components is exact, while the
% generalizability coefficient reported for splits data is a lower bound
% (the tutorial explains why; with no item main effect in the generator the
% bound is tight). At the observed design (n_s splits, harmonic-mean n-bar_i
% items per split):
%
%   D = sigma_p^2 / (sigma_p^2 + (sigma_ps^2 + sigma_s^2) / n_s
%                    + sigma_item^2 / (n_s * n-bar_i))
%   G = the same without the sigma_s^2 / n_s term
%
% (documentation/splits_observed_design_formulas.md). The manual's truth
% table is computed from these numbers and the committed file's weights.
%
% REPRODUCIBILITY. rng(SEED, 'twister'); measurements rounded to 6 decimals
% before writing (exact CSV round trip). Standalone base MATLAB.
%
% SEED SELECTION (disclosed). The committed seed is the first in the fixed
% sequence 20260829, 20260830, ... whose realization met the manual's
% pre-registered realization criteria (MANUAL_PRODUCTION.md, section 5:
% converged with zero divergent transitions and dependability within .03 of
% the truth above). The criteria were written before the first regenerated
% fit. The committed seed, 20260829, is position 1 in the sequence: the
% adopted 16-split design sampled with zero divergences at the first seed
% tried (MANUAL_PRODUCTION.md, section 5).

% ---- fixed design parameters ------------------------------------------------
SEED       = 20260829;            % committed generator seed (see SEED SELECTION)
NSUB       = 36;                  % participants, ids 301..336
NSPLITS    = 16;                  % split means per participant (see design note)
NITEMS     = [4 40];              % items per split: randi range (unequal on purpose)
MU         = -3.5;
SIGMA_P    = 2.0;
SIGMA_S    = 0.5;                 % split main effect (shared across persons)
SIGMA_PS   = 0.7;                 % person x split interaction
SIGMA_ITEM = 5.0;

% ---- design constants, returned for the manual's truth tables ---------------
P = struct('seed', SEED, 'nsub', NSUB, 'nsplits', NSPLITS, 'nitems', NITEMS, ...
    'mu', MU, 'sigma_p', SIGMA_P, 'sigma_s', SIGMA_S, 'sigma_ps', SIGMA_PS, ...
    'sigma_item', SIGMA_ITEM);

% ---- output location --------------------------------------------------------
if nargin < 1 || isempty(outfile)
    helpersdir = fileparts(mfilename('fullpath'));        % tests/helpers
    repo = fileparts(fileparts(helpersdir));              % repo root
    outfile = fullfile(repo, 'test_data', 'manual_ern_splits.csv');
end

% ---- simulate ---------------------------------------------------------------
if nargin < 2 || isempty(seed)
    seed = SEED;                               % the committed realization
end
P.seed = seed;
rng(seed, 'twister');

% Split main effects: one draw per split index, shared by every participant
% at that index (persons crossed with splits), drawn before the participant
% loop exactly as the single-session generator draws its trial main effect.
seff = SIGMA_S * randn(NSPLITS, 1);

id = []; meas = []; weight = [];
for p = 1:NSUB
    pid = 300 + p;
    peff = SIGMA_P * randn;                    % this person's true effect

    n_i = randi(NITEMS, NSPLITS, 1);           % items behind each split mean
    pseff = SIGMA_PS * randn(NSPLITS, 1);      % this person's split interactions
    % A split mean of n_i items has residual SD sigma_item / sqrt(n_i).
    splitmeans = MU + peff + seff + pseff + ...
        (SIGMA_ITEM ./ sqrt(n_i)) .* randn(NSPLITS, 1);

    id     = [id;     repmat(pid, NSPLITS, 1)];    %#ok<AGROW>
    meas   = [meas;   splitmeans];                 %#ok<AGROW>
    weight = [weight; n_i];                        %#ok<AGROW>
end

meas = round(meas, 6);   % exact CSV round trip (see single-session generator)

T = table(id, meas, weight, 'VariableNames', {'id', 'meas', 'weight'});

writetable(T, outfile);
fprintf('psyrat_make_manual_splits: wrote %d rows to %s\n', height(T), outfile);
end
