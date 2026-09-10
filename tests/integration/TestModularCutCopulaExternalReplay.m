classdef TestModularCutCopulaExternalReplay < PsyRATTestBase
    % EXTERNAL replay of the owner's completed modular cut-copula reference
    % fits (analysis 19, one-facet, person-specific dynamic concurrent
    % difference reliability) through the MATLAB stage-2/stage-3 machinery.
    %
    % WHAT THIS IS. The reference bundle
    %   ~/Desktop/stan_chisquare/person_specific_dynamic_concurrent_modular
    % contains two COMPLETED 151-person one-facet runs (observed data, 4 chains
    % x 10,000 draws, completion manifests present):
    %   theta_minus_alpha_fla_cor_trial_one_dimension_...   (K = 1)
    %   theta_minus_alpha_fla_cor_trial_two_dimensions_...  (K = 3, with the
    %                                                        Dim1*Dim2 term)
    % This test pushes each run's own posterior through PsyRAT's committed
    % stage 2 (PSYRAT_GAMMA_CUT_COPULA_STAGE) and stage 3
    % (PSYRAT_SSREL_DIFFDYNREL_RESCOR_GAMMA_LS) and compares against the
    % reference's OWN outputs on disk. Both sides therefore consume the SAME
    % posterior draws, so any disagreement is arithmetic, not Monte Carlo -
    % the same design the 8c external validation used
    % (documentation/gamma_locationscale_rollout.md, increment 8c).
    %
    % It validates the ONE-FACET path only. The bundle's two-facet runs were
    % prepared but never sampled (zero draw rows, no completion manifest), so
    % there is nothing to replay for analysis 20; its arithmetic is pinned by
    % the frozen fixtures in tests/baselines/modular_cut_copula/ instead.
    %
    % THE THREE RECORDED TRAPS, and where each is handled here:
    %  (1) mu_offset / log_sigma_offset are DATA-DERIVED centering constants;
    %      the fitted intercepts are deviations from them. Stage 2 consumes
    %      them inside PSYRAT_GAMMA_MARGIN_PREDICTORS via standata; stage 3
    %      receives ABSOLUTE intercepts, so the offsets are added back when the
    %      b / b_sigma / mu_ss / er_ss arrays are built (localStageThreeArgs).
    %  (2) Posterior means are insufficient (Jensen): every quantity is
    %      computed per draw and only then summarized, exactly as the
    %      reference's compute_person_draw does.
    %  (3) The draw subset must be ALIGNED: the copula correlation exists only
    %      for the 500 retained stage-1 draws. testDrawSelection pins the
    %      retained set against the reference's copula_draw_selection.csv, and
    %      the loader reads exactly those chain rows.
    %
    % WHAT IS AND IS NOT COMPARED IN STAGE 2. Given a marginal draw, everything
    % on the rho grid is deterministic: the conditional mean/sd/median/
    % quantiles/mode, the boundary mass, and every PIT diagnostic. Those are
    % compared PER DRAW. rho_e itself is a single inverse-CDF draw whose
    % uniform comes from the platform RNG (R's Mersenne stream seeded
    % seed+700001 there, MATLAB's mt19937ar here), so it CANNOT match per draw
    % and is compared distributionally only. Stage 3 then consumes the
    % REFERENCE'S OWN rho_e draws (modular_copula_draws.csv), which keeps the
    % per-participant comparison purely arithmetic.
    %
    % SUMMARY CONVENTIONS. PsyRAT's point estimates are means over draws, which
    % match the reference summary's `mean` column sharply. Interval endpoints
    % use MATLAB's quantile method against R's default type 7, so the ll/ul
    % comparisons carry a small method difference at 500 draws and get a looser
    % ceiling; the mean comparisons are the exact ones.
    %
    % BUNDLE-GATED: skips (Incomplete, not Failed) when the bundle is absent.
    % Override the bundle location with the MPSDC_REFERENCE_DIR environment
    % variable, mirroring make_modular_cut_copula_fixtures.R.
    %
    % See also PSYRAT_GAMMA_CUT_COPULA_STAGE, PSYRAT_GAMMA_MARGIN_PREDICTORS,
    % PSYRAT_SSREL_DIFFDYNREL_RESCOR_GAMMA_LS.

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

    properties (TestParameter)
        % The two completed one-facet reference runs. 'one_dimension' fits a
        % single standardized predictor (K = 1); 'two_dimensions' fits two plus
        % their product (K = 3), which exercises the event-crossed slope layout.
        run = {'one_dimension', 'two_dimensions'}
    end

    methods (Test)

        function testProductionSelectorImplementsTheDocumentedRule(testCase)
            % Trap (3), the half that needs NO bundle: the PRODUCTION selector,
            % executed on a 40,000-draw synthetic (the reference configuration)
            % with trivial data, must return exactly the documented rule - even
            % spacing with round-half-to-even. This is the assertion that
            % actually pins the code, and it is deliberately UNGATED: before
            % the pre-merge correction it sat behind the bundle assume below
            % and executed only where the owner's reference bundle exists.
            nd_total = 40000;
            expected = unique(localRoundHalfToEven(linspace(1, nd_total, 500)))';

            % All-zero draws over a two-observation dataset: the margins are
            % valid (mu = sigma = 1), so the stage runs end to end while
            % costing nothing, and out.sel is the production selection.
            synth = struct( ...
                'alpha_mu', zeros(nd_total, 2), 'alpha_sigma', zeros(nd_total, 2), ...
                'beta_mu', zeros(nd_total, 2, 1), 'beta_sigma', zeros(nd_total, 2, 1), ...
                'u_p', zeros(nd_total, 1, 4), 'u_i', zeros(nd_total, 1, 2));
            sdata = struct('N', 2, 'K', 1, 'Z', zeros(2, 1), ...
                'y_1', [1; 1], 'y_2', [1; 1], 'p_id', [1; 1], 'i_id', [1; 1], ...
                'mu_offset', [0 0], 'log_sigma_offset', [0 0]);
            out = psyrat_gamma_cut_copula_stage(synth, sdata, 12345);
            testCase.verifyEqual(out.sel, expected, ...
                ['psyrat_gamma_cut_copula_stage''s deterministic draw ', ...
                'selection does not implement even spacing with ', ...
                'round-half-to-even at the reference draw count.']);
        end

        function testDrawSelectionMatchesReference(testCase, run)
            % Trap (3), the bundle half: the reference's own retained-draw CSV
            % must match the same documented rule, so the loader's alignment
            % assumption holds. The production selector itself is pinned by the
            % ungated test above.
            info = localRunInfo(run);
            testCase.assumeTrue(~isempty(info), localSkipMessage());

            sel_csv = localReadSelection(info);
            nd_total = max(sel_csv.draw);
            testCase.verifyEqual(nd_total, 40000, ...
                'The reference run should hold 4 chains x 10,000 draws.');

            % The documented rule, replicated independently of the production
            % code so the two implementations check each other.
            expected = unique(localRoundHalfToEven(linspace(1, nd_total, 500)))';
            testCase.verifyEqual(sel_csv.draw, expected, ...
                ['The reference retained-draw set does not match even ', ...
                'spacing with round-half-to-even; the loader''s alignment ', ...
                'assumption would be wrong.']);
        end

        function testStageTwoDeterministicPosteriorMatchesReference(testCase, run)
            % Stage 2 on the reference's own retained marginal draws. Every
            % deterministic per-draw quantity must reproduce the reference's
            % modular_copula_draws.csv; rho_e (the one random quantity) is
            % compared distributionally against its summary.
            info = localRunInfo(run);
            testCase.assumeTrue(~isempty(info), localSkipMessage());
            data = localLoadRun(run);

            t0 = tic;
            out = psyrat_gamma_cut_copula_stage(data.draws, data.standata, 12345);
            fprintf('[%s] stage 2 over %d retained draws: %.1f s\n', ...
                run, numel(out.sel), toc(t0));

            % Passing exactly the 500 retained draws means the stage retains
            % all of them; anything else is a loader defect.
            testCase.assertEqual(numel(out.sel), height(data.copula), ...
                'The stage did not retain every supplied draw.');

            % PsyRAT field name -> reference CSV column. Components 1/2 are
            % theta/alpha in the reference's naming.
            pairs = { ...
                'rho_mean',           'rho_conditional_mean'; ...
                'rho_sd',             'rho_conditional_sd'; ...
                'rho_median',         'rho_conditional_median'; ...
                'rho_q025',           'rho_conditional_q025'; ...
                'rho_q975',           'rho_conditional_q975'; ...
                'rho_mode',           'rho_conditional_mode'; ...
                'rho_boundary_mass',  'rho_boundary_mass'; ...
                'pit_z_1_mean',       'pit_z_theta_mean'; ...
                'pit_z_1_sd',         'pit_z_theta_sd'; ...
                'pit_z_2_mean',       'pit_z_alpha_mean'; ...
                'pit_z_2_sd',         'pit_z_alpha_sd'; ...
                'pit_normal_score_correlation', 'pit_normal_score_correlation'; ...
                'pit_1_abs_gt_1_96',  'pit_theta_abs_gt_1_96'; ...
                'pit_2_abs_gt_1_96',  'pit_alpha_abs_gt_1_96'; ...
                'pit_joint_upper_1_96', 'pit_joint_upper_1_96'; ...
                'pit_joint_lower_1_96', 'pit_joint_lower_1_96'; ...
                'pit_1_lower_clipped', 'pit_theta_lower_clipped'; ...
                'pit_1_upper_clipped', 'pit_theta_upper_clipped'; ...
                'pit_2_lower_clipped', 'pit_alpha_lower_clipped'; ...
                'pit_2_upper_clipped', 'pit_alpha_upper_clipped'; ...
                'pit_any_clip_fraction', 'pit_any_clip_fraction'};

            % Report first, assert after, so a failed ceiling still leaves the
            % full measurement on record (the measure-then-pin protocol).
            fprintf('[%s] stage-2 per-draw agreement (max |MATLAB - reference|):\n', run);
            worst = zeros(size(pairs, 1), 1);
            for k = 1:size(pairs, 1)
                got = out.(pairs{k, 1});
                want = data.copula.(pairs{k, 2});
                %max() omits NaN by default, so without this a NaN on either
                %side would vanish from the pooled ceiling and the comparison
                %would silently pass on the remaining quantities.
                testCase.verifyTrue(all(isfinite(got(:))) && all(isfinite(want(:))), ...
                    sprintf('%s: non-finite value in the stage-2 comparison.', pairs{k, 1}));
                worst(k) = max(abs(got - want));
                fprintf('  %-30s %.3e\n', pairs{k, 1}, worst(k));
            end

            % MEASURED 2026-08-10, both runs: grid-valued outputs (median,
            % quantiles, mode) agreed to <= 2e-16 - the grid indices matched
            % exactly, with no one-step drift - and the worst quantity overall
            % was pit_z_2_sd at 3.2e-9, consistent with the two platforms'
            % incomplete-gamma and inverse-normal tails differing in their
            % last digits (the same rationale the frozen fixtures record for
            % their ~1e-7 normal-score tolerance). The ceiling sits ~30x above
            % that measured worst; it is a regression pin on the measured
            % agreement, not a claim about the smallest detectable defect.
            testCase.verifyLessThan(max(worst), 1e-7, ...
                'A stage-2 deterministic quantity exceeds the measured-agreement ceiling.');

            % The two flags are deterministic functions of the quantities
            % above, far from their thresholds in these runs; exact equality.
            testCase.verifyEqual(out.invalidity_flag, ...
                localLogical(data.copula.copula_stage_invalidity_flag), ...
                'invalidity_flag disagrees with the reference.');
            testCase.verifyEqual(out.calibration_flag, ...
                localLogical(data.copula.pit_calibration_review_flag), ...
                'calibration_flag disagrees with the reference.');

            % rho_e: platform RNG differs, so only the distribution can agree.
            % Both samples share the conditional-posterior mixture, so their
            % means differ by two independent draws of the inverse-CDF noise;
            % six standard errors is generous without being vacuous.
            se = std(out.rho_e) / sqrt(numel(out.rho_e));
            ref_mean = data.copula_summary.mean(strcmp(data.copula_summary.quantity, 'rho_e'));
            fprintf('[%s] rho_e mean: MATLAB %.6f vs reference %.6f (se %.2e)\n', ...
                run, mean(out.rho_e), ref_mean, se);
            testCase.verifyLessThan(abs(mean(out.rho_e) - ref_mean), 6 * se, ...
                'The rho_e draw distribution is off beyond RNG noise.');
        end

        function testStageThreeMatchesReferencePerParticipant(testCase, run)
            % Stage 3 on the reference's own retained draws AND its own rho_e
            % draws, so both sides run identical inputs; comparison against
            % person_reliability_summary.csv at both D-study designs.
            info = localRunInfo(run);
            testCase.assumeTrue(~isempty(info), localSkipMessage());
            data = localLoadRun(run);

            args = localStageThreeArgs(data);

            % The reference evaluates two designs per participant: 'actual'
            % (the participant's own trial count) and 'standardized' (a common
            % n_i-prime, 341 in both runs). Same machinery, different idtable.
            designs = {'actual', 'standardized'};
            for di = 1:2
                design = designs{di};
                idtable = localIdTable(data, design);

                t0 = tic;
                [gen, dep, diagnostic] = psyrat_ssrel_diffdynrel_rescor_gamma_ls( ...
                    args{:}, 'idtable', idtable, 'ndim', data.ndim, 'CI', 0.95);
                fprintf('[%s/%s] stage 3 over %d participants x %d draws: %.1f s\n', ...
                    run, design, height(idtable), size(data.draws.alpha_mu, 1), toc(t0));

                % PsyRAT column <- source table, against reference quantity;
                % 'rel'/'abs' picks the difference scale (bounded coefficients
                % compare absolutely; score-scale quantities relatively).
                checks = { ...
                    gen.rel_pt,             'G_delta',                 'abs'; ...
                    dep.rel_pt,             'D_delta',                 'abs'; ...
                    gen.icc_pt,             'ICC_relative_delta',      'abs'; ...
                    dep.icc_pt,             'ICC_absolute_delta',      'abs'; ...
                    gen.sem_pt,             'SEM_relative_delta',      'rel'; ...
                    dep.sem_pt,             'SEM_absolute_delta',      'rel'; ...
                    gen.bp_var,             'var_delta_p',             'rel'; ...
                    gen.ss_errvar,          'relative_error_delta',    'rel'; ...
                    diagnostic.res_cov_pt,  'residual_covariance_person', 'rel'; ...
                    diagnostic.res_cor_pt,  'residual_correlation_person', 'abs'; ...
                    diagnostic.rho_e_pt,    'rho_e_copula',            'abs'};

                fprintf('[%s/%s] per-participant agreement against the reference mean:\n', run, design);
                worst_abs = 0; worst_rel = 0;
                for k = 1:size(checks, 1)
                    got = checks{k, 1};
                    want = localReferenceColumn(data, checks{k, 2}, design, 'mean');
                    %max() omits NaN, so a NaN comparand (e.g. a NaN SEM from a
                    %negative error-variance draw) would silently drop out of
                    %the pooled ceilings below. Fail loudly instead.
                    testCase.verifyTrue(all(isfinite(got(:))) && all(isfinite(want(:))), ...
                        sprintf('%s/%s: non-finite value in the %s comparison.', ...
                        run, design, checks{k, 2}));
                    d = got - want;
                    if strcmp(checks{k, 3}, 'rel')
                        d = d ./ max(abs(want), realmin);
                        worst_rel = max(worst_rel, max(abs(d)));
                    else
                        worst_abs = max(worst_abs, max(abs(d)));
                    end
                    fprintf('  %-28s mean %+.3e  max|.| %.3e  rms %.3e (%s)\n', ...
                        checks{k, 2}, mean(d), max(abs(d)), sqrt(mean(d.^2)), checks{k, 3});
                end

                % Rank agreement across participants: the per-participant
                % ORDERING is the clinically consequential output.
                rho_s = localSpearman(gen.rel_pt, ...
                    localReferenceColumn(data, 'G_delta', design, 'mean'));
                fprintf('[%s/%s] Spearman(G) across participants: %.6f\n', run, design, rho_s);

                % Interval endpoints: R type-7 vs MATLAB quantile at 500
                % draws; report and hold loosely.
                ci_worst = 0;
                ci_pairs = {gen.rel_ll, 'G_delta', 'q025'; gen.rel_ul, 'G_delta', 'q975'; ...
                            dep.rel_ll, 'D_delta', 'q025'; dep.rel_ul, 'D_delta', 'q975'};
                for k = 1:size(ci_pairs, 1)
                    want = localReferenceColumn(data, ci_pairs{k, 2}, design, ci_pairs{k, 3});
                    testCase.verifyTrue(all(isfinite(ci_pairs{k, 1}(:))) && all(isfinite(want(:))), ...
                        sprintf('%s/%s: non-finite %s %s interval endpoint.', ...
                        run, design, ci_pairs{k, 2}, ci_pairs{k, 3}));
                    ci_worst = max(ci_worst, max(abs(ci_pairs{k, 1} - want)));
                end
                fprintf('[%s/%s] worst G/D interval endpoint difference: %.3e\n', ...
                    run, design, ci_worst);

                % MEASURED 2026-08-10 across all four run x design cells:
                % bounded coefficients agreed to <= 6.4e-14 (worst:
                % residual_correlation_person), score-scale quantities to
                % <= 7.4e-13 relative (worst: residual_covariance_person,
                % which folds the quadrature), Spearman was exactly 1.000000
                % everywhere, and the worst interval endpoint difference was
                % 3.3e-3 (R type-7 vs MATLAB quantile interpolation at 500
                % draws - a method difference, not an arithmetic one). The
                % ceilings sit two to three orders above those measured
                % values; like the stage-2 ceiling they are regression pins on
                % the measured agreement, not detectability claims.
                testCase.verifyLessThan(worst_abs, 1e-11, ...
                    sprintf('%s/%s: a bounded coefficient differs beyond the measured-agreement ceiling.', run, design));
                testCase.verifyLessThan(worst_rel, 1e-10, ...
                    sprintf('%s/%s: a score-scale quantity differs beyond the measured-agreement ceiling.', run, design));
                testCase.verifyGreaterThan(rho_s, 0.999999, ...
                    sprintf('%s/%s: participant ordering disagrees with the reference.', run, design));
                testCase.verifyLessThan(ci_worst, 0.01, ...
                    sprintf('%s/%s: an interval endpoint differs beyond quantile-method slack.', run, design));
            end
        end

    end
end

% ---------------------------------------------------------------------------
% Bundle discovery
% ---------------------------------------------------------------------------

function root = localBundleRoot()
%LOCALBUNDLEROOT Reference bundle location, mirroring the fixture generator:
%   the MPSDC_REFERENCE_DIR environment variable, else the Desktop default.
root = getenv('MPSDC_REFERENCE_DIR');
if isempty(root)
    root = fullfile(localHome(), 'Desktop', 'stan_chisquare', ...
        'person_specific_dynamic_concurrent_modular');
end
end

function home = localHome()
%LOCALHOME The user home directory without the Statistics/OS toolboxes.
if ispc
    home = getenv('USERPROFILE');
else
    home = getenv('HOME');
end
end

function msg = localSkipMessage()
%LOCALSKIPMESSAGE One place for the assume text, so every gate reads the same.
msg = ['The modular cut-copula reference bundle is not present (set ', ...
    'MPSDC_REFERENCE_DIR or place it at ~/Desktop/stan_chisquare/', ...
    'person_specific_dynamic_concurrent_modular). External replay skipped.'];
end

function sel = localReadSelection(info)
%LOCALREADSELECTION The retained-draw table with sane column names. The CSV
%   headers are posterior-package names (.chain/.iteration/.draw), which MATLAB
%   would otherwise mangle version-dependently; preserve then rename.
raw = readtable(info.selection_file, 'VariableNamingRule', 'preserve');
sel = table(raw.selected_row, raw.('.chain'), raw.('.iteration'), raw.('.draw'), ...
    'VariableNames', {'selected_row', 'chain', 'iteration', 'draw'});
end

function flag = localLogical(x)
%LOCALLOGICAL Robust logical parse: R writes TRUE/FALSE unquoted, which
%   readtable may deliver as logical, numeric, cellstr or string.
if islogical(x)
    flag = x;
elseif isnumeric(x)
    flag = x ~= 0;
else
    flag = strcmpi(strtrim(cellstr(x)), 'TRUE');
end
end

function info = localRunInfo(run)
%LOCALRUNINFO Paths for one completed reference run, or [] when absent.
%   A run counts as present only with its COMPLETION manifest: the bundle also
%   carries two prepared-but-never-sampled two-facet runs whose chain files
%   hold zero draw rows, and those must never be mistaken for replayable fits.
prefix = sprintf(['theta_minus_alpha_fla_cor_trial_%s_threaded_modular_', ...
    'observed_full_person_specific_logsigma_cut_copula_chisq_v1'], run);
base = fullfile(localBundleRoot(), ...
    'person_specific_dynamic_concurrent_modular_reliability_cut_copula_v1', prefix);

info = struct();
info.prefix = prefix;
info.base = base;
info.completion_file = fullfile(base, [prefix '_completion_manifest.csv']);
info.standata_file = fullfile(base, [prefix '_standata.json']);
info.selection_file = fullfile(base, 'marginal_fit', [prefix '_copula_draw_selection.csv']);
info.copula_file = fullfile(base, 'modular_copula', [prefix '_modular_copula_draws.csv']);
info.copula_summary_file = fullfile(base, 'modular_copula', [prefix '_modular_copula_summary.csv']);
info.person_file = fullfile(base, 'analysis_data', [prefix '_person_table.csv']);
info.summary_file = fullfile(base, 'person_specific_reliability', ...
    [prefix '_person_reliability_summary.csv']);

needed = {info.completion_file, info.standata_file, info.selection_file, ...
    info.copula_file, info.copula_summary_file, info.person_file, info.summary_file};
if ~all(cellfun(@(f) exist(f, 'file') == 2, needed))
    info = [];
end
end

% ---------------------------------------------------------------------------
% Reference-run loading (cached: the chain read is the expensive part)
% ---------------------------------------------------------------------------

function data = localLoadRun(run)
%LOCALLOADRUN Load one reference run: standata, the RETAINED chain draws, the
%   copula outputs, the person table and the reliability summary. Cached per
%   run label because two test methods need the same ~40 MB of parsing.
persistent cache
if isempty(cache)
    cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end
if isKey(cache, run)
    data = cache(run);
    return
end

info = localRunInfo(run);
if isempty(info)
    error('TestModularCutCopulaExternalReplay:bundle', ...
        'localLoadRun called with the bundle absent; the assume gate failed.');
end

data = struct();
data.info = info;
data.ndim = 1 + strcmp(run, 'two_dimensions');

% ---- standata (trap 1 carrier: the offsets ride in here) ----------------
raw = jsondecode(fileread(info.standata_file));
data.standata = struct('N', raw.N, 'K', raw.K, 'Z', raw.Z, ...
    'y_1', raw.y_theta, 'y_2', raw.y_alpha, ...
    'p_id', raw.p_id, 'i_id', raw.i_id, ...
    'mu_offset', raw.mu_offset(:)', 'log_sigma_offset', raw.log_sigma_offset(:)');
if data.standata.K == 1
    data.standata.Z = data.standata.Z(:);
end
data.P = raw.P;
data.I = raw.I;

% ---- retained-draw alignment (trap 3) -----------------------------------
sel_csv = localReadSelection(info);
data.sel = sel_csv;

% ---- chain CSVs, retained rows only -------------------------------------
[mat, names] = localReadRetainedChains(info, sel_csv);
data.draws = localDrawStruct(mat, names, data.standata.K, data.P, data.I);

% ---- reference stage-2/3 outputs ----------------------------------------
data.copula = readtable(info.copula_file, 'VariableNamingRule', 'preserve');
% The copula draws must sit on the SAME retained set, in the same order.
if ~isequal(data.copula.('.draw'), sel_csv.draw)
    error('TestModularCutCopulaExternalReplay:alignment', ...
        ['The reference modular_copula_draws rows are not the retained ', ...
        'stage-1 draws in order; the replay cannot pair rho with margins.']);
end
data.copula_summary = readtable(info.copula_summary_file);
data.person = readtable(info.person_file);
summary = readtable(info.summary_file);
% Long format: one row per (participant, design, quantity). Keep as a table
% and index on demand; localReferenceColumn does the ordering.
data.summary = summary;

cache(run) = data;
end

function [mat, names] = localReadRetainedChains(info, sel_csv)
%LOCALREADRETAINEDCHAINS Parse the four CmdStan chain CSVs and keep only the
%   retained iterations, stacked chain-major so row k corresponds to
%   sel_csv row k (.draw ascending). Only the parameter columns stage 2/3
%   consume are read ('%f'); everything else is skipped ('%*f') to halve the
%   transient memory of a 10,000 x 2,892 parse.
first_chain = localChainFile(info, 1);
fid = fopen(first_chain, 'r');
if fid < 0
    error('TestModularCutCopulaExternalReplay:chains', ...
        'Cannot open chain CSV %s.', first_chain);
end
cleanup = onCleanup(@() fclose(fid));
line = fgetl(fid);
while ischar(line) && startsWith(line, '#')
    line = fgetl(fid);
end
allnames = strsplit(line, ',');
clear cleanup

wanted = localWantedColumn(allnames);
names = allnames(wanted);
fmt = repmat({'%*f'}, 1, numel(allnames));
fmt(wanted) = {'%f'};
fmt = strjoin(fmt, '');

blocks = cell(4, 1);
for c = 1:4
    fid = fopen(localChainFile(info, c), 'r');
    if fid < 0
        error('TestModularCutCopulaExternalReplay:chains', ...
            'Cannot open chain %d CSV.', c);
    end
    cleanup = onCleanup(@() fclose(fid));
    % Skip the leading comments and the header, then parse; trailing comment
    % blocks (adaptation, timing) are eaten by CommentStyle.
    line = fgetl(fid);
    while ischar(line) && startsWith(line, '#')
        line = fgetl(fid);
    end
    parsed = textscan(fid, fmt, 'Delimiter', ',', 'CommentStyle', '#', ...
        'CollectOutput', true);
    clear cleanup
    chain = parsed{1};

    iters = sel_csv.iteration(sel_csv.chain == c);
    if max(iters) > size(chain, 1)
        error('TestModularCutCopulaExternalReplay:chains', ...
            ['Chain %d holds %d draws but the retained set references ', ...
            'iteration %d; the chain files do not match the selection.'], ...
            c, size(chain, 1), max(iters));
    end
    blocks{c} = chain(iters, :);
end
mat = vertcat(blocks{:});
end

function file = localChainFile(info, c)
%LOCALCHAINFILE The cmdstanr chain CSV path for chain c. The per-chain suffix
%   is cmdstanr's internal counter and varies (chain1-1 but chain2-2), so the
%   file is globbed rather than assumed; exactly one CSV must match.
folder = fullfile(info.base, 'chains', sprintf('chain%d', c));
hits = dir(fullfile(folder, sprintf('%s_chain%d-*.csv', info.prefix, c)));
if numel(hits) ~= 1
    error('TestModularCutCopulaExternalReplay:chains', ...
        'Expected exactly one chain-%d CSV under %s; found %d.', ...
        c, folder, numel(hits));
end
file = fullfile(folder, hits(1).name);
end

function keep = localWantedColumn(allnames)
%LOCALWANTEDCOLUMN Logical mask over CSV columns for the parameters the
%   replay consumes. CmdStan writes array elements as name.i or name.i.j.
prefixes = {'alpha_mu.', 'alpha_sigma.', 'beta_mu.', 'beta_sigma.', ...
    'sd_p_joint.', 'R_p_joint.', 'sd_i.', 'u_p_joint.', 'u_i.'};
keep = strcmp(allnames, 'cor_i');
for k = 1:numel(prefixes)
    keep = keep | startsWith(allnames, prefixes{k});
end
end

function draws = localDrawStruct(mat, names, K, P, I)
%LOCALDRAWSTRUCT Assemble the stage-2/3 draw arrays from named columns.
%   Every element is fetched BY NAME (a one-off map), so the assembly is
%   independent of the CSV's column order and fails loudly on a missing name.
idx = containers.Map(names, 1:numel(names));
col = @(name) localColumn(mat, idx, name);
nd = size(mat, 1);

draws = struct();
draws.alpha_mu = [col('alpha_mu.1'), col('alpha_mu.2')];
draws.alpha_sigma = [col('alpha_sigma.1'), col('alpha_sigma.2')];

draws.beta_mu = zeros(nd, 2, K);
draws.beta_sigma = zeros(nd, 2, K);
for m = 1:2
    for k = 1:K
        draws.beta_mu(:, m, k) = col(sprintf('beta_mu.%d.%d', m, k));
        draws.beta_sigma(:, m, k) = col(sprintf('beta_sigma.%d.%d', m, k));
    end
end

draws.sd_p_joint = zeros(nd, 4);
for k = 1:4
    draws.sd_p_joint(:, k) = col(sprintf('sd_p_joint.%d', k));
end
draws.R_p_joint = zeros(nd, 4, 4);
for i = 1:4
    for j = 1:4
        draws.R_p_joint(:, i, j) = col(sprintf('R_p_joint.%d.%d', i, j));
    end
end
draws.sd_i = [col('sd_i.1'), col('sd_i.2')];
draws.cor_i = col('cor_i');

draws.u_p = zeros(nd, P, 4);
for c = 1:4
    for p = 1:P
        draws.u_p(:, p, c) = col(sprintf('u_p_joint.%d.%d', p, c));
    end
end
draws.u_i = zeros(nd, I, 2);
for c = 1:2
    for i = 1:I
        draws.u_i(:, i, c) = col(sprintf('u_i.%d.%d', i, c));
    end
end
end

function v = localColumn(mat, idx, name)
%LOCALCOLUMN One named column, with a diagnosable failure on absence.
if ~isKey(idx, name)
    error('TestModularCutCopulaExternalReplay:column', ...
        ['Chain CSV column ''%s'' was not found; the parameter layout ', ...
        'differs from the documented reference model.'], name);
end
v = mat(:, idx(name));
end

% ---------------------------------------------------------------------------
% Stage-3 input assembly
% ---------------------------------------------------------------------------

function args = localStageThreeArgs(data)
%LOCALSTAGETHREEARGS Name/value arguments for the per-participant calculator,
%   built from the retained reference draws plus the REFERENCE'S OWN rho_e.
%
%   Trap (1) lives here: the calculator takes ABSOLUTE intercepts and
%   per-subject absolutes, so mu_offset / log_sigma_offset are added back to
%   every log-mean and log-residual-SD quantity.
d = data.draws;
mu_off = data.standata.mu_offset;
ls_off = data.standata.log_sigma_offset;
K = data.standata.K;
nd = size(d.alpha_mu, 1);

b = d.alpha_mu + repmat(mu_off, nd, 1);
b_sigma = d.alpha_sigma + repmat(ls_off, nd, 1);

% Event-crossed slope layout: odd columns event 1, even columns event 2,
% dimension terms in the reference's Z order (Dim1[, Dim2, Dim1*Dim2]).
b_dim = zeros(nd, 2 * K);
b_sigma_dim = zeros(nd, 2 * K);
for k = 1:K
    b_dim(:, 2*k - 1) = d.beta_mu(:, 1, k);
    b_dim(:, 2*k) = d.beta_mu(:, 2, k);
    b_sigma_dim(:, 2*k - 1) = d.beta_sigma(:, 1, k);
    b_sigma_dim(:, 2*k) = d.beta_sigma(:, 2, k);
end

% Per-subject ABSOLUTE log mean and log residual SD at z = 0: intercept plus
% the participant's own effect, offsets restored. Person-block column order is
% [mean_1, mean_2, logsigma_1, logsigma_2], the same order sd_id uses.
mu_ss1 = mu_off(1) + d.alpha_mu(:, 1) + squeeze(d.u_p(:, :, 1));
mu_ss2 = mu_off(2) + d.alpha_mu(:, 2) + squeeze(d.u_p(:, :, 2));
er_ss1 = ls_off(1) + d.alpha_sigma(:, 1) + squeeze(d.u_p(:, :, 3));
er_ss2 = ls_off(2) + d.alpha_sigma(:, 2) + squeeze(d.u_p(:, :, 4));

% Trap (2)/(3) closure: the rho fed forward is the REFERENCE'S own draw per
% retained marginal draw, so both implementations run identical inputs and
% the comparison is arithmetic, not Monte Carlo.
rho = data.copula.rho_e;

args = { ...
    'b', b, 'b_sigma', b_sigma, 'b_dim', b_dim, 'b_sigma_dim', b_sigma_dim, ...
    'sd_id', d.sd_p_joint, 'sd_trl', d.sd_i, ...
    'cor_id', d.R_p_joint, 'cor_i', d.cor_i, ...
    'er_ss1', er_ss1, 'er_ss2', er_ss2, ...
    'mu_ss1', mu_ss1, 'mu_ss2', mu_ss2, 'rho', rho};
end

function idtable = localIdTable(data, design)
%LOCALIDTABLE The ssinfo table for one D-study design. 'actual' evaluates each
%   participant at their own trial count; 'standardized' at the reference's
%   common n_i-prime (the summary carries it, identical for every participant).
p = data.person;
n = height(p);
if strcmp(design, 'actual')
    trls = p.n_i_actual;
else
    rows = strcmp(data.summary.design, 'standardized');
    trls = repmat(data.summary.n_i_prime(find(rows, 1)), n, 1);
end
% The concurrent design scores both events from the SAME trials, so the two
% per-event counts are one number.
if data.ndim == 1
    idtable = table(p.p_index, p.p_index, p.Dim1_z, trls, trls, trls, ...
        'VariableNames', {'id', 'id2', 'z1', 'trls1', 'trls2', 'trls'});
else
    idtable = table(p.p_index, p.p_index, p.Dim1_z, p.Dim2_z, trls, trls, trls, ...
        'VariableNames', {'id', 'id2', 'z1', 'z2', 'trls1', 'trls2', 'trls'});
end
end

% ---------------------------------------------------------------------------
% Reference summary access and small numerics
% ---------------------------------------------------------------------------

function v = localReferenceColumn(data, quantity, design, stat)
%LOCALREFERENCECOLUMN One summary statistic for one quantity at one design,
%   ordered by p_index to match the idtable/person-table row order.
rows = strcmp(data.summary.quantity, quantity) & strcmp(data.summary.design, design);
sub = data.summary(rows, :);
if height(sub) ~= height(data.person)
    error('TestModularCutCopulaExternalReplay:summary', ...
        'Expected one ''%s''/%s row per participant; found %d of %d.', ...
        quantity, design, height(sub), height(data.person));
end
[~, order] = sort(sub.p_index);
sub = sub(order, :);
v = sub.(stat);
end

function rho = localSpearman(a, b)
%LOCALSPEARMAN Spearman rank correlation without the Statistics Toolbox.
r = corrcoef(localTiedRank(a), localTiedRank(b));
rho = r(1, 2);
end

function r = localTiedRank(x)
%LOCALTIEDRANK Average ranks, ties shared (matches the usual definition).
x = x(:);
n = numel(x);
[s, ord] = sort(x);
r = zeros(n, 1);
i = 1;
while i <= n
    j = i;
    while j < n && s(j + 1) == s(j)
        j = j + 1;
    end
    r(ord(i:j)) = (i + j) / 2;
    i = j + 1;
end
end

function r = localRoundHalfToEven(x)
%LOCALROUNDHALFTOEVEN Round to nearest, ties to even (R's rounding), the same
%   rule psyrat_gamma_cut_copula_stage documents for its retained-draw set.
r = round(x);
half = abs(x - fix(x)) == 0.5;
r(half) = 2 * round(x(half) / 2);
end
