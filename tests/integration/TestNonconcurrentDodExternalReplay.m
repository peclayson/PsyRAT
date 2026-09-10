classdef TestNonconcurrentDodExternalReplay < PsyRATTestBase
    % EXTERNAL replay of the owner's completed nonconcurrent DoD reference
    % fits (analysis 28, one-facet, person-specific dynamic nonconcurrent
    % difference-of-differences) through the MATLAB stage-3 machinery.
    %
    % WHAT THIS IS. The reference bundle
    %   person_specific_dynamic_nonconcurrent_dod_modular
    % contains two COMPLETED 151-person one-facet runs (observed data, gates
    % passed, six-pair residual-independence audit passed, 500 retained
    % reliability draws):
    %   ..._trial_one_dimension_..._observed_full_...   (K = 1)
    %   ..._trial_two_dimensions_..._observed_full_...  (K = 3)
    % This test pushes each run's own retained posterior draws through
    % PsyRAT's committed stage 3 (PSYRAT_SSREL_DODIFFDYNREL_GAMMA_LS) and
    % compares against the reference's own person_reliability_summary. Both
    % sides consume the SAME draws, so any disagreement is arithmetic, not
    % Monte Carlo - the same design as TestModularCutCopulaExternalReplay,
    % minus its stage 2: NO dependence module exists in this design, so the
    % replay is stage 3 alone.
    %
    % It validates the ONE-FACET path only: the bundle's two-facet runs are
    % synthetic-occasion development artifacts (publishable = FALSE), so
    % analysis 29's external grounding is the frozen two-facet fixture lane
    % plus live CmdStan recovery.
    %
    % THE RECORDED TRAPS, and where each is handled:
    %  (1) mu_offset / log_sigma_offset are DATA-DERIVED centres; the fitted
    %      alpha_mu / alpha_sigma are deviations. Stage 3 receives ABSOLUTE
    %      intercepts, so the offsets are added back in localStageThreeArgs.
    %  (2) Per-draw computation, summarized only afterwards (Jensen).
    %  (3) The retained-draw set comes from the reference's own
    %      reliability_draw_selection.csv; the loader reads exactly those
    %      chain rows, in that order.
    %  (4) OWNER RULING (2026-08-12): the toolbox implements the HARMONIC
    %      count rule only. The reference summary also carries exact-shared
    %      (*_exact_shared) rows; those quantities are deliberately not
    %      compared.
    %
    % SUMMARY CONVENTIONS. Point estimates are means over draws and match
    % the reference's `mean` column sharply (the comparisons that matter).
    % q025/q975 use MATLAB's quantile against R's type 7 at 500 draws, so
    % they get a looser ceiling. The reference's own runs carry small
    % posterior_invalid fractions (0.24% / 0.013%); reproducing the flagged
    % draws is part of the arithmetic and needs no special handling here
    % because flags ride per draw.
    %
    % BUNDLE-GATED: skips (Incomplete, not Failed) when the bundle is
    % absent. Override the location with MPSDOD_REFERENCE_DIR, mirroring
    % make_nonconcurrent_dod_fixtures.R.

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
        run = {'one_dimension', 'two_dimensions'}
    end

    methods (Test)

        function testStageThreeMatchesReferencePerParticipant(testCase, run)
            info = localRunInfo(run);
            testCase.assumeTrue(~isempty(info), localSkipMessage());

            data = localLoadRun(run);
            args = localStageThreeArgs(data);
            [~, dl] = psyrat_ssrel_dodiffdynrel_gamma_ls(args{:});

            % Comparison map: drawlevel field (or per-cell slice) against the
            % reference summary quantity, per design. Means sharp; intervals
            % looser (quantile-method difference at 500 draws).
            meanTol = 5e-7;    % relative; CSV round-trip + assembly order
            ciTol = 0.02;      % absolute on [0,1]-ish quantities at 500 draws

            P = data.P;
            for design = {'actual', 'standardized'}
                dsn = design{1};
                if strcmp(dsn, 'actual')
                    G = dl.G; D = dl.D; relv = dl.rel; absv = dl.abs_;
                else
                    G = dl.G_std; D = dl.D_std; relv = dl.rel_std; absv = dl.abs_std;
                end
                localCompare(testCase, data, dsn, 'G_dod', G, meanTol, ciTol, P);
                localCompare(testCase, data, dsn, 'D_dod', D, meanTol, ciTol, P);
                localCompare(testCase, data, dsn, 'relative_error_dod', relv, ...
                    meanTol, inf, P);
                localCompare(testCase, data, dsn, 'absolute_error_dod', absv, ...
                    meanTol, inf, P);
                localCompare(testCase, data, dsn, 'universe_dod', dl.u, ...
                    meanTol, inf, P);
            end
            % Count-free quantities (identical across designs; compare once)
            localCompare(testCase, data, 'actual', ...
                'ICC_relative_matched_facet_dod', dl.iccg, meanTol, ciTol, P);
            localCompare(testCase, data, 'actual', ...
                'ICC_absolute_matched_facet_dod', dl.iccd, meanTol, ciTol, P);
            % Per-cell person quantities: residual SDs and person-conditional
            % nu (the person-specific estimand's defining outputs)
            for q = 1:4
                localCompare(testCase, data, 'actual', ...
                    sprintf('residual_sd_%s_person', data.labels{q}), ...
                    squeeze(dl.sigma_person(:, q, :)), meanTol, inf, P);
                localCompare(testCase, data, 'actual', ...
                    sprintf('nu_%s_at_person_expected_mean', data.labels{q}), ...
                    squeeze(dl.nu_person(:, q, :)), meanTol, inf, P);
            end
            s4 = squeeze(sum(dl.sigma_person.^2, 2));
            localCompare(testCase, data, 'actual', 'residual_var_dod_person', ...
                s4, meanTol, inf, P);

            % The invalid fractions must reproduce the reference's flagged
            % draws exactly: flags are 0/1 per draw, so their means are
            % rationals with denominator 500 and match to fp round-off.
            localCompare(testCase, data, 'actual', ...
                'posterior_invalid_primary_any', dl.invalid, 1e-9, inf, P);
        end

    end
end

% ---------------------------------------------------------------------------

function localCompare(testCase, data, design, quantity, drawmat, meanTol, ciTol, P)
%DRAWMAT is ndraws x P in p_index order. The reference summary is long format
%keyed by (ID, design, quantity) with mean/q025/q975 columns.
ref = data.summary(strcmp(data.summary.design, design) & ...
    strcmp(data.summary.quantity, quantity), :);
testCase.assertEqual(height(ref), P, sprintf( ...
    '%s (%s): expected %d reference rows, found %d.', quantity, design, ...
    P, height(ref)));
%order reference rows by p_index so column p of drawmat pairs with row p
ref = sortrows(ref, 'p_index');
got_mean = mean(drawmat, 1, 'omitnan')';
scale = max(1, abs(ref.mean));
testCase.verifyLessThanOrEqual(max(abs(got_mean - ref.mean) ./ scale), ...
    meanTol, sprintf(['%s (%s): posterior means diverge from the ', ...
    'reference beyond the arithmetic ceiling.'], quantity, design));
if isfinite(ciTol)
    got_lo = quantile(drawmat, 0.025, 1)';
    got_hi = quantile(drawmat, 0.975, 1)';
    testCase.verifyLessThanOrEqual(max(abs(got_lo - ref.q025)), ciTol, ...
        sprintf('%s (%s): q025 beyond the quantile-method ceiling.', ...
        quantity, design));
    testCase.verifyLessThanOrEqual(max(abs(got_hi - ref.q975)), ciTol, ...
        sprintf('%s (%s): q975 beyond the quantile-method ceiling.', ...
        quantity, design));
end
end

function info = localRunInfo(run)
%Locate one completed observed-data one-facet run; [] when absent.
root = localBundleRoot();
resultsdir = fullfile(root, ...
    'person_specific_dynamic_nonconcurrent_dod_modular_reliability_zero_residual_covariance_v1');
stem = sprintf(['fla_error_correct_theta_alpha_dod_trial_%s_threaded_', ...
    'modular_observed_full_person_specific_logsigma_zero_residual_', ...
    'covariance_chisq_v1'], run);
rundir = fullfile(resultsdir, stem);
info = [];
if ~isfolder(rundir)
    return;
end
cand = struct();
cand.rundir = rundir;
cand.stem = stem;
cand.standata_file = fullfile(rundir, sprintf('%s_standata.json', stem));
cand.selection_file = fullfile(rundir, 'marginal_fit', ...
    sprintf('%s_reliability_draw_selection.csv', stem));
cand.person_file = fullfile(rundir, 'analysis_data', ...
    sprintf('%s_person_table.csv', stem));
cand.summary_file = fullfile(rundir, 'person_specific_reliability', ...
    sprintf('%s_person_reliability_summary.csv', stem));
need = {cand.standata_file, cand.selection_file, cand.person_file, ...
    cand.summary_file};
for i = 1:numel(need)
    if ~isfile(need{i})
        return;
    end
end
for c = 1:4
    if isempty(localChainFile(cand, c))
        return;
    end
end
info = cand;
end

function data = localLoadRun(run)
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
    error('TestNonconcurrentDodExternalReplay:bundle', ...
        'localLoadRun called with the bundle absent; the assume gate failed.');
end

data = struct();
data.info = info;
data.ndim = 1 + strcmp(run, 'two_dimensions');
data.labels = {'error_theta', 'correct_theta', 'error_alpha', 'correct_alpha'};

raw = jsondecode(fileread(info.standata_file));
data.K = raw.K;
data.P = raw.P;
data.mu_offset = raw.mu_offset(:)';
data.log_sigma_offset = raw.log_sigma_offset(:)';

sel = readtable(info.selection_file, 'VariableNamingRule', 'preserve');
data.sel = sel;

[mat, names] = localReadRetainedChains(info, sel);
data.draws = localDrawStruct(mat, names, data.K, data.P);

data.person = readtable(info.person_file, 'VariableNamingRule', 'preserve');
data.summary = readtable(info.summary_file, 'VariableNamingRule', 'preserve');

cache(run) = data;
end

function [mat, names] = localReadRetainedChains(info, sel)
%Parse the four cmdstanr chain CSVs and keep only the retained iterations,
%stacked in .draw order. Only stage-3 parameter columns are read.
first_chain = localChainFile(info, 1);
fid = fopen(first_chain, 'r');
assert(fid > 0, 'Cannot open chain CSV %s.', first_chain);
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
    assert(fid > 0, 'Cannot open chain %d CSV.', c);
    cleanup = onCleanup(@() fclose(fid));
    line = fgetl(fid);
    while ischar(line) && startsWith(line, '#')
        line = fgetl(fid);
    end
    parsed = textscan(fid, fmt, 'Delimiter', ',', 'CommentStyle', '#', ...
        'CollectOutput', true);
    clear cleanup
    chain = parsed{1};
    iters = sel.('.iteration')(sel.('.chain') == c);
    assert(max(iters) <= size(chain, 1), ...
        'Chain %d shorter than the retained selection.', c);
    blocks{c} = chain(iters, :);
end
mat = vertcat(blocks{:});
%rows must land in .draw order: sel is chain-major ascending by construction;
%verify rather than assume
[~, order] = sortrows([sel.('.chain'), sel.('.iteration')]);
assert(issorted(sel.('.draw')(order)), ...
    'Retained selection is not chain-major; the stacking would misalign.');
end

function keep = localWantedColumn(allnames)
%Stage-3 parameters only: alpha_mu, alpha_sigma, beta_mu, beta_sigma,
%sd_p_joint, R_p_joint, u_p_joint, sd_i, R_i.
prefixes = {'alpha_mu.', 'alpha_sigma.', 'beta_mu.', 'beta_sigma.', ...
    'sd_p_joint.', 'R_p_joint.', 'u_p_joint.', 'sd_i.', 'R_i.'};
keep = false(1, numel(allnames));
for i = 1:numel(allnames)
    for p = 1:numel(prefixes)
        if startsWith(allnames{i}, prefixes{p})
            keep(i) = true;
            break;
        end
    end
end
end

function draws = localDrawStruct(mat, names, K, P)
%Reassemble the named arrays from the retained chain matrix.
idx = containers.Map(names, num2cell(1:numel(names)));
nd = size(mat, 1);
g1 = @(nm) mat(:, idx(nm));
draws = struct();
draws.alpha_mu = zeros(nd, 4);
draws.alpha_sigma = zeros(nd, 4);
for q = 1:4
    draws.alpha_mu(:, q) = g1(sprintf('alpha_mu.%d', q));
    draws.alpha_sigma(:, q) = g1(sprintf('alpha_sigma.%d', q));
end
draws.beta_mu = zeros(nd, 4, K);
draws.beta_sigma = zeros(nd, 4, K);
for q = 1:4
    for k = 1:K
        draws.beta_mu(:, q, k) = g1(sprintf('beta_mu.%d.%d', q, k));
        draws.beta_sigma(:, q, k) = g1(sprintf('beta_sigma.%d.%d', q, k));
    end
end
draws.sd_p_joint = zeros(nd, 8);
for j = 1:8
    draws.sd_p_joint(:, j) = g1(sprintf('sd_p_joint.%d', j));
end
draws.R_p_joint = zeros(nd, 8, 8);
for i = 1:8
    for j = 1:8
        draws.R_p_joint(:, i, j) = g1(sprintf('R_p_joint.%d.%d', i, j));
    end
end
draws.u_p_joint = zeros(nd, P, 8);
for p = 1:P
    for j = 1:8
        draws.u_p_joint(:, p, j) = g1(sprintf('u_p_joint.%d.%d', p, j));
    end
end
draws.sd_i = zeros(nd, 4);
for q = 1:4
    draws.sd_i(:, q) = g1(sprintf('sd_i.%d', q));
end
draws.R_i = zeros(nd, 4, 4);
for i = 1:4
    for j = 1:4
        draws.R_i(:, i, j) = g1(sprintf('R_i.%d.%d', i, j));
    end
end
end

function args = localStageThreeArgs(data)
%Build the calculator inputs from the retained reference draws (trap 1: the
%offsets are added back so every intercept is ABSOLUTE).
d = data.draws;
b = data.mu_offset + d.alpha_mu;
b_sigma = data.log_sigma_offset + d.alpha_sigma;
mu_ss = cell(1, 4); er_ss = cell(1, 4);
for q = 1:4
    mu_ss{q} = b(:, q) + d.u_p_joint(:, :, q);
    er_ss{q} = b_sigma(:, q) + d.u_p_joint(:, :, 4 + q);
end

pt = data.person;
if data.ndim == 2
    idtable = table(pt.ID, pt.p_index, pt.Dim1_z, pt.Dim2_z, ...
        pt.n_i_error_actual, pt.n_i_correct_actual, ...
        pt.n_i_error_actual, pt.n_i_correct_actual, ...
        'VariableNames', {'id','id2','z1','z2','trls1','trls2','trls3','trls4'});
else
    idtable = table(pt.ID, pt.p_index, pt.Dim1_z, ...
        pt.n_i_error_actual, pt.n_i_correct_actual, ...
        pt.n_i_error_actual, pt.n_i_correct_actual, ...
        'VariableNames', {'id','id2','z1','trls1','trls2','trls3','trls4'});
end

%standardized n' from the reference's own standardized rows (settings
%defaults 49/341 unless overridden at run time)
stdrows = data.summary(strcmp(data.summary.design, 'standardized'), :);
nprime = [stdrows.n_i_error(1), stdrows.n_i_correct(1), ...
    stdrows.n_i_error(1), stdrows.n_i_correct(1)];

args = {'b', b, 'b_sigma', b_sigma, 'b_dim', d.beta_mu, ...
    'b_sigma_dim', d.beta_sigma, 'sd_id', d.sd_p_joint, ...
    'cor_id', d.R_p_joint, 'sd_trl', d.sd_i, 'cor_trl', d.R_i, ...
    'mu_ss1', mu_ss{1}, 'mu_ss2', mu_ss{2}, 'mu_ss3', mu_ss{3}, ...
    'mu_ss4', mu_ss{4}, 'er_ss1', er_ss{1}, 'er_ss2', er_ss{2}, ...
    'er_ss3', er_ss{3}, 'er_ss4', er_ss{4}, 'idtable', idtable, ...
    'nprime', nprime, 'ndim', data.ndim, 'CI', 0.95};
end

function file = localChainFile(info, c)
%The cmdstanr chain CSV (globbed; the per-chain suffix varies).
d = fullfile(info.rundir, 'chains', sprintf('chain%d', c));
hits = dir(fullfile(d, sprintf('%s_chain%d-*.csv', info.stem, c)));
if numel(hits) ~= 1
    file = '';
else
    file = fullfile(hits(1).folder, hits(1).name);
end
end

function root = localBundleRoot()
root = getenv('MPSDOD_REFERENCE_DIR');
if isempty(root)
    root = fullfile(localHome(), 'Desktop', 'stan_chisquare', ...
        'person_specific_dynamic_nonconcurrent_dod_modular');
end
end

function home = localHome()
if ispc
    home = getenv('USERPROFILE');
else
    home = getenv('HOME');
end
end

function msg = localSkipMessage()
msg = ['The nonconcurrent DoD reference bundle (with completed observed ', ...
    'one-facet runs) is not present (set MPSDOD_REFERENCE_DIR or place it ', ...
    'at ~/Desktop/stan_chisquare/person_specific_dynamic_nonconcurrent_', ...
    'dod_modular). External replay skipped.'];
end
