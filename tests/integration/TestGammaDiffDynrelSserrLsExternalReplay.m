classdef TestGammaDiffDynrelSserrLsExternalReplay < PsyRATTestBase
    % EXTERNAL replay of the owner's completed dynamic_nonconcurrent
    % location-scale reference fit (one dimension, observed data, 151 persons)
    % through PsyRAT's committed analysis-13 converter
    % (PSYRAT_SSREL_DIFFDYNREL_GAMMA_LS), compared against the reference's OWN
    % per-person outputs on disk. Both sides consume the SAME posterior draws,
    % so any disagreement is arithmetic, not Monte Carlo - the design the 8c
    % external validation established
    % (documentation/gamma_locationscale_rollout.md) and
    % TestModularCutCopulaExternalReplay follows for analysis 19.
    %
    % WHAT THIS PINS. The S17-FIX acceptance result (2026-08-11, recorded in
    % PSYRAT_AUDIT_FINDINGS.md, S17 resolution bullet): with the pooled
    % per-participant residual, PsyRAT's per-person G agrees with the
    % reference implementation's stored G_delta up to the ONE term the two
    % deliberately treat differently - the cross-event mean-surface covariance
    % the reference retains in its per-person pi_e slot and PsyRAT zeroes per
    % the S14 owner decision. Measured on this bundle: mean gap +1.23e-06, max
    % |gap| 3.69e-06, Spearman 1.000000. The PRE-fix conditional-only residual
    % produced mean +3.0e-05 / max +9.8e-05, so the ceilings below (pooled
    % regression pins ~5-8x above measured, ~3-5x BELOW the pre-fix defect)
    % fail loudly on a revert while tolerating platform noise.
    %
    % THE RECORDED TRAPS, and where each is handled here:
    %  (1) mu_offset / log_sigma_offset are DATA-DERIVED centering constants;
    %      the fitted intercepts are deviations from them. Read from the
    %      run's standata JSON and added back when b / b_sigma are built.
    %  (2) Posterior means are insufficient (Jensen): every PsyRAT quantity is
    %      computed per draw by the production converter and only then
    %      summarized, exactly as the reference's compute_person_draw does.
    %  (3) The draw subset must be ALIGNED: the reference's person estimand
    %      uses unique(round(seq(1, ndraws, length.out = 1000))) over the
    %      chain-major post-warmup draws (R/reliability.R, the keep rule;
    %      verified 2026-08-11 to reproduce the stored posterior_row set
    %      exactly). This test derives the same set with linspace and then
    %      PROVES the full alignment - columns, offsets, rows and z - by
    %      recomputing the reference's own stored per-person conditional
    %      residual means from the loaded draws (the mapping gate). A wrong
    %      column, offset, row set or z ordering cannot pass that gate.
    %
    % The stored per-person summary means are over the same 1,000 draws, so
    % the mean-vs-mean comparisons are exact up to the S14 term; interval
    % columns are NOT compared (R type-7 vs MATLAB quantile method).
    %
    % BUNDLE-GATED: skips (Incomplete, not Failed) when the bundle is absent.
    % Override the bundle location with the DNC_REFERENCE_DIR environment
    % variable, mirroring MPSDC_REFERENCE_DIR on the analysis-19 replay.
    %
    % Budget: ~2-4 minutes, dominated by reading the four chain CSVs
    % (selected columns only). No CmdStan and no R required.
    %
    % See also PSYRAT_SSREL_DIFFDYNREL_GAMMA_LS, PSYRAT_GAMMA_VARCOMPS_LS,
    % PSYRAT_GAMMA_CROSSCOV, TESTMODULARCUTCOPULAEXTERNALREPLAY.

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

    methods (Test)

        function testPooledPerParticipantReplayMatchesReference(testCase)
            info = localRunInfo();
            testCase.assumeTrue(~isempty(info), localSkipMessage());

            fprintf('\n[S17-FIX external replay] converter: %s\n', ...
                which('psyrat_ssrel_diffdynrel_gamma_ls'));

            % ---- reference per-person summary (design == 'actual') --------
            ref = localReadPersonSummary(info);
            P = height(ref.persons);
            testCase.assertEqual(P, 151, ...
                'The reference person table should carry 151 participants.');

            % ---- offsets ---------------------------------------------------
            standata = jsondecode(fileread(info.standata_file));
            testCase.assertEqual(double(standata.P), 151);
            mu_off = standata.mu_offset(:)';           % 1 x 2, [error correct]
            ls_off = standata.log_sigma_offset(:)';    % 1 x 2

            % ---- parameter draws on the reference's own subsample ----------
            [par, up3, up4] = localLoadDraws(info, P);
            keep = unique(round(linspace(1, par.ntotal, 1000)));
            fprintf('[S17-FIX external replay] %d post-warmup draws, %d kept\n', ...
                par.ntotal, numel(keep));

            b           = [mu_off(1) + par.alpha_mu(keep,1), ...
                           mu_off(2) + par.alpha_mu(keep,2)];
            b_sigma     = [ls_off(1) + par.alpha_sigma(keep,1), ...
                           ls_off(2) + par.alpha_sigma(keep,2)];
            b_dim       = par.beta_mu(keep,:);         % [event1 event2], K = 1
            b_sigma_dim = par.beta_sigma(keep,:);
            sd_id       = par.sd_p_joint(keep,:);      % [mu_e mu_c sig_e sig_c]
            sd_trl      = par.sd_i(keep,:);
            cor_p       = par.cor_p(keep);
            cor_i       = par.cor_i(keep);
            er_ss1      = b_sigma(:,1) + up3(keep,:);  % z = 0 per-person log-SD
            er_ss2      = b_sigma(:,2) + up4(keep,:);

            % ---- THE MAPPING GATE (trap 3) --------------------------------
            % Recompute the reference's own stored per-person conditional
            % residual means from the loaded draws. Any wrong column, offset,
            % row subset or z assignment fails here, before anything is
            % compared. Ceiling 1e-10: the identical arithmetic agrees to
            % ~1e-15; the slack is float-format headroom, not tolerance for a
            % mapping error (the nearest wrong mapping is orders off).
            z = ref.persons.Dim1_z(:)';
            condErr = mean(exp(2*(er_ss1 + b_sigma_dim(:,1) * z)), 1)';
            condCor = mean(exp(2*(er_ss2 + b_sigma_dim(:,2) * z)), 1)';
            gateErr = max(abs(condErr ./ ref.persons.resid_err - 1));
            gateCor = max(abs(condCor ./ ref.persons.resid_cor - 1));
            fprintf(['[S17-FIX external replay] mapping gate: max rel diff ' ...
                '%.3e (error ev), %.3e (correct ev)\n'], gateErr, gateCor);
            testCase.assertLessThan(max(gateErr, gateCor), 1e-10, ...
                ['Mapping gate failed: the loaded draws do not reproduce the ' ...
                'reference''s own stored per-person conditional residuals, so ' ...
                'a column, offset, row-subset or z mapping is wrong and no ' ...
                'downstream comparison is meaningful.']);

            % ---- signal gate: var_delta_p is parameterization-shared -------
            % Per draw through the PRODUCTION helpers, then averaged; the
            % reference's var_delta_p is the same arithmetic (8c: exact).
            sigU = zeros(P,1);
            for k = 1:P
                aA = b(:,1) + b_dim(:,1) * z(k);
                aB = b(:,2) + b_dim(:,2) * z(k);
                vcA = psyrat_gamma_varcomps_ls(aA, sd_id(:,1), sd_trl(:,1), 0, 0);
                vcB = psyrat_gamma_varcomps_ls(aB, sd_id(:,2), sd_trl(:,2), 0, 0);
                cc  = psyrat_gamma_crosscov(aA, sd_id(:,1), sd_trl(:,1), ...
                    aB, sd_id(:,2), sd_trl(:,2), cor_p, cor_i);
                sigU(k) = mean(vcA.sigma_p2 + vcB.sigma_p2 - 2*cc.cov_p_obs);
            end
            sigGate = max(abs(sigU ./ ref.persons.var_delta_p - 1));
            fprintf('[S17-FIX external replay] signal gate: max rel diff %.3e\n', ...
                sigGate);
            testCase.verifyLessThan(sigGate, 1e-9, ...
                ['The observed-scale signal must match the reference exactly ' ...
                '(shared log-mean submodel; 8c measured ~1e-15).']);

            % ---- the replay: production converter, per draw ----------------
            ss = table(ref.persons.ID, ref.persons.p_index, ref.persons.Dim1_z, ...
                ref.persons.n_i_error, ref.persons.n_i_correct, ...
                round((ref.persons.n_i_error + ref.persons.n_i_correct)/2), ...
                'VariableNames', {'id','id2','z1','trls1','trls2','trls'});
            cor_id = zeros(numel(keep), 4, 4);
            cor_id(:,1,2) = cor_p; cor_id(:,2,1) = cor_p;

            [gen, ~] = psyrat_ssrel_diffdynrel_gamma_ls( ...
                'b',b, 'b_sigma',b_sigma, 'b_dim',b_dim, ...
                'b_sigma_dim',b_sigma_dim, 'sd_id',sd_id, 'sd_trl',sd_trl, ...
                'cor_id',cor_id, 'cor_i',cor_i, ...
                'er_ss1',er_ss1, 'er_ss2',er_ss2, ...
                'idtable',ss, 'ndim',1, 'CI',0.95);

            gap = gen.rel_pt - ref.persons.G_delta;
            rho = localSpearman(gen.rel_pt, ref.persons.G_delta);
            fprintf(['[S17-FIX external replay] G vs reference: mean %+0.3e  ' ...
                'max|.| %.3e  RMS %.3e  Spearman %.6f\n'], ...
                mean(gap), max(abs(gap)), sqrt(mean(gap.^2)), rho);

            % Ceilings: measured mean +1.23e-06 / max 3.69e-06 (the S14 share
            % on this bundle); pre-fix defect mean +3.0e-05 / max +9.8e-05.
            % 1e-05 / 2e-05 sit ~8x/5x above measured and ~3x/5x below the
            % defect, so a conditional-only revert fails BOTH.
            testCase.verifyLessThan(abs(mean(gap)), 1e-5, ...
                ['Per-person G departs from the reference by more than the ' ...
                'S14 cross-term share - the pooled residual (S17-FIX) has ' ...
                'regressed.']);
            testCase.verifyLessThan(max(abs(gap)), 2e-5, ...
                ['At least one participant''s G departs from the reference ' ...
                'by more than the S14 cross-term share.']);
            testCase.verifyGreaterThan(rho, 0.99999, ...
                'Per-person ordering must match the reference.');
        end

    end
end

% ---------------------------------------------------------------------------
% Bundle discovery
% ---------------------------------------------------------------------------

function info = localRunInfo()
%LOCALRUNINFO Paths for the completed one-dimension LS observed-data run, or
%   [] when the bundle (or any required file) is absent.
root = getenv('DNC_REFERENCE_DIR');
if isempty(root)
    root = fullfile(localHome(), 'Desktop', 'stan_chisquare', ...
        'dynamic_nonconcurrent');
end
stem = ['fla_err_theta_minus_fla_cor_theta_person_trial_one_dimension_' ...
    'threaded_logsigma_absolute_sd_nonconcurrent_chisq_model_implied_v2_' ...
    'observed_full_fixed_direct_sd_v2_20260804'];
base = fullfile(root, 'dynamic_nonconcurrent_reliability', stem);

info = struct();
info.standata_file = fullfile(base, [stem '_standata.json']);
info.summary_file  = fullfile(base, 'reliability', ...
    [stem '_reliability_summary.csv']);
info.chain_files = cell(1,4);
for k = 1:4
    info.chain_files{k} = fullfile(base, 'chains', sprintf('chain%d', k), ...
        sprintf('%s_chain%d-%d.csv', stem, k, k));
end
required = [{info.standata_file, info.summary_file}, info.chain_files];
if ~all(cellfun(@(f) exist(f, 'file') == 2, required))
    info = [];
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
msg = ['The dynamic_nonconcurrent reference bundle is not present (set ', ...
    'DNC_REFERENCE_DIR or place it at ~/Desktop/stan_chisquare/', ...
    'dynamic_nonconcurrent). External replay skipped.'];
end

% ---------------------------------------------------------------------------
% Reference readers
% ---------------------------------------------------------------------------

function ref = localReadPersonSummary(info)
%LOCALREADPERSONSUMMARY Per-person reference quantities from the long-format
%   summary CSV, design == 'actual' (the participant's own trial counts - the
%   design PsyRAT's per-participant table uses). The stored `mean` column is
%   the mean over the SAME 1,000-draw subsample the parameter export uses.
raw = readtable(info.summary_file, 'VariableNamingRule', 'preserve');
raw = raw(strcmp(raw.design, 'actual'), :);

persons = unique(raw(:, {'p_index', 'ID', 'Dim1_z', 'n_i_error', ...
    'n_i_correct'}), 'rows');
persons = sortrows(persons, 'p_index');
persons.ID = string(persons.ID);

    function col = pick(qname)
        sub = raw(strcmp(raw.quantity, qname), {'p_index', 'mean'});
        sub = sortrows(sub, 'p_index');
        assert(isequal(sub.p_index, persons.p_index), ...
            'quantity %s does not cover every person exactly once', qname);
        col = sub.mean;
    end

persons.G_delta     = pick('G_delta');
persons.var_delta_p = pick('var_delta_p');
persons.resid_err   = pick('residual_var_error_person');
persons.resid_cor   = pick('residual_var_correct_person');

ref = struct('persons', persons);
end

function [par, up3, up4] = localLoadDraws(info, P)
%LOCALLOADDRAWS The needed CmdStan columns from all four chains, bound
%   chain-major (chain 1 rows first) - the same order posterior::as_draws_df
%   gives the reference's keep rule. Only the required 318 columns are read.
scalarNames = {'alpha_mu.1','alpha_mu.2','alpha_sigma.1','alpha_sigma.2', ...
    'beta_mu.1.1','beta_mu.2.1','beta_sigma.1.1','beta_sigma.2.1', ...
    'sd_p_joint.1','sd_p_joint.2','sd_p_joint.3','sd_p_joint.4', ...
    'sd_i.1','sd_i.2','R_p_joint.1.2','cor_i'};
upNames3 = arrayfun(@(p) sprintf('u_p_joint.%d.3', p), 1:P, 'UniformOutput', false);
upNames4 = arrayfun(@(p) sprintf('u_p_joint.%d.4', p), 1:P, 'UniformOutput', false);
wanted = [scalarNames, upNames3, upNames4];

blocks = cell(1, numel(info.chain_files));
for k = 1:numel(info.chain_files)
    opts = detectImportOptions(info.chain_files{k}, 'FileType', 'text', ...
        'CommentStyle', '#', 'VariableNamingRule', 'preserve');
    assert(all(ismember(wanted, opts.VariableNames)), ...
        'chain %d is missing expected parameter columns', k);
    opts.SelectedVariableNames = wanted;
    blocks{k} = readtable(info.chain_files{k}, opts);
end
tbl = vertcat(blocks{:});

par = struct();
par.ntotal      = height(tbl);
par.alpha_mu    = [tbl.('alpha_mu.1'),    tbl.('alpha_mu.2')];
par.alpha_sigma = [tbl.('alpha_sigma.1'), tbl.('alpha_sigma.2')];
par.beta_mu     = [tbl.('beta_mu.1.1'),   tbl.('beta_mu.2.1')];
par.beta_sigma  = [tbl.('beta_sigma.1.1'), tbl.('beta_sigma.2.1')];
par.sd_p_joint  = [tbl.('sd_p_joint.1'), tbl.('sd_p_joint.2'), ...
                   tbl.('sd_p_joint.3'), tbl.('sd_p_joint.4')];
par.sd_i        = [tbl.('sd_i.1'), tbl.('sd_i.2')];
par.cor_p       = tbl.('R_p_joint.1.2');
par.cor_i       = tbl.('cor_i');

up3 = zeros(par.ntotal, P);
up4 = zeros(par.ntotal, P);
for p = 1:P
    up3(:,p) = tbl.(sprintf('u_p_joint.%d.3', p));
    up4(:,p) = tbl.(sprintf('u_p_joint.%d.4', p));
end
end

function r = localSpearman(a, b)
%LOCALSPEARMAN Rank correlation without the Statistics Toolbox.
ra = localRanks(a(:)); rb = localRanks(b(:));
ra = ra - mean(ra); rb = rb - mean(rb);
r = (ra' * rb) / sqrt((ra' * ra) * (rb' * rb));
end

function r = localRanks(x)
[~, idx] = sort(x); r = zeros(size(x)); r(idx) = 1:numel(x);
end
