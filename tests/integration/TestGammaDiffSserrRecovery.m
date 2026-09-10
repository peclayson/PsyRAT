classdef TestGammaDiffSserrRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square
    % SUBJECT-LEVEL two-event DIFFERENCE-score design (analysis 8,
    % ic_diff_sserrvar, NON-concurrent). The gamma case-8 model is the bivariate
    % one-facet gamma difference model (case 7); the subject-level distinction is
    % the per-participant CONDITIONAL residual per event
    %   sigma_pi_e2(e,p) = 2*exp(2*(b[e] + u_pe) + 2*sd_trl[e]^2 - (b_nu[e] + v_pe))
    % (psyrat_gamma_extract_diff_sserr), fed with the group person cross-covariance
    % through psyrat_ssrel_diff. This test simulates a bivariate person block
    % (correlated event means + per-event person-VARYING dispersion) and confirms:
    %   (1) the person log-nu SD s_v is recovered per event (the subject-level
    %       mechanism);
    %   (2) the observed-scale cross-event person covariance (the difference
    %       numerator's off-diagonal) is recovered;
    %   (3) the recovered per-subject difference residual TRACKS the truth across
    %       participants (subject-level differentiation); and
    %   (4) the true per-subject difference generalizability coefficients fall
    %       inside the recovered 95% credible intervals (coverage), with the
    %       sample-mean coefficient recovered.
    %
    % Ground truth uses the closed forms validated offline (no CmdStan) in
    % TestGammaDiffSserrReliability + TestGammaSserrConversion. CmdStan-gated
    % (skips Incomplete, not Failed, when CmdStan/MatlabStan absent).

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end

    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            localCleanupArtifacts(localCmdStanOutputDir(testCase.projectRoot()));
        end
    end

    methods (Test)
        function testSubjectLevelDifferenceReliabilityRecovered(testCase)
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required to draw continuous-df scaled chi-square.');

            % ---- Known truth (log expected-score scale), per event ----
            b     = [log(5.5), log(6.0)];    % per-event log-mean intercepts
            b_nu  = [log(15),  log(13)];     % per-event log-nu intercepts
            s_p   = [0.40, 0.42];            % per-event person log-mean SDs
            s_v   = [0.55, 0.55];            % per-event person log-nu SDs
            s_i   = [0.20, 0.22];            % per-event trial log-mean SDs
            rho_mean = 0.5;                  % cross-event person-mean correlation (log scale)
            rho_pv   = -0.3;                 % within-event mean<->log-nu correlation

            nsub = 60; ntrl = 40;
            [dataTbl, u_p, v_p] = localSimulateGammaDiffSserr( ...
                nsub, ntrl, b, b_nu, s_p, s_v, s_i, rho_mean, rho_pv);

            % ---- Closed-form observed-scale truth ----
            % numerator: per-event person universe variance (estimand-invariant) and
            % the observed cross-event person covariance (difference off-diagonal)
            sigma_p2 = zeros(1,2);
            for e = 1:2
                cov_pv_e = rho_pv * s_p(e) * s_v(e);
                vc = psyrat_gamma_varcomps(b(e), s_p(e), s_i(e), exp(b_nu(e)), s_v(e), cov_pv_e);
                sigma_p2(e) = vc.sigma_p2;
            end
            cc = psyrat_gamma_crosscov(b(1), s_p(1), s_i(1), b(2), s_p(2), s_i(2), rho_mean, 0);
            cov_p_obs_true = cc.cov_p_obs;
            var_pdiff_true = sigma_p2(1) + sigma_p2(2) - 2*cov_p_obs_true;

            % denominator: each participant's CONDITIONAL per-event residual, summed
            % over events (non-concurrent) and projected over the per-event n'
            trueResid = zeros(nsub,1);
            for e = 1:2
                trueResid = trueResid + (2 .* exp(2*(b(e) + u_p(:,e)) + 2*s_i(e)^2 ...
                    - (b_nu(e) + v_p(:,e)))) ./ ntrl;
            end
            trueG = var_pdiff_true ./ (var_pdiff_true + trueResid);   % generalizability, nsub x 1

            % ---- Fit (case 8 gamma, non-concurrent) ----
            rel = localRunDiffSserrGamma(testCase, dataTbl, 'recov_gamma_diff_sserr');

            testCase.verifyEqual(rel.analysis, 'ic_diff_sserrvar');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % ---- Recovered per-subject difference reliability (mirrors psyrat_relsummary) ----
            idvar   = cell2mat(rel.out.id_varcov{1});    % nd x 2 x 2
            trlvar  = cell2mat(rel.out.trl_varcov{1});   % nd x 2 x 2
            bsigma  = cell2mat(rel.out.b_sigma{1});      % nd x 2
            er1     = rel.out.er_var_ss1{1};             % nd x nsub
            er2     = rel.out.er_var_ss2{1};             % nd x nsub
            wpss    = rel.out.wp_cov_ss{1};              % nd x nsub (0)
            id_match = rel.out.id_matches{1};

            idtable = id_match;
            idtable.trls1 = repmat(ntrl, height(idtable), 1);
            idtable.trls2 = repmat(ntrl, height(idtable), 1);
            idtable.trls  = min(idtable.trls1, idtable.trls2);

            diff_table = psyrat_ssrel_diff('bp',idvar,'bt',trlvar,'er_var',bsigma,...
                'wp_cov',zeros(size(bsigma,1),1),'idtable',idtable,'CI',0.95,...
                'est','gen','er_var_ss1',er1,'er_var_ss2',er2,'wp_cov_ss',wpss);

            % align recovered per-subject quantities to truth via the id string ('S###')
            nrows = height(diff_table);
            recG = diff_table.rel_pt; recLL = diff_table.rel_ll; recUL = diff_table.rel_ul;
            recResid = diff_table.ss_errvar;
            tG = zeros(nrows,1); tResid = zeros(nrows,1);
            for r = 1:nrows
                idx = sscanf(char(string(diff_table.id(r))), 'S%d');
                tG(r) = trueG(idx);
                tResid(r) = trueResid(idx);
            end

            % (1) s_v recovered per event (gamma_disp.s_v is [nd x 2])
            svHat = mean(rel.out.gamma_disp.s_v, 1);
            % (2) observed cross-event person covariance recovered (numerator off-diagonal)
            covHat = mean(idvar(:,1,2));
            fprintf(['[diff sserr recovery] s_v hat=[%.3f %.3f] (true %.2f)  ', ...
                'cov_p obs hat=%.3f (true %.3f)  corr(resid)=%.3f  coverage=%.2f  ', ...
                'mean recG=%.3f (true %.3f)\n'], svHat(1), svHat(2), s_v(1), ...
                covHat, cov_p_obs_true, localPearson(recResid, tResid), ...
                mean(tG >= recLL & tG <= recUL), mean(recG), mean(tG));

            testCase.verifyGreaterThan(min(svHat), 0.35, sprintf( ...
                'Recovered person log-nu SD [%.3f %.3f] far below true %.2f.', ...
                svHat(1), svHat(2), s_v(1)));
            testCase.verifyLessThan(max(svHat), 0.85);

            % (2) cross-event person covariance recovered with the right sign + scale
            testCase.verifyGreaterThan(covHat, 0, ...
                'Recovered cross-event person covariance should be positive.');
            testCase.verifyLessThan(abs(covHat - cov_p_obs_true), 0.5*abs(cov_p_obs_true) + 1, ...
                sprintf('Cross-event person covariance off: recovered %.3f vs true %.3f.', ...
                covHat, cov_p_obs_true));

            % (3) subject-level DIFFERENTIATION: recovered per-subject difference
            % residual tracks the truth across participants.
            rResid = localPearson(recResid, tResid);
            testCase.verifyGreaterThan(rResid, 0.5, sprintf( ...
                ['Recovered per-subject difference residual does not track the truth ', ...
                '(Pearson r=%.3f); subject-level differentiation not recovered.'], rResid));

            % (4) COVERAGE: most true per-subject coefficients inside their 95%% CIs,
            % and the sample-mean coefficient recovered.
            coverage = mean(tG >= recLL & tG <= recUL);
            testCase.verifyGreaterThanOrEqual(coverage, 0.80, sprintf( ...
                'Only %.0f%% of true per-subject coefficients inside the 95%% CIs.', ...
                100*coverage));
            testCase.verifyLessThan(abs(mean(recG) - mean(tG)), 0.05, sprintf( ...
                'Sample-mean per-subject difference reliability off: recovered %.3f vs true %.3f.', ...
                mean(recG), mean(tG)));
        end
    end
end

function [tbl, u_p, v_p] = localSimulateGammaDiffSserr(nsub, ntrl, b, b_nu, s_p, s_v, s_i, rho_mean, rho_pv)
% Two-event Persons x Trials long table from a bivariate log-linked mean
% (correlated event means + crossed per-event trial main effects) and a
% person-VARYING per-event scaled-chi-square dispersion. Returns per-subject
% mean effects u_p [nsub x 2] and log-nu effects v_p [nsub x 2] so the caller
% can form the CONDITIONAL per-subject difference residual truth.
% Y = gamrnd(nu_pe/2, 2*mu/nu_pe) has E[Y]=mu, Var[Y]=2*mu^2/nu_pe.
rng(12345);

% 4-D person block [u1, u2, v1, v2] with SDs [s_p1 s_p2 s_v1 s_v2] and
% correlations: mean1<->mean2 = rho_mean, mean_e<->lognu_e = rho_pv (per event).
R = eye(4);
R(1,2) = rho_mean; R(2,1) = rho_mean;
R(1,3) = rho_pv;   R(3,1) = rho_pv;
R(2,4) = rho_pv;   R(4,2) = rho_pv;
sd = [s_p(1), s_p(2), s_v(1), s_v(2)];
L = chol(R, 'lower');
Z = randn(4, nsub);
blk = (diag(sd) * L * Z)';          % nsub x 4
u_p = blk(:,1:2);                   % person mean effects, per event
v_p = blk(:,3:4);                   % person log-nu effects, per event

% crossed per-event trial main effects (independent across events)
u_t = [s_i(1) * randn(ntrl,1), s_i(2) * randn(ntrl,1)];

nrow = nsub * 2 * ntrl;
ids  = cell(nrow,1);
evs  = cell(nrow,1);
meas = zeros(nrow,1);
r = 0;
evnames = {'E1','E2'};
for p = 1:nsub
    for e = 1:2
        nu_pe = exp(b_nu(e) + v_p(p,e));
        for t = 1:ntrl
            r = r + 1;
            ids{r} = sprintf('S%03d', p);
            evs{r} = evnames{e};
            mu = exp(b(e) + u_p(p,e) + u_t(t,e));
            meas(r) = gamrnd(nu_pe/2, 2*mu/nu_pe);
        end
    end
end
tbl = table(ids, meas, evs, 'VariableNames', {'id', 'meas', 'event'});
end

function r = localPearson(x, y)
% Pearson correlation via base MATLAB (no Statistics Toolbox).
x = x(:) - mean(x(:)); y = y(:) - mean(y(:));
r = (x' * y) / (sqrt(x' * x) * sqrt(y' * y));
end

function rel = localRunDiffSserrGamma(testCase, dataTbl, tag)
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    ... %pin log-nu: this design is movable, and since rulings 33-36 the
    ... %engine default is location-scale. Without this the test would
    ... %silently change estimand and check log-nu targets against an
    ... %LS fit. The location-scale twin has its own *LsRecovery test.
    'gammascale', 1, ...
    'diffest', 2, ...
    'sserrvar', 2, ...
    'chains', 2, ...
    'warmup', 750, ...
    'sampling', 750, ...
    'seed', 12345, ...
    'verbose', 0, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

function localCleanupArtifacts(baseDir)
p = fullfile(baseDir, 'recov_gamma_diff_sserr');
if isfolder(p)
    try
        rmdir(p, 's');
    catch
        warning('tests:cleanupDirFailed', ...
            'Could not remove recovery artifact directory ''%s''.', p);
    end
end
pats = {'cmdstan*.stan', 'cmdstan*.hpp', 'cmdstan18-*', 'output-*.csv', 'temp.data.R'};
for i = 1:numel(pats)
    files = dir(fullfile(baseDir, pats{i}));
    for j = 1:numel(files)
        if ~files(j).isdir
            delete(fullfile(baseDir, files(j).name));
        end
    end
end
dirs = dir(fullfile(baseDir, 'cmdstan_*'));
dirs = dirs([dirs.isdir]);
dirs = dirs(~ismember({dirs.name}, {'.', '..'}));
for i = 1:numel(dirs)
    try
        rmdir(fullfile(baseDir, dirs(i).name), 's');
    catch
        % best-effort cleanup
    end
end
end
