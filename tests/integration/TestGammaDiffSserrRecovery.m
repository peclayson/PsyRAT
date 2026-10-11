classdef TestGammaDiffSserrRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square
    % SUBJECT-LEVEL two-event DIFFERENCE-score design (analysis 8,
    % ic_diff_sserrvar, NON-concurrent). The gamma case-8 model is the bivariate
    % one-facet gamma difference model (case 7); the subject-level distinction is
    % the per-participant POOLED residual per event (S20; psyrat_gamma_extract_diff_sserr)
    %   sigma_pi_e2(e,p) = M_e + 2*exp(2*(b[e] + u_pe) + 2*sd_trl[e]^2 - (b_nu[e] + v_pe)),
    % M_e the POPULATION person x trial mean interaction of event e and the
    % exponential the person's OWN trial-averaged observation variance, fed with
    % the group person cross-covariance through psyrat_ssrel_diff. This test
    % simulates a bivariate person block (correlated event means + per-event
    % person-VARYING dispersion) and confirms:
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
    % SCORE, ESTIMAND, CONDITIONING, COUNTS. The score is the difference between
    % the two events' ntrl-trial means of single-trial gamma observations. The
    % estimand is the per-person generalizability of that difference,
    %   G_p = S / (S + E_{1,p}/ntrl + E_{2,p}/ntrl),
    %   S   = P_1 + P_2 - 2 CovP12,     E_{e,p} = M_e + C_{e,p},
    % relative error. S, P_e, CovP12 and M_e are POPULATION quantities at the
    % event intercepts b(e); only C_{e,p} is conditioned on the person's own
    % (u_pe, v_pe). The person covariance belongs to the SIGNAL and is not
    % subtracted from the residual; the design is non-concurrent, so the
    % residual covariance is zero. Allocation: nsub = 60 persons x 2 events x
    % ntrl = 40 crossed trials per event, n' = 40 for both events.
    %
    % Ground truth is the independent 3A reference
    % (tests/helpers/PsyRAT3ARecoveryReference.m, pinned offline in
    % tests/TestGamma3ARecoveryReference.m); no production converter is called
    % for the truth. Before this correction the truth took P_e/CovP12 from
    % psyrat_gamma_varcomps/psyrat_gamma_crosscov and omitted M_e from the
    % per-event error; earlier logs remain historical evidence for that stale
    % target. Design, seed, sampler settings and every assertion tolerance are
    % unchanged. CmdStan-gated (skips Incomplete, not Failed, when
    % CmdStan/MatlabStan absent). When PSYRAT_2C_RECOVERY_CAPTURE_DIR is set,
    % data and truth are captured BEFORE the fit, a fit failure still captures
    % whatever CmdStan wrote, the complete run is captured BEFORE any assertion,
    % and the class teardown skips its scratch cleanup.

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end

    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            % When evidence capture is requested, leave the raw scratch untouched
            % so it survives even a failed capture or a failed test; root's
            % external artifact guard archives/restores the owned scratch.
            if ~isempty(PsyRAT2CRecoveryCapture.captureDir())
                fprintf('[diff sserr recovery] %s set; scratch cleanup skipped.\n', ...
                    PsyRAT2CRecoveryCapture.EnvVar);
                return;
            end
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

            % ---- Independent closed-form observed-scale truth (3A reference) ----
            % SIGNAL (population, estimand-invariant): per event, B_e = exp(2 b(e)
            % + s_p(e)^2 + s_i(e)^2), P_e = B_e expm1(s_p(e)^2), and the observed
            % cross-event person covariance CovP12 = m_1 m_2 expm1(rho_mean s_p(1) s_p(2))
            % with m_e = exp(b(e) + (s_p(e)^2 + s_i(e)^2)/2). S = P_1 + P_2 - 2 CovP12.
            % The dispersion parameters (s_v, rho_pv) shape the planted (u_p, v_p)
            % only; none of these population components depends on them.
            P_true = zeros(1,2); I_true = zeros(1,2); M_true = zeros(1,2); m_true = zeros(1,2);
            for e = 1:2
                comp = PsyRAT3ARecoveryReference.oneFacetMeanComponents(b(e), s_p(e), s_i(e));
                P_true(e) = comp.P; I_true(e) = comp.I; M_true(e) = comp.M;
                m_true(e) = PsyRAT3ARecoveryReference.eventMean(b(e), s_p(e), s_i(e));
            end
            cov_p_obs_true = PsyRAT3ARecoveryReference.personCrossCovariance( ...
                m_true(1), m_true(2), rho_mean * s_p(1) * s_p(2));
            var_pdiff_true = P_true(1) + P_true(2) - 2*cov_p_obs_true;

            % RESIDUAL: each event's POOLED own residual E_{e,p} = M_e + C_{e,p}
            % (M_e = B_e expm1(s_p(e)^2) expm1(s_i(e)^2) the population person x
            % trial mean interaction at b(e), added exactly once; C_{e,p} =
            % 2 exp(2 (b(e) + u_pe) + 2 s_i(e)^2 - b_nu(e) - v_pe) the person's own
            % trial-averaged observation variance), summed over events at the
            % per-event count: error_p = E_{1,p}/ntrl + E_{2,p}/ntrl, residual
            % covariance zero (non-concurrent). G_p = S / (S + error_p).
            C_true = zeros(nsub,2); E_true = zeros(nsub,2);
            for e = 1:2
                C_true(:,e) = PsyRAT3ARecoveryReference.ownObservationVariance( ...
                    b(e), u_p(:,e), s_i(e), b_nu(e), v_p(:,e));
                E_true(:,e) = PsyRAT3ARecoveryReference.pooledOwnResidual(M_true(e), C_true(:,e));
            end
            trueResid = PsyRAT3ARecoveryReference.differenceError( ...
                E_true(:,1), ntrl, E_true(:,2), ntrl, 0);                       % nsub x 1
            trueG = PsyRAT3ARecoveryReference.differenceCoefficient( ...
                var_pdiff_true, trueResid);                                     % generalizability, nsub x 1
            truth = struct('b', b, 'b_nu', b_nu, 's_p', s_p, 's_v', s_v, 's_i', s_i, ...
                'rho_mean', rho_mean, 'rho_pv', rho_pv, 'nsub', nsub, 'ntrl', ntrl, ...
                'P', P_true, 'I', I_true, 'M', M_true, 'm', m_true, ...
                'cov_p_obs', cov_p_obs_true, 'signal', var_pdiff_true, ...
                'u_p', u_p, 'v_p', v_p, 'C_p', C_true, 'E_p', E_true, ...
                'error_p', trueResid, 'G_p', trueG, ...
                'residual_definition', ['error_p = (M_1 + C_{1,p})/ntrl + (M_2 + C_{2,p})/ntrl ', ...
                '(population per-event mean interaction plus own observation variance, ', ...
                'residual covariance zero); G_p = S/(S + error_p), S = P_1 + P_2 - 2 CovP12'], ...
                'generator', 'continuous gamrnd(nu_pe/2, 2 mu/nu_pe), seed 12345');

            % ---- Fit (case 8 gamma, non-concurrent), with evidence capture ----
            % (a no-op unless the capture root is set; see localCaptureStage)
            fitTag = 'recov_gamma_diff_sserr';
            outDir = fullfile(localCmdStanOutputDir(testCase.projectRoot()), fitTag);
            captureMeta = struct('test', 'TestGammaDiffSserrRecovery/testSubjectLevelDifferenceReliabilityRecovered', ...
                'analysis', 'ic_diff_sserrvar', 'gammascale', 1, 'diffest', 2, 'chains', 2, ...
                'warmup', 750, 'sampling', 750, 'seed', 12345, 'nsub', nsub, 'ntrl', ntrl, ...
                'design', 'two-event non-concurrent difference, persons x trials, subject-level');
            localCaptureStage('static8_lognu_diff_sserr_prefit', [], truth, dataTbl, outDir, ...
                captureMeta, 'prefit');
            try
                rel = localRunDiffSserrGamma(testCase, dataTbl, fitTag);
            catch fitErr
                localCaptureStage('static8_lognu_diff_sserr_fitfailed', ...
                    struct('fit_error', getReport(fitErr, 'extended', 'hyperlinks', 'off')), ...
                    truth, dataTbl, outDir, captureMeta, 'fitfailed');
                rethrow(fitErr);
            end
            localCaptureStage('static8_lognu_diff_sserr', rel, truth, dataTbl, outDir, ...
                captureMeta, 'postfit');

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

function localCaptureStage(tag, rel, truth, dataTbl, outDir, meta, stage)
% Evidence capture through the existing PsyRAT2CRecoveryCapture API (a no-op
% unless PSYRAT_2C_RECOVERY_CAPTURE_DIR is set). Three stages, three tags:
%   prefit    - data + truth only (no REL, no CmdStan output yet). The helper
%               writes the MAT/CSV/JSON first and then reports the raw capture
%               as INCOMPLETE because no chain CSV exists; that is the expected
%               outcome here and is tolerated, so data and truth are on disk
%               before CmdStan starts.
%   fitfailed - the fit threw: preserve the error report with whatever CmdStan
%               wrote (possibly partial or no chain CSVs); an incomplete raw
%               capture is again tolerated, and the caller rethrows the fit error.
%   postfit   - the complete REL + raw output, BEFORE any assertion; a capture
%               failure here is a real failure and propagates (sibling pattern).
% Any other capture error (non-absolute root, overlap with the scratch, an
% existing destination) propagates at every stage. Nothing is deleted.
if isempty(PsyRAT2CRecoveryCapture.captureDir())
    return;
end
meta.stage = stage;
try
    PsyRAT2CRecoveryCapture.capture(tag, rel, truth, dataTbl, outDir, 'meta', meta);
catch capErr
    if strcmp(stage, 'postfit') || ...
            ~strcmp(capErr.identifier, 'PsyRAT2CRecoveryCapture:rawCaptureIncomplete')
        rethrow(capErr);
    end
    fprintf(['[3A recovery capture] %s (%s): data/truth preserved; raw CmdStan capture ', ...
        'incomplete as expected at this stage (%s)\n'], tag, stage, capErr.message);
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
