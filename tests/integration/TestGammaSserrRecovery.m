classdef TestGammaSserrRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square
    % SUBJECT-LEVEL design (sserr, analysis 6). The gamma subject-level model IS
    % the one-facet gamma model: the person-varying log-dispersion log nu_p =
    % beta + v_p already gives each participant their own scaled-chi-square
    % observation variance, so per-subject reliability comes from the POOLED
    % per-subject residual (S20; psyrat_gamma_extract_sserr)
    %   sigma_pi_e2(p) = M + 2*exp(2*(alpha + u_p) + 2*s_i^2 - (beta + v_p)),
    % where M is the POPULATION person x trial mean-interaction component and the
    % exponential is the person's OWN trial-averaged observation variance, fed
    % through the Gaussian subject-level path (psyrat_ssrel). This test simulates
    % a person-VARYING dispersion (so the per-subject residuals genuinely differ)
    % and confirms:
    %   (1) the fit recovers the person log-nu SD s_v (the mechanism that makes
    %       the design subject-level);
    %   (2) the recovered per-subject residual DIFFERENTIATES participants (it
    %       tracks the true per-subject residual across the sample); and
    %   (3) the true per-subject generalizability coefficients fall inside the
    %       recovered 95% credible intervals (coverage), with the sample-mean
    %       coefficient recovered.
    %
    % SCORE, ESTIMAND, CONDITIONING, COUNTS. The score is a single-trial gamma
    % observation with a log-linked mean; the reliability is that of each
    % person's ntrl-trial MEAN. The estimand is the per-person generalizability
    %   G_p = P / (P + E_p/ntrl),   E_p = M + C_p,
    % relative error (gcoeff 2, the trial main effect I excluded). P and M are
    % POPULATION quantities at the grand intercept alpha (not recomputed at
    % alpha + u_p); only C_p is conditioned on the person's own (u_p, v_p).
    % Allocation: nsub = 60 persons x ntrl = 40 crossed trials, one event, one
    % occasion; every person is projected at n' = ntrl = 40.
    %
    % Ground truth is the independent 3A reference
    % (tests/helpers/PsyRAT3ARecoveryReference.m, pinned offline in
    % tests/TestGamma3ARecoveryReference.m). No production converter is called
    % for the truth. Before this correction the truth took P from
    % psyrat_gamma_varcomps and used the observation-only residual (no M);
    % earlier logs remain historical evidence for that stale target. Design,
    % seed, sampler settings and every assertion tolerance are unchanged.
    % CmdStan-gated (skips Incomplete, not Failed, when CmdStan/MatlabStan
    % absent). When PSYRAT_2C_RECOVERY_CAPTURE_DIR is set, the simulated data
    % and planted truth are captured BEFORE the fit, a fit failure still
    % captures whatever CmdStan wrote, and the complete REL/raw output is
    % captured after the fit and BEFORE any assertion; the class teardown then
    % skips its scratch cleanup so the raw output survives.

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
                fprintf('[sserr recovery] %s set; scratch cleanup skipped.\n', ...
                    PsyRAT2CRecoveryCapture.EnvVar);
                return;
            end
            localCleanupArtifacts(localCmdStanOutputDir(testCase.projectRoot()));
        end
    end

    methods (Test)
        function testSubjectLevelReliabilityRecovered(testCase)
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required to draw continuous-df scaled chi-square.');

            % ---- Known truth (log expected-score scale) ----
            alpha = log(5.5);   % log-mean grand intercept
            s_p   = 0.40;       % person log-mean SD
            s_i   = 0.20;       % trial  log-mean SD
            beta  = log(15);    % population log-nu intercept
            s_v   = 0.60;       % person log-nu SD (dispersion heterogeneity)
            rho   = -0.30;      % person mean <-> log-nu correlation

            nsub = 60; ntrl = 40;
            [dataTbl, u_p, v_p] = localSimulateGammaSserr( ...
                nsub, ntrl, alpha, s_p, s_i, beta, s_v, rho);

            % ---- Independent closed-form observed-scale truth (3A reference) ----
            % P (signal) and M (person x trial mean interaction) are POPULATION
            % components at the grand intercept: B = exp(2 alpha + s_p^2 + s_i^2),
            % P = B expm1(s_p^2), M = B expm1(s_p^2) expm1(s_i^2). C_p is each
            % participant's OWN trial-averaged observation variance
            % 2 exp(2 (alpha + u_p) + 2 s_i^2 - beta - v_p). The pooled residual
            % E_p = M + C_p enters exactly once, and G_p = P / (P + E_p/ntrl) is
            % the relative-error coefficient psyrat_ssrel forms under gcoeff 2.
            % The dispersion parameters (s_v, rho) shape the planted (u_p, v_p)
            % only; P and M do not depend on them.
            comp = PsyRAT3ARecoveryReference.oneFacetMeanComponents(alpha, s_p, s_i);
            sigma_p2_true = comp.P;
            trueObs   = PsyRAT3ARecoveryReference.ownObservationVariance( ...
                alpha, u_p, s_i, beta, v_p);                                    % nsub x 1
            trueResid = PsyRAT3ARecoveryReference.pooledOwnResidual(comp.M, trueObs); % nsub x 1
            trueG = PsyRAT3ARecoveryReference.oneFacetCoefficients( ...
                sigma_p2_true, comp.I, trueResid, ntrl);                        % generalizability
            truth = struct('alpha', alpha, 's_p', s_p, 's_i', s_i, 'beta', beta, ...
                's_v', s_v, 'rho', rho, 'nsub', nsub, 'ntrl', ntrl, ...
                'B', comp.B, 'P', comp.P, 'I', comp.I, 'M', comp.M, ...
                'u_p', u_p, 'v_p', v_p, 'C_p', trueObs, 'E_p', trueResid, 'G_p', trueG, ...
                'residual_definition', ['E_p = M + C_p (population person x trial mean ', ...
                'interaction at alpha plus own observation variance); ', ...
                'G_p = P/(P + E_p/ntrl), relative error'], ...
                'generator', 'continuous gamrnd(nu_p/2, 2 mu/nu_p), seed 12345');

            % ---- Fit, with evidence capture (a no-op unless the capture root is set) ----
            % Stage 1 preserves the data and truth before CmdStan starts; stage 2
            % preserves a failed fit's error and partial output; stage 3 is the
            % sibling-pattern capture of the complete run before any assertion.
            fitTag = 'recov_gamma_sserr';
            outDir = fullfile(localCmdStanOutputDir(testCase.projectRoot()), fitTag);
            captureMeta = struct('test', 'TestGammaSserrRecovery/testSubjectLevelReliabilityRecovered', ...
                'analysis', 'ic_sserrvar', 'gammascale', 1, 'chains', 2, 'warmup', 750, ...
                'sampling', 750, 'seed', 12345, 'nsub', nsub, 'ntrl', ntrl, ...
                'design', 'one-facet persons x trials, subject-level');
            localCaptureStage('static6_lognu_sserr_prefit', [], truth, dataTbl, outDir, ...
                captureMeta, 'prefit');
            try
                rel = localRunSserrGamma(testCase, dataTbl, fitTag);
            catch fitErr
                localCaptureStage('static6_lognu_sserr_fitfailed', ...
                    struct('fit_error', getReport(fitErr, 'extended', 'hyperlinks', 'off')), ...
                    truth, dataTbl, outDir, captureMeta, 'fitfailed');
                rethrow(fitErr);
            end
            localCaptureStage('static6_lognu_sserr', rel, truth, dataTbl, outDir, ...
                captureMeta, 'postfit');

            testCase.verifyEqual(rel.analysis, 'ic_sserrvar');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % ---- Recovered per-subject reliability (mirrors psyrat_relsummary) ----
            idtable = rel.out.id_matches{1};
            idtable.trls = repmat(ntrl, height(idtable), 1);
            ssrel = psyrat_ssrel( ...
                'bp',     rel.out.gro_sds{1}(:,1), ...
                'wp_pop', rel.out.pop_sdlog(:,1), ...
                'wp_ss',  rel.out.ind_sdlog{1}, ...
                'i',      rel.out.sig_trl(:,1), ...
                'gcoeff', 2, ...
                'idtable', idtable, ...
                'CI', 0.95);

            % align the recovered per-subject quantities to the truth by parsing
            % the original id string ('S###' -> subject index)
            nrows = height(ssrel);
            recResid = zeros(nrows,1); tG = zeros(nrows,1);
            recG = ssrel.dep_pt; recLL = ssrel.dep_ll; recUL = ssrel.dep_ul;
            for r = 1:nrows
                idx = sscanf(char(string(ssrel.id(r))), 'S%d');
                recResid(r) = ssrel.ss_errvar(r);
                tG(r) = trueG(idx);
            end
            trueResidAligned = arrayfun(@(r) trueResid(sscanf(char(string(ssrel.id(r))),'S%d')), (1:nrows)');

            % (1) s_v recovered (person log-nu SD is gro_sds column 2, provenance)
            svHat = mean(cell2mat(rel.out.gro_sds{1}(:,2)));
            fprintf(['[sserr recovery] s_v hat=%.3f (true %.2f)  corr(resid)=%.3f  ', ...
                'coverage=%.2f  mean recG=%.3f (true %.3f)\n'], ...
                svHat, s_v, localPearson(recResid, trueResidAligned), ...
                mean(tG >= recLL & tG <= recUL), mean(recG), mean(tG));
            testCase.verifyGreaterThan(svHat, 0.40, sprintf( ...
                'Recovered person log-nu SD %.3f far below the true %.2f.', svHat, s_v));
            testCase.verifyLessThan(svHat, 0.85);

            % (2) subject-level DIFFERENTIATION: recovered per-subject residual
            % tracks the true per-subject residual across participants.
            rResid = localPearson(recResid, trueResidAligned);
            testCase.verifyGreaterThan(rResid, 0.5, sprintf( ...
                ['Recovered per-subject residual does not track the truth ', ...
                '(Pearson r=%.3f); subject-level differentiation not recovered.'], rResid));

            % (3) COVERAGE: most true per-subject coefficients fall inside their
            % recovered 95%% CIs, and the sample-mean coefficient is recovered.
            coverage = mean(tG >= recLL & tG <= recUL);
            testCase.verifyGreaterThanOrEqual(coverage, 0.80, sprintf( ...
                'Only %.0f%% of true per-subject coefficients inside the 95%% CIs.', ...
                100*coverage));
            testCase.verifyLessThan(abs(mean(recG) - mean(tG)), 0.05, sprintf( ...
                'Sample-mean per-subject reliability off: recovered %.3f vs true %.3f.', ...
                mean(recG), mean(tG)));
        end
    end
end

function [tbl, u_p, v_p] = localSimulateGammaSserr(nsub, ntrl, alpha, s_p, s_i, beta, s_v, rho)
% One-facet Persons x Trials long table from a log-linked mean (person + crossed
% trial main effects) and a person-VARYING scaled-chi-square dispersion. Returns
% the per-subject (u_p, v_p) so the caller can form the CONDITIONAL per-subject
% residual truth. Y = gamrnd(nu_p/2, 2*mu/nu_p) has E[Y]=mu, Var[Y]=2*mu^2/nu_p.
rng(12345);
u_t = s_i * randn(ntrl, 1);                          % crossed trial main effects
z1 = randn(nsub, 1); z2 = randn(nsub, 1);
u_p = s_p * z1;                                      % person mean effect
v_p = s_v * (rho * z1 + sqrt(1 - rho^2) * z2);       % person log-nu effect (corr rho)
nu_p = exp(beta + v_p);

nrow = nsub * ntrl;
ids  = cell(nrow, 1);
meas = zeros(nrow, 1);
r = 0;
for p = 1:nsub
    for t = 1:ntrl
        r = r + 1;
        ids{r}  = sprintf('S%03d', p);
        mu = exp(alpha + u_p(p) + u_t(t));
        meas(r) = gamrnd(nu_p(p) / 2, 2 * mu / nu_p(p));
    end
end
tbl = table(ids, meas, 'VariableNames', {'id', 'meas'});
end

function r = localPearson(x, y)
% Pearson correlation via base MATLAB (no Statistics Toolbox).
x = x(:) - mean(x(:)); y = y(:) - mean(y(:));
r = (x' * y) / (sqrt(x' * x) * sqrt(y' * y));
end

function rel = localRunSserrGamma(testCase, dataTbl, tag)
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
p = fullfile(baseDir, 'recov_gamma_sserr');
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
