classdef TestTrtSserrRecovery < PsyRATTestBase
    % Parameter recovery for the subject-level crossed test-retest design
    % (trt_sserrvar, analysis 25).
    %
    % WHY THIS EXISTS. Analysis 25 shipped with golden Stan text, a routing
    % assertion in the end-to-end workflow test, and offline fixtures -- but no
    % numeric check of any kind. Its recovery was asserted only in the body of
    % the commit that added it ("live CmdStan recovery (per-subject residual
    % corr 0.98; sig_dp 0.44 vs 0.40)"), which is not a repeatable test. RC-24
    % then exposed the design through the GUI, so the numeric claim needs to
    % stand on something re-runnable.
    %
    % WHAT IT PINS. The design is Persons x Trials x Occasions with a
    % PER-SUBJECT residual. The interesting failure mode is the same one
    % TestSubjectLevelTrialRecovery pins for the one-facet case: a per-subject
    % residual is a greedy parameter, and if the crossed facet effects are not
    % separately identified it will absorb them. So the test checks that the
    % five crossed facet SDs and the per-subject residual recover SEPARATELY,
    % not just that the fit completes.
    %
    % Like the other CmdStan integration tests, dependency setup uses
    % assumptions, so this reports Incomplete (not Failed) when CmdStan or
    % MatlabStan is unavailable.
    %
    % FIRST LIVE RUN, 2026-07-28 (CmdStan 2.38.0, bundled MatlabStan, seed
    % 12345, 2 chains x 500/500, ~2 min). Recorded so a future change has
    % something concrete to diff against:
    %
    %   parameter              planted   posterior mean   90% CI
    %   sigma_p  gro_sds[1]      3.0          3.0         [2.5 , 3.7 ]
    %   sigma_i  sig_trl         1.2          1.2         [0.86, 1.6 ]
    %   sigma_pi sig_trlxid      1.0          1.0         [0.96, 1.1 ]
    %   sigma_po sig_occxid      0.8          0.91        [0.76, 1.1 ]
    %   sigma_io sig_trlxocc     0.6          0.59        [0.47, 0.74]
    %   sigma_o  sig_occ         0.9          0.63        [0.02, 2.2 ]  <- see (b)
    %   mu       Intercept       5.0          4.8         [3.5 , 6.0 ]
    %
    % All R_hat = 1.00. Sampler notes worth knowing before reading a future
    % failure as a code regression: ~0.4% divergent transitions, treedepth
    % pegged at 9 (511 leapfrog steps per draw; Stan's default ceiling is 10, so
    % not saturating), and 27 sporadic lkj_corr_cholesky_lpdf boundary
    % exceptions. The posterior geometry here is genuinely hard -- crossed
    % variance components plus a per-person residual scale -- and none of these
    % indicate bias on their own.

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

        function testCrossedFacetsAndPerSubjectResidualRecoverSeparately(testCase)
            % Known truth for a crossed p x t x o design with a per-subject
            % residual:
            %   meas_pto = mu + alpha_p + beta_t + gamma_o
            %            + (ab)_pt + (ag)_po + (bg)_to + e_pto
            % with e_pto ~ N(0, sigPOE(p)^2), i.e. heteroscedastic BY PERSON.
            truth = struct('sigP',3.0, 'sigI',1.2, 'sigO',0.9, ...
                'sigPI',1.0, 'sigPO',0.8, 'sigIO',0.6, 'mu',5);
            nsub = 40; ntrl = 20; nocc = 3;

            rng(24680);
            alpha = truth.sigP  * randn(nsub, 1);
            beta  = truth.sigI  * randn(ntrl, 1);
            gamma = truth.sigO  * randn(nocc, 1);
            ab    = truth.sigPI * randn(nsub, ntrl);
            ag    = truth.sigPO * randn(nsub, nocc);
            bg    = truth.sigIO * randn(ntrl, nocc);
            % per-subject residual SDs in [0.8, 2.6]
            sigPOE = 0.8 + 1.8 * rand(nsub, 1);

            dataTbl = localSimulate(truth.mu, alpha, beta, gamma, ab, ag, bg, sigPOE);
            rel = localRunTrtSserr(testCase, dataTbl, 'trtss_recov');

            % ---- (a) the five crossed facet SDs are estimated and finite ----
            facets = {'sig_trl','sigI'; 'sig_occ','sigO'; ...
                      'sig_trlxid','sigPI'; 'sig_occxid','sigPO'; ...
                      'sig_trlxocc','sigIO'};
            for k = 1:size(facets,1)
                fld = facets{k,1};
                testCase.verifyTrue(isfield(rel.out, fld), ...
                    sprintf('rel.out.%s missing', fld));
                testCase.verifyTrue(all(isfinite(rel.out.(fld)(:))), ...
                    sprintf('rel.out.%s has non-finite draws', fld));
            end

            % ---- (b) each facet SD recovers its planted value -----------------
            % Wide-but-real bands: these are variance components from a single
            % finite sample, so the bar is "recovers", not "matches to 2 d.p.".
            %
            % sig_occ is EXCLUDED from the point-estimate band and checked by
            % interval coverage instead. It is the occasion main effect, and this
            % design has only nocc = 3 occasions: a variance component estimated
            % from three levels is weakly identified, so its posterior is diffuse
            % and its mean is shrunk toward zero. Observed on the first live run:
            % planted 0.9, posterior mean 0.63, 90% CI [0.02, 2.2] -- the widest
            % interval and the lowest ESS in the model. Demanding a tight point
            % estimate there would be asking for more than three levels can
            % support, and would make this test seed-fragile. Coverage is the
            % statistically appropriate check; the shrinkage itself is correct
            % behavior, not a defect.
            for k = 1:size(facets,1)
                fld = facets{k,1};
                if strcmp(fld,'sig_occ'); continue; end
                want = truth.(facets{k,2});
                got  = mean(rel.out.(fld)(:));
                testCase.verifyGreaterThan(got, 0.5*want, sprintf( ...
                    '%s = %.3f is far below the planted %.3f', fld, got, want));
                testCase.verifyLessThan(got, 2.0*want, sprintf( ...
                    '%s = %.3f is far above the planted %.3f', fld, got, want));
            end

            % sig_occ: the 90% credible interval must COVER the planted value.
            occDraws = rel.out.sig_occ(:);
            occLo = quantile(occDraws, 0.05);
            occHi = quantile(occDraws, 0.95);
            testCase.verifyLessThanOrEqual(occLo, truth.sigO, sprintf( ...
                'sig_occ 90%% CI [%.3f, %.3f] excludes the planted %.3f from below', ...
                occLo, occHi, truth.sigO));
            testCase.verifyGreaterThanOrEqual(occHi, truth.sigO, sprintf( ...
                'sig_occ 90%% CI [%.3f, %.3f] excludes the planted %.3f from above', ...
                occLo, occHi, truth.sigO));

            % ---- (c) person SD recovers -------------------------------------
            gro = cell2mat(rel.out.gro_sds{1});   % [draws x 2]
            sigPhat = mean(gro(:,1));
            testCase.verifyGreaterThan(sigPhat, 0.5*truth.sigP);
            testCase.verifyLessThan(sigPhat, 2.0*truth.sigP);

            % ---- (d) THE SEPARATION CHECK -----------------------------------
            % The per-subject residual must track sigPOE(p) and must NOT have
            % absorbed the crossed facet effects. If it had, every subject's
            % residual would be inflated toward the lumped value
            % sqrt(sigPOE^2 + sigPI^2 + sigPO^2 + sigIO^2).
            ind_sdlog = cell2mat(rel.out.ind_sdlog{1});   % [draws x nsub]
            pop_sdlog = rel.out.pop_sdlog(:,1);
            resid_sd  = mean(exp(pop_sdlog + ind_sdlog), 1)';   % nsub x 1

            relTarget  = mean(sigPOE);
            lumpTarget = mean(sqrt(sigPOE.^2 + truth.sigPI^2 + ...
                truth.sigPO^2 + truth.sigIO^2));
            testCase.verifyLessThan(mean(resid_sd), 0.5*(relTarget + lumpTarget), ...
                sprintf(['the per-subject residual (%.3f) sits closer to the ' ...
                'LUMPED target (%.3f) than to the true residual (%.3f): the ' ...
                'crossed facets are being absorbed'], mean(resid_sd), ...
                lumpTarget, relTarget));

            % and it must track the per-person truth, not just its average
            % (base corrcoef, so no Statistics Toolbox dependency)
            cc = corrcoef(resid_sd, sigPOE);
            testCase.verifyGreaterThan(cc(1,2), 0.4, ...
                'recovered per-subject residuals do not track the planted ones');

            % ---- (e) the design tag is what the viewer routes on -------------
            testCase.verifyEqual(rel.analysis, 'trt_sserrvar');
        end

    end
end

function dataTbl = localSimulate(mu, alpha, beta, gamma, ab, ag, bg, sigPOE)
% Balanced crossed p x t x o long table in the toolbox I/O contract
% (id, meas, time), one row per single-trial score.
nsub = numel(alpha); ntrl = numel(beta); nocc = numel(gamma);
n = nsub*ntrl*nocc;
ids  = cell(n,1);
tims = cell(n,1);
meas = zeros(n,1);
r = 0;
for pp = 1:nsub
    for oo = 1:nocc
        for tt = 1:ntrl
            r = r + 1;
            ids{r}  = sprintf('S%02d', pp);
            tims{r} = sprintf('t%d', oo);
            meas(r) = mu + alpha(pp) + beta(tt) + gamma(oo) ...
                + ab(pp,tt) + ag(pp,oo) + bg(tt,oo) ...
                + sigPOE(pp)*randn;
        end
    end
end
dataTbl = table(ids, meas, tims, 'VariableNames', {'id','meas','time'});
end

function rel = localRunTrtSserr(testCase, dataTbl, tag)
% Subject-level error variance (sserrvar = 2) WITH an occasion facet, i.e.
% analysis 25. Fixed seed, per the reproducibility contract.
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'sserrvar', 2, ...
    'chains', 2, ...
    'warmup', 500, ...
    'sampling', 500, ...
    'seed', 12345, ...
    'verbose', 1, ...
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
% Remove the per-run output directory and any stray CmdStan build/output files.
p = fullfile(baseDir, 'trtss_recov');
if isfolder(p)
    try
        rmdir(p, 's');
    catch
        warning('tests:cleanupDirFailed', ...
            ['Could not remove recovery artifact directory ''%s''. ', ...
            'How to fix: close any process using it and delete it manually.'], p);
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
end
