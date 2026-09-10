classdef TestGammaDynrelSserrRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square
    % SUBJECT-LEVEL dynamic/conditional reliability designs: analysis 26
    % (ic_dynrel_sserrvar, one facet) and analysis 27 (ic_dynrel_sserrvar_trt,
    % trial + occasion), both with family='gamma'.
    %
    % WHY THIS EXISTS. The offline tests prove the calculators
    % (TestGammaDynrelSserrReliability) and the summary wiring
    % (TestDynrelSubjNondiffReliability), but neither exercises the ESTIMATOR:
    % that psyrat_computevarcomp initializes the gamma REL.out slots for
    % ind_bs/ind_nu/ssinfo, extracts those per-person effects from the fit, and
    % appends ssinfo per stratum. Before this work a gamma run of analysis 26/27
    % ERRORED on the first stratum ("Reference to non-existent field 'ssinfo'"),
    % because the {end+1} append was guarded on the analysis but not on the
    % family. This test is what pins that plumbing.
    %
    % These are also the FIRST live recovery tests for analyses 26/27 in either
    % family: PSYRAT_CODEX_TODOS records that their Gaussian recovery passed
    % during development but was deliberately not re-run after the salvage.
    %
    % WHAT IS ASSERTED. The planted per-person dispersion ordering must come back:
    % the recovered per-participant conditional coefficient must correlate with
    % the truth computed from the planted (u_p, v_p) at each participant's own z.
    % That is the assertion that fails if ind_bs/ind_nu are misaligned against
    % id2, or if the conditional conversion is wrong. The population surface is
    % checked for finiteness and monotone sanity only - it is the marginal
    % estimand and is covered by TestGammaDynrelRecovery.
    %
    % Dependency setup and the R-hat gate both use assumptions, so a missing
    % CmdStan or a non-converged fit reports Incomplete, never a false PASS.

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

        function testOneFacetSubjectLevelGammaDynrelRecovers(testCase)
            % ---- Known truth (log expected-score scale) ----
            alpha0 = log(5.0);  b    = 0.25;    % log-mean intercept + dimension slope
            lognu0 = log(20);   b_nu = -0.35;   % log-nu intercept + dimension slope
            sP = 0.40;   % person log-mean SD
            sI = 0.25;   % trial main-effect log-mean SD
            sV = 0.30;   % person log-nu SD
            rho = 0.35;  % corr(person log-mean, person log-nu)

            nsub = 45; ntrl = 30;

            [dataTbl, truth] = localSimulateGammaDynrelSserr(nsub, ntrl, ...
                alpha0, b, lognu0, b_nu, sP, sI, sV, rho);

            rel = localRunDynrelSserrGamma(testCase, dataTbl, 'recov_gamma_dynrelss');

            testCase.verifyEqual(rel.analysis, 'ic_dynrel_sserrvar');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');

            [maxRhat, nConvParams] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[gamma dynrel sserr recovery] max R-hat=%.4f over %d params; converged=%d\n', ...
                maxRhat, nConvParams, rel.out.conv.converged);
            testCase.assumeGreaterThan(nConvParams, 0, ...
                'Convergence table empty; R-hat unextractable, so recovery cannot be validated.');
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf(['Gamma ', ...
                'subject-level dynrel fit non-converged (max R-hat=%.4f >= 1.1 / ', ...
                'n_eff gate); recovery not validated in this environment.'], maxRhat));

            % ---- The estimator stored what the subject-level path needs ----
            % (this is precisely what was missing before: a gamma analysis-26 run
            % had no ind_bs/ind_nu/ssinfo at all)
            for f = {'ind_bs','ind_nu','ssinfo'}
                testCase.verifyTrue(isfield(rel.out, f{1}), sprintf( ...
                    ['REL.out.%s missing on a gamma analysis-26 run - the ', ...
                    'subject-level slots are not being initialized/extracted.'], f{1}));
            end
            indbs = cell2mat(rel.out.ind_bs{1});   % draws x NSUB
            indnu = cell2mat(rel.out.ind_nu{1});
            testCase.verifySize(indbs, [size(rel.out.mu,1) nsub]);
            testCase.verifySize(indnu, [size(rel.out.mu,1) nsub]);
            testCase.verifyEqual(height(rel.out.ssinfo{1}), nsub);

            % ---- The per-participant table is produced and is per-participant ----
            summ = psyrat_dynrel_summary(rel, 'CI', 0.95);
            T = summ.strata(1).ssrel_table;
            testCase.verifyNotEmpty(T, ...
                'No per-participant table: analysis 26 degraded to analysis 11.');
            testCase.verifyEqual(height(T), nsub);
            testCase.verifyTrue(all(isfinite(T.gen_pt)));

            % ---- Recovery: the per-participant ordering matches the truth ----
            r = corr_(T.gen_pt, truth.G_true(:));
            fprintf('[gamma dynrel sserr recovery] per-participant G corr = %.3f\n', r);
            testCase.verifyGreaterThan(r, 0.70, sprintf( ...
                ['Recovered per-participant coefficients correlate only %.3f with ', ...
                'the planted truth. Either the conditional conversion is wrong or ', ...
                'ind_bs/ind_nu are misaligned against id2.'], r));

            % SCALE. Asserting on the coefficient itself is nearly useless here:
            % G saturates near 1 at these trial counts, so an upper bound can never
            % fire and a lower bound tolerates a large residual error. Assert on the
            % implied RESIDUAL instead, which is unsaturated and is what a wrong
            % conditional conversion actually gets wrong.
            %
            % Invert the coefficient for each participant: from
            % G = bp / (bp + err/n) we get err/n = bp*(1-G)/G, so the ratio of
            % implied to planted error variance is comparable without needing bp.
            errRatio = localImpliedErrRatio(T.gen_pt, truth.G_true(:));
            fprintf('[gamma dynrel sserr recovery] implied error-variance ratio = %.3f\n', errRatio);
            testCase.verifyGreaterThan(errRatio, 0.60, sprintf( ...
                ['Implied per-participant error variance is %.2fx the planted value ', ...
                '- the conditional residual is too small.'], errRatio));
            testCase.verifyLessThan(errRatio, 1.70, sprintf( ...
                ['Implied per-participant error variance is %.2fx the planted value ', ...
                '- the conditional residual is too large (a person-facet leak into ', ...
                'V_w would show up here).'], errRatio));

            % genuine between-person spread (else it has collapsed to the surface)
            testCase.verifyGreaterThan(std(T.gen_pt), 1e-3, ...
                ['No between-person spread in the per-participant coefficients; ', ...
                'the conditional conversion has degenerated to the marginal one.']);
        end

        function testTwoFacetSubjectLevelGammaDynrelRecovers(testCase)
            % The analysis-27 twin. Smaller design (a crossed two-facet gamma fit
            % is expensive) - the assertion is the same: the plumbing exists and
            % the per-participant ordering is recovered.
            alpha0 = log(5.0);  b    = 0.20;
            lognu0 = log(20);   b_nu = -0.30;
            sP = 0.35; sO = 0.20; sT = 0.15; sPO = 0.12; sPT = 0.10; sOT = 0.08;
            sV = 0.30; rho = 0.35;

            nsub = 30; nocc = 3; ntrl = 10;

            [dataTbl, truth] = localSimulateGammaDynrelSserrTrt(nsub, nocc, ntrl, ...
                alpha0, b, lognu0, b_nu, sP, sO, sT, sPO, sPT, sOT, sV, rho);

            rel = localRunDynrelSserrTrtGamma(testCase, dataTbl, 'recov_gamma_dynrelsstrt');

            testCase.verifyEqual(rel.analysis, 'ic_dynrel_sserrvar_trt');
            testCase.verifyEqual(rel.family, 'gamma');

            [maxRhat, nConvParams] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[gamma dynrel trt sserr recovery] max R-hat=%.4f over %d params; converged=%d\n', ...
                maxRhat, nConvParams, rel.out.conv.converged);
            testCase.assumeGreaterThan(nConvParams, 0, ...
                'Convergence table empty; R-hat unextractable.');
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf(['Gamma ', ...
                'two-facet subject-level dynrel fit non-converged (max R-hat=%.4f); ', ...
                'recovery not validated in this environment.'], maxRhat));

            for f = {'ind_bs','ind_nu','ssinfo'}
                testCase.verifyTrue(isfield(rel.out, f{1}), sprintf( ...
                    'REL.out.%s missing on a gamma analysis-27 run.', f{1}));
            end

            summ = psyrat_dynrel_summary(rel, 'CI', 0.95);
            T = summ.strata(1).ssrel_table;
            testCase.verifyNotEmpty(T, ...
                'No per-participant table: analysis 27 degraded to analysis 14.');
            testCase.verifyEqual(height(T), nsub);
            testCase.verifyTrue(all(isfinite(T.gen_pt)));

            r = corr_(T.gen_pt, truth.G_true(:));
            fprintf('[gamma dynrel trt sserr recovery] per-participant G corr = %.3f\n', r);
            testCase.verifyGreaterThan(r, 0.60, sprintf( ...
                ['Recovered per-participant coefficients correlate only %.3f with ', ...
                'the planted truth.'], r));

            % SCALE, on the unsaturated implied residual rather than the saturated
            % coefficient - the correlation above is rank-invariant and cannot see a
            % constant multiplicative error in V_w.
            errRatio = localImpliedErrRatio(T.gen_pt, truth.G_true(:));
            fprintf('[gamma dynrel trt sserr recovery] implied error-variance ratio = %.3f\n', errRatio);
            testCase.verifyGreaterThan(errRatio, 0.55, sprintf( ...
                'Implied per-participant error variance is %.2fx the planted value.', errRatio));
            testCase.verifyLessThan(errRatio, 1.80, sprintf( ...
                ['Implied per-participant error variance is %.2fx the planted value ', ...
                '- a person-facet leak into V_w would show up here.'], errRatio));

            testCase.verifyGreaterThan(std(T.gen_pt), 1e-3);
        end

    end
end

function [tbl, truth] = localSimulateGammaDynrelSserr(nsub, ntrl, ...
    alpha0, b, lognu0, b_nu, sP, sI, sV, rho)
% One-facet dynamic gamma table. The dimension value z is a BETWEEN-person
% covariate (one per participant, the Rast & Clayson design), carried in a column
% the toolbox standardizes. The person block (u_p, v_p) is correlated; nu_p is
% subject-constant AND z-conditioned via b_nu.
rng(12345);
z1 = randn(nsub,1);
z2 = randn(nsub,1);
u_p = sP * z1;
v_p = sV * (rho * z1 + sqrt(1 - rho^2) * z2);

% raw (unstandardized) dimension values; psyrat_zscore_subjectlevel standardizes
zraw = randn(nsub,1);
zstd = (zraw - mean(zraw)) ./ std(zraw);

u_t = sI * randn(ntrl,1);

nrow = nsub * ntrl;
ids  = cell(nrow,1);
meas = zeros(nrow,1);
dim1 = zeros(nrow,1);
r = 0;
for p = 1:nsub
    lognu_p = lognu0 + b_nu * zstd(p) + v_p(p);
    nuDraw = max(2, round(exp(lognu_p)));
    for t = 1:ntrl
        r = r + 1;
        ids{r}  = sprintf('S%02d', p);
        dim1(r) = zraw(p);
        logmu = alpha0 + b * zstd(p) + u_p(p) + u_t(t);
        mu = exp(logmu);
        meas(r) = (mu / nuDraw) * sum(randn(nuDraw,1).^2);
    end
end
tbl = table(ids, meas, dim1, 'VariableNames', {'id','meas','dim1'});

% ---- closed-form per-participant truth, at each participant's own z ----
G = zeros(nsub,1);
for p = 1:nsub
    alpha_z = alpha0 + b * zstd(p);
    lognu_z = lognu0 + b_nu * zstd(p);
    vc = psyrat_gamma_varcomps(alpha_z, sP, sI, exp(lognu_z), 0, 0);
    sigE2 = 2 * exp(2*(alpha_z + u_p(p)) + 2*sI^2 - (lognu_z + v_p(p)));
    G(p) = vc.sigma_p2 / (vc.sigma_p2 + sigE2 / ntrl);
end
truth.G_true = G;
truth.u_p = u_p; truth.v_p = v_p; truth.zstd = zstd;
end

function [tbl, truth] = localSimulateGammaDynrelSserrTrt(nsub, nocc, ntrl, ...
    alpha0, b, lognu0, b_nu, sP, sO, sT, sPO, sPT, sOT, sV, rho)
% Two-facet (trial + occasion) dynamic gamma table, same construction as the
% one-facet simulator plus the five crossed log-mean facets.
rng(12345);
a1 = randn(nsub,1); a2 = randn(nsub,1);
u_p = sP * a1;
v_p = sV * (rho * a1 + sqrt(1 - rho^2) * a2);

zraw = randn(nsub,1);
zstd = (zraw - mean(zraw)) ./ std(zraw);

u_o  = sO  * randn(nocc,1);
u_t  = sT  * randn(ntrl,1);
u_po = sPO * randn(nsub,nocc);
u_pt = sPT * randn(nsub,ntrl);
u_ot = sOT * randn(ntrl,nocc);

nrow = nsub * nocc * ntrl;
ids  = cell(nrow,1);
time = zeros(nrow,1);
meas = zeros(nrow,1);
dim1 = zeros(nrow,1);
r = 0;
for p = 1:nsub
    lognu_p = lognu0 + b_nu * zstd(p) + v_p(p);
    nuDraw = max(2, round(exp(lognu_p)));
    for o = 1:nocc
        for t = 1:ntrl
            r = r + 1;
            ids{r}  = sprintf('S%02d', p);
            time(r) = o;
            dim1(r) = zraw(p);
            logmu = alpha0 + b * zstd(p) + u_p(p) + u_o(o) + u_t(t) ...
                + u_po(p,o) + u_pt(p,t) + u_ot(t,o);
            mu = exp(logmu);
            meas(r) = (mu / nuDraw) * sum(randn(nuDraw,1).^2);
        end
    end
end
tbl = table(ids, meas, time, dim1, 'VariableNames', {'id','meas','time','dim1'});

% ---- closed-form per-participant truth (nocc = 1, the summary default) ----
V_w = sO^2 + sT^2 + sPO^2 + sPT^2 + sOT^2;
G = zeros(nsub,1);
for p = 1:nsub
    alpha_z = alpha0 + b * zstd(p);
    lognu_z = lognu0 + b_nu * zstd(p);
    vc = psyrat_gamma_varcomps_trt(alpha_z, sP, sO, sT, sPO, sPT, sOT, ...
        exp(lognu_z), 0, 0);
    errSs = sqrt(2 * exp(2*(alpha_z + u_p(p)) + 2*V_w - (lognu_z + v_p(p))));
    [~, G(p), ~] = psyrat_rel_trt('gcoeff',2,'reltype',3, ...
        'bp',sqrt(vc.sigma_p2),'bo',sqrt(vc.sigma_o2),'bt',sqrt(vc.sigma_t2), ...
        'txp',sqrt(vc.sigma_pt2),'oxp',sqrt(vc.sigma_po2), ...
        'txo',sqrt(vc.sigma_ot2),'err',errSs, ...
        'obs',ntrl,'nocc',1,'CI',0.95);
end
truth.G_true = G;
truth.zstd = zstd;
end

function ratio = localImpliedErrRatio(Ghat, Gtrue)
% Ratio of implied to planted per-participant error variance, recovered by
% inverting the coefficient. From G = bp/(bp + err/n), err/n = bp*(1-G)/G, so
% (1-G)/G is proportional to the error variance at fixed bp and n. Taking the mean
% ratio across participants gives a scale check that does NOT saturate as G -> 1,
% unlike a bound on G itself.
oddsHat  = (1 - Ghat(:))  ./ Ghat(:);
oddsTrue = (1 - Gtrue(:)) ./ Gtrue(:);
ratio = mean(oddsHat) / mean(oddsTrue);
end

function r = corr_(a, b)
% Pearson correlation without the Statistics Toolbox.
a = a(:) - mean(a(:));
b = b(:) - mean(b(:));
r = (a' * b) / sqrt((a' * a) * (b' * b));
end

function rel = localRunDynrelSserrGamma(testCase, dataTbl, tag)
% dynrel=2 (a Dimension 1 column) + sserrvar=2 routes to analysis 26.
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
    'dynrel', 2, ...
    'sserrvar', 2, ...
    'chains', 2, ...
    'warmup', 750, ...
    'sampling', 750, ...
    'seed', 12345, ...
    'verbose', 1, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
end

function rel = localRunDynrelSserrTrtGamma(testCase, dataTbl, tag)
% dynrel=2 + a time facet + sserrvar=2 routes to analysis 27.
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
    'dynrel', 2, ...
    'sserrvar', 2, ...
    'chains', 2, ...
    'warmup', 750, ...
    'sampling', 750, ...
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
% Remove the per-run output directories and any stray CmdStan build/output files.
for tag = {'recov_gamma_dynrelss','recov_gamma_dynrelsstrt'}
    p = fullfile(baseDir, tag{1});
    if isfolder(p)
        try
            rmdir(p, 's');
        catch
            warning('tests:cleanupDirFailed', ...
                ['Could not remove recovery artifact directory ''%s''. ', ...
                'How to fix: close any process using it and delete it manually.'], p);
        end
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
