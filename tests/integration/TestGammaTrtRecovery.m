classdef TestGammaTrtRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square
    % test-retest design (analysis 5, family='gamma').
    %
    % Data are simulated from the crossed Persons x Occasions x Trials model the
    % family assumes: a LOG-LINKED mean carrying six independent normal random
    % effects (person, occasion, trial, and the three two-way interactions) and a
    % scaled-chi-square observation Y = (mu/nu) * chi2(nu). The known log-scale
    % SDs are mapped to observed-score-scale variance components by the same
    % conversion the estimator uses (psyrat_gamma_varcomps_trt), and the true
    % generalizability coefficients (coefficients of equivalence, stability, and
    % equivalence-and-stability) are computed in closed form via psyrat_rel_trt.
    % The test fits the model through psyrat_computevarcomp and confirms each true
    % coefficient lands inside the recovered 95% credible interval.
    %
    % The three GENERALIZABILITY coefficients are asserted (not the dependability
    % coefficients): they depend on the person, person x occasion, person x trial,
    % and residual components but NOT on the occasion or trial MAIN effects, so
    % they are well identified even with only a handful of occasions. The occasion
    % main effect is weakly identified at NOCC = 5 (the same few-levels limit the
    % Gaussian test-retest carries), so its dependability coefficients are left
    % out of the strict checks.
    %
    % Like the other CmdStan integration tests, the dependency setup uses
    % assumptions, so the test is skipped cleanly (reported Incomplete, not Failed)
    % when CmdStan/MatlabStan are unavailable (e.g. on CI runners).

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
        function testGeneralizabilityCoefficientsRecovered(testCase)
            % ---- Known truth (log expected-score scale) ----
            alpha = log(5.5);   % log-mean grand intercept
            nu    = 18;         % population scaled-chi-square df (integer for the
                                % toolbox-free chi-square draw below)
            sP  = 0.35;  % person
            sO  = 0.22;  % occasion
            sT  = 0.15;  % trial
            sPO = 0.12;  % person x occasion
            sPT = 0.10;  % person x trial
            sOT = 0.08;  % trial x occasion

            nsub = 40; nocc = 5; ntrl = 12;

            [dataTbl] = localSimulateGammaTrt(nsub, nocc, ntrl, alpha, nu, ...
                sP, sO, sT, sPO, sPT, sOT);

            % ---- Closed-form observed-scale truth and true coefficients ----
            vc = psyrat_gamma_varcomps_trt(alpha, sP, sO, sT, sPO, sPT, sOT, nu);
            truthSD = struct( ...
                'bp',  sqrt(vc.sigma_p2),  'bo',  sqrt(vc.sigma_o2), ...
                'bt',  sqrt(vc.sigma_t2),  'txp', sqrt(vc.sigma_pt2), ...
                'oxp', sqrt(vc.sigma_po2), 'txo', sqrt(vc.sigma_ot2), ...
                'err', sqrt(vc.sigma_res2));

            % ---- Fit ----
            rel = localRunTrtGamma(testCase, dataTbl, 'recov_gamma_trt');

            % Provenance: the gamma test-retest path actually ran under CmdStan.
            testCase.verifyEqual(rel.analysis, 'trt');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % All converted SD-equivalent components are present and finite.
            for f = {'sig_id','sig_occ','sig_trl','sig_trlxid','sig_occxid', ...
                    'sig_trlxocc','sig_err'}
                testCase.verifyTrue(isfield(rel.out, f{1}), ...
                    sprintf('Missing observed-scale component %s', f{1}));
                testCase.verifyTrue(all(isfinite(rel.out.(f{1})(:))), ...
                    sprintf('Non-finite draws in %s', f{1}));
            end

            % Person component sanity: the well-identified observed-scale person
            % SD-equivalent should sit near its closed-form truth (loose bounds to
            % tolerate CmdStan/toolchain drift and finite-sample realization).
            sigIdHat = mean(rel.out.sig_id(:));
            testCase.verifyGreaterThan(sigIdHat, 0.55 * truthSD.bp);
            testCase.verifyLessThan(sigIdHat,    1.60 * truthSD.bp);

            % Headline recovery: each true generalizability coefficient must fall
            % inside the recovered 95% CI, at the design's trial/occasion counts.
            for reltype = 1:3
                trueG = localCoeff(truthSD, 2, reltype, ntrl, nocc);
                [ll, ~, ul] = psyrat_rel_trt('gcoeff', 2, 'reltype', reltype, ...
                    'bp',  rel.out.sig_id(:),      'bo',  rel.out.sig_occ(:), ...
                    'bt',  rel.out.sig_trl(:),     'txp', rel.out.sig_trlxid(:), ...
                    'oxp', rel.out.sig_occxid(:),  'txo', rel.out.sig_trlxocc(:), ...
                    'err', rel.out.sig_err(:), 'obs', ntrl, 'nocc', nocc, 'CI', 0.95);
                testCase.verifyGreaterThanOrEqual(trueG, ll, sprintf( ...
                    ['True generalizability coefficient (reltype %d) %.3f below ', ...
                    'the recovered 95%% CI lower limit %.3f.'], reltype, trueG, ll));
                testCase.verifyLessThanOrEqual(trueG, ul, sprintf( ...
                    ['True generalizability coefficient (reltype %d) %.3f above ', ...
                    'the recovered 95%% CI upper limit %.3f.'], reltype, trueG, ul));
            end
        end
    end
end

function tbl = localSimulateGammaTrt(nsub, nocc, ntrl, alpha, nu, ...
    sP, sO, sT, sPO, sPT, sOT)
% Build a balanced crossed Persons x Occasions x Trials long table from a
% log-linked mean with six independent normal random effects and a scaled-
% chi-square observation. Trial index is nominal and crossed across occasions
% (trial t shares its main effect u_t and person x trial effect u_pt across
% occasions), matching the test-retest model's facet structure. Uses only base
% MATLAB randn; the chi-square draw is the sum of nu (integer) squared standard
% normals, so no Statistics Toolbox is required.
rng(12345);
u_p  = sP  * randn(nsub, 1);
u_o  = sO  * randn(nocc, 1);
u_t  = sT  * randn(ntrl, 1);
u_po = sPO * randn(nsub, nocc);
u_pt = sPT * randn(nsub, ntrl);
u_ot = sOT * randn(ntrl, nocc);

nrow = nsub * nocc * ntrl;
ids  = cell(nrow, 1);
time = zeros(nrow, 1);
meas = zeros(nrow, 1);
r = 0;
for p = 1:nsub
    for o = 1:nocc
        for t = 1:ntrl
            r = r + 1;
            ids{r}  = sprintf('S%02d', p);
            time(r) = o;
            logmu = alpha + u_p(p) + u_o(o) + u_t(t) ...
                + u_po(p, o) + u_pt(p, t) + u_ot(t, o);
            mu = exp(logmu);
            chi2 = sum(randn(nu, 1).^2);   % chi-square(nu), nu integer
            meas(r) = (mu / nu) * chi2;     % scaled chi-square: E = mu, Var = 2 mu^2 / nu
        end
    end
end
tbl = table(ids, meas, time, 'VariableNames', {'id', 'meas', 'time'});
end

function g = localCoeff(sd, gcoeff, reltype, ntrl, nocc)
% Closed-form coefficient at scalar (true) SD-equivalents: psyrat_rel_trt returns
% the single reliability value as its point estimate when the inputs are scalars.
[~, g, ~] = psyrat_rel_trt('gcoeff', gcoeff, 'reltype', reltype, ...
    'bp', sd.bp, 'bo', sd.bo, 'bt', sd.bt, 'txp', sd.txp, ...
    'oxp', sd.oxp, 'txo', sd.txo, 'err', sd.err, ...
    'obs', ntrl, 'nocc', nocc, 'CI', 0.95);
end

function rel = localRunTrtGamma(testCase, dataTbl, tag)
% Run the test-retest pipeline under the Gamma family with a fixed seed.
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'family', 'gamma', ...
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
% Remove the per-run output directory and any stray CmdStan build/output files.
p = fullfile(baseDir, 'recov_gamma_trt');
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
