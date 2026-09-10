classdef TestGammaDiffTrtRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square
    % two-facet (test-retest) two-event difference-score design (analysis 10,
    % family='gamma', non-concurrent).
    %
    % Two conditions (cells) are simulated from the crossed Persons x Occasions x
    % Trials scaled-chi-square model. Persons, occasions, and the person x occasion
    % interaction are SHARED across conditions (correlated log-scale effects), so
    % their cross-cell covariances are non-zero; the trial, person x trial, and
    % trial x occasion effects are condition-specific (independent), matching the
    % non-concurrent policy that zeros the trial-indexed cross-covariances. The
    % known log-scale parameters are mapped to observed-scale difference-score
    % variance/covariance components with the same helpers the estimator uses
    % (psyrat_gamma_varcomps_trt + psyrat_gamma_crosscov_trt) and the true
    % generalizability difference-score coefficients (equivalence, stability,
    % equivalence-and-stability) are computed in closed form via psyrat_diffrel_trt.
    % The test fits the model through psyrat_computevarcomp and confirms each true
    % coefficient lands inside the recovered 95% credible interval.
    %
    % The GENERALIZABILITY coefficients are asserted: they depend on the person,
    % person x occasion, person x trial, and residual components but not on the
    % occasion or trial MAIN effects, so they are well identified with only a few
    % occasions.
    %
    % Skipped cleanly (Incomplete, not Failed) when CmdStan/MatlabStan are
    % unavailable, like the other CmdStan integration tests.

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
        function testDiffGeneralizabilityCoefficientsRecovered(testCase)
            % ---- Known truth (log expected-score scale), two conditions f/s ----
            aF = log(5.5); aS = log(5.0);
            nu = 18;                 % population df (integer for the chi-square draw)
            % per-condition log-scale SDs: [person occ trial p*o p*t t*o]
            sF = struct('p',0.36,'o',0.22,'t',0.16,'po',0.13,'pt',0.11,'ot',0.09);
            sS = struct('p',0.32,'o',0.20,'t',0.14,'po',0.12,'pt',0.10,'ot',0.08);
            % cross-condition log-scale correlations. Shared facets (person,
            % occasion, person x occasion) are correlated; the condition-specific
            % (non-concurrent) trial-indexed facets are independent (0).
            cor = struct('p',0.6,'o',0.5,'po',0.45,'t',0,'pt',0,'ot',0);

            nsub = 35; nocc = 4; ntrl = 8;

            dataTbl = localSimulateGammaDiffTrt(nsub, nocc, ntrl, aF, aS, nu, sF, sS, cor);

            % ---- Closed-form observed-scale truth (per cell + cross-cell) ----
            vcF = psyrat_gamma_varcomps_trt(aF, sF.p, sF.o, sF.t, sF.po, sF.pt, sF.ot, nu);
            vcS = psyrat_gamma_varcomps_trt(aS, sS.p, sS.o, sS.t, sS.po, sS.pt, sS.ot, nu);
            cc  = psyrat_gamma_crosscov_trt(aF, sF.p, sF.o, sF.t, sF.po, sF.pt, sF.ot, ...
                aS, sS.p, sS.o, sS.t, sS.po, sS.pt, sS.ot, ...
                cor.p, cor.o, cor.t, cor.po, cor.pt, cor.ot);

            % Truth varcov arrays (non-concurrent policy: person/occ/person x occ
            % off-diagonals kept, trial-indexed zeroed), matching the extractor.
            truth = struct();
            truth.bp  = localVc(vcF.sigma_p2,  vcS.sigma_p2,  cc.cov_p);   % person
            truth.bo  = localVc(vcF.sigma_o2,  vcS.sigma_o2,  cc.cov_o);   % occasion
            truth.bpo = localVc(vcF.sigma_po2, vcS.sigma_po2, cc.cov_po);  % person x occasion
            truth.bt  = localVc(vcF.sigma_t2,  vcS.sigma_t2,  0);          % trial
            truth.bpi = localVc(vcF.sigma_pt2, vcS.sigma_pt2, 0);          % person x trial
            truth.boi = localVc(vcF.sigma_ot2, vcS.sigma_ot2, 0);          % trial x occasion
            truth.er  = [0.5*log(vcF.sigma_res2), 0.5*log(vcS.sigma_res2)];

            % ---- Fit ----
            rel = localRunDiffTrtGamma(testCase, dataTbl, 'recov_gamma_diff_trt');

            testCase.verifyEqual(rel.analysis, 'trt_diff');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % Reconstruct the estimated observed-scale (co)variance draws.
            est = struct();
            est.bp  = localCell2Arr(rel.out.id_varcov{1});
            est.bo  = localCell2Arr(rel.out.occ_varcov{1});
            est.bpo = localCell2Arr(rel.out.oid_varcov{1});
            est.bt  = localCell2Arr(rel.out.trl_varcov{1});
            est.bpi = localCell2Arr(rel.out.tid_varcov{1});
            est.boi = localCell2Arr(rel.out.to_varcov{1});
            est.er  = localCell2Arr(rel.out.b_sigma{1});

            for f = {'bp','bo','bpo','bt','bpi','boi','er'}
                testCase.verifyTrue(all(isfinite(est.(f{1})(:))), ...
                    sprintf('Non-finite draws in %s', f{1}));
            end

            % Headline recovery: each true generalizability difference-score
            % coefficient must fall inside the recovered 95% CI.
            for reltype = 1:3
                trueDS = localDiffCoeff(truth, reltype, ntrl, nocc);
                estDS  = localDiffCoeff(est,   reltype, ntrl, nocc);
                testCase.verifyGreaterThanOrEqual(trueDS.pt, estDS.ll, sprintf( ...
                    ['True difference-score G (reltype %d) %.3f below the ', ...
                    'recovered 95%% CI lower limit %.3f.'], reltype, trueDS.pt, estDS.ll));
                testCase.verifyLessThanOrEqual(trueDS.pt, estDS.ul, sprintf( ...
                    ['True difference-score G (reltype %d) %.3f above the ', ...
                    'recovered 95%% CI upper limit %.3f.'], reltype, trueDS.pt, estDS.ul));
            end
        end
    end
end

function ds = localDiffCoeff(v, reltype, ntrl, nocc)
% Difference-score generalizability coefficient via psyrat_diffrel_trt, using the
% same field->argument mapping as psyrat_relsummary (bp/bpi/bpo/bt/bo/boi/er_var,
% er_cov=0). Scalar (truth) inputs give a point coefficient; draw inputs give a CI.
ds = psyrat_diffrel_trt('bp',v.bp,'bpi',v.bpi,'bpo',v.bpo, ...
    'bt',v.bt,'bo',v.bo,'boi',v.boi,'er_var',v.er,'er_cov',0, ...
    'obs',[ntrl ntrl],'nocc',[nocc nocc],'reltype',reltype,'est','gen','CI',0.95);
end

function a = localVc(var_f, var_s, cov_fs)
% Build a [1 x 2 x 2] cross-condition (co)variance array for scalar truth inputs.
a = zeros(1,2,2);
a(1,1,1) = var_f; a(1,2,2) = var_s;
a(1,1,2) = cov_fs; a(1,2,1) = cov_fs;
end

function a = localCell2Arr(c)
% Reconstruct a numeric array from the {num2cell(...)} storage used in REL.out.
a = cell2mat(c);
end

function tbl = localSimulateGammaDiffTrt(nsub, nocc, ntrl, aF, aS, nu, sF, sS, cor)
% Balanced two-condition crossed Persons x Occasions x Trials scaled-chi-square
% long table. Person, occasion, and person x occasion effects are shared across
% conditions (correlated log-scale effects); trial, person x trial, and
% trial x occasion effects are condition-specific (independent) - the
% non-concurrent structure. Base-MATLAB randn only; chi-square(nu) is a sum of nu
% (integer) squared standard normals.
rng(12345);

% Correlated (shared) effects across conditions.
[uP_f,  uP_s]  = local_pair(nsub,        sF.p,  sS.p,  cor.p);
[uO_f,  uO_s]  = local_pair(nocc,        sF.o,  sS.o,  cor.o);
[uPO_f, uPO_s] = local_pair(nsub*nocc,   sF.po, sS.po, cor.po);
uPO_f = reshape(uPO_f, nsub, nocc); uPO_s = reshape(uPO_s, nsub, nocc);

% Condition-specific (independent) trial-indexed effects.
uT_f  = sF.t  * randn(ntrl,1);       uT_s  = sS.t  * randn(ntrl,1);
uPT_f = sF.pt * randn(nsub,ntrl);    uPT_s = sS.pt * randn(nsub,ntrl);
uOT_f = sF.ot * randn(ntrl,nocc);    uOT_s = sS.ot * randn(ntrl,nocc);

conds = {'condF','condS'};
nrow = 2 * nsub * nocc * ntrl;
ids  = cell(nrow,1); evs = cell(nrow,1);
time = zeros(nrow,1); meas = zeros(nrow,1);
r = 0;
for c = 1:2
    for p = 1:nsub
        for o = 1:nocc
            for t = 1:ntrl
                r = r + 1;
                ids{r} = sprintf('S%02d', p);
                evs{r} = conds{c};
                time(r) = o;
                if c == 1
                    logmu = aF + uP_f(p) + uO_f(o) + uT_f(t) ...
                        + uPO_f(p,o) + uPT_f(p,t) + uOT_f(t,o);
                else
                    logmu = aS + uP_s(p) + uO_s(o) + uT_s(t) ...
                        + uPO_s(p,o) + uPT_s(p,t) + uOT_s(t,o);
                end
                mu = exp(logmu);
                chi2 = sum(randn(nu,1).^2);
                meas(r) = (mu / nu) * chi2;
            end
        end
    end
end
tbl = table(ids, meas, evs, time, 'VariableNames', {'id','meas','event','time'});
end

function [ef, es] = local_pair(n, sf, ss, c)
% Correlated (c ~= 0) or independent bivariate normal effect pair, base randn.
z1 = randn(n,1); z2 = randn(n,1);
ef = sf * z1;
es = ss * (c * z1 + sqrt(1 - c^2) * z2);
end

function rel = localRunDiffTrtGamma(testCase, dataTbl, tag)
% Fit the two-facet difference-score pipeline under the Gamma family, fixed seed.
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
p = fullfile(baseDir, 'recov_gamma_diff_trt');
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
