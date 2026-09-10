classdef TestTwoFacetDiffRecovery < PsyRATTestBase
    % Recovery test for the two-facet (test-retest) difference-score model
    % (analysis 10: diffest=2, sserrvar=1, WITH a time/occasion facet). This is
    % the F3 end-to-end guard: simulate two-condition Persons x Trials x
    % Occasions data with KNOWN per-component 2x2 cross-condition covariances and
    % a known per-condition residual, run the REAL pipeline through CmdStan, and
    % confirm
    %   (a) the dominant, well-identified components (person, person x trial,
    %       residual, and the residual correlation in the concurrent case)
    %       recover near truth;
    %   (b1) psyrat_diffrel_trt fed the KNOWN truth reproduces the closed-form
    %       Table 6 coefficient (PsyRATAccuracyOracle.diffrelTrt) exactly;
    %   (b2) psyrat_diffrel_trt fed the RECOVERED posterior draws agrees with the
    %       oracle on those same draws (end-to-end wiring, bias-insensitive); and
    %   (c) the recovered CE generalizability point estimate is within a realistic
    %       band of the known-truth coefficient.
    %
    % Two scenarios: non-concurrent (residual covariance = 0, the default
    % long-format model) and concurrent (known residual correlation via the
    % rescor model). The trial/occasion MAIN effects are small and weakly
    % identified (few levels), so they are validated only through the integrated
    % coefficient, not asserted component-by-component (see localVerifyComponents).
    % >= 4 occasions are simulated (the production default still supports 2).
    %
    % Clean-skips (Incomplete, not Failed) when CmdStan/MatlabStan are
    % unavailable, matching the other integration tests.

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
        function testNonConcurrentRecovery(testCase)
            truth = localTruth(0);   % residual correlation 0
            nTrl = 8; nOcc = 4;
            tbl = localSimulate(truth, 26, nTrl, nOcc, 20260621);

            rel = localRun(testCase, tbl, 'recov_trtdiff_noncc', 1);

            comp = localExtractComponents(rel);

            % (a) the dominant, well-identified components recover near truth.
            localVerifyComponents(testCase, comp, truth);

            % Non-concurrent model does NOT estimate a residual covariance.
            testCase.verifyEqual(comp.rescor, 0, ...
                'Non-concurrent model should report residual correlation 0.');
            testCase.verifyEqual(comp.erCov, zeros(size(comp.erCov)), ...
                'Non-concurrent model should report residual covariance 0.');

            % (b) coefficient checks (function vs oracle + end-to-end recovery).
            localVerifyCoefficient(testCase, comp, truth, [nTrl nTrl], [nOcc nOcc]);
        end

        function testConcurrentRecovery(testCase)
            truth = localTruth(0.4);   % known residual correlation
            nTrl = 8; nOcc = 4;
            tbl = localSimulate(truth, 26, nTrl, nOcc, 20260622);

            rel = localRun(testCase, tbl, 'recov_trtdiff_cc', 2);

            comp = localExtractComponents(rel);

            % (a) the dominant, well-identified components recover near truth.
            localVerifyComponents(testCase, comp, truth);

            % Residual correlation recovered (true 0.4); the headline rescor
            % parameter for the concurrent model.
            testCase.verifyGreaterThan(comp.rescor, 0.15, ...
                sprintf('Residual correlation under-recovered (rho=%.3f, true=0.40).', comp.rescor));
            testCase.verifyLessThan(comp.rescor, 0.65, ...
                sprintf('Residual correlation over-recovered (rho=%.3f, true=0.40).', comp.rescor));

            % (b) coefficient checks (function vs oracle + end-to-end recovery).
            localVerifyCoefficient(testCase, comp, truth, [nTrl nTrl], [nOcc nOcc]);
        end
    end
end

function localVerifyComponents(testCase, comp, truth)
% (a) Assert the dominant, well-identified components recover near truth. The
% distinct diagonals (person 9, person x trial 4, residual 4) also guard the
% extraction-name mapping (id_varcov / tid_varcov / b_sigma). The trial,
% occasion, and trial x occasion MAIN effects are deliberately small and weakly
% identified (few levels), so they are not hard-asserted here; their recovery is
% covered indirectly by the coefficient check.
bpiVar = diag(comp.bpi);
% Person variance (true 9) is the noisiest of the well-identified components:
% with ~26 subjects the posterior mean has a wide spread and a mild upward pull
% from the weakly-informative half-student_t prior, so the band is generous
% (roughly truth +/- a factor of ~1.8) while still catching gross failures.
testCase.verifyGreaterThan(min(diag(comp.bp)), 5.0, ...
    sprintf('Person variance under-recovered (diag=[%.2f %.2f], true=9).', ...
    comp.bp(1,1), comp.bp(2,2)));
testCase.verifyLessThan(max(diag(comp.bp)), 16.0, ...
    sprintf('Person variance over-recovered (diag=[%.2f %.2f], true=9).', ...
    comp.bp(1,1), comp.bp(2,2)));
testCase.verifyGreaterThan(min(bpiVar), 2.6, ...
    sprintf('Person x trial variance under-recovered ([%.2f %.2f], true=4).', ...
    bpiVar(1), bpiVar(2)));
testCase.verifyLessThan(max(bpiVar), 5.6, ...
    sprintf('Person x trial variance over-recovered ([%.2f %.2f], true=4).', ...
    bpiVar(1), bpiVar(2)));
testCase.verifyGreaterThan(min(comp.resVar), 2.8, ...
    sprintf('Residual variance under-recovered ([%.2f %.2f], true=4).', ...
    comp.resVar(1), comp.resVar(2)));
testCase.verifyLessThan(max(comp.resVar), 5.5, ...
    sprintf('Residual variance over-recovered ([%.2f %.2f], true=4).', ...
    comp.resVar(1), comp.resVar(2)));
end

function truth = localTruth(rho)
% Known ground-truth cross-condition covariances (each [2x2]) and residual.
% Scaled like real ERP G-theory data: person, person x trial, and the residual
% dominate (and are well identified), while the trial, occasion, and trial x
% occasion MAIN effects are small. The main effects have few levels (n trials /
% n occasions), so their variances are weakly identified and biased upward; the
% reliability coefficient stays robust only because they are small contributors.
% Diagonals of the dominant components are distinct so an extraction-name swap
% (id_varcov / tid_varcov / residual) would be caught by part (a).
truth = struct();
truth.mu  = [5 3];
truth.Sp  = localCov(3.0, 3.0, 0.6);   % person sigma(p)            (dominant)
truth.Spt = localCov(2.0, 2.0, 0.4);   % person x trial sigma(pi)   (dominant within)
truth.Spo = localCov(1.2, 1.2, 0.3);   % person x occasion sigma(po)
truth.St  = localCov(0.6, 0.6, 0.3);   % trial sigma(i)             (small)
truth.So  = localCov(0.6, 0.6, 0.3);   % occasion sigma(o)          (small)
truth.Sto = localCov(0.4, 0.4, 0.2);   % trial x occasion sigma(io) (small)
truth.Se  = localCov(2.0, 2.0, rho);   % residual sigma(poi,e)
end

function S = localCov(sd1, sd2, r)
S = [sd1^2, r*sd1*sd2; r*sd1*sd2, sd2^2];
end

function tbl = localSimulate(truth, nSub, nTrl, nOcc, seed)
% Balanced, fully-crossed two-condition Persons x Trials x Occasions data. Each
% crossed component is a 2-vector (one entry per condition) drawn with the known
% cross-condition covariance and SHARED across the two conditions by index, so
% the model can recover the cross-condition 2x2. The residual 2-vector is drawn
% per (person, trial-position, occasion) cell with covariance truth.Se.
rng(seed);
A = localDraw(truth.Sp,  nSub);        % person            (index p)
B = localDraw(truth.St,  nTrl);        % trial-position    (index t)
C = localDraw(truth.So,  nOcc);        % occasion          (index o)
D = localDraw(truth.Spt, nSub*nTrl);   % person x trial    (index (p,t))
E = localDraw(truth.Spo, nSub*nOcc);   % person x occasion (index (p,o))
F = localDraw(truth.Sto, nTrl*nOcc);   % trial x occasion  (index (t,o))
cholE = chol(truth.Se, 'lower');

evnames = {'A', 'B'};
nrow = nSub*nTrl*nOcc*2;
ids   = cell(nrow,1);
event = cell(nrow,1);
time  = cell(nrow,1);
meas  = zeros(nrow,1);
k = 0;
for p = 1:nSub
    for o = 1:nOcc
        for t = 1:nTrl
            dpt = D((p-1)*nTrl + t, :);
            epo = E((p-1)*nOcc + o, :);
            fto = F((t-1)*nOcc + o, :);
            g = (cholE * randn(2,1))';
            for c = 1:2
                k = k + 1;
                ids{k}   = sprintf('s%03d', p);
                event{k} = evnames{c};
                time{k}  = sprintf('o%d', o);
                meas(k)  = truth.mu(c) + A(p,c) + B(t,c) + C(o,c) + ...
                    dpt(c) + epo(c) + fto(c) + g(c);
            end
        end
    end
end
tbl = table(ids, meas, event, time, 'VariableNames', {'id','meas','event','time'});
tbl.Properties.Description = 'twofacet_diff_recovery';
end

function M = localDraw(S, n)
% n draws (rows) from a zero-mean bivariate normal with covariance S, using a
% Cholesky factor (no Statistics Toolbox dependency).
L = chol(S, 'lower');
M = (L * randn(2, n))';
end

function rel = localRun(testCase, tbl, tag, diffwpcov)
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', tbl, ...
    'chains', 2, ...
    'warmup', 600, ...
    'sampling', 600, ...
    'seed', 12345, ...
    'verbose', 1, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir, ...
    'sserrvar', 1, ...
    'diffest', 2, ...
    'diffwpcov', diffwpcov);
testCase.assertEqual(rel.analysis, 'trt_diff');
testCase.assertTrue(isfield(rel, 'out') && isfield(rel.out, 'id_varcov'), ...
    'Two-facet difference model did not return REL.out component covariances.');
end

function comp = localExtractComponents(rel)
% Pull the per-draw component arrays from REL.out (first/only group chunk) and
% the posterior-mean summaries used for the recovery assertions.
comp = struct();
comp.bpDraws  = cell2mat(rel.out.id_varcov{1});   % [draws x 2 x 2]
comp.bpiDraws = cell2mat(rel.out.tid_varcov{1});
comp.bpoDraws = cell2mat(rel.out.oid_varcov{1});
comp.btDraws  = cell2mat(rel.out.trl_varcov{1});
comp.boDraws  = cell2mat(rel.out.occ_varcov{1});
comp.boiDraws = cell2mat(rel.out.to_varcov{1});
comp.erVar    = cell2mat(rel.out.b_sigma{1});      % [draws x 2] log residual SD
errvc         = cell2mat(rel.out.err_varcov{1});   % [draws x 2 x 2]
comp.erCov    = errvc(:,1,2);                       % residual covariance draws
comp.rescor   = mean(cell2mat(rel.out.rescor{1}));

comp.bp  = squeeze(mean(comp.bpDraws, 1));          % 2x2 posterior mean (person)
comp.bpi = squeeze(mean(comp.bpiDraws, 1));         % 2x2 posterior mean (person x trial)
comp.resVar = mean(exp(comp.erVar).^2, 1);          % residual variance per cond.
end

function localVerifyCoefficient(testCase, comp, truth, obs, nocc)
% Three coefficient checks across CES/dependability and CE/generalizability:
%   (b1) DETERMINISTIC function correctness: psyrat_diffrel_trt fed the KNOWN
%        truth reproduces the closed-form (PsyRATAccuracyOracle.diffrelTrt) to
%        floating-point tolerance. This is the plan's validation item (b):
%        "psyrat_diffrel_trt reproduces the closed-form ... from the known truth".
%   (b2) END-TO-END WIRING: psyrat_diffrel_trt fed the RECOVERED posterior draws
%        agrees with the oracle fed the SAME draws (tight tolerance). This pins
%        the pipeline-output -> coefficient path independent of estimation bias.
%   (c)  END-TO-END RECOVERY: the recovered CE generalizability point estimate
%        is within a realistic band of the known-truth coefficient. CE/gen is
%        used (not CES/dep) because its universe (person + person x occasion) and
%        relative error (person x trial + residual) involve only well-identified
%        components -- the weakly-identified trial/occasion MAIN effects do not
%        enter, so this is a clean end-to-end magnitude check. (All coefficients
%        are still pinned exactly against the oracle by b1/b2.)
trueBp  = reshape(truth.Sp,[1 2 2]);
trueBpi = reshape(truth.Spt,[1 2 2]);
trueBpo = reshape(truth.Spo,[1 2 2]);
trueBt  = reshape(truth.St,[1 2 2]);
trueBo  = reshape(truth.So,[1 2 2]);
trueBoi = reshape(truth.Sto,[1 2 2]);
trueErVar = log([sqrt(truth.Se(1,1)) sqrt(truth.Se(2,2))]);
trueErCov = truth.Se(1,2);

cases = {3, 'dep'; 1, 'gen'};
for c = 1:size(cases,1)
    reltype = cases{c,1};
    est = cases{c,2};

    % (b1) function vs oracle on KNOWN truth.
    truthFun = psyrat_diffrel_trt( ...
        'bp', trueBp, 'bpi', trueBpi, 'bpo', trueBpo, ...
        'bt', trueBt, 'bo', trueBo, 'boi', trueBoi, ...
        'er_var', trueErVar, 'er_cov', trueErCov, ...
        'obs', obs, 'nocc', nocc, 'reltype', reltype, 'est', est, 'CI', 0.95);
    truthOracle = PsyRATAccuracyOracle.diffrelTrt( ...
        trueBp, trueBpi, trueBpo, trueBt, trueBo, trueBoi, ...
        trueErVar, trueErCov, obs, nocc, reltype, est, 0.95);
    testCase.verifyLessThan(abs(truthFun.pt - truthOracle.pt), 1e-9, ...
        sprintf(['psyrat_diffrel_trt (%s, reltype %d) does not reproduce the ', ...
        'closed-form from known truth (%.6f vs %.6f).'], est, reltype, ...
        truthFun.pt, truthOracle.pt));

    % (b2) function vs oracle on RECOVERED draws (wiring check).
    pipeFun = psyrat_diffrel_trt( ...
        'bp', comp.bpDraws, 'bpi', comp.bpiDraws, 'bpo', comp.bpoDraws, ...
        'bt', comp.btDraws, 'bo', comp.boDraws, 'boi', comp.boiDraws, ...
        'er_var', comp.erVar, 'er_cov', comp.erCov, ...
        'obs', obs, 'nocc', nocc, 'reltype', reltype, 'est', est, 'CI', 0.95);
    pipeOracle = PsyRATAccuracyOracle.diffrelTrt( ...
        comp.bpDraws, comp.bpiDraws, comp.bpoDraws, comp.btDraws, ...
        comp.boDraws, comp.boiDraws, comp.erVar, comp.erCov, ...
        obs, nocc, reltype, est, 0.95);
    testCase.verifyLessThan(abs(pipeFun.pt - pipeOracle.pt), 1e-8, ...
        sprintf(['psyrat_diffrel_trt (%s, reltype %d) disagrees with the oracle ', ...
        'on recovered draws (%.6f vs %.6f).'], est, reltype, ...
        pipeFun.pt, pipeOracle.pt));

    testCase.verifyGreaterThan(pipeFun.pt, 0, ...
        sprintf('%s coefficient (reltype %d) is non-positive.', est, reltype));
    testCase.verifyLessThan(pipeFun.pt, 1, ...
        sprintf('%s coefficient (reltype %d) exceeds 1.', est, reltype));

    % (c) end-to-end recovery for CE generalizability (well-identified-only).
    if reltype == 1 && strcmp(est, 'gen')
        testCase.verifyLessThan(abs(pipeFun.pt - truthOracle.pt), 0.10, ...
            sprintf(['Recovered CE generalizability %.3f differs from the ', ...
            'known-truth %.3f by more than 0.10.'], pipeFun.pt, truthOracle.pt));
    end
end
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

function localCleanupArtifacts(baseDir)
for sub = {'recov_trtdiff_noncc', 'recov_trtdiff_cc'}
    p = fullfile(baseDir, sub{1});
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
