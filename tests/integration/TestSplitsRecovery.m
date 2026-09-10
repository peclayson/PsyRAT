classdef TestSplitsRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the data-splits pipeline. Covers both
    % the PARALLEL path (splits=2: equal items-per-split, routed through the
    % standard one-facet / test-retest machinery, cases 1/5) and the NONPARALLEL
    % observed-design path (splits=3: unequal items-per-split, cases 23/24 ->
    % psyrat_splits_summary; Rocha Tables 4/5). Each test simulates split-mean
    % data with KNOWN variance components. The NONPARALLEL cases use REAL items-
    % per-split (n_i) variation -- that variation is what identifies the per-row
    % interaction and the n_i-weighted per-item residual (near-equal n_i collapses
    % them, per the increment-2 finding) -- and then confirm (a) the well-powered
    % components are recovered and (b) the end-to-end dependability/generalizability
    % from psyrat_splits_summary match the closed-form coefficients computed from
    % the TRUE components at the observed design. The PARALLEL case holds n_i equal
    % by design (the split mean is the score; n_i plays no role for cases 1/5) and
    % checks routing plus recovery of the standard one-facet components.
    %
    % Like the other CmdStan integration tests, the dependency setup uses
    % assumptions, so the test is skipped cleanly (reported Incomplete, not
    % Failed) when CmdStan/MatlabStan are unavailable (e.g. on CI runners). It is
    % slow (compiles + samples three models).

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

        function testParallelSingleRecoversTruth(testCase)
            % splits=2 parallel one-facet: items-per-split (n_i) are equal, so the
            % split mean is the score and the run must route through the STANDARD
            % one-facet model (case 1, analysis 'ic') -- NOT the nonparallel
            % ic_splits path -- carrying REL.splits=2 for downstream labeling.
            % Confirms the parallel path fits end-to-end under CmdStan and recovers
            % the well-powered person/residual components. The reliability calc the
            % viewer then runs is the standard one-facet path: validated by the
            % existing one-facet suite and proven byte-identical to a single-trial
            % run by TestSplitsParallel, so it is not re-derived here.
            [tbl, tr] = localSimParallelSingle(2025);
            rel = localRunSplitsParallel(testCase, tbl, 'recov_splits_parallel');

            testCase.verifyEqual(rel.analysis, 'ic');
            testCase.verifyEqual(rel.splits, 2);

            sp = mean(rel.out.sig_u(:));
            se = mean(rel.out.sig_e(:));
            testCase.verifyGreaterThanOrEqual(sp, tr.sigP * 0.7);
            testCase.verifyLessThanOrEqual(sp,  tr.sigP * 1.3);
            testCase.verifyGreaterThanOrEqual(se, tr.sigE * 0.75);
            testCase.verifyLessThanOrEqual(se,  tr.sigE * 1.25);

            % The split main effect (sig_trl = sigma_i) must be recovered, not
            % collapsed into the residual: a regression that dropped the split
            % facet would drive sig_trl toward zero and leak its variance into
            % sig_e. The band is deliberately wide -- the split SD is identified
            % by only ksplit=12 levels, so it is weakly powered (like the
            % nonparallel person x split component, checked at [0.5, 1.8]x) -- but
            % the lower bound still excludes a collapse.
            st = mean(rel.out.sig_trl(:));
            testCase.verifyGreaterThanOrEqual(st, tr.sigS * 0.4);
            testCase.verifyLessThanOrEqual(st,  tr.sigS * 2.2);
        end

        function testSingleOccasionRecoversTruth(testCase)
            % p x (i:s), single occasion (Rocha Table 4).
            [tbl, tr] = localSimSingle(2024);
            rel = localRunSplits(testCase, tbl, 'recov_splits_ic');
            testCase.verifyEqual(rel.analysis, 'ic_splits');

            % well-powered components recovered (person, per-item residual,
            % person x split). Loose multiplicative bounds tolerate toolchain
            % drift while catching gross breakage.
            sp  = mean(rel.out.sig_u(:));
            se  = mean(rel.out.sig_err(:));
            sps = mean(rel.out.sig_splitxid(:));
            testCase.verifyGreaterThanOrEqual(sp, tr.sigP * 0.7);
            testCase.verifyLessThanOrEqual(sp,  tr.sigP * 1.3);
            testCase.verifyGreaterThanOrEqual(se, tr.sigE * 0.75);
            testCase.verifyLessThanOrEqual(se,  tr.sigE * 1.25);
            testCase.verifyGreaterThanOrEqual(sps, tr.sigPS * 0.5);
            testCase.verifyLessThanOrEqual(sps,  tr.sigPS * 1.8);

            % end-to-end: summary D/G vs the TRUE closed-form coefficients at the
            % observed design (n_s = ksplit, n_i = harmonic mean of weights).
            summary = psyrat_splits_summary(rel, 'CI', .95);
            nbar = numel(tbl.weight) / sum(1 ./ tbl.weight);
            [~, Dtrue, ~] = PsyRATAccuracyOracle.splitsSingle(1, tr.sigP, ...
                tr.sigS, tr.sigPS, tr.sigE, tr.ksplit, nbar, .95);
            [~, Gtrue, ~] = PsyRATAccuracyOracle.splitsSingle(2, tr.sigP, ...
                tr.sigS, tr.sigPS, tr.sigE, tr.ksplit, nbar, .95);
            c = summary.strata(1).coef(1);
            testCase.verifyEqual(summary.strata(1).nsplit, tr.ksplit);
            testCase.verifyEqual(c.D(2), Dtrue, 'AbsTol', 0.06);
            testCase.verifyEqual(c.G(2), Gtrue, 'AbsTol', 0.06);
            testCase.verifyLessThanOrEqual(c.D(2), c.G(2) + 1e-9);
        end

        function testMultiOccasionRecoversTruth(testCase)
            % p x o x (i:s), multi occasion (Rocha Table 5).
            [tbl, tr] = localSimTrt(909);
            rel = localRunSplits(testCase, tbl, 'recov_splits_trt');
            testCase.verifyEqual(rel.analysis, 'trt_splits');

            sp = mean(rel.out.sig_id(:));
            se = mean(rel.out.sig_err(:));
            testCase.verifyGreaterThanOrEqual(sp, tr.sigP * 0.7);
            testCase.verifyLessThanOrEqual(sp,  tr.sigP * 1.3);
            testCase.verifyGreaterThanOrEqual(se, tr.sigE * 0.75);
            testCase.verifyLessThanOrEqual(se,  tr.sigE * 1.25);

            % end-to-end: summary D/G vs TRUE closed form (CE and CES). The
            % occasion-involving facets are weakly identified with two occasions,
            % so the coefficient tolerance is wider than the single-occasion case.
            summary = psyrat_splits_summary(rel, 'CI', .95);
            nbar = numel(tbl.weight) / sum(1 ./ tbl.weight);
            testCase.verifyEqual(summary.strata(1).nsplit, tr.ksplit);
            testCase.verifyEqual(summary.strata(1).nocc, tr.nocc);
            for r = [1 3]
                [~, Dtrue, ~] = PsyRATAccuracyOracle.splitsTrt(1, r, tr.sigP, ...
                    tr.sigO, tr.sigS, tr.sigPS, tr.sigPO, tr.sigOS, tr.sigPOS, ...
                    tr.sigE, tr.ksplit, tr.nocc, nbar, .95);
                [~, Gtrue, ~] = PsyRATAccuracyOracle.splitsTrt(2, r, tr.sigP, ...
                    tr.sigO, tr.sigS, tr.sigPS, tr.sigPO, tr.sigOS, tr.sigPOS, ...
                    tr.sigE, tr.ksplit, tr.nocc, nbar, .95);
                c = summary.strata(1).coef(r);
                testCase.verifyEqual(c.D(2), Dtrue, 'AbsTol', 0.10);
                testCase.verifyEqual(c.G(2), Gtrue, 'AbsTol', 0.10);
                testCase.verifyLessThanOrEqual(c.D(2), c.G(2) + 1e-9);
            end
        end

    end
end

% =============================== local helpers ===============================

function [tbl, tr] = localSimParallelSingle(seed)
% Balanced PARALLEL one-facet split-mean design with EQUAL items-per-split (n_i).
% Equal n_i is what makes the splits parallel: the split mean is the score and
% splits=2 routes to the standard one-facet model (case 1), which partitions the
% split-mean scores into person (sig_u), split main effect (sig_trl), and a
% scalar residual (sig_e). Truth is a clean person / split-main / split-mean-
% residual decomposition; the per-item count is held constant so it plays no
% role (unlike the nonparallel design, where n_i variation identifies the
% weighted residual).
rng(seed);
nsub = 80; ksplit = 12; mu = 5; ni = 30;
tr.sigP = 2.5; tr.sigS = 0.8; tr.sigE = 1.2; tr.ksplit = ksplit; tr.ni = ni;
u = tr.sigP * randn(nsub,1);
t = tr.sigS * randn(ksplit,1);
nrow = nsub*ksplit;
ids = cell(nrow,1); meas = zeros(nrow,1); weight = zeros(nrow,1);
r = 0;
for p = 1:nsub
    for s = 1:ksplit
        r = r + 1;
        ids{r} = sprintf('S%02d',p);
        weight(r) = ni;              % EQUAL n_i across all splits (parallel)
        ebar = tr.sigE * randn;      % split-mean residual (scalar sig_e)
        meas(r) = mu + u(p) + t(s) + ebar;
    end
end
tbl = table(ids, meas, weight, 'VariableNames', {'id','meas','weight'});
end

function rel = localRunSplitsParallel(testCase, tbl, tag)
% Run the PARALLEL-splits estimation (splits = 2) with a fixed seed. Parallel
% splits route through the standard one-facet machinery (case 1), so this is a
% normal internal-consistency fit on split-mean scores with the required (equal)
% weight column present.
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', tbl, ...
    'splits', 2, ...
    'chains', 4, ...
    'warmup', 1000, ...
    'sampling', 1000, ...
    'seed', 12345, ...
    'verbose', 1, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
end

function [tbl, tr] = localSimSingle(seed)
% Balanced p x (i:s) split-mean design with a known partition and real n_i
% variation. Split levels are shared across persons by position (the model's
% within-person split index), so the split main effect is identified.
rng(seed);
nsub = 80; ksplit = 12; mu = 5;
tr.sigP = 2.5; tr.sigS = 1.0; tr.sigPS = 1.2; tr.sigE = 4.0; tr.ksplit = ksplit;
u = tr.sigP * randn(nsub,1);
t = tr.sigS * randn(ksplit,1);
nrow = nsub*ksplit;
ids = cell(nrow,1); meas = zeros(nrow,1); weight = zeros(nrow,1);
r = 0;
for p = 1:nsub
    for s = 1:ksplit
        r = r + 1;
        ids{r} = sprintf('S%02d',p);
        ni = randi([8 100]);              % real n_i variation (wide => identifies sig_err)
        weight(r) = ni;
        ps = tr.sigPS * randn;            % person x split
        ebar = (tr.sigE / sqrt(ni)) * randn;   % split-mean residual = sig_e/sqrt(n_i)
        meas(r) = mu + u(p) + t(s) + ps + ebar;
    end
end
tbl = table(ids, meas, weight, 'VariableNames', {'id','meas','weight'});
end

function [tbl, tr] = localSimTrt(seed)
% Balanced p x o x (i:s) split-mean design with a known partition and real n_i
% variation.
rng(seed);
% Five occasions, not the realistic test-retest two: the occasion variance
% component (sig_occ) carries only ~1 df with two occasion levels, so it is
% weakly identified and prior-sensitive (under the wide half-Cauchy prior its
% posterior mean is tail-dominated). Five occasions make the TRUTH identifiable
% so the end-to-end calc can be validated; the observed-design summary itself is
% correct at any n_o. This few-level imprecision is inherent to any two-occasion
% random-occasion design (incl. the existing test-retest), not splits-specific;
% see documentation/splits_observed_design_formulas.md section 6a.
nsub = 80; nocc = 5; ksplit = 6; mu = 5;
tr.sigP = 2.5; tr.sigO = 0.6; tr.sigS = 0.8; tr.sigPS = 1.0; tr.sigPO = 0.9;
tr.sigOS = 0.4; tr.sigPOS = 0.7; tr.sigE = 4.0; tr.ksplit = ksplit; tr.nocc = nocc;
u  = tr.sigP  * randn(nsub,1);
oc = tr.sigO  * randn(nocc,1);
t  = tr.sigS  * randn(ksplit,1);
ps = tr.sigPS * randn(nsub,ksplit);
po = tr.sigPO * randn(nsub,nocc);
os = tr.sigOS * randn(nocc,ksplit);
nrow = nsub*nocc*ksplit;
ids = cell(nrow,1); meas = zeros(nrow,1); weight = zeros(nrow,1); tm = cell(nrow,1);
r = 0;
for p = 1:nsub
    for o = 1:nocc
        for s = 1:ksplit
            r = r + 1;
            ids{r} = sprintf('S%02d',p);
            tm{r}  = sprintf('t%d',o);
            ni = randi([8 100]);
            weight(r) = ni;
            pos = tr.sigPOS * randn;       % person x occasion x split
            ebar = (tr.sigE / sqrt(ni)) * randn;
            meas(r) = mu + u(p) + oc(o) + t(s) + ps(p,s) + po(p,o) + ...
                os(o,s) + pos + ebar;
        end
    end
end
tbl = table(ids, meas, weight, tm, ...
    'VariableNames', {'id','meas','weight','time'});
end

function rel = localRunSplits(testCase, tbl, tag)
% Run the nonparallel-splits estimation (splits = 3) with a fixed seed.
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', tbl, ...
    'splits', 3, ...
    'chains', 4, ...
    'warmup', 1000, ...
    'sampling', 1000, ...
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
for sub = {'recov_splits_ic', 'recov_splits_trt', 'recov_splits_parallel'}
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
    rmdir(fullfile(baseDir, dirs(i).name), 's');
end
end
