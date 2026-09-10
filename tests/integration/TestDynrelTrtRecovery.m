classdef TestDynrelTrtRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the trial + occasion two-facet
    % dynamic/conditional reliability model (Rast & Clayson, analysis 14;
    % dynrel = 2 with a time/occasion facet).
    %
    % A balanced person x occasion x trial design is simulated with a fully
    % crossed location structure (person, occasion, trial, person x trial,
    % person x occasion, trial x occasion) and a dimension-conditioned,
    % person-specific log-residual
    %   log sigma_pto,e = alpha + beta_sigma * z + delta_p
    % plus a subject-level standardized dimension z. With one observation per
    % (person, occasion, trial) cell the residual is the highest-order person x
    % trial x occasion term confounded with pure error, exactly what the model
    % estimates. The model must recover the well-identified components
    % (sigma_p, the person x trial / person x occasion SDs, the residual scale
    % alpha and slope beta_sigma, sigma_delta_p) and the dynamic two-facet
    % G(z) surface. The occasion / trial / trial x occasion MAIN effects have few
    % levels (few occasions/trial positions), so they are checked only with broad
    % sanity bounds, and the absolute D(z) surface (which depends on them) is
    % checked with a looser tolerance than the relative G(z) surface.
    %
    % Like the other CmdStan integration tests, the dependency setup uses
    % assumptions, so this skips cleanly (Incomplete, not Failed) when
    % CmdStan/MatlabStan are unavailable. It is slow (compiles + samples) and
    % lives in tests/integration so the unit lane (runAllUnitTests) excludes it.

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
        function testOneDimensionRecoversTruth(testCase)
            % ---- known truth ----
            rng(20260622, 'twister');
            nsub = 50; nocc = 4; ntrl = 6;
            grand = 15;
            sigP = 4; sigDP = 0.25;
            sigOcc = 1.5; sigTrl = 1.5; sigTxP = 1.0; sigOxP = 1.0; sigTxO = 0.8;
            sigLog0 = log(3); bMean = 1.5; bSigma = 0.3;

            zraw = randn(nsub, 1);
            zsub = (zraw - mean(zraw)) / std(zraw);    % standardize (re-std = identity)
            nu_p = sigP * randn(nsub, 1);              % person location
            del_p = sigDP * randn(nsub, 1);            % person scale
            occ_main = sigOcc * randn(nocc, 1);        % occasion main effect
            trl_main = sigTrl * randn(ntrl, 1);        % trial main effect
            pt = sigTxP * randn(nsub, ntrl);           % person x trial
            po = sigOxP * randn(nsub, nocc);           % person x occasion
            to = sigTxO * randn(ntrl, nocc);           % trial x occasion

            n = nsub * nocc * ntrl;
            ids = cell(n, 1); tm = cell(n, 1);
            meas = zeros(n, 1); dim1 = zeros(n, 1);
            r = 0;
            for pp = 1:nsub
                for oo = 1:nocc
                    for tt = 1:ntrl
                        r = r + 1;
                        mu = grand + bMean * zsub(pp) + nu_p(pp) ...
                            + occ_main(oo) + trl_main(tt) ...
                            + pt(pp, tt) + po(pp, oo) + to(tt, oo);
                        sg = exp(sigLog0 + bSigma * zsub(pp) + del_p(pp));
                        ids{r} = sprintf('S%02d', pp);
                        tm{r} = sprintf('t%d', oo);
                        meas(r) = mu + sg * randn;
                        dim1(r) = zsub(pp);
                    end
                end
            end
            dataTbl = table(ids, meas, tm, dim1, ...
                'VariableNames', {'id', 'meas', 'time', 'dim1'});

            rel = localRunDynrelTrt(testCase, dataTbl, 'dynreltrt_1dim');

            % ---- posterior means ----
            gro = cell2mat(rel.out.gro_sds{1});        % draws x 2 (sigma_p, sigma_dp)
            log0 = rel.out.pop_sdlog(:, 1);
            bsig = cell2mat(rel.out.b_sigma{1});       % draws x 1
            sOcc = rel.out.sig_occ(:, 1);
            sTrl = rel.out.sig_trl(:, 1);
            sTxP = rel.out.sig_trlxid(:, 1);
            sOxP = rel.out.sig_occxid(:, 1);
            sTxO = rel.out.sig_trlxocc(:, 1);

            % ---- recovery: well-identified components ----
            testCase.verifyEqual(mean(gro(:, 1)), sigP, 'AbsTol', 1.0);
            testCase.verifyEqual(mean(log0), sigLog0, 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsig(:, 1)), bSigma, 'AbsTol', 0.15);
            testCase.verifyEqual(mean(gro(:, 2)), sigDP, 'AbsTol', 0.30);
            testCase.verifyEqual(mean(sTxP), sigTxP, 'AbsTol', 0.6);
            testCase.verifyEqual(mean(sOxP), sigOxP, 'AbsTol', 0.6);

            % ---- few-level facets (occasion / trial / trial x occasion): broad
            % sanity bounds only (4 occasions, 6 trial positions -> weakly
            % identified variance components) ----
            testCase.verifyGreaterThan(mean(sOcc), 0);
            testCase.verifyLessThan(mean(sOcc), sigOcc + 3);
            testCase.verifyGreaterThan(mean(sTrl), 0);
            testCase.verifyLessThan(mean(sTrl), sigTrl + 3);
            testCase.verifyGreaterThan(mean(sTxO), 0);
            testCase.verifyLessThan(mean(sTxO), sigTxO + 3);

            % ---- two-facet G(z)/D(z) recovered vs truth closed form (reltype 3,
            % both facets random) ----
            zg = [-1.5 0 1.5]; nt = ntrl; no = nocc;
            out = psyrat_rel_dynrel_trt('sig_p', gro(:, 1), 'sig_log0', log0, ...
                'b_sigma', bsig(:, 1), 'sig_occ', sOcc, 'sig_trl', sTrl, ...
                'sig_trlxid', sTxP, 'sig_occxid', sOxP, 'sig_trlxocc', sTxO, ...
                'ndim', 1, 'z1', zg, 'obs', nt, 'nocc', no, 'reltype', 3, ...
                'CI', .95);
            for a = 1:numel(zg)
                sige = exp(sigLog0 + bSigma * zg(a));
                uni = sigP^2;
                gerr = sigTxP^2/nt + sigOxP^2/no + sige^2/(nt*no);
                derr = gerr + sigTrl^2/nt + sigOcc^2/no + sigTxO^2/(nt*no);
                Gtr = uni / (uni + gerr);
                Dtr = uni / (uni + derr);
                testCase.verifyEqual(out.G.pt(a), Gtr, 'AbsTol', 0.08);
                testCase.verifyEqual(out.D.pt(a), Dtr, 'AbsTol', 0.12);
            end
        end
    end
end

function rel = localRunDynrelTrt(testCase, dataTbl, tag)
% Run the two-facet dynamic-reliability pipeline (dynrel = 2 + time facet) with a
% fixed seed.
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'dynrel', 2, ...
    'chains', 2, ...
    'warmup', 700, ...
    'sampling', 700, ...
    'seed', 12345, ...
    'verbose', 1, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
testCase.verifyEqual(rel.analysis, 'ic_dynrel_trt');
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

function localCleanupArtifacts(baseDir)
% Remove the per-run output directory and any stray CmdStan build/output files.
p = fullfile(baseDir, 'dynreltrt_1dim');
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
