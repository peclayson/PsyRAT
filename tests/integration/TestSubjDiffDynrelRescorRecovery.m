classdef TestSubjDiffDynrelRescorRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the SUBJECT-LEVEL CONCURRENT (rescor)
    % dynamic difference reliability models (Rast & Clayson, Stage 4d): one-facet
    % (analysis 19; dynrel=2 + diffest=2 + sserrvar=2 + diffwpcov=2) and trial +
    % occasion two-facet (analysis 20; + a time facet).
    %
    % This combines the subject-level non-concurrent design (cases 13/16: a
    % per-person event-specific scale random effect delta_{p,e}) with the concurrent
    % design (a true POPULATION residual correlation rho between the two events):
    %   log sigma_e,event(z,p) = a_event + b_sigma_event * z_p + delta_{p,event}
    %   [e1; e2] ~ MVN(0, diag([s1,s2]) * [1 rho; rho 1] * diag([s1,s2]))
    % The model must recover the population residual correlation rescor, the
    % difference universe-score variance sigma^2_p_delta (from the location
    % id_varcov), the between-person residual-scale SDs sd_id[3:4], the per-event
    % residual scale slopes, and the per-participant phi (each subject's own
    % z-conditioned SDs with the shared rescor) vs the closed form.
    %
    % Sample size >= 70 (the LKJ/N lesson). Skips cleanly when CmdStan/MatlabStan are
    % unavailable. Slow (compiles + samples two models); lives in tests/integration.

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
        function testOneFacetSubjectConcurrentRecovers(testCase)
            rng(60606060, 'twister');
            nsub = 72; ntrl = 20;
            Sp = [16 6; 6 9];
            mu = [5 2]; bmean = [1.0 0.5];
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10];
            sds = [0.5 0.4];                    % between-person residual-scale SDs
            rho = 0.4;                          % TRUE population residual correlation
            trlSD = [1.0 0.8];
            uniTrue = Sp(1,1) + Sp(2,2) - 2*Sp(1,2);

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            Lp = chol(Sp, 'lower');
            nu_p = (Lp * randn(2, nsub))';
            delta = [randn(nsub,1)*sds(1), randn(nsub,1)*sds(2)];
            trlE = [randn(ntrl,1)*trlSD(1), randn(ntrl,1)*trlSD(2)];

            ids = {}; meas = []; ev = {}; dim1 = []; events = {'gain','loss'};
            for p = 1:nsub
                for t = 1:ntrl
                    s1 = exp(a(1) + bsig(1)*zsub(p) + delta(p,1));
                    s2 = exp(a(2) + bsig(2)*zsub(p) + delta(p,2));
                    Lr = [s1 0; rho*s2 sqrt(1-rho^2)*s2];
                    e = Lr * randn(2,1);
                    m1 = mu(1) + bmean(1)*zsub(p) + nu_p(p,1) + trlE(t,1) + e(1);
                    m2 = mu(2) + bmean(2)*zsub(p) + nu_p(p,2) + trlE(t,2) + e(2);
                    ids = [ids; {sprintf('S%02d',p)}; {sprintf('S%02d',p)}]; %#ok<AGROW>
                    ev  = [ev; events(1); events(2)];                        %#ok<AGROW>
                    meas = [meas; m1; m2];                                   %#ok<AGROW>
                    dim1 = [dim1; zsub(p); zsub(p)];                         %#ok<AGROW>
                end
            end
            T = table(ids, meas, ev, dim1, ...
                'VariableNames', {'id','meas','event','dim1'});

            rel = localRun(testCase, T, 'subjdiffrescor1f');
            testCase.assertEqual(rel.analysis, 'ic_diff_dynrel_sserrvar_rescor');

            idvc = cell2mat(rel.out.id_varcov{1});
            bsg  = cell2mat(rel.out.b_sigma{1});
            bsd  = cell2mat(rel.out.b_sigma_dim{1});
            sdid = cell2mat(rel.out.sd_id{1});
            rc   = cell2mat(rel.out.rescor{1});

            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', 4.5);
            testCase.verifyEqual(mean(rc), rho, 'AbsTol', 0.15);
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsd(:,1)), bsig(1), 'AbsTol', 0.18);
            testCase.verifyEqual(mean(bsd(:,2)), bsig(2), 'AbsTol', 0.18);
            testCase.verifyEqual(mean(sdid(:,3)), sds(1), 'AbsTol', 0.25);
            testCase.verifyEqual(mean(sdid(:,4)), sds(2), 'AbsTol', 0.25);

            % per-participant table: mean phi vs the closed form (each subject's own
            % z + scale RE, with the shared rescor)
            summ = psyrat_dynrel_summary(rel, 'CI', .95);
            tab = summ.strata(1).ssrel_table;
            testCase.assertEqual(height(tab), nsub);
            genTrue = zeros(nsub,1);
            for p = 1:nsub
                s1 = exp(a(1) + bsig(1)*zsub(p) + delta(p,1));
                s2 = exp(a(2) + bsig(2)*zsub(p) + delta(p,2));
                relerr = (s1^2 + s2^2 - 2*rho*s1*s2)/ntrl;
                genTrue(p) = uniTrue / (uniTrue + relerr);
            end
            testCase.verifyEqual(mean(tab.gen_pt), mean(genTrue), 'AbsTol', 0.10);
            testCase.verifyTrue(all(tab.gen_pt > 0 & tab.gen_pt < 1));
            testCase.verifyTrue(all(tab.dep_pt <= tab.gen_pt + 1e-9));
        end

        function testTwoFacetSubjectConcurrentRecovers(testCase)
            rng(70707070, 'twister');
            nsub = 72; nocc = 4; ntrl = 5;
            Sp = [16 6; 6 9];
            mu = [5 2]; bmean = [1.0 0.5];
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10];
            sds = [0.5 0.4]; rho = 0.4;
            occSD = [1.0 0.8]; trlSD = [1.0 0.8];
            ptSD  = [1.0 0.9]; poSD  = [0.9 0.8]; toSD = [0.5 0.4];
            uniTrue = Sp(1,1) + Sp(2,2) - 2*Sp(1,2);
            bpiTrue = ptSD(1)^2 + ptSD(2)^2;
            bpoTrue = poSD(1)^2 + poSD(2)^2;

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            Lp = chol(Sp, 'lower');
            nu_p = (Lp * randn(2, nsub))';
            delta = [randn(nsub,1)*sds(1), randn(nsub,1)*sds(2)];
            occE = [randn(nocc,1)*occSD(1), randn(nocc,1)*occSD(2)];
            trlE = [randn(ntrl,1)*trlSD(1), randn(ntrl,1)*trlSD(2)];
            ptE  = cat(3, randn(nsub,ntrl)*ptSD(1), randn(nsub,ntrl)*ptSD(2));
            poE  = cat(3, randn(nsub,nocc)*poSD(1), randn(nsub,nocc)*poSD(2));
            toE  = cat(3, randn(ntrl,nocc)*toSD(1), randn(ntrl,nocc)*toSD(2));

            ids = {}; meas = []; ev = {}; tm = {}; dim1 = []; events = {'gain','loss'};
            for p = 1:nsub
                for o = 1:nocc
                    for t = 1:ntrl
                        s1 = exp(a(1) + bsig(1)*zsub(p) + delta(p,1));
                        s2 = exp(a(2) + bsig(2)*zsub(p) + delta(p,2));
                        Lr = [s1 0; rho*s2 sqrt(1-rho^2)*s2];
                        e = Lr * randn(2,1);
                        base1 = mu(1) + bmean(1)*zsub(p) + nu_p(p,1) + occE(o,1) ...
                            + trlE(t,1) + ptE(p,t,1) + poE(p,o,1) + toE(t,o,1);
                        base2 = mu(2) + bmean(2)*zsub(p) + nu_p(p,2) + occE(o,2) ...
                            + trlE(t,2) + ptE(p,t,2) + poE(p,o,2) + toE(t,o,2);
                        ids = [ids; {sprintf('S%02d',p)}; {sprintf('S%02d',p)}]; %#ok<AGROW>
                        ev  = [ev; events(1); events(2)];                        %#ok<AGROW>
                        tm  = [tm; {sprintf('t%d',o)}; {sprintf('t%d',o)}];       %#ok<AGROW>
                        meas = [meas; base1+e(1); base2+e(2)];                   %#ok<AGROW>
                        dim1 = [dim1; zsub(p); zsub(p)];                         %#ok<AGROW>
                    end
                end
            end
            T = table(ids, meas, ev, tm, dim1, ...
                'VariableNames', {'id','meas','event','time','dim1'});

            rel = localRun(testCase, T, 'subjdiffrescor2f');
            testCase.assertEqual(rel.analysis, 'ic_diff_dynrel_sserrvar_trt_rescor');

            idvc = cell2mat(rel.out.id_varcov{1});
            bsg  = cell2mat(rel.out.b_sigma{1});
            sdid = cell2mat(rel.out.sd_id{1});
            rc   = cell2mat(rel.out.rescor{1});

            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', 4.5);
            testCase.verifyEqual(mean(rc), rho, 'AbsTol', 0.15);
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.25);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.25);
            testCase.verifyEqual(mean(sdid(:,3)), sds(1), 'AbsTol', 0.25);
            testCase.verifyEqual(mean(sdid(:,4)), sds(2), 'AbsTol', 0.25);

            % per-participant table: mean phi vs closed form (reltype 3, observed
            % occasions; each subject's own z + scale RE + shared rescor)
            summ = psyrat_dynrel_summary(rel, 'CI', .95, 'reltype', 3, ...
                'nocc', 'observed');
            tab = summ.strata(1).ssrel_table;
            testCase.assertEqual(height(tab), nsub);
            genTrue = zeros(nsub,1);
            for p = 1:nsub
                s1 = exp(a(1) + bsig(1)*zsub(p) + delta(p,1));
                s2 = exp(a(2) + bsig(2)*zsub(p) + delta(p,2));
                res = s1^2 + s2^2 - 2*rho*s1*s2;
                relerr = bpiTrue/ntrl + bpoTrue/nocc + res/(ntrl*nocc);
                genTrue(p) = uniTrue / (uniTrue + relerr);
            end
            testCase.verifyEqual(mean(tab.gen_pt), mean(genTrue), 'AbsTol', 0.10);
            testCase.verifyTrue(all(tab.gen_pt > 0 & tab.gen_pt < 1));
        end
    end
end

function rel = localRun(testCase, dataTbl, tag)
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'dynrel', 2, ...
    'diffest', 2, ...
    'sserrvar', 2, ...
    'diffwpcov', 2, ...
    'chains', 2, ...
    'warmup', 800, ...
    'sampling', 800, ...
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
for tag = {'subjdiffrescor1f','subjdiffrescor2f'}
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
