classdef TestDiffDynrelRescorRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the CONCURRENT / residual-correlation
    % (rescor) GROUP-LEVEL dynamic difference reliability models (Rast & Clayson,
    % Stage 4c): one-facet (analysis 17; dynrel=2 + diffest=2 + sserrvar=1 +
    % diffwpcov=2) and trial + occasion two-facet (analysis 18; + a time facet).
    %
    % The two events are simulated CONCURRENTLY: at each (person, [occasion,] trial)
    % cell a bivariate residual is drawn with a TRUE residual correlation rho and
    % heteroscedastic per-event SDs that depend on a subject-level standardized
    % dimension:
    %   log sigma_e,event(z,p) = a_event + b_sigma_event * z_p
    %   [e1; e2] ~ MVN(0, diag([s1,s2]) * [1 rho; rho 1] * diag([s1,s2]))
    % The concurrent bivariate model must recover the estimated residual correlation
    % rescor, the difference universe-score variance sigma^2_p_delta (from the
    % location id_varcov), the per-event residual scale intercepts/slopes, and the
    % dynamic difference reliability surface G_delta(z) vs the closed form (which now
    % includes the residual covariance term -2*rho*s1(z)*s2(z)).
    %
    % Sample size >= 70 (the Stage-4a LKJ/N lesson for the difference universe
    % variance). Skips cleanly when CmdStan/MatlabStan are unavailable. Slow
    % (compiles + samples two models); lives in tests/integration.

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
        function testOneFacetConcurrentRecovers(testCase)
            rng(40404040, 'twister');
            nsub = 72; ntrl = 20;
            Sp = [16 6; 6 9];                  % person event var-cov (gain,loss)
            mu = [5 2]; bmean = [1.0 0.5];
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10];
            rho = 0.4;                          % TRUE residual correlation
            trlSD = [1.0 0.8];
            uniTrue = Sp(1,1) + Sp(2,2) - 2*Sp(1,2); % = 13

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            Lp = chol(Sp, 'lower');
            nu_p = (Lp * randn(2, nsub))';
            trlE = [randn(ntrl,1)*trlSD(1), randn(ntrl,1)*trlSD(2)];

            ids = {}; meas = []; ev = {}; dim1 = []; events = {'gain','loss'};
            for p = 1:nsub
                for t = 1:ntrl
                    s1 = exp(a(1) + bsig(1)*zsub(p));
                    s2 = exp(a(2) + bsig(2)*zsub(p));
                    Lr = [s1 0; rho*s2 sqrt(1-rho^2)*s2]; % chol of residual cov
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

            rel = localRun(testCase, T, 'diffdynrescor1f', false);
            testCase.assertEqual(rel.analysis, 'ic_diff_dynrel_rescor');

            idvc = cell2mat(rel.out.id_varcov{1});
            bsg  = cell2mat(rel.out.b_sigma{1});
            bsd  = cell2mat(rel.out.b_sigma_dim{1});
            rc   = cell2mat(rel.out.rescor{1});

            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', 4.5);
            testCase.verifyEqual(mean(rc), rho, 'AbsTol', 0.15);
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.30);
            % interleaved b_sigma_dim: col 1 = event 1, col 2 = event 2
            testCase.verifyEqual(mean(bsd(:,1)), bsig(1), 'AbsTol', 0.18);
            testCase.verifyEqual(mean(bsd(:,2)), bsig(2), 'AbsTol', 0.18);

            % dynamic surface G_delta(z) vs the closed form WITH the residual
            % covariance term (reltype is N/A for one-facet)
            summ = psyrat_dynrel_summary(rel, 'CI', .95, 'obs', ntrl);
            st = summ.strata(1);
            for idx = [1 numel(st.z1)]
                zv = st.z1(idx);
                s1 = exp(a(1) + bsig(1)*zv); s2 = exp(a(2) + bsig(2)*zv);
                res = s1^2 + s2^2 - 2*rho*s1*s2;       % difference residual var
                Gtr = uniTrue / (uniTrue + res/ntrl);
                testCase.verifyEqual(st.G.pt(idx), Gtr, 'AbsTol', 0.10);
            end
            testCase.verifyTrue(all(st.G.pt(:) > 0 & st.G.pt(:) < 1));
        end

        function testTwoFacetConcurrentRecovers(testCase)
            rng(50505050, 'twister');
            nsub = 72; nocc = 4; ntrl = 5;
            Sp = [16 6; 6 9];
            mu = [5 2]; bmean = [1.0 0.5];
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10];
            rho = 0.4;
            occSD = [1.0 0.8]; trlSD = [1.0 0.8];
            ptSD  = [1.0 0.9]; poSD  = [0.9 0.8]; toSD = [0.5 0.4];
            uniTrue = Sp(1,1) + Sp(2,2) - 2*Sp(1,2);

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            Lp = chol(Sp, 'lower');
            nu_p = (Lp * randn(2, nsub))';
            occE = [randn(nocc,1)*occSD(1), randn(nocc,1)*occSD(2)];
            trlE = [randn(ntrl,1)*trlSD(1), randn(ntrl,1)*trlSD(2)];
            ptE  = cat(3, randn(nsub,ntrl)*ptSD(1), randn(nsub,ntrl)*ptSD(2));
            poE  = cat(3, randn(nsub,nocc)*poSD(1), randn(nsub,nocc)*poSD(2));
            toE  = cat(3, randn(ntrl,nocc)*toSD(1), randn(ntrl,nocc)*toSD(2));

            ids = {}; meas = []; ev = {}; tm = {}; dim1 = []; events = {'gain','loss'};
            for p = 1:nsub
                for o = 1:nocc
                    for t = 1:ntrl
                        s1 = exp(a(1) + bsig(1)*zsub(p));
                        s2 = exp(a(2) + bsig(2)*zsub(p));
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

            rel = localRun(testCase, T, 'diffdynrescor2f', true);
            testCase.assertEqual(rel.analysis, 'ic_diff_dynrel_trt_rescor');

            idvc = cell2mat(rel.out.id_varcov{1});
            bsg  = cell2mat(rel.out.b_sigma{1});
            bsd  = cell2mat(rel.out.b_sigma_dim{1});
            rc   = cell2mat(rel.out.rescor{1});

            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', 4.5);
            testCase.verifyEqual(mean(rc), rho, 'AbsTol', 0.15);
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.25);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.25);
            testCase.verifyEqual(mean(bsd(:,1)), bsig(1), 'AbsTol', 0.16);
            testCase.verifyEqual(mean(bsd(:,2)), bsig(2), 'AbsTol', 0.16);

            % two-facet difference surface G_delta(z) vs closed form (reltype 3,
            % both facets random; observed occasions)
            bpiTrue = ptSD(1)^2 + ptSD(2)^2;
            bpoTrue = poSD(1)^2 + poSD(2)^2;
            summ = psyrat_dynrel_summary(rel, 'CI', .95, 'reltype', 3, ...
                'nocc', 'observed');
            st = summ.strata(1);
            for idx = [1 numel(st.z1)]
                zv = st.z1(idx);
                s1 = exp(a(1) + bsig(1)*zv); s2 = exp(a(2) + bsig(2)*zv);
                res = s1^2 + s2^2 - 2*rho*s1*s2;
                Gtr = uniTrue / (uniTrue + bpiTrue/ntrl + bpoTrue/nocc ...
                    + res/(ntrl*nocc));
                testCase.verifyEqual(st.G.pt(idx), Gtr, 'AbsTol', 0.10);
            end
        end
    end
end

function rel = localRun(testCase, dataTbl, tag, includeTime) %#ok<INUSD>
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'dynrel', 2, ...
    'diffest', 2, ...
    'sserrvar', 1, ...
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
for tag = {'diffdynrescor1f','diffdynrescor2f'}
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
