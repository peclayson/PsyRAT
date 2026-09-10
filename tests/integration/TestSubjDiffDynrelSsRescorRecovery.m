classdef TestSubjDiffDynrelSsRescorRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the SUBJECT-LEVEL CONCURRENT dynamic
    % difference reliability models with a PER-SUBJECT residual correlation (Rast &
    % Clayson, Stage 4e): one-facet (analysis 21; dynrel=2 + diffest=2 + sserrvar=2 +
    % diffwpcov=2 + ssrescor=2) and trial + occasion two-facet (analysis 22; + time).
    %
    % This extends the Stage-4d population-rescor design by making the residual
    % correlation PER-SUBJECT: each participant's correlation is drawn from the
    % population distribution rho[p] = tanh(rescor_mu + sd_rescor * z_rescor[p]):
    %   log sigma_e,event(z,p) = a_event + b_sigma_event * z_p + delta_{p,event}
    %   [e1; e2]_p ~ MVN(0, diag([s1,s2]) * [1 rho_p; rho_p 1] * diag([s1,s2]))
    % The model must recover the population location rescor_mu (and average
    % correlation), the between-person spread sd_rescor (identifiability-fragile -
    % asserted with a generous band), the genuine per-subject variation in rho, the
    % difference universe-score variance sigma^2_p_delta, the per-event residual
    % scale slopes, and the per-participant phi (each subject's own z-conditioned SDs
    % with that subject's OWN rho) vs the closed form.
    %
    % Sample size >= 72 (the LKJ/N lesson; per-subject rho needs many subjects AND
    % enough trials/subject to identify the spread). Skips cleanly when
    % CmdStan/MatlabStan are unavailable. Slow (compiles + samples two models).

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
        function testOneFacetPerSubjectRhoRecovers(testCase)
            rng(40404040, 'twister');
            nsub = 80; ntrl = 24;
            Sp = [16 6; 6 9];
            mu = [5 2]; bmean = [1.0 0.5];
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10];
            sds = [0.5 0.4];                    % between-person residual-scale SDs
            rescorMu = atanh(0.4);              % TRUE population location (avg rho ~ 0.4)
            sdRescor = 0.5;                     % TRUE between-person SD of atanh(rho)
            trlSD = [1.0 0.8];
            uniTrue = Sp(1,1) + Sp(2,2) - 2*Sp(1,2);

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            zres = randn(nsub,1);                       % per-person atanh effect
            rhoTrue = tanh(rescorMu + sdRescor*zres);   % TRUE per-subject correlation
            Lp = chol(Sp, 'lower');
            nu_p = (Lp * randn(2, nsub))';
            delta = [randn(nsub,1)*sds(1), randn(nsub,1)*sds(2)];
            trlE = [randn(ntrl,1)*trlSD(1), randn(ntrl,1)*trlSD(2)];

            ids = {}; meas = []; ev = {}; dim1 = []; events = {'gain','loss'};
            for p = 1:nsub
                rp = rhoTrue(p);
                for t = 1:ntrl
                    s1 = exp(a(1) + bsig(1)*zsub(p) + delta(p,1));
                    s2 = exp(a(2) + bsig(2)*zsub(p) + delta(p,2));
                    Lr = [s1 0; rp*s2 sqrt(1-rp^2)*s2];
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

            rel = localRun(testCase, T, 'subjdiffssrho1f');
            testCase.assertEqual(rel.analysis, 'ic_diff_dynrel_sserrvar_rho');

            idvc = cell2mat(rel.out.id_varcov{1});
            bsg  = cell2mat(rel.out.b_sigma{1});
            bsd  = cell2mat(rel.out.b_sigma_dim{1});
            sdid = cell2mat(rel.out.sd_id{1});
            rmu  = cell2mat(rel.out.rescor_mu{1});
            sdr  = cell2mat(rel.out.sd_rescor{1});
            rhoss = cell2mat(rel.out.rho_ss{1});   % draws x NSUB

            % difference universe-score variance var1+var2-2cov. The model recovers
            % the SAMPLE realization of the person location block, not the population
            % parameters: for a finite nsub the empirical difference-of-variances
            % departs from the population value by sampling error, and LKJ(2) on the
            % joint person block additionally shrinks the location covariance (the
            % documented 4a/4b/4d difference-contrast effect, NOT a wiring bug). So
            % compare the recovered value to the EMPIRICAL difference variance of the
            % generated person effects (what the data actually contain).
            Cnu = cov(nu_p);
            empUni = Cnu(1,1) + Cnu(2,2) - 2*Cnu(1,2);
            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, empUni, 'AbsTol', 4.5);
            % population location (atanh scale) and the implied average correlation
            testCase.verifyEqual(mean(rmu), rescorMu, 'AbsTol', 0.30);
            testCase.verifyEqual(mean(mean(rhoss,2)), mean(rhoTrue), 'AbsTol', 0.15);
            % between-person spread is identifiability-fragile: assert it is
            % identified (not collapsed to 0, not exploded) within a generous band
            testCase.verifyTrue(mean(sdr) > 0.10 && mean(sdr) < 1.20, ...
                sprintf('sd_rescor recovered = %.3f outside the broad band [0.10,1.20]', mean(sdr)));
            % the recovered per-subject correlations genuinely vary across subjects
            testCase.verifyTrue(std(mean(rhoss,1)) > 0.02, ...
                'recovered per-subject rho shows no between-subject variation');
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsd(:,1)), bsig(1), 'AbsTol', 0.18);
            testCase.verifyEqual(mean(bsd(:,2)), bsig(2), 'AbsTol', 0.18);
            testCase.verifyEqual(mean(sdid(:,3)), sds(1), 'AbsTol', 0.25);
            testCase.verifyEqual(mean(sdid(:,4)), sds(2), 'AbsTol', 0.25);

            % per-participant table: mean phi vs the closed form (each subject's own
            % z + scale RE, with that subject's OWN true rho)
            summ = psyrat_dynrel_summary(rel, 'CI', .95);
            tab = summ.strata(1).ssrel_table;
            testCase.assertEqual(height(tab), nsub);
            genTrue = zeros(nsub,1);
            for p = 1:nsub
                rp = rhoTrue(p);
                s1 = exp(a(1) + bsig(1)*zsub(p) + delta(p,1));
                s2 = exp(a(2) + bsig(2)*zsub(p) + delta(p,2));
                relerr = (s1^2 + s2^2 - 2*rp*s1*s2)/ntrl;
                genTrue(p) = uniTrue / (uniTrue + relerr);
            end
            testCase.verifyEqual(mean(tab.gen_pt), mean(genTrue), 'AbsTol', 0.10);
            testCase.verifyTrue(all(tab.gen_pt > 0 & tab.gen_pt < 1));
            testCase.verifyTrue(all(tab.dep_pt <= tab.gen_pt + 1e-9));
        end

        function testTwoFacetPerSubjectRhoRecovers(testCase)
            rng(50505050, 'twister');
            nsub = 72; nocc = 4; ntrl = 6;
            Sp = [16 6; 6 9];
            mu = [5 2]; bmean = [1.0 0.5];
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10];
            sds = [0.5 0.4];
            rescorMu = atanh(0.4); sdRescor = 0.5;
            occSD = [1.0 0.8]; trlSD = [1.0 0.8];
            ptSD  = [1.0 0.9]; poSD  = [0.9 0.8]; toSD = [0.5 0.4];
            uniTrue = Sp(1,1) + Sp(2,2) - 2*Sp(1,2);
            bpiTrue = ptSD(1)^2 + ptSD(2)^2;
            bpoTrue = poSD(1)^2 + poSD(2)^2;

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            zres = randn(nsub,1);
            rhoTrue = tanh(rescorMu + sdRescor*zres);
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
                rp = rhoTrue(p);
                for o = 1:nocc
                    for t = 1:ntrl
                        s1 = exp(a(1) + bsig(1)*zsub(p) + delta(p,1));
                        s2 = exp(a(2) + bsig(2)*zsub(p) + delta(p,2));
                        Lr = [s1 0; rp*s2 sqrt(1-rp^2)*s2];
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

            rel = localRun(testCase, T, 'subjdiffssrho2f');
            testCase.assertEqual(rel.analysis, 'ic_diff_dynrel_sserrvar_trt_rho');

            idvc = cell2mat(rel.out.id_varcov{1});
            bsg  = cell2mat(rel.out.b_sigma{1});
            sdid = cell2mat(rel.out.sd_id{1});
            rmu  = cell2mat(rel.out.rescor_mu{1});
            sdr  = cell2mat(rel.out.sd_rescor{1});
            rhoss = cell2mat(rel.out.rho_ss{1});

            % difference universe-score variance var1+var2-2cov. The model recovers
            % the SAMPLE realization of the person location block, not the population
            % parameters: for a finite nsub the empirical difference-of-variances
            % departs from the population value by sampling error, and LKJ(2) on the
            % joint person block additionally shrinks the location covariance (the
            % documented 4a/4b/4d difference-contrast effect, NOT a wiring bug). So
            % compare the recovered value to the EMPIRICAL difference variance of the
            % generated person effects (what the data actually contain).
            Cnu = cov(nu_p);
            empUni = Cnu(1,1) + Cnu(2,2) - 2*Cnu(1,2);
            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, empUni, 'AbsTol', 4.5);
            testCase.verifyEqual(mean(rmu), rescorMu, 'AbsTol', 0.30);
            testCase.verifyEqual(mean(mean(rhoss,2)), mean(rhoTrue), 'AbsTol', 0.15);
            testCase.verifyTrue(mean(sdr) > 0.10 && mean(sdr) < 1.20, ...
                sprintf('sd_rescor recovered = %.3f outside the broad band [0.10,1.20]', mean(sdr)));
            testCase.verifyTrue(std(mean(rhoss,1)) > 0.02, ...
                'recovered per-subject rho shows no between-subject variation');
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.25);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.25);
            testCase.verifyEqual(mean(sdid(:,3)), sds(1), 'AbsTol', 0.25);
            testCase.verifyEqual(mean(sdid(:,4)), sds(2), 'AbsTol', 0.25);

            % per-participant table: mean phi vs closed form (reltype 3, observed
            % occasions; each subject's own z + scale RE + own true rho)
            summ = psyrat_dynrel_summary(rel, 'CI', .95, 'reltype', 3, ...
                'nocc', 'observed');
            tab = summ.strata(1).ssrel_table;
            testCase.assertEqual(height(tab), nsub);
            genTrue = zeros(nsub,1);
            for p = 1:nsub
                rp = rhoTrue(p);
                s1 = exp(a(1) + bsig(1)*zsub(p) + delta(p,1));
                s2 = exp(a(2) + bsig(2)*zsub(p) + delta(p,2));
                res = s1^2 + s2^2 - 2*rp*s1*s2;
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
    'ssrescor', 2, ...
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
for tag = {'subjdiffssrho1f','subjdiffssrho2f'}
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
