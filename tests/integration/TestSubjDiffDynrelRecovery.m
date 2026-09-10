classdef TestSubjDiffDynrelRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the SUBJECT-LEVEL dynamic difference
    % reliability model (Rast & Clayson, analysis 13; dynrel=2 + diffest=2 +
    % sserrvar=2).
    %
    % Joint gain/loss data are simulated with correlated person/trial LOCATION
    % effects and a PER-PERSON event-specific residual whose log-SD depends on a
    % subject-level standardized dimension AND a person-specific scale random
    % effect:
    %   log sigma_e,event(z,p) = a_event + b_sigma_event * z_p + delta_{p,event}
    %   delta_{p,.} ~ N(0, diag(sd_sigma)^2)
    % The model (joint 4-coefficient (0 + event | p | ID) on the location and the
    % log-residual) must recover the difference universe-score variance
    % sigma^2_p_delta (from the location id_varcov), the per-event residual
    % intercepts a_event, the event-specific scale slopes b_sigma_dim, and the
    % between-person residual-scale SDs sd_id[3:4]. It must also produce a
    % per-participant reliability table whose mean phi matches the closed form.
    %
    % Skips cleanly (assumeTrue) when CmdStan/MatlabStan are unavailable. Slow
    % (compiles + samples); lives in tests/integration so the unit lane excludes
    % it.

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
        function testSubjectDifferenceRecoversTruth(testCase)
            rng(515151, 'twister');
            nsub = 50; ntrl = 30;
            Sp = [16 6; 6 9];                  % person event var-cov (gain,loss)
            mu = [5 2]; bmean = [1.0 0.5];      % cell means + mean dim slopes
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10]; % resid log-SD intercept + slope
            sds = [0.5 0.4];                   % between-person residual-scale SDs
            st = [1.5 1.2];                    % trial SDs (independent)
            uniTrue = Sp(1,1) + Sp(2,2) - 2*Sp(1,2); % = 13

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            Lp = chol(Sp, 'lower');
            nu_p = (Lp * randn(2, nsub))';
            delta = [randn(nsub,1)*sds(1), randn(nsub,1)*sds(2)]; % person scale RE
            beta_t = [randn(ntrl,1)*st(1), randn(ntrl,1)*st(2)];

            ids = {}; meas = []; ev = {}; dim1 = []; events = {'gain', 'loss'};
            for p = 1:nsub
                for e = 1:2
                    for t = 1:ntrl
                        mn = mu(e) + bmean(e)*zsub(p) + nu_p(p,e) + beta_t(t,e);
                        sg = exp(a(e) + bsig(e)*zsub(p) + delta(p,e));
                        ids{end+1,1} = sprintf('S%02d', p); %#ok<AGROW>
                        ev{end+1,1} = events{e}; %#ok<AGROW>
                        meas(end+1,1) = mn + randn*sg; %#ok<AGROW>
                        dim1(end+1,1) = zsub(p); %#ok<AGROW>
                    end
                end
            end
            T = table(ids, meas, ev, dim1, ...
                'VariableNames', {'id', 'meas', 'event', 'dim1'});

            rel = localRunSubjDiffDynrel(testCase, T, 'subjdiffdynrecov');
            testCase.assertEqual(rel.analysis, 'ic_diff_dynrel_sserrvar');

            idvc = cell2mat(rel.out.id_varcov{1});
            bsg = cell2mat(rel.out.b_sigma{1});
            bsd = cell2mat(rel.out.b_sigma_dim{1});
            sdid = cell2mat(rel.out.sd_id{1});

            % difference universe-score variance (from the location block)
            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', 4.0);

            % per-event residual intercepts and event-specific scale slopes
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsd(:,1)), bsig(1), 'AbsTol', 0.18);
            testCase.verifyEqual(mean(bsd(:,2)), bsig(2), 'AbsTol', 0.18);

            % between-person residual-scale SDs (sd_id columns 3-4)
            testCase.verifyEqual(mean(sdid(:,3)), sds(1), 'AbsTol', 0.22);
            testCase.verifyEqual(mean(sdid(:,4)), sds(2), 'AbsTol', 0.22);

            % per-participant table: mean estimated phi vs the closed-form truth
            summ = psyrat_dynrel_summary(rel, 'CI', .95);
            tab = summ.strata(1).ssrel_table;
            testCase.assertEqual(height(tab), nsub);

            genTrue = zeros(nsub,1);
            for p = 1:nsub
                wp1 = exp(a(1) + bsig(1)*zsub(p) + delta(p,1))^2;
                wp2 = exp(a(2) + bsig(2)*zsub(p) + delta(p,2))^2;
                relerr = wp1/ntrl + wp2/ntrl;
                genTrue(p) = uniTrue / (uniTrue + relerr);
            end
            testCase.verifyEqual(mean(tab.gen_pt), mean(genTrue), 'AbsTol', 0.08);
            % every per-subject coefficient is a proper reliability in (0,1)
            testCase.verifyTrue(all(tab.gen_pt > 0 & tab.gen_pt < 1));
            testCase.verifyTrue(all(tab.dep_pt <= tab.gen_pt + 1e-9));
        end
    end
end

function rel = localRunSubjDiffDynrel(testCase, dataTbl, tag)
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'dynrel', 2, ...
    'diffest', 2, ...
    'sserrvar', 2, ...
    'chains', 2, ...
    'warmup', 700, ...
    'sampling', 700, ...
    'seed', 12345, ...
    'verbose', 1, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
testCase.verifyEqual(rel.analysis, 'ic_diff_dynrel_sserrvar');
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

function localCleanupArtifacts(baseDir)
p = fullfile(baseDir, 'subjdiffdynrecov');
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
