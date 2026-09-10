classdef TestSubjDiffDynrelTrtRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the SUBJECT-LEVEL trial + occasion
    % two-facet dynamic difference reliability model (Rast & Clayson, analysis 16 /
    % Stage 4b; dynrel=2 + diffest=2 + sserrvar=2 + a time/occasion facet).
    %
    % This composes the group-level two-facet difference design (analysis 15) with
    % the per-person event-specific scale random effect of the one-facet
    % subject-level model (analysis 13). A balanced person x event x occasion x
    % trial design is simulated with correlated person event LOCATION effects, the
    % event-specific crossed facet effects (occasion / trial / person x trial /
    % person x occasion / trial x occasion), and a per-event residual whose log-SD
    % depends on a subject-level standardized dimension AND a per-person scale RE:
    %   log sigma_pto,e,event(z,p) = a_event + b_sigma_event * z_p + delta_{p,event}
    %   delta_{p,.} ~ N(0, diag(sd_sigma)^2)
    % The model (joint 4-coefficient (0 + event | p | ID) on the location and the
    % log-residual, the other five facets location-only) must recover the
    % difference universe-score variance sigma^2_p_delta (from the location
    % id_varcov), the per-event residual intercepts a_event and scale slopes
    % b_sigma_dim, and the between-person residual-scale SDs sd_id[3:4]. It must
    % also produce a per-participant reliability table whose mean phi matches the
    % closed form, and a typical-person reference surface G_delta(z) matching its
    % closed form.
    %
    % Sample size: the difference universe-score variance is sensitive to the person
    % event CORRELATION, which LKJ(2) regularizes toward 0 under small N (biasing
    % the contrast upward); 72 subjects let the data override that prior (the Stage
    % 4a lesson). Skips cleanly when CmdStan/MatlabStan are unavailable. Slow
    % (compiles + samples); lives in tests/integration so the unit lane excludes it.

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
        function testSubjectTwoFacetDifferenceRecoversTruth(testCase)
            rng(20260623, 'twister');
            nsub = 72; nocc = 4; ntrl = 5;
            Sp = [16 6; 6 9];                  % person event var-cov (gain,loss)
            mu = [5 2]; bmean = [1.0 0.5];      % cell means + mean dim slopes
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10]; % resid log-SD intercept + slope
            sds = [0.5 0.4];                    % between-person residual-scale SDs
            occSD = [1.0 0.8]; trlSD = [1.0 0.8];
            ptSD  = [1.0 0.9]; poSD  = [0.9 0.8];
            toSD  = [0.5 0.4];
            uniTrue = Sp(1,1) + Sp(2,2) - 2*Sp(1,2); % = 13
            bpiTrue = ptSD(1)^2 + ptSD(2)^2;          % person x trial contrast var
            bpoTrue = poSD(1)^2 + poSD(2)^2;          % person x occasion contrast var

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            Lp = chol(Sp, 'lower');
            nu_p = (Lp * randn(2, nsub))';                 % nsub x 2 person location
            delta = [randn(nsub,1)*sds(1), randn(nsub,1)*sds(2)]; % person scale RE
            occE = [randn(nocc,1)*occSD(1), randn(nocc,1)*occSD(2)];
            trlE = [randn(ntrl,1)*trlSD(1), randn(ntrl,1)*trlSD(2)];
            ptE  = cat(3, randn(nsub,ntrl)*ptSD(1), randn(nsub,ntrl)*ptSD(2));
            poE  = cat(3, randn(nsub,nocc)*poSD(1), randn(nsub,nocc)*poSD(2));
            toE  = cat(3, randn(ntrl,nocc)*toSD(1), randn(ntrl,nocc)*toSD(2));

            ids = {}; meas = []; ev = {}; tm = {}; dim1 = []; events = {'gain','loss'};
            for p = 1:nsub
                for e = 1:2
                    for o = 1:nocc
                        for t = 1:ntrl
                            mn = mu(e) + bmean(e)*zsub(p) + nu_p(p,e) ...
                                + occE(o,e) + trlE(t,e) ...
                                + ptE(p,t,e) + poE(p,o,e) + toE(t,o,e);
                            sg = exp(a(e) + bsig(e)*zsub(p) + delta(p,e));
                            ids{end+1,1} = sprintf('S%02d', p); %#ok<AGROW>
                            ev{end+1,1} = events{e};            %#ok<AGROW>
                            tm{end+1,1} = sprintf('t%d', o);    %#ok<AGROW>
                            meas(end+1,1) = mn + randn*sg;      %#ok<AGROW>
                            dim1(end+1,1) = zsub(p);            %#ok<AGROW>
                        end
                    end
                end
            end
            T = table(ids, meas, ev, tm, dim1, ...
                'VariableNames', {'id', 'meas', 'event', 'time', 'dim1'});

            rel = localRunSubjDiffDynrelTrt(testCase, T, 'subjdiffdyntrtrecov');
            testCase.assertEqual(rel.analysis, 'ic_diff_dynrel_sserrvar_trt');

            idvc  = cell2mat(rel.out.id_varcov{1});
            tidvc = cell2mat(rel.out.tid_varcov{1});
            oidvc = cell2mat(rel.out.oid_varcov{1});
            bsg = cell2mat(rel.out.b_sigma{1});
            bsd = cell2mat(rel.out.b_sigma_dim{1});
            sdid = cell2mat(rel.out.sd_id{1});

            % ---- well-identified components ----
            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', 4.5);
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.25);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.25);
            testCase.verifyEqual(mean(bsd(:,1)), bsig(1), 'AbsTol', 0.16);
            testCase.verifyEqual(mean(bsd(:,2)), bsig(2), 'AbsTol', 0.16);

            % ---- the subject-level addition: between-person residual-scale SDs
            % (sd_id columns 3-4), the quantity the 4-coefficient person block adds.
            testCase.verifyEqual(mean(sdid(:,3)), sds(1), 'AbsTol', 0.25);
            testCase.verifyEqual(mean(sdid(:,4)), sds(2), 'AbsTol', 0.25);

            % ---- typical-person reference surface G_delta(z) vs closed form
            % (reltype 3, both facets random; observed occasions) ----
            zg = [-1.5 0 1.5];
            out = psyrat_rel_diffdynrel_trt('id_varcov', idvc, ...
                'trl_varcov', cell2mat(rel.out.trl_varcov{1}), ...
                'occ_varcov', cell2mat(rel.out.occ_varcov{1}), ...
                'tid_varcov', tidvc, 'oid_varcov', oidvc, ...
                'to_varcov', cell2mat(rel.out.to_varcov{1}), ...
                'b_sigma', bsg, 'b_sigma_dim', bsd, 'ndim', 1, 'z1', zg, ...
                'obs', ntrl, 'nocc', nocc, 'reltype', 3, 'CI', .95);
            for k = 1:numel(zg)
                wp1 = exp(a(1) + bsig(1)*zg(k))^2;   % typical person: delta = 0
                wp2 = exp(a(2) + bsig(2)*zg(k))^2;
                res = wp1 + wp2;
                Gtr = uniTrue / (uniTrue + bpiTrue/ntrl + bpoTrue/nocc ...
                    + res/(ntrl*nocc));
                testCase.verifyEqual(out.G.pt(k), Gtr, 'AbsTol', 0.10);
            end

            % ---- per-participant table: mean estimated phi vs the closed-form
            % truth at each subject's own z AND scale RE (observed occasions) ----
            summ = psyrat_dynrel_summary(rel, 'CI', .95, 'reltype', 3, ...
                'nocc', 'observed');
            tab = summ.strata(1).ssrel_table;
            testCase.assertEqual(height(tab), nsub);

            genTrue = zeros(nsub,1);
            for p = 1:nsub
                wp1 = exp(a(1) + bsig(1)*zsub(p) + delta(p,1))^2;
                wp2 = exp(a(2) + bsig(2)*zsub(p) + delta(p,2))^2;
                relerr = bpiTrue/ntrl + bpoTrue/nocc + (wp1+wp2)/(ntrl*nocc);
                genTrue(p) = uniTrue / (uniTrue + relerr);
            end
            testCase.verifyEqual(mean(tab.gen_pt), mean(genTrue), 'AbsTol', 0.10);
            % every per-subject coefficient is a proper reliability in (0,1)
            testCase.verifyTrue(all(tab.gen_pt > 0 & tab.gen_pt < 1));
            testCase.verifyTrue(all(tab.dep_pt <= tab.gen_pt + 1e-9));
        end
    end
end

function rel = localRunSubjDiffDynrelTrt(testCase, dataTbl, tag)
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'dynrel', 2, ...
    'diffest', 2, ...
    'sserrvar', 2, ...
    'chains', 2, ...
    'warmup', 800, ...
    'sampling', 800, ...
    'seed', 12345, ...
    'verbose', 1, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
testCase.verifyEqual(rel.analysis, 'ic_diff_dynrel_sserrvar_trt');
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

function localCleanupArtifacts(baseDir)
p = fullfile(baseDir, 'subjdiffdyntrtrecov');
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
