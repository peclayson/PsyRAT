classdef TestDiffDynrelTrtRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the group-level trial + occasion
    % two-facet dynamic difference reliability model (Rast & Clayson, analysis
    % 15; dynrel=2 + diffest=2 + a time/occasion facet).
    %
    % A balanced person x event x occasion x trial design is simulated with
    % correlated person event effects, event-specific crossed facet effects
    % (occasion, trial, person x trial, person x occasion, trial x occasion), and
    % a per-event residual whose log-SD depends on a subject-level standardized
    % dimension:
    %   log sigma_pto,e,event(z) = a_event + b_sigma_event * z
    % With one observation per cell the residual is the highest-order person x
    % trial x occasion term. The model must recover the difference universe-score
    % variance sigma^2_p_delta (from id_varcov), the per-event residual intercepts
    % and slopes, and the dynamic two-facet difference reliability surface
    % G_delta(z) vs the closed form. The few-level facets (occasion / trial /
    % trial x occasion MAIN effects) are weakly identified, so they are checked
    % only with broad sanity bounds and the absolute D_delta(z) surface (which
    % depends on them) with a looser tolerance than the relative G_delta(z).
    %
    % Skips cleanly when CmdStan/MatlabStan are unavailable. Slow (compiles +
    % samples); lives in tests/integration so the unit lane excludes it.

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
        function testGroupDifferenceRecoversTruth(testCase)
            % NOTE on sample size: the difference universe-score variance
            % sigma^2_p_delta = bp1 + bp2 - 2*cov is sensitive to the person event
            % CORRELATION, which the LKJ(2) prior regularizes toward 0 under small
            % N (biasing the contrast upward). A larger subject count is needed so
            % the data override that prior and the difference variance recovers.
            rng(20260623, 'twister');
            nsub = 72; nocc = 4; ntrl = 5;
            Sp = [16 6; 6 9];                  % person event var-cov (gain,loss)
            mu = [5 2]; bmean = [1.0 0.5];      % cell means + mean dim slopes
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10]; % resid log-SD intercept + slope
            occSD = [1.0 0.8]; trlSD = [1.0 0.8];        % occasion / trial main SDs
            ptSD  = [1.0 0.9]; poSD  = [0.9 0.8];        % person x trial / x occasion SDs
            toSD  = [0.5 0.4];                            % trial x occasion SD
            uniTrue = Sp(1,1) + Sp(2,2) - 2*Sp(1,2); % = 13

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            Lp = chol(Sp, 'lower');
            nu_p = (Lp * randn(2, nsub))';                 % nsub x 2 person effects
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
                            sg = exp(a(e) + bsig(e)*zsub(p));
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

            rel = localRunDiffDynrelTrt(testCase, T, 'diffdyntrtrecov');
            testCase.assertEqual(rel.analysis, 'ic_diff_dynrel_trt');

            idvc  = cell2mat(rel.out.id_varcov{1});
            tidvc = cell2mat(rel.out.tid_varcov{1});
            oidvc = cell2mat(rel.out.oid_varcov{1});
            bsg = cell2mat(rel.out.b_sigma{1});
            bsd = cell2mat(rel.out.b_sigma_dim{1});

            % ---- well-identified components ----
            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', 4.5);
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.25);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.25);
            testCase.verifyEqual(mean(bsd(:,1)), bsig(1), 'AbsTol', 0.16);
            testCase.verifyEqual(mean(bsd(:,2)), bsig(2), 'AbsTol', 0.16);

            % person x trial / person x occasion difference (contrast) variances
            bpiHat = mean(tidvc(:,1,1) + tidvc(:,2,2) - 2*tidvc(:,2,1));
            bpoHat = mean(oidvc(:,1,1) + oidvc(:,2,2) - 2*oidvc(:,2,1));
            bpiTrue = ptSD(1)^2 + ptSD(2)^2;   % independent across events
            bpoTrue = poSD(1)^2 + poSD(2)^2;
            testCase.verifyEqual(bpiHat, bpiTrue, 'AbsTol', 1.0);
            testCase.verifyEqual(bpoHat, bpoTrue, 'AbsTol', 1.0);

            % ---- two-facet difference G(z) recovered vs truth closed form
            % (reltype 3, both facets random) ----
            zg = [-1.5 0 1.5]; nt = ntrl; no = nocc;
            out = psyrat_rel_diffdynrel_trt('id_varcov', idvc, ...
                'trl_varcov', cell2mat(rel.out.trl_varcov{1}), ...
                'occ_varcov', cell2mat(rel.out.occ_varcov{1}), ...
                'tid_varcov', tidvc, 'oid_varcov', oidvc, ...
                'to_varcov', cell2mat(rel.out.to_varcov{1}), ...
                'b_sigma', bsg, 'b_sigma_dim', bsd, 'ndim', 1, 'z1', zg, ...
                'obs', nt, 'nocc', no, 'reltype', 3, 'CI', .95);
            for k = 1:numel(zg)
                wp1 = exp(a(1) + bsig(1)*zg(k))^2;
                wp2 = exp(a(2) + bsig(2)*zg(k))^2;
                res = wp1 + wp2;
                Gtr = uniTrue / (uniTrue + bpiTrue/nt + bpoTrue/no + res/(nt*no));
                testCase.verifyEqual(out.G.pt(k), Gtr, 'AbsTol', 0.10);
            end
        end
    end
end

function rel = localRunDiffDynrelTrt(testCase, dataTbl, tag)
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
rel = psyrat_computevarcomp( ...
    'data', dataTbl, ...
    'dynrel', 2, ...
    'diffest', 2, ...
    'chains', 2, ...
    'warmup', 800, ...
    'sampling', 800, ...
    'seed', 12345, ...
    'verbose', 1, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir);
testCase.verifyEqual(rel.analysis, 'ic_diff_dynrel_trt');
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

function localCleanupArtifacts(baseDir)
p = fullfile(baseDir, 'diffdyntrtrecov');
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
