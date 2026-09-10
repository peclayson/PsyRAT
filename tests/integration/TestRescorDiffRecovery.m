classdef TestRescorDiffRecovery < PsyRATTestBase
    % Recovery test for the two-event difference-score residual-covariance
    % ("rescor") model: analysis 7 with diffwpcov=2 / diffrescor=2.
    %
    % This guards the B10 fix. Before B10 the rescor branch crashed in data prep
    % before any Stan was generated (psyrat_pairdiffrows produced wide rows that
    % psyrat_build_diff_pairs could not consume), and the paired-data trial index
    % was globally unique per pair, which confounded the trial facet with the
    % residual. This test simulates concurrent paired data with a KNOWN residual
    % correlation and known variance components, runs the model through CmdStan,
    % and asserts the residual correlation and the variance decomposition are
    % recovered (in particular that residual variance is no longer absorbed by an
    % over-parameterized trial facet).
    %
    % The data are concurrent by construction (each trial yields both event
    % measurements), which is the regime the rescor model assumes (see the S3/S4
    % positional-pairing caveat documented in psyrat_computevarcomp).
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
            outDir = localCmdStanOutputDir(testCase.projectRoot());
            pats = {'cmdstan*.stan', 'cmdstan*.hpp', 'cmdstan18-*', ...
                'output-*.csv', 'temp.data.R'};
            for i = 1:numel(pats)
                files = dir(fullfile(outDir, pats{i}));
                for j = 1:numel(files)
                    if ~files(j).isdir
                        delete(fullfile(outDir, files(j).name));
                    end
                end
            end
            dirs = dir(fullfile(outDir, 'cmdstan_*'));
            dirs = dirs([dirs.isdir]);
            dirs = dirs(~ismember({dirs.name}, {'.', '..'}));
            for i = 1:numel(dirs)
                try
                    rmdir(fullfile(outDir, dirs(i).name), 's');
                catch
                    warning('tests:cleanupDirFailed', ...
                        ['Could not remove CmdStan artifact directory ''%s''.'], ...
                        fullfile(outDir, dirs(i).name));
                end
            end
        end
    end

    methods (Test)
        function testRescorRecoversKnownResidualCorrelation(testCase)
            % Known ground truth.
            truth = struct('rho', 0.4, 'sd_id', [3 3], 'sigma', [2 2], ...
                'sd_trl', [0 0], 'b', [5 3]);

            tbl = localSimulateConcurrentPairs(truth, 28, 18, 20260619);

            outDir = localCmdStanOutputDir(testCase.projectRoot());
            if ~isfolder(outDir)
                mkdir(outDir);
            end

            REL = psyrat_computevarcomp( ...
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
                'diffwpcov', 2);

            testCase.assertTrue(isfield(REL, 'out') && isfield(REL.out, 'rescor'), ...
                'Rescor model did not return REL.out.rescor.');

            rhoHat = mean(cell2mat(REL.out.rescor{1}));
            sdId = mean(cell2mat(REL.out.sd_id{1}), 1);
            sdTrl = mean(cell2mat(REL.out.sd_trl{1}), 1);
            sigmaHat = mean(exp(cell2mat(REL.out.b_sigma{1})), 1);

            % Residual correlation recovered (the headline rescor parameter).
            testCase.verifyGreaterThan(rhoHat, 0.20, ...
                sprintf('Residual correlation under-recovered (rho_hat=%.3f, true=0.40).', rhoHat));
            testCase.verifyLessThan(rhoHat, 0.60, ...
                sprintf('Residual correlation over-recovered (rho_hat=%.3f, true=0.40).', rhoHat));

            % Subject SDs recovered (true 3 per event).
            testCase.verifyGreaterThan(min(sdId), 2.0, ...
                sprintf('Subject SD under-recovered (sd_id=[%.2f %.2f], true=3).', sdId(1), sdId(2)));
            testCase.verifyLessThan(max(sdId), 4.0, ...
                sprintf('Subject SD over-recovered (sd_id=[%.2f %.2f], true=3).', sdId(1), sdId(2)));

            % Residual SDs recovered (true 2 per event). This is the term that
            % was previously absorbed by the confounded trial facet.
            testCase.verifyGreaterThan(min(sigmaHat), 1.4, ...
                sprintf('Residual SD under-recovered (sigma=[%.2f %.2f], true=2).', sigmaHat(1), sigmaHat(2)));
            testCase.verifyLessThan(max(sigmaHat), 2.6, ...
                sprintf('Residual SD over-recovered (sigma=[%.2f %.2f], true=2).', sigmaHat(1), sigmaHat(2)));

            % Trial facet should now be small (true 0), i.e. not absorbing the
            % residual. Lenient upper bound for finite-sample shrinkage.
            testCase.verifyLessThan(max(sdTrl), 1.2, ...
                sprintf('Trial SD inflated (sd_trl=[%.2f %.2f], true=0); trial facet may be confounding with residual.', sdTrl(1), sdTrl(2)));
        end
    end
end

function tbl = localSimulateConcurrentPairs(truth, nSub, nTrial, seed)
% Deterministic concurrent two-event data: each (subject,trial) yields paired
% event-A/event-B measurements with subject random intercepts and bivariate
% residuals correlated at truth.rho. No trial-position effect (truth.sd_trl=0).
rng(seed);
sigma = truth.sigma;
covE = [sigma(1)^2, truth.rho*sigma(1)*sigma(2); ...
        truth.rho*sigma(1)*sigma(2), sigma(2)^2];
cholE = chol(covE, 'lower');

ids = cell(nSub*nTrial*2, 1);
evs = cell(nSub*nTrial*2, 1);
meas = zeros(nSub*nTrial*2, 1);
k = 0;
for s = 1:nSub
    a = truth.sd_id(:) .* randn(2, 1);
    for t = 1:nTrial
        e = cholE * randn(2, 1);
        k = k + 1;
        ids{k} = sprintf('s%02d', s); evs{k} = 'A';
        meas(k) = truth.b(1) + a(1) + e(1);
        k = k + 1;
        ids{k} = sprintf('s%02d', s); evs{k} = 'B';
        meas(k) = truth.b(2) + a(2) + e(2);
    end
end
tbl = table(ids, meas, evs, 'VariableNames', {'id', 'meas', 'event'});
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
end
