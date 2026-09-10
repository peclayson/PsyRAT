classdef TestSsrelDiffSemDecisionType < PsyRATTestBase
    % Regression coverage for finding F6(e): the per-subject difference-score
    % SEM in psyrat_ssrel_diff must use the error variance that matches the
    % requested decision type. With est='gen' (relative/generalizability) the
    % SEM is sqrt(rel_err); with est='dep' (absolute/dependability) it is
    % sqrt(abs_err). Both calls share identical synthetic inputs so the only
    % difference is 'est'.
    %
    % A single MCMC draw is used so sem_pt collapses to sqrt(max(err,0))
    % exactly (no Jensen gap from averaging sqrt over draws), which lets the
    % SEM^2 == err equalities be checked to tolerance.

    methods (Test)
        function testGenAndDepSemUseMatchingErrorVariance(testCase)
            [bp, bt, erVar, wpCov, idtable] = localSingleDrawInputs();

            outDep = psyrat_ssrel_diff( ...
                'bp', bp, 'bt', bt, 'er_var', erVar, ...
                'wp_cov', wpCov, 'idtable', idtable, 'CI', 0.95, ...
                'est', 'dep');

            outGen = psyrat_ssrel_diff( ...
                'bp', bp, 'bt', bt, 'er_var', erVar, ...
                'wp_cov', wpCov, 'idtable', idtable, 'CI', 0.95, ...
                'est', 'gen');

            wp1 = exp(erVar(1, 1))^2;
            wp2 = exp(erVar(1, 2))^2;
            bt1 = bt(1, 1, 1);
            bt2 = bt(1, 2, 2);
            btCov = bt(1, 2, 1);

            for r = 1:height(idtable)
                obs = [idtable.trls1(r) idtable.trls2(r)];
                hm = psyrat_harmmean(obs);

                relErr = wp1 / obs(1) + wp2 / obs(2) ...
                    - (2 * wpCov) / hm;
                absErr = relErr + bt1 / obs(1) + bt2 / obs(2) ...
                    - (2 * btCov) / hm;

                % Guard the test's own assumptions: both error variances are
                % positive (so max(.,0) does not clamp) and the absolute error
                % strictly exceeds the relative error for these inputs.
                testCase.assertGreaterThan(relErr, 0);
                testCase.assertGreaterThan(absErr, relErr);

                % gen SEM^2 == rel_err, dep SEM^2 == abs_err
                testCase.verifyEqual(outGen.sem_pt(r)^2, relErr, 'AbsTol', 1e-12);
                testCase.verifyEqual(outDep.sem_pt(r)^2, absErr, 'AbsTol', 1e-12);

                % gen SEM <= dep SEM (relative error never exceeds absolute)
                testCase.verifyLessThanOrEqual(outGen.sem_pt(r), outDep.sem_pt(r));
            end
        end
    end
end

function [bp, bt, erVar, wpCov, idtable] = localSingleDrawInputs()
% Deterministic single-draw difference-score inputs. One MCMC draw keeps the
% posterior-mean SEM equal to the per-draw SEM so the SEM^2 == err equalities
% are exact. Two subjects with asymmetric trial counts exercise the
% per-subject loop and the harmonic-mean cross term.
bp = zeros(1, 2, 2);
bp(1, 1, 1) = 0.60;
bp(1, 2, 2) = 0.50;
bp(1, 2, 1) = 0.10;
bp(1, 1, 2) = bp(1, 2, 1);

bt = zeros(1, 2, 2);
bt(1, 1, 1) = 0.25;
bt(1, 2, 2) = 0.20;
bt(1, 2, 1) = 0.04;
bt(1, 1, 2) = bt(1, 2, 1);

% log residual SDs so that exp(er)^2 recovers the residual variances
erVar = log(sqrt([0.50 0.40]));

wpCov = 0.05;

idtable = table({'s1'; 's2'}, [8; 12], [10; 9], ...
    'VariableNames', {'id', 'trls1', 'trls2'});
idtable.trls = min(idtable.trls1, idtable.trls2);
end
