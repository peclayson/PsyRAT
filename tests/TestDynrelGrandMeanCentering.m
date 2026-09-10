classdef TestDynrelGrandMeanCentering < matlab.unittest.TestCase
    % Pins the maintainer-confirmed estimand for multi-group dynamic
    % reliability: the subject-level dimension is standardized ONCE across the
    % whole sample (grand-mean centering), NOT separately within each group.
    %
    % Confirmed design decision (2026-06-24): groups share one common predictor
    % scale so the per-group conditional surfaces G(z)/D(z) are comparable at the
    % same z, and "z = 0" is the overall-sample mean of the dimension. Scores are
    % deliberately NOT centered within group. A future switch to within-group
    % centering (which would set each group's mean z to ~0 and break cross-group
    % comparability) must fail this test.
    %
    % Deterministic and CmdStan-free: it exercises psyrat_zscore_subjectlevel
    % directly, so it lives in the unit lane (runAllUnitTests) and always runs.

    methods (Test)
        function testStandardizationIsGrandMeanNotWithinGroup(testCase)
            % Two groups with deliberately DIFFERENT dimension means, so grand-
            % mean and within-group centering give measurably different z.
            rng(424242, 'twister');
            nA = 30; nB = 30; ntrl = 5;

            % subject-level raw dimension values, offset between groups (e.g. a
            % clinical group scoring higher on a questionnaire than controls).
            rawA = 40 + 8 * randn(nA, 1);
            rawB = 60 + 8 * randn(nB, 1);
            raw  = [rawA; rawB];
            grp  = [repmat({'A'}, nA, 1); repmat({'B'}, nB, 1)];

            % long table: ntrl rows per subject; the dimension is constant within
            % a subject (a between-person covariate). meas is irrelevant to the
            % standardizer and is just filler.
            nsub = nA + nB;
            id    = cell(nsub * ntrl, 1);
            group = cell(nsub * ntrl, 1);
            dim1  = zeros(nsub * ntrl, 1);
            meas  = randn(nsub * ntrl, 1);
            r = 0;
            for s = 1:nsub
                for t = 1:ntrl
                    r = r + 1;
                    id{r}    = sprintf('S%03d', s);
                    group{r} = grp{s};
                    dim1(r)  = raw(s);
                end
            end
            T = table(id, meas, group, dim1, ...
                'VariableNames', {'id', 'meas', 'group', 'dim1'});

            z = psyrat_zscore_subjectlevel(T, 'dim1');

            % collapse to one z per subject (constant within subject)
            [~, ia] = unique(T.id, 'stable');
            zsub = z(ia);
            gsub = T.group(ia);

            % (1) GRAND standardization: across ALL subjects the standardized
            %     dimension has mean 0 and (sample, N-1) SD 1, by construction.
            testCase.verifyEqual(mean(zsub), 0, 'AbsTol', 1e-9, ...
                'Pooled standardized dimension must have overall mean 0.');
            testCase.verifyEqual(std(zsub), 1, 'AbsTol', 1e-9, ...
                'Pooled standardized dimension must have overall SD 1.');

            % (2) The between-group offset is PRESERVED, proving the dimension was
            %     NOT centered within group. Group A (low raw mean) sits clearly
            %     below 0 on the common scale; group B (high raw mean) clearly
            %     above. Within-group centering would force both group means to
            %     ~0, so these checks would then fail.
            mA = mean(zsub(strcmp(gsub, 'A')));
            mB = mean(zsub(strcmp(gsub, 'B')));
            testCase.verifyLessThan(mA, -0.3, ...
                'Low-mean group must have negative mean z under grand centering.');
            testCase.verifyGreaterThan(mB, 0.3, ...
                'High-mean group must have positive mean z under grand centering.');
            testCase.verifyGreaterThan(mB - mA, 0.8, ...
                ['Groups must stay separated on the common z scale (within-group ', ...
                'centering would collapse this separation to ~0).']);
        end
    end
end
