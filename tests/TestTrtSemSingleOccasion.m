classdef TestTrtSemSingleOccasion < PsyRATTestBase
    %Regression coverage for the test-retest (TRT) SEM single-occasion fix.
    %
    %The default TRT estimand is a single-occasion score (n'_o = 1), so the
    %reported G/D coefficients are produced by psyrat_rel_trt with nocc left
    %at its default of 1. The variance-table SEM must be drawn from the SAME
    %error variance as that coefficient, i.e. SEM = sqrt(coefficient error
    %variance). The bug (scientific audit finding F2 / sec.14) recomputed the
    %SEM error term with nocc = numel(time), so with >1 occasion the reported
    %SEM no longer matched the coefficient (empirically about half).
    %
    %These tests construct a TRT result with >1 occasion, run the variance
    %table, and confirm the reported SEM equals the nocc = 1 value (which is
    %sqrt of the reported coefficient's error variance) and is invariant to
    %the number of occasions.

    methods (Test)

        function testTrtSemMatchesSingleOccasionCoefficientErrorVariance(testCase)
            %With 4 occasions, the reported SEM must equal the single-occasion
            %(nocc = 1) value - i.e. sqrt of the reported coefficient's error
            %variance - NOT the nocc = numel(time) value the bug produced.
            obs = 6;
            nocc = 4;
            gcoeff = 1;   %dependability
            reltype = 1;  %coefficient of equivalence
            d = testCase.sampleDraws();

            psyrat_data = testCase.makeTrtData(d, obs, nocc, gcoeff, reltype);

            tbl = psyrat_variancet('psyrat_data', psyrat_data, 'gui', 0);
            reportedSem = testCase.parsePointEstimate(testCase.semString(tbl));

            %Reported SEM should match the single-occasion error variance...
            expectedSemNocc1 = TestTrtSemSingleOccasion.semPointEstimate( ...
                d, obs, 1, gcoeff, reltype);
            testCase.verifyEqual(reportedSem, round(expectedSemNocc1, 2), ...
                'AbsTol', 0.011, ...
                ['Reported TRT SEM should equal sqrt(error variance) of the ' ...
                 'single-occasion (nocc = 1) coefficient.']);

            %...and must NOT match the buggy nocc = numel(time) value, which is
            %materially different (roughly half) for these draws.
            buggySemNocc4 = TestTrtSemSingleOccasion.semPointEstimate( ...
                d, obs, nocc, gcoeff, reltype);
            testCase.verifyGreaterThan( ...
                abs(round(expectedSemNocc1, 2) - round(buggySemNocc4, 2)), 0.5, ...
                ['Test draws must make the nocc = 1 and nocc = numel(time) ' ...
                 'SEMs clearly distinguishable.']);
            testCase.verifyGreaterThan(abs(reportedSem - round(buggySemNocc4, 2)), 0.5, ...
                'Reported SEM must not use nocc = numel(time).');
        end

        function testTrtSemInvariantToNumberOfOccasions(testCase)
            %The SEM must depend only on the (single-occasion) estimand, not on
            %how many occasions were collected. Across every gcoeff/reltype
            %branch, a 1-occasion and a 4-occasion result built from identical
            %variance-component draws must report an identical SEM string.
            obs = 6;
            d = testCase.sampleDraws();

            for gcoeff = 1:2
                for reltype = 1:3
                    singleOcc = testCase.makeTrtData(d, obs, 1, gcoeff, reltype);
                    multiOcc = testCase.makeTrtData(d, obs, 4, gcoeff, reltype);

                    tblSingle = psyrat_variancet('psyrat_data', singleOcc, 'gui', 0);
                    tblMulti = psyrat_variancet('psyrat_data', multiOcc, 'gui', 0);

                    msg = sprintf(['TRT SEM must be invariant to the number of ' ...
                        'occasions (gcoeff=%d, reltype=%d).'], gcoeff, reltype);
                    testCase.verifyEqual(testCase.semString(tblMulti), ...
                        testCase.semString(tblSingle), msg);

                    %And the reported SEM equals the analytic nocc = 1 value for
                    %this branch (sqrt of the coefficient error variance).
                    expected = TestTrtSemSingleOccasion.semPointEstimate( ...
                        d, obs, 1, gcoeff, reltype);
                    reported = testCase.parsePointEstimate(testCase.semString(tblMulti));
                    testCase.verifyEqual(reported, round(expected, 2), ...
                        'AbsTol', 0.011, msg);
                end
            end
        end

    end

    methods (Access = private)

        function d = sampleDraws(~)
            %Posterior draws (standard-deviation units) for one group/event.
            %A large residual relative to the other components makes the
            %nocc-dependent terms dominate, so the nocc = 1 and nocc = 4 SEMs
            %are clearly distinguishable.
            d = struct();
            d.sig_id     = [2.00; 2.10; 1.90; 2.05; 1.95];   %between-person
            d.sig_occ    = [0.30; 0.35; 0.28; 0.32; 0.31];   %between-occasion
            d.sig_trl    = [0.20; 0.25; 0.23; 0.21; 0.22];   %between-trial
            d.sig_trlxid = [0.20; 0.22; 0.18; 0.21; 0.19];   %trial x person
            d.sig_occxid = [0.40; 0.42; 0.38; 0.41; 0.39];   %occasion x person
            d.sig_trlxocc = [0.09; 0.10; 0.11; 0.09; 0.10];  %trial x occasion
            d.sig_err    = [4.90; 5.00; 5.10; 4.80; 5.20];   %residual
        end

        function psyrat_data = makeTrtData(testCase, d, obs, nocc, gcoeff, reltype)
            %Build a minimal TRT results struct. nocc only controls how many
            %occasion labels populate rel.time (the bug trigger); the variance
            %draws and trial count are held fixed.
            psyrat_data = struct();
            psyrat_data.rel = struct('groups', 'none', 'events', 'none', ...
                'analysis', 'trt');
            psyrat_data.rel.time = arrayfun(@(k) sprintf('t%d', k), ...
                1:nocc, 'UniformOutput', false);

            %Reported coefficient is computed with nocc left at its default of
            %1, matching psyrat_relsummary, so the table's reported coefficient
            %genuinely is the single-occasion estimand.
            [iccLl, iccPt, iccUl] = psyrat_rel_trt('gcoeff', gcoeff, ...
                'reltype', reltype, 'bp', d.sig_id, 'bo', d.sig_occ, ...
                'bt', d.sig_trl, 'txp', d.sig_trlxid, 'oxp', d.sig_occxid, ...
                'txo', d.sig_trlxocc, 'err', d.sig_err, 'obs', obs, 'CI', 0.95);

            witsd = mean(d.sig_trlxid);
            betsd = mean(d.sig_id);

            ev = struct();
            ev.name = 'measure';
            ev.icc = struct('m', iccPt, 'll', iccLl, 'ul', iccUl);
            ev.betsd = struct('m', betsd, 'll', betsd, 'ul', betsd);
            ev.witsd = struct('m', witsd, 'll', witsd, 'ul', witsd);
            ev.rel = struct('m', iccPt, 'll', iccLl, 'ul', iccUl, 'meas', obs);
            ev.trlinfo = struct('min', obs, 'max', obs, 'mean', obs, ...
                'med', obs, 'std', 0);

            psyrat_data.relsummary = struct();
            psyrat_data.relsummary.gcoeff = gcoeff;
            psyrat_data.relsummary.reltype = reltype;
            %keep the name consistent with the numeric selector; the fixture
            %sweeps gcoeff 1:2 and psyrat_variancet labels ICC/SEM from this
            gnames = {'dep', 'gen'};
            psyrat_data.relsummary.gcoeff_name = gnames{gcoeff};
            psyrat_data.relsummary.reltype_name = 'ic';
            psyrat_data.relsummary.ciperc = 0.95;
            psyrat_data.relsummary.group = struct();
            psyrat_data.relsummary.group(1).name = '';
            psyrat_data.relsummary.group(1).event = ev;

            psyrat_data.relsummary.data = struct();
            psyrat_data.relsummary.data.g(1).e(1).sig_id.raw = d.sig_id;
            psyrat_data.relsummary.data.g(1).e(1).sig_occ.raw = d.sig_occ;
            psyrat_data.relsummary.data.g(1).e(1).sig_trl.raw = d.sig_trl;
            psyrat_data.relsummary.data.g(1).e(1).sig_trlxid.raw = d.sig_trlxid;
            psyrat_data.relsummary.data.g(1).e(1).sig_occxid.raw = d.sig_occxid;
            psyrat_data.relsummary.data.g(1).e(1).sig_trlxocc.raw = d.sig_trlxocc;
            psyrat_data.relsummary.data.g(1).e(1).sig_err.raw = d.sig_err;

            testCase.assertEqual(numel(psyrat_data.rel.time), nocc);
        end

        function pt = parsePointEstimate(~, semStr)
            %Pull the leading point estimate from ' %0.2f CI [%0.2f %0.2f]'.
            tok = regexp(semStr, '[-+]?\d*\.?\d+', 'match');
            assert(~isempty(tok), 'SEM string had no numeric point estimate: %s', semStr);
            pt = str2double(tok{1});
        end

        function s = semString(~, tbl)
            %Resolve the SEM column by prefix. The SEM is estimand-dependent, so
            %psyrat_variancet names it SEM_Dep or SEM_Gen (RC-11); these tests
            %check the VALUE and should not be coupled to which one is selected.
            vn = tbl.Properties.VariableNames;
            hit = vn(startsWith(vn, 'SEM'));
            assert(isscalar(hit), 'expected exactly one SEM column, found %d', numel(hit));
            s = tbl.(hit{1}){1};
        end

    end

    methods (Static, Access = private)

        function sem = semPointEstimate(d, obs, nocc, gcoeff, reltype)
            %Analytic SEM point estimate: sqrt of the per-draw error variance
            %(identical err_term to psyrat_rel_trt / psyrat_sem_stats),
            %averaged over draws. With nocc = 1 this is sqrt(error variance) of
            %the reported coefficient.
            txp = d.sig_trlxid.^2;
            oxp = d.sig_occxid.^2;
            bt = d.sig_trl.^2;
            bo = d.sig_occ.^2;
            txo = d.sig_trlxocc.^2;
            err = d.sig_err.^2;

            if gcoeff == 1 && reltype == 1
                errTerm = (txp ./ obs) + (err ./ (obs*nocc)) + ...
                    (bt ./ obs) + (txo ./ (obs*nocc));
            elseif gcoeff == 1 && reltype == 2
                errTerm = (oxp ./ nocc) + (err ./ (obs*nocc)) + ...
                    (bo ./ nocc) + (txo ./ (obs*nocc));
            elseif gcoeff == 1 && reltype == 3
                errTerm = (txp ./ obs) + (oxp ./ nocc) + ...
                    (err ./ (obs*nocc)) + (bt ./ obs) + ...
                    (bo ./ nocc) + (txo ./ (obs*nocc));
            elseif gcoeff == 2 && reltype == 1
                errTerm = (txp ./ obs) + (err ./ (obs*nocc));
            elseif gcoeff == 2 && reltype == 2
                errTerm = (oxp ./ nocc) + (err ./ (obs*nocc));
            else %gcoeff == 2 && reltype == 3
                errTerm = (txp ./ obs) + (oxp ./ nocc) + (err ./ (obs*nocc));
            end

            sem = mean(sqrt(errTerm));
        end

    end

end
