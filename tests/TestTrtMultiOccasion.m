classdef TestTrtMultiOccasion < PsyRATTestBase
    %Coverage for the user-selectable multi-occasion composite test-retest
    %(TRT) reliability option (author estimand convention, scientific audit
    %sec.14 / finding F2).
    %
    %The default TRT estimand is a single-occasion score (n'_o = 1): the
    %reliability of a score from n'_i trials at one occasion, generalizing over
    %the occasion facet. A multi-occasion composite score (n'_o = k) is the
    %reliability of a score averaged across k occasions. The SAME n'_o must be
    %threaded into BOTH the reported coefficient (psyrat_rel_trt) and its SEM
    %(psyrat_variancet): single- and multi-occasion coefficients are never
    %paired with each other's SEM.
    %
    %These tests exercise the production WIRING (psyrat_relsummary ->
    %psyrat_variancet), not only the formula function, because the recently
    %fixed one-facet cut-score bug (F1) escaped the suite precisely because
    %only the function - not the wiring - was tested.

    methods (Test)

        % ---- relsummary-level wiring tests -------------------------------

        function testDefaultIsSingleOccasionAndUnchanged(testCase)
            %With no occasion option supplied, relsummary must produce the
            %single-occasion (n'_o = 1) estimand - byte-identical to the prior
            %behavior - and label it as such.
            d = testCase.sampleDraws();
            psyrat_data = testCase.makeTrtRel(d, 8, {'t1','t2','t3','t4'});

            out = testCase.runTrt(psyrat_data, 1, 1, {});

            testCase.verifyEqual(out.relsummary.noccmode, 1, ...
                'Default occasion mode must be single-occasion.');
            testCase.verifyEqual(out.relsummary.nocc, 1, ...
                'Default effective n''_o must be 1.');
            testCase.verifyEqual(out.relsummary.nobservedocc, 4, ...
                'Observed occasion count should come from numel(rel.time).');

            %the reported coefficient is psyrat_rel_trt at the relsummary obs
            %with n'_o = 1
            obs = out.relsummary.group(1).event(1).rel.meas;
            [~, expPt, ~] = psyrat_rel_trt('gcoeff', 1, 'reltype', 1, ...
                'bp', d.sig_id, 'bo', d.sig_occ, 'bt', d.sig_trl, ...
                'txp', d.sig_trlxid, 'oxp', d.sig_occxid, ...
                'txo', d.sig_trlxocc, 'err', d.sig_err, ...
                'obs', obs, 'nocc', 1, 'CI', 0.95);
            testCase.verifyEqual(out.relsummary.group(1).event(1).rel.m, ...
                expPt, 'AbsTol', 1e-9, ...
                'Default coefficient must equal the single-occasion value.');

            %and an explicit single-occasion request gives the identical result
            outExplicit = testCase.runTrt(psyrat_data, 1, 1, {'noccmode', 1});
            testCase.verifyEqual(outExplicit.relsummary.group(1).event(1).rel.m, ...
                out.relsummary.group(1).event(1).rel.m, 'AbsTol', 1e-12, ...
                'Explicit noccmode=1 must match the default.');
        end

        function testMultiOccasionCoefficientAndSemShareNocc(testCase)
            %The central guarantee: with multi-occasion (k = observed
            %occasions), the reported coefficient equals psyrat_rel_trt(...,
            %'nocc',k) AND the reported SEM equals sqrt(that coefficient's error
            %variance) at the SAME n'_o = k. They are never mismatched.
            d = testCase.sampleDraws();
            occ = {'t1','t2','t3','t4'};
            k = numel(occ);
            psyrat_data = testCase.makeTrtRel(d, 8, occ);

            for gcoeff = 1:2
                for reltype = 1:3
                    out = testCase.runTrt(psyrat_data, gcoeff, reltype, ...
                        {'noccmode', 2});

                    testCase.verifyEqual(out.relsummary.noccmode, 2);
                    testCase.verifyEqual(out.relsummary.nocc, k, ...
                        'Multi-occasion k should default to observed occasions.');

                    obs = out.relsummary.group(1).event(1).rel.meas;

                    %(1) reported coefficient == psyrat_rel_trt at n'_o = k
                    [~, expPt, ~] = psyrat_rel_trt('gcoeff', gcoeff, ...
                        'reltype', reltype, 'bp', d.sig_id, 'bo', d.sig_occ, ...
                        'bt', d.sig_trl, 'txp', d.sig_trlxid, ...
                        'oxp', d.sig_occxid, 'txo', d.sig_trlxocc, ...
                        'err', d.sig_err, 'obs', obs, 'nocc', k, 'CI', 0.95);
                    msg = sprintf('gcoeff=%d reltype=%d', gcoeff, reltype);
                    testCase.verifyEqual( ...
                        out.relsummary.group(1).event(1).rel.m, expPt, ...
                        'AbsTol', 1e-9, ...
                        ['Multi-occasion coefficient must use n''_o = k. ' msg]);

                    %(2) reported SEM == sqrt(coefficient error variance) at k
                    tbl = psyrat_variancet('psyrat_data', out, 'gui', 0);
                    reportedSem = testCase.parsePointEstimate(testCase.semString(tbl));
                    expSem = TestTrtMultiOccasion.semPointEstimate( ...
                        d, obs, k, gcoeff, reltype);
                    testCase.verifyEqual(reportedSem, round(expSem, 2), ...
                        'AbsTol', 0.011, ...
                        ['Multi-occasion SEM must be sqrt(error var) at ' ...
                         'n''_o = k. ' msg]);

                    %the multi-occasion SEM must DIFFER from the single-occasion
                    %SEM (otherwise the option does nothing for these draws)
                    expSemSingle = TestTrtMultiOccasion.semPointEstimate( ...
                        d, obs, 1, gcoeff, reltype);
                    testCase.verifyGreaterThan( ...
                        abs(round(expSem, 2) - round(expSemSingle, 2)), 0.1, ...
                        ['Draws must make single- vs multi-occasion SEM ' ...
                         'clearly distinguishable. ' msg]);
                end
            end
        end

        function testMultiOccasionDiffersFromSingleOccasion(testCase)
            %A multi-occasion composite (k>1) must report a different
            %coefficient AND SEM than the single-occasion default for the same
            %fabricated draws, confirming the option is actually wired.
            d = testCase.sampleDraws();
            psyrat_data = testCase.makeTrtRel(d, 8, {'t1','t2','t3','t4'});

            single = testCase.runTrt(psyrat_data, 1, 1, {});
            multi = testCase.runTrt(psyrat_data, 1, 1, {'noccmode', 2});

            testCase.verifyNotEqual(multi.relsummary.group(1).event(1).rel.m, ...
                single.relsummary.group(1).event(1).rel.m);

            semSingle = testCase.parsePointEstimate(testCase.semString( ...
                psyrat_variancet('psyrat_data', single, 'gui', 0)));
            semMulti = testCase.parsePointEstimate(testCase.semString( ...
                psyrat_variancet('psyrat_data', multi, 'gui', 0)));
            testCase.verifyNotEqual(semMulti, semSingle, ...
                'Multi-occasion SEM must differ from single-occasion SEM.');
        end

        function testExplicitKProjectionIsHonored(testCase)
            %An explicit k (D-study projection) overrides the observed-occasion
            %default and is used consistently in the coefficient and SEM.
            d = testCase.sampleDraws();
            psyrat_data = testCase.makeTrtRel(d, 8, {'t1','t2','t3','t4'});
            kProj = 2;

            out = testCase.runTrt(psyrat_data, 1, 1, {'noccmode', 2, 'nocc', kProj});

            testCase.verifyEqual(out.relsummary.nocc, kProj, ...
                'Explicit k must override the observed-occasion default.');

            obs = out.relsummary.group(1).event(1).rel.meas;
            [~, expPt, ~] = psyrat_rel_trt('gcoeff', 1, 'reltype', 1, ...
                'bp', d.sig_id, 'bo', d.sig_occ, 'bt', d.sig_trl, ...
                'txp', d.sig_trlxid, 'oxp', d.sig_occxid, ...
                'txo', d.sig_trlxocc, 'err', d.sig_err, ...
                'obs', obs, 'nocc', kProj, 'CI', 0.95);
            testCase.verifyEqual(out.relsummary.group(1).event(1).rel.m, ...
                expPt, 'AbsTol', 1e-9);

            tbl = psyrat_variancet('psyrat_data', out, 'gui', 0);
            reportedSem = testCase.parsePointEstimate(testCase.semString(tbl));
            expSem = TestTrtMultiOccasion.semPointEstimate(d, obs, kProj, 1, 1);
            testCase.verifyEqual(reportedSem, round(expSem, 2), 'AbsTol', 0.011);
        end

        function testInvalidNoccRejected(testCase)
            %relsummary must reject non-positive / non-integer composite sizes.
            d = testCase.sampleDraws();
            psyrat_data = testCase.makeTrtRel(d, 8, {'t1','t2','t3','t4'});

            testCase.verifyError(@() testCase.runTrt( ...
                psyrat_data, 1, 1, {'noccmode', 2, 'nocc', 0}), ...
                'varargin:nocc');
            testCase.verifyError(@() testCase.runTrt( ...
                psyrat_data, 1, 1, {'noccmode', 2, 'nocc', 2.5}), ...
                'varargin:nocc');
            testCase.verifyError(@() testCase.runTrt( ...
                psyrat_data, 1, 1, {'noccmode', 3}), ...
                'varargin:noccmode');
        end

        % ---- function-level test -----------------------------------------

        function testRelTrtNoccScalesErrorTerm(testCase)
            %Sanity at the function level: increasing n'_o monotonically
            %lowers the (shared) occasion-scaled error and matches the manual
            %formula for the coefficient of stability (dependability).
            d = testCase.sampleDraws();
            obs = 8;

            bp = d.sig_id.^2; txp = d.sig_trlxid.^2;
            oxp = d.sig_occxid.^2; bo = d.sig_occ.^2;
            txo = d.sig_trlxocc.^2; err = d.sig_err.^2;

            for k = [1 2 4]
                [~, pt, ~] = psyrat_rel_trt('gcoeff', 1, 'reltype', 2, ...
                    'bp', d.sig_id, 'bo', d.sig_occ, 'bt', d.sig_trl, ...
                    'txp', d.sig_trlxid, 'oxp', d.sig_occxid, ...
                    'txo', d.sig_trlxocc, 'err', d.sig_err, ...
                    'obs', obs, 'nocc', k, 'CI', 0.95);

                uni = bp + (txp ./ obs);
                errTerm = (oxp ./ k) + (err ./ (obs*k)) + ...
                    (bo ./ k) + (txo ./ (obs*k));
                manual = mean(uni ./ (uni + errTerm));
                testCase.verifyEqual(pt, manual, 'AbsTol', 1e-12, ...
                    sprintf('psyrat_rel_trt must match manual formula at k=%d', k));
            end
        end

        function testTrialsPlotHonorsNoccAndMatchesHeadlessCurve(testCase)
            %RC-10. The reliability-vs-trials FIGURE passed no n'_o, so it drew
            %the single-occasion curve while the table beside it reported the
            %multi-occasion coefficient - under an identical "Dependability"
            %header. The user reads a required trial count off where the curve
            %crosses the cutoff line, which is the same question the table's
            %cutoff column answers, so the two artifacts disagreed.
            %
            %Driving the figure headlessly and reading its YData back is the
            %only way to catch this: the plot is not reachable from the calc
            %layer, and the sibling headless curve had the same omission, so
            %comparing the two would have agreed while both were wrong.
            psyrat_data = PsyRATTestDataFactory.makeTRTSummaryData();
            ntrials = 20;

            y1 = testCase.plotCurveY(psyrat_data, ntrials);

            multi = psyrat_data;
            multi.relsummary.noccmode = 2;
            multi.relsummary.nocc = 4;
            y4 = testCase.plotCurveY(multi, ntrials);

            %the curve must actually move; under the shipped default this fix
            %is a strict no-op, so an unmoved curve means nothing was tested
            testCase.verifyGreaterThan(max(abs(y4 - y1)), 1e-3, ...
                'The plotted curve must honor the selected n''_o.');

            %and it must be the SAME quantity the headless *_curve export
            %writes - the parity psyrat_report's local_curve_table promises.
            %Until RC-10 both sides matched only because both omitted nocc.
            if psyrat_data.relsummary.gcoeff == 2
                ptcol = 'Generalizability_PointEstimate';
            else
                ptcol = 'Dependability_PointEstimate';
            end
            report = psyrat_report(multi, 'ntrials', ntrials);
            testCase.verifyEqual(round(y4(:), 4), ...
                report.tables.curve.(ptcol), 'AbsTol', 1e-12, ...
                'The figure and the headless curve must be the same estimand.');
        end

        function testTrialsPlotNamesTheOccasionEstimand(testCase)
            %RC-10. The axis label reads only "Dependability", so without the
            %occasion estimand in the window name a saved screenshot of the
            %single-occasion curve is indistinguishable from the composite one.
            psyrat_data = PsyRATTestDataFactory.makeTRTSummaryData();

            f = psyrat_trt_relvtrialsplot('psyrat_data', psyrat_data, ...
                'trials', [1 20], 'relcutoff', 0.7, 'relline', 2, 'CI', 0.95);
            closer = onCleanup(@() close(f));
            testCase.verifySubstring(get(f, 'Name'), ...
                'Single-occasion TRT reliability, n''_o = 1');
            clear closer;

            multi = psyrat_data;
            multi.relsummary.noccmode = 2;
            multi.relsummary.nocc = 4;
            f2 = psyrat_trt_relvtrialsplot('psyrat_data', multi, ...
                'trials', [1 20], 'relcutoff', 0.7, 'relline', 2, 'CI', 0.95);
            closer2 = onCleanup(@() close(f2));
            testCase.verifySubstring(get(f2, 'Name'), ...
                'Multi-occasion composite TRT reliability, n''_o = 4');
        end

    end

    methods (Access = private)

        function y = plotCurveY(~, psyrat_data, ntrials)
            %Draw the reliability-vs-trials figure and read the plotted curve
            %back. The figure also holds a flat cutoff reference line, so the
            %curve is identified by its length rather than by draw order.
            f = psyrat_trt_relvtrialsplot('psyrat_data', psyrat_data, ...
                'trials', [1 ntrials], 'relcutoff', 0.7, 'relline', 2, ...
                'CI', 0.95);
            closer = onCleanup(@() close(f));
            y = [];
            lines = findall(f, 'Type', 'line');
            for k = 1:numel(lines)
                yy = get(lines(k), 'YData');
                if numel(yy) == ntrials
                    y = yy(:);
                end
            end
            assert(~isempty(y), 'Could not find the plotted reliability curve.');
        end

        function d = sampleDraws(~)
            %Posterior draws (standard-deviation units) with a large residual,
            %so the occasion-scaled (n'_o-dependent) terms dominate and the
            %single- vs multi-occasion coefficients/SEMs are clearly different.
            d = struct();
            d.sig_id     = [2.00; 2.10; 1.90; 2.05; 1.95];
            d.sig_occ    = [0.30; 0.35; 0.28; 0.32; 0.31];
            d.sig_trl    = [0.20; 0.25; 0.23; 0.21; 0.22];
            d.sig_trlxid = [0.20; 0.22; 0.18; 0.21; 0.19];
            d.sig_occxid = [0.40; 0.42; 0.38; 0.41; 0.39];
            d.sig_trlxocc = [0.09; 0.10; 0.11; 0.09; 0.10];
            d.sig_err    = [4.90; 5.00; 5.10; 4.80; 5.20];
        end

        function psyrat_data = makeTrtRel(~, d, ntrl, occ)
            %Build a minimal PRE-relsummary TRT results struct (rel.out draws +
            %a rel.data trial table), the input that psyrat_relsummary consumes.
            subs = {'s1','s2','s3'};
            ids = {}; times = {}; meas = [];
            for s = 1:numel(subs)
                for o = 1:numel(occ)
                    for t = 1:ntrl
                        ids{end+1,1} = subs{s}; %#ok<AGROW>
                        times{end+1,1} = occ{o}; %#ok<AGROW>
                        meas(end+1,1) = numel(meas) + 1; %#ok<AGROW>
                    end
                end
            end

            rel = struct();
            rel.analysis = 'trt';
            rel.groups = 'none';
            rel.events = 'none';
            rel.filename = 'trt_multiocc_fixture.psyrat';
            rel.nchains = 2;
            rel.niter = 100;
            rel.time = occ(:);
            rel.data = table(ids, meas, times, ...
                'VariableNames', {'id','meas','time'});

            rel.out = struct();
            rel.out.labels = {'measure'};
            rel.out.mu = [0.10; 0.12; 0.11; 0.09; 0.10];
            rel.out.sig_id = d.sig_id;
            rel.out.sig_occ = d.sig_occ;
            rel.out.sig_trl = d.sig_trl;
            rel.out.sig_trlxid = d.sig_trlxid;
            rel.out.sig_occxid = d.sig_occxid;
            rel.out.sig_trlxocc = d.sig_trlxocc;
            rel.out.sig_err = d.sig_err;

            psyrat_data = struct('rel', rel);
        end

        function out = runTrt(~, psyrat_data, gcoeff, reltype, extra)
            %Drive psyrat_relsummary for a test-retest analysis.
            args = {'psyrat_data', psyrat_data, 'analysis', 'trt', ...
                'gcoeff', gcoeff, 'reltype', reltype, ...
                'relcutoff', 0.2, 'meascutoff', 2, 'relcentmeas', 1, ...
                'CI', 0.95};
            out = psyrat_relsummary(args{:}, extra{:});
        end

        function pt = parsePointEstimate(~, semStr)
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
            %averaged over draws. This is sqrt(the coefficient's error variance)
            %at the supplied n'_o.
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
