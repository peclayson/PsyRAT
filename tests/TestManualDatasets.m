classdef TestManualDatasets < PsyRATTestBase
    % Welds the manual's tutorial datasets to their generators, and both to
    % the toolbox's loader.
    %
    % The user manual prints expected numbers produced from the committed
    % CSVs in test_data/ (manual_ern_singlesession.csv,
    % manual_ern_testretest.csv, manual_ern_splits.csv, manual_ern_dynrel.csv,
    % manual_power_concurrent.csv) plus the legacy-header autodetect fixture
    % (legacy_headers_fixture.csv). Three
    % properties have to hold or the manual silently rots:
    %
    %   1. IDENTITY: each committed CSV is exactly what its committed, seeded
    %      generator in tests/helpers/ produces. Nobody can hand-edit a CSV
    %      (or edit a generator) without this test forcing the pair, and the
    %      manual's expected values, back into agreement.
    %   2. CONTRACT: each CSV round-trips through psyrat_loadfile under the
    %      canonical column names, with the design the manual describes
    %      (between-person groups, complete cells, positive-integer weights).
    %   3. SHAPE: the headline design facts the manual states (participant
    %      counts, events, occasions) are pinned here, so a parameter change
    %      in a generator fails loudly instead of drifting past the prose.
    %
    % Numeric comparisons use a 1e-9 absolute tolerance rather than byte
    % equality so the weld survives MATLAB-version float-printing quirks;
    % the generators round measurements to 6 decimals before writing, so in
    % practice the round trip is exact.

    methods (Test)

        % ---- 1. identity: committed CSV == regenerated table ----------------

        function testSinglesessionRegenerationMatchesCommitted(testCase)
            localVerifyRegen(testCase, @psyrat_make_manual_singlesession, ...
                'manual_ern_singlesession.csv');
        end

        function testTestretestRegenerationMatchesCommitted(testCase)
            localVerifyRegen(testCase, @psyrat_make_manual_testretest, ...
                'manual_ern_testretest.csv');
        end

        function testSplitsRegenerationMatchesCommitted(testCase)
            localVerifyRegen(testCase, @psyrat_make_manual_splits, ...
                'manual_ern_splits.csv');
        end

        function testDynrelRegenerationMatchesCommitted(testCase)
            localVerifyRegen(testCase, @psyrat_make_manual_dynrel, ...
                'manual_ern_dynrel.csv');
        end

        function testPowerConcurrentRegenerationMatchesCommitted(testCase)
            localVerifyRegen(testCase, @psyrat_make_manual_power_concurrent, ...
                'manual_power_concurrent.csv');
        end

        function testLegacyFixtureRegenerationMatchesCommitted(testCase)
            localVerifyRegen(testCase, @psyrat_make_legacy_headers_fixture, ...
                'legacy_headers_fixture.csv');
        end

        % ---- 2 + 3. loader round trip and pinned design ---------------------

        function testSinglesessionLoadsAndMatchesTheDocumentedDesign(testCase)
            f = localDataFile(testCase, 'manual_ern_singlesession.csv');
            d = psyrat_loadfile('file', f, 'idcol', 'id', 'meascol', 'meas', ...
                'groupcol', 'group', 'eventcol', 'event');

            testCase.verifyEqual(height(d), 8908);
            testCase.verifyEqual(numel(unique(d.id)), 80);
            testCase.verifyEqual(sort(unique(d.group)), {'GroupA'; 'GroupB'});
            testCase.verifyEqual(numel(unique(d.event)), 4);
            testCase.verifyFalse(any(ismissing(d.meas)));

            % Groups are BETWEEN-person: every participant appears in exactly
            % one group. (The retired testdata.csv violated this; the manual
            % teaches group as a between-person stratum, so pin it.)
            g = varfun(@(x) numel(unique(x)), d, ...
                'InputVariables', 'group', 'GroupingVariables', 'id');
            testCase.verifyTrue(all(g.Fun_group == 1), ...
                'a participant appears in more than one group');

            % Every participant has every event (the loader's completeness
            % requirement; a missing cell would abort a subset analysis).
            e = varfun(@(x) numel(unique(x)), d, ...
                'InputVariables', 'event', 'GroupingVariables', 'id');
            testCase.verifyTrue(all(e.Fun_event == 4), ...
                'a participant is missing an event cell');

            % Trial counts stay inside the documented ranges the tutorials
            % state: errors randi([12 36]), corrects randi([24 40]). The .70
            % walkthrough depends on the error range reaching 36.
            iserr = contains(d.event, '_err');
            ce = groupsummary(d(iserr, :), {'id', 'event'});
            testCase.verifyTrue(all(ce.GroupCount >= 12 & ce.GroupCount <= 36), ...
                'error-trial counts left the documented randi([12 36]) range');
            cc = groupsummary(d(~iserr, :), {'id', 'event'});
            testCase.verifyTrue(all(cc.GroupCount >= 24 & cc.GroupCount <= 40), ...
                'correct-trial counts left the documented randi([24 40]) range');

            % Directional sanity on the simulated ERN: error trials are more
            % negative than correct trials on average. This is a gross-defect
            % tripwire (sign flips, swapped labels), not a statistical test.
            testCase.verifyLessThan(mean(d.meas(iserr)), ...
                mean(d.meas(~iserr)) - 1);
        end

        function testTestretestLoadsAndMatchesTheDocumentedDesign(testCase)
            f = localDataFile(testCase, 'manual_ern_testretest.csv');
            d = psyrat_loadfile('file', f, 'idcol', 'id', 'meascol', 'meas', ...
                'eventcol', 'event', 'timecol', 'time');

            testCase.verifyEqual(height(d), 3840);
            testCase.verifyEqual(numel(unique(d.id)), 30);
            testCase.verifyEqual(sort(unique(d.time)), {'t1'; 't2'});
            testCase.verifyEqual(sort(unique(d.event)), {'cor'; 'err'});

            % Balanced by design: 24 err + 40 cor per participant per occasion.
            c = groupsummary(d, {'id', 'time', 'event'});
            testCase.verifyTrue(all(ismember(c.GroupCount, [24 40])), ...
                'trial counts are not the documented 24/40 per cell');
        end

        function testSplitsLoadsAndMatchesTheDocumentedDesign(testCase)
            f = localDataFile(testCase, 'manual_ern_splits.csv');
            d = psyrat_loadfile('file', f, 'idcol', 'id', 'meascol', 'meas', ...
                'weightcol', 'weight');

            testCase.verifyEqual(height(d), 576);
            testCase.verifyEqual(numel(unique(d.id)), 36);
            testCase.verifyTrue(all(d.weight == round(d.weight) & d.weight > 0), ...
                'weights are not positive integers');
            testCase.verifyTrue(all(d.weight >= 4 & d.weight <= 40), ...
                'weights left the documented randi([4 40]) range');
        end

        function testDynrelLoadsAndMatchesTheDocumentedDesign(testCase)
            f = localDataFile(testCase, 'manual_ern_dynrel.csv');
            d = psyrat_loadfile('file', f, 'idcol', 'id', 'meascol', 'meas', ...
                'groupcol', 'group', 'eventcol', 'event', 'dim1col', 'worry');

            testCase.verifyEqual(height(d), 4577);
            testCase.verifyEqual(numel(unique(d.id)), 80);
            testCase.verifyEqual(sort(unique(d.group)), {'GroupA'; 'GroupB'});
            testCase.verifyEqual(sort(unique(d.event)), {'inc_cor'; 'inc_err'});
            testCase.verifyFalse(any(ismissing(d.meas)));

            % Groups are BETWEEN-person and every participant has both events
            % (the loader's completeness requirement).
            g = varfun(@(x) numel(unique(x)), d, ...
                'InputVariables', 'group', 'GroupingVariables', 'id');
            testCase.verifyTrue(all(g.Fun_group == 1), ...
                'a participant appears in more than one group');
            e = varfun(@(x) numel(unique(x)), d, ...
                'InputVariables', 'event', 'GroupingVariables', 'id');
            testCase.verifyTrue(all(e.Fun_event == 2), ...
                'a participant is missing an event cell');

            % The covariate arrives under the loader's canonical name dim1,
            % numeric, and BETWEEN-person: one value per participant (the
            % standardizer refuses a column that varies within participant),
            % an integer on the documented 16 to 80 scale, with between-person
            % variance to standardize.
            testCase.verifyTrue(isnumeric(d.dim1), 'dim1 is not numeric');
            w = varfun(@(x) numel(unique(x)), d, ...
                'InputVariables', 'dim1', 'GroupingVariables', 'id');
            testCase.verifyTrue(all(w.Fun_dim1 == 1), ...
                'worry varies within a participant');
            w1 = varfun(@(x) x(1), d, ...
                'InputVariables', 'dim1', 'GroupingVariables', 'id');
            testCase.verifyTrue(all(w1.Fun_dim1 == round(w1.Fun_dim1)), ...
                'worry is not integer-valued');
            testCase.verifyTrue(all(w1.Fun_dim1 >= 16 & w1.Fun_dim1 <= 80), ...
                'worry left the documented 16 to 80 range');
            testCase.verifyGreaterThan(std(w1.Fun_dim1), 0, ...
                'worry has no between-participant variance');

            % Trial counts stay inside the documented ranges.
            iserr = contains(d.event, '_err');
            ce = groupsummary(d(iserr, :), {'id', 'event'});
            testCase.verifyTrue(all(ce.GroupCount >= 12 & ce.GroupCount <= 36), ...
                'error-trial counts left the documented randi([12 36]) range');
            cc = groupsummary(d(~iserr, :), {'id', 'event'});
            testCase.verifyTrue(all(cc.GroupCount >= 24 & cc.GroupCount <= 40), ...
                'correct-trial counts left the documented randi([24 40]) range');

            % Directional tripwire on the simulated ERN, as for the
            % single-session file.
            testCase.verifyLessThan(mean(d.meas(iserr)), ...
                mean(d.meas(~iserr)) - 1);
        end

        function testPowerConcurrentLoadsAndMatchesTheDocumentedDesign(testCase)
            f = localDataFile(testCase, 'manual_power_concurrent.csv');
            d = psyrat_loadfile('file', f, 'idcol', 'id', 'meascol', 'meas', ...
                'eventcol', 'event', 'dim1col', 'worry');

            testCase.verifyEqual(height(d), 2454);
            testCase.verifyEqual(numel(unique(d.id)), 40);
            testCase.verifyEqual(sort(unique(d.event)), {'alpha'; 'theta'});
            testCase.verifyFalse(any(ismissing(d.meas)));

            % The gamma family's precondition: every score strictly positive.
            testCase.verifyTrue(all(d.meas > 0), ...
                'a power score is not strictly positive');

            % CONCURRENT design: the two bands are scored from the same trials,
            % so their counts are equal within every participant (the pairing
            % the residual-covariance models rely on) and inside the
            % documented randi([20 40]) range.
            c = groupsummary(d, {'id', 'event'});
            nt = c.GroupCount(strcmp(c.event, 'theta'));
            na = c.GroupCount(strcmp(c.event, 'alpha'));
            testCase.verifyEqual(nt, na, ...
                'theta and alpha counts differ within a participant');
            testCase.verifyTrue(all(nt >= 20 & nt <= 40), ...
                'trial counts left the documented randi([20 40]) range');

            % The covariate, as for the dynrel file.
            testCase.verifyTrue(isnumeric(d.dim1), 'dim1 is not numeric');
            w = varfun(@(x) numel(unique(x)), d, ...
                'InputVariables', 'dim1', 'GroupingVariables', 'id');
            testCase.verifyTrue(all(w.Fun_dim1 == 1), ...
                'worry varies within a participant');
            w1 = varfun(@(x) x(1), d, ...
                'InputVariables', 'dim1', 'GroupingVariables', 'id');
            testCase.verifyTrue(all(w1.Fun_dim1 >= 16 & w1.Fun_dim1 <= 80), ...
                'worry left the documented 16 to 80 range');
            testCase.verifyGreaterThan(std(w1.Fun_dim1), 0, ...
                'worry has no between-participant variance');

            % Directional tripwire: theta power is generated larger than alpha.
            testCase.verifyGreaterThan(mean(d.meas(strcmp(d.event, 'theta'))), ...
                mean(d.meas(strcmp(d.event, 'alpha'))));
        end

        function testLegacyFixtureKeepsItsNonCanonicalHeaders(testCase)
            % The whole point of this fixture is its header row: none of the
            % canonical names, and a score column ('ern') that matches no
            % Measurement candidate literal, so the autodetect fallback is
            % exercised. Pin the headers so a well-meaning rename does not
            % quietly defeat TestStartprocColumnAutodetect.
            f = localDataFile(testCase, 'legacy_headers_fixture.csv');
            tbl = psyrat_readtable('file', f);
            testCase.verifyEqual(tbl.Properties.VariableNames, ...
                {'subjid', 'group', 'event', 'ern'});
            testCase.verifyEqual(height(tbl), 80);
        end

    end
end

% ---- local helpers ----------------------------------------------------------

function f = localDataFile(testCase, name)
f = fullfile(testCase.projectRoot(), 'test_data', name);
end

function localVerifyRegen(testCase, generator, committedName)
%Regenerate into a temp file and require agreement with the committed CSV.
committed = localDataFile(testCase, committedName);
testCase.assertTrue(isfile(committed), sprintf( ...
    'committed dataset missing: %s (run %s to create it)', ...
    committed, func2str(generator)));

tmp = [tempname '.csv'];
cleaner = onCleanup(@() delete(tmp)); %#ok<NASGU>
generator(tmp);

a = readtable(committed, 'TextType', 'char');
b = readtable(tmp, 'TextType', 'char');

testCase.verifyEqual(b.Properties.VariableNames, a.Properties.VariableNames, ...
    'regenerated column names differ from the committed CSV');
testCase.verifyEqual(height(b), height(a), ...
    'regenerated row count differs from the committed CSV');

for c = 1:width(a)
    col = a.Properties.VariableNames{c};
    if isnumeric(a.(col))
        testCase.verifyEqual(b.(col), a.(col), 'AbsTol', 1e-9, sprintf( ...
            'numeric column %s differs between generator and committed CSV', col));
    else
        testCase.verifyEqual(b.(col), a.(col), sprintf( ...
            'column %s differs between generator and committed CSV', col));
    end
end
end
