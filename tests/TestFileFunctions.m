classdef TestFileFunctions < PsyRATTestBase

    methods (Test)
        function testDefaultsReturnsExpectedPreferenceFields(testCase)
            prefs = psyrat_defaults();

            testCase.verifyTrue(isfield(prefs, 'proc'));
            testCase.verifyTrue(isfield(prefs, 'view'));
            testCase.verifyEqual(prefs.proc.nchains, 4);
            testCase.verifyEqual(prefs.proc.nwarmup, 5000);
            testCase.verifyEqual(prefs.proc.nsampling, 5000);
            testCase.verifyEqual(prefs.proc.niter, 10000);
            testCase.verifyEqual(prefs.proc.niter, prefs.proc.nwarmup + prefs.proc.nsampling);
            testCase.verifyEqual(prefs.proc.diffwpcov, 1);
            testCase.verifyEqual(prefs.view.depvalue, 0.8, 'AbsTol', 1e-12);
            testCase.verifyTrue(isfield(prefs.view, 'criterioncutoff'));
            testCase.verifyTrue(isfield(prefs.view, 'plotcriterion'));
            testCase.verifyTrue(isfield(prefs.view, 'tablecriterion'));
            testCase.verifyEqual(prefs.view.plotcriterion, 1);
            testCase.verifyEqual(prefs.view.tablecriterion, 1);
        end

        function testDefineVersionIsSemanticVersionString(testCase)
            % The version may carry a pre-release suffix ('0.1.0-beta', B4/W7),
            % so the anchor allows an optional '-<tag>' after x.y.z. The LEADING
            % x.y.z is still required and still anchored, because that is what
            % local_compare_semver in psyrat_checkgithubreleases parses with
            % sscanf('%d.%d.%d') and rejects with 'psyrat:versionparse' when it
            % does not yield three parts. A suffix is safe; a malformed leading
            % triple is not, and this pin still catches that.
            ver = psyrat_defineversion();
            testCase.verifyClass(ver, 'char');
            testCase.verifyNotEmpty(regexp(ver, '^\d+\.\d+\.\d+(-[0-9A-Za-z.]+)?$', 'once'));

            % Positive control for the half that actually matters downstream:
            % the leading triple must survive the semver parse the release
            % checker performs.
            parts = sscanf(ver, '%d.%d.%d')';
            testCase.verifyNumElements(parts, 3, ...
                'version must expose a parseable x.y.z to the release checker');
        end

        function testFindPrefsDataFromVararginAndFallback(testCase)
            prefsIn = struct('a', 1);
            dataIn = struct('b', 2);
            [prefsOut, dataOut] = psyrat_findprefsdata({'psyrat_prefs', prefsIn, 'psyrat_data', dataIn});
            testCase.verifyEqual(prefsOut, prefsIn);
            testCase.verifyEqual(dataOut, dataIn);

            [prefsOut2, dataOut2] = psyrat_findprefsdata({'x', 1});
            testCase.verifyEqual(prefsOut2, '');
            testCase.verifyEqual(dataOut2, '');
        end

        function testReadTableLoadsCsv(testCase)
            file = localWriteCsv(table({'s1'; 's2'}, [1.1; 2.2], ...
                'VariableNames', {'id', 'measurement'}));

            out = psyrat_readtable('file', file);
            testCase.verifyEqual(width(out), 2);
            testCase.verifyEqual(height(out), 2);
            testCase.verifyTrue(any(strcmpi(out.Properties.VariableNames, 'id')));
        end

        function testReadTableValidationAndUnsupportedExtension(testCase)
            testCase.verifyError(@() psyrat_readtable('idcol', ','), 'varargin:nofile');
            testCase.verifyError(@() psyrat_readtable('file', 'dummy.unsupported'), 'ext:filetype');
            testCase.verifyError(@() psyrat_readtable('file'), 'varargin:incomplete');
        end

        function testLoadFileBasicDefaultColumns(testCase)
            file = localWriteCsv(table({'s1'; 's1'; 's2'}, [1.1; 1.2; 2.1], ...
                'VariableNames', {'id', 'measurement'}));

            out = psyrat_loadfile('file', file);
            testCase.verifyEqual(out.Properties.VariableNames, {'id', 'meas'});
            testCase.verifyEqual(height(out), 3);
            testCase.verifyClass(out.id, 'cell');
        end

        function testLoadFileConvertsNumericIdsToStrings(testCase)
            file = localWriteCsv(table([1; 1; 2], [1.1; 1.2; 2.1], ...
                'VariableNames', {'id', 'measurement'}));

            out = psyrat_loadfile('file', file);
            testCase.verifyClass(out.id, 'cell');
            testCase.verifyEqual(out.id{1}, '1');
            testCase.verifyEqual(out.id{end}, '2');
        end

        function testLoadFileValidationErrors(testCase)
            missingMeas = localWriteCsv(table({'s1'; 's2'}, [1; 2], ...
                'VariableNames', {'id', 'score'}));
            testCase.verifyError(@() psyrat_loadfile('file', missingMeas), 'varargin:colheaders');

            nonNumericMeas = localWriteCsv(table({'s1'; 's2'}, {'x'; 'y'}, ...
                'VariableNames', {'id', 'measurement'}));
            testCase.verifyError(@() psyrat_loadfile('file', nonNumericMeas), 'meas:notnumeric');

            testCase.verifyError(@() psyrat_loadfile('idcol', 'id'), 'varargin:nofile');
            testCase.verifyError(@() psyrat_loadfile('file'), 'varargin:incomplete');
        end

        function testLoadFileFiltersGroupsAndEventsWhenUsingToolboxStruct(testCase)
            raw = table( ...
                {'s1'; 's1'; 's2'; 's2'; 's3'}, ...
                [1.0; 1.1; 2.0; 2.1; 3.0], ...
                {'g1'; 'g1'; 'g2'; 'g2'; 'g1'}, ...
                {'A'; 'B'; 'A'; 'B'; 'A'}, ...
                'VariableNames', {'id', 'measurement', 'group', 'event'});

            psyrat_data = struct();
            psyrat_data.raw.filename = localWriteCsv(raw);
            psyrat_data.raw.data = raw;
            psyrat_data.proc.idheader = 'id';
            psyrat_data.proc.measheader = 'measurement';
            psyrat_data.proc.groupheader = 'group';
            psyrat_data.proc.whichgroups = {'g1'};
            psyrat_data.proc.eventheader = 'event';
            psyrat_data.proc.whichevents = {'A'};
            psyrat_data.proc.timeheader = '';
            psyrat_data.proc.whichtimes = '';

            out = psyrat_loadfile('psyrat_prefs', struct(), 'psyrat_data', psyrat_data);

            testCase.verifyEqual(unique(out.group), {'g1'});
            testCase.verifyEqual(unique(out.event), {'A'});
            testCase.verifyEqual(height(out), 2);
        end

        function testLoadFileFiltersGroupsAndEventsInDirectArgPath(testCase)
            % B6: programmatic which* filtering must subset rows in the
            % direct-argument path, not only the GUI/struct path. Previously the
            % parsed whichgroups/whichevents/whichtimes locals were ignored
            % because the filter blocks only consulted psyrat_data.proc.*.
            raw = table( ...
                {'s1'; 's1'; 's2'; 's2'; 's3'}, ...
                [1.0; 1.1; 2.0; 2.1; 3.0], ...
                {'g1'; 'g1'; 'g2'; 'g2'; 'g1'}, ...
                {'A'; 'B'; 'A'; 'B'; 'A'}, ...
                'VariableNames', {'id', 'measurement', 'group', 'event'});
            file = localWriteCsv(raw);

            out = psyrat_loadfile('file', file, ...
                'groupcol', 'group', 'whichgroups', {'g1'}, ...
                'eventcol', 'event', 'whichevents', {'A'});

            testCase.verifyEqual(unique(out.group), {'g1'});
            testCase.verifyEqual(unique(out.event), {'A'});
            testCase.verifyEqual(height(out), 2);

            % A missing requested group still errors in the direct-arg path.
            testCase.verifyError(@() psyrat_loadfile('file', file, ...
                'groupcol', 'group', 'whichgroups', {'g3'}), ...
                'groups:groupmismatch');
        end

        function testLoadFileErrorsOnMissingRequestedGroup(testCase)
            raw = table( ...
                {'s1'; 's1'; 's2'; 's2'}, ...
                [1.0; 1.1; 2.0; 2.1], ...
                {'g1'; 'g1'; 'g2'; 'g2'}, ...
                'VariableNames', {'id', 'measurement', 'group'});

            psyrat_data = struct();
            psyrat_data.raw.filename = localWriteCsv(raw);
            psyrat_data.raw.data = raw;
            psyrat_data.proc.idheader = 'id';
            psyrat_data.proc.measheader = 'measurement';
            psyrat_data.proc.groupheader = 'group';
            psyrat_data.proc.whichgroups = {'g3'};
            psyrat_data.proc.eventheader = '';
            psyrat_data.proc.whichevents = '';
            psyrat_data.proc.timeheader = '';
            psyrat_data.proc.whichtimes = '';

            testCase.verifyError(@() psyrat_loadfile( ...
                'psyrat_prefs', struct(), 'psyrat_data', psyrat_data), ...
                'groups:groupmismatch');
        end

        function testLoadFilePopulatesOccasionColumnFromTimecol(testCase)
            % Regression (B1): the programmatic 'timecol' argument was
            % ignored because the lookup searched for 'eventcol', and the
            % direct-argument path then dereferenced an empty psyrat_data
            % during filtering. The occasion column should now populate.
            file = localWriteCsv(table( ...
                {'s1'; 's1'; 's2'; 's2'}, ...
                [1.1; 1.2; 2.1; 2.2], ...
                {'t1'; 't2'; 't1'; 't2'}, ...
                'VariableNames', {'id', 'measurement', 'occ'}));

            out = psyrat_loadfile('file', file, 'timecol', 'occ');
            testCase.verifyTrue(any(strcmp(out.Properties.VariableNames, 'time')));
            testCase.verifyEqual(height(out), 4);
            testCase.verifyEqual(unique(out.time), {'t1'; 't2'});
        end

        function testDefaultPriorsExposesPopulationLevelScales(testCase)
            % Guards the SCI-003 prior data model and the SCI-004 widened
            % occasion-prior defaults, and that the priors are exposed on the
            % preference struct.
            p = psyrat_default_priors();
            testCase.verifyEqual(p.single.mu, 100);
            testCase.verifyEqual(p.single.sig_u, 40);
            testCase.verifyEqual(p.single.sig_e, 40);
            testCase.verifyEqual(p.trt.sig_id, 10);
            testCase.verifyEqual(p.trt.sig_occ, 5);      % SCI-004 widened
            testCase.verifyEqual(p.trt.sig_trlxocc, 1);  % SCI-004 widened
            testCase.verifyEqual(p.diff.b, 5);
            testCase.verifyEqual(p.diff.sd_id_sigma, 2.5);

            prefs = psyrat_defaults();
            testCase.verifyTrue(isfield(prefs.proc, 'priors'));
            testCase.verifyEqual(prefs.proc.priors.trt.sig_occ, 5);
        end

        function testLoadFileReportsAllMissingRequiredHeaders(testCase)
            % Regression (B2): a headerror/headererror typo meant only one
            % missing required header was reported. With both the id and
            % measurement columns absent, both names should be listed.
            bad = localWriteCsv(table({'s1'; 's2'}, [1; 2], ...
                'VariableNames', {'subj', 'score'}));

            err = [];
            try
                psyrat_loadfile('file', bad);
            catch err %#ok<NASGU>
            end
            testCase.verifyNotEmpty(err);
            testCase.verifyEqual(err.identifier, 'varargin:colheaders');
            testCase.verifyTrue(contains(err.message, 'Subject ID'));
            testCase.verifyTrue(contains(err.message, 'Measurement'));
        end
    end
end

function file = localWriteCsv(tbl)
file = [tempname '.csv'];
writetable(tbl, file);
end
