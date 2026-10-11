classdef TestReadtableWarningState < PsyRATTestBase
    %G90: psyrat_readtable turns MATLAB:table:ModifiedVarnames off around its
    %reads, but it discards the state that warning('off',...) returns and
    %then runs an unconditional warning('on',...) once the reads finish. The
    %caller's own setting is therefore lost on two routes:
    %  (a) every exit that reaches the reset leaves the warning ON, even for
    %      a caller that had it off: a successful read, and the "no
    %      delimiter worked" error, which is raised after the reset;
    %  (b) readtable itself throwing skips the reset, so the warning stays
    %      OFF, even for a caller that had it on.
    %The intended fix keeps the returned state and restores it through
    %onCleanup, which covers every exit. Re-enabling the warning in a catch
    %would repair (b) and keep (a), so each exit is tested from BOTH initial
    %states:
    %
    %  exit (how it is reached)                  caller off   caller on
    %  successful read (two-column CSV)          RED          control
    %  readtable throws (the file is missing)    control      RED
    %  no delimiter worked ('@'-delimited file)  RED          control
    %
    %RED and control are the expectations on the unfixed code; nothing here
    %was executed by the author. The controls pass before the fix and must
    %still pass after it.
    %
    %NOT pinned here: that the read still silences the warning at all.
    %Deleting both warning calls would pass all six tests, so the fix's diff
    %has to be read for that.
    %
    %Hygiene: the WHOLE warning state is captured before each test and
    %restored after it, whatever the outcome, and every file lives in a
    %fresh temporary folder that is removed afterwards. No network, no fit.

    properties (Constant)
        WarnId = 'MATLAB:table:ModifiedVarnames';
        %psyrat_readtable's own error identifiers. The erroring-read tests
        %must fail inside readtable, after the warning was turned off, not
        %in one of these checks (all but the last fire before the warning
        %call; the last fires after the reset).
        OwnIds = {'varargin:incomplete', 'varargin:nofile', 'ext:filetype', ...
            'loadingfile:delimiter'};
    end

    methods (TestMethodSetup)
        function preserveWarningState(testCase)
            %the whole state, not just the one identifier, so a failing
            %test cannot leave any warning setting behind for the next class
            prior = warning;
            testCase.addTeardown(@() warning(prior));
        end
    end

    methods (Test)

        % ---- exit 1: a successful read -----------------------------------

        function testSuccessfulReadKeepsCallerOff(testCase)
            file = testCase.writeFile('two_columns.csv', ...
                sprintf('id,meas\ns1,1.1\ns2,2.2\n'));
            testCase.setCallerState('off');
            out = psyrat_readtable('file', file);
            testCase.assertEqual(width(out), 2, 'premise: the read must succeed');
            testCase.verifyEqual(testCase.callerState(), 'off', ...
                ['G90: a successful read turned the caller''s OFF ' ...
                'MATLAB:table:ModifiedVarnames back on']);
        end

        function testSuccessfulReadKeepsCallerOn(testCase)
            %control: a fix that simply leaves the warning off fails here
            file = testCase.writeFile('two_columns.csv', ...
                sprintf('id,meas\ns1,1.1\ns2,2.2\n'));
            testCase.setCallerState('on');
            out = psyrat_readtable('file', file);
            testCase.assertEqual(width(out), 2, 'premise: the read must succeed');
            testCase.verifyEqual(testCase.callerState(), 'on', ...
                'a successful read must leave the caller''s ON state on');
        end

        % ---- exit 2: readtable itself throws -----------------------------

        function testReadtableErrorKeepsCallerOn(testCase)
            testCase.setCallerState('on');
            testCase.readMissingFile();
            testCase.verifyEqual(testCase.callerState(), 'on', ...
                ['G90: readtable threw before the reset, leaving the ' ...
                'caller''s ON MATLAB:table:ModifiedVarnames switched off']);
        end

        function testReadtableErrorKeepsCallerOff(testCase)
            %control: re-enabling the warning in a catch fails here
            testCase.setCallerState('off');
            testCase.readMissingFile();
            testCase.verifyEqual(testCase.callerState(), 'off', ...
                'an erroring read must leave the caller''s OFF state off');
        end

        % ---- exit 3: no delimiter worked (error raised after the reset) --

        function testNoDelimiterErrorKeepsCallerOff(testCase)
            testCase.setCallerState('off');
            testCase.readAtDelimitedFile();
            testCase.verifyEqual(testCase.callerState(), 'off', ...
                ['G90: the "no delimiter worked" error is raised after the ' ...
                'unconditional reset, which turned the caller''s OFF state on']);
        end

        function testNoDelimiterErrorKeepsCallerOn(testCase)
            %control
            testCase.setCallerState('on');
            testCase.readAtDelimitedFile();
            testCase.verifyEqual(testCase.callerState(), 'on', ...
                'the "no delimiter worked" error must leave the caller''s ON state on');
        end

    end

    methods (Access = private)

        function tdir = tempFolder(testCase)
            %a fresh folder, removed after the test whatever its outcome
            tdir = tempname;
            mkdir(tdir);
            testCase.addTeardown(@() rmdir(tdir, 's'));
        end

        function file = writeFile(testCase, name, contents)
            %write CONTENTS verbatim to NAME in a fresh temporary folder
            file = fullfile(testCase.tempFolder(), name);
            fid = fopen(file, 'w');
            testCase.assertGreaterThan(fid, 0, sprintf('could not create %s', file));
            fprintf(fid, '%s', contents);
            fclose(fid);
        end

        function setCallerState(testCase, state)
            %the state a caller holds before calling psyrat_readtable
            warning(state, testCase.WarnId);
            testCase.assertEqual(testCase.callerState(), state, ...
                'premise: the caller''s warning state was not set');
        end

        function state = callerState(testCase)
            s = warning('query', testCase.WarnId);
            state = s.state;
        end

        function readMissingFile(testCase)
            %A genuine erroring read: the extension passes psyrat_readtable's
            %own checks, the warning is turned off, and readtable then throws
            %because the file does not exist. Deterministic on every
            %platform; the identifier is not pinned (it is MATLAB's, not ours).
            missing = fullfile(testCase.tempFolder(), 'no_such_file.csv');
            testCase.assertFalse(isfile(missing), 'premise: the file must not exist');
            err = [];
            try
                psyrat_readtable('file', missing);
            catch err
            end
            testCase.assertNotEmpty(err, 'premise: reading a missing file must throw');
            testCase.assertFalse(ismember(err.identifier, testCase.OwnIds), sprintf( ...
                ['premise: the error must come from readtable, after the ' ...
                'warning was turned off, not from psyrat_readtable''s own ' ...
                'checks (got ''%s'')'], err.identifier));
        end

        function readAtDelimitedFile(testCase)
            %Neither readtable's auto-detection nor the cycling list splits
            %on '@' (the premise TestReadtableDelimiterOption pins), so every
            %attempt yields one column and psyrat_readtable raises its own
            %error after the reset.
            file = testCase.writeFile('at_delimited.txt', ...
                sprintf('id@meas\ns1@1.1\ns2@2.2\n'));
            testCase.assertError(@() psyrat_readtable('file', file), ...
                'loadingfile:delimiter', ...
                'premise: the @-delimited file must exhaust the delimiter cycle');
        end

    end
end
