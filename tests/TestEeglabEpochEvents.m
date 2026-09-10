classdef TestEeglabEpochEvents < PsyRATTestBase
    % Unit tests for psyrat_eeglab_epoch_events, the EEGLAB time-locking-event
    % resolver. EEGLAB is not installed in this environment and the adapter's
    % read path requires pop_loadset, so the event-labeling logic is isolated
    % in this pure helper and exercised here with synthetic epoch structs.
    %
    % Behavior under test (owner decisions 2026-06-19): the time-locking event
    % is the eventtype whose eventlatency is closest to zero; when an epoch has
    % multiple events but no usable eventlatency, fall back to the first listed
    % event and count it as unresolved so the adapter can warn.

    methods (Test)

        function testPicksTimeLockingEventNotFirst(testCase)
            % A fixation event precedes the time-locking stimulus inside the
            % epoch window. The first listed event is 'fix' (legacy bug), but
            % the latency-0 event is 'stim'.
            ep = struct();
            ep.eventtype = {'fix','stim'};
            ep.eventlatency = {-150, 0};
            [event,nUn] = psyrat_eeglab_epoch_events(ep,1);
            testCase.verifyEqual(event{1},'stim');
            testCase.verifyEqual(nUn,0);
        end

        function testClosestToZeroNotExactZero(testCase)
            % Resampling leaves the time-locking event at a small non-zero
            % latency; closest-to-zero still selects it over a distant event.
            ep = struct();
            ep.eventtype = {'fix','stim'};
            ep.eventlatency = {-300, 0.4};
            [event,nUn] = psyrat_eeglab_epoch_events(ep,1);
            testCase.verifyEqual(event{1},'stim');
            testCase.verifyEqual(nUn,0);
        end

        function testSingleEventEpochUnambiguous(testCase)
            % One event in the epoch: taken directly, never counted unresolved.
            ep = struct();
            ep.eventtype = {'stim'};
            ep.eventlatency = {0};
            [event,nUn] = psyrat_eeglab_epoch_events(ep,1);
            testCase.verifyEqual(event{1},'stim');
            testCase.verifyEqual(nUn,0);
        end

        function testNonCellEventType(testCase)
            % EEGLAB sometimes stores a lone event type as a plain char.
            ep = struct();
            ep.eventtype = 'stim';
            ep.eventlatency = 0;
            [event,nUn] = psyrat_eeglab_epoch_events(ep,1);
            testCase.verifyEqual(event{1},'stim');
            testCase.verifyEqual(nUn,0);
        end

        function testMissingEventLatencyFallsBackAndCounts(testCase)
            % Multiple events but no eventlatency field at all: keep the legacy
            % first-event result and flag the epoch as unresolved.
            ep = struct('eventtype',{{'fix','stim'}}); % no eventlatency field
            testCase.verifyFalse(isfield(ep,'eventlatency'));
            [event,nUn] = psyrat_eeglab_epoch_events(ep,1);
            testCase.verifyEqual(event{1},'fix');
            testCase.verifyEqual(nUn,1);
        end

        function testMismatchedLatencyLengthFallsBack(testCase)
            % eventlatency length does not match eventtype: treated as unusable.
            ep = struct();
            ep.eventtype = {'a','b','c'};
            ep.eventlatency = {0, 100};
            [event,nUn] = psyrat_eeglab_epoch_events(ep,1);
            testCase.verifyEqual(event{1},'a');
            testCase.verifyEqual(nUn,1);
        end

        function testEmptyEventTypeYieldsBlank(testCase)
            % No event type recorded for the epoch -> empty event label.
            ep = struct();
            ep.eventtype = {};
            [event,nUn] = psyrat_eeglab_epoch_events(ep,1);
            testCase.verifyEqual(event{1},'');
            testCase.verifyEqual(nUn,0);
        end

        function testMultipleEpochsAggregateUnresolvedCount(testCase)
            % Mixed batch: one resolvable epoch and two without usable latency
            % ([] and empty-cell). Counts and labels are per epoch.
            ep(1).eventtype = {'fix','stim'}; ep(1).eventlatency = {-150, 0};
            ep(2).eventtype = {'a','b'};      ep(2).eventlatency = [];
            ep(3).eventtype = {'c','d'};      ep(3).eventlatency = {};
            [event,nUn] = psyrat_eeglab_epoch_events(ep,3);
            testCase.verifyEqual(event,{'stim';'a';'c'});
            testCase.verifyEqual(nUn,2);
        end

    end
end
