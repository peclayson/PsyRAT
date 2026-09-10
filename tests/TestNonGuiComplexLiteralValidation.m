classdef TestNonGuiComplexLiteralValidation < PsyRATTestBase
    % Regression guard for the complex-literal sweep OUTSIDE the viewers.
    %
    % Sibling of TestViewerComplexLiteralValidation, which covers the fourteen
    % GUI conditions. This file covers the four sites that sweep left alone
    % pending owner rulings, all ruled FIX on 2026-08-12:
    %
    %   psyrat_run                    adapt_delta, max_treedepth
    %   psyrat_computevarcomp         seven sampler guards in the validation
    %                                 prologue (chains, legacy iter, warmup,
    %                                 sampling, seed, adapt_delta, max_treedepth)
    %   psyrat_scoring_spec_create    window bounds, neighborhood_samples,
    %                                 area_fraction
    %   psyrat_scoring_spec_from_inputs   three builder parses (neighborhood,
    %                                 area_fraction, channel selector)
    %   psyrat_trt_reloverallt        the recalculation trial count n'
    %
    % The legacy-iter and channel-selector guards were MISSED by the first pass
    % and found by an adversarial review of it; both are recorded here because
    % "the sweep is complete" is exactly the claim that keeps being wrong.
    %
    % WHY THEY ALL FAILED THE SAME WAY. See the sibling's header for the full
    % idiom table. The short version: str2double parses complex literals, and
    % isnan is false for one, isfinite is TRUE, round/floor act on BOTH parts so
    % round(x) == x and floor(x) == x hold, max/min compare by MAGNITUDE, and
    % every ordering comparison reads the REAL part alone. mod() is the lone
    % idiom that rejects a complex, and it does so by ERRORING
    % (MATLAB:math:mustBeReal), which is why ~isreal has to lead or sit early in
    % every condition rather than trail it.
    %
    % THREE OF THESE PATHS WERE SILENT rather than merely ugly:
    %
    %   adapt_delta      psyrat_run's guard, psyrat_computevarcomp's guard AND
    %                    MatlabStan's validateattributes all passed a complex
    %                    value. Measured 2026-08-12: the unpatched engine
    %                    completed a CmdStan fit instead of erroring. What
    %                    CmdStan made of the token delta=0.5+0.5i was NOT
    %                    established, and no claim here rests on it.
    %   legacy iter      cleared the only guard it gets, then produced complex
    %                    warmup/sampling counts via max(round(iter),2) that the
    %                    warmup/sampling guards never re-check because those sit
    %                    in the branch not taken.
    %   scoring window   psyrat_scoring_spec_create BUILT a spec from
    %                    window_start_ms = 100+5i, stored the complex value and
    %                    folded it into spec_hash (measured). Read but not
    %                    measured: psyrat_project_score's resolver is
    %                    find(time_ms >= start_ms & time_ms <= end_ms) -- also
    %                    real-part only -- so scoring would proceed on the
    %                    real-part window, i.e. scores recorded under a window
    %                    that is not the window stored.
    %
    % The other guards did fail before the fix, but inside mod() with a MATLAB
    % internal message instead of the toolbox's own "How to fix" text.
    %
    % STRATEGY, same as the sibling: every check is DIFFERENTIAL rather than
    % message-matching. Driving a field with a complex value must produce the
    % same observable as driving it with plainly invalid input, which is the
    % contract these guards already document. Each differential pair carries a
    % POSITIVE CONTROL asserting a valid entry produces the OPPOSITE observable
    % -- without one, a differential assertion can pass for the wrong reason.
    %
    % NO CmdStan IS INVOLVED, by two different mechanisms. The ENGINE tests pass
    % a table with no 'meas' column, so a run that clears the sampler prologue
    % stops at varargin:colheaders instead of proceeding to a fit -- and that
    % identifier IS their positive control, since it proves execution reached
    % data validation, which sits after the guards under test. The psyrat_run
    % tests are protected differently: psyrat_run raises psyrat_run:savepath
    % before it examines the table at all, which is why THEIR positive control
    % asserts that exact identifier.

    methods (Test)

        %------------------------------------------------------------------
        % psyrat_run -- the CLI twin of the silent adapt_delta path
        %------------------------------------------------------------------

        function testRunAdaptDeltaRejectsComplexLikeOutOfRange(testCase)
            T = localMinimalTable();

            complexId = localErrId(@() psyrat_run('data',T,'adapt_delta',0.5+0.5i));
            invalidId = localErrId(@() psyrat_run('data',T,'adapt_delta',5));
            validId   = localErrId(@() psyrat_run('data',T,'adapt_delta',0.95));

            %POSITIVE CONTROL: the observable must be able to say "accepted".
            %Assert the EXACT downstream identifier rather than merely "not
            %this field's" -- assertNotEqual would also pass for a failure
            %raised BEFORE the sampler block (a bad data source, an inputParser
            %rejection), which would prove nothing about this guard.
            testCase.assertEqual(validId, 'psyrat_run:savepath', ...
                ['positive control: a valid adapt_delta must clear the ' ...
                'sampler-override block and stop later, at the savepath check.']);

            testCase.verifyEqual(complexId, invalidId, ...
                ['a complex adapt_delta must be rejected exactly like an ' ...
                'out-of-range one.']);
            testCase.verifyEqual(complexId, 'psyrat_run:adapt_delta', ...
                ['0.5+0.5i must raise psyrat_run:adapt_delta. Before the fix ' ...
                'it cleared this guard entirely and reached CmdStan as the ' ...
                'token delta=0.5+0.5i.']);
        end

        function testRunMaxTreedepthRejectsComplexLikeNonInteger(testCase)
            T = localMinimalTable();

            complexId = localErrId(@() psyrat_run('data',T,'max_treedepth',12+1i));
            invalidId = localErrId(@() psyrat_run('data',T,'max_treedepth',12.5));
            validId   = localErrId(@() psyrat_run('data',T,'max_treedepth',12));

            testCase.assertEqual(validId, 'psyrat_run:savepath', ...
                ['positive control: a valid max_treedepth must clear this ' ...
                'block and stop later, at the savepath check.']);

            testCase.verifyEqual(complexId, invalidId, ...
                ['a complex max_treedepth must be rejected exactly like a ' ...
                'non-integer one.']);
            testCase.verifyEqual(complexId, 'psyrat_run:max_treedepth', ...
                ['12+1i must raise psyrat_run:max_treedepth rather than dying ' ...
                'inside mod() with MATLAB:math:mustBeReal.']);
        end

        %------------------------------------------------------------------
        % psyrat_computevarcomp -- all six guards in the validation prologue
        %------------------------------------------------------------------

        function testEngineSamplerGuardsRejectComplexLikeInvalid(testCase)
            %One case per guard: {field, complex value, plainly-invalid value,
            %expected identifier}. The invalid comparator is chosen to fail the
            %SAME guard for an ordinary reason.
            cases = { ...
                'chains',        4+1i,     0,    'varargin:nchains'; ...
                'warmup',        100+1i,  -1,    'varargin:warmup'; ...
                'sampling',      100+1i,   0,    'varargin:sampling'; ...
                'seed',          5+1i,    -1,    'varargin:seed'; ...
                'adapt_delta',   0.5+0.5i, 5,    'varargin:adapt_delta'; ...
                'max_treedepth', 12+1i,    12.5, 'varargin:max_treedepth'};

            %POSITIVE CONTROL FIRST. Every sampler setting valid, but the table
            %has no 'meas' column: execution must clear all six guards and stop
            %at data validation. If this ever changes, the differentials below
            %stop meaning "the guard fired" and the file needs revisiting.
            validId = localEngineErrId('adapt_delta',0.95,'max_treedepth',12,'seed',7);
            testCase.assertEqual(validId, 'varargin:colheaders', ...
                ['positive control: with all sampler settings valid the run ' ...
                'must reach data validation, proving these guards can pass.']);

            for k = 1:size(cases,1)
                field = cases{k,1};
                complexId = localEngineErrId(field, cases{k,2});
                invalidId = localEngineErrId(field, cases{k,3});

                testCase.verifyEqual(complexId, invalidId, ...
                    sprintf(['a complex ''%s'' must be rejected exactly like ' ...
                    'a plainly invalid one.'], field));
                testCase.verifyEqual(complexId, cases{k,4}, ...
                    sprintf(['a complex ''%s'' must raise %s. Before the fix ' ...
                    'this was MATLAB:math:mustBeReal out of mod() for five of ' ...
                    'the six, and for adapt_delta no error at all.'], ...
                    field, cases{k,4}));
            end
        end

        function testEngineLegacyIterRejectsComplexLikeNonPositive(testCase)
            %The SEVENTH guard in the prologue, and the second silent one. It
            %sits in the legacy 'iter' branch, which is taken INSTEAD of the
            %warmup/sampling branch -- so those two guards never see these
            %values and this is the only check they get. A complex iter then
            %flows through niter = max(round(iter),2), where round acts on both
            %parts and max compares by magnitude, and both derived counts reach
            %MatlabStan complex.
            %
            %Missed by the first pass of this increment and found by adversarial
            %review; the arguments differ from the other engine cases because
            %this branch is only reachable when warmup and sampling are BOTH
            %absent.
            bad = table({'1';'2'}, 'VariableNames', {'id'});
            complexId = localErrId(@() psyrat_computevarcomp( ...
                'data',bad,'chains',4,'iter',10000+5i));
            invalidId = localErrId(@() psyrat_computevarcomp( ...
                'data',bad,'chains',4,'iter',0));
            validId   = localErrId(@() psyrat_computevarcomp( ...
                'data',bad,'chains',4,'iter',10000));

            testCase.assertEqual(validId, 'varargin:colheaders', ...
                ['positive control: a valid legacy iter must clear this guard ' ...
                'and reach data validation.']);

            testCase.verifyEqual(complexId, invalidId, ...
                ['a complex legacy ''iter'' must be rejected exactly like a ' ...
                'non-positive one.']);
            testCase.verifyEqual(complexId, 'varargin:niter', ...
                ['10000+5i must raise varargin:niter. Before the fix it ' ...
                'cleared this guard entirely -- isfinite is TRUE and > 0 reads ' ...
                'the real part -- and produced complex warmup and sampling ' ...
                'counts that no later guard re-checks.']);
        end

        %------------------------------------------------------------------
        % psyrat_scoring_spec_create -- the second silent path
        %------------------------------------------------------------------

        function testScoringWindowRejectsComplexLikeNaN(testCase)
            base = localScoringBase('time_window_mean_amplitude');

            complexId = localErrId(@() psyrat_scoring_spec_create(base{:}, ...
                'window_start_ms',100+5i,'window_end_ms',300));
            invalidId = localErrId(@() psyrat_scoring_spec_create(base{:}, ...
                'window_start_ms',NaN,'window_end_ms',300));
            validId   = localErrId(@() psyrat_scoring_spec_create(base{:}, ...
                'window_start_ms',100,'window_end_ms',300));

            testCase.assertEmpty(validId, ...
                'positive control: a real 100..300 ms window must be accepted.');

            testCase.verifyEqual(complexId, invalidId, ...
                'a complex window bound must be rejected exactly like NaN.');
            testCase.verifyEqual(complexId, 'psyrat_scoring:window', ...
                ['100+5i must raise psyrat_scoring:window. Before the fix the ' ...
                'spec was BUILT with the complex bound, scoring ran on its ' ...
                'real part, and the value reached spec_hash and the ' ...
                'provenance header.']);
        end

        function testScoringNeighborhoodRejectsComplexLikeZero(testCase)
            base = localScoringBase('local_peak_amplitude');
            base = [base {'window_start_ms',100,'window_end_ms',300}];

            complexId = localErrId(@() psyrat_scoring_spec_create(base{:}, ...
                'method_params',struct('neighborhood_samples',3+4i)));
            invalidId = localErrId(@() psyrat_scoring_spec_create(base{:}, ...
                'method_params',struct('neighborhood_samples',0)));
            validId   = localErrId(@() psyrat_scoring_spec_create(base{:}, ...
                'method_params',struct('neighborhood_samples',3)));

            testCase.assertEmpty(validId, ...
                'positive control: a neighborhood of 3 samples must be accepted.');

            testCase.verifyEqual(complexId, invalidId, ...
                ['a complex neighborhood_samples must be rejected exactly ' ...
                'like a non-positive one.']);
            testCase.verifyEqual(complexId, 'psyrat_scoring:localPeakParams', ...
                ['3+4i must be rejected: floor() acts on both parts, so the ' ...
                'floor(n) ~= n integrality test passed it unchanged.']);
        end

        function testScoringAreaFractionRejectsComplexLikeOutOfRange(testCase)
            base = localScoringBase('fractional_area_latency');
            base = [base {'window_start_ms',100,'window_end_ms',300}];

            complexId = localErrId(@() psyrat_scoring_spec_create(base{:}, ...
                'method_params',struct('area_fraction',0.5+0.1i)));
            invalidId = localErrId(@() psyrat_scoring_spec_create(base{:}, ...
                'method_params',struct('area_fraction',1.5)));
            validId   = localErrId(@() psyrat_scoring_spec_create(base{:}, ...
                'method_params',struct('area_fraction',0.5)));

            testCase.assertEmpty(validId, ...
                'positive control: an area fraction of 0.5 must be accepted.');

            testCase.verifyEqual(complexId, invalidId, ...
                ['a complex area_fraction must be rejected exactly like an ' ...
                'out-of-range one.']);
            testCase.verifyEqual(complexId, 'psyrat_scoring:areaLatencyParams', ...
                '0.5+0.1i must raise psyrat_scoring:areaLatencyParams.');
        end

        function testScoringBuilderFallsBackOnComplexLikeGarbage(testCase)
            %The GUI-facing builder has its own clamps ahead of the validator,
            %and both were blind. Their documented contract is a silent fallback
            %to the conventional default, so the observable here is the VALUE
            %the spec ends up with, not an error.
            %The complex literals here have a real part that DIFFERS from the
            %fallback default (7 against the default 3; 0.7 against 0.5). That
            %is deliberate and load-bearing: an earlier draft used '3-4i' and
            %'0.5+0.1i', whose real parts ARE the defaults, so an implementation
            %that silently coerced with real(n) instead of rejecting would have
            %produced identical observables and both assertions would have
            %passed. Real-part coercion is the most likely wrong fix for this
            %bug class, so the fixture has to be able to see it.
            complexSpec = psyrat_scoring_spec_from_inputs( ...
                'local_peak_amplitude', localBuilderInputs('neighborhood','7-4i'));
            garbageSpec = psyrat_scoring_spec_from_inputs( ...
                'local_peak_amplitude', localBuilderInputs('neighborhood','abc'));
            validSpec   = psyrat_scoring_spec_from_inputs( ...
                'local_peak_amplitude', localBuilderInputs('neighborhood','5'));

            testCase.assertEqual(validSpec.method_params.neighborhood_samples, 5, ...
                ['positive control: a valid neighborhood must be carried ' ...
                'through, not replaced by the default.']);

            testCase.verifyEqual(complexSpec.method_params.neighborhood_samples, ...
                garbageSpec.method_params.neighborhood_samples, ...
                ['complex neighborhood text must fall back exactly like any ' ...
                'other unparseable text.']);
            testCase.verifyEqual(complexSpec.method_params.neighborhood_samples, 3, ...
                ['''7-4i'' must fall back to the default 3, NOT be coerced to ' ...
                'its real part 7.']);
            testCase.verifyTrue(isreal(complexSpec.method_params.neighborhood_samples), ...
                ['the stored neighborhood must be real. Before the fix ''7-4i'' ' ...
                'survived round(), isfinite() and n < 1 and reached the spec.']);

            complexArea = psyrat_scoring_spec_from_inputs( ...
                'fractional_area_latency', localBuilderInputs('area_fraction','0.7+0.1i'));
            garbageArea = psyrat_scoring_spec_from_inputs( ...
                'fractional_area_latency', localBuilderInputs('area_fraction','abc'));
            validArea   = psyrat_scoring_spec_from_inputs( ...
                'fractional_area_latency', localBuilderInputs('area_fraction','0.7'));

            testCase.assertEqual(validArea.method_params.area_fraction, 0.7, ...
                'positive control: a valid area fraction must be carried through.');

            testCase.verifyEqual(complexArea.method_params.area_fraction, ...
                garbageArea.method_params.area_fraction, ...
                'complex area-fraction text must fall back like any other.');
            testCase.verifyEqual(complexArea.method_params.area_fraction, 0.5, ...
                ['''0.7+0.1i'' must take the 0.5 default, NOT be coerced to ' ...
                'its real part 0.7.']);
            testCase.verifyTrue(isreal(complexArea.method_params.area_fraction), ...
                'the stored area fraction must be real.');
        end

        function testScoringChannelSelectorRejectsComplexLikeLabels(testCase)
            %A third numeric parse in the same builder, found by adversarial
            %review after the first pass swept only the two method-parameter
            %clamps. A complex channel INDEX would have been carried through as
            %the selector, and psyrat_scoring_spec_create only checks that the
            %selector is non-empty. The documented fallback for unparseable
            %channel text is to treat it as channel LABELS, so that is the
            %observable.
            complexSpec = psyrat_scoring_spec_from_inputs( ...
                'time_window_mean_amplitude', localBuilderInputs('channel','1+2i'));
            garbageSpec = psyrat_scoring_spec_from_inputs( ...
                'time_window_mean_amplitude', localBuilderInputs('channel','Cz'));
            validSpec = psyrat_scoring_spec_from_inputs( ...
                'time_window_mean_amplitude', localBuilderInputs('channel','3'));

            testCase.assertEqual(validSpec.channel_selector, 3, ...
                ['positive control: a real numeric channel must stay a numeric ' ...
                'index, which is what makes "became labels" a distinguishable ' ...
                'outcome at all.']);

            testCase.verifyEqual(class(complexSpec.channel_selector), ...
                class(garbageSpec.channel_selector), ...
                ['a complex channel index must land in the same branch as any ' ...
                'other unparseable channel text.']);
            testCase.verifyFalse(isnumeric(complexSpec.channel_selector), ...
                ['''1+2i'' must not survive as a numeric channel index: ' ...
                'psyrat_scoring_spec_create only checks the selector is ' ...
                'non-empty, so a complex index would reach scoring.']);
        end

        %------------------------------------------------------------------
        % psyrat_trt_reloverallt -- the recalculation trial count n'
        %------------------------------------------------------------------

        function testTrtRecalcTrialCountRejectsComplexLikeGarbage(testCase)
            [complexShown, complexThrew] = localRunTrtNtrls(testCase, '3-4i');
            garbageShown = localRunTrtNtrls(testCase, 'abc');
            [validShown,  validThrew]    = localRunTrtNtrls(testCase, '30');

            %POSITIVE CONTROL: a valid count must NOT raise the dialog, and must
            %get far enough to hit the fixture's own limit. Asserting the throw
            %as well as the absence of the dialog is what proves the accept path
            %actually ran rather than dying somewhere earlier.
            testCase.assertFalse(validShown, ...
                'positive control: a valid trial count must not raise the dialog.');
            testCase.assertEqual(validThrew, 'MATLAB:nonExistentField', ...
                ['positive control: an accepted count must proceed past the ' ...
                'guard and reach REL.out, which this summary fixture lacks.']);

            testCase.verifyEqual(complexShown, garbageShown, ...
                ['a complex trial count must raise the same dialog as any ' ...
                'other invalid text.']);
            testCase.verifyTrue(complexShown, ...
                ['''3-4i'' must be refused. Before this guard existed it ' ...
                'reached psyrat_rel_trt as the error-term divisor and crashed ' ...
                'in quantile(), in a file with no try/catch.']);
            testCase.verifyEmpty(complexThrew, ...
                ['a refused count must return cleanly, not throw. An empty ' ...
                'throw is what distinguishes "the guard rejected it" from ' ...
                '"the callback died before raising the dialog".']);
        end

        function testTrtRecalcTrialCountRefusesNonPositiveAndFractional(testCase)
            %The ruled policy, not merely the isreal addition: n' is a count of
            %measurement replications entering a mean, so zero, negative and
            %fractional values are all outside its domain.
            [validShown, validThrew] = localRunTrtNtrls(testCase, '30');
            testCase.assertFalse(validShown, ...
                'positive control: a valid trial count must not raise the dialog.');
            testCase.assertEqual(validThrew, 'MATLAB:nonExistentField', ...
                'positive control: an accepted count must proceed past the guard.');

            %Each case asserts BOTH that the dialog was raised and that nothing
            %was thrown, so a callback that died before errordlg cannot pass.
            localVerifyRefused(testCase, '', ...
                'a blank trial count must be refused, not silently read as NaN.');
            localVerifyRefused(testCase, '0', ...
                ['zero must be refused. For reltype 1 and 3 it makes the error ' ...
                'term Inf and the coefficient prints as exactly 0.00; for ' ...
                'reltype 2 obs also divides the universe score, giving NaN. ' ...
                'Both are wrong, and the first is quietly plausible.']);
            localVerifyRefused(testCase, '-5', ...
                ['a negative count must be refused: it can push the ' ...
                'coefficient outside [0,1].']);
            localVerifyRefused(testCase, '2.5', ...
                ['a fractional count must be refused: you cannot average 2.5 ' ...
                'trials, and the figure title formats it with %d.']);

            %The floor is 1. The vs-trials PLOT fields reject 1 because their
            %value is a plot RANGE endpoint and a single point is a degenerate
            %line; obs here is a scalar, so a single-trial composite is a
            %legitimate D-study target. Asserting the throw as well as the
            %absent dialog is what makes this an ACCEPT rather than "nothing
            %visible happened".
            [oneShown, oneThrew] = localRunTrtNtrls(testCase, '1');
            testCase.verifyFalse(oneShown, ...
                'a single-trial composite must be accepted; the floor is 1.');
            testCase.verifyEqual(oneThrew, 'MATLAB:nonExistentField', ...
                'a single-trial composite must proceed past the guard.');
        end

    end
end


%==========================================================================
% helpers
%==========================================================================

function localVerifyRefused(testCase, text, why)
%A refused entry must raise the dialog AND return cleanly. Asserting only the
%dialog would let a callback that throws before errordlg pass as a rejection.
[shown, threw] = localRunTrtNtrls(testCase, text);
testCase.verifyTrue(shown, why);
testCase.verifyEmpty(threw, ...
    sprintf('refusing ''%s'' must return cleanly, not throw.', text));
end


function T = localMinimalTable()
%A well-formed table for the psyrat_run tests. Its CONTENTS are never examined:
%psyrat_run validates the sampler overrides and then the savepath, both before
%any column checking, which is why these tests cannot start a fit either.
T = table({'1';'1';'2';'2'}, [1.1;1.2;2.1;2.2], 'VariableNames', {'id','meas'});
end


function id = localEngineErrId(varargin)
%Call psyrat_computevarcomp with valid sampler defaults, overridden by the
%name/value pairs passed in, and return the error identifier.
%
%The table deliberately has NO 'meas' column: a call that clears the sampler
%prologue then stops at varargin:colheaders instead of starting a fit, so these
%tests never touch CmdStan. The argument list is expanded through a local
%variable because MATLAB cannot index the result of a function call.
bad = table({'1';'2'}, 'VariableNames', {'id'});
args = {'data',bad,'chains',4,'warmup',100,'sampling',100};
for i = 1:2:numel(varargin)
    ind = find(strcmp(varargin{i},args),1);
    if isempty(ind)
        args = [args varargin(i) varargin(i+1)]; %#ok<AGROW>
    else
        args{ind+1} = varargin{i+1};
    end
end
id = localErrId(@() psyrat_computevarcomp(args{:}));
end


function base = localScoringBase(method)
base = {'name','fixture','method',method, ...
    'channel_selector',{'Cz'},'polarity','positive'};
end


function in = localBuilderInputs(field,value)
%Raw GUI-style string inputs for psyrat_scoring_spec_from_inputs. Everything
%except the field under test is valid, so the spec builds.
in = struct('name','fixture','window_start_ms','100','window_end_ms','300', ...
    'channel','Cz','polarity','positive');
in.(field) = value;
end


function id = localErrId(fh)
%Run fh and return the error identifier, or '' when it does not error.
id = '';
try
    fh();
catch ME
    id = ME.identifier;
end
end


function [shown, threwId] = localRunTrtNtrls(testCase, text)
%Drive the "Estimate New Reliability Coefficient" flow in
%psyrat_trt_reloverallt and report whether the error dialog was raised.
%
%Two steps: the overall-table screen carries a button that opens a "Specify
%Inputs" screen, and that screen carries the trial-count edit box and its own
%estimate button, whose callback is the local psyrat_shownewrel where the guard
%lives.
%
%On the ACCEPT path this fixture cannot complete the recalculation -- it has no
%REL.out -- so the callback throws MATLAB:nonExistentField well past the guard.
%That throw is caught rather than swallowed: it is REPORTED back as threwId so
%the caller can assert on it. That distinction matters. Swallowing it silently
%would let a reject-path bug that throws the SAME identifier before reaching
%errordlg (say, a bad field reference while building the message) masquerade as
%"no dialog, therefore accepted", and the floor-of-1 assertion would then report
%acceptance for a path that actually rejected. Callers assert threwId is EMPTY
%on every reject case, which closes that hole.
close('all','force');
cleanup = onCleanup(@() close('all','force')); %#ok<NASGU>

pd = PsyRATTestDataFactory.makeTRTSummaryData();
%The recalculation button is gated on the test-retest reliability types, and
%the shared fixture is a one-facet 'ic' summary; these two fields are what the
%"Specify Inputs" screen prints in its header.
pd.relsummary.reltype_name = 'trt';
pd.rel.filename = 'fixture.psyrat';
pd.proc.measheader = 'meas';

psyrat_trt_reloverallt('psyrat_data', pd, 'gui', 1);

openbtn = localEstimateButton(testCase, []);
localFire(openbtn);

edits = findall(groot,'Style','edit');
testCase.assertNotEmpty(edits, 'the Specify Inputs screen has no edit field.');
set(edits(1), 'String', text);

estbtn = localEstimateButton(testCase, openbtn);
threwId = '';
try
    localFire(estbtn);
catch ME
    if ~strcmp(ME.identifier,'MATLAB:nonExistentField')
        rethrow(ME);
    end
    threwId = ME.identifier;
end

%B21 stage 1: the rejection dialog now carries its purpose title and must be
%modal (the untitled non-modal form was invisible behind the Specify Inputs
%figure -- the census passed while the live user saw nothing, which is why
%this census alone could not catch B21).
dlg = findall(groot,'Type','figure','Name','Invalid trial count');
shown = ~isempty(dlg);
if shown
    testCase.verifyEqual(dlg(1).WindowStyle, 'modal', ...
        'the trial-count rejection dialog must be modal (B21 stage 1)');
end
end


function btn = localEstimateButton(testCase, exclude)
%Both screens label their button {'Estimate New';'Reliability Coefficient'};
%pass the first one as `exclude` to reach the second.
btns = findall(groot,'Style','pushbutton');
btn = [];
for k = 1:numel(btns)
    s = btns(k).String;
    if iscell(s); s = strjoin(s,' '); end
    if ~contains(char(s),'Estimate New'); continue; end
    if ~isempty(exclude) && isequal(btns(k),exclude); continue; end
    btn = btns(k);
    return;
end
testCase.assertNotEmpty(btn, 'the "Estimate New" button was not found.');
end


function localFire(h)
%These buttons wire their callbacks as cell arrays ({@fcn, psyrat_data, ...}),
%so the trailing entries have to be unpacked as extra arguments.
cb = h.Callback;
if iscell(cb)
    fcn = cb{1};
    fcn(h, [], cb{2:end});
else
    cb(h, []);
end
end
