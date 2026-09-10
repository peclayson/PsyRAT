classdef TestRelsummaryInteractiveOption < PsyRATTestBase
    %G3 regression pins (owner ruling 2026-09-05, pulled forward from the
    %2026-08-24 council's post-release wave 1). psyrat_relsummary resolved
    %its interactive flag from psyrat_data.proc.interactive alone, and a
    %.psyrat saved by the GUI carried no such field (psyrat_startproc never
    %wrote it; only psyrat_run did), so the flag failed OPEN: a scripted
    %loop over GUI-produced files defaulted to interactive and could stall on
    %a modal dialog, once per dataset, with no argument that could say
    %otherwise. The fix adds an explicit 'interactive' name-value option that
    %wins over the stored field; a GUI-saved file still resolves to
    %interactive when the option is absent, as chapter 14 documents.
    %
    %WHY THESE TESTS DO NOT OBSERVE THE DIALOG. The threshold-failure dialog
    %("Data do not meet reliability threshold") is raised only by the ten
    %inline relerr.nogooddata guards, and since the G39 strict fold every one
    %of them is expected dead for nevents >= 1 (the profiler record inside
    %psyrat_relsummary, "SUPERSEDED 2026-08-29 (G39 strict fold)" and the
    %paragraphs around it, works this through: the fold killed the six guards
    %that were live, two were already dead by arm, and two joined the G54
    %never-firing backstop class): when nothing clears the cutoff the group
    %badids hold every participant, the fallback table is non-empty, and the
    %guard never fires.
    %No committed fixture reaches a live notice, so a test that asserted "no
    %dialog appeared" would pass whether or not the option worked (the
    %differential-assertion rule). The observable pinned here is therefore the
    %RESOLUTION of the flag, read back from the function's own validation and
    %from a source-level pin of the precedence chain.
    %
    %RED against the pre-G3 code: the validation test (no 'interactive'
    %option existed, so a bad value was silently ignored and no error with the
    %pinned identifier was raised) and the precedence pin (the option was not
    %consulted). The accepted-values test is the positive control and passes
    %on both builds.

    methods (Test)

        function testRejectsANonScalarOrNonLogicalOption(testCase)
            %The option is validated the way the file's other options are: a
            %char, a vector, a complex or a NaN value is refused with the
            %named identifier rather than silently coerced.
            pd = localGuiShapedFixture();
            bad = {'yes', [true false], 1+1i, NaN};
            for k = 1:numel(bad)
                testCase.verifyError(@() psyrat_relsummary('psyrat_data', pd, ...
                    'analysis', 'sing', 'depcutoff', 0.30, 'meascutoff', 2, ...
                    'depcentmeas', 1, 'CI', 0.95, 'interactive', bad{k}), ...
                    'varargin:interactive');
            end
        end

        function testAcceptedOptionValuesRunToCompletion(testCase)
            %Positive control for the validation test: logical and numeric
            %scalars are accepted in both directions and the summary runs.
            pd = localGuiShapedFixture();
            good = {true, false, 1, 0};
            for k = 1:numel(good)
                [out, relerr] = psyrat_relsummary('psyrat_data', pd, ...
                    'analysis', 'sing', 'depcutoff', 0.30, 'meascutoff', 2, ...
                    'depcentmeas', 1, 'CI', 0.95, 'interactive', good{k});
                testCase.verifyTrue(isfield(out, 'relsummary'));
                testCase.verifyEqual(relerr.nogooddata, 0);
            end
        end

        function testResolutionOrderIsPinnedInSource(testCase)
            %Source-level pin of the precedence chain: explicit option, then
            %psyrat_data.proc.interactive, then the interactive default. The
            %three arms are matched on their own code lines (continuations
            %joined first, per the repo's regexp lesson), so a reorder or a
            %dropped arm fails here even though no fixture can make the
            %dialog observable.
            src = fileread(which('psyrat_relsummary'));
            src = regexprep(src, '\.\.\.\s*[\r\n]+\s*', ' ');
            optIdx = regexp(src, 'find\(strcmpi\(''interactive'',varargin\),1\)', 'once');
            fldIdx = regexp(src, 'isfield\(psyrat_data\.proc,''interactive''\)', 'once');
            defIdx = regexp(src, 'else\s+interactive\s*=\s*true;', 'once');
            testCase.assertNotEmpty(optIdx, 'the explicit option must be read from varargin');
            testCase.assertNotEmpty(fldIdx, 'the stored proc.interactive field must still be consulted');
            testCase.assertNotEmpty(defIdx, 'the interactive default must be the last resort');
            testCase.verifyLessThan(optIdx, fldIdx, ...
                'G3: the explicit option must be resolved before the stored field');
            testCase.verifyLessThan(fldIdx, defIdx, ...
                'the stored field must be resolved before the default');
            %the stored-field arm must be guarded so it never overrides an
            %explicit option
            testCase.verifyTrue(~isempty(regexp(src, ...
                'if\s+isempty\(interactive\)\s+if\s+isfield\(psyrat_data,''proc''\)', 'once')), ...
                'the stored field may only fill an unset option');
        end

    end
end

function pd = localGuiShapedFixture()
%The multi-group x event one-facet fixture (relsummary 'sing' case 4), with
%any proc.interactive field removed so the struct is shaped like a .psyrat the
%GUI saved before G3. The factory's own docstring records the label contract
%this fixture pins; nothing about it is changed here beyond the field removal.
pd = PsyRATTestDataFactory.makeICRelDataMultiGroupEvent();
if isfield(pd, 'proc') && isstruct(pd.proc) && isfield(pd.proc, 'interactive')
    pd.proc = rmfield(pd.proc, 'interactive');
end
end
