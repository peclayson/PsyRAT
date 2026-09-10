classdef TestCiCoverageConsistency < PsyRATTestBase
    %The credible-interval width a user requests must govern EVERY interval in
    %the same output, and the table must say which width that is.
    %
    %THE DEFECT (council 2026-08-24, numerical lane finding 2, widened here).
    %The manual's scripting chapter (documentation/manual/14_scripting.md,
    %section 14.2) promises "CI (default .95) is accepted by every route". That
    %promise lived in documentation/scripted_use.md when this class was written
    %and moved to the manual on 2026-09-04. psyrat_relsummary does parse it and thread it into the
    %reliability coefficients as 'CI',ciperc. But the between-/within-person
    %standard deviations were read off literal .025/.975 quantiles at ~25 sites,
    %and psyrat_variancet's figure title hardcoded the string "95%". Meanwhile
    %the SEM (psyrat_sem_stats) and the ICC both DID honor ciperc. So a scripted
    %run with 'CI',.99 produced one exported table, under a title asserting 95%,
    %whose three column families carried two different coverages -- SDs at 95%,
    %SEMs and ICCs at 99%, and the title wrong for two of the three.
    %
    %Not reachable from the GUI: psyrat_relfigures contains no 'CI' argument, so
    %the GUI path always takes the .95 default. It is reachable from the
    %documented psyrat_run -> psyrat_relsummary -> psyrat_report workflow.
    %
    %WHY THE FIXTURE CARRIES 400 DRAWS. With a handful of draws, quantile()
    %clamps everything below the first plotting position to the minimum, so
    %quantile(x,.005) and quantile(x,.025) return the SAME number and the
    %widening assertion below could not fail no matter what the code did. 400
    %strictly increasing draws put both probabilities inside the interpolation
    %range. testTheFixtureCanActuallyDetectTheDefect is the guard on that.

    methods (Test)

        function testComponentSdsUseTheToolboxWideEdgeConvention(testCase)
            %The regression half, in two parts.
            %
            %PART 1 -- exact. The component SDs must now derive their edge the
            %same way every other interval in the toolbox does
            %(psyrat_rel_sing.m:109, psyrat_cutscore_sing.m:114 and
            %PsyRATAccuracyOracle all spell it (1-ci)/2). This is an equality,
            %not a tolerance: the whole point is that one convention governs.
            %
            %PART 2 -- bounded drift, stated rather than hidden. (1-.95)/2 is 6
            %ULP above the literal .025 the shipped code used, so default output
            %is NOT bit-identical; a bound can move in the last place. Pinning
            %the SIZE of that movement is what makes the change reviewable -- an
            %assertion of exact equality with the old literal would simply be
            %false, and no assertion at all would let real drift through.
            pd  = localFixture();
            out = localSummarize(pd, 0.95);
            ev  = out.relsummary.group(1).event(1);

            sig_u = pd.rel.out.sig_u;
            sig_e = pd.rel.out.sig_e;
            edge  = (1-0.95)/2;

            testCase.verifyEqual(ev.betsd.ll, quantile(sig_u, edge), ...
                ['Between-person lower bound does not match the toolbox edge '...
                'convention (1-ci)/2 that psyrat_rel_sing already uses.']);
            testCase.verifyEqual(ev.betsd.ul, quantile(sig_u, 1-edge), ...
                'Between-person upper bound does not match the convention.');
            testCase.verifyEqual(ev.witsd.ll, quantile(sig_e, edge), ...
                'Within-person lower bound does not match the convention.');
            testCase.verifyEqual(ev.witsd.ul, quantile(sig_e, 1-edge), ...
                'Within-person upper bound does not match the convention.');

            %Part 2: within a few ULP of the literals this replaced.
            legacy = { ...
                ev.betsd.ll, quantile(sig_u, .025), 'betsd.ll'
                ev.betsd.ul, quantile(sig_u, .975), 'betsd.ul'
                ev.witsd.ll, quantile(sig_e, .025), 'witsd.ll'
                ev.witsd.ul, quantile(sig_e, .975), 'witsd.ul'};
            for k = 1:size(legacy,1)
                got = legacy{k,1}; was = legacy{k,2};
                testCase.verifyLessThanOrEqual(abs(got-was), 4*eps(was), ...
                    sprintf(['%s moved more than 4 ULP from the value the '...
                    'literal .025/.975 quantiles produced. The default path '...
                    'is allowed to drift in the last place, not to change.'], ...
                    legacy{k,3}));
            end
        end

        function testComponentSdsFollowTheRequestedCi(testCase)
            %The differential half: the SAME fixture at two CI widths. A wider
            %requested interval must widen the component SDs, exactly as it
            %already widened the coefficients.
            pd = localFixture();

            narrow = localSummarize(pd, 0.95);
            wide   = localSummarize(pd, 0.99);

            evN = narrow.relsummary.group(1).event(1);
            evW = wide.relsummary.group(1).event(1);

            testCase.verifyGreaterThan(localWidth(evW.betsd), localWidth(evN.betsd), ...
                ['A 99% request must produce a WIDER between-person SD '...
                'interval than a 95% request. Equality means the component '...
                'quantiles are still pinned to a fixed 95% while the '...
                'coefficients follow the user (council numerical finding 2).']);
            testCase.verifyGreaterThan(localWidth(evW.witsd), localWidth(evN.witsd), ...
                'A 99% request must widen the within-person SD interval too.');
        end

        function testTheFixtureCanActuallyDetectTheDefect(testCase)
            %POSITIVE CONTROL. The assertions above are differential, so they
            %are only meaningful if this fixture is capable of showing a
            %difference at all. Two independent ways it could go vacuous:
            %
            %  (a) too few draws, so both probabilities clamp to the same
            %      order statistic and NOTHING widens; and
            %  (b) a fixture that never reaches the component-SD code.
            %
            %Guard (a) by asserting the raw draws themselves widen between the
            %two probabilities -- that is a property of the fixture alone, with
            %no production code involved. Guard (b) by asserting the
            %coefficient interval (which was never broken) widens too: if the
            %CI argument were being dropped on the floor entirely, this fails
            %and points at the harness rather than at the fix.
            pd = localFixture();
            sig_u = pd.rel.out.sig_u;

            rawNarrow = quantile(sig_u, .975) - quantile(sig_u, .025);
            rawWide   = quantile(sig_u, .995) - quantile(sig_u, .005);
            testCase.assertGreaterThan(rawWide, rawNarrow, ...
                ['The fixture cannot distinguish 95% from 99%: its draws '...
                'clamp to the same order statistic at both probabilities, so '...
                'the widening tests above would pass vacuously. Add draws.']);

            narrow = localSummarize(pd, 0.95);
            wide   = localSummarize(pd, 0.99);
            testCase.assertGreaterThan( ...
                localWidth(wide.relsummary.group(1).event(1).dep), ...
                localWidth(narrow.relsummary.group(1).event(1).dep), ...
                ['The reliability coefficient interval did not widen at 99%. '...
                'That path was never part of this defect, so a failure here '...
                'means the CI argument is not reaching psyrat_relsummary at '...
                'all and the other tests prove nothing.']);
        end

        function testVarianceTableTitleStatesTheRequestedCoverage(testCase)
            %The user-facing half. psyrat_variancet's title used to read "95%"
            %unconditionally; the number shown must be the coverage actually
            %computed. Built with gui,1 because the title is only assembled on
            %the render path -- the headless (gui,0) call returns the vartable
            %without ever forming the string.
            pd = localFixture();
            pd.rel.filename = 'ci_coverage_fixture.psyrat';

            testCase.verifyEqual(localTableTitleAtCi(testCase, pd, 0.95), ...
                'Point and 95% Interval Estimates for the Between- and Within-Person Standard Deviations, SEMs, and ICCs');

            %the differential: same data, different requested width
            testCase.verifyEqual(localTableTitleAtCi(testCase, pd, 0.99), ...
                'Point and 99% Interval Estimates for the Between- and Within-Person Standard Deviations, SEMs, and ICCs');
        end

    end
end

function pd = localFixture()
%A one-facet fixture with enough draws that .5%/99.5% and 2.5%/97.5% are
%genuinely different quantiles. Deterministic on purpose: linspace rather than
%rand, so there is no seed to carry and the expected values are reproducible by
%inspection.
pd = PsyRATTestDataFactory.makeRawSingRelData();

ndraw = 400;
pd.rel.out.mu    = linspace(0.05, 0.15, ndraw)';
pd.rel.out.sig_u = linspace(0.90, 1.30, ndraw)';
pd.rel.out.sig_e = linspace(0.80, 1.20, ndraw)';
end

function out = localSummarize(pd, ciperc)
%The documented scripted entry point, varying only the requested CI.
[out, ~] = psyrat_relsummary('psyrat_data', pd, ...
    'analysis', 'sing', ...
    'depcutoff', 0.70, ...
    'meascutoff', 2, ...
    'depcentmeas', 1, ...
    'CI', ciperc);
end

function w = localWidth(s)
w = s.ul - s.ll;
end

function name = localTableTitleAtCi(testCase, pd, ciperc)
%Render the standard-deviation table and read the title off the figure it
%builds, then close it. Reading the live figure rather than scanning the source
%keeps this a behavioral assertion.
out = localSummarize(pd, ciperc);

fig = [];
cleanup = onCleanup(@() localCloseIfValid(fig));
psyrat_variancet('psyrat_data', out, 'gui', 1);

fig = findobj('Type', 'figure', 'Tag', 'psyrat_output');
testCase.assertNotEmpty(fig, ...
    'psyrat_variancet(gui,1) did not build a figure, so no title can be read.');
name = get(fig(1), 'Name');
end

function localCloseIfValid(fig)
if ~isempty(fig)
    delete(fig(ishandle(fig)));
end
end
