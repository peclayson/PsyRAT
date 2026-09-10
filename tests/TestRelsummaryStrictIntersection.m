classdef TestRelsummaryStrictIntersection < PsyRATTestBase
    %G39 ruling (2026-08-29): group-level inclusion is STRICT -- a participant
    %is retained only by meeting the cutoff in EVERY event, and an event that
    %retained nobody EMPTIES the group's retained set instead of dropping out
    %of the intersection. Before the ruling, a dead event was silently skipped,
    %which made the retained sample depend on which events happened to die:
    %tightening the threshold could ENLARGE the group (killing an event removed
    %a constraint), a dead event's row printed n Included 0 while its siblings
    %printed the group figure, and n Excluded complemented neither.
    %
    %These are the FIRST tests of the mixed state (one event dead, another
    %alive at the same cutoff) -- no committed test exercised it in any family,
    %which is how the skip survived. The fixture starts from the factory's
    %case-4 (two groups x two events) plain-IC struct and separates the two
    %events: 'cor' clears the cutoffs used here within a trial or two, 'err' is
    %hopeless at the observed counts, with s1 alone carrying enough err trials
    %to survive a lax cutoff.

    methods (Test)

        function testDeadEventZeroesTheGroupAndCountsComplement(testCase)
            pd = localMixedFixture();
            [out, relerr] = localSummary(testCase, pd, 0.60);

            %the run must NOT abort: per-event rows stay informative
            testCase.verifyEqual(relerr.nogooddata, 0, ...
                'a dead event must not abort the whole summary');

            for gloc = 1:2
                g = out.relsummary.group(gloc);
                testCase.verifyTrue(iscell(g.goodids) && isempty(g.goodids), ...
                    sprintf(['group %d: a dead event must empty the ' ...
                    'group''s retained set'], gloc));
                testCase.verifyEqual(sort(g.badids), sort({'s1';'s2';'s3'}), ...
                    sprintf(['group %d: every participant must land in ' ...
                    'badids, so included + excluded complement'], gloc));

                for eloc = 1:2
                    testCase.verifyEqual(g.event(eloc).goodn, 0, sprintf( ...
                        ['group %d event %d: goodn is the GROUP-level ' ...
                        'retained count and must read 0 on every row'], ...
                        gloc, eloc));
                end

                %the dead event's reliability-at-cutoff carries the -1
                %sentinel (its cutoff extrapolates past the observed counts)
                errloc = find(strcmp({g.event.name}, 'err'), 1);
                corloc = find(strcmp({g.event.name}, 'cor'), 1);
                testCase.assertNotEmpty(errloc);
                testCase.assertNotEmpty(corloc);
                testCase.verifyEqual(g.event(errloc).rel.m, -1, ...
                    'the dead event''s cutoff row must carry the sentinel');
                testCase.verifyGreaterThan(g.event(corloc).rel.m, 0, ...
                    'the surviving event''s cutoff row must stay informative');
            end
            testCase.verifyEqual(relerr.trlmax, 1, ...
                'the extrapolated err cutoff must set the trlmax flag');
        end

        function testRetentionIsMonotoneInTheCutoff(testCase)
            %The pre-G39 artifact, staged directly: at a lax cutoff the err
            %event retains only s1, so the intersection is {s1}; at a strict
            %cutoff err retains nobody. Under the skip, the strict run's group
            %LOOSENED to cor's full roster (3 > 1, retention increasing in the
            %cutoff); under the strict fold it must empty.
            pd = localMixedFixture();

            outLax = localSummary(testCase, pd, 0.15);
            nLax = numel(outLax.relsummary.group(1).goodids);
            testCase.verifyEqual(nLax, 1, ...
                ['anti-vacuity: the lax cutoff must genuinely constrain ' ...
                'through the err event (only s1 carries enough err trials)']);
            testCase.verifyEqual(outLax.relsummary.group(1).goodids, {'s1'});

            outStrict = localSummary(testCase, pd, 0.60);
            nStrict = numel(outStrict.relsummary.group(1).goodids);
            testCase.verifyEqual(nStrict, 0, ...
                'killing the err event must empty the group, not loosen it');
            testCase.verifyLessThanOrEqual(nStrict, nLax, ...
                'retention must be monotone nonincreasing in the cutoff');
        end

    end
end

%--------------------------------------------------------------------------
function [out, relerr] = localSummary(testCase, pd, depcutoff) %#ok<INUSD>
[out, relerr] = psyrat_relsummary( ...
    'psyrat_data', pd, ...
    'analysis', 'sing', ...
    'depcutoff', depcutoff, ...
    'meascutoff', 2, ...
    'depcentmeas', 1, ...
    'gcoeff', 1, ...
    'CI', 0.95);
end

function pd = localMixedFixture()
%Factory case-4 fixture, reshaped so the two events separate at a cutoff:
%'cor' (sig_e ~ 0.3 against sig_u ~ 0.45) clears the thresholds used here
%within a trial or two (at .60 the GroupA cor cutoff is 2 trials, not 1 -- the
%factory's sig_trl ~ 0.22 counts toward absolute error); 'err' (sig_e ~ 5)
%needs ~168 trials at .60 and ~20 at .15, against observed counts of 6
%(s2, s3) and 30 (s1). Headless (interactive false) so the strict-fold console
%notice never blocks.
pd = PsyRATTestDataFactory.makeICRelDataMultiGroupEvent();
for col = 1:numel(pd.rel.out.labels)
    if startsWith(pd.rel.out.labels{col}, 'err_;_')
        pd.rel.out.sig_e(:, col) = 5.0 + 0.01 * (1:4)';
    else
        pd.rel.out.sig_e(:, col) = 0.30 + 0.01 * (1:4)';
    end
end

%s1 gets 24 extra err trials in each group (30 total vs 6 for s2/s3)
extraIds = {}; extraGroup = {}; extraEvent = {}; extraMeas = [];
for g = {'GroupA', 'GroupB'}
    for t = 1:24
        extraIds{end+1, 1} = 's1'; %#ok<AGROW>
        extraGroup{end+1, 1} = g{1}; %#ok<AGROW>
        extraEvent{end+1, 1} = 'err'; %#ok<AGROW>
        extraMeas(end+1, 1) = 0.05 * t; %#ok<AGROW>
    end
end
extra = table(extraIds, extraGroup, extraEvent, extraMeas, ...
    'VariableNames', {'id', 'group', 'event', 'meas'});
pd.rel.data{2} = [pd.rel.data{2}; extra];

pd.proc = struct('interactive', false);
end
