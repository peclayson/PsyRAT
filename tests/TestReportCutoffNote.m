classdef TestReportCutoffNote < PsyRATTestBase
    %G1 regression pins (owner ruling 2026-09-05, pulled forward from the
    %2026-08-24 council's post-release wave 1). psyrat_relsummary records two
    %D-study failures per cell without an error: no trial count in the
    %projected range reached the threshold (trlcutoff = -1 and the coefficient
    %slots = -1, "Cutoff not calculable"), or the cutoff it found exceeds every
    %participant's observed count (the cutoff is kept and the coefficient slots
    %are overwritten with -1, "Extrapolation beyond data"). The GUI raises a
    %dialog for each from relerr; the headless report had neither the dialog
    %nor relerr, so an exported .cutoff table carried bare -1 cells with no
    %explanation. psyrat_report now returns report.cutoff_note and writes it
    %into every exported header.
    %
    %Fixtures are copies of TestUnreachableCutoffSentinel's G54 fixtures (three
    %participants, five trials each, sig_u about 2 against sig_e about 5), so
    %the sentinel arithmetic is the one that test already pins: at 0.9999 no
    %count up to max + 1000 reaches the threshold; at 0.30 the point estimate
    %crosses by the third trial. The extrapolation state cannot be reached on
    %the no-group, no-event route: sing case 1 has no bad-ids fallback, so
    %when nobody meets a found cutoff its trial-count table is empty and the
    %trlcutoff > trlinfo.max test is vacuously false (refutation wave,
    %2026-09-05). Cases 2-4 fall back to the excluded participants' counts, so
    %the extrapolation test uses TestRelsummaryStrictIntersection's mixed
    %two-group x two-event fixture, whose 'err' event needs well over a hundred
    %trials at .60 (the exact cutoff differs by group) against 6 to 30 observed
    %while 'cor' clears the cutoff. Every test
    %asserts the fixture is in the state it claims before reading the note, so
    %a fixture drift fails loudly rather than letting an assertion pass for
    %the wrong reason.
    %
    %RED against the pre-G1 code: the note field does not exist, so the two
    %sentinel tests and the header test fail; the clean-run test passes on both
    %builds (positive control: the note must be able to be empty).

    methods (Test)

        function testSingNotCalculableProducesTheNote(testCase)
            pd = localSingRelData();
            [out, relerr] = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'sing', 'depcutoff', 0.9999, 'meascutoff', 2, ...
                'depcentmeas', 1, 'CI', 0.95);
            ev = out.relsummary.group(1).event(1);
            testCase.assertEqual(relerr.trlcutoff, 1, 'fixture must be in the not-calculable state');
            testCase.assertEqual(ev.trlcutoff, -1, 'fixture must carry the -1 sentinel');

            report = psyrat_report(out);
            testCase.assertTrue(isfield(report, 'cutoff_note'), ...
                'G1: psyrat_report must return cutoff_note for the one-facet family');
            testCase.verifyTrue(contains(report.cutoff_note, 'could not be calculated'), ...
                'the note must say the cutoff could not be calculated');
            testCase.verifyTrue(contains(report.cutoff_note, 'the measurement'), ...
                'a run with no groups or events names the lone cell "the measurement"');
            testCase.verifyTrue(contains(report.cutoff_note, 'print -1'), ...
                'the note must tell the reader what the -1 cells mean');
            testCase.verifyTrue(contains(report.cutoff_note, 'dependability'), ...
                'the note names the coefficient the summary was built for');
        end

        function testSingExtrapolationProducesTheNote(testCase)
            pd = localMixedGroupEventFixture();
            [out, relerr] = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'sing', 'depcutoff', 0.60, 'meascutoff', 2, ...
                'depcentmeas', 1, 'gcoeff', 1, 'CI', 0.95);
            %fixture guards: every cutoff was found (no -1 trlcutoff), the two
            %'err' cells extrapolate (cutoff beyond the observed counts, the
            %coefficient slots overwritten with -1) and the two 'cor' cells
            %stay clean, so the note must name exactly the former
            testCase.assertEqual(relerr.trlcutoff, 0, 'every cutoff must have been found');
            testCase.assertEqual(relerr.trlmax, 1, 'fixture must be in the extrapolated state');
            errcut = zeros(1, 2);
            for gloc = 1:2
                g = out.relsummary.group(gloc);
                errloc = find(strcmp({g.event.name}, 'err'), 1);
                corloc = find(strcmp({g.event.name}, 'cor'), 1);
                testCase.assertNotEmpty(errloc);
                testCase.assertNotEmpty(corloc);
                testCase.assertGreaterThan(g.event(errloc).trlcutoff, g.event(errloc).trlinfo.max, ...
                    'the err cutoff must exceed every observed count');
                testCase.assertEqual(g.event(errloc).rel.m, -1, ...
                    'the err coefficient slot must carry the -1 sentinel');
                testCase.assertGreaterThan(g.event(corloc).rel.m, 0, ...
                    'the cor cell must stay clean');
                errcut(gloc) = g.event(errloc).trlcutoff;
            end

            report = psyrat_report(out);
            testCase.assertTrue(isfield(report, 'cutoff_note'));
            testCase.verifyTrue(contains(report.cutoff_note, 'extrapolation beyond the data'), ...
                'the note must say the cutoff is an extrapolation');
            testCase.verifyTrue(contains(report.cutoff_note, sprintf('GroupA - err (cutoff %d)', errcut(1))), ...
                'the note must name the GroupA err cell with its extrapolated cutoff');
            testCase.verifyTrue(contains(report.cutoff_note, sprintf('GroupB - err (cutoff %d)', errcut(2))), ...
                'the note must name the GroupB err cell with its extrapolated cutoff');
            testCase.verifyFalse(contains(report.cutoff_note, '- cor'), ...
                'the clean cor cells must not be named');
            testCase.verifyFalse(contains(report.cutoff_note, 'could not be calculated'), ...
                'an extrapolated cutoff is not a not-calculable one');
        end

        function testCleanRunHasAnEmptyNote(testCase)
            %Positive control: the same fixture at a reachable threshold gives
            %no note, so the two tests above are driven by the sentinel state.
            pd = localSingRelData();
            [out, relerr] = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'sing', 'depcutoff', 0.30, 'meascutoff', 2, ...
                'depcentmeas', 1, 'CI', 0.95);
            ev = out.relsummary.group(1).event(1);
            testCase.assertEqual([relerr.trlcutoff relerr.trlmax], [0 0], 'fixture must be clean');
            testCase.assertTrue(ev.trlcutoff > 0 && ev.trlcutoff <= 5 && ev.rel.m ~= -1, ...
                'the clean fixture must carry a real cutoff and coefficient');

            report = psyrat_report(out);
            testCase.verifyTrue(~isfield(report, 'cutoff_note') || isempty(report.cutoff_note), ...
                'a clean run must not carry a cutoff note');
        end

        function testNoteReachesEveryExportedHeader(testCase)
            pd = localSingRelData();
            [out, ~] = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'sing', 'depcutoff', 0.9999, 'meascutoff', 2, ...
                'depcentmeas', 1, 'CI', 0.95);
            testCase.assertEqual(out.relsummary.group(1).event(1).trlcutoff, -1);
            tdir = tempname;
            mkdir(tdir);
            testCase.addTeardown(@() rmdir(tdir, 's'));

            report = psyrat_report(out, 'outdir', tdir, 'format', '.csv');
            testCase.assertNotEmpty(report.files, 'the report must have written files');
            for k = 1:numel(report.files)
                txt = fileread(report.files{k});
                testCase.verifyTrue(contains(txt, 'could not be calculated'), ...
                    sprintf('the sentinel note must be in the header of %s', report.files{k}));
            end
            %and the -1 cells themselves are unchanged: the note explains, it
            %does not rewrite the table the GUI export also produces
            cutoffFile = report.files{find(contains(report.files, 'cutoff'), 1)};
            testCase.verifyTrue(contains(fileread(cutoffFile), '-1'), ...
                'the cutoff table must still print the -1 the GUI export prints');
        end

        function testTrtNotCalculableProducesTheNote(testCase)
            %The test-retest families keep the coefficient under relcutoff, not
            %rel; the note reader must find the sentinel there too.
            pd = localTrtRelData();
            [out, relerr] = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'trt', 'gcoeff', 1, 'reltype', 3, ...
                'relcutoff', 0.9999, 'meascutoff', 2, 'relcentmeas', 1, ...
                'CI', 0.95);
            ev = out.relsummary.group(1).event(1);
            testCase.assertEqual(relerr.trlcutoff, 1, 'fixture must be in the not-calculable state');
            testCase.assertEqual(ev.trlcutoff, -1);
            testCase.assertEqual(ev.relcutoff.m, -1, 'the trt sentinel lives under relcutoff');

            report = psyrat_report(out);
            testCase.assertTrue(isfield(report, 'cutoff_note'), ...
                'G1: psyrat_report must return cutoff_note for the test-retest family');
            testCase.verifyTrue(contains(report.cutoff_note, 'could not be calculated'));
        end

    end
end

function psyrat_data = localSingRelData()
%Copy of TestUnreachableCutoffSentinel's localSingRelData('table') (G54): a
%minimal pre-relsummary one-facet struct routed to sing analysis 1. Three
%participants x five trials; sig_u about 2 against sig_e about 5, so the
%posterior-mean dependability is about .44 at the five observed trials and
%about .56 at eight projected trials. Kept as a copy so this file stays
%self-contained; if the G54 fixture changes, change both.
subs = {'s1','s2','s3'};
ntrl = 5;
ids = {}; meas = [];
for s = 1:numel(subs)
    for t = 1:ntrl
        ids{end+1,1} = subs{s}; %#ok<AGROW>
        meas(end+1,1) = numel(meas) + 1; %#ok<AGROW>
    end
end

rel = struct();
rel.analysis = 'ic';
rel.groups = 'none';
rel.events = 'none';
rel.filename = 'g1_sing_cutoff_note_fixture.psyrat';
rel.nchains = 2;
rel.niter = 100;
rel.data = table(ids, meas, 'VariableNames', {'id','meas'});

rel.out = struct();
rel.out.labels = {'measure'};
rel.out.mu    = [0.10; 0.12; 0.11; 0.09; 0.10];
rel.out.sig_u = [2.00; 2.10; 1.90; 2.05; 1.95];
rel.out.sig_e = [4.90; 5.00; 5.10; 4.80; 5.20];
rel.out.sig_trl = [0.20; 0.25; 0.23; 0.21; 0.22];

psyrat_data = struct('rel', rel);
%headless: never block on a modal
psyrat_data.proc = struct('interactive', false, 'measheader', 'meas');
end

function psyrat_data = localTrtRelData()
%Copy of TestUnreachableCutoffSentinel's localTrtRelData('table') (G54): a
%minimal pre-relsummary test-retest struct routed to trt analysis 1. Three
%participants x two occasions x five trials; at reltype 3 with the
%single-occasion estimand the occasion terms put a hard asymptote near .94,
%so 0.9999 is unreachable at any trial count.
subs = {'s1','s2','s3'};
occ = {'t1','t2'};
ntrl = 5;
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
rel.filename = 'g1_trt_cutoff_note_fixture.psyrat';
rel.nchains = 2;
rel.niter = 100;
rel.time = occ(:);
rel.data = table(ids, meas, times, 'VariableNames', {'id','meas','time'});

rel.out = struct();
rel.out.labels = {'measure'};
rel.out.mu = [0.10; 0.12; 0.11; 0.09; 0.10];
rel.out.sig_id      = [2.00; 2.10; 1.90; 2.05; 1.95];
rel.out.sig_occ     = [0.30; 0.35; 0.28; 0.32; 0.31];
rel.out.sig_trl     = [0.20; 0.25; 0.23; 0.21; 0.22];
rel.out.sig_trlxid  = [0.20; 0.22; 0.18; 0.21; 0.19];
rel.out.sig_occxid  = [0.40; 0.42; 0.38; 0.41; 0.39];
rel.out.sig_trlxocc = [0.09; 0.10; 0.11; 0.09; 0.10];
rel.out.sig_err     = [4.90; 5.00; 5.10; 4.80; 5.20];

psyrat_data = struct('rel', rel);
psyrat_data.proc = struct('interactive', false, 'measheader', 'meas');
end

function pd = localMixedGroupEventFixture()
%Copy of TestRelsummaryStrictIntersection's localMixedFixture: the factory's
%two-group x two-event plain-IC struct (relsummary 'sing' case 4) reshaped so
%the two events separate at a cutoff. 'cor' (sig_e about 0.3 against sig_u
%about 0.45) clears .60 within a trial or two; 'err' (sig_e about 5) needs
%well over a hundred trials (the exact cutoff differs by group), against
%observed counts of 6 (s2, s3) and 30 (s1), so its
%cutoff is found but extrapolates past every participant and the case-4
%bad-ids fallback marks the extrapolation. Kept as a copy so this file stays
%self-contained; if that fixture changes, change both.
pd = PsyRATTestDataFactory.makeICRelDataMultiGroupEvent();
for col = 1:numel(pd.rel.out.labels)
    if startsWith(pd.rel.out.labels{col}, 'err_;_')
        pd.rel.out.sig_e(:, col) = 5.0 + 0.01 * (1:4)';
    else
        pd.rel.out.sig_e(:, col) = 0.30 + 0.01 * (1:4)';
    end
end
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
%headless: never block on a modal
pd.proc = struct('interactive', false, 'measheader', 'meas');
end
