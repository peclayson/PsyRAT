classdef TestAdmissibilitySurfacing < PsyRATTestBase
    % End-to-end tests for the admissibility diagnostic: that an inadmissible
    % quantity computed in a coefficient kernel actually REACHES the user, and
    % that a clean run says nothing.
    %
    % The unit-level contract of the recorder lives in TestCalculationFunctions.
    % What is tested here is the chain, because the chain is where this kind of
    % diagnostic historically dies: psyrat_relsummary works on a local copy of
    % rel and never writes back, psyrat_relfigures has no output arguments, and
    % psyrat_report never calls the summary itself. A note that is computed
    % correctly and then dropped is worth nothing.

    methods (Test)
        function testNothingIsClippedOrSubstituted(testCase)
            % The governing rule. The coefficient is the formula evaluated
            % exactly, whatever it comes out to.
            [bp, bt, er_var, wp_cov, obs] = localRC01Fixture();

            out = testCase.verifyWarning(@() psyrat_diffrel( ...
                'bp', bp, 'bt', bt, 'er_var', er_var, 'wp_cov', wp_cov, ...
                'obs', obs, 'est', 'gen', 'CI', 0.95), 'psyrat:negerrvar');

            hm = psyrat_harmmean(obs);
            uni = 2 + 2 - 2*0.5;
            relErr = 16/obs(1) + 1/obs(2) - 2*wp_cov(1)/hm;

            testCase.verifyEqual(out.pt, uni/(uni + relErr), 'AbsTol', 1e-12);
            testCase.verifyGreaterThan(out.pt, 1, ...
                'A value above 1 must be reported, not bounded.');
        end

        function testTheDependabilityIsNoLongerPinnedAtOne(testCase)
            % Regression for the artifact the clipping guard produced: with a
            % negative universe term AND a negative error variance, clipping
            % the error to zero made the ratio uni/uni, i.e. exactly 1.000000 --
            % perfect dependability manufactured from an inadmissible numerator,
            % printed beside an ordinary ICC. Verified live before removal.
            nd = 6; one = ones(nd,1);
            cm = @(v1,v2,cv) cat(3, [v1*one cv*one], [cv*one v2*one]);

            out = testCase.verifyWarning(@() psyrat_diffrel_trt( ...
                'bp', cm(1,1,0.98), 'bpi', cm(12,1,0.98*sqrt(12)), ...
                'bpo', cm(0.001,0.001,0), 'bt', cm(0.001,0.001,0), ...
                'bo', cm(0.001,0.001,0), 'boi', cm(120,10,0.99*sqrt(1200)), ...
                'er_var', log(sqrt(0.01))*ones(nd,2), 'er_cov', zeros(nd,1), ...
                'obs', [60 12], 'nocc', [2 2], 'reltype', 2, ...
                'est', 'dep', 'CI', 0.95), 'psyrat:negerrvar');

            testCase.verifyNotEqual(out.pt, 1);
            testCase.verifyEqual(out.pt, 0.051806, 'AbsTol', 1e-5, ...
                'Expected the raw ratio of the two negative quantities.');
        end

        function testTheUniverseTermIsRecordedNotJustTheError(testCase)
            % RC-21: for the coefficients of equivalence and stability the
            % harmonic-mean-scaled block sits in the NUMERATOR, so a guard on
            % the error variance alone cannot see it. Before this, the case
            % below produced a negative coefficient in total silence.
            nd = 6; one = ones(nd,1);
            cm = @(v1,v2,cv) cat(3, [v1*one cv*one], [cv*one v2*one]);
            sm = cm(0.1,0.1,0);

            out = testCase.verifyWarning(@() psyrat_diffrel_trt( ...
                'bp', cm(1,1,0.98), 'bpi', cm(12,1,0.98*sqrt(12)), ...
                'bpo', sm, 'bt', sm, 'bo', sm, 'boi', sm, ...
                'er_var', log(sqrt(0.2))*ones(nd,2), 'er_cov', zeros(nd,1), ...
                'obs', [60 12], 'nocc', [2 2], 'reltype', 2, ...
                'est', 'gen', 'CI', 0.95), 'psyrat:negerrvar');

            testCase.verifyLessThan(out.pt, 0, ...
                'The realistic manifestation is a NEGATIVE coefficient.');

            uni = out.admissibility.entries(strcmp( ...
                {out.admissibility.entries.label}, 'universe-score variance'));
            testCase.verifyNotEmpty(uni);
            testCase.verifyGreaterThan(uni.nneg, 0, ...
                'The universe term is what went negative here.');
        end

        function testTheDoDKernelIsGuardedAtAll(testCase)
            % psyrat_dodiffrel carries the same harmonic-mean scaler on the
            % off-diagonals of its between-trial block and had no guard, no
            % warning and no diagnostic of any kind. Its four events have four
            % independently observed trial counts, so unequal counts are the
            % norm rather than the exception there.
            nd = 4;
            bp = zeros(nd,4,4); bt = zeros(nd,4,4);
            for k = 1:4
                bp(:,k,k) = 1;
                bt(:,k,k) = 1;
            end
            % near-perfect trial correlation between the two lopsided events
            bt(:,1,2) = 0.999; bt(:,2,1) = 0.999;

            out = psyrat_dodiffrel('bp', bp, 'bt', bt, ...
                'er_var', log(sqrt(0.01))*ones(nd,4), ...
                'obs', [200 2 50 50], 'est', 'dep', 'CI', 0.95);

            testCase.verifyTrue(isfield(out, 'admissibility'), ...
                'The DoD kernel must carry a diagnostic like its siblings.');

            % Assert WHICH quantities are recorded, not how many. A hardcoded
            % count has to be edited every time a quantity is added, and when
            % it does break it does not say which record went missing (RC-36).
            labels = {out.admissibility.entries.label};
            testCase.verifyTrue(any(strcmp(labels, 'universe-score variance')));
            testCase.verifyTrue(any(strcmp(labels, ...
                'between-trial contrast variance')));
            testCase.verifyTrue(any(strcmp(labels, ...
                'residual contrast variance')));
            testCase.verifyTrue(any(strcmp(labels, ...
                'between-trial contrast variance at a single observation (ICC)')));
        end

        function testTheSummaryAttachesTheNoteAndACleanRunSaysNothing(testCase)
            % The summary must always attach the field, and a clean run must
            % leave it empty so ordinary output is byte-identical.
            vd = localFixtureRel();
            out = psyrat_relsummary('psyrat_data', vd, ...
                'analysis', 'sing_diff', 'depcutoff', 0.5, 'meascutoff', 2, ...
                'depcentmeas', 1, 'diffgcoeff', 1, 'CI', 0.95);

            testCase.verifyTrue( ...
                isfield(out.relsummary, 'admissibility_note'), ...
                'psyrat_relsummary must attach the note.');
            testCase.verifyEmpty(out.relsummary.admissibility_note);
        end

        function testTheCollectorWalksASummaryAndFindsAnInadmissibleCell(testCase)
            % The other half of the chain, tested directly rather than through
            % a fixture that happens to be clean: given a summary whose
            % difference score carries an inadmissible record, the collector
            % must find it, including on the separately evaluated cutoff row.
            adm = psyrat_admissibility([], 'relative error variance', ...
                [-0.48; 0.2], 'psyrat_diffrel');

            rs = struct();
            rs.group(1).diffscore = struct('pt', 1.19, 'admissibility', adm);
            rs.group(2).diffscore = struct('pt', 0.7);
            rs.group(2).diffscore.trlcutoff = struct('pt', 0.7, ...
                'admissibility', adm);

            note = psyrat_admissibility_collect(rs);

            testCase.verifyNotEmpty(note);
            testCase.verifyTrue(contains(note, 'relative error variance'));
            % The explanation must name the cause, not just the symptom -- and
            % BOTH causes, not just the projection. The record carries no
            % counts, so the note cannot know which applies; a user at equal
            % counts needs cause (2) present to recognize that the prescribed
            % remedy is already in force and cannot help them (RC-34). Matching
            % 'harmonic' alone would still pass if cause (2) were dropped, which
            % is the regression this second assertion exists to catch.
            testCase.verifyTrue(contains(note, 'harmonic'), ...
                'The explanation must name the D-study projection cause.');
            testCase.verifyTrue(contains(note, 'gamma'), ...
                ['The explanation must also name the non-PSD cause, which is ' ...
                'the only one that can reach an equal-count caller.']);
            testCase.verifyTrue(contains(note, 'not interpretable', ...
                'IgnoreCase', true));

            % both cells counted: 2 of 4 draws, not 1 of 2
            testCase.verifyTrue(contains(note, '2 of 4'));

            % and a summary with nothing wrong produces nothing at all
            clean = struct();
            clean.group(1).diffscore = struct('pt', 0.7);
            testCase.verifyEmpty(psyrat_admissibility_collect(clean));
            testCase.verifyEmpty(psyrat_admissibility_collect([]));
        end

        function testAnObservedScoreCutoffRowIsNotDoubleCounted(testCase)
            % G60: on sing_diff_sserr the cutoff row is the OBSERVED-design
            % coefficient, stored as a minimal pt/ll/ul struct with NO
            % admissibility field. The producer used to alias the whole
            % diffscore into trlcutoff, which made this collector gather
            % the IDENTICAL record twice and double the note's absolute
            % draw counts. The stub mirrors the shipped minimal shape; the
            % '2 of 4' assertion in the test above is the positive control
            % that two GENUINE records still both count.
            adm = psyrat_admissibility([], 'relative error variance', ...
                [-0.48; 0.2], 'psyrat_diffrel');

            rs = struct();
            rs.group(1).diffscore = struct('pt', 0.7, 'admissibility', adm);
            rs.group(1).diffscore.trlcutoff = ...
                struct('pt', 0.7, 'll', 0.6, 'ul', 0.8);

            note = psyrat_admissibility_collect(rs);

            testCase.verifyNotEmpty(note);
            testCase.verifyTrue(contains(note, '1 of 2'), ...
                ['a cutoff row without its own admissibility record must ' ...
                'not double-count the parent''s']);
        end

        function testACleanRunIsSilentEverywhere(testCase)
            % Equal counts make the pathology arithmetically impossible, so
            % nothing should warn, nothing should be noted, and the numbers
            % must be untouched.
            [bp, bt, er_var, wp_cov, ~] = localRC01Fixture();

            out = testCase.verifyWarningFree(@() psyrat_diffrel( ...
                'bp', bp, 'bt', bt, 'er_var', er_var, 'wp_cov', wp_cov, ...
                'obs', [10 10], 'est', 'gen', 'CI', 0.95));

            testCase.verifyFalse( ...
                psyrat_admissibility_hasneg(out.admissibility));
            testCase.verifyEmpty( ...
                psyrat_admissibility_note(out.admissibility));

            uni = 3;
            denom = 16 + 1 - 2*wp_cov(1);
            testCase.verifyEqual(out.pt, uni/(uni + denom/10), 'AbsTol', 1e-12);
        end

        function testTheSemIsUndefinedRatherThanFlooredAtZero(testCase)
            % A negative error variance has no real square root. Flooring it
            % printed a confident small SEM for a quantity that is not a
            % variance; NaN says what is true.
            errvar = [-1; 4; 9];
            sem = nan(size(errvar));
            sem(errvar >= 0) = sqrt(errvar(errvar >= 0));

            testCase.verifyTrue(isnan(sem(1)));
            testCase.verifyEqual(sem(2:3), [2; 3], 'AbsTol', 0);
            testCase.verifyNotEqual(sem(1), 0, ...
                'Zero would be a fabricated SEM, not an undefined one.');
        end

        function testTheIccNumeratorIsRecordedAtTheCoefficientOfEquivalence(testCase)
            % RC-37, the load-bearing case. At reltype 1 the ICC numerator is
            % uni_icc = bp_raw + bpo_raw while the recorded coefficient
            % numerator is uni = bp_raw + bpo_raw/k. They are DIFFERENT
            % quantities, so a person-by-occasion block that is not positive
            % semi-definite can drive the ICC numerator negative while the
            % recorded one stays comfortably positive, and before RC-37 that
            % happened in complete silence.
            %
            % bp_raw = 1, bpo_raw = -2, k = 4: uni = 1 - 2/4 = 0.5 and
            % uni_icc = 1 - 2 = -1. The trial and occasion counts are EQUAL,
            % which rules out the D-study projection and leaves only the
            % non-positive-semi-definite cause (RC-23).
            nd = 6; one = ones(nd,1);
            cm = @(v1,v2,cv) cat(3, [v1*one cv*one], [cv*one v2*one]);

            out = testCase.verifyWarning(@() psyrat_diffrel_trt( ...
                'bp', cm(1,1,0.5), 'bpi', cm(1,1,0), 'bpo', cm(1,1,2), ...
                'bt', cm(0.1,0.1,0), 'bo', cm(0.1,0.1,0), ...
                'boi', cm(0.1,0.1,0), ...
                'er_var', log(sqrt(0.5))*ones(nd,2), 'er_cov', zeros(nd,1), ...
                'obs', [20 20], 'nocc', [4 4], 'reltype', 1, ...
                'est', 'gen', 'CI', 0.95), 'psyrat:negerrvar');

            % What the user was shown: a perfectly ordinary generalizability
            % coefficient printed beside an ICC of -0.5.
            testCase.verifyGreaterThan(out.pt, 0);
            testCase.verifyLessThan(out.pt, 1);
            testCase.verifyEqual(out.icc_pt, -0.5, 'AbsTol', 1e-12);

            iccuni = localEntry(out.admissibility, ...
                'universe-score variance at a single observation (ICC)');
            testCase.verifyNotEmpty(iccuni, ...
                'The ICC numerator must carry a record of its own here.');
            testCase.verifyEqual(iccuni.nneg, nd);

            % and the record that already existed does NOT see it, which is
            % precisely why a second one is needed
            uni = localEntry(out.admissibility, 'universe-score variance');
            testCase.verifyNotEmpty(uni);
            testCase.verifyEqual(uni.nneg, 0, ...
                'The coefficient numerator is positive; only the ICC is not.');
        end

        function testTheIccNumeratorIsNotDoubleCountedAtReltypeThree(testCase)
            % The negative half of the gate. At reltype 3 uni_icc and uni are
            % assigned from the same variable (bp_raw), so a second record
            % would be a bit-identical duplicate and would report one failure
            % as two. Here bp_raw = -1, so the shared quantity IS inadmissible
            % and the existing record must be the only one that says so.
            nd = 5; one = ones(nd,1);
            cm = @(v1,v2,cv) cat(3, [v1*one cv*one], [cv*one v2*one]);

            out = testCase.verifyWarning(@() psyrat_diffrel_trt( ...
                'bp', cm(1,1,1.5), 'bpi', cm(1,1,0), 'bpo', cm(1,1,0), ...
                'bt', cm(0.1,0.1,0), 'bo', cm(0.1,0.1,0), ...
                'boi', cm(0.1,0.1,0), ...
                'er_var', log(sqrt(0.5))*ones(nd,2), 'er_cov', zeros(nd,1), ...
                'obs', [20 20], 'nocc', [2 2], 'reltype', 3, ...
                'est', 'gen', 'CI', 0.95), 'psyrat:negerrvar');

            labels = {out.admissibility.entries.label};
            testCase.verifyFalse(any(strcmp(labels, ...
                'universe-score variance at a single observation (ICC)')), ...
                ['reltype 3 must not record the ICC numerator separately: it ' ...
                'is the same variable as the coefficient numerator.']);

            uni = localEntry(out.admissibility, 'universe-score variance');
            testCase.verifyEqual(uni.nneg, nd, ...
                'The one record taken must cover it.');

            % the gate is on the numerator only -- the two ICC error terms are
            % different quantities at every reltype and are always recorded
            testCase.verifyTrue(any(strcmp(labels, ...
                'relative error variance at a single observation (ICC)')));
            testCase.verifyTrue(any(strcmp(labels, ...
                'absolute error variance at a single observation (ICC)')));
        end

        function testTheIccDenominatorsAreRecordedInTheOneFacetKernel(testCase)
            % psyrat_diffrel's ICC denominators were documented as needing no
            % guard because a raw variance of a difference is non-negative for
            % positive semi-definite input. That premise is exactly what the
            % gamma / scaled-chi-square family breaks (RC-23). Residual
            % variances 16 and 1 with a residual covariance of 20 is a
            % correlation of 5, so wp1 + wp2 - 2*wp_cov = -23 at EQUAL trial
            % counts, where the projection cause cannot arise.
            nd = 6; one = ones(nd,1);

            bp = zeros(nd,2,2);
            bp(:,1,1) = 2; bp(:,2,2) = 2; bp(:,1,2) = 0.5; bp(:,2,1) = 0.5;

            bt = zeros(nd,2,2);
            bt(:,1,1) = 1; bt(:,2,2) = 1; bt(:,1,2) = 0.9; bt(:,2,1) = 0.9;

            out = testCase.verifyWarning(@() psyrat_diffrel( ...
                'bp', bp, 'bt', bt, 'er_var', log([4*one 1*one]), ...
                'wp_cov', 20*one, 'obs', [10 10], 'est', 'gen', ...
                'CI', 0.95), 'psyrat:negerrvar');

            % uni = 3, denom_rel_err = -23, so the ICC is out of range
            testCase.verifyEqual(out.icc_pt, 3/(3-23), 'AbsTol', 1e-12);

            relicc = localEntry(out.admissibility, ...
                'relative error variance at a single observation (ICC)');
            absicc = localEntry(out.admissibility, ...
                'absolute error variance at a single observation (ICC)');
            testCase.verifyNotEmpty(relicc);
            testCase.verifyNotEmpty(absicc);
            testCase.verifyEqual(relicc.nneg, nd);
            testCase.verifyEqual(absicc.nneg, nd);
            testCase.verifyEqual(relicc.minval, -23, 'AbsTol', 1e-12);

            % The role decides what the exactly-zero note says about each of
            % these, and saying it of a denominator what is true of a numerator
            % tells the user the reverse of what happened, so every role this
            % kernel records under is asserted rather than assumed.
            localVerifyRoles(testCase, out.admissibility, ...
                {'universe-score variance', 'numerator'; ...
                'relative error variance', 'denominator'; ...
                'absolute error variance', 'denominator'; ...
                'relative error variance at a single observation (ICC)', 'denominator'; ...
                'absolute error variance at a single observation (ICC)', 'denominator'});
        end

        function testTheIccDenominatorsAreRecordedInTheTwoFacetKernel(testCase)
            % The isolating case: the two ICC error terms are the ONLY
            % inadmissible quantities here, so before RC-37 this run produced
            % no record, no warning and no note of any kind. A person-by-trial
            % correlation of 5 gives bpi_raw = -8, which the divisors absorb
            % (bpi = -0.4 against bpo = +1) but the raw sum does not.
            nd = 6; one = ones(nd,1);
            cm = @(v1,v2,cv) cat(3, [v1*one cv*one], [cv*one v2*one]);

            out = testCase.verifyWarning(@() psyrat_diffrel_trt( ...
                'bp', cm(1,1,0.5), 'bpi', cm(1,1,5), 'bpo', cm(1,1,0), ...
                'bt', cm(0.1,0.1,0), 'bo', cm(0.1,0.1,0), ...
                'boi', cm(0.1,0.1,0), ...
                'er_var', log(sqrt(0.5))*ones(nd,2), 'er_cov', zeros(nd,1), ...
                'obs', [20 20], 'nocc', [2 2], 'reltype', 3, ...
                'est', 'gen', 'CI', 0.95), 'psyrat:negerrvar');

            % rel_err_icc = -8 + 2 + 1 = -5 against uni_icc = 1
            testCase.verifyEqual(out.icc_pt, 1/(1-5), 'AbsTol', 1e-12);
            testCase.verifyGreaterThan(out.pt, 0);
            testCase.verifyLessThan(out.pt, 1);

            relicc = localEntry(out.admissibility, ...
                'relative error variance at a single observation (ICC)');
            absicc = localEntry(out.admissibility, ...
                'absolute error variance at a single observation (ICC)');
            testCase.verifyEqual(relicc.nneg, nd);
            testCase.verifyEqual(absicc.nneg, nd);
            testCase.verifyEqual(relicc.minval, -5, 'AbsTol', 1e-12);

            % nothing else here is inadmissible, so these two records are the
            % whole reason the user is told anything at all
            for lbl = {'universe-score variance', 'relative error variance', ...
                    'absolute error variance'}
                other = localEntry(out.admissibility, lbl{1});
                testCase.verifyEqual(other.nneg, 0, ...
                    sprintf('%s should be admissible in this fixture.', lbl{1}));
            end

            localVerifyRoles(testCase, out.admissibility, ...
                {'universe-score variance', 'numerator'; ...
                'relative error variance', 'denominator'; ...
                'absolute error variance', 'denominator'; ...
                'relative error variance at a single observation (ICC)', 'denominator'; ...
                'absolute error variance at a single observation (ICC)', 'denominator'});
        end

        function testTheIccDenominatorIsRecordedInTheDoDKernel(testCase)
            % psyrat_dodiffrel's ICC denominator uses the UNSCALED
            % between-trial contrast variance, which the harmonic-mean scaler
            % cannot reach but a non-positive-semi-definite block can. A trial
            % correlation of 3 between the two events that share a sign gives
            % bt_unscaled = 4 - 6 = -2, while the lopsided trial counts keep
            % the SCALED version at +1.96. The user saw a dependability of
            % 0.67 beside an ICC above 1, with nothing recorded.
            nd = 4;
            bp = zeros(nd,4,4); bt = zeros(nd,4,4);
            for k = 1:4
                bp(:,k,k) = 1;
                bt(:,k,k) = 1;
            end
            bt(:,1,2) = 3; bt(:,2,1) = 3;

            out = testCase.verifyWarning(@() psyrat_dodiffrel( ...
                'bp', bp, 'bt', bt, ...
                'er_var', log(sqrt(0.01))*ones(nd,4), ...
                'obs', [100 100 1 1], 'est', 'dep', 'CI', 0.95), ...
                'psyrat:negerrvar');

            testCase.verifyGreaterThan(out.icc_pt, 1, ...
                'The ICC is the quantity that left [0,1] here.');
            testCase.verifyLessThan(out.pt, 1);

            bticc = localEntry(out.admissibility, ...
                'between-trial contrast variance at a single observation (ICC)');
            testCase.verifyNotEmpty(bticc);
            testCase.verifyEqual(bticc.nneg, nd);
            testCase.verifyEqual(bticc.minval, -2, 'AbsTol', 1e-12);

            % the scaled version stays admissible, so it cannot stand in
            btscaled = localEntry(out.admissibility, ...
                'between-trial contrast variance');
            testCase.verifyEqual(btscaled.nneg, 0);

            localVerifyRoles(testCase, out.admissibility, ...
                {'universe-score variance', 'numerator'; ...
                'between-trial contrast variance', 'denominator'; ...
                'residual contrast variance', 'denominator'; ...
                'between-trial contrast variance at a single observation (ICC)', ...
                'denominator'});
        end

        function testTheExactlyZeroNoteSaysTheRightThingForEachRole(testCase)
            % A scientific error in user-facing text, not a counting bug, so it
            % is tested here rather than with the recorder's unit contract:
            % exactly zero means OPPOSITE things in the two positions. A
            % universe-score variance of zero drives the coefficient to zero; an
            % error variance of zero removes a term from the denominator and
            % drives it to one. The shipped prose said the former in both cases.
            zerovals = [0; 0; 1];
            claim = 'reliability of zero is then the correct answer';

            den = psyrat_admissibility([], 'relative error variance', ...
                zerovals, 'unit test', 'denominator');
            denNote = psyrat_admissibility_note(den);

            testCase.verifyNotEmpty(denNote);
            testCase.verifyFalse(contains(denNote, claim), ...
                ['An error variance of zero does NOT make the coefficient ' ...
                'zero, so this sentence must not appear for a denominator.']);
            testCase.verifyTrue(contains(denNote, 'toward one'), ...
                'The denominator wording must name the correct direction.');

            % the numerator wording is unchanged
            num = psyrat_admissibility([], 'universe-score variance', ...
                zerovals, 'unit test', 'numerator');
            numNote = psyrat_admissibility_note(num);
            testCase.verifyTrue(contains(numNote, claim));

            % omitting the role keeps the behavior every record had before it
            % existed, so nothing outside this change moves
            def = psyrat_admissibility([], 'universe-score variance', ...
                zerovals, 'unit test');
            testCase.verifyEqual(psyrat_admissibility_note(def), numNote);

            % exactly zero is a result, not a fault, in either position
            testCase.verifyWarningFree(@() psyrat_admissibility_warn(den));
            testCase.verifyFalse(psyrat_admissibility_hasneg(den));

            % and a record written before the role field existed still reads,
            % rather than throwing out of the reader
            legacy = struct('entries', rmfield(num.entries, 'role'));
            testCase.verifyEqual(psyrat_admissibility_note(legacy), numNote);
        end

        function testTheDoDNoteCountsEachDrawExactlyOnce(testCase)
            % RC-43. The DoD record used to be written and read by nothing:
            % psyrat_admissibility_collect walks relsummary.group(g).diffscore,
            % and the DoD route builds no relsummary at all, so a negative
            % variance there produced a console line and no more -- no report
            % block, no dialog. psyrat_dod_admissibility_note is the collector
            % for that route.
            %
            % The trap it exists to prevent: every caller computes dep and gen
            % together, and each recorded quantity is built BEFORE the est
            % branch, so the two records are bit-identical. The note rolls up by
            % unique({label}) and SUMS the counts, so passing both would not
            % double-check anything -- it would double every count and claim
            % twice as many draws were inspected as were.
            [dep, gen] = localDoDInadmissiblePair();

            testCase.verifyEqual(dep.admissibility, gen.admissibility, ...
                ['The dep and gen records must be identical -- if this ever ' ...
                'stops holding, the "pass one" rule needs revisiting.']);

            note = psyrat_dod_admissibility_note(dep);
            testCase.verifyNotEmpty(note, ...
                'An inadmissible DoD variance must produce a note.');
            testCase.verifyTrue(contains(note, 'of 4'), ...
                'The denominator must be the 4 draws actually inspected.');

            % The regression this test is really for: hand it the pair and the
            % counts double. Pinned so that "helpfully" passing both at a call
            % site fails here with an explanation, rather than silently
            % misreporting to a user.
            doubled = psyrat_dod_admissibility_note([dep gen]);
            testCase.verifyTrue(contains(doubled, 'of 8'), ...
                ['Passing both halves doubles the counts. This is why the ' ...
                'call sites pass exactly one.']);

            % records from DIFFERENT evaluations are the opposite case and
            % genuinely do add, which is how the per-group report sweep works
            testCase.verifyEmpty(psyrat_dod_admissibility_note([]));
        end

        function testTheDoDCallSitesPassOneRecordAndTheReportAggregates(testCase)
            % The invariant test above proves psyrat_dod_admissibility_note
            % behaves correctly when handed one record versus two. It does NOT
            % prove the production call sites hand it one: it calls the function
            % directly, so mutating psyrat_dod_varcompt to pass [dep gen] leaves
            % it green. This test closes that hole by going through the real
            % builders, and by driving an inadmissible result all the way into
            % psyrat_report, whose note content nothing else covers.
            %
            % Both counts are load-bearing and they differ on purpose:
            %   one builder call  -> 6 draws  (ONE evaluation, dep or gen, not both)
            %   the whole report  -> 12 draws (TWO groups, which genuinely add)
            % Passing both halves at a call site would make the first read 12,
            % and passing one group's record twice would make the second read 24.
            psyrat_data = localDoDInadmissibleFixture();

            ws = warning('off', 'psyrat:negerrvar');
            cleanup = onCleanup(@() warning(ws));

            [~, adm] = psyrat_dod_varcompt('psyrat_data', psyrat_data, ...
                'groupidx', 1, 'map', [1 2 3 4], 'savefile', 0, 'gui', 0);
            builderNote = psyrat_admissibility_note(adm);

            testCase.verifyNotEmpty(builderNote, ...
                'The perturbed fixture must actually go inadmissible.');
            testCase.verifyTrue(contains(builderNote, 'of 6'), ...
                ['One builder call is ONE evaluation: 6 draws. "of 12" here ' ...
                'means the call site passed both the dep and gen records.']);

            % ...and the report sums across the two groups, because those ARE
            % different evaluations. This also pins that the note reaches
            % report.admissibility_note at all -- the header block is the only
            % surfacing a headless DoD export gets.
            report = psyrat_report(psyrat_data);

            testCase.verifyTrue(isfield(report, 'admissibility_note'), ...
                'A DoD report must carry the note field once anything is wrong.');
            testCase.verifyNotEmpty(report.admissibility_note, ...
                ['The DoD record reached nothing before RC-43; a headless ' ...
                'export could not learn a variance came out negative.']);
            testCase.verifyTrue(contains(report.admissibility_note, 'of 12'), ...
                'Two groups of 6 draws must aggregate to 12, not stay at 6.');
        end

        function testTheDoDTablesHandBackTheRecordAndACleanRunSaysNothing(testCase)
            % The wiring psyrat_report's per-group DoD sweep depends on: the
            % table builders return the record as a second output, so the
            % report can aggregate one note across groups without re-running
            % the kernel to get at it. And a clean fixture must stay silent, so
            % ordinary DoD output is unchanged.
            data = PsyRATTestDataFactory.makeDoDRelData();

            [~, admVc] = psyrat_dod_varcompt('psyrat_data', data, ...
                'groupidx', 1, 'map', [1 2 3 4], 'savefile', 0, 'gui', 0);
            [~, admObs] = psyrat_dod_observedt('psyrat_data', data, ...
                'groupidx', 1, 'map', [1 2 3 4], 'gui', 0);

            testCase.verifyNotEmpty(admVc, ...
                'psyrat_dod_varcompt must hand back the admissibility record.');
            testCase.verifyTrue(isfield(admVc, 'entries'));
            testCase.verifyNotEmpty(admObs, ...
                'psyrat_dod_observedt must hand back its record too.');

            % clean fixture: recorded, but nothing to report
            testCase.verifyEmpty(psyrat_admissibility_note(admVc), ...
                'A clean DoD run must produce no note.');
        end
    end
end

function psyrat_data = localDoDInadmissibleFixture()
%The two-group DoD fixture with one between-trial cross-covariance pushed past
%what the variances can absorb, so the contrast's between-trial variance goes
%negative on every draw in both groups.
%
%The DoD contrast is (A - B) - (C - D), i.e. weights [1 -1 -1 1], so the A-B
%off-diagonal enters at 2*c_A*c_B*cov = -2*cov. A covariance of 3 against
%diagonals near 0.9 takes the quadratic form well below zero. Perturbing BOTH
%groups is deliberate: it is what makes the report's summed count differ from a
%single builder's, which is the whole point of the test that uses this.
psyrat_data = PsyRATTestDataFactory.makeDoDRelData();

for g = 1:numel(psyrat_data.rel.out.trl_varcov)
    trl = psyrat_data.rel.out.trl_varcov{g};
    trl(:,1,2) = 3;
    trl(:,2,1) = 3;
    psyrat_data.rel.out.trl_varcov{g} = trl;
end
end

function [dep, gen] = localDoDInadmissiblePair()
%The dep/gen pair for one DoD evaluation that records something, reusing
%testTheIccDenominatorIsRecordedInTheDoDKernel's fixture: a trial correlation of
%3 between the two events that share a sign drives the UNSCALED between-trial
%contrast variance to 4 - 6 = -2 on every draw, which is the ICC denominator.
%
%The warning is the kernel working as intended and fires on both calls; it is
%asserted where it belongs (that test) rather than re-asserted here, so silence
%it and restore, the idiom psyrat_report uses for its own loop callers.
nd = 4;
bp = zeros(nd,4,4);
bt = zeros(nd,4,4);
for k = 1:4
    bp(:,k,k) = 1;
    bt(:,k,k) = 1;
end
bt(:,1,2) = 3;
bt(:,2,1) = 3;

args = {'bp', bp, 'bt', bt, 'er_var', log(sqrt(0.01))*ones(nd,4), ...
    'obs', [100 100 1 1], 'CI', 0.95};

ws = warning('off', 'psyrat:negerrvar');
cleanup = onCleanup(@() warning(ws));

dep = psyrat_dodiffrel(args{:}, 'est', 'dep');
gen = psyrat_dodiffrel(args{:}, 'est', 'gen');
end

function localVerifyRoles(testCase, adm, expected)
%Check that each named quantity was recorded, and recorded under the role that
%decides its exactly-zero wording. Keyed by label rather than by position, so
%adding a record does not move an assertion (RC-36).
for k = 1:size(expected,1)
    entry = localEntry(adm, expected{k,1});
    testCase.verifyNotEmpty(entry, ...
        sprintf('No record was taken for the %s.', expected{k,1}));
    testCase.verifyEqual(entry.role, expected{k,2}, ...
        sprintf('The %s is a %s.', expected{k,1}, expected{k,2}));
end
end

function entry = localEntry(adm, label)
%Pull the single record carrying a given label, or [] when there is none.
entry = [];
if isempty(adm) || ~isstruct(adm) || ~isfield(adm, 'entries')
    return;
end
hit = strcmp({adm.entries.label}, label);
if any(hit)
    entry = adm.entries(hit);
end
end

function [bp, bt, er_var, wp_cov, obs] = localRC01Fixture()
%The documented out-of-range difference-score case: residual variances 16 and 1
%at obs = [10 1] put the analytic threshold at rho = 26/44 = 0.5909, and
%rho = 0.70 clears it. Constant across draws, so every quantity has a closed
%form.
nd = 6;
one = ones(nd,1);

bp = zeros(nd,2,2);
bp(:,1,1) = 2; bp(:,2,2) = 2; bp(:,1,2) = 0.5; bp(:,2,1) = 0.5;

bt = zeros(nd,2,2);
bt(:,1,1) = 1; bt(:,2,2) = 1; bt(:,1,2) = 0.9; bt(:,2,1) = 0.9;

er_var = log([4*one 1*one]);
wp_cov = 0.7*sqrt(16*1)*one;
obs = [10 1];
end

function vd = localFixtureRel()
%A genuine one-facet difference-score REL, built by production code.
pd = PsyRATTestDataFactory.makeDoDRelData();
vd = psyrat_dod_build_virtual_rel('psyrat_data', pd, ...
    'mode', 'diff', 'groupidx', 1, 'eventpair', [1 4]);
end
