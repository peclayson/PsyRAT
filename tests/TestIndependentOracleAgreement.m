classdef TestIndependentOracleAgreement < matlab.unittest.TestCase
    %TestIndependentOracleAgreement  Check production against a PAPER-SOURCED oracle.
    %
    %WHAT MAKES THIS DIFFERENT FROM TestCalculationAccuracyOracle
    %   TestCalculationAccuracyOracle compares production against PsyRATAccuracyOracle,
    %   whose expressions are the same algebra as the production functions. That is a
    %   REGRESSION test: it proves the formulas have not changed. It cannot prove they are
    %   right, because both sides would have to be wrong in the same way for it to notice.
    %
    %   This class compares production against PsyRATIndependentOracle, which was
    %   transcribed from the published tables of Rocha et al. (2026) by a reader working
    %   from the page alone, with no access to the PsyRAT source. Agreement here is
    %   EXTERNAL corroboration of the formulas, not merely evidence that nothing moved.
    %
    %IF A TEST HERE FAILS
    %   Do NOT edit PsyRATIndependentOracle to match production. That would convert it into
    %   a second mirror and destroy the only independent check in the repository. A failure
    %   means one of three things, and which one is a scientific judgment for the
    %   maintainer: (a) a deliberate estimand choice PsyRAT makes that the table does not
    %   specify, (b) a transcription error in the oracle, or (c) a defect in production.
    %
    %TEST DESIGN
    %   Comparisons run at several bases, not one. The headline basis uses DISTINCT PRIMES
    %   for every component so that any omitted term, duplicated term, or wrong divisor
    %   changes the result -- with round numbers, two different term sets can coincide by
    %   accident and a real error can hide. n'_i and n'_o are deliberately different so a
    %   divisor swap is detectable; at n'_i == n'_o that mutation is invisible.
    %
    %   Production takes STANDARD DEVIATIONS and squares them internally; the oracle takes
    %   VARIANCES. Every call site below converts explicitly, so a unit error cannot pass
    %   silently between the two.
    %
    %   Copyright (C) 2016-2025 Peter E. Clayson
    %
    %   This program is free software: you can redistribute it and/or modify it under the
    %   terms of the GNU General Public License as published by the Free Software
    %   Foundation, either version 3 of the License, or any later version.
    %
    %   This program is distributed in the hope that it will be useful, but WITHOUT ANY
    %   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    %   PARTICULAR PURPOSE. See the GNU General Public License for more details.
    %
    %   You should have received a copy of the GNU General Public License along with this
    %   program (gpl.txt). If not, see <http://www.gnu.org/licenses/>.

    properties (Constant)
        %Machine-precision agreement is the right bar: both sides evaluate a closed form on
        %identical numbers, so anything above rounding is a real disagreement.
        TOL = 1e-12
        NDRAW = 200      %degenerate draws (all identical) => posterior mean is exact
    end

    methods (Static)
        function b = bases()
            %Several parameter sets, so agreement is not an artifact of one lucky point.
            %  prime      distinct primes; any dropped/duplicated term or wrong divisor moves it
            %  sentinel   D and G land on exact decimals (0.5 and 0.8) at n' = 1
            %  lopsided   large facet effects relative to person variance
            %  tiny_i     sigma2_i near zero -- the LEGITIMATE case where D and G coincide
            b = struct( ...
                'prime',    struct('p',2,'o',3,'i',5,'po',7,'pi',11,'oi',13,'poi_e',17,'ni',2,'no',3), ...
                'sentinel', struct('p',1,'o',0.2,'i',0.75,'po',0.1,'pi',0.15,'oi',0.05,'poi_e',0.25,'ni',1,'no',1), ...
                'lopsided', struct('p',0.5,'o',9,'i',11,'po',0.25,'pi',0.3,'oi',7,'poi_e',13,'ni',7,'no',5), ...
                'tiny_i',   struct('p',3,'o',2,'i',1e-9,'po',1.5,'pi',2.5,'oi',1e-9,'poi_e',4,'ni',3,'no',4));
        end

        function b = diffBases()
            %Bases for the Rocha Table 6 difference-score coefficients.
            %
            %Difference scores need a per-MEASURE parameter set, not a per-facet one: each
            %component appears three times (X's variance, Y's variance, their covariance)
            %and each facet has two replicate counts (one per measure). Field names carry
            %the measure suffix so a transposition is visible when a value is read back.
            %
            %EVERY BASIS USES UNEQUAL REPLICATE COUNTS, and that is the point. At
            %niX == niY the harmonic, geometric and arithmetic means all coincide, so a
            %wrong covariance divisor is arithmetically invisible -- which is exactly the
            %blindness RC-16 found in the existing external difference-score fixture. Only
            %unequal counts can anchor the divisor at all.
            %
            %  prime      distinct primes for the variances, so any dropped, duplicated or
            %             misplaced term moves the result; small distinct covariances keep
            %             every printed block positive
            %  persondom  person variance dominates the facet effects, i.e. the regime a
            %             real ERP difference score is estimated in, with coefficients high
            %             enough that a small error would still leave a plausible number
            %  zerocov    every covariance zero, all seven families -- the NON-CONCURRENT
            %             default, where PsyRAT omits wp_cov/er_cov entirely. Included
            %             even though it is BLIND to the divisor (0/anything == 0): its job
            %             is to check that the rest of the formula is right when the
            %             covariance terms drop out, not to pin the harmonic mean.
            b = struct( ...
                'prime', struct( ...
                    'pX',2,   'pY',3,   'cov_p',0.5, ...
                    'poX',5,  'poY',7,  'cov_po',0.11, ...
                    'piX',11, 'piY',13, 'cov_pi',0.17, ...
                    'poiX_e',17, 'poiY_e',19, 'cov_poi_e',0.23, ...
                    'oX',23,  'oY',29,  'cov_o',0.29, ...
                    'iX',31,  'iY',37,  'cov_i',0.31, ...
                    'oiX',41, 'oiY',43, 'cov_oi',0.37, ...
                    'niX',13, 'niY',29, 'noX',2, 'noY',5), ...
                'persondom', struct( ...
                    'pX',20,  'pY',24,  'cov_p',8, ...
                    'poX',1,  'poY',1.5,'cov_po',0.3, ...
                    'piX',2,  'piY',2.5,'cov_pi',0.4, ...
                    'poiX_e',3, 'poiY_e',3.5, 'cov_poi_e',0.5, ...
                    'oX',0.5, 'oY',0.6, 'cov_o',0.1, ...
                    'iX',0.7, 'iY',0.8, 'cov_i',0.15, ...
                    'oiX',0.9,'oiY',1.1,'cov_oi',0.2, ...
                    'niX',20, 'niY',40, 'noX',2, 'noY',4), ...
                'zerocov', struct( ...
                    'pX',6,   'pY',9,   'cov_p',0, ...
                    'poX',1.5,'poY',2.5,'cov_po',0, ...
                    'piX',3.5,'piY',4.5,'cov_pi',0, ...
                    'poiX_e',5.5, 'poiY_e',6.5, 'cov_poi_e',0, ...
                    'oX',0.8, 'oY',1.2, 'cov_o',0, ...
                    'iX',1.4, 'iY',1.6, 'cov_i',0, ...
                    'oiX',1.8,'oiY',2.2,'cov_oi',0, ...
                    'niX',7,  'niY',11, 'noX',3, 'noY',5));
        end
    end

    methods (Test)

        function testOneFacetGlobalCoefficients(tc)
            %Rocha Table 2: the G/D pair. This is the core absolute-vs-relative claim --
            %same numerator, and dependability's denominator carries sigma2_i/n' where
            %generalizability's does not.
            bs = TestIndependentOracleAgreement.bases();
            for f = fieldnames(bs)'
                b = bs.(f{1});
                [bp, wp, ii] = tc.oneFacetSDs(b);

                expG = PsyRATIndependentOracle.gCoefficientOneFacet(b.p, b.poi_e, b.ni);
                expD = PsyRATIndependentOracle.dCoefficientOneFacet(b.p, b.poi_e, b.i, b.ni);

                [~, gotG, ~] = psyrat_rel_sing('gcoeff',2,'metric','global', ...
                    'bp',bp,'wp',wp,'i',ii,'obs',b.ni,'CI',0.95);
                [~, gotD, ~] = psyrat_rel_sing('gcoeff',1,'metric','global', ...
                    'bp',bp,'wp',wp,'i',ii,'obs',b.ni,'CI',0.95);

                tc.verifyEqual(gotG, expG, 'AbsTol', tc.TOL, ...
                    sprintf('one-facet G disagrees with the paper at basis "%s"', f{1}));
                tc.verifyEqual(gotD, expD, 'AbsTol', tc.TOL, ...
                    sprintf('one-facet D disagrees with the paper at basis "%s"', f{1}));

                %Ordering must hold wherever sigma2_i > 0, and the "tiny_i" basis is the
                %documented legitimate coincidence, kept here so the suite records that
                %D == G is EXPECTED there rather than treating it as a failure.
                if b.i > 1e-6
                    tc.verifyLessThan(gotD, gotG, ...
                        sprintf('D must be strictly below G when sigma2_i > 0 (basis "%s")', f{1}));
                end
            end
        end

        function testOneFacetIccs(tc)
            %The ICC rows carry no divisor in the printed table: they are the
            %single-observation forms, not the global coefficients evaluated at n'.
            bs = TestIndependentOracleAgreement.bases();
            for f = fieldnames(bs)'
                b = bs.(f{1});
                [bp, wp, ii] = tc.oneFacetSDs(b);
                [~, gotRel, ~] = psyrat_rel_sing('gcoeff',2,'metric','icc', ...
                    'bp',bp,'wp',wp,'i',ii,'CI',0.95);
                [~, gotAbs, ~] = psyrat_rel_sing('gcoeff',1,'metric','icc', ...
                    'bp',bp,'wp',wp,'i',ii,'CI',0.95);
                tc.verifyEqual(gotRel, ...
                    PsyRATIndependentOracle.iccRelativeOneFacet(b.p, b.poi_e), ...
                    'AbsTol', tc.TOL, sprintf('one-facet relative ICC, basis "%s"', f{1}));
                tc.verifyEqual(gotAbs, ...
                    PsyRATIndependentOracle.iccAbsoluteOneFacet(b.p, b.poi_e, b.i), ...
                    'AbsTol', tc.TOL, sprintf('one-facet absolute ICC, basis "%s"', f{1}));
            end
        end

        function testTwoFacetAllSixCoefficients(tc)
            %Rocha Table 3: CE / CS / CES x {G, D}. The three coefficient types differ in
            %which facet is fixed, and the two D-coefficients that have a fixed facet drop
            %that facet's main effect from the absolute error (CE omits sigma2_o, CS omits
            %sigma2_i). Those omissions are printed in the table and are the sharpest
            %available test that production implements the right partition.
            bs = TestIndependentOracleAgreement.bases();
            for f = fieldnames(bs)'
                b = bs.(f{1});
                a = tc.twoFacetArgs(b);

                exp2 = { ...
                    1, 2, PsyRATIndependentOracle.ceG(b.p,b.po,b.pi,b.poi_e,b.ni,b.no), 'CE G'; ...
                    1, 1, PsyRATIndependentOracle.ceD(b.p,b.po,b.pi,b.poi_e,b.i,b.oi,b.ni,b.no), 'CE D'; ...
                    2, 2, PsyRATIndependentOracle.csG(b.p,b.pi,b.po,b.poi_e,b.ni,b.no), 'CS G'; ...
                    2, 1, PsyRATIndependentOracle.csD(b.p,b.pi,b.po,b.poi_e,b.o,b.oi,b.ni,b.no), 'CS D'; ...
                    3, 2, PsyRATIndependentOracle.cesG(b.p,b.po,b.pi,b.poi_e,b.ni,b.no), 'CES G'; ...
                    3, 1, PsyRATIndependentOracle.cesD(b.p,b.po,b.pi,b.poi_e,b.o,b.i,b.oi,b.ni,b.no), 'CES D'};

                for k = 1:size(exp2,1)
                    [~, got, ~] = psyrat_rel_trt('gcoeff',exp2{k,2},'reltype',exp2{k,1}, ...
                        a{:}, 'obs',b.ni,'nocc',b.no,'CI',0.95);
                    tc.verifyEqual(got, exp2{k,3}, 'AbsTol', tc.TOL, ...
                        sprintf('%s disagrees with the paper at basis "%s"', exp2{k,4}, f{1}));
                end
            end
        end

        function testTwoFacetIccs(tc)
            %The printed ICC rows are the single-observation forms, so production's
            %analogue is psyrat_rel_trt evaluated at obs = 1, nocc = 1.
            bs = TestIndependentOracleAgreement.bases();
            for f = fieldnames(bs)'
                b = bs.(f{1});
                a = tc.twoFacetArgs(b);
                expI = { ...
                    1, 2, PsyRATIndependentOracle.ceIccRelative(b.p,b.po,b.pi,b.poi_e), 'CE ICC rel'; ...
                    1, 1, PsyRATIndependentOracle.ceIccAbsolute(b.p,b.po,b.pi,b.poi_e,b.i,b.oi), 'CE ICC abs'; ...
                    2, 2, PsyRATIndependentOracle.csIccRelative(b.p,b.pi,b.po,b.poi_e), 'CS ICC rel'; ...
                    2, 1, PsyRATIndependentOracle.csIccAbsolute(b.p,b.pi,b.po,b.poi_e,b.o,b.oi), 'CS ICC abs'; ...
                    3, 2, PsyRATIndependentOracle.cesIccRelative(b.p,b.po,b.pi,b.poi_e), 'CES ICC rel'; ...
                    3, 1, PsyRATIndependentOracle.cesIccAbsolute(b.p,b.po,b.pi,b.poi_e,b.o,b.i,b.oi), 'CES ICC abs'};
                for k = 1:size(expI,1)
                    [~, got, ~] = psyrat_rel_trt('gcoeff',expI{k,2},'reltype',expI{k,1}, ...
                        a{:}, 'obs',1,'nocc',1,'CI',0.95);
                    tc.verifyEqual(got, expI{k,3}, 'AbsTol', tc.TOL, ...
                        sprintf('%s disagrees with the paper at basis "%s"', expI{k,4}, f{1}));
                end
            end
        end

        function testSubjectLevelTwoFacetIccsMatchPaper(tc)
            %RC-06. The SUBJECT-level two-facet ICC (psyrat_ssrel_trt) must be the
            %same Rocha Table 3 ICC row as the group-level one, with that
            %participant's residual variance standing in for sigma2_poi,e. Before
            %the fix it routed through psyrat_icc, which sees only sigma2_p and the
            %residual, so all six (reltype x gcoeff) cells returned one identical
            %number and the "relative" and "absolute" columns could not differ.
            %
            %The fixture mirrors the real caller (psyrat_relsummary's
            %local_ssrel_trt_call): ONLY the residual is per-subject, while
            %txp_ss/oxp_ss broadcast the population person x trial / person x
            %occasion SDs across subjects.
            %
            %Two properties are pinned here that a formula check alone would miss:
            %  1. nocc is deliberately passed as 3, and the expectation stays the
            %     n'_o = 1 form -- the ICC must NOT scale with a multi-occasion
            %     composite (psyrat_relsummary.m:5082-5086).
            %  2. s3 carries s1's residual with a different trial count, so its ICC
            %     must equal s1's -- the ICC must NOT scale with n'_i either.
            n = tc.NDRAW;
            ssPoiE = [0.36; 2.25; 0.36];    %per-subject residual VARIANCES
            trls   = [4; 9; 20];
            idtable = table({'s1';'s2';'s3'}, (1:3)', trls, ...
                'VariableNames', {'id','id2','trls'});
            errSs = [repmat(sqrt(ssPoiE(1)),n,1), ...
                     repmat(sqrt(ssPoiE(2)),n,1), ...
                     repmat(sqrt(ssPoiE(3)),n,1)];

            bs = TestIndependentOracleAgreement.bases();
            for f = fieldnames(bs)'
                b = bs.(f{1});
                a = tc.twoFacetArgs(b);
                expI = { ...
                    1, 2, @(e) PsyRATIndependentOracle.ceIccRelative(b.p,b.po,b.pi,e), 'CE ICC rel'; ...
                    1, 1, @(e) PsyRATIndependentOracle.ceIccAbsolute(b.p,b.po,b.pi,e,b.i,b.oi), 'CE ICC abs'; ...
                    2, 2, @(e) PsyRATIndependentOracle.csIccRelative(b.p,b.pi,b.po,e), 'CS ICC rel'; ...
                    2, 1, @(e) PsyRATIndependentOracle.csIccAbsolute(b.p,b.pi,b.po,e,b.o,b.oi), 'CS ICC abs'; ...
                    3, 2, @(e) PsyRATIndependentOracle.cesIccRelative(b.p,b.po,b.pi,e), 'CES ICC rel'; ...
                    3, 1, @(e) PsyRATIndependentOracle.cesIccAbsolute(b.p,b.po,b.pi,e,b.o,b.i,b.oi), 'CES ICC abs'};

                for k = 1:size(expI,1)
                    out = psyrat_ssrel_trt('gcoeff',expI{k,2},'reltype',expI{k,1}, ...
                        a{:}, 'idtable',idtable,'CI',0.95,'nocc',3, ...
                        'txp_ss',repmat(sqrt(b.pi),n,3), ...
                        'oxp_ss',repmat(sqrt(b.po),n,3), ...
                        'err_ss',errSs);

                    for r = 1:height(idtable)
                        tc.verifyEqual(out.icc_pt(r), expI{k,3}(ssPoiE(r)), ...
                            'AbsTol', tc.TOL, sprintf( ...
                            'subject-level %s, subject %d, disagrees with the paper at basis "%s"', ...
                            expI{k,4}, r, f{1}));
                    end

                    tc.verifyEqual(out.icc_pt(3), out.icc_pt(1), 'AbsTol', tc.TOL, ...
                        sprintf(['subject-level %s scales with the trial count at basis ' ...
                        '"%s"; an ICC is the single-observation coefficient'], expI{k,4}, f{1}));
                end

                %The two coefficient types must actually separate whenever some
                %facet-only variance is present, which is the user-visible symptom
                %RC-06 reported. tiny_i is excluded: sigma2_i = sigma2_oi = 1e-9 is
                %the LEGITIMATE case where relative and absolute coincide.
                if ~strcmp(f{1},'tiny_i')
                    for rt = 1:3
                        argsRt = {a{:}, 'idtable',idtable,'CI',0.95,'nocc',3, ...
                            'txp_ss',repmat(sqrt(b.pi),n,3), ...
                            'oxp_ss',repmat(sqrt(b.po),n,3), ...
                            'err_ss',errSs}; %#ok<CCAT>
                        oRel = psyrat_ssrel_trt('gcoeff',2,'reltype',rt,argsRt{:});
                        oAbs = psyrat_ssrel_trt('gcoeff',1,'reltype',rt,argsRt{:});
                        tc.verifyNotEqual(oAbs.icc_pt, oRel.icc_pt, sprintf( ...
                            ['relative and absolute subject-level ICCs are identical at ' ...
                            'reltype %d, basis "%s"'], rt, f{1}));
                    end
                end
            end
        end

        function testTwoFacetRangeBranchMatchesScalarBranchAndPaper(tc)
            %psyrat_rel_trt implements the six Table 3 formulas TWICE: a scalar branch
            %(fires when obs is a scalar) and a hand-duplicated range branch (fires when
            %obs is a 2-element [lo hi] vector). They are currently identical, but they are
            %maintained separately, so a correction applied to one and not the other would
            %be invisible to any check that only ever passes a scalar.
            %
            %That blind spot was real: the external corroboration this class is built on
            %compares against closed-form table values, which requires scalar inputs, so it
            %exercised only the scalar copy. The range branch is nonetheless what production
            %uses to build the reliability-vs-trials curve (psyrat_relsummary.m:3294 and
            %siblings), i.e. a user-visible surface.
            %
            %This test closes that gap: every column of the range output is checked against
            %the paper-sourced oracle at its own n'_i, and against the scalar branch called
            %at the same n'_i.
            b = TestIndependentOracleAgreement.bases().prime;
            a = tc.twoFacetArgs(b);
            lo = 1; hi = 6;

            spec = { ...
                1, 2, @(n) PsyRATIndependentOracle.ceG(b.p,b.po,b.pi,b.poi_e,n,b.no),  'CE G'; ...
                1, 1, @(n) PsyRATIndependentOracle.ceD(b.p,b.po,b.pi,b.poi_e,b.i,b.oi,n,b.no), 'CE D'; ...
                2, 2, @(n) PsyRATIndependentOracle.csG(b.p,b.pi,b.po,b.poi_e,n,b.no),  'CS G'; ...
                2, 1, @(n) PsyRATIndependentOracle.csD(b.p,b.pi,b.po,b.poi_e,b.o,b.oi,n,b.no), 'CS D'; ...
                3, 2, @(n) PsyRATIndependentOracle.cesG(b.p,b.po,b.pi,b.poi_e,n,b.no), 'CES G'; ...
                3, 1, @(n) PsyRATIndependentOracle.cesD(b.p,b.po,b.pi,b.poi_e,b.o,b.i,b.oi,n,b.no), 'CES D'};

            for k = 1:size(spec,1)
                rt = spec{k,1}; gc = spec{k,2}; oracleFn = spec{k,3}; name = spec{k,4};

                [~, gotRange, ~] = psyrat_rel_trt('gcoeff',gc,'reltype',rt, a{:}, ...
                    'obs',[lo hi],'nocc',b.no,'CI',0.95);

                tc.verifyNumElements(gotRange, hi-lo+1, sprintf( ...
                    '%s: range output must have one column per trial count in obs(1):obs(2)', name));

                for n = lo:hi
                    idx = n - lo + 1;   %columns align to obs(1):obs(2), not to 1:hi
                    tc.verifyEqual(gotRange(idx), oracleFn(n), 'AbsTol', tc.TOL, ...
                        sprintf('%s range branch disagrees with the paper at n''_i = %d', name, n));

                    [~, gotScalar, ~] = psyrat_rel_trt('gcoeff',gc,'reltype',rt, a{:}, ...
                        'obs',n,'nocc',b.no,'CI',0.95);
                    tc.verifyEqual(gotRange(idx), gotScalar, 'AbsTol', tc.TOL, sprintf( ...
                        ['%s: the range branch and the scalar branch disagree at n''_i = %d. ' ...
                         'These are two hand-maintained copies of the same formula; a fix ' ...
                         'applied to only one of them produces exactly this failure.'], name, n));
                end
            end
        end

        function testOneFacetRangeBranchMatchesScalarBranchAndPaper(tc)
            %psyrat_rel_sing has the same two-copy structure as psyrat_rel_trt (a scalar
            %path and an obs = [lo hi] loop), and the same blind spot applies. The range
            %branch drives the one-facet reliability-vs-trials curve.
            b = TestIndependentOracleAgreement.bases().prime;
            [bp, wp, ii] = tc.oneFacetSDs(b);
            lo = 2; hi = 7;   %lo > 1 deliberately: the column indexing must align to
                              %obs(1):obs(2), and an off-by-one is only visible when lo ~= 1.

            for gc = [1 2]
                [~, gotRange, ~] = psyrat_rel_sing('gcoeff',gc,'metric','global', ...
                    'bp',bp,'wp',wp,'i',ii,'obs',[lo hi],'CI',0.95);
                tc.verifyNumElements(gotRange, hi-lo+1, ...
                    'one-facet range output must have one column per requested trial count');

                for n = lo:hi
                    idx = n - lo + 1;
                    if gc == 1
                        expected = PsyRATIndependentOracle.dCoefficientOneFacet(b.p, b.poi_e, b.i, n);
                    else
                        expected = PsyRATIndependentOracle.gCoefficientOneFacet(b.p, b.poi_e, n);
                    end
                    tc.verifyEqual(gotRange(idx), expected, 'AbsTol', tc.TOL, sprintf( ...
                        'one-facet gcoeff=%d range branch disagrees with the paper at n''_i = %d', gc, n));

                    [~, gotScalar, ~] = psyrat_rel_sing('gcoeff',gc,'metric','global', ...
                        'bp',bp,'wp',wp,'i',ii,'obs',n,'CI',0.95);
                    tc.verifyEqual(gotRange(idx), gotScalar, 'AbsTol', tc.TOL, sprintf( ...
                        'one-facet gcoeff=%d: range and scalar branches disagree at n''_i = %d', gc, n));
                end
            end
        end

        function testCutScoreRangeBranchMatchesScalarBranchAndPaper(tc)
            %The cut-score functions carry the same duplicated range branch. Since a
            %cut-score must always be the absolute-error quantity, a divergence between the
            %two copies could silently make the curve relative while the reported scalar
            %stays absolute.
            b = TestIndependentOracleAgreement.bases().prime;
            [bp, wp, ~] = tc.oneFacetSDs(b);
            mu = 1.3; C = 4.7;
            muD = repmat(mu, tc.NDRAW, 1);
            lo = 2; hi = 6;

            [~, gotRange, ~] = psyrat_cutscore_sing('gcoeff',1,'bp',bp,'wp',wp, ...
                'mu',muD,'cut',C,'obs',[lo hi],'CI',0.95);
            tc.verifyNumElements(gotRange, hi-lo+1, ...
                'cut-score range output must have one column per requested trial count');

            for n = lo:hi
                idx = n - lo + 1;
                tc.verifyEqual(gotRange(idx), ...
                    PsyRATIndependentOracle.cutScoreOneFacet(b.p,b.poi_e,b.i,n,mu,C), ...
                    'AbsTol', tc.TOL, sprintf( ...
                    'cut-score range branch disagrees with the paper at n''_i = %d', n));

                [~, gotScalar, ~] = psyrat_cutscore_sing('gcoeff',1,'bp',bp,'wp',wp, ...
                    'mu',muD,'cut',C,'obs',n,'CI',0.95);
                tc.verifyEqual(gotRange(idx), gotScalar, 'AbsTol', tc.TOL, sprintf( ...
                    'cut-score range and scalar branches disagree at n''_i = %d', n));
            end
        end

        function testCutScoresAreAbsoluteErrorAndMatchThePaper(tc)
            %The invariant the maintainer asked to have pinned: a cut-score is always a
            %dependability (absolute-error) quantity. Rocha prints exactly one cut-score
            %schema, instantiated once per univariate design, with no relative variant --
            %in every table the relative cell is blank.
            bs = TestIndependentOracleAgreement.bases();
            for f = fieldnames(bs)'
                b = bs.(f{1});
                mu = 1.3; C = 4.7;
                [bp, wp, ~] = tc.oneFacetSDs(b);
                muD = repmat(mu, tc.NDRAW, 1);

                [~, got1, ~] = psyrat_cutscore_sing('gcoeff',1,'bp',bp,'wp',wp, ...
                    'mu',muD,'cut',C,'obs',b.ni,'CI',0.95);
                tc.verifyEqual(got1, ...
                    PsyRATIndependentOracle.cutScoreOneFacet(b.p,b.poi_e,b.i,b.ni,mu,C), ...
                    'AbsTol', tc.TOL, sprintf('one-facet cut-score, basis "%s"', f{1}));

                a = tc.twoFacetArgs(b);
                [~, got2, ~] = psyrat_cutscore_trt('gcoeff',1,'reltype',3, a{:}, ...
                    'mu',muD,'cut',C,'obs',b.ni,'nocc',b.no,'CI',0.95);
                tc.verifyEqual(got2, ...
                    PsyRATIndependentOracle.cutScoreTwoFacet(b.p,b.po,b.pi,b.poi_e, ...
                        b.o,b.i,b.oi,b.ni,b.no,mu,C), ...
                    'AbsTol', tc.TOL, sprintf('two-facet CES cut-score, basis "%s"', f{1}));
            end
        end

        function testCutScoreRejectsRelativeSelector(tc)
            %A relative-error cut score is not a quantity the paper defines, so requesting
            %one must FAIL LOUDLY rather than be silently coerced to absolute -- a silent
            %coercion would hand the user a number for a question that has no answer.
            b = TestIndependentOracleAgreement.bases().prime;
            [bp, wp, ~] = tc.oneFacetSDs(b);
            muD = repmat(1.3, tc.NDRAW, 1);
            tc.verifyError(@() psyrat_cutscore_sing('gcoeff',2,'bp',bp,'wp',wp, ...
                'mu',muD,'cut',4.7,'obs',b.ni,'CI',0.95), 'varargin:gcoeff');
            a = tc.twoFacetArgs(b);
            tc.verifyError(@() psyrat_cutscore_trt('gcoeff',2,'reltype',3, a{:}, ...
                'mu',muD,'cut',4.7,'obs',b.ni,'nocc',b.no,'CI',0.95), 'varargin:gcoeff');
        end

        function testCutScoreAtTheMeanReducesToDependability(tc)
            %Rocha states this in prose on p.7: when the cut score matches the scale mean,
            %the cut-score coefficient equals the global dependability coefficient. It is a
            %strong check because it must hold EXACTLY, and it fails if the cut-score error
            %term is the relative one -- so it discriminates D from G rather than merely
            %confirming a number.
            bs = TestIndependentOracleAgreement.bases();
            for f = fieldnames(bs)'
                b = bs.(f{1});
                mu = 2.75;
                [bp, wp, ii] = tc.oneFacetSDs(b);
                muD = repmat(mu, tc.NDRAW, 1);

                [~, atMean, ~] = psyrat_cutscore_sing('gcoeff',1,'bp',bp,'wp',wp, ...
                    'mu',muD,'cut',mu,'obs',b.ni,'CI',0.95);
                [~, globalD, ~] = psyrat_rel_sing('gcoeff',1,'metric','global', ...
                    'bp',bp,'wp',wp,'i',ii,'obs',b.ni,'CI',0.95);
                [~, globalG, ~] = psyrat_rel_sing('gcoeff',2,'metric','global', ...
                    'bp',bp,'wp',wp,'i',ii,'obs',b.ni,'CI',0.95);

                tc.verifyEqual(atMean, globalD, 'AbsTol', tc.TOL, ...
                    sprintf('cut==mu must equal global dependability (basis "%s")', f{1}));
                if b.i > 1e-6
                    tc.verifyNotEqual(atMean, globalG, ...
                        'cut==mu must NOT equal global generalizability when sigma2_i > 0');
                end
            end
        end

        function testAbsoluteOnlyComponentsReachDependabilityOnly(tc)
            %The component-sensitivity fingerprint, which tests the ESTIMAND rather than a
            %single number: perturb one component that belongs to absolute error only, and
            %dependability must move while generalizability stays bitwise identical.
            %
            %Run at nocc = 3 deliberately. At nocc = 1 the occasion main effect and the
            %occasion-by-trial interaction are not identified and D approaches G by
            %accident, so a selector that was accepted and then ignored would look correct.
            b = TestIndependentOracleAgreement.bases().prime;
            a = tc.twoFacetArgs(b);
            bump = @(fld, mult) tc.twoFacetArgs(setfield(b, fld, b.(fld)*mult));

            %reltype -> which facet-only components its ABSOLUTE error carries.
            %CE fixes occasion (so sigma2_o is absent); CS fixes trial (sigma2_i absent).
            spec = {1, {'o',false}, {'i',true},  {'oi',true}; ...
                    2, {'o',true},  {'i',false}, {'oi',true}; ...
                    3, {'o',true},  {'i',true},  {'oi',true}};

            for r = 1:size(spec,1)
                rt = spec{r,1};
                [~, baseD, ~] = psyrat_rel_trt('gcoeff',1,'reltype',rt, a{:}, ...
                    'obs',b.ni,'nocc',b.no,'CI',0.95);
                [~, baseG, ~] = psyrat_rel_trt('gcoeff',2,'reltype',rt, a{:}, ...
                    'obs',b.ni,'nocc',b.no,'CI',0.95);

                for c = 2:4
                    fld = spec{r,c}{1}; mustMove = spec{r,c}{2};
                    a2 = bump(fld, 2);
                    [~, pD, ~] = psyrat_rel_trt('gcoeff',1,'reltype',rt, a2{:}, ...
                        'obs',b.ni,'nocc',b.no,'CI',0.95);
                    [~, pG, ~] = psyrat_rel_trt('gcoeff',2,'reltype',rt, a2{:}, ...
                        'obs',b.ni,'nocc',b.no,'CI',0.95);

                    if mustMove
                        tc.verifyNotEqual(pD, baseD, sprintf( ...
                            ['reltype %d dependability must respond to sigma2_%s; ' ...
                             'no response means the term is missing from its absolute error'], rt, fld));
                    else
                        tc.verifyEqual(pD, baseD, 'AbsTol', tc.TOL, sprintf( ...
                            ['reltype %d dependability must NOT respond to sigma2_%s ' ...
                             '(that facet is fixed); a response means an extra term leaked in'], rt, fld));
                    end
                    %Generalizability is never a function of a facet main effect, in any reltype.
                    tc.verifyEqual(pG, baseG, 'AbsTol', tc.TOL, sprintf( ...
                        ['reltype %d generalizability must NOT respond to sigma2_%s -- ' ...
                         'relative error contains only person-carrying components'], rt, fld));
                end
            end
        end

        %% ---------------- Rocha Table 6: difference scores ----------------------
        %
        %  WHY THIS GROUP EXISTS (finding RC-17). Difference-score reliability is where
        %  PsyRAT's users are most exposed: the coefficient is a ratio of small differences
        %  of correlated quantities, so a sign or divisor error is easy to make and hard to
        %  eyeball. Before these tests, no paper-sourced closed form anchored ANY
        %  difference-score coefficient. The suite did contain difference-score checks, but
        %  they were mirrors -- the expected value was the production algebra retyped -- so
        %  a review probe could flip the sign on the between-person covariance, move the
        %  coefficient from 0.9108 to 0.9543, and have eleven regression tests stay green.
        %
        %  These tests are expected to PASS on arrival: production already reproduces
        %  Table 6. Their value is prospective. The covariance divisor they pin is the same
        %  quantity findings RC-01 / S13 are still adjudicating, so if the toolbox later
        %  moves off the harmonic mean, these tests fail LOUDLY and the deviation from the
        %  published formula becomes an explicit, recorded decision instead of a silent one.
        %  When that happens, record the deviation -- do not edit the oracle.

        function testOneFacetDifferenceCoefficientsMatchThePaper(tc)
            %Rocha Table 6, one-facet block (p x i): the G/D pair for a difference score.
            %
            %The unit boundary here is worse than anywhere else in this class, and each of
            %the three conversions below is explicit for that reason. psyrat_diffrel takes
            %its person and trial terms as VARIANCES AND COVARIANCES packed into a
            %[draw x 2 x 2] array, but takes the residual as a LOG STANDARD DEVIATION
            %(psyrat_diffrel.m squares exp(er_var)), while the residual COVARIANCE that sits
            %alongside it is on the plain covariance scale. Three scales, three arguments.
            bs = TestIndependentOracleAgreement.diffBases();
            for f = fieldnames(bs)'
                b = bs.(f{1});

                bp = tc.diffCovArray(b.pX, b.pY, b.cov_p);           %variances + covariance
                bt = tc.diffCovArray(b.iX, b.iY, b.cov_i);           %trial main effect
                er_var = tc.diffLogSD(b.poiX_e, b.poiY_e);           %LOG SD, not variance
                wp_cov = repmat(b.cov_poi_e, tc.NDRAW, 1);           %plain covariance
                obs = [b.niX b.niY];

                expG = PsyRATIndependentOracle.diffGOneFacet( ...
                    b.pX, b.pY, b.cov_p, b.poiX_e, b.poiY_e, b.cov_poi_e, b.niX, b.niY);
                expD = PsyRATIndependentOracle.diffDOneFacet( ...
                    b.pX, b.pY, b.cov_p, b.poiX_e, b.poiY_e, b.cov_poi_e, ...
                    b.iX, b.iY, b.cov_i, b.niX, b.niY);

                gotG = psyrat_diffrel('bp',bp,'bt',bt,'er_var',er_var,'wp_cov',wp_cov, ...
                    'obs',obs,'est','gen','CI',0.95);
                gotD = psyrat_diffrel('bp',bp,'bt',bt,'er_var',er_var,'wp_cov',wp_cov, ...
                    'obs',obs,'est','dep','CI',0.95);

                tc.verifyEqual(gotG.pt, expG, 'AbsTol', tc.TOL, sprintf( ...
                    'one-facet difference G disagrees with Rocha Table 6 at basis "%s"', f{1}));
                tc.verifyEqual(gotD.pt, expD, 'AbsTol', tc.TOL, sprintf( ...
                    'one-facet difference D disagrees with Rocha Table 6 at basis "%s"', f{1}));

                %Sanity floor: BOTH coefficients bounded on BOTH sides. A one-sided pair of
                %checks would miss D > 1, which is the specific pathology the RC-01 residue
                %is about -- a large positive trial-level covariance can drive the absolute
                %error negative, and the equality assertions above would then faithfully
                %confirm the oracle on a coefficient nobody could report.
                tc.verifyGreaterThan([gotG.pt gotD.pt], 0, ...
                    sprintf('basis "%s" produced a non-positive coefficient', f{1}));
                tc.verifyLessThan([gotG.pt gotD.pt], 1, ...
                    sprintf('basis "%s" produced a coefficient above 1', f{1}));
            end
        end

        function testTwoFacetCesDifferenceCoefficientsMatchThePaper(tc)
            %Rocha Table 6, two-facet block (p x i x o). The printed formulas are CES only
            %(both facets random), per the table's own Note; CE and CS are NOT anchored --
            %see the Table 6 header in PsyRATIndependentOracle for why.
            %
            %Occasion counts are unequal as well as trial counts, so the trial and occasion
            %harmonic means are themselves different numbers (prime: 17.95 vs 2.86;
            %persondom: 26.67 vs 2.67) and swapping one for the other changes the result.
            %The condition for that swap to be invisible is hmi == hmo -- NOT noX == noY,
            %which would still leave hmi differing from hmo whenever niX != niY. The zerocov
            %basis is blind to it either way, since its covariance terms are zero.
            bs = TestIndependentOracleAgreement.diffBases();
            for f = fieldnames(bs)'
                b = bs.(f{1});
                [c, n] = tc.diffTwoFacetOracleArgs(b);
                a = tc.diffTwoFacetProductionArgs(b);

                expG = PsyRATIndependentOracle.diffCesGTwoFacet(c, n);
                expD = PsyRATIndependentOracle.diffCesDTwoFacet(c, n);

                gotG = psyrat_diffrel_trt(a{:}, 'reltype',3, 'est','gen', 'CI',0.95);
                gotD = psyrat_diffrel_trt(a{:}, 'reltype',3, 'est','dep', 'CI',0.95);

                tc.verifyEqual(gotG.pt, expG, 'AbsTol', tc.TOL, sprintf( ...
                    'two-facet CES difference G disagrees with Rocha Table 6 at basis "%s"', f{1}));
                tc.verifyEqual(gotD.pt, expD, 'AbsTol', tc.TOL, sprintf( ...
                    'two-facet CES difference D disagrees with Rocha Table 6 at basis "%s"', f{1}));

                %Both coefficients, both bounds -- see the one-facet test for why a
                %one-sided guard would miss the out-of-range case that matters.
                tc.verifyGreaterThan([gotG.pt gotD.pt], 0, ...
                    sprintf('basis "%s" produced a non-positive coefficient', f{1}));
                tc.verifyLessThan([gotG.pt gotD.pt], 1, ...
                    sprintf('basis "%s" produced a coefficient above 1', f{1}));
            end
        end

        function testDifferenceCovarianceDivisorIsTheHarmonicMean(tc)
            %A negative control for the two tests above, and the reason they are worth
            %having at all.
            %
            %An equality test proves agreement but not SENSITIVITY: if the covariance terms
            %happened to contribute nothing, the tests would pass under any divisor and the
            %anchor would be decorative. So recompute the paper's closed form three more
            %times with the harmonic mean replaced by each plausible wrong divisor -- X's
            %own count, Y's own count, and the ARITHMETIC mean, which is the specific
            %alternative the RC-01 / S13 feasibility argument raises -- and require
            %production to disagree with all three.
            %
            %The zerocov basis is excluded on purpose: with every covariance at zero the
            %divisor cannot matter, so demanding disagreement there would be demanding that
            %0/a differ from 0/b.
            b = TestIndependentOracleAgreement.diffBases().prime;

            %VERIFY, not assume. An assumption failure is reported as Incomplete and
            %runAllUnitTests.m:110 treats assumption-filtered skips as EXPECTED, so a basis
            %edited to equal trial counts would silently retire this control while the suite
            %still reported green. That is precisely the failure mode RC-16 documented for
            %the fixture whose coverage was incidental. Rebalancing the basis must BREAK the
            %build, not quietly remove the check.
            tc.verifyNotEqual(b.niX, b.niY, ...
                ['this control needs unequal trial counts -- at niX == niY the harmonic, ' ...
                 'geometric and arithmetic means coincide and it can distinguish nothing']);

            bp = tc.diffCovArray(b.pX, b.pY, b.cov_p);
            bt = tc.diffCovArray(b.iX, b.iY, b.cov_i);
            er_var = tc.diffLogSD(b.poiX_e, b.poiY_e);
            wp_cov = repmat(b.cov_poi_e, tc.NDRAW, 1);

            got = psyrat_diffrel('bp',bp,'bt',bt,'er_var',er_var,'wp_cov',wp_cov, ...
                'obs',[b.niX b.niY],'est','dep','CI',0.95);

            %Rebuild Table 6's one-facet D with a substituted covariance divisor. Written
            %out here rather than parameterized into the oracle, so the oracle keeps exactly
            %one divisor and this test owns the wrong ones.
            uni = b.pX + b.pY - 2*b.cov_p;
            dWith = @(divCov) uni ./ (uni + ...
                b.poiX_e/b.niX + b.poiY_e/b.niY - 2*b.cov_poi_e/divCov + ...
                b.iX/b.niX     + b.iY/b.niY     - 2*b.cov_i/divCov);

            %The gap must exceed the tolerance the sibling equality tests run at, not merely
            %be nonzero. A bare inequality would certify "these divisors are
            %distinguishable" even when they differ by 1e-16 -- and if the gap were below
            %TOL, the equality assertions in the two tests above would themselves pass under
            %a WRONG divisor while this control still reported the anchor as discriminating.
            %That is false assurance, so require a difference this suite can actually resolve.
            wrong = {'n''_iX alone', b.niX; ...
                     'n''_iY alone', b.niY; ...
                     'the arithmetic mean', (b.niX + b.niY)/2};
            for w = 1:size(wrong,1)
                tc.verifyGreaterThan(abs(got.pt - dWith(wrong{w,2})), tc.TOL, sprintf( ...
                    ['a coefficient built on %s as the covariance divisor is not ' ...
                     'separated from production by more than the suite''s own tolerance, ' ...
                     'so the Table 6 tests cannot distinguish divisors and are not ' ...
                     'anchoring anything'], wrong{w,1}));
            end

            %And confirm the positive half of the control: the harmonic mean DOES match.
            tc.verifyEqual(got.pt, dWith(2/((1/b.niX) + (1/b.niY))), 'AbsTol', tc.TOL, ...
                'production does not use the harmonic mean Rocha Table 6''s Note specifies');
        end

    end

    methods (Access = private)
        function [bp, wp, ii] = oneFacetSDs(tc, b)
            %Production wants SDs and, for the one-facet path, a TOTAL within-person SD
            %(wp^2 = sigma2_pi,e + sigma2_i) plus the item SD separately, because it
            %recovers the relative residual by subtraction. Convert explicitly here so the
            %variance/SD boundary is visible at the call site.
            n  = tc.NDRAW;
            bp = repmat(sqrt(b.p), n, 1);
            wp = repmat(sqrt(b.poi_e + b.i), n, 1);
            ii = repmat(sqrt(b.i), n, 1);
        end

        function a = twoFacetArgs(tc, b)
            n = tc.NDRAW;
            a = {'bp',  repmat(sqrt(b.p),     n, 1), ...
                 'bo',  repmat(sqrt(b.o),     n, 1), ...
                 'bt',  repmat(sqrt(b.i),     n, 1), ...
                 'txp', repmat(sqrt(b.pi),    n, 1), ...
                 'oxp', repmat(sqrt(b.po),    n, 1), ...
                 'txo', repmat(sqrt(b.oi),    n, 1), ...
                 'err', repmat(sqrt(b.poi_e), n, 1)};
        end

        function A = diffCovArray(tc, varX, varY, cov)
            %Pack one difference-score component into the [draw x 2 x 2] layout the
            %difference functions parse: X's variance on (1,1), Y's on (2,2), the covariance
            %on (2,1). Filled symmetrically even though production reads only the lower
            %triangle, because an asymmetric array would be a covariance matrix in name only
            %and would quietly hide a future change of which triangle is read.
            %
            %These stay on the VARIANCE scale. Unlike the one-facet and two-facet univariate
            %paths, psyrat_diffrel does not square these arguments.
            n = tc.NDRAW;
            A = zeros(n, 2, 2);
            A(:,1,1) = varX;
            A(:,2,2) = varY;
            A(:,2,1) = cov;
            A(:,1,2) = cov;
        end

        function E = diffLogSD(tc, varX, varY)
            %The residual argument is the odd one out: psyrat_diffrel and
            %psyrat_diffrel_trt both compute exp(er_var).^2, so er_var must be the LOG of a
            %standard deviation. Converting variance -> SD -> log here, in one place, is the
            %whole reason this helper exists: written inline at four call sites, a missing
            %sqrt would give a residual that is off by a square and still look plausible.
            E = [repmat(log(sqrt(varX)), tc.NDRAW, 1), ...
                 repmat(log(sqrt(varY)), tc.NDRAW, 1)];
        end

        function [c, n] = diffTwoFacetOracleArgs(~, b)
            %Map a basis onto the oracle's Table 6 struct fields. The oracle takes structs
            %rather than 25 positional arguments precisely so this mapping is checkable by
            %eye, one paper symbol per line.
            c = struct( ...
                'sigma2_pX',     b.pX,     'sigma2_pY',     b.pY,     'sigmaXY_p',     b.cov_p, ...
                'sigma2_poX',    b.poX,    'sigma2_poY',    b.poY,    'sigmaXY_po',    b.cov_po, ...
                'sigma2_piX',    b.piX,    'sigma2_piY',    b.piY,    'sigmaXY_pi',    b.cov_pi, ...
                'sigma2_poiX_e', b.poiX_e, 'sigma2_poiY_e', b.poiY_e, 'sigmaXY_poi_e', b.cov_poi_e, ...
                'sigma2_oX',     b.oX,     'sigma2_oY',     b.oY,     'sigmaXY_o',     b.cov_o, ...
                'sigma2_iX',     b.iX,     'sigma2_iY',     b.iY,     'sigmaXY_i',     b.cov_i, ...
                'sigma2_oiX',    b.oiX,    'sigma2_oiY',    b.oiY,    'sigmaXY_oi',    b.cov_oi);
            n = struct('iX', b.niX, 'iY', b.niY, 'oX', b.noX, 'oY', b.noY);
        end

        function a = diffTwoFacetProductionArgs(tc, b)
            %The same basis in psyrat_diffrel_trt's argument shapes. Six [draw x 2 x 2]
            %arrays, plus the log-SD residual and its plain-scale covariance.
            %
            %Note which production name carries which paper component: 'bt' is the TRIAL
            %main effect (not the person x trial interaction, which is 'bpi'), and the
            %three-way residual arrives split across 'er_var' and 'er_cov' rather than as a
            %seventh array.
            a = {'bp',  tc.diffCovArray(b.pX,  b.pY,  b.cov_p), ...
                 'bpi', tc.diffCovArray(b.piX, b.piY, b.cov_pi), ...
                 'bpo', tc.diffCovArray(b.poX, b.poY, b.cov_po), ...
                 'bt',  tc.diffCovArray(b.iX,  b.iY,  b.cov_i), ...
                 'bo',  tc.diffCovArray(b.oX,  b.oY,  b.cov_o), ...
                 'boi', tc.diffCovArray(b.oiX, b.oiY, b.cov_oi), ...
                 'er_var', tc.diffLogSD(b.poiX_e, b.poiY_e), ...
                 'er_cov', repmat(b.cov_poi_e, tc.NDRAW, 1), ...
                 'obs',  [b.niX b.niY], ...
                 'nocc', [b.noX b.noY]};
        end
    end
end
