classdef TestCalculationFunctions < PsyRATTestBase

    methods (Test)
        function testDepSingleObservationMatchesManualFormula(testCase)
            bp = [0.4; 0.6; 0.8; 1.0];
            wp = [0.9; 1.0; 1.1; 1.2];
            obs = 5;
            ci = 0.95;

            [ll, pt, ul] = psyrat_dep('bp', bp, 'wp', wp, 'obs', obs, 'CI', ci);

            dep = bp.^2 ./ (bp.^2 + (wp.^2 ./ obs));
            ciedge = (1 - ci) / 2;

            testCase.verifyEqual(ll, quantile(dep, ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(pt, mean(dep), 'AbsTol', 1e-12);
            testCase.verifyEqual(ul, quantile(dep, 1 - ciedge), 'AbsTol', 1e-12);
        end

        function testDepRangeObservationProducesMonotonicPointEstimate(testCase)
            bp = [0.5; 0.55; 0.60; 0.65];
            wp = [1.0; 1.1; 0.9; 1.2];
            obs = [2 6];

            [~, pt, ~] = psyrat_dep('bp', bp, 'wp', wp, 'obs', obs, 'CI', 0.95);

            % pt aligns to obs(1):obs(2) with one column per trial count and
            % no leading zeros when obs(1) > 1 (B8 fix). Dependability rises
            % monotonically with the number of trials.
            testCase.verifyEqual(numel(pt), obs(2) - obs(1) + 1);
            testCase.verifyGreaterThanOrEqual(diff(pt), 0);
        end

        function testDepRangeOutputAlignsToRequestedTrialCounts(testCase)
            % Regression (B8): the range form previously preallocated
            % zeros(0,N) and indexed by the absolute trial number, leaving
            % leading zeros when obs(1) > 1. Outputs must now align to
            % obs(1):obs(2), and remain identical to the scalar form when
            % obs(1) == 1.
            bp = [0.7; 0.8; 0.9];
            wp = [1.0; 1.1; 1.2];

            % obs(1) == 1: one column per trial count, value at column k
            % equals the scalar dependability at k trials.
            [~, ptFull, ~] = psyrat_dep('bp', bp, 'wp', wp, 'obs', [1 10], 'CI', 0.95);
            testCase.verifySize(ptFull, [1 10]);
            for k = [1 4 10]
                depK = bp.^2 ./ (bp.^2 + wp.^2 ./ k);
                testCase.verifyEqual(ptFull(k), mean(depK), 'AbsTol', 1e-12);
            end

            % obs(1) > 1: length equals the range, column idx maps to
            % trials obs(1)+idx-1, with no leading zeros.
            [~, ptOffset, ~] = psyrat_dep('bp', bp, 'wp', wp, 'obs', [5 8], 'CI', 0.95);
            testCase.verifySize(ptOffset, [1 4]);
            for idx = 1:4
                trials = 4 + idx; % 5,6,7,8
                depK = bp.^2 ./ (bp.^2 + wp.^2 ./ trials);
                testCase.verifyEqual(ptOffset(idx), mean(depK), 'AbsTol', 1e-12);
            end
        end

        function testDepWarningsForSmallCI(testCase)
            testCase.verifyWarning(@() localCallDepWithSmallCI(), 'ci:width');
        end

        function testDepInputValidation(testCase)
            testCase.verifyError(@() psyrat_dep('bp', 1, 'wp'), 'varargin:incomplete');
            testCase.verifyError(@() psyrat_dep('bp', 1), 'varargin:wp');
            testCase.verifyError(@() psyrat_dep('wp', 1, 'obs', 1, 'CI', 0.95), 'varargin:bp');
            testCase.verifyError(@() psyrat_dep('bp', 1, 'obs', 1, 'CI', 0.95), 'varargin:wp');
            testCase.verifyError(@() psyrat_dep('bp', 1, 'wp', 1, 'CI', 0.95), 'varargin:obs');
            testCase.verifyError(@() psyrat_dep('bp', 1, 'wp', 1, 'obs', 1, 'CI', 1.2), 'varargin:ci');
        end

        function testICCMatchesManualFormula(testCase)
            bp = [0.7; 0.8; 0.9; 1.0];
            wp = [0.5; 0.6; 0.7; 0.8];
            ci = 0.95;

            [ll, pt, ul] = psyrat_icc('bp', bp, 'wp', wp, 'CI', ci);

            icc = bp.^2 ./ (bp.^2 + wp.^2);
            ciedge = (1 - ci) / 2;

            testCase.verifyEqual(ll, quantile(icc, ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(pt, mean(icc), 'AbsTol', 1e-12);
            testCase.verifyEqual(ul, quantile(icc, 1 - ciedge), 'AbsTol', 1e-12);
        end

        function testICCInputValidation(testCase)
            testCase.verifyError(@() psyrat_icc('bp', 1, 'CI', 0.95), 'varargin:wp');
            testCase.verifyError(@() psyrat_icc('bp', 1, 'wp', 1, 'CI', -0.1), 'varargin:ci');
            testCase.verifyWarning(@() localCallICCWithSmallCI(), 'ci:width');
        end
        function testCoreOneFacetPointEstimatesStayWithinUnitInterval(testCase)
            bp = [0.4; 0.6; 0.8; 1.0];
            wp = [0.9; 1.0; 1.1; 1.2];

            [~, depPt, ~] = psyrat_dep('bp', bp, 'wp', wp, 'obs', 4, 'CI', 0.95);
            [~, iccPt, ~] = psyrat_icc('bp', bp, 'wp', wp, 'CI', 0.95);

            testCase.verifyGreaterThanOrEqual(depPt, 0);
            testCase.verifyLessThanOrEqual(depPt, 1);
            testCase.verifyGreaterThanOrEqual(iccPt, 0);
            testCase.verifyLessThanOrEqual(iccPt, 1);
        end
        
        
        function testRelSingMatchesAbsoluteFormulas(testCase)
            bp = [0.45; 0.55; 0.65; 0.75];
            wp = [0.90; 1.00; 1.10; 1.20];
            obs = 8;
            ci = 0.95;
            
            [ll1, pt1, ul1] = psyrat_rel_sing( ...
                'gcoeff', 1, 'metric', 'global', ...
                'bp', bp, 'wp', wp, 'obs', obs, 'CI', ci);
            [ll2, pt2, ul2] = psyrat_dep('bp', bp, 'wp', wp, 'obs', obs, 'CI', ci);
            
            testCase.verifyEqual(ll1, ll2, 'AbsTol', 1e-12);
            testCase.verifyEqual(pt1, pt2, 'AbsTol', 1e-12);
            testCase.verifyEqual(ul1, ul2, 'AbsTol', 1e-12);
            
            [ll3, pt3, ul3] = psyrat_rel_sing( ...
                'gcoeff', 1, 'metric', 'icc', ...
                'bp', bp, 'wp', wp, 'CI', ci);
            [ll4, pt4, ul4] = psyrat_icc('bp', bp, 'wp', wp, 'CI', ci);
            
            testCase.verifyEqual(ll3, ll4, 'AbsTol', 1e-12);
            testCase.verifyEqual(pt3, pt4, 'AbsTol', 1e-12);
            testCase.verifyEqual(ul3, ul4, 'AbsTol', 1e-12);
        end
        
        function testRelSingRelativeUsesItemVariance(testCase)
            bp = [0.45; 0.55; 0.65; 0.75];
            wp = [0.90; 1.00; 1.10; 1.20];
            i = [0.20; 0.25; 0.30; 0.35];
            obs = 8;
            ci = 0.95;
            ciedge = (1 - ci) / 2;
            
            bpv = bp.^2;
            wpv = wp.^2;
            iv = i.^2;
            piv = wpv - iv;
            
            relGlobal = bpv ./ (bpv + (piv ./ obs));
            [llg, ptg, ulg] = psyrat_rel_sing( ...
                'gcoeff', 2, 'metric', 'global', ...
                'bp', bp, 'wp', wp, 'i', i, 'obs', obs, 'CI', ci);
            
            testCase.verifyEqual(llg, quantile(relGlobal, ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(ptg, mean(relGlobal), 'AbsTol', 1e-12);
            testCase.verifyEqual(ulg, quantile(relGlobal, 1 - ciedge), 'AbsTol', 1e-12);
            
            relICC = bpv ./ (bpv + piv);
            [lli, pti, uli] = psyrat_rel_sing( ...
                'gcoeff', 2, 'metric', 'icc', ...
                'bp', bp, 'wp', wp, 'i', i, 'CI', ci);
            
            testCase.verifyEqual(lli, quantile(relICC, ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(pti, mean(relICC), 'AbsTol', 1e-12);
            testCase.verifyEqual(uli, quantile(relICC, 1 - ciedge), 'AbsTol', 1e-12);
        end
        
        function testRelSingValidation(testCase)
            testCase.verifyError(@() psyrat_rel_sing('gcoeff',1,'metric'), 'varargin:incomplete');
            testCase.verifyError(@() psyrat_rel_sing('gcoeff',1,'metric','bad','bp',1,'wp',1,'CI',0.95), ...
                'varargin:metric');
            testCase.verifyError(@() psyrat_rel_sing('gcoeff',1,'metric','global','bp',1,'wp',1,'CI',0.95), ...
                'varargin:obs');
        end

        function testRelSingOneFacetTrialSplit(testCase)
            %One-facet Rocha Table 2 split (S10). The relsummary wiring passes
            %wp = sqrt(sig_e^2 + sig_trl^2) (total within) and i = sig_trl
            %(trial main effect). Then the generalizability coefficient
            %(gcoeff=2) must EXCLUDE the trial main effect (error = sig_e^2/n)
            %while dependability (gcoeff=1) INCLUDES it (error =
            %(sig_e^2+sig_trl^2)/n). Consequently G >= D, with equality when
            %sig_trl = 0.
            bp      = [0.45; 0.55; 0.65; 0.75];
            sig_e   = [0.90; 1.00; 1.10; 1.20];
            sig_trl = [0.40; 0.50; 0.30; 0.35];
            wp      = sqrt(sig_e.^2 + sig_trl.^2);
            obs = 8; ci = 0.95;

            %Generalizability equals dependability computed on the residual-only
            %within SD (the trial main effect is removed from relative error).
            [llG, ptG, ulG] = psyrat_rel_sing('gcoeff',2,'metric','global', ...
                'bp',bp,'wp',wp,'i',sig_trl,'obs',obs,'CI',ci);
            [llGref, ptGref, ulGref] = PsyRATAccuracyOracle.dep(bp,sig_e,obs,ci);
            testCase.verifyEqual(llG, llGref, 'AbsTol', 1e-12);
            testCase.verifyEqual(ptG, ptGref, 'AbsTol', 1e-12);
            testCase.verifyEqual(ulG, ulGref, 'AbsTol', 1e-12);

            %Dependability equals dependability computed on the total within SD.
            [~, ptD] = psyrat_rel_sing('gcoeff',1,'metric','global', ...
                'bp',bp,'wp',wp,'i',sig_trl,'obs',obs,'CI',ci);
            [~, ptDref] = PsyRATAccuracyOracle.dep(bp,wp,obs,ci);
            testCase.verifyEqual(ptD, ptDref, 'AbsTol', 1e-12);

            %The trial main effect inflates the absolute error only, so G >= D.
            testCase.verifyGreaterThanOrEqual(ptG, ptD);

            %With no trial main effect (sig_trl = 0) the two coincide.
            zero_trl = zeros(size(sig_e));
            [~, ptG0] = psyrat_rel_sing('gcoeff',2,'metric','global', ...
                'bp',bp,'wp',sig_e,'i',zero_trl,'obs',obs,'CI',ci);
            [~, ptD0] = psyrat_rel_sing('gcoeff',1,'metric','global', ...
                'bp',bp,'wp',sig_e,'i',zero_trl,'obs',obs,'CI',ci);
            testCase.verifyEqual(ptG0, ptD0, 'AbsTol', 1e-12);
        end

        function testCutScoreSingMatchesTable2Formulas(testCase)
            % Cut-scores are absolute-error only (Rocha 2026, Table 2): the
            % only valid one-facet cut-score is the dependability (gcoeff = 1)
            % form, whose error keeps the trial main effect via wp^2.
            bp = [0.40; 0.60; 0.80; 1.00];
            wp = [0.90; 1.00; 1.10; 1.20];
            mu = [1.10; 1.00; 0.90; 1.20];
            cut = 1.00;
            obs = 5;
            ci = 0.95;
            ciedge = (1 - ci) / 2;

            bpv = bp.^2;
            wpv = wp.^2;
            offsq = (mu - cut).^2;

            depRel = (bpv + offsq) ./ (bpv + offsq + (wpv ./ obs));

            [lld, ptd, uld] = psyrat_cutscore_sing( ...
                'gcoeff', 1, 'bp', bp, 'wp', wp, 'mu', mu, ...
                'cut', cut, 'obs', obs, 'CI', ci);

            testCase.verifyEqual(lld, quantile(depRel, ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(ptd, mean(depRel), 'AbsTol', 1e-12);
            testCase.verifyEqual(uld, quantile(depRel, 1 - ciedge), 'AbsTol', 1e-12);

            % A relative-error cut-score is not a valid quantity: gcoeff = 2
            % is rejected.
            testCase.verifyError(@() psyrat_cutscore_sing( ...
                'gcoeff', 2, 'bp', bp, 'wp', wp, 'mu', mu, ...
                'cut', cut, 'obs', obs, 'CI', ci), ...
                'varargin:gcoeff');
        end

        function testCutScoreSingValidation(testCase)
            testCase.verifyError(@() psyrat_cutscore_sing('gcoeff',1,'bp',1), 'varargin:wp');
            testCase.verifyError(@() psyrat_cutscore_sing( ...
                'gcoeff',1,'bp',1,'wp',1,'mu',0,'cut',0,'obs',2,'CI',1.1), ...
                'varargin:ci');
            testCase.verifyError(@() psyrat_cutscore_sing( ...
                'gcoeff',7,'bp',1,'wp',1,'mu',0,'cut',0,'obs',2,'CI',0.95), ...
                'varargin:gcoeff');
            % Cut-scores are absolute-only: the relative decision (gcoeff = 2)
            % must be rejected, not silently computed.
            testCase.verifyError(@() psyrat_cutscore_sing( ...
                'gcoeff',2,'bp',1,'wp',1,'mu',0,'cut',0,'obs',2,'CI',0.95), ...
                'varargin:gcoeff');
        end
        function testRelTrtMatchesManualFormulaForDepEquivalence(testCase)
            bp = [0.50; 0.60; 0.55; 0.65];
            bo = [0.10; 0.12; 0.11; 0.10];
            bt = [0.20; 0.22; 0.24; 0.21];
            txp = [0.18; 0.17; 0.16; 0.19];
            oxp = [0.13; 0.14; 0.15; 0.14];
            txo = [0.09; 0.10; 0.11; 0.10];
            err = [0.70; 0.75; 0.72; 0.71];
            obs = 6;
            ci = 0.95;

            [ll, pt, ul] = psyrat_rel_trt( ...
                'gcoeff', 1, 'reltype', 1, ...
                'bp', bp, 'bo', bo, 'bt', bt, ...
                'txp', txp, 'oxp', oxp, 'txo', txo, 'err', err, ...
                'obs', obs, 'CI', ci);

            bp2 = bp.^2;
            bt2 = bt.^2;
            txp2 = txp.^2;
            oxp2 = oxp.^2;
            txo2 = txo.^2;
            err2 = err.^2;

            uni = bp2 + oxp2;
            absErr = (txp2 ./ obs) + (err2 ./ obs) + (bt2 ./ obs) + (txo2 ./ obs);
            rel = uni ./ (uni + absErr);

            ciedge = (1 - ci) / 2;
            testCase.verifyEqual(ll, quantile(rel, ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(pt, mean(rel), 'AbsTol', 1e-12);
            testCase.verifyEqual(ul, quantile(rel, 1 - ciedge), 'AbsTol', 1e-12);
        end

        function testRelTrtRangeObservationOutputShape(testCase)
            args = { ...
                'gcoeff', 2, 'reltype', 2, ...
                'bp', [0.5; 0.6; 0.7], ...
                'bo', [0.1; 0.1; 0.2], ...
                'bt', [0.2; 0.2; 0.3], ...
                'txp', [0.1; 0.2; 0.3], ...
                'oxp', [0.1; 0.2; 0.1], ...
                'txo', [0.05; 0.05; 0.08], ...
                'err', [0.7; 0.8; 0.9], ...
                'CI', 0.95};

            % Range output aligns to obs(1):obs(2): one column per requested
            % trial count with no leading zeros when obs(1) > 1 (B7). The
            % obs(1) > 1 result must equal the matching slice of the obs(1) == 1
            % result (inputs are fixed draws, so the values are identical).
            [~, ptRange, ~] = psyrat_rel_trt(args{:}, 'obs', [2 5]);
            testCase.verifyEqual(numel(ptRange), 4);
            testCase.verifyGreaterThanOrEqual(ptRange, zeros(1, 4));
            testCase.verifyLessThanOrEqual(ptRange, ones(1, 4));

            [~, ptFull, ~] = psyrat_rel_trt(args{:}, 'obs', [1 5]);
            testCase.verifyEqual(numel(ptFull), 5);
            testCase.verifyEqual(ptRange, ptFull(2:5), 'AbsTol', 1e-12);
        end
        
        function testRelTrtMatchesManualFormulaForDepCES(testCase)
            bp = [0.50; 0.60; 0.55; 0.65];
            bo = [0.10; 0.12; 0.11; 0.10];
            bt = [0.20; 0.22; 0.24; 0.21];
            txp = [0.18; 0.17; 0.16; 0.19];
            oxp = [0.13; 0.14; 0.15; 0.14];
            txo = [0.09; 0.10; 0.11; 0.10];
            err = [0.70; 0.75; 0.72; 0.71];
            obs = 6;
            nocc = 2;
            ci = 0.95;
            
            [ll, pt, ul] = psyrat_rel_trt( ...
                'gcoeff', 1, 'reltype', 3, ...
                'bp', bp, 'bo', bo, 'bt', bt, ...
                'txp', txp, 'oxp', oxp, 'txo', txo, 'err', err, ...
                'obs', obs, 'nocc', nocc, 'CI', ci);
            
            bp2 = bp.^2;
            bo2 = bo.^2;
            bt2 = bt.^2;
            txp2 = txp.^2;
            oxp2 = oxp.^2;
            txo2 = txo.^2;
            err2 = err.^2;
            
            uni = bp2;
            absErr = (txp2 ./ obs) + (oxp2 ./ nocc) + (err2 ./ (obs * nocc)) + ...
                (bt2 ./ obs) + (bo2 ./ nocc) + (txo2 ./ (obs * nocc));
            rel = uni ./ (uni + absErr);
            
            ciedge = (1 - ci) / 2;
            testCase.verifyEqual(ll, quantile(rel, ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(pt, mean(rel), 'AbsTol', 1e-12);
            testCase.verifyEqual(ul, quantile(rel, 1 - ciedge), 'AbsTol', 1e-12);
        end
        
        function testRelTrtMatchesManualFormulaForGenCES(testCase)
            bp = [0.50; 0.60; 0.55; 0.65];
            bo = [0.10; 0.12; 0.11; 0.10];
            bt = [0.20; 0.22; 0.24; 0.21];
            txp = [0.18; 0.17; 0.16; 0.19];
            oxp = [0.13; 0.14; 0.15; 0.14];
            txo = [0.09; 0.10; 0.11; 0.10];
            err = [0.70; 0.75; 0.72; 0.71];
            obs = 6;
            nocc = 2;
            ci = 0.95;
            
            [ll, pt, ul] = psyrat_rel_trt( ...
                'gcoeff', 2, 'reltype', 3, ...
                'bp', bp, 'bo', bo, 'bt', bt, ...
                'txp', txp, 'oxp', oxp, 'txo', txo, 'err', err, ...
                'obs', obs, 'nocc', nocc, 'CI', ci);
            
            bp2 = bp.^2;
            txp2 = txp.^2;
            oxp2 = oxp.^2;
            err2 = err.^2;
            
            relErr = (txp2 ./ obs) + (oxp2 ./ nocc) + (err2 ./ (obs * nocc));
            rel = bp2 ./ (bp2 + relErr);
            
            ciedge = (1 - ci) / 2;
            testCase.verifyEqual(ll, quantile(rel, ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(pt, mean(rel), 'AbsTol', 1e-12);
            testCase.verifyEqual(ul, quantile(rel, 1 - ciedge), 'AbsTol', 1e-12);
        end
        
        function testRelTrtInputValidation(testCase)
            testCase.verifyError(@() psyrat_rel_trt('gcoeff', 1, 'reltype'), 'varargin:incomplete');
            testCase.verifyError(@() psyrat_rel_trt('gcoeff', 1), 'varargin:reltype');
            testCase.verifyError(@() psyrat_rel_trt( ...
                'gcoeff', 1, 'bp', 1, 'bo', 1, 'bt', 1, ...
                'txp', 1, 'oxp', 1, 'txo', 1, 'err', 1, 'obs', 2, 'CI', 0.95), ...
                'varargin:reltype');
            testCase.verifyError(@() psyrat_rel_trt( ...
                'gcoeff', 1, 'reltype', 1, 'bp', 1, 'bo', 1, 'bt', 1, ...
                'txp', 1, 'oxp', 1, 'txo', 1, 'err', 1, 'obs', 2, 'CI', 1.5), ...
                'varargin:ci');
        end
        
        function testCutScoreDepMatchesManualFormula(testCase)
            bp = [0.4; 0.6; 0.8; 1.0];
            wp = [0.9; 1.0; 1.1; 1.2];
            i = [0.2; 0.2; 0.25; 0.3];
            mu = [1.1; 1.0; 0.9; 1.2];
            cut = 1.0;
            obs = 5;
            ci = 0.95;
            
            [ll, pt, ul] = psyrat_cutscore_dep( ...
                'bp', bp, 'wp', wp, 'i', i, ...
                'mu', mu, 'cut', cut, 'obs', obs, 'CI', ci, 'est', 'dep');
            
            bp2 = bp.^2;
            wp2 = wp.^2;
            i2 = i.^2;
            offsq = (mu - cut).^2;
            rel = (bp2 + offsq) ./ (bp2 + offsq + (wp2 ./ obs) + (i2 ./ obs));
            ciedge = (1 - ci) / 2;
            
            testCase.verifyEqual(ll, quantile(rel, ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(pt, mean(rel), 'AbsTol', 1e-12);
            testCase.verifyEqual(ul, quantile(rel, 1 - ciedge), 'AbsTol', 1e-12);
        end
        
        function testCutScoreDepValidation(testCase)
            % Cut-scores are absolute-error only, so the relative/gen
            % cut-score case is no longer exercised (a relative-error cut-score
            % is not a valid quantity). Only input validation is checked here.
            bp = [0.4; 0.6; 0.8; 1.0];
            wp = [0.9; 1.0; 1.1; 1.2];

            testCase.verifyError(@() psyrat_cutscore_dep('bp', bp, 'wp', wp, 'obs', 1, 'CI', .95), ...
                'varargin:mu');
        end
        
        function testCutScoreTrtMatchesManualFormulaForDepCES(testCase)
            bp = [0.50; 0.60; 0.55; 0.65];
            bo = [0.10; 0.12; 0.11; 0.10];
            bt = [0.20; 0.22; 0.24; 0.21];
            txp = [0.18; 0.17; 0.16; 0.19];
            oxp = [0.13; 0.14; 0.15; 0.14];
            txo = [0.09; 0.10; 0.11; 0.10];
            err = [0.70; 0.75; 0.72; 0.71];
            mu = [1.2; 1.0; 0.9; 1.1];
            cut = 1.0;
            obs = 6;
            nocc = 2;
            ci = 0.95;
            
            [ll, pt, ul] = psyrat_cutscore_trt( ...
                'gcoeff', 1, 'reltype', 3, ...
                'bp', bp, 'bo', bo, 'bt', bt, ...
                'txp', txp, 'oxp', oxp, 'txo', txo, 'err', err, ...
                'mu', mu, 'cut', cut, 'obs', obs, 'nocc', nocc, 'CI', ci);
            
            bp2 = bp.^2;
            bo2 = bo.^2;
            bt2 = bt.^2;
            txp2 = txp.^2;
            oxp2 = oxp.^2;
            txo2 = txo.^2;
            err2 = err.^2;
            offsq = (mu - cut).^2;
            
            absErr = (txp2 ./ obs) + (oxp2 ./ nocc) + (err2 ./ (obs * nocc)) + ...
                (bt2 ./ obs) + (bo2 ./ nocc) + (txo2 ./ (obs * nocc));
            rel = (bp2 + offsq) ./ (bp2 + offsq + absErr);
            
            ciedge = (1 - ci) / 2;
            testCase.verifyEqual(ll, quantile(rel, ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(pt, mean(rel), 'AbsTol', 1e-12);
            testCase.verifyEqual(ul, quantile(rel, 1 - ciedge), 'AbsTol', 1e-12);
            
            testCase.verifyError(@() psyrat_cutscore_trt( ...
                'gcoeff', 1, 'reltype', 4, ...
                'bp', bp, 'bo', bo, 'bt', bt, ...
                'txp', txp, 'oxp', oxp, 'txo', txo, 'err', err, ...
                'mu', mu, 'cut', cut, 'obs', obs, 'CI', ci), ...
                'varargin:reltype');

            % Cut-scores are absolute-only: a relative-error (gcoeff = 2)
            % cut-score is not a valid quantity and must be rejected.
            testCase.verifyError(@() psyrat_cutscore_trt( ...
                'gcoeff', 2, 'reltype', 3, ...
                'bp', bp, 'bo', bo, 'bt', bt, ...
                'txp', txp, 'oxp', oxp, 'txo', txo, 'err', err, ...
                'mu', mu, 'cut', cut, 'obs', obs, 'nocc', nocc, 'CI', ci), ...
                'varargin:gcoeff');
        end

        function testCutScoreRangeOutputsAlignToObsRange(testCase)
            % B7: psyrat_cutscore_{sing,dep,trt} range outputs align to
            % obs(1):obs(2) (one column per trial count, no leading zeros when
            % obs(1) > 1), matching psyrat_dep/psyrat_rel_sing. Inputs are fixed
            % draws, so the obs(1) > 1 result must equal the matching slice of
            % the obs(1) == 1 result.
            bp = [0.50; 0.60; 0.70];
            wp = [0.90; 1.00; 1.10];
            i = [0.20; 0.20; 0.25];
            mu = [1.10; 1.00; 0.90];
            cut = 1.0; ci = 0.95;

            [~, ptR, ~] = psyrat_cutscore_sing('gcoeff', 1, 'bp', bp, 'wp', wp, ...
                'mu', mu, 'cut', cut, 'obs', [2 5], 'CI', ci);
            [~, ptF, ~] = psyrat_cutscore_sing('gcoeff', 1, 'bp', bp, 'wp', wp, ...
                'mu', mu, 'cut', cut, 'obs', [1 5], 'CI', ci);
            testCase.verifyEqual(numel(ptR), 4);
            testCase.verifyEqual(numel(ptF), 5);
            testCase.verifyEqual(ptR, ptF(2:5), 'AbsTol', 1e-12);

            [~, ptR, ~] = psyrat_cutscore_dep('bp', bp, 'wp', wp, 'i', i, ...
                'mu', mu, 'cut', cut, 'obs', [2 5], 'CI', ci, 'est', 'dep');
            [~, ptF, ~] = psyrat_cutscore_dep('bp', bp, 'wp', wp, 'i', i, ...
                'mu', mu, 'cut', cut, 'obs', [1 5], 'CI', ci, 'est', 'dep');
            testCase.verifyEqual(numel(ptR), 4);
            testCase.verifyEqual(ptR, ptF(2:5), 'AbsTol', 1e-12);

            bo = [0.10; 0.12; 0.11]; bt = [0.20; 0.22; 0.24];
            txp = [0.18; 0.17; 0.16]; oxp = [0.13; 0.14; 0.15];
            txo = [0.09; 0.10; 0.11]; err = [0.70; 0.75; 0.72];
            trtArgs = {'gcoeff', 1, 'reltype', 3, 'bp', bp, 'bo', bo, 'bt', bt, ...
                'txp', txp, 'oxp', oxp, 'txo', txo, 'err', err, 'mu', mu, ...
                'cut', cut, 'nocc', 2, 'CI', ci};
            [~, ptR, ~] = psyrat_cutscore_trt(trtArgs{:}, 'obs', [2 5]);
            [~, ptF, ~] = psyrat_cutscore_trt(trtArgs{:}, 'obs', [1 5]);
            testCase.verifyEqual(numel(ptR), 4);
            testCase.verifyEqual(ptR, ptF(2:5), 'AbsTol', 1e-12);
        end

        function testSubjectSpecificCutScoreReturnsAugmentedTable(testCase)
            bp = { [0.50; 0.60; 0.55] };
            wp_pop = log([0.70; 0.75; 0.80]);
            wp_ss = { [0.02 0.03 0.01; 0.01 0.02 0.03; 0.03 0.01 0.02] };
            mu = [1.00 1.20 1.10; 1.10 1.00 1.20; 1.05 1.15 1.05];
            idtable = table({'s1'; 's2'; 's3'}, [1; 2; 3], [5; 6; 7], ...
                'VariableNames', {'id', 'id2', 'trls'});
            cut = 1.0;
            ci = 0.95;
            
            out = psyrat_sscutscore( ...
                'bp', bp, ...
                'wp_pop', wp_pop, ...
                'wp_ss', wp_ss, ...
                'mu', mu, ...
                'cut', cut, ...
                'idtable', idtable, ...
                'CI', ci, ...
                'est', 'dep');
            
            newVars = {'cutdep_pt', 'cutdep_ll', 'cutdep_ul'};
            for i = 1:numel(newVars)
                testCase.verifyTrue(any(strcmp(out.Properties.VariableNames, newVars{i})));
            end
            
            bpVar = cell2mat(bp).^2;
            wpVar = exp(wp_pop + cell2mat(wp_ss)).^2;
            offsq = (mu - cut).^2;
            rel = (bpVar + offsq) ./ (bpVar + offsq + (wpVar ./ idtable.trls(:)'));
            ciedge = (1 - ci) / 2;
            
            testCase.verifyEqual(out.cutdep_ll', quantile(rel, ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.cutdep_pt', mean(rel), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.cutdep_ul', quantile(rel, 1 - ciedge), 'AbsTol', 1e-12);
            testCase.verifyError(@() psyrat_sscutscore('bp', bp, 'wp_pop', wp_pop, 'wp_ss', wp_ss, ...
                'mu', mu, 'cut', cut, 'CI', ci), ...
                'varargin:idtable');
        end

        function testSubjectSpecificCutScoreTrialMainEffectLowersDependability(testCase)
            % F6a: psyrat_sscutscore is retained (not called from
            % subroutines/) on the strength of its numeric pinning, so the
            % sigma_i^2 (trial main effect) half of its absolute error
            % err = sigma_pi,e^2(s)/n' + sigma_i^2/n' must be exercised with
            % a nonzero sigma_i, not only the i = 0 default the sibling lane
            % uses. sigma_i draws are SDs (squared internally); the input is
            % accepted when it carries one element or one per posterior draw
            % (validated by element count after i(:) flattening).
            bp = { [0.50; 0.60; 0.55] };
            wp_pop = log([0.70; 0.75; 0.80]);
            wp_ss = { [0.02 0.03 0.01; 0.01 0.02 0.03; 0.03 0.01 0.02] };
            mu = [1.00 1.20 1.10; 1.10 1.00 1.20; 1.05 1.15 1.05];
            idtable = table({'s1'; 's2'; 's3'}, [1; 2; 3], [5; 6; 7], ...
                'VariableNames', {'id', 'id2', 'trls'});
            cut = 1.0;
            ci = 0.95;
            sigI = [0.40; 0.35; 0.45]; % sigma_i draws (one per posterior draw)

            outWithI = psyrat_sscutscore( ...
                'bp', bp, ...
                'wp_pop', wp_pop, ...
                'wp_ss', wp_ss, ...
                'mu', mu, ...
                'cut', cut, ...
                'idtable', idtable, ...
                'CI', ci, ...
                'est', 'dep', ...
                'i', sigI);
            outNoI = psyrat_sscutscore( ...
                'bp', bp, ...
                'wp_pop', wp_pop, ...
                'wp_ss', wp_ss, ...
                'mu', mu, ...
                'cut', cut, ...
                'idtable', idtable, ...
                'CI', ci, ...
                'est', 'dep');

            bpVar = cell2mat(bp).^2;
            wpVar = exp(wp_pop + cell2mat(wp_ss)).^2;
            offsq = (mu - cut).^2;
            errAbs = (wpVar ./ idtable.trls(:)') + (sigI.^2 ./ idtable.trls(:)');
            rel = (bpVar + offsq) ./ (bpVar + offsq + errAbs);
            ciedge = (1 - ci) / 2;

            testCase.verifyEqual(outWithI.cutdep_ll', quantile(rel, ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(outWithI.cutdep_pt', mean(rel), 'AbsTol', 1e-12);
            testCase.verifyEqual(outWithI.cutdep_ul', quantile(rel, 1 - ciedge), 'AbsTol', 1e-12);

            % Nonzero sigma_i adds sigma_i^2/n' to the absolute error, so it
            % must strictly lower the dependability point estimate for every
            % subject relative to the i = 0 call on the same inputs.
            testCase.verifyLessThan(outWithI.cutdep_pt, outNoI.cutdep_pt);
        end

        function testDiffRelReturnsExpectedFieldsForDependability(testCase)
            [bp, bt, erVar] = localDiffRelInputs();
            out = psyrat_diffrel( ...
                'bp', bp, ...
                'bt', bt, ...
                'er_var', erVar, ...
                'obs', [8 10], ...
                'est', 'dep', ...
                'CI', 0.95);

            expectedFields = { ...
                'est_type', 'pt', 'll', 'ul', ...
                'bp_cov_pt', 'bp_cov_ll', 'bp_cov_ul', ...
                'bt_cov_pt', 'bt_cov_ll', 'bt_cov_ul', ...
                'wp_cov', 'icc_pt', 'icc_ll', 'icc_ul'};
            for i = 1:numel(expectedFields)
                testCase.verifyTrue(isfield(out, expectedFields{i}));
            end
            testCase.verifyEqual(out.est_type, 'dependability');
            testCase.verifyGreaterThanOrEqual(out.pt, 0);
            testCase.verifyLessThanOrEqual(out.pt, 1);
        end

        function testDiffRelGeneralizabilityAndValidation(testCase)
            [bp, bt, erVar] = localDiffRelInputs();
            out = psyrat_diffrel( ...
                'bp', bp, ...
                'bt', bt, ...
                'er_var', erVar, ...
                'obs', [8 10], ...
                'est', 'gen', ...
                'CI', 0.95);
            testCase.verifyEqual(out.est_type, 'generalizability');
            testCase.verifyError(@() psyrat_diffrel( ...
                'bp', bp, 'bt', bt, 'er_var', erVar, 'obs', [8 10], 'est', 'gen', 'CI', 1.2), ...
                'varargin:ci');
        end
        
        function testDiffRelDefaultsWpCovToZero(testCase)
            [bp, bt, erVar] = localDiffRelInputs();
            outDefault = psyrat_diffrel( ...
                'bp', bp, ...
                'bt', bt, ...
                'er_var', erVar, ...
                'obs', [8 10], ...
                'est', 'dep', ...
                'CI', 0.95);
            outZero = psyrat_diffrel( ...
                'bp', bp, ...
                'bt', bt, ...
                'er_var', erVar, ...
                'wp_cov', zeros(size(erVar,1),1), ...
                'obs', [8 10], ...
                'est', 'dep', ...
                'CI', 0.95);
            
            testCase.verifyEqual(outDefault.pt, outZero.pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(outDefault.ll, outZero.ll, 'AbsTol', 1e-12);
            testCase.verifyEqual(outDefault.ul, outZero.ul, 'AbsTol', 1e-12);
            testCase.verifyEqual(outDefault.wp_cov, 0, 'AbsTol', 1e-12);
        end
        
        function testDiffRelUsesProvidedWpCov(testCase)
            [bp, bt, erVar] = localDiffRelInputs();
            wpCov = 0.10 * ones(size(erVar,1),1);
            outZero = psyrat_diffrel( ...
                'bp', bp, ...
                'bt', bt, ...
                'er_var', erVar, ...
                'wp_cov', zeros(size(erVar,1),1), ...
                'obs', [8 10], ...
                'est', 'dep', ...
                'CI', 0.95);
            outCov = psyrat_diffrel( ...
                'bp', bp, ...
                'bt', bt, ...
                'er_var', erVar, ...
                'wp_cov', wpCov, ...
                'obs', [8 10], ...
                'est', 'dep', ...
                'CI', 0.95);
            
            testCase.verifyNotEqual(outCov.pt, outZero.pt);
            testCase.verifyEqual(outCov.wp_cov, mean(wpCov), 'AbsTol', 1e-12);
        end

        function testDiffRelMatchesTable6OneFacetDependabilityWithResidualCov(testCase)
            bp = localCovmat(...
                [0.50; 0.55; 0.60; 0.65],...
                [0.45; 0.50; 0.55; 0.60],...
                [0.10; 0.11; 0.10; 0.12]);
            bt = localCovmat(...
                [0.25; 0.27; 0.26; 0.28],...
                [0.20; 0.22; 0.24; 0.23],...
                [0.04; 0.05; 0.04; 0.06]);
            
            sd1 = [0.75; 0.78; 0.79; 0.77];
            sd2 = [0.82; 0.80; 0.81; 0.83];
            er_var = log([sd1 sd2]);
            wp1 = exp(er_var(:,1)) .^ 2;
            wp2 = exp(er_var(:,2)) .^ 2;
            wp_cov = [0.03; 0.02; 0.03; 0.02];
            obs = [8 10];
            ci = 0.95;
            
            out = psyrat_diffrel(...
                'bp', bp, ...
                'bt', bt, ...
                'er_var', er_var, ...
                'wp_cov', wp_cov, ...
                'obs', obs, ...
                'est', 'dep', ...
                'CI', ci);
            
            bp1 = bp(:,1,1);
            bp2 = bp(:,2,2);
            bp_cov = bp(:,2,1);
            
            bt1 = bt(:,1,1);
            bt2 = bt(:,2,2);
            bt_cov = bt(:,2,1);
            
            uni = bp1 + bp2 - (2 .* bp_cov);
            rel_err = (wp1 ./ obs(1)) + (wp2 ./ obs(2)) - ...
                ((2 .* wp_cov) ./ psyrat_harmmean(obs));
            abs_err = rel_err + (bt1 ./ obs(1)) + (bt2 ./ obs(2)) - ...
                ((2 .* bt_cov) ./ psyrat_harmmean(obs));
            
            rel = uni ./ (uni + abs_err);
            
            denom_rel_err = wp1 + wp2 - (2 .* wp_cov);
            denom_abs_err = denom_rel_err + bt1 + bt2 - (2 .* bt_cov);
            icc = uni ./ (uni + denom_abs_err);
            
            ciedge = (1-ci)/2;
            testCase.verifyEqual(out.ll, quantile(rel,ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.pt, mean(rel), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.ul, quantile(rel,1-ciedge), 'AbsTol', 1e-12);
            
            testCase.verifyEqual(out.icc_ll, quantile(icc,ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.icc_pt, mean(icc), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.icc_ul, quantile(icc,1-ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.wp_cov, mean(wp_cov), 'AbsTol', 1e-12);
            testCase.verifyTrue(out.wp_cov_est);
        end
        
        function testDiffRelDefaultResidualCovarianceIsZero(testCase)
            bp = localCovmat(...
                [0.50; 0.55; 0.60; 0.65],...
                [0.45; 0.50; 0.55; 0.60],...
                [0.10; 0.11; 0.10; 0.12]);
            bt = localCovmat(...
                [0.25; 0.27; 0.26; 0.28],...
                [0.20; 0.22; 0.24; 0.23],...
                [0.04; 0.05; 0.04; 0.06]);
            er_var = log([...
                0.75 0.82; ...
                0.78 0.80; ...
                0.79 0.81; ...
                0.77 0.83]);
            
            out1 = psyrat_diffrel(...
                'bp', bp, ...
                'bt', bt, ...
                'er_var', er_var, ...
                'obs', [8 10], ...
                'est', 'gen', ...
                'CI', 0.95);
            
            out2 = psyrat_diffrel(...
                'bp', bp, ...
                'bt', bt, ...
                'er_var', er_var, ...
                'wp_cov', 0, ...
                'obs', [8 10], ...
                'est', 'gen', ...
                'CI', 0.95);
            
            testCase.verifyFalse(out1.wp_cov_est);
            testCase.verifyEqual(out1.wp_cov, 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(out1.wp_cov_ll, 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(out1.wp_cov_ul, 0, 'AbsTol', 1e-12);
            testCase.verifyTrue(out2.wp_cov_est);
            
            testCase.verifyEqual(out1.pt, out2.pt, 'AbsTol', 1e-12);
            testCase.verifyEqual(out1.ll, out2.ll, 'AbsTol', 1e-12);
            testCase.verifyEqual(out1.ul, out2.ul, 'AbsTol', 1e-12);
            testCase.verifyEqual(out1.icc_pt, out2.icc_pt, 'AbsTol', 1e-12);
        end
        
        function testDiffRelTrtMatchesTable6CESDependability(testCase)
            bp = localCovmat(...
                [0.55; 0.58; 0.60; 0.62],...
                [0.50; 0.52; 0.54; 0.56],...
                [0.14; 0.15; 0.16; 0.15]);
            bpi = localCovmat(...
                [0.25; 0.26; 0.27; 0.28],...
                [0.22; 0.23; 0.24; 0.25],...
                [0.06; 0.06; 0.07; 0.06]);
            bpo = localCovmat(...
                [0.20; 0.21; 0.22; 0.23],...
                [0.18; 0.19; 0.20; 0.21],...
                [0.05; 0.05; 0.06; 0.05]);
            bt = localCovmat(...
                [0.18; 0.19; 0.20; 0.21],...
                [0.16; 0.17; 0.18; 0.19],...
                [0.04; 0.04; 0.05; 0.04]);
            bo = localCovmat(...
                [0.10; 0.11; 0.12; 0.13],...
                [0.09; 0.10; 0.11; 0.12],...
                [0.02; 0.02; 0.03; 0.02]);
            boi = localCovmat(...
                [0.08; 0.09; 0.10; 0.11],...
                [0.07; 0.08; 0.09; 0.10],...
                [0.02; 0.02; 0.02; 0.02]);
            
            sd1 = [0.70; 0.72; 0.74; 0.76];
            sd2 = [0.68; 0.70; 0.72; 0.74];
            er_var = log([sd1 sd2]);
            er_cov = [0.03; 0.03; 0.04; 0.04];
            
            obs = [8 10];
            nocc = [2 3];
            ci = 0.95;
            
            out = psyrat_diffrel_trt(...
                'bp', bp, ...
                'bpi', bpi, ...
                'bpo', bpo, ...
                'bt', bt, ...
                'bo', bo, ...
                'boi', boi, ...
                'er_var', er_var, ...
                'er_cov', er_cov, ...
                'obs', obs, ...
                'nocc', nocc, ...
                'reltype', 3, ...
                'est', 'dep', ...
                'CI', ci);
            
            rel = localTable6TrtRel(...
                bp,bpi,bpo,bt,bo,boi,exp(er_var(:,1)).^2,exp(er_var(:,2)).^2,...
                er_cov,obs,nocc,3,'dep');
            icc = localTable6TrtRel(...
                bp,bpi,bpo,bt,bo,boi,exp(er_var(:,1)).^2,exp(er_var(:,2)).^2,...
                er_cov,[1 1],[1 1],3,'dep');
            
            ciedge = (1-ci)/2;
            testCase.verifyEqual(out.ll, quantile(rel,ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.pt, mean(rel), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.ul, quantile(rel,1-ciedge), 'AbsTol', 1e-12);
            
            testCase.verifyEqual(out.icc_ll, quantile(icc,ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.icc_pt, mean(icc), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.icc_ul, quantile(icc,1-ciedge), 'AbsTol', 1e-12);
            testCase.verifyTrue(out.wp_cov_est);
        end
        
        function testDiffRelTrtMatchesTable6CEGeneralizability(testCase)
            bp = localCovmat(...
                [0.55; 0.58; 0.60; 0.62],...
                [0.50; 0.52; 0.54; 0.56],...
                [0.14; 0.15; 0.16; 0.15]);
            bpi = localCovmat(...
                [0.25; 0.26; 0.27; 0.28],...
                [0.22; 0.23; 0.24; 0.25],...
                [0.06; 0.06; 0.07; 0.06]);
            bpo = localCovmat(...
                [0.20; 0.21; 0.22; 0.23],...
                [0.18; 0.19; 0.20; 0.21],...
                [0.05; 0.05; 0.06; 0.05]);
            bt = localCovmat(...
                [0.18; 0.19; 0.20; 0.21],...
                [0.16; 0.17; 0.18; 0.19],...
                [0.04; 0.04; 0.05; 0.04]);
            bo = localCovmat(...
                [0.10; 0.11; 0.12; 0.13],...
                [0.09; 0.10; 0.11; 0.12],...
                [0.02; 0.02; 0.03; 0.02]);
            boi = localCovmat(...
                [0.08; 0.09; 0.10; 0.11],...
                [0.07; 0.08; 0.09; 0.10],...
                [0.02; 0.02; 0.02; 0.02]);
            
            sd1 = [0.70; 0.72; 0.74; 0.76];
            sd2 = [0.68; 0.70; 0.72; 0.74];
            er_var = log([sd1 sd2]);
            er_cov = [0.03; 0.03; 0.04; 0.04];
            
            obs = [8 10];
            nocc = [2 3];
            ci = 0.95;
            
            out = psyrat_diffrel_trt(...
                'bp', bp, ...
                'bpi', bpi, ...
                'bpo', bpo, ...
                'bt', bt, ...
                'bo', bo, ...
                'boi', boi, ...
                'er_var', er_var, ...
                'er_cov', er_cov, ...
                'obs', obs, ...
                'nocc', nocc, ...
                'reltype', 1, ...
                'est', 'gen', ...
                'CI', ci);
            
            rel = localTable6TrtRel(...
                bp,bpi,bpo,bt,bo,boi,exp(er_var(:,1)).^2,exp(er_var(:,2)).^2,...
                er_cov,obs,nocc,1,'gen');
            
            ciedge = (1-ci)/2;
            testCase.verifyEqual(out.ll, quantile(rel,ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.pt, mean(rel), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.ul, quantile(rel,1-ciedge), 'AbsTol', 1e-12);
        end
        
        function testDiffRelTrtInputValidation(testCase)
            bp = localCovmat([0.5; 0.6],[0.4; 0.5],[0.1; 0.1]);
            bpi = bp; bpo = bp; bt = bp; bo = bp; boi = bp;
            er_var = log([0.7 0.8; 0.75 0.85]);
            
            testCase.verifyError(@() psyrat_diffrel_trt(...
                'bp', bp, 'bpi', bpi, 'bpo', bpo, 'bt', bt, ...
                'bo', bo, 'boi', boi, 'er_var', er_var, ...
                'obs', [8 10], 'nocc', [2 2], 'reltype', 4, ...
                'est', 'dep', 'CI', 0.95), 'varargin:reltype');
            
            testCase.verifyError(@() psyrat_diffrel_trt(...
                'bp', bp, 'bpi', bpi, 'bpo', bpo, 'bt', bt, ...
                'bo', bo, 'boi', boi, 'er_var', er_var, ...
                'obs', [8 10], 'nocc', [2 2], 'reltype', 3, ...
                'est', 'dep', 'CI', 1.2), 'varargin:ci');
        end

        function testSubjectSpecificReliabilityReturnsAugmentedTable(testCase)
            bp = { [0.50; 0.60; 0.55] };
            wp_pop = log([0.70; 0.75; 0.80]);
            wp_ss = { [0.02 0.03 0.01; 0.01 0.02 0.03; 0.03 0.01 0.02] };
            idtable = table({'s1'; 's2'; 's3'}, [1; 2; 3], [5; 6; 7], ...
                'VariableNames', {'id', 'id2', 'trls'});

            out = psyrat_ssrel( ...
                'bp', bp, ...
                'wp_pop', wp_pop, ...
                'wp_ss', wp_ss, ...
                'idtable', idtable, ...
                'CI', 0.95);

            newVars = {'dep_pt', 'dep_ll', 'dep_ul', 'icc_pt', 'icc_ll', 'icc_ul', ...
                'sem_pt', 'sem_ll', 'sem_ul', ...
                'bp_var', 'ss_errvar', 'trl_var', 'pop_errvar'};
            for i = 1:numel(newVars)
                testCase.verifyTrue(any(strcmp(out.Properties.VariableNames, newVars{i})));
            end
            testCase.verifyEqual(height(out), height(idtable));

            %Without a trial main effect supplied, trl_var is zero and the
            %absolute (default) and relative coefficients coincide.
            testCase.verifyEqual(out.trl_var, zeros(height(idtable), 1), 'AbsTol', 1e-12);
            relOut = psyrat_ssrel('bp', bp, 'wp_pop', wp_pop, 'wp_ss', wp_ss, ...
                'idtable', idtable, 'CI', 0.95, 'gcoeff', 2);
            testCase.verifyEqual(out.dep_pt, relOut.dep_pt, 'AbsTol', 1e-12);
        end

        function testSubjectSpecificReliabilityValidation(testCase)
            testCase.verifyError(@() psyrat_ssrel('bp', {1}, 'wp_pop', 1, 'wp_ss', {1}, 'CI', 0.95), ...
                'varargin:idtable');
            %invalid decision type is rejected
            idtable = table({'s1'}, 1, 5, 'VariableNames', {'id', 'id2', 'trls'});
            testCase.verifyError(@() psyrat_ssrel('bp', {1}, 'wp_pop', 1, 'wp_ss', {1}, ...
                'idtable', idtable, 'CI', 0.95, 'gcoeff', 9), 'varargin:gcoeff');
        end

        function testSubjectSpecificTrtReliabilityReturnsAugmentedTable(testCase)
            idtable = table({'s1';'s2';'s3'},[5;6;7], ...
                'VariableNames',{'id','trls'});
            bp = [0.45;0.50;0.55];
            bo = [0.10;0.11;0.12];
            bt = [0.20;0.21;0.22];
            txp = [0.12;0.13;0.14];
            oxp = [0.09;0.10;0.11];
            txo = [0.07;0.08;0.09];
            err = [0.70;0.72;0.74];

            out = psyrat_ssrel_trt( ...
                'gcoeff',1,'reltype',1,...
                'bp',bp,'bo',bo,'bt',bt,'txp',txp,'oxp',oxp,'txo',txo,'err',err,...
                'idtable',idtable,'CI',0.95);

            req = {'rel_pt','rel_ll','rel_ul','dep_pt','dep_ll','dep_ul',...
                'icc_pt','icc_ll','icc_ul','sem_pt','sem_ll','sem_ul'};
            for i = 1:numel(req)
                testCase.verifyTrue(any(strcmp(out.Properties.VariableNames,req{i})));
            end
            testCase.verifyEqual(height(out),height(idtable));
        end
        
        function testSubjectSpecificTrtReliabilityCanUsePerSubjectComponents(testCase)
            idtable = table({'s1';'s2'},[1;2],[6;6], ...
                'VariableNames',{'id','id2','trls'});
            bp = [0.45;0.50;0.55];
            bo = [0.10;0.11;0.12];
            bt = [0.20;0.21;0.22];
            txp = [0.18;0.18;0.18];
            oxp = [0.15;0.15;0.15];
            txo = [0.07;0.08;0.09];
            err = [0.75;0.75;0.75];
            txp_ss = [0.08 0.35; 0.09 0.36; 0.10 0.37];
            oxp_ss = [0.07 0.30; 0.08 0.31; 0.09 0.32];
            err_ss = [0.45 1.05; 0.46 1.06; 0.47 1.07];
            
            out = psyrat_ssrel_trt( ...
                'gcoeff',1,'reltype',1,...
                'bp',bp,'bo',bo,'bt',bt,'txp',txp,'oxp',oxp,'txo',txo,'err',err,...
                'idtable',idtable,'CI',0.95,...
                'txp_ss',txp_ss,'oxp_ss',oxp_ss,'err_ss',err_ss);
            
            testCase.verifyNotEqual(out.rel_pt(1),out.rel_pt(2));
            testCase.verifyNotEqual(out.sem_pt(1),out.sem_pt(2));
        end

        function testDodResidualUsesSquaredContrastWeights(testCase)
            % Regression: psyrat_dodiffrel summed residual variances without the
            % contrast weights, correct only for +/-1 contrasts. The residual
            % term of a contrast c over independent residuals must be
            % sum_i c_i^2 * sigma^2_err_i (divided by obs_i when scaled).
            bp = reshape(eye(4), [1 4 4]);    % between-person cov, 1 draw
            bt = reshape(zeros(4), [1 4 4]);  % no trial cov (isolate residual)
            v = [1 2 3 4];                    % target residual variances
            er = 0.5 * log(v);                % log residual SD so exp(er)^2 = v
            obs = [10 10 10 10];

            % contrast that excludes events 3 and 4
            out = psyrat_dodiffrel('bp', bp, 'bt', bt, 'er_var', er, ...
                'obs', obs, 'est', 'dep', 'CI', 0.95, 'cvec', [1 -1 0 0]);
            testCase.verifyEqual(out.er_contrast_var_pt, v(1) + v(2), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.er_scaled_var_pt, v(1)/10 + v(2)/10, 'AbsTol', 1e-12);

            % a non-unit weight scales by its square
            out2 = psyrat_dodiffrel('bp', bp, 'bt', bt, 'er_var', er, ...
                'obs', obs, 'est', 'dep', 'CI', 0.95, 'cvec', [2 0 0 0]);
            testCase.verifyEqual(out2.er_contrast_var_pt, 4 * v(1), 'AbsTol', 1e-12);

            % default +/-1 contrast is unchanged (c_i^2 = 1 for all events)
            out3 = psyrat_dodiffrel('bp', bp, 'bt', bt, 'er_var', er, ...
                'obs', obs, 'est', 'dep', 'CI', 0.95);
            testCase.verifyEqual(out3.er_contrast_var_pt, sum(v), 'AbsTol', 1e-12);
        end

        function testHarmMeanMatchesLiteralValues(testCase)
            % Finding RC-16: psyrat_harmmean was almost invisible to this suite. It is
            % called by 603 of 604 tests only indirectly, and PsyRATAccuracyOracle CALLS it
            % rather than reimplementing it -- so the oracle agreed with production by
            % construction for this one helper, and a wrong harmonic mean would have been
            % invisible to the accuracy lane. The single test that did constrain it
            % (TestSubjDiffDynrelReliability.m:34) does so INCIDENTALLY, by inlining
            % 2/(1/n1 + 1/n2) as a hand reference at unequal trial counts; rebalance that
            % fixture and the coverage silently drops to zero.
            %
            % The fix is literal expected VALUES, not a reimplementation. A reimplementation
            % is just a second copy of the same expression and fails in the same direction;
            % a hardcoded 24 cannot. These three stay mutation-sensitive permanently.
            %
            % This matters most right now because findings RC-01 / S13 are adjudicating
            % whether the harmonic mean is the right covariance divisor at all. Pinning the
            % ARITHMETIC first means that when the estimand question is settled, any change
            % shows up as a deliberate edit to a test with a literal in it.
            testCase.verifyEqual(psyrat_harmmean([20 30]), 24, 'AbsTol', 1e-12, ...
                'harmonic mean of 20 and 30 is 2/(1/20+1/30) = 24');
            testCase.verifyEqual(psyrat_harmmean([1 1]), 1, 'AbsTol', 1e-12, ...
                'harmonic mean of equal values is that value');
            testCase.verifyEqual(psyrat_harmmean([2 4 8]), 24/7, 'AbsTol', 1e-12, ...
                'harmonic mean of 2, 4 and 8 is 3/(1/2+1/4+1/8) = 24/7');

            % Order invariance, and the property that separates a harmonic mean from the
            % arithmetic mean it is most likely to be mistaken for: HM < AM whenever the
            % inputs differ. If someone swapped the implementation for mean(x), the three
            % literals above would catch it -- this makes the reason explicit.
            testCase.verifyEqual(psyrat_harmmean([30 20]), 24, 'AbsTol', 1e-12);
            testCase.verifyLessThan(psyrat_harmmean([20 30]), mean([20 30]));
        end

        function testDodDiffTrialCovarianceUsesHarmonicMeanDivisor(testCase)
            % Finding RC-16, second half. psyrat_dodiffrel.m:141 does NOT call
            % psyrat_harmmean -- it recomputes 2/((1/obs(i))+(1/obs(j))) inline. Three
            % independent harmonic-mean implementations exist in the toolbox
            % (psyrat_harmmean.m:35, psyrat_dodiffrel.m:141, psyrat_splits_summary.m:147)
            % and they are deliberately NOT consolidated (minimal-diff rule). The
            % consequence is that a future edit to psyrat_harmmean would not propagate here,
            % so this copy needs its own pin -- and it needs it BEFORE the RC-01 scaler work,
            % because this is one of the sites that work would have to change.
            %
            % The probe isolates a single off-diagonal element so exactly one divisor
            % survives into the answer. bt below is not a valid covariance matrix (it is not
            % positive semi-definite); that is deliberate and harmless, because
            % psyrat_dodiffrel validates only shape and the quantity under test is a
            % quadratic form, not a sampled covariance.
            bp = reshape(eye(4), [1 4 4]);     % u = cvec'*I*cvec = sum(c.^2) = 4
            obs = [10 20 30 40];               % UNEQUAL: at equal counts every candidate
                                               % divisor coincides and the test is vacuous
            er = 0.5 * log([1 1 1 1]);         % residual variance 1 per event

            % Only bt(1,2) = bt(2,1) = 1 is nonzero, so with the default contrast
            % [1 -1 -1 1] the scaled trial contrast is 2*c1*c2/nharm(1,2) = -2/nharm(1,2).
            btOff = zeros(4); btOff(1,2) = 1; btOff(2,1) = 1;
            outOff = psyrat_dodiffrel('bp', bp, 'bt', reshape(btOff, [1 4 4]), ...
                'er_var', er, 'obs', obs, 'est', 'dep', 'CI', 0.95);

            nharm12 = 2 / ((1/obs(1)) + (1/obs(2)));      % 2/(1/10+1/20) = 40/3
            testCase.verifyEqual(outOff.bt_scaled_var_pt, -2 / nharm12, 'AbsTol', 1e-12, ...
                ['the off-diagonal trial covariance must be divided by the harmonic mean ' ...
                 'of the two events'' trial counts']);

            % WHAT THIS TEST CAN AND CANNOT DISTINGUISH -- established by mutating the
            % production line and observing the result, not by inspection.
            %
            % It CANNOT distinguish the harmonic mean from "divide each off-diagonal element
            % by its own row's count" (nharm = obs(i)). Those are the SAME NUMBER here, and
            % not by coincidence: because 2/HM(a,b) == 1/a + 1/b, and because a covariance
            % matrix is symmetric so the (i,j) and (j,i) elements both enter the quadratic
            % form, replacing the harmonic mean with each element's own row count sums to an
            % identical total. Mutating psyrat_dodiffrel.m:141 that way leaves this test
            % green. That is a property of the algebra, not a gap in the test -- the two
            % implementations are equivalent for any symmetric input, so there is nothing to
            % catch. Worth knowing before anyone "fixes" one into the other.
            %
            % It CAN distinguish the harmonic mean from the ARITHMETIC mean, which is the
            % alternative the RC-01 / S13 discussion actually turns on: mutating the line to
            % (obs(i)+obs(j))/2 fails this test. It also catches a single shared count.
            for wrong = [obs(1), obs(2), (obs(1)+obs(2))/2]
                testCase.verifyNotEqual(outOff.bt_scaled_var_pt, -2 / wrong, ...
                    sprintf('divisor %g must not reproduce the harmonic-mean result', wrong));
            end

            % And the other half of the i == j branch: a DIAGONAL element is divided by that
            % event's own count, never by a harmonic mean. Without this, a mutation that used
            % the harmonic mean everywhere would pass the check above.
            btDiag = zeros(4); btDiag(1,1) = 1;
            outDiag = psyrat_dodiffrel('bp', bp, 'bt', reshape(btDiag, [1 4 4]), ...
                'er_var', er, 'obs', obs, 'est', 'dep', 'CI', 0.95);
            testCase.verifyEqual(outDiag.bt_scaled_var_pt, 1 / obs(1), 'AbsTol', 1e-12, ...
                'a diagonal trial variance is divided by that event''s own trial count');
        end

        function testSubjectSpecificDiffReliabilityReturnsAugmentedTable(testCase)
            [bp, bt, erVar] = localDiffRelInputs();
            idtable = table({'s1';'s2';'s3'},[6;7;8],[6;7;8], ...
                'VariableNames',{'id','trls1','trls2'});
            idtable.trls = min(idtable.trls1,idtable.trls2);

            out = psyrat_ssrel_diff( ...
                'bp',bp,...
                'bt',bt,...
                'er_var',erVar,...
                'idtable',idtable,...
                'CI',0.95,...
                'est','dep');

            req = {'rel_pt','rel_ll','rel_ul','dep_pt','dep_ll','dep_ul',...
                'icc_pt','icc_ll','icc_ul','sem_pt','sem_ll','sem_ul'};
            for i = 1:numel(req)
                testCase.verifyTrue(any(strcmp(out.Properties.VariableNames,req{i})));
            end
            testCase.verifyEqual(height(out),height(idtable));
        end
        
        function testSubjectSpecificDiffReliabilityCanUsePerSubjectComponents(testCase)
            [bp, bt, erVar] = localDiffRelInputs();
            idtable = table({'s1';'s2'},[1;2],[6;6],[6;6], ...
                'VariableNames',{'id','id2','trls1','trls2'});
            idtable.trls = min(idtable.trls1,idtable.trls2);
            er_ss1 = [log(0.35) log(0.90); log(0.36) log(0.92); ...
                log(0.34) log(0.95); log(0.33) log(0.96); log(0.32) log(0.98)];
            er_ss2 = [log(0.30) log(0.80); log(0.31) log(0.82); ...
                log(0.29) log(0.84); log(0.28) log(0.86); log(0.27) log(0.88)];
            wp_cov_ss = [0.00 0.02; 0.00 0.02; 0.00 0.02; 0.00 0.02; 0.00 0.02];
            
            out = psyrat_ssrel_diff( ...
                'bp',bp,...
                'bt',bt,...
                'er_var',erVar,...
                'idtable',idtable,...
                'CI',0.95,...
                'est','dep',...
                'er_var_ss1',er_ss1,...
                'er_var_ss2',er_ss2,...
                'wp_cov_ss',wp_cov_ss);
            
            testCase.verifyNotEqual(out.rel_pt(1),out.rel_pt(2));
            testCase.verifyNotEqual(out.sem_pt(1),out.sem_pt(2));
        end

        function testRelSummarySubjectLevelDiffBuildsDiffTable(testCase)
            psyrat_data = PsyRATTestDataFactory.makeICDiffSSErrRelData();

            [outData, relerr] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'sing_diff_sserr', ...
                'depcutoff', 0.5, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'diffgcoeff', 1, ...
                'CI', 0.95);

            testCase.verifyEqual(relerr.nogooddata, 0);
            testCase.verifyEqual(outData.rel.analysis, 'ic_diff_sserrvar');
            testCase.verifyTrue(isfield(outData.relsummary.group(1), 'diffscore'));
            testCase.verifyTrue(isfield(outData.relsummary.group(1).diffscore, 'ssrel_table'));
            ss = outData.relsummary.group(1).diffscore.ssrel_table;
            req = {'id','id2','trls1','trls2','rel_pt','icc_pt','sem_pt','ind2include','ind2exclude'};
            for i = 1:numel(req)
                testCase.verifyTrue(any(strcmp(ss.Properties.VariableNames, req{i})));
            end
            testCase.verifyEqual(height(ss), 3);
            testCase.verifyNotEqual(ss.rel_pt(1), ss.rel_pt(3));
        end

        function testRelSummarySserrGroupOnlyDoesNotPoolTrials(testCase)
            %S12 regression: subject-level (sserrvar=2), group-only design
            %(analysis 6 'ic_sserrvar', relsummary 'sing_sserr' case 2). Two
            %groups reuse the same subject ids with different, known
            %per-group trial counts; the pre-fix bug summed them across
            %groups (id 's1': 5+9=14 pooled into both groups).
            psyrat_data = PsyRATTestDataFactory.makeICSSErrRelDataMultiGroup();

            [outData, relerr] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'sing_sserr', ...
                'depcutoff', 0.01, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'gcoeff', 1, ...
                'CI', 0.95);

            testCase.verifyEqual(relerr.nogooddata, 0);

            ssA = outData.relsummary.group(1).event(1).ssrel_table;
            ssB = outData.relsummary.group(2).event(1).ssrel_table;
            trlsA = containers.Map(cellstr(string(ssA.id)), num2cell(ssA.trls));
            trlsB = containers.Map(cellstr(string(ssB.id)), num2cell(ssB.trls));

            testCase.verifyEqual(trlsA('s1'), 5);
            testCase.verifyEqual(trlsA('s2'), 7);
            testCase.verifyEqual(trlsB('s1'), 9);
            testCase.verifyEqual(trlsB('s2'), 11);

            %independent numeric cross-check: dep_pt against the accuracy
            %oracle fed the SAME known-correct per-group counts, not merely
            %"differs across groups". wp_ss columns are reordered by ssA.id2
            %since innerjoin does not guarantee ssA's row order matches the
            %original per-subject column order in ind_sdlog.
            wpSsAll = cell2mat(psyrat_data.rel.out.ind_sdlog{1});
            wpSsOrdered = num2cell(wpSsAll(:,ssA.id2));
            oracleA = PsyRATAccuracyOracle.ssrel( ...
                psyrat_data.rel.out.gro_sds{1}, psyrat_data.rel.out.pop_sdlog(:,1), ...
                wpSsOrdered, ssA, 0.95, psyrat_data.rel.out.sig_trl(:,1), 1);
            testCase.verifyEqual(ssA.dep_pt, oracleA.dep_pt, 'AbsTol', 1e-10);
        end

        function testRelSummarySserrGroupEventDoesNotPoolTrials(testCase)
            %S12 regression: subject-level, group+event design (analysis 6,
            %relsummary 'sing_sserr' case 4). Two groups x two events reuse
            %the same subject ids with distinct known counts in every cell;
            %the pre-fix bug filtered by event only, pooling across groups
            %within each event.
            psyrat_data = PsyRATTestDataFactory.makeICSSErrRelDataMultiGroupEvent();

            [outData, relerr] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'sing_sserr', ...
                'depcutoff', 0.01, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'gcoeff', 1, ...
                'CI', 0.95);

            testCase.verifyEqual(relerr.nogooddata, 0);

            expected = struct( ...
                'GroupA', struct('cor', [5 7], 'err', [4 6]), ...
                'GroupB', struct('cor', [9 11], 'err', [8 10]));
            groupNames = {'GroupA','GroupB'}; eventNames = {'cor','err'};
            for g = 1:2
                for e = 1:2
                    ss = outData.relsummary.group(g).event(e).ssrel_table;
                    trlsById = containers.Map(cellstr(string(ss.id)), num2cell(ss.trls));
                    exp2 = expected.(groupNames{g}).(eventNames{e});
                    testCase.verifyEqual(trlsById('s1'), exp2(1), ...
                        sprintf('%s/%s id s1', groupNames{g}, eventNames{e}));
                    testCase.verifyEqual(trlsById('s2'), exp2(2), ...
                        sprintf('%s/%s id s2', groupNames{g}, eventNames{e}));
                end
            end
        end

        function testRelSummaryTrtSserrGroupDoesNotPoolTrials(testCase)
            %S12 regression: subject-level crossed test-retest (trt_sserrvar,
            %analysis 25, relsummary 'trt_sserr'). Two groups reuse the same
            %subject ids with different known trials-per-occasion; the
            %pre-fix bug never filtered by group at all.
            psyrat_data = PsyRATTestDataFactory.makeTRTSSErrRelDataMultiGroup();

            [outData, relerr] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'trt_sserr', ...
                'gcoeff', 1, ...
                'reltype', 3, ...
                'depcutoff', 0.01, ...
                'meascutoff', 2, ...
                'CI', 0.95);

            testCase.verifyEqual(relerr.nogooddata, 0);

            ssA = outData.relsummary.group(1).event(1).ssrel_table;
            ssB = outData.relsummary.group(2).event(1).ssrel_table;
            trlsA = containers.Map(cellstr(string(ssA.id)), num2cell(ssA.trls));
            trlsB = containers.Map(cellstr(string(ssB.id)), num2cell(ssB.trls));

            testCase.verifyEqual(trlsA('s1'), 4);
            testCase.verifyEqual(trlsA('s2'), 6);
            testCase.verifyEqual(trlsB('s1'), 8);
            testCase.verifyEqual(trlsB('s2'), 10);
        end

        function testRelSummarySingGroupEventLabelParsingMatchesEstimator(testCase)
            %Label-contract guard: the plain-IC group+event estimation path
            %(psyrat_computevarcomp case 4) emits REL.out.labels as
            %event_;_group (e.g. 'cor_;_GroupA'), so the 'sing' case-4 parse
            %reads lblstr(1)=event, lblstr(2)=group. This pins that each
            %(group,event) cell is labelled correctly; it fails if the parse
            %is swapped to group-first (both ismember lookups would return
            %empty and error on real event_;_group labels).
            psyrat_data = PsyRATTestDataFactory.makeICRelDataMultiGroupEvent();

            [outData, relerr] = psyrat_relsummary( ...
                'psyrat_data', psyrat_data, ...
                'analysis', 'sing', ...
                'depcutoff', 0.01, ...
                'meascutoff', 2, ...
                'depcentmeas', 1, ...
                'gcoeff', 1, ...
                'CI', 0.95);

            testCase.verifyEqual(relerr.nogooddata, 0);
            testCase.verifyEqual(outData.relsummary.group(1).name, 'GroupA');
            testCase.verifyEqual(outData.relsummary.group(2).name, 'GroupB');
            testCase.verifyEqual(outData.relsummary.group(1).event(1).name, 'cor');
            testCase.verifyEqual(outData.relsummary.group(1).event(2).name, 'err');
            testCase.verifyEqual(outData.relsummary.group(2).event(1).name, 'cor');
            testCase.verifyEqual(outData.relsummary.group(2).event(2).name, 'err');
        end

        function testRelsummaryPublishesOneGroupPerRelGroup(testCase)
            %RC-44 publish-contract guard. psyrat_relsummary hands the summary
            %off in exactly one place ('%store relsummary information' near the
            %end of the function), and its consumers split on how they size
            %their group loop. Most size from rel.groups and then index
            %relsummary.group(gloc) or relsummary.data.g(gloc) with no bounds
            %check -- 23 such loops across 13 files, the largest counts in
            %psyrat_variancet (5), psyrat_report (3) and psyrat_ssrelplot (2)
            %-- so a summary holding fewer group elements than rel.groups
            %has entries is an index-out-of-bounds throw there. A few size from
            %the published array instead and are deliberately fail-open
            %(psyrat_admissibility_collect:57, psyrat_startview_criterion:347),
            %so the same short array would drop a group from those surfaces
            %without any notice. The failure is loud in some places and silent
            %in others, which is why the one-element-per-group contract is
            %pinned here rather than left implicit. If a future design
            %legitimately needs to skip a group, this test failing is the
            %signal that those consumer loops must be changed first; do not
            %relax the assertion to make the new behavior pass.
            %
            %Three fixtures cover the branches that build group arrays
            %differently: subject-level group-only (sing_sserr case 2),
            %subject-level group+event (sing_sserr case 4), and plain IC
            %group+event (sing case 4).
            cases = { ...
                'makeICSSErrRelDataMultiGroup',      'sing_sserr'; ...
                'makeICSSErrRelDataMultiGroupEvent', 'sing_sserr'; ...
                'makeICRelDataMultiGroupEvent',      'sing'};

            for c = 1:size(cases, 1)
                fixture = cases{c, 1};
                analysis = cases{c, 2};
                psyrat_data = PsyRATTestDataFactory.(fixture)();

                %expected count derived exactly as the rel.groups-sized
                %consumers derive it (cf. psyrat_ssrelplot's group check):
                %'none' collapses to one implicit group, otherwise one element
                %per rel.groups entry.
                if strcmpi(psyrat_data.rel.groups, 'none')
                    expNGroups = 1;
                else
                    expNGroups = length(psyrat_data.rel.groups);
                end

                %all three fixtures are genuinely multi-group; assert that so
                %the check cannot quietly degrade into a trivial 1 == 1 if a
                %fixture is ever changed.
                testCase.verifyEqual(expNGroups, 2, ...
                    sprintf('%s must remain a multi-group fixture', fixture));

                [outData, relerr] = psyrat_relsummary( ...
                    'psyrat_data', psyrat_data, ...
                    'analysis', analysis, ...
                    'depcutoff', 0.01, ...
                    'meascutoff', 2, ...
                    'depcentmeas', 1, ...
                    'gcoeff', 1, ...
                    'CI', 0.95);

                %a nogooddata run returns before the publish site, so the
                %contract below is only meaningful on a run that published.
                testCase.verifyEqual(relerr.nogooddata, 0, ...
                    sprintf('%s: fixture must yield usable data', fixture));

                testCase.verifyEqual(numel(outData.relsummary.group), expNGroups, ...
                    sprintf('%s: relsummary.group must hold one element per rel.groups entry', ...
                    fixture));
                testCase.verifyEqual(numel(outData.relsummary.data.g), expNGroups, ...
                    sprintf('%s: relsummary.data.g must hold one element per rel.groups entry', ...
                    fixture));
            end
        end

        function testCheckConvConvergedAndNotConverged(testCase)
            rel = struct();
            rel.nchains = 4;
            rel.out.conv.data = {{ ...
                'name', 'n_eff', 'Rhat'; ...
                'mu', 120, 1.00; ...
                'sig', 140, 1.01}};

            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 1);

            rel.out.conv.data = {{ ...
                'name', 'n_eff', 'Rhat'; ...
                'mu', 120, 1.20; ...
                'sig', 140, 1.01}};
            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 0);
        end

        function testCheckConvAcceptsNonNestedFormat(testCase)
            rel = struct();
            rel.nchains = 3;
            rel.out.conv.data = { ...
                'name', 'n_eff', 'Rhat'; ...
                'mu', 90, 1.00; ...
                'sig', 95, 1.01};

            out = psyrat_checkconv(rel);
            testCase.verifyEqual(out.out.conv.converged, 1);
        end

        function testDoDDiffRelMatchesManualFormula(testCase)
            bp = localDoDCovDraws();
            bt = bp .* 0.8;
            erVar = log([ ...
                0.72 0.80 0.76 0.82; ...
                0.74 0.81 0.78 0.83; ...
                0.73 0.79 0.77 0.84; ...
                0.75 0.82 0.79 0.85; ...
                0.76 0.83 0.80 0.86]);
            obs = [8 10 9 11];
            ci = 0.95;

            outDep = psyrat_dodiffrel( ...
                'bp', bp, 'bt', bt, 'er_var', erVar, ...
                'obs', obs, 'est', 'dep', 'CI', ci);
            outGen = psyrat_dodiffrel( ...
                'bp', bp, 'bt', bt, 'er_var', erVar, ...
                'obs', obs, 'est', 'gen', 'CI', ci);

            [depDraws, depIcc, genDraws, genIcc] = localDoDManualRel(bp, bt, erVar, obs);
            ciedge = (1 - ci) / 2;

            testCase.verifyEqual(outDep.pt, mean(depDraws), 'AbsTol', 1e-12);
            testCase.verifyEqual(outDep.ll, quantile(depDraws, ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(outDep.ul, quantile(depDraws, 1 - ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(outDep.icc_pt, mean(depIcc), 'AbsTol', 1e-12);

            testCase.verifyEqual(outGen.pt, mean(genDraws), 'AbsTol', 1e-12);
            testCase.verifyEqual(outGen.ll, quantile(genDraws, ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(outGen.ul, quantile(genDraws, 1 - ciedge), 'AbsTol', 1e-12);
            testCase.verifyEqual(outGen.icc_pt, mean(genIcc), 'AbsTol', 1e-12);
        end

        function testDoDVirtualRelBuilderSingleAndDiffModes(testCase)
            psyrat_data = PsyRATTestDataFactory.makeDoDRelData();

            outSingle = psyrat_dod_build_virtual_rel( ...
                'psyrat_data', psyrat_data, ...
                'mode', 'single', ...
                'groupidx', 2, ...
                'eventidx', 3);

            testCase.verifyEqual(outSingle.rel.analysis, 'ic');
            testCase.verifyEqual(outSingle.rel.groups, 'none');
            testCase.verifyEqual(outSingle.rel.events, {'ERP_C'});
            testCase.verifyEqual(unique(string(outSingle.rel.data.group)), "Clinical");
            testCase.verifyEqual(unique(string(outSingle.rel.data.event)), "ERP_C");
            testCase.verifyTrue(isfield(outSingle.rel.out, 'sig_u'));
            testCase.verifyTrue(isfield(outSingle.rel.out, 'sig_e'));

            outDiff = psyrat_dod_build_virtual_rel( ...
                'psyrat_data', psyrat_data, ...
                'mode', 'diff', ...
                'groupidx', 1, ...
                'eventpair', [1 4]);

            testCase.verifyEqual(outDiff.rel.analysis, 'ic_diff');
            testCase.verifyEqual(outDiff.rel.events, {'ERP_A'; 'ERP_D'});
            testCase.verifyEqual(unique(string(outDiff.rel.data.group)), "Control");
            testCase.verifyTrue(isfield(outDiff.rel.out, 'id_varcov'));
            testCase.verifyTrue(isfield(outDiff.rel.out, 'b_sigma'));
        end

        function testDoDVirtualRelSingleCarriesTrialMainEffectAndMean(testCase)
            %RC-15 + RC-18. The DoD -> single-event view builds a virtual
            %one-facet REL. It used to drop the trial main effect (sigma_i) and
            %hardcode the universe-score mean to zero, so psyrat_relsummary
            %backfilled sigma_i with zeros and every criterion coefficient on the
            %route was evaluated as though the mean were exactly 0.
            %
            %Expected values are hand-derived from the fixture rather than
            %recomputed with the production expression, so an indexing error
            %(off-diagonal element, wrong group, wrong event) fails here instead
            %of being mirrored into the expectation.
            psyrat_data = PsyRATTestDataFactory.makeDoDRelData();

            out = psyrat_dod_build_virtual_rel( ...
                'psyrat_data', psyrat_data, ...
                'mode', 'single', ...
                'groupidx', 2, ...
                'eventidx', 3);

            %sigma_i is the SD of the SELECTED cell's trial variance: the (3,3)
            %diagonal of group 2's trial covariance. Fixture: baseTrl(3,3) = 0.88,
            %plus the group-2 shift 0.03, plus the per-draw drift 0.005*(d-1).
            expVar = 0.88 + 0.03 + (0:5)' * 0.005;
            testCase.verifyEqual(out.rel.out.sig_trl, sqrt(expVar), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.rel.out.sig_trl_source, 'estimated');

            %The mean is the selected cell's intercept, not a placeholder. Fixture:
            %baseCellMean(3) = 2.40, plus the group-2 shift 0.05, plus 0.01*(d-1).
            expMu = 2.40 + 0.05 + (0:5)' * 0.01;
            testCase.verifyEqual(out.rel.out.mu, expMu, 'AbsTol', 1e-12);
            testCase.verifyEqual(out.rel.out.mu_source, 'estimated');

            %Guard the specific failure mode both findings describe: a mean that
            %is identically zero looks plausible and is not otherwise detectable.
            testCase.verifyTrue(all(out.rel.out.mu ~= 0), ...
                'mu must be the estimated cell mean, not the zeros placeholder.');
        end

        function testDoDVirtualRelSingleFlagsUnestimatedMean(testCase)
            %The gamma DoD model parameterizes a per-cell LOG-mean rather than an
            %observed-scale mean, and result structs saved before the extraction
            %existed carry no b_cell either. Both must fall back to the legacy
            %zeros AND say so, because the criterion button is gated on the flag:
            %silently returning zeros is what made RC-18 invisible.
            psyrat_data = PsyRATTestDataFactory.makeDoDRelData();
            psyrat_data.rel.out = rmfield(psyrat_data.rel.out, 'b_cell');

            out = psyrat_dod_build_virtual_rel( ...
                'psyrat_data', psyrat_data, ...
                'mode', 'single', ...
                'groupidx', 1, ...
                'eventidx', 2);

            testCase.verifyEqual(out.rel.out.mu, zeros(6, 1));
            testCase.verifyEqual(out.rel.out.mu_source, 'unavailable');

            %sigma_i does not depend on b_cell, so it must still be propagated.
            testCase.verifyEqual(out.rel.out.sig_trl_source, 'estimated');
            testCase.verifyTrue(all(out.rel.out.sig_trl > 0));
        end

        function testDoDVirtualRelSingleSeparatesDependabilityFromGeneralizability(testCase)
            %The consequence that made RC-15 undetectable by inspection: with
            %sigma_i zeroed, absolute error equals relative error, so the two
            %coefficients are bit-identical across the whole view and no number a
            %user could look at reveals the omission. Assert they now differ, and
            %in the correct direction (absolute error is larger, so dependability
            %is the smaller coefficient).
            psyrat_data = PsyRATTestDataFactory.makeDoDRelData();
            out = psyrat_dod_build_virtual_rel( ...
                'psyrat_data', psyrat_data, ...
                'mode', 'single', ...
                'groupidx', 2, ...
                'eventidx', 3);
            o = out.rel.out;

            wp = sqrt(o.sig_e.^2 + o.sig_trl.^2);
            [~, dep] = psyrat_rel_sing('gcoeff', 1, 'metric', 'global', ...
                'bp', o.sig_u, 'wp', wp, 'i', o.sig_trl, 'obs', 20, 'CI', 0.95);
            [~, gen] = psyrat_rel_sing('gcoeff', 2, 'metric', 'global', ...
                'bp', o.sig_u, 'wp', wp, 'i', o.sig_trl, 'obs', 20, 'CI', 0.95);

            testCase.verifyLessThan(dep, gen, ...
                'Dependability must be strictly below generalizability once sigma_i enters absolute error.');
            testCase.verifyGreaterThan(gen - dep, 1e-6, ...
                'A vanishing gap means sigma_i is not reaching the coefficient.');
        end

        function testDoDVirtualRelSingleOnGammaLayoutInheritsFixAndFlagsMean(testCase)
            %The gamma family validation for this tier. psyrat_dod_build_virtual_rel
            %contains no reference to 'family' or 'gamma' -- it is family-blind by
            %construction, and the gamma DoD extractor deliberately writes its
            %observed-scale components into the same REL.out slots so the whole
            %downstream stack runs unchanged. That means gamma inherits both the
            %defect and the fix, and it must be checked rather than assumed.
            %
            %Two things differ under gamma and are asserted separately:
            % (1) trl_varcov IS populated, so sigma_i propagates exactly as in the
            %     Gaussian case -- the RC-15 fix is inherited for free.
            % (2) there is NO b_cell. The gamma DoD model parameterizes a per-cell
            %     LOG-mean, and the extraction branch does not store even that, so
            %     an observed-scale mean needs its own derivation. The mean must
            %     therefore be reported as unavailable, not silently zeroed.
            pd = PsyRATTestDataFactory.makeDoDRelData();

            %Reshape REL.out to the exact field set the gamma DoD branch writes
            %(psyrat_gamma_extract_dodiff -> REL.out): no b_cell, plus err_var,
            %the observed-scale correlations, and the dispersion diagnostics.
            g = pd.rel.out;
            ndraw = size(g.b_sigma{1}, 1);
            gammaOut = struct();
            gammaOut.b_sigma    = g.b_sigma;
            gammaOut.id_varcov  = g.id_varcov;
            gammaOut.trl_varcov = g.trl_varcov;
            gammaOut.err_var    = {exp(g.b_sigma{1}).^2, exp(g.b_sigma{2}).^2};
            gammaOut.cor_id     = g.id_varcov;
            gammaOut.cor_trl    = g.trl_varcov;
            gammaOut.gamma_disp = struct('s_v', 0.4 * ones(ndraw,1), ...
                'cov_pv', zeros(ndraw,1), 'disp_factor', ones(ndraw,1));
            gammaOut.conv = g.conv;

            pd.rel.out = gammaOut;
            pd.rel.family = 'gamma';

            out = psyrat_dod_build_virtual_rel( ...
                'psyrat_data', pd, ...
                'mode', 'single', ...
                'groupidx', 2, ...
                'eventidx', 3);

            %(1) Same hand-derived value as the Gaussian case: the fix is inherited.
            expVar = 0.88 + 0.03 + (0:5)' * 0.005;
            testCase.verifyEqual(out.rel.out.sig_trl, sqrt(expVar), 'AbsTol', 1e-12);
            testCase.verifyEqual(out.rel.out.sig_trl_source, 'estimated');

            %(2) The mean is honestly reported as unavailable rather than zeroed
            %behind the user's back. The GUI gates the criterion button on this.
            testCase.verifyEqual(out.rel.out.mu_source, 'unavailable');
            testCase.verifyEqual(out.rel.out.mu, zeros(ndraw, 1));

            %And the coefficients still separate, which is the point of (1).
            o = out.rel.out;
            wp = sqrt(o.sig_e.^2 + o.sig_trl.^2);
            [~, dep] = psyrat_rel_sing('gcoeff', 1, 'metric', 'global', ...
                'bp', o.sig_u, 'wp', wp, 'i', o.sig_trl, 'obs', 20, 'CI', 0.95);
            [~, gen] = psyrat_rel_sing('gcoeff', 2, 'metric', 'global', ...
                'bp', o.sig_u, 'wp', wp, 'i', o.sig_trl, 'obs', 20, 'CI', 0.95);
            testCase.verifyLessThan(dep, gen);
        end

        function testDoDVirtualRelBuilderValidation(testCase)
            psyrat_data = PsyRATTestDataFactory.makeDoDRelData();

            testCase.verifyError(@() psyrat_dod_build_virtual_rel( ...
                'psyrat_data', psyrat_data, ...
                'mode', 'diff', ...
                'groupidx', 1, ...
                'eventpair', [2 2]), ...
                'varargin:eventpair');
        end

        function testDiffRelReportsTheRawValueAndFlagsItAsInadmissible(testCase)
            % RC-01. At unequal trial counts the covariance term is scaled by
            % the harmonic mean of the two counts, which is smaller than
            % either, so the subtraction can exceed the variances it is taken
            % from. This fixture is the case documented in the findings
            % ledger: residual variances 16 and 1 at obs = [10 1] put the
            % analytic threshold at rho = 26/44 = 0.5909, and rho = 0.70
            % drives BOTH error variances negative.
            [bp, bt, er_var, wp_cov, obs] = localRC01Fixture();

            % the fixture really is pathological: check the unclipped algebra
            % before checking that the guard contains it. Without this the
            % test could pass on a benign fixture and prove nothing.
            hm = psyrat_harmmean(obs);
            relErrRaw = 16/obs(1) + 1/obs(2) - 2*wp_cov(1)/hm;
            absErrRaw = relErrRaw + 1/obs(1) + 1/obs(2) - 2*bt(1,2,1)/hm;
            testCase.verifyLessThan(relErrRaw, 0);
            testCase.verifyLessThan(absErrRaw, 0);
            uni = bp(1,1,1) + bp(1,2,2) - 2*bp(1,2,1);
            testCase.verifyGreaterThan(uni/(uni + relErrRaw), 1);
            testCase.verifyGreaterThan(uni/(uni + absErrRaw), 1);

            gen = testCase.verifyWarning(@() psyrat_diffrel( ...
                'bp', bp, 'bt', bt, 'er_var', er_var, 'wp_cov', wp_cov, ...
                'obs', obs, 'est', 'gen', 'CI', 0.95), 'psyrat:negerrvar');
            dep = testCase.verifyWarning(@() psyrat_diffrel( ...
                'bp', bp, 'bt', bt, 'er_var', er_var, 'wp_cov', wp_cov, ...
                'obs', obs, 'est', 'dep', 'CI', 0.95), 'psyrat:negerrvar');

            % NOTHING IS CLIPPED. The reported value is the formula evaluated
            % exactly, so it equals the raw algebra above and lies OUTSIDE
            % [0,1]. Asserting the bound would be asserting the old, hiding
            % behaviour; asserting the raw value is what pins the new contract.
            testCase.verifyEqual(gen.pt, uni/(uni + relErrRaw), 'AbsTol', 1e-12);
            testCase.verifyEqual(dep.pt, uni/(uni + absErrRaw), 'AbsTol', 1e-12);
            testCase.verifyGreaterThan(gen.pt, 1);
            testCase.verifyGreaterThan(dep.pt, 1);

            % and the diagnostic that must travel with it says so
            for out = [gen dep]
                testCase.verifyTrue(isfield(out, 'admissibility'));
                testCase.verifyTrue( ...
                    psyrat_admissibility_hasneg(out.admissibility));
                testCase.verifyNotEmpty( ...
                    psyrat_admissibility_note(out.admissibility));
            end
        end

        function testDiffRelClipDoesNotMoveIcc(testCase)
            % The ICC denominators are the same quantities at a single
            % observation, where the covariance cannot outrun the variances,
            % so they are never clipped. This pins the findings ledger's
            % correction that the ICC does NOT move under RC-01.
            [bp, bt, er_var, wp_cov, obs] = localRC01Fixture();

            gen = testCase.verifyWarning(@() psyrat_diffrel( ...
                'bp', bp, 'bt', bt, 'er_var', er_var, 'wp_cov', wp_cov, ...
                'obs', obs, 'est', 'gen', 'CI', 0.95), 'psyrat:negerrvar');

            uni = bp(:,1,1) + bp(:,2,2) - 2*bp(:,2,1);
            denomRelErr = 16 + 1 - 2*wp_cov;
            iccExpected = uni ./ (uni + denomRelErr);

            testCase.verifyEqual(gen.icc_pt, mean(iccExpected), 'AbsTol', 1e-12);
            testCase.verifyGreaterThan(gen.icc_pt, 0);
            testCase.verifyLessThan(gen.icc_pt, 1);
        end

        function testDiffRelLeavesNonNegativeErrorVarianceUntouched(testCase)
            % Regression: the guard must not perturb an ordinary result. Same
            % components as the pathological fixture, but at a common trial
            % count, where the harmonic, geometric and arithmetic means all
            % coincide and the error variance reduces to the single-
            % observation variance of the difference divided by n.
            [bp, bt, er_var, wp_cov, ~] = localRC01Fixture();
            n = 10;

            out = testCase.verifyWarningFree(@() psyrat_diffrel( ...
                'bp', bp, 'bt', bt, 'er_var', er_var, 'wp_cov', wp_cov, ...
                'obs', [n n], 'est', 'gen', 'CI', 0.95));

            uni = bp(:,1,1) + bp(:,2,2) - 2*bp(:,2,1);
            denomRelErr = 16 + 1 - 2*wp_cov;
            expected = uni ./ (uni + denomRelErr./n);

            testCase.verifyEqual(out.pt, mean(expected), 'AbsTol', 1e-12);
            testCase.verifyGreaterThan(out.pt, 0);
            testCase.verifyLessThan(out.pt, 1);
        end

        function testDiffRelTrtFlagsInadmissibleWithoutModifying(testCase)
            % The two-facet twin carries the same harmonic-mean scaler on
            % every D-study adjusted component, so it needs the same guard.
            % Unequal trial counts plus a large person-by-trial covariance
            % drive the relative error negative for the coefficient of trial
            % equivalence and stability.
            %
            % reltype 3 specifically: there uni = bp_raw, so guarding the
            % error variance is enough to bound the coefficient. For reltypes
            % 1 and 2 the universe term itself carries a harmonic-mean-scaled
            % block and is NOT guarded, so no [0,1] claim is made for them
            % here. That gap is filed as a finding, not fixed by this guard.
            nd = 8;
            one = ones(nd,1);
            cov2 = @(v1,v2,cv) localCovmat(v1*one, v2*one, cv*one);

            bp   = cov2(2.0, 2.0, 0.5);
            bpi  = cov2(16.0, 1.0, 3.6);   % rho = 0.90 at variances 16 and 1
            bpo  = cov2(0.2, 0.2, 0.05);
            bt   = cov2(0.3, 0.3, 0.10);
            bo   = cov2(0.2, 0.2, 0.05);
            boi  = cov2(0.1, 0.1, 0.02);
            er_var = log([0.30*one 0.30*one]);
            er_cov = 0.02*one;
            obs = [10 1];
            nocc = [2 2];

            out = testCase.verifyWarning(@() psyrat_diffrel_trt( ...
                'bp', bp, 'bpi', bpi, 'bpo', bpo, 'bt', bt, 'bo', bo, ...
                'boi', boi, 'er_var', er_var, 'er_cov', er_cov, ...
                'obs', obs, 'nocc', nocc, 'reltype', 3, 'est', 'gen', ...
                'CI', 0.95), 'psyrat:negerrvar');

            % Nothing is clipped, so the assertion is on the DIAGNOSTIC, not
            % on a bound. reltype 3 is used here because its universe term is
            % a raw variance of a difference; for reltypes 1 and 2 the
            % universe term carries a harmonic-mean-scaled block of its own
            % and can go negative too, which is covered separately.
            testCase.verifyTrue(isfield(out, 'admissibility'));
            testCase.verifyTrue(psyrat_admissibility_hasneg(out.admissibility));

            labels = {out.admissibility.entries.label};
            testCase.verifyTrue(any(strcmp(labels, 'relative error variance')));
            testCase.verifyTrue(any(strcmp(labels, 'universe-score variance')));
        end

        function testAdmissibilityRecordsWithoutModifyingAnything(testCase)
            % The whole point of the recorder: it counts, and it does not
            % touch the values. Nothing in this design may repair a quantity.
            v = [-2; -1e-9; 0; 0.5; 3];

            adm = psyrat_admissibility([], 'relative error variance', v, ...
                'unit test');

            testCase.verifyEqual(adm.entries(1).ndraw, 5);
            testCase.verifyEqual(adm.entries(1).nneg, 2);
            testCase.verifyEqual(adm.entries(1).nzero, 1);
            testCase.verifyEqual(adm.entries(1).minval, -2, 'AbsTol', 0);
            testCase.verifyTrue(psyrat_admissibility_hasneg(adm));

            % NaN is counted as neither: a quantity that was never computed is
            % a different problem from one computed and found impossible.
            admn = psyrat_admissibility([], 'x', [NaN; -1; 2], 'unit test');
            testCase.verifyEqual(admn.entries(1).nneg, 1);
            testCase.verifyEqual(admn.entries(1).ndraw, 3);

            % An all-positive quantity is silent and produces no note, so a
            % clean run reads exactly as it did before this existed.
            clean = psyrat_admissibility([], 'x', [0.5; 3], 'unit test');
            testCase.verifyFalse(psyrat_admissibility_hasneg(clean));
            testCase.verifyEmpty(psyrat_admissibility_note(clean));
            testCase.verifyWarningFree(@() psyrat_admissibility_warn(clean));

            % Exactly zero is reported, but as a real result rather than a
            % fault: it means the difference score has no between-person
            % variance. It must NOT raise the warning.
            zeroadm = psyrat_admissibility([], 'universe-score variance', ...
                [0; 0; 1], 'unit test');
            testCase.verifyFalse(psyrat_admissibility_hasneg(zeroadm));
            testCase.verifyNotEmpty(psyrat_admissibility_note(zeroadm));
            testCase.verifyWarningFree(@() psyrat_admissibility_warn(zeroadm));

            % Readers fail open on anything malformed or legacy.
            testCase.verifyEmpty(psyrat_admissibility_note([]));
            testCase.verifyFalse(psyrat_admissibility_hasneg(struct('a',1)));
        end
    end
end

function [bp, bt, er_var, wp_cov, obs] = localRC01Fixture()
%The out-of-range difference-score case documented in the findings ledger
%(RC-01). Residual SDs of 4 and 1 give residual variances of 16 and 1, and at
%obs = [10 1] the analytic threshold for a negative relative error variance is
%rho = (16/10 + 1)/(1.1*sqrt(16)) = 26/44 = 0.5909. rho = 0.70 clears it.
%Components are constant across draws so every quantity has a closed form.
nd = 6;
one = ones(nd,1);

bp = localCovmat(2*one, 2*one, 0.5*one);      % uni = 2 + 2 - 2*0.5 = 3
bt = localCovmat(1*one, 1*one, 0.9*one);
er_var = log([4*one 1*one]);                  % residual variances 16 and 1
wp_cov = 0.7*sqrt(16*1)*one;                  % rho_e = 0.70
obs = [10 1];
end

function localCallDepWithSmallCI()
psyrat_dep('bp', [1; 2; 3], 'wp', [1; 1; 1], 'obs', 3, 'CI', 0.40);
end

function localCallICCWithSmallCI()
psyrat_icc('bp', [1; 2; 3], 'wp', [1; 1; 1], 'CI', 0.40);
end

function [bp, bt, erVar] = localDiffRelInputs()
n = 5;
bp = zeros(n, 2, 2);
bp(:, 1, 1) = [0.50; 0.55; 0.60; 0.65; 0.70];
bp(:, 2, 2) = [0.45; 0.50; 0.55; 0.60; 0.65];
bp(:, 2, 1) = [0.10; 0.11; 0.10; 0.12; 0.11];
bp(:, 1, 2) = bp(:, 2, 1);

bt = zeros(n, 2, 2);
bt(:, 1, 1) = [0.25; 0.27; 0.26; 0.28; 0.30];
bt(:, 2, 2) = [0.20; 0.22; 0.24; 0.23; 0.25];
bt(:, 2, 1) = [0.04; 0.05; 0.04; 0.06; 0.05];
bt(:, 1, 2) = bt(:, 2, 1);

erVar = log([ ...
    0.75 0.82; ...
    0.78 0.80; ...
    0.79 0.81; ...
    0.77 0.83; ...
    0.76 0.84]);
end

function bp = localDoDCovDraws()
n = 5;
bp = zeros(n, 4, 4);
for d = 1:n
    m = [0.60 0.10 0.06 0.05; ...
        0.10 0.58 0.07 0.04; ...
        0.06 0.07 0.55 0.03; ...
        0.05 0.04 0.03 0.57] + eye(4) * (d - 1) * 0.01;
    bp(d, :, :) = m;
end
end

function [depDraws, depIcc, genDraws, genIcc] = localDoDManualRel(bp, bt, erVarLog, obs)
nd = size(bp, 1);
c = [1; -1; -1; 1];
erVar = exp(erVarLog).^2;

u = zeros(nd, 1);
btScaled = zeros(nd, 1);
btUnscaled = zeros(nd, 1);
erScaled = sum(erVar .* (1 ./ obs), 2);
erUnscaled = sum(erVar, 2);

for d = 1:nd
    sbp = squeeze(bp(d, :, :));
    sbt = squeeze(bt(d, :, :));
    u(d) = c' * sbp * c;
    btUnscaled(d) = c' * sbt * c;

    sbtScaled = zeros(4, 4);
    for i = 1:4
        for j = 1:4
            if i == j
                sbtScaled(i, j) = sbt(i, j) / obs(i);
            else
                nharm = 2 / ((1 / obs(i)) + (1 / obs(j)));
                sbtScaled(i, j) = sbt(i, j) / nharm;
            end
        end
    end
    btScaled(d) = c' * sbtScaled * c;
end

depDraws = u ./ (u + btScaled + erScaled);
genDraws = u ./ (u + erScaled);
depIcc = u ./ (u + btUnscaled + erUnscaled);
genIcc = u ./ (u + erUnscaled);
end

function m = localCovmat(v1,v2,cv)
n = length(v1);
m = zeros(n,2,2);
m(:,1,1) = v1(:);
m(:,2,2) = v2(:);
m(:,2,1) = cv(:);
m(:,1,2) = cv(:);
end

function rel = localTable6TrtRel(bp,bpi,bpo,bt,bo,boi,wp1,wp2,er_cov,obs,nocc,reltype,est)
hmi = psyrat_harmmean(obs);
hmo = psyrat_harmmean(nocc);

bp_raw = bp(:,1,1) + bp(:,2,2) - (2 .* bp(:,2,1));
bpi_raw = bpi(:,1,1) + bpi(:,2,2) - (2 .* bpi(:,2,1));
bpo_raw = bpo(:,1,1) + bpo(:,2,2) - (2 .* bpo(:,2,1));
bt_raw = bt(:,1,1) + bt(:,2,2) - (2 .* bt(:,2,1));
bo_raw = bo(:,1,1) + bo(:,2,2) - (2 .* bo(:,2,1));
boi_raw = boi(:,1,1) + boi(:,2,2) - (2 .* boi(:,2,1));
bpoi_raw = wp1 + wp2 - (2 .* er_cov);

bpi_adj = (bpi(:,1,1) ./ obs(1)) + (bpi(:,2,2) ./ obs(2)) - ((2 .* bpi(:,2,1)) ./ hmi);
bpo_adj = (bpo(:,1,1) ./ nocc(1)) + (bpo(:,2,2) ./ nocc(2)) - ((2 .* bpo(:,2,1)) ./ hmo);
bt_adj = (bt(:,1,1) ./ obs(1)) + (bt(:,2,2) ./ obs(2)) - ((2 .* bt(:,2,1)) ./ hmi);
bo_adj = (bo(:,1,1) ./ nocc(1)) + (bo(:,2,2) ./ nocc(2)) - ((2 .* bo(:,2,1)) ./ hmo);
boi_adj = (boi(:,1,1) ./ (obs(1)*nocc(1))) + (boi(:,2,2) ./ (obs(2)*nocc(2))) - ...
    ((2 .* boi(:,2,1)) ./ (hmi*hmo));
bpoi_adj = (wp1 ./ (obs(1)*nocc(1))) + (wp2 ./ (obs(2)*nocc(2))) - ...
    ((2 .* er_cov) ./ (hmi*hmo));

switch reltype
    case 1
        uni = bp_raw + bpo_adj;
        uni_raw = bp_raw + bpo_raw;
        rel_err = bpi_adj + bpoi_adj;
        rel_err_raw = bpi_raw + bpoi_raw;
        abs_err = rel_err + bt_adj + boi_adj;
        abs_err_raw = rel_err_raw + bt_raw + boi_raw;
    case 2
        uni = bp_raw + bpi_adj;
        uni_raw = bp_raw + bpi_raw;
        rel_err = bpo_adj + bpoi_adj;
        rel_err_raw = bpo_raw + bpoi_raw;
        abs_err = rel_err + bo_adj + boi_adj;
        abs_err_raw = rel_err_raw + bo_raw + boi_raw;
    case 3
        uni = bp_raw;
        uni_raw = bp_raw;
        rel_err = bpi_adj + bpo_adj + bpoi_adj;
        rel_err_raw = bpi_raw + bpo_raw + bpoi_raw;
        abs_err = rel_err + bt_adj + bo_adj + boi_adj;
        abs_err_raw = rel_err_raw + bt_raw + bo_raw + boi_raw;
end

switch lower(est)
    case 'dep'
        rel = uni ./ (uni + abs_err);
    case 'gen'
        rel = uni ./ (uni + rel_err);
end

if isequal(obs,[1 1]) && isequal(nocc,[1 1])
    switch lower(est)
        case 'dep'
            rel = uni_raw ./ (uni_raw + abs_err_raw);
        case 'gen'
            rel = uni_raw ./ (uni_raw + rel_err_raw);
    end
end
end
