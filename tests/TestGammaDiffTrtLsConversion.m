classdef TestGammaDiffTrtLsConversion < matlab.unittest.TestCase
    %TESTGAMMADIFFTRTLSCONVERSION Offline checks for the LOCATION-SCALE gamma
    %TWO-FACET (test-retest) two-event DIFFERENCE conversion (analysis 10,
    %gammascale = 2, non-concurrent).
    %
    %WHY THIS FILE EXISTS. psyrat_gamma_extract_diff_trt_ls is a LOCAL function of
    %psyrat_computevarcomp: it cannot be called from a test and it only runs
    %against a live CmdStan fit. This class therefore pins the algebra it rests
    %on, through the callable helpers it delegates to
    %(psyrat_gamma_varcomps_trt_ls, psyrat_gamma_crosscov_trt) plus the REAL
    %downstream consumer (psyrat_diffrel_trt).
    %
    %WHAT IS DIFFERENT FROM THE ONE-FACET SIBLING (TestGammaDiffLsConversion).
    %Analysis 10 is GROUP level and has no subject-level twin, so there are no
    %per-subject encoding tests here. What replaces them is the two-facet
    %question: the log-nu residual is an expectation over the five within-person
    %facets and therefore carries exp(2*V); the location-scale residual is a
    %single per-participant sigma with no facet terms and carries none. That
    %asymmetry is the reason a psyrat_gamma_varcomps_trt_ls had to be written at
    %all, and it is the single most likely place for a future edit to go wrong
    %"for symmetry".
    %
    %A green run here does NOT establish that estimation recovers truth; that is
    %the live-recovery lane's job (TestGammaDiffTrtLsRecovery). It establishes
    %that the algebra the extractor rests on is what this repo believes it is.
    %
    %See documentation/gamma_scale_submodel_decision.md section H and finding S16.

    % Copyright (C) 2016-2025 Peter E. Clayson
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program (gpl.txt). If not, see
    %     <http://www.gnu.org/licenses/>.

    methods (Test)

        function testSixGroupComponentsAreParameterizationInvariant(testCase)
            % THE LOAD-BEARING ASSUMPTION of psyrat_gamma_extract_diff_trt_ls. It
            % builds all six crossed varcov arrays from per-event
            % psyrat_gamma_varcomps_trt_ls diagonals plus psyrat_gamma_crosscov_trt
            % off-diagonals, exactly as its log-nu twin does with
            % psyrat_gamma_varcomps_trt. If the per-event diagonals disagreed
            % between the two converters, the location-scale difference coefficient
            % would differ from the log-nu one for a reason having nothing to do
            % with the residual - and the "only the residual changes" claim in the
            % builder header, which the golden diff appears to confirm, would be
            % false at the conversion layer.
            %
            % Two events with DIFFERENT amplitudes and different scale parameters,
            % because a symmetric setup could hide an amplitude-dependent error.
            c = localCells();

            for e = 1:2
                vcLs = psyrat_gamma_varcomps_trt_ls(c.alpha(e), c.s_p(e), ...
                    c.s_o(e), c.s_t(e), c.s_po(e), c.s_pt(e), c.s_ot(e), ...
                    c.logsig(e), c.s_sd(e));
                % Deliberately UNRELATED log-nu inputs: the point is that whatever
                % the dispersion side is doing, the mean-surface fields do not move.
                vcNu = psyrat_gamma_varcomps_trt(c.alpha(e), c.s_p(e), ...
                    c.s_o(e), c.s_t(e), c.s_po(e), c.s_pt(e), c.s_ot(e), ...
                    c.nu(e), c.s_v(e), c.covpv(e));

                for f = {'sigma_p2','sigma_o2','sigma_t2', ...
                        'sigma_po2','sigma_pt2','sigma_ot2','var_mu','mu_bar'}
                    testCase.verifyEqual(vcLs.(f{1}), vcNu.(f{1}), 'AbsTol', 0, ...
                        sprintf(['Event %d: %s is built from the log-mean submodel '...
                        'alone and must not depend on the scale parameterization.'], ...
                        e, f{1}));
                end

                % ... and the residual MUST differ, or every assertion above is
                % vacuous (they would pass just as well if the two converters were
                % the same function). This is the non-vacuity guard.
                testCase.verifyNotEqual(vcLs.sigma_res2, vcNu.sigma_res2, ...
                    sprintf(['Event %d: the two parameterizations must give '...
                    'DIFFERENT residuals; if they agree the invariance '...
                    'assertions above prove nothing.'], e));
            end
        end

        function testAssembledCrossCellMatricesAreParameterizationInvariant(testCase)
            % psyrat_gamma_crosscov_trt is called with IDENTICAL arguments on both
            % paths - it takes no dispersion input at all - so the off-diagonals of
            % all six varcov arrays are parameterization-invariant by construction.
            % This pins that the fully ASSEMBLED matrices agree, not just the
            % diagonals, which is what the extractor header claims and what
            % justifies reusing psyrat_gamma_crosscov_trt unchanged.
            %
            % The non-concurrent policy is applied here exactly as the extractors
            % apply it: person/occasion/person-x-occasion keep their cross-cell
            % covariance; the three trial-indexed arrays are zeroed.
            c = localCells();

            [lsF, lsS] = localVarcompPair(c, true);
            [nuF, nuS] = localVarcompPair(c, false);

            cc = psyrat_gamma_crosscov_trt(c.alpha(1), c.s_p(1), c.s_o(1), ...
                c.s_t(1), c.s_po(1), c.s_pt(1), c.s_ot(1), ...
                c.alpha(2), c.s_p(2), c.s_o(2), c.s_t(2), c.s_po(2), ...
                c.s_pt(2), c.s_ot(2), ...
                c.cor_p, c.cor_o, c.cor_t, c.cor_po, c.cor_pt, c.cor_ot);

            % {component field, cross-cell covariance used by the non-concurrent
            % policy}. Zero means the trial-indexed cross-cov is dropped.
            spec = { ...
                'sigma_p2',  cc.cov_p; ...
                'sigma_o2',  cc.cov_o; ...
                'sigma_po2', cc.cov_po; ...
                'sigma_t2',  0; ...
                'sigma_pt2', 0; ...
                'sigma_ot2', 0};

            for k = 1:size(spec,1)
                f = spec{k,1};
                mLs = localVarcov(lsF.(f), lsS.(f), spec{k,2});
                mNu = localVarcov(nuF.(f), nuS.(f), spec{k,2});
                testCase.verifyEqual(mLs, mNu, 'AbsTol', 0, ...
                    sprintf(['The assembled %s cross-cell matrix must be '...
                    'parameterization-invariant.'], f));
            end
        end

        function testResidualCarriesNoWithinPersonFacetTerm(testCase)
            % THE SCIENTIFIC PIN of this increment. It is the two-facet analogue of
            % the one-facet no-trial-averaging-factor pin, and the reason
            % psyrat_gamma_varcomps_trt_ls exists rather than reusing the log-nu
            % converter with substituted arguments.
            %
            % The log-nu conditional variance 2*mu^2/nu is MEAN-COUPLED, so its
            % expectation over the design picks up exp(2*V) where V sums ALL six
            % log-mean variances: e_cond_var = (2/nu)*exp(2*alpha + 2*V)*F. Under
            % location-scale sigma_p is a single per-participant quantity with no
            % facet terms, so e_cond_var = exp(2*log_sigma + 2*s_sd^2): no alpha,
            % no V, no cov_pv.
            %
            % Asserted as a SENSITIVITY, which is the property that would actually
            % break: moving the five within-person facet SDs must leave the
            % location-scale residual term untouched, while it demonstrably moves
            % the log-nu one. A future edit that reintroduces an exp(2*V) factor
            % fails the first half; an edit that accidentally made the log-nu path
            % facet-free would fail the second, so neither half can pass vacuously.
            c = localCells();
            e = 1;

            base = psyrat_gamma_varcomps_trt_ls(c.alpha(e), c.s_p(e), c.s_o(e), ...
                c.s_t(e), c.s_po(e), c.s_pt(e), c.s_ot(e), c.logsig(e), c.s_sd(e));
            % Same person SD and same scale submodel; every OTHER facet doubled.
            moved = psyrat_gamma_varcomps_trt_ls(c.alpha(e), c.s_p(e), ...
                2*c.s_o(e), 2*c.s_t(e), 2*c.s_po(e), 2*c.s_pt(e), 2*c.s_ot(e), ...
                c.logsig(e), c.s_sd(e));

            testCase.verifyEqual(moved.e_cond_var, base.e_cond_var, 'AbsTol', 0, ...
                ['The location-scale observation variance must not depend on any '...
                'facet SD. If this fails, an exp(2*V) mean-coupling factor has '...
                'been reintroduced, which is exactly what this parameterization '...
                'exists to remove.']);

            % It is also exactly the closed form, not merely facet-free.
            testCase.verifyEqual(base.e_cond_var, ...
                exp(2*c.logsig(e) + 2*c.s_sd(e)^2), 'RelTol', 1e-12, ...
                'e_cond_var must be exp(2*log_sigma + 2*s_sd^2).');

            % Non-vacuity: the log-nu twin DOES move under the same perturbation.
            baseNu = psyrat_gamma_varcomps_trt(c.alpha(e), c.s_p(e), c.s_o(e), ...
                c.s_t(e), c.s_po(e), c.s_pt(e), c.s_ot(e), c.nu(e), c.s_v(e), c.covpv(e));
            movedNu = psyrat_gamma_varcomps_trt(c.alpha(e), c.s_p(e), ...
                2*c.s_o(e), 2*c.s_t(e), 2*c.s_po(e), 2*c.s_pt(e), 2*c.s_ot(e), ...
                c.nu(e), c.s_v(e), c.covpv(e));
            testCase.verifyGreaterThan(movedNu.e_cond_var, baseNu.e_cond_var, ...
                ['The log-nu observation variance MUST move with the facet SDs; '...
                'if it does not, the contrast this test draws is meaningless.']);
        end

        function testReducesToTheOneFacetLocationScaleConverter(testCase)
            % The log-nu pair satisfies psyrat_gamma_varcomps_trt -> ...varcomps
            % when the occasion and interaction SDs vanish. The location-scale pair
            % must satisfy the same reduction, or the two-facet converter has
            % introduced an algebra difference that is not the scale submodel.
            %
            % THE RESIDUALS ARE NOT THE SAME QUANTITY, and the reduction has to say
            % so. The one-facet sigma_pi_e2 is "P x trial interaction + observation
            % dispersion" - the one-facet design cannot separate them. The two-facet
            % sigma_res2 is "P x O x T + dispersion", with P x trial broken out as
            % its own sigma_pt2. So the identity is
            %   sigma_res2 + sigma_pt2 = sigma_pi_e2,
            % NOT sigma_res2 = sigma_pi_e2. Note sigma_pt2 is NONZERO here even
            % though s_pt = 0: by inclusion-exclusion the interaction component is
            % B*(exp(vp)-1)*(exp(vt)-1) on the observed scale, which the log link
            % makes positive whenever both main effects are. Asserting the naive
            % equality instead would fail by exactly that component (measured 0.187
            % at this design point, ~3% of the residual), which is a real quantity
            % and not a tolerance problem.
            alpha = 1.70; s_p = 0.35; s_t = 0.20; logsig = 0.61; s_sd = 0.50;

            trt = psyrat_gamma_varcomps_trt_ls(alpha, s_p, 0, s_t, 0, 0, 0, ...
                logsig, s_sd);
            one = psyrat_gamma_varcomps_ls(alpha, s_p, s_t, logsig, s_sd);

            testCase.verifyEqual(trt.mu_bar, one.mu_bar, 'RelTol', 1e-12, ...
                'Reduction: mu_bar must match the one-facet converter.');
            testCase.verifyEqual(trt.sigma_p2, one.sigma_p2, 'RelTol', 1e-12, ...
                'Reduction: person component must match the one-facet converter.');
            testCase.verifyEqual(trt.sigma_t2, one.sigma_i2, 'RelTol', 1e-12, ...
                'Reduction: trial component must match the one-facet converter.');
            testCase.verifyEqual(trt.var_mu, one.var_mu, 'RelTol', 1e-12, ...
                'Reduction: mean-surface variance must match.');
            testCase.verifyEqual(trt.e_cond_var, one.e_cond_var, 'AbsTol', 0, ...
                'Reduction: the scale submodel is identical in both converters.');
            testCase.verifyEqual(trt.sigma_res2 + trt.sigma_pt2, one.sigma_pi_e2, ...
                'RelTol', 1e-12, ...
                ['Reduction: the two-facet residual PLUS the person x trial '...
                'component must equal the one-facet residual, which confounds '...
                'the two.']);

            % The confounded component must be materially nonzero, or the identity
            % above is indistinguishable from the naive one and proves less.
            testCase.verifyGreaterThan(trt.sigma_pt2, 0.01 * one.sigma_pi_e2, ...
                ['The P x trial component must be a real share of the one-facet '...
                'residual, or this test would pass under the wrong identity too.']);

            % The genuinely occasion-indexed components must vanish, so the
            % reduction is not passing because everything collapsed to a degenerate
            % point. sigma_pt2 is deliberately NOT in this list.
            for f = {'sigma_o2','sigma_po2','sigma_ot2'}
                testCase.verifyEqual(trt.(f{1}), 0, 'AbsTol', 1e-12, ...
                    sprintf('Reduction: %s must vanish with no occasion facet.', f{1}));
            end
            testCase.verifyGreaterThan(trt.sigma_p2, 0, ...
                'Reduction must be non-degenerate: the person component is positive.');

            % The SAME identity must hold for the log-nu pair, which shows it is a
            % property of the observed-scale decomposition rather than something
            % the location-scale branch introduced.
            trtNu = psyrat_gamma_varcomps_trt(alpha, s_p, 0, s_t, 0, 0, 0, ...
                exp(2.89), 0.70, 0.05);
            oneNu = psyrat_gamma_varcomps(alpha, s_p, s_t, exp(2.89), 0.70, 0.05);
            testCase.verifyEqual(trtNu.sigma_res2 + trtNu.sigma_pt2, ...
                oneNu.sigma_pi_e2, 'RelTol', 1e-12, ...
                ['The same reduction identity must hold under log-nu; it is a '...
                'property of the decomposition, not of this increment.']);

            % s_pt = 0 IS REQUIRED and is not implied by zeroing the three
            % occasion terms. The one-facet model carries no log-scale person x
            % trial effect at all, so a nonzero s_pt inflates V and has no
            % counterpart on the other side of the identity. Pinned as a NEGATIVE
            % case, because the reduction assertions above set s_pt = 0 and would
            % pass just as well under a header that forgot to say so - which is
            % exactly the wording error this pins against.
            %
            % e_cond_var still agrees, because the scale submodel never sees any
            % facet; it is only the mean-surface decomposition that diverges.
            sptBreak = 0.25;
            trtBreak = psyrat_gamma_varcomps_trt_ls(alpha, s_p, 0, s_t, 0, ...
                sptBreak, 0, logsig, s_sd);
            testCase.verifyEqual(trtBreak.e_cond_var, one.e_cond_var, 'AbsTol', 0, ...
                ['The scale submodel must stay identical even when the mean '...
                'submodel gains a person x trial effect.']);
            relErr = (trtBreak.sigma_res2 + trtBreak.sigma_pt2) / one.sigma_pi_e2 - 1;
            testCase.verifyGreaterThan(relErr, 0.40, sprintf( ...
                ['With s_pt = %.2f the reduction identity must FAIL, and visibly '...
                '(measured %+.1f%%). If this passes, either s_pt has stopped '...
                'entering the mean surface or the identity has been restated '...
                'too loosely.'], sptBreak, 100*relErr));
        end

        function testResidualEncodingRoundTripsThroughTheConsumer(testCase)
            % psyrat_diffrel_trt reads er_var as a log-scale SD [draws x 2] and
            % squares it. The extractor therefore stores
            % b_sigma = 0.5*log(sigma_res2), the same encoding the log-nu twin uses.
            % This pins the round trip and then drives the REAL consumer with
            % location-scale-shaped inputs, checking the difference-score
            % generalizability against a hand reference - so it tests the assembled
            % contract rather than restating the formula.
            c = localCells();
            nd = 40;
            rng(20260804, 'twister');   % local, deterministic; no hidden state

            [vcF, vcS] = localVarcompPair(c, true);

            % Per-draw jitter so the consumer sees genuine draw vectors.
            resF = vcF.sigma_res2 .* exp(0.05*randn(nd,1));
            resS = vcS.sigma_res2 .* exp(0.05*randn(nd,1));
            er_var = [0.5*log(resF), 0.5*log(resS)];

            testCase.verifyEqual(exp(er_var(:,1)).^2, resF, 'RelTol', 1e-12, ...
                'The 0.5*log encoding must round-trip to the residual variance.');
            testCase.verifyEqual(exp(er_var(:,2)).^2, resS, 'RelTol', 1e-12, ...
                'The 0.5*log encoding must round-trip to the residual variance.');

            % Observed-scale crossed components, shaped as the extractor assembles
            % them. Person and occasion families carry a cross-cell covariance; the
            % three trial-indexed families are zeroed (non-concurrent policy).
            bp  = localVarcov(4.0*ones(nd,1), 5.0*ones(nd,1), -0.6*ones(nd,1));
            bpo = localVarcov(0.30*ones(nd,1), 0.22*ones(nd,1), 0.05*ones(nd,1));
            bo  = localVarcov(0.40*ones(nd,1), 0.35*ones(nd,1), 0.08*ones(nd,1));
            bpi = localVarcov(0.25*ones(nd,1), 0.18*ones(nd,1), zeros(nd,1));
            bt  = localVarcov(0.15*ones(nd,1), 0.12*ones(nd,1), zeros(nd,1));
            boi = localVarcov(0.10*ones(nd,1), 0.09*ones(nd,1), zeros(nd,1));

            % Deliberately UNEQUAL trial counts, so the harmonic-mean covariance
            % divisor psyrat_diffrel_trt applies is actually exercised rather than
            % collapsing to the arithmetic count.
            obs = [20 25]; nocc = [2 2];
            out = psyrat_diffrel_trt('bp',bp,'bpi',bpi,'bpo',bpo,'bt',bt, ...
                'bo',bo,'boi',boi,'er_var',er_var,'obs',obs,'nocc',nocc, ...
                'reltype',3,'est','gen','CI',0.95,'er_cov',0);

            % Hand reference for the coefficient of equivalence and stability with a
            % structurally zero residual covariance. Each D-study adjusted block
            % divides its two variances by their OWN counts and the cross-condition
            % covariance by the harmonic mean of those counts - the Rocha Table 6
            % scaler (see psyrat_harmmean and the RC-01/S13 note in
            % psyrat_diffrel_trt; that finding is not this increment's business,
            % but the reference must match the formula as implemented).
            hmi = psyrat_harmmean(obs); hmo = psyrat_harmmean(nocc);
            varP = bp(:,1,1) + bp(:,2,2) - 2*bp(:,1,2);
            relErr = (bpi(:,1,1)./obs(1) + bpi(:,2,2)./obs(2) - 2*bpi(:,1,2)./hmi) ...
                + (bpo(:,1,1)./nocc(1) + bpo(:,2,2)./nocc(2) - 2*bpo(:,1,2)./hmo) ...
                + (resF./(obs(1)*nocc(1)) + resS./(obs(2)*nocc(2)));
            Gdraw = varP ./ (varP + relErr);

            testCase.verifyEqual(out.pt, mean(Gdraw), 'RelTol', 1e-8, ...
                ['The non-concurrent two-facet difference residual must enter as '...
                'sigma^2_res,1/(n_t1*n_o1) + sigma^2_res,2/(n_t2*n_o2) with zero '...
                'covariance.']);
        end

        function testNonconcurrentResidualCovarianceIsStructurallyZero(testCase)
            % Owner ruling H: for the non-concurrent difference the residual
            % covariance is structurally zero, NOT merely assumed. Correlating the
            % two events' person SCALE effects correlates the residual VARIANCES,
            % not the residuals, since E[e1*e2] = E[E[e1|p]*E[e2|p]] = 0.
            %
            % This is the configuration that would break if someone "fixed" the zero
            % covariance by deriving one from the fitted person scale correlation:
            % strongly correlated scale effects, and a coefficient that must be
            % identical to the uncorrelated case because er_cov stays 0.
            nd = 30;
            rng(20260804, 'twister');

            w1 = 0.5*randn(nd,1);
            wLow  = 0.5*randn(nd,1);                          % rho ~ 0
            wHigh = 0.9*w1 + sqrt(1-0.9^2)*0.5*randn(nd,1);   % rho = 0.9

            bp  = localVarcov(4.0*ones(nd,1), 5.0*ones(nd,1), -0.6*ones(nd,1));
            bpo = localVarcov(0.30*ones(nd,1), 0.22*ones(nd,1), 0.05*ones(nd,1));
            bo  = localVarcov(0.40*ones(nd,1), 0.35*ones(nd,1), 0.08*ones(nd,1));
            bpi = localVarcov(0.25*ones(nd,1), 0.18*ones(nd,1), zeros(nd,1));
            bt  = localVarcov(0.15*ones(nd,1), 0.12*ones(nd,1), zeros(nd,1));
            boi = localVarcov(0.10*ones(nd,1), 0.09*ones(nd,1), zeros(nd,1));

            args = {'bp',bp,'bpi',bpi,'bpo',bpo,'bt',bt,'bo',bo,'boi',boi, ...
                'obs',[20 25],'nocc',[2 2],'reltype',3,'est','gen','CI',0.95, ...
                'er_cov',0};

            outLow  = psyrat_diffrel_trt(args{:}, 'er_var', [0.61 + w1, 0.48 + wLow]);
            outHigh = psyrat_diffrel_trt(args{:}, 'er_var', [0.61 + w1, 0.48 + wHigh]);

            % The two runs have different residual DRAWS, so the coefficients differ;
            % what must hold is that neither run derived a residual covariance from
            % the scale correlation. Verified by reconstructing each coefficient
            % with er_cov = 0 and requiring an exact match.
            hmi = psyrat_harmmean([20 25]); hmo = psyrat_harmmean([2 2]);
            for pair = {{outLow, wLow}, {outHigh, wHigh}}
                o = pair{1}{1}; w2 = pair{1}{2};
                res1 = exp(0.61 + w1).^2;
                res2 = exp(0.48 + w2).^2;
                varP = bp(:,1,1) + bp(:,2,2) - 2*bp(:,1,2);
                relErr = (bpi(:,1,1)./20 + bpi(:,2,2)./25 - 2*bpi(:,1,2)./hmi) ...
                    + (bpo(:,1,1)./2 + bpo(:,2,2)./2 - 2*bpo(:,1,2)./hmo) ...
                    + res1./(20*2) + res2./(25*2);
                testCase.verifyEqual(o.pt, mean(varP ./ (varP + relErr)), ...
                    'RelTol', 1e-8, ...
                    ['The difference residual must be the SUM of the two event '...
                    'residuals regardless of how strongly the person scale '...
                    'effects are correlated.']);
            end

            % Non-vacuity: the two runs must actually differ, or the loop above
            % would be checking the same arithmetic twice.
            testCase.verifyNotEqual(outLow.pt, outHigh.pt, ...
                ['The two scale-correlation settings must produce different '...
                'residual draws, or this test compares nothing.']);
        end

    end
end

function c = localCells()
%Shared two-cell design point. Deliberately ASYMMETRIC: event 2 has the lower
%amplitude and the higher residual SD, so both terms push the implied log-nu
%contrast the same way and the amplitude/precision blend is maximal. Same
%principle as the increment-5 recovery design point.
c = struct();
c.alpha  = [1.70; 1.40];    % per-event log-mean intercept
c.s_p    = [0.35; 0.42];    % person
c.s_o    = [0.18; 0.15];    % occasion
c.s_t    = [0.20; 0.11];    % trial
c.s_po   = [0.22; 0.19];    % person x occasion
c.s_pt   = [0.16; 0.13];    % person x trial
c.s_ot   = [0.09; 0.07];    % occasion x trial
% Location-scale scale submodel
c.logsig = [0.61; 0.75];
c.s_sd   = [0.50; 0.33];
% Deliberately UNRELATED log-nu inputs
c.nu     = exp([2.89; 2.70]);
c.s_v    = [0.70; 0.40];
c.covpv  = [0.05; -0.02];
% Cross-event correlations (mean submodel only, shared by both paths)
c.cor_p = 0.55; c.cor_o = 0.30; c.cor_t = 0.20;
c.cor_po = 0.25; c.cor_pt = 0.15; c.cor_ot = 0.10;
end

function [vcF, vcS] = localVarcompPair(c, isLs)
%Per-cell components for both events, under the requested parameterization.
if isLs
    vcF = psyrat_gamma_varcomps_trt_ls(c.alpha(1), c.s_p(1), c.s_o(1), ...
        c.s_t(1), c.s_po(1), c.s_pt(1), c.s_ot(1), c.logsig(1), c.s_sd(1));
    vcS = psyrat_gamma_varcomps_trt_ls(c.alpha(2), c.s_p(2), c.s_o(2), ...
        c.s_t(2), c.s_po(2), c.s_pt(2), c.s_ot(2), c.logsig(2), c.s_sd(2));
else
    vcF = psyrat_gamma_varcomps_trt(c.alpha(1), c.s_p(1), c.s_o(1), ...
        c.s_t(1), c.s_po(1), c.s_pt(1), c.s_ot(1), c.nu(1), c.s_v(1), c.covpv(1));
    vcS = psyrat_gamma_varcomps_trt(c.alpha(2), c.s_p(2), c.s_o(2), ...
        c.s_t(2), c.s_po(2), c.s_pt(2), c.s_ot(2), c.nu(2), c.s_v(2), c.covpv(2));
end
end

function m = localVarcov(v1, v2, cv)
%Assemble the [draws x 2 x 2] observed-scale var-cov exactly as the extractors
%do, so comparisons are of the assembled matrix rather than of loose components.
n = max([numel(v1) numel(v2) numel(cv)]);
m = zeros(n,2,2);
m(:,1,1) = v1(:) .* ones(n,1);
m(:,2,2) = v2(:) .* ones(n,1);
m(:,1,2) = cv(:) .* ones(n,1);
m(:,2,1) = cv(:) .* ones(n,1);
end
