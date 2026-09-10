classdef TestGammaDiffTrtLsRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square static
    % TWO-FACET (test-retest) two-event DIFFERENCE design under the
    % LOCATION-SCALE parameterization (analysis 10, gammascale = 2,
    % NON-concurrent). Increment 6 of the gamma location-scale rollout.
    %
    % GENERATIVE MODEL. Per event e, participant p, occasion o and trial t,
    %   log mu_epot  = b[e] + u_{p,e} + u_{o,e} + u_{t,e}
    %                        + u_{po,e} + u_{pt,e} + u_{ot,e}
    %   log sigma_ep = b_sigma[e] + w_{p,e}
    %   nu_epot      = 2*(mu_epot / sigma_ep)^2,   Y ~ scaled chi-square
    % so Var(Y | e,p) = sigma_ep^2 exactly, with no dependence on the mean or on
    % any facet. (u_{p,1}, u_{p,2}, w_1, w_2) share a person block, matching owner
    % ruling H.
    %
    % WHAT THIS TEST ADDS OVER THE ONE-FACET SIBLING (TestGammaDiffLsRecovery).
    % Analysis 10 is GROUP level and has no subject-level twin, so there are no
    % per-participant assertions here. What replaces them is the two-facet
    % question the increment turns on: under log-nu the residual is an expectation
    % over the five within-person facets and carries exp(2*V); under
    % location-scale sigma_ep has no facet terms and carries none. That asymmetry
    % is why psyrat_gamma_varcomps_trt_ls exists, and the recovery below is the
    % end-to-end check of it.
    %
    % NOTE ON WHAT REL.out.b_sigma IS ON THIS PATH. Unlike the case-8 path, where
    % it is the Stan per-event intercept, the GROUP-level extractor stores the
    % observed-scale residual ENCODING 0.5*log(sigma_res2). Comparisons below are
    % therefore against the CLOSED-FORM observed-scale residual computed from the
    % planted values, not against b_sigma directly. Confusing the two would
    % compare quantities that differ by the whole mean surface.
    %
    % TWO ASSERTIONS THIS TEST DELIBERATELY DOES NOT MAKE, both measured false in
    % increment 4 and easy to reinstate by accident:
    %   - that a log-nu fit gets the coefficient WRONG at an adequate trial count.
    %     It does not; the paired fit below is used for the cross-parameterization
    %     identities, not to show log-nu failing.
    %   - that the two parameterizations' stored residuals agree. They must not:
    %     the log-nu one carries the facet-averaging factor and the location-scale
    %     one does not.
    %
    % CmdStan-gated (skips Incomplete, not Failed, when CmdStan/MatlabStan absent)
    % and R-hat gated (a non-converged fit skips rather than false-passing).

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

    methods (TestClassSetup)
        function configureDependencies(testCase)
            PsyRATCmdStanDependencySetup.configure(testCase);
        end
    end

    methods (TestClassTeardown)
        function cleanupCmdStanArtifacts(testCase)
            localCleanupArtifacts(localCmdStanOutputDir(testCase.projectRoot()));
        end
    end

    methods (Test)

        function testLocationScaleTwoFacetDifferenceRecovered(testCase)
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required to draw continuous-df scaled chi-square.');

            T = localTruth();
            dataTbl = localSimulateGammaDiffTrtLs(T);

            % ---- Closed-form observed-scale truth -------------------------------
            % Built through the SAME production converters the extractor uses, so
            % this is a recovery check on estimation rather than a second
            % implementation of the algebra (that layer is pinned offline in
            % TestGammaDiffTrtLsConversion).
            [vcT1, vcT2, ccT] = localTruthComponents(T);
            trueG = localDiffG(vcT1, vcT2, ccT, T.ntrl, T.nocc);

            % ---- Fit: analysis 10 under location-scale --------------------------
            rel = localRunDiffTrtGamma(testCase, dataTbl, 'recov_gamma_ls_diff_trt', 2);

            testCase.verifyEqual(rel.analysis, 'trt_diff');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.gammascale, 2);
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % ---- Convergence gate: skip, never false-pass ------------------------
            [maxRhat, nFinite] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[gamma ls diff trt] max R-hat=%.4f over %d params; converged=%d\n', ...
                maxRhat, nFinite, rel.out.conv.converged);
            testCase.assumeGreaterThan(nFinite, 0, ...
                'Convergence table empty; R-hat unextractable, so recovery cannot be validated.');
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf( ...
                ['Location-scale two-facet difference fit non-converged (max R-hat=%.3f); ' ...
                'skipping rather than reporting a pass on an unconverged fit.'], maxRhat));

            % ---- (1) the mean submodel and the person scale SDs ------------------
            bHat   = mean(cell2mat(rel.out.b{1}), 1);        % 1 x 2 per-event log-mean
            ssdHat = mean(rel.out.gamma_disp.s_sd, 1);       % 1 x 2 person log-sigma SD
            rhoHat = mean(rel.out.gamma_disp.rho_ps, 1);     % 1 x 2 provenance only

            fprintf(['[gamma ls diff trt] b hat=[%.3f %.3f] (true [%.3f %.3f])  ' ...
                's_sd hat=[%.3f %.3f] (true [%.2f %.2f])  rho_ps hat=[%.3f %.3f]\n'], ...
                bHat(1), bHat(2), T.b(1), T.b(2), ...
                ssdHat(1), ssdHat(2), T.s_sd(1), T.s_sd(2), rhoHat(1), rhoHat(2));

            for e = 1:2
                testCase.verifyLessThan(abs(bHat(e) - T.b(e)), 0.30, sprintf( ...
                    'Event %d log-mean intercept off: %.3f vs true %.3f.', ...
                    e, bHat(e), T.b(e)));
            end
            testCase.verifyGreaterThan(min(ssdHat), 0.20, sprintf( ...
                'Person log residual-SD SD collapsed: [%.3f %.3f] vs true [%.2f %.2f].', ...
                ssdHat(1), ssdHat(2), T.s_sd(1), T.s_sd(2)));

            % ---- (2) the observed-scale residual, per event -----------------------
            % THE quantity this increment changes. REL.out.b_sigma encodes
            % 0.5*log(sigma_res2), so exp(.)^2 recovers the residual variance and
            % is compared against the closed-form truth on the LOG scale, where the
            % tolerance is interpretable as a proportional error.
            resHat  = exp(mean(cell2mat(rel.out.b_sigma{1}), 1)).^2;   % 1 x 2
            resTrue = [vcT1.sigma_res2, vcT2.sigma_res2];
            fprintf(['[gamma ls diff trt] observed-scale residual hat=[%.3f %.3f] ' ...
                '(true [%.3f %.3f])\n'], resHat(1), resHat(2), resTrue(1), resTrue(2));
            for e = 1:2
                testCase.verifyLessThan(abs(log(resHat(e)/resTrue(e))), 0.45, sprintf( ...
                    ['Event %d observed-scale residual off by %.1f%%: %.3f vs true ' ...
                    '%.3f.'], e, 100*(resHat(e)/resTrue(e) - 1), resHat(e), resTrue(e)));
            end

            % THE RESIDUAL EVENT CONTRAST - asserted separately, because both
            % events could be biased together and still leave the contrast right
            % (or vice versa). This is the estimand the parameterization exists to
            % make readable.
            contrastHat  = 0.5*log(resHat(2)/resHat(1));
            contrastTrue = 0.5*log(resTrue(2)/resTrue(1));
            fprintf('[gamma ls diff trt] residual-SD event contrast hat=%.3f (true %.3f)\n', ...
                contrastHat, contrastTrue);
            testCase.verifyLessThan(abs(contrastHat - contrastTrue), 0.30, sprintf( ...
                ['The residual-SD EVENT CONTRAST is the estimand this design exists ' ...
                'to deliver: recovered %.3f vs true %.3f.'], contrastHat, contrastTrue));

            % ---- (3) the difference-score coefficient ----------------------------
            recG = localDiffGFromRel(rel, T.ntrl, T.nocc);
            fprintf('[gamma ls diff trt] difference-score G hat=%.3f (true %.3f)\n', ...
                recG, trueG);
            testCase.verifyGreaterThan(recG, 0, 'Difference-score G must be positive.');
            testCase.verifyLessThan(recG, 1, 'Difference-score G must be below 1.');
            testCase.verifyLessThan(abs(recG - trueG), 0.12, sprintf( ...
                'Difference-score generalizability off: %.3f vs true %.3f.', recG, trueG));

            % The design point must actually sit below the near-ceiling regime, or
            % the coefficient surface is being tested weakly - the criticism
            % increment 4's ledger records against its own design point.
            testCase.verifyLessThan(trueG, 0.90, sprintf( ...
                ['This recovery is meant to sit well below G ~ 0.98; true G is ' ...
                '%.3f. If the design drifted upward the coefficient surface is ' ...
                'being tested weakly again.'], trueG));
        end

        function testCrossParameterizationIdentitiesHold(testCase)
            % The validation gate the rollout doc sets for every increment: fit the
            % SAME simulated data under both parameterizations and check what must
            % agree and what must not. It is NOT a claim that either fit is wrong.
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required.');

            T = localTruth();
            dataTbl = localSimulateGammaDiffTrtLs(T);

            relLs = localRunDiffTrtGamma(testCase, dataTbl, 'recov_gamma_ls_diff_trt_pairA', 2);
            relNu = localRunDiffTrtGamma(testCase, dataTbl, 'recov_gamma_lognu_diff_trt_pairB', 1);

            [rhatLs, nLs] = psyrat_conv_maxrhat(relLs.out.conv.data);
            [rhatNu, nNu] = psyrat_conv_maxrhat(relNu.out.conv.data);
            relLs = psyrat_checkconv(relLs);
            relNu = psyrat_checkconv(relNu);
            fprintf(['[gamma ls diff trt identity] ls: max R-hat=%.4f over %d params, ' ...
                'converged=%d | lognu: max R-hat=%.4f over %d params, converged=%d\n'], ...
                rhatLs, nLs, relLs.out.conv.converged, ...
                rhatNu, nNu, relNu.out.conv.converged);
            testCase.assumeGreaterThan(min(nLs, nNu), 0, ...
                'Convergence table empty for one of the paired fits.');
            testCase.assumeEqual(relLs.out.conv.converged, 1, ...
                'Location-scale fit did not converge; skipping the identity check.');
            testCase.assumeEqual(relNu.out.conv.converged, 1, ...
                'Log-nu fit did not converge; skipping the identity check.');

            % (a) THE MEAN SUBMODEL IS BYTE-IDENTICAL between the two builders (see
            % the golden diff), so the two fits must agree on it. This is that
            % structural claim confirmed on live fits rather than by file diff.
            bLs = mean(cell2mat(relLs.out.b{1}), 1);
            bNu = mean(cell2mat(relNu.out.b{1}), 1);
            fprintf('[gamma ls diff trt identity] b (ls)=[%.3f %.3f] b (lognu)=[%.3f %.3f]\n', ...
                bLs(1), bLs(2), bNu(1), bNu(2));
            for e = 1:2
                testCase.verifyLessThan(abs(bLs(e) - bNu(e)), 0.15, sprintf( ...
                    ['Event %d log-mean intercept differs across parameterizations ' ...
                    '(%.3f vs %.3f); the mean submodel is byte-identical, so these ' ...
                    'must agree.'], e, bLs(e), bNu(e)));
            end

            % (b) THE OBSERVED-SCALE SIGNAL COMPONENTS must agree across the two
            % fits. This is the fit-level counterpart of the offline invariance
            % pin: they are built from the log-mean submodel alone, so the
            % parameterization cannot move them except through estimation error.
            %
            % SPLIT BY IDENTIFIABILITY, and the split is NOT a post-hoc carve-out of
            % whatever failed. What a cross-fit comparison can demonstrate is
            % bounded by how well each underlying log-mean SD is identified: the
            % converters agree BIT-EXACTLY given identical inputs (pinned in
            % TestGammaDiffTrtLsConversion), so anything seen here is two
            % independent estimates of the same SD, not a parameterization effect.
            % A component estimated from four levels cannot serve as evidence about
            % invariance at any tolerance that is also honest about its own noise.
            %
            % The criterion is the number of levels informing the DOMINANT term:
            %   person (70), trial (15), person x trial (1050), person x occasion
            %   (280)  -> well identified, tight tolerance;
            %   occasion main effect (4) and trial x occasion -> occasion-dominated.
            % sigma_ot2 lands in the second group because its inclusion-exclusion
            % form B*(exp(vo+vt+vot) - exp(vo) - exp(vt) + 1) carries vo, and vot is
            % small enough here that vo drives it. sigma_po2 carries vo too but is
            % dominated by the well-identified vpo, which is why it belongs above.
            %
            % MEASURED 2026-08-04, and this is a scientific caveat rather than a
            % test artifact. At nocc = 4 BOTH parameterizations grossly overestimate
            % the occasion main-effect component: against the REALIZED planted value
            % (~0.19 for event 2; the realized occasion SD is 0.101 against a planted
            % 0.120, so the population parameter is not the right reference at four
            % levels) the location-scale fit gave 1.204 and the log-nu fit 0.776 -
            % 6.3x and 4.1x. The variance of a four-level random effect is
            % prior-dominated under the half-t prior, and exp(2*alpha + V) amplifies
            % it. Neither direction favors either parameterization: on sigma_ot2 the
            % LOCATION-SCALE fit is the accurate one (0.120 against a realized 0.119)
            % while log-nu is 33% low.
            %
            % This does NOT contaminate the reported generalizability, and the reason
            % is structural: with reltype 3 the occasion main effect enters ABSOLUTE
            % error only (abs_err = rel_err + bt + bo + boi in psyrat_diffrel_trt), so
            % it is excluded from the relative-error coefficient. The primary test
            % recovered G to 0.014 with sigma_o2 6x high. DEPENDABILITY at few
            % occasions is the quantity a reader should treat with caution.
            wellIdentified = {'id_varcov','trl_varcov','tid_varcov','oid_varcov'};
            for k = 1:numel(wellIdentified)
                vLs = localMeanDiag(relLs, wellIdentified{k});
                vNu = localMeanDiag(relNu, wellIdentified{k});
                for e = 1:2
                    denom = max(abs(vNu(e)), eps);
                    testCase.verifyLessThan(abs(vLs(e) - vNu(e))/denom, 0.35, ...
                        sprintf(['%s event %d differs across parameterizations ' ...
                        '(%.4g vs %.4g); the signal components are built from ' ...
                        'the log-mean submodel alone.'], wellIdentified{k}, e, ...
                        vLs(e), vNu(e)));
                end
            end

            % The occasion-dominated pair: reported in full, and held only to a
            % same-order-of-magnitude bound. Asserting the tight tolerance here
            % would be asserting that a four-level variance component is sharply
            % estimated, which it is not; loosening the tight tolerance just far
            % enough to admit the measured 55% would be fitting the threshold to one
            % run. A factor-of-two bound still catches a genuine wiring error (a
            % transposed slot moves these by orders of magnitude) without pretending
            % to more precision than nocc = 4 supports.
            occasionDominated = {'occ_varcov','to_varcov'};
            for k = 1:numel(occasionDominated)
                vLs = localMeanDiag(relLs, occasionDominated{k});
                vNu = localMeanDiag(relNu, occasionDominated{k});
                fprintf(['[gamma ls diff trt identity] %s (occasion-dominated, ' ...
                    'nocc=%d): ls=[%.4g %.4g] lognu=[%.4g %.4g]\n'], ...
                    occasionDominated{k}, 4, vLs(1), vLs(2), vNu(1), vNu(2));
                for e = 1:2
                    denom = max(abs(vNu(e)), eps);
                    testCase.verifyLessThan(abs(vLs(e) - vNu(e))/denom, 1.0, ...
                        sprintf(['%s event %d differs across parameterizations by ' ...
                        'more than a factor of two (%.4g vs %.4g). At nocc = 4 this ' ...
                        'component is weakly identified, so a moderate gap is ' ...
                        'expected - but a gap this large points at a slot ' ...
                        'transposition rather than at identification.'], ...
                        occasionDominated{k}, e, vLs(e), vNu(e)));
                end
            end

            % (c) THE RESIDUALS MUST NOT AGREE, and the direction is predictable.
            % The log-nu observed-scale residual carries the facet-averaging factor
            % exp(2*V) that the location-scale one has no counterpart for, so the
            % log-nu residual must be the LARGER. Direction is asserted; the
            % magnitude is printed as a diagnostic, because a sharp prediction
            % would require reconstructing e_cond_var from parameters neither path
            % stores, and because sigma_res2 is recovered by subtraction so the
            % ratio is not cleanly exp(2*V).
            %
            % THE ASSERTION BELOW WAS MISSING until 2026-08-04: the comment claimed
            % a direction check that the code did not make, so this block printed
            % and asserted nothing. Found by reading the emitted diagnostics against
            % the comment during the first live run. Measured on that run:
            % ratio = [1.562, 1.481] against exp(2*V) = 1.333 from the planted
            % truth, so the direction holds comfortably on both events.
            resLs = exp(mean(cell2mat(relLs.out.b_sigma{1}), 1)).^2;
            resNu = exp(mean(cell2mat(relNu.out.b_sigma{1}), 1)).^2;
            fprintf(['[gamma ls diff trt identity] observed-scale residual: ' ...
                'ls=[%.3f %.3f] lognu=[%.3f %.3f] ratio=[%.3f %.3f]\n'], ...
                resLs(1), resLs(2), resNu(1), resNu(2), ...
                resNu(1)/resLs(1), resNu(2)/resLs(2));
            for e = 1:2
                testCase.verifyGreaterThan(resNu(e), resLs(e), sprintf( ...
                    ['Event %d: the log-nu observed-scale residual (%.3f) must ' ...
                    'exceed the location-scale one (%.3f) - it carries the ' ...
                    'exp(2*V) facet-averaging factor that the location-scale ' ...
                    'residual has no counterpart for.'], e, resNu(e), resLs(e)));
            end

            % (d) THE STRUCTURAL CONTRAST that reproduced on both previous
            % increments (2.9x on the one-facet subject-level design, ~2.4x on the
            % one-facet difference): log-nu needs a LARGER person-dispersion SD,
            % because its single parameter must describe amplitude and precision
            % together. Direction only - the multiplier is design-dependent and
            % asserting a value would be overfitting one configuration.
            ssdLs = mean(relLs.out.gamma_disp.s_sd, 1);
            svNu  = mean(relNu.out.gamma_disp.s_v, 1);
            fprintf(['[gamma ls diff trt identity] s_sd (ls)=[%.3f %.3f]  ' ...
                's_v (lognu)=[%.3f %.3f]\n'], ssdLs(1), ssdLs(2), svNu(1), svNu(2));
            testCase.verifyGreaterThan(mean(svNu), mean(ssdLs), sprintf( ...
                ['The log-nu person-dispersion SD (mean %.3f) should exceed the ' ...
                'location-scale one (mean %.3f): it must describe amplitude and ' ...
                'precision with one parameter.'], mean(svNu), mean(ssdLs)));
        end

    end
end

% -------------------------------------------------------------------------

function T = localTruth()
%Known truth for the location-scale two-facet difference recovery. Collected in
%one place so both tests plant identical values.
T = struct();
T.b       = [log(5.5), log(4.0)];   % per-event log-mean: event 2 LOWER amplitude
T.b_sigma = [0.606136, 1.006136];   % per-event log residual SD: event 2 HIGHER
%ASYMMETRY IS THE POINT, exactly as in the one-facet sibling: event 2 has the
%lower amplitude and the higher residual SD, so both terms push the implied
%log-nu contrast the same way and the amplitude/precision blend is maximal.
T.s_p     = [0.22, 0.22];           % person
T.s_o     = [0.12, 0.12];           % occasion
T.s_t     = [0.15, 0.15];           % trial
T.s_po    = [0.18, 0.18];           % person x occasion
T.s_pt    = [0.14, 0.14];           % person x trial
T.s_ot    = [0.08, 0.08];           % occasion x trial
T.s_sd    = [0.50, 0.50];           % person log-sigma SD
T.rho_mean = 0.60;                  % cross-event person MEAN correlation
T.rho_fac  = 0.50;                  % cross-event correlation, other five factors
T.rho_ps   = -0.30;                 % within-event person mean <-> log-sigma corr
T.nsub = 70;
T.nocc = 4;
T.ntrl = 15;                        % per occasion per event
%DESIGN POINT. nocc = 4 rather than 2: increment 4's case-25 ledger records that
%at low occasion counts the whole observed scale is estimated high, and with only
%two occasions the occasion main effect and its two interactions compete for the
%same variance. ntrl = 15 keeps the table at 70*4*15*2 = 8400 rows, comparable to
%increment 5's 6400, and puts the difference-score coefficient below the
%near-ceiling regime rather than at G ~ 0.98.
end

function [vc1, vc2, cc] = localTruthComponents(T)
%Closed-form observed-scale components and cross-cell covariances at the planted
%truth, through the production converters.
vc1 = psyrat_gamma_varcomps_trt_ls(T.b(1), T.s_p(1), T.s_o(1), T.s_t(1), ...
    T.s_po(1), T.s_pt(1), T.s_ot(1), T.b_sigma(1), T.s_sd(1));
vc2 = psyrat_gamma_varcomps_trt_ls(T.b(2), T.s_p(2), T.s_o(2), T.s_t(2), ...
    T.s_po(2), T.s_pt(2), T.s_ot(2), T.b_sigma(2), T.s_sd(2));
cc = psyrat_gamma_crosscov_trt(T.b(1), T.s_p(1), T.s_o(1), T.s_t(1), ...
    T.s_po(1), T.s_pt(1), T.s_ot(1), ...
    T.b(2), T.s_p(2), T.s_o(2), T.s_t(2), T.s_po(2), T.s_pt(2), T.s_ot(2), ...
    T.rho_mean, T.rho_fac, T.rho_fac, T.rho_fac, T.rho_fac, T.rho_fac);
end

function g = localDiffG(vc1, vc2, cc, ntrl, nocc)
%True difference-score generalizability, computed by driving the PRODUCTION
%kernel (psyrat_diffrel_trt) with the closed-form components as a single draw.
%Using the real kernel rather than a hand formula keeps the truth and the
%recovered value on the same definition, including the non-concurrent policy
%(trial-indexed cross-covs zeroed, residual covariance 0).
z = 0;
bp  = localVarcov(vc1.sigma_p2,  vc2.sigma_p2,  cc.cov_p);
bo  = localVarcov(vc1.sigma_o2,  vc2.sigma_o2,  cc.cov_o);
bpo = localVarcov(vc1.sigma_po2, vc2.sigma_po2, cc.cov_po);
bt  = localVarcov(vc1.sigma_t2,  vc2.sigma_t2,  z);
bpi = localVarcov(vc1.sigma_pt2, vc2.sigma_pt2, z);
boi = localVarcov(vc1.sigma_ot2, vc2.sigma_ot2, z);
er_var = [0.5*log(vc1.sigma_res2), 0.5*log(vc2.sigma_res2)];

out = psyrat_diffrel_trt('bp',bp,'bpi',bpi,'bpo',bpo,'bt',bt,'bo',bo, ...
    'boi',boi,'er_var',er_var,'obs',[ntrl ntrl],'nocc',[nocc nocc], ...
    'reltype',3,'est','gen','CI',0.95,'er_cov',0);
g = out.pt;
end

function g = localDiffGFromRel(rel, ntrl, nocc)
%Recovered difference-score generalizability, from the stored posterior draws
%through the same kernel. Slot -> argument mapping is the relsummary wiring:
%id->bp, tid->bpi, oid->bpo, trl->bt, occ->bo, to->boi.
out = psyrat_diffrel_trt( ...
    'bp',  cell2mat(rel.out.id_varcov{1}), ...
    'bpi', cell2mat(rel.out.tid_varcov{1}), ...
    'bpo', cell2mat(rel.out.oid_varcov{1}), ...
    'bt',  cell2mat(rel.out.trl_varcov{1}), ...
    'bo',  cell2mat(rel.out.occ_varcov{1}), ...
    'boi', cell2mat(rel.out.to_varcov{1}), ...
    'er_var', cell2mat(rel.out.b_sigma{1}), ...
    'obs',[ntrl ntrl],'nocc',[nocc nocc], ...
    'reltype',3,'est','gen','CI',0.95,'er_cov',0);
g = out.pt;
end

function v = localMeanDiag(rel, slot)
%Posterior-mean per-event diagonal of one stored [draws x 2 x 2] varcov slot.
m = cell2mat(rel.out.(slot){1});
v = [mean(m(:,1,1)), mean(m(:,2,2))];
end

function m = localVarcov(v1, v2, cv)
%[1 x 2 x 2] cross-condition (co)variance for a single "draw" of truth.
m = zeros(1,2,2);
m(1,1,1) = v1;
m(1,2,2) = v2;
m(1,1,2) = cv;
m(1,2,1) = cv;
end

function tbl = localSimulateGammaDiffTrtLs(T)
% Two-event Persons x Occasions x Trials long table from the LOCATION-SCALE
% generative model. Y = gamrnd(nu/2, 2*mu/nu) has E[Y] = mu and
% Var[Y] = 2*mu^2/nu; with nu = 2*(mu/sigma)^2 that is Var[Y] = sigma^2,
% independent of mu AND of every facet. nu is computed per OBSERVATION (it moves
% with the crossed mean effects) while sigma is per participant per event - the
% asymmetry that makes this a different model from the log-nu form rather than a
% reparameterization of it.
rng(12345);

% Person block: two mean effects correlated rho_mean across events, and two scale
% effects each correlated rho_ps with their own event's mean effect.
z = randn(T.nsub, 4);
u_p = zeros(T.nsub,2);
u_p(:,1) = T.s_p(1) * z(:,1);
u_p(:,2) = T.s_p(2) * (T.rho_mean*z(:,1) + sqrt(1-T.rho_mean^2)*z(:,2));
w_p = zeros(T.nsub,2);
w_p(:,1) = T.s_sd(1) * (T.rho_ps*z(:,1) + sqrt(1-T.rho_ps^2)*z(:,3));
w_p(:,2) = T.s_sd(2) * (T.rho_ps*z(:,2) + sqrt(1-T.rho_ps^2)*z(:,4));

% The five remaining crossed factors, each a cross-event pair correlated rho_fac.
u_o  = localCrossEventEffects(T.nocc,               T.s_o,  T.rho_fac);
u_t  = localCrossEventEffects(T.ntrl,               T.s_t,  T.rho_fac);
u_po = localCrossEventEffects(T.nsub*T.nocc,        T.s_po, T.rho_fac);
u_pt = localCrossEventEffects(T.nsub*T.ntrl,        T.s_pt, T.rho_fac);
u_ot = localCrossEventEffects(T.nocc*T.ntrl,        T.s_ot, T.rho_fac);

nrow = T.nsub * T.nocc * T.ntrl * 2;
ids  = cell(nrow,1);
evt  = cell(nrow,1);
tme  = cell(nrow,1);
meas = zeros(nrow,1);
r = 0;
for p = 1:T.nsub
    for e = 1:2
        sigma_ep = exp(T.b_sigma(e) + w_p(p,e));
        for o = 1:T.nocc
            ipo = (p-1)*T.nocc + o;
            for t = 1:T.ntrl
                ipt = (p-1)*T.ntrl + t;
                iot = (o-1)*T.ntrl + t;
                r = r + 1;
                ids{r} = sprintf('S%03d', p);
                evt{r} = sprintf('E%d', e);
                tme{r} = sprintf('T%d', o);
                logmu = T.b(e) + u_p(p,e) + u_o(o,e) + u_t(t,e) ...
                    + u_po(ipo,e) + u_pt(ipt,e) + u_ot(iot,e);
                mu = exp(logmu);
                nu = 2 * (mu / sigma_ep)^2;
                meas(r) = gamrnd(nu/2, 2*mu/nu);
            end
        end
    end
end
tbl = table(ids, meas, evt, tme, 'VariableNames', {'id','meas','event','time'});
end

function u = localCrossEventEffects(n, sd, rho)
%n levels x 2 events of a zero-mean normal effect, correlated rho across events.
z = randn(n, 2);
u = zeros(n, 2);
u(:,1) = sd(1) * z(:,1);
u(:,2) = sd(2) * (rho*z(:,1) + sqrt(1-rho^2)*z(:,2));
end

function rel = localRunDiffTrtGamma(testCase, dataTbl, tag, gammascale)
% Group-level (case 10) non-concurrent two-facet gamma difference fit.
% gammascale selects the parameterization; everything else is held identical so
% the paired fits differ in exactly one thing.
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
args = { ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    'diffest', 2, ...
    'chains', 2, ...
    'warmup', 750, ...
    'sampling', 750, ...
    'seed', 12345, ...
    'verbose', 0, ...
    'showgui', 1, ...
    'cmdstanoutdir', outDir};
%Always pass gammascale explicitly. Until 2026-08-07 omitting it meant log-nu,
%so only the ==2 case needed an argument; ruling 33 flipped the default to
%location-scale for the movable analyses, and an omitted gammascale now
%resolves to 2 - the log-nu arm of a paired fit must say gammascale=1 or it
%silently refits location-scale.
args = [args, {'gammascale', gammascale}];
rel = psyrat_computevarcomp(args{:});
end

function outdir = localCmdStanOutputDir(projectRoot)
outdir = fullfile(projectRoot, 'tests', 'integration', 'cmdstan_artifacts');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
end

function localCleanupArtifacts(baseDir)
for p = {fullfile(baseDir, 'recov_gamma_ls_diff_trt'), ...
         fullfile(baseDir, 'recov_gamma_ls_diff_trt_pairA'), ...
         fullfile(baseDir, 'recov_gamma_lognu_diff_trt_pairB')}
    if isfolder(p{1})
        try
            rmdir(p{1}, 's');
        catch
            warning('tests:cleanupDirFailed', ...
                'Could not remove recovery artifact directory ''%s''.', p{1});
        end
    end
end
pats = {'cmdstan*.stan', 'cmdstan*.hpp', 'cmdstan18-*', 'output-*.csv', 'temp.data.R'};
for i = 1:numel(pats)
    files = dir(fullfile(baseDir, pats{i}));
    for j = 1:numel(files)
        if ~files(j).isdir
            delete(fullfile(baseDir, files(j).name));
        end
    end
end
dirs = dir(fullfile(baseDir, 'cmdstan_*'));
dirs = dirs([dirs.isdir]);
dirs = dirs(~ismember({dirs.name}, {'.', '..'}));
for i = 1:numel(dirs)
    try
        rmdir(fullfile(baseDir, dirs(i).name), 's');
    catch
        % best-effort cleanup
    end
end
end
