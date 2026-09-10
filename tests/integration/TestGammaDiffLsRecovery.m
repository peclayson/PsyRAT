classdef TestGammaDiffLsRecovery < PsyRATTestBase
    % Live CmdStan known-truth recovery for the Gamma / scaled-chi-square static
    % two-event DIFFERENCE designs under the LOCATION-SCALE parameterization
    % (analysis 7 group-level and analysis 8 subject-level, gammascale = 2,
    % NON-concurrent). Increment 5 of the gamma location-scale rollout.
    %
    % GENERATIVE MODEL. Per event e and participant p,
    %   log mu_ept     = b[e] + u_{p,e} + u_{t,e}
    %   log sigma_ep   = b_sigma[e] + w_{p,e}
    %   nu_ept         = 2*(mu_ept / sigma_ep)^2,   Y ~ scaled chi-square
    % so Var(Y | e,p) = sigma_ep^2 exactly, with no dependence on the mean.
    % (u_1, u_2, w_1, w_2) share a person block, matching owner ruling H.
    %
    % TWO DESIGN CHOICES THAT ARE THE POINT OF THIS TEST, not incidental:
    %
    % (1) b_sigma IS ASYMMETRIC ACROSS EVENTS (0.606 vs 1.006). The whole reason
    %     the difference designs are in scope is the EVENT contrast: since
    %     log sigma_e = log mu_e + 0.5*log(2) - 0.5*log(nu_e), a log-nu event
    %     contrast blends amplitude with residual SD. A symmetric design would
    %     validate the plumbing and leave that contrast - the one quantity the
    %     increment exists for - untested. The event MEANS are asymmetric too, and
    %     in the direction that makes the blend maximal: event 2 has the LOWER
    %     amplitude and the HIGHER residual SD, so both terms push the implied
    %     log-nu contrast the same way.
    %
    % (2) THE DESIGN POINT IS DELIBERATELY LOWER-RELIABILITY than increment 4's.
    %     That run sat at G ~ 0.98 and the rollout ledger records it as a
    %     self-acknowledged weak test of the coefficient surface. This one targets
    %     a difference-score G near 0.7 - a realistic ERP value, and the regime
    %     where the decision record predicts the parameterization matters most
    %     (person dispersion less sharply identified).
    %
    % TWO ASSERTIONS THIS TEST DELIBERATELY DOES NOT MAKE. Both were measured
    % false during increment 4 and are easy to reinstate by accident:
    %   - that a log-nu fit gets per-participant reliability WRONG at an adequate
    %     trial count. It does not; ind_nu is free per person and absorbs
    %     v_p = const + 2*u_p - 2*w_p. The paired fit below is used for the
    %     cross-parameterization IDENTITIES, not to show log-nu failing.
    %   - that partial corr(planted u_p, coefficient | planted w_p) is ~0. The
    %     coefficient is nonlinear in w_p while partial correlation removes only
    %     the linear part, so the true value is nonzero and design-dependent.
    %     Where that quantity appears it is compared against the TRUTH's own
    %     value, never against zero.
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

        function testLocationScaleDifferenceRecovered(testCase)
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required to draw continuous-df scaled chi-square.');

            T = localTruth();
            [dataTbl, u_p, w_p] = localSimulateGammaDiffLs(T);

            % ---- Closed-form observed-scale truth -------------------------------
            % The numerator is a POPULATION quantity built from the log-MEAN
            % submodel alone, so psyrat_gamma_varcomps_ls and psyrat_gamma_crosscov
            % give it exactly as the log-nu path would (pinned offline in
            % TestGammaDiffLsConversion).
            sigma_p2 = zeros(1,2);
            for e = 1:2
                vc = psyrat_gamma_varcomps_ls(T.b(e), T.s_p(e), T.s_i(e), ...
                    T.b_sigma(e), T.s_sd(e));
                sigma_p2(e) = vc.sigma_p2;
            end
            cc = psyrat_gamma_crosscov(T.b(1), T.s_p(1), T.s_i(1), ...
                T.b(2), T.s_p(2), T.s_i(2), T.rho_mean, 0);
            var_pdiff_true = sigma_p2(1) + sigma_p2(2) - 2*cc.cov_p_obs;

            % Each participant's per-event residual is POOLED (S17-FIX): the
            % population mean-surface person x trial term exp(2*b_e)*K_e (K
            % written out from the section-J identity, an independent oracle)
            % plus their own conditional variance sigma_ep^2 - no amplitude
            % term and no trial-averaging factor on the conditional part - and
            % the non-concurrent difference residual is the SUM over events
            % (ruling H).
            trueResid = zeros(T.nsub,1);
            pi_true = zeros(1,2);
            for e = 1:2
                Ke = exp(T.s_p(e)^2 + T.s_i(e)^2) * ...
                    (exp(T.s_p(e)^2) - 1) * (exp(T.s_i(e)^2) - 1);
                pi_true(e) = exp(2*T.b(e)) * Ke;
                trueResid = trueResid + ...
                    (pi_true(e) + exp(2*(T.b_sigma(e) + w_p(:,e)))) ./ T.ntrl;
            end
            trueG = var_pdiff_true ./ (var_pdiff_true + trueResid);

            % ---- Fit: case 8 (subject-level) under location-scale ---------------
            % Case 8 is fitted rather than case 7 because it emits BOTH: the
            % group-level components (inherited from the case-7 extraction) and the
            % per-subject residuals. One fit therefore exercises both extractions.
            rel = localRunDiffGamma(testCase, dataTbl, 'recov_gamma_ls_diff', 2);

            testCase.verifyEqual(rel.analysis, 'ic_diff_sserrvar');
            testCase.verifyEqual(rel.family, 'gamma');
            testCase.verifyEqual(rel.gammascale, 2);
            testCase.verifyEqual(rel.engine, 'cmdstan');

            % ---- Convergence gate: skip, never false-pass ------------------------
            [maxRhat, nFinite] = psyrat_conv_maxrhat(rel.out.conv.data);
            rel = psyrat_checkconv(rel);
            fprintf('[gamma ls diff] max R-hat=%.4f over %d params; converged=%d\n', ...
                maxRhat, nFinite, rel.out.conv.converged);
            testCase.assumeGreaterThan(nFinite, 0, ...
                'Convergence table empty; R-hat unextractable, so recovery cannot be validated.');
            testCase.assumeEqual(rel.out.conv.converged, 1, sprintf( ...
                ['Location-scale difference fit non-converged (max R-hat=%.3f); ' ...
                'skipping rather than reporting a pass on an unconverged fit.'], maxRhat));

            % ---- (1) the per-event scale submodel is recovered -------------------
            % Since S17-FIX, REL.out.b_sigma on the case-8 path stores the POOLED
            % typical-subject reference 0.5*log(pi_mean_e + exp(2*b_sigma_e)),
            % not the bare Stan per-event b_sigma. Its truth reference pools the
            % same way; at this design point the pooled shift is small (+0.008 /
            % +0.002 on the two events) but the reference must be the right one.
            bsigTrue = 0.5 .* log(pi_true + exp(2 .* T.b_sigma));  % 1 x 2
            bsigHat = mean(cell2mat(rel.out.b_sigma{1}), 1);       % 1 x 2
            ssdHat  = mean(rel.out.gamma_disp.s_sd, 1);            % 1 x 2

            fprintf(['[gamma ls diff] pooled b_sigma hat=[%.3f %.3f] (true [%.3f %.3f])  ' ...
                's_sd hat=[%.3f %.3f] (true [%.2f %.2f])\n'], ...
                bsigHat(1), bsigHat(2), bsigTrue(1), bsigTrue(2), ...
                ssdHat(1), ssdHat(2), T.s_sd(1), T.s_sd(2));

            for e = 1:2
                testCase.verifyLessThan(abs(bsigHat(e) - bsigTrue(e)), 0.30, ...
                    sprintf('Event %d pooled residual log-SD off: %.3f vs true %.3f.', ...
                    e, bsigHat(e), bsigTrue(e)));
            end

            % THE EVENT CONTRAST - the quantity this increment exists to make
            % readable. Asserted on its own, because both intercepts could be
            % biased together and still leave the contrast right (or vice versa).
            % The reference is the contrast of the POOLED per-event references.
            contrastHat  = bsigHat(2) - bsigHat(1);
            contrastTrue = bsigTrue(2) - bsigTrue(1);
            fprintf('[gamma ls diff] pooled b_sigma contrast hat=%.3f (true %.3f)\n', ...
                contrastHat, contrastTrue);
            testCase.verifyLessThan(abs(contrastHat - contrastTrue), 0.25, sprintf( ...
                ['The residual-SD EVENT CONTRAST is the estimand this design exists ' ...
                'to deliver: recovered %.3f vs true %.3f.'], contrastHat, contrastTrue));

            testCase.verifyGreaterThan(min(ssdHat), 0.20, sprintf( ...
                'Person log residual-SD SD collapsed: [%.3f %.3f] vs true [%.2f %.2f].', ...
                ssdHat(1), ssdHat(2), T.s_sd(1), T.s_sd(2)));

            % ---- (2) per-subject residuals track the planted truth ---------------
            er1 = rel.out.er_var_ss1{1};    % nd x nsub
            er2 = rel.out.er_var_ss2{1};
            wpss = rel.out.wp_cov_ss{1};
            testCase.verifyEqual(max(abs(wpss(:))), 0, ...
                'Non-concurrent runs must carry a structurally zero residual covariance.');

            idvar  = cell2mat(rel.out.id_varcov{1});
            trlvar = cell2mat(rel.out.trl_varcov{1});
            bsigma = cell2mat(rel.out.b_sigma{1});
            id_match = rel.out.id_matches{1};
            idtable = id_match;
            idtable.trls1 = repmat(T.ntrl, height(idtable), 1);
            idtable.trls2 = repmat(T.ntrl, height(idtable), 1);
            idtable.trls  = min(idtable.trls1, idtable.trls2);

            diff_table = psyrat_ssrel_diff('bp',idvar,'bt',trlvar,'er_var',bsigma, ...
                'wp_cov',zeros(size(bsigma,1),1),'idtable',idtable,'CI',0.95, ...
                'est','gen','er_var_ss1',er1,'er_var_ss2',er2,'wp_cov_ss',wpss);

            nrows = height(diff_table);
            idx = arrayfun(@(r) sscanf(char(string(diff_table.id(r))),'S%d'), (1:nrows)');
            recG     = diff_table.rel_pt;
            recResid = diff_table.ss_errvar;
            tG       = trueG(idx);
            tResid   = trueResid(idx);
            plantedW = w_p(idx,:);

            rResid = localPearson(recResid, tResid);
            fprintf(['[gamma ls diff] corr(per-subject resid, truth)=%.3f  ' ...
                'mean recG=%.3f (true %.3f)  G range [%.3f %.3f]\n'], ...
                rResid, mean(recG), mean(tG), min(recG), max(recG));

            testCase.verifyGreaterThan(rResid, 0.5, sprintf( ...
                ['Recovered per-subject difference residual does not track the truth ' ...
                '(r=%.3f); subject-level differentiation not recovered.'], rResid));
            testCase.verifyTrue(all(recG > 0 & recG < 1), ...
                'Every per-subject generalizability must lie in (0,1).');
            testCase.verifyLessThan(abs(mean(recG) - mean(tG)), 0.12, sprintf( ...
                'Mean per-subject generalizability off: %.3f vs true %.3f.', ...
                mean(recG), mean(tG)));

            % The design point must actually BE lower-reliability, or the
            % "stronger test of the coefficient surface" claim above is empty.
            testCase.verifyLessThan(mean(tG), 0.90, sprintf( ...
                ['This recovery is meant to sit well below increment 4''s G ~ 0.98; ' ...
                'true mean G is %.3f. If the design drifted upward, the coefficient ' ...
                'surface is being tested weakly again.'], mean(tG)));

            % Higher person residual SD must mean lower reliability. Summed over
            % events because the difference residual is the sum.
            rGW = localPearson(sum(plantedW,2), recG);
            fprintf('[gamma ls diff] corr(planted w_p sum, recovered G)=%.3f\n', rGW);
            testCase.verifyLessThan(rGW, -0.4, sprintf( ...
                ['Participants with larger planted residual SD must receive lower ' ...
                'reliability (r=%.3f).'], rGW));
        end

        function testCrossParameterizationIdentitiesHold(testCase)
            % The validation gate the rollout doc sets for every increment: fit the
            % SAME simulated data under both parameterizations and check the
            % algebraic identities that relate them.
            %
            %   b_sigma,e = b_e + 0.5*log(2) - 0.5*b_nu,e     (per-event intercepts)
            %
            % This is a statement about the two FITS, so it needs both. It is NOT a
            % claim that either fit is wrong - see the class header.
            %
            % Note the identity is checked through what the two paths actually STORE,
            % which is not the same quantity on both sides: the log-nu typical-subject
            % residual carries an exp(2*s_i^2) trial-averaging factor and the
            % location-scale one does not. The predicted ratio is therefore
            % exp(2*s_i^2), NOT 1. See the derivation at the comparison below.
            testCase.assumeTrue(logical(exist('gamrnd','file')), ...
                'Statistics Toolbox (gamrnd) required.');

            T = localTruth();
            dataTbl = localSimulateGammaDiffLs(T);

            relLs = localRunDiffGamma(testCase, dataTbl, 'recov_gamma_ls_diff_pairA', 2);
            relNu = localRunDiffGamma(testCase, dataTbl, 'recov_gamma_lognu_diff_pairB', 1);

            [rhatLs, nLs] = psyrat_conv_maxrhat(relLs.out.conv.data);
            [rhatNu, nNu] = psyrat_conv_maxrhat(relNu.out.conv.data);
            relLs = psyrat_checkconv(relLs);
            relNu = psyrat_checkconv(relNu);
            fprintf(['[gamma ls diff identity] ls: max R-hat=%.4f over %d params, ' ...
                'converged=%d | lognu: max R-hat=%.4f over %d params, converged=%d\n'], ...
                rhatLs, nLs, relLs.out.conv.converged, ...
                rhatNu, nNu, relNu.out.conv.converged);
            testCase.assumeGreaterThan(min(nLs, nNu), 0, ...
                'Convergence table empty for one of the paired fits.');
            testCase.assumeEqual(relLs.out.conv.converged, 1, ...
                'Location-scale fit did not converge; skipping the identity check.');
            testCase.assumeEqual(relNu.out.conv.converged, 1, ...
                'Log-nu fit did not converge; skipping the identity check.');

            % Per-event log-mean intercepts are shared by both models (the mean
            % submodel is byte-identical), so either fit supplies b.
            bLs  = mean(cell2mat(relLs.out.b{1}), 1);          % 1 x 2
            bNu  = mean(cell2mat(relNu.out.b{1}), 1);          % 1 x 2
            bsig = mean(cell2mat(relLs.out.b_sigma{1}), 1);    % 1 x 2, pooled b_sigma slot (S17-FIX)

            % THE TWO STORED RESIDUALS ARE NOT THE SAME QUANTITY, and getting this
            % reference right is the whole point of the check. Both paths store a
            % TYPICAL-SUBJECT per-event residual in REL.out.b_sigma, but:
            %
            %   log-nu:  exp(2*b_sigma) = 2*exp(2*b_e + 2*s_i^2 - b_nu,e)
            %   LS:      exp(2*b_sigma) = pi_mean_e + sigma_e^2     (S17-FIX pools)
            %
            % and ruling H's identity b_sigma,e = b_e + 0.5*log(2) - 0.5*b_nu,e says
            % sigma_e^2 = 2*exp(2*b_e - b_nu,e). Substituting,
            %
            %   resNu / resLs = exp(2*s_i^2) / (1 + pi_mean_e/sigma_e^2),
            %
            % i.e. the pre-S17-FIX prediction exp(2*s_i^2) divided by the pooling
            % correction (the log-nu subject-level path is NOT pooled - S17-FIX
            % is scoped to the location-scale designs). The exp(2*s_i^2) factor
            % is the trial-averaging term the log-nu residual carries because its
            % conditional variance is mean-coupled, and which the location-scale
            % conditional term has no counterpart for. Asserting these two
            % "agree" would be asserting the identity against the WRONG
            % reference, and would pass only because the tolerance is loose
            % enough to swallow the very asymmetry this increment turns on.
            resLs = exp(2*bsig);          % pooled typical-subject residual (S17-FIX)
            bsigNuSlot = mean(cell2mat(relNu.out.b_sigma{1}), 1);
            resNu = exp(2*bsigNuSlot);                         % carries exp(2*s_i^2)

            % pooled correction from the PLANTED truth: pi_mean_e/sigma_e^2, with
            % K written out from the section-J identity (independent oracle)
            Kvec = exp(T.s_p.^2 + T.s_i.^2) .* (exp(T.s_p.^2) - 1) .* (exp(T.s_i.^2) - 1);
            piOverCond = exp(2 .* T.b) .* Kvec ./ exp(2 .* T.b_sigma);
            predictedLogRatio = 2 * T.s_i.^2 - log(1 + piOverCond);   % 1 x 2

            fprintf(['[gamma ls diff identity] b (ls)=[%.3f %.3f] b (lognu)=[%.3f %.3f]\n' ...
                '                          typical-subject residual: ls=[%.3f %.3f] ' ...
                'lognu=[%.3f %.3f]\n'], bLs(1), bLs(2), bNu(1), bNu(2), ...
                resLs(1), resLs(2), resNu(1), resNu(2));

            % (a) the shared mean submodel must agree across the two fits
            for e = 1:2
                testCase.verifyLessThan(abs(bLs(e) - bNu(e)), 0.15, sprintf( ...
                    ['Event %d log-mean intercept differs across parameterizations ' ...
                    '(%.3f vs %.3f); the mean submodel is byte-identical, so these ' ...
                    'must agree.'], e, bLs(e), bNu(e)));
            end

            % (b) the residual identity, referenced to exp(2*s_i^2) rather than to 1.
            %
            % HONEST LIMIT OF THIS ASSERTION, stated so the number is not over-read.
            % The predicted log-ratio is 2*s_i^2 = 0.08 at this design point, while
            % log(resNu/resLs) = 2*(b_sigma_nu - b_sigma_ls) and the per-intercept
            % recovery error measured on the primary test was ~0.07. So the SIGNAL is
            % comparable to unpaired estimation noise. What makes the check
            % informative anyway is that the two fits use the SAME simulated data and
            % the SAME seed, so their errors are correlated and the DIFFERENCE is far
            % better determined than either intercept alone. The measured value is
            % printed either way, and the tolerance is set to what a paired
            % comparison at this design point can actually support - it is a
            % consistency check on the reference point, NOT a sharp discrimination of
            % exp(2*s_i^2) from 1. A sharp version would need a materially larger
            % s_i, which would change the design point of the primary test.
            for e = 1:2
                logRatio = log(resNu(e) / resLs(e));
                fprintf(['[gamma ls diff identity] event %d log(resNu/resLs)=%.4f ' ...
                    '(predicted 2*s_i^2 - log(1+pi/cond) = %.4f; ratio %.3f vs ' ...
                    'predicted %.3f)\n'], ...
                    e, logRatio, predictedLogRatio(e), ...
                    exp(logRatio), exp(predictedLogRatio(e)));
                testCase.verifyLessThan(abs(logRatio - predictedLogRatio(e)), 0.25, ...
                    sprintf(['Event %d: the two parameterizations'' typical-subject ' ...
                    'residuals do not stand in the relation the cross-parameterization ' ...
                    'identity predicts. log(resNu/resLs)=%.4f, predicted %.4f ' ...
                    '(ls %.3f, log-nu %.3f).'], e, logRatio, predictedLogRatio(e), ...
                    resLs(e), resNu(e)));

                % Direction, asserted separately because it is the qualitative
                % content and survives even where the magnitude is noise-limited:
                % the log-nu typical-subject residual must EXCEED the pooled
                % location-scale one at this design point - the trial-averaging
                % factor it carries (log 0.08 per event) dominates the pooling
                % correction on the LS side (log 0.016 / 0.004). A design point
                % where the pooled term rivaled exp(2*s_i^2) would flip this;
                % re-derive before re-tuning localTruth.
                testCase.verifyGreaterThan(resNu(e), resLs(e), sprintf( ...
                    ['Event %d: the log-nu typical-subject residual (%.3f) must exceed ' ...
                    'the pooled location-scale one (%.3f) at this design point - the ' ...
                    'exp(2*s_i^2) trial-averaging factor outweighs the S17-FIX ' ...
                    'pooling correction here.'], e, resNu(e), resLs(e)));
            end

            % (c) THE STRUCTURAL CONTRAST that survives (increment 4 measured 2.9x
            % on the one-facet design): log-nu needs a LARGER person-dispersion SD
            % than location-scale, because its parameter carries the amplitude term
            % as well as the precision term. Direction only - the multiplier is
            % design-dependent and asserting a value here would be overfitting.
            ssdLs = mean(relLs.out.gamma_disp.s_sd, 1);
            svNu  = mean(relNu.out.gamma_disp.s_v, 1);
            fprintf('[gamma ls diff identity] s_sd (ls)=[%.3f %.3f]  s_v (lognu)=[%.3f %.3f]\n', ...
                ssdLs(1), ssdLs(2), svNu(1), svNu(2));
            testCase.verifyGreaterThan(mean(svNu), mean(ssdLs), sprintf( ...
                ['The log-nu person-dispersion SD (mean %.3f) should exceed the ' ...
                'location-scale one (mean %.3f): it must describe amplitude and ' ...
                'precision with one parameter.'], mean(svNu), mean(ssdLs)));
        end

    end
end

% -------------------------------------------------------------------------

function T = localTruth()
%Known truth for the location-scale two-event difference recovery. Collected in
%one place so the two tests plant identical values.
T = struct();
T.b       = [log(5.5), log(4.0)];   % per-event log-mean: event 2 LOWER amplitude
T.b_sigma = [0.606136, 1.006136];   % per-event log residual SD: event 2 HIGHER
T.s_p     = [0.20, 0.20];           % person log-mean SD (small -> lower reliability)
T.s_i     = [0.20, 0.20];           % trial  log-mean SD
T.s_sd    = [0.50, 0.50];           % person log-sigma SD
T.rho_mean = 0.60;                  % cross-event person MEAN correlation
T.rho_ps   = -0.30;                 % within-event person mean <-> log-sigma corr
T.nsub = 80;
T.ntrl = 40;                        % per event
%DESIGN POINT, chosen deliberately and worth not re-tuning casually. Computed
%from the closed form above: true mean difference-score G ~ 0.67. Increment 4 sat
%at G ~ 0.98 and its own ledger records that as a weak test of the coefficient
%surface, so this is a large step down. It is not lower still because of a real
%trade-off: reducing ntrl further (ntrl = 25 gives G ~ 0.55) pushes the person
%scale effects toward the weakly-identified regime, where recovery becomes
%noise-limited and the assertions turn flaky rather than more informative. ~0.67
%is both a realistic ERP difference-score reliability and adequately identified.
end

function [tbl, u_p, w_p] = localSimulateGammaDiffLs(T)
% Two-event Persons x Trials long table from the LOCATION-SCALE generative model.
% Y = gamrnd(nu/2, 2*mu/nu) has E[Y] = mu and Var[Y] = 2*mu^2/nu; with
% nu = 2*(mu/sigma)^2 that is Var[Y] = sigma^2, independent of mu. nu is computed
% per OBSERVATION (it moves with the trial effect) while sigma is per participant
% per event -- the asymmetry that makes this a different model from the log-nu
% form rather than a reparameterization of it.
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

u_t = [T.s_i(1)*randn(T.ntrl,1), T.s_i(2)*randn(T.ntrl,1)];

nrow = T.nsub * T.ntrl * 2;
ids  = cell(nrow,1);
evt  = cell(nrow,1);
meas = zeros(nrow,1);
r = 0;
for p = 1:T.nsub
    for e = 1:2
        sigma_ep = exp(T.b_sigma(e) + w_p(p,e));
        for t = 1:T.ntrl
            r = r + 1;
            ids{r} = sprintf('S%03d', p);
            evt{r} = sprintf('E%d', e);
            mu = exp(T.b(e) + u_p(p,e) + u_t(t,e));
            nu = 2 * (mu / sigma_ep)^2;
            meas(r) = gamrnd(nu/2, 2*mu/nu);
        end
    end
end
tbl = table(ids, meas, evt, 'VariableNames', {'id','meas','event'});
end

function r = localPearson(x, y)
% Pearson correlation via base MATLAB (no Statistics Toolbox).
x = x(:) - mean(x(:)); y = y(:) - mean(y(:));
r = (x' * y) / (sqrt(x' * x) * sqrt(y' * y));
end

function rel = localRunDiffGamma(testCase, dataTbl, tag, gammascale)
% Subject-level (case 8) non-concurrent gamma difference fit. gammascale selects
% the parameterization; everything else is held identical so the paired fits
% differ in exactly one thing.
baseDir = localCmdStanOutputDir(testCase.projectRoot());
outDir = fullfile(baseDir, tag);
if ~exist(outDir, 'dir'), mkdir(outDir); end
args = { ...
    'data', dataTbl, ...
    'family', 'gamma', ...
    'diffest', 2, ...
    'sserrvar', 2, ...
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
for p = {fullfile(baseDir, 'recov_gamma_ls_diff'), ...
         fullfile(baseDir, 'recov_gamma_ls_diff_pairA'), ...
         fullfile(baseDir, 'recov_gamma_lognu_diff_pairB')}
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
