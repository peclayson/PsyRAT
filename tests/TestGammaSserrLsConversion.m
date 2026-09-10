classdef TestGammaSserrLsConversion < matlab.unittest.TestCase
    %TESTGAMMASSERRLSCONVERSION Offline checks for the LOCATION-SCALE gamma
    %subject-level conversion (analyses 6 and 25, gammascale = 2).
    %
    %WHY THIS FILE EXISTS. psyrat_gamma_extract_sserr_ls and its two-facet twin
    %are LOCAL functions of psyrat_computevarcomp and cannot be called from a
    %test, and they only execute against a live CmdStan fit. This class therefore
    %works on RESTATEMENTS of the extractors' algebra, not on the extractors:
    %  (a) the invariance that lets them pass NaN for the dispersion arguments,
    %  (b) the log-SD encoding psyrat_ssrel depends on, and
    %  (c) the contrast between the two estimands - that under location-scale a
    %      participant's own amplitude no longer moves their coefficient.
    %
    %WHAT A GREEN RUN HERE DOES NOT ESTABLISH, stated because an earlier version of
    %this header overstated (c) as "the scientific property the whole increment
    %exists to deliver", which implied production coverage this class cannot have.
    %It does not establish that estimation recovers truth - that is the
    %live-recovery lane's job - and it does not pin the extractors themselves. A
    %regression inside psyrat_gamma_extract_sserr_ls, including re-coupling
    %ind_sdlog to ind_bs (the S16 failure mode on analysis 6), leaves this class
    %green. What it establishes is that the algebra the extractors rest on is what
    %this repo believes it is.
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

        function testOneFacetGroupComponentsAreParameterizationInvariant(testCase)
            % psyrat_gamma_extract_sserr_ls reads sigma_p2 (the reliability
            % numerator) and sigma_i2 (the absolute-error term) from
            % psyrat_gamma_varcomps_ls, while its log-nu twin reads them from
            % psyrat_gamma_varcomps. If those two disagreed, the location-scale
            % coefficient would differ from the log-nu one for a reason that has
            % nothing to do with the residual, and the comparison the decision
            % record makes would be meaningless.
            alpha = [1.70; 1.55];  s_p = [0.35; 0.42];  s_i = [0.20; 0.11];
            lognu = [0.61; 0.55];  s_sd = [0.50; 0.33];

            vcLs = psyrat_gamma_varcomps_ls(alpha, s_p, s_i, lognu, s_sd);
            vcNu = psyrat_gamma_varcomps(alpha, s_p, s_i, exp([2.89; 2.70]), ...
                [0.7; 0.4], [0.05; -0.02]);

            testCase.verifyEqual(vcLs.sigma_p2, vcNu.sigma_p2, ...
                'Numerator sigma_p2 must not depend on the scale parameterization.');
            testCase.verifyEqual(vcLs.sigma_i2, vcNu.sigma_i2, ...
                'Trial component sigma_i2 must not depend on the scale parameterization.');
            testCase.verifyEqual(vcLs.mu_bar, vcNu.mu_bar);
        end

        function testTrtGroupComponentsIgnoreDispersionArguments(testCase)
            % THE LOAD-BEARING ASSUMPTION of psyrat_gamma_extract_trt_sserr_ls.
            % That extractor calls psyrat_gamma_varcomps_trt with NaN for
            % (nu, s_v, cov_pv) because the six group components it consumes are
            % built from alpha and the log-MEAN SDs alone. This test is what makes
            % that safe rather than lucky: if a future edit routes a dispersion
            % argument into any of the six, this fails immediately.
            a = [1.70; 1.55]; sp = [0.35; 0.42]; so = [0.18; 0.22];
            st = [0.20; 0.11]; spo = [0.09; 0.14]; spt = [0.12; 0.07];
            sot = [0.05; 0.08];

            lo = psyrat_gamma_varcomps_trt(a,sp,so,st,spo,spt,sot, ...
                exp([2.89; 2.70]), [0.7; 0.4], [0.05; -0.02]);
            hi = psyrat_gamma_varcomps_trt(a,sp,so,st,spo,spt,sot, ...
                exp([0.10; 5.00]), [2.5; 0.0], [-0.9; 0.9]);

            fields = {'sigma_p2','sigma_o2','sigma_t2','sigma_po2', ...
                'sigma_pt2','sigma_ot2','mu_bar'};
            for k = 1:numel(fields)
                testCase.verifyEqual(lo.(fields{k}), hi.(fields{k}), ...
                    sprintf(['%s changed when only the dispersion arguments ' ...
                    'changed. psyrat_gamma_extract_trt_sserr_ls passes NaN for ' ...
                    'those arguments precisely because this cannot happen.'], ...
                    fields{k}));
            end
        end

        function testNaNDispersionArgumentsLeaveGroupComponentsFinite(testCase)
            % The NaN sentinel must not contaminate what the extractor reads.
            % Pinning this separately from the invariance above, because "does
            % not depend on" and "survives NaN" are different claims: a
            % 0*NaN anywhere in the expression would satisfy the first and fail
            % the second.
            a = 1.70; sp = 0.35; so = 0.18; st = 0.20;
            spo = 0.09; spt = 0.12; sot = 0.05;

            ref = psyrat_gamma_varcomps_trt(a,sp,so,st,spo,spt,sot, exp(2.89), 0.7, 0.05);
            nanArgs = psyrat_gamma_varcomps_trt(a,sp,so,st,spo,spt,sot, NaN, NaN, NaN);

            fields = {'sigma_p2','sigma_o2','sigma_t2','sigma_po2', ...
                'sigma_pt2','sigma_ot2','mu_bar'};
            for k = 1:numel(fields)
                testCase.verifyTrue(isfinite(nanArgs.(fields{k})), ...
                    sprintf('%s went non-finite under NaN dispersion arguments.', fields{k}));
                testCase.verifyEqual(nanArgs.(fields{k}), ref.(fields{k}), ...
                    sprintf('%s differs under the NaN sentinel.', fields{k}));
            end
        end

        function testLogSdEncodingRoundTripsExactly(testCase)
            % psyrat_ssrel forms wp = exp(wp_pop + wp_ss)^2. Since S17-FIX the
            % analysis-6 extractor manufactures that encoding from the pooled
            % residual (0.5*log of pi_mean plus the conditional variance), as
            % the log-nu path always did for its own residual; the two-facet
            % analysis-25 extractor still passes the fitted Intercept_sigma and
            % ind_sd straight through.
            %
            % WHAT REMAINS HERE, AND WHAT WAS REMOVED. This method used to assert
            % the PASS-THROUGH as well - "there is no arithmetic between the fitted
            % parameters and the slots", called "the real content" - by assigning
            % pop_sdlog = logsig and ind_sdlog = ind_sd and then verifying each
            % against the variable it had just been assigned from, at AbsTol 0.
            % Those two assertions could not fail under any change to the toolbox,
            % and the pass-through they named lives in a local function of
            % psyrat_computevarcomp.m that no test in this file can call. Both are
            % gone, along with the production comment that credited them.
            %
            % What survives is the RECONSTRUCTION check, which is real: exp(a+b)^2
            % and exp(2*(a+b)) are algebraically identical but are different
            % floating-point operation sequences, so they differ in the last ULP.
            % psyrat_ssrel consumes the first spelling, so the two must stay within
            % machine precision. AbsTol 0 would assert a property of IEEE
            % association rather than of this code, which is why it is RelTol.
            logsig = [0.61; 0.55; 0.70];
            ind_sd = [0.10 -0.25 0.40; 0.05 -0.10 0.22; -0.30 0.15 0.02];

            wp       = exp(logsig + ind_sd).^2;
            expected = exp(2.*(logsig + ind_sd));
            testCase.verifyEqual(wp, expected, 'RelTol', 1e-15, ...
                'The location-scale encoding must reconstruct sigma_p^2 to machine precision.');
        end

        function testSubjectCoefficientIsIndependentOfOwnAmplitude(testCase)
            % THE POINT OF THE WHOLE INCREMENT (finding S16, owner ruling H).
            %
            % Under log-nu the per-participant residual is
            %   2*exp(2*(alpha + u_p) + 2*s_i^2 - (beta + v_p)),
            % evaluated at that participant's OWN log-mean, while the numerator is
            % a population quantity. So at equal dispersion a higher-amplitude
            % participant is scored as LESS reliable. Under location-scale the
            % residual is exp(2*(Intercept_sigma + w_p)) and never touches mu, so
            % amplitude drops out entirely.
            %
            % SCOPE, corrected on two counts. (a) Both branches are computed here
            % with hand-written arithmetic that mirrors psyrat_ssrel.m; psyrat_ssrel
            % itself is NEVER called, and an earlier version of this header said
            % "the SAME psyrat_ssrel kernel", which reads as though it were.
            % (b) This does not pin the production extractor:
            % psyrat_gamma_extract_sserr_ls is a local function of
            % psyrat_computevarcomp.m, so a regression re-coupling ind_sdlog to
            % ind_bs there - the S16 failure mode on analysis 6 - leaves this class
            % green. The analogous property IS pinned against production code for
            % analysis 26 by TestGammaDynrelLsSserrReliability; the static
            % analysis-6 path has no such coverage.
            %
            % What this method genuinely establishes is the CONTRAST between the two
            % estimands under a shared design point, which is the S16 argument.
            alpha = 1.70; s_i = 0.20; beta = 2.89; s_p = 0.35;
            u_p   = [-0.5 0 0.5];        % three amplitudes, identical dispersion
            v_p   = [0 0 0];
            ntrl  = 30;

            vc = psyrat_gamma_varcomps(alpha, s_p, s_i, exp(beta));
            bp = sqrt(vc.sigma_p2);       % population numerator, common to both
            i2 = vc.sigma_i2;

            % log-nu: residual carries u_p
            piMean = max(vc.var_mu - vc.sigma_p2 - vc.sigma_i2, 0);   % pooled person x trial term (S17-FIX / S20)
            wpNu = piMean + 2 .* exp(2.*(alpha + u_p) + 2.*s_i.^2 - (beta + v_p));
            depNu = bp.^2 ./ (bp.^2 + (wpNu + i2)./ntrl);

            % location-scale: residual is read off the person scale effect alone.
            % Give every participant the SAME scale effect, so the only thing
            % varying across them is amplitude.
            w_p    = [0 0 0];
            logsig = 0.5*log(2) + alpha + s_i.^2 - 0.5*beta;   % matched at u_p = 0
            wpLs   = piMean + exp(2.*(logsig + w_p));
            depLs  = bp.^2 ./ (bp.^2 + (wpLs + i2)./ntrl);

            % POSITIVE CONTROL for the flatness assertion below. On its own, that
            % assertion compares three copies of one number and so cannot fail; what
            % makes the design point meaningful is that it is SENSITIVE to the
            % failure mode. Re-couple the residual to the participant's own log-mean
            % - the S16 regression - and the coefficient must move materially across
            % the same u_p grid. Without this control the flatness check would be
            % indistinguishable from a degenerate design point.
            wpLsCoupled  = exp(2.*(logsig + w_p + u_p));
            depLsCoupled = bp.^2 ./ (bp.^2 + (wpLsCoupled + i2)./ntrl);
            testCase.verifyGreaterThan(abs(depLsCoupled(1) - depLsCoupled(3)), 0.01, ...
                ['Positive control failed: at this design point even a residual ' ...
                'coupled to u_p would not move the coefficient, so the flatness ' ...
                'asserted below would prove nothing. Re-choose the design point.']);

            % log-nu: amplitude moves the coefficient, and downward.
            testCase.verifyTrue(depNu(1) > depNu(2) && depNu(2) > depNu(3), ...
                ['Under log-nu a higher-amplitude participant must receive a ' ...
                'LOWER coefficient at equal dispersion. If this stops holding, ' ...
                'finding S16 has been invalidated and the ruling needs revisiting.']);
            testCase.verifyGreaterThan(depNu(1) - depNu(3), 0.01, ...
                'The log-nu amplitude effect must be materially large, not epsilon.');

            % location-scale: it does not.
            testCase.verifyEqual(depLs(1), depLs(3), 'AbsTol', 0, ...
                ['Under location-scale the per-participant coefficient must be ' ...
                'EXACTLY independent of the participant''s own log-mean.']);

            % and the two agree at the reference point, confirming the contrast
            % above is about amplitude rather than an accidental level shift.
            testCase.verifyEqual(depLs(2), depNu(2), 'RelTol', 1e-12, ...
                'The two parameterizations must agree at u_p = 0 by construction.');
        end

        function testCrossParameterizationPersonIdentity(testCase)
            % Owner ruling H's validation requirement, person-level half:
            %   u_sigma = u_mu - 0.5*u_nu
            % Derived by taking 0.5*log of the log-nu residual and splitting into
            % population and person parts. Pinned numerically so a sign or factor
            % slip in either extractor is caught offline.
            alpha = 1.70; s_i = 0.20; beta = 2.89;
            u_p = [-0.40 0.15 0.62];
            v_p = [ 0.30 -0.55 0.11];

            sigNu = 2 .* exp(2.*(alpha + u_p) + 2.*s_i.^2 - (beta + v_p));
            logSigmaTotal = 0.5 .* log(sigNu);

            popPart    = 0.5*log(2) + alpha + s_i.^2 - 0.5*beta;  % u_p = v_p = 0
            personPart = logSigmaTotal - popPart;

            testCase.verifyEqual(personPart, u_p - 0.5.*v_p, 'RelTol', 1e-12, ...
                'The person-level identity u_sigma = u_mu - 0.5*u_nu must hold exactly.');

            % Intercept half of the identity, and the ONE documented residue.
            % Ruling H gives b_sigma = b + 0.5*log(2) - 0.5*b_nu. The population
            % part above equals that PLUS s_i^2, and that surplus is exactly the
            % trial-averaging factor exp(2*s_i^2) the decision record names as the
            % reason the two forms are different models rather than a
            % reparameterization. Pinning it here keeps that difference visible
            % instead of folded into a tolerance.
            rulingIntercept = alpha + 0.5*log(2) - 0.5*beta;
            testCase.verifyEqual(popPart - rulingIntercept, s_i.^2, 'RelTol', 1e-12, ...
                ['The only departure from the ruling''s intercept identity must ' ...
                'be the trial-averaging term s_i^2.']);
        end

        function testPersonScaleEffectIsNotCoupledToAmplitudeInSource(testCase)
            % SOURCE-SCANNING GUARD, and the offline lane's ONLY protection for the
            % estimand claim of this whole increment (finding S16, owner ruling H).
            %
            % THE PROPERTY. Under location-scale the participant-specific part of
            % the residual is read off the person SCALE effect alone, with no
            % term in ind_bs (the person log-MEAN effect, u_p). That decoupling
            % is what makes the per-participant coefficient an ABSOLUTE-precision
            % quantity rather than precision confounded with signal amplitude.
            % It is the difference between the two estimands, not an
            % implementation detail. S17-FIX did not change it: the pooled
            % person x trial term is built from the POPULATION alpha, s_p and
            % s_i, so ind_bs still must not appear in any residual-channel line.
            %
            % WHY A SOURCE SCAN. The two extractors that carry it -
            % psyrat_gamma_extract_sserr_ls (analysis 6) and
            % psyrat_gamma_extract_trt_sserr_ls (analysis 25) - are LOCAL functions
            % of psyrat_computevarcomp.m, so no offline test can call them, and
            % they execute only against a live CmdStan fit. This is the tactic
            % B13/B14 use for exactly that situation.
            %
            % WHY THE BEHAVIORAL TESTS IN THIS CLASS CANNOT DO IT. They restate the
            % algebra rather than running the extractor, so they compute whatever
            % the test itself writes down. An earlier version of
            % testSubjectCoefficientIsIndependentOfOwnAmplitude asserted this
            % property against a design point built with u_p already absent - true
            % by construction, and green under the very regression it named.
            %
            % WHAT A VIOLATION WOULD DO. Re-coupling ind_sdlog to ind_bs
            % reintroduces the log-nu behavior: at equal dispersion a
            % higher-amplitude participant is scored as LESS reliable. Every
            % recovery test would still recover, every coefficient would still lie
            % in [0,1], and the contamination would surface only as a spurious
            % association between participants' amplitude and their estimated
            % reliability. No gate in this repo sees that.
            src = fileread(localEngineSourcePath());

            % Per-extractor STRICT expected form of the ind_sdlog assignment.
            % Analysis 6 pools (S17-FIX): the remainder over the pooled
            % typical-subject reference, with ind_sd inside the conditional term
            % and pi_mean built from POPULATION quantities. Analysis 25 (two
            % facets, no pooling analog) keeps the pure pass-through. Any term
            % in ind_bs fails both patterns.
            extractors = { ...
                'psyrat_gamma_extract_sserr_ls', ...
                ['^\s*ind_sdlog\s*=\s*0\.5\.\*log\(pi_mean\s*\+\s*exp\(' ...
                 '2\.\*\(logsig\(:\)\s*\+\s*ind_sd\)\)\)\s*-\s*pop_sdlog' ...
                 '\s*;\s*(%.*)?$'], ...
                ['must assign the S17-FIX pooled remainder ind_sdlog = ' ...
                 '0.5.*log(pi_mean + exp(2.*(logsig(:) + ind_sd))) - pop_sdlog']; ...
                'psyrat_gamma_extract_trt_sserr_ls', ...
                '^\s*ind_sdlog\s*=\s*ind_sd\s*;\s*(%.*)?$', ...
                'must assign ind_sdlog = ind_sd and nothing else'};

            for k = 1:size(extractors,1)
                body = localFunctionBody(testCase, src, extractors{k,1});
                bodyLines = regexp(body, '\r\n|\n|\r', 'split');

                % The expected assignment, in STRICT form. A drifted variant -
                % `... + ind_bs`, a scaled copy, a rename - fails here. A
                % trailing comment is allowed.
                %
                % NOTE the assertion is on the match COUNT. If this pattern ever
                % stops matching - a rename, a reformat - the test FAILS rather
                % than passing vacuously. That is deliberate: MATLAB spells the
                % word boundary \< and \>, NOT \b (which is a literal backspace,
                % char 8), and a scanner that silently matches nothing is the
                % failure mode B14 was written to avoid.
                isExpected = ~cellfun(@isempty, regexp(bodyLines, ...
                    extractors{k,2}, 'once'));
                testCase.verifyEqual(nnz(isExpected), 1, ...
                    sprintf(['%s %s. Exactly one such line is required and %d ' ...
                    'were found. If a term in ind_bs was added, that is finding ' ...
                    'S16 reappearing: the per-participant coefficient would once ' ...
                    'again track the participant''s own amplitude, which is the ' ...
                    'behavior the location-scale parameterization exists to ' ...
                    'remove.'], ...
                    extractors{k,1}, extractors{k,3}, nnz(isExpected)));

                % S16 x S17-FIX interaction guard: NO residual-channel line may
                % reference ind_bs, including the pooled term's ingredients. The
                % pooled pi_mean must stay population-built.
                residualLines = ~cellfun(@isempty, regexp(bodyLines, ...
                    '^\s*(pi_mean|pop_sdlog|ind_sdlog)\s*=', 'once'));
                touchesIndBs = ~cellfun(@isempty, ...
                    regexp(bodyLines, '\<ind_bs\>', 'once'));
                testCase.verifyEqual(nnz(residualLines & touchesIndBs), 0, ...
                    sprintf(['A residual-channel assignment in %s references ' ...
                    'ind_bs - the S16 coupling the location-scale ' ...
                    'parameterization exists to remove.'], extractors{k,1}));

                % POSITIVE CONTROL. The assertion above is only meaningful if
                % ind_bs is actually IN SCOPE at that point - otherwise "it is not
                % coupled to ind_bs" would be true merely because the variable does
                % not exist, and the guard would be measuring nothing. It is in
                % scope: both extractors pull it from the fit before the assignment
                % and store it for provenance. If this control ever fails, the
                % strict check above has stopped being evidence of anything and
                % must be re-derived rather than trusted.
                hasIndBs = ~cellfun(@isempty, regexp(bodyLines, ...
                    '^\s*ind_bs\s*=', 'once'));
                testCase.verifyGreaterThanOrEqual(nnz(hasIndBs), 1, ...
                    sprintf(['Positive control failed: ind_bs is no longer bound ' ...
                    'in %s, so the decoupling asserted above is vacuous rather ' ...
                    'than a property of the code.'], extractors{k,1}));
            end
        end

    end
end

function p = localEngineSourcePath()
%Absolute path to the estimator core, derived from THIS file's location so the
%scan does not depend on pwd or on which copy of the toolbox is first on the
%MATLAB path (a .claude/worktrees copy shadowing the real one has bitten this
%repo before).
here = fileparts(mfilename('fullpath'));          % <root>/tests
p = fullfile(fileparts(here), 'subroutines', 'estimation', 'psyrat_computevarcomp.m');
end

function body = localFunctionBody(testCase, src, fname)
%Return the source text of local function FNAME, from its `function` line to the
%line before the next one. Scanning the BODY rather than the whole file is what
%keeps this guard specific: psyrat_computevarcomp.m contains other `ind_sdlog =`
%assignments, and the log-nu extractors legitimately MANUFACTURE the encoding
%(`0.5 .* log(sig_pi_e2) - pop_sdlog`) instead of passing it through. A file-wide
%scan would either match those too or have to special-case them.
lines = regexp(src, '\r\n|\n|\r', 'split');
isFuncStart = ~cellfun(@isempty, regexp(lines, '^function\s', 'once'));
startHit = find(~cellfun(@isempty, ...
    regexp(lines, ['^function\s.*\<' fname '\>\s*\('], 'once')));

testCase.assertEqual(numel(startHit), 1, sprintf( ...
    ['Expected exactly one definition of %s, found %d. It was renamed, ' ...
    'duplicated, or promoted to its own file - fix this scanner rather than ' ...
    'deleting it.'], fname, numel(startHit)));

laterStarts = find(isFuncStart);
laterStarts = laterStarts(laterStarts > startHit);
if isempty(laterStarts)
    stopAt = numel(lines);
else
    stopAt = laterStarts(1) - 1;
end
body = strjoin(lines(startHit:stopAt), newline);
end
