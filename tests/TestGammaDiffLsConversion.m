classdef TestGammaDiffLsConversion < matlab.unittest.TestCase
    %TESTGAMMADIFFLSCONVERSION Offline checks for the LOCATION-SCALE gamma
    %two-event DIFFERENCE conversion (analyses 7 and 8, gammascale = 2).
    %
    %WHY THIS FILE EXISTS. psyrat_gamma_extract_diff_ls and
    %psyrat_gamma_extract_diff_sserr_ls are LOCAL functions of
    %psyrat_computevarcomp: they cannot be called from a test and they only run
    %against a live CmdStan fit. This class therefore pins the algebra they rest
    %on, using the callable helpers they delegate to, plus the REAL downstream
    %consumer (psyrat_ssrel_diff).
    %
    %Three layers, mirroring TestGammaSserrLsConversion:
    %  (a) the numerator invariance that lets psyrat_gamma_crosscov and the
    %      id_varcov/trl_varcov assembly be reused unchanged;
    %  (b) the per-subject log-SD encoding round-trip psyrat_ssrel_diff depends on;
    %  (c) the difference-specific property - sigma^2_Delta,p = sigma^2_1,p +
    %      sigma^2_2,p with a structurally zero residual covariance - and the
    %      absent trial-averaging factor that distinguishes this parameterization
    %      from its log-nu twin.
    %
    %A green run here does NOT establish that estimation recovers truth; that is
    %the live-recovery lane's job (TestGammaDiffLsRecovery). It establishes that
    %the algebra the extractors rest on is what this repo believes it is.
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

        function testPerEventGroupComponentsAreParameterizationInvariant(testCase)
            % THE LOAD-BEARING ASSUMPTION of psyrat_gamma_extract_diff_ls. It
            % builds id_varcov/trl_varcov from per-event psyrat_gamma_varcomps_ls
            % diagonals plus psyrat_gamma_crosscov off-diagonals, exactly as its
            % log-nu twin does with psyrat_gamma_varcomps. If the per-event
            % diagonals disagreed between the two converters, the location-scale
            % difference coefficient would differ from the log-nu one for a reason
            % having nothing to do with the residual - and the whole "only the
            % residual changes" argument in the builder header would be false.
            %
            % Two events with DIFFERENT amplitudes and different scale parameters,
            % because a symmetric setup could hide an amplitude-dependent error.
            alpha = [1.70; 1.40];   % per-event log-mean (event 1, event 2)
            s_p   = [0.35; 0.42];
            s_i   = [0.20; 0.11];
            % Location-scale scale inputs
            logsig = [0.61; 0.48];
            s_sd   = [0.50; 0.33];
            % Deliberately UNRELATED log-nu inputs: the point is that whatever the
            % dispersion side is doing, the mean-surface fields do not move.
            nu    = exp([2.89; 2.70]);
            s_v   = [0.70; 0.40];
            covpv = [0.05; -0.02];

            vcLs = psyrat_gamma_varcomps_ls(alpha, s_p, s_i, logsig, s_sd);
            vcNu = psyrat_gamma_varcomps(alpha, s_p, s_i, nu, s_v, covpv);

            testCase.verifyEqual(vcLs.sigma_p2, vcNu.sigma_p2, ...
                'Per-event numerator sigma_p2 must not depend on the scale parameterization.');
            testCase.verifyEqual(vcLs.sigma_i2, vcNu.sigma_i2, ...
                'Per-event trial component sigma_i2 must not depend on the scale parameterization.');
            testCase.verifyEqual(vcLs.mu_bar, vcNu.mu_bar, ...
                'Per-event mean surface must not depend on the scale parameterization.');

            % ... and the residual MUST differ, or the test above would be vacuous
            % (it would pass just as well if the two converters were the same
            % function). This is the non-vacuity guard.
            testCase.verifyNotEqual(vcLs.sigma_pi_e2, vcNu.sigma_pi_e2, ...
                ['The two parameterizations must give DIFFERENT residuals; if they '...
                'agree here the invariance assertions above prove nothing.']);
        end

        function testCrossEventCovariancesAreReusedUnchanged(testCase)
            % psyrat_gamma_crosscov is called with IDENTICAL arguments on both
            % paths - it takes no dispersion input at all - so the off-diagonals
            % of id_varcov/trl_varcov are parameterization-invariant by
            % construction. This pins that the full assembled matrices agree, not
            % just the diagonals, which is what the extractor header claims and
            % what justifies having no psyrat_gamma_varcomps_diff_ls.
            alpha = [1.70; 1.40]; s_p = [0.35; 0.42]; s_i = [0.20; 0.11];
            cor_p = 0.55; cor_i = 0.20;

            vcLs_1 = psyrat_gamma_varcomps_ls(alpha(1), s_p(1), s_i(1), 0.61, 0.50);
            vcLs_2 = psyrat_gamma_varcomps_ls(alpha(2), s_p(2), s_i(2), 0.48, 0.33);
            vcNu_1 = psyrat_gamma_varcomps(alpha(1), s_p(1), s_i(1), exp(2.89), 0.70, 0.05);
            vcNu_2 = psyrat_gamma_varcomps(alpha(2), s_p(2), s_i(2), exp(2.70), 0.40, -0.02);

            cc = psyrat_gamma_crosscov(alpha(1), s_p(1), s_i(1), ...
                alpha(2), s_p(2), s_i(2), cor_p, cor_i);

            idLs = localVarcov(vcLs_1.sigma_p2, vcLs_2.sigma_p2, cc.cov_p_obs);
            idNu = localVarcov(vcNu_1.sigma_p2, vcNu_2.sigma_p2, cc.cov_p_obs);
            trlLs = localVarcov(vcLs_1.sigma_i2, vcLs_2.sigma_i2, cc.cov_i_obs);
            trlNu = localVarcov(vcNu_1.sigma_i2, vcNu_2.sigma_i2, cc.cov_i_obs);

            testCase.verifyEqual(idLs, idNu, ...
                'The assembled person var-cov must be parameterization-invariant.');
            testCase.verifyEqual(trlLs, trlNu, ...
                'The assembled trial var-cov must be parameterization-invariant.');
        end

        function testPerSubjectResidualEncodingRoundTripsExactly(testCase)
            % psyrat_ssrel_diff reads the per-subject residual as
            % exp(er_var_ss).^2. Since S17-FIX the extractor stores the POOLED
            % residual log-SD, er_var_ss = 0.5*log(pi_mean_e +
            % exp(2*(b_sigma[e] + w_{p,e}))) (psyrat_computevarcomp.m,
            % psyrat_gamma_extract_diff_sserr_ls); the CONDITIONAL term inside
            % it is still b_sigma[e] + w_{p,e} with no conversion, which is what
            % the check below exercises.
            %
            % SCOPE, corrected. This exercises NO production code: the extractor
            % is a local function this file cannot call, so nothing here can
            % detect a change to it (the strict source scan in
            % testPooledResidualEncodingIsInSource is what pins the stored
            % form). What survives below is an IEEE-association check on
            % exp(x)^2 vs exp(2x) over 50 draws x 4 subjects x 2 events - worth
            % keeping, because psyrat_ssrel_diff consumes exp(er_var_ss).^2 and
            % the two spellings must not drift apart at the tolerances used
            % downstream. A behavioral regression in the encoding itself is
            % caught only by the CmdStan-gated
            % tests/integration/TestGammaDiffLsRecovery.
            nd = 50; NSUB = 4;
            rng(20260804, 'twister');   % local, deterministic; no hidden state
            b_sigma = [0.61 + 0.02*randn(nd,1), 0.48 + 0.02*randn(nd,1)];
            w = {0.5*randn(nd,NSUB), 0.5*randn(nd,NSUB)};

            for e = 1:2
                er_var_ss = b_sigma(:,e) + w{e};      % the stored encoding

                % REMOVED: verifyEqual(exp(er_var_ss).^2, sigma_pe.^2, 'AbsTol', 0)
                % with sigma_pe = exp(er_var_ss) assigned on the line above -
                % the same expression compared against itself, which could not
                % fail under any change to the toolbox.
                testCase.verifyEqual(exp(er_var_ss).^2, exp(2.*(b_sigma(:,e) + w{e})), ...
                    'RelTol', 1e-15, ...
                    'exp(b_sigma + w)^2 must agree with exp(2*(b_sigma + w)).');
            end
        end

        function testSubjectResidualCarriesNoTrialAveragingFactor(testCase)
            % WHAT THIS TEST IS AND IS NOT. It DOCUMENTS the algebra below and
            % quantifies the error a "symmetry" edit would introduce. It does NOT
            % pin the production extractor, and an earlier version of this header
            % claimed it did ("If someone adds the term back, this fails with a
            % number that names the size of the error"). That was false: the
            % extractor is a local function of psyrat_computevarcomp.m, so no test
            % in this file can call it, and restoring exp(2*s_i^2) there leaves
            % every assertion here green. The claim was removed rather than
            % weakened, because a test that names itself THE SCIENTIFIC PIN and is
            % not one is worse than no test.
            %
            % What DOES cover the property: nothing in the offline lane. The live
            % tests/integration/TestGammaDiffLsRecovery reads b_sigma rather than
            % er_var_ss and its per-subject tolerances are an order of magnitude
            % looser than the ~1.083 inflation this factor would cause, so it would
            % not reliably catch it either. Treat the property as UNPINNED.
            %
            % The log-nu per-subject residual is an expectation over the trial
            % facet of a MEAN-COUPLED conditional variance, so it picks up the
            % lognormal second moment of the trial effect sitting INSIDE mu:
            %   sigma_pi_e2 = 2*exp(2*(b + u_pe) + 2*s_i^2 - (b_nu + v_pe)).
            % Under location-scale Var(Y_pt | mu, sigma_{e,p}) = sigma^2_{e,p} has
            % no t in it, so there is no trial average to take and the factor has
            % no counterpart on the CONDITIONAL term - which since S17-FIX is
            % one summand of the stored pooled residual, not the whole of it.
            % (s_i does legitimately reach the residual through the POOLED
            % person x trial term pi_mean_e, where it is a mean-surface
            % quantity, not a scaling of the conditional variance.)
            % psyrat_ssrel_diff divides er_var_ss by the participant's TRIAL
            % COUNT, so what is stored must be the SINGLE-TRIAL residual
            % variance.
            %
            % The two assertions below state the consequence quantitatively:
            % restoring the factor would inflate the conditional term of every
            % per-subject residual by exactly exp(2*s_i^2), a constant > 1
            % whenever there is any trial variance, which would depress every
            % per-subject dependability. Both are algebraically forced by the
            % two lines that construct them, so they are documentation with a
            % worked number, not guards.
            s_i = 0.20;
            b_sigma = 0.61;
            w_pe = [-0.4; 0.0; 0.4];

            lsResidual = exp(2.*(b_sigma + w_pe));     % the conditional term
            withFactor = lsResidual .* exp(2.*s_i.^2); % the wrong version

            testCase.verifyEqual(withFactor ./ lsResidual, ...
                repmat(exp(2*s_i^2), size(w_pe)), 'RelTol', 1e-12, ...
                'The erroneous trial factor would be a constant exp(2*s_i^2) inflation.');
            testCase.verifyGreaterThan(exp(2*s_i^2), 1, ...
                ['With any trial variance the factor exceeds 1, so restoring it '...
                'would systematically DEPRESS subject-level dependability.']);

            % REMOVED: an assertion that the residual "must not depend on s_i at
            % all", comparing lsResidual against a second variable assigned the
            % textually identical expression exp(2.*(b_sigma + w_pe)). s_i appeared
            % in neither side, so nothing was varied and it could not fail. The
            % property it named is real and is stated in the header; restating it
            % as a self-comparison only made the class look better covered than it
            % is. Pinning it for real needs the production encoding at
            % psyrat_computevarcomp.m to be reachable, which a local function is
            % not - a source scan is the repo's usual answer (see B13/B14 in
            % TestGammaScaleOption) and is a maintainer call, not a mechanical fix.
        end

        function testMarginalResidualIsTheClosedFormSubjectAverage(testCase)
            % The deliberate estimand gap between REL.out.b_sigma (the
            % TYPICAL-subject reference, w = 0 - since S17-FIX the POOLED
            % typical-subject residual) and REL.out.err_varcov (the
            % subject-AVERAGED population component) is documented as
            % "do not reconcile them". Under location-scale the gap ON THE
            % CONDITIONAL-VARIANCE TERM is available in closed form, so it can
            % be CHECKED rather than merely asserted:
            %   E_p[exp(2*(b_sigma + w_p))] = exp(2*b_sigma + 2*s_sd^2),
            % which is exactly psyrat_gamma_varcomps_ls's e_cond_var. That makes
            % the marginalization factor exp(2*s_sd^2) - the F the provenance
            % reports - a testable quantity rather than a claim. (The pooled
            % pi_mean term is common to both sides of that gap and drops out.)
            b_sigma = 0.61; s_sd = 0.50;
            NSUB = 200000;
            rng(20260804, 'twister');
            w_p = s_sd .* randn(NSUB,1);

            empirical = mean(exp(2.*(b_sigma + w_p)));
            closedForm = exp(2*b_sigma + 2*s_sd^2);

            % Monte-Carlo over a lognormal is heavy-tailed, so the tolerance is
            % loose on purpose; the point is that the two agree in kind, not that
            % this sample size pins three decimals.
            testCase.verifyEqual(empirical, closedForm, 'RelTol', 0.05, ...
                'Subject-averaged residual must match the closed-form marginal.');

            % The closed form must be exactly what psyrat_gamma_varcomps_ls uses
            % for e_cond_var, since the extractor's b_sigma comes from that field.
            vc = psyrat_gamma_varcomps_ls(1.70, 0.35, 0.20, b_sigma, s_sd);
            testCase.verifyEqual(vc.e_cond_var, closedForm, 'RelTol', 1e-12, ...
                'e_cond_var must be exp(2*log_sigma + 2*s_sd^2).');

            % ... and the typical-subject reference is strictly smaller, which is
            % WHY the two must not be reconciled: under the log link the typical
            % subject is the median, not the mean.
            testCase.verifyLessThan(exp(2*b_sigma), closedForm, ...
                'The typical-subject residual must sit below the subject-averaged one.');
        end

        function testNonconcurrentDifferenceResidualIsTheSumOfEventResiduals(testCase)
            % Owner ruling H: for the non-concurrent difference,
            %   sigma^2_{Delta,p} = sigma^2_{1,p} + sigma^2_{2,p},
            % the residual covariance being structurally zero. This drives the REAL
            % consumer (psyrat_ssrel_diff) with location-scale-shaped inputs and
            % checks the resulting per-subject generalizability against a hand
            % reference, so it tests the assembled contract rather than restating
            % the formula.
            %
            % Structurally zero, not merely assumed: correlating the two events'
            % person SCALE effects correlates the residual VARIANCES, not the
            % residuals, since E[e1*e2] = E[E[e1|p]*E[e2|p]] = 0. The inputs below
            % therefore use CORRELATED per-subject residual variances with
            % wp_cov_ss = 0, which is the configuration that would break if
            % someone "fixed" the zero covariance.
            nd = 200; NSUB = 3;
            rng(20260804, 'twister');

            % observed-scale person (co)variance across the two events
            sp1 = 4.0 + 0.2*randn(nd,1);
            sp2 = 5.0 + 0.2*randn(nd,1);
            cov_p = -0.6 + 0.1*randn(nd,1);
            si1 = 0.15 + 0.02*randn(nd,1);
            si2 = 0.12 + 0.02*randn(nd,1);

            bp = zeros(nd,2,2);
            bp(:,1,1) = sp1; bp(:,2,2) = sp2; bp(:,1,2) = cov_p; bp(:,2,1) = cov_p;
            bt = zeros(nd,2,2);
            bt(:,1,1) = si1; bt(:,2,2) = si2;   % non-concurrent trial cross-cov = 0

            % LOCATION-SCALE per-subject residual log-SDs: b_sigma[e] + w_{p,e},
            % with the two events' scale effects deliberately CORRELATED (w2 built
            % from w1) to exercise the structurally-zero covariance claim.
            b_sigma_stan = [0.61 + 0.02*randn(nd,1), 0.48 + 0.02*randn(nd,1)];
            w1 = 0.5*randn(nd,NSUB);
            w2 = 0.7*w1 + sqrt(1-0.7^2)*0.5*randn(nd,NSUB);
            er_var_ss1 = b_sigma_stan(:,1) + w1;
            er_var_ss2 = b_sigma_stan(:,2) + w2;

            id_match = table((1:NSUB)', (1:NSUB)', 'VariableNames', {'id','id2'});
            idtable = id_match;
            idtable.trls  = repmat(20, NSUB, 1);
            idtable.trls1 = repmat(20, NSUB, 1);   % o1
            idtable.trls2 = repmat(25, NSUB, 1);   % o2

            out = psyrat_ssrel_diff('bp',bp,'bt',bt,'er_var',b_sigma_stan,...
                'wp_cov',zeros(nd,1),'idtable',idtable,'CI',0.95,'est','gen',...
                'er_var_ss1',er_var_ss1,'er_var_ss2',er_var_ss2,...
                'wp_cov_ss',zeros(nd,NSUB));

            for r = 1:NSUB
                o1 = idtable.trls1(r); o2 = idtable.trls2(r);
                wp1 = exp(er_var_ss1(:,r)).^2;    % sigma^2_{1,p}
                wp2 = exp(er_var_ss2(:,r)).^2;    % sigma^2_{2,p}
                var_pdiff = sp1 + sp2 - 2*cov_p;
                % the ruling's difference residual, per trial count, wp_cov = 0
                rel_err = wp1./o1 + wp2./o2;
                Gdraw = var_pdiff ./ (var_pdiff + rel_err);
                testCase.verifyEqual(out.rel_pt(r), mean(Gdraw), 'RelTol', 1e-10, ...
                    sprintf(['Subject %d: the non-concurrent difference residual '...
                    'must be sigma^2_1,p/o1 + sigma^2_2,p/o2 with zero covariance.'], r));
            end
        end

        function testPooledResidualEncodingIsInSource(testCase)
            % SOURCE-SCANNING GUARD for the analysis-8 extractor,
            % psyrat_gamma_extract_diff_sserr_ls - a LOCAL function of
            % psyrat_computevarcomp.m that no offline test can call (the tactic
            % B13/B14 use for exactly that situation; the analysis-6 twin is
            % scanned the same way in TestGammaSserrLsConversion). Three lines
            % are pinned in STRICT form, so a drifted variant fails rather than
            % passing vacuously - the assertions are on match COUNTS, and a
            % rename or reformat FAILS loudly by design:
            %  (1) the person scale effect comes from the SCALE slots of r_id
            %      (e+2), not the mean slots - the S16 decoupling;
            %  (2) er_var_ss stores the POOLED residual log-SD (S17-FIX):
            %      0.5*log(pi_mean_e + exp(2*(b_sigma[e] + w_pe)));
            %  (3) the population reference b_sigma stores the pooled
            %      typical-subject (w = 0) form of the same quantity.
            src = fileread(localDiffEngineSourcePath());
            body = localDiffFunctionBody(testCase, src, ...
                'psyrat_gamma_extract_diff_sserr_ls');
            bodyLines = regexp(body, '\r\n|\n|\r', 'split');

            pins = { ...
                '^\s*w_pe\s*=\s*r_id\(:,:,e\+2\)\s*;\s*(%.*)?$', ...
                'the scale-slot read w_pe = r_id(:,:,e+2)'; ...
                ['^\s*er_var_ss\{e\}\s*=\s*0\.5\.\*log\(pi_mean_e\s*\+\s*' ...
                 'exp\(2\.\*\(bsig\(:,e\)\s*\+\s*w_pe\)\)\)\s*;\s*(%.*)?$'], ...
                'the pooled per-subject encoding (S17-FIX)'; ...
                ['^\s*b_sigma\(:,e\)\s*=\s*0\.5\.\*log\(pi_mean_e\s*\+\s*' ...
                 'exp\(2\.\*bsig\(:,e\)\)\)\s*;\s*(%.*)?$'], ...
                'the pooled typical-subject reference (S17-FIX)'};

            for k = 1:size(pins,1)
                hits = ~cellfun(@isempty, regexp(bodyLines, pins{k,1}, 'once'));
                testCase.verifyEqual(nnz(hits), 1, sprintf( ...
                    ['Expected exactly one line matching %s in ' ...
                    'psyrat_gamma_extract_diff_sserr_ls and found %d. If the ' ...
                    'code is right and only the spelling moved, update this ' ...
                    'scanner deliberately rather than deleting it.'], ...
                    pins{k,2}, nnz(hits)));
            end
        end

    end
end

function m = localVarcov(v1, v2, c)
%Assemble the [draws x 2 x 2] observed-scale var-cov exactly as the extractors do,
%so the comparison is of the assembled matrix rather than of loose components.
n = numel(v1);
m = zeros(n,2,2);
m(:,1,1) = v1(:);
m(:,2,2) = v2(:);
m(:,1,2) = c(:);
m(:,2,1) = c(:);
end

function p = localDiffEngineSourcePath()
%Absolute path to the estimator core, derived from THIS file's location so the
%scan does not depend on pwd or on which copy of the toolbox is first on the
%MATLAB path (a .claude/worktrees copy shadowing the real one has bitten this
%repo before). Same helper as TestGammaSserrLsConversion's, duplicated because
%test-local functions cannot be shared across classdef files.
here = fileparts(mfilename('fullpath'));          % <root>/tests
p = fullfile(fileparts(here), 'subroutines', 'estimation', 'psyrat_computevarcomp.m');
end

function body = localDiffFunctionBody(testCase, src, fname)
%Return the source text of local function FNAME, from its `function` line to the
%line before the next one, failing loudly if it is missing or duplicated.
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
