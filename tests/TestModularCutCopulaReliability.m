classdef TestModularCutCopulaReliability < PsyRATTestBase
    %TESTMODULARCUTCOPULARELIABILITY Stage 3 of the modular cut-copula workflow:
    % the per-participant concurrent difference-score calculator (analysis 19)
    % and the cross-covariance wrapper it depends on.
    %
    % No CmdStan: every test drives the calculators from synthetic draws, so what
    % is pinned is the conversion arithmetic and the estimand choices, not that
    % the model recovers truth. Live recovery and the external-oracle replay are
    % what check correctness end to end.
    %
    % Organized around the three claims this design makes that its non-concurrent
    % sibling (analysis 13) does not:
    %   1 the residuals are COUPLED, through a quadrature-derived covariance
    %   2 the residual channel is POOLED, per Rast & Clayson
    %   3 the participant's OWN location effect enters - through the error, via
    %     the point at which the copula covariance is evaluated - while the
    %     signal stays a population quantity

    methods (Test)

        %% -- the reference oracle --------------------------------------------

        function testPerParticipantCoefficientsMatchTheReferenceOracle(testCase)
            %Every reported coefficient, per participant, against the
            %reference's own reliability_onefacet, from the frozen in-repo
            %stage-3 fixture (the one-facet mirror of the two-facet lane in
            %TestModularCutCopulaTrtReliability). The one-facet path also has
            %the external replay against the owner's completed fits, but that
            %check is BUNDLE-GATED and skips wherever the bundle is absent (CI
            %included); this lane is the external pin that always runs. The
            %fixture's standardized n_i' (7) differs from its simulated
            %design's I (5), so the D-study divisor is exercised rather than
            %coinciding with the actual counts.
            fx = localLoadOnefacetStage3();
            [gen,dep] = psyrat_ssrel_diffdynrel_rescor_gamma_ls( ...
                fx.args{:}, 'CI', 0.95);

            for p = 1:height(fx.people)
                testCase.verifyEqual(gen.rel_pt(p), ...
                    fx.meanOver('G_delta', p), 'RelTol', 1e-7, ...
                    sprintf('participant %d: G_delta', p));
                testCase.verifyEqual(dep.rel_pt(p), ...
                    fx.meanOver('D_delta', p), 'RelTol', 1e-7, ...
                    sprintf('participant %d: D_delta', p));
                testCase.verifyEqual(gen.icc_pt(p), ...
                    fx.meanOver('ICC_relative_delta', p), 'RelTol', 1e-7, ...
                    sprintf('participant %d: relative ICC', p));
                testCase.verifyEqual(dep.icc_pt(p), ...
                    fx.meanOver('ICC_absolute_delta', p), 'RelTol', 1e-7, ...
                    sprintf('participant %d: absolute ICC', p));
                testCase.verifyEqual(gen.sem_pt(p), ...
                    fx.meanOver('SEM_relative_delta', p), 'RelTol', 1e-7, ...
                    sprintf('participant %d: relative SEM', p));
                testCase.verifyEqual(dep.sem_pt(p), ...
                    fx.meanOver('SEM_absolute_delta', p), 'RelTol', 1e-7, ...
                    sprintf('participant %d: absolute SEM', p));
            end
        end

        function testDiagnosticMatchesTheReferenceResidualCovariance(testCase)
            %The per-participant quadrature provenance must describe the same
            %quantity the reference computed, with the observed-scale residual
            %correlation reported separately from the copula correlation.
            fx = localLoadOnefacetStage3();
            [~,~,diagnostic] = psyrat_ssrel_diffdynrel_rescor_gamma_ls( ...
                fx.args{:}, 'CI', 0.95);

            for p = 1:height(fx.people)
                testCase.verifyEqual(diagnostic.res_cov_pt(p), ...
                    fx.meanOver('residual_covariance_person', p), 'RelTol', 1e-7, ...
                    sprintf('participant %d: residual covariance', p));
                testCase.verifyEqual(diagnostic.res_cor_pt(p), ...
                    fx.meanOver('residual_correlation_person', p), 'RelTol', 1e-7, ...
                    sprintf('participant %d: residual correlation', p));
            end
            testCase.verifyEqual(diagnostic.quad_unconverged_n, ...
                zeros(height(fx.people),1), ...
                'The reference fixture should not need the non-strict fallback.');
            testCase.verifyNotEqual(diagnostic.res_cor_pt, diagnostic.rho_e_pt, ...
                'Observed-scale residual correlation is NOT the copula correlation.');
        end

        %% -- the cross-covariance wrapper ------------------------------------

        function testWrapperReproducesTheShippedTermsExactly(testCase)
            %The person and trial cross-covariances must be bit-identical to the
            %shipped non-concurrent helper, or the two designs would silently
            %disagree on components that are supposed to be shared.
            [a1,p1,i1,a2,p2,i2,cp,ci] = localCrossArgs();
            base = psyrat_gamma_crosscov(a1,p1,i1,a2,p2,i2,cp,ci);
            got = psyrat_gamma_crosscov_rescor(a1,p1,i1,a2,p2,i2,cp,ci,0);

            testCase.verifyEqual(got.cov_p_obs, base.cov_p_obs, 'AbsTol', 0);
            testCase.verifyEqual(got.cov_i_obs, base.cov_i_obs, 'AbsTol', 0);
        end

        function testMeanSurfaceCrossTermEqualsTheSubtractionItReplaces(testCase)
            %The product form is used because the subtraction cancels. Both must
            %agree at ordinary magnitudes; the product is simply the accurate one
            %when the log covariances get small.
            [a1,p1,i1,a2,p2,i2,cp,ci] = localCrossArgs();
            got = psyrat_gamma_crosscov_rescor(a1,p1,i1,a2,p2,i2,cp,ci,0);

            subtraction = got.cov_mean_total - got.cov_p_obs - got.cov_i_obs;
            testCase.verifyEqual(got.cov_pimean_obs, subtraction, 'RelTol', 1e-9);

            %And the pooled channel is that term plus whatever residual
            %covariance the caller supplied.
            cov_e = 0.37 * ones(size(got.cov_pimean_obs));
            withE = psyrat_gamma_crosscov_rescor(a1,p1,i1,a2,p2,i2,cp,ci,cov_e);
            testCase.verifyEqual(withE.cov_pi_e_obs, ...
                withE.cov_pimean_obs + cov_e, 'AbsTol', 1e-12);
        end

        function testWrapperRejectsAMismatchedResidualLength(testCase)
            [a1,p1,i1,a2,p2,i2,cp,ci] = localCrossArgs();
            testCase.verifyError( ...
                @() psyrat_gamma_crosscov_rescor(a1,p1,i1,a2,p2,i2,cp,ci,[1;2]), ...
                'psyrat_gamma_crosscov_rescor:size');
        end

        %% -- claim 2: the pooled residual ------------------------------------

        function testPooledResidualIsTheMeanSurfaceTermPlusThePersonVariance(testCase)
            %The identity the calculator relies on: passing a participant's own
            %log residual SD as log_sigma with s_sd = 0 makes sigma_pi_e2 the
            %POOLED channel - the mean-surface person x trial term built from the
            %population alpha, plus that participant's conditional variance.
            %That is what Rast & Clayson's residual is, and it is the whole
            %difference from analysis 13's residual.
            alpha = log(6.1); sp = 0.33; si = 0.19; ss = log(1.7);

            vc = psyrat_gamma_varcomps_ls(alpha, sp, si, ss, 0);

            %Independent reconstruction: the mean-surface residual is the total
            %lognormal variance minus the two named components.
            v = sp^2 + si^2;
            var_mu = exp(2*alpha + v) * expm1(v);
            p_term = exp(2*alpha + v) * expm1(sp^2);
            i_term = exp(2*alpha + v) * expm1(si^2);
            pi_mean = var_mu - p_term - i_term;

            testCase.verifyEqual(vc.sigma_pi_e2, pi_mean + exp(2*ss), ...
                'RelTol', 1e-10, ...
                'sigma_pi_e2 must be the pooled residual channel.');
            %And the mean-surface term is strictly positive, so pooling always
            %ADDS to the residual - the direction S17 records.
            testCase.verifyGreaterThan(pi_mean, 0);
        end

        %% -- claims 1 and 3: the calculator ----------------------------------

        function testProducesOneAdmissibleRowPerParticipant(testCase)
            f = localFixture();
            [gen,dep,diagnostic] = localRun(f);

            testCase.verifyEqual(height(gen), height(f.idtable));
            testCase.verifyEqual(height(dep), height(f.idtable));
            testCase.verifyEqual(height(diagnostic), height(f.idtable));
            %row order and the real id2 are preserved, not the remapped 1 used
            %per call
            testCase.verifyEqual(gen.id2(:)', f.idtable.id2(:)');
            testCase.verifyTrue(all(gen.rel_pt > 0 & gen.rel_pt < 1));
            %absolute error can only exceed relative error, so dependability can
            %never exceed generalizability for the same participant
            testCase.verifyTrue(all(dep.rel_pt <= gen.rel_pt + 1e-12));
        end

        function testCopulaCorrelationChangesTheAnswer(testCase)
            %CLAIM 1. If rho did nothing, the concurrent design would be its
            %non-concurrent sibling wearing a different name. A positive residual
            %correlation REDUCES the difference score's error (the two events'
            %errors partly cancel), so reliability must rise relative to rho = 0.
            f = localFixture();
            zero = localRun(f, 'rho', zeros(f.ndraws,1));
            pos  = localRun(f, 'rho', 0.6*ones(f.ndraws,1));
            neg  = localRun(f, 'rho', -0.6*ones(f.ndraws,1));

            testCase.verifyTrue(all(pos.rel_pt > zero.rel_pt), ...
                'A positive residual correlation must raise difference reliability.');
            testCase.verifyTrue(all(neg.rel_pt < zero.rel_pt), ...
                'A negative residual correlation must lower it.');
        end

        function testAtZeroCorrelationOnlyTheMeanSurfaceCrossTermSurvives(testCase)
            %The residual covariance channel must collapse to the mean-surface
            %term alone when the copula says the residuals are independent -
            %not to zero, because the interaction cross term is still there.
            f = localFixture();
            [~,~,diagnostic] = localRun(f, 'rho', zeros(f.ndraws,1));
            testCase.verifyEqual(diagnostic.res_cov_pt, ...
                zeros(height(diagnostic),1), 'AbsTol', 1e-12);
            testCase.verifyEqual(diagnostic.quad_nodes_max, ...
                zeros(height(diagnostic),1), ...
                'The quadrature must short-circuit rather than integrate at rho = 0.');
        end

        function testTheParticipantsOwnLocationEffectMovesTheErrorOnly(testCase)
            %CLAIM 3, and the reason mu_ss1/2 are stored at estimation time.
            %Shifting one participant's own location effect must move THEIR
            %coefficient - because it moves the point at which the copula
            %covariance is evaluated - and must leave every other participant
            %untouched.
            f = localFixture();
            base = localRun(f);

            bumped = f;
            bumped.mu_ss1(:,2) = bumped.mu_ss1(:,2) + 0.8;
            moved = localRun(bumped);

            testCase.verifyNotEqual(moved.rel_pt(2), base.rel_pt(2));
            others = [1 3:height(base)];
            testCase.verifyEqual(moved.rel_pt(others), base.rel_pt(others), ...
                'AbsTol', 0, ...
                'One participant''s location effect must not touch another''s.');

            %...and it is inert when the residuals are uncoupled, which proves it
            %acts THROUGH the copula term rather than through the signal.
            b0 = localRun(f, 'rho', zeros(f.ndraws,1));
            m0 = localRun(bumped, 'rho', zeros(f.ndraws,1));
            testCase.verifyEqual(m0.rel_pt, b0.rel_pt, 'AbsTol', 0, ...
                ['With no residual coupling the participant''s own location '...
                'effect must not enter at all - the signal is a population '...
                'quantity.']);
        end

        function testDrawMisalignmentIsRejected(testCase)
            %The copula correlation exists only for the retained draws. Recycling
            %a short rho against full-length margins would pair unrelated draws
            %and produce a plausible number from an incoherent posterior.
            f = localFixture();
            testCase.verifyError(@() localRun(f, 'rho', zeros(f.ndraws-1,1)), ...
                'varargin:rho');
            testCase.verifyError(@() localRun(f, 'rho', 1.5*ones(f.ndraws,1)), ...
                'varargin:rho');

            short = f;
            short.mu_ss1 = short.mu_ss1(:,1:2);
            testCase.verifyError(@() localRun(short), 'varargin:mu_ss');
        end

        function testQuadratureNonConvergenceIsCountedNotFatal(testCase)
            %A tail draw can put the two components in regimes where no node pair
            %meets tolerance. That must be recorded per participant rather than
            %destroying the analysis - but it must never be silent.
            out = psyrat_gamma_copula_rescov(1.4e6, 1.6e-4, 0.0644, [], [], false);
            testCase.verifyFalse(out.converged);
            testCase.verifyTrue(isfinite(out.value));
            testCase.verifyGreaterThan(out.error, 1e-6);

            %strict mode (the default) still errors, so no existing caller
            %silently changed behavior
            testCase.verifyError( ...
                @() psyrat_gamma_copula_rescov(1.4e6, 1.6e-4, 0.0644), ...
                'psyrat_gamma_copula_rescov:convergence');
            %and a well-behaved argument reports converged under both modes
            ok = psyrat_gamma_copula_rescov(18, 30, 0.3, [], [], false);
            testCase.verifyTrue(ok.converged);
        end

        %% -- the summary disclosure (analysis 19) -------------------------

        function testOneFacetSummaryCarriesTheCutPosteriorDisclosure(testCase)
            %Same contract as the two-facet sibling test in
            %TestModularCutCopulaTrtReliability: the stage's provenance labels,
            %its diagnostic flag counts, and the copula correlation must reach
            %the summary. Analysis 19 is the design that ships live, so its
            %disclosure path gets its own pin rather than an inference from the
            %two-facet one. A clean stage-2 (no flag fields) must not warn.
            REL = localFakeOnefacetREL();
            s = testCase.verifyWarningFree(...
                @() psyrat_dynrel_summary(REL,'CI',0.95,'ngrid',4));

            testCase.verifyEqual(s.copula_labels.inference_framework, ...
                'modular_cut_two_stage_gamma_margins_gaussian_copula_v1');
            testCase.verifyEqual(s.copula_labels.joint_posterior_status, ...
                'not_a_full_joint_bayesian_posterior');
            testCase.verifyEqual(s.strata(1).copula_flags.n_draws, 6);
            testCase.verifyEqual(s.strata(1).copula_flags.n_invalid, 0);

            vc = s.strata(1).varcomp;
            row = strcmp(vc.symbol,'rho_e,12(copula)');
            testCase.verifyEqual(sum(row), 1, ...
                'Exactly one copula-correlation row in the varcomp table.');
            rc = REL.out.copula{1}.rho_e;
            testCase.verifyEqual(vc.estimate(row), mean(rc), 'AbsTol', 1e-12);
            testCase.verifyEqual(vc.component{row}, ...
                'Copula residual correlation (events)');

            %and the per-participant table still comes out alongside it
            testCase.verifyNotEmpty(s.strata(1).ssrel_table, ...
                'analysis 19 must still produce its per-participant table.');

            %invalidity flags warn on this path too
            REL2 = localFakeOnefacetREL();
            REL2.out.copula{1}.invalidity_flag = true(6,1);
            testCase.verifyWarning(...
                @() psyrat_dynrel_summary(REL2,'CI',0.95,'ngrid',4), ...
                'dynrel:copulainvalidity');
        end

        function testSharedFlagCounterContract(testCase)
            %The ONE reader both disclosure carriers share (2026-08-11 code
            %review replaced two inline copies that already disagreed on
            %missing-field handling). recorded=false marks a stratum nothing
            %could be summarized from; flags_complete=false marks best-effort
            %zero counts, so an affirmative "none flagged" is never built
            %from absent evidence.
            full = struct('sel',(1:4)','invalidity_flag',logical([1;0;1;0]),...
                'calibration_flag',false(4,1));
            f = psyrat_copula_flag_counts(full);
            testCase.verifyTrue(f.recorded && f.flags_complete);
            testCase.verifyEqual([f.n_draws f.n_invalid f.n_calibration], [4 2 0]);

            noflags = struct('sel',(1:4)');
            f2 = psyrat_copula_flag_counts(noflags);
            testCase.verifyTrue(f2.recorded);
            testCase.verifyFalse(f2.flags_complete, ...
                'Absent flag fields must be distinguishable from zero flagged.');
            testCase.verifyEqual([f2.n_draws f2.n_invalid f2.n_calibration], [4 0 0]);

            f3 = psyrat_copula_flag_counts(struct('rho_e',0.3));
            testCase.verifyFalse(f3.recorded, ...
                'A sel-less stored struct is unrecorded, not zero draws of evidence.');
            f4 = psyrat_copula_flag_counts([]);
            testCase.verifyFalse(f4.recorded, ...
                'A non-struct entry must not error (best-effort contract).');
        end

        function testSharedLabelLookupContract(testCase)
            %psyrat_copula_labels is the single label locator for the summary
            %and the provenance text, requiring inference_framework so the two
            %carriers can never disagree about whether labels exist.
            lbl = struct('inference_framework','tok','copula_feedback_to_margins',false,...
                'joint_posterior_status','stat');
            testCase.verifyEqual(psyrat_copula_labels({struct('sel',1), ...
                struct('sel',1,'labels',lbl)}).inference_framework, 'tok', ...
                'The first stratum CARRYING labels wins, not the first stratum.');
            testCase.verifyEmpty(psyrat_copula_labels({struct('sel',1)}));
            testCase.verifyEmpty(psyrat_copula_labels({struct('labels',...
                struct('joint_posterior_status','stat'))}), ...
                'Labels without inference_framework are not well-formed.');
            testCase.verifyEmpty(psyrat_copula_labels('not a cell'));
        end

    end
end

% ---------------------------------------------------------------------------
% Local helpers
% ---------------------------------------------------------------------------

function [a1,p1,i1,a2,p2,i2,cp,ci] = localCrossArgs()
%LOCALCROSSARGS Draw-vector arguments for the cross-covariance wrapper.
n = 40;
rng(3141, 'twister');
a1 = log(6.2) + 0.05*randn(n,1);  p1 = 0.32 + 0.01*rand(n,1);  i1 = 0.18 + 0.01*rand(n,1);
a2 = log(5.7) + 0.05*randn(n,1);  p2 = 0.29 + 0.01*rand(n,1);  i2 = 0.15 + 0.01*rand(n,1);
cp = 0.45 + 0.02*rand(n,1);
ci = 0.25 + 0.02*rand(n,1);
end

function f = localFixture()
%LOCALFIXTURE Synthetic stage-1 draws plus a matching ssinfo table.
%   Values are deliberately distinct per participant and per event so an
%   indexing error cannot pass by symmetry.
rng(2718, 'twister');
nd = 12;
nsub = 4;
f.ndraws = nd;
f.b = [log(6.1) + 0.02*randn(nd,1), log(5.6) + 0.02*randn(nd,1)];
f.b_sigma = [log(1.8) + 0.02*randn(nd,1), log(1.6) + 0.02*randn(nd,1)];
f.b_dim = [0.10 + 0.01*randn(nd,1), -0.06 + 0.01*randn(nd,1)];      % KDIM = 2
f.b_sigma_dim = [0.05 + 0.01*randn(nd,1), -0.03 + 0.01*randn(nd,1)];
f.sd_id = [0.30*ones(nd,1), 0.34*ones(nd,1), 0.20*ones(nd,1), 0.22*ones(nd,1)];
f.sd_trl = [0.12*ones(nd,1), 0.14*ones(nd,1)];
f.cor_id = repmat(reshape(eye(4),1,4,4), nd, 1, 1);
f.cor_id(:,1,2) = 0.45; f.cor_id(:,2,1) = 0.45;
f.cor_i = 0.25*ones(nd,1);
f.er_ss1 = log(1.8) + 0.15*repmat(1:nsub, nd, 1)/nsub;
f.er_ss2 = log(1.6) + 0.12*repmat(1:nsub, nd, 1)/nsub;
f.mu_ss1 = log(6.1) + 0.20*repmat(1:nsub, nd, 1)/nsub;
f.mu_ss2 = log(5.6) + 0.18*repmat(1:nsub, nd, 1)/nsub;
f.idtable = table((1:nsub)', (1:nsub)', linspace(-1.2, 1.2, nsub)', ...
    repmat(9, nsub, 1), repmat(9, nsub, 1), repmat(9, nsub, 1), ...
    'VariableNames', {'id','id2','z1','trls1','trls2','trls'});
f.rho = 0.3*ones(nd,1);
end

function [gen,dep,diagnostic] = localRun(f, varargin)
%LOCALRUN Calls the calculator, with optional per-call overrides.
rho = f.rho;
for k = 1:2:numel(varargin)
    if strcmp(varargin{k}, 'rho'); rho = varargin{k+1}; end
end
[gen,dep,diagnostic] = psyrat_ssrel_diffdynrel_rescor_gamma_ls( ...
    'b', f.b, 'b_sigma', f.b_sigma, 'b_dim', f.b_dim, ...
    'b_sigma_dim', f.b_sigma_dim, 'sd_id', f.sd_id, 'sd_trl', f.sd_trl, ...
    'cor_id', f.cor_id, 'cor_i', f.cor_i, ...
    'er_ss1', f.er_ss1, 'er_ss2', f.er_ss2, ...
    'mu_ss1', f.mu_ss1, 'mu_ss2', f.mu_ss2, 'rho', rho, ...
    'idtable', f.idtable, 'ndim', 1, 'CI', 0.95);
end

function fx = localLoadOnefacetStage3()
%LOCALLOADONEFACETSTAGE3 Rebuilds the one-facet calculator's argument list from
%the frozen reference fixtures, as DRAW VECTORS - the one-facet mirror of
%localLoadTwofacetStage3 in TestModularCutCopulaTrtReliability, and the same
%facet-name mapping applies to the one facet this design has:
%   reference i (trial) -> sd_trl
%The person block is the reference's sd_p_joint / cor_p_mu; the person SCALE
%SDs (sd_id columns 3/4) are zero here because the calculator consumes each
%participant's own realized log residual SD through er_ss instead.

here = fileparts(mfilename('fullpath'));
rd = @(f) readtable(fullfile(here,'baselines','modular_cut_copula',f), ...
    'VariableNamingRule','preserve');

meta   = rd('modular_cut_copula_onefacet_meta.csv');
long   = rd('modular_cut_copula_onefacet_draws.csv');
people = rd('modular_cut_copula_onefacet_persontable.csv');
pars   = rd('modular_cut_copula_onefacet_stage3_params.csv');
exp3   = rd('modular_cut_copula_onefacet_stage3_expected.csv');

%The STANDARDIZED design applies one n_i to everybody, which is what the
%calculator's shared idtable trial counts represent.
std3 = exp3(strcmp(exp3.design,'standardized'), :);
n_i_prime = std3.n_i_prime(1);

draws = unique(long.draw);
nd = numel(draws);
K = meta.K(1);   % one_dimension: a single slope per event
mu_off  = [meta.mu_offset_1(1) meta.mu_offset_2(1)];
sig_off = [meta.log_sigma_offset_1(1) meta.log_sigma_offset_2(1)];
nsub = height(people);

%population log-scale parameters, per draw
b = zeros(nd,2); bsig = zeros(nd,2);
bdim = zeros(nd,2*K); bsigdim = zeros(nd,2*K);
for d = 1:nd
    for m = 1:2
        b(d,m)    = mu_off(m)  + localLookup(long,'alpha_mu',m,NaN,draws(d));
        bsig(d,m) = sig_off(m) + localLookup(long,'alpha_sigma',m,NaN,draws(d));
        for k = 1:K
            bdim(d,2*(k-1)+m)    = localLookup(long,'beta_mu',m,k,draws(d));
            bsigdim(d,2*(k-1)+m) = localLookup(long,'beta_sigma',m,k,draws(d));
        end
    end
end

pop = @(nm,ev) arrayfun(@(d) localLookupPar(pars,nm,ev,d), draws);
sd_id = [pop('sd_p_joint',1) pop('sd_p_joint',2) zeros(nd,1) zeros(nd,1)];
corid = zeros(nd,4,4);
cp = pop('cor_p_mu',NaN);
for d = 1:nd
    corid(d,:,:) = eye(4);
    corid(d,1,2) = cp(d); corid(d,2,1) = cp(d);
end

%per-subject realized effects, as draws x NSUB
mu1 = zeros(nd,nsub); mu2 = zeros(nd,nsub);
ss1 = zeros(nd,nsub); ss2 = zeros(nd,nsub);
for d = 1:nd
    for p = 1:nsub
        pidx = people.p_index(p);
        mu1(d,p) = b(d,1)    + localLookup(long,'u_p_joint',pidx,1,draws(d));
        mu2(d,p) = b(d,2)    + localLookup(long,'u_p_joint',pidx,2,draws(d));
        ss1(d,p) = bsig(d,1) + localLookup(long,'u_p_joint',pidx,3,draws(d));
        ss2(d,p) = bsig(d,2) + localLookup(long,'u_p_joint',pidx,4,draws(d));
    end
end

%The idtable: one shared trial n' per event (the design is CONCURRENT, so both
%events are measured on the same trials and the reference carries a single n_i).
ssinfo = table(people.ID, (1:nsub)', people.Dim1_z, ...
    repmat(n_i_prime,nsub,1), repmat(n_i_prime,nsub,1), repmat(n_i_prime,nsub,1), ...
    'VariableNames',{'id','id2','z1','trls1','trls2','trls'});

fx.args = {'b',b,'b_sigma',bsig,'b_dim',bdim,'b_sigma_dim',bsigdim, ...
    'sd_id',sd_id,'cor_id',corid, ...
    'sd_trl',[pop('sd_i',1) pop('sd_i',2)], ...
    'cor_i',pop('cor_i',NaN), ...
    'er_ss1',ss1,'er_ss2',ss2,'mu_ss1',mu1,'mu_ss2',mu2, ...
    'rho',pop('rho_e',NaN), 'idtable',ssinfo, 'ndim',1};

fx.people = people;

%mean over draws of one expected column for one participant, matching the way
%psyrat_ssrel_diff summarizes (mean of the per-draw ratio).
fx.meanOver = @(col,p) mean(std3.(col)(std3.p_index == people.p_index(p)));
end

function v = localLookup(long, name, i1, i2, d)
sel = strcmp(long.param,name) & long.draw == d & long.idx1 == i1;
if ~isnan(i2); sel = sel & long.idx2 == i2; end
v = long.value(sel);
assert(isscalar(v), 'localLookup: %s[%g,%g] draw %g matched %d rows', ...
    name, i1, i2, d, numel(v));
end

function v = localLookupPar(pars, name, ev, d)
sel = strcmp(pars.param,name) & pars.draw == d;
if ~isnan(ev); sel = sel & pars.event == ev; end
v = pars.value(sel);
assert(isscalar(v), 'localLookupPar: %s (event %g) draw %g matched %d rows', ...
    name, ev, d, numel(v));
end

function REL = localFakeOnefacetREL()
%LOCALFAKEONEFACETREL A minimal analysis-19 gamma result with the REL.out layout
%psyrat_computevarcomp's case-19 gamma arm produces, built from localFixture so
%the draw values stay well conditioned. What is exercised is the summary routing
%and the disclosure contract, not the numbers. The copula stub carries the
%stage's labels verbatim and deliberately NO flag fields, so the absent-field
%counting path (older stored results) is what the disclosure test pins.
f = localFixture();
REL = struct();
REL.analysis = 'ic_diff_dynrel_sserrvar_rescor';
REL.family = 'gamma';
REL.gammascale = 2;
REL.ndim = 1;
REL.dim_names = {'dim1'};
REL.out = struct();
REL.out.labels = {'none'};
REL.out.b           = {num2cell(f.b)};
REL.out.b_sigma     = {num2cell(f.b_sigma)};
REL.out.b_dim       = {num2cell(f.b_dim)};
REL.out.b_sigma_dim = {num2cell(f.b_sigma_dim)};
REL.out.sd_id       = {num2cell(f.sd_id)};
REL.out.sd_trl      = {num2cell(f.sd_trl)};
REL.out.cor_id      = {num2cell(f.cor_id)};
REL.out.cor_i       = {num2cell(f.cor_i)};
REL.out.mu_ss1      = {num2cell(f.mu_ss1)};
REL.out.mu_ss2      = {num2cell(f.mu_ss2)};
REL.out.er_var_ss1  = {num2cell(f.er_ss1)};
REL.out.er_var_ss2  = {num2cell(f.er_ss2)};
sel = (1:2:f.ndraws)';
REL.out.copula      = {struct(...
    'sel',sel,...
    'rho_e',0.30 + (0:numel(sel)-1)' * 0.01,...
    'labels',struct(...
        'inference_framework','modular_cut_two_stage_gamma_margins_gaussian_copula_v1',...
        'copula_feedback_to_margins',false,...
        'joint_posterior_status','not_a_full_joint_bayesian_posterior'))};
REL.out.ssinfo   = {f.idtable};
REL.out.dimz     = {f.idtable.z1};
REL.out.ntrials  = 9;
REL.out.elabels  = {'e1';'e2'};
REL.out.glabels  = {'none'};
end
