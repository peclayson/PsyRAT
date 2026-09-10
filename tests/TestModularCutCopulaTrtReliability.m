classdef TestModularCutCopulaTrtReliability < PsyRATTestBase
    %TESTMODULARCUTCOPULATRTRELIABILITY Stage 3 of the modular cut-copula
    % workflow for the TWO-FACET design: the per-participant concurrent
    % difference-score calculator (analysis 20) and the psyrat_dynrel_summary
    % dispatch that reaches it.
    %
    % No CmdStan. The calculator is driven from the frozen reference fixtures and
    % from synthetic draws, so what is pinned is the conversion arithmetic, the
    % estimand choices and the wiring - not that the model recovers truth.
    %
    % THE FIXTURE COMPARISON IS THE ONLY EXTERNAL CHECK THIS DESIGN HAS. The
    % one-facet modular design can be replayed against the owner's completed
    % 151-person fits; the two-facet cannot, because the reference bundle's
    % two-facet runs were prepared against a simulated second occasion and were
    % never sampled. See TestGammaDiffDynrelTrtLsConversion's header.
    %
    % The fixture's 'standardized' D-study rows are the ones compared, not
    % 'actual': psyrat_ssrel_diff_trt applies ONE occasion n' to every
    % participant (matching the population surface), which is the standardized
    % design's contract. The actual-design rows use each person's own n_o and
    % have no counterpart in this calculator's signature.

    methods (Test)

        %% -- the reference oracle ---------------------------------------------

        function testPerParticipantCoefficientsMatchTheReferenceOracle(testCase)
            %Every reported coefficient, per participant, against the reference's
            %own reliability_twofacet. Both sides average a per-draw ratio over
            %draws, so the comparison is mean-to-mean rather than ratio-of-means.
            fx = localLoadTwofacetStage3();

            for rel = [3 1 2]   % CES, CE, CS
                [gen,dep] = psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls( ...
                    fx.args{:}, 'reltype', rel, 'nocc', fx.n_o_prime, 'CI', 0.95);

                switch rel
                    case 3; gname = 'G_CES_delta';
                    case 1; gname = 'G_CE_delta';
                    case 2; gname = 'G_CS_delta';
                end
                for p = 1:height(fx.people)
                    want = fx.meanOver(gname, p);
                    testCase.verifyEqual(gen.rel_pt(p), want, 'RelTol', 1e-7, ...
                        sprintf('participant %d, reltype %d: %s', p, rel, gname));
                end

                if rel == 3
                    %the CES rows are the only ones the reference reports a
                    %dependability, ICC and SEM for, so those are checked once.
                    for p = 1:height(fx.people)
                        testCase.verifyEqual(dep.rel_pt(p), ...
                            fx.meanOver('D_CES_delta', p), 'RelTol', 1e-7, ...
                            sprintf('participant %d: D_CES_delta', p));
                        testCase.verifyEqual(gen.icc_pt(p), ...
                            fx.meanOver('ICC_relative_CES_delta', p), 'RelTol', 1e-7, ...
                            sprintf('participant %d: relative ICC', p));
                        testCase.verifyEqual(dep.icc_pt(p), ...
                            fx.meanOver('ICC_absolute_CES_delta', p), 'RelTol', 1e-7, ...
                            sprintf('participant %d: absolute ICC', p));
                        testCase.verifyEqual(gen.sem_pt(p), ...
                            fx.meanOver('SEM_relative_delta', p), 'RelTol', 1e-7, ...
                            sprintf('participant %d: relative SEM', p));
                        testCase.verifyEqual(dep.sem_pt(p), ...
                            fx.meanOver('SEM_absolute_delta', p), 'RelTol', 1e-7, ...
                            sprintf('participant %d: absolute SEM', p));
                    end
                end
            end
        end

        function testDiagnosticReportsTheReferenceResidualCovariance(testCase)
            %The per-participant quadrature provenance must describe the same
            %quantity the reference computed, and must report the observed-scale
            %residual CORRELATION separately from the copula correlation - they
            %are different numbers and the table exists so they are not confused.
            fx = localLoadTwofacetStage3();
            [~,~,diag] = psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls( ...
                fx.args{:}, 'reltype', 3, 'nocc', fx.n_o_prime, 'CI', 0.95);

            for p = 1:height(fx.people)
                testCase.verifyEqual(diag.res_cov_pt(p), ...
                    fx.meanOver('residual_covariance_person', p), 'RelTol', 1e-7, ...
                    sprintf('participant %d: residual covariance', p));
                testCase.verifyEqual(diag.res_cor_pt(p), ...
                    fx.meanOver('residual_correlation_person', p), 'RelTol', 1e-7, ...
                    sprintf('participant %d: residual correlation', p));
            end
            testCase.verifyEqual(diag.quad_unconverged_n, ...
                zeros(height(fx.people),1), ...
                'The reference fixture should not need the non-strict fallback.');
            testCase.verifyNotEqual(diag.res_cor_pt, diag.rho_e_pt, ...
                'Observed-scale residual correlation is NOT the copula correlation.');
        end

        %% -- the estimand claims ------------------------------------------------

        function testPersonByTrialStaysOutOfTheResidualChannel(testCase)
            %THE CENTRAL CLAIM. The residual channel handed to the coefficient is
            %the three-way term plus dispersion, NOT the one-facet's pooled
            %person x trial + error. Verified against the reference's own
            %var_theta_pio_e, which excludes person x trial by construction
            %(R/reliability.R:246, 266).
            fx = localLoadTwofacetStage3();
            s = fx.perDraw(1,1);   % draw 1, participant 1

            vc = psyrat_gamma_varcomps_trt_ls(s.alpha_a, s.sp_a, s.so_a, ...
                s.st_a, s.spo_a, s.spt_a, s.sot_a, log(s.sigma_a), 0);

            testCase.verifyEqual(vc.sigma_res2, s.var_theta_pio_e, 'RelTol', 1e-9, ...
                'The residual channel must be the reference pio_e.');
            testCase.verifyGreaterThan(vc.sigma_pt2, 0, ...
                'Person x trial must be a separate, nonzero component here.');
            testCase.verifyEqual(vc.sigma_pt2, s.var_theta_pi, 'RelTol', 1e-9, ...
                'and must equal the reference person-by-trial component.');
        end

        function testTheParticipantsOwnLocationEffectEntersThroughTheError(testCase)
            %mu_ss carries each participant's realized location effect and sets
            %the point at which the copula covariance is evaluated. It must move
            %the result when rho is nonzero, and go INERT at rho = 0 - which is
            %what shows it enters through the error rather than the signal.
            a = localSyntheticArgs();

            shifted = a;
            i1 = find(strcmp(shifted,'mu_ss1'),1);
            shifted{i1+1} = shifted{i1+1} + 0.4;

            base = psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls(a{:});
            moved = psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls(shifted{:});
            testCase.verifyNotEqual(moved.rel_pt, base.rel_pt, ...
                'A shifted person location must move the coefficient at rho ~= 0.');

            nd = localNDraws();
            z0 = localReplace(a, 'rho', zeros(nd,1));
            z0shift = localReplace(shifted, 'rho', zeros(nd,1));
            testCase.verifyEqual( ...
                psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls(z0shift{:}).rel_pt, ...
                psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls(z0{:}).rel_pt, ...
                'AbsTol', 0, ...
                'At rho = 0 the quadrature vanishes, so mu_ss must be inert.');
        end

        %% -- input contract -----------------------------------------------------

        function testRejectsMisalignedOrInadmissibleDraws(testCase)
            %Draw alignment is the caller's job: rho exists only for the retained
            %stage-1 draws. Recycling a short rho against full-length margins
            %would produce a plausible number from an incoherent posterior.
            a = localSyntheticArgs();
            nd = localNDraws();

            short = localReplace(a, 'rho', 0.2*ones(nd-1,1));
            testCase.verifyError( ...
                @() psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls(short{:}), ...
                'varargin:rho');

            atOne = localReplace(a, 'rho', ones(nd,1));
            testCase.verifyError( ...
                @() psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls(atOne{:}), ...
                'varargin:rho');
        end

        function testRejectsShortPerSubjectMatrices(testCase)
            %psyrat_ssrel_diff_trt silently falls back to the POPULATION residual
            %for an id2 it cannot index, which would turn a subject-level table
            %into a population one without any error. Fail loudly instead.
            a = localSyntheticArgs();
            ss1 = localArg(a,'er_ss1');
            bad1 = localReplace(a, 'er_ss1', ss1(:,1));
            testCase.verifyError( ...
                @() psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls(bad1{:}), ...
                'varargin:er_ss');

            mu1 = localArg(a,'mu_ss1');
            bad2 = localReplace(a, 'mu_ss1', mu1(:,1));
            testCase.verifyError( ...
                @() psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls(bad2{:}), ...
                'varargin:mu_ss');
        end

        function testRejectsAPerEventFacetGivenOnce(testCase)
            a = localSyntheticArgs();
            bad = localReplace(a, 'sd_oid', 0.07*ones(localNDraws(),1));
            testCase.verifyError( ...
                @() psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls(bad{:}), ...
                'varargin:facetsd');
        end

        %% -- the summary dispatch -----------------------------------------------

        function testSummaryReachesTheTwoFacetGammaBranch(testCase)
            %REGRESSION GUARD for the highest-risk wiring here. Until this
            %increment psyrat_dynrel_summary refused analysis 20 under gamma with
            %dynrel:modularcutpending. With the guard gone, the branch order is
            %what keeps a gamma result out of the Gaussian two-facet branch,
            %which reads *_varcov / rescor / err_varcov slots the gamma arm never
            %fills. Producing a real surface AND a per-participant table is only
            %reachable through the new branch.
            REL = localFakeTwofacetREL();

            s = psyrat_dynrel_summary(REL,'CI',0.95,'ngrid',4);

            testCase.verifyEqual(numel(s.strata), 1);
            testCase.verifySize(s.strata(1).G.pt, [1 4], ...
                'Two-facet gamma surface has the wrong grid shape.');
            testCase.verifyTrue(all(isfinite(s.strata(1).G.pt)) && ...
                all(s.strata(1).G.pt > 0 & s.strata(1).G.pt < 1), ...
                'Two-facet gamma G(z) is not a coefficient.');
            testCase.verifyNotEmpty(s.strata(1).varcomp, ...
                'Two-facet gamma variance-component table is empty.');
            testCase.verifyNotEmpty(s.strata(1).ssrel_table, ...
                'analysis 20 must produce a per-participant table.');
            testCase.verifyNotEmpty(s.strata(1).quaddiag, ...
                'analysis 20 must surface the quadrature provenance.');
            testCase.verifyEqual(height(s.strata(1).ssrel_table), 3, ...
                'One row per participant.');
        end

        function testSummarySubsamplesToTheRetainedCopulaDraws(testCase)
            %Every quantity must be computed on the retained draw set, surface
            %included - the owner's decision, taken so the surface and the table
            %stay per-draw consistent instead of the surface silently using draws
            %the copula never saw. Shrinking sel must therefore change the
            %surface, not just the table.
            REL = localFakeTwofacetREL();
            %Constant rho for THIS test only: with the stub's varying rho_e,
            %truncating the vector changes the surface through rho alone and
            %the assertion stops isolating the margin-draw subsetting it
            %exists to pin (2026-08-11 code review).
            REL.out.copula{1}.rho_e(:) = 0.3;
            full = psyrat_dynrel_summary(REL,'CI',0.95,'ngrid',4);

            REL2 = REL;
            REL2.out.copula{1}.sel = REL.out.copula{1}.sel(1:2);
            REL2.out.copula{1}.rho_e = REL.out.copula{1}.rho_e(1:2);
            cut = psyrat_dynrel_summary(REL2,'CI',0.95,'ngrid',4);

            testCase.verifyNotEqual(cut.strata(1).G.pt, full.strata(1).G.pt, ...
                'The SURFACE must be recomputed on the retained draws too.');
        end

        function testSummaryRejectsBarredTwoFacetGammaSiblings(testCase)
            %The two-facet gamma location-scale branch supports only analysis
            %20; its barred siblings satisfy the same branch condition and
            %would otherwise be summarized through the analysis-20 layout -
            %a missing-field error at best, numbers computed from the wrong
            %quantities at worst (guard hardening, 2026-08-11 post-merge
            %follow-up). Each must fail by name.
            barred = {'ic_diff_dynrel_trt', 'ic_diff_dynrel_sserrvar_trt', ...
                'ic_diff_dynrel_trt_rescor', 'ic_diff_dynrel_sserrvar_trt_rho'};
            for k = 1:numel(barred)
                REL = localFakeTwofacetREL();
                REL.analysis = barred{k};
                testCase.verifyError( ...
                    @() psyrat_dynrel_summary(REL,'CI',0.95,'ngrid',4), ...
                    'dynrel:gammatrtdesign', barred{k});
            end
        end

        function testSummaryStillFailsLoudlyWithoutTheCopulaStage(testCase)
            %The pending guard is gone, but a gamma two-facet result with no
            %stage-2 output is still not summarizable. It must fail rather than
            %fall through to a branch that reads slots this arm never fills.
            REL = localFakeTwofacetREL();
            REL.out = rmfield(REL.out,'copula');
            testCase.verifyError(@() psyrat_dynrel_summary(REL,'CI',0.95,'ngrid',4), ...
                'MATLAB:nonExistentField');
        end

        function testSummaryCarriesTheCutPosteriorDisclosure(testCase)
            %The stage's provenance labels, its diagnostic flag counts, and the
            %copula correlation itself must reach the summary - before this
            %increment all three died in REL.out (pre-merge review finding: the
            %disclosure chain was broken while the committed docs claimed the
            %labels rode through every output surface). A clean stage-2 (no
            %flag fields at all - the best-effort path) must not warn.
            REL = localFakeTwofacetREL();
            s = testCase.verifyWarningFree(...
                @() psyrat_dynrel_summary(REL,'CI',0.95,'ngrid',4));

            %labels ride verbatim from the stage output (single source of truth)
            testCase.verifyEqual(s.copula_labels.inference_framework, ...
                'modular_cut_two_stage_gamma_margins_gaussian_copula_v1');
            testCase.verifyFalse(s.copula_labels.copula_feedback_to_margins);
            testCase.verifyEqual(s.copula_labels.joint_posterior_status, ...
                'not_a_full_joint_bayesian_posterior');

            %flag counts: 4 retained draws, absent flag fields count as zero
            testCase.verifyEqual(s.strata(1).copula_flags.n_draws, 4);
            testCase.verifyEqual(s.strata(1).copula_flags.n_invalid, 0);
            testCase.verifyEqual(s.strata(1).copula_flags.n_calibration, 0);

            %rho_e row: exactly one, summarizing the RETAINED draws
            vc = s.strata(1).varcomp;
            row = strcmp(vc.symbol,'rho_e,12(copula)');
            testCase.verifyEqual(sum(row), 1, ...
                'Exactly one copula-correlation row in the varcomp table.');
            rc = REL.out.copula{1}.rho_e;
            testCase.verifyEqual(vc.estimate(row), mean(rc), 'AbsTol', 1e-12);
            testCase.verifyEqual(vc.ci_lower(row), quantile(rc,0.025), ...
                'AbsTol', 1e-12);
            testCase.verifyEqual(vc.ci_upper(row), quantile(rc,0.975), ...
                'AbsTol', 1e-12);
            testCase.verifyEqual(vc.component{row}, ...
                'Copula residual correlation (events)');
        end

        function testSummaryWarnsOnStageInvalidityFlags(testCase)
            %The stage records per-draw numerical invalidity but never blocks
            %(diagnose-not-hide); the summary is where the user hears about it.
            %Calibration flags are descriptive model checks and must be counted
            %WITHOUT warning.
            REL = localFakeTwofacetREL();
            REL.out.copula{1}.invalidity_flag = logical([1;0;1;0]);
            s = testCase.verifyWarning(...
                @() psyrat_dynrel_summary(REL,'CI',0.95,'ngrid',4), ...
                'dynrel:copulainvalidity');
            testCase.verifyEqual(s.strata(1).copula_flags.n_invalid, 2);

            REL2 = localFakeTwofacetREL();
            REL2.out.copula{1}.calibration_flag = logical([0;1;0;0]);
            s2 = testCase.verifyWarningFree(...
                @() psyrat_dynrel_summary(REL2,'CI',0.95,'ngrid',4));
            testCase.verifyEqual(s2.strata(1).copula_flags.n_calibration, 1);
            testCase.verifyEqual(s2.strata(1).copula_flags.n_invalid, 0);
        end

    end
end

%% ------------------------------------------------------------------------
%% helpers
%% ------------------------------------------------------------------------

function n = localNDraws()
n = 6;
end

function a = localSyntheticArgs()
%LOCALSYNTHETICARGS A small two-facet calculator call: 3 participants, one
%dimension, all facet SDs distinct across events so an argument slip cannot pass
%by symmetry.
n = localNDraws();
nsub = 3;
step = @(v) v + (0:n-1)' * (v/40);
rep = @(v) repmat(v, 1, nsub) + (0:nsub-1) * 0.01;
corid = zeros(n,4,4);
for q = 1:n
    corid(q,:,:) = eye(4);
    corid(q,1,2) = 0.45; corid(q,2,1) = 0.45;
end
ssinfo = table((1:nsub)', (1:nsub)', linspace(-1,1,nsub)', ...
    repmat(20,nsub,1), repmat(20,nsub,1), repmat(20,nsub,1), ...
    repmat(2,nsub,1), repmat(2,nsub,1), ...
    'VariableNames',{'id','id2','z1','trls1','trls2','trls','occs1','occs2'});

a = {'b',[step(1.62) step(1.44)], 'b_sigma',[step(-0.5) step(-0.62)], ...
    'b_dim',[step(0.10) step(0.08)], 'b_sigma_dim',[step(-0.05) step(-0.04)], ...
    'sd_id',[step(0.301) step(0.283) step(0.201) step(0.184)], ...
    'sd_trl',[step(0.121) step(0.113)], 'sd_occ',[step(0.104) step(0.096)], ...
    'sd_tid',[step(0.093) step(0.086)], 'sd_oid',[step(0.082) step(0.074)], ...
    'sd_to',[step(0.061) step(0.052)], 'cor_id',corid, ...
    'cor_t',step(0.31), 'cor_o',step(0.24), 'cor_pt',step(0.19), ...
    'cor_po',step(-0.14), 'cor_ot',step(0.27), ...
    'er_ss1',rep(step(-0.5)), 'er_ss2',rep(step(-0.62)), ...
    'mu_ss1',rep(step(1.62)), 'mu_ss2',rep(step(1.44)), ...
    'rho',0.3*ones(n,1), 'idtable',ssinfo, 'ndim',1, ...
    'reltype',3, 'nocc',2, 'CI',.95};
end

function REL = localFakeTwofacetREL()
%LOCALFAKETWOFACETREL A minimal analysis-20 gamma result with the exact REL.out
%layout psyrat_computevarcomp's case-20 gamma arm produces. Draw values are
%arbitrary but well conditioned; what is exercised is the routing and the field
%contract, not the numbers.
n = 8;
nsub = 3;
step = @(v) v + (0:n-1)' * (v/40);
rep = @(v) repmat(v, 1, nsub) + (0:nsub-1) * 0.01;
corid = zeros(n,4,4);
for q = 1:n
    corid(q,:,:) = eye(4);
    corid(q,1,2) = 0.45; corid(q,2,1) = 0.45;
end

REL = struct();
REL.analysis = 'ic_diff_dynrel_sserrvar_trt_rescor';
REL.family = 'gamma';
REL.gammascale = 2;
REL.ndim = 1;
REL.dim_names = {'dim1'};
REL.out = struct();
REL.out.labels = {'none'};
REL.out.b           = {num2cell([step(1.62) step(1.44)])};
REL.out.b_sigma     = {num2cell([step(-0.5) step(-0.62)])};
REL.out.b_dim       = {num2cell([step(0.10) step(0.08)])};
REL.out.b_sigma_dim = {num2cell([step(-0.05) step(-0.04)])};
REL.out.sd_id       = {num2cell([step(0.301) step(0.283) step(0.201) step(0.184)])};
REL.out.cor_id      = {num2cell(corid)};
REL.out.sd_trl      = {num2cell([step(0.121) step(0.113)])};
REL.out.sd_occ      = {num2cell([step(0.104) step(0.096)])};
REL.out.sd_tid      = {num2cell([step(0.093) step(0.086)])};
REL.out.sd_oid      = {num2cell([step(0.082) step(0.074)])};
REL.out.sd_to       = {num2cell([step(0.061) step(0.052)])};
REL.out.cor_t       = {num2cell(step(0.31))};
REL.out.cor_o       = {num2cell(step(0.24))};
REL.out.cor_pt      = {num2cell(step(0.19))};
REL.out.cor_po      = {num2cell(step(-0.14))};
REL.out.cor_ot      = {num2cell(step(0.27))};
REL.out.mu_ss1      = {num2cell(rep(step(1.62)))};
REL.out.mu_ss2      = {num2cell(rep(step(1.44)))};
REL.out.er_var_ss1  = {num2cell(rep(step(-0.5)))};
REL.out.er_var_ss2  = {num2cell(rep(step(-0.62)))};
%rho_e varies across the retained draws so the table row's mean and CI are
%non-vacuous; no flag fields by default (the counting is best-effort and the
%absent-field path must read as zero flagged, not error), labels verbatim as
%the stage archives them.
REL.out.copula      = {struct(...
    'sel',(1:2:n)',...
    'rho_e',0.30 + (0:numel(1:2:n)-1)' * 0.01,...
    'labels',struct(...
        'inference_framework','modular_cut_two_stage_gamma_margins_gaussian_copula_v1',...
        'copula_feedback_to_margins',false,...
        'joint_posterior_status','not_a_full_joint_bayesian_posterior'))};
REL.out.ssinfo      = {table((1:nsub)', (1:nsub)', linspace(-1,1,nsub)', ...
    repmat(20,nsub,1), repmat(20,nsub,1), repmat(20,nsub,1), ...
    repmat(2,nsub,1), repmat(2,nsub,1), ...
    'VariableNames',{'id','id2','z1','trls1','trls2','trls','occs1','occs2'})};
REL.out.dimz     = {linspace(-1,1,nsub)'};
REL.out.ntrials  = 20;
REL.out.nocc     = 2;
REL.out.elabels  = {'e1';'e2'};
REL.out.glabels  = {'none'};
end

function v = localArg(args, name)
idx = find(strcmp(args, name), 1);
assert(~isempty(idx), 'localArg: %s not present', name);
v = args{idx+1};
end

function args = localReplace(args, name, value)
idx = find(strcmp(args, name), 1);
assert(~isempty(idx), 'localReplace: %s not present', name);
args{idx+1} = value;
end

function fx = localLoadTwofacetStage3()
%LOCALLOADTWOFACETSTAGE3 Rebuilds the calculator's full argument list from the
%frozen reference fixtures, as DRAW VECTORS (the calculator consumes all draws at
%once, where TestGammaDiffDynrelTrtLsConversion consumes one draw at a time).
%
%FACET NAMES CROSS OVER between the reference and the toolbox, and the fixture
%values are deliberately all-distinct so a mis-mapping fails rather than passing
%by symmetry:
%   reference i  (trial)             -> sd_trl
%   reference o  (occasion)          -> sd_occ
%   reference pi (person x trial)    -> sd_tid
%   reference po (person x occasion) -> sd_oid
%   reference io (trial x occasion)  -> sd_to

here = fileparts(mfilename('fullpath'));
rd = @(f) readtable(fullfile(here,'baselines','modular_cut_copula',f), ...
    'VariableNamingRule','preserve');

meta   = rd('modular_cut_copula_twofacet_meta.csv');
long   = rd('modular_cut_copula_twofacet_draws.csv');
people = rd('modular_cut_copula_twofacet_persontable.csv');
pars   = rd('modular_cut_copula_twofacet_stage3_params.csv');
exp3   = rd('modular_cut_copula_twofacet_stage3_expected.csv');

%The STANDARDIZED design applies one (n_i, n_o) to everybody, which is what
%psyrat_ssrel_diff_trt's shared-nocc signature can represent. See the header.
std3 = exp3(strcmp(exp3.design,'standardized'), :);
n_i_prime = std3.n_i_prime(1);
n_o_prime = std3.n_o_prime(1);

draws = unique(long.draw);
nd = numel(draws);
K = meta.K(1);
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
ssinfo = table(people.ID, (1:nsub)', people.Dim1_z, people.Dim2_z, ...
    repmat(n_i_prime,nsub,1), repmat(n_i_prime,nsub,1), repmat(n_i_prime,nsub,1), ...
    repmat(n_o_prime,nsub,1), repmat(n_o_prime,nsub,1), ...
    'VariableNames',{'id','id2','z1','z2','trls1','trls2','trls','occs1','occs2'});

fx.args = {'b',b,'b_sigma',bsig,'b_dim',bdim,'b_sigma_dim',bsigdim, ...
    'sd_id',sd_id,'cor_id',corid, ...
    'sd_trl',[pop('sd_i',1) pop('sd_i',2)], ...
    'sd_occ',[pop('sd_o',1) pop('sd_o',2)], ...
    'sd_tid',[pop('sd_pi',1) pop('sd_pi',2)], ...
    'sd_oid',[pop('sd_po',1) pop('sd_po',2)], ...
    'sd_to', [pop('sd_io',1) pop('sd_io',2)], ...
    'cor_t',pop('cor_i',NaN), 'cor_o',pop('cor_o',NaN), ...
    'cor_pt',pop('cor_pi',NaN), 'cor_po',pop('cor_po',NaN), ...
    'cor_ot',pop('cor_io',NaN), ...
    'er_ss1',ss1,'er_ss2',ss2,'mu_ss1',mu1,'mu_ss2',mu2, ...
    'rho',pop('rho_e',NaN), 'idtable',ssinfo, 'ndim',2};

fx.people = people;
fx.n_i_prime = n_i_prime;
fx.n_o_prime = n_o_prime;
fx.draws = draws;

%mean over draws of one expected column for one participant, matching the way
%psyrat_diffrel_trt summarizes (mean of the per-draw ratio).
fx.meanOver = @(col,p) mean(std3.(col)(std3.p_index == people.p_index(p)));

%one (draw, participant) slice, for the component-level checks
fx.perDraw = @(d,p) localSlice(std3, people, pars, d, p);
end

function s = localSlice(std3, people, pars, d, p)
%LOCALSLICE The reference's own component values for one draw and participant,
%plus the population SDs needed to recompute them.
row = std3(std3.posterior_row == d & std3.p_index == people.p_index(p), :);
assert(height(row) == 1, 'localSlice matched %d rows', height(row));
s = struct();
s.sp_a  = localLookupPar(pars,'sd_p_joint',1,d);
s.st_a  = localLookupPar(pars,'sd_i',1,d);
s.so_a  = localLookupPar(pars,'sd_o',1,d);
s.spt_a = localLookupPar(pars,'sd_pi',1,d);
s.spo_a = localLookupPar(pars,'sd_po',1,d);
s.sot_a = localLookupPar(pars,'sd_io',1,d);
s.sigma_a = row.residual_sd_theta_person;
%recover the population alpha at this participant's z from the reference's own
%population expected value, so this slice does not restate the linear predictor.
v_total = s.sp_a^2 + s.st_a^2 + s.so_a^2 + s.spt_a^2 + s.spo_a^2 + s.sot_a^2;
s.alpha_a = log(row.expected_theta_population) - 0.5*v_total;
s.var_theta_pio_e = row.var_theta_pio_e;
s.var_theta_pi = row.var_theta_pi;
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
