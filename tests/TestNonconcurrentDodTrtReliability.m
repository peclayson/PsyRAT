classdef TestNonconcurrentDodTrtReliability < matlab.unittest.TestCase
    %TESTNONCONCURRENTDODTRTRELIABILITY The two-facet per-participant stage-3
    % calculator (analysis 29, PSYRAT_SSREL_DODIFFDYNREL_TRT_GAMMA_LS) against
    % the FROZEN EXTERNAL ORACLE: the reference bundle's own stage-3 output on
    % the deterministic two-facet lane (synthetic occasions; the reference's
    % own two-facet fits were never sampled on real data, so this frozen lane
    % plus live CmdStan recovery are the two-facet path's external checks).
    %
    % Beyond the one-facet class this pins THE DECOMPOSITION MAPPING: the
    % reference's CES/CE/CS families against reltype 3/1/2, on all-distinct
    % facet values so a swapped facet or decomposition fails rather than
    % passing by symmetry - including the property that CE/CS universes carry
    % their scaled facet term and therefore differ between the actual and
    % standardized designs.
    %
    % Tolerances as the one-facet class: RelTol 1e-9, AbsTol 1e-12 floor;
    % flags exact.

    methods (Test)

        function testDecompositionsMatchFrozenReference(testCase)
            fx = localLoadTrtLane();
            [~, dl] = localRunTrt(fx, 3);   %drawlevel carries all three
            act = fx.actual; stdz = fx.standardized;
            flat = @(m) m(:);

            pairs = {'ces','CES'; 'ce','CE'; 'cs','CS'};
            for i = 1:3
                lo = pairs{i,1}; up = pairs{i,2};
                testCase.verifyEqual(flat(dl.(sprintf('u_%s',lo))'), ...
                    act.(sprintf('universe_%s_dod',up)), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-12, ...
                    sprintf('universe %s (actual)', up));
                testCase.verifyEqual(flat(dl.(sprintf('rel_%s',lo))'), ...
                    act.(sprintf('relative_error_%s_dod',up)), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-12, ...
                    sprintf('relative error %s (actual)', up));
                testCase.verifyEqual(flat(dl.(sprintf('abs_%s',lo))'), ...
                    act.(sprintf('absolute_error_%s_dod',up)), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-12, ...
                    sprintf('absolute error %s (actual)', up));
                testCase.verifyEqual(flat(dl.(sprintf('G_%s',lo))'), ...
                    act.(sprintf('G_%s_dod',up)), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-12, sprintf('G %s (actual)', up));
                testCase.verifyEqual(flat(dl.(sprintf('D_%s',lo))'), ...
                    act.(sprintf('D_%s_dod',up)), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-12, sprintf('D %s (actual)', up));
                testCase.verifyEqual(flat(dl.(sprintf('G_%s_std',lo))'), ...
                    stdz.(sprintf('G_%s_dod',up)), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-12, ...
                    sprintf('G %s (standardized)', up));
                testCase.verifyEqual(flat(dl.(sprintf('D_%s_std',lo))'), ...
                    stdz.(sprintf('D_%s_dod',up)), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-12, ...
                    sprintf('D %s (standardized)', up));
                testCase.verifyEqual(flat(dl.(sprintf('iccg_%s',lo))'), ...
                    act.(sprintf('ICC_relative_%s_matched_facet_dod',up)), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-12, sprintf('relative ICC %s', up));
                testCase.verifyEqual(flat(dl.(sprintf('iccd_%s',lo))'), ...
                    act.(sprintf('ICC_absolute_%s_matched_facet_dod',up)), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-12, sprintf('absolute ICC %s', up));
            end

            testCase.verifyEqual(flat(dl.invalid'), ...
                act.posterior_invalid_primary_any, 'AbsTol', 0, ...
                'invalid flag (actual)');
            testCase.verifyEqual(flat(dl.invalid_std'), ...
                stdz.posterior_invalid_primary_any, 'AbsTol', 0, ...
                'invalid flag (standardized)');

            %CE/CS universes carry a scaled term, so they must genuinely
            %differ between designs on this fixture (n' primes != actual);
            %the CES universe, pure person variance, must NOT move.
            t1 = localRunTrt(fx, 1);
            testCase.verifyTrue(all(t1.uni_pt ~= t1.uni_pt_std), ...
                'The CE universe must move with the D-study design.');
            t3 = localRunTrt(fx, 3);
            testCase.verifyEqual(t3.uni_pt, t3.uni_pt_std, 'AbsTol', 0, ...
                'The CES universe is design-invariant by construction.');
        end

        function testReltypeSelectsTheMatchingDecomposition(testCase)
            %The reltype argument must surface exactly the mapped
            %decomposition in the primary table columns: 1 = CE, 2 = CS,
            %3 = CES (the toolbox reltype convention of psyrat_diffrel_trt).
            fx = localLoadTrtLane();
            [t3, dl] = localRunTrt(fx, 3);
            t1 = localRunTrt(fx, 1);
            t2 = localRunTrt(fx, 2);
            testCase.verifyEqual(t3.gen_pt, mean(dl.G_ces, 1, 'omitnan')', ...
                'AbsTol', 1e-12, 'reltype 3 must surface CES.');
            testCase.verifyEqual(t1.gen_pt, mean(dl.G_ce, 1, 'omitnan')', ...
                'AbsTol', 1e-12, 'reltype 1 must surface CE.');
            testCase.verifyEqual(t2.gen_pt, mean(dl.G_cs, 1, 'omitnan')', ...
                'AbsTol', 1e-12, 'reltype 2 must surface CS.');
            testCase.verifyTrue(all(t1.gen_pt ~= t2.gen_pt), ...
                'CE and CS must differ on all-distinct facet values.');
        end

        function testPersonQuantitiesMatchFrozenReference(testCase)
            fx = localLoadTrtLane();
            [~, dl] = localRunTrt(fx, 3);
            act = fx.actual;
            labels = fx.labels;
            for q = 1:4
                testCase.verifyEqual(localFlatCell(dl.expected_person, q), ...
                    act.(sprintf('expected_%s_person',labels{q})), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-12, ...
                    sprintf('expected_%s_person', labels{q}));
                testCase.verifyEqual(localFlatCell(dl.expected_pop, q), ...
                    act.(sprintf('expected_%s_population',labels{q})), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-12, ...
                    sprintf('expected_%s_population', labels{q}));
                testCase.verifyEqual(localFlatCell(dl.sigma_person, q), ...
                    act.(sprintf('residual_sd_%s_person',labels{q})), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-12, ...
                    sprintf('residual_sd_%s_person', labels{q}));
                testCase.verifyEqual(localFlatCell(dl.nu_person, q), ...
                    act.(sprintf('nu_%s_at_person_expected_mean',labels{q})), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-12, sprintf('nu_%s', labels{q}));
            end
            s4 = localFlatCell(dl.sigma_person.^2, 1) + ...
                localFlatCell(dl.sigma_person.^2, 2) + ...
                localFlatCell(dl.sigma_person.^2, 3) + ...
                localFlatCell(dl.sigma_person.^2, 4);
            testCase.verifyEqual(s4, act.residual_var_dod_person, ...
                'RelTol', 1e-9, 'AbsTol', 1e-12, 'residual sum of four');
        end

        function testCouplingVocabularyRejected(testCase)
            for coupled = {'rho','rho_e','cov_e_obs','cov_res_e_obs','er_cov','rescor','copula'}
                testCase.verifyError(@() psyrat_ssrel_dodiffdynrel_trt_gamma_ls( ...
                    coupled{1}, 0.3), ...
                    'psyrat_ssrel_dodiffdynrel_trt_gamma_ls:noresidualcoupling');
                testCase.verifyError(@() psyrat_rel_dodiffdynrel_trt_gamma_ls( ...
                    coupled{1}, 0.3), ...
                    'psyrat_rel_dodiffdynrel_trt_gamma_ls:noresidualcoupling');
            end
        end

        function testSummaryDispatch(testCase)
            %A synthetic analysis-29 REL through psyrat_dynrel_summary: the
            %two-facet DoD branch must route with the per-cell trial AND
            %occasion n' vectors, honor the reltype selector, and carry the
            %labels/flags.
            fx = localLoadTrtLane();
            REL = localSyntheticTrtRel(fx);
            w = warning('off','dynrel:dodinvalidity');
            cleaner = onCleanup(@() warning(w));
            summ = psyrat_dynrel_summary(REL,'ngrid',4,'reltype',1);

            st = summ.strata(1);
            testCase.verifyEqual(numel(st.obs_ev), 4);
            testCase.verifyEqual(st.nocc, 1, ...
                ['Default occasion n'' is 1 (the owner decision every trt '...
                'dynrel variant follows and what the GUI''s "1 occasion" '...
                'choice means).']);
            summObs = psyrat_dynrel_summary(REL,'ngrid',4,'reltype',1,...
                'nocc','observed');
            testCase.verifyEqual(summObs.strata(1).nocc, 2, ...
                '''observed'' must resolve to the stored per-cell medians.');
            testCase.verifyError(@() psyrat_dynrel_summary(REL,'ngrid',4,...
                'reltype',1,'nocc',NaN), 'dynrel:dodnocc', ...
                'A NaN occasion override must be rejected, not floored to 1.');
            testCase.verifyEqual(height(st.ssrel_table), fx.P);
            testCase.verifyEqual(unique(st.ssrel_table.reltype), 1, ...
                'The reltype selector must reach the per-participant table.');
            testCase.verifyEqual(summ.dod_labels.estimand_version, ...
                'signed_four_component_nonconcurrent_dod_v1');
            testCase.verifyTrue(st.dod_flags.recorded);

            %reltype must move the surface (CE vs CES differ on all-distinct
            %facet values)
            summ3 = psyrat_dynrel_summary(REL,'ngrid',4,'reltype',3);
            testCase.verifyTrue(any(st.G.pt(:) ~= summ3.strata(1).G.pt(:)), ...
                'reltype 1 and 3 must produce different surfaces.');
        end

        function testSurfaceMatchesZeroEffectParticipant(testCase)
            %The typical-person surface at grid (z1,z2) must equal the
            %standardized-design coefficients of a zero-person-effect
            %participant at the same point, for a universe-moving reltype (CE)
            %and the pure-person one (CES) alike.
            fx = localLoadTrtLane();
            z1p = 0.45; z2p = -0.7;
            for reltype = [1 3]
                surf = psyrat_rel_dodiffdynrel_trt_gamma_ls('b', fx.b, ...
                    'b_sigma', fx.b_sigma, 'b_dim', fx.b_dim, ...
                    'b_sigma_dim', fx.b_sigma_dim, 'sd_id', fx.sd_id, ...
                    'cor_id', fx.cor_id, 'sd_trl', fx.sd_trl, ...
                    'cor_trl', fx.cor_trl, 'sd_occ', fx.sd_occ, ...
                    'cor_occ', fx.cor_occ, 'sd_tid', fx.sd_tid, ...
                    'cor_tid', fx.cor_tid, 'sd_oid', fx.sd_oid, ...
                    'cor_oid', fx.cor_oid, 'sd_to', fx.sd_to, ...
                    'cor_to', fx.cor_to, 'ndim', 2, 'z1', z1p, 'z2', z2p, ...
                    'obs', fx.nprime, 'nocc', fx.nprime_o, ...
                    'reltype', reltype, 'CI', 0.95, 'cvec', fx.cvec);

                fake = table(1, 1, z1p, z2p, ...
                    fx.idtable.trls1(1), fx.idtable.trls2(1), ...
                    fx.idtable.trls3(1), fx.idtable.trls4(1), ...
                    fx.idtable.occs1(1), fx.idtable.occs2(1), ...
                    fx.idtable.occs3(1), fx.idtable.occs4(1), ...
                    'VariableNames', {'id','id2','z1','z2','trls1','trls2', ...
                    'trls3','trls4','occs1','occs2','occs3','occs4'});
                zero = struct();
                for q = 1:4
                    zero.(sprintf('mu_ss%d',q)) = fx.b(:,q);
                    zero.(sprintf('er_ss%d',q)) = fx.b_sigma(:,q);
                end
                ss = psyrat_ssrel_dodiffdynrel_trt_gamma_ls('b', fx.b, ...
                    'b_sigma', fx.b_sigma, 'b_dim', fx.b_dim, ...
                    'b_sigma_dim', fx.b_sigma_dim, 'sd_id', fx.sd_id, ...
                    'cor_id', fx.cor_id, 'sd_trl', fx.sd_trl, ...
                    'cor_trl', fx.cor_trl, 'sd_occ', fx.sd_occ, ...
                    'cor_occ', fx.cor_occ, 'sd_tid', fx.sd_tid, ...
                    'cor_tid', fx.cor_tid, 'sd_oid', fx.sd_oid, ...
                    'cor_oid', fx.cor_oid, 'sd_to', fx.sd_to, ...
                    'cor_to', fx.cor_to, ...
                    'mu_ss1', zero.mu_ss1, 'mu_ss2', zero.mu_ss2, ...
                    'mu_ss3', zero.mu_ss3, 'mu_ss4', zero.mu_ss4, ...
                    'er_ss1', zero.er_ss1, 'er_ss2', zero.er_ss2, ...
                    'er_ss3', zero.er_ss3, 'er_ss4', zero.er_ss4, ...
                    'idtable', fake, 'nprime', fx.nprime, ...
                    'nprime_o', fx.nprime_o, 'reltype', reltype, ...
                    'ndim', 2, 'CI', 0.95, 'cvec', fx.cvec);

                testCase.verifyEqual(ss.gen_pt_std, surf.G.pt, 'AbsTol', 1e-12, ...
                    sprintf('Surface G (reltype %d) must match.', reltype));
                testCase.verifyEqual(ss.dep_pt_std, surf.D.pt, 'AbsTol', 1e-12, ...
                    sprintf('Surface D (reltype %d) must match.', reltype));
                testCase.verifyEqual(ss.iccg_pt, surf.ICCg.pt, 'AbsTol', 1e-12, ...
                    sprintf('Surface relative ICC (reltype %d) must match.', reltype));
            end
        end

        function testZeroCountCellIsFlooredNotFatal(testCase)
            %B16, analysis-29 arm: per-cell TRIAL counts and per-cell
            %OCCASION counts both reach the harmonic composites, so a zero
            %in either vector tripped psyrat_dod_composite's n > 0 guard and
            %killed the whole summary. Floored at 1 (mirroring
            %psyrat_ssrel_diff); the output table keeps the true zeros.
            fx = localLoadTrtLane();
            fz = fx;
            fz.idtable.trls3(2) = 0;
            fz.idtable.occs1(1) = 0;
            ssrel = localRunTrt(fz, 3);
            testCase.verifyTrue(all(isfinite(ssrel.gen_pt)) && ...
                all(isfinite(ssrel.dep_pt)), ...
                'zero trial/occasion cells must yield finite floored coefficients');

            %floor semantics: exactly the arithmetic of explicit ones
            f1 = fx;
            f1.idtable.trls3(2) = 1;
            f1.idtable.occs1(1) = 1;
            want = localRunTrt(f1, 3);
            testCase.verifyEqual(ssrel.gen_pt, want.gen_pt, ...
                'the floored zeros must compute exactly as n = 1');
            testCase.verifyEqual(ssrel.dep_pt, want.dep_pt, ...
                'the floored zeros must compute exactly as n = 1');

            %the floor applies to the composite inputs only
            testCase.verifyEqual(ssrel.trls3(2), 0, ...
                'the output table must keep the true zero trial count');
            testCase.verifyEqual(ssrel.occs1(1), 0, ...
                'the output table must keep the true zero occasion count');
        end

        function testFlooredRowIsDisclosedByTwoNFlooredColumns(testCase)
            %B16 residual (2026-08-18), analysis-29 arm. This design floors
            %BOTH facets, so one combined counter would tell the reader that
            %something was substituted without saying whether it was a trial
            %cell or an occasion cell -- different problems with different
            %fixes. Two columns keep them distinguishable.
            fx = localLoadTrtLane();

            %negative control first
            clean = localRunTrt(fx, 3);
            vn = clean.Properties.VariableNames;
            testCase.verifyTrue(ismember('n_floored', vn) && ...
                ismember('n_floored_occ', vn), ...
                'both disclosure columns must exist on every row');
            testCase.verifyEqual(unique(clean.n_floored), 0, ...
                'no empty trial cell, so nothing may be reported as floored');
            testCase.verifyEqual(unique(clean.n_floored_occ), 0, ...
                'no empty occasion cell, so nothing may be reported as floored');

            %one empty TRIAL cell on participant 2, one empty OCCASION cell on
            %participant 1 -- the differential: each must land in its own
            %column and on its own row, or a single shared counter would pass
            fz = fx;
            fz.idtable.trls3(2) = 0;
            fz.idtable.occs1(1) = 0;
            ssrel = localRunTrt(fz, 3);

            testCase.verifyEqual(ssrel.n_floored(2), 1, ...
                'the empty trial cell must be disclosed on participant 2');
            testCase.verifyEqual(ssrel.n_floored_occ(2), 0, ...
                'participant 2 had no empty occasion cell');
            testCase.verifyEqual(ssrel.n_floored_occ(1), 1, ...
                'the empty occasion cell must be disclosed on participant 1');
            testCase.verifyEqual(ssrel.n_floored(1), 0, ...
                'participant 1 had no empty trial cell');

            %counts stay true
            testCase.verifyEqual(ssrel.trls3(2), 0, ...
                'the count column must still carry the true zero');
            testCase.verifyEqual(ssrel.occs1(1), 0, ...
                'the occasion column must still carry the true zero');
        end

        function testFloorCoversExactlyTheEmptyCellZero(testCase)
            %B16 hardening (refuter-driven, 2026-08-17), analysis-29 arm:
            %the floor's envelope must be EXACTLY the real zero on BOTH
            %count vectors. max(1,NaN) is silently 1 (MATLAB's max omits
            %NaN) and negatives floor the same way, so a max(1,...) floor
            %hides corrupt counts from psyrat_dod_composite's n > 0 guard.
            %See the one-facet twin for the full rationale.
            fx = localLoadTrtLane();

            fnan = fx;
            fnan.idtable.occs2(1) = NaN;
            testCase.verifyError(@() localRunTrt(fnan, 3), 'varargin:n', ...
                'a NaN occasion count must still raise the composite guard');

            fneg = fx;
            fneg.idtable.trls1(2) = -2;
            testCase.verifyError(@() localRunTrt(fneg, 3), 'varargin:n', ...
                'a negative trial count must still raise the composite guard');
        end

    end
end

% ---------------------------------------------------------------------------

function fx = localLoadTrtLane()
here = fileparts(mfilename('fullpath'));
root = fullfile(here, 'baselines', 'nonconcurrent_dod');
rd = @(n) readtable(fullfile(root, sprintf('nonconcurrent_dod_twofacet_%s.csv', ...
    n)), 'VariableNamingRule', 'preserve');

meta = rd('meta');
draws = rd('draws');
fx.persontable = rd('persontable');
expected = rd('stage3_expected');
fx.actual = sortrows(expected(strcmp(expected.design,'actual'),:), ...
    {'posterior_row','p_index'});
fx.standardized = sortrows(expected(strcmp(expected.design,'standardized'),:), ...
    {'posterior_row','p_index'});

fx.nd = meta.n_draws;
fx.K = meta.K;
fx.P = meta.P;
fx.labels = {char(meta.outcome_1), char(meta.outcome_2), ...
    char(meta.outcome_3), char(meta.outcome_4)};
fx.cvec = [meta.w1 meta.w2 meta.w3 meta.w4];
fx.mu_offset = [meta.mu_offset_1 meta.mu_offset_2 meta.mu_offset_3 meta.mu_offset_4];
fx.log_sigma_offset = [meta.log_sigma_offset_1 meta.log_sigma_offset_2 ...
    meta.log_sigma_offset_3 meta.log_sigma_offset_4];
fx.nprime = [meta.n_i_error_std meta.n_i_correct_std ...
    meta.n_i_error_std meta.n_i_correct_std];
fx.nprime_o = [meta.n_o_error_std meta.n_o_correct_std ...
    meta.n_o_error_std meta.n_o_correct_std];

fx.b = fx.mu_offset + localParam(draws, 'alpha_mu');
fx.b_sigma = fx.log_sigma_offset + localParam(draws, 'alpha_sigma');
fx.b_dim = localParam3(draws, 'beta_mu', fx.K);
fx.b_sigma_dim = localParam3(draws, 'beta_sigma', fx.K);
fx.sd_id = localParamWide(draws, 'sd_p_joint', 8);
fx.cor_id = localParamMat(draws, 'R_p_joint', 8);
%reference facet names -> toolbox blocks: i=trl, o=occ, pi=tid, po=oid, io=to
fx.sd_trl = localParamWide(draws, 'sd_i', 4);
fx.cor_trl = localParamMat(draws, 'R_i', 4);
fx.sd_occ = localParamWide(draws, 'sd_o', 4);
fx.cor_occ = localParamMat(draws, 'R_o', 4);
fx.sd_tid = localParamWide(draws, 'sd_pi', 4);
fx.cor_tid = localParamMat(draws, 'R_pi', 4);
fx.sd_oid = localParamWide(draws, 'sd_po', 4);
fx.cor_oid = localParamMat(draws, 'R_po', 4);
fx.sd_to = localParamWide(draws, 'sd_io', 4);
fx.cor_to = localParamMat(draws, 'R_io', 4);
u_p = localParamMat(draws, 'u_p_joint', [fx.P 8]);
for q = 1:4
    fx.(sprintf('mu_ss%d',q)) = fx.b(:,q) + u_p(:,:,q);
    fx.(sprintf('er_ss%d',q)) = fx.b_sigma(:,q) + u_p(:,:,4+q);
end

pt = fx.persontable;
fx.idtable = table(pt.ID, pt.p_index, pt.Dim1_z, pt.Dim2_z, ...
    pt.n_i_error_actual, pt.n_i_correct_actual, ...
    pt.n_i_error_actual, pt.n_i_correct_actual, ...
    pt.n_o_error_actual, pt.n_o_correct_actual, ...
    pt.n_o_error_actual, pt.n_o_correct_actual, ...
    'VariableNames', {'id','id2','z1','z2','trls1','trls2','trls3','trls4', ...
    'occs1','occs2','occs3','occs4'});
fx.ndim = 2;
end

function REL = localSyntheticTrtRel(fx)
%A synthetic analysis-29 REL in the engine's storage idiom.
REL = struct();
REL.analysis = 'ic_dodiff_dynrel_sserrvar_trt';
REL.family = 'gamma';
REL.gammascale = 2;
REL.ndim = 2;
REL.dim_names = {'dim1','dim2'};
out = struct();
out.labels = {'none'};
out.b = {num2cell(fx.b)};
out.b_sigma = {num2cell(fx.b_sigma)};
out.b_dim = {num2cell(fx.b_dim)};
out.b_sigma_dim = {num2cell(fx.b_sigma_dim)};
out.sd_id = {num2cell(fx.sd_id)};
out.sd_trl = {num2cell(fx.sd_trl)};
out.sd_occ = {num2cell(fx.sd_occ)};
out.sd_tid = {num2cell(fx.sd_tid)};
out.sd_oid = {num2cell(fx.sd_oid)};
out.sd_to = {num2cell(fx.sd_to)};
out.cor_id = {num2cell(fx.cor_id)};
out.cor_trl = {num2cell(fx.cor_trl)};
out.cor_occ = {num2cell(fx.cor_occ)};
out.cor_tid = {num2cell(fx.cor_tid)};
out.cor_oid = {num2cell(fx.cor_oid)};
out.cor_to = {num2cell(fx.cor_to)};
for q = 1:4
    out.(sprintf('mu_ss%d',q)) = {num2cell(fx.(sprintf('mu_ss%d',q)))};
    out.(sprintf('er_var_ss%d',q)) = {num2cell(fx.(sprintf('er_ss%d',q)))};
end
out.ssinfo = {fx.idtable};
out.dimz = {[fx.idtable.z1, fx.idtable.z2]};
out.ntrials = [5;7;5;7];
out.nocc = [2;2;2;2];
out.elabels = fx.labels(:);
out.glabels = {'none'};
out.dod_labels = struct(...
    'inference_framework',...
    'modular_four_gamma_margins_all_residual_covariances_fixed_zero_v1',...
    'cut_posterior_status',...
    'no_dependence_module_exists_because_all_constituents_are_assumed_nonconcurrent',...
    'residual_dependence_structure',...
    'all six pairwise conditional residual covariances fixed exactly to zero',...
    'scale_parameterization','person_specific_log_absolute_residual_sd',...
    'estimand_version','signed_four_component_nonconcurrent_dod_v1');
REL.out = out;
end

function [ssrel, dl] = localRunTrt(fx, reltype)
if nargout > 1
    [ssrel, dl] = local_call(fx, reltype);
else
    ssrel = local_call(fx, reltype);
end
end

function varargout = local_call(fx, reltype)
[varargout{1:nargout}] = psyrat_ssrel_dodiffdynrel_trt_gamma_ls('b', fx.b, ...
    'b_sigma', fx.b_sigma, 'b_dim', fx.b_dim, ...
    'b_sigma_dim', fx.b_sigma_dim, 'sd_id', fx.sd_id, 'cor_id', fx.cor_id, ...
    'sd_trl', fx.sd_trl, 'cor_trl', fx.cor_trl, ...
    'sd_occ', fx.sd_occ, 'cor_occ', fx.cor_occ, ...
    'sd_tid', fx.sd_tid, 'cor_tid', fx.cor_tid, ...
    'sd_oid', fx.sd_oid, 'cor_oid', fx.cor_oid, ...
    'sd_to', fx.sd_to, 'cor_to', fx.cor_to, ...
    'mu_ss1', fx.mu_ss1, 'mu_ss2', fx.mu_ss2, ...
    'mu_ss3', fx.mu_ss3, 'mu_ss4', fx.mu_ss4, ...
    'er_ss1', fx.er_ss1, 'er_ss2', fx.er_ss2, ...
    'er_ss3', fx.er_ss3, 'er_ss4', fx.er_ss4, ...
    'idtable', fx.idtable, 'nprime', fx.nprime, 'nprime_o', fx.nprime_o, ...
    'reltype', reltype, 'ndim', fx.ndim, 'CI', 0.95, 'cvec', fx.cvec);
end

function v = localFlatCell(arr3, q)
x = squeeze(arr3(:, q, :));
xt = x';
v = xt(:);
end

function v = localParam(draws, nm)
sub = draws(strcmp(draws.param, nm), :);
nd = max(sub.draw);
v = zeros(nd, 4);
for i = 1:height(sub)
    v(sub.draw(i), sub.idx1(i)) = sub.value(i);
end
end

function v = localParamWide(draws, nm, k)
sub = draws(strcmp(draws.param, nm), :);
nd = max(sub.draw);
v = zeros(nd, k);
for i = 1:height(sub)
    v(sub.draw(i), sub.idx1(i)) = sub.value(i);
end
end

function v = localParam3(draws, nm, K)
sub = draws(strcmp(draws.param, nm), :);
nd = max(sub.draw);
v = zeros(nd, 4, K);
for i = 1:height(sub)
    v(sub.draw(i), sub.idx1(i), sub.idx2(i)) = sub.value(i);
end
end

function v = localParamMat(draws, nm, k)
if isscalar(k), k = [k k]; end
sub = draws(strcmp(draws.param, nm), :);
nd = max(sub.draw);
v = zeros(nd, k(1), k(2));
for i = 1:height(sub)
    v(sub.draw(i), sub.idx1(i), sub.idx2(i)) = sub.value(i);
end
end
