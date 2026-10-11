classdef TestGammaDynrelSurfaceTarget < PsyRATTestBase
    %TESTGAMMADYNRELSURFACETARGET Q9/GUI-1 (owner rulings 9-12, 2026-09-13):
    % the explicit gamma SURFACE target for the non-difference dynamic designs
    % (11/26/14/27, both scale parameterizations, one and two dimensions), the
    % Q8 legacy routing of unstamped analysis-27 results, and every consumer
    % that must agree about the resolved target (summary metadata, component
    % table, plot caption/UserData, headless report, viewer control/cache).
    %
    % INDEPENDENT REFERENCES. The frozen 2C Q8/Q9 expectation set (local
    % development evidence, not shipped; 70-digit Decimal; support-union products;
    % produced without running production code), transcribed to double.
    %
    %   two-facet LS  (s_p .6, s_o = s_i .5, sigma .2, n_i 3, n_o 2, CES):
    %     typical G .79974834803861694  D .547157189479854
    %     primary (s_w .45) G .799665194144230765  D .5471182657569164
    %   two-facet log-nu (nu 18):
    %     typical G .76965090115835433  D .5328997943609791
    %     primary (s_v .5, cov .2) G .75366888065707248  D .5251886610478044
    %     primary (s_v 0)          G .740124734017914 (NOT the typical value)
    %     equality control: s_p .6, s_v 1.2, cov .72 -> Var(2U - V) = 0
    %   one-facet LS  (n_i 3): typical G .9129552957512670 D .7611352468937906
    %                          primary (s_w .45) G .9126770738475140 D .7609418551313894
    %   one-facet log-nu:      typical G .8537938983160641 D .7195663459702325
    %                          primary (s_v .5, cov .2) G .8236158308287681 D .6980113945461412
    %
    % Nothing here was executed by the author; root runs the class.

    methods (Static, Access = private)

        function d = rep(v, n)
            d = repmat(v, n, 1);
        end

        function args = localTrtArgs(param, nd, z1, extra)
            %Two-facet main-fixture calculator arguments at zero slopes, with
            %the observation-variance heterogeneity given in extra
            %(sig_v/rho_pv for log-nu, sig_sd for LS).
            rp = @(x) TestGammaDynrelSurfaceTarget.rep(x, nd);
            base = {'alpha0', rp(log(5)), 'b', rp(0), ...
                'sig_p', rp(.6), 'sig_occ', rp(.5), 'sig_trl', rp(.5), ...
                'sig_occxid', rp(0), 'sig_trlxid', rp(0), 'sig_trlxocc', rp(0), ...
                'ndim', 1, 'z1', z1, 'obs', 3, 'nocc', 2, 'reltype', 3, 'CI', .95};
            if strcmp(param,'ls')
                args = [base, {'logsig0', rp(log(.2)), 'b_sigma', rp(0), ...
                    'sig_sd', rp(extra.sig_sd)}];
            else
                args = [base, {'lognu0', rp(log(18)), 'b_nu', rp(0), ...
                    'sig_v', rp(extra.sig_v), 'rho_pv', rp(extra.rho_pv)}];
            end
        end

        function rel = localTrtRel(param, stamped, stored_target)
            %Synthetic analysis-27 result: two strata with distinct z, counts
            %and person effects, subject rows in a PERMUTED id2 order.
            %stored_target (optional) is the stamp's saved surface target.
            if nargin < 3; stored_target = 'population_marginal'; end
            nd = 4;
            rp = @(x) repmat(x, nd, 1);
            rel = struct();
            rel.analysis = 'ic_dynrel_sserrvar_trt';
            rel.family = 'gamma';
            rel.ndim = 1; rel.dim_names = {'dim1'};
            rel.filename = 'gamma_trt_target_fixture.psyrat';
            o = struct();
            o.labels = {'A','B'};
            zA = [.5; -.8; 1.1]; zB = [-.3; .9];
            o.dimz = {zA, zB};
            o.ntrials = [8 10];
            o.nocc = [2 2];
            o.mu = [rp(log(5)) rp(log(4))];
            o.sig_occ = [rp(.5) rp(.4)];
            o.sig_trl = [rp(.5) rp(.45)];
            o.sig_trlxid = [rp(0) rp(.2)];
            o.sig_occxid = [rp(0) rp(.3)];
            o.sig_trlxocc = [rp(0) rp(.25)];
            o.b = {num2cell(rp(.2)), num2cell(rp(-.1))};
            o.ssinfo = { ...
                table((101:103)', [3;1;2], zA, [9;4;6], [2;2;2], ...
                'VariableNames', {'id','id2','z1','trls','occs'}), ...
                table((201:202)', [2;1], zB, [12;5], [2;2], ...
                'VariableNames', {'id','id2','z1','trls','occs'})};
            if strcmp(param,'ls')
                rel.gammascale = 2;
                rho = .3; L = [1 0; rho sqrt(1-rho^2)];
                o.pop_sdlog = [rp(log(.2)) rp(log(.3))];
                o.gro_sds = {num2cell(rp([.6 .45])), num2cell(rp([.5 .3]))};
                o.b_sigma = {num2cell(rp(.15)), num2cell(rp(-.05))};
                o.ind_sdlog = {num2cell(rp([.3 -.15 .05])), num2cell(rp([-.2 .1]))};
            else
                rel.gammascale = 1;
                rho = 2/3; L = [1 0; rho sqrt(1-rho^2)];
                o.pop_lognu = [rp(log(18)) rp(log(15))];
                o.gro_sds = {num2cell(rp([.6 .5])), num2cell(rp([.5 .4]))};
                o.b_nu = {num2cell(rp(-.25)), num2cell(rp(.1))};
                o.ind_bs = {num2cell(rp([.25 -.10 .05])), num2cell(rp([-.15 .2]))};
                o.ind_nu = {num2cell(rp([-.20 .30 0])), num2cell(rp([.1 -.1]))};
            end
            ch = zeros(1,2,2); ch(1,:,:) = L;
            o.chol_corrmat = {num2cell(repmat(ch, nd, 1)), num2cell(repmat(ch, nd, 1))};
            if stamped
                if strcmp(param,'ls'); gp = 'ls'; else; gp = 'lognu'; end
                o.gamma_target = struct('gamma_target_version', 2, ...
                    'surface_target', stored_target, ...
                    'reference_target', 'not_applicable', ...
                    'subject_target', 'own_observation_variance_with_population_mean_components_at_own_z', ...
                    'remainder_scheme', 'population_highest_order_plus_person_observation_v2', ...
                    'parameterization', gp, ...
                    'count_source', 'median_trials_per_occupied_occasion', ...
                    'surface_count_source', 'default_stratum_median', ...
                    'analysis', 'ic_dynrel_sserrvar_trt');
            end
            rel.out = o;
        end

        function rel = localOneFacetRel(id2)
            %Synthetic analysis-26 log-nu result, one stratum, id2 as given.
            nd = 4;
            rp = @(x) repmat(x, nd, 1);
            n = numel(id2);
            z = linspace(-1, 1, n)';
            rel = struct();
            rel.analysis = 'ic_dynrel_sserrvar';
            rel.family = 'gamma';
            rel.ndim = 1; rel.dim_names = {'dim1'};
            rel.filename = 'gamma_onefacet_target_fixture.psyrat';
            rho = 2/3; L = [1 0; rho sqrt(1-rho^2)];
            ch = zeros(1,2,2); ch(1,:,:) = L;
            o = struct();
            o.labels = {'none'};
            o.dimz = {z};
            o.ntrials = 3;
            o.mu = rp(log(5));
            o.pop_lognu = rp(log(18));
            o.gro_sds = {num2cell(rp([.6 .5]))};
            o.sig_trl = rp(.5);
            o.b = {num2cell(rp(.2))};
            o.b_nu = {num2cell(rp(-.25))};
            o.chol_corrmat = {num2cell(repmat(ch, nd, 1))};
            o.ind_bs = {num2cell(rp(linspace(-.2, .2, n)))};
            o.ind_nu = {num2cell(rp(linspace(.1, -.1, n)))};
            o.ssinfo = {table((1:n)' + 300, id2(:), z, repmat(3, n, 1), ...
                'VariableNames', {'id','id2','z1','trls'})};
            rel.out = o;
        end

        function pd = localData(rel)
            pd = struct('rel', rel, 'ver', 'test');
            pd.proc = struct('savename', rel.filename);
        end

        function localVerifyLayout(testCase, fig, label)
            %Every uicontrol inside the figure, and the whole window inside a
            %1728x1117 display (the same invariant TestTrtSserrvarViewerRoute
            %pins for the test-retest viewer).
            u = findall(fig, 'Type', 'uicontrol');
            testCase.assertNotEmpty(u, sprintf('%s: viewer holds no uicontrols', label));
            ys = arrayfun(@(h) h.Position(2), u);
            hs = arrayfun(@(h) h.Position(4), u);
            testCase.verifyGreaterThanOrEqual(min(ys), 0, sprintf( ...
                '%s: a control starts below the figure bottom (y = %.2f)', label, min(ys)));
            testCase.verifyLessThanOrEqual(max(ys + hs), fig.Position(4), sprintf( ...
                '%s: a control ends above the figure top (y+h = %.2f, height %.2f)', ...
                label, max(ys + hs), fig.Position(4)));
            testCase.verifyLessThanOrEqual(fig.Position(2) + fig.Position(4), 1117, sprintf( ...
                '%s: the window must fit a 1728x1117 display (origin %.0f + height %.0f)', ...
                label, fig.Position(2), fig.Position(4)));
        end

        function fig = localViewerFig(testCase, pd, prefs)
            %prefs (optional) sets e.g. guis.fsize; default psyrat_defaults()
            if nargin < 3; prefs = psyrat_defaults(); end
            name = 'Dynamic Reliability: View Results';
            close(findall(0, 'Type', 'figure', 'Name', name));
            psyrat_startview_dynrel('psyrat_prefs', prefs, 'psyrat_data', pd);
            fig = findall(0, 'Type', 'figure', 'Name', name);
            testCase.assertNotEmpty(fig, 'Dynrel viewer figure not created');
            fig = fig(1);
        end

        function out = localFireButton(testCase, fig, label)
            %Fire the push button whose label contains LABEL and return the
            %newest 'psyrat_output' figure it produced.
            before = findall(0, 'Type', 'figure', 'Tag', 'psyrat_output');
            btns = findall(fig, 'Style', 'pushbutton');
            fired = false;
            for b = btns(:)'
                s = b.String; if iscell(s), s = strjoin(s, ' '); end
                if contains(string(s), label)
                    cb = b.Callback;
                    cb(b, []);
                    fired = true;
                    break;
                end
            end
            testCase.assertTrue(fired, sprintf('button "%s" not found', label));
            after = findall(0, 'Type', 'figure', 'Tag', 'psyrat_output');
            out = setdiff(after, before);
            testCase.assertEqual(numel(out), 1, 'expected exactly one new output figure');
        end

    end

    methods (Test)

        % ---------- validation ----------

        function testInvalidTargetIsRejectedNeverCoerced(testCase)
            nd = 4;
            argsL = TestGammaDynrelSurfaceTarget.localTrtArgs('ls', nd, [0 .5], struct('sig_sd', .45));
            argsN = TestGammaDynrelSurfaceTarget.localTrtArgs('lognu', nd, [0 .5], struct('sig_v', .5, 'rho_pv', 2/3));
            testCase.verifyError(@() psyrat_rel_dynrel_trt_gamma_ls(argsL{:}, 'surface_target', 'typical'), ...
                'dynrel:badsurfacetarget');
            testCase.verifyError(@() psyrat_rel_dynrel_trt_gamma(argsN{:}, 'surface_target', 'marginal'), ...
                'dynrel:badsurfacetarget');
            testCase.verifyError(@() psyrat_gamma_surface_target(1), 'dynrel:badsurfacetarget');
            testCase.verifyEqual(psyrat_gamma_surface_target([]), 'population_marginal');
            testCase.verifyEqual(psyrat_gamma_surface_target("typical_person"), 'typical_person');
            % summary: bad value, and any value on a Gaussian run
            rel = TestGammaDynrelSurfaceTarget.localTrtRel('lognu', true);
            testCase.verifyError(@() psyrat_dynrel_summary(rel, 'surface_target', 'population'), ...
                'dynrel:badsurfacetarget');
            gauss = PsyRATTestDataFactory.makeDynrelSubjNondiffSummaryData();
            testCase.verifyError(@() psyrat_dynrel_summary(gauss.rel, 'surface_target', 'typical_person'), ...
                'dynrel:surfacetargetunsupported');
            testCase.verifyError(@() psyrat_dynrel_summary(gauss.rel, 'remainder_scheme', 'legacy_observation_only_v1'), ...
                'dynrel:remainderschemeunsupported');
            % subject calculators: bad scheme
            testCase.verifyError(@() psyrat_ssrel_dynrel_trt_gamma_ls('alpha0', 1, 'b', 0, ...
                'logsig0', 0, 'b_sigma', 0, 'sig_p', .5, 'sig_occ', .1, 'sig_trl', .1, ...
                'sig_occxid', 0, 'sig_trlxid', 0, 'sig_trlxocc', 0, 'ind_sd', 0, ...
                'idtable', table(1, 1, 0, 3, 2, 'VariableNames', {'id','id2','z1','trls','occs'}), ...
                'ndim', 1, 'reltype', 3, 'nocc', 2, 'CI', .95, 'remainder_scheme', 'v3'), ...
                'dynrel:badremainderscheme');
        end

        % ---------- two-facet calculators against frozen anchors ----------

        function testTwoFacetLocationScaleTargets(testCase)
            nd = 4; z1 = [0 .7];
            args = TestGammaDynrelSurfaceTarget.localTrtArgs('ls', nd, z1, struct('sig_sd', .45));
            def = psyrat_rel_dynrel_trt_gamma_ls(args{:});
            pop = psyrat_rel_dynrel_trt_gamma_ls(args{:}, 'surface_target', 'population_marginal');
            typ = psyrat_rel_dynrel_trt_gamma_ls(args{:}, 'surface_target', 'typical_person');
            testCase.verifyEqual(def.surface_target, 'population_marginal');
            testCase.verifyEqual(def.G, pop.G, 'the default must be the population-marginal primary');
            testCase.verifyEqual(def.D, pop.D);
            % zero slopes: every z equals the frozen anchor
            testCase.verifyEqual(pop.G.pt, [.799665194144230765 .799665194144230765], 'RelTol', 1e-12);
            testCase.verifyEqual(pop.D.pt, [.5471182657569164 .5471182657569164], 'RelTol', 1e-12);
            testCase.verifyEqual(typ.G.pt, [.79974834803861694 .79974834803861694], 'RelTol', 1e-12);
            testCase.verifyEqual(typ.D.pt, [.547157189479854 .547157189479854], 'RelTol', 1e-12);
            testCase.verifyEqual(typ.sige_pt, sqrt([2.1052177285562118 2.1052177285562118]), 'RelTol', 1e-12);
            testCase.verifyEqual(typ.surface_target, 'typical_person');
            % s_sd = 0: the two targets coincide (to the subtraction form's rounding)
            args0 = TestGammaDynrelSurfaceTarget.localTrtArgs('ls', nd, z1, struct('sig_sd', 0));
            pop0 = psyrat_rel_dynrel_trt_gamma_ls(args0{:});
            typ0 = psyrat_rel_dynrel_trt_gamma_ls(args0{:}, 'surface_target', 'typical_person');
            testCase.verifyEqual(pop0.G.pt, typ0.G.pt, 'RelTol', 1e-12);
            testCase.verifyEqual(pop0.G.pt(1), .79974834803861694, 'RelTol', 1e-12);
        end

        function testTwoFacetLogNuTargetsAndEqualityControls(testCase)
            nd = 4; z1 = [0 -.4];
            args = TestGammaDynrelSurfaceTarget.localTrtArgs('lognu', nd, z1, struct('sig_v', .5, 'rho_pv', 2/3));
            pop = psyrat_rel_dynrel_trt_gamma(args{:});
            typ = psyrat_rel_dynrel_trt_gamma(args{:}, 'surface_target', 'typical_person');
            testCase.verifyEqual(pop.G.pt(1), .75366888065707248, 'RelTol', 1e-12, 'log-nu primary G');
            testCase.verifyEqual(pop.D.pt(1), .5251886610478044, 'RelTol', 1e-12, 'log-nu primary D');
            testCase.verifyEqual(typ.G.pt, [.76965090115835433 .76965090115835433], 'RelTol', 1e-12, 'log-nu typical G');
            testCase.verifyEqual(typ.D.pt(1), .5328997943609791, 'RelTol', 1e-12, 'log-nu typical D');
            testCase.verifyEqual(typ.sige_pt(1), sqrt(9.6160005853868930), 'RelTol', 1e-12);
            % s_v = 0 is NOT the companion: the primary still integrates u_p
            args0 = TestGammaDynrelSurfaceTarget.localTrtArgs('lognu', nd, z1, struct('sig_v', 0, 'rho_pv', 0));
            pop0 = psyrat_rel_dynrel_trt_gamma(args0{:});
            typ0 = psyrat_rel_dynrel_trt_gamma(args0{:}, 'surface_target', 'typical_person');
            testCase.verifyEqual(pop0.G.pt(1), .740124734017914, 'RelTol', 1e-12, 's_v = 0 primary');
            testCase.verifyEqual(typ0.G.pt(1), .76965090115835433, 'RelTol', 1e-12);
            testCase.verifyGreaterThan(abs(pop0.G.pt(1) - typ0.G.pt(1)), 1e-3, ...
                'zero scale variance alone must not reproduce the typical-person surface');
            % nonvacuous equality: Var(2U - V) = 4 v_P + v_V - 4 cov = 0 at
            % s_p .6, s_v 1.2, cov .72 (rho = 1), with nonzero population signal
            argsE = TestGammaDynrelSurfaceTarget.localTrtArgs('lognu', nd, z1, struct('sig_v', 1.2, 'rho_pv', 1));
            popE = psyrat_rel_dynrel_trt_gamma(argsE{:});
            typE = psyrat_rel_dynrel_trt_gamma(argsE{:}, 'surface_target', 'typical_person');
            testCase.verifyEqual(popE.G.pt, typE.G.pt, 'RelTol', 1e-12, 'PSD-boundary equality control');
            testCase.verifyEqual(popE.D.pt, typE.D.pt, 'RelTol', 1e-12);
            testCase.verifyGreaterThan(popE.G.pt(1), .5, 'the equality control must keep nonzero signal');
        end

        % ---------- one-facet calculators ----------

        function testOneFacetTargetsMatchFrozen(testCase)
            nd = 4; z1 = [0 .3];
            rp = @(x) TestGammaDynrelSurfaceTarget.rep(x, nd);
            argsL = {'alpha0', rp(log(5)), 'b', rp(0), 'logsig0', rp(log(.2)), 'b_sigma', rp(0), ...
                'sig_p', rp(.6), 'sig_sd', rp(.45), 'sig_i', rp(.5), ...
                'ndim', 1, 'z1', z1, 'obs', 3, 'CI', .95};
            popL = psyrat_rel_dynrel_gamma_ls(argsL{:});
            typL = psyrat_rel_dynrel_gamma_ls(argsL{:}, 'surface_target', 'typical_person');
            testCase.verifyEqual(popL.G.pt(1), .9126770738475140, 'RelTol', 1e-12, 'one-facet LS primary G');
            testCase.verifyEqual(popL.D.pt(1), .7609418551313894, 'RelTol', 1e-12);
            testCase.verifyEqual(typL.G.pt, [.9129552957512670 .9129552957512670], 'RelTol', 1e-12, 'one-facet LS typical G');
            testCase.verifyEqual(typL.D.pt(1), .7611352468937906, 'RelTol', 1e-12);
            testCase.verifyEqual(typL.sige_pt(1), sqrt(5.7028494835756152), 'RelTol', 1e-12);
            testCase.verifyEqual(popL.sige_pt(1), sqrt(5.7228215835778859), 'RelTol', 1e-10);
            testCase.verifyError(@() psyrat_rel_dynrel_gamma_ls(argsL{:}, 'surface_target', 'x'), 'dynrel:badsurfacetarget');

            argsN = {'alpha0', rp(log(5)), 'b', rp(0), 'lognu0', rp(log(18)), 'b_nu', rp(0), ...
                'sig_p', rp(.6), 'sig_v', rp(.5), 'rho_pv', rp(2/3), 'sig_i', rp(.5), ...
                'ndim', 1, 'z1', z1, 'obs', 3, 'CI', .95};
            popN = psyrat_rel_dynrel_gamma(argsN{:});
            typN = psyrat_rel_dynrel_gamma(argsN{:}, 'surface_target', 'typical_person');
            testCase.verifyEqual(popN.G.pt(1), .8236158308287681, 'RelTol', 1e-12, 'one-facet log-nu primary G');
            testCase.verifyEqual(popN.D.pt(1), .6980113945461412, 'RelTol', 1e-12);
            testCase.verifyEqual(typN.G.pt, [.8537938983160641 .8537938983160641], 'RelTol', 1e-12, 'one-facet log-nu typical G');
            testCase.verifyEqual(typN.D.pt(1), .7195663459702325, 'RelTol', 1e-12);
            testCase.verifyEqual(typN.sige_pt(1), sqrt(10.242630791075971), 'RelTol', 1e-12);
            testCase.verifyError(@() psyrat_rel_dynrel_gamma(argsN{:}, 'surface_target', 'x'), 'dynrel:badsurfacetarget');
        end

        function testTwoDimensionalSurfacesHonorTheTarget(testCase)
            nd = 4;
            rp = @(x) TestGammaDynrelSurfaceTarget.rep(x, nd);
            z1 = [-1 0 1]; z2 = [-.5 .5];
            argsN = {'alpha0', rp(log(5)), 'b', repmat([0 0 0], nd, 1), ...
                'lognu0', rp(log(18)), 'b_nu', repmat([0 0 0], nd, 1), ...
                'sig_p', rp(.6), 'sig_v', rp(.5), 'rho_pv', rp(2/3), ...
                'sig_occ', rp(.5), 'sig_trl', rp(.5), 'sig_occxid', rp(0), ...
                'sig_trlxid', rp(0), 'sig_trlxocc', rp(0), ...
                'ndim', 2, 'z1', z1, 'z2', z2, 'obs', 3, 'nocc', 2, 'reltype', 3, 'CI', .95};
            pop = psyrat_rel_dynrel_trt_gamma(argsN{:});
            typ = psyrat_rel_dynrel_trt_gamma(argsN{:}, 'surface_target', 'typical_person');
            testCase.verifySize(typ.G.pt, [3 2]);
            testCase.verifyEqual(typ.G.pt, repmat(.76965090115835433, 3, 2), 'RelTol', 1e-12);
            testCase.verifyEqual(pop.G.pt, repmat(.75366888065707248, 3, 2), 'RelTol', 1e-12);
            argsL = {'alpha0', rp(log(5)), 'b', repmat([0 0 0], nd, 1), ...
                'logsig0', rp(log(.2)), 'b_sigma', repmat([0 0 0], nd, 1), ...
                'sig_p', rp(.6), 'sig_sd', rp(.45), 'sig_i', rp(.5), ...
                'ndim', 2, 'z1', z1, 'z2', z2, 'obs', 3, 'CI', .95};
            typL = psyrat_rel_dynrel_gamma_ls(argsL{:}, 'surface_target', 'typical_person');
            testCase.verifyEqual(typL.G.pt, repmat(.9129552957512670, 3, 2), 'RelTol', 1e-12);
        end

        % ---------- summary: selector, subjects, counts, metadata ----------

        function testSummarySelectorChangesSurfaceOnly(testCase)
            for param = {'lognu', 'ls'}
                rel = TestGammaDynrelSurfaceTarget.localTrtRel(param{1}, true);
                sA  = psyrat_dynrel_summary(rel, 'CI', .95);
                sA2 = psyrat_dynrel_summary(rel, 'CI', .95, 'surface_target', 'population_marginal');
                sB  = psyrat_dynrel_summary(rel, 'CI', .95, 'surface_target', 'typical_person');
                sO  = psyrat_dynrel_summary(rel, 'CI', .95, 'surface_target', 'typical_person', 'obs', 30);
                testCase.verifyEqual(sA.marginal, 0, 'the historical Boolean stays 0 for gamma');
                testCase.verifyEqual(sA.surface_target, 'population_marginal');
                testCase.verifyEqual(sB.surface_target, 'typical_person');
                gt = sA.gamma_target;
                testCase.verifyEqual(gt.gamma_target_version, 2);
                testCase.verifyEqual(gt.remainder_scheme, 'population_highest_order_plus_person_observation_v2');
                testCase.verifyEqual(gt.subject_scheme_source, 'stamped_v2');
                testCase.verifyTrue(gt.stamped);
                testCase.verifyEqual(gt.parameterization, param{1});
                testCase.verifyEqual(gt.count_source, 'median_trials_per_occupied_occasion');
                testCase.verifyEqual(gt.surface_count_source, 'default_stratum_median');
                testCase.verifyEqual(sO.gamma_target.surface_count_source, 'user_override');
                testCase.verifyEqual(gt.analysis, 'ic_dynrel_sserrvar_trt');
                for s = 1:2
                    testCase.verifyEqual(sA.strata(s).G, sA2.strata(s).G, 'default == explicit primary');
                    testCase.verifyGreaterThan(max(abs(sA.strata(s).G.pt(:) - sB.strata(s).G.pt(:))), 1e-6, ...
                        'the companion must differ from the primary when the scale SD is nonzero');
                    testCase.verifyGreaterThan(max(abs(sB.strata(s).G.pt(:) - sO.strata(s).G.pt(:))), 1e-6, ...
                        'the surface n'' override must move the surface');
                    % subject table: fixed across the surface target AND the surface count
                    testCase.verifyTrue(isequal(sB.strata(s).ssrel_table, sA.strata(s).ssrel_table), ...
                        'the surface target must not touch the per-participant table');
                    testCase.verifyTrue(isequal(sO.strata(s).ssrel_table, sA.strata(s).ssrel_table), ...
                        'the surface n'' override must not touch the per-participant table');
                    testCase.verifyEqual(sA.strata(s).ssrel_table.trls, rel.out.ssinfo{s}.trls, ...
                        'subject counts stay the stored median per occupied occasion');
                    testCase.verifyEqual(sA.strata(s).ssrel_table.id2, rel.out.ssinfo{s}.id2, ...
                        'permuted id2 order preserved');
                    testCase.verifyEqual(sO.strata(s).obs, 30);
                    % residual row of the component table names the target
                    lab = string(sA.strata(s).varcomp{:,1});
                    testCase.verifyTrue(any(contains(lab, 'population-marginal surface')));
                    labB = string(sB.strata(s).varcomp{:,1});
                    testCase.verifyTrue(any(contains(labB, 'typical-person surface')));
                    rowA = sA.strata(s).varcomp(contains(lab, 'Residual variance'), :);
                    rowB = sB.strata(s).varcomp(contains(labB, 'Residual variance'), :);
                    testCase.verifyGreaterThan(rowA{1,3}, rowB{1,3}, ...
                        'the marginal residual exceeds the typical-person residual');
                end
                % the Gaussian Boolean is ignored under gamma with the documented warning,
                % and still resolves to the population surface
                sM = testCase.verifyWarning(@() psyrat_dynrel_summary(rel, 'CI', .95, 'marginal', 1), ...
                    'dynrel:marginalgamma');
                testCase.verifyEqual(sM.surface_target, 'population_marginal');
                testCase.verifyEqual(sM.strata(1).G, sA.strata(1).G);
            end
        end

        function testSubjectPointsSitOnTheCompanionAtEqualCounts(testCase)
            % a zero-effect participant with the surface's counts sits on the
            % typical-person companion, not on the primary (Q8 A + Q9)
            rel = TestGammaDynrelSurfaceTarget.localTrtRel('ls', true);
            % make subject 1 of stratum A (row 1, draw column id2(1) = 3)
            % zero-effect at z = 0 with the surface's own counts
            rel.out.ssinfo{1}.z1(1) = 0;
            rel.out.ssinfo{1}.trls(1) = rel.out.ntrials(1);
            col = rel.out.ssinfo{1}.id2(1);
            m = cell2mat(rel.out.ind_sdlog{1}); m(:,col) = 0; rel.out.ind_sdlog{1} = num2cell(m);
            sB = psyrat_dynrel_summary(rel, 'CI', .95, 'surface_target', 'typical_person', 'ngrid', 3);
            sA = psyrat_dynrel_summary(rel, 'CI', .95, 'ngrid', 3);
            % evaluate both surfaces at z = 0 through the calculator on the same grid point
            o = rel.out; gro = cell2mat(o.gro_sds{1});
            base = {'alpha0', o.mu(:,1), 'b', cell2mat(o.b{1}), 'logsig0', o.pop_sdlog(:,1), ...
                'b_sigma', cell2mat(o.b_sigma{1}), 'sig_p', gro(:,1), 'sig_sd', gro(:,2), ...
                'sig_occ', o.sig_occ(:,1), 'sig_trl', o.sig_trl(:,1), 'sig_trlxid', o.sig_trlxid(:,1), ...
                'sig_occxid', o.sig_occxid(:,1), 'sig_trlxocc', o.sig_trlxocc(:,1), ...
                'ndim', 1, 'z1', 0, 'obs', sA.strata(1).obs, 'nocc', sA.strata(1).nocc, 'reltype', 3, 'CI', .95};
            typ = psyrat_rel_dynrel_trt_gamma_ls(base{:}, 'surface_target', 'typical_person');
            pop = psyrat_rel_dynrel_trt_gamma_ls(base{:});
            g_subject = sB.strata(1).ssrel_table.gen_pt(1);
            testCase.verifyEqual(g_subject, typ.G.pt, 'RelTol', 1e-12, ...
                'an equal-count zero-effect participant must sit on the typical-person companion');
            testCase.verifyGreaterThan(abs(g_subject - pop.G.pt), 1e-6, ...
                'and not on the population-marginal primary (s_sd > 0)');
        end

        function testLegacyUnstampedTwoFacetKeepsObservationOnlyScheme(testCase)
            for param = {'lognu', 'ls'}
                fresh = TestGammaDynrelSurfaceTarget.localTrtRel(param{1}, true);
                old   = TestGammaDynrelSurfaceTarget.localTrtRel(param{1}, false);
                sF = psyrat_dynrel_summary(fresh, 'CI', .95);
                sL = psyrat_dynrel_summary(old, 'CI', .95);
                testCase.verifyEqual(sL.gamma_target.remainder_scheme, 'legacy_observation_only_v1');
                testCase.verifyEqual(sL.gamma_target.subject_scheme_source, 'legacy_unstamped');
                testCase.verifyFalse(sL.gamma_target.stamped);
                % the surface itself is the same raw-draw computation either way
                testCase.verifyEqual(sL.strata(1).G, sF.strata(1).G, 'surface does not depend on the stamp');
                % but the subject coefficients are the OLD (observation-only) ones,
                % strictly higher than v2 because R is missing
                testCase.verifyTrue(all(sL.strata(1).ssrel_table.gen_pt > sF.strata(1).ssrel_table.gen_pt + 1e-9), ...
                    'legacy (no R) coefficients must exceed the v2 ones');
                % and equal a direct legacy calculator call
                o = old.out; gro = cell2mat(o.gro_sds{1});
                common = {'alpha0', o.mu(:,1), 'b', cell2mat(o.b{1}), 'sig_p', gro(:,1), ...
                    'sig_occ', o.sig_occ(:,1), 'sig_trl', o.sig_trl(:,1), 'sig_occxid', o.sig_occxid(:,1), ...
                    'sig_trlxid', o.sig_trlxid(:,1), 'sig_trlxocc', o.sig_trlxocc(:,1), ...
                    'idtable', o.ssinfo{1}, 'ndim', 1, 'reltype', 3, 'nocc', 1, 'CI', .95, ...
                    'remainder_scheme', 'legacy_observation_only_v1'};
                if strcmp(param{1}, 'ls')
                    [gen, ~] = psyrat_ssrel_dynrel_trt_gamma_ls(common{:}, 'logsig0', o.pop_sdlog(:,1), ...
                        'b_sigma', cell2mat(o.b_sigma{1}), 'ind_sd', cell2mat(o.ind_sdlog{1}));
                else
                    [gen, ~] = psyrat_ssrel_dynrel_trt_gamma(common{:}, 'lognu0', o.pop_lognu(:,1), ...
                        'b_nu', cell2mat(o.b_nu{1}), 'ind_bs', cell2mat(o.ind_bs{1}), 'ind_nu', cell2mat(o.ind_nu{1}));
                end
                testCase.verifyEqual(sL.strata(1).ssrel_table.gen_pt, gen.dep_pt, 'AbsTol', 1e-12);
                % explicit v2 on unstamped draws is refused (migration's job);
                % explicit legacy on stamped draws is allowed and labeled
                testCase.verifyError(@() psyrat_dynrel_summary(old, 'CI', .95, ...
                    'remainder_scheme', 'population_highest_order_plus_person_observation_v2'), ...
                    'dynrel:remainderschemeunstamped');
                sX = psyrat_dynrel_summary(fresh, 'CI', .95, 'remainder_scheme', 'legacy_observation_only_v1');
                testCase.verifyEqual(sX.gamma_target.subject_scheme_source, 'explicit_override');
                testCase.verifyEqual(sX.strata(1).ssrel_table.gen_pt, sL.strata(1).ssrel_table.gen_pt, 'AbsTol', 1e-12);
                testCase.verifyError(@() psyrat_dynrel_summary(fresh, 'CI', .95, 'remainder_scheme', 'v9'), ...
                    'dynrel:badremainderscheme');
                % the surface target still switches on legacy draws
                sLt = psyrat_dynrel_summary(old, 'CI', .95, 'surface_target', 'typical_person');
                testCase.verifyTrue(isequal(sLt.strata(1).ssrel_table, sL.strata(1).ssrel_table));
                testCase.verifyGreaterThan(max(abs(sLt.strata(1).G.pt(:) - sL.strata(1).G.pt(:))), 1e-6);
            end
        end

        function testOneFacetAcceptsPermutedIdOrderButNotInvalidIds(testCase)
            rel = TestGammaDynrelSurfaceTarget.localOneFacetRel([3 1 2]);
            s = psyrat_dynrel_summary(rel, 'CI', .95);
            testCase.verifyEqual(s.gamma_target.remainder_scheme, 'population_mean_remainder_plus_person_observation_v1');
            testCase.verifyEqual(s.gamma_target.subject_scheme_source, 'unstamped_s20_convention');
            testCase.verifyEqual(s.strata(1).ssrel_table.id2, [3;1;2]);
            o = rel.out; gro = cell2mat(o.gro_sds{1});
            [gen, ~] = psyrat_ssrel_dynrel_gamma('alpha0', o.mu, 'b', cell2mat(o.b{1}), ...
                'lognu0', o.pop_lognu, 'b_nu', cell2mat(o.b_nu{1}), 'sig_p', gro(:,1), ...
                'sig_i', o.sig_trl, 'ind_bs', cell2mat(o.ind_bs{1}), 'ind_nu', cell2mat(o.ind_nu{1}), ...
                'idtable', o.ssinfo{1}, 'ndim', 1, 'CI', .95);
            testCase.verifyEqual(s.strata(1).ssrel_table.gen_pt, gen.dep_pt, 'AbsTol', 1e-12);
            % the selector works on 26 too and leaves the subjects alone
            sT = psyrat_dynrel_summary(rel, 'CI', .95, 'surface_target', 'typical_person');
            testCase.verifyTrue(isequal(sT.strata(1).ssrel_table, s.strata(1).ssrel_table));
            testCase.verifyGreaterThan(max(abs(sT.strata(1).G.pt(:) - s.strata(1).G.pt(:))), 1e-6);
            % invalid id2 (a missing column) is still refused
            bad = TestGammaDynrelSurfaceTarget.localOneFacetRel([1 2 4]);
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), ...
                'psyrat_dynrel_summary:ssinfoindex');
        end

        % ---------- report, plot, viewer ----------

        function testHeadlessReportForwardsTargetAndCaptions(testCase)
            rel = TestGammaDynrelSurfaceTarget.localTrtRel('lognu', true);
            pd = TestGammaDynrelSurfaceTarget.localData(rel);
            rA = psyrat_report(pd);
            rB = psyrat_report(pd, 'surface_target', 'typical_person');
            testCase.verifyEqual(rA.gamma_target.surface_target, 'population_marginal');
            testCase.verifyEqual(rB.gamma_target.surface_target, 'typical_person');
            testCase.verifyTrue(contains(rA.estimand_note, 'POPULATION-MARGINAL'));
            testCase.verifyTrue(contains(rB.estimand_note, 'TYPICAL PERSON'));
            testCase.verifyTrue(contains(rB.estimand_note, 'surface_target = typical_person'));
            testCase.verifyTrue(contains(rB.estimand_note, 'own log-mean and log-nu'));
            testCase.verifyTrue(contains(rB.estimand_note, 'scheme v2'));
            testCase.verifyFalse(contains(rB.estimand_note, 'delta_p'));
            testCase.verifyEqual(rB.tables.ssrel.gen_pt, rA.tables.ssrel.gen_pt, ...
                'the report''s per-participant table must not move with the target');
            testCase.verifyNotEqual(rB.tables.surface.Generalizability, rA.tables.surface.Generalizability);
            % legacy result: the note and the metadata name the LEGACY scheme
            % (an explanatory comparison against scheme v2 is allowed; what is
            % forbidden is claiming the v2 residual for the old file)
            old = TestGammaDynrelSurfaceTarget.localTrtRel('ls', false);
            rL = psyrat_report(TestGammaDynrelSurfaceTarget.localData(old));
            testCase.verifyTrue(contains(rL.estimand_note, 'LEGACY observation-only residual'));
            testCase.verifyTrue(contains(rL.estimand_note, 'remainder_scheme = legacy_observation_only_v1'));
            testCase.verifyTrue(contains(rL.estimand_note, 'has not been migrated'));
            testCase.verifyFalse(contains(rL.estimand_note, 'residual = pooled population highest-order'));
            testCase.verifyFalse(contains(rL.estimand_note, 'observation variance (scheme v2)'));
            testCase.verifyTrue(contains(rL.estimand_note, 'own log residual SD'));
            testCase.verifyEqual(rL.gamma_target.remainder_scheme, 'legacy_observation_only_v1');
            testCase.verifyEqual(rL.gamma_target.subject_scheme_source, 'legacy_unstamped');
            % and the provenance lines say the same from the REL alone
            L = string(psyrat_provenance_lines(TestGammaDynrelSurfaceTarget.localData(old)));
            testCase.verifyTrue(any(contains(L, 'LEGACY observation-only')));
            testCase.verifyFalse(any(contains(L, 'POOLED, scheme')));
            L2 = string(psyrat_provenance_lines(pd));
            testCase.verifyTrue(any(contains(L2, 'POOLED, scheme population_highest_order')));
            % a Gaussian run refuses the option outright
            gauss = PsyRATTestDataFactory.makeDynrelSubjNondiffSummaryData();
            testCase.verifyError(@() psyrat_report(gauss, 'surface_target', 'typical_person'), ...
                'dynrel:surfacetargetunsupported');
        end

        function testPlotCaptionAndUserDataRetainTheTarget(testCase)
            rel = TestGammaDynrelSurfaceTarget.localTrtRel('ls', true);
            sB = psyrat_dynrel_summary(rel, 'CI', .95, 'surface_target', 'typical_person', 'ngrid', 5);
            fh = psyrat_dynrelplot(sB, 'gcoeff', 2);
            c = onCleanup(@() close(fh, 'force'));
            ud = fh.UserData;
            testCase.assertTrue(isstruct(ud) && isfield(ud, 'dynrel'));
            testCase.verifyEqual(ud.dynrel.surface_target, 'typical_person');
            testCase.verifyEqual(ud.dynrel.gamma_target.remainder_scheme, ...
                'population_highest_order_plus_person_observation_v2');
            testCase.verifyEqual(ud.dynrel.gcoeff, 2);
            testCase.verifyEqual([ud.dynrel.strata.obs], [sB.strata.obs]);
            cap = findall(fh, 'Tag', 'psyrat_dynrel_target_caption');
            testCase.assertNumElements(cap, 1);
            % the caption must be an EXPORTABLE graphics annotation: exportgraphics
            % drops UI components, and a uicontrol caption was absent from the
            % saved PNG (root's image gate)
            testCase.verifyFalse(isa(cap, 'matlab.ui.control.UIControl'), ...
                'the caption must not be a uicontrol (excluded from PNG/PDF export)');
            testCase.verifyTrue(isa(cap, 'matlab.graphics.shape.TextBox'), ...
                'the caption must be a textbox annotation');
            txt = localCapText(cap);
            testCase.verifyTrue(contains(txt, 'TYPICAL-PERSON'));
            testCase.verifyTrue(contains(txt, 'population_highest_order'));
            testCase.verifyTrue(contains(txt, 'mixed population/person target'));
            testCase.verifyTrue(contains(txt, 'own effects in the observation variance only'));
            % the caption band sits below every axes (no overlap with x-labels)
            axs = findobj(fh, 'Type', 'axes');
            for a = axs(:)'
                p = get(a, 'Position');
                testCase.verifyGreaterThanOrEqual(p(2), cap.Position(2) + cap.Position(4) + 0.05, ...
                    'axes must be lifted clear of the caption band');
            end
            % the native save button is still a uicontrol (hidden during export)
            btn = findall(fh, 'Style', 'pushbutton');
            testCase.verifyTrue(any(contains(string(get(btn, 'String')), 'Save figure')));
            % legacy result: caption says LEGACY
            old = TestGammaDynrelSurfaceTarget.localTrtRel('ls', false);
            sL = psyrat_dynrel_summary(old, 'CI', .95, 'ngrid', 5);
            fhL = psyrat_dynrelplot(sL, 'gcoeff', 1);
            cL = onCleanup(@() close(fhL, 'force'));
            capL = findall(fhL, 'Tag', 'psyrat_dynrel_target_caption');
            testCase.verifyTrue(contains(localCapText(capL), 'LEGACY'));
            testCase.verifyTrue(contains(localCapText(capL), 'POPULATION-MARGINAL'));
            % Gaussian figure: no caption, no UserData
            gauss = PsyRATTestDataFactory.makeDynrelSubjNondiffSummaryData();
            sG = psyrat_dynrel_summary(gauss.rel, 'CI', .95, 'ngrid', 5);
            fhG = psyrat_dynrelplot(sG, 'gcoeff', 2);
            cG = onCleanup(@() close(fhG, 'force'));
            testCase.verifyEmpty(findall(fhG, 'Tag', 'psyrat_dynrel_target_caption'));
            testCase.verifyEmpty(fhG.UserData);
        end

        function testViewerControlSelectsTargetAndCacheRoundTrips(testCase)
            rel = TestGammaDynrelSurfaceTarget.localTrtRel('lognu', true);
            pd = TestGammaDynrelSurfaceTarget.localData(rel);
            fig = TestGammaDynrelSurfaceTarget.localViewerFig(testCase, pd);
            c = onCleanup(@() close(findall(0, 'Type', 'figure'), 'force'));
            pop = findobj(fig, 'Tag', 'psyrat_dynrel_estimand');
            testCase.assertNumElements(pop, 1);
            testCase.verifyEqual(pop.Value, 1, 'the primary is the default');
            testCase.verifyTrue(contains(pop.String{1}, 'Population-marginal'));
            testCase.verifyTrue(contains(pop.String{2}, 'Typical-person'));
            testCase.verifyEqual(pop.Enable, 'on');
            % A: primary
            fA = TestGammaDynrelSurfaceTarget.localFireButton(testCase, fig, 'Plot reliability surface');
            testCase.verifyEqual(fA.UserData.dynrel.surface_target, 'population_marginal');
            yA = get(findobj(fA, 'Type', 'line'), 'YData');
            % B: companion
            pop.Value = 2;
            fB = TestGammaDynrelSurfaceTarget.localFireButton(testCase, fig, 'Plot reliability surface');
            testCase.verifyEqual(fB.UserData.dynrel.surface_target, 'typical_person');
            yB = get(findobj(fB, 'Type', 'line'), 'YData');
            testCase.verifyNotEqual(yA, yB, 'the companion must draw a different curve');
            % A again: the cache must return A's surface, not B's
            pop.Value = 1;
            fA2 = TestGammaDynrelSurfaceTarget.localFireButton(testCase, fig, 'Plot reliability surface');
            testCase.verifyEqual(fA2.UserData.dynrel.surface_target, 'population_marginal');
            yA2 = get(findobj(fA2, 'Type', 'line'), 'YData');
            testCase.verifyEqual(yA2, yA, 'A -> B -> A must reproduce A exactly');
            % the per-participant window carries the on-screen scheme note and the
            % same table under both targets
            pop.Value = 2;
            fT = TestGammaDynrelSurfaceTarget.localFireButton(testCase, fig, 'Show per-participant table');
            note = findobj(fT, 'Tag', 'psyrat_dynrel_subject_target_note');
            testCase.assertNumElements(note, 1);
            testCase.verifyTrue(contains(string(note.String), 'scheme v2'));
            % Gaussian viewer keeps its own strings and default
            gauss = PsyRATTestDataFactory.makeDynrelSubjNondiffSummaryData();
            figG = TestGammaDynrelSurfaceTarget.localViewerFig(testCase, gauss);
            popG = findobj(figG, 'Tag', 'psyrat_dynrel_estimand');
            testCase.verifyTrue(contains(popG.String{1}, 'delta_p'));
            testCase.verifyEqual(popG.Value, 1);
        end

        % ---------- saved explicit target: precedence, reopen, refusal ----------

        function testStoredExplicitTargetReopensAndExplicitOptionOverrides(testCase)
            for param = {'lognu', 'ls'}
                fresh = TestGammaDynrelSurfaceTarget.localTrtRel(param{1}, true);                    % stamp: primary
                saved = TestGammaDynrelSurfaceTarget.localTrtRel(param{1}, true, 'typical_person'); % stamp: companion
                old   = TestGammaDynrelSurfaceTarget.localTrtRel(param{1}, false);                   % no stamp
                % resolver precedence: explicit > stored > historical
                [t, g] = psyrat_gamma_target_resolve(saved, [], [], []);
                testCase.verifyEqual(t, 'typical_person');
                testCase.verifyEqual(g.surface_target_source, 'stored_stamp');
                [t, g] = psyrat_gamma_target_resolve(saved, 'population_marginal', [], []);
                testCase.verifyEqual(t, 'population_marginal');
                testCase.verifyEqual(g.surface_target_source, 'explicit_option');
                [t, g] = psyrat_gamma_target_resolve(fresh, [], [], []);
                testCase.verifyEqual(t, 'population_marginal');
                testCase.verifyEqual(g.surface_target_source, 'stored_stamp');
                [t, g] = psyrat_gamma_target_resolve(old, [], [], []);
                testCase.verifyEqual(t, 'population_marginal');
                testCase.verifyEqual(g.surface_target_source, 'historical_default');
                testCase.verifyFalse(g.stamped);
                % summary default follows the stored choice; explicit primary overrides;
                % the per-participant table is identical under all of them
                sS = psyrat_dynrel_summary(saved, 'CI', .95);
                sP = psyrat_dynrel_summary(saved, 'CI', .95, 'surface_target', 'population_marginal');
                sF = psyrat_dynrel_summary(fresh, 'CI', .95);
                sT = psyrat_dynrel_summary(fresh, 'CI', .95, 'surface_target', 'typical_person');
                testCase.verifyEqual(sS.surface_target, 'typical_person');
                testCase.verifyEqual(sS.gamma_target.surface_target_source, 'stored_stamp');
                testCase.verifyEqual(sP.surface_target, 'population_marginal');
                testCase.verifyEqual(sF.surface_target, 'population_marginal');
                for s = 1:2
                    testCase.verifyEqual(sS.strata(s).G, sT.strata(s).G, 'stored companion == explicit companion');
                    testCase.verifyEqual(sP.strata(s).G, sF.strata(s).G, 'explicit primary == fresh primary');
                    testCase.verifyGreaterThan(max(abs(sS.strata(s).G.pt(:) - sP.strata(s).G.pt(:))), 1e-6);
                    testCase.verifyTrue(isequal(sS.strata(s).ssrel_table, sP.strata(s).ssrel_table));
                    testCase.verifyTrue(isequal(sS.strata(s).ssrel_table, sF.strata(s).ssrel_table));
                end
                % the legacy marginal Boolean is never read as a target: an unstamped
                % result opens as the population surface whatever marginal says
                sM = testCase.verifyWarning(@() psyrat_dynrel_summary(old, 'CI', .95, 'marginal', 1), ...
                    'dynrel:marginalgamma');
                testCase.verifyEqual(sM.surface_target, 'population_marginal');
                % headless report follows the stored choice and names it
                rS = psyrat_report(TestGammaDynrelSurfaceTarget.localData(saved));
                testCase.verifyEqual(rS.gamma_target.surface_target, 'typical_person');
                testCase.verifyTrue(contains(rS.estimand_note, 'TYPICAL PERSON'));
                rP = psyrat_report(TestGammaDynrelSurfaceTarget.localData(saved), 'surface_target', 'population_marginal');
                testCase.verifyEqual(rP.gamma_target.surface_target, 'population_marginal');
                testCase.verifyEqual(rP.tables.ssrel.gen_pt, rS.tables.ssrel.gen_pt);
                % viewer: reopens on choice 2 for the saved companion, choice 1 for
                % the fresh primary and for an unstamped legacy result; the first
                % plot honors the reopened choice
                c = onCleanup(@() close(findall(0, 'Type', 'figure'), 'force')); %#ok<NASGU>
                figS = TestGammaDynrelSurfaceTarget.localViewerFig(testCase, TestGammaDynrelSurfaceTarget.localData(saved));
                popS = findobj(figS, 'Tag', 'psyrat_dynrel_estimand');
                testCase.verifyEqual(popS.Value, 2, 'a saved typical_person must reopen as choice 2');
                fS = TestGammaDynrelSurfaceTarget.localFireButton(testCase, figS, 'Plot reliability surface');
                testCase.verifyEqual(fS.UserData.dynrel.surface_target, 'typical_person');
                figF = TestGammaDynrelSurfaceTarget.localViewerFig(testCase, TestGammaDynrelSurfaceTarget.localData(fresh));
                popF = findobj(figF, 'Tag', 'psyrat_dynrel_estimand');
                testCase.verifyEqual(popF.Value, 1, 'a fresh (stamped primary) result opens as choice 1');
                figO = TestGammaDynrelSurfaceTarget.localViewerFig(testCase, TestGammaDynrelSurfaceTarget.localData(old));
                popO = findobj(figO, 'Tag', 'psyrat_dynrel_estimand');
                testCase.verifyEqual(popO.Value, 1, 'an unstamped legacy result opens as choice 1');
                close(findall(0, 'Type', 'figure'), 'force');
            end
        end

        function testUnsupportedOrIncoherentStampIsRefusedNotDowngraded(testCase)
            base = TestGammaDynrelSurfaceTarget.localTrtRel('lognu', true);
            % unsupported version
            bad = base; bad.out.gamma_target.gamma_target_version = 99;
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetversion');
            bad = base; bad.out.gamma_target = rmfield(bad.out.gamma_target, 'gamma_target_version');
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetversion');
            % contradictory identifiers
            bad = base; bad.out.gamma_target.parameterization = 'ls';
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetincoherent');
            bad = base; bad.out.gamma_target.remainder_scheme = 'legacy_observation_only_v1';
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetincoherent');
            bad = base; bad.out.gamma_target.subject_target = 'not_applicable';
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetincoherent');
            bad = base; bad.out.gamma_target.count_source = 'mean_trials_per_occupied_occasion';
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetincoherent');
            bad = base; bad.out.gamma_target.analysis = 'ic_dynrel_trt';
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetincoherent');
            % a static-25-shaped stamp on a dynamic result
            bad = base; bad.out.gamma_target.surface_target = 'not_applicable';
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetincoherent');
            % a partial stamp
            bad = base; bad.out.gamma_target = rmfield(bad.out.gamma_target, 'surface_target');
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetincoherent');
            % and the report / provenance refuse the same way (no silent legacy)
            bad = base; bad.out.gamma_target.gamma_target_version = 99;
            testCase.verifyError(@() psyrat_report(TestGammaDynrelSurfaceTarget.localData(bad)), ...
                'dynrel:gammatargetversion');
            % positive control: the untouched stamp resolves
            [t, g] = psyrat_gamma_target_resolve(base, [], [], []);
            testCase.verifyEqual(t, 'population_marginal');
            testCase.verifyTrue(g.stamped);
        end

        function testFigureUserDataSummaryRoundTripsThroughSavedFig(testCase)
            rel = TestGammaDynrelSurfaceTarget.localTrtRel('ls', true);
            tmp = [tempname '.fig'];
            cleanFile = onCleanup(@() localDeleteIfExists(tmp));
            for target = {'population_marginal', 'typical_person'}
                summ = psyrat_dynrel_summary(rel, 'CI', .95, 'surface_target', target{1}, 'ngrid', 5);
                fh = psyrat_dynrelplot(summ, 'gcoeff', 2);
                testCase.assertTrue(isfield(fh.UserData, 'summary'));
                testCase.verifyTrue(isequal(fh.UserData.summary, summ), ...
                    'the figure must carry the exact summary it was drawn from');
                savefig(fh, tmp);
                close(fh, 'force');
                fh2 = openfig(tmp, 'invisible');
                c2 = onCleanup(@() close(fh2, 'force')); %#ok<NASGU>
                ud = fh2.UserData;
                testCase.assertTrue(isstruct(ud) && isfield(ud, 'summary') && isfield(ud, 'dynrel'));
                testCase.verifyEqual(ud.summary.gamma_target.surface_target, target{1});
                testCase.verifyEqual(ud.dynrel.surface_target, target{1});
                % the exportable caption and its tag survive the FIG round trip
                cap2 = findall(fh2, 'Tag', 'psyrat_dynrel_target_caption');
                testCase.assertNumElements(cap2, 1, 'caption annotation lost on reopen');
                testCase.verifyTrue(isa(cap2, 'matlab.graphics.shape.TextBox'));
                if strcmp(target{1}, 'typical_person')
                    testCase.verifyTrue(contains(localCapText(cap2), 'TYPICAL-PERSON'));
                else
                    testCase.verifyTrue(contains(localCapText(cap2), 'POPULATION-MARGINAL'));
                end
                testCase.verifyEqual(ud.summary.gamma_target.remainder_scheme, ...
                    'population_highest_order_plus_person_observation_v2');
                for s = 1:2
                    testCase.verifyEqual(ud.summary.strata(s).G.pt, summ.strata(s).G.pt);
                    testCase.verifyEqual(ud.summary.strata(s).D.pt, summ.strata(s).D.pt);
                    testCase.verifyEqual(ud.summary.strata(s).obs, summ.strata(s).obs);
                    testCase.verifyEqual(ud.summary.strata(s).nocc, summ.strata(s).nocc);
                    testCase.verifyTrue(isequal(ud.summary.strata(s).ssrel_table, summ.strata(s).ssrel_table));
                    testCase.verifyTrue(isequal(ud.summary.strata(s).varcomp, summ.strata(s).varcomp));
                end
                % the curves on the reopened figure are the saved target's curves
                y = get(findobj(fh2, 'Type', 'line'), 'YData');
                if ~iscell(y); y = {y}; end
                for s = 1:2
                    want = summ.strata(s).G.pt(:)';
                    hit = any(cellfun(@(v) isequal(size(v), size(want)) && ...
                        max(abs(v - want)) <= 1e-12 * max(1, max(abs(want))), y));
                    testCase.verifyTrue(hit, sprintf('stratum %d curve not found on the reopened figure', s));
                end
                clear c2
            end
            % the summary carries no raw draws
            testCase.verifyFalse(isfield(summ, 'out'));
        end

        function testMalformedMarkerAndGammascaleAreRefused(testCase)
            base = TestGammaDynrelSurfaceTarget.localTrtRel('lognu', true);
            old  = TestGammaDynrelSurfaceTarget.localTrtRel('lognu', false);
            % a PRESENT but empty / non-struct / non-scalar marker is refused,
            % never read as an absent (historical) marker
            bad = base; bad.out.gamma_target = [];
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetmalformed');
            bad = base; bad.out.gamma_target = 'population_marginal';
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetmalformed');
            bad = base; bad.out.gamma_target = [base.out.gamma_target base.out.gamma_target];
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetmalformed');
            bad = base; bad.out.gamma_target = struct();
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetversion');
            % the two identifiers the first cut did not validate
            bad = base; bad.out.gamma_target = rmfield(bad.out.gamma_target, 'reference_target');
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetincoherent');
            bad = base; bad.out.gamma_target = rmfield(bad.out.gamma_target, 'surface_count_source');
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetincoherent');
            bad = base; bad.out.gamma_target.reference_target = 'typical_person_zero_person_effects';
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetincoherent');
            bad = base; bad.out.gamma_target.surface_count_source = 'user_override';
            testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:gammatargetincoherent');
            % a summary-time n' override is reported as such without touching the
            % stored marker's validation
            s = psyrat_dynrel_summary(base, 'CI', .95, 'obs', 20);
            testCase.verifyEqual(s.gamma_target.surface_count_source, 'user_override');
            testCase.verifyTrue(s.gamma_target.stamped);
            % gammascale: present must be exactly 1 or 2 (stamped and unstamped)
            for gs = {3, NaN, [], '2', [1 2], 1i}
                bad = base; bad.gammascale = gs{1};
                testCase.verifyError(@() psyrat_dynrel_summary(bad, 'CI', .95), 'dynrel:badgammascale');
                badold = old; badold.gammascale = gs{1};
                testCase.verifyError(@() psyrat_dynrel_summary(badold, 'CI', .95), 'dynrel:badgammascale');
            end
            % absent gammascale on an unstamped result is the documented log-nu fallback
            oldnogs = rmfield(old, 'gammascale');
            sN = psyrat_dynrel_summary(oldnogs, 'CI', .95);
            testCase.verifyEqual(sN.gamma_target.parameterization, 'lognu');
            testCase.verifyEqual(sN.gamma_target.subject_scheme_source, 'legacy_unstamped');
            % but a stamped marker must stay coherent with it: an LS stamp on a
            % result without gammascale is contradictory, not silently log-nu
            badls = TestGammaDynrelSurfaceTarget.localTrtRel('ls', true);
            badls = rmfield(badls, 'gammascale');
            testCase.verifyError(@() psyrat_dynrel_summary(badls, 'CI', .95), 'dynrel:gammatargetincoherent');
            % the validator itself reports the same verdicts without throwing
            v = psyrat_gamma_target_validate(base);
            testCase.verifyEqual(v.status, 'valid');
            v = psyrat_gamma_target_validate(old);
            testCase.verifyEqual(v.status, 'absent');
            bad = base; bad.out.gamma_target = [];
            v = psyrat_gamma_target_validate(bad);
            testCase.verifyEqual(v.status, 'invalid');
            testCase.verifyEqual(v.errid, 'dynrel:gammatargetmalformed');
            % the viewer opens (choice 1) and the first plot raises the refusal
            pd = TestGammaDynrelSurfaceTarget.localData(bad);
            fig = TestGammaDynrelSurfaceTarget.localViewerFig(testCase, pd);
            c = onCleanup(@() close(findall(0, 'Type', 'figure'), 'force')); %#ok<NASGU>
            pop = findobj(fig, 'Tag', 'psyrat_dynrel_estimand');
            testCase.verifyEqual(pop.Value, 1);
            btns = findall(fig, 'Style', 'pushbutton');
            for b = btns(:)'
                sb = b.String; if iscell(sb), sb = strjoin(sb, ' '); end
                if contains(string(sb), 'Plot reliability surface')
                    cb = b.Callback;
                    testCase.verifyError(@() cb(b, []), 'dynrel:gammatargetmalformed');
                end
            end
        end

        function testProvenanceRefusesInvalidMarkerDirectly(testCase)
            % DIRECT psyrat_provenance_lines on static 25 (both scales) and
            % dynamic 27: a marker the validator rejects yields an explicit
            % uninterpretable-marker disclosure, never POOLED and never LEGACY
            nd = 4;
            fixtures = {};
            for param = {'lognu', 'ls'}
                rel = struct('analysis', 'trt_sserrvar', 'family', 'gamma', ...
                    'filename', 'static25_marker_fixture.psyrat');
                if strcmp(param{1}, 'ls')
                    rel.gammascale = 2;
                    rel.out.gamma_disp = struct('is_ls', true, 's_sd', repmat(.45, nd, 1), ...
                        'rho_ps', repmat(.3, nd, 1), 'disp_factor', repmat(exp(2*.45^2), nd, 1));
                else
                    rel.gammascale = 1;
                    rel.out.gamma_disp = struct('s_v', repmat(.5, nd, 1), ...
                        'rho_pv', repmat(2/3, nd, 1), 'disp_factor', repmat(exp(.125 - .4), nd, 1));
                end
                valid = struct('gamma_target_version', 2, ...
                    'surface_target', 'not_applicable', ...
                    'reference_target', 'typical_person_zero_person_effects', ...
                    'subject_target', 'own_observation_variance_with_population_mean_components_at_own_z', ...
                    'remainder_scheme', 'population_highest_order_plus_person_observation_v2', ...
                    'parameterization', param{1}, ...
                    'count_source', 'mean_trials_per_occupied_occasion', ...
                    'surface_count_source', 'not_applicable', ...
                    'analysis', 'trt_sserrvar');
                fixtures{end+1} = struct('rel', rel, 'valid', valid); %#ok<AGROW>
            end
            dyn = TestGammaDynrelSurfaceTarget.localTrtRel('lognu', true);
            dynrel = dyn; dynrel.out = rmfield(dyn.out, 'gamma_target');
            fixtures{end+1} = struct('rel', dynrel, 'valid', dyn.out.gamma_target);

            for k = 1:numel(fixtures)
                rel = fixtures{k}.rel;
                % bare version field: first cut printed POOLED
                bad = rel; bad.out.gamma_target = struct('gamma_target_version', 2);
                L = string(psyrat_provenance_lines(TestGammaDynrelSurfaceTarget.localData(bad)));
                testCase.verifyTrue(any(contains(L, 'INVALID or uninterpretable')), sprintf('fixture %d: bare version not flagged', k));
                testCase.verifyTrue(any(contains(L, 'marker present but invalid')));
                testCase.verifyFalse(any(contains(L, 'POOLED')));
                testCase.verifyFalse(any(contains(L, 'LEGACY observation-only')));
                testCase.verifyFalse(any(contains(L, 'LEGACY, unstamped')));
                % unsupported version: first cut printed LEGACY / no stamp
                bad = rel; bad.out.gamma_target = fixtures{k}.valid; bad.out.gamma_target.gamma_target_version = 99;
                L = string(psyrat_provenance_lines(TestGammaDynrelSurfaceTarget.localData(bad)));
                testCase.verifyTrue(any(contains(L, 'INVALID or uninterpretable')));
                testCase.verifyFalse(any(contains(L, 'LEGACY')));
                testCase.verifyFalse(any(contains(L, 'POOLED')));
                % present-but-empty marker
                bad = rel; bad.out.gamma_target = [];
                L = string(psyrat_provenance_lines(TestGammaDynrelSurfaceTarget.localData(bad)));
                testCase.verifyTrue(any(contains(L, 'INVALID or uninterpretable')));
                testCase.verifyFalse(any(contains(L, 'LEGACY')));
                % contradictory identifier (wrong reference_target / count source)
                bad = rel; bad.out.gamma_target = fixtures{k}.valid;
                bad.out.gamma_target.reference_target = 'something_else';
                L = string(psyrat_provenance_lines(TestGammaDynrelSurfaceTarget.localData(bad)));
                testCase.verifyTrue(any(contains(L, 'INVALID or uninterpretable')));
                testCase.verifyFalse(any(contains(L, 'POOLED')));
                % positive controls: valid -> POOLED with coherent identifiers;
                % absent -> LEGACY with legacy identifiers
                good = rel; good.out.gamma_target = fixtures{k}.valid;
                L = string(psyrat_provenance_lines(TestGammaDynrelSurfaceTarget.localData(good)));
                testCase.verifyTrue(any(contains(L, 'POOLED, scheme population_highest_order')));
                id = L(contains(L, 'Stored-result contract'));
                testCase.assertNumElements(id, 1);
                testCase.verifyTrue(contains(id, ['reference_target = ' fixtures{k}.valid.reference_target]));
                testCase.verifyTrue(contains(id, ['count_source = ' fixtures{k}.valid.count_source]));
                testCase.verifyTrue(contains(id, ['surface_target = ' fixtures{k}.valid.surface_target]));
                testCase.verifyFalse(any(contains(L, 'INVALID')));
                L = string(psyrat_provenance_lines(TestGammaDynrelSurfaceTarget.localData(rel)));
                testCase.verifyTrue(any(contains(L, 'LEGACY observation-only')));
                testCase.verifyTrue(any(contains(L, 'LEGACY, unstamped')));
                testCase.verifyFalse(any(contains(L, 'INVALID')));
            end
        end

        % ---------- active vs stored contract lines ----------

        function testActiveContractLineCarriesAllIdentifiersOnEveryCoveredDesign(testCase)
            % line 2 of the shared writer is the literal active contract for
            % 11/14/26/27, subject and non-subject outputs alike; line 3 stays
            % the per-participant note where one exists
            ids = {'gamma_target_version = 2', 'surface_target = ', 'subject_target = ', ...
                'remainder_scheme = ', 'parameterization = ', 'count_source = ', ...
                'surface_count_source = ', 'reference_target = ', 'analysis = '};
            % 27 (subject, two-facet) and 26 (subject, one-facet)
            r27 = TestGammaDynrelSurfaceTarget.localTrtRel('lognu', true);
            r26 = TestGammaDynrelSurfaceTarget.localOneFacetRel([1 2 3]);
            % 14 and 11 (group-level) derived from the same draws with the
            % group-level marker
            r14 = r27; r14.analysis = 'ic_dynrel_trt';
            r14.out = rmfield(r14.out, 'ssinfo');
            r14.out.gamma_target.analysis = 'ic_dynrel_trt';
            r14.out.gamma_target.subject_target = 'not_applicable';
            r14.out.gamma_target.remainder_scheme = 'population_mean_remainder_plus_selected_surface_observation_v1';
            r14.out.gamma_target.count_source = 'not_applicable';
            r11 = r26; r11.analysis = 'ic_dynrel';
            r11.out = rmfield(r11.out, 'ssinfo');
            fixtures = {r27, 'ic_dynrel_sserrvar_trt', 'population_highest_order_plus_person_observation_v2', ...
                    'median_trials_per_occupied_occasion', 'own_observation_variance_with_population_mean_components_at_own_z', true; ...
                r26, 'ic_dynrel_sserrvar', 'population_mean_remainder_plus_person_observation_v1', ...
                    'median_trials_per_occupied_occasion', 'own_observation_variance_with_population_mean_components_at_own_z', true; ...
                r14, 'ic_dynrel_trt', 'population_mean_remainder_plus_selected_surface_observation_v1', ...
                    'not_applicable', 'not_applicable', false; ...
                r11, 'ic_dynrel', 'population_mean_remainder_plus_selected_surface_observation_v1', ...
                    'not_applicable', 'not_applicable', false};
            for k = 1:size(fixtures,1)
                s = psyrat_dynrel_summary(fixtures{k,1}, 'CI', .95, 'ngrid', 3);
                for perpart = [false true]
                    L = psyrat_gamma_target_lines(s, perpart);
                    testCase.assertGreaterThanOrEqual(numel(L), 2, sprintf('%s: fewer than two lines', fixtures{k,2}));
                    act = L{2};
                    testCase.verifyTrue(startsWith(act, 'Active summary contract:'), ...
                        sprintf('%s: line 2 must be the active contract', fixtures{k,2}));
                    for j = 1:numel(ids)
                        testCase.verifyTrue(contains(act, ids{j}), ...
                            sprintf('%s: identifier %s missing from the active contract', fixtures{k,2}, ids{j}));
                    end
                    testCase.verifyTrue(contains(act, ['analysis = ' fixtures{k,2}]));
                    testCase.verifyTrue(contains(act, ['remainder_scheme = ' fixtures{k,3}]));
                    testCase.verifyTrue(contains(act, ['count_source = ' fixtures{k,4} ';']));
                    testCase.verifyTrue(contains(act, ['subject_target = ' fixtures{k,5} ';']));
                    testCase.verifyTrue(contains(act, 'surface_target = population_marginal'));
                    testCase.verifyTrue(contains(act, 'surface_count_source = default_stratum_median'));
                    testCase.verifyTrue(contains(act, 'reference_target = not_applicable'));
                    testCase.verifyTrue(contains(act, 'parameterization = lognu'));
                    if perpart && fixtures{k,6}
                        testCase.verifyTrue(startsWith(L{3}, 'Estimand: per-participant table'), ...
                            'line 3 must remain the per-participant note');
                    else
                        testCase.verifyEqual(numel(L), 2);
                    end
                end
                % the report's note (variance/surface outputs too) carries the same line
                rep = psyrat_report(TestGammaDynrelSurfaceTarget.localData(fixtures{k,1}));
                testCase.verifyTrue(contains(rep.estimand_note, 'Active summary contract:'));
                testCase.verifyTrue(contains(rep.estimand_note, ['remainder_scheme = ' fixtures{k,3}]));
                % and the stored-result contract is a separate line family in the provenance
                P = string(rep.provenance);
                if fixtures{k,6} && strcmp(fixtures{k,2}, 'ic_dynrel_sserrvar_trt')
                    testCase.verifyTrue(any(startsWith(P, 'Stored-result contract')));
                end
            end
        end

        function testActiveAndStoredContractsAreDistinguishedWhenTheyDiffer(testCase)
            % (a) stored primary, view-time typical: the active line says
            % typical_person, the stored-result line still says population_marginal
            rel = TestGammaDynrelSurfaceTarget.localTrtRel('ls', true);   % stamp: population_marginal
            pd = TestGammaDynrelSurfaceTarget.localData(rel);
            rep = psyrat_report(pd, 'surface_target', 'typical_person');
            note = rep.estimand_note;
            testCase.verifyTrue(contains(note, 'Active summary contract:'));
            testCase.verifyTrue(contains(note, 'surface_target = typical_person'));
            testCase.verifyTrue(contains(note, 'surface_target_source = explicit_option'));
            testCase.verifyFalse(contains(note, 'surface_target = population_marginal'));
            P = string(rep.provenance);
            stored = P(startsWith(P, 'Stored-result contract'));
            testCase.assertNumElements(stored, 1);
            testCase.verifyTrue(contains(stored, 'surface_target = population_marginal'));
            testCase.verifyTrue(contains(stored, 'before any view-time override'));
            testCase.verifyFalse(contains(stored, 'typical_person'));
            % the viewer's export header carries both families in the same order
            % (active from the summary, stored from the provenance) - checked on
            % the assembled lines rather than through the save dialog
            summ = psyrat_dynrel_summary(rel, 'CI', .95, 'surface_target', 'typical_person', 'ngrid', 3);
            act = psyrat_gamma_target_lines(summ, true);
            testCase.verifyTrue(contains(act{2}, 'surface_target = typical_person'));
            testCase.verifyTrue(contains(act{2}, 'remainder_scheme = population_highest_order_plus_person_observation_v2'));
            testCase.verifyTrue(contains(act{2}, 'subject_scheme_source = stamped_v2'));
            % (b) stored v2, explicit legacy subject scheme requested: the active
            % line says legacy + explicit_override, the stored line still says v2
            summX = psyrat_dynrel_summary(rel, 'CI', .95, 'ngrid', 3, ...
                'remainder_scheme', 'legacy_observation_only_v1');
            actX = psyrat_gamma_target_lines(summX, true);
            testCase.verifyTrue(contains(actX{2}, 'remainder_scheme = legacy_observation_only_v1'));
            testCase.verifyTrue(contains(actX{2}, 'subject_scheme_source = explicit_override'));
            testCase.verifyTrue(contains(actX{3}, 'explicitly requested'));
            L = string(psyrat_provenance_lines(pd));
            storedX = L(startsWith(L, 'Stored-result contract'));
            testCase.assertNumElements(storedX, 1);
            testCase.verifyTrue(contains(storedX, 'remainder_scheme = population_highest_order_plus_person_observation_v2'));
            testCase.verifyFalse(contains(storedX, 'legacy'));
            % (c) unstamped legacy: stored line is the legacy contract with the
            % historical primary, and it does not claim the active choice
            old = TestGammaDynrelSurfaceTarget.localTrtRel('ls', false);
            Lo = string(psyrat_provenance_lines(TestGammaDynrelSurfaceTarget.localData(old)));
            so = Lo(startsWith(Lo, 'Stored-result contract'));
            testCase.assertNumElements(so, 1);
            testCase.verifyTrue(contains(so, 'LEGACY, unstamped'));
            testCase.verifyTrue(contains(so, 'surface_target = population_marginal'));
            testCase.verifyTrue(contains(so, 'surface_count_source = default_stratum_median'));
            testCase.verifyTrue(contains(so, 'remainder_scheme = legacy_observation_only_v1'));
            so_act = psyrat_gamma_target_lines(psyrat_dynrel_summary(old, 'CI', .95, 'ngrid', 3, ...
                'surface_target', 'typical_person'), true);
            testCase.verifyTrue(contains(so_act{2}, 'surface_target = typical_person'));
            testCase.verifyTrue(contains(so_act{2}, 'subject_scheme_source = legacy_unstamped'));
        end

        function testTypicalPersonRelationWordingIsParameterizationSpecific(testCase)
            relL = TestGammaDynrelSurfaceTarget.localTrtRel('ls', true);
            relN = TestGammaDynrelSurfaceTarget.localTrtRel('lognu', true);
            lL = psyrat_gamma_target_lines(psyrat_dynrel_summary(relL, 'CI', .95, 'ngrid', 3, ...
                'surface_target', 'typical_person'), false);
            lN = psyrat_gamma_target_lines(psyrat_dynrel_summary(relN, 'CI', .95, 'ngrid', 3, ...
                'surface_target', 'typical_person'), false);
            % LS: w_p = 0 and the zero-scale-heterogeneity equality control; no
            % claim that zeroing the scale variance fails to give the companion
            testCase.verifyTrue(contains(lL{1}, 'w_p = 0'));
            testCase.verifyTrue(contains(lL{1}, 's_sd is zero'));
            testCase.verifyTrue(contains(lL{1}, 'equality control'));
            testCase.verifyFalse(contains(lL{1}, 'zeroing the person'));
            testCase.verifyFalse(contains(lL{1}, 'log-nu'));
            % log-nu: both effects zero; s_v = 0 alone is not the companion;
            % the exact equality condition is stated
            testCase.verifyTrue(contains(lN{1}, 'u_p = v_p = 0'));
            testCase.verifyTrue(contains(lN{1}, 'zeroing the person log-nu variance alone'));
            testCase.verifyTrue(contains(lN{1}, '4 s_p^2 + s_v^2 - 4 cov_pv = 0'));
            testCase.verifyFalse(contains(lN{1}, 's_sd'));
        end

        % ---------- caption wording and static-25 identifiers ----------

        function testCaptionWordingNamesJointDistributionAndMixedTarget(testCase)
            relN = TestGammaDynrelSurfaceTarget.localTrtRel('lognu', true);
            relL = TestGammaDynrelSurfaceTarget.localTrtRel('ls', true);
            sN = psyrat_dynrel_summary(relN, 'CI', .95, 'ngrid', 3);
            sL = psyrat_dynrel_summary(relL, 'CI', .95, 'ngrid', 3);
            lN = psyrat_gamma_target_lines(sN, true);
            lL = psyrat_gamma_target_lines(sL, true);
            % primary: what the observation variance is marginalized over
            testCase.verifyTrue(contains(lN{1}, 'JOINT person location'));
            testCase.verifyTrue(contains(lN{1}, 'log-nu distribution'));
            testCase.verifyFalse(contains(lN{1}, 'person scale population'));
            testCase.verifyTrue(contains(lL{1}, 'log residual-SD distribution'));
            testCase.verifyFalse(contains(lL{1}, 'person scale population'));
            % retained components named, not "all components plus remainder"
            testCase.verifyTrue(contains(lN{1}, 'six population mean components'));
            testCase.verifyTrue(contains(lN{1}, 'pooled population highest-order'));
            testCase.verifyFalse(contains(lN{1}, 'all population mean components'));
            % subject caption: mixed target, own effects only in the observation variance
            testCase.verifyTrue(contains(lN{3}, 'MIXED population/person target'));
            testCase.verifyTrue(contains(lN{3}, 'enter ONLY the observation variance'));
            testCase.verifyTrue(contains(lN{3}, 'not a fully person-conditioned'));
            testCase.verifyTrue(contains(lL{3}, 'their own log residual SD'));
            % one-facet names the person and trial mains plus the pooled PI term
            rel1 = TestGammaDynrelSurfaceTarget.localOneFacetRel([1 2 3]);
            l1 = psyrat_gamma_target_lines(psyrat_dynrel_summary(rel1, 'CI', .95, 'ngrid', 3), true);
            testCase.verifyTrue(contains(l1{1}, 'person and trial main-effect components'));
            testCase.verifyTrue(contains(l1{1}, 'pooled population person x trial'));
            testCase.verifyTrue(contains(l1{3}, 'scheme v1'));
        end

        function testLegacyCaptionStatesOmissionNotBias(testCase)
            old = TestGammaDynrelSurfaceTarget.localTrtRel('ls', false);
            fresh = TestGammaDynrelSurfaceTarget.localTrtRel('ls', true);
            sL = psyrat_dynrel_summary(old, 'CI', .95, 'ngrid', 3);
            sX = psyrat_dynrel_summary(fresh, 'CI', .95, 'ngrid', 3, ...
                'remainder_scheme', 'legacy_observation_only_v1');
            lL = psyrat_gamma_target_lines(sL, true);
            lX = psyrat_gamma_target_lines(sX, true);
            % unstamped: R omitted, "can be larger", equality only at R = 0
            testCase.verifyTrue(contains(lL{3}, 'OMITTED'));
            testCase.verifyTrue(contains(lL{3}, 'can be larger'));
            testCase.verifyTrue(contains(lL{3}, 'coincide when R is zero'));
            testCase.verifyFalse(contains(lL{3}, 'coincide only'), ...
                'zero signal also equalizes the two schemes, so "only" overstates');
            testCase.verifyTrue(contains(lL{3}, 'has not been migrated'));
            testCase.verifyFalse(contains(lL{3}, 'optimistically biased'));
            testCase.verifyFalse(contains(lL{3}, 'do not sit on either'));
            % explicit override on stamped data must not claim the file predates v2
            testCase.verifyTrue(contains(lX{3}, 'explicitly requested'));
            testCase.verifyFalse(contains(lX{3}, 'has not been migrated'));
            testCase.verifyFalse(contains(lX{3}, 'predates'));
            testCase.verifyTrue(contains(lX{3}, 'OMITTED'));
            % the figure caption agrees
            fh = psyrat_dynrelplot(sX, 'gcoeff', 2);
            c = onCleanup(@() close(fh, 'force')); %#ok<NASGU>
            cap = findall(fh, 'Tag', 'psyrat_dynrel_target_caption');
            testCase.verifyTrue(contains(localCapText(cap), 'explicitly requested'));
            testCase.verifyFalse(contains(localCapText(cap), 'not migrated'));
            fhL = psyrat_dynrelplot(sL, 'gcoeff', 2);
            cL = onCleanup(@() close(fhL, 'force')); %#ok<NASGU>
            capL = findall(fhL, 'Tag', 'psyrat_dynrel_target_caption');
            testCase.verifyTrue(contains(localCapText(capL), 'not migrated'));
            % provenance (REL-level) legacy line: omission wording, no bias claim
            L = string(psyrat_provenance_lines(TestGammaDynrelSurfaceTarget.localData(old)));
            leg = L(contains(L, 'LEGACY observation-only'));
            testCase.assertNumElements(leg, 1);
            testCase.verifyTrue(contains(leg, 'OMITTED'));
            testCase.verifyTrue(contains(leg, 'can be larger'));
            testCase.verifyFalse(contains(leg, 'optimistically biased'));
        end

        function testStaticTwoFacetProvenanceCarriesTargetIdentifiers(testCase)
            % static analysis 25 (both scales, stamped and unstamped): the
            % provenance lines - which every static CSV/XLSX export header and
            % the headless report's .provenance include - carry the explicit
            % target identifiers and the count source, not merely a POOLED word
            nd = 4;
            for param = {'lognu', 'ls'}
                rel = struct('analysis', 'trt_sserrvar', 'family', 'gamma', ...
                    'filename', 'static25_fixture.psyrat');
                if strcmp(param{1}, 'ls')
                    rel.gammascale = 2;
                    rel.out.gamma_disp = struct('is_ls', true, 's_sd', repmat(.45, nd, 1), ...
                        'rho_ps', repmat(.3, nd, 1), 'disp_factor', repmat(exp(2*.45^2), nd, 1));
                else
                    rel.gammascale = 1;
                    rel.out.gamma_disp = struct('s_v', repmat(.5, nd, 1), ...
                        'rho_pv', repmat(2/3, nd, 1), 'disp_factor', repmat(exp(.125 - .4), nd, 1));
                end
                % unstamped: legacy identifiers, marked unstamped
                L0 = string(psyrat_provenance_lines(TestGammaDynrelSurfaceTarget.localData(rel)));
                id0 = L0(contains(L0, 'Stored-result contract'));
                testCase.assertNumElements(id0, 1, sprintf('%s: identifier line missing (unstamped)', param{1}));
                testCase.verifyTrue(contains(id0, 'LEGACY, unstamped'));
                testCase.verifyTrue(contains(id0, 'surface_target = not_applicable'));
                testCase.verifyTrue(contains(id0, 'reference_target = typical_person_zero_person_effects'));
                testCase.verifyTrue(contains(id0, 'subject_target = own_observation_variance_with_population_mean_components_at_own_z'));
                testCase.verifyTrue(contains(id0, 'remainder_scheme = legacy_observation_only_v1'));
                testCase.verifyTrue(contains(id0, ['parameterization = ' param{1}]));
                testCase.verifyTrue(contains(id0, 'count_source = mean_trials_per_occupied_occasion'));
                testCase.verifyTrue(any(contains(L0, 'LEGACY observation-only')));
                testCase.verifyFalse(any(contains(L0, 'POOLED, scheme')));
                % stamped with the case-25 contract: identifiers read from the stamp
                rel.out.gamma_target = struct('gamma_target_version', 2, ...
                    'surface_target', 'not_applicable', ...
                    'reference_target', 'typical_person_zero_person_effects', ...
                    'subject_target', 'own_observation_variance_with_population_mean_components_at_own_z', ...
                    'remainder_scheme', 'population_highest_order_plus_person_observation_v2', ...
                    'parameterization', param{1}, ...
                    'count_source', 'mean_trials_per_occupied_occasion', ...
                    'surface_count_source', 'not_applicable', ...
                    'analysis', 'trt_sserrvar');
                pd = TestGammaDynrelSurfaceTarget.localData(rel);
                L2 = string(psyrat_provenance_lines(pd));
                id2 = L2(contains(L2, 'Stored-result contract'));
                testCase.assertNumElements(id2, 1, sprintf('%s: identifier line missing (stamped)', param{1}));
                testCase.verifyTrue(contains(id2, 'gamma_target_version 2'));
                testCase.verifyFalse(contains(id2, 'LEGACY'));
                testCase.verifyTrue(contains(id2, 'surface_target = not_applicable'));
                testCase.verifyTrue(contains(id2, 'reference_target = typical_person_zero_person_effects'));
                testCase.verifyTrue(contains(id2, 'remainder_scheme = population_highest_order_plus_person_observation_v2'));
                testCase.verifyTrue(contains(id2, ['parameterization = ' param{1}]));
                testCase.verifyTrue(contains(id2, 'count_source = mean_trials_per_occupied_occasion'));
                testCase.verifyTrue(contains(id2, 'surface_count_source = not_applicable'));
                testCase.verifyTrue(contains(id2, 'analysis = trt_sserrvar'));
                testCase.verifyTrue(contains(id2, 'before any view-time override'));
                testCase.verifyTrue(contains(id0, 'surface_count_source = not_applicable'));
                pooled = L2(contains(L2, 'POOLED, scheme population_highest_order'));
                testCase.assertNumElements(pooled, 1);
                testCase.verifyTrue(contains(pooled, 'MIXED population/person target'));
                testCase.verifyTrue(contains(pooled, 'six population mean components'));
                % a stamp missing a contract field is INVALID: nothing is invented
                % and no scheme is asserted
                rel2 = rel; rel2.out.gamma_target = rmfield(rel2.out.gamma_target, 'reference_target');
                L3 = string(psyrat_provenance_lines(TestGammaDynrelSurfaceTarget.localData(rel2)));
                testCase.verifyTrue(any(contains(L3, 'INVALID or uninterpretable')));
                testCase.verifyFalse(any(contains(L3, 'POOLED')));
                testCase.verifyFalse(any(contains(L3, 'LEGACY')));
            end
        end

        function testViewerShowsRestrictedNoteForGammaTwoFacetDesignsOnly(testCase)
            % Queue row 3a (owner rulings 2-4, 2026-09-21): the gamma two-facet
            % designs (14 and 27 here; 25 has its own viewer) are RESTRICTED
            % and must say so ON SCREEN, in the viewer's main window, not only
            % in the export header. The words are the shared
            % psyrat_gamma_twofacet_restriction text. The Gaussian arms and
            % the one-facet gamma designs (11, 26; ruling 4) stay unlabeled.
            % The layout invariant matches TestTrtSserrvarViewerRoute: every
            % control inside the figure and the window inside a 1728x1117
            % display, so the extra row cannot push the buttons off screen.
            R = psyrat_gamma_twofacet_restriction();
            tag = 'psyrat_dynrel_restriction_note';
            c = onCleanup(@() close(findall(0, 'Type', 'figure'), 'force')); %#ok<NASGU>

            % positives: 27 under both parameterizations, and 14 under gamma
            % (the 27 fixture with the analysis label swapped; the main
            % window reads no subject-level slot at build time, and the stamp
            % resolver sits in a try/catch)
            rel27ls = TestGammaDynrelSurfaceTarget.localTrtRel('ls', true);
            rel27nu = TestGammaDynrelSurfaceTarget.localTrtRel('lognu', true);
            rel14 = TestGammaDynrelSurfaceTarget.localTrtRel('ls', false);
            rel14.analysis = 'ic_dynrel_trt';
            pos = {rel27ls, rel27nu, rel14};
            for k = 1:numel(pos)
                pd = TestGammaDynrelSurfaceTarget.localData(pos{k});
                fig = TestGammaDynrelSurfaceTarget.localViewerFig(testCase, pd);
                note = findobj(fig, 'Tag', tag);
                testCase.verifyNumElements(note, 1, sprintf( ...
                    '%s gammascale %d: expected one on-screen RESTRICTED note', ...
                    pos{k}.analysis, pos{k}.gammascale));
                if ~isempty(note)
                    txt = string(note(1).String);
                    if numel(txt) > 1; txt = strjoin(txt, ' '); end
                    testCase.verifyEqual(char(txt), R.screen, ...
                        'the note must be the shared on-screen text');
                    testCase.verifyTrue(contains(txt, 'Restricted'));
                end
                TestGammaDynrelSurfaceTarget.localVerifyLayout(testCase, fig, pos{k}.analysis);
                close(fig, 'force');
            end

            % negatives: the Gaussian 14 arm (Supported) and the gamma
            % one-facet subject-level design 26 (outside the label, ruling 4)
            pdG = PsyRATTestDataFactory.makeDynrelTrtSummaryData();
            rel26 = TestGammaDynrelSurfaceTarget.localOneFacetRel([2; 1; 3]);
            neg = {pdG, TestGammaDynrelSurfaceTarget.localData(rel26)};
            for k = 1:numel(neg)
                fig = TestGammaDynrelSurfaceTarget.localViewerFig(testCase, neg{k});
                testCase.verifyEmpty(findobj(fig, 'Tag', tag), sprintf( ...
                    '%s must not carry the gamma two-facet RESTRICTED note', ...
                    neg{k}.rel.analysis));
                TestGammaDynrelSurfaceTarget.localVerifyLayout(testCase, fig, neg{k}.rel.analysis);
                close(fig, 'force');
            end
        end

        function testPerParticipantTargetNoteFitsItsBox(testCase)
            % G111 (owner ruling 2026-09-23: blocks 0.1.1-beta; the repair is
            % a taller note area). The per-participant window's target note
            % must show its whole text on screen. It wraps to 6-9 lines: the
            % analysis-27 LEGACY scheme text is the longest, and fsize-2 runs
            % 9-11 pt over psyrat_guifsize's 11-13 range. The fixed box it had
            % (about 53 px) showed about 3.5 lines in the Item 14 native look.
            % The line count is the note's own textwrap count; under -batch
            % that call reproduced the native look's line breaks exactly (6
            % and 7 lines, 2026-09-23). The line HEIGHT is not Extent's: a
            % text control's Extent adds one em per line, while the rendered
            % lines sit about 1.23 em apart (exportapp renders at 9 and 11
            % pt, 2026-09-23: 14.8 and 18.0 px against ems of 12 and 14.7
            % px), and that spacing reproduces the native 3.5 visible lines
            % where Extent predicts 4. So each line counts as 1.25 em, and
            % the box must hold the lines plus half a line. The em is read
            % from Extent, in a hidden text control with the note's font.
            % The note must also sit wholly below the table, with the table
            % no shorter than its former 311 px and the window inside a
            % 1728x1117 display. A window without the note keeps its
            % [200 200 920 380] geometry. Given equal line counts and ems,
            % the fit check follows from the production sizing rule, so it
            % guards that rule; the independent parts are the anchors (the
            % em is the font size at 96 px/in, and textwrap reports at least
            % the lines the unwrapped text's width needs) and the scheme
            % identifiers. Renders and the native look check the pixels.
            tag = 'psyrat_dynrel_subject_target_note';
            c = onCleanup(@() close(findall(0, 'Type', 'figure'), 'force')); %#ok<NASGU>
            v1 = 'remainder_scheme = population_mean_remainder_plus_person_observation_v1';
            v2 = 'remainder_scheme = population_highest_order_plus_person_observation_v2';
            lg = 'remainder_scheme = legacy_observation_only_v1';
            V = {'26 log-nu', TestGammaDynrelSurfaceTarget.localOneFacetRel([2; 1; 3]), v1; ...
                '27 v2 LS', TestGammaDynrelSurfaceTarget.localTrtRel('ls', true), v2; ...
                '27 v2 log-nu', TestGammaDynrelSurfaceTarget.localTrtRel('lognu', true), v2; ...
                '27 LEGACY LS', TestGammaDynrelSurfaceTarget.localTrtRel('ls', false), lg; ...
                '27 LEGACY log-nu', TestGammaDynrelSurfaceTarget.localTrtRel('lognu', false), lg};
            for fs = [11 13]
                prefs = psyrat_defaults();
                prefs.guis.fsize = fs;
                for k = 1:size(V, 1)
                    lbl = sprintf('%s, fsize %d', V{k, 1}, fs);
                    pd = TestGammaDynrelSurfaceTarget.localData(V{k, 2});
                    fig = TestGammaDynrelSurfaceTarget.localViewerFig(testCase, pd, prefs);
                    fT = TestGammaDynrelSurfaceTarget.localFireButton(testCase, fig, 'Show per-participant table');
                    note = findobj(fT, 'Tag', tag);
                    testCase.assertNumElements(note, 1, sprintf('%s: expected one target note', lbl));
                    % the text and font are unchanged: one char row holding
                    % the whole per-participant target, through its last
                    % sentence (the part the fixed box hid)
                    s = note.String;
                    testCase.verifyTrue(ischar(s) && isrow(s), sprintf('%s: the note must stay one char row', lbl));
                    testCase.verifyTrue(startsWith(s, 'Estimand: per-participant table = MIXED'), lbl);
                    testCase.verifyTrue(contains(s, V{k, 3}), lbl);
                    testCase.verifyTrue(endsWith(s, 'selected occasion count as shown per stratum'), lbl);
                    testCase.verifyEqual(note.FontSize, fs - 2, lbl);
                    % the fit: the wrapped lines at the rendered line
                    % spacing (1.25 em), plus half a line
                    np = getpixelposition(note);
                    n = numel(textwrap(note, {s}));
                    em = localTextHeight(note, {'X'; 'X'}, np(3)) - localTextHeight(note, {'X'}, np(3));
                    need = (n + 0.5) * 1.25 * em;
                    % anchors outside the production rule: a degenerate
                    % Extent (em near 0) or a textwrap undercount would
                    % otherwise shrink the box and the need together
                    fontpx = note.FontSize * 96/72;
                    testCase.verifyGreaterThanOrEqual(em, 0.9 * fontpx, sprintf( ...
                        '%s: em %.2f px is not the %g pt font size (%.2f px)', lbl, em, note.FontSize, fontpx));
                    testCase.verifyLessThanOrEqual(em, 1.1 * fontpx, sprintf( ...
                        '%s: em %.2f px is not the %g pt font size (%.2f px)', lbl, em, note.FontSize, fontpx));
                    wide = localTextExtent(note, s, np(3));
                    testCase.verifyGreaterThanOrEqual(n, ceil(wide(3) / np(3)), sprintf( ...
                        '%s: textwrap reports %d lines, but the unwrapped text is %.0f px in a %.0f px box', ...
                        lbl, n, wide(3), np(3)));
                    testCase.verifyGreaterThanOrEqual(np(4), need, sprintf( ...
                        ['%s: %d wrapped lines at 1.25 em (em = %.2f px) plus half a line ' ...
                        'need %.1f px, but the box is %.1f px'], lbl, n, em, need, np(4)));
                    % the layout
                    tb = findobj(fT, 'Type', 'uitable');
                    testCase.assertNumElements(tb, 1, lbl);
                    tp = getpixelposition(tb);
                    fp = fT.Position;
                    testCase.verifyGreaterThanOrEqual(np(2), 0, sprintf('%s: the note starts below the window', lbl));
                    testCase.verifyLessThanOrEqual(np(2) + np(4), tp(2), sprintf('%s: the note overlaps the table', lbl));
                    testCase.verifyLessThanOrEqual(tp(2) + tp(4), fp(4), sprintf('%s: the table ends above the window', lbl));
                    testCase.verifyGreaterThanOrEqual(tp(4), 311, sprintf('%s: the table shrank to %.1f px', lbl, tp(4)));
                    testCase.verifyLessThanOrEqual(fp(2) + fp(4), 1117, sprintf( ...
                        '%s: the window must fit a 1728x1117 display (origin %.0f + height %.0f)', lbl, fp(2), fp(4)));
                    close(fT, 'force');
                    close(fig, 'force');
                end
            end
            % negative: a per-participant window without the note (the
            % Gaussian analysis-26 arm carries no gamma target) is unchanged
            pdG = PsyRATTestDataFactory.makeDynrelSubjNondiffSummaryData();
            figG = TestGammaDynrelSurfaceTarget.localViewerFig(testCase, pdG);
            fG = TestGammaDynrelSurfaceTarget.localFireButton(testCase, figG, 'Show per-participant table');
            testCase.verifyEmpty(findobj(fG, 'Tag', tag), 'the Gaussian per-participant window carries no target note');
            testCase.verifyEqual(fG.Position, [200 200 920 380], 'a window without the note keeps its geometry');
        end

    end
end

function localDeleteIfExists(f)
if exist(f, 'file') == 2
    delete(f);
end
end

function t = localCapText(h)
%The caption's text as one string (its String is a cellstr, one line each).
t = strjoin(string(h.String), ' ');
end

function px = localTextHeight(ref, lines, width)
%Pixel height of LINES (a cellstr, one entry per line) in a hidden text
%control with REF's font (see localTextExtent).
e = localTextExtent(ref, lines, width);
px = e(4);
end

function e = localTextExtent(ref, str, width)
%Extent of STR (a cellstr of lines, or one char row, which Extent measures
%unwrapped) set in a hidden text control on REF's figure with REF's font,
%WIDTH px wide. The explicit empty CreateFcn keeps any figure-level autofit
%(psyrat_autofitbuttontext) out of the measurement, as in
%TestUitableWidthSlack.
c = uicontrol(ancestor(ref, 'figure'), 'Style', 'text', 'Units', 'pixels', ...
    'Position', [1 1 width 4000], 'HorizontalAlignment', 'left', ...
    'FontName', ref.FontName, 'FontUnits', ref.FontUnits, 'FontSize', ref.FontSize, ...
    'String', str, 'Visible', 'off', 'CreateFcn', '');
e = get(c, 'Extent');
delete(c);
end
