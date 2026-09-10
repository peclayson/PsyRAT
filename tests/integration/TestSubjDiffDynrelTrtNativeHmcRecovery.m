classdef TestSubjDiffDynrelTrtNativeHmcRecovery < PsyRATTestBase
    % Native HMC (engine = 3) known-truth recovery for the SUBJECT-LEVEL trial +
    % occasion two-facet dynamic difference reliability model (Rast & Clayson,
    % analysis 16; dynrel=2 + diffest=2 + sserrvar=2 + a time/occasion facet).
    % This is the group-level two-facet dynamic difference (case 15) with the
    % person (id) factor widened to a 4-D joint location + log-residual block (a
    % 4x4 LKJ Cholesky, as case 13), so each participant has their own per-event
    % residual; the other five crossed factors stay 2-D. Matches
    % psyrat_build_stan_diffdynrel_sserr_trt.
    %
    % The CmdStan path is validated by TestSubjDiffDynrelTrtRecovery (nsub=72). To
    % bound this high-dimensional native model's runtime, the generator here uses
    % INDEPENDENT person/facet event LOCATION effects (so the LKJ(2) prior agrees
    % with the truth and the difference universe variance recovers at modest N),
    % while keeping the per-person residual-SCALE random effect that defines the
    % subject-level variant. The well-identified components, the between-person
    % residual-scale SDs (sd_id columns 3-4), and the relative G_delta(z) surface
    % are asserted.
    %
    % Requires NO CmdStan (pure MATLAB hmcSampler + dlarray AD gradients); skips
    % cleanly (Incomplete) when either toolbox is absent. It is the most
    % parameter-rich native model (4-D person block + five crossed factors), so it
    % is slow and lives in tests/integration.

    methods (Test)
        function testRecoversTruth(testCase)
            av = psyrat_estimation_engines_available();
            testCase.assumeTrue(av.hmc, ...
                'Skipped: Statistics and Machine Learning Toolbox (hmcSampler) not available.');
            testCase.assumeTrue(av.ad, ...
                'Skipped: Deep Learning Toolbox (dlarray AD gradients) not available.');

            % ---- known truth (INDEPENDENT person/facet event location effects;
            % per-person residual-scale random effect) ----
            rng(20260624, 'twister');
            nsub = 25; nocc = 3; ntrl = 4;
            mu = [5 2]; bmean = [1.0 0.5];
            a = [log(3) log(2.5)]; bsig = [0.25 -0.10];
            sdP = [3 2.5];                    % person event LOCATION SDs (independent)
            sds = [0.5 0.4];                  % between-person residual-SCALE SDs
            occSD = [1.0 0.8]; trlSD = [1.0 0.8];
            ptSD  = [1.0 0.9]; poSD  = [0.9 0.8]; toSD = [0.5 0.4];
            uniTrue = sdP(1)^2 + sdP(2)^2;    % sigma^2_p_delta (independent events)
            bpiTrue = ptSD(1)^2 + ptSD(2)^2;
            bpoTrue = poSD(1)^2 + poSD(2)^2;

            zraw = randn(nsub,1); zsub = (zraw - mean(zraw)) / std(zraw);
            nu_p = [randn(nsub,1)*sdP(1), randn(nsub,1)*sdP(2)];
            delta = [randn(nsub,1)*sds(1), randn(nsub,1)*sds(2)]; % person scale RE
            occE = [randn(nocc,1)*occSD(1), randn(nocc,1)*occSD(2)];
            trlE = [randn(ntrl,1)*trlSD(1), randn(ntrl,1)*trlSD(2)];
            ptE  = cat(3, randn(nsub,ntrl)*ptSD(1), randn(nsub,ntrl)*ptSD(2));
            poE  = cat(3, randn(nsub,nocc)*poSD(1), randn(nsub,nocc)*poSD(2));
            toE  = cat(3, randn(ntrl,nocc)*toSD(1), randn(ntrl,nocc)*toSD(2));

            ids = {}; meas = []; ev = {}; tm = {}; dim1 = []; events = {'gain','loss'};
            for p = 1:nsub
                for e = 1:2
                    for o = 1:nocc
                        for t = 1:ntrl
                            mn = mu(e) + bmean(e)*zsub(p) + nu_p(p,e) ...
                                + occE(o,e) + trlE(t,e) ...
                                + ptE(p,t,e) + poE(p,o,e) + toE(t,o,e);
                            sg = exp(a(e) + bsig(e)*zsub(p) + delta(p,e));
                            ids{end+1,1} = sprintf('S%02d', p); %#ok<AGROW>
                            ev{end+1,1} = events{e};            %#ok<AGROW>
                            tm{end+1,1} = sprintf('t%d', o);    %#ok<AGROW>
                            meas(end+1,1) = mn + randn*sg;      %#ok<AGROW>
                            dim1(end+1,1) = zsub(p);            %#ok<AGROW>
                        end
                    end
                end
            end
            T = table(ids, meas, ev, tm, dim1, ...
                'VariableNames', {'id', 'meas', 'event', 'time', 'dim1'});

            % ---- run native HMC through the full pipeline ----
            rel = psyrat_computevarcomp( ...
                'data', T, ...
                'dynrel', 2, ...
                'diffest', 2, ...
                'sserrvar', 2, ...
                'engine', 3, ...
                'chains', 2, ...
                'warmup', 250, ...
                'sampling', 250, ...
                'seed', 12345, ...
                'verbose', 1, ...
                'showgui', 1);

            testCase.verifyEqual(rel.analysis, 'ic_diff_dynrel_sserrvar_trt');
            testCase.verifyEqual(rel.engine, 'hmc');

            % ---- extracted draws ----
            idvc  = cell2mat(rel.out.id_varcov{1});
            tidvc = cell2mat(rel.out.tid_varcov{1});
            oidvc = cell2mat(rel.out.oid_varcov{1});
            bsg = cell2mat(rel.out.b_sigma{1});
            bsd = cell2mat(rel.out.b_sigma_dim{1});
            sdid = cell2mat(rel.out.sd_id{1});
            nd = size(bsg, 1);
            testCase.verifyEqual(size(sdid), [nd 4]);
            testCase.verifyEqual(size(idvc), [nd 2 2]);

            % ---- well-identified components (generous bands) ----
            uniHat = mean(idvc(:,1,1) + idvc(:,2,2) - 2*idvc(:,2,1));
            testCase.verifyEqual(uniHat, uniTrue, 'AbsTol', max(8.0, 0.4*uniTrue));
            testCase.verifyEqual(mean(bsg(:,1)), a(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsg(:,2)), a(2), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(bsd(:,1)), bsig(1), 'AbsTol', 0.20);
            testCase.verifyEqual(mean(bsd(:,2)), bsig(2), 'AbsTol', 0.20);

            % ---- the subject-level addition: between-person residual-scale SDs
            % (sd_id columns 3-4) -- the defining quantity of this variant ----
            testCase.verifyEqual(mean(sdid(:,3)), sds(1), 'AbsTol', 0.30);
            testCase.verifyEqual(mean(sdid(:,4)), sds(2), 'AbsTol', 0.30);

            % ---- typical-person reference surface G_delta(z) vs closed form ----
            zg = [-1.5 0 1.5];
            out = psyrat_rel_diffdynrel_trt('id_varcov', idvc, ...
                'trl_varcov', cell2mat(rel.out.trl_varcov{1}), ...
                'occ_varcov', cell2mat(rel.out.occ_varcov{1}), ...
                'tid_varcov', tidvc, 'oid_varcov', oidvc, ...
                'to_varcov', cell2mat(rel.out.to_varcov{1}), ...
                'b_sigma', bsg, 'b_sigma_dim', bsd, 'ndim', 1, 'z1', zg, ...
                'obs', ntrl, 'nocc', nocc, 'reltype', 3, 'CI', .95);
            for k = 1:numel(zg)
                wp1 = exp(a(1) + bsig(1)*zg(k))^2;   % typical person: delta = 0
                wp2 = exp(a(2) + bsig(2)*zg(k))^2;
                res = wp1 + wp2;
                Gtr = uniTrue / (uniTrue + bpiTrue/ntrl + bpoTrue/nocc ...
                    + res/(ntrl*nocc));
                testCase.verifyEqual(out.G.pt(k), Gtr, 'AbsTol', 0.12);
            end
        end
    end
end
