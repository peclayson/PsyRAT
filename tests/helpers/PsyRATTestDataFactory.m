classdef PsyRATTestDataFactory
    % Reusable deterministic test fixtures for PsyRAT toolbox tests.

    methods (Static)
        function psyrat_data = makeICSummaryData()
            ids = {'s1'; 's2'; 's3'; 's4'};
            id2 = (1:4)';
            trls = [3; 5; 6; 8];
            dep_pt = [0.65; 0.78; 0.84; 0.90];
            dep_ll = dep_pt - 0.06;
            dep_ul = dep_pt + 0.05;
            icc_pt = [0.40; 0.50; 0.62; 0.70];
            icc_ll = icc_pt - 0.07;
            icc_ul = icc_pt + 0.06;
            bp_var = repmat(0.25, size(id2));
            ss_errvar = [0.55; 0.40; 0.32; 0.28];
            pop_errvar = repmat(0.50, size(id2));

            ssrel_table = table(ids, id2, trls, ...
                dep_pt, dep_ll, dep_ul, ...
                icc_pt, icc_ll, icc_ul, ...
                bp_var, ss_errvar, pop_errvar);

            ev = struct();
            ev.name = 'measure';
            ev.trlcutoff = 4;
            ev.rel = struct('m', 0.82, 'll', 0.73, 'ul', 0.90);
            ev.dep = struct('m', 0.85, 'll', 0.76, 'ul', 0.92, 'meas', mean(trls));
            ev.trlinfo = struct( ...
                'min', min(trls), ...
                'max', max(trls), ...
                'mean', mean(trls), ...
                'med', median(trls), ...
                'std', std(trls));
            ev.goodn = 3;
            ev.icc = struct('m', 0.56, 'll', 0.43, 'ul', 0.68);
            ev.betsd = struct('m', 0.51, 'll', 0.41, 'ul', 0.62);
            ev.witsd = struct('m', 1.02, 'll', 0.91, 'ul', 1.13);
            ev.ssrel_table = ssrel_table;

            psyrat_data = struct();
            psyrat_data.rel = struct('groups', 'none', 'events', 'none', 'analysis', 'ic');

            psyrat_data.relsummary = struct();
            psyrat_data.relsummary.gcoeff = 1;
            psyrat_data.relsummary.gcoeff_name = 'dep';
            psyrat_data.relsummary.depcutoff = 0.80;
            psyrat_data.relsummary.meascutoff = 'Point Estimate';
            psyrat_data.relsummary.ciperc = 0.95;
            psyrat_data.relsummary.group = struct();
            psyrat_data.relsummary.group(1).name = '';
            psyrat_data.relsummary.group(1).event = ev;
            psyrat_data.relsummary.group(1).goodids = ids(2:4);
            psyrat_data.relsummary.group(1).badids = ids(1);

            psyrat_data.relsummary.data = struct();
            psyrat_data.relsummary.data.g(1).e(1).mu.raw = [0.08; 0.10; 0.12; 0.11];
            psyrat_data.relsummary.data.g(1).e(1).sig_u.raw = [0.40; 0.50; 0.55; 0.60];
            psyrat_data.relsummary.data.g(1).e(1).sig_e.raw = [0.90; 1.00; 1.10; 0.95];
        end

        function psyrat_data = makeICSummaryDataWithTrial()
            %One-facet summary fixture carrying a nonzero trial main effect
            %(sig_trl = sigma_i). Used to exercise the criterion (cut-score)
            %wiring so the dependability absolute error includes sigma_i^2
            %per Rocha 2026 Table 2. Built from makeICSummaryData so the only
            %difference is the added sig_trl draw and a small trial count
            %(which amplifies the per-observation sigma_i^2 contribution).
            psyrat_data = PsyRATTestDataFactory.makeICSummaryData();
            psyrat_data.relsummary.data.g(1).e(1).sig_trl.raw = ...
                [0.50; 0.60; 0.55; 0.45];

            trls = [3; 4; 5; 4];
            psyrat_data.relsummary.group(1).event(1).trlinfo.min = min(trls);
            psyrat_data.relsummary.group(1).event(1).trlinfo.max = max(trls);
            psyrat_data.relsummary.group(1).event(1).trlinfo.mean = mean(trls);
            psyrat_data.relsummary.group(1).event(1).trlinfo.med = median(trls);
            psyrat_data.relsummary.group(1).event(1).trlinfo.std = std(trls);
        end

        function psyrat_data = makeICDiffSummaryData()
            %One-facet difference-score summary fixture (ic_diff, one group x
            %two events), shaped like the psyrat_relsummary 'sing_diff' output
            %psyrat_depoverallt consumes: per-event dep/trlinfo/goodn rows plus
            %the group-level diffscore whose overall-table row carries '---'
            %placeholder cells in the badn and trial-summary columns (the G47
            %CSV-export surface).
            %trlsB has an even spread so its median is FRACTIONAL (7.5):
            %the CSV loop's %d conversion renders that in scientific
            %notation, a reachable rendering the byte-identical pins must
            %cover. evB.goodn differs from evA.goodn on purpose, so the
            %diff row's goodn cell (read from event(2) by the builder) is
            %a real differential against a wrong event index.
            trlsA = [4; 6; 6; 8];
            trlsB = [5; 7; 8; 9];

            evA = struct();
            evA.name = 'loss';
            evA.trlcutoff = 4;
            evA.dep = struct('m', 0.85, 'll', 0.76, 'ul', 0.92);
            evA.trlinfo = struct('min', min(trlsA), 'max', max(trlsA), ...
                'mean', mean(trlsA), 'med', median(trlsA), 'std', std(trlsA));
            evA.goodn = 3;

            evB = evA;
            evB.name = 'gain';
            evB.dep = struct('m', 0.88, 'll', 0.80, 'ul', 0.94);
            evB.trlinfo = struct('min', min(trlsB), 'max', max(trlsB), ...
                'mean', mean(trlsB), 'med', median(trlsB), 'std', std(trlsB));
            evB.goodn = 2;

            psyrat_data = struct();
            psyrat_data.rel = struct('groups', 'none', ...
                'events', {{'loss'; 'gain'}}, 'analysis', 'ic_diff');

            psyrat_data.relsummary = struct();
            psyrat_data.relsummary.gcoeff = 1;
            psyrat_data.relsummary.gcoeff_name = 'dep';
            psyrat_data.relsummary.depcutoff = 0.80;
            psyrat_data.relsummary.meascutoff = 'Point Estimate';
            psyrat_data.relsummary.ciperc = 0.95;
            psyrat_data.relsummary.group = struct();
            psyrat_data.relsummary.group(1).name = '';
            psyrat_data.relsummary.group(1).event = [evA, evB];
            psyrat_data.relsummary.group(1).goodids = {'s2'; 's3'; 's4'};
            psyrat_data.relsummary.group(1).badids = {'s1'};
            %the draw-based difference-score coefficient (a real value; the
            %diff row keeps this cell numeric-looking while dashing its counts)
            psyrat_data.relsummary.group(1).diffscore = ...
                struct('pt', 0.72, 'll', 0.61, 'ul', 0.81);
        end

        function psyrat_data = makeTRTSummaryData()
            psyrat_data = struct();
            psyrat_data.rel = struct('groups', 'none', 'events', 'none', 'analysis', 'trt');

            ev = struct();
            ev.name = 'measure';
            ev.trlcutoff = 5;
            ev.relcutoff = struct('m', 0.81, 'll', 0.71, 'ul', 0.89);
            ev.rel = struct('m', 0.79, 'll', 0.69, 'ul', 0.88);
            ev.trlinfo = struct('min', 4, 'max', 9, 'mean', 6.2, 'med', 6, 'std', 1.9);
            ev.goodn = 8;
            ev.icc = struct('m', 0.52, 'll', 0.40, 'ul', 0.65);
            ev.betsd = struct('m', 0.62, 'll', 0.51, 'ul', 0.73);
            ev.witsd = struct('m', 1.10, 'll', 0.96, 'ul', 1.24);

            psyrat_data.relsummary = struct();
            psyrat_data.relsummary.gcoeff = 1;
            psyrat_data.relsummary.reltype = 1;
            psyrat_data.relsummary.gcoeff_name = 'dep';
            psyrat_data.relsummary.reltype_name = 'ic';
            psyrat_data.relsummary.relcutoff = 0.80;
            psyrat_data.relsummary.meascutoff = 'Point Estimate';
            psyrat_data.relsummary.ciperc = 0.95;
            psyrat_data.relsummary.group = struct();
            psyrat_data.relsummary.group(1).name = '';
            psyrat_data.relsummary.group(1).event = ev;
            psyrat_data.relsummary.group(1).goodids = {'s1'; 's2'; 's3'; 's4'; 's5'; 's6'; 's7'; 's8'};
            psyrat_data.relsummary.group(1).badids = {'s9'; 's10'};

            psyrat_data.relsummary.data = struct();
            psyrat_data.relsummary.data.g(1).e(1).mu.raw = [0.10; 0.12; 0.11; 0.09];
            psyrat_data.relsummary.data.g(1).e(1).sig_id.raw = [0.40; 0.44; 0.48; 0.52];
            psyrat_data.relsummary.data.g(1).e(1).sig_occ.raw = [0.10; 0.12; 0.11; 0.10];
            psyrat_data.relsummary.data.g(1).e(1).sig_trl.raw = [0.20; 0.25; 0.23; 0.21];
            psyrat_data.relsummary.data.g(1).e(1).sig_trlxid.raw = [0.15; 0.16; 0.17; 0.15];
            psyrat_data.relsummary.data.g(1).e(1).sig_occxid.raw = [0.12; 0.13; 0.14; 0.12];
            psyrat_data.relsummary.data.g(1).e(1).sig_trlxocc.raw = [0.09; 0.10; 0.11; 0.09];
            psyrat_data.relsummary.data.g(1).e(1).sig_err.raw = [0.70; 0.75; 0.72; 0.68];
        end

        function psyrat_data = makeICDiffSSErrRelData()
            ndraw = 5;
            ids = {'s1'; 's2'; 's3'};
            id2 = (1:3)';
            events = {'loss'; 'gain'};

            idVarcov = zeros(ndraw, 2, 2);
            trlVarcov = zeros(ndraw, 2, 2);
            bSigma = zeros(ndraw, 2);
            errVarcov = zeros(ndraw, 2, 2);
            erVarSs1 = zeros(ndraw, 3);
            erVarSs2 = zeros(ndraw, 3);
            wpCovSs = zeros(ndraw, 3);

            for d = 1:ndraw
                drift = (d - 1) * 0.01;
                idVarcov(d,:,:) = [0.60 + drift 0.12; 0.12 0.54 + drift];
                trlVarcov(d,:,:) = [0.30 + drift 0.05; 0.05 0.28 + drift];
                bSigma(d,:) = log([0.70 + drift 0.76 + drift]);
                errVarcov(d,:,:) = [0.49 + drift 0.08; 0.08 0.58 + drift];
                erVarSs1(d,:) = log([0.50 + drift 0.75 + drift 0.95 + drift]);
                erVarSs2(d,:) = log([0.45 + drift 0.65 + drift 0.85 + drift]);
                wpCovSs(d,:) = [0.01 0.04 0.08] + drift;
            end

            rowIds = [repmat(ids(1), 11, 1); repmat(ids(2), 13, 1); repmat(ids(3), 15, 1)];
            rowEvents = [repmat(events(1), 6, 1); repmat(events(2), 5, 1); ...
                repmat(events(1), 7, 1); repmat(events(2), 6, 1); ...
                repmat(events(1), 8, 1); repmat(events(2), 7, 1)];
            meas = (1:numel(rowIds))' ./ 10;

            rel = struct();
            rel.analysis = 'ic_diff_sserrvar';
            rel.groups = 'none';
            rel.events = events;
            rel.diff_names = events;
            rel.diffwpcov = 2;
            rel.data = table(rowIds, meas, rowEvents, ...
                'VariableNames', {'id','meas','event'});
            rel.out = struct();
            rel.out.b_sigma(:,1) = {num2cell(bSigma)};
            rel.out.id_varcov(:,1) = {num2cell(idVarcov)};
            rel.out.trl_varcov(:,1) = {num2cell(trlVarcov)};
            rel.out.err_varcov(:,1) = {num2cell(errVarcov)};
            rel.out.er_var_ss1(:,1) = {erVarSs1};
            rel.out.er_var_ss2(:,1) = {erVarSs2};
            rel.out.wp_cov_ss(:,1) = {wpCovSs};
            rel.out.id_matches(:,1) = {table(ids, id2, 'VariableNames', {'id','id2'})};
            rel.out.elabels(:,1) = {events};
            rel.out.glabels(:,1) = {'none'};
            rel.out.conv.data = {{'name','n_eff','Rhat'; 'b_sigma', 200, 1.00}};

            psyrat_data = struct('rel', rel);
        end

        function psyrat_data = makeTRTSSErrRelData()
            %Subject-level crossed test-retest fixture (trt_sserrvar, analysis 25):
            %a synthetic REL with the sserr person block (cell-wrapped gro_sds /
            %ind_sdlog) plus the five crossed facet SDs and a crossed id/meas/time
            %data table. Per-subject log-residual offsets differ across subjects so
            %the per-participant coefficient is not constant. Consumed by running
            %psyrat_relsummary(..., 'analysis', 'trt_sserr', ...).
            ndraw = 8; NSUB = 3; NOCC = 2; NTRL = 5;
            ids = {'s1'; 's2'; 's3'};
            id2 = (1:3)';

            gro = zeros(ndraw, 2); popsd = zeros(ndraw, 1); indsd = zeros(ndraw, NSUB);
            sOcc = zeros(ndraw, 1); sTrl = zeros(ndraw, 1); sTID = zeros(ndraw, 1);
            sOID = zeros(ndraw, 1); sTO = zeros(ndraw, 1); mu = zeros(ndraw, 1);
            baseIndsd = [-0.3 0.0 0.3];  % per-subject log-residual offsets
            for d = 1:ndraw
                dr = (d - 1) * 0.005;
                gro(d, :) = [2.0 + dr, 0.4 + dr];   % [person-loc SD, person-resid log-SD spread]
                popsd(d) = 0.2 + dr;                % population log-residual
                indsd(d, :) = baseIndsd + dr;       % per-subject log-residual
                sOcc(d) = 0.5 + dr; sTrl(d) = 0.5 + dr; sTID(d) = 0.4 + dr;
                sOID(d) = 0.3 + dr; sTO(d) = 0.2 + dr; mu(d) = 5.0 + dr;
            end

            rowIds = {}; rowMeas = []; rowTime = {};
            for p = 1:NSUB
                for o = 1:NOCC
                    for t = 1:NTRL
                        rowIds{end+1, 1} = ids{p}; %#ok<AGROW>
                        rowTime{end+1, 1} = sprintf('t%d', o); %#ok<AGROW>
                        rowMeas(end+1, 1) = 5 + 0.1 * t + 0.2 * p; %#ok<AGROW>
                    end
                end
            end

            rel = struct();
            rel.analysis = 'trt_sserrvar';
            rel.groups = 'none';
            rel.events = 'none';
            rel.time = {'t1'; 't2'};
            rel.data = table(rowIds, rowMeas, rowTime, ...
                'VariableNames', {'id', 'meas', 'time'});
            rel.out = struct();
            rel.out.labels = {'none'};
            rel.out.mu = mu;
            rel.out.gro_sds(:, 1) = {num2cell(gro)};
            rel.out.pop_sdlog = popsd;
            rel.out.ind_sdlog(:, 1) = {num2cell(indsd)};
            rel.out.sig_occ = sOcc;
            rel.out.sig_trl = sTrl;
            rel.out.sig_trlxid = sTID;
            rel.out.sig_occxid = sOID;
            rel.out.sig_trlxocc = sTO;
            rel.out.id_matches(:, 1) = {table(ids, id2, 'VariableNames', {'id', 'id2'})};
            rel.out.conv.data = {{'name', 'n_eff', 'r_hat'; 'Intercept', 200, 1.00}};

            psyrat_data = struct('rel', rel);
        end

        function psyrat_data = makeICSSErrRelDataMultiGroup()
            %Subject-level single-occasion fixture (ic_sserrvar, analysis 6,
            %group-only / relsummary 'sing_sserr' case 2) with TWO groups that
            %deliberately reuse the SAME subject ids, each with a DIFFERENT
            %known trial count per (id,group). Regression fixture for S12: the
            %case-2 trial-count join previously grouped by 'id' alone, so a
            %repeated id's counts were summed across groups (s1: 5+9=14, s2:
            %7+11=18) and that pooled value was reused for BOTH groups. The
            %fix must instead recover each group's true per-id count exactly.
            ndraw = 6;
            ids = {'s1'; 's2'};
            id2 = (1:2)';
            groupNames = {'GroupA'; 'GroupB'};
            ngroup = 2;
            trueCounts = [5 7; 9 11]; % rows = group, cols = id (s1,s2)

            gro = cell(1, ngroup); indsd = cell(1, ngroup); idMatches = cell(1, ngroup);
            indBs = cell(1, ngroup);
            popsd = zeros(ndraw, ngroup); sTrl = zeros(ndraw, ngroup); mu = zeros(ndraw, ngroup);

            rowIds = {}; rowGroup = {}; rowMeas = [];
            for g = 1:ngroup
                groDraws = zeros(ndraw, 1);
                indDraws = zeros(ndraw, numel(ids));
                for d = 1:ndraw
                    drift = (d - 1) * 0.01;
                    groDraws(d) = 0.60 + drift + (g - 1) * 0.05;
                    popsd(d, g) = log(0.50 + drift + (g - 1) * 0.03);
                    indDraws(d, :) = log([0.45 0.55] + drift + (g - 1) * 0.02);
                    sTrl(d, g) = 0.30 + drift;
                    mu(d, g) = 5.0 + drift + (g - 1) * 0.2;
                end
                gro{g} = num2cell(groDraws);
                indsd{g} = num2cell(indDraws);
                indBs{g} = num2cell(zeros(ndraw, numel(ids)));
                idMatches{g} = table(ids, id2, 'VariableNames', {'id', 'id2'});

                for idi = 1:numel(ids)
                    n = trueCounts(g, idi);
                    for t = 1:n
                        rowIds{end+1, 1} = ids{idi}; %#ok<AGROW>
                        rowGroup{end+1, 1} = groupNames{g}; %#ok<AGROW>
                        rowMeas(end+1, 1) = 0.1 * t; %#ok<AGROW>
                    end
                end
            end

            rel = struct();
            rel.analysis = 'ic_sserrvar';
            rel.groups = groupNames;
            rel.events = 'none';
            rel.data = table(rowIds, rowGroup, rowMeas, ...
                'VariableNames', {'id', 'group', 'meas'});
            rel.out = struct();
            rel.out.labels = groupNames;
            rel.out.mu = mu;
            rel.out.ind_bs = indBs;
            rel.out.gro_sds = gro;
            rel.out.pop_sdlog = popsd;
            rel.out.ind_sdlog = indsd;
            rel.out.sig_trl = sTrl;
            rel.out.id_matches = idMatches;
            rel.out.conv.data = repmat({{'name', 'n_eff', 'Rhat'; 'Intercept', 200, 1.00}}, 1, ngroup);

            psyrat_data = struct('rel', rel);
            psyrat_data.expectedTrls = struct('GroupA', [5 7], 'GroupB', [9 11]);
        end

        function psyrat_data = makeICSSErrRelDataMultiGroupEvent()
            %Subject-level single-occasion fixture (ic_sserrvar, analysis 6,
            %group+event / relsummary 'sing_sserr' case 4) with TWO groups x
            %TWO events, deliberately reusing the SAME subject ids in every
            %(group,event) cell with a DIFFERENT known trial count each.
            %Regression fixture for S12: the case-4 trial-count join
            %previously filtered by event only, so a repeated id's counts
            %were summed across groups within each event.
            ndraw = 6;
            ids = {'s1'; 's2'};
            id2 = (1:2)';
            groupNames = {'GroupA'; 'GroupB'};
            eventNames = {'cor'; 'err'};
            ngroup = 2; nevent = 2;
            %true per-(group,event,id) counts, deliberately all distinct
            trueCounts = struct( ...
                'GroupA', struct('cor', [5 7], 'err', [4 6]), ...
                'GroupB', struct('cor', [9 11], 'err', [8 10]));

            ncol = ngroup * nevent;
            gro = cell(1, ncol); indsd = cell(1, ncol); idMatches = cell(1, ncol);
            indBs = cell(1, ncol);
            popsd = zeros(ndraw, ncol); sTrl = zeros(ndraw, ncol); mu = zeros(ndraw, ncol);
            labels = cell(ncol, 1);

            rowIds = {}; rowGroup = {}; rowEvent = {}; rowMeas = [];
            col = 0;
            for g = 1:ngroup
                for e = 1:nevent
                    col = col + 1;
                    labels{col} = strcat(groupNames{g}, '_;_', eventNames{e});

                    groDraws = zeros(ndraw, 1);
                    indDraws = zeros(ndraw, numel(ids));
                    for d = 1:ndraw
                        drift = (d - 1) * 0.01;
                        groDraws(d) = 0.60 + drift + (g - 1) * 0.05 + (e - 1) * 0.02;
                        popsd(d, col) = log(0.50 + drift + (g - 1) * 0.03);
                        indDraws(d, :) = log([0.45 0.55] + drift + (g - 1) * 0.02);
                        sTrl(d, col) = 0.30 + drift;
                        mu(d, col) = 5.0 + drift + (g - 1) * 0.2 + (e - 1) * 0.1;
                    end
                    gro{col} = num2cell(groDraws);
                    indsd{col} = num2cell(indDraws);
                    indBs{col} = num2cell(zeros(ndraw, numel(ids)));
                    idMatches{col} = table(ids, id2, 'VariableNames', {'id', 'id2'});

                    counts = trueCounts.(groupNames{g}).(eventNames{e});
                    for idi = 1:numel(ids)
                        n = counts(idi);
                        for t = 1:n
                            rowIds{end+1, 1} = ids{idi}; %#ok<AGROW>
                            rowGroup{end+1, 1} = groupNames{g}; %#ok<AGROW>
                            rowEvent{end+1, 1} = eventNames{e}; %#ok<AGROW>
                            rowMeas(end+1, 1) = 0.1 * t; %#ok<AGROW>
                        end
                    end
                end
            end

            rel = struct();
            rel.analysis = 'ic_sserrvar';
            rel.groups = groupNames;
            rel.events = eventNames;
            rel.data = table(rowIds, rowGroup, rowEvent, rowMeas, ...
                'VariableNames', {'id', 'group', 'event', 'meas'});
            rel.out = struct();
            rel.out.labels = labels;
            rel.out.mu = mu;
            rel.out.ind_bs = indBs;
            rel.out.gro_sds = gro;
            rel.out.pop_sdlog = popsd;
            rel.out.ind_sdlog = indsd;
            rel.out.sig_trl = sTrl;
            rel.out.id_matches = idMatches;
            rel.out.conv.data = repmat({{'name', 'n_eff', 'Rhat'; 'Intercept', 200, 1.00}}, 1, ncol);

            psyrat_data = struct('rel', rel);
            psyrat_data.expectedTrls = trueCounts;
        end

        function psyrat_data = makeTRTSSErrRelDataMultiGroup()
            %Subject-level crossed test-retest fixture (trt_sserrvar, analysis
            %25, relsummary 'trt_sserr') with TWO groups that deliberately
            %reuse the SAME subject ids, each with a DIFFERENT known
            %trials-per-occasion count per (id,group). Regression fixture for
            %S12: the trt_sserr trial-count join never included 'group' at
            %all, so repeated ids were pooled across groups exactly as in the
            %sing_sserr case-2 bug, feeding psyrat_ssrel_trt directly.
            ndraw = 6; NOCC = 2;
            ids = {'s1'; 's2'};
            id2 = (1:2)';
            groupNames = {'GroupA'; 'GroupB'};
            ngroup = 2;
            %true trials-PER-OCCASION per (group,id); NOCC occasions each
            trueTrlsPerOcc = struct('GroupA', [4 6], 'GroupB', [8 10]);

            gro = cell(1, ngroup); indsd = cell(1, ngroup); idMatches = cell(1, ngroup);
            popsd = zeros(ndraw, ngroup); sOcc = zeros(ndraw, ngroup);
            sTrl = zeros(ndraw, ngroup); sTID = zeros(ndraw, ngroup);
            sOID = zeros(ndraw, ngroup); sTO = zeros(ndraw, ngroup);

            rowIds = {}; rowGroup = {}; rowTime = {}; rowMeas = [];
            for g = 1:ngroup
                groDraws = zeros(ndraw, 2);
                indDraws = zeros(ndraw, numel(ids));
                for d = 1:ndraw
                    dr = (d - 1) * 0.005;
                    groDraws(d, :) = [2.0 + dr, 0.4 + dr] + (g - 1) * 0.1;
                    popsd(d, g) = 0.2 + dr;
                    indDraws(d, :) = [-0.3 0.3] + dr + (g - 1) * 0.05;
                    sOcc(d, g) = 0.5 + dr; sTrl(d, g) = 0.5 + dr; sTID(d, g) = 0.4 + dr;
                    sOID(d, g) = 0.3 + dr; sTO(d, g) = 0.2 + dr;
                end
                gro{g} = num2cell(groDraws);
                indsd{g} = num2cell(indDraws);
                idMatches{g} = table(ids, id2, 'VariableNames', {'id', 'id2'});

                counts = trueTrlsPerOcc.(groupNames{g});
                for idi = 1:numel(ids)
                    n = counts(idi);
                    for o = 1:NOCC
                        for t = 1:n
                            rowIds{end+1, 1} = ids{idi}; %#ok<AGROW>
                            rowGroup{end+1, 1} = groupNames{g}; %#ok<AGROW>
                            rowTime{end+1, 1} = sprintf('t%d', o); %#ok<AGROW>
                            rowMeas(end+1, 1) = 5 + 0.1 * t; %#ok<AGROW>
                        end
                    end
                end
            end

            rel = struct();
            rel.analysis = 'trt_sserrvar';
            rel.groups = groupNames;
            rel.events = 'none';
            rel.time = {'t1'; 't2'};
            rel.data = table(rowIds, rowGroup, rowMeas, rowTime, ...
                'VariableNames', {'id', 'group', 'meas', 'time'});
            rel.out = struct();
            rel.out.labels = groupNames;
            rel.out.gro_sds = gro;
            rel.out.pop_sdlog = popsd;
            rel.out.ind_sdlog = indsd;
            rel.out.sig_occ = sOcc;
            rel.out.sig_trl = sTrl;
            rel.out.sig_trlxid = sTID;
            rel.out.sig_occxid = sOID;
            rel.out.sig_trlxocc = sTO;
            rel.out.id_matches = idMatches;
            rel.out.conv.data = repmat({{'name', 'n_eff', 'r_hat'; 'Intercept', 200, 1.00}}, 1, ngroup);

            psyrat_data = struct('rel', rel);
            psyrat_data.expectedTrlsPerOcc = trueTrlsPerOcc;
        end

        function psyrat_data = makeICRelDataMultiGroupEvent()
            %Plain (non-subject-level) internal-consistency fixture with TWO
            %groups x TWO events (analysis 4, relsummary 'sing' case 4).
            %Correctness guard for the case-4 label contract. The plain-IC
            %group+event estimation path (psyrat_computevarcomp case 4)
            %emits REL.out.labels as event_;_group (e.g. 'cor_;_GroupA';
            %see the docstring example 'Error_;_Controls'), so the 'sing'
            %case-4 parse must read lblstr(1) as the EVENT and lblstr(2) as
            %the GROUP. This fixture mirrors that event_;_group order; the
            %test asserts each (group,event) cell is labelled correctly and
            %fails if the parse is swapped to group-first. REL.data is a
            %1 x nevent CELL array (one cell per event, each containing
            %every group's rows for that event) -- the eventarray.data
            %layout the 'ic' group+event estimation case (4) produces.
            ids = {'s1'; 's2'; 's3'};
            groupNames = {'GroupA'; 'GroupB'};
            eventNames = {'cor'; 'err'};
            ngroup = 2; nevent = 2;
            ncol = ngroup * nevent;

            labels = cell(ncol, 1);
            mu = zeros(4, ncol); sig_u = zeros(4, ncol);
            sig_e = zeros(4, ncol); sig_trl = zeros(4, ncol);
            col = 0;
            for g = 1:ngroup
                for e = 1:nevent
                    col = col + 1;
                    labels{col} = strcat(eventNames{e}, '_;_', groupNames{g});
                    mu(:, col) = 0.1 + 0.01 * (1:4)' + 0.02 * (g - 1) + 0.01 * (e - 1);
                    sig_u(:, col) = 0.45 + 0.01 * (1:4)' + 0.03 * (g - 1);
                    sig_e(:, col) = 0.95 + 0.01 * (1:4)' + 0.02 * (g - 1);
                    sig_trl(:, col) = 0.20 + 0.01 * (1:4)';
                end
            end

            dataByEvent = cell(1, nevent);
            for e = 1:nevent
                rowIds = {}; rowGroup = {}; rowEvent = {}; rowMeas = [];
                for g = 1:ngroup
                    for idi = 1:numel(ids)
                        for t = 1:6
                            rowIds{end+1, 1} = ids{idi}; %#ok<AGROW>
                            rowGroup{end+1, 1} = groupNames{g}; %#ok<AGROW>
                            rowEvent{end+1, 1} = eventNames{e}; %#ok<AGROW>
                            rowMeas(end+1, 1) = 0.1 * t; %#ok<AGROW>
                        end
                    end
                end
                dataByEvent{e} = table(rowIds, rowGroup, rowEvent, rowMeas, ...
                    'VariableNames', {'id', 'group', 'event', 'meas'});
            end

            rel = struct();
            rel.analysis = 'ic';
            rel.groups = groupNames;
            rel.events = eventNames;
            rel.data = dataByEvent;
            rel.out = struct();
            rel.out.labels = labels;
            rel.out.mu = mu;
            rel.out.sig_u = sig_u;
            rel.out.sig_e = sig_e;
            rel.out.sig_trl = sig_trl;

            psyrat_data = struct('rel', rel);
        end

        function psyrat_data = makeRawSingRelData()
            ids = [repmat({'s1'}, 4, 1); repmat({'s2'}, 5, 1); repmat({'s3'}, 6, 1)];
            meas = [1.2; 1.0; 1.1; 1.3; 0.7; 0.9; 1.0; 0.8; 1.1; 1.4; 1.5; 1.3; 1.2; 1.4; 1.3];

            rel = struct();
            rel.analysis = 'ic';
            rel.groups = 'none';
            rel.events = 'none';
            rel.out = struct();
            rel.out.labels = {'measure'};
            rel.out.mu = [0.1; 0.2; 0.15; 0.05];
            rel.out.sig_u = [0.45; 0.50; 0.47; 0.52];
            rel.out.sig_e = [0.95; 1.00; 1.05; 0.98];
            rel.data = table(ids, meas, 'VariableNames', {'id', 'meas'});

            psyrat_data = struct('rel', rel);
        end

        function psyrat_data = makeSingleEventViewData()
            psyrat_data = PsyRATTestDataFactory.makeRawSingRelData();
            psyrat_data.rel.analysis = 'ic';
            psyrat_data.rel.events = {'Reward'};
            psyrat_data.rel.groups = 'none';
            psyrat_data.rel.filename = 'single_event_fixture.psyrat';
            psyrat_data.rel.nchains = 4;
            psyrat_data.rel.niter = 1000;
            psyrat_data.proc = struct('measheader', 'Score');
            psyrat_data.ver = 'test';
        end

        function psyrat_data = makeDoDRelData()
            ndraw = 6;
            nevent = 4;
            ngroup = 2;

            eventNames = {'ERP_A'; 'ERP_B'; 'ERP_C'; 'ERP_D'};
            groupNames = {'Control'; 'Clinical'};

            idVarcov = cell(1, ngroup);
            trlVarcov = cell(1, ngroup);
            bSigma = cell(1, ngroup);
            bCell = cell(1, ngroup);
            convData = cell(1, ngroup);

            % Per-cell mean intercepts. Deliberately NOT centered on zero and
            % spread across both signs, in the microvolt range an ERP DoD design
            % actually produces: a fixture centered at 0 cannot distinguish an
            % extracted mean from the zeros placeholder the DoD -> single-event
            % view used to substitute, so the criterion coefficient's (mu - cut)^2
            % offset would test the same either way.
            baseCellMean = [-5.20 -3.10 2.40 6.70];

            baseId = [0.62 0.08 0.05 0.04; ...
                0.08 0.58 0.06 0.03; ...
                0.05 0.06 0.54 0.02; ...
                0.04 0.03 0.02 0.56];
            baseTrl = [0.95 0.10 0.06 0.05; ...
                0.10 0.90 0.07 0.04; ...
                0.06 0.07 0.88 0.03; ...
                0.05 0.04 0.03 0.92];

            for g = 1:ngroup
                idDraws = zeros(ndraw, nevent, nevent);
                trlDraws = zeros(ndraw, nevent, nevent);
                sigmaDraws = zeros(ndraw, nevent);
                cellMeanDraws = zeros(ndraw, nevent);
                for d = 1:ndraw
                    drift = (d - 1) * 0.005;
                    gshift = (g - 1) * 0.02;
                    idDraws(d, :, :) = baseId + eye(nevent) .* (drift + gshift);
                    trlDraws(d, :, :) = baseTrl + eye(nevent) .* (drift + (g - 1) * 0.03);
                    sigmaDraws(d, :) = log([0.70 0.76 0.74 0.79] + drift + (g - 1) * 0.015);
                    cellMeanDraws(d, :) = baseCellMean + (d - 1) * 0.01 + (g - 1) * 0.05;
                end
                idVarcov{g} = idDraws;
                trlVarcov{g} = trlDraws;
                bSigma{g} = sigmaDraws;
                bCell{g} = cellMeanDraws;
                convData{g} = {'param', 'n_eff', 'Rhat'; 'mu', 200, 1.00};
            end

            rows = {};
            baseCounts = [10 12 9 11; 8 7 10 6];
            gidNames = {'s1', 's2', 's3'; 's4', 's5', 's6'};
            for g = 1:ngroup
                for idi = 1:3
                    for e = 1:nevent
                        ntrials = baseCounts(g, e) + mod(idi + e, 2);
                        for t = 1:ntrials
                            rows(end+1, :) = {gidNames{g, idi}, groupNames{g}, eventNames{e}, 0.1 * t}; %#ok<AGROW>
                        end
                    end
                end
            end

            rel = struct();
            rel.analysis = 'ic_dodiff';
            rel.groups = groupNames;
            rel.events = eventNames;
            rel.filename = 'dod_fixture.psyrat';
            rel.nchains = 4;
            rel.niter = 2000;
            rel.data = cell2table(rows, 'VariableNames', {'id', 'group', 'event', 'meas'});
            rel.dod_map = {'ERP_A', 'ERP_B', 'ERP_C', 'ERP_D'};
            rel.out = struct();
            rel.out.id_varcov = idVarcov;
            rel.out.trl_varcov = trlVarcov;
            rel.out.b_sigma = bSigma;
            rel.out.b_cell = bCell;
            rel.out.conv = struct('data', {convData});

            psyrat_data = struct();
            psyrat_data.rel = rel;
            psyrat_data.proc = struct('measheader', 'Score');
            psyrat_data.ver = 'test';
        end

        function psyrat_data = makeSplitsSummaryData()
            % Single-occasion nonparallel data-splits result (ic_splits, Rocha
            % Table 4). Mirrors the localSingleREL fixture in
            % TestSplitsNonparallel: one stratum, 3 persons x 4 splits with
            % varying items-per-split (n_i). Scalar component "draws" collapse
            % each [ll pt ul] to the point estimate, which is all the headless
            % report needs to exercise.
            id = repelem((1:3)', 4);
            meas = (1:12)';                                  % arbitrary; unused by the calc
            weight = [10;20;30;40; 15;25;35;45; 12;22;32;42]; % n_i variation

            rel = struct();
            rel.analysis = 'ic_splits';
            rel.splits = 3;                                  % nonparallel splits
            rel.filename = 'splits_fixture.psyrat';
            rel.data = table(id, meas, weight);
            rel.groups = 'none';
            rel.events = 'none';
            rel.out = struct();
            rel.out.sig_u        = 2.0;   % sigma_p
            rel.out.sig_trl      = 1.0;   % sigma_s (split main)
            rel.out.sig_splitxid = 1.5;   % sigma_ps (person x split)
            rel.out.sig_err      = 3.0;   % per-item residual
            rel.out.labels = 'measure';

            psyrat_data = struct();
            psyrat_data.rel = rel;
            psyrat_data.proc = struct('savename', 'splits_fixture.psyrat');
            psyrat_data.ver = 'test';
        end

        function psyrat_data = makeSplitsTRTSummaryData()
            % Multi-occasion nonparallel data-splits result (trt_splits, Rocha
            % Table 5). Mirrors the localTrtREL fixture in TestSplitsNonparallel:
            % one stratum, 3 persons x 2 occasions x 4 splits, n_i varying across
            % rows. Reports the three test-retest reltypes (CE/CS/CES) per
            % stratum.
            id = repelem((1:3)', 8);
            time = repmat(repelem((1:2)', 4), 3, 1);
            meas = (1:24)';
            weight = (10:33)';                               % n_i variation across rows

            rel = struct();
            rel.analysis = 'trt_splits';
            rel.splits = 3;                                  % nonparallel splits
            rel.filename = 'splits_trt_fixture.psyrat';
            rel.data = table(id, meas, weight, time);
            rel.groups = 'none';
            rel.events = 'none';
            rel.time = {'t1','t2'};
            rel.out = struct();
            rel.out.sig_id       = 2.0;
            rel.out.sig_occ      = 0.8;
            rel.out.sig_trl      = 1.0;
            rel.out.sig_trlxid   = 1.5;
            rel.out.sig_occxid   = 1.2;
            rel.out.sig_trlxocc  = 0.5;
            rel.out.sig_posxid   = 0.9;
            rel.out.sig_err      = 3.0;
            rel.out.labels = {'none'};

            psyrat_data = struct();
            psyrat_data.rel = rel;
            psyrat_data.proc = struct('savename', 'splits_trt_fixture.psyrat');
            psyrat_data.ver = 'test';
        end

        function psyrat_data = makeDynrelDiffSummaryData()
            % Group-level two-event dynamic difference reliability
            % (ic_diff_dynrel, Rast & Clayson; one dimension, no per-participant
            % table). Mirrors the localFakeDiffREL fixture in
            % TestDiffDynrelReliability: constant posterior "draws" so the
            % credible interval collapses to the point estimate -- all the
            % headless report needs to flatten the surface / variance components.
            nd = 50;
            idvc = zeros(nd,2,2); trvc = zeros(nd,2,2);
            for k = 1:nd
                idvc(k,1,1) = 9; idvc(k,2,2) = 7; idvc(k,1,2) = 2; idvc(k,2,1) = 2;
                trvc(k,1,1) = 1.5; trvc(k,2,2) = 1.0; trvc(k,1,2) = 0.3; trvc(k,2,1) = 0.3;
            end
            o = struct();
            o.id_varcov   = {num2cell(idvc)};
            o.trl_varcov  = {num2cell(trvc)};
            o.b_sigma     = {num2cell(repmat([log(3) log(2.5)], nd, 1))};
            o.b_sigma_dim = {num2cell(repmat([0.2 -0.1], nd, 1))};
            o.b_dim       = {num2cell(repmat([1 0.5], nd, 1))};
            o.dimz        = {linspace(-1.5, 1.5, 20)'};
            o.ntrials     = 20;
            o.labels      = {'none'};

            rel = struct();
            rel.analysis = 'ic_diff_dynrel';
            rel.ndim = 1; rel.dim_names = {'dim1'};
            rel.filename = 'dynrel_diff_fixture.psyrat';
            rel.out = o;

            psyrat_data = struct();
            psyrat_data.rel = rel;
            psyrat_data.proc = struct('savename', 'dynrel_diff_fixture.psyrat');
            psyrat_data.ver = 'test';
        end

        function psyrat_data = makeDynrelSubjSummaryData()
            % Subject-level two-event dynamic difference reliability
            % (ic_diff_dynrel_sserrvar, Rast & Clayson; one dimension). Mirrors
            % the localFakeSubjREL fixture in TestSubjDiffDynrelReliability:
            % 3 subjects, identity person correlation, constant draws. This
            % variant DOES produce a per-participant (ssrel) table.
            nd = 30; nsub = 3;
            idvc1 = zeros(1,2,2);
            idvc1(1,1,1) = 9; idvc1(1,2,2) = 7; idvc1(1,1,2) = 2; idvc1(1,2,1) = 2;
            trvc1 = zeros(1,2,2);
            trvc1(1,1,1) = 1.5; trvc1(1,2,2) = 1.0; trvc1(1,1,2) = 0.3; trvc1(1,2,1) = 0.3;
            idvc = repmat(idvc1, nd, 1); trvc = repmat(trvc1, nd, 1);

            z = [-1; 0; 1];                         % per-subject standardized values
            trls1 = [20; 18; 22]; trls2 = [20; 16; 24];
            erss1 = [log(3); log(3.2); log(2.8)];   % per-subject z=0 log-residuals
            erss2 = [log(2.5); log(2.6); log(2.4)];
            sd = [2.5 2.0 0.4 0.3];                 % person joint-block SDs (loc1 loc2 sc1 sc2)
            Lid = zeros(nd, 4, 4);
            for k = 1:nd; Lid(k, :, :) = eye(4); end % identity correlation cholesky

            o = struct();
            o.id_varcov   = {num2cell(idvc)};
            o.trl_varcov  = {num2cell(trvc)};
            o.b_sigma     = {num2cell(repmat([log(3) log(2.5)], nd, 1))};
            o.b_sigma_dim = {num2cell(repmat([0.2 -0.1], nd, 1))};
            o.b_dim       = {num2cell(repmat([1 0.5], nd, 1))};
            o.er_var_ss1  = {num2cell(repmat(erss1', nd, 1))};
            o.er_var_ss2  = {num2cell(repmat(erss2', nd, 1))};
            o.sd_id       = {num2cell(repmat(sd, nd, 1))};
            o.L_id        = {num2cell(Lid)};
            id = (101:100+nsub)';
            ssinfo = table(id, (1:nsub)', z, trls1, trls2, round((trls1+trls2)/2), ...
                'VariableNames', {'id','id2','z1','trls1','trls2','trls'});
            o.ssinfo  = {ssinfo};
            o.dimz    = {z};
            o.ntrials = 20;
            o.labels  = {'none'};

            rel = struct();
            rel.analysis = 'ic_diff_dynrel_sserrvar';
            rel.ndim = 1; rel.dim_names = {'dim1'};
            rel.filename = 'dynrel_subj_fixture.psyrat';
            rel.out = o;

            psyrat_data = struct();
            psyrat_data.rel = rel;
            psyrat_data.proc = struct('savename', 'dynrel_subj_fixture.psyrat');
            psyrat_data.ver = 'test';
        end

        function psyrat_data = makeDynrelTrtSummaryData()
            % Trial + occasion two-facet dynamic reliability (ic_dynrel_trt,
            % Rast & Clayson; one dimension). Mirrors the localFabricateRel
            % fixture in TestDynrelTrtReliability: single (constant) draw, the
            % storage layout case 14 produces. Exercises the reltype/nocc
            % passthrough (the two-facet variant consumes them).
            rho = 0.3; L = [1 0; rho sqrt(1 - rho^2)];
            choldraw = zeros(1, 2, 2); choldraw(1, :, :) = L;

            o = struct();
            o.labels       = {'none'};
            o.dimz         = {(-1.5:1:1.5)'};       % 4 subject-level z values
            o.ntrials      = 8;                     % trial n'
            o.nocc         = 2;                     % observed occasion count
            o.gro_sds      = {num2cell([4, 0.5])};  % (sigma_p, sigma_delta_p)
            o.pop_sdlog    = log(3);                % Intercept_sigma
            o.sig_occ      = 1.5;
            o.sig_trl      = 2;
            o.sig_trlxid   = 1;
            o.sig_occxid   = 0.8;
            o.sig_trlxocc  = 0.5;
            o.b            = {num2cell(0.5)};       % mean dimension slope
            o.b_sigma      = {num2cell(0.25)};      % scale dimension slope
            o.chol_corrmat = {num2cell(choldraw)};

            rel = struct();
            rel.analysis = 'ic_dynrel_trt';
            rel.ndim = 1; rel.dim_names = {'dim1'};
            rel.filename = 'dynrel_trt_fixture.psyrat';
            rel.out = o;

            psyrat_data = struct();
            psyrat_data.rel = rel;
            psyrat_data.proc = struct('savename', 'dynrel_trt_fixture.psyrat');
            psyrat_data.ver = 'test';
        end

        function psyrat_data = makeDynrelTrt2DSummaryData()
            % Two-dimension variant of makeDynrelTrtSummaryData, to exercise the
            % ngrid^2 surface flatten (the 2D branch of local_dynrel_tables). The
            % design row is X(z) = [z1 z2 z1*z2], so the slopes are length 3 and
            % dimz is N x 2 (psyrat_dynrel_summary's 2D path reads dimz(:,1:2)).
            rho = 0.3; L = [1 0; rho sqrt(1 - rho^2)];
            choldraw = zeros(1, 2, 2); choldraw(1, :, :) = L;
            zcol = (-1.5:1:1.5)';

            o = struct();
            o.labels       = {'none'};
            o.dimz         = {[zcol zcol]};         % N x 2 standardized values
            o.ntrials      = 8;
            o.nocc         = 2;
            o.gro_sds      = {num2cell([4, 0.5])};
            o.pop_sdlog    = log(3);
            o.sig_occ      = 1.5;
            o.sig_trl      = 2;
            o.sig_trlxid   = 1;
            o.sig_occxid   = 0.8;
            o.sig_trlxocc  = 0.5;
            o.b            = {num2cell([0.5 0.25 0])};   % [z1 z2 z1*z2] mean slopes
            o.b_sigma      = {num2cell([0.25 0.125 0])}; % [z1 z2 z1*z2] scale slopes
            o.chol_corrmat = {num2cell(choldraw)};

            rel = struct();
            rel.analysis = 'ic_dynrel_trt';
            rel.ndim = 2; rel.dim_names = {'dim1','dim2'};
            rel.filename = 'dynrel_trt2d_fixture.psyrat';
            rel.out = o;

            psyrat_data = struct();
            psyrat_data.rel = rel;
            psyrat_data.proc = struct('savename', 'dynrel_trt2d_fixture.psyrat');
            psyrat_data.ver = 'test';
        end

        function psyrat_data = makeDynrelSubjNondiffSummaryData()
            % One-facet NON-difference subject-level dynamic reliability
            % (ic_dynrel_sserrvar, analysis 26): the same one-facet location-scale
            % layout case 11 produces, PLUS the per-subject residual (ind_sdlog)
            % and the per-participant lookup (ssinfo). 3 subjects, constant draws
            % (nd rows) so the point estimates are exact but psyrat_ssrel's
            % per-subject mean/quantile still reduce column-wise. Produces the
            % typical-person surface AND a per-participant (ssrel) table.
            nsub = 3; nd = 50;
            rho = 0.3; L = [1 0; rho sqrt(1 - rho^2)];
            choldraw = zeros(1, 2, 2); choldraw(1, :, :) = L;

            z     = [0.5; 0.0; -1.0];      % per-subject standardized values
            trls  = [10; 15; 20];          % per-subject trial counts (n')
            indsd = [-0.2, 0.0, 0.3];      % per-subject log-residual RE (1 x NSUB)
            id    = (101:100+nsub)';

            o = struct();
            o.mu           = zeros(nd,1);
            o.ind_bs       = {num2cell(zeros(nd,nsub))};
            o.gro_sds      = {num2cell(repmat([3, 0.5], nd, 1))}; % (sigma_p, sigma_delta_p)
            o.pop_sdlog    = repmat(log(2), nd, 1);               % Intercept_sigma
            o.ind_sdlog    = {num2cell(repmat(indsd, nd, 1))};    % nd x NSUB
            o.sig_trl      = repmat(1.5, nd, 1);                  % trial main effect (sigma_i)
            o.b            = {num2cell(repmat(0.5, nd, 1))};      % mean dimension slope
            o.b_sigma      = {num2cell(repmat(0.1, nd, 1))};      % scale dimension slope
            o.chol_corrmat = {num2cell(repmat(choldraw, nd, 1))};
            o.dimz         = {z};
            o.ntrials      = median(trls);
            o.ssinfo       = {table(id, (1:nsub)', z, trls, ...
                'VariableNames', {'id','id2','z1','trls'})};
            o.labels       = {'none'};

            rel = struct();
            rel.analysis = 'ic_dynrel_sserrvar';
            rel.ndim = 1; rel.dim_names = {'dim1'};
            rel.filename = 'dynrel_subj_nondiff_fixture.psyrat';
            rel.out = o;

            psyrat_data = struct();
            psyrat_data.rel = rel;
            psyrat_data.proc = struct('savename', 'dynrel_subj_nondiff_fixture.psyrat');
            psyrat_data.ver = 'test';
        end

        function psyrat_data = makeDynrelSubjNondiffTrtSummaryData()
            % Two-facet (trial + occasion) NON-difference subject-level dynamic
            % reliability (ic_dynrel_sserrvar_trt, analysis 27): the ic_dynrel_trt
            % layout PLUS the per-subject residual (ind_sdlog) and the
            % per-participant lookup (ssinfo with trials-per-occasion + occasion
            % count). Constant draws (nd rows). Produces the two-facet surface AND
            % a per-participant (ssrel) table.
            nsub = 3; nd = 50;
            rho = 0.3; L = [1 0; rho sqrt(1 - rho^2)];
            choldraw = zeros(1, 2, 2); choldraw(1, :, :) = L;

            z     = [0.5; 0.0; -1.0];
            trls  = [8; 10; 6];            % per-subject trials-per-occasion (n')
            occs  = [2; 2; 2];             % per-subject occasion counts
            indsd = [-0.2, 0.0, 0.3];
            id    = (201:200+nsub)';

            o = struct();
            o.labels       = {'none'};
            o.dimz         = {z};
            o.ntrials      = 8;
            o.nocc         = 2;
            o.gro_sds      = {num2cell(repmat([4, 0.5], nd, 1))};
            o.pop_sdlog    = repmat(log(3), nd, 1);
            o.ind_sdlog    = {num2cell(repmat(indsd, nd, 1))};    % nd x NSUB
            o.sig_occ      = repmat(1.5, nd, 1);
            o.sig_trl      = repmat(2, nd, 1);
            o.sig_trlxid   = repmat(1, nd, 1);
            o.sig_occxid   = repmat(0.8, nd, 1);
            o.sig_trlxocc  = repmat(0.5, nd, 1);
            o.b            = {num2cell(repmat(0.5, nd, 1))};      % mean dimension slope
            o.b_sigma      = {num2cell(repmat(0.25, nd, 1))};     % scale dimension slope
            o.chol_corrmat = {num2cell(repmat(choldraw, nd, 1))};
            o.ssinfo       = {table(id, (1:nsub)', z, trls, occs, ...
                'VariableNames', {'id','id2','z1','trls','occs'})};

            rel = struct();
            rel.analysis = 'ic_dynrel_sserrvar_trt';
            rel.ndim = 1; rel.dim_names = {'dim1'};
            rel.filename = 'dynrel_subj_nondiff_trt_fixture.psyrat';
            rel.out = o;

            psyrat_data = struct();
            psyrat_data.rel = rel;
            psyrat_data.proc = struct('savename', 'dynrel_subj_nondiff_trt_fixture.psyrat');
            psyrat_data.ver = 'test';
        end

        function psyrat_data = makeGammaDynrelSubjNondiffSummaryData()
            % GAMMA one-facet NON-difference subject-level dynamic reliability
            % (ic_dynrel_sserrvar, analysis 26, family='gamma'). Same shape as
            % makeDynrelSubjNondiffSummaryData, but the gamma REL.out layout: there
            % is no log-residual submodel, so pop_sdlog/b_sigma are absent and the
            % per-participant residual is derived from each person's own log-mean
            % (ind_bs) and log-nu (ind_nu). pop_lognu/b_nu are gamma-specific.
            % Constant draws (nd rows) so point estimates are exact.
            nsub = 3; nd = 50;
            rho = 0.3; L = [1 0; rho sqrt(1 - rho^2)];
            choldraw = zeros(1, 2, 2); choldraw(1, :, :) = L;

            z     = [0.5; 0.0; -1.0];      % per-subject standardized values
            trls  = [10; 15; 20];          % per-subject trial counts (n')
            indbs = [ 0.25, 0.00, -0.30];  % per-subject log-mean effect u_p
            indnu = [-0.20, 0.10,  0.35];  % per-subject log-nu  effect v_p
            id    = (301:300+nsub)';

            o = struct();
            o.mu           = repmat(log(5), nd, 1);               % Intercept (alpha0)
            o.pop_lognu    = repmat(log(16), nd, 1);              % Intercept_nu
            % s_v is deliberately DISTINCT from every other SD in the fixture: if
            % it equalled sig_trl, a wiring bug that sourced sig_i from gro_sds(:,2)
            % (the neighbouring column) would produce identical numbers and the
            % summary-wiring test could not see it.
            o.gro_sds      = {num2cell(repmat([0.40, 0.27], nd, 1))}; % (s_p, s_v)
            o.sig_trl      = repmat(0.30, nd, 1);                 % trial log-mean SD (s_i)
            o.b            = {num2cell(repmat(0.20, nd, 1))};     % log-mean slope
            o.b_nu         = {num2cell(repmat(-0.25, nd, 1))};    % log-nu slope
            o.chol_corrmat = {num2cell(repmat(choldraw, nd, 1))};
            o.ind_bs       = {num2cell(repmat(indbs, nd, 1))};    % nd x NSUB
            o.ind_nu       = {num2cell(repmat(indnu, nd, 1))};    % nd x NSUB
            o.dimz         = {z};
            o.ntrials      = median(trls);
            o.ssinfo       = {table(id, (1:nsub)', z, trls, ...
                'VariableNames', {'id','id2','z1','trls'})};
            o.labels       = {'none'};

            rel = struct();
            rel.analysis = 'ic_dynrel_sserrvar';
            rel.family = 'gamma';
            rel.ndim = 1; rel.dim_names = {'dim1'};
            rel.filename = 'gamma_dynrel_subj_nondiff_fixture.psyrat';
            rel.out = o;

            psyrat_data = struct();
            psyrat_data.rel = rel;
            psyrat_data.proc = struct('savename', 'gamma_dynrel_subj_nondiff_fixture.psyrat');
            psyrat_data.ver = 'test';
        end

        function psyrat_data = makeGammaLsDynrelSubjNondiffSummaryData()
            % LOCATION-SCALE GAMMA one-facet NON-difference subject-level dynamic
            % reliability (ic_dynrel_sserrvar, analysis 26, family='gamma',
            % gammascale=2). The location-scale twin of
            % makeGammaDynrelSubjNondiffSummaryData.
            %
            % The REL.out layout is deliberately the one psyrat_computevarcomp
            % writes under gammascale=2, which is the GAUSSIAN dynrel's slot
            % names: pop_sdlog / b_sigma / ind_sdlog, with NO pop_lognu, b_nu or
            % ind_nu. Two consequences a test can lean on: routing must come from
            % rel.gammascale rather than the family, and any consumer that gates
            % on pop_lognu existing (psyrat_provenance_lines) must fall silent.
            %
            % Constant draws (nd rows) so point estimates are exact.
            nsub = 3; nd = 50;
            rho = 0.3; L = [1 0; rho sqrt(1 - rho^2)];
            choldraw = zeros(1, 2, 2); choldraw(1, :, :) = L;

            z     = [0.5; 0.0; -1.0];      % per-subject standardized values
            trls  = [10; 15; 20];          % per-subject trial counts (n')
            indbs = [ 0.25, 0.00, -0.30];  % per-subject log-mean effect u_p
            indsd = [-0.20, 0.10,  0.35];  % per-subject log residual-SD effect w_p
            id    = (401:400+nsub)';

            o = struct();
            o.mu           = repmat(log(5), nd, 1);               % Intercept (alpha0)
            o.pop_sdlog    = repmat(log(1.4), nd, 1);             % Intercept_sigma
            % s_sd is deliberately DISTINCT from every other SD in the fixture,
            % for the same reason its log-nu twin keeps s_v distinct: if it
            % equalled sig_trl, a wiring bug sourcing sig_i from the neighbouring
            % gro_sds column would produce identical numbers and go unseen.
            o.gro_sds      = {num2cell(repmat([0.40, 0.27], nd, 1))}; % (s_p, s_sd)
            o.sig_trl      = repmat(0.30, nd, 1);                 % trial log-mean SD (s_i)
            o.b            = {num2cell(repmat(0.20, nd, 1))};     % log-mean slope
            o.b_sigma      = {num2cell(repmat(-0.25, nd, 1))};    % log residual-SD slope
            o.chol_corrmat = {num2cell(repmat(choldraw, nd, 1))};
            % ind_bs is stored by the estimator for parity and as a diagnostic;
            % the location-scale per-participant calculator never reads it. It is
            % kept here so the fixture matches what a real run produces.
            o.ind_bs       = {num2cell(repmat(indbs, nd, 1))};    % nd x NSUB
            o.ind_sdlog    = {num2cell(repmat(indsd, nd, 1))};    % nd x NSUB
            o.dimz         = {z};
            o.ntrials      = median(trls);
            o.ssinfo       = {table(id, (1:nsub)', z, trls, ...
                'VariableNames', {'id','id2','z1','trls'})};
            o.labels       = {'none'};

            rel = struct();
            rel.analysis = 'ic_dynrel_sserrvar';
            rel.family = 'gamma';
            rel.gammascale = 2;
            rel.ndim = 1; rel.dim_names = {'dim1'};
            rel.filename = 'gamma_ls_dynrel_subj_nondiff_fixture.psyrat';
            rel.out = o;

            psyrat_data = struct();
            psyrat_data.rel = rel;
            psyrat_data.proc = struct('savename', 'gamma_ls_dynrel_subj_nondiff_fixture.psyrat');
            psyrat_data.ver = 'test';
        end

        function psyrat_data = makeGammaDynrelSubjNondiffTrtSummaryData()
            % GAMMA two-facet (trial + occasion) NON-difference subject-level
            % dynamic reliability (ic_dynrel_sserrvar_trt, analysis 27,
            % family='gamma'). The gamma twin of
            % makeDynrelSubjNondiffTrtSummaryData: pop_lognu/b_nu replace
            % pop_sdlog/b_sigma, and ind_bs/ind_nu replace ind_sdlog.
            nsub = 3; nd = 50;
            rho = 0.3; L = [1 0; rho sqrt(1 - rho^2)];
            choldraw = zeros(1, 2, 2); choldraw(1, :, :) = L;

            z     = [0.5; 0.0; -1.0];
            trls  = [8; 10; 6];            % per-subject trials-per-occasion (n')
            occs  = [2; 2; 2];
            indbs = [ 0.25, 0.00, -0.30];
            indnu = [-0.20, 0.10,  0.35];
            id    = (401:400+nsub)';

            o = struct();
            o.labels       = {'none'};
            o.dimz         = {z};
            o.ntrials      = 8;
            o.nocc         = 2;
            o.mu           = repmat(log(5), nd, 1);
            o.pop_lognu    = repmat(log(18), nd, 1);
            o.gro_sds      = {num2cell(repmat([0.40, 0.30], nd, 1))};  % (s_p, s_v)
            o.sig_occ      = repmat(0.22, nd, 1);
            o.sig_trl      = repmat(0.18, nd, 1);
            o.sig_trlxid   = repmat(0.12, nd, 1);
            o.sig_occxid   = repmat(0.14, nd, 1);
            o.sig_trlxocc  = repmat(0.09, nd, 1);
            o.b            = {num2cell(repmat(0.20, nd, 1))};
            o.b_nu         = {num2cell(repmat(-0.25, nd, 1))};
            o.chol_corrmat = {num2cell(repmat(choldraw, nd, 1))};
            o.ind_bs       = {num2cell(repmat(indbs, nd, 1))};
            o.ind_nu       = {num2cell(repmat(indnu, nd, 1))};
            o.ssinfo       = {table(id, (1:nsub)', z, trls, occs, ...
                'VariableNames', {'id','id2','z1','trls','occs'})};

            rel = struct();
            rel.analysis = 'ic_dynrel_sserrvar_trt';
            rel.family = 'gamma';
            rel.ndim = 1; rel.dim_names = {'dim1'};
            rel.filename = 'gamma_dynrel_subj_nondiff_trt_fixture.psyrat';
            rel.out = o;

            psyrat_data = struct();
            psyrat_data.rel = rel;
            psyrat_data.proc = struct('savename', 'gamma_dynrel_subj_nondiff_trt_fixture.psyrat');
            psyrat_data.ver = 'test';
        end
    end
end
