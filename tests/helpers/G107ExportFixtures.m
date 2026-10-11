classdef G107ExportFixtures
    %Fixtures, an independent raw-draw reference, an RFC 4180 reader and
    %recording dialog and fopen doubles for TestSubjectDifferenceExport
    %(G107, analysis 8's per-participant difference rows in both subject
    %exports; G109, the CSV Save callbacks' final destination, text escaping
    %and write-failure cleanup).
    %
    %INDEPENDENCE. Nothing here calls PsyRATAccuracyOracle or a production
    %coefficient kernel. diffReference recomputes every per-participant
    %difference coefficient, ICC, SEM and error variance from the raw
    %posterior draws stored on rel.out and from trial counts tallied directly
    %from rel.data, with its own quantile (refQuantile). It encodes the
    %toolbox's DECLARED convention (the Rocha et al., 2026, Table 6 shape,
    %each cross-event covariance divided by the harmonic mean of the
    %participant's two counts), so it re-derives that convention
    %independently; it is not external ground truth for the convention.
    %
    %SCHEMA, ONE PLACE. guiDiffLabel, directionText, guiDiffColumns and
    %headlessDiffColumns hold the difference-section label and the column
    %names the tests require; a rename edits only these four methods.
    %
    %DIALOG DOUBLES. installShadows puts recording uiputfile, questdlg and
    %uiwait shadows first on the path. uiputfile answers from a queue (a 0
    %entry or an empty queue is Cancel), questdlg consents or declines from a
    %mode, and uiwait never blocks, so a consent prompt cannot hang a
    %headless run. Every call is recorded (dialogState); that is how a test
    %tells whether a prompt was issued and which file it named. A double
    %answers the CALLBACK; it says nothing about what the native macOS panel
    %returns. installErrordlgRecorder (opt-in, G109) records every errordlg
    %call so a test can count dialogs; installIoFaults and installFopenWatch
    %are the I/O fault seams.

    methods (Static)

        %% ---- fixtures ----------------------------------------------------

        function pd = baseRel()
            %The shared ic_diff_sserrvar draws (PsyRATTestDataFactory: 5
            %draws, participants s1-s3, loss/gain counts 6/5, 7/6 and 8/7,
            %so trls1 differs from trls2 for everyone) plus the provenance
            %fields every export header reads.
            pd = PsyRATTestDataFactory.makeICDiffSSErrRelData();
            pd = G107ExportFixtures.addHeaderFields(pd, 'g107_fixture.psyrat');
        end

        function pd = addHeaderFields(pd, filename)
            %Deterministic header and provenance fields. The save callbacks
            %read ver, filename and nchains/niter; psyrat_provenance_lines
            %reads seed, engine, family, priors and the convergence verdict.
            pd.ver = 'test';
            pd.rel.filename = filename;
            pd.rel.nchains = 1;
            pd.rel.niter = 5;
            pd.rel.seed = 12345;
            pd.rel.engine = 'cmdstan';
            pd.rel.family = 'gaussian';
            pd.rel.priors = psyrat_default_priors;
            pd.rel.out.conv.converged = 1;
            pd.rel.out.conv.ndivergent = 0;
        end

        function pd = groupedRel()
            %Two groups that share participant IDs with unequal sizes (3 and
            %2), input rows reversed, a PERMUTED id2 mapping in the second
            %group (s2 -> column 1, s1 -> column 2), different trial counts
            %for s2 (two first-event rows dropped) and different posterior
            %draws everywhere in the second group. A group mix-up, an id2
            %mis-mapping or a count mix-up therefore moves a value.
            pd = G107ExportFixtures.baseRel();
            rel = pd.rel;
            ev = cellstr(string(rel.events(:)));
            ndraw = size(rel.out.er_var_ss1{1}, 1);
            idB = zeros(ndraw, 2, 2);
            trB = zeros(ndraw, 2, 2);
            erB = zeros(ndraw, 2, 2);
            bsB = zeros(ndraw, 2);
            e1B = zeros(ndraw, 2);
            e2B = zeros(ndraw, 2);
            wcB = zeros(ndraw, 2);
            for d = 1:ndraw
                drift = (d - 1) * 0.02;
                idB(d, :, :) = [0.80 + drift, 0.10; 0.10, 0.70 + drift];
                trB(d, :, :) = [0.20 + drift, 0.03; 0.03, 0.25 + drift];
                erB(d, :, :) = [0.36 + drift, 0.05; 0.05, 0.44 + drift];
                bsB(d, :) = log([0.60 + drift, 0.66 + drift]);
                %columns follow the permuted id2 below: column 1 is s2,
                %column 2 is s1
                e1B(d, :) = log([0.90 + drift, 0.55 + drift]);
                e2B(d, :) = log([0.70 + drift, 0.50 + drift]);
                wcB(d, :) = [0.05, 0.02] + drift;
            end
            dataA = rel.data;
            dataA.group = repmat({'GrpAlpha'}, height(dataA), 1);
            dataB = rel.data(~strcmp(rel.data.id, 's3'), :);
            drop = find(strcmp(dataB.id, 's2') & strcmp(dataB.event, ev{1}), 2);
            dataB(drop, :) = [];
            dataB.group = repmat({'GrpBeta'}, height(dataB), 1);
            rel.data = [flipud(dataA); flipud(dataB)];
            rel.groups = {'GrpAlpha'; 'GrpBeta'};
            rel.out.id_varcov = [rel.out.id_varcov, {num2cell(idB)}];
            rel.out.trl_varcov = [rel.out.trl_varcov, {num2cell(trB)}];
            rel.out.err_varcov = [rel.out.err_varcov, {num2cell(erB)}];
            rel.out.b_sigma = [rel.out.b_sigma, {num2cell(bsB)}];
            rel.out.er_var_ss1 = [rel.out.er_var_ss1, {e1B}];
            rel.out.er_var_ss2 = [rel.out.er_var_ss2, {e2B}];
            rel.out.wp_cov_ss = [rel.out.wp_cov_ss, {wcB}];
            rel.out.id_matches = [rel.out.id_matches, ...
                {table({'s2'; 's1'}, [1; 2], 'VariableNames', {'id', 'id2'})}];
            rel.out.elabels = [rel.out.elabels, rel.out.elabels];
            rel.out.glabels = {'GrpAlpha', 'GrpBeta'};
            rel.out.conv.data = [rel.out.conv.data, rel.out.conv.data];
            pd.rel = rel;
        end

        function pd = reversedEventRel()
            %The base draws with the EVENT ORDER reversed, as a fit whose
            %data listed 'gain' first would store them: rel.events and
            %rel.diff_names become {'gain'; 'loss'}, every event-indexed
            %block is permuted to match, and the gain rows come first in the
            %data. The difference coefficient is symmetric in the two
            %events, so every per-participant value must equal the base run
            %while trls1/trls2 and the contrast labels swap.
            pd = G107ExportFixtures.baseRel();
            rel = pd.rel;
            idv = cell2mat(rel.out.id_varcov{1});
            trv = cell2mat(rel.out.trl_varcov{1});
            erv = cell2mat(rel.out.err_varcov{1});
            bsv = cell2mat(rel.out.b_sigma{1});
            rel.out.id_varcov = {num2cell(idv(:, [2 1], [2 1]))};
            rel.out.trl_varcov = {num2cell(trv(:, [2 1], [2 1]))};
            rel.out.err_varcov = {num2cell(erv(:, [2 1], [2 1]))};
            rel.out.b_sigma = {num2cell(bsv(:, [2 1]))};
            first = rel.out.er_var_ss1{1};
            rel.out.er_var_ss1 = {rel.out.er_var_ss2{1}};
            rel.out.er_var_ss2 = {first};
            rel.events = flipud(cellstr(string(rel.events(:))));
            rel.diff_names = rel.events;
            rel.out.elabels = {rel.events};
            [~, order] = sort(~strcmp(rel.data.event, rel.events{1}));
            rel.data = rel.data(order, :);
            pd.rel = rel;
        end

        function [old, new] = specialIdMap()
            %Participant IDs that need RFC 4180 quoting: a comma, a double
            %quote, and an embedded line feed.
            old = {'s1'; 's2'; 's3'};
            new = {'subject,one'; 'say "hi"'; sprintf('line one\nline two')};
        end

        function pd = specialIdRel()
            %The base draws with the special IDs injected at the DATA level
            %(rel.data.id and the id_matches table), ahead of the summary, so
            %the summary tables carry them the way the pipeline would. The
            %file loader is not exercised; whether it admits such IDs is a
            %separate question.
            pd = G107ExportFixtures.baseRel();
            [old, new] = G107ExportFixtures.specialIdMap();
            t = pd.rel.out.id_matches{1};
            for k = 1:numel(old)
                pd.rel.data.id(strcmp(pd.rel.data.id, old{k})) = new(k);
                t.id(strcmp(t.id, old{k})) = new(k);
            end
            pd.rel.out.id_matches{1} = t;
            pd.rel.filename = 'g107_special_ids.psyrat';
        end

        function pd = specialGroupRel()
            %groupedRel with group names that need quoting (a comma and a
            %double quote), so the difference rows' Group field and the
            %section labels exercise the escaping too.
            pd = G107ExportFixtures.groupedRel();
            oldg = {'GrpAlpha'; 'GrpBeta'};
            newg = {'Site, North'; 'Site "South"'};
            for k = 1:numel(oldg)
                pd.rel.data.group(strcmp(pd.rel.data.group, oldg{k})) = newg(k);
            end
            pd.rel.groups = newg;
            pd.rel.out.glabels = newg';
        end

        function pd = lineBreakGroupRel()
            %groupedRel with one group name holding a line feed and the other
            %a carriage return. A hand-built struct: psyrat_preflight_validate
            %refuses id, group, event and time labels containing a comma, a
            %double quote or a line break on every new fit, so only an older
            %result or a built struct can carry such names. The export must
            %still keep each such label in one field.
            pd = G107ExportFixtures.groupedRel();
            oldg = {'GrpAlpha'; 'GrpBeta'};
            newg = {sprintf('Site\nNorth'); sprintf('Site\rSouth')};
            for k = 1:numel(oldg)
                pd.rel.data.group(strcmp(pd.rel.data.group, oldg{k})) = newg(k);
            end
            pd.rel.groups = newg;
            pd.rel.out.glabels = newg';
        end

        function pd = icDiffViewData()
            %A group-level ic_diff summary that psyrat_variancet and its View
            %All Standard Deviations child (psyrat_viewdiffcomp) can render:
            %the TestUitableWidthSlack recipe (per-event variance structs plus
            %the difference components, deterministic values).
            pd = PsyRATTestDataFactory.makeICDiffSummaryData();
            pd = G107ExportFixtures.addHeaderFields(pd, 'g107_view_all.psyrat');
            for g = 1:numel(pd.relsummary.group)
                for e = 1:numel(pd.relsummary.group(g).event)
                    base = 1 + 0.1 * e;
                    names = {'icc', 'betsd', 'trlsd', 'witsd', 'sem'};
                    for k = 1:numel(names)
                        v = base * k / 5;
                        pd.relsummary.group(g).event(e).(names{k}) = ...
                            struct('m', v, 'll', v * 0.8, 'ul', v * 1.2);
                    end
                end
                ds = pd.relsummary.group(g).diffscore;
                names = {'bt1_var', 'bt2_var', 'bp_cov', 'bt_cov', 'icc'};
                vals = [2.1 1.7 0.9 0.4 0.55];
                for k = 1:numel(names)
                    ds.([names{k} '_pt']) = vals(k);
                    ds.([names{k} '_ll']) = vals(k) * 0.8;
                    ds.([names{k} '_ul']) = vals(k) * 1.2;
                end
                ds.wp_cov = 0.3;
                pd.relsummary.group(g).diffscore = ds;
            end
        end

        function pd = criterionSingData()
            %The one-facet criterion fixture the TestExportCallbacks
            %criterion tests use.
            pd = PsyRATTestDataFactory.makeICSummaryDataWithTrial();
            pd = G107ExportFixtures.addHeaderFields(pd, 'g107_criterion.psyrat');
        end

        function prefs = criterionPrefs(cutoff)
            %Criterion table only (no plot), as in TestExportCallbacks.
            prefs = psyrat_defaults();
            prefs.view.plotcriterion = 0;
            prefs.view.tablecriterion = 1;
            prefs.view.criterioncutoff = cutoff;
            prefs.view.reltype = 1;
            prefs.view.gcoeff = 2;
            prefs.view.depcentmeas = 1;
            prefs.view.relcentmeas = 1;
        end

        function [pd, relerr] = summarize(pd, gc, ci, depcutoff, meascutoff)
            %The real analysis-8 summary (psyrat_relsummary 'sing_diff_sserr');
            %no fit.
            [pd, relerr] = psyrat_relsummary('psyrat_data', pd, ...
                'analysis', 'sing_diff_sserr', 'depcutoff', depcutoff, ...
                'meascutoff', meascutoff, 'depcentmeas', 1, ...
                'diffgcoeff', gc, 'CI', ci);
        end

        function cut = mixedCutoffValue(pd0, gc, ci, meas)
            %A cutoff midway between the two smallest selected-bound
            %difference values pooled over groups, derived from a clean run
            %of the SAME draws, so the lowest participant is excluded and the
            %rest are included. The per-participant values do not depend on
            %the cutoff, so re-summarizing at this cutoff moves only flags.
            clean = G107ExportFixtures.summarize(pd0, gc, ci, 0.001, meas);
            vals = [];
            for g = 1:numel(clean.relsummary.group)
                T = clean.relsummary.group(g).diffscore.ssrel_table;
                vals = [vals; G107ExportFixtures.boundValues(T, meas)]; %#ok<AGROW>
            end
            vals = unique(vals(isfinite(vals)));
            if numel(vals) < 2
                error('G107ExportFixtures:noMixedCutoff', ...
                    'The fixture needs two distinct difference values to stage a mixed cutoff.');
            end
            cut = (vals(1) + vals(2)) / 2;
        end

        function cut = mixedCutoff(testCase, pd0, gc, ci, meas)
            %mixedCutoffValue with the failure reported as a failed weld.
            cut = [];
            try
                cut = G107ExportFixtures.mixedCutoffValue(pd0, gc, ci, meas);
            catch ME
                testCase.assertFail(['weld: ' ME.message]);
            end
        end

        function v = boundValues(T, meas)
            %The bound the inclusion rule compares (1 lower limit, 2 point,
            %3 upper limit), mirroring the meascutoff codes.
            switch meas
                case 1
                    v = T.dep_ll;
                case 2
                    v = T.dep_pt;
                otherwise
                    v = T.dep_ul;
            end
        end

        %% ---- independent reference ---------------------------------------

        function R = diffReference(pd, gloc, gc, ci, depcutoff, meascutoff)
            %Per-participant difference-score reference for group GLOC,
            %recomputed from raw draws. For participant p with own counts n1
            %and n2 (tallied from rel.data within the group), per draw:
            %  universe  u   = s2_p1 + s2_p2 - 2*s_p12            (id_varcov)
            %  relative  rel = w1/n1 + w2/n2 - 2*c/h, w = exp(er_var_ss)^2,
            %                  c = wp_cov_ss, h = 2/(1/n1 + 1/n2)
            %  absolute  abs = rel + t1/n1 + t2/n2 - 2*t12/h       (trl_varcov)
            %  coefficient u/(u + abs) for dependability (gc = 1) or
            %              u/(u + rel) for generalizability (gc = 2)
            %  ICC (n' = 1, derived) u/(u + e1) with e1 the same error at
            %              n1 = n2 = 1 (w1 + w2 - 2c, plus the trial block
            %              for dependability)
            %  SEM sqrt(error), undefined (NaN) for a negative draw
            %  errvar      mean of the RELATIVE error of the averaged
            %              contrast at the participant's own counts
            %Summaries: posterior mean, and equal-tailed limits at
            %(1 - ci)/2 and 1 - (1 - ci)/2 from refQuantile. The flag applies
            %the cutoff to the bound meascutoff names. A participant without
            %rows in both events has no difference score and is skipped,
            %never fabricated.
            rel = pd.rel;
            gnames = G107ExportFixtures.groupLabels(pd);
            enames = G107ExportFixtures.eventNames(pd);
            data = rel.data;
            if ~isempty(gnames{gloc})
                data = data(strcmp(cellstr(string(data.group)), gnames{gloc}), :);
            end
            ids = cellstr(string(data.id));
            evs = cellstr(string(data.event));
            idm = rel.out.id_matches{gloc};
            bp = cell2mat(rel.out.id_varcov{gloc});
            bt = cell2mat(rel.out.trl_varcov{gloc});
            er1 = rel.out.er_var_ss1{gloc};
            er2 = rel.out.er_var_ss2{gloc};
            wc = rel.out.wp_cov_ss{gloc};
            u = bp(:, 1, 1) + bp(:, 2, 2) - 2 .* bp(:, 2, 1);
            t1 = bt(:, 1, 1);
            t2 = bt(:, 2, 2);
            t12 = bt(:, 2, 1);
            lo = (1 - ci) / 2;
            hi = 1 - lo;
            q = @(x, p) G107ExportFixtures.refQuantile(x, p);
            idc = cellstr(string(idm.id));
            n = numel(idc);
            keep = false(n, 1);
            out = nan(n, 12);
            for k = 1:n
                n1 = sum(strcmp(ids, idc{k}) & strcmp(evs, enames{1}));
                n2 = sum(strcmp(ids, idc{k}) & strcmp(evs, enames{2}));
                if n1 == 0 || n2 == 0
                    continue;
                end
                keep(k) = true;
                s = idm.id2(k);
                w1 = exp(er1(:, s)) .^ 2;
                w2 = exp(er2(:, s)) .^ 2;
                c = wc(:, s);
                h = 2 / (1 / n1 + 1 / n2);
                relerr = w1 ./ n1 + w2 ./ n2 - 2 .* c ./ h;
                abserr = relerr + t1 ./ n1 + t2 ./ n2 - 2 .* t12 ./ h;
                icc1 = w1 + w2 - 2 .* c;
                if gc == 1
                    err = abserr;
                    iccerr = icc1 + t1 + t2 - 2 .* t12;
                else
                    err = relerr;
                    iccerr = icc1;
                end
                coef = u ./ (u + err);
                icc = u ./ (u + iccerr);
                sem = nan(size(err));
                ok = err >= 0;
                sem(ok) = sqrt(err(ok));
                out(k, :) = [n1, n2, mean(coef), q(coef, lo), q(coef, hi), ...
                    mean(icc), q(icc, lo), q(icc, hi), ...
                    mean(sem), q(sem, lo), q(sem, hi), mean(relerr)];
            end
            R = array2table(out(keep, :), 'VariableNames', {'trls1', 'trls2', ...
                'pt', 'll', 'ul', 'icc_pt', 'icc_ll', 'icc_ul', ...
                'sem_pt', 'sem_ll', 'sem_ul', 'errvar'});
            R = addvars(R, idc(keep), 'Before', 1, 'NewVariableNames', 'id');
            switch meascutoff
                case 1
                    sel = R.ll;
                case 2
                    sel = R.pt;
                otherwise
                    sel = R.ul;
            end
            R.include = sel >= depcutoff;
        end

        function q = refQuantile(x, p)
            %MATLAB's documented default quantile, re-implemented rather than
            %called: NaNs are dropped, the k-th sorted draw sits at
            %cumulative probability (k - 0.5)/n, values between are linearly
            %interpolated, and probabilities outside [0.5/n, 1 - 0.5/n] take
            %the extreme draw. With 5 draws an 80% or 95% interval therefore
            %lands on the smallest and largest draw exactly, and a 50%
            %interval interpolates.
            x = sort(x(~isnan(x)));
            n = numel(x);
            if n == 0
                q = NaN;
                return;
            end
            pos = n * p + 0.5;
            if pos <= 1
                q = x(1);
            elseif pos >= n
                q = x(n);
            else
                k = floor(pos);
                q = x(k) + (pos - k) * (x(k + 1) - x(k));
            end
        end

        %% ---- labels and the export schema -------------------------------

        function g = groupLabels(pd)
            %Stratum labels: {''} for an ungrouped run (the headless
            %sscoeffs Group convention), else the group names.
            rel = pd.rel;
            if ischar(rel.groups) && strcmpi(rel.groups, 'none')
                g = {''};
            else
                g = cellstr(string(rel.groups(:)))';
            end
        end

        function e = eventNames(pd)
            %The fitted events in model order (rel.events == rel.diff_names).
            e = cellstr(string(pd.rel.events(:)))';
        end

        function labels = eventSectionLabels(pd)
            %The event-section labels psyrat_savesscoeffs writes, in its
            %group-major loop order ('Measurement', group, event, or
            %'group - event').
            rel = pd.rel;
            g = G107ExportFixtures.groupLabels(pd);
            if ischar(rel.events) && strcmpi(rel.events, 'none')
                e = {''};
            else
                e = cellstr(string(rel.events(:)))';
            end
            labels = {};
            for gi = 1:numel(g)
                for ei = 1:numel(e)
                    if numel(g) == 1 && numel(e) == 1
                        labels{end+1} = 'Measurement'; %#ok<AGROW>
                    elseif numel(g) > 1 && numel(e) == 1
                        labels{end+1} = g{gi}; %#ok<AGROW>
                    elseif numel(g) == 1 && numel(e) > 1
                        labels{end+1} = e{ei}; %#ok<AGROW>
                    else
                        labels{end+1} = [g{gi} ' - ' e{ei}]; %#ok<AGROW>
                    end
                end
            end
        end

        function labels = diffSectionLabels(pd)
            %One difference-section label per group, in group
            %order, appended after every event section.
            g = G107ExportFixtures.groupLabels(pd);
            e = G107ExportFixtures.eventNames(pd);
            labels = cellfun(@(x) G107ExportFixtures.guiDiffLabel(x, e{1}, e{2}), ...
                g, 'UniformOutput', false);
        end

        function label = guiDiffLabel(gname, eA, eB)
            %'Difference Score (loss minus gain)', prefixed with
            %'<group> - ' in grouped runs (the event-section convention).
            base = sprintf('Difference Score (%s)', ...
                G107ExportFixtures.directionText(eA, eB));
            if isempty(gname)
                label = base;
            else
                label = [gname ' - ' base];
            end
        end

        function s = directionText(eA, eB)
            %Fitted event A minus fitted event B (rel.events{1}
            %minus rel.events{2}, the model's column order). The coefficient
            %is symmetric in the two events; the label states the contrast.
            s = sprintf('%s minus %s', eA, eB);
        end

        function s = coefLabel(gc)
            %Difference coefficient label from diffgcoeff.
            if isequal(gc, 2)
                s = 'Generalizability';
            else
                s = 'Dependability';
            end
        end

        function s = eventCoefLabel(pd)
            %The event header's coefficient label, read the way
            %psyrat_savesscoeffs reads it (relsummary.gcoeff_name).
            s = 'Dependability';
            if isfield(pd.relsummary, 'gcoeff_name') && ...
                    strcmpi(pd.relsummary.gcoeff_name, 'gen')
                s = 'Generalizability';
            end
        end

        function h = guiEventHeader(coef)
            %The shipped event-section column row (psyrat_variancet.m,
            %psyrat_savesscoeffs); it must not change.
            h = {'ID', 'Participant Meets Reliability Cutoff', 'Trial Counts', ...
                coef, [coef ' Lower Limit'], [coef ' Upper Limit'], ...
                'ICC', 'ICC Lower Limit', 'ICC Upper Limit', ...
                'SEM', 'SEM Lower Limit', 'SEM Upper Limit', 'Error Variance'};
        end

        function [keys, names] = guiDiffColumns()
            %GUI difference-section columns (semantic key, header).
            %The difference error variance is deliberately NOT called
            %'Error Variance': the event column of that name is a
            %trial-level residual variance, and this one is the relative
            %error variance of the averaged contrast at the participant's
            %own counts.
            spec = { ...
                'id',            'ID'; ...
                'group',         'Group'; ...
                'scoretype',     'Score Type'; ...
                'eventA',        'Event A'; ...
                'eventB',        'Event B'; ...
                'direction',     'Direction'; ...
                'coefficient',   'Coefficient'; ...
                'include',       'Difference Score Meets Reliability Cutoff'; ...
                'cutoff',        'Reliability Cutoff'; ...
                'bound',         'Cutoff Bound'; ...
                'trls1',         'Trial Counts (Event A)'; ...
                'trls2',         'Trial Counts (Event B)'; ...
                'pt',            'Coefficient Posterior Mean'; ...
                'll',            'Coefficient Lower Limit'; ...
                'ul',            'Coefficient Upper Limit'; ...
                'intervalkind',  'Interval'; ...
                'intervalwidth', 'Interval Width'; ...
                'icc_pt',        'ICC (n''=1; derived)'; ...
                'icc_ll',        'ICC Lower Limit'; ...
                'icc_ul',        'ICC Upper Limit'; ...
                'sem_pt',        'SEM'; ...
                'sem_ll',        'SEM Lower Limit'; ...
                'sem_ul',        'SEM Upper Limit'; ...
                'errvar',        'Relative Error Variance of the Mean Difference'};
            keys = spec(:, 1)';
            names = spec(:, 2)';
        end

        function names = headlessDiffColumns()
            %tables.sscoeffs_diff variables, Group-first like
            %tables.sscoeffs. dep_* hold the selected difference coefficient
            %(the sscoeffs convention); IntervalWidth is the proportion
            %(relsummary.ciperc); diff_rel_errvar is the difference
            %table's ss_errvar under a name that cannot be mistaken for the
            %event table's ss_errvar.
            names = {'Group', 'ScoreType', 'EventA', 'EventB', 'Direction', ...
                'Coefficient', 'id', 'ind2include', 'Cutoff', 'CutoffBound', ...
                'trls1', 'trls2', 'dep_pt', 'dep_ll', 'dep_ul', ...
                'IntervalKind', 'IntervalWidth', 'icc_pt', 'icc_ll', 'icc_ul', ...
                'sem_pt', 'sem_ll', 'sem_ul', 'diff_rel_errvar'};
        end

        function names = headlessEventColumns()
            %The shipped tables.sscoeffs variables (TestHeadlessReport:222).
            names = {'Group', 'Event', 'id', 'ind2include', 'trls', ...
                'dep_pt', 'dep_ll', 'dep_ul', 'icc_pt', 'icc_ll', 'icc_ul', ...
                'sem_pt', 'sem_ll', 'sem_ul', 'ss_errvar'};
        end

        function T = flattenEventTables(pd)
            %The shipped tables.sscoeffs contract, rebuilt from the stored
            %event tables: every group x event table's exported columns with
            %Group and Event prepended, stacked in group-major order.
            cols = {'id', 'ind2include', 'trls', 'dep_pt', 'dep_ll', 'dep_ul', ...
                'icc_pt', 'icc_ll', 'icc_ul', 'sem_pt', 'sem_ll', 'sem_ul', 'ss_errvar'};
            rel = pd.rel;
            g = G107ExportFixtures.groupLabels(pd);
            if ischar(rel.events) && strcmpi(rel.events, 'none')
                e = {''};
            else
                e = cellstr(string(rel.events(:)))';
            end
            parts = {};
            for gi = 1:numel(g)
                for ei = 1:numel(e)
                    t = pd.relsummary.group(gi).event(ei).ssrel_table(:, cols);
                    t = addvars(t, repmat(e(ei), height(t), 1), ...
                        'Before', 1, 'NewVariableNames', 'Event');
                    t = addvars(t, repmat(g(gi), height(t), 1), ...
                        'Before', 1, 'NewVariableNames', 'Group');
                    parts{end+1} = t; %#ok<AGROW>
                end
            end
            T = vertcat(parts{:});
        end

        %% ---- CSV reading --------------------------------------------------

        function recs = parseCsv(txt)
            %Independent RFC 4180 reader (not the toolbox writer and not
            %readtable). A double-quoted field may hold commas, doubled
            %quotes and line breaks; an unquoted field ends at a comma or a
            %line break (LF, CRLF, or a bare CR, as Python's csv reader and
            %spreadsheets treat it). A quote inside an unquoted field is kept
            %literally, as Python's csv reader keeps it. An unterminated
            %quoted field is an error, never a silent truncation. Returns a
            %cell of records, each a 1xN cell of char fields.
            recs = {};
            fields = {};
            field = '';
            inq = false;
            n = numel(txt);
            i = 1;
            while i <= n
                c = txt(i);
                if inq
                    if c == '"'
                        if i < n && txt(i + 1) == '"'
                            field(end+1) = '"'; %#ok<AGROW>
                            i = i + 2;
                        else
                            inq = false;
                            i = i + 1;
                        end
                    else
                        field(end+1) = c; %#ok<AGROW>
                        i = i + 1;
                    end
                elseif c == '"' && isempty(field)
                    inq = true;
                    i = i + 1;
                elseif c == ','
                    fields{end+1} = field; %#ok<AGROW>
                    field = '';
                    i = i + 1;
                elseif c == char(13) && i < n && txt(i + 1) == char(10)
                    fields{end+1} = field; %#ok<AGROW>
                    recs{end+1} = fields; %#ok<AGROW>
                    fields = {};
                    field = '';
                    i = i + 2;
                elseif c == char(10) || c == char(13)
                    fields{end+1} = field; %#ok<AGROW>
                    recs{end+1} = fields; %#ok<AGROW>
                    fields = {};
                    field = '';
                    i = i + 1;
                else
                    field(end+1) = c; %#ok<AGROW>
                    i = i + 1;
                end
            end
            if inq
                error('G107ExportFixtures:unterminatedQuote', ...
                    'The CSV ends inside a quoted field.');
            end
            if ~isempty(field) || ~isempty(fields)
                fields{end+1} = field;
                recs{end+1} = fields;
            end
        end

        function S = guiSections(recs, labels)
            %Split parsed records into labeled sections. A section starts at
            %a single-field record equal to one of LABELS and runs to the
            %next such record or the end of the file. notes: single-field
            %records before the section's 'ID' column row; data: multi-field
            %records after it; trailing: single-field records after it.
            %Whitespace-only spacer records (the ' ' lines the header block
            %already uses) are dropped from notes and trailing, so a spacer
            %is never mistaken for the stray first half of a split ID.
            isLabel = false(1, numel(recs));
            for r = 1:numel(recs)
                isLabel(r) = numel(recs{r}) == 1 && any(strcmp(recs{r}{1}, labels));
            end
            starts = find(isLabel);
            S = struct('label', {}, 'index', {}, 'header', {}, 'notes', {}, ...
                'data', {}, 'trailing', {});
            for k = 1:numel(starts)
                stop = numel(recs);
                if k < numel(starts)
                    stop = starts(k + 1) - 1;
                end
                span = recs(starts(k) + 1:stop);
                h = find(cellfun(@(x) numel(x) > 1 && strcmp(x{1}, 'ID'), span), 1);
                if isempty(h)
                    header = {};
                    notes = G107ExportFixtures.singles(span);
                    data = {};
                    trailing = {};
                else
                    header = span{h};
                    notes = G107ExportFixtures.singles(span(1:h - 1));
                    after = span(h + 1:end);
                    multi = cellfun(@numel, after) > 1;
                    data = after(multi);
                    trailing = G107ExportFixtures.singles(after(~multi));
                end
                notes = notes(~cellfun(@(x) isempty(strtrim(x)), notes));
                trailing = trailing(~cellfun(@(x) isempty(strtrim(x)), trailing));
                S(end+1) = struct('label', {recs{starts(k)}{1}}, 'index', starts(k), ...
                    'header', {header}, 'notes', {notes}, 'data', {data}, ...
                    'trailing', {trailing}); %#ok<AGROW>
            end
        end

        function out = singles(recs)
            %The text of every single-field record, in order (1xN cell).
            out = {};
            for r = 1:numel(recs)
                if numel(recs{r}) == 1
                    out{end+1} = recs{r}{1}; %#ok<AGROW>
                end
            end
        end

        function v = parseNum(x)
            %A number from a CSV field ('NaN', '' and a trailing % allowed)
            %or an XLSX cell (double, logical or missing).
            if iscell(x)
                if isempty(x)
                    v = NaN;
                    return;
                end
                x = x{1};
            end
            if isnumeric(x) || islogical(x)
                v = double(x);
                if isempty(v)
                    v = NaN;
                end
                return;
            end
            if isa(x, 'missing') || (isstring(x) && any(ismissing(x)))
                v = NaN;
                return;
            end
            s = strtrim(strrep(char(x), '%', ''));
            if isempty(s)
                v = NaN;
            else
                v = str2double(s);
            end
        end

        function v = parseFlag(x)
            %0/1 from '0'/'1', 'true'/'false', a logical or a number; NaN
            %for anything else.
            if iscell(x)
                x = x{1};
            end
            if islogical(x) || isnumeric(x)
                v = double(x);
                return;
            end
            if isa(x, 'missing')
                v = NaN;
                return;
            end
            switch lower(strtrim(char(string(x))))
                case {'1', 'true'}
                    v = 1;
                case {'0', 'false'}
                    v = 0;
                otherwise
                    v = NaN;
            end
        end

        function d = decimals(s)
            %Absolute decimal places a printed number carries: digits after
            %the point, less the exponent when one is printed ('1.234e-05'
            %carries 8).
            s = strtrim(char(s));
            d = 0;
            dot = find(s == '.', 1);
            if isempty(dot)
                return;
            end
            k = dot + 1;
            while k <= numel(s) && s(k) >= '0' && s(k) <= '9'
                d = d + 1;
                k = k + 1;
            end
            ex = regexp(s, '[eE]([-+]?\d+)$', 'tokens', 'once');
            if ~isempty(ex)
                d = d - str2double(ex{1});
            end
        end

        %% ---- dialog doubles -----------------------------------------------

        function d = installShadows()
            %Write the three recording doubles to a fresh folder and put it
            %first on the path. The caller removes it (removeShadows).
            d = tempname;
            mkdir(d);
            put = { ...
                'function [filename, pathname, filterindex] = uiputfile(varargin)'
                '%TEST SHADOW (G107/G109): answers from a queue; a 0 entry or an empty queue is Cancel.'
                'global G107_UIPUTFILE_QUEUE G107_UIPUTFILE_CALLS'
                'if isempty(G107_UIPUTFILE_CALLS), G107_UIPUTFILE_CALLS = {}; end'
                'G107_UIPUTFILE_CALLS{end+1} = varargin;'
                'filename = 0; pathname = 0; filterindex = 0;'
                'if isempty(G107_UIPUTFILE_QUEUE), return; end'
                'answer = G107_UIPUTFILE_QUEUE{1};'
                'G107_UIPUTFILE_QUEUE(1) = [];'
                'if isequal(answer, 0), return; end'
                '[p, n, e] = fileparts(answer);'
                'filename = [n e]; pathname = [p filesep]; filterindex = 1;'
                'end'};
            quest = { ...
                'function answer = questdlg(varargin)'
                '%TEST SHADOW (G109): consents or declines from a global mode; every call is recorded.'
                'global G107_QUESTDLG_MODE G107_QUESTDLG_CALLS'
                'if isempty(G107_QUESTDLG_CALLS), G107_QUESTDLG_CALLS = {}; end'
                'G107_QUESTDLG_CALLS{end+1} = varargin;'
                'btns = {};'
                'for k = 3:numel(varargin)'
                '    a = varargin{k};'
                '    if ischar(a) || (isstring(a) && isscalar(a)), btns{end+1} = char(a); end %#ok<AGROW>'
                'end'
                'if isempty(btns), btns = {''Yes'', ''No'', ''Cancel''}; end'
                'answer = '''';'
                'if isequal(G107_QUESTDLG_MODE, ''consent'')'
                '    hit = find(~cellfun(@isempty, regexpi(btns, ''(replace|overwrite|^yes)'', ''once'')), 1);'
                '    if ~isempty(hit), answer = btns{hit}; end'
                'end'
                'end'};
            wait = { ...
                'function uiwait(varargin)'
                '%TEST SHADOW (G109): never blocks, so an unexpected modal cannot hang a headless run.'
                'global G107_UIWAIT_CALLS'
                'if isempty(G107_UIWAIT_CALLS), G107_UIWAIT_CALLS = 0; end'
                'G107_UIWAIT_CALLS = G107_UIWAIT_CALLS + 1;'
                'end'};
            G107ExportFixtures.writeLines(fullfile(d, 'uiputfile.m'), put);
            G107ExportFixtures.writeLines(fullfile(d, 'questdlg.m'), quest);
            G107ExportFixtures.writeLines(fullfile(d, 'uiwait.m'), wait);
            addpath(d, '-begin');
            rehash;
        end

        function removeShadows(d)
            %Drop the doubles and their globals so later tests see the real
            %dialogs.
            clear global G107_UIPUTFILE_QUEUE G107_UIPUTFILE_CALLS G107_QUESTDLG_MODE G107_QUESTDLG_CALLS G107_UIWAIT_CALLS
            if ~isempty(d) && any(strcmp(strsplit(path, pathsep), d))
                rmpath(d);
            end
            clear('uiputfile', 'questdlg', 'uiwait');
            rehash;
            if ~isempty(d) && isfolder(d)
                rmdir(d, 's');
            end
        end

        function writeLines(file, lines)
            fid = fopen(file, 'w');
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
        end

        function d = installFopenWatch(folder, refuse)
            %Put a recording fopen first on the path for the write-failure
            %tests. Every call is forwarded to the built-in fopen, except an
            %open for writing or appending of a path under FOLDER: that open
            %is recorded (fopenWrites) and, when REFUSE is true, refused with
            %fid = -1, the way a folder that refuses new files answers. It is
            %a narrow seam for one test at a time; removeFopenWatch takes it
            %off the path. The caller checks which('fopen') before relying
            %on it.
            global G107_FOPEN_WATCH G107_FOPEN_REFUSE G107_FOPEN_WRITES %#ok<GVMIS>
            d = tempname;
            mkdir(d);
            watch = { ...
                'function varargout = fopen(varargin)'
                '%TEST SHADOW (G109): forwards to the built-in fopen; records, and in refusal mode refuses, write opens under one folder.'
                'global G107_FOPEN_WATCH G107_FOPEN_REFUSE G107_FOPEN_WRITES'
                'n = max(1, nargout);'
                'if numel(varargin) >= 2 && ~isempty(G107_FOPEN_WATCH) && ...'
                '        (ischar(varargin{1}) || isstring(varargin{1})) && ...'
                '        (ischar(varargin{2}) || isstring(varargin{2}))'
                '    p = char(varargin{1});'
                '    perm = lower(char(varargin{2}));'
                '    if ~isempty(perm) && any(perm(1) == ''wa'') && startsWith(p, G107_FOPEN_WATCH)'
                '        if isempty(G107_FOPEN_WRITES), G107_FOPEN_WRITES = {}; end'
                '        G107_FOPEN_WRITES{end+1} = p;'
                '        if isequal(G107_FOPEN_REFUSE, true)'
                '            refused = {-1, ''Permission denied (G109 test shadow)''};'
                '            varargout = refused(1:min(n, 2));'
                '            return;'
                '        end'
                '    end'
                'end'
                '[varargout{1:n}] = builtin(''fopen'', varargin{:});'
                'end'};
            G107ExportFixtures.writeLines(fullfile(d, 'fopen.m'), watch);
            G107_FOPEN_WATCH = char(folder);
            G107_FOPEN_REFUSE = logical(refuse);
            G107_FOPEN_WRITES = {};
            addpath(d, '-begin');
            rehash;
        end

        function writes = fopenWrites()
            %The write/append opens under the watched folder since the watch
            %was installed.
            global G107_FOPEN_WRITES %#ok<GVMIS>
            writes = G107_FOPEN_WRITES;
            if isempty(writes)
                writes = {};
            end
        end

        function ids = openFileIDs()
            %Identifiers of every file MATLAB has open, for the leak checks.
            %openedFiles is the supported enumeration. fopen('all') is a
            %staged removal that raises MATLAB:fopen:FopenAllStagedRemoval on
            %R2025b, so it is only the fallback for a release that has no
            %openedFiles. Returned as a row either way, so a before/after
            %setdiff compares like with like.
            if ismember(exist('openedFiles'), [2 3 5 6]) %#ok<EXIST>
                ids = openedFiles();
            else
                ids = fopen('all');
            end
            ids = reshape(double(ids), 1, []);
        end

        function removeFopenWatch(d)
            clear global G107_FOPEN_WATCH G107_FOPEN_REFUSE G107_FOPEN_WRITES
            if ~isempty(d) && any(strcmp(strsplit(path, pathsep), d))
                rmpath(d);
            end
            clear('fopen');
            rehash;
            if ~isempty(d) && isfolder(d)
                rmdir(d, 's');
            end
        end

        function d = installIoFaults(folder, fails)
            %Put recording fclose, ferror and movefile doubles first on the
            %path for the write-failure tests. Each forwards to its built-in.
            %For a file under FOLDER each records the call and, while its
            %failure budget in FAILS (fields fclose, ferror, movefile; optional
            %movefileThrows and delete) lasts, reports a failure instead:
            %  fclose   - closes the file through the built-in FIRST, then
            %             returns -1 (a close that failed to flush);
            %  ferror   - returns a non-empty message (a write error the
            %             stream recorded without fprintf throwing);
            %  movefile - moves nothing and returns status 0 with a message
            %             (throws that message when called with no outputs,
            %             or always when the optional field movefileThrows is
            %             true: a move that raises instead of returning its
            %             status, review F1);
            %  delete   - OPT-IN: written only when FAILS has a delete field.
            %             Deletes nothing and returns quietly for a file path
            %             under FOLDER (a temporary file that cannot be
            %             removed). Handle objects (figures) still dispatch to
            %             their class delete method, not to this double.
            %Budgets count down one per watched call, so the criterion
            %writer's own header close is the first fclose it meets. A budget
            %of 0 is the control: every call is forwarded unchanged and only
            %recorded. ioFaultLog reports what was seen and what failed;
            %removeIoFaults takes the doubles off the path.
            global G107_IO_WATCH G107_IO_FAILS G107_IO_LOG %#ok<GVMIS>
            d = tempname;
            mkdir(d);
            closeDouble = { ...
                'function status = fclose(fid, varargin)'
                '%TEST SHADOW (G109): forwards to the built-in fclose; a watched file is really closed first, then may report -1.'
                'global G107_IO_WATCH G107_IO_FAILS G107_IO_LOG'
                'name = '''';'
                'if ~isempty(G107_IO_WATCH) && isnumeric(fid) && isscalar(fid) && fid > 2'
                '    name = builtin(''fopen'', fid);'
                'end'
                'status = builtin(''fclose'', fid, varargin{:});'
                'if ~isempty(name) && startsWith(name, G107_IO_WATCH)'
                '    G107_IO_LOG.fcloseSeen{end+1} = name;'
                '    if G107_IO_FAILS.fclose > 0'
                '        G107_IO_FAILS.fclose = G107_IO_FAILS.fclose - 1;'
                '        G107_IO_LOG.fcloseFailed{end+1} = name;'
                '        status = -1;'
                '    end'
                'end'
                'end'};
            errorDouble = { ...
                'function [msg, errnum] = ferror(fid, varargin)'
                '%TEST SHADOW (G109): forwards to the built-in ferror; for a watched file it may report a write error.'
                'global G107_IO_WATCH G107_IO_FAILS G107_IO_LOG'
                'name = '''';'
                'if ~isempty(G107_IO_WATCH) && isnumeric(fid) && isscalar(fid) && fid > 2'
                '    name = builtin(''fopen'', fid);'
                'end'
                'if ~isempty(name) && startsWith(name, G107_IO_WATCH)'
                '    G107_IO_LOG.ferrorSeen{end+1} = name;'
                '    if G107_IO_FAILS.ferror > 0'
                '        G107_IO_FAILS.ferror = G107_IO_FAILS.ferror - 1;'
                '        G107_IO_LOG.ferrorFailed{end+1} = name;'
                '        msg = ''No space left on device (G109 test shadow)'';'
                '        errnum = 28;'
                '        return;'
                '    end'
                'end'
                '[msg, errnum] = builtin(''ferror'', fid, varargin{:});'
                'end'};
            moveDouble = { ...
                'function varargout = movefile(varargin)'
                '%TEST SHADOW (G109): forwards to the built-in movefile; a move of a watched file may be refused, moving nothing.'
                'global G107_IO_WATCH G107_IO_FAILS G107_IO_LOG'
                'src = '''';'
                'if ~isempty(varargin) && (ischar(varargin{1}) || isstring(varargin{1}))'
                '    src = char(varargin{1});'
                'end'
                'if ~isempty(G107_IO_WATCH) && ~isempty(src) && startsWith(src, G107_IO_WATCH)'
                '    G107_IO_LOG.movefileSeen{end+1} = src;'
                '    if G107_IO_FAILS.movefile > 0'
                '        G107_IO_FAILS.movefile = G107_IO_FAILS.movefile - 1;'
                '        G107_IO_LOG.movefileFailed{end+1} = src;'
                '        msg = ''Operation not permitted (G109 test shadow)'';'
                '        if nargout == 0 || G107_IO_FAILS.movefileThrows'
                '            error(''G109Test:movefileRefused'', ''%s'', msg);'
                '        end'
                '        refused = {false, msg, ''G109Test:movefileRefused''};'
                '        varargout = refused(1:nargout);'
                '        return;'
                '    end'
                'end'
                'if nargout == 0'
                '    builtin(''movefile'', varargin{:});'
                'else'
                '    [varargout{1:nargout}] = builtin(''movefile'', varargin{:});'
                'end'
                'end'};
            deleteDouble = { ...
                'function varargout = delete(varargin)'
                '%TEST SHADOW (G109): forwards to the built-in delete; a budgeted delete of a file path under one folder deletes nothing.'
                'global G107_IO_WATCH G107_IO_FAILS G107_IO_LOG'
                'p = '''';'
                'if numel(varargin) == 1 && (ischar(varargin{1}) || (isstring(varargin{1}) && isscalar(varargin{1})))'
                '    p = char(varargin{1});'
                'end'
                'if ~isempty(G107_IO_WATCH) && ~isempty(p) && startsWith(p, G107_IO_WATCH)'
                '    G107_IO_LOG.deleteSeen{end+1} = p;'
                '    if G107_IO_FAILS.delete > 0'
                '        G107_IO_FAILS.delete = G107_IO_FAILS.delete - 1;'
                '        G107_IO_LOG.deleteFailed{end+1} = p;'
                '        return;'
                '    end'
                'end'
                'if nargout == 0'
                '    builtin(''delete'', varargin{:});'
                'else'
                '    [varargout{1:nargout}] = builtin(''delete'', varargin{:});'
                'end'
                'end'};
            G107ExportFixtures.writeLines(fullfile(d, 'fclose.m'), closeDouble);
            G107ExportFixtures.writeLines(fullfile(d, 'ferror.m'), errorDouble);
            G107ExportFixtures.writeLines(fullfile(d, 'movefile.m'), moveDouble);
            if isfield(fails, 'delete')
                G107ExportFixtures.writeLines(fullfile(d, 'delete.m'), deleteDouble);
            end
            G107_IO_WATCH = char(folder);
            G107_IO_FAILS = struct('fclose', fails.fclose, 'ferror', fails.ferror, ...
                'movefile', fails.movefile, ...
                'movefileThrows', isfield(fails, 'movefileThrows') && logical(fails.movefileThrows), ...
                'delete', 0);
            if isfield(fails, 'delete')
                G107_IO_FAILS.delete = fails.delete;
            end
            G107_IO_LOG = G107ExportFixtures.emptyIoLog();
            addpath(d, '-begin');
            rehash;
        end

        function log = emptyIoLog()
            %The ioFaultLog record before any watched call.
            log = struct('fcloseSeen', {{}}, 'fcloseFailed', {{}}, ...
                'ferrorSeen', {{}}, 'ferrorFailed', {{}}, ...
                'movefileSeen', {{}}, 'movefileFailed', {{}}, ...
                'deleteSeen', {{}}, 'deleteFailed', {{}});
        end

        function log = ioFaultLog()
            %What the fclose, ferror, movefile and (when installed) delete
            %doubles saw (watched calls) and which of those they failed,
            %since installIoFaults.
            global G107_IO_LOG %#ok<GVMIS>
            log = G107_IO_LOG;
            if isempty(log)
                log = G107ExportFixtures.emptyIoLog();
            end
        end

        function removeIoFaults(d)
            clear global G107_IO_WATCH G107_IO_FAILS G107_IO_LOG
            if ~isempty(d) && any(strcmp(strsplit(path, pathsep), d))
                rmpath(d);
            end
            clear('fclose', 'ferror', 'movefile', 'delete');
            rehash;
            if ~isempty(d) && isfolder(d)
                rmdir(d, 's');
            end
        end

        function d = installErrordlgRecorder()
            %Put a recording errordlg first on the path (G109 one-dialog
            %checks, review F2). Every call is recorded (errordlgCalls) and
            %then shown through msgbox with the error icon, which is how
            %errordlg itself draws, so capture still reads the window's text.
            %CALLS are counted, not windows: errordlg's 'modal' mode is
            %documented to reuse an open dialog with the same title, so two
            %'File Not Saved' dialogs could leave one window holding only the
            %second text. removeErrordlgRecorder takes it off the path. The
            %caller checks which('errordlg') before relying on it.
            global G107_ERRORDLG_CALLS %#ok<GVMIS>
            d = tempname;
            mkdir(d);
            rec = { ...
                'function varargout = errordlg(varargin)'
                '%TEST SHADOW (G109): records every call, then shows it through msgbox with the error icon, as errordlg does.'
                'global G107_ERRORDLG_CALLS'
                'if isempty(G107_ERRORDLG_CALLS), G107_ERRORDLG_CALLS = {}; end'
                'G107_ERRORDLG_CALLS{end+1} = varargin;'
                'args = {'''', ''Error Dialog'', ''non-modal''};'
                'n = min(3, numel(varargin));'
                'args(1:n) = varargin(1:n);'
                'h = msgbox(args{1}, args{2}, ''error'', args{3});'
                'if nargout > 0, varargout{1} = h; end'
                'end'};
            G107ExportFixtures.writeLines(fullfile(d, 'errordlg.m'), rec);
            G107_ERRORDLG_CALLS = {};
            addpath(d, '-begin');
            rehash;
        end

        function resetErrordlgCalls()
            %Forget the errordlg calls recorded so far.
            global G107_ERRORDLG_CALLS %#ok<GVMIS>
            G107_ERRORDLG_CALLS = {};
        end

        function calls = errordlgCalls()
            %Every errordlg call since the recorder was installed or reset,
            %each as its argument cell ({message, title, mode}).
            global G107_ERRORDLG_CALLS %#ok<GVMIS>
            calls = G107_ERRORDLG_CALLS;
            if isempty(calls)
                calls = {};
            end
        end

        function txt = errordlgText(calls)
            %The messages of the recorded CALLS, lower case, each message's
            %lines joined by spaces (so a phrase is found whatever the
            %window's line wrapping) and the messages joined by newlines.
            parts = cell(1, numel(calls));
            for k = 1:numel(calls)
                msg = {};
                if ~isempty(calls{k})
                    msg = G107ExportFixtures.flattenText(calls{k}{1});
                end
                parts{k} = strjoin(msg, ' ');
            end
            txt = lower(strjoin(parts, newline));
        end

        function removeErrordlgRecorder(d)
            clear global G107_ERRORDLG_CALLS
            if ~isempty(d) && any(strcmp(strsplit(path, pathsep), d))
                rmpath(d);
            end
            clear('errordlg');
            rehash;
            if ~isempty(d) && isfolder(d)
                rmdir(d, 's');
            end
        end

        function verifyReplaceFailureText(testCase, txt, name, what)
            %TXT (errordlgText) states that the finished file could not
            %replace NAME and asks the user to check it, and says nothing the
            %write-failure (abort) dialog says or that nobody observed.
            testCase.verifyTrue(contains(txt, ['could not replace ' lower(name)]), ...
                sprintf('%s: the dialog must say the new file could not replace %s', what, name));
            testCase.verifyTrue(contains(txt, 'check it before relying on it'), ...
                sprintf('%s: the dialog must ask the user to check %s', what, name));
            testCase.verifyFalse(contains(txt, 'could not be written completely'), ...
                sprintf('%s: a replace failure must not be reported as an incomplete write', what));
            testCase.verifyFalse(contains(txt, 'did not replace any existing file'), ...
                sprintf('%s: after a failed replacement nobody observed that nothing was replaced', what));
            testCase.verifyFalse(contains(txt, 'is unchanged'), ...
                sprintf('%s: after a failed replacement the dialog must not claim the file is unchanged', what));
        end

        function verifyWriteFailureText(testCase, txt, what)
            %TXT (errordlgText) is the write-failure (abort) report: the new
            %file could not be written completely, and it is not reported as
            %a failed replacement.
            testCase.verifyTrue(contains(txt, 'could not be written completely'), ...
                sprintf('%s: the dialog must say the new file could not be written completely', what));
            testCase.verifyFalse(contains(txt, 'could not replace'), ...
                sprintf('%s: a write failure must not be reported as a failed replacement', what));
        end

        function decoys = plantStageDecoys(folder)
            %Three entries named the way psyrat_export_stage names its
            %temporary files, none of them this export's: a file holding
            %bytes (psyrat_partial_1.csv), an EMPTY file
            %(psyrat_partial_2.csv) and a FOLDER (psyrat_partial_3.csv, with a
            %sentinel inside). An export beside them must leave all three
            %byte- and tree-identical. Returns their names.
            decoys = {'psyrat_partial_1.csv', 'psyrat_partial_2.csv', 'psyrat_partial_3.csv'};
            G107ExportFixtures.plant(fullfile(folder, decoys{1}), ...
                sprintf('G109 DECOY: a file this export did not create\n'));
            G107ExportFixtures.plant(fullfile(folder, decoys{2}), '');
            mkdir(fullfile(folder, decoys{3}));
            G107ExportFixtures.plant(fullfile(folder, decoys{3}, 'inside.txt'), ...
                sprintf('G109 DECOY FOLDER SENTINEL\n'));
        end

        function resetDialogs(queue, mode)
            %QUEUE: the uiputfile answers in order (full paths, or 0 for
            %Cancel). MODE: 'consent' or 'decline' for questdlg.
            global G107_UIPUTFILE_QUEUE G107_UIPUTFILE_CALLS G107_QUESTDLG_MODE G107_QUESTDLG_CALLS G107_UIWAIT_CALLS %#ok<GVMIS>
            G107_UIPUTFILE_QUEUE = queue;
            G107_UIPUTFILE_CALLS = {};
            G107_QUESTDLG_MODE = mode;
            G107_QUESTDLG_CALLS = {};
            G107_UIWAIT_CALLS = 0;
        end

        function st = dialogState()
            %Every recorded dialog call since the last resetDialogs.
            global G107_UIPUTFILE_CALLS G107_QUESTDLG_CALLS G107_UIWAIT_CALLS %#ok<GVMIS>
            st = struct();
            st.uiputfileCalls = G107_UIPUTFILE_CALLS;
            if isempty(st.uiputfileCalls)
                st.uiputfileCalls = {};
            end
            st.questdlgCalls = G107_QUESTDLG_CALLS;
            if isempty(st.questdlgCalls)
                st.questdlgCalls = {};
            end
            st.uiwaitCalls = G107_UIWAIT_CALLS;
            if isempty(st.uiwaitCalls)
                st.uiwaitCalls = 0;
            end
        end

        function tf = prompted(st)
            %True when the callback asked anything after the first save
            %dialog: a questdlg, or a second uiputfile.
            tf = ~isempty(st.questdlgCalls) || numel(st.uiputfileCalls) > 1;
        end

        function txt = promptText(st)
            %All text the callback put in front of the user AFTER the first
            %save dialog (questdlg arguments and any re-prompt's uiputfile
            %arguments, including its default name).
            parts = {};
            for k = 1:numel(st.questdlgCalls)
                parts = [parts, G107ExportFixtures.flattenText(st.questdlgCalls{k})]; %#ok<AGROW>
            end
            for k = 2:numel(st.uiputfileCalls)
                parts = [parts, G107ExportFixtures.flattenText(st.uiputfileCalls{k})]; %#ok<AGROW>
            end
            txt = strjoin(parts, newline);
        end

        function parts = flattenText(x)
            parts = {};
            if ischar(x)
                parts = cellstr(x)';
            elseif isstring(x)
                parts = cellstr(x(:))';
            elseif iscell(x)
                for k = 1:numel(x)
                    parts = [parts, G107ExportFixtures.flattenText(x{k})]; %#ok<AGROW>
                end
            end
        end

        %% ---- GUI driving --------------------------------------------------

        function figs = openWindows(testCase, builder)
            %Run a table builder (gui = 1) and return the windows it opened.
            before = findall(groot, 'Type', 'figure');
            builder();
            after = findall(groot, 'Type', 'figure');
            figs = after(~ismember(after, before));
            testCase.assertNotEmpty(figs, 'the builder opened no window');
        end

        function button = findButton(testCase, figs, labelfrag)
            %The pushbutton whose label contains LABELFRAG.
            button = [];
            btns = findall(figs, 'Style', 'pushbutton');
            for b = btns(:)'
                s = get(b, 'String');
                if iscell(s)
                    s = strjoin(s, ' ');
                end
                if contains(string(s), labelfrag)
                    button = b;
                    return;
                end
            end
            testCase.assertFail(sprintf('no "%s" button on the window', labelfrag));
        end

        function win = openChild(testCase, figs, callbackName)
            %Fire the button whose Callback is CALLBACKNAME and return the
            %window it opened (the TestUitableWidthSlack route to the View All
            %window).
            win = [];
            btns = findall(figs, 'Style', 'pushbutton');
            for b = btns(:)'
                cb = get(b, 'Callback');
                if iscell(cb) && isa(cb{1}, 'function_handle') && ...
                        strcmp(func2str(cb{1}), callbackName)
                    before = findall(groot, 'Type', 'figure');
                    feval(cb{1}, b, [], cb{2:end});
                    after = findall(groot, 'Type', 'figure');
                    win = after(~ismember(after, before));
                    testCase.assertNumElements(win, 1, ...
                        sprintf('%s must open one window', callbackName));
                    return;
                end
            end
            testCase.assertFail(sprintf('no button on the window calls %s', callbackName));
        end

        function out = fire(button, varargin)
            %Fire BUTTON's own Callback (feval of the stored cell, the
            %callback-invocation route). Arguments after BUTTON, when given,
            %REPLACE the stored callback arguments; that is how a malformed
            %summary is injected into a window built from a valid one.
            cb = get(button, 'Callback');
            if iscell(cb)
                args = cb(2:end);
                if ~isempty(varargin)
                    args = varargin;
                end
                fn = @() feval(cb{1}, button, [], args{:});
            else
                fn = @() cb(button, []);
            end
            out = G107ExportFixtures.capture(fn);
        end

        function out = capture(fn)
            %Run FN and report what a user would have seen: an uncaught
            %error (a GUI callback would print it to the Command Window),
            %and the text of every non-output window it opened (errordlg,
            %warndlg, msgbox).
            before = findall(groot, 'Type', 'figure');
            out = struct('error', [], 'errorText', '', 'dialogText', '', ...
                'dialogs', []);
            try
                fn();
            catch ME
                out.error = ME;
                out.errorText = ME.message;
            end
            after = findall(groot, 'Type', 'figure');
            newf = after(~ismember(after, before));
            txt = {};
            for f = newf(:)'
                if strcmp(get(f, 'Tag'), 'psyrat_output')
                    continue;
                end
                txt{end+1} = get(f, 'Name'); %#ok<AGROW>
                for h = findall(f, 'Type', 'text')'
                    txt = [txt, cellstr(string(get(h, 'String')))']; %#ok<AGROW>
                end
                for h = findall(f, 'Type', 'uicontrol')'
                    s = get(h, 'String');
                    if ~isempty(s)
                        txt = [txt, cellstr(string(s))']; %#ok<AGROW>
                    end
                end
            end
            out.dialogText = strjoin(txt, newline);
            out.dialogs = newf;
        end

        %% ---- files --------------------------------------------------------

        function b = fileBytes(file)
            %Exact bytes of FILE ([] when it does not exist).
            b = [];
            if ~isfile(file)
                return;
            end
            fid = fopen(file, 'r');
            b = fread(fid, Inf, '*uint8');
            fclose(fid);
        end

        function names = folderFiles(folder)
            %Sorted names of the files (not folders) in FOLDER, as a column.
            listing = dir(folder);
            listing = listing(~[listing.isdir]);
            names = {listing.name};
            names = sort(names(:));
        end

        function plant(file, text)
            %Write TEXT to FILE (a sentinel whose bytes must survive).
            fid = fopen(file, 'w');
            fwrite(fid, text);
            fclose(fid);
        end

        function tree = folderTree(folder)
            %Every file and folder below FOLDER, recursively (paths relative
            %to FOLDER, folders ending in filesep, sorted), with each file's
            %bytes. Unlike folderFiles, this sees a partial file dropped
            %inside a subfolder, and a changed byte anywhere. Call it with no
            %I/O doubles installed.
            base = dir(folder);
            base = base(strcmp({base.name}, '.'));
            root = base(1).folder;
            listing = dir(fullfile(folder, '**', '*'));
            listing = listing(~ismember({listing.name}, {'.', '..'}));
            paths = cell(numel(listing), 1);
            bytes = cell(numel(listing), 1);
            for k = 1:numel(listing)
                full = fullfile(listing(k).folder, listing(k).name);
                paths{k} = full(numel(root) + 2:end);
                if listing(k).isdir
                    paths{k} = [paths{k} filesep];
                    bytes{k} = [];
                else
                    bytes{k} = G107ExportFixtures.fileBytes(full);
                end
            end
            [paths, order] = sort(paths);
            tree = struct('paths', {paths}, 'bytes', {bytes(order)});
        end

        function name = nameOfLength(n, ext)
            %An ASCII file name exactly N characters long, ending in EXT.
            name = [repmat('a', 1, n - numel(ext)) ext];
        end

        function tf = canCreateName(folder, name)
            %True when a file called NAME can be created in FOLDER; the probe
            %file is removed again. Call it with no I/O doubles installed.
            fid = fopen(fullfile(folder, name), 'w');
            tf = fid >= 0;
            if tf
                fclose(fid);
                delete(fullfile(folder, name));
            end
        end

        function n = longestCreatableName(folder, ext, lo, hi)
            %The longest ASCII name ending in EXT that FOLDER accepts, found by
            %bisection between LO (must be creatable) and HI (must not be).
            %It measures whatever limit applies in that folder, a name-length
            %or a full-path limit, instead of assuming one. NaN when LO cannot
            %be created or HI can, so the caller fails with a message instead
            %of testing against a limit that is not there.
            n = NaN;
            if ~G107ExportFixtures.canCreateName(folder, G107ExportFixtures.nameOfLength(lo, ext)) || ...
                    G107ExportFixtures.canCreateName(folder, G107ExportFixtures.nameOfLength(hi, ext))
                return;
            end
            while hi - lo > 1
                mid = floor((lo + hi) / 2);
                if G107ExportFixtures.canCreateName(folder, G107ExportFixtures.nameOfLength(mid, ext))
                    lo = mid;
                else
                    hi = mid;
                end
            end
            n = lo;
        end

    end
end
