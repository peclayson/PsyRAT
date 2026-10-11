classdef G108DiffFixture
    %G108DIFFFIXTURE Fixtures and independent reference arithmetic for G108 (the
    %analysis-8 difference caterpillar's reference line) and the RC-31 residual
    %(the analysis-7 per-event trials curve). Used by
    %TestDifferencePresentationContracts.
    %
    %INDEPENDENCE. No method in this class calls toolbox code. The builders only
    %assemble REL structures from plain arrays, and every expected value is
    %computed here from those arrays (or, for a saved result, from the raw draw
    %arrays in REL.out) with plain arithmetic. A match against production is
    %therefore a match against a separate derivation, not against a copy of
    %psyrat_diffrel, psyrat_rel_sing or PsyRATAccuracyOracle.
    %
    %THE FORMULAS REPRODUCED. The G108 contract adopts the group target the
    %summary already stores (group(g).diffscore.pt and .icc_pt), so the reference
    %arithmetic reproduces that target's definition rather than proposing one:
    %  universe-score variance of the contrast  u = s11 + s22 - 2*s12
    %  relative error at counts (n1,n2)         w1/n1 + w2/n2 - 2*c/h
    %  absolute error (dependability) adds      b11/n1 + b22/n2 - 2*b12/h
    %where h = 2/(1/n1 + 1/n2) is the harmonic mean of the two counts. That
    %divisor is the TOOLBOX covariance convention (RC-01/S13), not a claim about
    %Rocha et al. (2026) Table 6. The ICC is the same ratio at one observation
    %per event. Each coefficient is the posterior mean of the per-draw ratio with
    %an equal-tailed quantile interval.
    %
    %CONTRACT NAMES, kept here as constants so the tests name them in one place:
    %  META        the diffscore field recording the count basis of the stored
    %              group target, written beside the coefficient by the
    %              sing_diff_sserr helper of psyrat_relsummary. Subfields: obs
    %              (1x2 counts the coefficient was evaluated at), centmeas
    %              ('mean' | 'median'), population ('retained' | 'full_roster'),
    %              nretained (count).
    %  CAPTION_TAG the Tag of the exportable caption psyrat_ssrelplot draws on the
    %              analysis-8 difference caterpillar.

    properties (Constant)
        META = 'obsbasis';
        CAPTION_TAG = 'psyrat_ssrel_reference_caption';
        EVENTS = {'loss'; 'gain'};
    end

    methods (Static)

        % ------------------------------------------------------------------
        % G108 group specifications (analysis 8, ic_diff_sserrvar)
        % ------------------------------------------------------------------

        function spec = groupA(name)
            %Primary group. Six participants with unequal, skewed paired counts
            %(so the retained mean differs from the retained median, and both
            %differ from the full-roster values), a residual correlation of .6
            %at the population and participant levels, a nonzero trial
            %covariance block, and heterogeneous participant residual variances.
            %
            %Hand-checked while authoring, at diffgcoeff 1, depcutoff .50,
            %meascutoff 2 (point estimate), mean: a3, a4, a5 and a6 are retained,
            %at mean counts [20 15.25]; the stored target is about .695 and the
            %one-facet fallback line the plot drew before G108 about .650; the
            %two lines give different red masks on both the coefficient and the
            %ICC panels. The tests re-establish every one of these separations at
            %run time from the independent arithmetic; nothing depends on the
            %rounded numbers in this comment.
            if nargin < 1; name = ''; end
            spec = struct();
            spec.name = name;
            spec.ids  = {'a1'; 'a2'; 'a3'; 'a4'; 'a5'; 'a6'};
            spec.n1   = [4; 6; 8; 20; 40; 12];     % trials, event 1 (loss)
            spec.n2   = [3; 5; 12; 10; 30; 9];     % trials, event 2 (gain)
            spec.k    = [4; 1; 0.05; 2; 0.5; 1.2]; % participant residual-variance multipliers
            spec.s    = [0.4; 0.8; 1.2; 1.6];      % per-draw scale of the person block
            spec.V    = [1.0 0.5; 0.5 1.0];        % person (co)variance at scale 1
            spec.sd   = [2 2];                     % population residual SD per event
            spec.rpop = 0.6;                       % population residual correlation
            spec.Bt   = [2.0 0.4; 0.4 2.0];        % trial (co)variance
            spec.rsub = 0.6;                       % participant residual correlation
        end

        function spec = groupB(name)
            %Second group, deliberately low in reliability: small counts, a weak
            %residual correlation (.2) and UNEQUAL event variances (person,
            %residual and trial blocks all differ by event). At depcutoff .50
            %nobody is retained under either coefficient (hand-checked best
            %point estimate about .41 for dependability and .46 for
            %generalizability), so a two-group run exercises the full-roster
            %fallback in one panel and partial retention in the other.
            if nargin < 1; name = 'GB'; end
            spec = struct();
            spec.name = name;
            spec.ids  = {'b1'; 'b2'; 'b3'; 'b4'; 'b5'};
            spec.n1   = [3; 4; 6; 5; 8];
            spec.n2   = [4; 3; 5; 7; 6];
            spec.k    = [2; 1.5; 1; 1.2; 0.8];
            spec.s    = [0.4; 0.8; 1.2; 1.6];
            spec.V    = [1.2 0.4; 0.4 0.9];
            spec.sd   = [3 2.5];
            spec.rpop = 0.2;
            spec.Bt   = [1.5 0.3; 0.3 1.2];
            spec.rsub = 0.2;
        end

        function [pd, truth] = build(specs)
            %Assemble an ic_diff_sserrvar psyrat_data (REL only, no summary) from
            %one or more group specifications, plus the plain-array truth the
            %reference arithmetic uses. A single spec with an empty name builds
            %an ungrouped run (rel.groups = 'none', no group column).
            if isstruct(specs); specs = {specs}; end
            G = numel(specs);
            events = G108DiffFixture.EVENTS;
            multi = G > 1 || ~isempty(specs{1}.name);

            %REL.out blocks: one cell column per group, the shapes a fitted
            %analysis-8 result carries (see PsyRATTestDataFactory
            %.makeICDiffSSErrRelData for the single-group original)
            idVc = cell(1,G); trVc = cell(1,G); bSc = cell(1,G); eVc = cell(1,G);
            e1c = cell(1,G); e2c = cell(1,G); wcc = cell(1,G); idmc = cell(1,G);
            elc = cell(1,G); tc = cell(1,G);
            rowId = {}; rowEv = {}; rowGr = {}; rowMeas = [];

            for g = 1:G
                sp = specs{g};
                nd = numel(sp.s);
                n  = numel(sp.ids);

                %population residual covariance and the participant residual
                %variances/covariances (constant across draws; only the person
                %block varies, which is what gives the intervals their width)
                C  = sp.rpop * sp.sd(1) * sp.sd(2);
                w1 = sp.k(:)' .* sp.sd(1)^2;          % 1 x n, event 1
                w2 = sp.k(:)' .* sp.sd(2)^2;          % 1 x n, event 2
                c  = sp.rsub .* sqrt(w1 .* w2);       % 1 x n

                A = zeros(nd,2,2); B = zeros(nd,2,2); S = zeros(nd,2); E = zeros(nd,2,2);
                for d = 1:nd
                    A(d,:,:) = sp.s(d) .* sp.V;
                    B(d,:,:) = sp.Bt;
                    S(d,:)   = log(sp.sd);                        % log residual SD
                    E(d,:,:) = [sp.sd(1)^2 C; C sp.sd(2)^2];
                end
                idVc{g} = num2cell(A);
                trVc{g} = num2cell(B);
                bSc{g}  = num2cell(S);
                eVc{g}  = num2cell(E);
                e1c{g}  = repmat(0.5 .* log(w1), nd, 1);          % log residual SD per participant
                e2c{g}  = repmat(0.5 .* log(w2), nd, 1);
                wcc{g}  = repmat(c, nd, 1);
                idmc{g} = table(sp.ids(:), (1:n)', 'VariableNames', {'id','id2'});
                elc{g}  = events;

                %one row per single-trial score; the score values never enter a
                %reliability coefficient, only the row counts do
                for j = 1:n
                    for t = 1:sp.n1(j)
                        rowId{end+1,1} = sp.ids{j}; %#ok<AGROW>
                        rowEv{end+1,1} = events{1}; %#ok<AGROW>
                        rowGr{end+1,1} = sp.name; %#ok<AGROW>
                        rowMeas(end+1,1) = j + 0.1 * t; %#ok<AGROW>
                    end
                    for t = 1:sp.n2(j)
                        rowId{end+1,1} = sp.ids{j}; %#ok<AGROW>
                        rowEv{end+1,1} = events{2}; %#ok<AGROW>
                        rowGr{end+1,1} = sp.name; %#ok<AGROW>
                        rowMeas(end+1,1) = j + 0.2 * t; %#ok<AGROW>
                    end
                end

                %plain-array truth for the reference arithmetic
                tr = struct();
                tr.name = sp.name;
                tr.ids  = sp.ids(:);
                tr.n1   = sp.n1(:);
                tr.n2   = sp.n2(:);
                tr.u    = sp.s(:) .* (sp.V(1,1) + sp.V(2,2) - 2 * sp.V(2,1));
                tr.W1   = repmat(sp.sd(1)^2, nd, 1);
                tr.W2   = repmat(sp.sd(2)^2, nd, 1);
                tr.C    = repmat(C, nd, 1);
                tr.b11  = repmat(sp.Bt(1,1), nd, 1);
                tr.b22  = repmat(sp.Bt(2,2), nd, 1);
                tr.b12  = repmat(sp.Bt(2,1), nd, 1);
                tr.w1   = repmat(w1, nd, 1);
                tr.w2   = repmat(w2, nd, 1);
                tr.c    = repmat(c, nd, 1);
                %the fallback plot's pop_errvar: the half-sum of the two event
                %residual variances averaged over draws AND participants
                tr.popwp = mean((w1 + w2) ./ 2);
                tc{g} = tr;
            end

            rel = struct();
            rel.analysis = 'ic_diff_sserrvar';
            if multi
                rel.groups = cellfun(@(sp) sp.name, specs(:), 'UniformOutput', false);
                rel.data = table(rowId, rowGr, rowMeas, rowEv, ...
                    'VariableNames', {'id','group','meas','event'});
            else
                rel.groups = 'none';
                rel.data = table(rowId, rowMeas, rowEv, ...
                    'VariableNames', {'id','meas','event'});
            end
            rel.events = events;
            rel.diff_names = events;
            rel.diffwpcov = 2;
            rel.filename = 'g108_fixture.psyrat';
            rel.family = 'gaussian';
            rel.engine = 'stan';
            rel.seed = 12345;
            rel.nchains = 1;
            rel.niter = numel(specs{1}.s);
            rel.out = struct();
            rel.out.b_sigma = bSc;
            rel.out.id_varcov = idVc;
            rel.out.trl_varcov = trVc;
            rel.out.err_varcov = eVc;
            rel.out.er_var_ss1 = e1c;
            rel.out.er_var_ss2 = e2c;
            rel.out.wp_cov_ss = wcc;
            rel.out.id_matches = idmc;
            rel.out.elabels = elc;
            if multi
                rel.out.glabels = rel.groups;
            else
                rel.out.glabels = {'none'};
            end
            rel.out.conv.data = repmat({{'name','n_eff','Rhat'; 'b_sigma', 200, 1.00}}, 1, G);

            pd = struct('rel', rel);
            pd.ver = '0-test';
            pd.proc = struct('measheader', 'meas');
            truth = [tc{:}];
        end

        function truth = truthFromRel(rel)
            %Extract the plain-array truth from any ic_diff_sserrvar REL -- a
            %fixture or a saved fit -- reading only REL.out and REL.data. Counts
            %are recounted from the data rows; participants without scores in
            %both events are dropped, the paired-roster rule the difference
            %table applies.
            events = cellstr(string(rel.events(:)));
            G = size(rel.out.id_varcov, 2);
            multi = ~(ischar(rel.groups) && strcmpi(rel.groups, 'none'));
            tc = cell(1, G);
            for g = 1:G
                A = G108DiffFixture.unwrap(rel.out.id_varcov{1,g});
                B = G108DiffFixture.unwrap(rel.out.trl_varcov{1,g});
                S = G108DiffFixture.unwrap(rel.out.b_sigma{1,g});
                nd = size(S, 1);
                if isfield(rel.out, 'err_varcov') && size(rel.out.err_varcov, 2) >= g && ...
                        ~isempty(rel.out.err_varcov{1,g})
                    E = G108DiffFixture.unwrap(rel.out.err_varcov{1,g});
                    C = E(:,2,1);
                else
                    C = zeros(nd, 1);
                end
                e1 = G108DiffFixture.unwrap(rel.out.er_var_ss1{1,g});
                e2 = G108DiffFixture.unwrap(rel.out.er_var_ss2{1,g});
                wc = G108DiffFixture.unwrap(rel.out.wp_cov_ss{1,g});
                idm = rel.out.id_matches{1,g};

                data = rel.data;
                name = '';
                if multi
                    gn = G108DiffFixture.groupNames(rel);
                    name = gn{g};
                    if ismember('group', data.Properties.VariableNames)
                        data = data(strcmp(string(data.group), string(name)), :);
                    end
                end

                ids = cellstr(string(idm.id));
                n1 = zeros(numel(ids), 1);
                n2 = zeros(numel(ids), 1);
                for j = 1:numel(ids)
                    isid = strcmp(string(data.id), string(ids{j}));
                    n1(j) = sum(isid & strcmp(string(data.event), string(events{1})));
                    n2(j) = sum(isid & strcmp(string(data.event), string(events{2})));
                end
                keep = n1 > 0 & n2 > 0;
                id2 = idm.id2(keep);

                t = struct();
                t.name = name;
                t.ids  = ids(keep);
                t.n1   = n1(keep);
                t.n2   = n2(keep);
                t.u    = A(:,1,1) + A(:,2,2) - 2 .* A(:,2,1);
                t.W1   = exp(S(:,1)).^2;
                t.W2   = exp(S(:,2)).^2;
                t.C    = C;
                t.b11  = B(:,1,1);
                t.b22  = B(:,2,2);
                t.b12  = B(:,2,1);
                t.w1   = exp(e1(:,id2)).^2;
                t.w2   = exp(e2(:,id2)).^2;
                t.c    = wc(:,id2);
                t.popwp = mean((exp(e1(:)).^2 + exp(e2(:)).^2) ./ 2);
                tc{g} = t;
            end
            truth = [tc{:}];
        end

        function ex = expected(truth, diffgcoeff, depcutoff, meascutoff, depcentmeas, ciperc)
            %Independent expectation for every group: each participant's
            %coefficient and ICC with intervals, the inclusion decision, the
            %count basis of the group target, the target itself (point, interval,
            %ICC), and a set of WRONG-BUT-PLAUSIBLE alternatives the tests use to
            %prove the fixture can tell the targets apart.
            lo = (1 - ciperc) / 2;
            hi = 1 - lo;
            cents = {'mean', 'median'};
            ec = cell(1, numel(truth));
            for g = 1:numel(truth)
                t = truth(g);
                n = numel(t.n1);

                P = struct('pt', zeros(n,1), 'll', zeros(n,1), 'ul', zeros(n,1), ...
                    'icc_pt', zeros(n,1), 'icc_ll', zeros(n,1), 'icc_ul', zeros(n,1));
                for j = 1:n
                    o = [max(1, t.n1(j)) max(1, t.n2(j))];
                    [r, ri] = G108DiffFixture.ratios(t.u, t.w1(:,j), t.w2(:,j), ...
                        t.c(:,j), t.b11, t.b22, t.b12, o, diffgcoeff);
                    P.pt(j) = mean(r);
                    P.ll(j) = quantile(r, lo);
                    P.ul(j) = quantile(r, hi);
                    P.icc_pt(j) = mean(ri);
                    P.icc_ll(j) = quantile(ri, lo);
                    P.icc_ul(j) = quantile(ri, hi);
                end

                switch meascutoff
                    case 1
                        incl = P.ll >= depcutoff;
                    case 2
                        incl = P.pt >= depcutoff;
                    case 3
                        incl = P.ul >= depcutoff;
                end

                %count basis: the retained participants, or the full roster
                %when nobody is retained
                if any(incl)
                    pool = incl;
                    population = 'retained';
                else
                    pool = true(n, 1);
                    population = 'full_roster';
                end
                obs = G108DiffFixture.central(t.n1(pool), t.n2(pool), depcentmeas);

                [R, Ri] = G108DiffFixture.ratios(t.u, t.W1, t.W2, t.C, ...
                    t.b11, t.b22, t.b12, obs, diffgcoeff);
                T = struct('pt', mean(R), 'll', quantile(R, lo), 'ul', quantile(R, hi), ...
                    'icc_pt', mean(Ri), 'icc_ll', quantile(Ri, lo), 'icc_ul', quantile(Ri, hi));

                %alternatives a wrong implementation could plausibly draw
                alt = struct();
                alt.othercent = mean(G108DiffFixture.ratios(t.u, t.W1, t.W2, t.C, ...
                    t.b11, t.b22, t.b12, ...
                    G108DiffFixture.central(t.n1(pool), t.n2(pool), 3 - depcentmeas), diffgcoeff));
                alt.fullroster = mean(G108DiffFixture.ratios(t.u, t.W1, t.W2, t.C, ...
                    t.b11, t.b22, t.b12, ...
                    G108DiffFixture.central(t.n1, t.n2, depcentmeas), diffgcoeff));
                alt.nocov = mean(G108DiffFixture.ratios(t.u, t.W1, t.W2, zeros(size(t.C)), ...
                    t.b11, t.b22, zeros(size(t.b12)), obs, diffgcoeff));
                alt.notrial = mean(G108DiffFixture.ratios(t.u, t.W1, t.W2, t.C, ...
                    t.b11, t.b22, t.b12, obs, 2));
                m = mean(min(t.n1, t.n2));
                alt.mincount = mean(G108DiffFixture.ratios(t.u, t.W1, t.W2, t.C, ...
                    t.b11, t.b22, t.b12, [m m], diffgcoeff));
                %the one-facet fallback psyrat_ssrelplot drew on a difference
                %table before G108 (the table carries no trl_var)
                alt.prefix = mean(t.u) / (mean(t.u) + t.popwp / m);
                alt.prefix_icc = mean(t.u) / (mean(t.u) + t.popwp);
                %the overall table's subject coefficient: a different quantity
                if any(incl)
                    alt.subjmean = mean(P.pt(incl));
                else
                    alt.subjmean = NaN;
                end

                e = struct();
                e.name = t.name;
                e.ids = t.ids;
                e.n1 = t.n1;
                e.n2 = t.n2;
                e.participants = P;
                e.included = incl;
                e.obs = obs;
                e.centmeas = cents{depcentmeas};
                e.population = population;
                e.nretained = sum(incl);
                e.target = T;
                e.alt = alt;
                ec{g} = e;
            end
            ex = [ec{:}];
        end

        function [r, ri] = ratios(u, v1, v2, cv, b11, b22, b12, obs, gc)
            %Per-draw difference-score coefficient (r) and ICC (ri) at counts
            %obs = [n1 n2]. gc = 1 dependability (absolute error, trial block
            %included), gc = 2 generalizability (relative error).
            h = 2 / (1 / obs(1) + 1 / obs(2));
            relerr = v1 ./ obs(1) + v2 ./ obs(2) - 2 .* cv ./ h;
            iccrel = v1 + v2 - 2 .* cv;
            if gc == 1
                err = relerr + b11 ./ obs(1) + b22 ./ obs(2) - 2 .* b12 ./ h;
                iccerr = iccrel + b11 + b22 - 2 .* b12;
            else
                err = relerr;
                iccerr = iccrel;
            end
            r = u ./ (u + err);
            ri = u ./ (u + iccerr);
        end

        function obs = central(a, b, which)
            %Central trial counts of the two events: 1 = mean, 2 = median.
            if which == 1
                obs = [mean(a) mean(b)];
            else
                obs = [median(a) median(b)];
            end
        end

        function [red, ord] = redRanks(pt, ll, ul, line)
            %Plot positions (1 = lowest point estimate) of the participants
            %whose interval lies wholly above or below LINE, the red-marker
            %rule the caterpillar documents. Point estimates must be distinct
            %(the tests assert it), so this ascending sort matches the plot's
            %stable sortrows exactly.
            [~, ord] = sort(pt(:));
            llS = ll(ord);
            ulS = ul(ord);
            red = find(ulS < line | llS > line).';
        end

        function m = margin(ll, ul, line)
            %Smallest distance between LINE and any interval endpoint. A mask
            %comparison is only meaningful when this is clear of rounding.
            m = min(abs([ll(:); ul(:)] - line));
        end

        % ------------------------------------------------------------------
        % RC-31 residual (analysis 7, ic_diff)
        % ------------------------------------------------------------------

        function [pd, arrays] = rc31Build()
            %A one-group ic_diff psyrat_data (REL only) with a NONZERO trial SD
            %in both events, so dependability and generalizability curves differ
            %materially at every n. Field shapes mirror a fitted analysis-7 REL
            %(and psyrat_dod_build_virtual_rel's 'diff' output). ARRAYS carries
            %the per-event variance draws the independent curve uses.
            events = G108DiffFixture.EVENTS;
            nd = 6;
            drift = (0:nd-1)' .* 0.02;
            sdId  = [1.00 + drift, 0.80 + drift];   % person SD per event
            sdTrl = [0.90 + drift, 1.10 - drift];   % trial SD per event (nonzero)
            sdRes = [1.20 + drift, 1.40 + drift];   % residual SD per event
            rId = 0.30; rTrl = 0.20; rRes = 0.25;
            A = zeros(nd,2,2); B = zeros(nd,2,2); E = zeros(nd,2,2);
            for d = 1:nd
                A(d,:,:) = G108DiffFixture.cov2(sdId(d,:), rId);
                B(d,:,:) = G108DiffFixture.cov2(sdTrl(d,:), rTrl);
                E(d,:,:) = G108DiffFixture.cov2(sdRes(d,:), rRes);
            end

            ids = arrayfun(@(k) sprintf('p%02d', k), (1:8)', 'UniformOutput', false);
            n1 = [10; 12; 8; 15; 9; 11; 14; 7];
            n2 = [9; 13; 10; 12; 8; 14; 11; 6];
            rowId = {}; rowEv = {}; rowMeas = [];
            for j = 1:numel(ids)
                for t = 1:n1(j)
                    rowId{end+1,1} = ids{j}; %#ok<AGROW>
                    rowEv{end+1,1} = events{1}; %#ok<AGROW>
                    rowMeas(end+1,1) = j + 0.1 * t; %#ok<AGROW>
                end
                for t = 1:n2(j)
                    rowId{end+1,1} = ids{j}; %#ok<AGROW>
                    rowEv{end+1,1} = events{2}; %#ok<AGROW>
                    rowMeas(end+1,1) = j + 0.2 * t; %#ok<AGROW>
                end
            end

            rel = struct();
            rel.analysis = 'ic_diff';
            rel.groups = 'none';
            rel.events = events;
            rel.data = table(rowId, rowMeas, rowEv, 'VariableNames', {'id','meas','event'});
            rel.diff_names = events;
            rel.diffwpcov = 2;
            rel.diffrescor = 1;
            rel.filename = 'rc31_fixture.psyrat';
            rel.family = 'gaussian';
            rel.engine = 'stan';
            rel.seed = 12345;
            rel.out = struct();
            rel.out.b = {num2cell(zeros(nd, 2))};
            rel.out.b_sigma = {num2cell(log(sdRes))};
            rel.out.sd_id = {num2cell(sdId)};
            rel.out.sd_trl = {num2cell(sdTrl)};
            rel.out.id_varcov = {num2cell(A)};
            rel.out.trl_varcov = {num2cell(B)};
            rel.out.rescor = {num2cell(repmat(rRes, nd, 1))};
            rel.out.err_varcov = {num2cell(E)};
            rel.out.elabels = events;
            rel.out.glabels = {'none'};
            rel.out.conv.data = {{'name','n_eff','Rhat'; 'b_sigma', 200, 1.00}};

            pd = struct('rel', rel);
            pd.ver = '0-test';
            pd.proc = struct('measheader', 'meas');

            arrays = struct();
            arrays.events = events;
            arrays.su2 = sdId .^ 2;    % person variance draws, nd x 2
            arrays.se2 = sdRes .^ 2;   % residual variance draws, nd x 2
            arrays.si2 = sdTrl .^ 2;   % trial variance draws, nd x 2
        end

        function [yD, yG] = onefacetCurves(su2, se2, si2, N)
            %Per-event one-facet curve over n = 1..N from variance draws:
            %dependability su2/(su2 + (se2 + si2)/n), generalizability
            %su2/(su2 + se2/n); posterior mean of the per-draw ratio.
            yD = zeros(1, N);
            yG = zeros(1, N);
            for n = 1:N
                yD(n) = mean(su2 ./ (su2 + (se2 + si2) ./ n));
                yG(n) = mean(su2 ./ (su2 + se2 ./ n));
            end
        end

        % ------------------------------------------------------------------
        % small utilities
        % ------------------------------------------------------------------

        function gn = groupNames(rel)
            %Group labels as a cellstr column; {''} for an ungrouped run. A
            %single group stored as a char is one name, not a column of letters.
            if ischar(rel.groups) && strcmpi(rel.groups, 'none')
                gn = {''};
            elseif ischar(rel.groups)
                gn = {rel.groups};
            else
                gn = cellstr(string(rel.groups(:)));
            end
        end

        function a = unwrap(v)
            %Posterior blocks are stored either as numeric arrays or as
            %num2cell of them; return numeric either way.
            if iscell(v)
                a = cell2mat(v);
            else
                a = v;
            end
        end

        function M = cov2(sd, r)
            %2 x 2 covariance matrix from two SDs and a correlation.
            M = [sd(1)^2, r * sd(1) * sd(2); r * sd(1) * sd(2), sd(2)^2];
        end

        function name = coefName(gc)
            %Coefficient label the toolbox prints for gcoeff/diffgcoeff = gc.
            if gc == 1
                name = 'Dependability';
            else
                name = 'Generalizability';
            end
        end

    end
end
