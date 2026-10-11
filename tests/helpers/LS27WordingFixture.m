classdef LS27WordingFixture
    %LS27WORDINGFIXTURE Draw-backed fixtures and independent arithmetic for the
    %dynamic-migration (Q8) facts that the corrected wording states: the
    %migration leaves the raw draws unchanged, while the next summary of a
    %stamped analysis-27 result derives different per-participant
    %coefficients. Used by TestRestrictedDisclosureWording and by the local
    %before/after trace.
    %
    %THE FIXTURES ARE THE REPOSITORY'S EXISTING DETERMINISTIC DRAWS, not new
    %data: PsyRATTestDataFactory's gamma analysis-27 and analysis-26 results
    %(log-nu slots), converted to the location-scale parameterization the way
    %TestGammaResultMigration's dynamic fixture converts them (pop_sdlog =
    %log 1.4, b_sigma and ind_sdlog taken from the log-nu slots, the log-nu
    %slots and, for the two-facet design, ind_bs removed), with gammascale set
    %explicitly (1 log-nu, 2 location-scale). The factory draws are constant
    %(50 identical rows), so their intervals are degenerate; variant 'varied'
    %spreads the same draws deterministically (fixed linear factors, no random
    %numbers) so that intervals are not.
    %
    %THE SCORE AND THE TWO TARGETS (stated before any arithmetic).
    %  Score: each participant's conditional two-facet coefficient at their
    %    own dimension value z_s (analysis 27, ic_dynrel_sserrvar_trt, gamma
    %    family): G_s (relative error) and D_s (absolute error), coefficient of
    %    equivalence and stability (reltype 3), occasion count n'_o = 1 (the
    %    summary default), trial count n'_i = the participant's own stored
    %    trls (median trials per occupied occasion).
    %  Components: the six observed-score mean components P, O, I, PO, PI, OI
    %    and the highest-order remainder R of the log-linked mean surface at
    %    alpha(z_s) = mu + b z_s, from the support-union enumeration below.
    %  Old target (legacy_observation_only_v1, an unstamped result):
    %    E_s = C_s, the participant's own observation variance.
    %  New target (population_highest_order_plus_person_observation_v2, a
    %    stamped result): E_s = R(z_s) + C_s.
    %  C_s: log-nu, 2 exp(2(alpha(z_s) + u_s) + 2 V_w - (log nu(z_s) + v_s))
    %    with V_w the sum of the five within-person log-mean variances;
    %    location-scale, exp(2(log sigma(z_s) + w_s)).
    %  G_s = P / (P + PI/n'_i + PO/n'_o + E_s/(n'_i n'_o)); D_s adds
    %    I/n'_i + O/n'_o + OI/(n'_i n'_o) to the error. Each is the posterior
    %    mean of the per-draw ratio, with an equal-tailed quantile interval.
    %  Consequence of the two targets alone: R >= 0 is added to every
    %    residual draw, so G_s and D_s can only fall, and fall strictly when R
    %    is positive. A participant then meets a fixed cutoff c before the
    %    migration and not after exactly when new < c <= old. Analysis 27 has
    %    no inclusion output of its own; that classification is the reader's.
    %
    %INDEPENDENCE. expected27 reads only the raw draw arrays in REL.out and
    %calls no PsyRAT function; components is the same support-union
    %enumeration TestGammaResultMigration uses as its independent reference
    %(frozen anchor: alpha = log 5, s_p = .6, s_o = s_i = .5, two-way SDs 0
    %gives R = 2.0652177285562118).

    methods (Static)

        function pd = dyn27(param, variant)
            %Analysis 27 ('lognu' | 'ls'), unstamped; variant 'constant'
            %(the factory draws) or 'varied'.
            if nargin < 2; variant = 'constant'; end
            pd = PsyRATTestDataFactory.makeGammaDynrelSubjNondiffTrtSummaryData();
            o = pd.rel.out;
            nd = size(o.mu, 1);
            f = ones(nd, 1); g = ones(nd, 1);
            if strcmp(variant, 'varied')
                f = linspace(0.8, 1.2, nd)';   % rising factor
                g = flipud(f);                 % falling factor, so draws are not collinear
                gro = cell2mat(o.gro_sds{1});
                o.gro_sds = {num2cell([gro(:,1) .* f, gro(:,2) .* g])};
                o.sig_occ = o.sig_occ .* g;
                o.sig_trl = o.sig_trl .* f;
                o.sig_trlxid = o.sig_trlxid .* g;
                o.sig_occxid = o.sig_occxid .* f;
                o.sig_trlxocc = o.sig_trlxocc .* g;
                o.mu = o.mu + 0.05 .* (f - 1);
                o.pop_lognu = o.pop_lognu + 0.10 .* (g - 1);
                o.ind_bs = {num2cell(cell2mat(o.ind_bs{1}) .* f)};
                o.ind_nu = {num2cell(cell2mat(o.ind_nu{1}) .* g)};
            elseif ~strcmp(variant, 'constant')
                error('LS27WordingFixture:variant', 'variant must be ''constant'' or ''varied''.');
            end
            if strcmp(param, 'ls')
                o.pop_sdlog = log(1.4) + 0.10 .* (g - 1);
                o.b_sigma = o.b_nu;
                o.ind_sdlog = o.ind_nu;
                o = rmfield(o, {'pop_lognu', 'b_nu', 'ind_nu', 'ind_bs'});
                pd.rel.gammascale = 2;
            elseif strcmp(param, 'lognu')
                pd.rel.gammascale = 1;
            else
                error('LS27WordingFixture:param', 'param must be ''lognu'' or ''ls''.');
            end
            pd.rel.out = o;
        end

        function pd = dyn14(param)
            %Analysis 14 from the same draws: the group-level label, with the
            %per-participant slots removed.
            pd = LS27WordingFixture.dyn27(param, 'constant');
            pd.rel.analysis = 'ic_dynrel_trt';
            drop = intersect(fieldnames(pd.rel.out), {'ssinfo'; 'ind_bs'; 'ind_nu'; 'ind_sdlog'});
            if ~isempty(drop); pd.rel.out = rmfield(pd.rel.out, drop); end
        end

        function pd = dyn26(param)
            %Analysis 26 (one-facet, subject-level) from the factory, converted
            %as for analysis 27.
            pd = PsyRATTestDataFactory.makeGammaDynrelSubjNondiffSummaryData();
            o = pd.rel.out;
            if strcmp(param, 'ls')
                o.pop_sdlog = repmat(log(1.4), size(o.mu, 1), 1);
                o.b_sigma = o.b_nu;
                o.ind_sdlog = o.ind_nu;
                o = rmfield(o, {'pop_lognu', 'b_nu', 'ind_nu'});
                pd.rel.gammascale = 2;
            else
                pd.rel.gammascale = 1;
            end
            pd.rel.out = o;
        end

        function E = expected27(rel, target, ci)
            %Per-participant G_s and D_s for an analysis-27 REL under TARGET
            %'legacy' (E_s = C_s) or 'v2' (E_s = R + C_s), from the raw draws
            %only. Rows follow each stratum's ssinfo row order, strata in
            %order. Also returns each participant's mean R and C for the trace.
            if nargin < 3; ci = .95; end
            switch target
                case 'legacy'; addR = false;
                case 'v2';     addR = true;
                otherwise; error('LS27WordingFixture:target', 'target must be ''legacy'' or ''v2''.');
            end
            o = rel.out;
            ls = isfield(rel, 'gammascale') && isequal(rel.gammascale, 2);
            twodim = isfield(rel, 'ndim') && isequal(rel.ndim, 2);
            lo = (1 - ci) / 2;
            rows = cell(0, 12);
            for s = 1:numel(o.labels)
                gro = cell2mat(o.gro_sds{s});
                %log-mean SDs in the order P, O, I, PO, PI, OI
                sd = [gro(:,1), o.sig_occ(:,s), o.sig_trl(:,s), ...
                    o.sig_occxid(:,s), o.sig_trlxid(:,s), o.sig_trlxocc(:,s)];
                v = sd .^ 2;
                bmean = cell2mat(o.b{s});
                si = o.ssinfo{s};
                for r = 1:height(si)
                    if twodim
                        x = [si.z1(r); si.z2(r); si.z1(r) * si.z2(r)];
                    else
                        x = si.z1(r);
                    end
                    col = si.id2(r);
                    alpha = o.mu(:,s) + bmean * x;
                    K = LS27WordingFixture.components(alpha, v);
                    if ls
                        logsig = o.pop_sdlog(:,s) + cell2mat(o.b_sigma{s}) * x;
                        w = cell2mat(o.ind_sdlog{s});
                        C = exp(2 .* (logsig + w(:, col)));
                    else
                        lognu = o.pop_lognu(:,s) + cell2mat(o.b_nu{s}) * x;
                        u = cell2mat(o.ind_bs{s});
                        vn = cell2mat(o.ind_nu{s});
                        Vw = sum(v(:, 2:6), 2);
                        C = 2 .* exp(2 .* (alpha + u(:, col)) + 2 .* Vw - (lognu + vn(:, col)));
                    end
                    Eres = C;
                    if addR; Eres = K.R + C; end
                    ni = max(1, si.trls(r));
                    no = 1;
                    relerr = K.PI ./ ni + K.PO ./ no + Eres ./ (ni * no);
                    abserr = relerr + K.I ./ ni + K.O ./ no + K.OI ./ (ni * no);
                    G = K.P ./ (K.P + relerr);
                    D = K.P ./ (K.P + abserr);
                    rows(end + 1, :) = {s, si.id(r), ni, mean(G), quantile(G, lo), ...
                        quantile(G, 1 - lo), mean(D), quantile(D, lo), quantile(D, 1 - lo), ...
                        mean(K.R), mean(C), mean(Eres)}; %#ok<AGROW>
                end
            end
            E = cell2table(rows, 'VariableNames', {'stratum', 'id', 'n_i', 'gen_pt', ...
                'gen_ll', 'gen_ul', 'dep_pt', 'dep_ll', 'dep_ul', 'R_mean', 'C_mean', 'E_mean'});
        end

        function K = components(alpha, v)
            %Observed-score components of mu = exp(alpha + six independent
            %normal effects) with log-scale variances V (columns P, O, I, PO,
            %PI, OI). With q = expm1(v) and B = exp(2 alpha + sum v), each
            %component is B times the sum, over the subsets of effects whose
            %facet supports have that union, of the product of their q. Every
            %term is nonnegative, so nothing cancels.
            q = expm1(v);
            B = exp(2 .* alpha + sum(v, 2));
            supp = logical([1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1]);  % P O I PO PI OI -> {P,O,I}
            acc = zeros(size(v, 1), 8);
            for m = 1:63
                inT = logical(bitget(m, 1:6));
                u = any(supp(inT, :), 1);
                idx = 1 + u(1) + 2 * u(2) + 4 * u(3);
                acc(:, idx) = acc(:, idx) + prod(q(:, inT), 2);
            end
            K = struct('P', B .* acc(:,2), 'O', B .* acc(:,3), 'PO', B .* acc(:,4), ...
                'I', B .* acc(:,5), 'PI', B .* acc(:,6), 'OI', B .* acc(:,7), ...
                'R', B .* acc(:,8));
        end

    end
end
