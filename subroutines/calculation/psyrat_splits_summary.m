function summary = psyrat_splits_summary(REL,varargin)
%Summarize a nonparallel data-splits run (Rocha Tables 4/5, observed design) into
%per-stratum dependability/generalizability coefficients, SEMs, and variance
%components.
%
% summary = psyrat_splits_summary(REL,'CI',.95,'obs',[],'nocc',[])
%
%This is the single place that turns the per-component posterior draws stored by
%psyrat_computevarcomp (case 23 = ic_splits, single occasion; case 24 =
%trt_splits, multi-occasion) into the reported reliability coefficients at the
%OBSERVED score. Case24 equally weights available occasion means, then averages
%universe/error variances equally over persons before forming each draw's ratio.
%Its numeric count overrides retain a balanced reference-design projection.
%There is no D-study projection range and no cut-score search: split means
%lose the item-level partition, so a projection over hypothetical item counts is
%not identified. The split designs deliberately bypass psyrat_relsummary's
%trial-cutoff machinery (it does not fit an observed-design split summary), the
%same way psyrat_dynrel_summary handles the dynamic-reliability family.
%
%The full derivation and the component->coefficient mapping are recorded in
%the development splits-formula note (not shipped; the Golden-Rule-1 traceability
%record). For case23, the per-item residual sig_err lumps the item main effect
%with the person-involving item residual; that case retains its corrected
%person-averaged target and existing G convention.
%Case24 D is exact conditional on the specified score and fitted restricted
%independent-effects model (independent split-row residuals). Its conventional
%relative partition has no blanket bound or estimation-bias guarantee. The
%fitted split identity is the shared ordinal within each person/occasion cell,
%not a physical row label; this is not unrestricted Rocha Table 5.
%
%Inputs
% REL - the result struct from psyrat_computevarcomp with REL.analysis ==
%  'ic_splits' (REL.out carries sig_u/sig_trl/sig_splitxid/sig_err) or
%  'trt_splits' (REL.out carries, per group/event stratum, sig_id/sig_occ/sig_trl/
%  sig_trlxid/sig_occxid/sig_trlxocc/sig_posxid/sig_err). REL.data holds the long
%  table (id/meas/weight[/group/event/time]).
%
%Optional Inputs
% CI   - credible-interval width in (0,1], default .95.
% obs  - Numeric split-count override, or [] (default) for the observed score.
%        Case23 retains its corrected effective person-averaged design. For
%        case24 any numeric obs or nocc selects a balanced projection; an
%        unspecified split count resolves to the rounded median cell count.
% nocc - [] or 'observed' retains the observed score unless obs is numeric.
%        A numeric value selects a balanced projection. An unspecified count
%        in a projection resolves to the stratum occasion union. Case23 inert.
%
%Output (struct)
% summary.analysis, summary.ci, summary.gbias_note
% summary.strata(s) for each stratum s:
%   .label   - stratum name ('measure' for the single-occasion design)
%   .nsplit  - case24 reference split count (not an irregular effective count)
%   .nocc    - case24 reference occasion count (1 for ic_splits)
%   .nbar_i  - case23 effective item scale; case24 pooled harmonic reference
%   .coef    - struct array of reported coefficients (1 for ic_splits; the three
%              reltypes CE/CS/CES for trt_splits). Each has .label, .D, .G,
%              .SEM_abs, .SEM_rel (each [ll pt ul]); .gbias is true in case23
%              and false in case24 (no bound asserted, not an unbiasedness claim).
%   .design  - case24 observed incidence, score weights, actual counts and
%              projection provenance; .mean_factors is observed A/B/C/W,
%              .coefficient_factors is the design used for the coefficients.
%   .varcomp - table of reported variance components (Greek labels; ll/pt/ul)

% Copyright (C) 2016-2026 Peter E. Clayson
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

ciperc  = psyrat_opt(varargin,'CI',.95);
obs_in  = psyrat_opt(varargin,'obs',[]);
nocc_in = psyrat_opt(varargin,'nocc',[]);

%obs (n_s) and nocc (n_o) overrides are single observed-design points: the split
%design reports the observed design only (no D-study projection range), so a
%non-scalar override is rejected rather than silently malforming the [ll pt ul]
%coefficients (a 2-element obs would make psyrat_rel_sing/trt return row vectors).
if ~isempty(obs_in) && ~(isnumeric(obs_in) && isscalar(obs_in))
    error('splits:badobs',...
        ['''obs'' must be a numeric scalar (splits per person at the observed '...
        'design). A range is not supported - the split design reports the '...
        'observed design only, with no projection.']);
end
if ~isempty(nocc_in) && ~(isnumeric(nocc_in) && isscalar(nocc_in)) && ...
        ~((ischar(nocc_in) || isstring(nocc_in)) && strcmpi(char(nocc_in),'observed'))
    error('splits:badnocc',...
        '''nocc'' must be a numeric scalar or the string ''observed''.');
end

%accepted-set guard: psyrat_analysis_family is the single source of truth for
%which analysis strings are data-splits variants (shared with the GUI router and
%the headless report dispatch).
if ~isfield(REL,'analysis') || ...
        ~strcmp(psyrat_analysis_family(REL.analysis),'splits')
    error('splits:badinput',...
        ['psyrat_splits_summary expects a nonparallel data-splits result '...
        '(REL.analysis = ''ic_splits'' or ''trt_splits'').']);
end

issingle = strcmp(REL.analysis,'ic_splits');

%stratum labels: the single-occasion design is one stratum; the multi-occasion
%design is chunked per group/event by psyrat_computevarcomp (case 24) with one
%REL.out column per stratum, in REL.out.labels order.
if issingle
    labels = {'measure'};
else
    labels = REL.out.labels;
    if ischar(labels); labels = {labels}; end
end
nstrata = numel(labels);

%reltype labels for the multi-occasion design (psyrat_rel_trt reltype order)
reltype_label = {'Equivalence (CE)','Stability (CS)',...
    'Equivalence and stability (CES)'};

summary = struct();
summary.analysis = REL.analysis;
summary.ci = ciperc;
summary.gbias_note = ['Generalizability (relative) is a downward-biased '...
    'observed-design lower bound from split means (the per-item residual lumps '...
    'the item main effect into relative error); the relative SEM is '...
    'correspondingly an upper bound. Dependability (absolute) is exact for '...
    'the person-averaged error. When participants contribute different '...
    'numbers of splits the coefficients use each participant''s own counts, '...
    'and the exclusion of the split main effect from relative error is then '...
    'approximate.'];
summary.strata = struct('label',{},'nsplit',{},'nocc',{},'nbar_i',{},...
    'coef',{},'varcomp',{});
if ~issingle
    summary.gbias_note = ['Dependability is exact conditional on the stated score '...
        'and fitted restricted independent-effects model with independent split-row '...
        'residuals and shared ordinal split identities. Generalizability retains the '...
        'conventional relative partition; no blanket lower-bound, estimation-bias, '...
        'or interval-coverage guarantee is asserted. Irregular incidence makes the '...
        'relative partition approximate. In observed summaries, missing occasions target available occasion '...
        'means, not an unobserved complete schedule; missingness is not modeled.'];
    summary.design_note = ['Observed: equal-weight mean of available occasion means; '...
        'universe and error variances are averaged equally over participants before '...
        'ratios and SEMs, separately within each stratum. CE/CS retain each person''s '...
        'fixed-facet support and can have person-specific conditional targets. Any '...
        'numeric split/occasion override selects a balanced projection with pooled '...
        'harmonic items per split. Scalar counts are reference counts; actual '...
        'observed incidence and projection provenance are recorded in design metadata.'];
    summary.strata = struct('label',{},'nsplit',{},'nocc',{},'nbar_i',{},...
        'coef',{},'varcomp',{},'design',{});
end

for s = 1:nstrata
    %--- observed design for this stratum (recomputed from REL.data; the
    %estimation stores no n_s/n_i, and recomputing keeps that committed core
    %untouched) ---
    [rows,group_label,event_label] = local_stratum_rows(REL,labels{s},issingle);
    w = double(REL.data.weight(rows));
    %defensive: a stratum label is built from the same normalized group/event
    %columns it matches against, so it always selects >=1 row in the GUI
    %pipeline (this branch is unreachable there). Guard a hand-built or
    %inconsistent REL against a silent NaN coefficient (nbar_i = 0/0 below, and
    %median([]) further down) by failing loud instead.
    if isempty(w)
        error('psyrat_splits_summary:emptyStratum',...
            ['Stratum ''%s'' matched no observations in REL.data. How to fix: '...
            'ensure REL.data''s group/event labels are consistent with the '...
            'stratum labels.'], labels{s});
    end
    %harmonic mean of items-per-split: the exact residual scale for an
    %equal-weight mean of split means (owner decision 2026-06-25; see the
    %derivation note section 5). Pooled over the stratum it is exact for the
    %person-averaged error only for a common complete design. Case23 replaces
    %it by its effective item scale when counts differ; case24 retains it
    %as projection/reference metadata and uses incidence weights otherwise.
    nbar_i = numel(w) / sum(1 ./ w);

    idsub = REL.data.id(rows);
    if issingle
        [~,~,g] = unique(idsub);
        counts = accumarray(g,1);                       % splits per person
        nsplit_obs = median(counts);
        nocc_obs = 1;
        %PER-PERSON error terms averaged over persons (owner ruling 2026-09-10,
        %work queue Q7; post-beta review PBA-06). Person p's absolute error is
        %  (sigma_ps^2 + sigma_s^2)/S_p + sigma_e^2 * sum_s(1/n_ps) / S_p^2,
        %and its mean over persons is the Table 4 closed form evaluated at
        %  obs = 1/mean_p(1/S_p)   and   nbar_i = 1/kappa,
        %  kappa = obs * mean_p(sum_s(1/n_ps)/S_p^2),
        %so the coefficient code below is unchanged. With equal counts these
        %ARE the median count and the pooled harmonic mean, so a balanced design
        %keeps the arithmetic above bit for bit; with unequal counts the median
        %and the pooled mean misweighted the persons (the ledger's probe
        %reported D = 0.683 for a true 0.578). A user-supplied 'obs' keeps its
        %projection meaning with the pooled items per split. The blind
        %re-derivation of 2026-09-10 agreed with these expressions; see
        %the development splits-formula note, section 5.
        if isempty(obs_in) && any(counts ~= counts(1))
            sumInv = accumarray(g, 1 ./ w);             % sum_s(1/n_ps) per person
            nsplit_obs = 1 / mean(1 ./ counts);
            nbar_i = 1 / (nsplit_obs * mean(sumInv ./ counts.^2));
        end
    else
        timesub = REL.data.time(rows);
        design = local_trt_design(idsub,timesub,w,find(rows),...
            group_label,event_label);
        nsplit_obs = median(design.cell_split_counts);
        nocc_obs = design.occasion_union;
    end

    %the observed default is reported as computed: a whole number of splits
    %when the counts are equal, the effective (harmonic-mean) count when the
    %single-occasion counts differ; case 24 keeps its rounded per-cell median,
    %and a user override is a whole number of splits, as before
    if isempty(obs_in)
        nsplit = nsplit_obs;
        if ~issingle
            nsplit = max(1,round(nsplit));
        end
    else
        nsplit = max(1,round(obs_in));
    end
    nocc = local_resolve_nocc(nocc_in,nocc_obs);
    if ~issingle
        design.mode = 'observed_available_occasions';
        design.score_rule = 'Equal-weight mean of available occasion means';
        design.projection_nsplit = [];
        design.projection_nocc = [];
        design.coefficient_factors = design.mean_factors;
        design.obs_override = obs_in;
        design.nocc_override = [];
        if isnumeric(nocc_in); design.nocc_override = nocc_in; end
        if ~isempty(obs_in) || (~isempty(nocc_in) && isnumeric(nocc_in))
            design.mode = 'balanced_projection';
            design.score_rule = 'Balanced projection using pooled harmonic items per split';
            design.fixed_facet_target = ...
                'CE/CS retain common balanced reference support; observed missingness is not projected';
            design.projection_nsplit = nsplit;
            design.projection_nocc = nocc;
            design.coefficient_factors = [1/nsplit 1/nocc 1/(nsplit*nocc) ...
                1/(nsplit*nocc*nbar_i)];
        end
        design.reference_nsplit = nsplit;
        design.reference_nocc = nocc;
        design.pooled_harmonic_items = nbar_i;
    end

    %--- coefficients + variance components ---
    if issingle
        sig_p  = REL.out.sig_u(:);
        sig_s  = REL.out.sig_trl(:);
        sig_ps = REL.out.sig_splitxid(:);
        sig_e  = REL.out.sig_err(:);

        %reuse psyrat_rel_sing: total within-person SD folds person x split,
        %split main, and the n_i-scaled per-item residual; i = split main is the
        %term removed from relative error (Rocha Table 4 = exact for D, lower
        %bound for G).
        wp = sqrt(sig_ps.^2 + sig_s.^2 + sig_e.^2 ./ nbar_i);
        D = local_relsing(1,sig_p,wp,sig_s,nsplit,ciperc);
        G = local_relsing(2,sig_p,wp,sig_s,nsplit,ciperc);
        %error variances for the SEMs (same expressions psyrat_rel_sing uses
        %internally; see derivation note section 3)
        relvar = (sig_ps.^2 + sig_e.^2 ./ nbar_i) ./ nsplit;
        absvar = (sig_ps.^2 + sig_s.^2 + sig_e.^2 ./ nbar_i) ./ nsplit;
        coef = local_makecoef('Internal consistency (equivalence)',...
            D,G,relvar,absvar,ciperc);

        varcomp = local_varcomp_single(sig_p,sig_s,sig_ps,sig_e,ciperc);
    else
        sig_p   = psyrat_draw_column(REL.out.sig_id,s);
        sig_o   = psyrat_draw_column(REL.out.sig_occ,s);
        sig_s   = psyrat_draw_column(REL.out.sig_trl,s);
        sig_ps  = psyrat_draw_column(REL.out.sig_trlxid,s);
        sig_po  = psyrat_draw_column(REL.out.sig_occxid,s);
        sig_os  = psyrat_draw_column(REL.out.sig_trlxocc,s);
        sig_pos = psyrat_draw_column(REL.out.sig_posxid,s);
        sig_e   = psyrat_draw_column(REL.out.sig_err,s);

        %reuse psyrat_rel_trt: err folds the p:o:s 3-way and the n_i-scaled
        %per-item residual so the 3-way lands at /(n_s n_o) and the item residual
        %at /(n_s n_o n_i) for the complete balanced restricted model or a
        %balanced projection. Irregular incidence uses independent row weights.
        err = sqrt(sig_pos.^2 + sig_e.^2 ./ nbar_i);
        coef = struct('label',{},'D',{},'G',{},'SEM_abs',{},'SEM_rel',{},...
            'gbias',{});
        for r = 1:3
            if ~design.irregular || strcmp(design.mode,'balanced_projection')
                % Preserve the complete balanced and explicit projection
                % arithmetic, including the legacy sqrt/square evaluation order.
                D = local_reltrt(1,r,sig_p,sig_o,sig_s,sig_ps,sig_po,sig_os,err,...
                    nsplit,nocc,ciperc);
                G = local_reltrt(2,r,sig_p,sig_o,sig_s,sig_ps,sig_po,sig_os,err,...
                    nsplit,nocc,ciperc);
                [relvar,absvar] = local_trt_errvar(r,sig_o,sig_s,sig_ps,sig_po,...
                    sig_os,err,nsplit,nocc);
            else
                [u,relvar,absvar] = local_trt_incidence(r,design.mean_factors,...
                    sig_p,sig_o,sig_s,sig_ps,sig_po,sig_os,sig_pos,sig_e);
                D = local_draw_summary(u ./ (u + absvar),ciperc);
                G = local_draw_summary(u ./ (u + relvar),ciperc);
            end
            coef(r) = local_makecoef(reltype_label{r},D,G,relvar,absvar,ciperc);
        end

        for r = 1:3
            % False means no bound is asserted, not evidence of unbiased
            % estimation. The explicit interpretation fields are authoritative.
            coef(r).gbias = false;
            coef(r).relative_approximation = design.irregular && ...
                strcmp(design.mode,'observed_available_occasions');
            coef(r).relative_interpretation = ...
                'Conventional relative partition; no blanket bound or estimation-bias guarantee';
        end

        varcomp = local_varcomp_trt(sig_p,sig_o,sig_s,sig_ps,sig_po,sig_os,...
            sig_pos,sig_e,ciperc);
    end

    st = struct('label',labels{s},'nsplit',nsplit,'nocc',nocc,...
        'nbar_i',nbar_i,'coef',coef,'varcomp',varcomp);
    if ~issingle; st.design = design; end
    summary.strata(s) = st;
end

end

% ------------------------------------------------------------------------------

function D = local_relsing(gcoeff,bp,wp,ival,obs,ciperc)
%Thin wrapper returning [ll pt ul] from psyrat_rel_sing for the split-mean
%one-facet design (Rocha Table 4).
[ll,pt,ul] = psyrat_rel_sing('gcoeff',gcoeff,'metric','global',...
    'bp',bp,'wp',wp,'i',ival,'obs',obs,'CI',ciperc);
D = [ll pt ul];
end

function D = local_reltrt(gcoeff,reltype,bp,bo,bt,txp,oxp,txo,err,obs,nocc,ciperc)
%Thin wrapper returning [ll pt ul] from psyrat_rel_trt for the split-mean
%two-facet design (Rocha Table 5). bt = split main, txp = person x split,
%txo = occasion x split (the split facet maps onto psyrat_rel_trt's "trial").
[ll,pt,ul] = psyrat_rel_trt('gcoeff',gcoeff,'reltype',reltype,...
    'bp',bp,'bo',bo,'bt',bt,'txp',txp,'oxp',oxp,'txo',txo,'err',err,...
    'obs',obs,'nocc',nocc,'CI',ciperc);
D = [ll pt ul];
end

function [relvar,absvar] = local_trt_errvar(reltype,sig_o,sig_s,sig_ps,sig_po,...
    sig_os,err,obs,nocc)
%Per-draw relative/absolute error VARIANCE for the two-facet split design,
%matching psyrat_rel_trt's err_term expressions exactly (see
%psyrat_rel_trt.m:245-319). Used only for the SEMs; the coefficients come from
%psyrat_rel_trt itself, so the two share one definition of the error partition.
bo  = sig_o.^2;  bt  = sig_s.^2;  txp = sig_ps.^2;
oxp = sig_po.^2; txo = sig_os.^2; er  = err.^2;
switch reltype
    case 1 % coefficient of equivalence (occasion fixed)
        relvar = txp ./ obs + er ./ (obs*nocc);
        absvar = relvar + bt ./ obs + txo ./ (obs*nocc);
    case 2 % coefficient of stability (split fixed)
        relvar = oxp ./ nocc + er ./ (obs*nocc);
        absvar = relvar + bo ./ nocc + txo ./ (obs*nocc);
    case 3 % coefficient of equivalence and stability (both random)
        relvar = txp ./ obs + oxp ./ nocc + er ./ (obs*nocc);
        absvar = relvar + bt ./ obs + bo ./ nocc + txo ./ (obs*nocc);
end
end

function c = local_makecoef(label,D,G,relvar,absvar,ciperc)
%Assemble one reported coefficient: dependability D and generalizability G
%(each [ll pt ul] from the reuse wrappers) plus the absolute/relative SEMs
%(posterior sqrt of the error variance). G and the relative SEM are flagged
%biased (lower/upper bound) per the split-mean identifiability collapse.
ciedge = (1-ciperc)/2;
semabs = sqrt(absvar);
semrel = sqrt(relvar);
c.label   = label;
c.D       = D;
c.G       = G;
c.SEM_abs = [quantile(semabs,ciedge) mean(semabs) quantile(semabs,1-ciedge)];
c.SEM_rel = [quantile(semrel,ciedge) mean(semrel) quantile(semrel,1-ciedge)];
c.gbias   = true;
end

function T = local_varcomp_single(sig_p,sig_s,sig_ps,sig_e,ciperc)
%Variance-component table for the single-occasion split design (Rocha Table 4).
%The residual is the per-item value sig_err^2 (it enters the coefficients scaled
%by 1/n_i; see the derivation note).
[add,getT] = local_vc_accumulator(ciperc);
add('Between-person variance (universe score)','sigma^2_p', sig_p.^2);
add('Split main-effect variance','sigma^2_s', sig_s.^2);
add('Person x split variance','sigma^2_ps', sig_ps.^2);
add('Per-item residual variance (item:split lumped)','sigma^2_(i:s),e', sig_e.^2);
T = getT();
end

function T = local_varcomp_trt(sig_p,sig_o,sig_s,sig_ps,sig_po,sig_os,sig_pos,...
    sig_e,ciperc)
%Variance-component table for the fitted restricted multi-occasion split model.
[add,getT] = local_vc_accumulator(ciperc);
add('Between-person variance (universe score)','sigma^2_p', sig_p.^2);
add('Occasion main-effect variance','sigma^2_o', sig_o.^2);
add('Split main-effect variance','sigma^2_s', sig_s.^2);
add('Person x split variance','sigma^2_ps', sig_ps.^2);
add('Person x occasion variance','sigma^2_po', sig_po.^2);
add('Occasion x split variance','sigma^2_os', sig_os.^2);
add('Person x occasion x split variance','sigma^2_pos', sig_pos.^2);
add('Per-item residual variance (item:split lumped)','sigma^2_po(i:s),e', sig_e.^2);
T = getT();
end

function [add,getT] = local_vc_accumulator(ciperc)
%Shared accumulator for the variance-component tables (same idiom as
%psyrat_dynrel_summary): add(name,symbol,draws) appends one row -- the posterior
%mean plus the central credible interval -- and getT() returns the assembled
%five-column table {component,symbol,estimate,ci_lower,ci_upper}.
ciedge = (1-ciperc)/2;
name = {}; sym = {}; ptv = []; llv = []; ulv = [];
    function add_(n,sy,draws)
        name{end+1,1} = n; sym{end+1,1} = sy;
        ptv(end+1,1) = mean(draws);
        llv(end+1,1) = quantile(draws,ciedge);
        ulv(end+1,1) = quantile(draws,1-ciedge);
    end
    function T = getT_()
        T = table(name,sym,ptv,llv,ulv,'VariableNames',...
            {'component','symbol','estimate','ci_lower','ci_upper'});
    end
add = @add_;
getT = @getT_;
end

function nocc = local_resolve_nocc(nocc_in,nocc_obs)
%Resolve the occasion n' for the multi-occasion split summary. [] or 'observed'
%-> the stratum's observed number of occasions (the observed-design default; no
%projection); a numeric scalar overrides. Always at least 1.
if isempty(nocc_in)
    nocc = nocc_obs;
elseif (ischar(nocc_in) || isstring(nocc_in)) && strcmpi(char(nocc_in),'observed')
    nocc = nocc_obs;
elseif isnumeric(nocc_in) && isscalar(nocc_in)
    nocc = nocc_in;
else
    nocc = nocc_obs;
end
nocc = max(1,round(nocc));
end

function [rows,group_label,event_label] = local_stratum_rows(REL,label,issingle)
%Logical index into REL.data for one stratum. The single-occasion design is one
%stratum (all rows). The multi-occasion design is chunked per group/event with
%labels 'group_;_event' (group x event), or the bare group/event name, or 'none';
%this reconstructs the subset by matching REL.data.group / REL.data.event.
n = height(REL.data);
group_label = 'none'; event_label = 'none';
if issingle
    rows = true(n,1);
    return;
end
hasg = ~local_isnone(REL.groups);
hase = ~local_isnone(REL.events);
if ~hasg && ~hase
    rows = true(n,1);
elseif hasg && ~hase
    group_label = char(string(label));
    rows = strcmp(string(REL.data.group),string(label));
elseif ~hasg && hase
    event_label = char(string(label));
    rows = strcmp(string(REL.data.event),string(label));
else
    % Enumerate true tuples: a delimiter in a label is valid when the
    % joined display label still identifies exactly one observed tuple.
    tuples = unique([string(REL.data.group) string(REL.data.event)],'rows');
    matches = find(tuples(:,1) + "_;_" + tuples(:,2) == string(label));
    if numel(matches) > 1
        error('psyrat_splits_summary:ambiguousStratum',...
            ['Stratum ''%s'' identifies multiple group/event tuples. '...
            'Use unambiguous group/event labels before fitting.'],label);
    elseif isempty(matches)
        rows = false(n,1);
    else
        group_label = char(tuples(matches,1));
        event_label = char(tuples(matches,2));
        rows = string(REL.data.group) == tuples(matches,1) & ...
            string(REL.data.event) == tuples(matches,2);
    end
end
end

function tf = local_isnone(x)
%True when a REL.groups/REL.events field encodes "no facet" ('none'). The
%engine's facet-absent sentinel is ALWAYS a char/string SCALAR; a cell carries
%real labels by construction, so it is never the sentinel. Unwrapping a cell
%here (B26) misread label sets whose alphabetically-first entry is literally
%'none' ({'none','patients'}) as "facet absent" and pooled every stratum.
tf = (ischar(x) || (isstring(x) && isscalar(x))) && strcmpi(char(x),'none');
end

function d = local_trt_design(ids,times,w,row_indices,group_label,event_label)
% Reconstruct fitted shared ordinal split identity in each TRUE tuple cell.
% Numeric codes avoid delimiter collisions; find preserves original row order
% even when a hand-built REL interleaves a cell's rows.
[person_ids,~,pg] = unique(ids);
[occasion_ids,~,og] = unique(times);
[cells,~,cg] = unique([pg og],'rows');
counts = accumarray(cg,1);
np = numel(person_ids);
occasions = accumarray(cells(:,1),1,[np 1]);
a = 1 ./ (occasions(pg) .* counts(cg));
ordinal = zeros(size(w));
for k = 1:numel(counts)
    ri = find(cg == k);
    ordinal(ri) = (1:numel(ri))';
end
factors = zeros(np,4);
for p = 1:np
    ri = pg == p;
    split_weights = accumarray(ordinal(ri),a(ri));
    occasion_weights = accumarray(og(ri),a(ri));
    factors(p,:) = [sum(split_weights.^2) sum(occasion_weights.^2) ...
        sum(a(ri).^2) sum(a(ri).^2 ./ w(ri))];
end
d.group = group_label;
d.event = event_label;
d.participant_count = np;
d.person_ids = string(person_ids);
d.occasion_ids = string(occasion_ids);
d.occasion_union = numel(occasion_ids);
d.person_occasion_counts = occasions;
d.cell_person_index = cells(:,1);
d.cell_occasion_index = cells(:,2);
d.cell_split_counts = counts;
d.row_indices = row_indices;
d.row_person_index = pg;
d.row_occasion_index = og;
d.row_cell_index = cg;
d.row_split_index = ordinal;
d.row_items = w;
d.row_weights = a;
d.factor_order = {'A','B','C','W'};
d.person_factors = factors;
d.mean_factors = mean(factors,1);
d.missing_occasions = any(occasions ~= d.occasion_union);
d.irregular = d.missing_occasions || any(counts ~= counts(1));
d.aggregation = 'Equal-person mean universe/error variances before each draw ratio and SEM';
d.split_identity = 'Shared ordinal within true person/occasion cell in original relative row order';
d.row_residual_assumption = 'Independent split-row residual with variance sig_err^2 / items';
d.fixed_facet_target = 'CE/CS retain person-specific fixed-facet support; missingness is not modeled';
end

function [u,relvar,absvar] = local_trt_incidence(r,f,p,o,s,ps,po,os,pos,e)
% Linearity permits averaging independent effect-incidence factors first.
% This is equal-person U/error averaging, not an average of coefficients.
A=f(1); B=f(2); C=f(3); W=f(4);
u = p.^2;
rowvar = pos.^2 .* C + e.^2 .* W;
switch r
    case 1
        u = u + po.^2 .* B;
        relvar = ps.^2 .* A + rowvar;
        absvar = relvar + s.^2 .* A + os.^2 .* C;
    case 2
        u = u + ps.^2 .* A;
        relvar = po.^2 .* B + rowvar;
        absvar = relvar + o.^2 .* B + os.^2 .* C;
    case 3
        relvar = ps.^2 .* A + po.^2 .* B + rowvar;
        absvar = relvar + s.^2 .* A + o.^2 .* B + os.^2 .* C;
end
end

function v = local_draw_summary(draws,ciperc)
ciedge = (1-ciperc)/2;
v = [quantile(draws,ciedge) mean(draws) quantile(draws,1-ciedge)];
end
