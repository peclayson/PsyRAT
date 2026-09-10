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
%OBSERVED design (n_s splits, n_o occasions, and the harmonic-mean items-per-split
%n_i). There is no D-study projection range and no cut-score search: split means
%lose the item-level partition, so a projection over hypothetical item counts is
%not identified. The split designs deliberately bypass psyrat_relsummary's
%trial-cutoff machinery (it does not fit an observed-design split summary), the
%same way psyrat_dynrel_summary handles the dynamic-reliability family.
%
%The full derivation and the component->coefficient mapping are recorded in
%documentation/splits_observed_design_formulas.md (the Golden-Rule-1 traceability
%record). In brief, from split means the per-item residual sig_err lumps the item
%MAIN effect with the person-involving item residual, so:
%   - dependability (absolute, D) is computed EXACTLY, and
%   - generalizability (relative, G) is a downward-biased observed-design LOWER
%     BOUND (and the relative SEM a corresponding upper bound).
%Both are reported; every relative coefficient carries a gbias flag.
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
% obs  - n_s (splits per person) for the coefficients. Scalar applied to every
%        stratum, or [] (default) to use each stratum's observed median split
%        count per person (per person x occasion for trt_splits).
% nocc - n_o (occasions) for the multi-occasion coefficients. [] (default) or
%        'observed' uses each stratum's observed number of occasions; a numeric
%        scalar overrides. Inert for ic_splits.
%
%Output (struct)
% summary.analysis, summary.ci, summary.gbias_note
% summary.strata(s) for each stratum s:
%   .label   - stratum name ('measure' for the single-occasion design)
%   .nsplit  - n_s used
%   .nocc    - n_o used (1 for ic_splits)
%   .nbar_i  - harmonic-mean items-per-split used to scale the per-item residual
%   .coef    - struct array of reported coefficients (1 for ic_splits; the three
%              reltypes CE/CS/CES for trt_splits). Each has .label, .D, .G,
%              .SEM_abs, .SEM_rel (each [ll pt ul]) and .gbias (true).
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
    'correspondingly an upper bound. Dependability (absolute) is exact.'];
summary.strata = struct('label',{},'nsplit',{},'nocc',{},'nbar_i',{},...
    'coef',{},'varcomp',{});

for s = 1:nstrata
    %--- observed design for this stratum (recomputed from REL.data; the
    %estimation stores no n_s/n_i, and recomputing keeps that committed core
    %untouched) ---
    rows = local_stratum_rows(REL,labels{s},issingle);
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
    %derivation note section 5).
    nbar_i = numel(w) / sum(1 ./ w);

    idsub = REL.data.id(rows);
    if issingle
        nsplit_obs = median(local_cellcounts(idsub));   % splits per person
        nocc_obs = 1;
    else
        timesub = REL.data.time(rows);
        %splits per person x occasion cell (the case-12 lesson: a per-cell n')
        nsplit_obs = median(local_cellcounts(...
            strcat(string(idsub),'_;_',string(timesub))));
        nocc_obs = numel(unique(timesub));
    end

    if isempty(obs_in); nsplit = nsplit_obs; else; nsplit = obs_in; end
    nsplit = max(1,round(nsplit));
    nocc = local_resolve_nocc(nocc_in,nocc_obs);

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
        sig_p   = REL.out.sig_id(:,s);
        sig_o   = REL.out.sig_occ(:,s);
        sig_s   = REL.out.sig_trl(:,s);
        sig_ps  = REL.out.sig_trlxid(:,s);
        sig_po  = REL.out.sig_occxid(:,s);
        sig_os  = REL.out.sig_trlxocc(:,s);
        sig_pos = REL.out.sig_posxid(:,s);
        sig_e   = REL.out.sig_err(:,s);

        %reuse psyrat_rel_trt: err folds the p:o:s 3-way and the n_i-scaled
        %per-item residual so the 3-way lands at /(n_s n_o) and the item residual
        %at /(n_s n_o n_i) (Rocha Table 5 under the nested-items assumption; see
        %derivation note section 4).
        err = sqrt(sig_pos.^2 + sig_e.^2 ./ nbar_i);
        coef = struct('label',{},'D',{},'G',{},'SEM_abs',{},'SEM_rel',{},...
            'gbias',{});
        for r = 1:3
            D = local_reltrt(1,r,sig_p,sig_o,sig_s,sig_ps,sig_po,sig_os,err,...
                nsplit,nocc,ciperc);
            G = local_reltrt(2,r,sig_p,sig_o,sig_s,sig_ps,sig_po,sig_os,err,...
                nsplit,nocc,ciperc);
            [relvar,absvar] = local_trt_errvar(r,sig_o,sig_s,sig_ps,sig_po,...
                sig_os,err,nsplit,nocc);
            coef(r) = local_makecoef(reltype_label{r},D,G,relvar,absvar,ciperc);
        end

        varcomp = local_varcomp_trt(sig_p,sig_o,sig_s,sig_ps,sig_po,sig_os,...
            sig_pos,sig_e,ciperc);
    end

    st = struct('label',labels{s},'nsplit',nsplit,'nocc',nocc,...
        'nbar_i',nbar_i,'coef',coef,'varcomp',varcomp);
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
%Variance-component table for the multi-occasion split design (Rocha Table 5,
%estimated nested-items component set).
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

function rows = local_stratum_rows(REL,label,issingle)
%Logical index into REL.data for one stratum. The single-occasion design is one
%stratum (all rows). The multi-occasion design is chunked per group/event with
%labels 'group_;_event' (group x event), or the bare group/event name, or 'none';
%this reconstructs the subset by matching REL.data.group / REL.data.event.
n = height(REL.data);
if issingle
    rows = true(n,1);
    return;
end
hasg = ~local_isnone(REL.groups);
hase = ~local_isnone(REL.events);
if ~hasg && ~hase
    rows = true(n,1);
elseif hasg && ~hase
    rows = strcmp(string(REL.data.group),string(label));
elseif ~hasg && hase
    rows = strcmp(string(REL.data.event),string(label));
else
    parts = split(string(label),'_;_');   % [group, event]
    rows = strcmp(string(REL.data.group),parts(1)) & ...
        strcmp(string(REL.data.event),parts(2));
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

function counts = local_cellcounts(key)
%Number of rows per unique key value (splits per person, or per person x
%occasion when key encodes both).
[~,~,g] = unique(key);
counts = accumarray(g,1);
end
