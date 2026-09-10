function summary = psyrat_dynrel_summary(REL,varargin)
%Summarize a dynamic/conditional reliability run (Rast & Clayson, analysis 11)
%into per-stratum reliability surfaces and variance components.
%
% summary = psyrat_dynrel_summary(REL,'CI',.95,'ngrid',25,'obs',[],...
%   'marginal',0)
%
%This is the single place that turns the per-group/event posterior draws stored
%by psyrat_computevarcomp (case 11) into the reliability surface G(z)/D(z) and
%the reported variance components. The figure (psyrat_dynrelplot), the variance
%table, and the GUI viewer all consume this struct, so the surface is computed
%once and consistently. The dynamic-reliability paradigm (a continuous surface
%over the standardized dimensions) does not fit psyrat_relsummary's trial-cutoff
%machinery, so it is summarized here instead of bloating that fragile core.
%
%Inputs
% REL - the result struct from psyrat_computevarcomp with REL.analysis ==
%  'ic_dynrel' (REL.out carries gro_sds, pop_sdlog, b_sigma, b, sig_trl,
%  chol_corrmat, dimz, ntrials, labels per stratum; REL.ndim / REL.dim_names).
%
%Optional Inputs
% CI       - credible-interval width in (0,1], default .95.
% ngrid    - number of grid points per dimension axis, default 25.
% obs      - n' (trial count) for the surface. Scalar applied to every stratum,
%            or [] (default) to use each stratum's median trial count
%            (REL.out.ntrials).
% marginal - 0 (default): typical person, delta_p = 0. 1: lognormal
%            population-average over the person scale random effect.
%
%Output (struct)
% summary.ndim, summary.dim_names, summary.ci, summary.marginal
% summary.isgamma, summary.isgamma_ls - family/parameterization flags for the
%            viewer and headless report (estimand labeling)
% summary.copula_labels - modular cut-copula provenance labels (any design
%            storing REL.out.copula; copied from the stage output via
%            psyrat_copula_labels; [] for every other design). Programmatic
%            carrier only - the printed disclosure comes from
%            psyrat_provenance_lines reading the same fields.
% summary.dod_labels - nonconcurrent DoD estimand labels (any design storing
%            REL.out.dod_labels, via psyrat_dod_labels; [] otherwise). Same
%            carrier doctrine as copula_labels; also keys the viewer's
%            typical-person estimand line.
% summary.strata(s) for each group/event stratum s:
%   .label                 - stratum name
%   .obs                   - n' used
%   .obs_ev                - per-event/per-cell n': [n1 n2] for the gamma
%                            difference designs, a 1x4 per-cell vector for the
%                            DoD designs 28/29 (both families); [] otherwise
%   .nocc                  - occasion n' (two-facet variants only; [])
%   .z1 (, .z2)            - grid axes (observed range)
%   .G, .D, .ICCg, .ICCd  - surfaces (ll/pt/ul) from psyrat_rel_dynrel
%   .sige_pt              - relative residual SD over the grid
%   .varcomp              - table of reported variance components (Greek labels)
%   .ssrel_table          - per-participant table (subject-level variants; [])
%   .quaddiag             - per-participant quadrature provenance (modular
%                           cut-copula designs only; [])
%   .copula_flags         - stage-2 diagnostic counts, n_draws/n_invalid/
%                           n_calibration (modular cut-copula designs only; [])
%   .dod_flags            - signed-projection invalidity counts, n_evals/
%                           n_invalid/n_invalid_std (nonconcurrent DoD designs
%                           only; [])

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

ciperc   = psyrat_opt(varargin,'CI',.95);
ngrid    = psyrat_opt(varargin,'ngrid',25);
obs_in   = psyrat_opt(varargin,'obs',[]);
marginal = psyrat_opt(varargin,'marginal',0);
%occasion n' and reltype are only consumed by the trial + occasion two-facet
%variant (ic_dynrel_trt); inert for the one-facet / difference analyses. nocc_in
%[] -> default 1 (owner decision); 'observed' -> each stratum's number of
%occasions. reltype 1/2/3 selects the equivalence / stability / both-random
%coefficient (see psyrat_rel_trt).
nocc_in  = psyrat_opt(varargin,'nocc',[]);
reltype  = psyrat_opt(varargin,'reltype',3);

%accepted-set guard: psyrat_analysis_family is the single source of truth for
%which analysis strings are dynamic-reliability variants (shared with the GUI
%router and the headless report dispatch).
if ~isfield(REL,'analysis') || ...
        ~strcmp(psyrat_analysis_family(REL.analysis),'dynrel')
    error('dynrel:badinput',...
        ['psyrat_dynrel_summary expects a dynamic-reliability result '...
        '(REL.analysis = ''ic_dynrel'', ''ic_dynrel_sserrvar'', '...
        '''ic_diff_dynrel'', ''ic_diff_dynrel_sserrvar'', ''ic_dynrel_trt'', '...
        '''ic_dynrel_sserrvar_trt'', '...
        '''ic_diff_dynrel_trt'', ''ic_diff_dynrel_sserrvar_trt'', '...
        '''ic_diff_dynrel_rescor'', ''ic_diff_dynrel_trt_rescor'', '...
        '''ic_diff_dynrel_sserrvar_rescor'', '...
        '''ic_diff_dynrel_sserrvar_trt_rescor'', '...
        '''ic_diff_dynrel_sserrvar_rho'', '...
        '''ic_diff_dynrel_sserrvar_trt_rho'', '...
        '''ic_dodiff_dynrel_sserrvar'', or '...
        '''ic_dodiff_dynrel_sserrvar_trt'').']);
end

%the two-event difference-score variants have a different REL.out layout
%(id_varcov/trl_varcov/b_sigma per event/b_sigma_dim) and compute the surface
%with psyrat_rel_diffdynrel. The subject-level variant (issubj) ADDS a
%per-participant table (each subject's phi at their own z, via psyrat_ssrel_diff)
%while still reporting the typical-person (r_sigma = 0) surface as a population
%reference. All variants produce the same uniform summary struct so the figure
%and viewer consume them identically.
issubj = strcmp(REL.analysis,'ic_diff_dynrel_sserrvar');
%the one-facet group-level CONCURRENT difference variant (analysis 17) reuses the
%same psyrat_rel_diffdynrel surface but feeds an estimated residual correlation
%(REL.out.rescor) so the z-conditioned residual covariance enters the difference
%variance. Group-level, so no per-participant table.
isdiffrescor = strcmp(REL.analysis,'ic_diff_dynrel_rescor');
%the one-facet SUBJECT-LEVEL CONCURRENT variant (analysis 19) = case 13's joint
%person block + the population residual correlation: it has BOTH a per-participant
%table (per-person scale RE) AND the rescor term. Reuses the isdiff path.
issubjrescor = strcmp(REL.analysis,'ic_diff_dynrel_sserrvar_rescor');
%the one-facet SUBJECT-LEVEL CONCURRENT variant with a PER-SUBJECT residual
%correlation (analysis 21, Stage 4e) = analysis 19 with the shared population
%correlation replaced by a per-person rho. It reuses the isdiff path exactly like
%issubjrescor, but the per-participant table feeds each subject's own rho (a
%draws x NSUB matrix in REL.out.rho_ss) instead of the shared scalar, and the
%reference surface uses the population-average correlation (REL.out.rescor =
%mean(rho)).
issubjrho = strcmp(REL.analysis,'ic_diff_dynrel_sserrvar_rho');
isdiff = strcmp(REL.analysis,'ic_diff_dynrel') || issubj || isdiffrescor ...
    || issubjrescor || issubjrho;
%the ONE-FACET non-difference SUBJECT-LEVEL variant (analysis 26,
%ic_dynrel_sserrvar) is the same one-facet location-scale model as ic_dynrel
%(case 11), which already estimates a per-subject residual scale (ind_sdlog); it
%ADDS a per-participant table (each subject's phi_s(z) at their own z, via
%psyrat_ssrel - the same per-subject calc the static ic_sserrvar uses) while
%still reporting the typical-person population surface as a reference. It reuses
%the one-facet (else) surface path below.
issubjnondiff = strcmp(REL.analysis,'ic_dynrel_sserrvar');
%the TWO-FACET non-difference SUBJECT-LEVEL variant (analysis 27,
%ic_dynrel_sserrvar_trt) is the two-facet analogue: the same model as
%ic_dynrel_trt (case 14) with ind_sdlog additionally extracted. It reuses the
%istrt surface path plus a per-participant table via psyrat_ssrel_trt.
issubjnondiff_trt = strcmp(REL.analysis,'ic_dynrel_sserrvar_trt');
%the trial + occasion two-facet non-difference variant (analysis 14) reuses the
%test-retest psyrat_rel_trt via psyrat_rel_dynrel_trt and reports the occasion
%facet alongside the trial facet (its own REL.out layout: gro_sds/chol_corrmat +
%the five crossed facet SDs + nocc). The subject-level two-facet variant
%(issubjnondiff_trt) reuses this same surface path, adding a per-participant table.
istrt  = strcmp(REL.analysis,'ic_dynrel_trt') || issubjnondiff_trt;
%the trial + occasion two-facet GROUP-LEVEL difference variant (analysis 15)
%reuses the two-facet difference psyrat_diffrel_trt via psyrat_rel_diffdynrel_trt
%(its REL.out layout: the six event-specific 2x2 covariance blocks + per-event
%b_sigma/b_sigma_dim/b_dim + nocc). Fixed-per-event residual, so no per-subject
%table and no marginal estimand (as in Stage 2a).
%the trial + occasion two-facet group-level CONCURRENT difference variant
%(analysis 18) reuses the same psyrat_rel_diffdynrel_trt surface but feeds the
%estimated residual correlation (REL.out.rescor). Group-level, so no per-subject
%table.
isdifftrtrescor = strcmp(REL.analysis,'ic_diff_dynrel_trt_rescor');
%the two-facet SUBJECT-LEVEL CONCURRENT variant (analysis 20) = case 16's joint
%person block + the population residual correlation. Reuses the isdifftrt path.
issubjtrtrescor = strcmp(REL.analysis,'ic_diff_dynrel_sserrvar_trt_rescor');
%the two-facet SUBJECT-LEVEL CONCURRENT variant with a PER-SUBJECT residual
%correlation (analysis 22, Stage 4e) = analysis 20 with the per-person rho. Reuses
%the isdifftrt path like issubjtrtrescor; the per-participant table feeds each
%subject's own rho (REL.out.rho_ss).
issubjtrtrho = strcmp(REL.analysis,'ic_diff_dynrel_sserrvar_trt_rho');
isdifftrt = strcmp(REL.analysis,'ic_diff_dynrel_trt') || isdifftrtrescor ...
    || issubjtrtrescor || issubjtrtrho;
%the trial + occasion two-facet SUBJECT-LEVEL difference variant (analysis 16) is
%the group-level difttrt model with the case-13 joint 4-coefficient person block,
%so it ADDS a per-participant table (each subject's phi at their own z, via
%psyrat_ssrel_diff_trt) and the person scale-RE components, while still reporting
%the typical-person (scale RE = 0) reference surface. Its REL.out layout extends
%isdifftrt with sd_id(4) / L_id(4x4) / er_var_ss1 / er_var_ss2 / ssinfo. It reuses
%the isdifftrt surface + varcomp path below.
issubjtrt = strcmp(REL.analysis,'ic_diff_dynrel_sserrvar_trt');
%the person-specific dynamic NONCONCURRENT difference-of-differences (analyses
%28 one-facet / 29 two-facet): four margins fitted jointly, all six
%observation-level residual covariances fixed to zero by the estimand. Two
%families exist - gamma location-scale and, since the 2026-08-20 owner ruling,
%Gaussian (SCIENTIFIC_FORMULA_AUDIT.md section 26) - and the engine guards
%guarantee no other stamp can be produced; the branch below still fails loudly
%by name on any other arrival rather than misrouting - the exact
%silent-fall-through failure the case-13 note in the isdiff branch records.
isdod = strcmp(REL.analysis,'ic_dodiff_dynrel_sserrvar');
isdodtrt = strcmp(REL.analysis,'ic_dodiff_dynrel_sserrvar_trt');

%the population-average ('marginal') estimand applies only to the non-difference
%one-facet surface (psyrat_rel_dynrel). The difference variants use a fixed
%per-event residual (psyrat_rel_diffdynrel/_trt take no marginal argument), so a
%marginal request here cannot be honored; force it off and warn rather than
%silently advertising summary.marginal = 1 over a typical-person surface. (The
%GUI already disables the estimand control for difference analyses; this guards
%direct programmatic callers.)
if marginal == 1 && (isdiff || isdifftrt)
    warning('dynrel:marginalnotsupported',...
        ['The population-average (marginal) estimand is not defined for '...
        'difference-score dynamic reliability; using the typical-person '...
        'surface (marginal = 0).']);
    marginal = 0;
end

%the gamma / scaled-chi-square dynrel surfaces carry no Gaussian delta_p, so
%the typical-person / population-average toggle does not apply to any of them;
%force it off (and warn if explicitly requested) rather than record a
%misleading summary.marginal. WHICH estimand each surface reports differs by
%design and is labeled downstream: the non-difference and two-event designs
%report the marginal-over-person-scale surface (estimand #2; log-nu under
%gammascale = 1, residual-SD heterogeneity under 2), while the nonconcurrent
%DoD designs (28/29) report the TYPICAL-PERSON surface (all person effects at
%zero; see psyrat_rel_dodiffdynrel_gamma_ls) - do not restate the marginal
%claim as a universal.
isgamma = isfield(REL,'family') && strcmp(REL.family,'gamma');
%the LOCATION-SCALE gamma variant (gammascale = 2) replaces the log-nu scale
%submodel with a log residual SD, so it reads DIFFERENT REL.out slots and
%converts through a different function. The family alone does not identify it.
%
%Gate on the family AND gammascale, in that order and never on gammascale
%alone: psyrat_computevarcomp stamps REL.gammascale on EVERY run and defaults it
%to 1, so a bare check would fire on Gaussian analyses too. isfield is still
%required because a .psyrat saved before the option existed carries no such
%field at all, and must fall through to the log-nu path.
isgamma_ls = isgamma && isfield(REL,'gammascale') && isequal(REL.gammascale,2);
if marginal == 1 && isgamma
    warning('dynrel:marginalgamma',...
        ['The population-average (marginal) toggle does not apply to the gamma / '...
        'scaled-chi-square dynamic reliability. The non-difference and '...
        'two-event designs report the estimand-#2 (marginal-over-dispersion) '...
        'surface; the nonconcurrent difference-of-differences designs report '...
        'the typical-person surface (all person effects at zero).']);
    marginal = 0;
end
%the nonconcurrent DoD designs report the typical-person surface in EVERY
%family (owner rulings 2026-08-12 gamma / 2026-08-20 Gaussian), so the
%marginal toggle is not honored under the GAUSSIAN arm either - which would
%otherwise slip past the gamma-only force above.
if marginal == 1 && (isdod || isdodtrt)
    warning('dynrel:marginalnotsupported',...
        ['The population-average (marginal) estimand is not defined for '...
        'the nonconcurrent difference-of-differences; using the '...
        'typical-person surface (marginal = 0).']);
    marginal = 0;
end

ndim = REL.ndim;
nstrata = numel(REL.out.labels);

summary = struct();
summary.ndim = ndim;
summary.dim_names = REL.dim_names;
summary.ci = ciperc;
summary.marginal = marginal;
%isgamma is carried on the summary so the viewer and headless report can label the
%estimand correctly. The Gaussian marginal flag does not describe either gamma
%surface: neither parameterization has a delta_p to set to zero, so the surface is
%ALWAYS the marginal one regardless of summary.marginal, which is forced to 0
%above. Without this flag both exports print "typical person (delta_p = 0)" on a
%gamma run, which is wrong in both halves - and doubly misleading on the
%subject-level per-participant table, which is CONDITIONAL on each person's own
%scale effect.
summary.isgamma = isgamma;
%isgamma_ls rides along for the same reason: the viewer and the headless report
%name the per-participant table's conditioning set, and under the location-scale
%parameterization that set is each person's own log RESIDUAL SD - not the
%log-mean and log-nu the log-nu branch conditions on. The participant's own
%log-mean does not enter the coefficient at all here (see
%psyrat_ssrel_dynrel_gamma_ls), so carrying the log-nu wording over would be
%wrong in both halves.
summary.isgamma_ls = isgamma_ls;
%The modular cut-copula provenance labels ride on the summary as the
%PROGRAMMATIC carrier for API consumers and tests; the PRINTED disclosure the
%user sees (report + export headers + viewer note) is emitted by
%psyrat_provenance_lines from the same underlying REL.out.copula fields through
%the same shared helpers, so the two carriers cannot drift. Keyed on the
%presence of out.copula - the field only the cut-copula arms write - rather
%than on exact analysis strings, so a future design built on the same stage
%(the barred 21/22) inherits the disclosure instead of silently missing it.
%[] for every other design. The stage (psyrat_gamma_cut_copula_stage) stays the
%single source of truth for the strings.
summary.copula_labels = [];
if isfield(REL.out,'copula')
    summary.copula_labels = psyrat_copula_labels(REL.out.copula);
end
%The nonconcurrent DoD estimand labels ride on the summary the same way:
%PRESENCE-keyed on REL.out.dod_labels - the field only the analysis-28/29 arms
%write - never on analysis strings, and read through the shared helper so this
%carrier and the printed provenance cannot drift. [] for every other design.
summary.dod_labels = psyrat_dod_labels(REL.out);
%obs_ev carries the vector n': PER-EVENT [n1 n2] for the gamma difference
%designs, and a PER-CELL 1x4 vector for the DoD designs 28/29 in both
%families (resolved below). Elsewhere it is populated only by the gamma
%difference variant, whose designs are routinely unbalanced (e.g. error vs correct
%trials). obs stays the scalar reference count for every variant so the plot label
%and the variance-component table keep their existing contract.
%The field list here must match the st struct assembled per stratum EXACTLY, in
%name and in order: `summary.strata(s) = st` is a struct-array assignment and
%rejects a mismatch outright. quaddiag and copula_flags are empty for every
%design except the modular cut-copula ones.
summary.strata = struct('label',{},'obs',{},'obs_ev',{},'nocc',{},'z1',{},'z2',{},...
    'G',{},'D',{},'ICCg',{},'ICCd',{},'sige_pt',{},'varcomp',{},...
    'ssrel_table',{},'quaddiag',{},'copula_flags',{},'dod_flags',{});

for s = 1:nstrata
    %observed standardized dimension values for this stratum (subject-level)
    zsub = REL.out.dimz{s};
    if isvector(zsub); zsub = zsub(:); end

    %n' for the surface. The DoD designs store PER-CELL median counts as a
    %4 x nstrata matrix, so their scalar reference is the column mean; the
    %per-cell vector itself is resolved just below.
    if isempty(obs_in)
        if isdod || isdodtrt
            obs = mean(REL.out.ntrials(:,s));
        else
            obs = REL.out.ntrials(s);
        end
    else
        obs = obs_in;
    end
    %reject a non-finite n' BEFORE the max(1,...) floor: MATLAB's two-input max
    %omits NaN, so max(1,NaN) silently returns 1 and the surface would be
    %reported at n' = 1 with no diagnostic (and the downstream validators, which
    %only see the floored value, could never fire).
    if any(~isfinite(obs))
        error('dynrel:obsnonfinite',...
            ['WARNING: the requested n'' (obs) is not finite.\nHow to fix: pass a '...
            'finite positive trial count, or omit ''obs'' to use the observed '...
            'per-stratum count.\n']);
    end
    obs = max(1,round(obs));

    %PER-EVENT n' for the gamma difference variant. The scalar above is the
    %median per-event count derived by dividing each subject's total rows by
    %nevent, which is the true per-event n' only under balanced replication;
    %difference designs frequently are not (the reference ERP design is 49 error
    %vs 341 correct trials). psyrat_diffrel already scales each event's variance
    %by its own n', so pass the exact counts. A user-supplied scalar 'obs'
    %override is broadcast to both events; a 2-element override is used as given.
    obs_ev = [];
    if isdiff && isgamma
        if ~isempty(obs_in)
            obs_ev = max(1,round(obs_in));
            if isscalar(obs_ev); obs_ev = [obs_ev obs_ev]; end
        elseif isfield(REL.out,'ntrials_ev') && ~isempty(REL.out.ntrials_ev)
            obs_ev = max(1,round(REL.out.ntrials_ev(:,s)'));
        else
            obs_ev = [obs obs];
        end
        %keep the scalar reference count consistent with the per-event pair the
        %surface actually used, so the plot label and varcomp table do not claim
        %a count no event has.
        obs = round(mean(obs_ev));
    end

    %PER-CELL n' for the DoD designs: a 1 x 4 vector in dodmap constituent
    %order. A user-supplied scalar 'obs' is broadcast to all four cells; a
    %4-element override is used as given; anything else is rejected rather
    %than truncated.
    if isdod || isdodtrt
        if ~isempty(obs_in)
            obs_ev = max(1,round(obs_in(:)'));
            if isscalar(obs_ev); obs_ev = repmat(obs_ev,1,4); end
            if numel(obs_ev) ~= 4
                error('dynrel:dodobs',...
                    ['WARNING: the DoD n'' override must be a scalar or a '...
                    '4-element per-cell vector; got %d elements.\n'],...
                    numel(obs_in));
            end
        else
            obs_ev = max(1,round(REL.out.ntrials(:,s)'));
        end
        obs = round(mean(obs_ev));
    end

    %grid axes spanning the observed standardized range. local_axis guards
    %against a stratum whose subjects share one dimension value (min == max),
    %which would otherwise collapse the axis to a single point and render a
    %degenerate plot; it also warns, since the scale slope is not identified
    %within a stratum that has no dimension variation. Shared by both variants.
    z1 = local_axis(zsub(:,1),ngrid,REL.dim_names{1},REL.out.labels{s});
    if ndim == 2
        z2 = local_axis(zsub(:,2),ngrid,REL.dim_names{2},REL.out.labels{s});
    end

    ssrel = [];
    %per-participant quadrature provenance; only the modular cut-copula
    %designs fill it, and it stays empty everywhere else.
    quaddiag = [];
    %stage-2 diagnostic flag counts (invalidity/calibration); only the modular
    %cut-copula designs fill it, and it stays empty everywhere else.
    cutflags = [];
    %signed-projection invalidity counts; only the nonconcurrent DoD designs
    %fill it, and it stays empty everywhere else.
    dodflags = [];
    nocc_used = []; %only the trial + occasion two-facet variant sets this
    if isdod || isdodtrt
        %person-specific dynamic NONCONCURRENT difference-of-differences
        %(analyses 28/29), two families: gamma LOCATION-SCALE and Gaussian.
        %Fail loudly by name on any other arrival - a log-nu gamma stamp, a
        %missing family field (a damaged or hand-edited REL), anything else:
        %the engine cannot produce one, but a bad stamp must not misroute
        %into a two-event branch (the documented case-13 silent-fall-through
        %failure mode). Family FIRST, never gammascale alone: the engine
        %stamps gammascale on EVERY run and the Gaussian default leaves it
        %at 1 (the isgamma_ls comment above documents the trap).
        isgaussdod = isfield(REL,'family') && strcmp(REL.family,'gaussian');
        if ~isgamma_ls && ~isgaussdod
            error('dynrel:dodgammals',...
                ['The person-specific dynamic difference-of-differences '...
                '(analyses 28/29) exists only as a gamma location-scale '...
                'result (family = ''gamma'', gammascale = 2) or a Gaussian '...
                'result (family = ''gaussian''); this result carries a '...
                'different family/parameterization stamp. How to fix: '...
                're-estimate with family = ''gaussian'', or with family = '...
                '''gamma'' (gammascale defaults to 2 for this design).']);
        end
        bb      = cell2mat(REL.out.b{s});           % draws x 4 (gamma: ABSOLUTE log-mean; Gaussian: natural-scale mean)
        bsig    = cell2mat(REL.out.b_sigma{s});     % draws x 4 (ABSOLUTE log-sigma, both families)
        bdim    = cell2mat(REL.out.b_dim{s});       % draws x 4 x KDIM (mean slopes, on the family's mean scale)
        bsigdim = cell2mat(REL.out.b_sigma_dim{s}); % draws x 4 x KDIM (log-sigma slopes)
        sdid    = cell2mat(REL.out.sd_id{s});       % draws x 8 [loc_1..4, logsd_1..4]
        sdtrl   = cell2mat(REL.out.sd_trl{s});      % draws x 4 (per-cell trial SD)
        corid   = cell2mat(REL.out.cor_id{s});      % draws x 8 x 8 (joint person corr)
        cortrl  = cell2mat(REL.out.cor_trl{s});     % draws x 4 x 4 (cross-cell trial)
        ssinfo  = REL.out.ssinfo{s};
        mu_ss = cell(1,4); er_ss = cell(1,4);
        for q = 1:4
            mu_ss{q} = cell2mat(REL.out.(sprintf('mu_ss%d',q)){s});
            er_ss{q} = cell2mat(REL.out.(sprintf('er_var_ss%d',q)){s});
        end

        dargs = {'b',bb,'b_sigma',bsig,'b_dim',bdim,'b_sigma_dim',bsigdim,...
            'sd_id',sdid,'cor_id',corid,'sd_trl',sdtrl,'cor_trl',cortrl,...
            'ndim',ndim,'z1',z1,'obs',obs_ev,'CI',ciperc};
        if ndim == 2; dargs = [dargs, {'z2',z2}]; end %#ok<AGROW>
        ssargs = {'b',bb,'b_sigma',bsig,'b_dim',bdim,'b_sigma_dim',bsigdim,...
            'sd_id',sdid,'cor_id',corid,'sd_trl',sdtrl,'cor_trl',cortrl,...
            'mu_ss1',mu_ss{1},'mu_ss2',mu_ss{2},'mu_ss3',mu_ss{3},...
            'mu_ss4',mu_ss{4},'er_ss1',er_ss{1},'er_ss2',er_ss{2},...
            'er_ss3',er_ss{3},'er_ss4',er_ss{4},'idtable',ssinfo,...
            'nprime',obs_ev,'ndim',ndim,'CI',ciperc};

        if isdodtrt
            %per-cell occasion n': default 1 (the OWNER DECISION every other
            %trt dynrel variant follows, and what the GUI's "1 occasion"
            %choice means); 'observed' (char or string) = this stratum's
            %per-cell medians; a scalar broadcasts to four cells. Non-finite
            %overrides are rejected BEFORE the max(1,round(.)) floor - the
            %two-input max drops NaN, so max(1,NaN) would silently become 1
            %(the documented trap the obs path also guards).
            if isempty(nocc_in)
                nocc4 = ones(1,4);
            elseif (ischar(nocc_in) || isstring(nocc_in)) && ...
                    strcmpi(char(nocc_in),'observed')
                nocc4 = max(1,round(REL.out.nocc(:,s)'));
            elseif isnumeric(nocc_in) && isscalar(nocc_in) && isfinite(nocc_in)
                nocc4 = repmat(max(1,round(nocc_in)),1,4);
            elseif isnumeric(nocc_in) && numel(nocc_in) == 4 && ...
                    all(isfinite(nocc_in))
                nocc4 = max(1,round(nocc_in(:)'));
            else
                error('dynrel:dodnocc',...
                    ['WARNING: the DoD occasion n'' must be ''observed'', a '...
                    'finite scalar, or a 4-element per-cell vector of finite '...
                    'counts.\n']);
            end
            nocc_used = round(mean(nocc4));
            sdocc = cell2mat(REL.out.sd_occ{s});
            sdtid = cell2mat(REL.out.sd_tid{s});
            sdoid = cell2mat(REL.out.sd_oid{s});
            sdto  = cell2mat(REL.out.sd_to{s});
            corocc = cell2mat(REL.out.cor_occ{s});
            cortid = cell2mat(REL.out.cor_tid{s});
            coroid = cell2mat(REL.out.cor_oid{s});
            corto  = cell2mat(REL.out.cor_to{s});
            facetargs = {'sd_occ',sdocc,'cor_occ',corocc,'sd_tid',sdtid,...
                'cor_tid',cortid,'sd_oid',sdoid,'cor_oid',coroid,...
                'sd_to',sdto,'cor_to',corto};
            %family selects the calculator set; the argument contracts are
            %identical (the slot SEMANTICS differ - Gaussian b/mu_ss* are
            %natural-scale means - and each calculator carries the right
            %interpretation; see SCIENTIFIC_FORMULA_AUDIT.md section 26)
            if isgamma_ls
                rel = psyrat_rel_dodiffdynrel_trt_gamma_ls(dargs{:},facetargs{:},...
                    'nocc',nocc4,'reltype',reltype);
                ssrel = psyrat_ssrel_dodiffdynrel_trt_gamma_ls(ssargs{:},...
                    facetargs{:},'nprime_o',nocc4,'reltype',reltype);
                varcomp = local_gamma_dod_varcomp_table_ls(bb,bsig,bdim,bsigdim,...
                    sdid,corid,ndim,REL.dim_names,ciperc,...
                    {'Trial',sdtrl,cortrl;'Occasion',sdocc,corocc;...
                    'Person x trial',sdtid,cortid;'Person x occasion',sdoid,coroid;...
                    'Trial x occasion',sdto,corto});
            else
                rel = psyrat_rel_dodiffdynrel_trt(dargs{:},facetargs{:},...
                    'nocc',nocc4,'reltype',reltype);
                ssrel = psyrat_ssrel_dodiffdynrel_trt(ssargs{:},...
                    facetargs{:},'nprime_o',nocc4,'reltype',reltype);
                varcomp = local_dod_varcomp_table(bb,bsig,bdim,bsigdim,...
                    sdid,corid,ndim,REL.dim_names,ciperc,...
                    {'Trial',sdtrl,cortrl;'Occasion',sdocc,corocc;...
                    'Person x trial',sdtid,cortid;'Person x occasion',sdoid,coroid;...
                    'Trial x occasion',sdto,corto});
            end
        else
            if isgamma_ls
                rel = psyrat_rel_dodiffdynrel_gamma_ls(dargs{:});
                ssrel = psyrat_ssrel_dodiffdynrel_gamma_ls(ssargs{:});
                varcomp = local_gamma_dod_varcomp_table_ls(bb,bsig,bdim,bsigdim,...
                    sdid,corid,ndim,REL.dim_names,ciperc,...
                    {'Trial',sdtrl,cortrl});
            else
                rel = psyrat_rel_dodiffdynrel(dargs{:});
                ssrel = psyrat_ssrel_dodiffdynrel(ssargs{:});
                varcomp = local_dod_varcomp_table(bb,bsig,bdim,bsigdim,...
                    sdid,corid,ndim,REL.dim_names,ciperc,...
                    {'Trial',sdtrl,cortrl});
            end
        end
        sige_pt = rel.sige_pt;
        dodflags = psyrat_dod_flag_counts(ssrel,size(bb,1));
        local_warn_dod_invalidity(dodflags,REL.out.labels{s});
        %the SURFACE records invalidity too (record-never-repair applies to
        %every evaluation, not only the per-participant table); without this
        %warning the surface's flags would be computed and then die unread
        if any(rel.invalid_frac(:) > 0)
            warning('dynrel:dodsurfaceinvalidity',...
                ['Stratum ''%s'': the typical-person DoD surface carries '...
                'flagged draws at %d of %d grid points (max flagged '...
                'fraction %.3f). Values are reported unrepaired.'],...
                char(string(REL.out.labels{s})),...
                nnz(rel.invalid_frac(:) > 0),numel(rel.invalid_frac),...
                max(rel.invalid_frac(:)));
        end
    elseif isdiff && ~isgamma
        %two-event difference score: reuse psyrat_diffrel via psyrat_rel_diffdynrel.
        %For the subject-level variant this is the typical-person (r_sigma = 0)
        %reference surface; the per-participant table is added below.
        idvc  = cell2mat(REL.out.id_varcov{s});    % draws x 2 x 2 (location)
        trvc  = cell2mat(REL.out.trl_varcov{s});   % draws x 2 x 2
        bsig  = cell2mat(REL.out.b_sigma{s});      % draws x 2 (per-event log-resid SD)
        bsdim = cell2mat(REL.out.b_sigma_dim{s});  % draws x KDIM (scale slopes)
        bdim  = cell2mat(REL.out.b_dim{s});        % draws x KDIM (mean slopes)
        dargs = {'id_varcov',idvc,'trl_varcov',trvc,'b_sigma',bsig,...
            'b_sigma_dim',bsdim,'ndim',ndim,'z1',z1,'obs',obs,'CI',ciperc};
        if ndim == 2; dargs = [dargs, {'z2',z2}]; end %#ok<AGROW>
        %concurrent variants feed the estimated residual correlation so the
        %z-conditioned residual covariance enters the difference variance. The group
        %(isdiffrescor) and population subject (issubjrescor) variants use a single
        %correlation; the per-subject variant (issubjrho, Stage 4e) uses the
        %population-AVERAGE correlation (REL.out.rescor = mean(rho)) for this
        %typical-person reference surface. rescor stays [] otherwise.
        rescor = [];
        if isdiffrescor || issubjrescor || issubjrho
            rescor = cell2mat(REL.out.rescor{s});      % draws x 1
            errvc  = cell2mat(REL.out.err_varcov{s});  % draws x 2 x 2 (at z = 0)
            dargs = [dargs, {'rescor',rescor}]; %#ok<AGROW>
        end
        rel = psyrat_rel_diffdynrel(dargs{:});
        varcomp = local_diff_varcomp_table(idvc,trvc,bsig,bsdim,bdim,ndim,...
            REL.dim_names,ciperc);
        if isdiffrescor || issubjrescor || issubjrho
            varcomp = [varcomp; local_rescor_varcomp(rescor,errvc,ciperc)];
        end
        if issubj || issubjrescor || issubjrho
            %per-participant phi at each subject's own z (psyrat_ssrel_diff), plus
            %the person scale-RE components. The per-subject residual covariance uses
            %each subject's OWN rho (issubjrho, a draws x NSUB matrix in
            %REL.out.rho_ss); the shared population rescor (issubjrescor); or 0
            %(issubj, non-concurrent).
            if issubjrho
                ss_corr = cell2mat(REL.out.rho_ss{s});   % draws x NSUB
            else
                ss_corr = rescor;                        % draws x 1 or []
            end
            ssrel = local_subj_ssrel_table(REL,s,idvc,trvc,bsig,bsdim,ndim,ciperc,...
                ss_corr);
            sdid = cell2mat(REL.out.sd_id{s});     % draws x 4 (1-2 loc, 3-4 scale)
            Lid  = cell2mat(REL.out.L_id{s});      % draws x 4 x 4 (joint cholesky)
            varcomp = [varcomp; local_subj_scale_varcomp(sdid,Lid,ciperc)];
            if issubjrho
                %Stage 4e: report the between-person spread of the per-subject
                %residual correlation (the population distribution of rho[s]).
                varcomp = [varcomp; local_rho_spread_varcomp(REL,s,ciperc)];
            end
        end
        sige_pt = [];
    elseif isdiff && isgamma_ls
        %two-event difference score, GAMMA LOCATION-SCALE (analysis 12,
        %family = 'gamma', gammascale = 2).
        %
        %THIS MUST BE TESTED BEFORE isgamma. isgamma_ls is defined as
        %isgamma && gammascale == 2, so a location-scale run satisfies BOTH; an
        %isgamma-first chain would route it to the log-nu branch below, which
        %reads REL.out.b_nu / b_nu_dim - fields a location-scale fit never
        %creates. Same ordering requirement as the two-facet dynrel branch
        %further down; see the note there.
        %
        %The dispersion side is the per-event log residual SD (b_sigma) and its
        %dimension slopes (b_sigma_dim); sd_id slots 3-4 hold the person
        %log-sigma SDs rather than log-nu SDs. Everything on the MEAN side is
        %read identically to the log-nu branch.
        bb     = cell2mat(REL.out.b{s});           % draws x 2 (per-event log-mean)
        bsig   = cell2mat(REL.out.b_sigma{s});     % draws x 2 (per-event log-sigma)
        bdim   = cell2mat(REL.out.b_dim{s});       % draws x KDIM (log-mean slopes)
        bsigdim = cell2mat(REL.out.b_sigma_dim{s});% draws x KDIM (log-sigma slopes)
        sdid   = cell2mat(REL.out.sd_id{s});       % draws x 4 [m1,m2,lsig1,lsig2]
        sdtrl  = cell2mat(REL.out.sd_trl{s});      % draws x 2 (per-event trial SD)
        corid  = cell2mat(REL.out.cor_id{s});      % draws x 4 x 4 (person corr)
        cori   = cell2mat(REL.out.cor_i{s});       % draws x 1 (cross-event trial)
        %MODULAR CUT-COPULA (analysis 19). The copula correlation exists only
        %for the RETAINED stage-1 draws, so the whole surface is recomputed on
        %that subset rather than pairing a short rho against full-length
        %margins. Owner decision 2026-08-09: subsample everything, keeping every
        %quantity per-draw consistent, in preference to a plug-in at the
        %posterior mean. The cost is the quadrature - one call per grid point per
        %retained draw - which the viewer's own cache absorbs.
        cutrho = [];
        if issubjrescor
            cop = REL.out.copula{s};
            sel = cop.sel(:);
            cutrho = cop.rho_e(:);
            cutflags = psyrat_copula_flag_counts(cop);
            local_warn_copula_invalidity(cutflags,REL.out.labels{s});
            bb = bb(sel,:); bsig = bsig(sel,:);
            bdim = bdim(sel,:); bsigdim = bsigdim(sel,:);
            sdid = sdid(sel,:); sdtrl = sdtrl(sel,:);
            corid = corid(sel,:,:); cori = cori(sel,:);
        end
        dargs = {'b',bb,'b_sigma',bsig,'b_dim',bdim,'b_sigma_dim',bsigdim,...
            'sd_id',sdid,'sd_trl',sdtrl,'cor_id',corid,'cor_i',cori,...
            'ndim',ndim,'z1',z1,'obs',obs_ev,'CI',ciperc};
        if ndim == 2; dargs = [dargs, {'z2',z2}]; end %#ok<AGROW>
        if ~isempty(cutrho); dargs = [dargs, {'rho',cutrho}]; end %#ok<AGROW>
        rel = psyrat_rel_diffdynrel_gamma_ls(dargs{:});
        varcomp = local_gamma_diff_varcomp_table_ls(bb,bsig,bdim,bsigdim,sdid,...
            sdtrl,corid,cori,ndim,REL.dim_names,ciperc);
        sige_pt = rel.sige_pt;
        %issubjrescor (analysis 19) must be admitted alongside issubj (13). Both
        %are folded into isdiff at the top of this file, so without it a
        %concurrent gamma run reaches this branch, produces the GROUP surface and
        %reports no per-participant table at all - silently, because nothing
        %throws. That is the exact failure the comment below records having
        %already happened once for case 13.
        if issubj || issubjrescor
            %SUBJECT-LEVEL read-out (analysis 13). The surface above is the
            %typical-person reference curve and is computed identically for both
            %analyses; this block adds the per-participant table that is the whole
            %difference between 12 and 13.
            %
            %WITHOUT THIS BLOCK A CASE-13 GAMMA RUN FAILS SILENTLY, and that is
            %worth stating because nothing would have thrown. issubj is folded into
            %isdiff at the top of this file, so an analysis-13 location-scale run
            %already satisfied `isdiff && isgamma_ls` and fell into this branch -
            %the GROUP-level one. Its REL.out is a strict superset of the eight
            %fields read above, so every one of them resolved, the group surface
            %came out, and the run simply reported no per-participant table at all.
            ss1 = cell2mat(REL.out.er_var_ss1{s});  % draws x NSUB (z=0, event 1)
            ss2 = cell2mat(REL.out.er_var_ss2{s});  % draws x NSUB (z=0, event 2)
            ssinfo = REL.out.ssinfo{s};
            if issubjrescor
                %CONCURRENT (analysis 19): the residuals are coupled, the
                %residual channel is POOLED per Rast & Clayson, and the
                %quadrature point carries each participant's own location
                %effect - which is why mu_ss1/2 are needed and analysis 13's
                %calculator cannot serve. Draws are already subsampled to the
                %retained set above, so the per-subject columns are indexed with
                %the same sel.
                mu1 = cell2mat(REL.out.mu_ss1{s});
                mu2 = cell2mat(REL.out.mu_ss2{s});
                [gg,dd,qdiag] = psyrat_ssrel_diffdynrel_rescor_gamma_ls(...
                    'b',bb,'b_sigma',bsig,...
                    'b_dim',bdim,'b_sigma_dim',bsigdim,'sd_id',sdid,'sd_trl',sdtrl,...
                    'cor_id',corid,'cor_i',cori,...
                    'er_ss1',ss1(sel,:),'er_ss2',ss2(sel,:),...
                    'mu_ss1',mu1(sel,:),'mu_ss2',mu2(sel,:),'rho',cutrho,...
                    'idtable',ssinfo,'ndim',ndim,'CI',ciperc);
                quaddiag = qdiag;
            else
                [gg,dd] = psyrat_ssrel_diffdynrel_gamma_ls('b',bb,'b_sigma',bsig,...
                    'b_dim',bdim,'b_sigma_dim',bsigdim,'sd_id',sdid,'sd_trl',sdtrl,...
                    'cor_id',corid,'cor_i',cori,'er_ss1',ss1,'er_ss2',ss2,...
                    'idtable',ssinfo,'ndim',ndim,'CI',ciperc);
            end
            ssrel = local_pack_ssrel_diff(ssinfo,gg,dd);
            %the person scale block is what makes participants differ, so it is
            %surfaced alongside the variance components. Only the CROSS-EVENT scale
            %correlation is appended: local_gamma_diff_varcomp_table_ls above is
            %built from the same sdid/corid and already reports the two person scale
            %SDs and the two within-event location-scale correlations, so appending
            %the full block duplicated four rows under second names. Lid is not
            %stored under gamma, so the correlation comes from cor_id directly.
            varcomp = [varcomp; local_subj_scale_varcomp_gamma(corid,ciperc)];
            if issubjrescor
                %the copula correlation is the estimand this design adds; report
                %it in the table like the Gaussian concurrent variants do.
                varcomp = [varcomp; local_cutrho_varcomp(cutrho,ciperc)]; %#ok<AGROW>
            end
        end
    elseif isdiff && isgamma
        %two-event difference score, GAMMA / scaled-chi-square (analysis 12,
        %family = 'gamma'). MUST be a separate branch from the Gaussian isdiff
        %path above: that path reads REL.out.b_sigma / b_sigma_dim (a log-residual
        %submodel) and REL.out.id_varcov / trl_varcov (fixed across z), NONE of
        %which a scaled-chi-square fit produces. Here the dispersion is the
        %per-event log-nu (b_nu + b_nu_dim) and the observed-scale (co)variances
        %are z-DEPENDENT, so they are rebuilt per z inside
        %psyrat_rel_diffdynrel_gamma via psyrat_gamma_varcomps +
        %psyrat_gamma_crosscov. Estimand #2 (marginal over person dispersion).
        bb     = cell2mat(REL.out.b{s});         % draws x 2 (per-event log-mean)
        bnu    = cell2mat(REL.out.b_nu{s});      % draws x 2 (per-event log-nu)
        bdim   = cell2mat(REL.out.b_dim{s});     % draws x KDIM (log-mean slopes)
        bnudim = cell2mat(REL.out.b_nu_dim{s});  % draws x KDIM (log-nu slopes)
        sdid   = cell2mat(REL.out.sd_id{s});     % draws x 4 [m1,m2,lognu1,lognu2]
        sdtrl  = cell2mat(REL.out.sd_trl{s});    % draws x 2 (per-event trial SD)
        corid  = cell2mat(REL.out.cor_id{s});    % draws x 4 x 4 (person corr)
        cori   = cell2mat(REL.out.cor_i{s});     % draws x 1 (cross-event trial)
        dargs = {'b',bb,'b_nu',bnu,'b_dim',bdim,'b_nu_dim',bnudim,...
            'sd_id',sdid,'sd_trl',sdtrl,'cor_id',corid,'cor_i',cori,...
            'ndim',ndim,'z1',z1,'obs',obs_ev,'CI',ciperc};
        if ndim == 2; dargs = [dargs, {'z2',z2}]; end %#ok<AGROW>
        rel = psyrat_rel_diffdynrel_gamma(dargs{:});
        varcomp = local_gamma_diff_varcomp_table(bb,bnu,bdim,bnudim,sdid,...
            sdtrl,corid,cori,ndim,REL.dim_names,ciperc);
        sige_pt = rel.sige_pt;
    elseif istrt
        %trial + occasion two-facet non-difference (Stage 3): reuse the
        %test-retest two-facet coefficient via psyrat_rel_dynrel_trt, feeding the
        %z-conditioned residual as 'err' and the five crossed facet SDs as the
        %static facet arguments. nocc resolves to 1 by default (owner decision)
        %or this stratum's observed occasion count when requested. The five
        %crossed facet SDs and nocc are shared by the Gaussian and gamma branches.
        nocc_used = local_resolve_nocc(nocc_in,REL,s);
        s_occ = REL.out.sig_occ(:,s);
        s_trl = REL.out.sig_trl(:,s);
        s_txp = REL.out.sig_trlxid(:,s);
        s_oxp = REL.out.sig_occxid(:,s);
        s_txo = REL.out.sig_trlxocc(:,s);
        if isgamma_ls
            %gamma LOCATION-SCALE two-facet dynrel (analysis 14, family = 'gamma',
            %gammascale = 2). The scale submodel is the log residual SD, so the
            %slots are the Gaussian dynrel_trt's (pop_sdlog, b_sigma, ind_sdlog)
            %and the observed-score conversion is psyrat_gamma_varcomps_trt_ls,
            %which takes no covariance argument: the residual never touches mu.
            %
            %THIS MUST BE TESTED BEFORE isgamma. isgamma_ls is defined as
            %isgamma && gammascale == 2, so a location-scale run satisfies BOTH; an
            %isgamma-first chain would route it to the log-nu branch, which reads
            %pop_lognu/b_nu - fields the location-scale REL.out deliberately never
            %creates - and would error on the first stratum.
            %
            %WATCH gro_sds COLUMN 2. It is the SAME SLOT that holds the person
            %log-nu SD (s_v) under gammascale = 1 and holds the person log-sigma SD
            %(s_sd) here. Handing it to the log-nu branch's 'sig_v' produces
            %numbers, not an error, so the two branches must never share this read.
            %
            %The five crossed facet SDs (s_occ ... s_txo, read above) are shared
            %with both other branches unchanged: they live on the MEAN submodel,
            %which this parameterization leaves untouched.
            alpha0  = REL.out.mu(:,s);                   % log-mean intercept (Intercept)
            logsig0 = REL.out.pop_sdlog(:,s);            % log residual-SD intercept (Intercept_sigma)
            gro     = cell2mat(REL.out.gro_sds{s});      % draws x 2 (s_p, s_sd)
            sig_p   = gro(:,1);
            sig_sd  = gro(:,2);
            bmean   = cell2mat(REL.out.b{s});            % draws x KDIM (log-mean slopes)
            bsig    = cell2mat(REL.out.b_sigma{s});      % draws x KDIM (log residual-SD slopes)
            chol    = cell2mat(REL.out.chol_corrmat{s}); % draws x 2 x 2
            %person mean<->log-sigma correlation: for a 2x2 cholesky L, R(1,2)=L(2,1).
            %Reported in the component table, and deliberately NOT passed to the
            %converter - see psyrat_gamma_varcomps_trt_ls.
            rho_ps  = chol(:,2,1);
            args = {'alpha0',alpha0,'b',bmean,'logsig0',logsig0,'b_sigma',bsig,...
                'sig_p',sig_p,'sig_sd',sig_sd,...
                'sig_occ',s_occ,'sig_trl',s_trl,'sig_trlxid',s_txp,...
                'sig_occxid',s_oxp,'sig_trlxocc',s_txo,...
                'ndim',ndim,'z1',z1,'obs',obs,'nocc',nocc_used,...
                'reltype',reltype,'CI',ciperc};
            if ndim == 2; args = [args, {'z2',z2}]; end %#ok<AGROW>
            rel = psyrat_rel_dynrel_trt_gamma_ls(args{:});
            varcomp = local_gamma_trt_ls_varcomp_table(alpha0,logsig0,bmean,bsig,...
                sig_p,sig_sd,s_occ,s_trl,s_txp,s_oxp,s_txo,rho_ps,ndim,...
                REL.dim_names,ciperc);
            sige_pt = rel.sige_pt;
            if issubjnondiff_trt
                %per-participant CONDITIONAL phi_s(z)/G_s(z) (analysis 27 under
                %gammascale = 2). Each participant's residual is read straight off
                %ind_sdlog rather than derived from two person effects, and their
                %own log-mean does not enter at all.
                ssrel = local_subj_ssrel_table_nondiff_trt_gamma_ls(REL,s,alpha0,...
                    bmean,logsig0,bsig,sig_p,s_occ,s_trl,s_txp,s_oxp,s_txo,...
                    reltype,nocc_used,ciperc,ndim);
            end
        elseif isgamma
            %gamma / scaled-chi-square two-facet dynrel: no log-residual submodel;
            %the residual is the mean-coupled dispersion, so the log-mean and
            %log-nu submodels are dimension-conditioned and the observed-score
            %components are converted per z via psyrat_gamma_varcomps_trt
            %(psyrat_rel_dynrel_trt_gamma). Estimand #2 (marginal over dispersion).
            alpha0 = REL.out.mu(:,s);                   % log-mean intercept (Intercept)
            lognu0 = REL.out.pop_lognu(:,s);            % log-nu intercept (Intercept_nu)
            gro    = cell2mat(REL.out.gro_sds{s});      % draws x 2 (s_p, s_v)
            sig_p  = gro(:,1);
            sig_v  = gro(:,2);
            bmean  = cell2mat(REL.out.b{s});            % draws x KDIM (log-mean slopes)
            bnu    = cell2mat(REL.out.b_nu{s});         % draws x KDIM (log-nu slopes)
            chol   = cell2mat(REL.out.chol_corrmat{s}); % draws x 2 x 2
            rho_pv = chol(:,2,1);                        % person mean<->log-nu corr
            args = {'alpha0',alpha0,'b',bmean,'lognu0',lognu0,'b_nu',bnu,...
                'sig_p',sig_p,'sig_v',sig_v,'rho_pv',rho_pv,...
                'sig_occ',s_occ,'sig_trl',s_trl,'sig_trlxid',s_txp,...
                'sig_occxid',s_oxp,'sig_trlxocc',s_txo,...
                'ndim',ndim,'z1',z1,'obs',obs,'nocc',nocc_used,...
                'reltype',reltype,'CI',ciperc};
            if ndim == 2; args = [args, {'z2',z2}]; end %#ok<AGROW>
            rel = psyrat_rel_dynrel_trt_gamma(args{:});
            varcomp = local_gamma_trt_varcomp_table(alpha0,lognu0,bmean,bnu,...
                sig_p,sig_v,s_occ,s_trl,s_txp,s_oxp,s_txo,rho_pv,ndim,...
                REL.dim_names,ciperc);
            sige_pt = rel.sige_pt;
            if issubjnondiff_trt
                %per-participant CONDITIONAL phi_s(z)/G_s(z) for the gamma
                %two-facet subject-level design (analysis 27). The surface above
                %is the population reference and MARGINALIZES over person
                %dispersion (estimand #2); this table instead CONDITIONS on each
                %participant's own log-mean (ind_bs) and log-nu (ind_nu) at their
                %own z. Unlike the Gaussian branch below, the five crossed facet
                %SDs cannot be broadcast as population values: under gamma they
                %also move with z - see psyrat_ssrel_dynrel_trt_gamma.
                ssrel = local_subj_ssrel_table_nondiff_trt_gamma(REL,s,alpha0,...
                    bmean,lognu0,bnu,sig_p,s_occ,s_trl,s_txp,s_oxp,s_txo,...
                    reltype,nocc_used,ciperc,ndim);
            end
        else
            gro   = cell2mat(REL.out.gro_sds{s});      % draws x 2 (sigma_p, sigma_delta_p)
            sig_p = gro(:,1);
            sig_dp = gro(:,2);
            log0  = REL.out.pop_sdlog(:,s);            % Intercept_sigma
            bsig  = cell2mat(REL.out.b_sigma{s});      % draws x KDIM (scale slopes)
            bmean = cell2mat(REL.out.b{s});            % draws x KDIM (mean slopes)
            chol  = cell2mat(REL.out.chol_corrmat{s}); % draws x 2 x 2
            rho   = chol(:,2,1);                        % person location-scale corr
            args = {'sig_p',sig_p,'sig_log0',log0,'b_sigma',bsig,...
                'sig_occ',s_occ,'sig_trl',s_trl,'sig_trlxid',s_txp,...
                'sig_occxid',s_oxp,'sig_trlxocc',s_txo,...
                'ndim',ndim,'z1',z1,'obs',obs,'nocc',nocc_used,...
                'reltype',reltype,'CI',ciperc,'sig_dp',sig_dp,'marginal',marginal};
            if ndim == 2; args = [args, {'z2',z2}]; end %#ok<AGROW>
            rel = psyrat_rel_dynrel_trt(args{:});
            varcomp = local_trt_varcomp_table(sig_p,s_occ,s_trl,s_txp,s_oxp,s_txo,...
                log0,bsig,bmean,sig_dp,rho,ndim,REL.dim_names,ciperc);
            sige_pt = rel.sige_pt;
            if issubjnondiff_trt
                %per-participant phi_s(z) at each subject's own z from their own
                %residual (ind_sdlog); the five crossed facet SDs are population,
                %so psyrat_ssrel_trt broadcasts them across subjects (as the
                %static case-25 path). The population surface above is the
                %typical-person reference.
                ssrel = local_subj_ssrel_table_nondiff_trt(REL,s,sig_p,log0,bsig,...
                    s_occ,s_trl,s_txp,s_oxp,s_txo,reltype,nocc_used,ciperc,ndim);
            end
        end
    elseif (isdifftrt || issubjtrt) && isgamma_ls
        %trial + occasion two-facet difference, GAMMA LOCATION-SCALE. Today the
        %only design that reaches here is analysis 20, the two-facet MODULAR
        %CUT-COPULA (family = 'gamma', gammascale = 2, concurrent,
        %subject-level); its non-concurrent siblings under gamma are not enabled.
        %
        %THIS MUST BE TESTED BEFORE THE GAUSSIAN isdifftrt BRANCH BELOW, for the
        %reason recorded on the one-facet chain: the Gaussian branch reads
        %REL.out.id_varcov / trl_varcov / occ_varcov / tid_varcov / oid_varcov /
        %to_varcov / rescor / err_varcov, and the gamma arm fills NONE of them
        %(its observed-scale (co)variances are z-DEPENDENT, so they are rebuilt
        %per z downstream, and the residual correlation is not a parameter of the
        %model at all). Falling through would either error on a missing field or,
        %worse, read a slot that happens to exist and report a number computed
        %from the wrong quantity.
        %
        %FAIL-LOUD GUARD. The barred siblings - 15 (ic_diff_dynrel_trt), 16
        %(ic_diff_dynrel_sserrvar_trt), 18 (ic_diff_dynrel_trt_rescor) and 22
        %(ic_diff_dynrel_sserrvar_trt_rho) - satisfy this branch's condition
        %too, and nothing below distinguishes them: a stored result from a
        %lifted upstream guard or a future increment would be summarized
        %through the analysis-20 layout, erroring on a missing field at best
        %and reporting numbers computed from the wrong quantities at worst.
        %Error by name instead; an increment that enables a sibling under
        %gamma must extend this branch deliberately.
        if ~issubjtrtrescor
            error('dynrel:gammatrtdesign', ...
                ['Analysis ''%s'' reached the two-facet gamma ' ...
                'location-scale summary branch, which supports only the ' ...
                'modular cut-copula design (analysis 20, ' ...
                'ic_diff_dynrel_sserrvar_trt_rescor). This design is not ' ...
                'enabled under family = ''gamma''.'], REL.analysis);
        end
        bb      = cell2mat(REL.out.b{s});            % draws x 2 (per-event log-mean)
        bsig    = cell2mat(REL.out.b_sigma{s});      % draws x 2 (per-event log-sigma)
        bdim    = cell2mat(REL.out.b_dim{s});        % draws x KDIM (log-mean slopes)
        bsigdim = cell2mat(REL.out.b_sigma_dim{s});  % draws x KDIM (log-sigma slopes)
        sdid    = cell2mat(REL.out.sd_id{s});        % draws x 4 [m1,m2,lsig1,lsig2]
        corid   = cell2mat(REL.out.cor_id{s});       % draws x 4 x 4 (person corr)
        %the five crossed facets: per-event log-SDs and cross-event correlations.
        %sd_tid is trial x person and sd_oid is occasion x person - the names
        %cross over against the converter's s_pt / s_po arguments, and the
        %calculators below are where that mapping is stated in full.
        sdtrl = cell2mat(REL.out.sd_trl{s});
        sdocc = cell2mat(REL.out.sd_occ{s});
        sdtid = cell2mat(REL.out.sd_tid{s});
        sdoid = cell2mat(REL.out.sd_oid{s});
        sdto  = cell2mat(REL.out.sd_to{s});
        cort  = cell2mat(REL.out.cor_t{s});
        coro  = cell2mat(REL.out.cor_o{s});
        corpt = cell2mat(REL.out.cor_pt{s});
        corpo = cell2mat(REL.out.cor_po{s});
        corot = cell2mat(REL.out.cor_ot{s});
        nocc_used = local_resolve_nocc(nocc_in,REL,s);

        %MODULAR CUT-COPULA (analysis 20), exactly as the one-facet branch above:
        %the copula correlation exists only for the RETAINED stage-1 draws, so
        %the whole surface is recomputed on that subset rather than pairing a
        %short rho against full-length margins. Owner decision 2026-08-09.
        cutrho = [];
        if issubjtrtrescor
            cop = REL.out.copula{s};
            sel = cop.sel(:);
            cutrho = cop.rho_e(:);
            cutflags = psyrat_copula_flag_counts(cop);
            local_warn_copula_invalidity(cutflags,REL.out.labels{s});
            bb = bb(sel,:); bsig = bsig(sel,:);
            bdim = bdim(sel,:); bsigdim = bsigdim(sel,:);
            sdid = sdid(sel,:); corid = corid(sel,:,:);
            sdtrl = sdtrl(sel,:); sdocc = sdocc(sel,:); sdtid = sdtid(sel,:);
            sdoid = sdoid(sel,:); sdto = sdto(sel,:);
            cort = cort(sel,:); coro = coro(sel,:); corpt = corpt(sel,:);
            corpo = corpo(sel,:); corot = corot(sel,:);
        end

        facetargs = {'sd_trl',sdtrl,'sd_occ',sdocc,'sd_tid',sdtid,...
            'sd_oid',sdoid,'sd_to',sdto,'cor_t',cort,'cor_o',coro,...
            'cor_pt',corpt,'cor_po',corpo,'cor_ot',corot};

        %A SCALAR trial n' here, not the per-event obs_ev the one-facet gamma
        %difference branch builds, and that is a property of the design rather
        %than a gap. These events are CONCURRENT - both are scored from the same
        %trials - so their per-event trial counts cannot differ, which is also
        %why the reference implementation carries a single n_i for this design.
        %(The unbalanced 49-vs-341 case the obs_ev block exists for is a
        %NON-concurrent one.) Per-participant counts still vary and still reach
        %the per-participant table, through ssinfo's trls1/trls2.
        dargs = [{'b',bb,'b_sigma',bsig,'b_dim',bdim,'b_sigma_dim',bsigdim,...
            'sd_id',sdid,'cor_id',corid}, facetargs, ...
            {'ndim',ndim,'z1',z1,'obs',obs,'nocc',nocc_used,...
            'reltype',reltype,'CI',ciperc}];
        if ndim == 2; dargs = [dargs, {'z2',z2}]; end %#ok<AGROW>
        if ~isempty(cutrho); dargs = [dargs, {'rho',cutrho}]; end %#ok<AGROW>
        rel = psyrat_rel_diffdynrel_trt_gamma_ls(dargs{:});
        varcomp = local_gamma_difftrt_varcomp_table_ls(bb,bsig,bdim,bsigdim,...
            sdid,sdtrl,sdocc,sdtid,sdoid,sdto,corid,cort,coro,corpt,corpo,...
            corot,ndim,REL.dim_names,ciperc);
        sige_pt = rel.sige_pt;

        if issubjtrtrescor
            %SUBJECT-LEVEL read-out. The surface above is the typical-person
            %reference curve; this adds the per-participant table, which is the
            %whole point of analysis 20. Draws are already subsampled to the
            %retained set above, so the per-subject columns take the same sel.
            ss1 = cell2mat(REL.out.er_var_ss1{s});   % draws x NSUB (z=0, event 1)
            ss2 = cell2mat(REL.out.er_var_ss2{s});
            mu1 = cell2mat(REL.out.mu_ss1{s});
            mu2 = cell2mat(REL.out.mu_ss2{s});
            ssinfo = REL.out.ssinfo{s};
            [gg,dd,qdiag] = psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls(...
                'b',bb,'b_sigma',bsig,'b_dim',bdim,'b_sigma_dim',bsigdim,...
                'sd_id',sdid,'cor_id',corid, facetargs{:}, ...
                'er_ss1',ss1(sel,:),'er_ss2',ss2(sel,:),...
                'mu_ss1',mu1(sel,:),'mu_ss2',mu2(sel,:),'rho',cutrho,...
                'idtable',ssinfo,'ndim',ndim,'reltype',reltype,...
                'nocc',nocc_used,'CI',ciperc);
            quaddiag = qdiag;
            ssrel = local_pack_ssrel_diff(ssinfo,gg,dd);
            %the person scale block is what makes participants differ. Only the
            %CROSS-EVENT scale correlation corid(:,3,4) is appended; the varcomp
            %table above already reports the two person scale SDs and both
            %within-event location-scale correlations from the same sdid/corid,
            %so appending the full block would duplicate four rows under second
            %names. L_id is not stored under gamma, so this comes from cor_id
            %directly.
            varcomp = [varcomp; local_subj_scale_varcomp_gamma(corid,ciperc)];
            %the copula correlation is the estimand this design adds; report it
            %in the table like the Gaussian concurrent variants do.
            varcomp = [varcomp; local_cutrho_varcomp(cutrho,ciperc)]; %#ok<AGROW>
        end
    elseif isdifftrt || issubjtrt
        %trial + occasion two-facet difference. Group-level (Stage 4a, isdifftrt)
        %and subject-level (Stage 4b, issubjtrt) share the typical-person reference
        %surface and the static difference-component table: reuse the two-facet
        %difference coefficient via psyrat_rel_diffdynrel_trt, feeding the
        %z-conditioned per-event residual and the six event covariance blocks. The
        %subject-level variant adds a per-participant table (each subject's phi at
        %their own z) and the person scale-RE components below.
        idvc  = cell2mat(REL.out.id_varcov{s});   % draws x 2 x 2 (person)
        trvc  = cell2mat(REL.out.trl_varcov{s});  % draws x 2 x 2 (trial)
        ocvc  = cell2mat(REL.out.occ_varcov{s});  % draws x 2 x 2 (occasion)
        tidvc = cell2mat(REL.out.tid_varcov{s});  % draws x 2 x 2 (person x trial)
        oidvc = cell2mat(REL.out.oid_varcov{s});  % draws x 2 x 2 (person x occasion)
        tovc  = cell2mat(REL.out.to_varcov{s});   % draws x 2 x 2 (trial x occasion)
        bsig  = cell2mat(REL.out.b_sigma{s});     % draws x 2 (per-event log-resid SD)
        bsdim = cell2mat(REL.out.b_sigma_dim{s}); % draws x KDIM (scale slopes)
        bdim  = cell2mat(REL.out.b_dim{s});       % draws x KDIM (mean slopes)
        nocc_used = local_resolve_nocc(nocc_in,REL,s);
        dargs = {'id_varcov',idvc,'trl_varcov',trvc,'occ_varcov',ocvc,...
            'tid_varcov',tidvc,'oid_varcov',oidvc,'to_varcov',tovc,...
            'b_sigma',bsig,'b_sigma_dim',bsdim,'ndim',ndim,'z1',z1,...
            'obs',obs,'nocc',nocc_used,'reltype',reltype,'CI',ciperc};
        if ndim == 2; dargs = [dargs, {'z2',z2}]; end %#ok<AGROW>
        %concurrent variants feed the estimated residual correlation. Group
        %(isdifftrtrescor) and population subject (issubjtrtrescor) use a single
        %correlation; the per-subject variant (issubjtrtrho, Stage 4e) uses the
        %population-AVERAGE correlation (REL.out.rescor = mean(rho)) for this
        %typical-person reference surface. rescor stays [] otherwise.
        rescor = [];
        if isdifftrtrescor || issubjtrtrescor || issubjtrtrho
            rescor = cell2mat(REL.out.rescor{s});      % draws x 1
            errvc  = cell2mat(REL.out.err_varcov{s});  % draws x 2 x 2 (at z = 0)
            dargs = [dargs, {'rescor',rescor}]; %#ok<AGROW>
        end
        rel = psyrat_rel_diffdynrel_trt(dargs{:});
        varcomp = local_difftrt_varcomp_table(idvc,trvc,ocvc,tidvc,oidvc,tovc,...
            bsig,bsdim,bdim,ndim,REL.dim_names,ciperc);
        if isdifftrtrescor || issubjtrtrescor || issubjtrtrho
            varcomp = [varcomp; local_rescor_varcomp(rescor,errvc,ciperc)];
        end
        if issubjtrt || issubjtrtrescor || issubjtrtrho
            %per-participant phi at each subject's own z (psyrat_ssrel_diff_trt),
            %plus the joint 4x4 person scale-RE components. The per-subject residual
            %covariance uses each subject's OWN rho (issubjtrtrho, a draws x NSUB
            %matrix in REL.out.rho_ss); the shared population rescor (issubjtrtrescor);
            %or 0 (issubjtrt, non-concurrent).
            if issubjtrtrho
                ss_corr = cell2mat(REL.out.rho_ss{s});   % draws x NSUB
            else
                ss_corr = rescor;                        % draws x 1 or []
            end
            ssrel = local_subj_ssrel_table_trt(REL,s,idvc,trvc,ocvc,tidvc,...
                oidvc,tovc,bsig,bsdim,ndim,reltype,nocc_used,ciperc,ss_corr);
            sdid = cell2mat(REL.out.sd_id{s});     % draws x 4 (1-2 loc, 3-4 scale)
            Lid  = cell2mat(REL.out.L_id{s});      % draws x 4 x 4 (joint cholesky)
            varcomp = [varcomp; local_subj_scale_varcomp(sdid,Lid,ciperc)];
            if issubjtrtrho
                %Stage 4e: report the between-person spread of the per-subject
                %residual correlation (the population distribution of rho[s]).
                varcomp = [varcomp; local_rho_spread_varcomp(REL,s,ciperc)];
            end
        end
        sige_pt = [];
    elseif isgamma_ls
        %one-facet non-difference, LOCATION-SCALE GAMMA (Stage 1, family =
        %'gamma', gammascale = 2). The scale submodel is the log residual SD, so
        %the slots are the Gaussian dynrel's (pop_sdlog, b_sigma, ind_sdlog) and
        %the observed-score conversion is psyrat_gamma_varcomps_ls, which takes
        %no covariance argument: the residual never touches mu.
        %
        %WATCH gro_sds COLUMN 2. It is the SAME SLOT that holds the person log-nu
        %SD (s_v) under gammascale = 1 and holds the person log-sigma SD (s_sd)
        %here. Handing it to the log-nu branch's 'sig_v' produces numbers, not an
        %error, so the two branches must never share this read.
        alpha0  = REL.out.mu(:,s);                   % log-mean intercept (Intercept)
        logsig0 = REL.out.pop_sdlog(:,s);            % log residual-SD intercept (Intercept_sigma)
        gro     = cell2mat(REL.out.gro_sds{s});      % draws x 2 (s_p, s_sd)
        sig_p   = gro(:,1);
        sig_sd  = gro(:,2);
        bmean   = cell2mat(REL.out.b{s});            % draws x KDIM (log-mean slopes)
        bsig    = cell2mat(REL.out.b_sigma{s});      % draws x KDIM (log residual-SD slopes)
        sig_i   = REL.out.sig_trl(:,s);              % trial log-mean SD (s_i)
        chol    = cell2mat(REL.out.chol_corrmat{s}); % draws x 2 x 2
        %person mean<->log-sigma correlation: for a 2x2 cholesky L, R(1,2)=L(2,1).
        %Reported in the component table, and deliberately NOT passed to the
        %converter - see psyrat_gamma_varcomps_ls.
        rho_ps  = chol(:,2,1);
        args = {'alpha0',alpha0,'b',bmean,'logsig0',logsig0,'b_sigma',bsig,...
            'sig_p',sig_p,'sig_sd',sig_sd,'sig_i',sig_i,...
            'ndim',ndim,'z1',z1,'obs',obs,'CI',ciperc};
        if ndim == 2; args = [args, {'z2',z2}]; end %#ok<AGROW>
        rel = psyrat_rel_dynrel_gamma_ls(args{:});
        varcomp = local_gamma_ls_varcomp_table(alpha0,logsig0,bmean,bsig,...
            sig_p,sig_sd,sig_i,rho_ps,ndim,REL.dim_names,ciperc);
        sige_pt = rel.sige_pt;
        if issubjnondiff
            %per-participant CONDITIONAL phi_s(z)/G_s(z) (analysis 26 under
            %gammascale = 2). Each participant's residual is read straight off
            %ind_sdlog rather than derived from two person effects, and their own
            %log-mean does not enter at all.
            ssrel = local_subj_ssrel_table_nondiff_gamma_ls(REL,s,alpha0,bmean,...
                logsig0,bsig,sig_p,sig_i,ndim,ciperc);
        end
    elseif isgamma
        %one-facet non-difference, GAMMA / scaled-chi-square (Stage 1, family =
        %'gamma'). No log-residual submodel: the residual is the mean-coupled
        %scaled-chi-square dispersion 2*mu(z)^2/nu(z), so the log-mean and log-nu
        %submodels are both dimension-conditioned and the observed-score
        %components are converted per z via psyrat_gamma_varcomps
        %(psyrat_rel_dynrel_gamma). Estimand #2 (marginal over person dispersion).
        alpha0 = REL.out.mu(:,s);                   % log-mean intercept (Intercept)
        lognu0 = REL.out.pop_lognu(:,s);            % log-nu intercept (Intercept_nu)
        gro    = cell2mat(REL.out.gro_sds{s});      % draws x 2 (s_p, s_v)
        sig_p  = gro(:,1);
        sig_v  = gro(:,2);
        bmean  = cell2mat(REL.out.b{s});            % draws x KDIM (log-mean slopes)
        bnu    = cell2mat(REL.out.b_nu{s});         % draws x KDIM (log-nu slopes)
        sig_i  = REL.out.sig_trl(:,s);              % trial log-mean SD (s_i)
        chol   = cell2mat(REL.out.chol_corrmat{s}); % draws x 2 x 2
        %person mean<->log-nu correlation: for a 2x2 cholesky L, R(1,2)=L(2,1)
        rho_pv = chol(:,2,1);
        args = {'alpha0',alpha0,'b',bmean,'lognu0',lognu0,'b_nu',bnu,...
            'sig_p',sig_p,'sig_v',sig_v,'rho_pv',rho_pv,'sig_i',sig_i,...
            'ndim',ndim,'z1',z1,'obs',obs,'CI',ciperc};
        if ndim == 2; args = [args, {'z2',z2}]; end %#ok<AGROW>
        rel = psyrat_rel_dynrel_gamma(args{:});
        varcomp = local_gamma_varcomp_table(alpha0,lognu0,bmean,bnu,...
            sig_p,sig_v,sig_i,rho_pv,ndim,REL.dim_names,ciperc);
        sige_pt = rel.sige_pt;
        if issubjnondiff
            %per-participant CONDITIONAL phi_s(z)/G_s(z) for the gamma
            %subject-level design (analysis 26). The surface above is the
            %population reference and MARGINALIZES over person dispersion
            %(estimand #2); this table instead CONDITIONS on each participant's
            %own log-mean (ind_bs) and log-nu (ind_nu) at their own z. Unlike the
            %Gaussian branch below it cannot reuse psyrat_ssrel in a single call,
            %because under gamma sigma_p2 and sigma_i2 also move with z and every
            %participant sits at a different z - see psyrat_ssrel_dynrel_gamma.
            ssrel = local_subj_ssrel_table_nondiff_gamma(REL,s,alpha0,bmean,...
                lognu0,bnu,sig_p,sig_i,ndim,ciperc);
        end
    else
        %one-facet non-difference (Stage 1)
        gro   = cell2mat(REL.out.gro_sds{s});      % draws x 2 (sigma_p, sigma_delta_p)
        sig_p = gro(:,1);
        sig_dp = gro(:,2);
        log0  = REL.out.pop_sdlog(:,s);            % Intercept_sigma
        bsig  = cell2mat(REL.out.b_sigma{s});      % draws x KDIM
        bmean = cell2mat(REL.out.b{s});            % draws x KDIM
        sig_i = REL.out.sig_trl(:,s);
        chol  = cell2mat(REL.out.chol_corrmat{s}); % draws x 2 x 2
        %person location/log-scale correlation: for a 2x2 cholesky L, R(1,2)=L(2,1)
        rho   = chol(:,2,1);
        args = {'sig_p',sig_p,'sig_log0',log0,'b_sigma',bsig,'sig_i',sig_i,...
            'ndim',ndim,'z1',z1,'obs',obs,'CI',ciperc,...
            'sig_dp',sig_dp,'marginal',marginal};
        if ndim == 2; args = [args, {'z2',z2}]; end %#ok<AGROW>
        rel = psyrat_rel_dynrel(args{:});
        varcomp = local_varcomp_table(sig_p,sig_i,log0,bsig,bmean,sig_dp,rho,...
            ndim,REL.dim_names,ciperc);
        sige_pt = rel.sige_pt;
        if issubjnondiff
            %per-participant phi_s(z) at each subject's own z from their own
            %residual (ind_sdlog), reusing psyrat_ssrel (the same per-subject
            %calc the static ic_sserrvar uses). The population surface above is
            %the typical-person reference.
            ssrel = local_subj_ssrel_table_nondiff(REL,s,sig_p,log0,bsig,sig_i,...
                ndim,ciperc);
        end
    end

    st = struct();
    st.label = REL.out.labels{s};
    st.obs = obs;
    st.obs_ev = obs_ev; %[n1 n2] gamma difference / 1x4 per-cell for DoD 28/29 (both families) / [] otherwise
    st.nocc = nocc_used; %occasion n' (only the two-facet variant; [] otherwise)
    st.z1 = rel.z1;
    if ndim == 2; st.z2 = rel.z2; else; st.z2 = []; end
    st.G = rel.G; st.D = rel.D; st.ICCg = rel.ICCg; st.ICCd = rel.ICCd;
    st.sige_pt = sige_pt;
    st.varcomp = varcomp;
    st.ssrel_table = ssrel;
    %empty except for the modular cut-copula designs, where it records how
    %hard the copula quadrature had to work per participant.
    st.quaddiag = quaddiag;
    %empty except for the modular cut-copula designs, where it carries the
    %stage-2 diagnostic counts (n_draws/n_invalid/n_calibration) so the flags
    %the stage records are visible downstream instead of dying in REL.out.
    st.copula_flags = cutflags;
    %empty except for the nonconcurrent DoD designs, where it carries the
    %signed-projection invalidity counts (record-never-repair; the reference's
    %posterior_invalid_* convention) so flagged draws are visible downstream.
    st.dod_flags = dodflags;

    summary.strata(s) = st;
end

end

function T = local_diff_varcomp_table(idvc,trvc,bsig,bsdim,bdim,ndim,...
    dim_names,ciperc)
%Posterior summaries of the reported difference-score components (Greek
%notation). idvc/trvc are draws x 2 x 2 person/trial event variance-covariances;
%bsig is draws x 2 per-event log-residual SD at z=0; bsdim/bdim are draws x KDIM
%event-specific scale/mean dimension slopes (odd columns event 1, even event 2).
[add,getT] = local_vc_accumulator(ciperc);

add('Difference universe-score variance (conditional on dimensions)',...
    'sigma^2_p_delta', idvc(:,1,1)+idvc(:,2,2)-2*idvc(:,2,1));
add('Person variance, event 1','sigma^2_p,1', idvc(:,1,1));
add('Person variance, event 2','sigma^2_p,2', idvc(:,2,2));
add('Person covariance (events)','sigma_p,12', idvc(:,2,1));
add('Trial difference variance','sigma^2_i_delta',...
    trvc(:,1,1)+trvc(:,2,2)-2*trvc(:,2,1));
add('Residual variance at z = 0, event 1','sigma^2_e,1(0)', exp(bsig(:,1)).^2);
add('Residual variance at z = 0, event 2','sigma^2_e,2(0)', exp(bsig(:,2)).^2);

%dimension-term names for the event-specific slope labels
terms = local_slopelabels(ndim,dim_names);
for k = 1:size(bdim,2)
    ev = 2 - mod(k,2);            % odd column -> event 1, even -> event 2
    tm = terms{ceil(k/2)};
    add(sprintf('Mean slope, event %d, %s',ev,tm),...
        sprintf('beta_%d,%d',ev,ceil(k/2)), bdim(:,k));
end
for k = 1:size(bsdim,2)
    ev = 2 - mod(k,2);
    tm = terms{ceil(k/2)};
    add(sprintf('Log-residual slope, event %d, %s',ev,tm),...
        sprintf('beta_sigma_%d,%d',ev,ceil(k/2)), bsdim(:,k));
end

T = getT();
end

function T = local_rescor_varcomp(rescor,errvc,ciperc)
%Residual-correlation components for the CONCURRENT variants (analyses 17/18): the
%estimated between-event residual correlation rho_e,12 and the residual covariance
%at z = 0 (errvc(:,2,1)). The covariance is z-dependent through the per-event SDs;
%it is reported at z = 0 here, while the reliability surface uses the z-conditioned
%value rho_e,12 * sigma_e,1(z) * sigma_e,2(z).
ciedge = (1-ciperc)/2;
rc = rescor(:);
cov0 = errvc(:,2,1);
name = {'Residual correlation (events)';'Residual covariance at z = 0'};
sym  = {'rho_e,12';'sigma_e,12(0)'};
draws = {rc; cov0};
ptv = zeros(2,1); llv = zeros(2,1); ulv = zeros(2,1);
for k = 1:2
    ptv(k) = mean(draws{k});
    llv(k) = quantile(draws{k},ciedge);
    ulv(k) = quantile(draws{k},1-ciedge);
end
T = table(name,sym,ptv,llv,ulv,'VariableNames',...
    {'component','symbol','estimate','ci_lower','ci_upper'});
end

function T = local_cutrho_varcomp(cutrho,ciperc)
%Copula residual correlation for the MODULAR CUT-COPULA designs (analyses 19/20):
%the Gaussian-copula correlation rho_e of the two events' PIT-transformed
%residuals, summarized over the RETAINED stage-1 draws only (cutrho is already
%the retained set, keeping this row per-draw consistent with the surface and the
%rest of the table). This is the copula correlation on the latent normal scale,
%not the observed-scale residual correlation; the observed-scale residual
%covariance is deliberately NOT reported here for the same z = 0 half-visibility
%reason local_gamma_difftrt_varcomp_table_ls documents. The symbol is
%DELIBERATELY DISTINCT from the Gaussian concurrent rho_e,12 row: under the
%identity link that row's rescor IS the observed-scale residual correlation,
%while this one is the latent-normal copula correlation - two different
%estimands that would otherwise collide under one symbol in the concatenated
%export whenever the same analysis number is run under both families
%(2026-08-11 code review). The cut-posterior status of this estimate is stated
%by summary.copula_labels and the provenance lines, not restated per row.
ciedge = (1-ciperc)/2;
rc = cutrho(:);
name = {'Copula residual correlation (events)'};
sym  = {'rho_e,12(copula)'};
ptv = mean(rc);
llv = quantile(rc,ciedge);
ulv = quantile(rc,1-ciedge);
T = table(name,sym,ptv,llv,ulv,'VariableNames',...
    {'component','symbol','estimate','ci_lower','ci_upper'});
end

function local_warn_copula_invalidity(flags,label)
%Numerical-invalidity disclosure for the modular cut-copula stage. The stage's
%contract is diagnose-not-hide: reliability is still computed, so the user must
%be told when stage-2 draws were numerically suspect. Calibration flags are
%descriptive model checks and do not warn; they ride on summary.strata(s)
%.copula_flags and the provenance lines instead.
if isempty(flags) || flags.n_invalid == 0
    return
end
warning('dynrel:copulainvalidity',...
    ['Stratum ''%s'': %d of %d retained stage-2 draws are flagged '...
    'numerically invalid (excessive PIT clipping or copula posterior mass '...
    'at the rho support bound). Reliability is still reported, but rho_e '...
    'and the residual-coupling terms from those draws are numerically '...
    'suspect; inspect summary.strata(s).copula_flags and consider more '...
    'stage-1 draws or a check of the margin fit.'],...
    char(string(label)),flags.n_invalid,flags.n_draws);
end

function T = local_rho_spread_varcomp(REL,s,ciperc)
%Between-person spread of the PER-SUBJECT residual correlation (analyses 21/22,
%Stage 4e). Each participant has their own residual correlation
%rho[s] = tanh(rescor_mu + sd_rescor*z_rescor[s]); this row summarizes the standard
%deviation of those per-subject correlations ACROSS subjects (computed per posterior
%draw, then summarized). The population-AVERAGE correlation itself is reported by the
%companion residual-correlation row (local_rescor_varcomp, which is fed mean(rho)).
ciedge = (1-ciperc)/2;
rho_ss = cell2mat(REL.out.rho_ss{s});   % draws x NSUB
%std across a single subject is undefined; report NaN rather than a spuriously
%precise 0 so a one-subject stratum reads as unestimable, not zero-variance.
if size(rho_ss,2) < 2
    sd_draws = nan(size(rho_ss,1),1);
else
    sd_draws = std(rho_ss,0,2);         % between-person SD of rho per draw
end
name = {'Between-person SD of residual correlation'};
sym  = {'SD(rho_e,s)'};
ptv = mean(sd_draws);
llv = quantile(sd_draws,ciedge);
ulv = quantile(sd_draws,1-ciedge);
T = table(name,sym,ptv,llv,ulv,'VariableNames',...
    {'component','symbol','estimate','ci_lower','ci_upper'});
end

function ax = local_axis(zvals,ngrid,dimname,label)
%Build a grid axis over the observed range of a standardized dimension. When
%the stratum has (near) no variation in this dimension (min == max), widen the
%axis to a small symmetric window so the surface still renders, and warn that
%the dimension is effectively constant in this stratum (its scale slope is not
%identified there).
lo = min(zvals);
hi = max(zvals);
if ~isfinite(lo) || ~isfinite(hi) || (hi - lo) <= eps(max(1,abs(hi)))
    z0 = lo;
    if ~isfinite(z0); z0 = 0; end
    lo = z0 - 0.5;
    hi = z0 + 0.5;
    warning('dynrel:constantdimension',...
        ['Dimension ''%s'' has effectively no variation within stratum '...
        '''%s''; its reliability surface over this dimension is shown across '...
        'a nominal range and should be interpreted with caution (the scale '...
        'slope is not identified without within-stratum variation).'],...
        char(string(dimname)),char(string(label)));
end
ax = linspace(lo,hi,ngrid);
end

function T = local_varcomp_table(sig_p,sig_i,log0,bsig,bmean,sig_dp,rho,...
    ndim,dim_names,ciperc)
%Posterior summaries of the reported components, labeled in Greek notation.
[add,getT] = local_vc_accumulator(ciperc);

add('Between-person variance (universe score, conditional on dimensions)',...
    'sigma^2_p', sig_p.^2);
add('Trial main-effect variance','sigma^2_i', sig_i.^2);
add('Residual variance at z = 0','sigma^2_pi,e(0)', exp(log0).^2);
add('Person scale SD','sigma_delta_p', sig_dp);
add('Location-scale correlation','rho', rho);

%dimension slope blocks (mean b and scale b_sigma)
slopelabels = local_slopelabels(ndim,dim_names);
for k = 1:numel(slopelabels)
    add(['Mean slope, ' slopelabels{k}], sprintf('beta_%d', k), bmean(:,k));
end
for k = 1:numel(slopelabels)
    add(['Log-residual slope, ' slopelabels{k}], sprintf('beta_sigma_%d', k),...
        bsig(:,k));
end

T = getT();
end

function ssrel = local_subj_ssrel_table(REL,s,idvc,trvc,bsig,bsdim,ndim,ciperc,...
    rescor)
%Per-participant subject-level dynamic difference reliability for stratum s.
%The dimension predictors are between-person (one z per subject), so each
%participant's per-event log-residual is conditioned at THEIR OWN z and then
%psyrat_ssrel_diff is reused unchanged (called once for generalizability =
%phi_rho,s, once for dependability = phi_delta,s). er_var_ss1/er_var_ss2 are the
%z = 0 per-subject log-residuals (b_sigma[e] + r_1_s{e}) indexed by id2; the
%dimension contribution b_sigma_dim[e] . z_terms(z_s) is added per subject.
%
%rescor (optional): the residual correlation for the CONCURRENT subject-level
%variants. Two shapes are accepted: a draws x 1 column is a single POPULATION
%correlation shared across subjects (analysis 19); a draws x NSUB matrix is a
%PER-SUBJECT correlation (analysis 21, Stage 4e), one column per subject aligned to
%the er_var_ss columns. In both cases each subject's per-event residual covariance
%is corr .* exp(er_ss1_z) .* exp(er_ss2_z) at the subject's own z-conditioned SDs;
%when absent/empty the residual covariance is 0 (non-concurrent, byte-identical to
%case 13).
if nargin < 9; rescor = []; end
%one-facet difference coefficient: person (bp) + trial (bt) blocks, with each
%subject's residual covariance carried on the wp_cov_ss channel (wp_cov = 0). The
%shared z-conditioning, residual covariance, and table assembly are in
%local_subj_ssrel_core; this wrapper only binds the coefficient.
coeff = @(est,ssinfo,e1,e2,zc) psyrat_ssrel_diff('bp',idvc,'bt',trvc,...
    'er_var',bsig,'idtable',ssinfo,'CI',ciperc,'est',est,...
    'er_var_ss1',e1,'er_var_ss2',e2,'wp_cov_ss',zc,'wp_cov',0);
ssrel = local_subj_ssrel_core(REL,s,bsdim,ndim,rescor,coeff);
end

function ssrel = local_subj_ssrel_table_trt(REL,s,idvc,trvc,ocvc,tidvc,oidvc,...
    tovc,bsig,bsdim,ndim,reltype,nocc,ciperc,rescor)
%Per-participant subject-level trial + occasion two-facet dynamic difference
%reliability for stratum s (Stage 4b). The two-facet analogue of
%local_subj_ssrel_table: each participant's per-event log-residual is conditioned
%at THEIR OWN z, then psyrat_ssrel_diff_trt is called once for generalizability
%(phi_rho,s) and once for dependability (phi_delta,s). The six event covariance
%blocks map to the psyrat_diffrel_trt inputs (id->bp, tid->bpi, oid->bpo, trl->bt,
%occ->bo, to->boi). The occasion n' is the shared reference value (nocc), matching
%the typical-person surface; each subject's observed occasions are in ssinfo.
%
%rescor (optional): the residual correlation for the CONCURRENT subject-level
%variants. A draws x 1 column is a single POPULATION correlation (analysis 20); a
%draws x NSUB matrix is a PER-SUBJECT correlation (analysis 22, Stage 4e), one
%column per subject. Each subject's per-event residual covariance is
%corr .* exp(er_ss1_z) .* exp(er_ss2_z); absent/empty -> 0 (non-concurrent,
%byte-identical to case 16).
if nargin < 15; rescor = []; end
%two-facet difference coefficient: the six event (co)variance blocks + occasion
%n', with each subject's residual covariance carried on the er_cov_ss channel
%(er_cov = 0). The shared z-conditioning, residual covariance, and table assembly
%are in local_subj_ssrel_core; this wrapper only binds the coefficient.
blockargs = {'bp',idvc,'bpi',tidvc,'bpo',oidvc,'bt',trvc,'bo',ocvc,'boi',tovc};
coeff = @(est,ssinfo,e1,e2,zc) psyrat_ssrel_diff_trt(blockargs{:},...
    'er_var',bsig,'idtable',ssinfo,'reltype',reltype,'nocc',nocc,'CI',ciperc,...
    'est',est,'er_var_ss1',e1,'er_var_ss2',e2,'er_cov',0,'er_cov_ss',zc);
ssrel = local_subj_ssrel_core(REL,s,bsdim,ndim,rescor,coeff);
end

function ssrel = local_subj_ssrel_core(REL,s,bsdim,ndim,rescor,coefffun)
%Shared body for the subject-level dynamic difference reliability tables
%(local_subj_ssrel_table and local_subj_ssrel_table_trt). Conditions each
%participant's per-event log-residual at THEIR OWN standardized dimension value(s),
%builds the per-subject residual covariance from rescor, then calls a model-specific
%difference coefficient once for generalizability and once for dependability.
%coefffun(est,ssinfo,er_ss1_z,er_ss2_z,zcov) wraps that coefficient
%(psyrat_ssrel_diff for one-facet, psyrat_ssrel_diff_trt for the two-facet trial +
%occasion variant) with the fixed (co)variance blocks already bound; only the
%z-conditioned per-event residuals and covariance vary per call.
%
%rescor (optional): the residual correlation for the CONCURRENT subject-level
%variants. A draws x 1 column is a single POPULATION correlation shared across
%subjects; a draws x NSUB matrix is a PER-SUBJECT correlation (Stage 4e), one column
%per subject aligned to the er_var_ss columns. Each subject's per-event residual
%covariance is corr .* exp(er_ss1_z) .* exp(er_ss2_z) at the subject's own
%z-conditioned SDs; when absent/empty the residual covariance is 0 (non-concurrent,
%byte-identical to the non-concurrent cases 13/16).
er_ss1 = cell2mat(REL.out.er_var_ss1{s});   % draws x NSUB
er_ss2 = cell2mat(REL.out.er_var_ss2{s});   % draws x NSUB
ssinfo = REL.out.ssinfo{s};
ndraws = size(er_ss1,1);
nsub   = size(er_ss1,2);

%the loop below writes column ssinfo.id2(c) of the draws x NSUB residual
%matrices, so id2 must be a permutation of 1:nsub (one row per subject, each a
%valid column index). Assert to fail loudly rather than silently skip or
%overwrite a subject's residual if a future change breaks the seqid invariant.
assert(isequal(sort(ssinfo.id2(:)'),1:nsub),...
    'psyrat_dynrel_summary:ssinfoindex',...
    'ssinfo.id2 must be a permutation of 1:NSUB.');

%event-specific scale slopes: odd cols -> event 1, even -> event 2 (Xdim layout)
ev1 = bsdim(:,1:2:end);
ev2 = bsdim(:,2:2:end);

%z-condition each subject at their own standardized dimension value(s)
er_ss1_z = er_ss1;
er_ss2_z = er_ss2;
for c = 1:height(ssinfo)
    col = ssinfo.id2(c);
    if ndim == 1
        zterms = ssinfo.z1(c);
    else
        zterms = [ssinfo.z1(c); ssinfo.z2(c); ssinfo.z1(c)*ssinfo.z2(c)];
    end
    er_ss1_z(:,col) = er_ss1(:,col) + ev1 * zterms(:);
    er_ss2_z(:,col) = er_ss2(:,col) + ev2 * zterms(:);
end

%per-subject residual covariance: 0 (non-concurrent); the population rescor
%(draws x 1, broadcast across subjects); or each subject's own rho (draws x NSUB,
%per-subject variant) -- times each subject's own z-conditioned per-event SDs.
if isempty(rescor)
    zcov = zeros(ndraws,nsub);
elseif size(rescor,2) == 1
    zcov = rescor(:) .* exp(er_ss1_z) .* exp(er_ss2_z);
else
    zcov = rescor .* exp(er_ss1_z) .* exp(er_ss2_z);
end

%gen = generalizability phi_rho,s; dep = dependability phi_delta,s
gen = coefffun('gen',ssinfo,er_ss1_z,er_ss2_z,zcov);
dep = coefffun('dep',ssinfo,er_ss1_z,er_ss2_z,zcov);

%clean per-participant table: relative (g) / absolute (d) ICC and SEM.
ssrel = ssinfo;
ssrel.gen_pt = gen.rel_pt; ssrel.gen_ll = gen.rel_ll; ssrel.gen_ul = gen.rel_ul;
ssrel.dep_pt = dep.rel_pt; ssrel.dep_ll = dep.rel_ll; ssrel.dep_ul = dep.rel_ul;
ssrel.iccg_pt = gen.icc_pt; ssrel.iccg_ll = gen.icc_ll; ssrel.iccg_ul = gen.icc_ul;
ssrel.iccd_pt = dep.icc_pt; ssrel.iccd_ll = dep.icc_ll; ssrel.iccd_ul = dep.icc_ul;
ssrel.sem_gen_pt = gen.sem_pt; ssrel.sem_gen_ll = gen.sem_ll; ssrel.sem_gen_ul = gen.sem_ul;
ssrel.sem_dep_pt = dep.sem_pt; ssrel.sem_dep_ll = dep.sem_ll; ssrel.sem_dep_ul = dep.sem_ul;
end

function ssrel = local_subj_ssrel_table_nondiff(REL,s,sig_p,log0,bsig,sig_i,...
    ndim,ciperc)
%Per-participant subject-level ONE-FACET non-difference dynamic reliability for
%stratum s (analysis 26). The dimension predictors are between-person (one z per
%subject), so each participant's log-residual is conditioned at THEIR OWN z and
%psyrat_ssrel is reused unchanged: once for generalizability (G_s, gcoeff = 2)
%and once for dependability (phi_s, gcoeff = 1). The between-person SD (sig_p),
%population residual intercept (log0 = Intercept_sigma), and trial main effect
%(sig_i) are shared; only the residual varies per subject and z.
ssinfo = REL.out.ssinfo{s};
indsd  = cell2mat(REL.out.ind_sdlog{s});   % draws x NSUB (per-subject log-resid RE)
nsub   = size(indsd,2);
%psyrat_ssrel broadcasts idtable.trls positionally against the residual columns,
%so the id2 column index must line up with the ssinfo row order. local_dynrel_ssinfo
%builds id2 = sort(unique(id2)) = 1:NSUB, so row r <-> id2 r <-> column r. Assert
%to fail loudly if a future change breaks that invariant.
assert(isequal(ssinfo.id2(:)',1:nsub),...
    'psyrat_dynrel_summary:ssinfoindex',...
    'ssinfo.id2 must equal 1:NSUB in row order.');

%z-condition each subject's log-residual at their own standardized dimension
%value(s): wp_ss(z_s) = ind_sdlog[s] + b_sigma . z_terms(z_s). Folding the
%population scale slope b_sigma into wp_ss lets psyrat_ssrel (which forms
%wp = exp(wp_pop + wp_ss)^2 with wp_pop = log0) evaluate each subject's residual
%at their own z without changing that shared calc.
wp_ss_z = indsd;
for c = 1:height(ssinfo)
    col = ssinfo.id2(c);
    wp_ss_z(:,col) = indsd(:,col) + bsig * local_zterms(ssinfo,c,ndim);
end

%gen = generalizability G_s (gcoeff = 2); dep = dependability phi_s (gcoeff = 1).
%bp/wp_ss are single-cell-wrapped so psyrat_ssrel's cell2mat returns them
%unchanged (bp = sig_p, wp_ss = the z-conditioned per-subject log-residuals).
gen = psyrat_ssrel('bp',{sig_p},'wp_pop',log0,'wp_ss',{wp_ss_z},...
    'idtable',ssinfo,'CI',ciperc,'i',sig_i,'gcoeff',2);
dep = psyrat_ssrel('bp',{sig_p},'wp_pop',log0,'wp_ss',{wp_ss_z},...
    'idtable',ssinfo,'CI',ciperc,'i',sig_i,'gcoeff',1);

ssrel = local_pack_ssrel_nondiff(ssinfo,gen,dep);
end

function ssrel = local_subj_ssrel_table_nondiff_gamma(REL,s,alpha0,bmean,...
    lognu0,bnu,sig_p,sig_i,ndim,ciperc)
%Per-participant subject-level ONE-FACET non-difference dynamic reliability for
%stratum s under the GAMMA / scaled-chi-square family (analysis 26). The gamma
%twin of local_subj_ssrel_table_nondiff.
%
%The Gaussian helper folds the population scale slope into a per-subject
%log-residual offset and calls psyrat_ssrel ONCE, because sigma_p and sigma_i are
%z-invariant there. Under gamma every observed-scale component carries
%exp(2*alpha(z)), so the between-person and trial components move with z too and
%each participant sits at their own z; psyrat_ssrel_dynrel_gamma therefore
%evaluates all of them per participant. The per-participant residual is derived
%from each person's own (u_p, v_p) rather than read off a log-residual parameter,
%which the gamma model does not have.
ssinfo = REL.out.ssinfo{s};
ind_bs = cell2mat(REL.out.ind_bs{s});   % draws x NSUB (person log-mean u_p)
ind_nu = cell2mat(REL.out.ind_nu{s});   % draws x NSUB (person log-nu  v_p)
nsub   = size(ind_bs,2);
%local_dynrel_ssinfo builds id2 = sort(unique(id2)) = 1:NSUB, so row r <-> id2 r
%<-> draw column r. Assert to fail loudly if a future change breaks that
%invariant and starts assigning per-person values to the wrong participants.
assert(isequal(ssinfo.id2(:)',1:nsub),...
    'psyrat_dynrel_summary:ssinfoindex',...
    'ssinfo.id2 must equal 1:NSUB in row order.');

[gen,dep] = psyrat_ssrel_dynrel_gamma('alpha0',alpha0,'b',bmean,...
    'lognu0',lognu0,'b_nu',bnu,'sig_p',sig_p,'sig_i',sig_i,...
    'ind_bs',ind_bs,'ind_nu',ind_nu,'idtable',ssinfo,'ndim',ndim,'CI',ciperc);

ssrel = local_pack_ssrel_nondiff(ssinfo,gen,dep);
end

function ssrel = local_subj_ssrel_table_nondiff_gamma_ls(REL,s,alpha0,bmean,...
    logsig0,bsig,sig_p,sig_i,ndim,ciperc)
%Per-participant subject-level ONE-FACET non-difference dynamic reliability for
%stratum s under the LOCATION-SCALE GAMMA family (analysis 26, gammascale = 2).
%The location-scale twin of local_subj_ssrel_table_nondiff_gamma.
%
%The log-nu twin must read TWO per-subject blocks (ind_bs and ind_nu) and derive
%the residual from them. Here the residual is a fitted parameter: ind_sdlog
%holds each person's log residual-SD offset, exactly as in the Gaussian
%subject-level path. ind_bs is not read at all - the participant's own log-mean
%does not enter the coefficient under this parameterization (see
%psyrat_ssrel_dynrel_gamma_ls). The estimator still stores it, for parity with
%the log-nu layout and as a diagnostic.
ssinfo = REL.out.ssinfo{s};
ind_sd = cell2mat(REL.out.ind_sdlog{s});   % draws x NSUB (person log-sigma w_p)
nsub   = size(ind_sd,2);
%local_dynrel_ssinfo builds id2 = sort(unique(id2)) = 1:NSUB, so row r <-> id2 r
%<-> draw column r. Assert to fail loudly if a future change breaks that
%invariant and starts assigning per-person values to the wrong participants.
assert(isequal(ssinfo.id2(:)',1:nsub),...
    'psyrat_dynrel_summary:ssinfoindex',...
    'ssinfo.id2 must equal 1:NSUB in row order.');

[gen,dep] = psyrat_ssrel_dynrel_gamma_ls('alpha0',alpha0,'b',bmean,...
    'logsig0',logsig0,'b_sigma',bsig,'sig_p',sig_p,'sig_i',sig_i,...
    'ind_sd',ind_sd,'idtable',ssinfo,'ndim',ndim,'CI',ciperc);

ssrel = local_pack_ssrel_nondiff(ssinfo,gen,dep);
end

function ssrel = local_subj_ssrel_table_nondiff_trt_gamma(REL,s,alpha0,bmean,...
    lognu0,bnu,sig_p,s_occ,s_trl,s_txp,s_oxp,s_txo,reltype,nocc,ciperc,ndim)
%Per-participant subject-level TWO-FACET non-difference dynamic reliability for
%stratum s under the GAMMA / scaled-chi-square family (analysis 27). The gamma
%twin of local_subj_ssrel_table_nondiff_trt.
%
%The Gaussian helper broadcasts the five crossed facet SDs as population values
%and varies only the residual per subject, because those SDs are z-invariant
%there. Under gamma every observed-scale component carries exp(2*alpha(z)), so
%all six move with z and each participant sits at their own z;
%psyrat_ssrel_dynrel_trt_gamma therefore recomputes them per participant. The
%per-participant residual is derived from each person's own (u_p, v_p) rather
%than read off a log-residual parameter, which the gamma model does not have.
ssinfo = REL.out.ssinfo{s};
ind_bs = cell2mat(REL.out.ind_bs{s});   % draws x NSUB (person log-mean u_p)
ind_nu = cell2mat(REL.out.ind_nu{s});   % draws x NSUB (person log-nu  v_p)
nsub   = size(ind_bs,2);
%the per-participant lookup indexes the draw columns by id2, so any row order is
%tolerated; the assert guards the id2/column correspondence so a future change
%cannot silently assign one participant's dispersion to another.
assert(isequal(sort(ssinfo.id2(:)'),1:nsub),...
    'psyrat_dynrel_summary:ssinfoindex',...
    'ssinfo.id2 must be a permutation of 1:NSUB.');

[gen,dep] = psyrat_ssrel_dynrel_trt_gamma('alpha0',alpha0,'b',bmean,...
    'lognu0',lognu0,'b_nu',bnu,'sig_p',sig_p,...
    'sig_occ',s_occ,'sig_trl',s_trl,'sig_occxid',s_oxp,...
    'sig_trlxid',s_txp,'sig_trlxocc',s_txo,...
    'ind_bs',ind_bs,'ind_nu',ind_nu,'idtable',ssinfo,'ndim',ndim,...
    'reltype',reltype,'nocc',nocc,'CI',ciperc);

ssrel = local_pack_ssrel_nondiff(ssinfo,gen,dep);
end

function ssrel = local_subj_ssrel_table_nondiff_trt_gamma_ls(REL,s,alpha0,...
    bmean,logsig0,bsig,sig_p,s_occ,s_trl,s_txp,s_oxp,s_txo,reltype,nocc,...
    ciperc,ndim)
%Per-participant subject-level TWO-FACET non-difference dynamic reliability for
%stratum s under the LOCATION-SCALE GAMMA family (analysis 27, gammascale = 2).
%The location-scale twin of local_subj_ssrel_table_nondiff_trt_gamma, and the
%two-facet sibling of local_subj_ssrel_table_nondiff_gamma_ls.
%
%The log-nu twin must read TWO per-subject blocks (ind_bs and ind_nu) and derive
%the residual from them, carrying an exp(2*V_w) within-person factor. Here the
%residual is a fitted parameter: ind_sdlog holds each person's log residual-SD
%offset, exactly as in the Gaussian subject-level path. ind_bs is not read at all
%- the participant's own log-mean does not enter the coefficient under this
%parameterization (see psyrat_ssrel_dynrel_trt_gamma_ls), and for this design the
%estimator does not store it either.
%
%The five crossed facet SDs are still passed through and still recomputed per
%participant, for the same reason the log-nu twin gives: under gamma every
%observed-scale SIGNAL component carries exp(2*alpha(z)), so all six move with z
%and each participant sits at their own z. That part is unchanged by the
%parameterization; only the residual differs.
ssinfo = REL.out.ssinfo{s};
ind_sd = cell2mat(REL.out.ind_sdlog{s});   % draws x NSUB (person log-sigma w_p)
nsub   = size(ind_sd,2);
%the per-participant lookup indexes the draw columns by id2, so any row order is
%tolerated; the assert guards the id2/column correspondence so a future change
%cannot silently assign one participant's residual to another.
assert(isequal(sort(ssinfo.id2(:)'),1:nsub),...
    'psyrat_dynrel_summary:ssinfoindex',...
    'ssinfo.id2 must be a permutation of 1:NSUB.');

[gen,dep] = psyrat_ssrel_dynrel_trt_gamma_ls('alpha0',alpha0,'b',bmean,...
    'logsig0',logsig0,'b_sigma',bsig,'sig_p',sig_p,...
    'sig_occ',s_occ,'sig_trl',s_trl,'sig_occxid',s_oxp,...
    'sig_trlxid',s_txp,'sig_trlxocc',s_txo,...
    'ind_sd',ind_sd,'idtable',ssinfo,'ndim',ndim,...
    'reltype',reltype,'nocc',nocc,'CI',ciperc);

ssrel = local_pack_ssrel_nondiff(ssinfo,gen,dep);
end

function ssrel = local_subj_ssrel_table_nondiff_trt(REL,s,sig_p,log0,bsig,...
    s_occ,s_trl,s_txp,s_oxp,s_txo,reltype,nocc,ciperc,ndim)
%Per-participant subject-level TWO-FACET non-difference dynamic reliability for
%stratum s (analysis 27). The two-facet analogue of
%local_subj_ssrel_table_nondiff: each participant's residual SD is conditioned at
%THEIR OWN z, then psyrat_ssrel_trt is reused once for generalizability
%(gcoeff = 2) and once for dependability (gcoeff = 1). Only the residual is
%per-subject; the five crossed facet SDs are population, so txp_ss/oxp_ss
%broadcast them across subjects (psyrat_ssrel_trt's per-subject path requires all
%three *_ss inputs), exactly as the static case-25 local_ssrel_trt_call.
ssinfo = REL.out.ssinfo{s};
indsd  = cell2mat(REL.out.ind_sdlog{s});   % draws x NSUB (per-subject log-resid RE)
nsub   = size(indsd,2);
%psyrat_ssrel_trt indexes the per-subject residual columns by idtable.id2, so it
%tolerates any row order; the assert still guards the id2/column correspondence.
assert(isequal(sort(ssinfo.id2(:)'),1:nsub),...
    'psyrat_dynrel_summary:ssinfoindex',...
    'ssinfo.id2 must be a permutation of 1:NSUB.');

%z-condition each subject's residual SD at their own z:
%err_ss(z_s) = exp(log0 + ind_sdlog[s] + b_sigma . z_terms(z_s)).
err_ss = exp(log0 + indsd);   % draws x NSUB at z = 0; overwritten per subject below
for c = 1:height(ssinfo)
    col = ssinfo.id2(c);
    err_ss(:,col) = exp(log0 + indsd(:,col) + bsig * local_zterms(ssinfo,c,ndim));
end
nsub_e = size(err_ss,2);

%the five crossed facet SDs are population (bo/bt/txp/oxp/txo); err = z = 0
%population residual (only psyrat_ssrel_trt's pop_errvar / non-ss fallback, both
%overridden by the per-subject err_ss path). obs per subject = idtable.trls
%(trials per occasion); nocc is the shared reference occasion count.
base = {'reltype',reltype,'bp',sig_p,'bo',s_occ,'bt',s_trl,...
    'txp',s_txp,'oxp',s_oxp,'txo',s_txo,'err',exp(log0),...
    'txp_ss',repmat(s_txp,1,nsub_e),'oxp_ss',repmat(s_oxp,1,nsub_e),...
    'err_ss',err_ss,'idtable',ssinfo,'CI',ciperc,'nocc',nocc};
gen = psyrat_ssrel_trt('gcoeff',2,base{:});
dep = psyrat_ssrel_trt('gcoeff',1,base{:});

ssrel = local_pack_ssrel_nondiff(ssinfo,gen,dep);
end

function zterms = local_zterms(ssinfo,c,ndim)
%The dimension design column for subject row c: [z1] for one dimension or
%[z1; z2; z1*z2] for two, a KDIM x 1 vector so b_sigma (draws x KDIM) * zterms is
%the draws x 1 scale contribution at that subject's z (matching the Xdim layout
%the estimator builds).
if ndim == 1
    zterms = ssinfo.z1(c);
else
    zterms = [ssinfo.z1(c); ssinfo.z2(c); ssinfo.z1(c)*ssinfo.z2(c)];
end
end

function ssrel = local_pack_ssrel_nondiff(ssinfo,gen,dep)
%Assemble the per-participant table from the generalizability (gen) and
%dependability (dep) psyrat_ssrel / psyrat_ssrel_trt outputs into the shared
%ssrel_table schema consumed by the viewer, plot, and headless report
%(gen_*/dep_*/iccg_*/iccd_*/sem_gen_*/sem_dep_*), matching the difference
%subject-level tables (local_subj_ssrel_core). Both calcs name their selected
%coefficient dep_*/icc_*/sem_* regardless of gcoeff, so the gcoeff = 2 call is
%the generalizability column and the gcoeff = 1 call the dependability column.
ssrel = ssinfo;
ssrel.gen_pt = gen.dep_pt; ssrel.gen_ll = gen.dep_ll; ssrel.gen_ul = gen.dep_ul;
ssrel.dep_pt = dep.dep_pt; ssrel.dep_ll = dep.dep_ll; ssrel.dep_ul = dep.dep_ul;
ssrel.iccg_pt = gen.icc_pt; ssrel.iccg_ll = gen.icc_ll; ssrel.iccg_ul = gen.icc_ul;
ssrel.iccd_pt = dep.icc_pt; ssrel.iccd_ll = dep.icc_ll; ssrel.iccd_ul = dep.icc_ul;
ssrel.sem_gen_pt = gen.sem_pt; ssrel.sem_gen_ll = gen.sem_ll; ssrel.sem_gen_ul = gen.sem_ul;
ssrel.sem_dep_pt = dep.sem_pt; ssrel.sem_dep_ll = dep.sem_ll; ssrel.sem_dep_ul = dep.sem_ul;
end

function ssrel = local_pack_ssrel_diff(ssinfo,gen,dep)
%Assemble the per-participant table from the generalizability (gen) and
%dependability (dep) psyrat_ssrel_diff outputs into the shared ssrel_table schema
%(gen_*/dep_*/iccg_*/iccd_*/sem_gen_*/sem_dep_*).
%
%This is the DIFFERENCE-score counterpart of local_pack_ssrel_nondiff. It reads
%rel_* because psyrat_ssrel_diff names its selected coefficient rel_*, while
%psyrat_ssrel names its dep_* regardless of gcoeff.
%
%CORRECTION, and read this before merging the two packers. An earlier version of
%this comment justified keeping them separate by claiming that reading dep_* here
%"would silently pick up the wrong column - or error, if the est = 'gen' call never
%created it". Neither hazard exists: psyrat_ssrel_diff aliases dep_pt/dep_ll/dep_ul
%to rel_pt/rel_ll/rel_ul unconditionally, outside any est branch, so gen.dep_pt is
%an exact copy of gen.rel_pt and is always created. psyrat_ssrel_trt aliases the
%same way. So a single dep_*-reading packer would in fact be numerically correct
%for every callee of both packers, and the stated reason does not support the
%stated conclusion. Whether to merge them is a live question, not a closed one;
%it is left alone here because it is a refactor, not a fix.
%
%The inline packing at the end of local_subj_ssrel_core is the same mapping; this
%function exists so the gamma difference path does not have to route through that
%local's Gaussian-specific z-conditioning.
ssrel = ssinfo;
ssrel.gen_pt = gen.rel_pt; ssrel.gen_ll = gen.rel_ll; ssrel.gen_ul = gen.rel_ul;
ssrel.dep_pt = dep.rel_pt; ssrel.dep_ll = dep.rel_ll; ssrel.dep_ul = dep.rel_ul;
ssrel.iccg_pt = gen.icc_pt; ssrel.iccg_ll = gen.icc_ll; ssrel.iccg_ul = gen.icc_ul;
ssrel.iccd_pt = dep.icc_pt; ssrel.iccd_ll = dep.icc_ll; ssrel.iccd_ul = dep.icc_ul;
ssrel.sem_gen_pt = gen.sem_pt; ssrel.sem_gen_ll = gen.sem_ll; ssrel.sem_gen_ul = gen.sem_ul;
ssrel.sem_dep_pt = dep.sem_pt; ssrel.sem_dep_ll = dep.sem_ll; ssrel.sem_dep_ul = dep.sem_ul;
end

function T = local_subj_scale_varcomp_gamma(corid,ciperc)
%Person scale-RE variance components for the GAMMA subject-level dynamic
%difference (analysis 13, location-scale). The gamma twin of
%local_subj_scale_varcomp.
%
%Same reported quantities, different source. The Gaussian path stores the joint
%cholesky L_id and reconstructs the correlations from it; the gamma case-13
%extraction stores the assembled cor_id (draws x 4 x 4) instead and does not store
%L_id at all, so the correlations are read directly. Deriving them from a cholesky
%that is not in REL.out would fail loudly rather than silently, but only against a
%live fit - hence a separate function rather than a widened one.
%WHY ONLY ONE ROW. This block is APPENDED to the table built by
%local_gamma_diff_varcomp_table_ls, which is already built from the same sdid and
%corid and already reports four of the five quantities this function used to add:
%  sdid(:,3)    -> 'Person log residual-SD SD, event 1'  (s_sd,1)
%  sdid(:,4)    -> 'Person log residual-SD SD, event 2'  (s_sd,2)
%  corid(:,1,3) -> 'rho_ps,1'
%  corid(:,2,4) -> 'rho_ps,2'
%Emitting them again under second names (sigma_delta_p,*, rho_p*,sigma*) put every
%one of those quantities in the viewer uitable and the CSV/XLSX export TWICE, with
%bit-identical estimates and intervals but different labels, which reads as four
%extra parameters rather than four repeats. Nothing downstream deduplicates.
%
%The justification removed with them was that this mirrors the Gaussian
%subject-level path. It does not: local_diff_varcomp_table never receives sd_id and
%emits none of these rows, so there the twin adds five genuinely absent rows. Here
%only the cross-event scale correlation is absent, so only it is added.
[add,getT] = local_vc_accumulator(ciperc);

%The scale-scale coupling ACROSS events - the one person-scale quantity the
%location-scale difference table does not already carry. It is estimated and
%reported but does NOT enter the coefficient: under location-scale the residual
%never touches mu, so a person mean<->log-sigma correlation cannot affect the
%conditional variance. Reporting it is how a user sees whether the two events'
%precisions travel together in this sample.
add('Person scale correlation (events)','rho_sigma1,sigma2', reshape(corid(:,3,4),[],1));

T = getT();
end

function T = local_subj_scale_varcomp(sdid,Lid,ciperc)
%Person scale-RE variance components from the joint 4-coefficient person block
%(case 13). sdid is draws x 4 (1-2 location SDs, 3-4 scale/log-residual SDs); Lid
%is the draws x 4 x 4 joint cholesky factor. The reported correlations are the
%person location-scale couplings (the two-event generalization of the case-11
%single-event location-scale correlation) plus the scale-scale correlation across
%events.
[add,getT] = local_vc_accumulator(ciperc);

%per-event person scale (log-residual) SDs
add('Person scale SD, event 1','sigma_delta_p,1', reshape(sdid(:,3),[],1));
add('Person scale SD, event 2','sigma_delta_p,2', reshape(sdid(:,4),[],1));

%correlations from the joint cholesky (C = L*L'): location-scale within event and
%scale-scale across events. L is lower-triangular, so each entry is a short dot
%product of the relevant rows.
Lel = @(i,j) reshape(Lid(:,i,j),[],1);
cor13 = Lel(1,1).*Lel(3,1);
cor24 = Lel(2,1).*Lel(4,1) + Lel(2,2).*Lel(4,2);
cor34 = Lel(3,1).*Lel(4,1) + Lel(3,2).*Lel(4,2) + Lel(3,3).*Lel(4,3);
add('Location-scale correlation, event 1','rho_p1,sigma1', cor13);
add('Location-scale correlation, event 2','rho_p2,sigma2', cor24);
add('Person scale correlation (events)','rho_sigma1,sigma2', cor34);

T = getT();
end

function nocc = local_resolve_nocc(nocc_in,REL,s)
%Resolve the occasion n' for the two-facet surface. [] -> 1 (owner default);
%'observed' -> this stratum's number of occasions (REL.out.nocc(s)); a numeric
%scalar -> that value. Always at least 1.
if isempty(nocc_in)
    nocc = 1;
elseif (ischar(nocc_in) || isstring(nocc_in)) && ...
        strcmpi(char(nocc_in),'observed')
    nocc = REL.out.nocc(s);
elseif isnumeric(nocc_in) && isscalar(nocc_in)
    nocc = nocc_in;
else
    nocc = 1;
end
nocc = max(1,round(nocc));
end

function T = local_gamma_varcomp_table(alpha0,lognu0,bmean,bnu,sig_p,sig_v,...
    sig_i,rho_pv,ndim,dim_names,ciperc)
%Posterior summaries of the reported gamma / scaled-chi-square dynrel components
%(one-facet). The gamma family estimates LOG-scale parameters, so the log-mean /
%log-nu intercepts, the person (mean, log-nu) SDs and their correlation, the
%trial log-mean SD, and the two dimension-slope blocks are reported on the log
%scale, alongside the OBSERVED-score variance components at z = 0 (via
%psyrat_gamma_varcomps, estimand #2) for interpretability. Reporting s_v / rho_pv
%/ the log-nu slopes matters: the dispersion side is weakly identified at low
%trial counts (mean-coupled residual; see SCIENTIFIC_FORMULA_AUDIT.md section 18).
[add,getT] = local_vc_accumulator(ciperc);

alpha0 = alpha0(:); lognu0 = lognu0(:);
sig_p = sig_p(:); sig_v = sig_v(:); sig_i = sig_i(:); rho_pv = rho_pv(:);
cov_pv = rho_pv .* sig_p .* sig_v;

%observed-score variance components at z = 0 (estimand #2, marginal over person
%dispersion), the same converter used per z by psyrat_rel_dynrel_gamma.
vc0 = psyrat_gamma_varcomps(alpha0, sig_p, sig_i, exp(lognu0), sig_v, cov_pv);
add('Between-person variance at z = 0 (observed score)','sigma^2_p(0)', vc0.sigma_p2);
add('Trial main-effect variance at z = 0 (observed score)','sigma^2_i(0)', vc0.sigma_i2);
add('Residual variance at z = 0 (observed score)','sigma^2_pi,e(0)', vc0.sigma_pi_e2);

%log-scale parameters (what the model estimates)
add('Log-mean grand intercept (z = 0)','alpha(0)', alpha0);
add('Log-nu grand intercept (z = 0)','beta(0)', lognu0);
add('Population dispersion df at z = 0','nu(0)', exp(lognu0));
add('Person log-mean SD','s_p', sig_p);
add('Person log-nu SD','s_v', sig_v);
add('Person mean<->log-nu correlation','rho_pv', rho_pv);
add('Trial log-mean SD','s_i', sig_i);

%dimension slope blocks (log-mean b and log-nu b_nu)
slopelabels = local_slopelabels(ndim,dim_names);
for k = 1:numel(slopelabels)
    add(['Log-mean slope, ' slopelabels{k}], sprintf('b_%d', k), bmean(:,k));
end
for k = 1:numel(slopelabels)
    add(['Log-nu slope, ' slopelabels{k}], sprintf('b_nu_%d', k), bnu(:,k));
end

T = getT();
end

function T = local_gamma_ls_varcomp_table(alpha0,logsig0,bmean,bsig,sig_p,...
    sig_sd,sig_i,rho_ps,ndim,dim_names,ciperc)
%Posterior summaries of the reported LOCATION-SCALE gamma dynrel components
%(one-facet, gammascale = 2). The location-scale twin of
%local_gamma_varcomp_table.
%
%Every row that names the scale submodel changes, and that is the whole point of
%having a separate builder: gro_sds column 2 is s_sd (the person log residual-SD
%SD) here, not s_v; the scale intercept is a log SD, not a log nu; and the
%dimension slopes are b_sigma, not b_nu. Reusing the log-nu builder would print
%a table of correct numbers under wrong names - it would not error - which is
%why the two are separate rather than one function with a flag.
%
%rho_ps (the person mean <-> log-sigma correlation) IS reported, because it is
%estimated and interpretable. It is deliberately absent from the observed-scale
%conversion below: psyrat_gamma_varcomps_ls takes no covariance argument,
%because the residual never touches mu.
[add,getT] = local_vc_accumulator(ciperc);

alpha0 = alpha0(:); logsig0 = logsig0(:);
sig_p = sig_p(:); sig_sd = sig_sd(:); sig_i = sig_i(:); rho_ps = rho_ps(:);

%observed-score variance components at z = 0, marginal over the person
%residual-SD population - the same converter used per z by
%psyrat_rel_dynrel_gamma_ls.
vc0 = psyrat_gamma_varcomps_ls(alpha0, sig_p, sig_i, logsig0, sig_sd);
add('Between-person variance at z = 0 (observed score)','sigma^2_p(0)', vc0.sigma_p2);
add('Trial main-effect variance at z = 0 (observed score)','sigma^2_i(0)', vc0.sigma_i2);
add('Residual variance at z = 0 (observed score)','sigma^2_pi,e(0)', vc0.sigma_pi_e2);

%log-scale parameters (what the model estimates)
add('Log-mean grand intercept (z = 0)','alpha(0)', alpha0);
add('Log residual-SD grand intercept (z = 0)','log sigma(0)', logsig0);
add('Population residual SD at z = 0','sigma(0)', exp(logsig0));
add('Person log-mean SD','s_p', sig_p);
add('Person log residual-SD SD','s_sd', sig_sd);
add('Person mean<->log-scale correlation','rho_p,sigma', rho_ps);
add('Trial log-mean SD','s_i', sig_i);

%dimension slope blocks (log-mean b and log residual-SD b_sigma)
slopelabels = local_slopelabels(ndim,dim_names);
for k = 1:numel(slopelabels)
    add(['Log-mean slope, ' slopelabels{k}], sprintf('b_%d', k), bmean(:,k));
end
for k = 1:numel(slopelabels)
    add(['Log residual-SD slope, ' slopelabels{k}], sprintf('b_sigma_%d', k),...
        bsig(:,k));
end

T = getT();
end

function T = local_gamma_diff_varcomp_table(bb,bnu,bdim,bnudim,sdid,sdtrl,...
    corid,cori,ndim,dim_names,ciperc)
%Posterior summaries of the reported gamma / scaled-chi-square DYNAMIC DIFFERENCE
%components (analysis 12). Structure follows local_gamma_varcomp_table (log-scale
%parameters plus OBSERVED-score components at z = 0) but is reported PER EVENT
%plus the difference contrasts, mirroring local_diff_varcomp_table.
%
%bb/bnu are draws x 2 per-event log-mean / log-nu intercepts at z = 0; bdim/bnudim
%are draws x KDIM event-crossed slope blocks (ODD columns event 1, EVEN event 2);
%sdid is draws x 4 [mean_1, mean_2, lognu_1, lognu_2]; sdtrl draws x 2; corid is
%the draws x 4 x 4 person correlation; cori the cross-event trial correlation.
%
%NOTE the observed-score rows are evaluated AT z = 0 only. Unlike the Gaussian
%difference path, every one of them moves with z here (the log link makes each
%component carry exp(2*alpha(z) + V)), so these are a reference slice, not the
%whole story - the z-varying picture is the G/D surface itself.
[add,getT] = local_vc_accumulator(ciperc);

sp_a = sdid(:,1); sp_b = sdid(:,2);
sv_a = sdid(:,3); sv_b = sdid(:,4);
si_a = sdtrl(:,1); si_b = sdtrl(:,2);
cov_pv_a = corid(:,1,3) .* sp_a .* sv_a;
cov_pv_b = corid(:,2,4) .* sp_b .* sv_b;
cor_p = corid(:,1,2);
cori  = cori(:);

%observed-score components at z = 0 (estimand #2 per event), the same converters
%used per z by psyrat_rel_diffdynrel_gamma.
vc_a = psyrat_gamma_varcomps(bb(:,1), sp_a, si_a, exp(bnu(:,1)), sv_a, cov_pv_a);
vc_b = psyrat_gamma_varcomps(bb(:,2), sp_b, si_b, exp(bnu(:,2)), sv_b, cov_pv_b);
cc   = psyrat_gamma_crosscov(bb(:,1), sp_a, si_a, bb(:,2), sp_b, si_b, cor_p, cori);

add('Difference universe-score variance at z = 0 (observed score)',...
    'sigma^2_p_delta(0)', vc_a.sigma_p2 + vc_b.sigma_p2 - 2*cc.cov_p_obs);
add('Person variance at z = 0, event 1','sigma^2_p,1(0)', vc_a.sigma_p2);
add('Person variance at z = 0, event 2','sigma^2_p,2(0)', vc_b.sigma_p2);
add('Person covariance at z = 0 (events)','sigma_p,12(0)', cc.cov_p_obs);
add('Trial difference variance at z = 0','sigma^2_i_delta(0)',...
    vc_a.sigma_i2 + vc_b.sigma_i2 - 2*cc.cov_i_obs);
add('Trial variance at z = 0, event 1','sigma^2_i,1(0)', vc_a.sigma_i2);
add('Trial variance at z = 0, event 2','sigma^2_i,2(0)', vc_b.sigma_i2);
add('Trial covariance at z = 0 (events)','sigma_i,12(0)', cc.cov_i_obs);
add('Residual variance at z = 0, event 1','sigma^2_pi,e,1(0)', vc_a.sigma_pi_e2);
add('Residual variance at z = 0, event 2','sigma^2_pi,e,2(0)', vc_b.sigma_pi_e2);

%log-scale parameters (what the model estimates)
add('Log-mean intercept (z = 0), event 1','alpha_1(0)', bb(:,1));
add('Log-mean intercept (z = 0), event 2','alpha_2(0)', bb(:,2));
add('Log-nu intercept (z = 0), event 1','beta_1(0)', bnu(:,1));
add('Log-nu intercept (z = 0), event 2','beta_2(0)', bnu(:,2));
add('Dispersion df at z = 0, event 1','nu_1(0)', exp(bnu(:,1)));
add('Dispersion df at z = 0, event 2','nu_2(0)', exp(bnu(:,2)));
add('Person log-mean SD, event 1','s_p,1', sp_a);
add('Person log-mean SD, event 2','s_p,2', sp_b);
add('Person log-nu SD, event 1','s_v,1', sv_a);
add('Person log-nu SD, event 2','s_v,2', sv_b);
add('Person mean<->log-nu correlation, event 1','rho_pv,1', corid(:,1,3));
add('Person mean<->log-nu correlation, event 2','rho_pv,2', corid(:,2,4));
add('Cross-event person-mean correlation (log scale)','rho_p', cor_p);
add('Trial log-mean SD, event 1','s_i,1', si_a);
add('Trial log-mean SD, event 2','s_i,2', si_b);
add('Cross-event trial correlation (log scale)','rho_i', cori);

%event-crossed dimension slope blocks (odd columns event 1, even event 2)
slopelabels = local_slopelabels(ndim,dim_names);
for k = 1:numel(slopelabels)
    add(['Log-mean slope, event 1, ' slopelabels{k}], ...
        sprintf('b_dim,1_%d', k), bdim(:,2*k-1));
    add(['Log-mean slope, event 2, ' slopelabels{k}], ...
        sprintf('b_dim,2_%d', k), bdim(:,2*k));
end
for k = 1:numel(slopelabels)
    add(['Log-nu slope, event 1, ' slopelabels{k}], ...
        sprintf('b_nu_dim,1_%d', k), bnudim(:,2*k-1));
    add(['Log-nu slope, event 2, ' slopelabels{k}], ...
        sprintf('b_nu_dim,2_%d', k), bnudim(:,2*k));
end

T = getT();
end

function T = local_gamma_diff_varcomp_table_ls(bb,bsig,bdim,bsigdim,sdid,sdtrl,...
    corid,cori,ndim,dim_names,ciperc)
%Posterior summaries of the reported gamma LOCATION-SCALE DYNAMIC DIFFERENCE
%components (analysis 12, gammascale = 2). Location-scale sibling of
%local_gamma_diff_varcomp_table; same row structure, with the dispersion rows
%re-expressed on the residual-SD scale.
%
%bb/bsig are draws x 2 per-event log-mean / log residual-SD intercepts at z = 0;
%bdim/bsigdim are draws x KDIM event-crossed slope blocks (ODD columns event 1,
%EVEN event 2); sdid is draws x 4 [mean_1, mean_2, logsigma_1, logsigma_2];
%sdtrl draws x 2; corid the draws x 4 x 4 person correlation; cori the
%cross-event trial correlation.
%
%NOTE the observed-score rows are evaluated AT z = 0 only, for the same reason as
%the log-nu sibling: every component moves with z under the log link, so these
%are a reference slice rather than the whole story.
%
%TWO DELIBERATE DIFFERENCES FROM THE LOG-NU TABLE, both consequences of the
%decoupled residual. (a) There is no cov_pv: psyrat_gamma_varcomps_ls takes no
%such argument because the residual never touches mu. (b) The person
%mean<->scale correlations ARE still reported, because they are estimated and
%scientifically interesting, but the row label says explicitly that they do not
%enter the conversion - otherwise a reader would reasonably assume they do, as
%they do in the log-nu table.
[add,getT] = local_vc_accumulator(ciperc);

sp_a = sdid(:,1); sp_b = sdid(:,2);
ssd_a = sdid(:,3); ssd_b = sdid(:,4);
si_a = sdtrl(:,1); si_b = sdtrl(:,2);
cor_p = corid(:,1,2);
cori  = cori(:);

%observed-score components at z = 0, the same converters used per z by
%psyrat_rel_diffdynrel_gamma_ls. Note the asymmetry: the LOCATION-SCALE converter
%for the diagonals, the UNCHANGED cross-event helper for the off-diagonals,
%because all eight of its arguments are mean-submodel quantities.
vc_a = psyrat_gamma_varcomps_ls(bb(:,1), sp_a, si_a, bsig(:,1), ssd_a);
vc_b = psyrat_gamma_varcomps_ls(bb(:,2), sp_b, si_b, bsig(:,2), ssd_b);
cc   = psyrat_gamma_crosscov(bb(:,1), sp_a, si_a, bb(:,2), sp_b, si_b, cor_p, cori);

add('Difference universe-score variance at z = 0 (observed score)',...
    'sigma^2_p_delta(0)', vc_a.sigma_p2 + vc_b.sigma_p2 - 2*cc.cov_p_obs);
add('Person variance at z = 0, event 1','sigma^2_p,1(0)', vc_a.sigma_p2);
add('Person variance at z = 0, event 2','sigma^2_p,2(0)', vc_b.sigma_p2);
add('Person covariance at z = 0 (events)','sigma_p,12(0)', cc.cov_p_obs);
add('Trial difference variance at z = 0','sigma^2_i_delta(0)',...
    vc_a.sigma_i2 + vc_b.sigma_i2 - 2*cc.cov_i_obs);
add('Trial variance at z = 0, event 1','sigma^2_i,1(0)', vc_a.sigma_i2);
add('Trial variance at z = 0, event 2','sigma^2_i,2(0)', vc_b.sigma_i2);
add('Trial covariance at z = 0 (events)','sigma_i,12(0)', cc.cov_i_obs);
add('Residual variance at z = 0, event 1','sigma^2_pi,e,1(0)', vc_a.sigma_pi_e2);
add('Residual variance at z = 0, event 2','sigma^2_pi,e,2(0)', vc_b.sigma_pi_e2);

%log-scale parameters (what the model estimates)
add('Log-mean intercept (z = 0), event 1','alpha_1(0)', bb(:,1));
add('Log-mean intercept (z = 0), event 2','alpha_2(0)', bb(:,2));
add('Log residual-SD intercept (z = 0), event 1','logsigma_1(0)', bsig(:,1));
add('Log residual-SD intercept (z = 0), event 2','logsigma_2(0)', bsig(:,2));
add('Residual SD at z = 0, event 1','sigma_1(0)', exp(bsig(:,1)));
add('Residual SD at z = 0, event 2','sigma_2(0)', exp(bsig(:,2)));
add('Person log-mean SD, event 1','s_p,1', sp_a);
add('Person log-mean SD, event 2','s_p,2', sp_b);
add('Person log residual-SD SD, event 1','s_sd,1', ssd_a);
add('Person log residual-SD SD, event 2','s_sd,2', ssd_b);
add(['Person mean<->log-sigma correlation, event 1 '...
    '(reported; does not enter the conversion)'],'rho_ps,1', corid(:,1,3));
add(['Person mean<->log-sigma correlation, event 2 '...
    '(reported; does not enter the conversion)'],'rho_ps,2', corid(:,2,4));
add('Cross-event person-mean correlation (log scale)','rho_p', cor_p);
add('Trial log-mean SD, event 1','s_i,1', si_a);
add('Trial log-mean SD, event 2','s_i,2', si_b);
add('Cross-event trial correlation (log scale)','rho_i', cori);

%event-crossed dimension slope blocks (odd columns event 1, even event 2)
slopelabels = local_slopelabels(ndim,dim_names);
for k = 1:numel(slopelabels)
    add(['Log-mean slope, event 1, ' slopelabels{k}], ...
        sprintf('b_dim,1_%d', k), bdim(:,2*k-1));
    add(['Log-mean slope, event 2, ' slopelabels{k}], ...
        sprintf('b_dim,2_%d', k), bdim(:,2*k));
end
for k = 1:numel(slopelabels)
    add(['Log residual-SD slope, event 1, ' slopelabels{k}], ...
        sprintf('b_sigma_dim,1_%d', k), bsigdim(:,2*k-1));
    add(['Log residual-SD slope, event 2, ' slopelabels{k}], ...
        sprintf('b_sigma_dim,2_%d', k), bsigdim(:,2*k));
end

T = getT();
end

function T = local_gamma_difftrt_varcomp_table_ls(bb,bsig,bdim,bsigdim,sdid,...
    sdtrl,sdocc,sdtid,sdoid,sdto,corid,cort,coro,corpt,corpo,corot,...
    ndim,dim_names,ciperc)
%Posterior summaries of the reported gamma LOCATION-SCALE TWO-FACET DYNAMIC
%DIFFERENCE components (analysis 20, gammascale = 2). The two-facet twin of
%local_gamma_diff_varcomp_table_ls: same row structure, six crossed components
%per event instead of two, and the residual is the THREE-WAY term rather than
%the one-facet's pooled person x trial + error.
%
%NOTE the observed-score rows are evaluated AT z = 0 only, as in every gamma
%table here: every component moves with z under the log link, so these are a
%reference slice rather than the whole story.
%
%THE RESIDUAL ROWS NAME A DIFFERENT QUANTITY FROM THE ONE-FACET TABLE'S. There
%sigma^2_pi,e pools person x trial with measurement error because the design
%confounds them; here person x trial is its own reported component and the
%residual is sigma^2_pto,e, the three-way interaction plus unmodelled noise
%(Rast & Clayson, in press, page 4). Do not carry the one-facet label across.
%
%The residual cross-covariance is deliberately NOT reported here, and the reason
%is that at z = 0 this table can only see HALF of it. The pooled channel is the
%mean-surface three-way cross term PLUS the copula-derived observation term, and
%the latter needs a quadrature point that this table has no basis to choose: the
%surface evaluates it at the typical person, the per-participant calculator at
%each person's own expected means, and the two are different numbers. Emitting
%the mean-surface half alone under a "residual covariance" label would invite
%exactly the confusion the per-participant diagnostic table exists to prevent.
%The five reported cross-covariances have no such ambiguity. See
%psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls for the per-participant quantities.
[add,getT] = local_vc_accumulator(ciperc);

sp_a = sdid(:,1); sp_b = sdid(:,2);
ssd_a = sdid(:,3); ssd_b = sdid(:,4);
%facet name cross-over, stated once here and again at every call site: sdtid is
%trial x person and supplies psyrat_gamma_varcomps_trt_ls's s_pt argument, sdoid
%is occasion x person and supplies s_po.
st_a = sdtrl(:,1); st_b = sdtrl(:,2);
so_a = sdocc(:,1); so_b = sdocc(:,2);
spt_a = sdtid(:,1); spt_b = sdtid(:,2);
spo_a = sdoid(:,1); spo_b = sdoid(:,2);
sot_a = sdto(:,1);  sot_b = sdto(:,2);
cor_p = corid(:,1,2);
cort = cort(:); coro = coro(:); corpt = corpt(:);
corpo = corpo(:); corot = corot(:);

%observed-score components at z = 0, the same converters psyrat_rel_diffdynrel_
%trt_gamma_ls uses per z. The cross-event helper is the CONCURRENT wrapper with a
%zero residual covariance: at z = 0 the mean-surface cross terms are what this
%table reports, and the residual cross term is per-participant (see above).
vc_a = psyrat_gamma_varcomps_trt_ls(bb(:,1), sp_a, so_a, st_a, spo_a, spt_a, ...
    sot_a, bsig(:,1), ssd_a);
vc_b = psyrat_gamma_varcomps_trt_ls(bb(:,2), sp_b, so_b, st_b, spo_b, spt_b, ...
    sot_b, bsig(:,2), ssd_b);
cc = psyrat_gamma_crosscov_trt_rescor( ...
    bb(:,1), sp_a, so_a, st_a, spo_a, spt_a, sot_a, ...
    bb(:,2), sp_b, so_b, st_b, spo_b, spt_b, sot_b, ...
    cor_p, coro, cort, corpo, corpt, corot, 0);

add('Difference universe-score variance at z = 0 (observed score)',...
    'sigma^2_p_delta(0)', vc_a.sigma_p2 + vc_b.sigma_p2 - 2*cc.cov_p);
add('Person variance at z = 0, event 1','sigma^2_p,1(0)', vc_a.sigma_p2);
add('Person variance at z = 0, event 2','sigma^2_p,2(0)', vc_b.sigma_p2);
add('Person covariance at z = 0 (events)','sigma_p,12(0)', cc.cov_p);

vcpairs = { ...
    'Trial',            'i',    vc_a.sigma_t2,  vc_b.sigma_t2,  cc.cov_t; ...
    'Occasion',         'o',    vc_a.sigma_o2,  vc_b.sigma_o2,  cc.cov_o; ...
    'Trial x person',   'pt',   vc_a.sigma_pt2, vc_b.sigma_pt2, cc.cov_pt; ...
    'Occasion x person','po',   vc_a.sigma_po2, vc_b.sigma_po2, cc.cov_po; ...
    'Trial x occasion', 'io',   vc_a.sigma_ot2, vc_b.sigma_ot2, cc.cov_ot};
for k = 1:size(vcpairs,1)
    nm = vcpairs{k,1}; sym = vcpairs{k,2};
    add([nm ' difference variance at z = 0'], ...
        sprintf('sigma^2_%s_delta(0)',sym), ...
        vcpairs{k,3} + vcpairs{k,4} - 2*vcpairs{k,5});
    add([nm ' variance at z = 0, event 1'], sprintf('sigma^2_%s,1(0)',sym), vcpairs{k,3});
    add([nm ' variance at z = 0, event 2'], sprintf('sigma^2_%s,2(0)',sym), vcpairs{k,4});
    add([nm ' covariance at z = 0 (events)'], sprintf('sigma_%s,12(0)',sym), vcpairs{k,5});
end

add('Residual variance at z = 0, event 1 (three-way + dispersion)',...
    'sigma^2_pto,e,1(0)', vc_a.sigma_res2);
add('Residual variance at z = 0, event 2 (three-way + dispersion)',...
    'sigma^2_pto,e,2(0)', vc_b.sigma_res2);

%log-scale parameters (what the model estimates)
add('Log-mean intercept (z = 0), event 1','alpha_1(0)', bb(:,1));
add('Log-mean intercept (z = 0), event 2','alpha_2(0)', bb(:,2));
add('Log residual-SD intercept (z = 0), event 1','logsigma_1(0)', bsig(:,1));
add('Log residual-SD intercept (z = 0), event 2','logsigma_2(0)', bsig(:,2));
add('Residual SD at z = 0, event 1','sigma_1(0)', exp(bsig(:,1)));
add('Residual SD at z = 0, event 2','sigma_2(0)', exp(bsig(:,2)));
add('Person log-mean SD, event 1','s_p,1', sp_a);
add('Person log-mean SD, event 2','s_p,2', sp_b);
add('Person log residual-SD SD, event 1','s_sigma,1', ssd_a);
add('Person log residual-SD SD, event 2','s_sigma,2', ssd_b);
add('Cross-event person log-mean correlation','rho_p', cor_p);
%The WITHIN-event person mean<->log-sigma correlations. They are estimated (the
%LKJ person block) and scientifically interesting, but they enter NO component of
%this design: psyrat_gamma_varcomps_trt_ls takes no cov_pv argument because the
%residual never touches mu. Reported for the same reason the one-facet table
%reports them, with the label saying so - a reader would otherwise reasonably
%assume they feed the conversion, as they do on the log-nu path.
add('Person mean<->log-sigma correlation, event 1 (not in the conversion)',...
    'rho_ps,1', corid(:,1,3));
add('Person mean<->log-sigma correlation, event 2 (not in the conversion)',...
    'rho_ps,2', corid(:,2,4));

sdpairs = { ...
    'trial','i',st_a,st_b,cort; 'occasion','o',so_a,so_b,coro; ...
    'trial x person','pt',spt_a,spt_b,corpt; ...
    'occasion x person','po',spo_a,spo_b,corpo; ...
    'trial x occasion','io',sot_a,sot_b,corot};
for k = 1:size(sdpairs,1)
    nm = sdpairs{k,1}; sym = sdpairs{k,2};
    add(sprintf('Log-mean %s SD, event 1',nm), sprintf('s_%s,1',sym), sdpairs{k,3});
    add(sprintf('Log-mean %s SD, event 2',nm), sprintf('s_%s,2',sym), sdpairs{k,4});
    add(sprintf('Cross-event %s correlation (log scale)',nm), ...
        sprintf('rho_%s',sym), sdpairs{k,5});
end

%event-crossed dimension slope blocks (odd columns event 1, even event 2)
slopelabels = local_slopelabels(ndim,dim_names);
for k = 1:numel(slopelabels)
    add(['Log-mean slope, event 1, ' slopelabels{k}], ...
        sprintf('b_dim,1_%d', k), bdim(:,2*k-1));
    add(['Log-mean slope, event 2, ' slopelabels{k}], ...
        sprintf('b_dim,2_%d', k), bdim(:,2*k));
end
for k = 1:numel(slopelabels)
    add(['Log residual-SD slope, event 1, ' slopelabels{k}], ...
        sprintf('b_sigma_dim,1_%d', k), bsigdim(:,2*k-1));
    add(['Log residual-SD slope, event 2, ' slopelabels{k}], ...
        sprintf('b_sigma_dim,2_%d', k), bsigdim(:,2*k));
end

T = getT();
end

function T = local_gamma_trt_varcomp_table(alpha0,lognu0,bmean,bnu,sig_p,sig_v,...
    s_occ,s_trl,s_txp,s_oxp,s_txo,rho_pv,ndim,dim_names,ciperc)
%Posterior summaries of the reported gamma / scaled-chi-square dynrel components
%(trial + occasion two-facet). The gamma family estimates LOG-scale parameters,
%so the log-mean / log-nu intercepts, person (mean, log-nu) SDs and correlation,
%the five crossed log-mean SDs, and the two dimension-slope blocks are reported on
%the log scale, alongside the OBSERVED-score variance components at z = 0 (via
%psyrat_gamma_varcomps_trt, estimand #2). Reporting s_v / rho_pv / the log-nu
%slopes matters: the dispersion side is weakly identified at low trial counts
%(mean-coupled residual; see SCIENTIFIC_FORMULA_AUDIT.md section 18).
[add,getT] = local_vc_accumulator(ciperc);

alpha0 = alpha0(:); lognu0 = lognu0(:);
sig_p = sig_p(:); sig_v = sig_v(:);
s_occ = s_occ(:); s_trl = s_trl(:); s_txp = s_txp(:); s_oxp = s_oxp(:);
s_txo = s_txo(:); rho_pv = rho_pv(:);
cov_pv = rho_pv .* sig_p .* sig_v;

%observed-score components at z = 0 (estimand #2). Argument order follows
%psyrat_gamma_extract_trt: (alpha, s_p, s_o, s_t, s_po, s_pt, s_ot, nu, ...) with
%s_po = occasion x person (s_oxp), s_pt = trial x person (s_txp), s_ot = trial x
%occasion (s_txo).
vc0 = psyrat_gamma_varcomps_trt(alpha0, sig_p, s_occ, s_trl, s_oxp, s_txp, ...
    s_txo, exp(lognu0), sig_v, cov_pv);
add('Person variance at z = 0 (observed score)','sigma^2_p(0)', vc0.sigma_p2);
add('Occasion variance at z = 0 (observed score)','sigma^2_o(0)', vc0.sigma_o2);
add('Trial variance at z = 0 (observed score)','sigma^2_t(0)', vc0.sigma_t2);
add('Trial x person variance at z = 0','sigma^2_pt(0)', vc0.sigma_pt2);
add('Occasion x person variance at z = 0','sigma^2_po(0)', vc0.sigma_po2);
add('Trial x occasion variance at z = 0','sigma^2_ot(0)', vc0.sigma_ot2);
add('Residual variance at z = 0 (observed score)','sigma^2_pto,e(0)', vc0.sigma_res2);

%log-scale parameters (what the model estimates)
add('Log-mean grand intercept (z = 0)','alpha(0)', alpha0);
add('Log-nu grand intercept (z = 0)','beta(0)', lognu0);
add('Population dispersion df at z = 0','nu(0)', exp(lognu0));
add('Person log-mean SD','s_p', sig_p);
add('Person log-nu SD','s_v', sig_v);
add('Person mean<->log-nu correlation','rho_pv', rho_pv);
add('Occasion log-mean SD','s_o', s_occ);
add('Trial log-mean SD','s_t', s_trl);
add('Trial x person log-mean SD','s_pt', s_txp);
add('Occasion x person log-mean SD','s_po', s_oxp);
add('Trial x occasion log-mean SD','s_ot', s_txo);

%dimension slope blocks (log-mean b and log-nu b_nu)
slopelabels = local_slopelabels(ndim,dim_names);
for k = 1:numel(slopelabels)
    add(['Log-mean slope, ' slopelabels{k}], sprintf('b_%d', k), bmean(:,k));
end
for k = 1:numel(slopelabels)
    add(['Log-nu slope, ' slopelabels{k}], sprintf('b_nu_%d', k), bnu(:,k));
end

T = getT();
end

function T = local_gamma_trt_ls_varcomp_table(alpha0,logsig0,bmean,bsig,sig_p,...
    sig_sd,s_occ,s_trl,s_txp,s_oxp,s_txo,rho_ps,ndim,dim_names,ciperc)
%Posterior summaries of the reported LOCATION-SCALE gamma dynrel components
%(trial + occasion two-facet, gammascale = 2). The location-scale twin of
%local_gamma_trt_varcomp_table, and the two-facet sibling of
%local_gamma_ls_varcomp_table.
%
%Every row that names the scale submodel changes, and that is the whole point of
%having a separate builder: gro_sds column 2 is s_sd (the person log residual-SD
%SD) here, not s_v; the scale intercept is a log SD, not a log nu; and the
%dimension slopes are b_sigma, not b_nu. Reusing the log-nu builder would print a
%table of correct numbers under wrong names - it would not error - which is why
%the two are separate rather than one function with a flag.
%
%The five crossed log-mean SD rows are IDENTICAL to the log-nu builder's, and
%deliberately so: those facets live on the mean submodel, which this
%parameterization leaves untouched. They are what makes the six observed-score
%signal components parameterization-invariant.
%
%rho_ps (the person mean <-> log-sigma correlation) IS reported, because it is
%estimated and interpretable. It is deliberately absent from the observed-scale
%conversion below: psyrat_gamma_varcomps_trt_ls takes no covariance argument,
%because the residual never touches mu. Note there is consequently no cov_pv line
%here, where the log-nu builder must compute one.
[add,getT] = local_vc_accumulator(ciperc);

alpha0 = alpha0(:); logsig0 = logsig0(:);
sig_p = sig_p(:); sig_sd = sig_sd(:);
s_occ = s_occ(:); s_trl = s_trl(:); s_txp = s_txp(:); s_oxp = s_oxp(:);
s_txo = s_txo(:); rho_ps = rho_ps(:);

%observed-score components at z = 0, marginal over the person residual-SD
%population - the same converter used per z by psyrat_rel_dynrel_trt_gamma_ls.
%Argument order follows the log-nu builder: (alpha, s_p, s_o, s_t, s_po, s_pt,
%s_ot, ...) with s_po = occasion x person (s_oxp), s_pt = trial x person (s_txp).
vc0 = psyrat_gamma_varcomps_trt_ls(alpha0, sig_p, s_occ, s_trl, s_oxp, s_txp, ...
    s_txo, logsig0, sig_sd);
add('Person variance at z = 0 (observed score)','sigma^2_p(0)', vc0.sigma_p2);
add('Occasion variance at z = 0 (observed score)','sigma^2_o(0)', vc0.sigma_o2);
add('Trial variance at z = 0 (observed score)','sigma^2_t(0)', vc0.sigma_t2);
add('Trial x person variance at z = 0','sigma^2_pt(0)', vc0.sigma_pt2);
add('Occasion x person variance at z = 0','sigma^2_po(0)', vc0.sigma_po2);
add('Trial x occasion variance at z = 0','sigma^2_ot(0)', vc0.sigma_ot2);
add('Residual variance at z = 0 (observed score)','sigma^2_pto,e(0)', vc0.sigma_res2);

%log-scale parameters (what the model estimates)
add('Log-mean grand intercept (z = 0)','alpha(0)', alpha0);
add('Log residual-SD grand intercept (z = 0)','log sigma(0)', logsig0);
add('Population residual SD at z = 0','sigma(0)', exp(logsig0));
add('Person log-mean SD','s_p', sig_p);
add('Person log residual-SD SD','s_sd', sig_sd);
add('Person mean<->log-scale correlation','rho_p,sigma', rho_ps);
add('Occasion log-mean SD','s_o', s_occ);
add('Trial log-mean SD','s_t', s_trl);
add('Trial x person log-mean SD','s_pt', s_txp);
add('Occasion x person log-mean SD','s_po', s_oxp);
add('Trial x occasion log-mean SD','s_ot', s_txo);

%dimension slope blocks (log-mean b and log residual-SD b_sigma)
slopelabels = local_slopelabels(ndim,dim_names);
for k = 1:numel(slopelabels)
    add(['Log-mean slope, ' slopelabels{k}], sprintf('b_%d', k), bmean(:,k));
end
for k = 1:numel(slopelabels)
    add(['Log residual-SD slope, ' slopelabels{k}], sprintf('b_sigma_%d', k),...
        bsig(:,k));
end

T = getT();
end

function T = local_trt_varcomp_table(sig_p,s_occ,s_trl,s_txp,s_oxp,s_txo,...
    log0,bsig,bmean,sig_dp,rho,ndim,dim_names,ciperc)
%Posterior summaries of the reported trial + occasion two-facet components
%(Greek notation). The residual is the dimension-conditioned highest-order
%person x trial x occasion term, reported here at z = 0 (the OVERALL-sample mean
%of each dimension, pooled across any groups); the dimension dependence is shown
%by the b_sigma slope block and the
%reliability surface.
[add,getT] = local_vc_accumulator(ciperc);

add('Between-person variance (universe score, conditional on dimensions)',...
    'sigma^2_p', sig_p.^2);
add('Occasion main-effect variance','sigma^2_o', s_occ.^2);
add('Trial main-effect variance','sigma^2_i', s_trl.^2);
add('Person x trial variance','sigma^2_pt', s_txp.^2);
add('Person x occasion variance','sigma^2_po', s_oxp.^2);
add('Trial x occasion variance','sigma^2_to', s_txo.^2);
add('Residual variance at z = 0','sigma^2_pto,e(0)', exp(log0).^2);
add('Person scale SD','sigma_delta_p', sig_dp);
add('Location-scale correlation','rho', rho);

%dimension slope blocks (mean b and scale b_sigma)
slopelabels = local_slopelabels(ndim,dim_names);
for k = 1:numel(slopelabels)
    add(['Mean slope, ' slopelabels{k}], sprintf('beta_%d', k), bmean(:,k));
end
for k = 1:numel(slopelabels)
    add(['Log-residual slope, ' slopelabels{k}], sprintf('beta_sigma_%d', k),...
        bsig(:,k));
end

T = getT();
end

function T = local_difftrt_varcomp_table(idvc,trvc,ocvc,tidvc,oidvc,tovc,...
    bsig,bsdim,bdim,ndim,dim_names,ciperc)
%Posterior summaries of the reported trial + occasion two-facet difference-score
%components (Greek notation). Each crossed component is reported as its
%gain-minus-loss CONTRAST variance (X11 + X22 - 2*X12) from the event 2x2
%covariance block, plus the event-specific per-event residual at z = 0 and the
%event-specific mean/log-residual dimension-slope blocks. The contrast variances
%are conditional on the dimensions (the mean dimension slopes remove the
%dimension dependence from the level components).
[add,getT] = local_vc_accumulator(ciperc);

%gain-minus-loss contrast variance from a draws x 2 x 2 event covariance block
delta = @(M) M(:,1,1) + M(:,2,2) - 2*M(:,2,1);

add('Difference universe-score variance (conditional on dimensions)',...
    'sigma^2_p_delta', delta(idvc));
add('Person variance, event 1','sigma^2_p,1', idvc(:,1,1));
add('Person variance, event 2','sigma^2_p,2', idvc(:,2,2));
add('Person covariance (events)','sigma_p,12', idvc(:,2,1));
add('Trial difference variance','sigma^2_i_delta', delta(trvc));
add('Occasion difference variance','sigma^2_o_delta', delta(ocvc));
add('Person x trial difference variance','sigma^2_pi_delta', delta(tidvc));
add('Person x occasion difference variance','sigma^2_po_delta', delta(oidvc));
add('Trial x occasion difference variance','sigma^2_io_delta', delta(tovc));
add('Residual variance at z = 0, event 1','sigma^2_pto,e,1(0)', exp(bsig(:,1)).^2);
add('Residual variance at z = 0, event 2','sigma^2_pto,e,2(0)', exp(bsig(:,2)).^2);

%event-specific slope blocks (odd columns event 1, even columns event 2)
terms = local_slopelabels(ndim,dim_names);
for k = 1:size(bdim,2)
    ev = 2 - mod(k,2);            % odd column -> event 1, even -> event 2
    tm = terms{ceil(k/2)};
    add(sprintf('Mean slope, event %d, %s',ev,tm),...
        sprintf('beta_%d,%d',ev,ceil(k/2)), bdim(:,k));
end
for k = 1:size(bsdim,2)
    ev = 2 - mod(k,2);
    tm = terms{ceil(k/2)};
    add(sprintf('Log-residual slope, event %d, %s',ev,tm),...
        sprintf('beta_sigma_%d,%d',ev,ceil(k/2)), bsdim(:,k));
end

T = getT();
end

function [add,getT] = local_vc_accumulator(ciperc)
%Shared accumulator for the variance-component tables. Returns two function
%handles that share this invocation's workspace (the standard MATLAB
%nested-function idiom): add(name,symbol,draws) appends one row -- the posterior
%mean plus the central credible interval of draws -- and getT() returns the
%assembled five-column table. Factored out so every local_*_varcomp_table builder
%shares one definition of the per-row math and the
%{component,symbol,estimate,ci_lower,ci_upper} column layout instead of
%re-declaring an identical nested closure and table() call.
ciedge = (1-ciperc)/2;
name = {}; sym = {}; ptv = []; llv = []; ulv = [];
    function add_(n,s,draws)
        name{end+1,1} = n; sym{end+1,1} = s;
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

function labels = local_slopelabels(ndim,dim_names)
%Dimension-slope row labels shared by the varcomp-table builders: the single
%dimension name for one dimension, or {dim1, dim2, dim1 x dim2} for two (the
%main-effect and interaction slope order used by every dynamic-reliability model).
%The event-specific difference builders index this list with ceil(k/2).
if ndim == 1
    labels = {char(string(dim_names{1}))};
else
    labels = {char(string(dim_names{1})),char(string(dim_names{2})),...
        [char(string(dim_names{1})) ' x ' char(string(dim_names{2}))]};
end
end

function T = local_dod_varcomp_table(bb,bsig,bdim,bsigdim,sdid,corid,...
    ndim,dim_names,ciperc,facets)
%Posterior summaries of the nonconcurrent DoD components for the GAUSSIAN
%family (analyses 28/29, family = 'gaussian'), the sibling of
%local_gamma_dod_varcomp_table_ls below. Under the identity link the person
%and facet SDs ARE observed-scale quantities and are constant in z
%(SCIENTIFIC_FORMULA_AUDIT.md section 26), so - unlike the gamma table's
%log-scale disclaimer - these rows are the observed-scale components
%themselves. The residual row stays exp(bsig).^2: the scale submodel is
%log-linked in both families. Argument shapes as the gamma sibling; bb is a
%draws x 4 NATURAL-SCALE per-cell mean intercept.
[add,getT] = local_vc_accumulator(ciperc);

for q = 1:4
    add(sprintf('Mean at z = 0, cell %d',q),...
        sprintf('mu_%d(0)',q), bb(:,q));
end
for q = 1:4
    add(sprintf('Residual variance at z = 0, cell %d',q),...
        sprintf('sigma^2_e,%d(0)',q), exp(bsig(:,q)).^2);
end
for q = 1:4
    add(sprintf('Person location SD, cell %d',q),...
        sprintf('tau_p,%d',q), sdid(:,q));
end
for q = 1:4
    add(sprintf('Person log-residual-SD SD, cell %d',q),...
        sprintf('tau_sigma,%d',q), sdid(:,4+q));
end
%the six cross-cell location correlations drive the universe-score variance
%through the signed contrast; the four within-cell location-scale
%correlations and the six cross-cell scale correlations are estimated by the
%joint block and reported for completeness (none enters the composite; the
%residual channel is exactly diagonal under this family).
for a = 1:3
    for b = (a+1):4
        add(sprintf('Person location correlation, cells %d-%d',a,b),...
            sprintf('rho_p,%d%d',a,b), corid(:,a,b));
    end
end
for q = 1:4
    add(sprintf('Person location-scale correlation, cell %d',q),...
        sprintf('rho_ps,%d',q), corid(:,q,4+q));
end
for a = 1:3
    for b = (a+1):4
        add(sprintf('Person scale correlation, cells %d-%d',a,b),...
            sprintf('rho_s,%d%d',a,b), corid(:,4+a,4+b));
    end
end
for f = 1:size(facets,1)
    flabel = facets{f,1};
    fsd = facets{f,2};
    fcor = facets{f,3};
    for q = 1:4
        add(sprintf('%s SD, cell %d',flabel,q),...
            sprintf('tau_%s,%d',lower(strrep(flabel,' ','')),q), fsd(:,q));
    end
    for a = 1:3
        for b = (a+1):4
            add(sprintf('%s correlation, cells %d-%d',flabel,a,b),...
                sprintf('rho_%s,%d%d',lower(strrep(flabel,' ','')),a,b),...
                fcor(:,a,b));
        end
    end
end
terms = local_slopelabels(ndim,dim_names);
for q = 1:4
    for k = 1:size(bdim,3)
        add(sprintf('Mean slope, cell %d, %s',q,terms{k}),...
            sprintf('beta_%d,%d',q,k), bdim(:,q,k));
    end
end
for q = 1:4
    for k = 1:size(bsigdim,3)
        add(sprintf('Log-residual slope, cell %d, %s',q,terms{k}),...
            sprintf('beta_sigma_%d,%d',q,k), bsigdim(:,q,k));
    end
end

T = getT();
end

function T = local_gamma_dod_varcomp_table_ls(bb,bsig,bdim,bsigdim,sdid,corid,...
    ndim,dim_names,ciperc,facets)
%Posterior summaries of the nonconcurrent DoD components (analyses 28/29),
%LOG-scale quantities: the observed-scale components are z-dependent under the
%log link (every one carries exp(2*alpha(z))), so the z-conditioned surfaces
%and the per-participant table are where observed-scale values live. bb/bsig
%are draws x 4 ABSOLUTE per-cell intercepts; bdim/bsigdim draws x 4 x KDIM;
%sdid draws x 8 [loc_1..4, logsd_1..4]; corid draws x 8 x 8; facets is an
%N x 3 cell {label, sd draws x 4, cor draws x 4 x 4} - the trial block for
%analysis 28, plus occasion/interaction blocks for 29.
[add,getT] = local_vc_accumulator(ciperc);

for q = 1:4
    add(sprintf('Log mean at z = 0, cell %d',q),...
        sprintf('alpha_%d(0)',q), bb(:,q));
end
for q = 1:4
    add(sprintf('Residual variance at z = 0, cell %d',q),...
        sprintf('sigma^2_e,%d(0)',q), exp(bsig(:,q)).^2);
end
for q = 1:4
    add(sprintf('Person location SD (log scale), cell %d',q),...
        sprintf('tau_p,%d',q), sdid(:,q));
end
for q = 1:4
    add(sprintf('Person log-residual-SD SD, cell %d',q),...
        sprintf('tau_sigma,%d',q), sdid(:,4+q));
end
%the six cross-cell location correlations drive the universe-score variance
%through the signed contrast; the four within-cell location-scale
%correlations and the six cross-cell scale correlations are estimated by the
%joint block and reported for completeness (none enters the conversion; the
%residual never touches the mean surface under this parameterization).
for a = 1:3
    for b = (a+1):4
        add(sprintf('Person location correlation, cells %d-%d',a,b),...
            sprintf('rho_p,%d%d',a,b), corid(:,a,b));
    end
end
for q = 1:4
    add(sprintf('Person location-scale correlation, cell %d',q),...
        sprintf('rho_ps,%d',q), corid(:,q,4+q));
end
for a = 1:3
    for b = (a+1):4
        add(sprintf('Person scale correlation, cells %d-%d',a,b),...
            sprintf('rho_s,%d%d',a,b), corid(:,4+a,4+b));
    end
end
for f = 1:size(facets,1)
    flabel = facets{f,1};
    fsd = facets{f,2};
    fcor = facets{f,3};
    for q = 1:4
        add(sprintf('%s SD (log scale), cell %d',flabel,q),...
            sprintf('tau_%s,%d',lower(strrep(flabel,' ','')),q), fsd(:,q));
    end
    for a = 1:3
        for b = (a+1):4
            add(sprintf('%s correlation, cells %d-%d',flabel,a,b),...
                sprintf('rho_%s,%d%d',lower(strrep(flabel,' ','')),a,b),...
                fcor(:,a,b));
        end
    end
end
terms = local_slopelabels(ndim,dim_names);
for q = 1:4
    for k = 1:size(bdim,3)
        add(sprintf('Mean slope, cell %d, %s',q,terms{k}),...
            sprintf('beta_%d,%d',q,k), bdim(:,q,k));
    end
end
for q = 1:4
    for k = 1:size(bsigdim,3)
        add(sprintf('Log-residual slope, cell %d, %s',q,terms{k}),...
            sprintf('beta_sigma_%d,%d',q,k), bsigdim(:,q,k));
    end
end

T = getT();
end

function local_warn_dod_invalidity(flags,label)
%Signed-projection invalidity disclosure for the nonconcurrent DoD designs.
%The convention is record-never-repair: a draw whose signed composite is
%materially negative or whose coefficient leaves [0,1] is reported as-is with
%a flag (the reference's posterior_invalid_* columns), so the user must be
%told when any evaluation was flagged. Small nonzero fractions are expected
%on real data - the reference bundle's own completed observed-data runs carry
%0.24%% and 0.013%% - so this warns rather than errors.
if isempty(flags) || ~flags.recorded
    return
end
if ~flags.flags_complete
    %absence of evidence must stay distinguishable from evidence of absence:
    %a table missing one design's flag column cannot earn an all-clear
    warning('dynrel:dodinvalidity',...
        ['Stratum ''%s'': the signed-projection invalidity flags were not '...
        'fully recorded (one D-study design''s column is missing); the '...
        'counted values below understate. %d flagged of %d evaluations '...
        'under the recorded design(s).'],...
        char(string(label)),flags.n_invalid + flags.n_invalid_std,...
        flags.n_evals);
    return
end
if flags.n_invalid == 0 && flags.n_invalid_std == 0
    return
end
warning('dynrel:dodinvalidity',...
    ['Stratum ''%s'': %d of %d per-participant draw evaluations are flagged '...
    'invalid under the actual D-study design (%d under the standardized '...
    'one): a signed composite was materially negative or a coefficient '...
    'left [0,1]. Values are reported unrepaired; inspect '...
    'summary.strata(s).dod_flags and the per-participant invalid_frac '...
    'columns.'],...
    char(string(label)),flags.n_invalid,flags.n_evals,flags.n_invalid_std);
end
