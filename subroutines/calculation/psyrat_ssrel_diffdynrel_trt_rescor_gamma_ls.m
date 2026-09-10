function [gen,dep,diagnostic] = psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls(varargin)
%Per-participant CONDITIONAL dynamic/conditional CONCURRENT DIFFERENCE-score
%reliability for the TWO-FACET (trial x occasion) gamma location-scale MODULAR
%CUT-COPULA design (analysis 20, ic_diff_dynrel_sserrvar_trt_rescor,
%family = 'gamma', gammascale = 2). Stage 3 of the modular workflow.
%
%[gen,dep,diagnostic] = psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls('b',b,...)
%
%WHAT THIS IS. The two-facet twin of PSYRAT_SSREL_DIFFDYNREL_RESCOR_GAMMA_LS.
%Read that file first: the per-participant loop, the z conditioning, the one-row
%coefficient call, the id2 remap and the draw-alignment contract are taken from
%it unchanged, as are two of its three scientific claims - the residuals are
%COUPLED through a quadrature-derived covariance, and the participant's own
%location effect enters through the ERROR while the signal stays a population
%quantity.
%
%ITS THIRD CLAIM DOES NOT CARRY ACROSS, and neither does its marginalization.
%Both differences are below, and both are consequences of the design rather than
%of plumbing.
%
%THE RESIDUAL SLOT IS NOT THE ONE-FACET RESIDUAL SLOT. In the one-facet design
%person x trial is confounded with measurement error - there is one observation
%per person x trial cell - so the toolbox reports them together and the pooled
%channel is sigma_pi_e2. Here person x trial is SEPARATELY ESTIMABLE, is reported
%as its own component, and PSYRAT_DIFFREL_TRT divides it by n_i where it divides
%the residual by n_i*n_o. Folding it in would change the D-study arithmetic, not
%just a label. The quantity pooled here is the highest-order term instead: the
%three-way mean-surface remainder plus the participant's own conditional
%variance, which is PSYRAT_GAMMA_VARCOMPS_TRT_LS's sigma_res2 evaluated at their
%own log residual SD with s_sd = 0.
%
%The reference implementation splits it the same way: R/reliability.R:241-246
%selects residual_name = "pio_e" and mean_name = "pio_mean" for the two-facet
%branch where the one-facet branch selects "pi_e"/"pi_mean", and :266-268 form
%pio_e = pio_mean + sigma^2 with the matching covariance line. DO NOT
%"reconcile" this with the one-facet calculator by pooling person x trial; the
%two designs disagree on purpose and both follow the reference.
%
%THE MARGINALIZATION SUMS ALL FIVE NON-PERSON FACETS. The expected observed score
%that sets the quadrature point is exp(mu_person + 0.5*v_other) with v_other the
%sum of the trial, occasion, trial x person, occasion x person and trial x
%occasion log-scale variances - not the trial variance alone, which is the
%one-facet rule (R/reliability.R:250-256). Everything else about the estimand is
%as in the one-facet sibling: the SIGNAL stays a population variance evaluated at
%the participant's z, while the participant's own realized location effect enters
%through the ERROR, by setting where the copula covariance is evaluated.
%
%THE ESTIMAND OF THE RESIDUAL COVARIANCE. Evaluated at expected means, it is a
%MEAN-MATCHED DIAGNOSTIC quantity, NOT integrated over every realized
%location-facet effect, and must not be described as a fully facet-marginalized
%residual covariance. See documentation/modular_cut_copula_workflow.md.
%
%DEVELOPMENT-ONLY. Unlike the one-facet design this configuration has no
%empirical anchor: the reference bundle's two-facet runs were prepared against a
%SIMULATED second occasion and were never sampled, so no completed external fit
%exists to replay against. The arithmetic here is pinned against the reference's
%own R implementation through frozen fixtures, not against results.
%
%DRAW ALIGNMENT IS THE CALLER'S JOB, exactly as in the one-facet sibling: the
%copula correlation exists only for the retained (thinned) stage-1 draws, so
%every array passed here must already be subsampled to that same set. A
%full-length marginal array against a short rho is rejected rather than recycled.
%
%INPUTS (name-value; every SD is on the LOG expected-score scale)
%   b, b_sigma, b_dim, b_sigma_dim, sd_id, cor_id, er_ss1, er_ss2,
%   mu_ss1, mu_ss2, rho, idtable, ndim, CI
%                 - as PSYRAT_SSREL_DIFFDYNREL_RESCOR_GAMMA_LS; see that file
%   sd_trl, sd_occ, sd_tid, sd_oid, sd_to
%                 - per-event log-SDs of the five non-person crossed facets
%                   (draws x 2 each): trial, occasion, trial x person,
%                   occasion x person, trial x occasion
%   cor_t, cor_o, cor_pt, cor_po, cor_ot
%                 - their cross-event correlations (draws x 1 each)
%   reltype       - 1 (CE), 2 (CS) or 3 (CES); see PSYRAT_DIFFREL_TRT
%   nocc          - occasion n' for the coefficient (scalar, shared across
%                   participants, matching the reference surface). Each
%                   participant's OBSERVED occasion counts ride along in the
%                   idtable (occs1/occs2) for reference only.
%
%OUTPUTS
%   gen, dep   - stacked one-row-per-participant tables, as psyrat_ssrel_diff_trt
%                returns for est = 'gen' and est = 'dep'
%   diagnostic - per-participant quadrature provenance, identical in shape to the
%                one-facet calculator's
%
%See also PSYRAT_SSREL_DIFFDYNREL_RESCOR_GAMMA_LS, PSYRAT_GAMMA_COPULA_RESCOV,
%PSYRAT_GAMMA_CROSSCOV_TRT_RESCOR, PSYRAT_GAMMA_VARCOMPS_TRT_LS,
%PSYRAT_SSREL_DIFF_TRT.

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
%

b           = psyrat_req(varargin,'b','Please input b (per-event log-mean intercepts).');
b_sigma     = psyrat_req(varargin,'b_sigma','Please input b_sigma (per-event log residual-SD intercepts).');
b_dim       = psyrat_req(varargin,'b_dim','Please input b_dim (log-mean dimension slopes).');
b_sigma_dim = psyrat_req(varargin,'b_sigma_dim','Please input b_sigma_dim (log residual-SD dimension slopes).');
sd_id       = psyrat_req(varargin,'sd_id','Please input sd_id (person log-scale SDs, draws x 4).');
sd_trl      = psyrat_req(varargin,'sd_trl','Please input sd_trl (per-event trial log-SDs).');
sd_occ      = psyrat_req(varargin,'sd_occ','Please input sd_occ (per-event occasion log-SDs).');
sd_tid      = psyrat_req(varargin,'sd_tid','Please input sd_tid (per-event trial x person log-SDs).');
sd_oid      = psyrat_req(varargin,'sd_oid','Please input sd_oid (per-event occasion x person log-SDs).');
sd_to       = psyrat_req(varargin,'sd_to','Please input sd_to (per-event trial x occasion log-SDs).');
cor_id      = psyrat_req(varargin,'cor_id','Please input cor_id (person 4x4 correlation draws).');
cor_t       = psyrat_req(varargin,'cor_t','Please input cor_t (cross-event trial correlation).');
cor_o       = psyrat_req(varargin,'cor_o','Please input cor_o (cross-event occasion correlation).');
cor_pt      = psyrat_req(varargin,'cor_pt','Please input cor_pt (cross-event trial x person correlation).');
cor_po      = psyrat_req(varargin,'cor_po','Please input cor_po (cross-event occasion x person correlation).');
cor_ot      = psyrat_req(varargin,'cor_ot','Please input cor_ot (cross-event trial x occasion correlation).');
er_ss1      = psyrat_req(varargin,'er_ss1','Please input er_ss1 (per-subject log residual SD at z=0, event 1).');
er_ss2      = psyrat_req(varargin,'er_ss2','Please input er_ss2 (per-subject log residual SD at z=0, event 2).');
mu_ss1      = psyrat_req(varargin,'mu_ss1','Please input mu_ss1 (per-subject log mean at z=0, event 1).');
mu_ss2      = psyrat_req(varargin,'mu_ss2','Please input mu_ss2 (per-subject log mean at z=0, event 2).');
rho         = psyrat_req(varargin,'rho','Please input rho (copula correlation per retained draw).');
idtable     = psyrat_req(varargin,'idtable','Please input idtable (ssinfo).');
ndim        = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
nocc        = psyrat_req(varargin,'nocc','Please input nocc (occasion n'').');
reltype     = psyrat_opt(varargin,'reltype',3);
ciperc      = psyrat_req(varargin,'CI','Please input CI.');

if ~any(ndim == [1 2])
    error('varargin:ndim','WARNING: ndim must be 1 or 2.\n');
end
if ~any(reltype == [1 2 3])
    error('varargin:reltype','WARNING: reltype must be 1, 2, or 3.\n');
end
if ~isscalar(nocc) || ~isfinite(nocc) || nocc < 1
    error('varargin:nocc','WARNING: nocc must be a finite positive scalar.\n');
end
if size(sd_id,2) ~= 4
    error('varargin:sd_id',...
        ['WARNING: sd_id must have 4 columns, ordered '...
        '[mean_1, mean_2, logsigma_1, logsigma_2]; got %d.\n'],size(sd_id,2));
end
if size(b,2) ~= 2 || size(b_sigma,2) ~= 2
    error('varargin:events',...
        ['WARNING: b and b_sigma must each have 2 columns, one per '...
        'difference-score event.\n']);
end
%Every crossed facet is per-event; a single column would broadcast silently and
%give both events the same facet SD, which is a different model.
facetnames = {'sd_trl','sd_occ','sd_tid','sd_oid','sd_to'};
facetvals  = {sd_trl,sd_occ,sd_tid,sd_oid,sd_to};
for f = 1:numel(facetnames)
    if size(facetvals{f},2) ~= 2
        error('varargin:facetsd',...
            ['WARNING: %s must have 2 columns, one per difference-score event; '...
            'got %d.\n'],facetnames{f},size(facetvals{f},2));
    end
end
kdim = size(b_dim,2);
if size(b_sigma_dim,2) ~= kdim || mod(kdim,2) ~= 0 || kdim ~= 2*(2*ndim-1)
    error('varargin:b_dim',...
        ['WARNING: b_dim and b_sigma_dim must both have %d columns for '...
        'ndim = %d (2 per dimension term, event-crossed); got %d and %d.\n'],...
        2*(2*ndim-1),ndim,kdim,size(b_sigma_dim,2));
end

ndraws = size(b,1);
rho = rho(:);
%DRAW ALIGNMENT. Every array must already be on the retained-draw set; see the
%header. Recycling a short rho would pair unrelated draws.
if numel(rho) ~= ndraws
    error('varargin:rho',...
        ['WARNING: rho has %d entries but the marginal draws have %d rows. '...
        'The copula correlation exists only for the RETAINED stage-1 draws, so '...
        'the caller must subsample every marginal array to that same set '...
        'before calling. How to fix: index the draw arrays by REL.out.copula '...
        '.sel.\n'],numel(rho),ndraws);
end
if any(~isfinite(rho)) || any(abs(rho) >= 1)
    error('varargin:rho',...
        'WARNING: every copula correlation must be finite with magnitude below one.\n');
end

nsub = height(idtable);
if size(er_ss1,2) < max(idtable.id2) || size(er_ss2,2) < max(idtable.id2)
    %FAIL LOUDLY: psyrat_ssrel_diff_trt's own range check silently falls back to
    %the POPULATION residual for an id2 it cannot index, turning a subject-level
    %table into a population one.
    error('varargin:er_ss',...
        ['WARNING: er_ss1/er_ss2 have %d and %d subject columns but idtable '...
        'references id2 up to %d.\n'],...
        size(er_ss1,2),size(er_ss2,2),max(idtable.id2));
end
if size(mu_ss1,2) < max(idtable.id2) || size(mu_ss2,2) < max(idtable.id2)
    error('varargin:mu_ss',...
        ['WARNING: mu_ss1/mu_ss2 have %d and %d subject columns but idtable '...
        'references id2 up to %d. These carry each participant''s own location '...
        'effect and set the point at which the copula covariance is '...
        'evaluated, so a short matrix would silently move the estimand.\n'],...
        size(mu_ss1,2),size(mu_ss2,2),max(idtable.id2));
end

%event-crossed slopes: odd columns -> event 1, even -> event 2.
bdim1 = b_dim(:,1:2:end);
bdim2 = b_dim(:,2:2:end);
bsd1  = b_sigma_dim(:,1:2:end);
bsd2  = b_sigma_dim(:,2:2:end);

%per-event log-scale SDs. NOTE the toolbox's tid is "trial x person" and so
%supplies psyrat_gamma_varcomps_trt_ls's s_pt, while oid is "occasion x person"
%and supplies s_po - the names CROSS OVER.
%
%That cross-over has to be held consistent at SIX places in this file, most of
%them positional and so unnamed at the call site: this unpacking; the two
%per-participant psyrat_gamma_varcomps_trt_ls calls (slots 5/6); the two
%POPULATION-fallback calls to the same function; the
%psyrat_gamma_crosscov_trt_rescor call (slots 5/6 per event AND the cor_po /
%cor_pt pair); the varcov assembly (sigma_pt2 -> tid, sigma_po2 -> oid); and the
%blockargs list ('bpi' <- tid, 'bpo' <- oid). A swap at any one of them is
%SILENT, because both indices are valid and the result is still a number. Change
%one, check all six. The surface sibling carries the same warning for its five.
sp_a  = sd_id(:,1);  sp_b  = sd_id(:,2);
st_a  = sd_trl(:,1); st_b  = sd_trl(:,2);
so_a  = sd_occ(:,1); so_b  = sd_occ(:,2);
spt_a = sd_tid(:,1); spt_b = sd_tid(:,2);
spo_a = sd_oid(:,1); spo_b = sd_oid(:,2);
sot_a = sd_to(:,1);  sot_b = sd_to(:,2);

cor_p  = cor_id(:,1,2);
cor_t  = cor_t(:);
cor_o  = cor_o(:);
cor_pt = cor_pt(:);
cor_po = cor_po(:);
cor_ot = cor_ot(:);

%ALL FIVE non-person facets enter the marginalization that turns a participant's
%log-mean into an expected observed score - not the trial facet alone. None of
%them carries a dimension slope, so this is formed once.
v_other_a = st_a.^2 + so_a.^2 + spt_a.^2 + spo_a.^2 + sot_a.^2;
v_other_b = st_b.^2 + so_b.^2 + spt_b.^2 + spo_b.^2 + sot_b.^2;

id_varcov  = zeros(ndraws,2,2);
trl_varcov = zeros(ndraws,2,2);
occ_varcov = zeros(ndraws,2,2);
tid_varcov = zeros(ndraws,2,2);
oid_varcov = zeros(ndraws,2,2);
to_varcov  = zeros(ndraws,2,2);

genRows = cell(nsub,1);
depRows = cell(nsub,1);
diagRows = cell(nsub,1);

%One warning per table rather than one per participant; see the one-facet sibling.
adm_ws = warning('off','psyrat:negerrvar');
cleanupWs = onCleanup(@() warning(adm_ws));

for c = 1:nsub

    if ndim == 1
        zt = idtable.z1(c);
    else
        zt = [idtable.z1(c); idtable.z2(c); idtable.z1(c)*idtable.z2(c)];
    end
    zt = zt(:);
    col = idtable.id2(c);

    %POPULATION log-mean and log residual SD at this participant's z. The signal
    %components are built from these, not from the participant's own effects.
    alpha_a = b(:,1) + bdim1 * zt;
    alpha_b = b(:,2) + bdim2 * zt;
    lsig_a  = b_sigma(:,1) + bsd1 * zt;
    lsig_b  = b_sigma(:,2) + bsd2 * zt;

    %THIS PARTICIPANT'S OWN log residual SD and log mean at their z.
    ss_a = er_ss1(:,col) + bsd1 * zt;
    ss_b = er_ss2(:,col) + bsd2 * zt;
    mu_a = mu_ss1(:,col) + bdim1 * zt;
    mu_b = mu_ss2(:,col) + bdim2 * zt;

    %Expected observed score for this participant, marginalizing over all five
    %non-person facets but NOT over their own location effect - that effect is
    %realized, not a population to average across.
    E_a = exp(mu_a + 0.5.*v_other_a);
    E_b = exp(mu_b + 0.5.*v_other_b);
    sig_a = exp(ss_a);
    sig_b = exp(ss_b);
    nu_a = 2 .* (E_a ./ sig_a).^2;
    nu_b = 2 .* (E_b ./ sig_b).^2;

    %QUADRATURE, per draw, NON-STRICT - both for the reasons the one-facet
    %calculator documents at length: the converter's node escalation is
    %per-argument, and a single tail draw can put the two events in regimes where
    %no node pair meets tolerance. The best estimate is kept, each such draw is
    %COUNTED, and the count is surfaced per participant. Zero is never
    %substituted, which would bias toward independence.
    cov_e = zeros(ndraws,1);
    nodes = zeros(ndraws,1);
    qerr = zeros(ndraws,1);
    unconverged = 0;
    for d = 1:ndraws
        q = psyrat_gamma_copula_rescov(nu_a(d), nu_b(d), rho(d), [], [], false);
        cov_e(d) = E_a(d) .* E_b(d) .* q.value;
        nodes(d) = q.nodes;
        qerr(d) = q.error;
        if ~q.converged
            unconverged = unconverged + 1;
        end
    end

    %OBSERVED-scale components at this participant's z. Passing the PARTICIPANT'S
    %log residual SD as log_sigma (with s_sd = 0) makes sigma_res2 the pooled
    %channel: the three-way mean-surface remainder, built from the population
    %alpha at z, plus that participant's own conditional variance. Person x trial
    %is NOT in it - see the header.
    vc_a = psyrat_gamma_varcomps_trt_ls(alpha_a, sp_a, so_a, st_a, ...
        spo_a, spt_a, sot_a, ss_a, 0);
    vc_b = psyrat_gamma_varcomps_trt_ls(alpha_b, sp_b, so_b, st_b, ...
        spo_b, spt_b, sot_b, ss_b, 0);
    cc = psyrat_gamma_crosscov_trt_rescor( ...
        alpha_a, sp_a, so_a, st_a, spo_a, spt_a, sot_a, ...
        alpha_b, sp_b, so_b, st_b, spo_b, spt_b, sot_b, ...
        cor_p, cor_o, cor_t, cor_po, cor_pt, cor_ot, cov_e);

    id_varcov(:,1,1)  = vc_a.sigma_p2;   id_varcov(:,2,2)  = vc_b.sigma_p2;
    id_varcov(:,1,2)  = cc.cov_p;        id_varcov(:,2,1)  = cc.cov_p;
    trl_varcov(:,1,1) = vc_a.sigma_t2;   trl_varcov(:,2,2) = vc_b.sigma_t2;
    trl_varcov(:,1,2) = cc.cov_t;        trl_varcov(:,2,1) = cc.cov_t;
    occ_varcov(:,1,1) = vc_a.sigma_o2;   occ_varcov(:,2,2) = vc_b.sigma_o2;
    occ_varcov(:,1,2) = cc.cov_o;        occ_varcov(:,2,1) = cc.cov_o;
    tid_varcov(:,1,1) = vc_a.sigma_pt2;  tid_varcov(:,2,2) = vc_b.sigma_pt2;
    tid_varcov(:,1,2) = cc.cov_pt;       tid_varcov(:,2,1) = cc.cov_pt;
    oid_varcov(:,1,1) = vc_a.sigma_po2;  oid_varcov(:,2,2) = vc_b.sigma_po2;
    oid_varcov(:,1,2) = cc.cov_po;       oid_varcov(:,2,1) = cc.cov_po;
    to_varcov(:,1,1)  = vc_a.sigma_ot2;  to_varcov(:,2,2)  = vc_b.sigma_ot2;
    to_varcov(:,1,2)  = cc.cov_ot;       to_varcov(:,2,1)  = cc.cov_ot;

    %psyrat_ssrel_diff_trt reads both residual channels as LOG SDs and squares
    %them, so the pooled variance is handed over as 0.5*log(.). The covariance
    %channel is read directly on the variance scale, hence no transform there.
    er_ss_a = 0.5*log(vc_a.sigma_res2);
    er_ss_b = 0.5*log(vc_b.sigma_res2);
    vcp_a = psyrat_gamma_varcomps_trt_ls(alpha_a, sp_a, so_a, st_a, ...
        spo_a, spt_a, sot_a, lsig_a, 0);
    vcp_b = psyrat_gamma_varcomps_trt_ls(alpha_b, sp_b, so_b, st_b, ...
        spo_b, spt_b, sot_b, lsig_b, 0);
    er_var_z = [0.5*log(vcp_a.sigma_res2), 0.5*log(vcp_b.sigma_res2)];

    row = idtable(c,:);
    row.id2 = 1;

    blockargs = {'bp',id_varcov,'bpi',tid_varcov,'bpo',oid_varcov,...
        'bt',trl_varcov,'bo',occ_varcov,'boi',to_varcov};
    g = psyrat_ssrel_diff_trt(blockargs{:},'er_var',er_var_z,...
        'idtable',row,'reltype',reltype,'nocc',nocc,'CI',ciperc,'est','gen',...
        'er_var_ss1',er_ss_a,'er_var_ss2',er_ss_b,...
        'er_cov_ss',cc.cov_res_e_obs,'er_cov',0);
    d_ = psyrat_ssrel_diff_trt(blockargs{:},'er_var',er_var_z,...
        'idtable',row,'reltype',reltype,'nocc',nocc,'CI',ciperc,'est','dep',...
        'er_var_ss1',er_ss_a,'er_var_ss2',er_ss_b,...
        'er_cov_ss',cc.cov_res_e_obs,'er_cov',0);

    g.id2(:) = idtable.id2(c);
    d_.id2(:) = idtable.id2(c);
    genRows{c} = g;
    depRows{c} = d_;

    %Provenance for the quadrature, per participant. res_cor is the residual
    %CORRELATION on the observed scale, which is NOT the copula correlation and
    %is reported so the two are not confused.
    res_cor = cov_e ./ (sig_a .* sig_b);
    diagRows{c} = table(idtable.id(c), idtable.id2(c), ...
        max(nodes), max(qerr), unconverged, unconverged/ndraws, ...
        mean(cov_e), mean(res_cor), mean(rho), ...
        'VariableNames',{'id','id2','quad_nodes_max','quad_err_max',...
        'quad_unconverged_n','quad_unconverged_frac',...
        'res_cov_pt','res_cor_pt','rho_e_pt'});
end

gen = vertcat(genRows{:});
dep = vertcat(depRows{:});
diagnostic = vertcat(diagRows{:});

end
