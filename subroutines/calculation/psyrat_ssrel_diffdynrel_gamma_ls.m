function [gen,dep] = psyrat_ssrel_diffdynrel_gamma_ls(varargin)
%Per-participant CONDITIONAL dynamic/conditional DIFFERENCE-score reliability for
%the ONE-FACET LOCATION-SCALE Gamma / scaled-chi-square subject-level design
%(analysis 13, ic_diff_dynrel_sserrvar, family = 'gamma', gammascale = 2).
%
%[gen,dep] = psyrat_ssrel_diffdynrel_gamma_ls('b',b,'b_sigma',bs,...
%  'b_dim',bd,'b_sigma_dim',bsd,'sd_id',sd_id,'sd_trl',sd_trl,...
%  'cor_id',cor_id,'cor_i',cor_i,'er_ss1',e1,'er_ss2',e2,...
%  'idtable',ssinfo,'ndim',1,'CI',.95)
%
%WHAT THIS IS. The cross of two existing files: psyrat_rel_diffdynrel_gamma_ls
%(GROUP level, two events, evaluated over a z GRID) and
%psyrat_ssrel_dynrel_gamma_ls (SUBJECT level, one event, evaluated at each
%PARTICIPANT'S OWN z). It is the subject-level read-out of the case-12
%location-scale model, which it shares outright - the emitted Stan is byte
%identical (psyrat_build_stan_diffdynrel_sserr_chisq is a one-line delegation),
%so this file is the entire difference between analyses 12 and 13.
%
%WHY IT CANNOT REUSE local_subj_ssrel_table (the Gaussian path). That local binds
%psyrat_ssrel_diff ONCE with a single person and trial (co)variance block shared
%across participants, and varies only each participant's residual. That is valid
%under the Gaussian model, where sigma_p and sigma_i do not depend on z. Under
%gamma the MEAN-SURFACE components carry exp(2*alpha(z)) - sigma_p2, sigma_i2 and
%both cross-event covariances - and each participant sits at a DIFFERENT z, so
%the person and trial blocks must be rebuilt per participant.
%
%SAY "MEAN-SURFACE", NOT "EVERY". An earlier draft of this paragraph said every
%observed-scale component carries exp(2*alpha(z)). That is true of the LOG-NU
%sibling, from which the sentence was copied, and false here: the whole point of
%location-scale is that the conditional residual does NOT carry the amplitude
%(e_cond_var = exp(2*log_sigma + 2*s_sd^2) contains no alpha at all), and
%sigma_pi_e2 is exp(2*alpha)*K + e_cond_var, which is not proportional to
%exp(2*alpha) either. The z-dependence that forces the per-participant loop comes
%from the mean submodel; it is enough on its own.
%
%psyrat_ssrel_diff indexes bp/bt positionally as single [ndraws x 2 x 2] arrays
%and validates NOTHING about their shape, so a per-subject block cannot be passed
%in one call - it would be silently misread rather than rejected. The resolution
%keeps the vetted arithmetic: psyrat_ssrel_diff is called ONCE PER PARTICIPANT on
%a one-row idtable and the rows are stacked. Nothing in psyrat_ssrel_diff
%changes. This is the same resolution, for the same reason, that
%psyrat_ssrel_dynrel_gamma documents for the non-difference designs.
%
%THREE POINTS THAT ARE EASY TO GET WRONG.
%
% (a) THE SIGNAL IS NOT CONDITIONED ON THE PARTICIPANT. Both sigma_p2 and
%     sigma_i2 are evaluated at alpha_k(z_s), NOT at alpha_k(z_s) + u_mu_s, but
%     they are two different kinds of quantity and the reason differs.
%     sigma_p2 is a BETWEEN-person variance: a participant's own mean effect says
%     where they sit, not how spread the population is, so conditioning it on that
%     participant is a category error. sigma_i2 is the TRIAL main effect - a
%     population quantity that has nothing to condition on in the first place.
%     (The log-nu sibling states this correctly in the singular; an earlier draft
%     of this file called both of them between-person variances.) The coefficient
%     is deliberately a hybrid - population signal at the participant's own z over
%     that participant's own realized measurement error - and that hybrid is what
%     makes the numbers comparable across participants.
%
% (b) THE RESIDUAL CARRIES NO TRIAL-AVERAGING FACTOR. There is no exp(2*s_i^2)
%     here. That factor belongs to the log-nu path only, where the conditional
%     variance 2*mu^2/nu is MEAN-COUPLED so the per-person residual is an
%     expectation over the trial facet. Under location-scale
%     Var(Y | mu, sigma) = sigma^2 has no trial index, so there is no trial
%     average to take, and adding the dimension slope does not change that - it
%     shifts the log residual SD without making the conditional variance a
%     function of the trial. psyrat_ssrel_diff divides the residual by each
%     participant's TRIAL COUNT, so what is passed must be the SINGLE-TRIAL
%     conditional variance. sd_trl never reaches the two log-SD channels handed
%     to psyrat_ssrel_diff (er_var and er_var_ss1/2), and that absence is the
%     observable signature. It does legitimately reach the TRIAL block, via
%     psyrat_gamma_varcomps_ls and psyrat_gamma_crosscov - "unused on the residual
%     side" is the precise claim, not "unused".
%
% (c) THE PER-PARTICIPANT RESIDUAL IS POOLED (S17-FIX, owner-directed
%     2026-08-09). Each participant's error is sigma_pi_e2 evaluated at their
%     own log residual SD: the mean-surface person x trial term
%     exp(2*alpha(z))*K - built from the POPULATION alpha(z), s_p and s_i, the
%     same term the group surface carries - PLUS that participant's own
%     conditional variance exp(2*(log_sigma(z) + w_s)). That restores the
%     COMBINED sigma_pt,e^2 form Rast & Clayson (in press, p. 4: the residual
%     captures the person-by-task interaction together with remaining
%     unattributed variation) define at the population level; the
%     per-participant instantiation - population person x trial term plus the
%     participant's own conditional variance - is the owner's reference
%     implementation's convention, and the construction the concurrent sibling
%     (analysis 19, psyrat_ssrel_diffdynrel_rescor_gamma_ls) uses. Before
%     S17-FIX this file passed the participant's conditional variance ALONE,
%     which understated every participant's error - the pre-fix behavior and
%     its magnitude (3.9% to 67% of the retained residual across plausible
%     designs, 0.02% in the owner's reference ERP data) are finding S17; see
%     gamma_scale_submodel_decision.md section J.6.
%
%     TWO CAVEATS THAT SURVIVE THE POOLING.
%       - The pooled term is the POPULATION-AVERAGE person x trial term at the
%         participant's z, not "the participant's own person x trial term".
%         That is an estimand choice, not an identifiability constraint - the
%         per-person share V_s is computable from the fitted u_s (J.6 derives
%         it) - and it is the choice the owner's reference implementation and
%         the group surface both make, which keeps u_s out of the coefficient.
%       - Against the plotted GROUP CURVE the direction is STILL not fixed. The
%         curve uses the residual MARGINALIZED over the person residual-SD
%         population, K + exp(2*log_sigma(z) + 2*s_sd^2), and the MEDIAN
%         per-event n'; a participant uses their own w_s and their OWN counts.
%         Either difference alone flips the comparison, so per-subject points
%         legitimately fall on both sides of the curve.
%
%     The cross-event half of the mean-surface term (the per-person instance of
%     S14) stays at ZERO by separate owner decision - the pooling is the
%     DIAGONAL only, per event. The two-facet designs 25/27 are structurally
%     different (person x trial is separately estimable there) and take no
%     analog of this pooling.
%
%INPUTS (name-value; every SD is on the LOG expected-score scale)
%   b            - [ndraws x 2] per-event log-mean intercepts at z = 0
%   b_sigma      - [ndraws x 2] per-event log residual-SD intercepts at z = 0
%   b_dim        - [ndraws x KDIM] event-crossed log-mean dimension slopes
%   b_sigma_dim  - [ndraws x KDIM] event-crossed log residual-SD dimension slopes
%                  (odd columns event 1, even columns event 2, matching the
%                  case-12 Xdim construction)
%   sd_id        - [ndraws x 4] person SDs [mean_1, mean_2, logsigma_1, logsigma_2]
%   sd_trl       - [ndraws x 2] per-event trial log-SDs
%   cor_id       - [ndraws x 4 x 4] person correlation draws (only (1,2) is read)
%   cor_i        - [ndraws x 1] cross-event trial correlation
%   er_ss1/er_ss2- [ndraws x NSUB] per-subject log residual SD AT z = 0, per event
%                  (b_sigma[e] + r_id[:,e+2]); the dimension contribution is added
%                  here, per participant, matching the Gaussian storage convention
%   idtable      - ssinfo table: id, id2, z1[, z2], trls1, trls2, trls
%   ndim         - 1 or 2
%   CI           - credible-interval mass (e.g. .95)
%
%OUTPUTS
%   gen, dep - stacked one-row-per-participant tables, as returned by
%              psyrat_ssrel_diff for est = 'gen' and est = 'dep'
%
%See also PSYRAT_REL_DIFFDYNREL_GAMMA_LS, PSYRAT_SSREL_DYNREL_GAMMA_LS,
%PSYRAT_GAMMA_VARCOMPS_LS, PSYRAT_GAMMA_CROSSCOV, PSYRAT_SSREL_DIFF.

% Copyright (C) 2016-2026 Peter E. Clayson
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

b           = psyrat_req(varargin,'b','Please input b (per-event log-mean intercepts).');
b_sigma     = psyrat_req(varargin,'b_sigma','Please input b_sigma (per-event log residual-SD intercepts).');
b_dim       = psyrat_req(varargin,'b_dim','Please input b_dim (log-mean dimension slopes).');
b_sigma_dim = psyrat_req(varargin,'b_sigma_dim','Please input b_sigma_dim (log residual-SD dimension slopes).');
sd_id       = psyrat_req(varargin,'sd_id','Please input sd_id (person log-scale SDs, draws x 4).');
sd_trl      = psyrat_req(varargin,'sd_trl','Please input sd_trl (per-event trial log-SDs).');
cor_id      = psyrat_req(varargin,'cor_id','Please input cor_id (person 4x4 correlation draws).');
cor_i       = psyrat_req(varargin,'cor_i','Please input cor_i (cross-event trial correlation).');
er_ss1      = psyrat_req(varargin,'er_ss1','Please input er_ss1 (per-subject log residual SD at z=0, event 1).');
er_ss2      = psyrat_req(varargin,'er_ss2','Please input er_ss2 (per-subject log residual SD at z=0, event 2).');
idtable     = psyrat_req(varargin,'idtable','Please input idtable (ssinfo).');
ndim        = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
ciperc      = psyrat_req(varargin,'CI','Please input CI.');

if ~any(ndim == [1 2])
    error('varargin:ndim','WARNING: ndim must be 1 or 2.\n');
end
if size(sd_id,2) ~= 4
    error('varargin:sd_id',...
        ['WARNING: sd_id must have 4 columns, ordered '...
        '[mean_1, mean_2, logsigma_1, logsigma_2]; got %d.\n'],size(sd_id,2));
end
if size(b,2) ~= 2 || size(b_sigma,2) ~= 2 || size(sd_trl,2) ~= 2
    error('varargin:events',...
        ['WARNING: b, b_sigma and sd_trl must each have 2 columns, one per '...
        'difference-score event.\n']);
end
%KDIM is 2 per dimension term and must split evenly into the two events.
kdim = size(b_dim,2);
if size(b_sigma_dim,2) ~= kdim || mod(kdim,2) ~= 0 || kdim ~= 2*(2*ndim-1)
    error('varargin:b_dim',...
        ['WARNING: b_dim and b_sigma_dim must both have %d columns for '...
        'ndim = %d (2 per dimension term, event-crossed); got %d and %d.\n'],...
        2*(2*ndim-1),ndim,kdim,size(b_sigma_dim,2));
end

nsub = height(idtable);
if size(er_ss1,2) < max(idtable.id2) || size(er_ss2,2) < max(idtable.id2)
    %FAIL LOUDLY. psyrat_ssrel_diff's own per-row range check falls back to the
    %POPULATION residual for any id2 it cannot index, with no warning - so a short
    %er_ss matrix would silently turn a subject-level table into a population one.
    error('varargin:er_ss',...
        ['WARNING: er_ss1/er_ss2 have %d and %d subject columns but idtable '...
        'references id2 up to %d.\n'],...
        size(er_ss1,2),size(er_ss2,2),max(idtable.id2));
end

%event-crossed slopes: odd columns -> event 1, even columns -> event 2, matching
%the case-12 Xdim construction [meas1.*z1, meas2.*z1, meas1.*z2, meas2.*z2, ...].
bdim1 = b_dim(:,1:2:end);      %draws x (KDIM/2)
bdim2 = b_dim(:,2:2:end);
bsd1  = b_sigma_dim(:,1:2:end);
bsd2  = b_sigma_dim(:,2:2:end);

%per-event person mean SDs and trial SDs (log scale, z-independent). The person
%log-sigma SDs sd_id(:,3:4) are deliberately NOT read: this is the CONDITIONAL
%estimand, so each participant's own realized scale effect is used instead of
%marginalizing over their population - which is precisely what makes the
%coefficient subject-level rather than a population coefficient evaluated at z.
sp_a = sd_id(:,1); sp_b = sd_id(:,2);
si_a = sd_trl(:,1); si_b = sd_trl(:,2);

%The CROSS-event person-mean correlation drives the difference numerator. As in
%the group-level twin there is no within-event mean<->scale covariance to form:
%psyrat_gamma_varcomps_ls has no cov_pv argument because the residual never
%touches mu. Every cor_id entry involving slots 3-4 is deliberately unread.
cor_p = cor_id(:,1,2);
cor_i = cor_i(:);

ndraws = size(b,1);
id_varcov  = zeros(ndraws,2,2);
trl_varcov = zeros(ndraws,2,2);
%non-concurrent: the two events come from different physical trials, so the
%conditional draws are independent given mu and sigma. One column because each
%call carries exactly one participant.
zc = zeros(ndraws,1);

genRows = cell(nsub,1);
depRows = cell(nsub,1);

%One warning per table, not one per participant: psyrat_diffrel warns on every
%call and this loop calls it 2*nsub times with deliberately unequal per-event
%trial counts. Suppress here and let the caller's surface pass report it.
adm_ws = warning('off','psyrat:negerrvar');
cleanupWs = onCleanup(@() warning(adm_ws));

for c = 1:nsub

    %dimension design row for THIS participant's z (KDIM/2 x 1), matching the
    %Xdim layout: [z1] for one dimension, [z1; z2; z1*z2] for two.
    if ndim == 1
        zt = idtable.z1(c);
    else
        zt = [idtable.z1(c); idtable.z2(c); idtable.z1(c)*idtable.z2(c)];
    end
    zt = zt(:);

    %z-conditioned per-event log-mean intercepts and POPULATION log residual SDs
    %at this participant's z. Each event is evaluated at its OWN z-specific
    %predictors (ruling 13); asymmetric slopes on either submodel legitimately
    %move the coefficient and must not be suppressed (ruling 15).
    alpha_a = b(:,1) + bdim1 * zt;
    alpha_b = b(:,2) + bdim2 * zt;
    lsig_a  = b_sigma(:,1) + bsd1 * zt;
    lsig_b  = b_sigma(:,2) + bsd2 * zt;

    %THIS PARTICIPANT'S OWN log residual SD per event. er_ss holds the z = 0 value
    %(b_sigma[e] + r_id[:,e+2]); the dimension contribution is added here, which is
    %the same split the Gaussian path documents in local_subj_ssrel_table.
    col = idtable.id2(c);
    ss_a = er_ss1(:,col) + bsd1 * zt;
    ss_b = er_ss2(:,col) + bsd2 * zt;

    %OBSERVED-score-scale components at this participant's z. Passing the
    %PARTICIPANT'S OWN log residual SD as log_sigma (with s_sd = 0) makes
    %sigma_pi_e2 the POOLED per-participant residual this design reports
    %(S17-FIX): the mean-surface person x trial term, built from the population
    %alpha at z, plus that participant's own conditional variance - the same
    %construction the concurrent sibling uses
    %(psyrat_ssrel_diffdynrel_rescor_gamma_ls, analysis 19). sigma_p2 and
    %sigma_i2 do not depend on the log_sigma argument, so the signal blocks are
    %unchanged by the pooling. NOTE the asymmetry between the two calls: the
    %converter is the LOCATION-SCALE one, while the cross-event helper is the
    %SAME one the log-nu path uses, because its whole argument list is
    %mean-submodel.
    vc_a = psyrat_gamma_varcomps_ls(alpha_a, sp_a, si_a, ss_a, 0);
    vc_b = psyrat_gamma_varcomps_ls(alpha_b, sp_b, si_b, ss_b, 0);
    cc   = psyrat_gamma_crosscov(alpha_a, sp_a, si_a, ...
                                 alpha_b, sp_b, si_b, cor_p, cor_i);

    %assemble the 2x2 person and trial (co)variance arrays psyrat_ssrel_diff
    %expects. It reads (:,1,1), (:,2,2) and the LOWER triangle (:,2,1); the upper
    %triangle is filled too so the arrays are symmetric for any future reader.
    id_varcov(:,1,1) = vc_a.sigma_p2;  id_varcov(:,2,2) = vc_b.sigma_p2;
    id_varcov(:,1,2) = cc.cov_p_obs;   id_varcov(:,2,1) = cc.cov_p_obs;
    trl_varcov(:,1,1) = vc_a.sigma_i2; trl_varcov(:,2,2) = vc_b.sigma_i2;
    trl_varcov(:,1,2) = cc.cov_i_obs;  trl_varcov(:,2,1) = cc.cov_i_obs;

    %psyrat_ssrel_diff reads BOTH residual channels as log-scale SDs and squares
    %them (exp(.)^2), so the pooled variances are handed over as 0.5*log(.).
    %er_var is the pooled POPULATION (typical-subject, w = 0) reference at this
    %z; er_ss_a/er_ss_b are this participant's own pooled residuals and are what
    %the subject-level path actually uses.
    er_ss_a = 0.5*log(vc_a.sigma_pi_e2);
    er_ss_b = 0.5*log(vc_b.sigma_pi_e2);
    vcp_a = psyrat_gamma_varcomps_ls(alpha_a, sp_a, si_a, lsig_a, 0);
    vcp_b = psyrat_gamma_varcomps_ls(alpha_b, sp_b, si_b, lsig_b, 0);
    er_var_z = [0.5*log(vcp_a.sigma_pi_e2), 0.5*log(vcp_b.sigma_pi_e2)];

    %ONE-ROW call. psyrat_ssrel_diff indexes the er_var_ss columns by idtable.id2,
    %and here exactly one column is supplied, so id2 must be remapped to 1 for the
    %call. Leaving the participant's global id2 in place would send it out of
    %range, and that failure is SILENT - the row would fall back to the population
    %residual and look like a plausible subject-level number. The original id2 is
    %restored on the returned row so the output table keeps its real index.
    row = idtable(c,:);
    row.id2 = 1;

    g = psyrat_ssrel_diff('bp',id_varcov,'bt',trl_varcov,'er_var',er_var_z,...
        'idtable',row,'CI',ciperc,'est','gen',...
        'er_var_ss1',er_ss_a,'er_var_ss2',er_ss_b,'wp_cov_ss',zc,'wp_cov',0);
    d = psyrat_ssrel_diff('bp',id_varcov,'bt',trl_varcov,'er_var',er_var_z,...
        'idtable',row,'CI',ciperc,'est','dep',...
        'er_var_ss1',er_ss_a,'er_var_ss2',er_ss_b,'wp_cov_ss',zc,'wp_cov',0);

    g.id2(:) = idtable.id2(c);
    d.id2(:) = idtable.id2(c);
    genRows{c} = g;
    depRows{c} = d;
end

gen = vertcat(genRows{:});
dep = vertcat(depRows{:});

end
