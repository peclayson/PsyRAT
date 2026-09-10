function [gen,dep] = psyrat_ssrel_dynrel_gamma_ls(varargin)
%Per-participant CONDITIONAL dynamic/conditional reliability for the ONE-FACET
%LOCATION-SCALE Gamma / scaled-chi-square subject-level design (analysis 26,
%ic_dynrel_sserrvar, family = 'gamma', gammascale = 2).
%
%[gen,dep] = psyrat_ssrel_dynrel_gamma_ls('alpha0',a0,'b',b,'logsig0',ls0,...
%  'b_sigma',bs,'sig_p',sp,'sig_i',si,'ind_sd',w,'idtable',ssinfo,'ndim',1,...
%  'CI',.95)
%
%WHAT THIS IS. The subject-level sibling of psyrat_rel_dynrel_gamma_ls, and the
%location-scale twin of psyrat_ssrel_dynrel_gamma. The surface function
%evaluates a POPULATION coefficient over a z grid, marginalizing over the person
%residual-SD population; this function evaluates one coefficient PER
%PARTICIPANT, at that participant's own z, CONDITIONING on that participant's
%own log residual SD.
%
%WHY IT IS SIMPLER THAN ITS LOG-NU TWIN, AND WHY THAT IS THE POINT.
%Under gammascale = 1 the conditional variance is 2*mu^2/nu, so the twin must
%DERIVE each participant's residual from two person effects (ind_bs = u_p and
%ind_nu = v_p) and carry the trial-averaging factor exp(2*s_i^2). Here the model
%estimates the quantity directly:
%
%   log sigma_s(z) = logsig0 + b_sigma . x(z) + w_s        (w_s = ind_sd)
%   Var(Y | s)     = sigma_s(z)^2 = exp(2*(logsig0 + b_sigma . x(z) + w_s))
%
%with no mu anywhere, and no trial-averaging factor because the conditional
%variance does not depend on the trial. That conditional variance is read
%straight off ind_sd, the same way the Gaussian subject-level path reads
%ind_sdlog.
%
%THE REPORTED PER-PARTICIPANT RESIDUAL IS POOLED (S17-FIX, owner-directed
%2026-08-09). What each participant is charged is not their conditional
%variance alone but sigma_pi_e2 evaluated at their own log residual SD: the
%mean-surface person x trial term exp(2*alpha(z))*K - built from the POPULATION
%alpha(z), s_p and s_i, the same term the group surface carries - PLUS the
%conditional variance above. That restores the COMBINED sigma_pt,e^2 form Rast
%& Clayson (in press, p. 4) define at the population level; the
%per-participant instantiation - population person x trial term plus the
%participant's own conditional variance - follows the owner's reference
%implementation, and analyses 13 and 19 use the same construction. The
%pre-fix conditional-only behavior and its magnitude are finding S17; see
%gamma_scale_submodel_decision.md section J.6.
%
%CONSEQUENCE: ind_bs IS NOT AN INPUT, AND MUST NOT BECOME ONE.
%The participant's OWN log-mean u_s does not appear in this calculation at all:
%
% (a) The residual is exp(2*alpha(z_s))*K + exp(2*(log sigma_s)), and neither
%     term touches u_s: the pooled person x trial term is the POPULATION
%     AVERAGE at the participant's z, not their own share. That is an estimand
%     choice, not an identifiability constraint - the per-person share V_s is
%     computable from the fitted u_s (gamma_scale_submodel_decision.md J.6
%     derives it) - and it is the choice the owner's reference implementation
%     and the group surface both make, which keeps u_s out of the coefficient.
%     Under the log-nu twin the residual is 2*mu(z,u_s)^2/nu, so u_s DOES
%     enter there and a higher-mean participant receives a strictly lower
%     conditional coefficient (that twin's header records G 0.974 -> 0.932 for
%     u_s = 0 -> 0.5).
% (b) sigma_p^2 is a BETWEEN-person variance, so it is evaluated at alpha(z_s)
%     and NOT at alpha(z_s) + u_s - conditioning on a participant tells you
%     where they sit, not how spread the population is. This is the same rule
%     the log-nu twin states, and it is unchanged here.
%
%So under gammascale = 2 the per-participant coefficient is INDEPENDENT of the
%participant's own mean level. That decoupling is the reason this
%parameterization exists, and it means the per-participant spread behaves
%qualitatively unlike its log-nu twin. It is a property, not an omission: a
%future reader who adds ind_bs "for parity" would reintroduce exactly the
%mean-dispersion blending the owner ruled against. The estimator still STORES
%ind_bs on REL.out for parity and as a diagnostic; nothing here reads it.
%
%WHAT DOES STILL MOVE WITH z. Each participant sits at their own z_s, and both
%submodels are evaluated there: sigma_p^2(z_s) and sigma_i^2(z_s) carry
%exp(2*alpha(z_s)) while the residual carries exp(2*(logsig0 + b_sigma.x(z_s))).
%Unlike the log-nu twin, those two do NOT share a common factor, so the log-mean
%slope b does not cancel out of the coefficient here either (see
%psyrat_rel_dynrel_gamma_ls for the algebra). Both slopes matter.
%
%WHY NOT REUSE psyrat_ssrel IN ONE CALL, as the Gaussian analysis 26 does. That
%path can, because under the Gaussian location-scale dynrel sigma_p and sigma_i
%are z-invariant, so one shared bp and i serve every subject and only wp_ss
%varies. Under gamma every observed-scale component carries exp(2*alpha(z)) and
%each participant sits at a different z, so bp and i differ per participant.
%psyrat_ssrel takes bp and i as single columns of draws with a hard size guard,
%so it cannot carry a per-subject bp/i in one call. The resolution is the twin's:
%call psyrat_ssrel ONCE PER PARTICIPANT on a one-row idtable and stack the rows.
%Nothing in psyrat_ssrel changes (its RC-19 argument contract is untouched), and
%the per-subject coefficient, ICC and SEM formulas stay in exactly one place.
%
%THE RESIDUAL IS HANDED OVER AS A POPULATION/PERSON SPLIT, not folded into one
%term. psyrat_ssrel forms wp = exp(wp_pop + wp_ss)^2 and treats it as the
%RELATIVE residual sigma_pi,e^2(s), adding i^2 only for the absolute branch -
%which, since S17-FIX, is exactly what the channel carries. wp_pop is the
%pooled TYPICAL-SUBJECT (w = 0) residual log-SD at this participant's z,
%0.5*log(exp(2*alpha(z))*K + exp(2*log sigma(z))), so psyrat_ssrel's pop_errvar
%column keeps meaning what its name says; wp_ss is the remainder, so the two
%recompose to the participant's own pooled residual exactly. The log-nu twin
%folds everything into wp_pop with wp_ss = 0 instead.
%
%(s_sd) is passed to psyrat_gamma_varcomps_ls as ZERO in both calls. sigma_p2
%and sigma_i2 do not depend on it, and the residual read off vcs must be
%CONDITIONAL on the participant's own realized w_s rather than marginalized
%over the person residual-SD population - vcs's log_sigma argument already
%carries w_s, and adding s_sd on top would double-count the person scale
%heterogeneity. That conditioning is what makes this calculation subject-level
%rather than a population coefficient evaluated at z.
%
%Inputs (all draw vectors are columns of posterior draws)
% alpha0  - log-mean grand intercept at z = 0 (Intercept)          [ndraws x 1]
% b       - log-mean dimension slopes                           [ndraws x KDIM]
% logsig0 - log residual-SD grand intercept at z = 0 (Intercept_sigma)
%                                                                  [ndraws x 1]
% b_sigma - log residual-SD dimension slopes                    [ndraws x KDIM]
% sig_p   - person log-mean SD (s_p)                               [ndraws x 1]
% sig_i   - trial main-effect log-mean SD (s_i)                    [ndraws x 1]
% ind_sd  - per-participant log residual-SD effects w_p         [ndraws x NSUB]
% idtable - ssinfo: one row per participant with id, id2, z1[, z2], trls.
%           id2 indexes the ind_sd COLUMNS.
% ndim    - 1 or 2 dimensions
% CI      - credible interval width in decimal format (.95 = 95%)
%
%Outputs
% gen - psyrat_ssrel output table for generalizability (gcoeff = 2), one row per
%       participant in idtable row order
% dep - the same for dependability (gcoeff = 1)
%
%Both carry psyrat_ssrel's dep_*/icc_*/sem_* column names regardless of gcoeff;
%psyrat_dynrel_summary's local_pack_ssrel_nondiff renames them into the shared
%ssrel_table schema.
%
%See also PSYRAT_SSREL_DYNREL_GAMMA, PSYRAT_REL_DYNREL_GAMMA_LS,
%PSYRAT_GAMMA_VARCOMPS_LS, PSYRAT_SSREL.

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

alpha0  = psyrat_req(varargin,'alpha0','Please input alpha0 (Intercept draws).');
b       = psyrat_req(varargin,'b','Please input b (log-mean dimension slopes).');
logsig0 = psyrat_req(varargin,'logsig0','Please input logsig0 (Intercept_sigma draws).');
b_sigma = psyrat_req(varargin,'b_sigma','Please input b_sigma (log residual-SD dimension slopes).');
sig_p   = psyrat_req(varargin,'sig_p','Please input sig_p (person log-mean SD draws).');
sig_i   = psyrat_req(varargin,'sig_i','Please input sig_i (trial log-mean SD draws).');
ind_sd  = psyrat_req(varargin,'ind_sd','Please input ind_sd (per-subject log residual-SD effects).');
idtable = psyrat_req(varargin,'idtable','Please input idtable (ssinfo).');
ndim    = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
ciperc  = psyrat_req(varargin,'CI','Please input CI.');

if ~any(ndim == [1 2])
    error('varargin:ndim','WARNING: ndim must be 1 or 2.\n');
end

%column vectors of draws
alpha0  = alpha0(:);
logsig0 = logsig0(:);
sig_p   = sig_p(:);
sig_i   = sig_i(:);

%b/b_sigma must be draws x KDIM (KDIM = 1 for one dimension; 3 for two),
%matching the Xdim layout the estimator builds. Same shape check
%psyrat_rel_dynrel_gamma_ls applies, so a transposed input fails loudly rather
%than broadcasting into a wrong answer.
expectedk = (ndim == 1) * 1 + (ndim == 2) * 3;
if size(b,2) ~= expectedk
    error('varargin:b',...
        'WARNING: b has %d columns but ndim=%d expects %d.\n',...
        size(b,2),ndim,expectedk);
end
if size(b_sigma,2) ~= expectedk
    error('varargin:bsigma',...
        'WARNING: b_sigma has %d columns but ndim=%d expects %d.\n',...
        size(b_sigma,2),ndim,expectedk);
end

if any(idtable.id2 > size(ind_sd,2)) || any(idtable.id2 < 1)
    error('varargin:id2',...
        ['WARNING: idtable.id2 indexes outside the per-subject draw columns '...
        '(max id2 = %d, columns = %d).\n'],max(idtable.id2),size(ind_sd,2));
end

nsub = height(idtable);
genRows = cell(nsub,1);
depRows = cell(nsub,1);

for c = 1:nsub
    col = idtable.id2(c);

    %dimension design row for THIS participant's z (KDIM x 1), matching the
    %Xdim layout: [z1] for one dimension, [z1; z2; z1*z2] for two.
    if ndim == 1
        xrow = idtable.z1(c);
    else
        xrow = [idtable.z1(c); idtable.z2(c); idtable.z1(c)*idtable.z2(c)];
    end

    %z-conditioned log-mean intercept and population log residual SD at this
    %participant's z
    alpha_z  = alpha0  + b       * xrow(:);   % [ndraws x 1]
    logsig_z = logsig0 + b_sigma * xrow(:);   % [ndraws x 1]

    %observed-scale GROUP components at this participant's z (vc, evaluated at
    %the population log residual SD with s_sd = 0), and the POOLED
    %per-participant residual (vcs, evaluated at the PARTICIPANT'S OWN log
    %residual SD, S17-FIX): sigma_pi_e2 there is the mean-surface person x
    %trial term exp(2*alpha(z))*K plus this participant's own conditional
    %variance - the same construction analyses 13 and 19 use. sigma_p2 and
    %sigma_i2 do not depend on the log_sigma argument, so the signal terms are
    %read off vc unchanged.
    vc  = psyrat_gamma_varcomps_ls(alpha_z, sig_p, sig_i, logsig_z, 0);
    vcs = psyrat_gamma_varcomps_ls(alpha_z, sig_p, sig_i, ...
        logsig_z + ind_sd(:,col), 0);

    %per-participant POOLED residual, as the population/person split
    %psyrat_ssrel expects: wp = exp(wp_pop + wp_ss)^2 = sigma_pi_e2(s). wp_pop
    %carries the pooled typical-subject (w = 0) reference, so psyrat_ssrel's
    %pop_errvar column keeps its meaning; wp_ss carries the remainder. No
    %trial-averaging factor - see the header.
    wp_pop = 0.5*log(vc.sigma_pi_e2);         % [ndraws x 1]
    wp_ss  = 0.5*log(vcs.sigma_pi_e2) - wp_pop;   % [ndraws x 1]

    row = idtable(c,:);

    %gcoeff = 2 -> generalizability G_s (relative error, excludes sigma_i^2)
    %gcoeff = 1 -> dependability  phi_s (absolute error, adds sigma_i^2)
    genRows{c} = psyrat_ssrel('bp',{sqrt(vc.sigma_p2)},'wp_pop',wp_pop,...
        'wp_ss',{wp_ss},'idtable',row,'CI',ciperc,...
        'i',sqrt(vc.sigma_i2),'gcoeff',2);
    depRows{c} = psyrat_ssrel('bp',{sqrt(vc.sigma_p2)},'wp_pop',wp_pop,...
        'wp_ss',{wp_ss},'idtable',row,'CI',ciperc,...
        'i',sqrt(vc.sigma_i2),'gcoeff',1);
end

gen = vertcat(genRows{:});
dep = vertcat(depRows{:});

end
