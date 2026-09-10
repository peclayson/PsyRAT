function [gen,dep] = psyrat_ssrel_dynrel_trt_gamma(varargin)
%Per-participant CONDITIONAL dynamic/conditional reliability for the TWO-FACET
%(trial + occasion) Gamma / scaled-chi-square subject-level design (analysis 27,
%ic_dynrel_sserrvar_trt, family = 'gamma').
%
%[gen,dep] = psyrat_ssrel_dynrel_trt_gamma('alpha0',a0,'b',b,'lognu0',l0,...
%   'b_nu',bnu,'sig_p',sp,'sig_occ',so,'sig_trl',st,'sig_occxid',spo,...
%   'sig_trlxid',spt,'sig_trlxocc',sot,'ind_bs',u,'ind_nu',v,...
%   'idtable',ssinfo,'ndim',1,'reltype',3,'nocc',n,'CI',.95)
%
%WHAT THIS IS. The two-facet sibling of psyrat_ssrel_dynrel_gamma, and the
%subject-level sibling of psyrat_rel_dynrel_trt_gamma. It combines the two
%pieces the other designs establish separately:
%
%  - from analysis 25 (psyrat_gamma_extract_trt_sserr): the multi-facet
%    conditional residual, whose within-person mean term sums the FIVE crossed
%    log-mean variances and excludes the person main effect;
%  - from analysis 26 (psyrat_ssrel_dynrel_gamma): evaluation at each
%    participant's own z, with every observed-scale component recomputed there.
%
%For participant s at their own z_s:
%
%  alpha(z_s) = alpha0 + b . x(z_s)        lognu(z_s) = lognu0 + b_nu . x(z_s)
%  sigma_p2/o2/t2/pt2/po2/ot2 at z_s   <- psyrat_gamma_varcomps_trt
%  V_w = s_o^2 + s_t^2 + s_po^2 + s_pt^2 + s_ot^2
%  sigma_pot,e^2(s,z_s) = 2*exp( 2*(alpha(z_s) + u_s) + 2*V_w - (lognu(z_s) + v_s) )
%
%WHY psyrat_ssrel_trt IS CALLED PER PARTICIPANT. psyrat_ssrel_trt already
%supports per-subject txp/oxp/err, but it binds bp, bo, bt and txo ONCE, outside
%its subject loop. Under the Gaussian two-facet dynrel that is correct, because
%only the residual moves with z. Under gamma every observed-scale component
%carries exp(2*alpha(z)), so all six crossed components move with z as well - and
%each participant sits at a different z. Rather than change psyrat_ssrel_trt (a
%calculator RC-06 recently reworked), this function calls it once per participant
%on a one-row idtable with that participant's own components, and stacks the
%rows. The coefficient, ICC and SEM formulas - including the ICC pinned to
%obs = nocc = 1 for composite invariance - stay in exactly one place.
%
%psyrat_ssrel_trt's gro_icc/gro_rel reference columns are computed from whatever
%it is handed, so on a one-row call they are per-participant rather than group
%quantities. That is harmless here: psyrat_dynrel_summary's
%local_pack_ssrel_nondiff copies only the dep_*/icc_*/sem_* columns into the
%shared ssrel_table schema and drops those two.
%
%WHAT IS CONDITIONED ON WHAT (the qualification (a)-(d) below assumes). Every
%component is conditioned on the participant's own z. Only the RESIDUAL is
%additionally conditioned on their own (u_s, v_s); the six crossed mean components
%remain population quantities at that z, i.e. still integrated over the person
%distribution. That is deliberate and follows the case-6 precedent
%(psyrat_gamma_extract_sserr), where the group components are likewise population
%values and only the residual is per-person - so the per-participant error term
%mixes a person-conditional residual with person-marginal facet components. The
%owner ruling this file implements ("own-z, all components") is about the z axis,
%not about (u_s, v_s).
%
%ESTIMAND NOTES. Identical to psyrat_ssrel_dynrel_gamma, which states them in
%full:
% (a) the between-person numerator is evaluated at alpha(z_s), NOT at
%     alpha(z_s) + u_s - conditioning on a participant says where they sit, not
%     how spread the population is;
% (b) the residual drops the person facet variance sig_p^2 because u_s is
%     conditioned on rather than integrated out, and is EXACT under the ID-only
%     nu structure;
% (c) (s_v, cov_pv) are passed to the converter as ZERO, which is what turns
%     estimand #2 off. Only the six mean components are consumed and none of them
%     depends on those arguments, so it is numerically inconsequential - but it
%     states the conditional intent;
% (d) the log-mean SLOPE b still cancels exactly, as on the group surface (every
%     component and the residual carry exp(2*alpha(z))); what does NOT cancel is
%     the participant's OWN log-mean effect u_s, which enters the residual but not
%     the between-person numerator, so a higher-mean participant receives a
%     strictly lower conditional coefficient. Two different quantities - do not
%     restate this as "the slope cancellation breaks", which is false.
%
%Inputs (all draw vectors are columns of posterior draws)
% alpha0      - log-mean grand intercept at z = 0            [ndraws x 1]
% b           - log-mean dimension slopes                 [ndraws x KDIM]
% lognu0      - log-nu grand intercept at z = 0               [ndraws x 1]
% b_nu        - log-nu dimension slopes                   [ndraws x KDIM]
% sig_p       - person log-mean SD (s_p)                      [ndraws x 1]
% sig_occ     - occasion log-mean SD (s_o)                    [ndraws x 1]
% sig_trl     - trial log-mean SD (s_t)                       [ndraws x 1]
% sig_occxid  - person x occasion log-mean SD (s_po)          [ndraws x 1]
% sig_trlxid  - person x trial log-mean SD (s_pt)             [ndraws x 1]
% sig_trlxocc - trial x occasion log-mean SD (s_ot)           [ndraws x 1]
% ind_bs      - per-participant log-mean effects u_p       [ndraws x NSUB]
% ind_nu      - per-participant log-nu effects v_p         [ndraws x NSUB]
% idtable     - ssinfo: one row per participant with id, id2, z1[, z2], trls,
%               occs. id2 indexes the ind_bs/ind_nu COLUMNS.
% ndim        - 1 or 2 dimensions
% reltype     - 1 = equivalence, 2 = stability, 3 = equivalence-and-stability
% nocc        - occasion count n'_o used for the coefficient
% CI          - credible interval width in decimal format (.95 = 95%)
%
%Outputs
% gen - psyrat_ssrel_trt output table for generalizability (gcoeff = 2), one row
%       per participant in idtable row order
% dep - the same for dependability (gcoeff = 1)

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
lognu0  = psyrat_req(varargin,'lognu0','Please input lognu0 (Intercept_nu draws).');
b_nu    = psyrat_req(varargin,'b_nu','Please input b_nu (log-nu dimension slopes).');
sig_p   = psyrat_req(varargin,'sig_p','Please input sig_p (person log-mean SD draws).');
s_o     = psyrat_req(varargin,'sig_occ','Please input sig_occ (occasion log-mean SD draws).');
s_t     = psyrat_req(varargin,'sig_trl','Please input sig_trl (trial log-mean SD draws).');
s_po    = psyrat_req(varargin,'sig_occxid','Please input sig_occxid (person x occasion SD draws).');
s_pt    = psyrat_req(varargin,'sig_trlxid','Please input sig_trlxid (person x trial SD draws).');
s_ot    = psyrat_req(varargin,'sig_trlxocc','Please input sig_trlxocc (trial x occasion SD draws).');
ind_bs  = psyrat_req(varargin,'ind_bs','Please input ind_bs (per-subject log-mean effects).');
ind_nu  = psyrat_req(varargin,'ind_nu','Please input ind_nu (per-subject log-nu effects).');
idtable = psyrat_req(varargin,'idtable','Please input idtable (ssinfo).');
ndim    = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
reltype = psyrat_req(varargin,'reltype','Please input reltype (1, 2, or 3).');
nocc    = psyrat_req(varargin,'nocc','Please input nocc (occasion count).');
ciperc  = psyrat_req(varargin,'CI','Please input CI.');

if ~any(ndim == [1 2])
    error('varargin:ndim','WARNING: ndim must be 1 or 2.\n');
end

%column vectors of draws
alpha0 = alpha0(:);
lognu0 = lognu0(:);
sig_p  = sig_p(:);
s_o    = s_o(:);
s_t    = s_t(:);
s_po   = s_po(:);
s_pt   = s_pt(:);
s_ot   = s_ot(:);

%b/b_nu must be draws x KDIM (KDIM = 1 for one dimension; 3 for two), matching the
%Xdim layout the estimator builds. Same shape check psyrat_rel_dynrel_trt_gamma
%applies, so a transposed input fails loudly rather than broadcasting.
expectedk = (ndim == 1) * 1 + (ndim == 2) * 3;
if size(b,2) ~= expectedk
    error('varargin:b',...
        'WARNING: b has %d columns but ndim=%d expects %d.\n',...
        size(b,2),ndim,expectedk);
end
if size(b_nu,2) ~= expectedk
    error('varargin:bnu',...
        'WARNING: b_nu has %d columns but ndim=%d expects %d.\n',...
        size(b_nu,2),ndim,expectedk);
end

nsub = height(idtable);
if size(ind_bs,2) ~= size(ind_nu,2)
    error('varargin:indblocks',...
        'WARNING: ind_bs has %d columns but ind_nu has %d.\n',...
        size(ind_bs,2),size(ind_nu,2));
end
if any(idtable.id2 > size(ind_bs,2)) || any(idtable.id2 < 1)
    error('varargin:id2',...
        ['WARNING: idtable.id2 indexes outside the per-subject draw columns '...
        '(max id2 = %d, columns = %d).\n'],max(idtable.id2),size(ind_bs,2));
end

%WITHIN-person log-mean variance: every log-mean facet EXCEPT the person main
%effect, which is dropped because u_p is conditioned on rather than integrated
%out. This is the 2*s_i^2 of the one-facet case generalized to five facets, and
%is the same V_w psyrat_gamma_extract_trt_sserr forms for the static analysis 25.
V_w = s_o.^2 + s_t.^2 + s_po.^2 + s_pt.^2 + s_ot.^2;   % [ndraws x 1]

genRows = cell(nsub,1);
depRows = cell(nsub,1);

for c = 1:nsub
    col = idtable.id2(c);

    %dimension design row for THIS participant's z (KDIM x 1), matching Xdim:
    %[z1] for one dimension, [z1; z2; z1*z2] for two.
    if ndim == 1
        xrow = idtable.z1(c);
    else
        xrow = [idtable.z1(c); idtable.z2(c); idtable.z1(c)*idtable.z2(c)];
    end

    alpha_z = alpha0 + b * xrow(:);          % [ndraws x 1]
    lognu_z = lognu0 + b_nu * xrow(:);       % [ndraws x 1]
    nu_z    = exp(lognu_z);

    %observed-scale crossed components at this participant's z. (s_v, cov_pv) = 0:
    %the conditional estimand; the six mean components do not depend on them.
    vc = psyrat_gamma_varcomps_trt(alpha_z, sig_p, s_o, s_t, s_po, s_pt, s_ot, ...
        nu_z, 0, 0);

    %per-participant CONDITIONAL residual at their own z and their own (u,v)
    sig_pot_e2 = 2 .* exp(2.*(alpha_z + ind_bs(:,col)) + 2.*V_w ...
        - (lognu_z + ind_nu(:,col)));        % [ndraws x 1]
    err_ss = sqrt(sig_pot_e2);

    row = idtable(c,:);

    %Slot mapping follows psyrat_gamma_extract_trt: oxp <- sigma_po2 (person x
    %occasion), txp <- sigma_pt2 (person x trial). The varcomps names and the Stan
    %names are transposed relative to each other, so this pairing is deliberate.
    %
    %NOTE the *_ss arguments are deliberately NOT passed. psyrat_ssrel_trt's
    %per-subject path indexes them by idtable.id2, which on a one-row call is that
    %participant's ORIGINAL index (3, 7, ...) rather than 1, so the bounds check
    %would fail and it would silently fall back to the non-_ss values anyway.
    %Since every argument below already carries THIS participant's own components,
    %the fallback is the correct path - passing *_ss as well would add an
    %ignored-for-most-participants input that reads as if it did something.
    base = {'reltype',reltype,...
        'bp',sqrt(vc.sigma_p2),'bo',sqrt(vc.sigma_o2),'bt',sqrt(vc.sigma_t2),...
        'txp',sqrt(vc.sigma_pt2),'oxp',sqrt(vc.sigma_po2),...
        'txo',sqrt(vc.sigma_ot2),'err',err_ss,...
        'idtable',row,'CI',ciperc,'nocc',nocc};

    %gcoeff = 2 -> generalizability G_s; gcoeff = 1 -> dependability phi_s
    genRows{c} = psyrat_ssrel_trt('gcoeff',2,base{:});
    depRows{c} = psyrat_ssrel_trt('gcoeff',1,base{:});
end

gen = vertcat(genRows{:});
dep = vertcat(depRows{:});

end
