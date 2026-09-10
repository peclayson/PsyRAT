function [ssrel,drawlevel] = psyrat_ssrel_dodiffdynrel_trt(varargin)
%Per-participant person-specific dynamic NONCONCURRENT difference-of-
%differences reliability for the TWO-FACET (trials x occasions) GAUSSIAN
%design (analysis 29, ic_dodiff_dynrel_sserrvar_trt, family = 'gaussian').
%The Gaussian sibling of PSYRAT_SSREL_DODIFFDYNREL_TRT_GAMMA_LS. Derivation:
%SCIENTIFIC_FORMULA_AUDIT.md section 26.
%
%[ssrel,drawlevel] = psyrat_ssrel_dodiffdynrel_trt('b',b,...)
%
%WHAT THIS IS. The two-facet twin of PSYRAT_SSREL_DODIFFDYNREL - read that
%file first for the identity-link simplifications (components are the model
%blocks directly; expected scores are the additive predictor; no nu
%read-out). What the second facet adds is exactly what the gamma trt sibling
%documents:
%
% (1) SEVEN COMPONENT MATRICES (p, i, o, pi, po, io, pio_e), with the
%     residual pooled at the HIGHEST order only: under the identity link
%     pio_e is EXACTLY diag(sigma_qp^2) - the three-way mean-surface
%     remainder is identically zero. Person x trial and person x occasion
%     stay separate components with their own divisors.
%
% (2) THREE DECOMPOSITIONS, selected by reltype: reltype 1 = CE (the scaled
%     person x occasion term joins the universe), reltype 2 = CS (the scaled
%     person x trial term joins), reltype 3 = CES (pure person variance).
%     Because CE/CS universes carry a SCALED term, the table carries uni_pt
%     AND uni_pt_std.
%
% (3) THE INVALIDITY FLAG IS A PROPERTY OF THE DRAW, NOT OF THE SELECTED
%     DECOMPOSITION: a draw is flagged when ANY of the six scaled terms is
%     materially negative or ANY of the three decompositions'
%     errors/coefficients is out of range.
%
%INPUTS - as the one-facet Gaussian twin, plus:
% sd_occ/cor_occ - [ndraws x 4] / [ndraws x 4 x 4] occasion block
% sd_tid/cor_tid - person x trial block
% sd_oid/cor_oid - person x occasion block
% sd_to /cor_to  - trial x occasion block
% idtable        - additionally occs1..occs4 (per-cell occasion counts)
% nprime         - [1 x 4] standardized per-cell TRIAL counts n'_i
% nprime_o       - [1 x 4] standardized per-cell OCCASION counts n'_o
% reltype        - 1 (CE), 2 (CS), or 3 (CES)
%
%OUTPUTS - as the one-facet Gaussian twin (no nu columns); the coefficient
%columns are the selected reltype's. drawlevel additionally carries G/D for
%ALL THREE decompositions (G_ces/G_ce/G_cs and D_*). ssrel carries TWO
%floor-disclosure columns -- n_floored for the trial counts and
%n_floored_occ for the occasion counts -- because this design floors both
%facets and they are not interchangeable.
%
%See also PSYRAT_SSREL_DODIFFDYNREL, PSYRAT_DOD_COMPOSITE,
%PSYRAT_DOD_COMPONENTS, PSYRAT_SSREL_DODIFFDYNREL_TRT_GAMMA_LS.

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

for coupled = {'rho','rho_e','cov_e_obs','cov_res_e_obs','er_cov','rescor','copula'}
    if any(strcmpi(coupled{1},varargin(1:2:end)))
        error('psyrat_ssrel_dodiffdynrel_trt:noresidualcoupling',... %Error code and associated error
            ['''%s'' is not an input: all six observation-level residual '...
            'covariances are fixed to zero by the nonconcurrent estimand.'],...
            coupled{1});
    end
end

b           = psyrat_req(varargin,'b','Please input b (per-cell natural-scale mean intercepts).');
b_sigma     = psyrat_req(varargin,'b_sigma','Please input b_sigma (per-cell absolute log residual-SD intercepts).');
b_dim       = psyrat_req(varargin,'b_dim','Please input b_dim (per-cell mean dimension slopes).');
b_sigma_dim = psyrat_req(varargin,'b_sigma_dim','Please input b_sigma_dim (per-cell log residual-SD dimension slopes).');
sd_id       = psyrat_req(varargin,'sd_id','Please input sd_id (joint person SDs, draws x 8).');
cor_id      = psyrat_req(varargin,'cor_id','Please input cor_id (joint person correlations, draws x 8 x 8).');
sd_trl      = psyrat_req(varargin,'sd_trl','Please input sd_trl (per-cell trial SDs, draws x 4).');
cor_trl     = psyrat_req(varargin,'cor_trl','Please input cor_trl (cross-cell trial correlations).');
sd_occ      = psyrat_req(varargin,'sd_occ','Please input sd_occ (per-cell occasion SDs).');
cor_occ     = psyrat_req(varargin,'cor_occ','Please input cor_occ (cross-cell occasion correlations).');
sd_tid      = psyrat_req(varargin,'sd_tid','Please input sd_tid (person x trial SDs).');
cor_tid     = psyrat_req(varargin,'cor_tid','Please input cor_tid (person x trial correlations).');
sd_oid      = psyrat_req(varargin,'sd_oid','Please input sd_oid (person x occasion SDs).');
cor_oid     = psyrat_req(varargin,'cor_oid','Please input cor_oid (person x occasion correlations).');
sd_to       = psyrat_req(varargin,'sd_to','Please input sd_to (trial x occasion SDs).');
cor_to      = psyrat_req(varargin,'cor_to','Please input cor_to (trial x occasion correlations).');
mu_ss       = local_ss_block(varargin,'mu_ss','per-subject natural-scale means');
er_ss       = local_ss_block(varargin,'er_ss','per-subject absolute log residual SDs');
idtable     = psyrat_req(varargin,'idtable','Please input idtable (ssinfo).');
nprime      = psyrat_req(varargin,'nprime','Please input nprime (standardized per-cell trial counts).');
nprime_o    = psyrat_req(varargin,'nprime_o','Please input nprime_o (standardized per-cell occasion counts).');
reltype     = psyrat_req(varargin,'reltype','Please input reltype (1 = CE, 2 = CS, 3 = CES).');
ndim        = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
ciperc      = psyrat_req(varargin,'CI','Please input CI.');

ind = find(strcmpi('cvec',varargin),1);
if ~isempty(ind)
    cvec = varargin{ind+1};
else
    cvec = [1 -1 -1 1];
end

if ~any(ndim == [1 2])
    error('varargin:ndim','WARNING: ndim must be 1 or 2.\n');
end
if ~any(reltype == [1 2 3])
    error('varargin:reltype','WARNING: reltype must be 1 (CE), 2 (CS) or 3 (CES).\n');
end
ndraws = size(b,1);
if size(b,2) ~= 4 || size(b_sigma,2) ~= 4
    error('varargin:cells',... %Error code and associated error
        'WARNING: b and b_sigma must each have 4 columns, one per cell.\n');
end
if size(sd_id,2) ~= 8
    error('varargin:sd_id',... %Error code and associated error
        'WARNING: sd_id must have 8 columns [loc_1..4, logsd_1..4]; got %d.\n',...
        size(sd_id,2));
end
kdim = 2*ndim - 1;
if size(b_dim,2) ~= 4 || size(b_dim,3) ~= kdim || ...
        size(b_sigma_dim,2) ~= 4 || size(b_sigma_dim,3) ~= kdim
    error('varargin:b_dim',... %Error code and associated error
        'WARNING: b_dim and b_sigma_dim must be [ndraws x 4 x %d] for ndim = %d.\n',...
        kdim,ndim);
end
for counts = {{nprime,'nprime'},{nprime_o,'nprime_o'}}
    cn = counts{1}{1};
    if numel(cn) ~= 4 || any(cn <= 0) || any(isnan(cn))
        error(['varargin:' counts{1}{2}],... %Error code and associated error
            'WARNING: %s must be four positive per-cell counts.\n',counts{1}{2});
    end
end
nprime = nprime(:)'; nprime_o = nprime_o(:)';

nsub = height(idtable);
for q = 1:4
    if size(mu_ss{q},2) < max(idtable.id2) || size(er_ss{q},2) < max(idtable.id2)
        error('varargin:ss',... %Error code and associated error
            ['WARNING: mu_ss%d/er_ss%d have %d/%d subject columns but idtable '...
            'references id2 up to %d.\n'],q,q,...
            size(mu_ss{q},2),size(er_ss{q},2),max(idtable.id2));
    end
end

sd_p = sd_id(:,1:4);
cor_p = cor_id(:,1:4,1:4);

wantdraw = nargout > 1;
if wantdraw
    z1nan = nan(ndraws,nsub);
    drawlevel = struct('u_ces',z1nan,'G_ces',z1nan,'D_ces',z1nan,...
        'u_ce',z1nan,'G_ce',z1nan,'D_ce',z1nan,...
        'u_cs',z1nan,'G_cs',z1nan,'D_cs',z1nan,...
        'rel_ces',z1nan,'abs_ces',z1nan,'rel_ce',z1nan,'abs_ce',z1nan,...
        'rel_cs',z1nan,'abs_cs',z1nan,...
        'G_ces_std',z1nan,'D_ces_std',z1nan,'G_ce_std',z1nan,...
        'D_ce_std',z1nan,'G_cs_std',z1nan,'D_cs_std',z1nan,...
        'iccg_ces',z1nan,'iccd_ces',z1nan,'iccg_ce',z1nan,'iccd_ce',z1nan,...
        'iccg_cs',z1nan,'iccd_cs',z1nan,...
        'invalid',z1nan,'invalid_std',z1nan,...
        'expected_person',nan(ndraws,4,nsub),'expected_pop',nan(ndraws,4,nsub),...
        'sigma_person',nan(ndraws,4,nsub));
end

rows = cell(nsub,1);
ciedge = (1-ciperc)/2;

for c = 1:nsub

    if ndim == 1
        zt = idtable.z1(c);
    else
        zt = [idtable.z1(c); idtable.z2(c); idtable.z1(c)*idtable.z2(c)];
    end
    zt = zt(:);
    col = idtable.id2(c);

    alpha = zeros(ndraws,4);
    lsig_own = zeros(ndraws,4);
    mu_own = zeros(ndraws,4);
    for q = 1:4
        dimterm = reshape(b_dim(:,q,:),ndraws,kdim) * zt;
        sigterm = reshape(b_sigma_dim(:,q,:),ndraws,kdim) * zt;
        alpha(:,q) = b(:,q) + dimterm;
        mu_own(:,q) = mu_ss{q}(:,col) + dimterm;
        lsig_own(:,q) = er_ss{q}(:,col) + sigterm;
    end

    mats = psyrat_dod_components('facet','trial_occasion',...
        'alpha',alpha,'log_sigma',lsig_own,...
        'sd_p',sd_p,'cor_p',cor_p,'sd_trl',sd_trl,'cor_trl',cor_trl,...
        'sd_occ',sd_occ,'cor_occ',cor_occ,'sd_tid',sd_tid,'cor_tid',cor_tid,...
        'sd_oid',sd_oid,'cor_oid',cor_oid,'sd_to',sd_to,'cor_to',cor_to);

    %floor exactly the empty-cell ZERO (B16): a zero trial OR occasion
    %cell would trip psyrat_dod_composite's n > 0 guard. ONLY zero is
    %floored -- deliberately NOT max(1,...): the two-input max omits NaN,
    %so a max floor would silently convert corrupt counts to n = 1 and
    %hide them from the composite guard, which must keep erroring loudly
    %on them. The floor applies to the composite inputs only; ssrel keeps
    %the true counts. Trials and occasions are counted SEPARATELY -- they
    %are different facets and a reader needs to know which one was empty.
    ni_act = [idtable.trls1(c) idtable.trls2(c) idtable.trls3(c) idtable.trls4(c)];
    n_floored = nnz(ni_act == 0);
    ni_act(ni_act == 0) = 1;
    no_act = [idtable.occs1(c) idtable.occs2(c) idtable.occs3(c) idtable.occs4(c)];
    n_floored_occ = nnz(no_act == 0);
    no_act(no_act == 0) = 1;

    dA = local_design(mats,cvec,ni_act,no_act);
    dS = local_design(mats,cvec,nprime,nprime_o);
    singles = local_singles(mats,cvec);

    %identity link: the expected score is the additive predictor itself
    %(facet effects have mean zero additively - no +0.5*variance lognormal
    %correction), and no nu read-out exists for this family
    sig_own = exp(lsig_own);
    E_own = mu_own;
    expected_dod_pop = mats.expected * cvec(:);
    expected_dod_own = E_own * cvec(:);
    res_var_dod = sum(sig_own.^2, 2);

    if wantdraw
        for dn = {'ces','ce','cs'}
            drawlevel.(sprintf('u_%s',dn{1}))(:,c) = dA.(dn{1}).u;
            drawlevel.(sprintf('G_%s',dn{1}))(:,c) = dA.(dn{1}).G;
            drawlevel.(sprintf('D_%s',dn{1}))(:,c) = dA.(dn{1}).D;
            drawlevel.(sprintf('rel_%s',dn{1}))(:,c) = dA.(dn{1}).rel;
            drawlevel.(sprintf('abs_%s',dn{1}))(:,c) = dA.(dn{1}).abs_;
            drawlevel.(sprintf('G_%s_std',dn{1}))(:,c) = dS.(dn{1}).G;
            drawlevel.(sprintf('D_%s_std',dn{1}))(:,c) = dS.(dn{1}).D;
            drawlevel.(sprintf('iccg_%s',dn{1}))(:,c) = singles.(dn{1}).iccg;
            drawlevel.(sprintf('iccd_%s',dn{1}))(:,c) = singles.(dn{1}).iccd;
        end
        drawlevel.invalid(:,c) = dA.invalid;
        drawlevel.invalid_std(:,c) = dS.invalid;
        drawlevel.expected_person(:,:,c) = E_own;
        drawlevel.expected_pop(:,:,c) = mats.expected;
        drawlevel.sigma_person(:,:,c) = sig_own;
    end

    dn = local_reltype_name(reltype);
    row = idtable(c,:);
    row = local_summarize(row,'gen',dA.(dn).G,ciedge);
    row = local_summarize(row,'dep',dA.(dn).D,ciedge);
    row = local_summarize(row,'iccg',singles.(dn).iccg,ciedge);
    row = local_summarize(row,'iccd',singles.(dn).iccd,ciedge);
    row = local_summarize(row,'sem_gen',local_sem(dA.(dn).rel),ciedge);
    row = local_summarize(row,'sem_dep',local_sem(dA.(dn).abs_),ciedge);
    row = local_summarize(row,'gen',dS.(dn).G,ciedge,'_std');
    row = local_summarize(row,'dep',dS.(dn).D,ciedge,'_std');
    row = local_summarize(row,'sem_gen',local_sem(dS.(dn).rel),ciedge,'_std');
    row = local_summarize(row,'sem_dep',local_sem(dS.(dn).abs_),ciedge,'_std');
    row.uni_pt = mean(dA.(dn).u,'omitnan');
    row.uni_pt_std = mean(dS.(dn).u,'omitnan');
    row.relerr_pt = mean(dA.(dn).rel,'omitnan');
    row.abserr_pt = mean(dA.(dn).abs_,'omitnan');
    row.relerr_pt_std = mean(dS.(dn).rel,'omitnan');
    row.abserr_pt_std = mean(dS.(dn).abs_,'omitnan');
    row.expected_dod_pop_pt = mean(expected_dod_pop,'omitnan');
    row.expected_dod_person_pt = mean(expected_dod_own,'omitnan');
    row.residual_var_dod_person_pt = mean(res_var_dod,'omitnan');
    for q = 1:4
        row.(sprintf('ressd%d_pt',q)) = mean(sig_own(:,q),'omitnan');
    end
    row.invalid_frac = mean(dA.invalid);
    row.invalid_frac_std = mean(dS.invalid);
    row.n_floored = n_floored;
    row.n_floored_occ = n_floored_occ;
    row.reltype = reltype;
    rows{c} = row;
end

ssrel = vertcat(rows{:});

end

function ss = local_ss_block(args,stem,what)
ss = cell(1,4);
for q = 1:4
    ss{q} = psyrat_req(args,sprintf('%s%d',stem,q),...
        sprintf('Please input %s%d (%s, cell %d).',stem,q,what,q));
end
end

function d = local_design(mats,cvec,ni,no)
%One D-study design: the six harmonic-scaled terms with their facet-specific
%counts (trial terms / n_i, occasion terms / n_o, product terms / n_i*n_o),
%then the three decompositions. The invalid flag pools every term and every
%decomposition - a property of the draw, not of the selected reltype.
u = psyrat_dod_composite('varcov',mats.p,'cvec',cvec,'rule','unscaled');
t.pi  = psyrat_dod_composite('varcov',mats.pi,'cvec',cvec,'rule','harmonic','n',ni);
t.po  = psyrat_dod_composite('varcov',mats.po,'cvec',cvec,'rule','harmonic','n',no);
t.pio = psyrat_dod_composite('varcov',mats.pio_e,'cvec',cvec,'rule','harmonic','n',ni.*no);
t.i   = psyrat_dod_composite('varcov',mats.i,'cvec',cvec,'rule','harmonic','n',ni);
t.o   = psyrat_dod_composite('varcov',mats.o,'cvec',cvec,'rule','harmonic','n',no);
t.io  = psyrat_dod_composite('varcov',mats.io,'cvec',cvec,'rule','harmonic','n',ni.*no);

p = u.value;
d.ces.u = p;
d.ces.rel = t.pi.value + t.po.value + t.pio.value;
d.ces.abs_ = d.ces.rel + t.i.value + t.o.value + t.io.value;
d.ce.u = p + t.po.value;
d.ce.rel = t.pi.value + t.pio.value;
d.ce.abs_ = d.ce.rel + t.i.value + t.io.value;
d.cs.u = p + t.pi.value;
d.cs.rel = t.po.value + t.pio.value;
d.cs.abs_ = d.cs.rel + t.o.value + t.io.value;

tol = 1e-10;
badcoef = @(x) ~isfinite(x) | x < -tol | x > 1 + tol;
invalid = u.invalid_negative;
for f = {'pi','po','pio','i','o','io'}
    invalid = invalid | t.(f{1}).invalid_negative;
end
for dn = {'ces','ce','cs'}
    d.(dn{1}).G = local_ratio(d.(dn{1}).u, d.(dn{1}).u + d.(dn{1}).rel);
    d.(dn{1}).D = local_ratio(d.(dn{1}).u, d.(dn{1}).u + d.(dn{1}).abs_);
    invalid = invalid | d.(dn{1}).rel < 0 | d.(dn{1}).abs_ < 0 | ...
        badcoef(d.(dn{1}).G) | badcoef(d.(dn{1}).D);
end
d.invalid = double(invalid);
end

function s = local_singles(mats,cvec)
%Matched-facet ICCs at n' = 1: unscaled composites.
p   = psyrat_dod_composite('varcov',mats.p,'cvec',cvec,'rule','unscaled');
spi = psyrat_dod_composite('varcov',mats.pi,'cvec',cvec,'rule','unscaled');
spo = psyrat_dod_composite('varcov',mats.po,'cvec',cvec,'rule','unscaled');
spio= psyrat_dod_composite('varcov',mats.pio_e,'cvec',cvec,'rule','unscaled');
si  = psyrat_dod_composite('varcov',mats.i,'cvec',cvec,'rule','unscaled');
so  = psyrat_dod_composite('varcov',mats.o,'cvec',cvec,'rule','unscaled');
sio = psyrat_dod_composite('varcov',mats.io,'cvec',cvec,'rule','unscaled');
pv = p.value;
relsum = spi.value + spo.value + spio.value;
s.ces.iccg = local_ratio(pv, pv + relsum);
s.ces.iccd = local_ratio(pv, pv + relsum + si.value + so.value + sio.value);
s.ce.iccg = local_ratio(pv + spo.value, pv + spo.value + spi.value + spio.value);
s.ce.iccd = local_ratio(pv + spo.value, ...
    pv + spo.value + spi.value + spio.value + si.value + sio.value);
s.cs.iccg = local_ratio(pv + spi.value, pv + spi.value + spo.value + spio.value);
s.cs.iccd = local_ratio(pv + spi.value, ...
    pv + spi.value + spo.value + spio.value + so.value + sio.value);
end

function dn = local_reltype_name(reltype)
switch reltype
    case 1, dn = 'ce';
    case 2, dn = 'cs';
    case 3, dn = 'ces';
end
end

function r = local_ratio(num,den)
r = num ./ den;
r(~isfinite(num) | ~isfinite(den) | num < 0 | den <= 0) = NaN;
end

function s = local_sem(v)
s = nan(size(v));
ok = isfinite(v) & v >= 0;
s(ok) = sqrt(v(ok));
end

function row = local_summarize(row,stem,draws,ciedge,suffix)
if nargin < 5, suffix = ''; end
q = quantile(draws,[ciedge 1-ciedge]);
row.(sprintf('%s_pt%s',stem,suffix)) = mean(draws,'omitnan');
row.(sprintf('%s_ll%s',stem,suffix)) = q(1);
row.(sprintf('%s_ul%s',stem,suffix)) = q(2);
end
