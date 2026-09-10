function out = psyrat_rel_dodiffdynrel_trt_gamma_ls(varargin)
%Dynamic/conditional reliability SURFACE of the four-cell NONCONCURRENT
%difference-of-differences for the TWO-FACET (trials x occasions) gamma
%location-scale design (analysis 29, family = 'gamma', gammascale = 2).
%
% out = psyrat_rel_dodiffdynrel_trt_gamma_ls('b',b,...,'sd_occ',sd_occ,...
%   'nocc',nocc,'reltype',reltype,...)
%
%The two-facet twin of PSYRAT_REL_DODIFFDYNREL_GAMMA_LS - read that file
%first for the TYPICAL-PERSON estimand statement (all person effects at zero;
%owner-approved 2026-08-12; no counterpart in the reference bundle, which
%computes only per-participant tables). This twin adds the seven-component
%structure and the three decompositions of the per-participant twin
%PSYRAT_SSREL_DODIFFDYNREL_TRT_GAMMA_LS: reltype 1 = CE, 2 = CS, 3 = CES,
%with the reltype-invariant draw-level invalidity flag.
%
%Additional inputs beyond the one-facet surface:
% sd_occ/cor_occ, sd_tid/cor_tid, sd_oid/cor_oid, sd_to/cor_to - facet blocks
% nocc    - [1 x 4] per-cell occasion counts n'_o
% reltype - 1 (CE), 2 (CS), or 3 (CES)
%
%See also PSYRAT_REL_DODIFFDYNREL_GAMMA_LS,
%PSYRAT_SSREL_DODIFFDYNREL_TRT_GAMMA_LS, PSYRAT_DOD_COMPOSITE.

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
        error('psyrat_rel_dodiffdynrel_trt_gamma_ls:noresidualcoupling',... %Error code and associated error
            ['''%s'' is not an input: all six observation-level residual '...
            'covariances are fixed to zero by the nonconcurrent estimand.'],...
            coupled{1});
    end
end

b           = psyrat_req(varargin,'b','Please input b (per-cell absolute log-mean intercepts).');
b_sigma     = psyrat_req(varargin,'b_sigma','Please input b_sigma (per-cell absolute log residual-SD intercepts).');
b_dim       = psyrat_req(varargin,'b_dim','Please input b_dim (per-cell log-mean dimension slopes).');
b_sigma_dim = psyrat_req(varargin,'b_sigma_dim','Please input b_sigma_dim (per-cell log residual-SD dimension slopes).');
sd_id       = psyrat_req(varargin,'sd_id','Please input sd_id (joint person SDs, draws x 8).');
cor_id      = psyrat_req(varargin,'cor_id','Please input cor_id (joint person correlations, draws x 8 x 8).');
sd_trl      = psyrat_req(varargin,'sd_trl','Please input sd_trl (per-cell trial SDs).');
cor_trl     = psyrat_req(varargin,'cor_trl','Please input cor_trl (cross-cell trial correlations).');
sd_occ      = psyrat_req(varargin,'sd_occ','Please input sd_occ (per-cell occasion SDs).');
cor_occ     = psyrat_req(varargin,'cor_occ','Please input cor_occ (cross-cell occasion correlations).');
sd_tid      = psyrat_req(varargin,'sd_tid','Please input sd_tid (person x trial SDs).');
cor_tid     = psyrat_req(varargin,'cor_tid','Please input cor_tid (person x trial correlations).');
sd_oid      = psyrat_req(varargin,'sd_oid','Please input sd_oid (person x occasion SDs).');
cor_oid     = psyrat_req(varargin,'cor_oid','Please input cor_oid (person x occasion correlations).');
sd_to       = psyrat_req(varargin,'sd_to','Please input sd_to (trial x occasion SDs).');
cor_to      = psyrat_req(varargin,'cor_to','Please input cor_to (trial x occasion correlations).');
ndim        = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
z1          = psyrat_req(varargin,'z1','Please input z1 (grid axis).');
obs         = psyrat_req(varargin,'obs','Please input obs (per-cell trial counts n''_i).');
nocc        = psyrat_req(varargin,'nocc','Please input nocc (per-cell occasion counts n''_o).');
reltype     = psyrat_req(varargin,'reltype','Please input reltype (1 = CE, 2 = CS, 3 = CES).');
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
if ndim == 2
    z2 = psyrat_req(varargin,'z2','Please input z2 (second grid axis).');
end
ndraws = size(b,1);
kdim = 2*ndim - 1;
if size(b_dim,2) ~= 4 || size(b_dim,3) ~= kdim || ...
        size(b_sigma_dim,2) ~= 4 || size(b_sigma_dim,3) ~= kdim
    error('varargin:b_dim',... %Error code and associated error
        'WARNING: b_dim and b_sigma_dim must be [ndraws x 4 x %d] for ndim = %d.\n',...
        kdim,ndim);
end
for counts = {{obs,'obs'},{nocc,'nocc'}}
    cn = counts{1}{1};
    if numel(cn) ~= 4 || any(cn <= 0) || any(isnan(cn))
        error(['varargin:' counts{1}{2}],... %Error code and associated error
            'WARNING: %s must be four positive per-cell counts.\n',counts{1}{2});
    end
end
obs = obs(:)'; nocc = nocc(:)';

sd_p = sd_id(:,1:4);
cor_p = cor_id(:,1:4,1:4);

z1 = z1(:)'; nz1 = numel(z1);
if ndim == 2; z2 = z2(:)'; nz2 = numel(z2); else; nz2 = 1; end

gridops = psyrat_dynrel_grid();
G = gridops.empty(nz1,nz2);  D = gridops.empty(nz1,nz2);
ICCg = gridops.empty(nz1,nz2);  ICCd = gridops.empty(nz1,nz2);
sige_pt = nan(nz1,nz2);
invalid_frac = nan(nz1,nz2);
ciedge = (1-ciperc)/2;

for a = 1:nz1
    for c = 1:nz2
        if ndim == 1
            zt = z1(a);
        else
            zt = [z1(a); z2(c); z1(a)*z2(c)];
        end
        zt = zt(:);

        alpha = zeros(ndraws,4);
        lsig = zeros(ndraws,4);
        for q = 1:4
            alpha(:,q) = b(:,q) + reshape(b_dim(:,q,:),ndraws,kdim) * zt;
            lsig(:,q) = b_sigma(:,q) + reshape(b_sigma_dim(:,q,:),ndraws,kdim) * zt;
        end

        mats = psyrat_gamma_dod_components_ls('facet','trial_occasion',...
            'alpha',alpha,'log_sigma',lsig,...
            'sd_p',sd_p,'cor_p',cor_p,'sd_trl',sd_trl,'cor_trl',cor_trl,...
            'sd_occ',sd_occ,'cor_occ',cor_occ,'sd_tid',sd_tid,'cor_tid',cor_tid,...
            'sd_oid',sd_oid,'cor_oid',cor_oid,'sd_to',sd_to,'cor_to',cor_to);

        u = psyrat_dod_composite('varcov',mats.p,'cvec',cvec,'rule','unscaled');
        t.pi  = psyrat_dod_composite('varcov',mats.pi,'cvec',cvec,'rule','harmonic','n',obs);
        t.po  = psyrat_dod_composite('varcov',mats.po,'cvec',cvec,'rule','harmonic','n',nocc);
        t.pio = psyrat_dod_composite('varcov',mats.pio_e,'cvec',cvec,'rule','harmonic','n',obs.*nocc);
        t.i   = psyrat_dod_composite('varcov',mats.i,'cvec',cvec,'rule','harmonic','n',obs);
        t.o   = psyrat_dod_composite('varcov',mats.o,'cvec',cvec,'rule','harmonic','n',nocc);
        t.io  = psyrat_dod_composite('varcov',mats.io,'cvec',cvec,'rule','harmonic','n',obs.*nocc);

        p = u.value;
        %ALL THREE decompositions are formed so the invalidity flag can pool
        %them (reltype-invariant, exactly as the per-participant twin does) -
        %pooling only the selected one is invariant today by an algebraic
        %argument, but a decomposition-specific edit would silently break it
        %(pre-merge review finding, 2026-08-13).
        dall = struct();
        dall.ce.uni = p + t.po.value;
        dall.ce.rel = t.pi.value + t.pio.value;
        dall.ce.abs_ = dall.ce.rel + t.i.value + t.io.value;
        dall.cs.uni = p + t.pi.value;
        dall.cs.rel = t.po.value + t.pio.value;
        dall.cs.abs_ = dall.cs.rel + t.o.value + t.io.value;
        dall.ces.uni = p;
        dall.ces.rel = t.pi.value + t.po.value + t.pio.value;
        dall.ces.abs_ = dall.ces.rel + t.i.value + t.o.value + t.io.value;
        switch reltype
            case 1  %CE: the scaled person x occasion term joins the universe
                uni = p + t.po.value;
                rel = t.pi.value + t.pio.value;
                abs_ = rel + t.i.value + t.io.value;
            case 2  %CS: the scaled person x trial term joins the universe
                uni = p + t.pi.value;
                rel = t.po.value + t.pio.value;
                abs_ = rel + t.o.value + t.io.value;
            case 3  %CES
                uni = p;
                rel = t.pi.value + t.po.value + t.pio.value;
                abs_ = rel + t.i.value + t.o.value + t.io.value;
        end
        Gd = local_ratio(uni, uni + rel);
        Dd = local_ratio(uni, uni + abs_);

        %matched-facet ICCs at n' = 1 (unscaled singles), per reltype
        spi = psyrat_dod_composite('varcov',mats.pi,'cvec',cvec,'rule','unscaled');
        spo = psyrat_dod_composite('varcov',mats.po,'cvec',cvec,'rule','unscaled');
        spio= psyrat_dod_composite('varcov',mats.pio_e,'cvec',cvec,'rule','unscaled');
        si  = psyrat_dod_composite('varcov',mats.i,'cvec',cvec,'rule','unscaled');
        so  = psyrat_dod_composite('varcov',mats.o,'cvec',cvec,'rule','unscaled');
        sio = psyrat_dod_composite('varcov',mats.io,'cvec',cvec,'rule','unscaled');
        switch reltype
            case 1
                iuni = p + spo.value;
                irel = spi.value + spio.value;
                iabs = irel + si.value + sio.value;
            case 2
                iuni = p + spi.value;
                irel = spo.value + spio.value;
                iabs = irel + so.value + sio.value;
            case 3
                iuni = p;
                irel = spi.value + spo.value + spio.value;
                iabs = irel + si.value + so.value + sio.value;
        end
        ig = local_ratio(iuni, iuni + irel);
        id_ = local_ratio(iuni, iuni + iabs);

        %the draw-level invalidity flag pools every term and ALL THREE
        %decompositions' error/coefficient checks (reltype-invariant by
        %construction, matching the per-participant twin and the reference)
        tol = 1e-10;
        badcoef = @(x) ~isfinite(x) | x < -tol | x > 1 + tol;
        inval = u.invalid_negative;
        for f = {'pi','po','pio','i','o','io'}
            inval = inval | t.(f{1}).invalid_negative;
        end
        for dn = {'ce','cs','ces'}
            dd = dall.(dn{1});
            Gx = local_ratio(dd.uni, dd.uni + dd.rel);
            Dx = local_ratio(dd.uni, dd.uni + dd.abs_);
            inval = inval | dd.rel < 0 | dd.abs_ < 0 | ...
                badcoef(Gx) | badcoef(Dx);
        end

        [G.ll(a,c),G.pt(a,c),G.ul(a,c)] = local_sum(Gd,ciedge);
        [D.ll(a,c),D.pt(a,c),D.ul(a,c)] = local_sum(Dd,ciedge);
        [ICCg.ll(a,c),ICCg.pt(a,c),ICCg.ul(a,c)] = local_sum(ig,ciedge);
        [ICCd.ll(a,c),ICCd.pt(a,c),ICCd.ul(a,c)] = local_sum(id_,ciedge);
        sige_pt(a,c) = mean(sqrt(max(spio.value,0)),'omitnan');
        invalid_frac(a,c) = mean(inval);
    end
end

out = struct();
out.ndim = ndim; out.obs = obs; out.nocc = nocc; out.reltype = reltype;
out.ci = ciperc;
out.z1 = z1;
if ndim == 2; out.z2 = z2; end
out.sige_pt = gridops.squeeze1d(sige_pt,ndim);
out.invalid_frac = gridops.squeeze1d(invalid_frac,ndim);
out.G = gridops.squeeze(G,ndim);
out.D = gridops.squeeze(D,ndim);
out.ICCg = gridops.squeeze(ICCg,ndim);
out.ICCd = gridops.squeeze(ICCd,ndim);
out.estimand = 'typical person (all person effects at zero)';

end

function r = local_ratio(num,den)
r = num ./ den;
r(~isfinite(num) | ~isfinite(den) | num < 0 | den <= 0) = NaN;
end

function [ll,pt,ul] = local_sum(draws,ciedge)
q = quantile(draws,[ciedge 1-ciedge]);
pt = mean(draws,'omitnan');
ll = q(1); ul = q(2);
end
