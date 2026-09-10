function out = psyrat_rel_diffdynrel_trt(varargin)
%Dynamic/conditional reliability of a trial + occasion two-facet difference
%score as a function of the standardized dimension predictors (Rast & Clayson,
%group-level / fixed-per-event residual; non-concurrent).
%
% out = psyrat_rel_diffdynrel_trt('id_varcov',id_varcov,'trl_varcov',trl_varcov,...
%   'occ_varcov',occ_varcov,'tid_varcov',tid_varcov,'oid_varcov',oid_varcov,...
%   'to_varcov',to_varcov,'b_sigma',b_sigma,'b_sigma_dim',b_sigma_dim,...
%   'ndim',ndim,'z1',z1axis,'z2',z2axis,'obs',ntrl,'nocc',nocc,...
%   'reltype',reltype,'CI',.95)
%
%The group-level two-facet dynamic difference model estimates per-event
%log-residual SDs that depend on the standardized dimensions:
%
%   log sigma_pto,e,event(z) = b_sigma[event] + b_sigma_dim[event] . z-terms
%
%so the gain-minus-loss two-facet difference reliability varies with z. This is
%computed by reusing the validated two-facet psyrat_diffrel_trt with the
%z-conditioned per-event residual er_var(z) = [log sigma_e,1(z), log sigma_e,2(z)]
%(psyrat_diffrel_trt exponentiates and squares internally). The six event
%variance-covariance blocks map to the psyrat_diffrel_trt inputs:
%
%   id_varcov  -> bp   (person)             tid_varcov -> bpi  (person x trial)
%   oid_varcov -> bpo  (person x occasion)  trl_varcov -> bt   (trial)
%   occ_varcov -> bo   (occasion)           to_varcov  -> boi  (trial x occasion)
%
%which are NOT z-dependent (the mean dimension slopes make the difference
%universe-score variance conditional on the dimensions), so the surface varies
%through the residual only. Residual covariance between the two events is not
%modeled here (er_cov = 0; non-concurrent). psyrat_diffrel_trt's three reltypes
%give the equivalence (occasion fixed), stability (trial fixed), and
%equivalence-and-stability (both random) difference coefficients.
%
%Inputs (name-value)
% id_varcov  - person event variance-covariance draws (draws x 2 x 2).
% trl_varcov - trial event variance-covariance draws (draws x 2 x 2).
% occ_varcov - occasion event variance-covariance draws (draws x 2 x 2).
% tid_varcov - person x trial event variance-covariance draws (draws x 2 x 2).
% oid_varcov - person x occasion event variance-covariance draws (draws x 2 x 2).
% to_varcov  - trial x occasion event variance-covariance draws (draws x 2 x 2).
% b_sigma    - per-event log-residual SD draws at z = 0 (draws x 2).
% b_sigma_dim- event-specific log-residual dimension-slope draws (draws x KDIM;
%              KDIM = 2 for one dimension, 6 for two). Odd columns are event 1,
%              even columns event 2, matching the case-15 Xdim design.
% ndim       - number of dimensions: 1 or 2.
% z1         - grid axis for dimension 1.
% obs        - n' (trial count) per event for the surface (scalar).
% CI         - credible-interval width in (0,1].
%
%Optional Inputs
% z2       - grid axis for dimension 2 (required if ndim = 2).
% nocc     - n' (occasion count) per event for the surface (scalar); default 1.
% reltype  - 1, 2, or 3 (see above); default 3 (both facets random).
%
%Output (struct)
% out.ndim, out.obs, out.nocc, out.reltype, out.ci
% out.z1 (, out.z2)        - grid axes
% out.G, out.D             - difference-score surfaces (ll/pt/ul); ndim=1 ->
%                            1 x nz1; ndim=2 -> nz1 x nz2 (rows z1, cols z2).
% out.ICCg, out.ICCd       - relative/absolute single-observation difference
%                            coefficient (the same reltype at obs=1, nocc=1).

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

if mod(length(varargin),2)
    error('varargin:incomplete',...
        strcat('WARNING: Inputs are incomplete \n\n',...
        'Make sure each variable input is paired with a value \n',...
        'See help psyrat_rel_diffdynrel_trt for more information about inputs'));
end

id_varcov  = psyrat_req(varargin,'id_varcov','Please input id_varcov.');
trl_varcov = psyrat_req(varargin,'trl_varcov','Please input trl_varcov.');
occ_varcov = psyrat_req(varargin,'occ_varcov','Please input occ_varcov.');
tid_varcov = psyrat_req(varargin,'tid_varcov','Please input tid_varcov.');
oid_varcov = psyrat_req(varargin,'oid_varcov','Please input oid_varcov.');
to_varcov  = psyrat_req(varargin,'to_varcov','Please input to_varcov.');
b_sigma    = psyrat_req(varargin,'b_sigma','Please input b_sigma (per-event log-residual SD).');
b_sig_dim  = psyrat_req(varargin,'b_sigma_dim','Please input b_sigma_dim (scale dimension slopes).');
ndim       = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
z1         = psyrat_req(varargin,'z1','Please input z1.');
obs        = psyrat_req(varargin,'obs','Please input obs (n'').');
ciperc     = psyrat_req(varargin,'CI','Please input CI.');
z2         = psyrat_opt(varargin,'z2',[]);
nocc       = psyrat_opt(varargin,'nocc',1);
reltype    = psyrat_opt(varargin,'reltype',3);
%concurrent variant (analysis 18): estimated residual-correlation draws. When
%supplied, the z-conditioned residual covariance rescor .* sigma_1(z) .* sigma_2(z)
%is fed to psyrat_diffrel_trt as er_cov; empty (default) -> 0 (non-concurrent).
rescor     = psyrat_opt(varargin,'rescor',[]);
if ~isempty(rescor); rescor = rescor(:); end

if ~any(ndim == [1 2])
    error('varargin:ndim','WARNING: ndim must be 1 or 2.\n');
end
if ndim == 2 && isempty(z2)
    error('varargin:z2','WARNING: z2 is required when ndim = 2.\n');
end
if ~isscalar(obs) || obs < 1
    error('varargin:obs','WARNING: obs (n'') must be a positive scalar.\n');
end
if ~isscalar(nocc) || nocc < 1
    error('varargin:nocc','WARNING: nocc must be a positive scalar.\n');
end
if ~any(reltype == [1 2 3])
    error('varargin:reltype','WARNING: reltype must be 1, 2, or 3.\n');
end

%expected number of scale slope columns: 2 events x (1 or 3 dimension terms)
expectedk = (ndim == 1) * 2 + (ndim == 2) * 6;
if size(b_sig_dim,2) ~= expectedk
    error('varargin:bsigmadim',...
        'WARNING: b_sigma_dim has %d columns but ndim=%d expects %d.\n',...
        size(b_sig_dim,2),ndim,expectedk);
end

%event-specific scale slopes: odd columns -> event 1, even columns -> event 2
ev1 = b_sig_dim(:,1:2:end); %draws x (KDIM/2)
ev2 = b_sig_dim(:,2:2:end);

z1 = z1(:)'; nz1 = numel(z1);
if ndim == 2; z2 = z2(:)'; nz2 = numel(z2); else; nz2 = 1; end

gridops = psyrat_dynrel_grid();  %shared grid-shaping helpers (empty/squeeze)
G = gridops.empty(nz1,nz2);  D = gridops.empty(nz1,nz2);
ICCg = gridops.empty(nz1,nz2);  ICCd = gridops.empty(nz1,nz2);

%shared static covariance blocks for psyrat_diffrel_trt (id->bp, tid->bpi,
%oid->bpo, trl->bt, occ->bo, to->boi).
blockargs = {'bp',id_varcov,'bpi',tid_varcov,'bpo',oid_varcov,...
    'bt',trl_varcov,'bo',occ_varcov,'boi',to_varcov};

%One warning per surface, not one per grid point: the kernel warns on every
%call and this loop calls it nz1*nz2 times. Suppress here, aggregate, report
%once. The gamma difference variant deliberately passes UNEQUAL per-event trial
%counts at every point, so without this a single surface would emit hundreds of
%identical blocks.
adm_ws = warning('off','psyrat:negerrvar');
adm_all = [];

for a = 1:nz1
    for c = 1:nz2
        %dimension-term vector matching the case-15 design
        if ndim == 1
            zterms = z1(a);
        else
            zterms = [z1(a); z2(c); z1(a)*z2(c)];
        end

        %z-conditioned per-event log-residual SD (draws x 2)
        er1 = b_sigma(:,1) + ev1 * zterms(:);
        er2 = b_sigma(:,2) + ev2 * zterms(:);
        er_var_z = [er1 er2];

        %concurrent variant: z-conditioned residual covariance from the estimated
        %residual correlation; 0 when non-concurrent (rescor empty)
        if isempty(rescor)
            er_cov_z = 0;
        else
            er_cov_z = rescor .* exp(er1) .* exp(er2);
        end

        g = psyrat_diffrel_trt(blockargs{:},'er_var',er_var_z,...
            'obs',[obs obs],'nocc',[nocc nocc],'reltype',reltype,...
            'est','gen','er_cov',er_cov_z,'CI',ciperc);
        d = psyrat_diffrel_trt(blockargs{:},'er_var',er_var_z,...
            'obs',[obs obs],'nocc',[nocc nocc],'reltype',reltype,...
            'est','dep','er_cov',er_cov_z,'CI',ciperc);

        %single-observation difference ICC (one trial, one occasion) is exactly
        %psyrat_diffrel_trt's icc_* output: at obs=[1 1]/nocc=[1 1] every
        %D-study-adjusted error component reduces (by exact division by 1) to the
        %same _raw component the icc is built from, so g.icc_*/d.icc_* are
        %bit-identical to a separate obs=[1 1] call. Reuse them instead of two
        %extra psyrat_diffrel_trt calls (matches the one-facet psyrat_rel_diffdynrel).
        G.ll(a,c) = g.ll; G.pt(a,c) = g.pt; G.ul(a,c) = g.ul;
        D.ll(a,c) = d.ll; D.pt(a,c) = d.pt; D.ul(a,c) = d.ul;
        ICCg.ll(a,c) = g.icc_ll; ICCg.pt(a,c) = g.icc_pt; ICCg.ul(a,c) = g.icc_ul;
        ICCd.ll(a,c) = d.icc_ll; ICCd.pt(a,c) = d.icc_pt; ICCd.ul(a,c) = d.icc_ul;

        if isfield(g,'admissibility')
            adm_all = [adm_all, g.admissibility, d.admissibility]; %#ok<AGROW>
        end
    end
end

warning(adm_ws);
psyrat_admissibility_warn(adm_all);

out = struct();
out.ndim = ndim; out.obs = obs; out.nocc = nocc; out.reltype = reltype;
out.ci = ciperc;
out.z1 = z1;
if ndim == 2; out.z2 = z2; end
out.G = gridops.squeeze(G,ndim);
out.D = gridops.squeeze(D,ndim);
out.ICCg = gridops.squeeze(ICCg,ndim);
out.ICCd = gridops.squeeze(ICCd,ndim);

end
