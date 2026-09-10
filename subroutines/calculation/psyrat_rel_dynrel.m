function out = psyrat_rel_dynrel(varargin)
%Dynamic/conditional one-facet reliability as a function of the standardized
%dimension predictors (Rast & Clayson, analysis 11).
%
% out = psyrat_rel_dynrel('sig_p',sig_p,'sig_log0',sig_log0,'b_sigma',b_sigma,...
%   'sig_i',sig_i,'ndim',ndim,'z1',z1axis,'z2',z2axis,'obs',nprime,'CI',.95)
%
%The dynamic reliability model estimates a dimension-conditioned RELATIVE
%residual SD
%
%   sigma_pi,e(z) = exp(Intercept_sigma + b_sigma . X(z))      (delta_p = 0)
%
%where X(z) is the dimension design row: [z1] for one dimension; [z1 z2 z1*z2]
%for two dimensions. Reliability is then the one-facet coefficient evaluated at
%each z, reusing the validated psyrat_rel_sing with the production
%total-within convention wp = sqrt(sigma_pi,e(z)^2 + sigma_i^2), i = sigma_i:
%
%   G(z) = sigma_p^2 / (sigma_p^2 + sigma_pi,e(z)^2 / n')                (Rast 2.2.1)
%   D(z) = sigma_p^2 / (sigma_p^2 + (sigma_pi,e(z)^2 + sigma_i^2) / n')  (D adds sigma_i)
%
%so generalizability excludes, and dependability includes, the trial main
%effect, exactly as Rocha Table 2 / the static one-facet path.
%
%Inputs (name-value)
% sig_p    - between-person (universe-score) SD posterior draws (draws x 1).
%            This is gro_sds(:,1); per the owner decision it is CONDITIONAL on
%            the dimensions (they enter the mean submodel).
% sig_log0 - log-residual intercept draws Intercept_sigma (draws x 1): the log
%            relative residual SD at z = 0. z = 0 is the OVERALL-sample mean of
%            each dimension (standardized once across the whole sample, pooled
%            across any groups; see psyrat_zscore_subjectlevel), not a per-group
%            mean.
% b_sigma  - log-residual (scale) dimension-slope draws (draws x KDIM).
% sig_i    - trial main-effect SD draws sigma_i = sig_trl (draws x 1).
% ndim     - number of dimensions: 1 or 2.
% z1       - vector of standardized values for dimension 1 (the grid axis).
% obs      - n' (trial count) at which to evaluate the coefficients (scalar).
% CI       - credible-interval width in (0,1], e.g., .95.
%
%Optional Inputs
% z2       - vector of standardized values for dimension 2 (required if ndim=2).
% sig_dp   - person log-residual SD draws gro_sds(:,2). Required only when
%            'marginal' is 1.
% marginal - 0 (default): evaluate at the typical person (delta_p = 0).
%            1: lognormal population-average residual, multiplying the relative
%            residual SD by exp(sigma_delta_p^2) before forming the coefficient.
%
%Output (struct)
% out.ndim, out.obs, out.marginal, out.ci
% out.z1            - dimension-1 grid axis (1 x nz1)
% out.z2            - dimension-2 grid axis (1 x nz2; only when ndim=2)
% out.G, out.D      - structs with fields pt/ll/ul (global coefficients over the
%                     grid). ndim=1 -> 1 x nz1; ndim=2 -> nz1 x nz2 (rows index
%                     z1, columns index z2).
% out.ICCg, out.ICCd- relative/absolute single-observation ICC over the grid.
% out.sige_pt       - posterior-mean RELATIVE residual SD sigma_pi,e(z) over the
%                     grid (for reporting the scale surface).

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
    error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',...
        'Make sure each variable input is paired with a value \n',...
        'See help psyrat_rel_dynrel for more information about inputs'));
end

sig_p    = psyrat_req(varargin,'sig_p','Please input sig_p (between-person SD draws).');
sig_log0 = psyrat_req(varargin,'sig_log0','Please input sig_log0 (Intercept_sigma draws).');
b_sigma  = psyrat_req(varargin,'b_sigma','Please input b_sigma (scale dimension slopes).');
sig_i    = psyrat_req(varargin,'sig_i','Please input sig_i (trial main-effect SD draws).');
ndim     = psyrat_req(varargin,'ndim','Please input ndim (1 or 2).');
z1       = psyrat_req(varargin,'z1','Please input z1 (dimension-1 grid axis).');
obs      = psyrat_req(varargin,'obs','Please input obs (n'' trial count).');
ciperc   = psyrat_req(varargin,'CI','Please input CI.');

z2       = psyrat_opt(varargin,'z2',[]);
sig_dp   = psyrat_opt(varargin,'sig_dp',[]);
marginal = psyrat_opt(varargin,'marginal',0);

if ~any(ndim == [1 2])
    error('varargin:ndim','WARNING: ndim must be 1 or 2.\n');
end
if ndim == 2 && isempty(z2)
    error('varargin:z2','WARNING: z2 is required when ndim = 2.\n');
end
if ~isscalar(obs) || obs < 1
    error('varargin:obs','WARNING: obs (n'') must be a positive scalar.\n');
end
if marginal == 1 && isempty(sig_dp)
    error('varargin:sig_dp',...
        'WARNING: sig_dp (person log-residual SD draws) is required when marginal = 1.\n');
end

%column vectors of draws
sig_p    = sig_p(:);
sig_log0 = sig_log0(:);
sig_i    = sig_i(:);
if isempty(sig_dp); sig_dp = zeros(size(sig_log0)); else; sig_dp = sig_dp(:); end

%b_sigma must be draws x KDIM (KDIM = 1 for one dimension; 3 for two). The
%expectedk column-count check below catches a transposed/mis-shaped input.
kdim = size(b_sigma,2);
expectedk = (ndim == 1) * 1 + (ndim == 2) * 3;
if kdim ~= expectedk
    error('varargin:bsigma',...
        'WARNING: b_sigma has %d columns but ndim=%d expects %d.\n',...
        kdim,ndim,expectedk);
end

z1 = z1(:)';
nz1 = numel(z1);
if ndim == 2
    z2 = z2(:)';
    nz2 = numel(z2);
else
    nz2 = 1;
end

%preallocate the coefficient surfaces (rows index z1, columns index z2)
gridops = psyrat_dynrel_grid();  %shared grid-shaping helpers (empty/squeeze)
G = gridops.empty(nz1,nz2);
D = gridops.empty(nz1,nz2);
ICCg = gridops.empty(nz1,nz2);
ICCd = gridops.empty(nz1,nz2);
sige_pt = nan(nz1,nz2);

for a = 1:nz1
    for c = 1:nz2
        %dimension design row matching the Stan model design
        if ndim == 1
            xrow = z1(a);
        else
            xrow = [z1(a) z2(c) z1(a)*z2(c)];
        end

        %relative residual SD at this z (typical person, delta_p = 0)
        sige_z = exp(sig_log0 + b_sigma * xrow(:));

        %optional lognormal population-average over the person scale random
        %effect: marginal residual SD = exp(alpha + b_sigma z) * exp(sigma_dp^2)
        if marginal == 1
            sige_z = sige_z .* exp(sig_dp.^2);
        end

        sige_pt(a,c) = mean(sige_z);

        %total within = sqrt(sigma_pi,e(z)^2 + sigma_i^2); psyrat_rel_sing
        %subtracts i^2 to recover the relative residual for G and keeps it for D
        %(production one-facet convention, psyrat_relsummary local_total_within).
        wp_z = sqrt(sige_z.^2 + sig_i.^2);

        [G.ll(a,c),G.pt(a,c),G.ul(a,c)] = psyrat_rel_sing(...
            'gcoeff',2,'metric','global','bp',sig_p,'wp',wp_z,'i',sig_i,...
            'obs',obs,'CI',ciperc);
        [D.ll(a,c),D.pt(a,c),D.ul(a,c)] = psyrat_rel_sing(...
            'gcoeff',1,'metric','global','bp',sig_p,'wp',wp_z,'i',sig_i,...
            'obs',obs,'CI',ciperc);
        [ICCg.ll(a,c),ICCg.pt(a,c),ICCg.ul(a,c)] = psyrat_rel_sing(...
            'gcoeff',2,'metric','icc','bp',sig_p,'wp',wp_z,'i',sig_i,'CI',ciperc);
        [ICCd.ll(a,c),ICCd.pt(a,c),ICCd.ul(a,c)] = psyrat_rel_sing(...
            'gcoeff',1,'metric','icc','bp',sig_p,'wp',wp_z,'i',sig_i,'CI',ciperc);
    end
end

out = struct();
out.ndim = ndim;
out.obs = obs;
out.marginal = marginal;
out.ci = ciperc;
out.z1 = z1;
if ndim == 2
    out.z2 = z2;
end
out.sige_pt = gridops.squeeze1d(sige_pt,ndim);
out.G    = gridops.squeeze(G,ndim);
out.D    = gridops.squeeze(D,ndim);
out.ICCg = gridops.squeeze(ICCg,ndim);
out.ICCd = gridops.squeeze(ICCd,ndim);

end
