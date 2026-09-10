function dodiff = psyrat_dodiffrel(varargin)
%Calculate task-free 4-event difference-of-differences reliability
%
%dodiff = psyrat_dodiffrel('bp',bp_varcov,'bt',bt_varcov,...
%  'er_var',b_sigma,'obs',[n1 n2 n3 n4],'est','dep','CI',.95)
%
%Required Inputs:
% bp - between-person covariance draws (ndraws x 4 x 4)
% bt - between-trial covariance draws (ndraws x 4 x 4)
% er_var - log residual SD draws for each event (ndraws x 4)
% obs - number of trials for each event [n1 n2 n3 n4]
% est - 'dep' (dependability, absolute) or 'gen' (generalizability, relative)
% CI - credible interval width in [0,1]
%
%Optional Inputs:
% cvec - signed contrast vector (default [1 -1 -1 1]) corresponding to
%  (ERP_1 - ERP_2) - (ERP_3 - ERP_4)
%
%Output:
% dodiff - structure with reliability and ICC summaries
%
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

if mod(length(varargin),2)
    error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',...
        'Make sure each variable input is paired with a value \n',...
        'See help psyrat_dodiffrel for more information about inputs'));
end

bp = local_req(varargin,'bp','Between-person covariance draws (bp) were not specified.');
bt = local_req(varargin,'bt','Between-trial covariance draws (bt) were not specified.');
er_var = local_req(varargin,'er_var','Residual log-SD draws (er_var) were not specified.');
obs = local_req(varargin,'obs','Trial counts (obs) were not specified.');
est = local_req(varargin,'est','Reliability estimate type (est) was not specified.');
ciperc = local_req(varargin,'ci','Credible interval width (CI) was not specified.');

ind = find(strcmpi('cvec',varargin),1);
if ~isempty(ind)
    cvec = varargin{ind+1};
else
    cvec = [1 -1 -1 1];
end

if numel(obs) ~= 4
    error('varargin:obs',... %Error code and associated error
        'obs must include four event-level trial counts [n1 n2 n3 n4].');
end
obs = obs(:)';
if any(isnan(obs)) || any(obs <= 0)
    error('varargin:obs',... %Error code and associated error
        'All obs values must be numeric and > 0.');
end

if numel(cvec) ~= 4
    error('varargin:cvec',... %Error code and associated error
        'cvec must contain four signed contrast weights.');
end
cvec = cvec(:);

if ~any(strcmpi(est,{'dep','gen'}))
    error('varargin:est',... %Error code and associated error
        'est must be either ''dep'' or ''gen''.');
end

if ~isscalar(ciperc) || isnan(ciperc) || ciperc <= 0 || ciperc >= 1
    error('varargin:ci',... %Error code and associated error
        'CI must be a scalar in the open interval (0,1).');
end
ciedge = (1-ciperc)/2;

if ndims(bp) ~= 3 || size(bp,2) ~= 4 || size(bp,3) ~= 4
    error('varargin:bp',... %Error code and associated error
        'bp must be ndraws x 4 x 4 covariance draws.');
end
if ndims(bt) ~= 3 || size(bt,2) ~= 4 || size(bt,3) ~= 4
    error('varargin:bt',... %Error code and associated error
        'bt must be ndraws x 4 x 4 covariance draws.');
end

if size(bp,1) ~= size(bt,1)
    error('varargin:draws',... %Error code and associated error
        'bp and bt must have the same number of posterior draws.');
end

if size(er_var,2) ~= 4
    error('varargin:er_var',... %Error code and associated error
        'er_var must be ndraws x 4 log residual SD draws.');
end

if size(er_var,1) ~= size(bp,1)
    error('varargin:draws',... %Error code and associated error
        'er_var must have the same number of posterior draws as bp/bt.');
end

ndraws = size(bp,1);

%Residual variances for each event.
er_var = exp(er_var).^2;

%Universe score variance from between-person covariance matrix.
u = zeros(ndraws,1);
bt_unscaled = zeros(ndraws,1);
bt_scaled = zeros(ndraws,1);
%Residuals are independent across the four events, so the contrast's residual
%variance is sum_i c_i^2 * sigma^2_err_i (divided by obs_i for the mean score).
%The squared contrast weights make this correct for any contrast, not only the
%default +/-1 weights (for which c_i^2 = 1 and this reduces to a plain sum).
cw2 = cvec.^2;
er_scaled = er_var * (cw2 ./ obs(:));
er_unscaled = er_var * cw2;

%Build trial-scaled between-trial covariance matrix for each draw.
for d = 1:ndraws
    sbp = squeeze(bp(d,:,:));
    sbt = squeeze(bt(d,:,:));
    u(d) = cvec' * sbp * cvec;
    bt_unscaled(d) = cvec' * sbt * cvec;

    sbt_scaled = zeros(4,4);
    for i = 1:4
        for j = 1:4
            if i == j
                sbt_scaled(i,j) = sbt(i,j) / obs(i);
            else
                nharm = 2 / ((1/obs(i)) + (1/obs(j)));
                sbt_scaled(i,j) = sbt(i,j) / nharm;
            end
        end
    end
    bt_scaled(d) = cvec' * sbt_scaled * cvec;
end

dep_denom = u + bt_scaled + er_scaled;
gen_denom = u + er_scaled;
icc_dep_denom = u + bt_unscaled + er_unscaled;
icc_gen_denom = u + er_unscaled;

%Record whether the quantities that must be variances actually are. This kernel
%carries the same hazard as the two-measure difference kernels and is more
%exposed to it, not less: the off-diagonals of the between-trial block are
%divided by the harmonic mean of each PAIR of trial counts (above), and the
%four events here have four independently observed counts, so equal counts are
%the exception rather than the rule. u is a quadratic form in the unscaled
%person block and cannot be negative for positive semi-definite input; the
%trial-scaled contrast variance can be, and until now nothing said so.
%
%The ICC denominators (:151-152) use the UNSCALED contrast variances, and the
%two of them are NOT alike. bt_unscaled is a quadratic form in the RAW
%between-trial block, so the harmonic-mean scaler cannot reach it, but a raw
%block that is not positive semi-definite still can -- the gamma /
%scaled-chi-square family clamps variance diagonals at zero while leaving the
%matching cross-covariances alone (RC-23) -- and until RC-37 an out-of-range ICC
%was printed in silence. It is recorded.
%
%er_unscaled is NOT recorded, and must not be: er_var = exp(er_var).^2 is
%non-negative elementwise, cw2 = cvec.^2 is non-negative, and a non-negative
%matrix times a non-negative vector cannot be negative, so the record could
%never fire and would only add a line that says nothing. Do not "complete the
%set" here.
%
%er_scaled (the third record below) shares that non-negativity -- it is the same
%product with obs in the denominator, and obs > 0 is enforced at :65-67, so no
%subtraction enters the chain and its NEGATIVE branch is unreachable too. Its
%record is kept anyway, and the asymmetry with er_unscaled is deliberate on both
%counts (RC-39):
%
% - The ZERO branch is what earns it. er_scaled reaches exactly zero whenever
%   er_var does (a log-SD of -Inf, the log(0) encoding), when a contrast weight
%   is zero, or when obs is Inf, which :65 does not reject. Zero here is worth
%   saying out loud precisely because it is an ERROR variance sitting in a
%   denominator: it drives the coefficient toward ONE, not toward zero, which is
%   the opposite of the reading a user brings to "a variance came out zero".
%   That role-conditioned wording is RC-37's fix and is set explicitly per call
%   site, so this record inherits the right claim rather than a guessed one.
% - Recording er_unscaled ALONGSIDE it would not add a second check, because the
%   two vanish on exactly the same draws. It would report one failure as two.
%   Same reasoning that keeps uni_icc unrecorded at reltype 3, where it and uni
%   are assigned from the same variable and are bit-identical by construction.
%
%So: one record for the residual, taken on the scaled form that actually enters
%the coefficient denominator. u and bt_scaled earn their records on the negative
%branch instead, being genuine quadratic forms (:132, :146) whose blocks the
%pairwise harmonic-mean rescaling at :141-142 can push out of positive
%semi-definiteness.
%
%NOTHING IS MODIFIED. The coefficient is the formula evaluated exactly.
adm = psyrat_admissibility([],'universe-score variance',u,'psyrat_dodiffrel',...
    'numerator');
adm = psyrat_admissibility(adm,'between-trial contrast variance',bt_scaled,...
    'psyrat_dodiffrel','denominator');
adm = psyrat_admissibility(adm,'residual contrast variance',er_scaled,...
    'psyrat_dodiffrel','denominator');
adm = psyrat_admissibility(adm,...
    'between-trial contrast variance at a single observation (ICC)',bt_unscaled,...
    'psyrat_dodiffrel','denominator');

psyrat_admissibility_warn(adm);

if strcmpi(est,'dep')
    rel_draws = u ./ dep_denom;
    icc_draws = u ./ icc_dep_denom;
    est_type = 'dependability';
else
    rel_draws = u ./ gen_denom;
    icc_draws = u ./ icc_gen_denom;
    est_type = 'generalizability';
end

qrel = quantile(rel_draws,[ciedge 1-ciedge]);
qicc = quantile(icc_draws,[ciedge 1-ciedge]);

dodiff = struct;
dodiff.est_type = est_type;
dodiff.pt = mean(rel_draws);
dodiff.ll = qrel(1);
dodiff.ul = qrel(2);
dodiff.icc_pt = mean(icc_draws);
dodiff.icc_ll = qicc(1);
dodiff.icc_ul = qicc(2);

qu = quantile(u,[ciedge 1-ciedge]);
dodiff.uni_var_pt = mean(u);
dodiff.uni_var_ll = qu(1);
dodiff.uni_var_ul = qu(2);

qbt = quantile(bt_unscaled,[ciedge 1-ciedge]);
dodiff.bt_contrast_var_pt = mean(bt_unscaled);
dodiff.bt_contrast_var_ll = qbt(1);
dodiff.bt_contrast_var_ul = qbt(2);

qbt_scaled = quantile(bt_scaled,[ciedge 1-ciedge]);
dodiff.bt_scaled_var_pt = mean(bt_scaled);
dodiff.bt_scaled_var_ll = qbt_scaled(1);
dodiff.bt_scaled_var_ul = qbt_scaled(2);

qer = quantile(er_unscaled,[ciedge 1-ciedge]);
dodiff.er_contrast_var_pt = mean(er_unscaled);
dodiff.er_contrast_var_ll = qer(1);
dodiff.er_contrast_var_ul = qer(2);

qer_scaled = quantile(er_scaled,[ciedge 1-ciedge]);
dodiff.er_scaled_var_pt = mean(er_scaled);
dodiff.er_scaled_var_ll = qer_scaled(1);
dodiff.er_scaled_var_ul = qer_scaled(2);

dodiff.rel_draws = rel_draws;
dodiff.icc_draws = icc_draws;
dodiff.obs = obs;
dodiff.cvec = cvec';

%the admissibility record travels with the coefficient it describes
%
%CONSUMERS: pass ONE of these records to psyrat_admissibility_note per
%evaluation, never both halves of a dep/gen pair. Callers here always compute
%the two together from the same draws, but every recorded quantity -- u,
%bt_scaled, er_scaled, bt_unscaled -- is built above, BEFORE the est branch, so
%the dep record and the gen record are bit-identical. The note rolls up by
%unique({label}) and SUMS nneg and ndraw, so handing it both would not check
%anything twice; it would silently double every count inside one line and tell
%the user that twice as many draws were inspected as actually were. Take either
%one and drop the other (RC-43).
%
%Records from DIFFERENT evaluations -- one per group, say -- are a different
%case and SHOULD be concatenated, which is how psyrat_admissibility_collect
%aggregates across groups.
dodiff.admissibility = adm;

end

function val = local_req(args,name,msg)
ind = find(strcmpi(name,args),1);
if isempty(ind)
    error(['varargin:' lower(name)],... %Error code and associated error
        msg);
end
val = args{ind+1};
end
