function idtable = psyrat_ssrel_diff_trt(varargin)
%Calculate subject-level TWO-FACET (trial + occasion) difference-score reliability
%summaries (Rast & Clayson dynamic/conditional reliability, Stage 4b / analysis
%16). This is the two-facet analogue of psyrat_ssrel_diff (the one-facet
%subject-level difference calc): it reuses the validated two-facet difference
%coefficient psyrat_diffrel_trt with the six crossed (co)variance blocks, but
%substitutes each participant's OWN per-event residual so the gain-minus-loss
%difference reliability phi_delta,s is reported per subject.
%
%The universe-score and facet (co)variances are SHARED across participants
%(group-level, as in the brms location-scale parameterization); only the per-event
%residual is subject-specific, taken from the z-conditioned per-subject log-residual
%draws (er_var_ss1/er_var_ss2) passed in by psyrat_dynrel_summary. The occasion n'
%is a single value applied to every subject (matching the population reference
%surface); each subject's observed per-event occasion counts are carried in the
%idtable (occs1/occs2) for reference only.
%
%SEM is computed in-calc (local_diff_trt_errvar) by replicating the
%reltype-specific relative/absolute error partition of psyrat_diffrel_trt, because
%psyrat_diffrel_trt returns phi and ICC but not SEM. This mirrors how the one-facet
%psyrat_ssrel_diff computes SEM manually and keeps psyrat_diffrel_trt unchanged.
%
%Required name/value inputs
% bp,bpi,bpo,bt,bo,boi - the six [draws x 2 x 2] event (co)variance blocks:
%   person (id->bp), person x trial (tid->bpi), person x occasion (oid->bpo),
%   trial (trl->bt), occasion (occ->bo), trial x occasion (to->boi).
% er_var   - [draws x 2] population per-event log-residual SD; the fallback when a
%            participant has no matching er_var_ss column.
% idtable  - per-subject table aligned to id2 (the er_var_ss column index). Must
%            include trls1/trls2 (per-event trial n'); occs1/occs2 are optional and
%            reported only.
% reltype  - 1 (CE), 2 (CS), or 3 (CES); see psyrat_diffrel_trt.
% nocc     - occasion n' for the coefficient (scalar; applied to every subject).
% CI       - credible-interval width in (0,1): .95 = 95%.
%
%Optional name/value inputs
% est          - 'dep' (default, dependability) or 'gen' (generalizability).
% er_var_ss1   - [draws x NSUB] per-subject z-conditioned event-1 log-residuals.
% er_var_ss2   - [draws x NSUB] per-subject z-conditioned event-2 log-residuals.
% er_cov       - residual covariance (scalar or [draws x 1]); default 0
%                (non-concurrent events).
% er_cov_ss    - [draws x NSUB] per-subject residual covariance; default none.
%
%Output
% idtable - the input table with per-row reliability columns appended: rel_pt/ll/ul
%           (phi for the chosen est), icc_pt/ll/ul, sem_pt/ll/ul, plus dep_* aliases
%           (= rel_*) for compatibility with the existing per-subject table/plot.

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
    error('varargin:incomplete','Inputs are incomplete');
end

req = {'bp','bpi','bpo','bt','bo','boi','er_var','idtable','reltype','nocc','CI'};
vals = struct;
for i = 1:length(req)
    %use the shared required-argument getter (psyrat_req) rather than a
    %hand-rolled find/error, matching every sibling calc file.
    vals.(req{i}) = psyrat_req(varargin,req{i},...
        sprintf('Please input %s.',req{i}));
end

idtable = vals.idtable;
reltype = vals.reltype;
nocc = vals.nocc;

est = psyrat_opt(varargin,'est','dep');
er_var_ss1 = psyrat_opt(varargin,'er_var_ss1',[]);
er_var_ss2 = psyrat_opt(varargin,'er_var_ss2',[]);
er_cov = psyrat_opt(varargin,'er_cov',0);
er_cov_ss = psyrat_opt(varargin,'er_cov_ss',[]);

if ~any(reltype == [1 2 3])
    error('varargin:reltype',...
        'reltype is invalid. Valid values are 1 (CE), 2 (CS), or 3 (CES).');
end
if ~any(strcmpi(est,{'dep','gen'}))
    error('varargin:est','est must be ''dep'' or ''gen''.');
end

%per-event trial n' contract (mirror psyrat_ssrel_diff): trls1/trls2 are required;
%fall back to a single trls column if only it is present.
if ~any(strcmp(idtable.Properties.VariableNames,'trls1'))
    if any(strcmp(idtable.Properties.VariableNames,'trls'))
        idtable.trls1 = idtable.trls;
    else
        error('varargin:idtable','idtable must include trls1 (or trls).');
    end
end
if ~any(strcmp(idtable.Properties.VariableNames,'trls2'))
    if any(strcmp(idtable.Properties.VariableNames,'trls'))
        idtable.trls2 = idtable.trls;
    else
        error('varargin:idtable','idtable must include trls2 (or trls).');
    end
end

%subject-level residual draws are usable only when both per-event draw matrices and
%the id2 index column are present (matches psyrat_ssrel_diff).
has_ss = ~isempty(er_var_ss1) && ~isempty(er_var_ss2) && ...
    any(strcmp(idtable.Properties.VariableNames,'id2'));

ciedge = (1-vals.CI)/2;

idtable.rel_pt = zeros(height(idtable),1);
idtable.rel_ll = zeros(height(idtable),1);
idtable.rel_ul = zeros(height(idtable),1);
idtable.icc_pt = zeros(height(idtable),1);
idtable.icc_ll = zeros(height(idtable),1);
idtable.icc_ul = zeros(height(idtable),1);
idtable.sem_pt = zeros(height(idtable),1);
idtable.sem_ll = zeros(height(idtable),1);
idtable.sem_ul = zeros(height(idtable),1);

%One warning per table, not one per participant -- see psyrat_ssrel_diff. This
%is the two-facet twin and the single worst-exposed path: per-participant trial
%counts are essentially never equal, and for the coefficient of stability the
%UNIVERSE term carries the harmonic-mean-scaled block, not just the error.
adm_ws = warning('off','psyrat:negerrvar');
adm_all = [];

for r = 1:height(idtable)
    %per-event trial n' for this subject (occasion n' is the shared scalar nocc)
    obs = [max(1,idtable.trls1(r)) max(1,idtable.trls2(r))];

    %default residual = the population per-event log-residual; replace with this
    %subject's z-conditioned draws when available (indexed by id2 = Stan column).
    er_draw = vals.er_var;
    ercov_draw = er_cov;
    if has_ss
        sid = idtable.id2(r);
        if isnumeric(sid) && isfinite(sid) && sid >= 1 && ...
                sid <= size(er_var_ss1,2) && sid <= size(er_var_ss2,2)
            sid = round(sid);
            er_draw = [er_var_ss1(:,sid) er_var_ss2(:,sid)];
            %Concurrent subject-level variants supply a per-subject residual
            %covariance in er_cov_ss, already on this subject's z-conditioned SD
            %scale. The non-concurrent path passes er_cov = 0 with an empty
            %er_cov_ss, so ercov_draw stays 0 here -- consistent with the
            %subject-specific SDs. A nonzero SCALAR er_cov with an empty
            %er_cov_ss (a population-scale covariance against subject-scale SDs)
            %would be inconsistent, but no caller produces that combination.
            if ~isempty(er_cov_ss) && sid <= size(er_cov_ss,2)
                ercov_draw = er_cov_ss(:,sid);
            end
        end
    end

    %phi + ICC from the validated two-facet difference coefficient
    diffscore = psyrat_diffrel_trt(...
        'bp',vals.bp,'bpi',vals.bpi,'bpo',vals.bpo,...
        'bt',vals.bt,'bo',vals.bo,'boi',vals.boi,...
        'er_var',er_draw,'obs',obs,'nocc',[nocc nocc],...
        'reltype',reltype,'est',est,'CI',vals.CI,'er_cov',ercov_draw);

    %SEM: relative-error SD (gen) or absolute-error SD (dep), per draw, summarized.
    %The person/universe block (vals.bp) is intentionally NOT passed: error variance
    %excludes the universe-score component (it enters only the phi numerator, which
    %comes from psyrat_diffrel_trt above, not from this SEM helper).
    [rel_err,abs_err] = local_diff_trt_errvar(vals.bpi,vals.bpo,...
        vals.bt,vals.bo,vals.boi,er_draw,obs,[nocc nocc],reltype,ercov_draw);
    %NaN rather than a floor at zero, for the reason given in psyrat_ssrel_diff:
    %a negative error variance has no real square root, and flooring it would
    %print a confident small SEM for a quantity that is not a variance.
    if strcmpi(est,'gen')
        errvar_draw = rel_err;
    else
        errvar_draw = abs_err;
    end

    sem_draw = nan(size(errvar_draw));
    sem_draw(errvar_draw >= 0) = sqrt(errvar_draw(errvar_draw >= 0));

    idtable.rel_pt(r) = diffscore.pt;
    idtable.rel_ll(r) = diffscore.ll;
    idtable.rel_ul(r) = diffscore.ul;
    idtable.icc_pt(r) = diffscore.icc_pt;
    idtable.icc_ll(r) = diffscore.icc_ll;
    idtable.icc_ul(r) = diffscore.icc_ul;
    idtable.sem_pt(r) = mean(sem_draw);
    idtable.sem_ll(r) = quantile(sem_draw,ciedge);
    idtable.sem_ul(r) = quantile(sem_draw,1-ciedge);

    if isfield(diffscore,'admissibility')
        adm_all = [adm_all, diffscore.admissibility]; %#ok<AGROW> one per row
    end
end

warning(adm_ws);
psyrat_admissibility_warn(adm_all);

%aliases for compatibility with the existing per-subject plotting/table functions
idtable.dep_pt = idtable.rel_pt;
idtable.dep_ll = idtable.rel_ll;
idtable.dep_ul = idtable.rel_ul;

end

function [rel_err,abs_err] = local_diff_trt_errvar(bpi,bpo,bt,bo,boi,...
    er_var,obs,nocc,reltype,er_cov)
%Per-draw relative- and absolute-error variances for the two-facet difference
%score, replicating the error partition in psyrat_diffrel_trt (which returns
%phi/ICC but not SEM). The person/universe block is deliberately absent: it enters
%the universe-score numerator, not the error variances. bpi..boi are [draws x 2 x 2]
%event blocks; er_var is [draws x 2] log-SD; obs/nocc are 2-element trial/occasion
%counts; er_cov is the residual covariance.
[bpi1,bpi2,bpi_cov] = local_parsecov(bpi);
[bpo1,bpo2,bpo_cov] = local_parsecov(bpo);
[bt1,bt2,bt_cov]    = local_parsecov(bt);
[bo1,bo2,bo_cov]    = local_parsecov(bo);
[boi1,boi2,boi_cov] = local_parsecov(boi);

wp1 = exp(er_var(:,1)) .^ 2;
wp2 = exp(er_var(:,2)) .^ 2;

ndraw = max([numel(bpi1),numel(bpo1),numel(bt1),numel(bo1),numel(boi1),...
    numel(wp1),numel(er_cov)]);
bpi1 = local_expand(bpi1,ndraw); bpi2 = local_expand(bpi2,ndraw); bpi_cov = local_expand(bpi_cov,ndraw);
bpo1 = local_expand(bpo1,ndraw); bpo2 = local_expand(bpo2,ndraw); bpo_cov = local_expand(bpo_cov,ndraw);
bt1  = local_expand(bt1,ndraw);  bt2  = local_expand(bt2,ndraw);  bt_cov  = local_expand(bt_cov,ndraw);
bo1  = local_expand(bo1,ndraw);  bo2  = local_expand(bo2,ndraw);  bo_cov  = local_expand(bo_cov,ndraw);
boi1 = local_expand(boi1,ndraw); boi2 = local_expand(boi2,ndraw); boi_cov = local_expand(boi_cov,ndraw);
wp1  = local_expand(wp1,ndraw);  wp2  = local_expand(wp2,ndraw);
er_cov = local_expand(er_cov,ndraw);

obs1 = obs(1); obs2 = obs(2);
nocc1 = nocc(1); nocc2 = nocc(2);
hmi = psyrat_harmmean(obs);
hmo = psyrat_harmmean(nocc);

%D-study adjusted difference-score error components (psyrat_diffrel_trt:214-221)
bpi = (bpi1 ./ obs1) + (bpi2 ./ obs2) - ((2 .* bpi_cov) ./ hmi);
bpo = (bpo1 ./ nocc1) + (bpo2 ./ nocc2) - ((2 .* bpo_cov) ./ hmo);
bt  = (bt1 ./ obs1) + (bt2 ./ obs2) - ((2 .* bt_cov) ./ hmi);
bo  = (bo1 ./ nocc1) + (bo2 ./ nocc2) - ((2 .* bo_cov) ./ hmo);
boi = (boi1 ./ (obs1*nocc1)) + (boi2 ./ (obs2*nocc2)) - ((2 .* boi_cov) ./ (hmi*hmo));
bpoi = (wp1 ./ (obs1*nocc1)) + (wp2 ./ (obs2*nocc2)) - ((2 .* er_cov) ./ (hmi*hmo));

switch reltype
    case 1 %coefficient of equivalence (occasion fixed)
        rel_err = bpi + bpoi;
        abs_err = rel_err + bt + boi;
    case 2 %coefficient of stability (trial fixed)
        rel_err = bpo + bpoi;
        abs_err = rel_err + bo + boi;
    case 3 %coefficient of equivalence and stability (both random)
        rel_err = bpi + bpo + bpoi;
        abs_err = rel_err + bt + bo + boi;
end
end

function [v1,v2,cv] = local_parsecov(inp)
%Parse a [draws x 2 x 2] (or single 2x2) event block into event-1 variance,
%event-2 variance, and the cross-event covariance draw vectors.
if iscell(inp)
    inp = cell2mat(inp);
end
if ismatrix(inp) && all(size(inp) == [2 2])
    inp = reshape(inp,[1 2 2]);
end
v1 = inp(:,1,1); v1 = v1(:);
v2 = inp(:,2,2); v2 = v2(:);
cv = inp(:,2,1); cv = cv(:);
end

function out = local_expand(v,ndraw)
v = v(:);
if isscalar(v)
    out = repmat(v,ndraw,1);
else
    out = v;
end
end
