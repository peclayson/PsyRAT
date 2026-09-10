function diffscore = psyrat_diffrel_trt(varargin)
%Calculate two-facet difference score reliability from GT variance terms
%
%diffscore = psyrat_diffrel_trt('bp',bp,'bpi',bpi,'bpo',bpo,...
% 'bt',bt,'bo',bo,'boi',boi,'er_var',er_var,'obs',[20 20],...
% 'nocc',[2 2],'reltype',3,'est','dep','CI',.95,'er_cov',0)
%
%
%Inputs
% bp - between-person (co)variance estimates [draw x 2 x 2]
% bpi - person x trial (co)variance estimates [draw x 2 x 2]
% bpo - person x occasion (co)variance estimates [draw x 2 x 2]
% bt - between-trial (co)variance estimates [draw x 2 x 2]
% bo - between-occasion (co)variance estimates [draw x 2 x 2]
% boi - trial x occasion (co)variance estimates [draw x 2 x 2]
% er_var - residual standard deviation estimates on log scale [draw x 2]
% obs - number of trials for measure X and Y [nX nY]
% nocc - number of occasions for measure X and Y [nX nY]
% reltype - reliability coefficient
%  1 - coefficient of equivalence (CE)
%  2 - coefficient of stability (CS)
%  3 - coefficient of equivalence and stability (CES)
% est - type of estimate to be outputted
%  'dep' - dependability (D coefficient), 'gen' - generalizability (G)
% CI - size of credible interval in decimal format: .95 = 95%
%
%Optional Inputs
% er_cov - residual covariance term sigma_XY(poi,e). This should only be
%  supplied when X and Y are concurrently measured and covariance is
%  estimated in the model. Default: 0 (non-concurrent events).
%
%Outputs
% diffscore - structure containing reliability estimates and covariance
% summaries for Table 6 two-facet formulas.

% Copyright (C) 2016-2025 Peter E. Clayson
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

if ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)
        error('varargin:incomplete',... %Error code and associated error
            strcat('WARNING: Inputs are incomplete \n\n',...
            'Make sure each variable input is paired with a value \n',...
            'See help psyrat_diffrel_trt for more information about inputs'));
    end
    
    bp = psyrat_getreq(varargin,'bp', ...
        'Between-person (co)variances not specified');
    bpi = psyrat_getreq(varargin,'bpi', ...
        'Person x trial (co)variances not specified');
    bpo = psyrat_getreq(varargin,'bpo', ...
        'Person x occasion (co)variances not specified');
    bt = psyrat_getreq(varargin,'bt', ...
        'Between-trial (co)variances not specified');
    bo = psyrat_getreq(varargin,'bo', ...
        'Between-occasion (co)variances not specified');
    boi = psyrat_getreq(varargin,'boi', ...
        'Trial x occasion (co)variances not specified');
    er_var = psyrat_getreq(varargin,'er_var', ...
        'Residual variance not specified');
    obs = psyrat_getreq(varargin,'obs', ...
        'Please input number of trials for each measure');
    nocc = psyrat_getreq(varargin,'nocc', ...
        'Please input number of occasions for each measure');
    reltype = psyrat_getreq(varargin,'reltype', ...
        'Please input reltype');
    est = psyrat_getreq(varargin,'est', ...
        'Reliability estimate not specified');
    
    ind = find(strcmpi('er_cov',varargin),1);
    if ~isempty(ind)
        er_cov = varargin{ind+1};
        er_cov_est = true;
    else
        er_cov = 0;
        er_cov_est = false;
    end
    
    %check if CI was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('CI',varargin),1);
    if ~isempty(ind)
        ciperc = varargin{ind+1};
        if ciperc > 1 || ciperc < 0
            error('varargin:ci',... %Error code and associated error
                strcat('WARNING: Size of credible interval should ',...
                'be a value between 0 and 1\n',...
                'A value of ',sprintf(' %2.2f',ciperc),...
                ' is invalid\n',...
                'See help psyrat_diffrel_trt for more information \n'));
        end
    else
        error('varargin:ci',... %Error code and associated error
            strcat('WARNING: Size of credible interval not specified \n\n',...
            'Please input CI. See help psyrat_diffrel_trt for more information \n'));
    end
    
end

if ~any(reltype == [1 2 3])
    error('varargin:reltype',... %Error code and associated error
        strcat('WARNING: reltype is invalid\n',...
        'Valid values are 1 (CE), 2 (CS), or 3 (CES)\n'));
end

if ~any(strcmpi(est,{'dep','gen'}))
    error('varargin:est',... %Error code and associated error
        strcat('WARNING: reliability estimate not specified \n\n',...
        'Please input est as ''dep'' or ''gen''\n',...
        'See help psyrat_diffrel_trt for more information \n'));
end

obs = psyrat_paircount(obs,'obs');
nocc = psyrat_paircount(nocc,'nocc');

%make sure the user understood the the CI input is width not edges, if
%inputted incorrectly, provide a warning, but don't change
if ciperc < .5
    str = sprintf(' %2.0f%%',100*ciperc);
    warning('ci:width',... %Warning code and associated warning
        strcat('WARNING: Size of credible interval is small \n\n',...
        'User specified a credible interval width of ',str,'%\n',...
        'If this was intended, ignore this warning.\n'));
end

%calculate edges for CI
ciedge = (1-ciperc)/2;

%parse (co)variance components
[bp1,bp2,bp_cov] = psyrat_parsecov(bp,'bp');
[bpi1,bpi2,bpi_cov] = psyrat_parsecov(bpi,'bpi');
[bpo1,bpo2,bpo_cov] = psyrat_parsecov(bpo,'bpo');
[bt1,bt2,bt_cov] = psyrat_parsecov(bt,'bt');
[bo1,bo2,bo_cov] = psyrat_parsecov(bo,'bo');
[boi1,boi2,boi_cov] = psyrat_parsecov(boi,'boi');

if size(er_var,2) ~= 2
    error('varargin:er_var',... %Error code and associated error
        strcat('WARNING: er_var should be [posterior draws x 2]\n',...
        'for the X and Y residual terms.\n'));
end

wp1 = exp(er_var(:,1)) .^ 2;
wp2 = exp(er_var(:,2)) .^ 2;

%align posterior draw sizes across inputs
ndraw = max([length(bp1),length(bpi1),length(bpo1),length(bt1),...
    length(bo1),length(boi1),length(wp1),length(er_cov)]);

bp1 = psyrat_expanddraw(bp1,ndraw,'bp');
bp2 = psyrat_expanddraw(bp2,ndraw,'bp');
bp_cov = psyrat_expanddraw(bp_cov,ndraw,'bp');

bpi1 = psyrat_expanddraw(bpi1,ndraw,'bpi');
bpi2 = psyrat_expanddraw(bpi2,ndraw,'bpi');
bpi_cov = psyrat_expanddraw(bpi_cov,ndraw,'bpi');

bpo1 = psyrat_expanddraw(bpo1,ndraw,'bpo');
bpo2 = psyrat_expanddraw(bpo2,ndraw,'bpo');
bpo_cov = psyrat_expanddraw(bpo_cov,ndraw,'bpo');

bt1 = psyrat_expanddraw(bt1,ndraw,'bt');
bt2 = psyrat_expanddraw(bt2,ndraw,'bt');
bt_cov = psyrat_expanddraw(bt_cov,ndraw,'bt');

bo1 = psyrat_expanddraw(bo1,ndraw,'bo');
bo2 = psyrat_expanddraw(bo2,ndraw,'bo');
bo_cov = psyrat_expanddraw(bo_cov,ndraw,'bo');

boi1 = psyrat_expanddraw(boi1,ndraw,'boi');
boi2 = psyrat_expanddraw(boi2,ndraw,'boi');
boi_cov = psyrat_expanddraw(boi_cov,ndraw,'boi');

wp1 = psyrat_expanddraw(wp1,ndraw,'er_var');
wp2 = psyrat_expanddraw(wp2,ndraw,'er_var');
er_cov = psyrat_expanddraw(er_cov,ndraw,'er_cov');

obs1 = obs(1);
obs2 = obs(2);
nocc1 = nocc(1);
nocc2 = nocc(2);

hmi = psyrat_harmmean(obs);
hmo = psyrat_harmmean(nocc);

%difference-score component variances/covariances
bp_raw = bp1 + bp2 - (2 .* bp_cov);
bpi_raw = bpi1 + bpi2 - (2 .* bpi_cov);
bpo_raw = bpo1 + bpo2 - (2 .* bpo_cov);
bt_raw = bt1 + bt2 - (2 .* bt_cov);
bo_raw = bo1 + bo2 - (2 .* bo_cov);
boi_raw = boi1 + boi2 - (2 .* boi_cov);
bpoi_raw = wp1 + wp2 - (2 .* er_cov);

%D-study adjusted error components
bpi = (bpi1 ./ obs1) + (bpi2 ./ obs2) - ((2 .* bpi_cov) ./ hmi);
bpo = (bpo1 ./ nocc1) + (bpo2 ./ nocc2) - ((2 .* bpo_cov) ./ hmo);
bt = (bt1 ./ obs1) + (bt2 ./ obs2) - ((2 .* bt_cov) ./ hmi);
bo = (bo1 ./ nocc1) + (bo2 ./ nocc2) - ((2 .* bo_cov) ./ hmo);
boi = (boi1 ./ (obs1*nocc1)) + (boi2 ./ (obs2*nocc2)) -...
    ((2 .* boi_cov) ./ (hmi*hmo));
bpoi = (wp1 ./ (obs1*nocc1)) + (wp2 ./ (obs2*nocc2)) -...
    ((2 .* er_cov) ./ (hmi*hmo));

if reltype == 1 %coefficient of equivalence
    uni = bp_raw + bpo;
    rel_err = bpi + bpoi;
    abs_err = rel_err + bt + boi;
    
    uni_icc = bp_raw + bpo_raw;
    rel_err_icc = bpi_raw + bpoi_raw;
    abs_err_icc = rel_err_icc + bt_raw + boi_raw;
    
elseif reltype == 2 %coefficient of stability
    uni = bp_raw + bpi;
    rel_err = bpo + bpoi;
    abs_err = rel_err + bo + boi;
    
    uni_icc = bp_raw + bpi_raw;
    rel_err_icc = bpo_raw + bpoi_raw;
    abs_err_icc = rel_err_icc + bo_raw + boi_raw;
    
elseif reltype == 3 %coefficient of trial equivalence and stability
    uni = bp_raw;
    rel_err = bpi + bpo + bpoi;
    abs_err = rel_err + bt + bo + boi;
    
    uni_icc = bp_raw;
    rel_err_icc = bpi_raw + bpo_raw + bpoi_raw;
    abs_err_icc = rel_err_icc + bt_raw + bo_raw + boi_raw;
end

%Record whether the quantities that must be variances actually are. NOTHING IS
%MODIFIED: the coefficient below is the formula evaluated exactly, and the
%diagnostic travels with it. See psyrat_admissibility.
%
%The UNIVERSE term is recorded, not just the error terms, and for this kernel
%that is the point. Every D-study adjusted component above subtracts a
%covariance scaled by a harmonic mean (of the trial counts, the occasion
%counts, or their product), and two of the three reltype branches put one of
%those blocks in the NUMERATOR: uni = bp_raw + bpo for the coefficient of
%equivalence and uni = bp_raw + bpi for the coefficient of stability. Rocha
%et al. (2026) Table 6's note prescribes exactly that addition, so the formula
%is right, but it means the universe-score variance can itself go negative.
%
%TWO causes, not one, and reltype 3 is NOT exempt from the second. See
%psyrat_admissibility's help for the full statement; in brief: (1) at UNEQUAL
%counts the harmonic mean falls below the geometric mean, so the subtraction can
%exceed the variances it is taken from -- the PROJECTION leaving the space of
%variances, which reaches the two branches named above; (2) a raw cross-condition
%block that is not positive semi-definite, which the gamma/scaled-chi-square
%family can produce (RC-23), and which reaches EVERY branch at ANY counts.
%reltype 3's uni = bp_raw = bp1 + bp2 - 2*bp_cov (:205) is negative whenever
%bp_cov > sqrt(bp1*bp2), so calling it structurally safe -- as this comment did
%until RC-35 -- was false. The counterexample is already in the repo:
%TestTrtDiffCurveWiring sets a between-person block at r > 1 on an EQUAL-count
%curve and drives exactly this quantity negative.
%
%The _icc quantities are variances of a difference at a single observation, so
%cause (1) cannot touch them -- there is no projection to leave -- and they are
%non-negative for positive semi-definite input. Cause (2) DOES reach them, at
%any counts, so since RC-37 they are recorded too, under labels of their own so
%that an inadmissible ICC can be told apart from an inadmissible coefficient.
%
%The two error terms are recorded on all three branches. The universe term is
%recorded on reltype 1 and 2 ONLY, and the omission at reltype 3 is deliberate
%rather than an oversight to be tidied up later: there uni (:242) and uni_icc
%(:246) are assigned from the same variable, bp_raw, so a record here would be
%bit-identical to the one taken immediately below rather than merely equal in
%value. Adding it would count the same draws twice under a second label and
%report a single failure as two.
%
%So an ordinary ICC printed beside an inadmissible coefficient is evidence for
%the projection, not proof of it; an ICC that is itself out of range points at
%the raw components instead.
adm = psyrat_admissibility([],'universe-score variance',uni,'psyrat_diffrel_trt',...
    'numerator');
adm = psyrat_admissibility(adm,'relative error variance',rel_err,'psyrat_diffrel_trt',...
    'denominator');
adm = psyrat_admissibility(adm,'absolute error variance',abs_err,'psyrat_diffrel_trt',...
    'denominator');

if reltype ~= 3
    adm = psyrat_admissibility(adm,...
        'universe-score variance at a single observation (ICC)',uni_icc,...
        'psyrat_diffrel_trt','numerator');
end

adm = psyrat_admissibility(adm,...
    'relative error variance at a single observation (ICC)',rel_err_icc,...
    'psyrat_diffrel_trt','denominator');
adm = psyrat_admissibility(adm,...
    'absolute error variance at a single observation (ICC)',abs_err_icc,...
    'psyrat_diffrel_trt','denominator');

%One warning per call; loop callers suppress this identifier and aggregate.
psyrat_admissibility_warn(adm);

switch lower(est)
    case 'dep'
        rel = uni ./ (uni + abs_err);
        icc = uni_icc ./ (uni_icc + abs_err_icc);
        est_type = 'dependability';
    case 'gen'
        rel = uni ./ (uni + rel_err);
        icc = uni_icc ./ (uni_icc + rel_err_icc);
        est_type = 'generalizability';
end

temp = quantile(rel,[ciedge, 1-ciedge]);
ll = temp(1);
ul = temp(2);
pt = mean(rel);

temp = quantile(icc,[ciedge, 1-ciedge]);
icc_ll = temp(1);
icc_ul = temp(2);
icc_pt = mean(icc);

%create output structure array
diffscore = [];

diffscore.est_type = est_type;
diffscore.reltype = reltype;

diffscore.pt = pt;
diffscore.ll = ll;
diffscore.ul = ul;

diffscore.bp_cov_pt = mean(bp_cov);
diffscore.bp_cov_ll = quantile(bp_cov,ciedge);
diffscore.bp_cov_ul = quantile(bp_cov,1-ciedge);

diffscore.bpi_cov_pt = mean(bpi_cov);
diffscore.bpi_cov_ll = quantile(bpi_cov,ciedge);
diffscore.bpi_cov_ul = quantile(bpi_cov,1-ciedge);

diffscore.bpo_cov_pt = mean(bpo_cov);
diffscore.bpo_cov_ll = quantile(bpo_cov,ciedge);
diffscore.bpo_cov_ul = quantile(bpo_cov,1-ciedge);

diffscore.bt_cov_pt = mean(bt_cov);
diffscore.bt_cov_ll = quantile(bt_cov,ciedge);
diffscore.bt_cov_ul = quantile(bt_cov,1-ciedge);

diffscore.bo_cov_pt = mean(bo_cov);
diffscore.bo_cov_ll = quantile(bo_cov,ciedge);
diffscore.bo_cov_ul = quantile(bo_cov,1-ciedge);

diffscore.boi_cov_pt = mean(boi_cov);
diffscore.boi_cov_ll = quantile(boi_cov,ciedge);
diffscore.boi_cov_ul = quantile(boi_cov,1-ciedge);

diffscore.wp_cov = mean(er_cov);
diffscore.wp_cov_ll = quantile(er_cov,ciedge);
diffscore.wp_cov_ul = quantile(er_cov,1-ciedge);
diffscore.wp_cov_est = er_cov_est;

diffscore.icc_pt = icc_pt;
diffscore.icc_ll = icc_ll;
diffscore.icc_ul = icc_ul;

%the admissibility record travels with the coefficient it describes
diffscore.admissibility = adm;

end

function value = psyrat_getreq(args,name,msg)
ind = find(strcmpi(name,args),1);
if ~isempty(ind)
    value = args{ind+1};
else
    error(['varargin:' lower(name)],... %Error code and associated error
        strcat('WARNING: ',msg,' \n\n',...
        'Please input ',name,'. See help psyrat_diffrel_trt for more information \n'));
end
end

function counts = psyrat_paircount(inp,name)
if length(inp) == 1
    counts = [inp inp];
elseif length(inp) == 2
    counts = inp(:)';
else
    error(['varargin:' lower(name)],... %Error code and associated error
        strcat('WARNING: ',name,' should be a scalar or 2-element vector.\n'));
end

if any(counts <= 0)
    error(['varargin:' lower(name)],... %Error code and associated error
        strcat('WARNING: ',name,' values should be positive.\n'));
end
end

function [v1,v2,cv] = psyrat_parsecov(inp,name)
if iscell(inp)
    inp = cell2mat(inp);
end

if ismatrix(inp) && all(size(inp) == [2 2])
    inp = reshape(inp,[1 2 2]);
end

if ndims(inp) ~= 3 || size(inp,2) ~= 2 || size(inp,3) ~= 2
    error(['varargin:' lower(name)],... %Error code and associated error
        strcat('WARNING: ',name,' should be [posterior draws x 2 x 2].\n'));
end

v1 = inp(:,1,1);
v2 = inp(:,2,2);
cv = inp(:,2,1);

v1 = v1(:);
v2 = v2(:);
cv = cv(:);
end

function out = psyrat_expanddraw(v,ndraw,name)
v = v(:);
if length(v) == 1
    out = repmat(v,ndraw,1);
elseif length(v) == ndraw
    out = v;
else
    error(['varargin:' lower(name)],... %Error code and associated error
        strcat('WARNING: Posterior draw sizes do not align for ',name,'.\n'));
end
end
