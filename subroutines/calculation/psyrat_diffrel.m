function diffscore = psyrat_diffrel(varargin)
%Calculate difference score reliability of single session data
%
%[idtable] = psyrat_diffrel('bp',id_varcov,'bt',trl_varcov,...
%  'er_var',b_sigma,'obs',[20 40],'CI',.95,'est','dep','wp_cov',0)
%
%
%Inputs
% bp - between-person (co)variance estimates from CmdStan (id_varcov)
% bt - between-trial (co)variance estimates from CmdStan (trl_varcov)
% er_var - log residual SD draws from CmdStan (b_sigma); exponentiated and
%  squared internally
% obs - observations for measure 1 and 2 [trls trls]
% est - type of estimate to be outputted
%  'dep' - dependability, 'gen' - generalizability
% CI - size of the credible interval in decimal format: .95 = 95%
%
%Optional Inputs
% wp_cov - within-person covariance between X and Y. Pass sigma_XY(pi,e)
%  itself, NOT twice it: the factor of 2 in Table 6's error term is applied
%  internally, so the error variance carries 2*wp_cov/n-bar. This should only
%  be supplied when X and Y are concurrently measured and covariance is
%  estimated in the model. Default: 0 (non-concurrent events).
%
%Outputs
% diffscore - array, which contains
%  est_type - estimate type used ('dep' or 'gen')
%  ll - lower limit of credible interval for dependability/generalizability
%  pt - point estimate of credible interval for dependability/generalizability
%  ul - upper limit of credible interval for dependability/generalizability
%  bp_cov_ll - lower limit of credible interval for between-person covariance
%  bp_cov_pt - point estimate of credible interval for between-person covariance
%  bp_cov_ul - upper limit of credible interval for between-person covariance
%  bt1_var_ll/pt/ul - credible interval for between-trial variance of measure 1
%  bt2_var_ll/pt/ul - credible interval for between-trial variance of measure 2
%  bt_cov_ll - lower limit of credible interval for between-trial covariance
%  bt_cov_pt - point estimate of credible interval for between-trial covariance
%  bt_cov_ul - upper limit of credible interval for between-trial covariance
%  wp_cov - within-person (residual) covariance
%  wp_cov_ll/ul - lower/upper limit of credible interval for within-person covariance
%  wp_cov_est - within-person covariance estimate used
%  icc_ll - lower limit of credible interval for ICC
%  icc_pt - point estimate of credible interval for ICC
%  icc_ul - upper limit of credible interval for ICC
% 

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
            'See help psyrat_diffrel for more information about inputs'));
    end
    
    %check if bp was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('bp',varargin),1);
    if ~isempty(ind)
        bp = varargin{ind+1};
    else
        error('varargin:bp',... %Error code and associated error
            strcat('WARNING: Between-person (co)variances not specified \n\n',...
            'Please input bp. See help psyrat_diffrel for more information \n'));
    end
    
    %check if bt was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('bt',varargin),1);
    if ~isempty(ind)
        bt = varargin{ind+1};
    else
        error('varargin:bt',... %Error code and associated error
            strcat('WARNING: Between-trial (co)variances not specified \n\n',...
            'Please input bt. See help psyrat_diffrel for more information \n'));
    end
    
    %check if er_var was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('er_var',varargin),1);
    if ~isempty(ind)
        er_var = varargin{ind+1};
    else
        error('varargin:er_var',... %Error code and associated error
            strcat('WARNING: Error variance not specified \n\n',...
            'Please input er_var. See help psyrat_diffrel for more information \n'));
    end
    
    %check if wp_cov was specified.
    %If it is not found, default to no residual covariance.
    ind = find(strcmpi('wp_cov',varargin),1);
    if ~isempty(ind)
        wp_cov = varargin{ind+1};
        wp_cov_est = true;
    else
        wp_cov = [];
        wp_cov_est = false;
    end
    
    %check if est was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('obs',varargin),1);
    if ~isempty(ind)
        obs = varargin{ind+1};
    else
        error('varargin:obs',... %Error code and associated error
            strcat('WARNING: please input number of trials for each measure \n\n',...
            'Please input obs. See help psyrat_diffrel for more information \n'));
    end
    
    %check if est was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('est',varargin),1);
    if ~isempty(ind)
        est = varargin{ind+1};
    else
        error('varargin:est',... %Error code and associated error
            strcat('WARNING: reliability estimate not specified \n\n',...
            'Please input est. See help psyrat_diffrel for more information \n'));
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
                'See help psyrat_diffrel for more information \n'));
        end
    else
        error('varargin:ci',... %Error code and associated error
            strcat('WARNING: Size of credible interval not specified \n\n',...
            'Please input CI. See help psyrat_diffrel for more information \n'));
    end
    
end

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
bp1 = bp(:,1,1);
bp2 = bp(:,2,2);
bp_cov = bp(:,2,1);

bt1 = bt(:,1,1);
bt2 = bt(:,2,2);
bt_cov = bt(:,2,1);

wp1 = exp(er_var(:,1)) .^ 2;
wp2 = exp(er_var(:,2)) .^ 2;

if isempty(wp_cov)
    wp_cov = zeros(size(wp1));
else
    if iscell(wp_cov)
        wp_cov = cell2mat(wp_cov);
    end
    wp_cov = wp_cov(:);
    if length(wp_cov) == 1 && length(wp1) > 1
        wp_cov = repmat(wp_cov,length(wp1),1);
    end
    if length(wp_cov) ~= length(wp1)
        error('varargin:wp_cov',... %Error code and associated error
            strcat('WARNING: wp_cov should have the same number of rows as er_var\n',...
            'Please input wp_cov as posterior draws matching er_var.\n'));
    end
end

obs1 = obs(1);
obs2 = obs(2);

uni = bp1 + bp2 - (2*bp_cov);

rel_err = (wp1 ./ obs1) + (wp2 ./ obs2) - ((2 .* wp_cov) ./ psyrat_harmmean(obs));
abs_err = rel_err + (bt1 ./ obs1) + (bt2 ./ obs2) -...
    ((2 .* bt_cov) ./ psyrat_harmmean(obs));

%Record whether the quantities that must be variances actually are. Both error
%terms subtract a covariance scaled by the harmonic mean of the two trial
%counts, which falls below their geometric mean whenever the counts differ, so
%a legitimate draw can drive either below zero and take the coefficient outside
%[0,1]. NOTHING IS MODIFIED HERE: the coefficient below is the formula
%evaluated exactly, and the diagnostic travels with it so the user is told what
%happened rather than shown a silently repaired number. See
%psyrat_admissibility.
%
%uni is recorded too. It is a raw variance of a difference and cannot be
%negative for positive semi-definite input, but it CAN be exactly zero, which
%is a real result worth telling the user about: perfectly correlated conditions
%of equal variance leave the difference score with no between-person variance.
adm = psyrat_admissibility([],'universe-score variance',uni,'psyrat_diffrel',...
    'numerator');
adm = psyrat_admissibility(adm,'relative error variance',rel_err,'psyrat_diffrel',...
    'denominator');
adm = psyrat_admissibility(adm,'absolute error variance',abs_err,'psyrat_diffrel',...
    'denominator');

%The ICC denominators below are the same quantities at a single observation
%(n' = 1), where the harmonic mean equals the count and the covariance term
%cannot outrun the variances. For any positive semi-definite input,
%wp1 + wp2 - 2*wp_cov >= (sqrt(wp1) - sqrt(wp2))^2 >= 0, and likewise for the
%between-trial block.
%
%That is a conditional, and until RC-37 this comment read it as a guarantee and
%concluded these needed no guard. They do. The premise is exactly what fails in
%the gamma / scaled-chi-square family, whose variance diagonals are clamped at
%zero while the matching cross-covariances are not (RC-23): a block that is not
%positive semi-definite drives either of these negative at ANY count, including
%the equal counts where the projection cause cannot arise, and the out-of-range
%ICC was then printed in silence beside the coefficient. They are recorded under
%their own labels so a reader can tell an inadmissible ICC from an inadmissible
%coefficient. The ICC NUMERATOR needs no second record: it is uni, recorded
%above, unchanged.
denom_rel_err = wp1 + wp2 - (2 .* wp_cov);
denom_abs_err = denom_rel_err + bt1 + bt2 - (2 .* bt_cov);

adm = psyrat_admissibility(adm,...
    'relative error variance at a single observation (ICC)',denom_rel_err,...
    'psyrat_diffrel','denominator');
adm = psyrat_admissibility(adm,...
    'absolute error variance at a single observation (ICC)',denom_abs_err,...
    'psyrat_diffrel','denominator');

%One warning per call, raised once every quantity has been recorded so the one
%block covers the coefficient and the ICC together. Callers that invoke this
%kernel in a loop (per participant, or per point on a conditional-reliability
%surface) suppress this identifier around the loop and report the aggregate
%once instead.
psyrat_admissibility_warn(adm);

%Validate the DOMAIN of est before the switch. The switch has no otherwise branch, so a
%mis-cased or misspelled value silently assigns nothing and then fails with an opaque
%"Undefined variable ll" instead of naming the bad input. psyrat_diffrel_trt validates
%this way at its top; this mirrors it.
if ~any(strcmpi(est,{'dep','gen'}))
    error('varargin:est',... %Error code and associated error
        strcat('WARNING: reliability estimate not specified \n\n',...
        'Please input est as ''dep'' (absolute) or ''gen'' (relative)\n',...
        'See help psyrat_diffrel for more information \n'));
end

switch est
    case 'dep'
        
        temp = quantile(uni ./ (uni + abs_err),[ciedge, 1-ciedge]);
        ll = temp(1);
        ul = temp(2);
        pt = mean(uni ./ (uni + abs_err));
        
        temp = quantile(uni ./ (uni + denom_abs_err),[ciedge 1-ciedge]);
        icc_ll = temp(1);
        icc_ul = temp(2);
        icc_pt = mean(uni ./ (uni + denom_abs_err));
        
        est_type = 'dependability';
        
    case 'gen'
        
        temp = quantile(uni ./ (uni + rel_err),[ciedge 1-ciedge]);
        ll = temp(1);
        ul = temp(2);
        pt = mean(uni ./ (uni + rel_err));
        
        temp = quantile(uni ./ (uni + denom_rel_err),[ciedge 1-ciedge]);
        icc_ll = temp(1);
        icc_ul = temp(2);
        icc_pt = mean(uni ./ (uni + denom_rel_err));
        
        est_type = 'generalizability';
        
end

%Outputs
% diffscore - array, which contains
%  est_type - estimate type used ('dep' or 'gen')
%  ll - lower limit of credible interval for dependability/generalizability
%  pt - point estimate of credible interval for dependability/generalizability
%  ul - upper limit of credible interval for dependability/generalizability
%  bp_cov_ll - lower limit of credible interval for between-person covariance
%  bp_cov_pt - point estimate of credible interval for between-person covariance
%  bp_cov_ul - upper limit of credible interval for between-person covariance
%  bt1_var_ll/pt/ul - credible interval for between-trial variance of measure 1
%  bt2_var_ll/pt/ul - credible interval for between-trial variance of measure 2
%  bt_cov_ll - lower limit of credible interval for between-trial covariance
%  bt_cov_pt - point estimate of credible interval for between-trial covariance
%  bt_cov_ul - upper limit of credible interval for between-trial covariance
%  wp_cov - within-person (residual) covariance
%  wp_cov_ll/ul - lower/upper limit of credible interval for within-person covariance
%  wp_cov_est - within-person covariance estimate used
%  icc_ll - lower limit of credible interval for ICC
%  icc_pt - point estimate of credible interval for ICC
%  icc_ul - upper limit of credible interval for ICC
% 

%create output structure array
diffscore = [];

diffscore.est_type = est_type;

diffscore.pt = pt;
diffscore.ll = ll;
diffscore.ul = ul;

diffscore.bp_cov_pt = mean(bp_cov);
diffscore.bp_cov_ll = quantile(bp_cov,ciedge);
diffscore.bp_cov_ul = quantile(bp_cov,1-ciedge);

%output between-trial variances as well
diffscore.bt1_var_pt = mean(bt1);
diffscore.bt1_var_ll = quantile(bt1,ciedge);
diffscore.bt1_var_ul = quantile(bt1,1-ciedge);

diffscore.bt2_var_pt = mean(bt2);
diffscore.bt2_var_ll = quantile(bt2,ciedge);
diffscore.bt2_var_ul = quantile(bt2,1-ciedge);

diffscore.bt_cov_pt = mean(bt_cov);
diffscore.bt_cov_ll = quantile(bt_cov,ciedge);
diffscore.bt_cov_ul = quantile(bt_cov,1-ciedge);

diffscore.wp_cov = mean(wp_cov);
diffscore.wp_cov_ll = quantile(wp_cov,ciedge);
diffscore.wp_cov_ul = quantile(wp_cov,1-ciedge);
diffscore.wp_cov_est = wp_cov_est;

diffscore.icc_pt = icc_pt;
diffscore.icc_ll = icc_ll;
diffscore.icc_ul = icc_ul;

%the admissibility record travels with the coefficient it describes, so a
%caller can surface it wherever the number is shown
diffscore.admissibility = adm;

end
