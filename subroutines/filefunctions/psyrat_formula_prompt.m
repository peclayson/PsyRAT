function result = psyrat_formula_prompt(varargin)
%Interactive runner for formula catalog entries.
%
%result = psyrat_formula_prompt('formula_id','t2_global_d')
%
%Output fields:
% result.ok
% result.command
% result.summary
% result.ll
% result.pt
% result.ul

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

formula_id = '';

if ~isempty(varargin)
    if mod(length(varargin),2)
        error('varargin:incomplete',... %Error code and associated error
            strcat('WARNING: Inputs are incomplete \n\n',...
            'Make sure each variable input is paired with a value \n'));
    end
    
    ind = find(strcmpi('formula_id',varargin),1);
    if ~isempty(ind)
        formula_id = char(string(varargin{ind+1}));
    end
end

if isempty(formula_id)
    error('varargin:formula_id',... %Error code and associated error
        'WARNING: formula_id is required.');
end

[template,titleText] = psyrat_formula_template('formula_id',formula_id);

result = struct();
result.ok = false;
result.command = template;
result.summary = '';
result.ll = [];
result.pt = [];
result.ul = [];

switch lower(formula_id)
    case {'t2_global_d','t2_global_g'}
        if strcmpi(formula_id,'t2_global_d')
            prompts = {'Between-person SD (bp)','Within-person total SD (wp)',...
                'Observations (n)','Credible interval (CI)'};
            defs = {'0.50','1.00','20','0.95'};
        else
            prompts = {'Between-person SD (bp)','Within-person total SD (wp)',...
                'Item SD (i)','Observations (n)','Credible interval (CI)'};
            defs = {'0.50','1.00','0.20','20','0.95'};
        end
        vals = psyrat_getvals(prompts,defs,titleText);
        if isempty(vals)
            return;
        end
        
        if strcmpi(formula_id,'t2_global_d')
            [ll,pt,ul] = psyrat_rel_sing('gcoeff',1,'metric','global',...
                'bp',vals(1),'wp',vals(2),'obs',vals(3),'CI',vals(4));
        else
            [ll,pt,ul] = psyrat_rel_sing('gcoeff',2,'metric','global',...
                'bp',vals(1),'wp',vals(2),'i',vals(3),'obs',vals(4),'CI',vals(5));
        end
        
    case {'t2_icc_abs','t2_icc_rel'}
        if strcmpi(formula_id,'t2_icc_abs')
            prompts = {'Between-person SD (bp)','Within-person total SD (wp)',...
                'Credible interval (CI)'};
            defs = {'0.50','1.00','0.95'};
        else
            prompts = {'Between-person SD (bp)','Within-person total SD (wp)',...
                'Item SD (i)','Credible interval (CI)'};
            defs = {'0.50','1.00','0.20','0.95'};
        end
        vals = psyrat_getvals(prompts,defs,titleText);
        if isempty(vals)
            return;
        end
        
        if strcmpi(formula_id,'t2_icc_abs')
            [ll,pt,ul] = psyrat_rel_sing('gcoeff',1,'metric','icc',...
                'bp',vals(1),'wp',vals(2),'CI',vals(3));
        else
            [ll,pt,ul] = psyrat_rel_sing('gcoeff',2,'metric','icc',...
                'bp',vals(1),'wp',vals(2),'i',vals(3),'CI',vals(4));
        end
        
    case 't2_cut_d'
        prompts = {'Between-person SD (bp)','Within-person total SD (wp)',...
            'Universe score mean (mu)','Criterion/Cut-score (cut)',...
            'Observations (n)','Credible interval (CI)'};
        defs = {'0.50','1.00','0.10','0','20','0.95'};
        vals = psyrat_getvals(prompts,defs,titleText);
        if isempty(vals)
            return;
        end

        [ll,pt,ul] = psyrat_cutscore_sing('gcoeff',1,...
            'bp',vals(1),'wp',vals(2),'mu',vals(3),'cut',vals(4),...
            'obs',vals(5),'CI',vals(6));

    case {'t3_rel_d','t3_rel_g'}
        prompts = {'Between-person SD (bp)','Occasion SD (bo)','Trial SD (bt)',...
            'Trial x Person SD (txp)','Occasion x Person SD (oxp)',...
            'Trial x Occasion SD (txo)','Residual SD (err)',...
            'Reliability type (1=CE,2=CS,3=CES)',...
            'Observations (n)','Occasions (nocc)','Credible interval (CI)'};
        defs = {'0.50','0.10','0.20','0.15','0.12','0.09','0.70','1','20','2','0.95'};
        vals = psyrat_getvals(prompts,defs,titleText);
        if isempty(vals)
            return;
        end
        
        if strcmpi(formula_id,'t3_rel_d')
            gcoeff = 1;
        else
            gcoeff = 2;
        end
        
        [ll,pt,ul] = psyrat_rel_trt(...
            'gcoeff',gcoeff,...
            'reltype',vals(8),...
            'bp',vals(1),'bo',vals(2),'bt',vals(3),...
            'txp',vals(4),'oxp',vals(5),'txo',vals(6),'err',vals(7),...
            'obs',vals(9),'nocc',vals(10),'CI',vals(11));
        
    case {'t3_icc_abs','t3_icc_rel'}
        prompts = {'Between-person SD (bp)','Occasion SD (bo)','Trial SD (bt)',...
            'Trial x Person SD (txp)','Occasion x Person SD (oxp)',...
            'Trial x Occasion SD (txo)','Residual SD (err)',...
            'Reliability type (1=CE,2=CS,3=CES)','Credible interval (CI)'};
        defs = {'0.50','0.10','0.20','0.15','0.12','0.09','0.70','1','0.95'};
        vals = psyrat_getvals(prompts,defs,titleText);
        if isempty(vals)
            return;
        end
        
        if strcmpi(formula_id,'t3_icc_abs')
            gcoeff = 1;
        else
            gcoeff = 2;
        end
        
        [ll,pt,ul] = psyrat_rel_trt(...
            'gcoeff',gcoeff,...
            'reltype',vals(8),...
            'bp',vals(1),'bo',vals(2),'bt',vals(3),...
            'txp',vals(4),'oxp',vals(5),'txo',vals(6),'err',vals(7),...
            'obs',1,'nocc',1,'CI',vals(9));
        
    case 't3_cut_d'
        prompts = {'Between-person SD (bp)','Occasion SD (bo)','Trial SD (bt)',...
            'Trial x Person SD (txp)','Occasion x Person SD (oxp)',...
            'Trial x Occasion SD (txo)','Residual SD (err)',...
            'Reliability type (1=CE,2=CS,3=CES)',...
            'Universe score mean (mu)','Criterion/Cut-score (cut)',...
            'Observations (n)','Occasions (nocc)','Credible interval (CI)'};
        defs = {'0.50','0.10','0.20','0.15','0.12','0.09','0.70','1','0.10','0','20','2','0.95'};
        vals = psyrat_getvals(prompts,defs,titleText);
        if isempty(vals)
            return;
        end

        [ll,pt,ul] = psyrat_cutscore_trt(...
            'gcoeff',1,...
            'reltype',vals(8),...
            'bp',vals(1),'bo',vals(2),'bt',vals(3),...
            'txp',vals(4),'oxp',vals(5),'txo',vals(6),'err',vals(7),...
            'mu',vals(9),'cut',vals(10),...
            'obs',vals(11),'nocc',vals(12),'CI',vals(13));
        
    case {'t6_diff1_d','t6_diff1_g'}
        prompts = {'BP variance X','BP variance Y','BP covariance XY',...
            'BT variance X','BT variance Y','BT covariance XY',...
            'Residual variance X','Residual variance Y','Residual covariance XY',...
            'Observations X','Observations Y','Credible interval (CI)'};
        defs = {'0.50','0.45','0.10','0.25','0.20','0.04','0.55','0.62','0.02','20','20','0.95'};
        vals = psyrat_getvals(prompts,defs,titleText);
        if isempty(vals)
            return;
        end
        
        bp = psyrat_covmat(vals(1),vals(2),vals(3));
        bt = psyrat_covmat(vals(4),vals(5),vals(6));
        er_var = log([sqrt(vals(7)) sqrt(vals(8))]);
        
        if strcmpi(formula_id,'t6_diff1_d')
            est = 'dep';
        else
            est = 'gen';
        end
        
        out = psyrat_diffrel(...
            'bp',bp,'bt',bt,'er_var',er_var,'wp_cov',vals(9),...
            'obs',[vals(10) vals(11)],...
            'est',est,'CI',vals(12));
        
        ll = out.ll;
        pt = out.pt;
        ul = out.ul;
        
    case {'t6_diff2_d','t6_diff2_g'}
        prompts = {'BP variance X','BP variance Y','BP covariance XY',...
            'BPI variance X','BPI variance Y','BPI covariance XY',...
            'BPO variance X','BPO variance Y','BPO covariance XY',...
            'BT variance X','BT variance Y','BT covariance XY',...
            'BO variance X','BO variance Y','BO covariance XY',...
            'BOI variance X','BOI variance Y','BOI covariance XY',...
            'Residual variance X','Residual variance Y','Residual covariance XY',...
            'Observations X','Observations Y','Occasions X','Occasions Y',...
            'Reliability type (1=CE,2=CS,3=CES)','Credible interval (CI)'};
        defs = {'0.50','0.45','0.10','0.25','0.22','0.06','0.20','0.18','0.05',...
            '0.18','0.16','0.04','0.10','0.09','0.02','0.08','0.07','0.02',...
            '0.50','0.55','0.03','20','20','2','2','3','0.95'};
        vals = psyrat_getvals(prompts,defs,titleText);
        if isempty(vals)
            return;
        end
        
        bp = psyrat_covmat(vals(1),vals(2),vals(3));
        bpi = psyrat_covmat(vals(4),vals(5),vals(6));
        bpo = psyrat_covmat(vals(7),vals(8),vals(9));
        bt = psyrat_covmat(vals(10),vals(11),vals(12));
        bo = psyrat_covmat(vals(13),vals(14),vals(15));
        boi = psyrat_covmat(vals(16),vals(17),vals(18));
        er_var = log([sqrt(vals(19)) sqrt(vals(20))]);
        
        if strcmpi(formula_id,'t6_diff2_d')
            est = 'dep';
        else
            est = 'gen';
        end
        
        out = psyrat_diffrel_trt(...
            'bp',bp,'bpi',bpi,'bpo',bpo,'bt',bt,'bo',bo,'boi',boi,...
            'er_var',er_var,'er_cov',vals(21),...
            'obs',[vals(22) vals(23)],'nocc',[vals(24) vals(25)],...
            'reltype',vals(26),'est',est,'CI',vals(27));
        
        ll = out.ll;
        pt = out.pt;
        ul = out.ul;
        
    case {'dyn_sigma_z','dyn_g','dyn_d','dyn_icc_rel','dyn_icc_abs',...
            'dyn_diff_g','dyn_diff_d','dyn_rescor'}
        %Dynamic/conditional reliability is a posterior surface over the
        %standardized dimension, not a single closed-form number, so there is no
        %interactive calculator. Throwing here surfaces the guidance through the
        %Formula Guide's Run Formula errordlg (an out.ok=false return would be
        %swallowed silently). Show Usage shows the formula; for the coefficient
        %forms it is also a runnable point-estimate template.
        error('psyrat:dynrelNoCalculator',...
            ['Dynamic/conditional reliability is a surface over the standardized ',...
            'dimension and is estimated from single-trial data via Specify Inputs ',...
            '(Bayesian location-scale fit). Use ''Show Usage'' to see the formula; for ',...
            'the coefficient forms (G(z), D(z), ICC) it is also a runnable point-estimate ',...
            'template you can copy and run.']);
    otherwise
        error('varargin:formula_id',... %Error code and associated error
            strcat('WARNING: Unknown formula_id: ',formula_id,'\n'));
end

result.ok = true;
result.ll = ll;
result.pt = pt;
result.ul = ul;
result.summary = sprintf('%s\nll = %.4f\npt = %.4f\nul = %.4f\n\n%s',...
    titleText,ll,pt,ul,result.command);

end

function vals = psyrat_getvals(prompts,defs,titleText)
raw = inputdlg(prompts,titleText,1,defs);
if isempty(raw)
    vals = [];
    return;
end

vals = zeros(1,length(raw));
for i = 1:length(raw)
    vals(i) = str2double(raw{i});
    if isnan(vals(i))
        errordlg(sprintf('Input %d is not numeric.',i),'Formula Input Error','modal');
        vals = [];
        return;
    end
end
end

function out = psyrat_covmat(v1,v2,cv)
out = zeros(1,2,2);
out(1,1,1) = v1;
out(1,2,2) = v2;
out(1,2,1) = cv;
out(1,1,2) = cv;
end
