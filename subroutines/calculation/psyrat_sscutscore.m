function idtable = psyrat_sscutscore(varargin)
%Calculate subject-level cut-score (criterion-referenced) dependability from
%CmdStan draws.
%
%idtable = psyrat_sscutscore('bp',gro_sds,'wp_pop',pop_sdlog,...
% 'wp_ss',ind_sdlog,'mu',ind_mu,'cut',cutscore,'idtable',idtable,...
% 'CI',.95,'i',sig_trl)
%
%Cut-scores are ABSOLUTE-error only. A criterion-referenced (cut-score)
%coefficient is a dependability/absolute decision (Phi(lambda)) by definition
%(Rocha 2026, Table 2; generalizability theory); the absolute error keeps the
%trial main effect sigma_i^2. A relative-error cut-score is not a valid
%quantity, so the only accepted decision type is est = 'dep'; est = 'gen'
%(relative) is rejected. This mirrors the absolute-only rule that
%psyrat_cutscore_sing / psyrat_cutscore_trt adopt for the group level (audit
%finding F7).
%
%Inputs
% bp - between-person standard deviation draws from CmdStan (gro_sds)
% wp_pop - population estimates of within-person variance components from
%  CmdStan on log scale (pop_sdlog). With the trial main effect estimated on
%  the mean, exp(wp_pop + wp_ss)^2 is the RELATIVE residual sigma_pi,e^2(s);
%  the absolute cut-score error adds sigma_i^2 (the i input).
% wp_ss - subject-level estimates of within-person variance components from
%  CmdStan on log scale (ind_sdlog)
% mu - subject-level universe score means (posterior draws x subjects)
% cut - cut-score
% idtable - table with id, id2, and trls
% CI - size of the credible interval in decimal format: .95 = 95%
%
%Optional Inputs
% i - trial main-effect standard deviation draws (sigma_i; default: 0). Kept
%  in the absolute error so the cut-score includes sigma_i^2/n'.
% est - decision type. Must be 'dep' (dependability/absolute); the default is
%  'dep'. est = 'gen' (relative) is invalid for a cut-score and raises an error.
%
%Outputs
% idtable - idtable concatenated with cut-score reliability estimates:
%  cutdep_pt/cutdep_ll/cutdep_ul (absolute/dependability)
%
%STATUS (audit finding F6a; resolved 2026-08-11, owner-approved): documented
%as unused. No production code calls this function: nothing under
%subroutines/ references it, and the toolbox's dynamic-dispatch sites
%(feval/str2func/eval) all resolve to other function names. It is exercised
%only by tests/TestCalculationFunctions.m and
%tests/TestCalculationAccuracyOracle.m (against
%PsyRATAccuracyOracle.sscutscore), which pin it as a regression baseline. It
%is retained rather than retired because the implementation is absolute-only
%(F7) and oracle-tested, and its natural future surface, per-subject
%criterion (cut-score) classification for the sserr family alongside the
%group-level psyrat_criterionfigures path, would be a scoped feature
%increment of its own. See SCIENTIFIC_FORMULA_AUDIT.md (F6a) and the SS-CUT
%row of SCIENTIFIC_FORMULA_TRACEABILITY_TABLE.md.

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
            'See help psyrat_sscutscore for more information about inputs'));
    end
    
    ind = find(strcmpi('bp',varargin),1);
    if ~isempty(ind)
        bp = varargin{ind+1};
    else
        error('varargin:bp',... %Error code and associated error
            strcat('WARNING: Between-person variance not specified \n\n',...
            'Please input bp. See help psyrat_sscutscore for more information \n'));
    end
    
    ind = find(strcmpi('wp_pop',varargin),1);
    if ~isempty(ind)
        wp_pop = varargin{ind+1};
    else
        error('varargin:wp_pop',... %Error code and associated error
            strcat('WARNING: Within-person population variance not specified \n\n',...
            'Please input wp_pop. See help psyrat_sscutscore for more information \n'));
    end
    
    ind = find(strcmpi('wp_ss',varargin),1);
    if ~isempty(ind)
        wp_ss = varargin{ind+1};
    else
        error('varargin:wp_ss',... %Error code and associated error
            strcat('WARNING: Subject-level within-person variance not specified \n\n',...
            'Please input wp_ss. See help psyrat_sscutscore for more information \n'));
    end
    
    ind = find(strcmpi('mu',varargin),1);
    if ~isempty(ind)
        mu = varargin{ind+1};
    else
        error('varargin:mu',... %Error code and associated error
            strcat('WARNING: Subject-level universe score means not specified \n\n',...
            'Please input mu. See help psyrat_sscutscore for more information \n'));
    end
    
    ind = find(strcmpi('cut',varargin),1);
    if ~isempty(ind)
        cut = varargin{ind+1};
    else
        error('varargin:cut',... %Error code and associated error
            strcat('WARNING: Cut-score not specified \n\n',...
            'Please input cut. See help psyrat_sscutscore for more information \n'));
    end
    
    ind = find(strcmpi('idtable',varargin),1);
    if ~isempty(ind)
        idtable = varargin{ind+1};
    else
        error('varargin:idtable',... %Error code and associated error
            strcat('WARNING: idtable not specified \n\n',...
            'Please input idtable. See help psyrat_sscutscore for more information \n'));
    end
    
    ind = find(strcmpi('CI',varargin),1);
    if ~isempty(ind)
        ciperc = varargin{ind+1};
        if ciperc > 1 || ciperc < 0
            error('varargin:ci',... %Error code and associated error
                strcat('WARNING: Size of credible interval should ',...
                'be a value between 0 and 1\n',...
                'A value of ',sprintf(' %2.2f',ciperc),...
                ' is invalid\n',...
                'See help psyrat_sscutscore for more information \n'));
        end
    else
        error('varargin:ci',... %Error code and associated error
            strcat('WARNING: Size of credible interval not specified \n\n',...
            'Please input CI. See help psyrat_sscutscore for more information \n'));
    end
    
    ind = find(strcmpi('i',varargin),1);
    if ~isempty(ind)
        i = varargin{ind+1};
    else
        i = 0;
    end
    
    ind = find(strcmpi('est',varargin),1);
    if ~isempty(ind)
        est = varargin{ind+1};
        if ~strcmpi(est,'dep')
            error('varargin:est',... %Error code and associated error
                strcat('WARNING: Cut-scores are absolute-error only.\n\n',...
                'A criterion-referenced cut-score is a dependability ',...
                '(absolute)\ndecision; a relative-error cut-score is not a ',...
                'valid quantity.\nThe only accepted value is est = ''dep''. ',...
                'est = ''gen'' (relative) is\nnot permitted. See help ',...
                'psyrat_sscutscore for more information.\n'));
        end
    else
        est = 'dep';
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

%estimate between-person and subject-level within-person variance
bp = cell2mat(bp) .^ 2;
wp = exp(wp_pop + cell2mat(wp_ss)) .^ 2;

%reshape means to posterior draws x subjects
mu = psyrat_mumat(mu,size(wp,1),size(wp,2));

i = i(:).^2;
if length(i) == 1
    i = repmat(i,size(wp,1),1);
elseif length(i) ~= size(wp,1)
    error('varargin:drawsize',... %Error code and associated error
        strcat('WARNING: Posterior draw sizes do not align across inputs.\n'));
end

offsq = (mu - cut) .^ 2;

%Absolute (dependability) cut-score only: wp is the per-subject relative
%residual sigma_pi,e^2(s) and i is the population trial main effect sigma_i^2,
%so the absolute error keeps sigma_i^2 (Rocha 2026, Table 2; audit F7). The
%est='gen' branch was removed because a relative-error cut-score is not valid.
err_term = (wp ./ idtable.trls(:)') + (i ./ idtable.trls(:)');
var_prefix = 'cutdep';

rel = (bp + offsq) ./ (bp + offsq + err_term);

rel_ll = quantile(rel,ciedge);
rel_pt = mean(rel);
rel_ul = quantile(rel,1-ciedge);

idtable.([var_prefix '_pt']) = rel_pt';
idtable.([var_prefix '_ll']) = rel_ll';
idtable.([var_prefix '_ul']) = rel_ul';

end

function mu = psyrat_mumat(mu,nrow,ncol)
%Coerce mu into a [posterior draws x subjects] matrix.

if iscell(mu)
    mu = cell2mat(mu);
end

if isvector(mu)
    if length(mu) == nrow
        mu = repmat(mu(:),1,ncol);
    elseif length(mu) == ncol
        mu = repmat(mu(:)',nrow,1);
    else
        error('varargin:mu',... %Error code and associated error
            strcat('WARNING: Shape of mu is incompatible with wp_ss draws.\n'));
    end
elseif size(mu,1) == nrow && size(mu,2) == ncol
    %already in expected shape
elseif size(mu,2) == nrow && size(mu,1) == ncol
    mu = mu';
else
    error('varargin:mu',... %Error code and associated error
        strcat('WARNING: Shape of mu is incompatible with wp_ss draws.\n'));
end

end
