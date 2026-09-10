function idtable = psyrat_ssrel(varargin)
%Calculate subject-level dependability and ICCs using variance components  
% from CmdStan
%
%[idtable] = psyrat_ssrel('bp',gro_sds,'wp_pop',pop_sdlog,'wp_ss',ind_sdlog,...
%  'idtable',idtable,'CI',.95,'i',sig_trl,'gcoeff',1)
%
%
%Inputs
% bp - between-person standard deviation draws from CmdStan (gro_sds);
%  squared internally
% wp_pop - population log residual SD draws from CmdStan (pop_sdlog);
%  exponentiated and squared internally
% wp_ss - subject-level log residual SD offset draws from CmdStan
%  (ind_sdlog); exponentiated and squared internally
% idtable - table with information about original id, new id2 for CmdStan,
%  and the number of trials retained for each id (each row id should line
%  up exactly with each column for wp_ss)
% CI - size of the credible interval in decimal format: .95 = 95%
%
%Optional Inputs
% i - trial main-effect standard deviation draws (sigma_i, the crossed trial
%  facet from the sserr model's sig_trl). Default: 0. With trial-on-mean
%  estimated, the per-subject residual wp = exp(wp_pop + wp_ss)^2 is the
%  RELATIVE residual sigma_pi,e^2(s); the absolute coefficient adds sigma_i^2.
% gcoeff - decision type: 1 = dependability/absolute (default; error keeps
%  sigma_i^2), 2 = generalizability/relative (error excludes sigma_i^2). When
%  i = 0 the two coincide (the pre-trial-effect behavior).
%
%ARGUMENT CONTRACT -- READ BEFORE COMPARING THIS FUNCTION WITH psyrat_rel_sing
%(RC-19). The two take the SAME argument names with OPPOSITE meanings:
%
%  psyrat_ssrel     wp = the RELATIVE residual sigma_pi,e^2(s); i^2 is ADDED
%                        for the absolute branch
%  psyrat_rel_sing  wp = the TOTAL within-person SD; i^2 is SUBTRACTED to
%                        recover the relative residual
%
%They reach the same two error variances from opposite directions, so a caller
%that hands one function the other's wp gets a silently plausible number rather
%than an error. This collision is the shared root cause of the sigma_i^2
%findings (RC-02, RC-03, RC-15); psyrat_cutscore_sing and psyrat_cutscore_dep
%carry the same one. The names are deliberately NOT being changed -- the call
%sites are numerous and the risk of a partial rename exceeds the ambiguity --
%so the contract is stated here and in psyrat_rel_sing's header instead.
%
%Outputs
% tableout - a table that is concatenated with idtable input.
%   id - original id
%   id2 - CmdStan id
%   trls - number of retained trials
%   dep_pt - selected-coefficient (phi_s when gcoeff=1, G_s when gcoeff=2)
%     point estimate
%   dep_ll - lower limit of the credible interval specified by CI
%   dep_ul - upper limit of the credible interval specified by CI
%   icc_pt - ICC point estimate (absolute or relative per gcoeff)
%   icc_ll - lower limit of the credible interval specified by CI for
%     ICC
%   icc_ul - upper limit of the credible interval specified by CI for
%     ICC
%   sem_pt - subject-level standard error of measurement point estimate
%   sem_ll - lower limit of the credible interval specified by CI for SEM
%   sem_ul - upper limit of the credible interval specified by CI for SEM
%   bp_var - population between-person variance estimate (same for all ids)
%   ss_errvar - subject-level RELATIVE residual variance estimates
%     (sigma_pi,e^2(s))
%   trl_var - population trial main-effect variance estimate sigma_i^2 (same
%     for all ids; 0 when i was not supplied)
%   pop_errvar - population-level residual variance estimate (same for all
%     ids)
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
        'See help psyrat_ssrel for more information about inputs'));
    end
    
    %check if bp was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('bp',varargin),1);
    if ~isempty(ind)
        bp = varargin{ind+1}; 
    else 
        error('varargin:bp',... %Error code and associated error
        strcat('WARNING: Between-person variance not specified \n\n',... 
        'Please input bp. See help psyrat_ssrel for more information \n'));
    end
    
    %check if wp_pop was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('wp_pop',varargin),1);
    if ~isempty(ind)
        wp_pop = varargin{ind+1}; 
    else 
        error('varargin:wp_pop',... %Error code and associated error
        strcat('WARNING: Within-person popluation variance not specified \n\n',... 
        'Please input wp_pop. See help psyrat_ssrel for more information \n'));
    end
    
    %check if wp_ss was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('wp_ss',varargin),1);
    if ~isempty(ind)
        wp_ss = varargin{ind+1}; 
    else 
        error('varargin:wp_ss',... %Error code and associated error
        strcat('WARNING: Subject-level within-person variance not specified \n\n',... 
        'Please input wp_ss. See help psyrat_ssrel for more information \n'));
    end
    
    %check if idtable was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('idtable',varargin),1);
    if ~isempty(ind)
        idtable = varargin{ind+1}; 
    else 
        error('varargin:idtable',... %Error code and associated error
        strcat('WARNING: idtable not specified \n\n',... 
        'Please input idtable. See help psyrat_ssrel for more information \n'));
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
                'See help psyrat_ssrel for more information \n'));
        end
    else
        error('varargin:ci',... %Error code and associated error
            strcat('WARNING: Size of credible interval not specified \n\n',...
            'Please input CI. See help psyrat_ssrel for more information \n'));
    end

    %check if i (trial main-effect SD draws, sigma_i) was specified.
    %If not, default to 0 so absolute and relative coincide (pre-trial-effect
    %behavior).
    ind = find(strcmpi('i',varargin),1);
    if ~isempty(ind)
        i = varargin{ind+1};
    else
        i = 0;
    end

    %check if gcoeff (decision type) was specified.
    %If not, default to 1 (dependability/absolute), preserving prior behavior.
    ind = find(strcmpi('gcoeff',varargin),1);
    if ~isempty(ind)
        gcoeff = varargin{ind+1};
        if ~any(gcoeff == [1 2])
            error('varargin:gcoeff',... %Error code and associated error
                strcat('WARNING: gcoeff is invalid. Valid values are 1 ',...
                '(absolute/dependability) or 2 (relative/generalizability)\n'));
        end
    else
        gcoeff = 1;
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

%estimate between-person variance and the subject-level RELATIVE residual
%variance. With the trial main effect now estimated on the mean of the sserr
%model, wp = exp(wp_pop + wp_ss)^2 is the relative residual sigma_pi,e^2(s)
%(the trial mean has been removed), per subject (posterior draws x subjects).
bp = cell2mat(bp) .^ 2;
wp = exp(wp_pop + cell2mat(wp_ss)) .^ 2;

%trial main-effect variance sigma_i^2 (population/crossed; shared across
%subjects). A column vector of draws (or scalar 0 when not supplied) that
%broadcasts across the subject columns of wp.
i2 = i(:) .^ 2;
if ~isscalar(i2) && size(i2,1) ~= size(wp,1)
    error('varargin:i',... %Error code and associated error
        strcat('WARNING: Posterior draw sizes do not align between i ',...
        '(sigma_i) and the within-person draws.\n'));
end

%absolute-error metrics keep sigma_i^2; relative-error metrics exclude it. The
%SELECTION rule is the same as the group-level psyrat_rel_sing gcoeff switch,
%but the ARITHMETIC is its mirror image, because the two functions are handed
%wp on different scales (RC-19):
%
%  psyrat_rel_sing  wp = TOTAL within-person SD, so it SUBTRACTS  (pierr = wp^2 - i^2)
%  psyrat_ssrel     wp = RELATIVE residual,      so it ADDS       (errwp = wp   + i^2)
%
%Both arrive at the same two error variances; only the starting point differs.
%An earlier version of this comment claimed the switch worked "exactly as"
%psyrat_rel_sing's, which inverted the operation and misled a verification pass.
%Do not "align" the two by copying one branch into the other -- reversing either
%one silently double-counts or double-removes sigma_i^2.
if gcoeff == 1
    errwp = wp + i2; %absolute: sigma_pi,e^2(s) + sigma_i^2
else
    errwp = wp;      %relative: sigma_pi,e^2(s) only
end

%estimate subject-level reliability (dependability phi_s or generalizability
%G_s, per gcoeff)
dep_ll = quantile(bp ./ (bp + (errwp ./ idtable.trls(:)')),ciedge);
dep_pt = mean(bp ./ (bp + (errwp ./ idtable.trls(:)')));
dep_ul = quantile(bp ./ (bp + (errwp ./ idtable.trls(:)')),1-ciedge);

%estimate subject-level ICCs (absolute or relative per gcoeff)
icc_ll = quantile(bp ./ (bp + errwp),ciedge);
icc_pt = mean(bp ./ (bp + errwp));
icc_ul = quantile(bp ./ (bp + errwp),1-ciedge);

%estimate subject-level SEMs (standard error of measurement for averages)
sem = sqrt(errwp ./ idtable.trls(:)');
sem_ll = quantile(sem,ciedge);
sem_pt = mean(sem);
sem_ul = quantile(sem,1-ciedge);

%put information into output table
idtable.dep_pt = dep_pt';
idtable.dep_ll = dep_ll';
idtable.dep_ul = dep_ul';

idtable.icc_pt = icc_pt';
idtable.icc_ll = icc_ll';
idtable.icc_ul = icc_ul';

idtable.sem_pt = sem_pt';
idtable.sem_ll = sem_ll';
idtable.sem_ul = sem_ul';

idtable.bp_var = repmat(mean(bp), height(idtable), 1);
%ss_errvar is the per-subject RELATIVE residual variance sigma_pi,e^2(s); the
%trial main effect sigma_i^2 is reported separately in trl_var (shared across
%subjects).
idtable.ss_errvar = mean(wp,1)';
idtable.trl_var = repmat(mean(i2), height(idtable), 1);
idtable.pop_errvar = repmat(mean(exp(wp_pop).^2),...
    height(idtable), 1);


end
