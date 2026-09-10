function idtable = psyrat_ssrel_diff(varargin)
%Calculate subject-level difference-score reliability summaries.

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

if mod(length(varargin),2)
    error('varargin:incomplete','Inputs are incomplete');
end

req = {'bp','bt','er_var','idtable','CI'};
vals = struct;
for i = 1:length(req)
    ind = find(strcmpi(req{i},varargin),1);
    if isempty(ind)
        error(['varargin:' req{i}],['Missing input: ' req{i}]);
    end
    vals.(req{i}) = varargin{ind+1};
end

idtable = vals.idtable;

ind = find(strcmpi('wp_cov',varargin),1);
if isempty(ind)
    wp_cov = 0;
else
    wp_cov = varargin{ind+1};
end

ind = find(strcmpi('er_var_ss1',varargin),1);
if isempty(ind)
    er_var_ss1 = [];
else
    er_var_ss1 = varargin{ind+1};
end

ind = find(strcmpi('er_var_ss2',varargin),1);
if isempty(ind)
    er_var_ss2 = [];
else
    er_var_ss2 = varargin{ind+1};
end

ind = find(strcmpi('wp_cov_ss',varargin),1);
if isempty(ind)
    wp_cov_ss = [];
else
    wp_cov_ss = varargin{ind+1};
end

ind = find(strcmpi('est',varargin),1);
if isempty(ind)
    est = 'dep';
else
    est = varargin{ind+1};
end

if ~any(strcmp(idtable.Properties.VariableNames,'trls'))
    error('varargin:idtable','idtable must include trls');
end

if ~any(strcmp(idtable.Properties.VariableNames,'trls1'))
    idtable.trls1 = idtable.trls;
end
if ~any(strcmp(idtable.Properties.VariableNames,'trls2'))
    idtable.trls2 = idtable.trls;
end

idtable.rel_pt = zeros(height(idtable),1);
idtable.rel_ll = zeros(height(idtable),1);
idtable.rel_ul = zeros(height(idtable),1);
idtable.icc_pt = zeros(height(idtable),1);
idtable.icc_ll = zeros(height(idtable),1);
idtable.icc_ul = zeros(height(idtable),1);
idtable.sem_pt = zeros(height(idtable),1);
idtable.sem_ll = zeros(height(idtable),1);
idtable.sem_ul = zeros(height(idtable),1);
idtable.bp_var = zeros(height(idtable),1);
idtable.ss_errvar = zeros(height(idtable),1);
idtable.pop_errvar = zeros(height(idtable),1);

bp1 = vals.bp(:,1,1);
bp2 = vals.bp(:,2,2);
bp_cov = vals.bp(:,2,1);
bt1 = vals.bt(:,1,1);
bt2 = vals.bt(:,2,2);
bt_cov = vals.bt(:,2,1);
wp1 = exp(vals.er_var(:,1)).^2;
wp2 = exp(vals.er_var(:,2)).^2;
has_ss = ~isempty(er_var_ss1) && ~isempty(er_var_ss2) && ...
    ~isempty(wp_cov_ss) && any(strcmp(idtable.Properties.VariableNames,'id2'));

ciedge = (1-vals.CI)/2;
pop_wp = mean((wp1 + wp2)./2);
if has_ss
    pop_wp = mean((exp(er_var_ss1(:)).^2 + exp(er_var_ss2(:)).^2)./2);
end

%The kernel warns once per call, and this loop calls it once per participant,
%so a real table would emit one block per row. Suppress the identifier here,
%accumulate the per-row records, and report the aggregate once after the loop.
%Per-participant trial counts are essentially never equal between two
%conditions, so this path is the most exposed one in the toolbox.
adm_ws = warning('off','psyrat:negerrvar');
adm_all = [];

for r = 1:height(idtable)
    obs = [max(1,idtable.trls1(r)) max(1,idtable.trls2(r))];
    er_draw = vals.er_var;
    wp_cov_draw = wp_cov;
    wp1_draw = wp1;
    wp2_draw = wp2;
    
    if has_ss
        sid = idtable.id2(r);
        if isnumeric(sid) && isfinite(sid) && sid >= 1 && ...
                sid <= size(er_var_ss1,2) && sid <= size(er_var_ss2,2) && ...
                sid <= size(wp_cov_ss,2)
            sid = round(sid);
            er_draw = [er_var_ss1(:,sid) er_var_ss2(:,sid)];
            wp_cov_draw = wp_cov_ss(:,sid);
            wp1_draw = exp(er_draw(:,1)).^2;
            wp2_draw = exp(er_draw(:,2)).^2;
        end
    else
        if isscalar(wp_cov)
            wp_cov_draw = repmat(wp_cov,size(vals.bp,1),1);
        else
            wp_cov_draw = wp_cov(:);
        end
    end
    
    diffscore = psyrat_diffrel(...
        'bp',vals.bp,...
        'bt',vals.bt,...
        'er_var',er_draw,...
        'wp_cov',wp_cov_draw,...
        'obs',obs,...
        'est',est,...
        'CI',vals.CI);

    rel_err = (wp1_draw ./ obs(1)) + (wp2_draw ./ obs(2));
    rel_err = rel_err - ((2 .* wp_cov_draw) ./ psyrat_harmmean(obs));
    abs_err = rel_err + (bt1 ./ obs(1)) + (bt2 ./ obs(2)) -...
        ((2 .* bt_cov) ./ psyrat_harmmean(obs));
    %The SEM is the square root of an error variance, so a negative variance
    %has no real root and the SEM is genuinely undefined for that draw. It is
    %reported as NaN rather than floored at zero: max(...,0) would print a
    %confident, small SEM for a quantity that is not a variance at all, and
    %would hide the very thing the user needs to see. The error variance itself
    %is recorded below and reaches the user through the difference score's
    %admissibility note.
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
    idtable.bp_var(r) = mean(bp1 + bp2 - (2*bp_cov));
    idtable.ss_errvar(r) = mean(rel_err);
    idtable.pop_errvar(r) = pop_wp;

    if isfield(diffscore,'admissibility')
        adm_all = [adm_all, diffscore.admissibility]; %#ok<AGROW> one per row
    end
end

warning(adm_ws);
psyrat_admissibility_warn(adm_all);

%aliases for compatibility with existing plotting/table functions
idtable.dep_pt = idtable.rel_pt;
idtable.dep_ll = idtable.rel_ll;
idtable.dep_ul = idtable.rel_ul;

end
