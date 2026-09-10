function idtable = psyrat_ssrel_trt(varargin)
%Calculate subject-level TRT reliability summaries from posterior draws.
%
%ESTIMAND NOTE. Only the residual is genuinely per-subject here. The callers
%(psyrat_relsummary's local_ssrel_trt_call, psyrat_dynrel_summary) broadcast the
%POPULATION person x trial and person x occasion SDs across subjects, because the
%location-scale model estimates a per-person residual scale, not per-person
%interaction components. So the spread of the coefficient, ICC and SEM columns
%across participants reflects differing residual scales (and differing trial
%counts), NOT person-specific occasion or trial effects. Do not read it as the
%latter.

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

req = {'gcoeff','reltype','bp','bo','bt','txp','oxp','txo','err','idtable','CI'};
vals = struct;
for i = 1:length(req)
    ind = find(strcmpi(req{i},varargin),1);
    if isempty(ind)
        error(['varargin:' req{i}],['Missing input: ' req{i}]);
    end
    vals.(req{i}) = varargin{ind+1};
end

%Validate the DOMAIN of the two selectors, not just their presence. The SEM branch below
%is an if/elseif chain whose final else silently absorbs any unmatched combination into
%the CES relative expression, so an out-of-range gcoeff or reltype would return a
%plausible number for a coefficient the caller did not ask for.
if ~any(vals.gcoeff == [1 2])
    error('varargin:gcoeff',... %Error code and associated error
        strcat('WARNING: gcoeff is invalid. Valid values are 1 (absolute) ',...
        'or 2 (relative)\n'));
end
if ~any(vals.reltype == [1 2 3])
    error('varargin:reltype',... %Error code and associated error
        strcat('WARNING: reltype is invalid. Valid values are 1 (equivalence), ',...
        '2 (stability), or 3 (equivalence and stability)\n'));
end

idtable = vals.idtable;

ind = find(strcmpi('nocc',varargin),1);
if isempty(ind)
    nocc = 1;
else
    nocc = varargin{ind+1};
end

ind = find(strcmpi('txp_ss',varargin),1);
if isempty(ind)
    txp_ss = [];
else
    txp_ss = varargin{ind+1};
end

ind = find(strcmpi('oxp_ss',varargin),1);
if isempty(ind)
    oxp_ss = [];
else
    oxp_ss = varargin{ind+1};
end

ind = find(strcmpi('err_ss',varargin),1);
if isempty(ind)
    err_ss = [];
else
    err_ss = varargin{ind+1};
end

if ~any(strcmp(idtable.Properties.VariableNames,'trls'))
    error('varargin:idtable','idtable must include trls');
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

bpv = vals.bp.^2;
bov = vals.bo.^2;
btv = vals.bt.^2;
txov = vals.txo.^2;
pop_errvar = mean(vals.err.^2);
has_ss = ~isempty(txp_ss) && ~isempty(oxp_ss) && ~isempty(err_ss) && ...
    any(strcmp(idtable.Properties.VariableNames,'id2'));

if has_ss
    pop_errvar = mean(err_ss(:).^2);
end

ciedge = (1-vals.CI)/2;

%Group-level reference quantities for the subject-level caterpillar plots
%(psyrat_ssrelplot). That plot receives only this table, never the variance
%components, so the group reference line has to be carried on the table. Both use
%the POPULATION draws (vals.err, not the pooled per-subject residual), so the line
%is the same quantity the group-level tables report for this run.
%
% gro_icc - single-observation coefficient (one trial, one occasion), the group
%  twin of the per-subject icc_* columns computed in the loop below.
% gro_rel - the selected coefficient at the mean trial count. mean(trls) is the
%  divisor psyrat_ssrelplot already used for its reference line, so carrying it
%  here changes only the variance PARTITION on that plot, not its n' semantics.
gro_obs = max(1,mean(idtable.trls));

[~,gro_icc,~] = psyrat_rel_trt(...
    'gcoeff',vals.gcoeff,...
    'reltype',vals.reltype,...
    'bp',vals.bp,...
    'bo',vals.bo,...
    'bt',vals.bt,...
    'txp',vals.txp,...
    'oxp',vals.oxp,...
    'txo',vals.txo,...
    'err',vals.err,...
    'obs',1,...
    'nocc',1,...
    'CI',vals.CI);

[~,gro_rel,~] = psyrat_rel_trt(...
    'gcoeff',vals.gcoeff,...
    'reltype',vals.reltype,...
    'bp',vals.bp,...
    'bo',vals.bo,...
    'bt',vals.bt,...
    'txp',vals.txp,...
    'oxp',vals.oxp,...
    'txo',vals.txo,...
    'err',vals.err,...
    'obs',gro_obs,...
    'nocc',nocc,...
    'CI',vals.CI);

idtable.gro_icc = repmat(gro_icc,height(idtable),1);
idtable.gro_rel = repmat(gro_rel,height(idtable),1);

for r = 1:height(idtable)
    obs = max(1,idtable.trls(r));
    txp_draw = vals.txp;
    oxp_draw = vals.oxp;
    err_draw = vals.err;
    
    if has_ss
        sid = idtable.id2(r);
        if isnumeric(sid) && isfinite(sid) && sid >= 1 && ...
                sid <= size(txp_ss,2) && sid <= size(oxp_ss,2) && ...
                sid <= size(err_ss,2)
            sid = round(sid);
            txp_draw = txp_ss(:,sid);
            oxp_draw = oxp_ss(:,sid);
            err_draw = err_ss(:,sid);
        end
    end
    
    [rll,rpt,rul] = psyrat_rel_trt(...
        'gcoeff',vals.gcoeff,...
        'reltype',vals.reltype,...
        'bp',vals.bp,...
        'bo',vals.bo,...
        'bt',vals.bt,...
        'txp',txp_draw,...
        'oxp',oxp_draw,...
        'txo',vals.txo,...
        'err',err_draw,...
        'obs',obs,...
        'nocc',nocc,...
        'CI',vals.CI);
    
    %ICC = the single-observation coefficient (one trial at one occasion), so it
    %uses the SAME two-facet partition as the coefficient above, evaluated at
    %obs = 1. Routing it through psyrat_icc instead would pass only bp and the
    %residual, discarding bo/bt/txp/oxp/txo and asserting they are zero -- which
    %returned one identical value for all six (gcoeff x reltype) combinations.
    %psyrat_rel_trt at obs = nocc = 1 reproduces Rocha (2026) Table 3's printed
    %ICC rows exactly, for every one of those six cells; that equivalence is
    %pinned against the paper by TestIndependentOracleAgreement/testTwoFacetIccs.
    %
    %nocc is PINNED to 1 rather than the run's nocc, for the composite-invariance
    %reason written out at psyrat_relsummary.m:5082-5086: an ICC must not scale
    %with a multi-occasion composite, or the same condition's ICC differs between
    %a trt run and a trt_diff run. Same construction as the group-level dynamic
    %two-facet ICC at psyrat_rel_dynrel_trt.m:204-210. Do not "fix" it to nocc.
    [icc_ll,icc_pt,icc_ul] = psyrat_rel_trt(...
        'gcoeff',vals.gcoeff,...
        'reltype',vals.reltype,...
        'bp',vals.bp,...
        'bo',vals.bo,...
        'bt',vals.bt,...
        'txp',txp_draw,...
        'oxp',oxp_draw,...
        'txo',vals.txo,...
        'err',err_draw,...
        'obs',1,...
        'nocc',1,...
        'CI',vals.CI);
    
    txpv = txp_draw.^2;
    oxpv = oxp_draw.^2;
    errv = err_draw.^2;

    if vals.gcoeff == 1 && vals.reltype == 1
        err_term = (txpv ./ obs) + (errv ./ (obs*nocc)) + (btv ./ obs) + (txov ./ (obs*nocc));
    elseif vals.gcoeff == 1 && vals.reltype == 2
        err_term = (oxpv ./ nocc) + (errv ./ (obs*nocc)) + (bov ./ nocc) + (txov ./ (obs*nocc));
    elseif vals.gcoeff == 1 && vals.reltype == 3
        err_term = (txpv ./ obs) + (oxpv ./ nocc) + (errv ./ (obs*nocc)) + (btv ./ obs) + (bov ./ nocc) + (txov ./ (obs*nocc));
    elseif vals.gcoeff == 2 && vals.reltype == 1
        err_term = (txpv ./ obs) + (errv ./ (obs*nocc));
    elseif vals.gcoeff == 2 && vals.reltype == 2
        err_term = (oxpv ./ nocc) + (errv ./ (obs*nocc));
    else
        err_term = (txpv ./ obs) + (oxpv ./ nocc) + (errv ./ (obs*nocc));
    end
    sem_draw = sqrt(err_term);

    idtable.rel_pt(r) = rpt;
    idtable.rel_ll(r) = rll;
    idtable.rel_ul(r) = rul;
    idtable.icc_pt(r) = icc_pt;
    idtable.icc_ll(r) = icc_ll;
    idtable.icc_ul(r) = icc_ul;
    idtable.sem_pt(r) = mean(sem_draw);
    idtable.sem_ll(r) = quantile(sem_draw,ciedge);
    idtable.sem_ul(r) = quantile(sem_draw,1-ciedge);
    idtable.bp_var(r) = mean(bpv);
    idtable.ss_errvar(r) = mean(errv);
    idtable.pop_errvar(r) = pop_errvar;
end

%aliases for compatibility with existing single-session subject-level code
idtable.dep_pt = idtable.rel_pt;
idtable.dep_ll = idtable.rel_ll;
idtable.dep_ul = idtable.rel_ul;

end
