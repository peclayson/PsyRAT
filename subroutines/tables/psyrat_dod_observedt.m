function [observedtable, adm] = psyrat_dod_observedt(varargin)
%Observed DoD reliability table using average recorded trials
%
%observedtable = psyrat_dod_observedt('psyrat_data',psyrat_data,...
%  'groupidx',1,'map',[1 2 3 4],'gui',1)
%
%Required Inputs:
% psyrat_data - PsyRAT results structure with rel.analysis='ic_dodiff'
% groupidx - selected group index (1 if groups='none')
% map - 1x4 event-index mapping in ERP order [ERP_1 ERP_2 ERP_3 ERP_4]
%
%Optional Inputs:
% gui - 1 show GUI (default), 0 return table only
% CI - credible interval width (default .95)
%
%Output:
% observedtable - table with observed DoD reliability for selected contrast
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
        ['Inputs are incomplete. Provide name/value pairs, for example ',...
        '''psyrat_dod_observedt(''psyrat_data'',psyrat_data,''groupidx'',1,', ...
        '''map'',[1 2 3 4],''gui'',1)''.']);
end

psyrat_data = local_req(varargin,'psyrat_data',...
    'psyrat_data was not provided.');
groupidx = local_req(varargin,'groupidx',...
    'groupidx was not provided.');
map = local_req(varargin,'map',...
    'map was not provided.');
gui = local_opt(varargin,'gui',1);
ciperc = local_opt(varargin,'CI',.95);

if ~isfield(psyrat_data,'rel') || ~isfield(psyrat_data.rel,'analysis') || ...
        ~strcmpi(psyrat_data.rel.analysis,'ic_dodiff')
    error('varargin:analysis',... %Error code and associated error
        'Observed DoD table requires rel.analysis = ''ic_dodiff''.');
end

if ~isnumeric(groupidx) || ~isscalar(groupidx) || isnan(groupidx) || ...
        groupidx < 1 || abs(round(groupidx)-groupidx) > 0
    error('varargin:groupidx',... %Error code and associated error
        'groupidx must be a positive integer scalar.');
end
groupidx = round(groupidx);

if ~isnumeric(map) || numel(map) ~= 4
    error('varargin:map',... %Error code and associated error
        'map must be a numeric vector with four event indices [ERP_1 ERP_2 ERP_3 ERP_4].');
end
map = round(map(:)');
if numel(unique(map)) ~= 4
    error('varargin:map',... %Error code and associated error
        'map must include four unique event indices.');
end

rel = psyrat_data.rel;
[groupnames,hasgroup] = local_groupnames(rel);
eventnames = cellstr(string(rel.events(:)));

if groupidx > numel(groupnames)
    error('varargin:groupidx',... %Error code and associated error
        'groupidx exceeds available groups.');
end
if any(map < 1) || any(map > numel(eventnames))
    error('varargin:map',... %Error code and associated error
        'map indices exceed the number of available events.');
end

[bp,bt,erlog] = local_groupdraws(rel,groupidx);
bp = bp(:,map,map);
bt = bt(:,map,map);
erlog = erlog(:,map);

[obs,ninc,nexc] = local_observed_trials(rel.data,hasgroup,groupnames{groupidx},eventnames(map));

dep = local_nan_summary();
gen = local_nan_summary();
if all(obs > 0)
    dep = psyrat_dodiffrel(...
        'bp',bp,...
        'bt',bt,...
        'er_var',erlog,...
        'obs',obs,...
        'est','dep',...
        'CI',ciperc);
    gen = psyrat_dodiffrel(...
        'bp',bp,...
        'bt',bt,...
        'er_var',erlog,...
        'obs',obs,...
        'est','gen',...
        'CI',ciperc);
end

contrast_label = sprintf('(%s - %s) - (%s - %s)',...
    eventnames{map(1)},eventnames{map(2)},eventnames{map(3)},eventnames{map(4)});

observedtable = table(...
    {contrast_label},...
    ninc,...
    nexc,...
    {local_ci_str(dep.pt,dep.ll,dep.ul)},...
    {local_ci_str(gen.pt,gen.ll,gen.ul)},...
    {local_ci_str(dep.icc_pt,dep.icc_ll,dep.icc_ul)},...
    {local_ci_str(gen.icc_pt,gen.icc_ll,gen.icc_ul)},...
    obs(1),obs(2),obs(3),obs(4),...
    'VariableNames',{'Label','n_Included','n_Excluded',...
    'Dependability','Generalizability','ICC_Dep','ICC_Gen',...
    'MeanTrials_ERP_1','MeanTrials_ERP_2','MeanTrials_ERP_3','MeanTrials_ERP_4'});

meta = struct;
meta.group = groupnames{groupidx};
meta.event_labels = eventnames(map);
meta.contrast = contrast_label;
meta.obs = obs;

%Hand the admissibility record back for the same reason psyrat_dod_varcompt
%does: psyrat_report sweeps this per group and aggregates one note (RC-43). dep
%only, because the dep and gen records are bit-identical. Empty when the zero-
%trial path substituted local_nan_summary, which carries no record.
adm = [];
if isstruct(dep) && isfield(dep,'admissibility')
    adm = dep.admissibility;
end

if gui == 1
    local_show_gui(psyrat_data,observedtable,meta);

    adm_note = psyrat_dod_admissibility_note(dep);
    if ~isempty(adm_note)
        warndlg(adm_note,'Coefficient is not interpretable as a reliability');
    end
end

end

function local_show_gui(psyrat_data,observedtable,meta)
figwidth = 1320;
figheight = 320;
rowspace = 25;
row = figheight - rowspace*2;

fig = psyrat_newfigure('unit','pix',...
    'position',[1120 560 figwidth figheight],...
    'menub','no',...
    'Tag','psyrat_output',...
    'name','Observed Difference-of-Differences Reliability',...
    'numbertitle','off',...
    'resize','off');

uicontrol(fig,'Style','text','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Observed Difference-of-Differences Reliability (Mean Recorded Trials)',...
    'Position',[0 row figwidth 24]);

uicontrol(fig,'Style','text','fontsize',12,...
    'HorizontalAlignment','left',...
    'String',['Group: ' meta.group],...
    'Position',[20 row-28 420 22]);

uicontrol(fig,'Style','text','fontsize',12,...
    'HorizontalAlignment','left',...
    'String',['Contrast: ' meta.contrast],...
    'Position',[20 row-52 figwidth-40 22]);

t = uitable('Parent',fig,'Position',[20 82 figwidth-40 130],...
    'Data',table2cell(observedtable));
set(t,'ColumnName',observedtable.Properties.VariableNames);
set(t,'ColumnWidth',{320 80 80 140 140 130 130 95 95 95 95});
set(t,'RowName',[]);

uicontrol(fig,'Style','push','fontsize',13,...
    'HorizontalAlignment','center',...
    'String','Save Table',...
    'Position',[figwidth/8 18 figwidth/4 44],...
    'Callback',{@local_save_table,psyrat_data,observedtable,meta});
end

function local_save_table(varargin)
psyrat_data = varargin{3};
observedtable = varargin{4};
meta = varargin{5};

%XLSX now writes cross-platform via writecell/writetable, so both formats are
%offered on every platform (including macOS).
[savename, savepath] = uiputfile(...
    {'*.xlsx','Excel File (.xlsx)';'*.csv','Comma-Separated Value File (.csv)'},...
    'Where would you like to save table?');

if isequal(savename,0) || isequal(savepath,0)
    return;
end

outfile = fullfile(savepath,savename);
[~,~,ext] = fileparts(outfile);

head = {};
head{end+1} = 'Table Generated on';
head{end+1} = datestr(clock);
head{end+1} = '';
head{end+1} = sprintf('PsyRAT Toolbox v%s',psyrat_data.ver);
head{end+1} = '';
head{end+1} = sprintf('Dataset: %s',psyrat_data.rel.filename);
head{end+1} = psyrat_iterline(psyrat_data.rel);
head{end+1} = sprintf('Group: %s',meta.group);
head{end+1} = sprintf('Contrast: %s',meta.contrast);
head{end+1} = sprintf('ERP_1: %s',meta.event_labels{1});
head{end+1} = sprintf('ERP_2: %s',meta.event_labels{2});
head{end+1} = sprintf('ERP_3: %s',meta.event_labels{3});
head{end+1} = sprintf('ERP_4: %s',meta.event_labels{4});
head{end+1} = sprintf('Observed trials used: [%0.3f %0.3f %0.3f %0.3f]',...
    meta.obs(1),meta.obs(2),meta.obs(3),meta.obs(4));
%run provenance (seed/engine/priors/convergence); fail-open
head = [head, psyrat_provenance_lines(psyrat_data)'];
head{end+1} = '';
head{end+1} = '';

if strcmpi(ext,'.xlsx')
    writecell(head',outfile);
    writetable(observedtable,outfile,'Range',strcat('A',num2str(length(head))));
elseif strcmpi(ext,'.csv')
    fid = fopen(outfile,'w');
    %quote per the CSV convention (B19): a comma-bearing header line must
    %stay ONE spreadsheet field; comma-free lines stay byte-identical
    for i = 1:numel(head)
        fprintf(fid,'%s\n',psyrat_quote_csv_header(head{i}));
    end

    vars = observedtable.Properties.VariableNames;
    fprintf(fid,'%s\n',strjoin(vars,','));
    row = table2cell(observedtable);
    out = cell(size(row));
    for c = 1:numel(row)
        out{c} = local_csv_scalar(row{c});
    end
    fprintf(fid,'%s\n',strjoin(out,','));
    fclose(fid);
end
end

function s = local_ci_str(pt,ll,ul)
if any(~isfinite([pt ll ul]))
    s = 'NaN';
else
    s = sprintf('%0.3f CI [%0.3f %0.3f]',pt,ll,ul);
end
end

function out = local_nan_summary
out = struct('pt',NaN,'ll',NaN,'ul',NaN,...
    'icc_pt',NaN,'icc_ll',NaN,'icc_ul',NaN);
end

function [obs,ninc,nexc] = local_observed_trials(datatable,hasgroup,gname,event_labels)
if hasgroup
    gmask = string(datatable.group) == string(gname);
    groupsub = datatable(gmask,:);
else
    groupsub = datatable;
end

obs = zeros(1,numel(event_labels));
event_ids = cell(1,numel(event_labels));
for e = 1:numel(event_labels)
    emask = string(groupsub.event) == string(event_labels{e});
    subset = groupsub(emask,:);
    if isempty(subset)
        obs(e) = 0;
        event_ids{e} = string.empty(0,1);
    else
        trltable = varfun(@length,subset,'GroupingVariables',{'id'});
        obs(e) = mean(trltable.GroupCount);
        event_ids{e} = unique(string(subset.id));
    end
end

allids = unique(string(groupsub.id));
if isempty(event_ids) || any(cellfun(@isempty,event_ids))
    include_ids = string.empty(0,1);
else
    include_ids = event_ids{1};
    for e = 2:numel(event_ids)
        include_ids = intersect(include_ids,event_ids{e});
    end
end

ninc = numel(include_ids);
nexc = numel(setdiff(allids,include_ids));
end

function [groupnames,hasgroup] = local_groupnames(rel)
if ischar(rel.groups) && strcmpi(rel.groups,'none')
    groupnames = {'none'};
    hasgroup = false;
else
    groupnames = cellstr(string(rel.groups(:)));
    hasgroup = true;
end
end

function [bp,bt,erlog] = local_groupdraws(rel,gidx)
bp = rel.out.id_varcov{1,gidx};
bt = rel.out.trl_varcov{1,gidx};
erlog = rel.out.b_sigma{1,gidx};
if iscell(bp)
    bp = cell2mat(bp);
end
if iscell(bt)
    bt = cell2mat(bt);
end
if iscell(erlog)
    erlog = cell2mat(erlog);
end
end

function val = local_req(args,name,msg)
ind = find(strcmpi(name,args),1);
if isempty(ind)
    error(['varargin:' lower(name)],... %Error code and associated error
        msg);
end
val = args{ind+1};
end

function val = local_opt(args,name,default_val)
ind = find(strcmpi(name,args),1);
if isempty(ind)
    val = default_val;
else
    val = args{ind+1};
end
end

function out = local_csv_scalar(v)
if ischar(v)
    out = v;
elseif isstring(v)
    out = char(v);
elseif isnumeric(v) || islogical(v)
    if isscalar(v)
        out = num2str(v,'%0.10g');
    else
        out = mat2str(v);
    end
else
    out = char(string(v));
end

out = strrep(out,'"','""');
if contains(out,',') || contains(out,'"') || contains(out,newline)
    out = ['"' out '"'];
end
end
