function [varcomptable, adm] = psyrat_dod_varcompt(varargin)
%Export summary variance components for a DoD contrast
%
%varcomptable = psyrat_dod_varcompt('psyrat_data',psyrat_data,...
%  'groupidx',1,'map',[1 2 3 4])
%
%Required Inputs:
% psyrat_data - PsyRAT results structure with rel.analysis='ic_dodiff'
% groupidx - selected group index
% map - [ERP_1 ERP_2 ERP_3 ERP_4] event indices
%
%Optional Inputs:
% obs - [n1 n2 n3 n4] trial counts (default: observed group means)
% CI - credible interval width (default .95)
% savefile - 1 save to file (default), 0 do not prompt for save
% gui - 1 show GUI table, 0 no GUI (default)
%
%Output:
% varcomptable - table with DoD variance components and reliability summaries
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
        '''psyrat_dod_varcompt(''psyrat_data'',psyrat_data,''groupidx'',1,', ...
        '''map'',[1 2 3 4])''.']);
end

psyrat_data = local_req(varargin,'psyrat_data',...
    'psyrat_data was not provided.');
groupidx = local_req(varargin,'groupidx',...
    'groupidx was not provided.');
map = local_req(varargin,'map',...
    'map was not provided.');
obs = local_opt(varargin,'obs',[]);
ciperc = local_opt(varargin,'CI',.95);
savefile = local_opt(varargin,'savefile',1);
gui = local_opt(varargin,'gui',0);

if ~isfield(psyrat_data,'rel') || ~isfield(psyrat_data.rel,'analysis') || ...
        ~strcmpi(psyrat_data.rel.analysis,'ic_dodiff')
    error('varargin:analysis',... %Error code and associated error
        'DoD variance export requires rel.analysis = ''ic_dodiff''.');
end

if ~isnumeric(groupidx) || ~isscalar(groupidx) || isnan(groupidx) || ...
        groupidx < 1 || abs(round(groupidx)-groupidx) > 0
    error('varargin:groupidx',... %Error code and associated error
        'groupidx must be a positive integer scalar.');
end
groupidx = round(groupidx);

if ~isnumeric(map) || numel(map) ~= 4
    error('varargin:map',... %Error code and associated error
        'map must be [ERP_1 ERP_2 ERP_3 ERP_4].');
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
        'map indices exceed available event count.');
end

if isempty(obs)
    obs = local_observed_trials(rel.data,hasgroup,groupnames{groupidx},eventnames(map));
end

if ~isnumeric(obs) || numel(obs) ~= 4 || any(~isfinite(obs)) || any(obs <= 0)
    error('varargin:obs',... %Error code and associated error
        ['obs must be four positive values. ',...
        'How to fix: choose events with recorded trials, or pass ''obs'',[n1 n2 n3 n4].']);
end
obs = obs(:)';

[bp,bt,erlog] = local_groupdraws(rel,groupidx);
bp = bp(:,map,map);
bt = bt(:,map,map);
erlog = erlog(:,map);

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

%The two ICC rows are labeled n'=1 because psyrat_dodiffrel builds their
%denominators from the UNSCALED contrast variances (its :151-152), so unlike the
%Dependability/Generalizability rows above them they do not depend on the
%entered trial counts. Nothing here changes what is computed.
rows = {...
    'Universe variance of DoD contrast'
    'Between-trial contrast variance (unscaled)'
    'Between-trial contrast variance (scaled by obs)'
    'Error contrast variance (unscaled)'
    'Error contrast variance (scaled by obs)'
    'Dependability'
    'Generalizability'
    'ICC (Dependability, n''=1)'
    'ICC (Generalizability, n''=1)'
    };

pt = [...
    dep.uni_var_pt
    dep.bt_contrast_var_pt
    dep.bt_scaled_var_pt
    dep.er_contrast_var_pt
    dep.er_scaled_var_pt
    dep.pt
    gen.pt
    dep.icc_pt
    gen.icc_pt
    ];

ll = [...
    dep.uni_var_ll
    dep.bt_contrast_var_ll
    dep.bt_scaled_var_ll
    dep.er_contrast_var_ll
    dep.er_scaled_var_ll
    dep.ll
    gen.ll
    dep.icc_ll
    gen.icc_ll
    ];

ul = [...
    dep.uni_var_ul
    dep.bt_contrast_var_ul
    dep.bt_scaled_var_ul
    dep.er_contrast_var_ul
    dep.er_scaled_var_ul
    dep.ul
    gen.ul
    dep.icc_ul
    gen.icc_ul
    ];

varcomptable = table(rows,pt,ll,ul,...
    'VariableNames',{'Component','Point','Low_CrI','High_CrI'});

meta = struct;
meta.group = groupnames{groupidx};
meta.event_labels = eventnames(map);
meta.contrast = sprintf('(%s - %s) - (%s - %s)',...
    eventnames{map(1)},eventnames{map(2)},eventnames{map(3)},eventnames{map(4)});
meta.obs = obs;

%Hand the admissibility record back so a caller that is assembling several
%evaluations -- psyrat_report's per-group DoD sweep -- can aggregate them into
%one note instead of re-running the kernel to get at it (RC-43). dep only: the
%dep and gen records are bit-identical, so returning both would double every
%count downstream. See psyrat_dod_admissibility_note.
adm = dep.admissibility;

%Build the note once, for whichever surface this call is headed to. dep only,
%because the dep and gen records are bit-identical -- see
%psyrat_dod_admissibility_note. Empty for a clean run, so both branches below
%stay exactly as they were before this existed.
adm_note = psyrat_dod_admissibility_note(dep);
meta.adm_note = adm_note;

if gui == 1
    local_show_gui(psyrat_data,varcomptable,meta);

    %This table prints the variance components themselves alongside the two
    %coefficients, so it is exactly where a negative one is visible and needs
    %explaining.
    if ~isempty(adm_note)
        warndlg(adm_note,'Coefficient is not interpretable as a reliability');
    end
elseif savefile == 1
    %The branch the app actually takes. psyrat_startview_dodiff's "Export DoD
    %Var Components" button calls this with 'savefile',1,'gui',0, and
    %psyrat_report passes 'gui',0,'savefile',0 -- so nothing in production
    %reaches the warndlg above, and gating the note on `gui` alone would have
    %left the one route that produces a file the user KEEPS with no warning on
    %it. A console warning does not survive into an exported file, which is the
    %same reason psyrat_report folds its note into the report header. Without
    %this, the button wrote a negative between-trial contrast variance and a
    %dependability outside [0,1] into a file with nothing attached (RC-43).
    local_save_table(psyrat_data,varcomptable,meta);
end

end

function local_show_gui(psyrat_data,varcomptable,meta)
figwidth = 980;
figheight = 420;
rowspace = 25;
row = figheight - rowspace*2;

fig = psyrat_newfigure('unit','pix',...
    'position',[1120 520 figwidth figheight],...
    'menub','no',...
    'Tag','psyrat_output',...
    'name','DoD Variance Components',...
    'numbertitle','off',...
    'resize','off');

uicontrol(fig,'Style','text','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Difference-of-Differences Variance Components (Summary + CrI)',...
    'Position',[0 row figwidth 24]);

uicontrol(fig,'Style','text','fontsize',12,...
    'HorizontalAlignment','left',...
    'String',['Group: ' meta.group],...
    'Position',[20 row-28 300 22]);
uicontrol(fig,'Style','text','fontsize',12,...
    'HorizontalAlignment','left',...
    'String',['Contrast: ' meta.contrast],...
    'Position',[20 row-52 figwidth-30 22]);

t = uitable('Parent',fig,'Position',[20 82 figwidth-40 238],...
    'Data',table2cell(varcomptable));
set(t,'ColumnName',varcomptable.Properties.VariableNames);
set(t,'ColumnWidth',{420 120 120 120});
set(t,'RowName',[]);

uicontrol(fig,'Style','push','fontsize',13,...
    'HorizontalAlignment','center',...
    'String','Save Table',...
    'Position',[figwidth/8 20 figwidth/4 44],...
    'Callback',{@local_save_callback,psyrat_data,varcomptable,meta});
end

function local_save_callback(varargin)
psyrat_data = varargin{3};
varcomptable = varargin{4};
meta = varargin{5};
local_save_table(psyrat_data,varcomptable,meta);
end

function local_save_table(psyrat_data,varcomptable,meta)
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
head{end+1} = 'Contrast form: (ERP_1 - ERP_2) - (ERP_3 - ERP_4)';
head{end+1} = sprintf('Selected contrast: %s',meta.contrast);
head{end+1} = sprintf('ERP_1: %s',meta.event_labels{1});
head{end+1} = sprintf('ERP_2: %s',meta.event_labels{2});
head{end+1} = sprintf('ERP_3: %s',meta.event_labels{3});
head{end+1} = sprintf('ERP_4: %s',meta.event_labels{4});
head{end+1} = sprintf('Obs vector used: [%0.3f %0.3f %0.3f %0.3f]',...
    meta.obs(1),meta.obs(2),meta.obs(3),meta.obs(4));
%run provenance (seed/engine/priors/convergence); fail-open
head = [head, psyrat_provenance_lines(psyrat_data)'];

%Admissibility note, when a quantity that had to be a variance was not. This is
%the only surfacing this route gets: it writes a file and raises no dialog, so
%without the note a user keeps a spreadsheet showing a negative variance and an
%out-of-range coefficient with nothing to explain them (RC-43). Split on
%newline because the note is newline-joined and each header line is written as
%its own record. Absent entirely for a clean run.
if isfield(meta,'adm_note') && ~isempty(meta.adm_note)
    head{end+1} = '';
    notelines = strsplit(meta.adm_note, newline);
    head = [head, notelines];
end

head{end+1} = '';
head{end+1} = '';

if strcmpi(ext,'.xlsx')
    writecell(head',outfile);
    writetable(varcomptable,outfile,'Range',strcat('A',num2str(length(head))));
elseif strcmpi(ext,'.csv')
    fid = fopen(outfile,'w');
    %quote per the CSV convention (B19): a comma-bearing header line must
    %stay ONE spreadsheet field; comma-free lines stay byte-identical
    for i = 1:numel(head)
        fprintf(fid,'%s\n',psyrat_quote_csv_header(head{i}));
    end
    fprintf(fid,'Component,Point,Low CrI,High CrI\n');
    for i = 1:height(varcomptable)
        fprintf(fid,'%s,%0.10g,%0.10g,%0.10g\n',...
            local_csv_scalar(varcomptable.Component{i}),...
            varcomptable.Point(i),...
            varcomptable.Low_CrI(i),...
            varcomptable.High_CrI(i));
    end
    fclose(fid);
end
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

function obs = local_observed_trials(datatable,hasgroup,gname,event_labels)
if hasgroup
    gmask = string(datatable.group) == string(gname);
    groupsub = datatable(gmask,:);
else
    groupsub = datatable;
end

obs = zeros(1,numel(event_labels));
for e = 1:numel(event_labels)
    emask = string(groupsub.event) == string(event_labels{e});
    subset = groupsub(emask,:);
    if isempty(subset)
        obs(e) = 0;
    else
        trltable = varfun(@length,subset,'GroupingVariables',{'id'});
        obs(e) = mean(trltable.GroupCount);
    end
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
else
    out = char(string(v));
end
out = strrep(out,'"','""');
if contains(out,',') || contains(out,'"') || contains(out,newline)
    out = ['"' out '"'];
end
end
