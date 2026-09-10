function psyrat_startview_criterion(varargin)
%Prepare criterion-score-specific outputs for one-facet or TRT analyses.
%
%psyrat_startview_criterion('psyrat_prefs',psyrat_prefs,...
%    'psyrat_data',psyrat_data,'analysis','sing')

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

[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%analysis-unit label fragments (trial vs split). The criterion view is reachable
%for a parallel-splits one-facet run, so its trial-labeled controls switch to
%"split"; psyrat_unitlabels defaults to "trial" when no splits flag is present,
%so non-splits output is byte-identical.
usplit = 1;
if isfield(psyrat_data,'rel') && isfield(psyrat_data.rel,'splits')
    usplit = psyrat_data.rel.splits;
end
L = psyrat_unitlabels(usplit);

if isempty(psyrat_prefs) || ~isstruct(psyrat_prefs)
    psyrat_prefs = psyrat_defaults;
end
if ~isfield(psyrat_prefs,'view')
    psyrat_prefs.view = struct;
end
if ~isfield(psyrat_prefs.view,'gcoeff')
    psyrat_prefs.view.gcoeff = 1;
end
if ~isfield(psyrat_prefs.view,'reltype')
    psyrat_prefs.view.reltype = 1;
end
if ~isfield(psyrat_prefs.view,'plotdepline')
    psyrat_prefs.view.plotdepline = 2;
end
if ~isfield(psyrat_prefs.view,'plotrelline')
    psyrat_prefs.view.plotrelline = 2;
end
if ~isfield(psyrat_prefs.view,'plotcriterion')
    psyrat_prefs.view.plotcriterion = 1;
end
if ~isfield(psyrat_prefs.view,'tablecriterion')
    psyrat_prefs.view.tablecriterion = 1;
end
if ~isfield(psyrat_prefs.view,'criterioncutoff')
    psyrat_prefs.view.criterioncutoff = psyrat_criterion_defaultcut(psyrat_data);
end

analysis = 'sing';
ind = find(strcmpi('analysis',varargin),1);
if ~isempty(ind)
    analysis = lower(char(string(varargin{ind+1})));
end
if ~any(strcmp(analysis,{'sing','trt'}))
    analysis = 'sing';
end

if ~isfield(psyrat_prefs,'guis') || ~isfield(psyrat_prefs.guis,'fsize')
    psyrat_prefs.guis.fsize = psyrat_guifsize;
end

%check if psyrat_gui is open.
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end

figwidth = 620;
figheight = 560;
rowspace = 35;
row = figheight - rowspace*2;
lcol = 30;
rcol = 350;

psyrat_gui = psyrat_newfigure('unit','pix',...
  'position',[450 420 figwidth figheight],...
  'menub','no',...
  'name','Criterion Score Outputs Setup',...
  'numbertitle','off',...
  'resize','off');

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',['Dataset:  ' psyrat_data.rel.filename],...
    'Tooltip','Dataset that was used',...
    'Position',[0 row figwidth 25]);

row = row - (rowspace*.45);

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',['Measurement:  ' psyrat_data.proc.measheader],...
    'Tooltip','Measurement that was analyzed',...
    'Position',[0 row figwidth 25]);

row = row - (rowspace*1.35);

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize+1,...
    'HorizontalAlignment','left',...
    'String','Criterion Score Settings',...
    'Position',[lcol row figwidth-lcol 24]);

row = row - rowspace;

str = sprintf(['Criterion score to evaluate cut-score reliability formulas.\n'...
    'This is a score threshold on the measurement scale (not a reliability threshold).']);
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Criterion Cutoff:',...
    'Tooltip',str,...
    'Position',[lcol row figwidth/3 25]);

inputs.cutoff = uicontrol(psyrat_gui,...
    'Style','edit',...
    'fontsize',psyrat_prefs.guis.fsize,...
    'String',num2str(psyrat_prefs.view.criterioncutoff),...
    'Position',[rcol row+6 figwidth/4 25]);

inputs.cut_error = uicontrol(psyrat_gui,'Style','text',...
    'fontsize',max(psyrat_prefs.guis.fsize-1,8),...
    'HorizontalAlignment','left',...
    'ForegroundColor',[0.70 0.00 0.00],...
    'String','Criterion cutoff must be numeric',...
    'Visible','off',...
    'Position',[lcol row-14 figwidth-lcol 16]);

row = row - rowspace;

str = sprintf(['Criterion-score outputs use dependability formulas\n'...
    '(absolute decisions) by design.']);
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Decision Coefficient: Dependability (Absolute)',...
    'Tooltip',str,...
    'Position',[lcol row figwidth-lcol 25]);

row = row - rowspace;

muval = psyrat_criterion_mu_display(psyrat_data);
str = sprintf(['Average score (mu) across currently loaded data used by\n'...
    'criterion dependability calculations.']);
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',sprintf('Average Score (mu): %0.4f',muval),...
    'Tooltip',str,...
    'Position',[lcol row figwidth-lcol 25]);

if strcmp(analysis,'trt')
    row = row - rowspace;
    
    str = sprintf(['For internal consistency, use coefficient of equivalence\n'...
        'For test-retest reliability, use coefficient of stability\n'...
        'For combined consistency/stability, use equivalence + stability']);
    uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
        'HorizontalAlignment','left',...
        'String','Type of Reliability Coefficient:',...
        'Tooltip',str,...
        'Position',[lcol row figwidth/3 25]);

    inputs.reltype = uicontrol(psyrat_gui,...
        'Style','pop',...
        'fontsize',psyrat_prefs.guis.fsize,...
        'String',{'Equivalence' 'Stability' 'Equivalence + Stability'},...
        'Value',psyrat_prefs.view.reltype,...
        'Position',[rcol row figwidth/4 25]);
end

row = row - (rowspace+6);

uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Output Selection',...
    'Position',[rcol row figwidth/4 24]);

uicontrol(psyrat_gui,'Style','text',...
    'fontsize',max(psyrat_prefs.guis.fsize-1,8),...
    'HorizontalAlignment','center',...
    'String','(checked = include)',...
    'Position',[rcol row-14 figwidth/4 20]);

row = row - (rowspace+14);

str = sprintf(['Display a plot showing reliability at the criterion cutoff\n'...
    'as a function of the ' L.numberof ' retained for averaging']);
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String',['Plot: ' L.Units ' vs Reliability at Criterion'],...
    'Tooltip',str,...
    'Position',[lcol row figwidth/2 40]);

inputs.plotcriterion = uicontrol(psyrat_gui,'Style','checkbox',...
    'Value',psyrat_prefs.view.plotcriterion,...
    'Position',[rcol + (figwidth/8) - 12 row+20 24 25]);

row = row - rowspace;

str = sprintf(['Display a table summarizing reliability at the criterion cutoff\n'...
    'using the selected ' L.unit ' summary settings']);
uicontrol(psyrat_gui,'Style','text','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','left',...
    'String','Table: Overall Reliability Summary at Criterion',...
    'Tooltip',str,...
    'Position',[lcol row figwidth/2 40]);

inputs.tablecriterion = uicontrol(psyrat_gui,'Style','checkbox',...
    'Value',psyrat_prefs.view.tablecriterion,...
    'Position',[rcol + (figwidth/8) - 12 row+20 24 25]);

row = row - rowspace*1.6;

btn_margin = 18;
btn_gap = 12;
btn_height = 50;
btn_width = floor((figwidth - (2*btn_margin) - (2*btn_gap)) / 3);
btn_x = btn_margin;

uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String','Back to Results Setup',...
    'Position',[btn_x row btn_width btn_height],...
    'Callback',{@psyrat_criterion_back,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,'analysis',analysis});
btn_x = btn_x + btn_width + btn_gap;

uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',{'Generate';'Criterion Outputs'},...
    'Position',[btn_x row btn_width btn_height],...
    'Tooltip','Generate criterion-specific plots/tables',...
    'Callback',{@psyrat_criterion_run,'psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,'analysis',analysis,'inputs',inputs});
btn_x = btn_x + btn_width + btn_gap;

uicontrol(psyrat_gui,'Style','push','fontsize',psyrat_prefs.guis.fsize,...
    'HorizontalAlignment','center',...
    'String',{'Close PsyRAT';'Outputs'},...
    'Position',[btn_x row btn_width btn_height],...
    'Tooltip','Close only PsyRAT-generated output figures and tables',...
    'Callback',{@psyrat_criterion_closefigs});

psyrat_gui.Tag = 'psyrat_gui';

end

function psyrat_criterion_back(varargin)
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);
ind = find(strcmpi('analysis',varargin),1);
analysis = varargin{ind+1};

psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end

if strcmpi(analysis,'trt')
    psyrat_startview_trt('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
else
    psyrat_startview_sing('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data);
end

end

function psyrat_criterion_run(varargin)
[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);
ind = find(strcmpi('analysis',varargin),1);
analysis = varargin{ind+1};
ind = find(strcmpi('inputs',varargin),1);
inputs = varargin{ind+1};

cutval = str2double(inputs.cutoff.String);
if psyrat_criterion_cutcheck(cutval) ~= 0
    if isfield(inputs,'cut_error') && ishandle(inputs.cut_error)
        inputs.cut_error.Visible = 'on';
    end
    if ishghandle(inputs.cutoff)
        uicontrol(inputs.cutoff);
    end
    return;
end

if isfield(inputs,'cut_error') && ishandle(inputs.cut_error)
    inputs.cut_error.Visible = 'off';
end

psyrat_prefs.view.criterioncutoff = cutval;
psyrat_prefs.view.plotcriterion = inputs.plotcriterion.Value;
psyrat_prefs.view.tablecriterion = inputs.tablecriterion.Value;
if isfield(inputs,'reltype') && ishghandle(inputs.reltype)
    psyrat_prefs.view.reltype = inputs.reltype.Value;
end

psyrat_criterionfigures('psyrat_prefs',psyrat_prefs,...
    'psyrat_data',psyrat_data,...
    'analysis',analysis);

end

function psyrat_criterion_closefigs(varargin)
outfigs = findall(0,'Type','figure','Tag','psyrat_output');
if ~isempty(outfigs)
    close(outfigs);
end
end

function checkout = psyrat_criterion_cutcheck(cutval)
%~isreal joins isnan because str2double parses complex literals like '2-3i',
%for which isnan is false. Such a cutoff was accepted here, passed the isnan-only
%guard in psyrat_criterionfigures, and reached offsq = (mu - cut).^2 in
%psyrat_cutscore_sing, where the posterior summary threw "First argument must be
%an array of real values" from quantile. The cutoff is also written to screen and
%into the exported provenance header with %0.4f, which prints only the real part,
%so an accepted complex cutoff would have been recorded as a different number
%than the one used.
if ~isreal(cutval) || isnan(cutval)
    checkout = 1;
else
    checkout = 0;
end
end

function cutval = psyrat_criterion_defaultcut(psyrat_data)
cutval = 0;
try
    if isfield(psyrat_data,'rel') && isfield(psyrat_data.rel,'out') && ...
            isfield(psyrat_data.rel.out,'mu') && ~isempty(psyrat_data.rel.out.mu)
        cutval = mean(psyrat_data.rel.out.mu(:));
    end
catch
    cutval = 0;
end
end

function muval = psyrat_criterion_mu_display(psyrat_data)
%Estimate the currently loaded average score (mu) for display in setup.
muval = NaN;

try
    if isfield(psyrat_data,'relsummary') && isfield(psyrat_data.relsummary,'data')
        muvals = [];
        for g = 1:length(psyrat_data.relsummary.data.g)
            for e = 1:length(psyrat_data.relsummary.data.g(g).e)
                d = psyrat_data.relsummary.data.g(g).e(e);
                if isfield(d,'mu') && isfield(d.mu,'raw') && ~isempty(d.mu.raw)
                    muvals(end+1,1) = mean(d.mu.raw(:)); %#ok<AGROW>
                end
            end
        end
        if ~isempty(muvals)
            muval = mean(muvals);
        end
    end
    
    if isnan(muval) && isfield(psyrat_data,'rel') && isfield(psyrat_data.rel,'out') && ...
            isfield(psyrat_data.rel.out,'mu') && ~isempty(psyrat_data.rel.out.mu)
        muval = mean(psyrat_data.rel.out.mu(:));
    end
catch
    muval = NaN;
end

if isnan(muval)
    muval = 0;
end
end
