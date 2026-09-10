function psyrat_startview_dynrel(varargin)
%Viewer for dynamic/conditional reliability results (Rast & Clayson, the dynamic
%reliability family, analyses 11-22 and 26-29: non-difference one-facet/two-facet
%at the group or subject level, group-level and subject-level difference scores,
%the concurrent residual-correlation variants, and the person-specific
%NONCONCURRENT difference-of-differences designs). Lets the user choose the
%coefficient (generalizability G(z) or dependability D(z)), the trial count n',
%the credible-interval width, an optional reference cutoff, and (for the
%non-difference variants only) the typical-person vs population-average estimand,
%then view the reliability surface and the variance-component table or export the
%table. The two-facet variants add a reliability-coefficient
%(equivalence/stability/both) and occasion-n' selector.
%For the subject-level variants (the difference analyses 13/16/19/20/21/22, the
%non-difference analyses 26/27, and the DoD analyses 28/29) it also exposes a
%per-participant reliability table (each subject's phi at their own standardized
%dimension value) to view or export.
%
%psyrat_startview_dynrel('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data)
%
%Required Input
% psyrat_data - PsyRAT Toolbox data structure array (REL.analysis = 'ic_dynrel'
%   or another dynrel-family analysis string)
% psyrat_prefs - PsyRAT Toolbox preferences structure array
%
%Output
% No variables are returned. Figures/tables are created on demand and can be
%  saved from the viewer.

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

[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

if ~isfield(psyrat_prefs,'guis') || ~isfield(psyrat_prefs.guis,'fsize')
    psyrat_prefs.guis.fsize = psyrat_guifsize;
end
fsize = psyrat_prefs.guis.fsize;

REL = psyrat_data.rel;
ndim = REL.ndim;

%the subject-level dynamic difference variants additionally expose a
%per-participant reliability table (each subject's phi at their own z):
%non-concurrent one-facet (13) / two-facet (16), CONCURRENT population-rescor
%one-facet (19) / two-facet (20), and CONCURRENT per-subject-rho one-facet (21) /
%two-facet (22).
issubjtrt = strcmp(REL.analysis,'ic_diff_dynrel_sserrvar_trt');
issubjrescor = strcmp(REL.analysis,'ic_diff_dynrel_sserrvar_rescor');
issubjtrtrescor = strcmp(REL.analysis,'ic_diff_dynrel_sserrvar_trt_rescor');
issubjrho = strcmp(REL.analysis,'ic_diff_dynrel_sserrvar_rho');
issubjtrtrho = strcmp(REL.analysis,'ic_diff_dynrel_sserrvar_trt_rho');
%the NON-difference subject-level dynamic variants (one-facet analysis 26 =
%'ic_dynrel_sserrvar', two-facet analysis 27 = 'ic_dynrel_sserrvar_trt') ALSO
%expose a per-participant table (each subject's phi_s(z)); unlike every
%difference subject-level variant they are non-difference, so they also keep the
%estimand selector below (isnondiff).
issubjnondiff = strcmp(REL.analysis,'ic_dynrel_sserrvar');
issubjnondiff_trt = strcmp(REL.analysis,'ic_dynrel_sserrvar_trt');
%the person-specific dynamic NONCONCURRENT difference-of-differences (analyses
%28/29) are subject-level by construction (each participant's own residual SD
%is the estimand), so both expose the per-participant table.
isdod = strcmp(REL.analysis,'ic_dodiff_dynrel_sserrvar');
isdodtrt = strcmp(REL.analysis,'ic_dodiff_dynrel_sserrvar_trt');
issubj = strcmp(REL.analysis,'ic_diff_dynrel_sserrvar') || issubjtrt ...
    || issubjrescor || issubjtrtrescor || issubjrho || issubjtrtrho ...
    || issubjnondiff || issubjnondiff_trt || isdod || isdodtrt;

%the trial + occasion two-facet variants (non-difference analysis 14, group-level
%difference analysis 15, subject-level difference analysis 16) additionally expose
%a reliability-coefficient selector (the three test-retest reltypes) and an
%occasion-n' selector (1 occasion by default, or the observed number).
istrt = strcmp(REL.analysis,'ic_dynrel_trt');
isdifftrt = strcmp(REL.analysis,'ic_diff_dynrel_trt');
%the concurrent (rescor) two-facet group-level difference variant (analysis 18)
%also uses the two-facet selectors; the one-facet concurrent variant (analysis 17)
%does not. The two-facet subject-level concurrent variants (population-rescor
%analysis 20, per-subject-rho analysis 22) also use the selectors.
isdifftrtrescor = strcmp(REL.analysis,'ic_diff_dynrel_trt_rescor');
istwofacet = istrt || isdifftrt || issubjtrt || isdifftrtrescor ...
    || issubjtrtrescor || issubjtrtrho || issubjnondiff_trt ...
    || isdodtrt; %two-facet variants use the selectors

%the estimand selector (hmarg: typical person vs population average) averages
%over the person-scale random effect and is meaningful only for the NON-difference
%variants (one-facet analysis 11 = 'ic_dynrel', trial+occasion analysis 14 =
%'ic_dynrel_trt'). For every difference variant the reported surface is the
%typical-person (delta_p = 0) reference, so the selector is disabled there to
%avoid implying a population-average difference estimand that is not computed.
isnondiff = strcmp(REL.analysis,'ic_dynrel') || strcmp(REL.analysis,'ic_dynrel_trt') ...
    || issubjnondiff || issubjnondiff_trt;

%close the setup gui if open
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end

figwidth = 560;
rowspace = 35;
figheight = 470;
if issubj
    figheight = figheight + rowspace*1.1; %room for the per-participant buttons
end
if istwofacet
    figheight = figheight + rowspace*2; %room for the reltype + occasion-n' rows
end
figheight = figheight + rowspace*1.6; %room for the navigation row (B18; 50-px buttons, sibling height)
row = figheight - rowspace*1.6;
lcol = 30;
rcol = (figwidth/9)*5;

psyrat_gui = psyrat_newfigure('unit','pix',...
    'position',[400 400 figwidth figheight],...
    'menub','no','name','Dynamic Reliability: View Results',...
    'numbertitle','off','resize','off');

%TAG REQUIRED (B48): psyrat_start's teardown matches on a 'psyrat_' Tag
%prefix, so an untagged viewer would survive it. The Tag is DELIBERATELY NOT
%'psyrat_gui' -- that name means the setup screen, which this function looks
%up and closes above, and reusing it would make that lookup ambiguous.
psyrat_gui.Tag = 'psyrat_gui_dynrelview';

uicontrol(psyrat_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String',['Dataset:  ' REL.filename],...
    'Position',[0 row figwidth 25]);
row = row - rowspace*.8;

dimstr = strjoin(cellfun(@(x) char(string(x)),REL.dim_names,...
    'UniformOutput',false),', ');
uicontrol(psyrat_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String',sprintf('Dimensions (standardized): %s',dimstr),...
    'Position',[0 row figwidth 25]);
row = row - rowspace;

%coefficient (G vs D)
uicontrol(psyrat_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left','String','Coefficient:',...
    'Position',[lcol row 220 25]);
hcoeff = uicontrol(psyrat_gui,'Style','pop','fontsize',fsize,...
    'String',{'Generalizability G(z)','Dependability D(z)'},'Value',1,...
    'Position',[rcol row 200 28]);
row = row - rowspace;

%trial count n'
uicontrol(psyrat_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String','Trials n'' (blank = per-group median):',...
    'Position',[lcol row 320 25]);
hobs = uicontrol(psyrat_gui,'Style','edit','fontsize',fsize,...
    'String','','Position',[rcol row 200 28]);
row = row - rowspace;

%credible interval
uicontrol(psyrat_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left','String','Credible interval (0-1):',...
    'Position',[lcol row 220 25]);
hci = uicontrol(psyrat_gui,'Style','edit','fontsize',fsize,...
    'String','0.95','Position',[rcol row 200 28]);
row = row - rowspace;

%reference cutoff
uicontrol(psyrat_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left','String','Reference cutoff (blank = none):',...
    'Position',[lcol row 300 25]);
hcut = uicontrol(psyrat_gui,'Style','edit','fontsize',fsize,...
    'String','','Position',[rcol row 200 28]);
row = row - rowspace;

%estimand: typical person vs population average (disabled for difference
%variants, where only the typical-person reference surface is computed).
if isnondiff; marg_enable = 'on'; else; marg_enable = 'off'; end
uicontrol(psyrat_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left','String','Estimand:','Enable',marg_enable,...
    'Position',[lcol row 220 25]);
hmarg = uicontrol(psyrat_gui,'Style','pop','fontsize',fsize,...
    'String',{'Typical person (delta_p = 0)','Population average (lognormal)'},...
    'Value',1,'Enable',marg_enable,'Position',[rcol row 230 28]);

%trial + occasion two-facet: coefficient reltype + occasion n' selectors. The
%popup Value (1/2/3) maps directly to the psyrat_rel_trt reltype; the default is
%3 (both facets random), the brms-notes dynamic coefficient.
hreltype = [];
hnocc = [];
if istwofacet
    row = row - rowspace;
    uicontrol(psyrat_gui,'Style','text','fontsize',fsize,...
        'HorizontalAlignment','left','String','Reliability coefficient:',...
        'Position',[lcol row 220 25]);
    hreltype = uicontrol(psyrat_gui,'Style','pop','fontsize',fsize,...
        'String',{'Equivalence (occasion fixed)',...
        'Stability (trial fixed)',...
        'Equivalence + stability (both random)'},'Value',3,...
        'Position',[rcol row 260 28]);
    row = row - rowspace;
    uicontrol(psyrat_gui,'Style','text','fontsize',fsize,...
        'HorizontalAlignment','left','String','Occasions n'':',...
        'Position',[lcol row 220 25]);
    hnocc = uicontrol(psyrat_gui,'Style','pop','fontsize',fsize,...
        'String',{'1 occasion','Observed occasions'},'Value',1,...
        'Position',[rcol row 230 28]);
end
row = row - rowspace*1.4;

%buttons
bw = (figwidth - 2*lcol - 20)/2;
uicontrol(psyrat_gui,'Style','push','fontsize',fsize,...
    'String','Plot reliability surface',...
    'Position',[lcol row bw 32],'Callback',@plotSurface);
uicontrol(psyrat_gui,'Style','push','fontsize',fsize,...
    'String','Show variance components',...
    'Position',[lcol+bw+20 row bw 32],'Callback',@showVarcomp);
row = row - rowspace*1.1;
uicontrol(psyrat_gui,'Style','push','fontsize',fsize,...
    'String','Save variance table (CSV)',...
    'Position',[lcol row bw 32],'Callback',@saveVarcomp);
uicontrol(psyrat_gui,'Style','push','fontsize',fsize,...
    'String','Done',...
    'Position',[lcol+bw+20 row bw 32],'Callback',@(~,~) close(psyrat_gui));

%subject-level variant: per-participant reliability table view/export
if issubj
    row = row - rowspace*1.1;
    uicontrol(psyrat_gui,'Style','push','fontsize',fsize,...
        'String','Show per-participant table',...
        'Position',[lcol row bw 32],'Callback',@showSubjTable);
    uicontrol(psyrat_gui,'Style','push','fontsize',fsize,...
        'String','Save per-participant table (CSV)',...
        'Position',[lcol+bw+20 row bw 32],'Callback',@saveSubjTable);
end

%navigation row (B18): every sibling viewer carries these two buttons; the
%setup GUI is already closed on entry, so without them Done strands the user
%at a bare prompt with no way back or to another stored result. Same
%local-callback idiom as psyrat_startview_sing/_trt. The two-line label needs
%the siblings' 50-px button height (btn_height = 50 in sing/trt/dodiff); 32
%clips its second line, so the row step and figheight reservation are 1.6.
row = row - rowspace*1.6;
uicontrol(psyrat_gui,'Style','push','fontsize',fsize,...
    'String','Back to Home',...
    'Position',[lcol row bw 50],'Callback',{@psyrat_svb,psyrat_gui});
uicontrol(psyrat_gui,'Style','push','fontsize',fsize,...
    'String',{'Open Another';'.psyrat File'},...
    'Position',[lcol+bw+20 row bw 50],...
    'Tooltip','Open a different processed PsyRAT (.psyrat) file',...
    'Callback',{@psyrat_view_loadnewfile,psyrat_gui});

% ---- nested callbacks (share REL/fsize) -----------------------------------

%cache the computed surface so repeated clicks (plot / show table / save) do not
%recompute the full per-stratum posterior surface; recompute only when the
%CI / n' / estimand inputs change.
cachedSumm = [];
cachedKey = '';

    function summ = buildSummary()
        ciperc = local_parse_pos(get(hci,'String'),0.95);
        if ciperc <= 0 || ciperc >= 1; ciperc = 0.95; end
        obsstr = strtrim(get(hobs,'String'));
        if isempty(obsstr); obs = []; else; obs = local_parse_pos(obsstr,[]); end
        marginal = get(hmarg,'Value') - 1; %1->0 typical, 2->1 marginal

        %trial + occasion two-facet selectors (reltype + occasion n'); inert for
        %the other dynrel analyses (handles are empty there).
        extra = {};
        relKey = '';
        if istwofacet
            reltype = get(hreltype,'Value');     %1/2/3 -> psyrat_rel_trt reltype
            if get(hnocc,'Value') == 2
                noccArg = 'observed';
            else
                noccArg = [];                    %default: 1 occasion
            end
            extra = {'reltype',reltype,'nocc',noccArg};
            if isempty(noccArg); noccLbl = '1'; else; noccLbl = noccArg; end
            relKey = sprintf('|rt%d|nocc%s',reltype,noccLbl);
        end

        %num2str at DEFAULT precision keeps only about six significant digits,
        %so two different requests could share a key while the summary resolved
        %them to different counts: it prints 20.499999 and 20.5 both as '20.5',
        %but max(1,round(obs)) turns those into n' = 20 vs 21. The second save
        %then reused the first surface, and the plot, table and export header
        %all reported the stale n'. 17 digits round-trips any double, so
        %distinct requests now get distinct keys.
        if isempty(obs); obsKey = 'auto'; else; obsKey = num2str(obs,17); end
        key = sprintf('%.6f|%s|%d%s',ciperc,obsKey,marginal,relKey);
        if strcmp(key,cachedKey) && ~isempty(cachedSumm)
            summ = cachedSumm;
            return;
        end

        summ = psyrat_dynrel_summary(REL,'CI',ciperc,'obs',obs,...
            'marginal',marginal,extra{:});
        %record the n' REQUEST alongside the resolved per-stratum counts so
        %the export header can label their source: [] = blank field (each
        %stratum's median default), a scalar = user override. Viewer-local
        %annotation; the summary contract stores only the RESOLVED counts.
        summ.obs_request = obs;
        cachedSumm = summ;
        cachedKey = key;
    end

    function plotSurface(~,~)
        summ = buildSummary();
        gcoeff = 3 - get(hcoeff,'Value'); %1->G(2), 2->D(1)
        cutoff = local_parse_pos(get(hcut,'String'),[]);
        psyrat_dynrelplot(summ,'gcoeff',gcoeff,'cutoff',cutoff);
    end

    function T = concatVarcomp(summ)
        %Stack the per-stratum variance-component tables via the shared builder
        %(also used by the headless report), prepending the stratum label and its
        %observed trial count.
        T = psyrat_dynrel_varcomp_concat(summ);
    end

    function showVarcomp(~,~)
        summ = buildSummary();
        T = concatVarcomp(summ);
        f = figure('Name','Dynamic reliability: variance components',...
            'NumberTitle','off','Tag','psyrat_output',...
            'Position',[200 200 870 360]);
        %Modular cut-copula runs must disclose ON SCREEN, not only in the
        %export header: without this note a user reads rho_e,12(copula) and
        %its CI off the table with nothing saying it is a cut-posterior
        %conditional interval (2026-08-11 code review). The table is shortened
        %to make room only when the note is present, so every other design's
        %window is unchanged.
        tablepos = [0.02 0.02 0.96 0.96];
        if ~isempty(summ.copula_labels)
            tablepos = [0.02 0.13 0.96 0.85];
            uicontrol(f,'Style','text','Units','normalized',...
                'Position',[0.02 0.01 0.96 0.11],...
                'HorizontalAlignment','left',...
                'String',['Modular CUT posterior: rho_e,12(copula) is '...
                'estimated conditionally on the retained margin draws with '...
                'no copula feedback; intervals are not full joint-posterior '...
                'credible intervals. See the exported header for the full '...
                'provenance.']);
        elseif isfield(summ,'dod_labels') && ~isempty(summ.dod_labels)
            %The nonconcurrent DoD estimand must disclose ON SCREEN for the
            %same reason the modular-cut note does: without it a user reads
            %the components off the table with nothing saying the four
            %constituents were TREATED as measured on distinct trials. Same
            %table-shortening mechanics; every other design's window is
            %byte-unchanged. The family noun follows summ.isgamma - both
            %families reach here since the 2026-08-20 Gaussian arm, and a
            %Gaussian run must not print "Gamma margins" (the adversarial
            %pre-merge review caught exactly that hardcoding).
            if isfield(summ,'isgamma') && summ.isgamma
                famnoun = 'Gamma';
            else
                famnoun = 'Gaussian';
            end
            tablepos = [0.02 0.13 0.96 0.85];
            uicontrol(f,'Style','text','Units','normalized',...
                'Position',[0.02 0.01 0.96 0.11],...
                'HorizontalAlignment','left',...
                'String',['Nonconcurrent counterfactual estimand: the four '...
                'constituents are treated as if measured on distinct '...
                'trials; all six conditional residual covariances are '...
                'fixed to zero by assumption. This is the full Bayesian '...
                'posterior of four conditionally independent ' famnoun ...
                ' margins (no dependence module exists) - not a concurrent '...
                'reliability, and not a cut posterior. See the exported '...
                'header for the full provenance.']);
        end
        %Component column widened to fit the longest label
        %("Between-person SD of residual correlation") used by the per-subject-rho
        %variants (cases 21/22).
        uitable(f,'Data',[cellstr(string(T.stratum)) cellstr(string(T.component)) ...
            cellstr(string(T.symbol)) num2cell(round(T.estimate,4)) ...
            num2cell(round(T.ci_lower,4)) num2cell(round(T.ci_upper,4))],...
            'ColumnName',{'Stratum','Component','Symbol','Estimate',...
            'CI lower','CI upper'},...
            'Units','normalized','Position',tablepos,...
            'ColumnWidth',{120 320 90 80 80 80});
    end

    function saveVarcomp(~,~)
        summ = buildSummary();
        T = concatVarcomp(summ);
        exportDynrelTable(T,summ,'dynrel_variance_components.xlsx',...
            'variance components','Save variance components',false);
    end

    function T = concatSsrel(summ)
        %stack the per-stratum per-participant tables with a stratum label via the
        %shared builder (also used by the headless report). Returns [] when no
        %stratum has a per-participant table; the callers' isempty guards handle
        %both [] and an empty table identically.
        T = psyrat_dynrel_ssrel_concat(summ);
    end

    function showSubjTable(~,~)
        summ = buildSummary();
        T = concatSsrel(summ);
        if isempty(T)
            warndlg('No per-participant table available.','Subject-level');
            return;
        end
        [data,cols,widths] = local_ssrel_display(T,ndim);
        f = figure('Name',...
            'Dynamic reliability: per-participant',...
            'NumberTitle','off','Tag','psyrat_output',...
            'Position',[200 200 920 380]);
        uitable(f,'Data',data,'ColumnName',cols,...
            'Units','normalized','Position',[0.02 0.02 0.96 0.96],...
            'ColumnWidth',widths);
    end

    function saveSubjTable(~,~)
        summ = buildSummary();
        T = concatSsrel(summ);
        if isempty(T)
            warndlg('No per-participant table available.','Subject-level');
            return;
        end
        exportDynrelTable(T,summ,'dynrel_diff_subjectlevel.xlsx',...
            'per-participant reliability','Save per-participant reliability',true);
    end

    function exportDynrelTable(T,summ,defaultfn,label,dlgtitle,isperpart)
        %Write a dynamic-reliability table to XLSX or CSV with a metadata header,
        %mirroring the Mac-safe writecell/writetable pattern used by the standard
        %viewers (see psyrat_dod_varcompt.m). XLSX is offered first so the default
        %extension matches the standard viewers; CSV remains available.
        %
        %isperpart says whether T is the PER-PARTICIPANT table (true) or the
        %variance-component table (false). It gates the conditional-estimand
        %header line, which describes the per-participant table specifically and
        %must not be written onto the variance-component export. It is a property
        %of the table being written, NOT of the run - gating it on a run-level
        %flag such as issubj puts the line on both exports.
        if nargin < 6; isperpart = false; end
        [fn,pth] = uiputfile(...
            {'*.xlsx','Excel file (*.xlsx)';'*.csv','CSV file (*.csv)'},...
            dlgtitle,defaultfn);
        if isequal(fn,0); return; end
        outfile = fullfile(pth,fn);
        [~,~,ext] = fileparts(outfile);

        %metadata header (one row-cell per line); writecell(head',...) writes each
        %entry down its own spreadsheet row. Two trailing blanks for spacing.
        head = {};
        head{end+1} = sprintf('PsyRAT dynamic reliability: %s',label);
        head{end+1} = sprintf('Generated: %s',datestr(clock));
        if isfield(psyrat_data,'ver')
            head{end+1} = sprintf('PsyRAT Toolbox v%s',psyrat_data.ver);
        end
        head{end+1} = sprintf('Dataset: %s',REL.filename);
        head{end+1} = sprintf('Dimensions (standardized): %s',dimstr);
        if isfield(REL,'nchains') && isfield(REL,'niter')
            head{end+1} = psyrat_iterline(REL);
        end
        head{end+1} = sprintf('Credible interval: %.3g',summ.ci);
        %the nonconcurrent DoD estimand lines apply in EVERY family (gamma
        %location-scale and, since the 2026-08-20 owner ruling, Gaussian), so
        %the DoD tests run BEFORE the family branch: a Gaussian DoD result
        %(summ.isgamma false) would otherwise fall through to the generic
        %typical-person line and silently drop the "reference curve only" and
        %per-participant conditional sentences. The wording is identical
        %across families - the per-cell log residual SD conditioning
        %statement is true for both, because the scale submodel is log-linked
        %in both (SCIENTIFIC_FORMULA_AUDIT.md section 26).
        %
        %The nonconcurrent DoD surfaces (analyses 28/29) are the
        %TYPICAL-PERSON estimand - all person effects, location AND scale, at
        %zero. Presence-keyed on the dod_labels carrier, matching the
        %disclosure doctrine.
        isdodsum = isfield(summ,'dod_labels') && ~isempty(summ.dod_labels);
        %B15 (2026-08-16): the DoD designs are identified by the ANALYSIS
        %string, not by the labels carrier. A DoD result whose dod_labels
        %is missing (damaged or pre-release stored file) must NOT fall
        %through to a family branch's estimand sentence - decline to assert
        %an estimand instead, mirroring the string-keyed DEGRADED fallback
        %psyrat_provenance_lines carries for the identical condition.
        %Mirrors psyrat_report.
        isdodrun = any(strcmp(REL.analysis, ...
            {'ic_dodiff_dynrel_sserrvar','ic_dodiff_dynrel_sserrvar_trt'}));
        if isdodsum
            head{end+1} = ['Estimand: surface = typical person (all '...
                'person effects at zero); reference curve only - the '...
                'per-participant table is the primary output'];
            if isperpart
                head{end+1} = ['Estimand: per-participant table = '...
                    'conditional on each person''s own per-cell log '...
                    'residual SDs, at their own dimension value(s)'];
            end
        elseif isdodrun
            head{end+1} = ['Estimand: cannot be determined from this '...
                'stored result - the dod_labels carrier is missing '...
                '(damaged or pre-release file); see the provenance '...
                'lines for the string-keyed disclosure'];
            if isperpart
                %same degraded arm as the surface line: without the
                %carrier, the per-participant estimand sentence would
                %otherwise claim a family branch's conditional wording
                head{end+1} = ['Estimand: per-participant table = '...
                    'undetermined without the dod_labels carrier (see '...
                    'the surface line above)'];
            end
        elseif isfield(summ,'isgamma') && summ.isgamma
            %Neither gamma parameterization has a delta_p, so the Gaussian
            %typical-person label never applies: the surface is always marginal
            %over the person scale population, and the per-participant table
            %exported alongside it is the complementary CONDITIONAL read-out.
            %Say which is which.
            %
            %WHAT the two marginalize and condition on DIFFERS, so the wording
            %below is branched rather than shared. Under gammascale = 1 there is
            %no log-residual submodel at all and the person quantity is the
            %log-nu dispersion. Under gammascale = 2 there IS one - that is the
            %whole change - and the per-participant table conditions on each
            %person's own log RESIDUAL SD. Their own log-mean does not enter it,
            %because the residual is decoupled from the mean there
            %(psyrat_ssrel_dynrel_gamma_ls). Mirrors psyrat_report.
            isls = isfield(summ,'isgamma_ls') && summ.isgamma_ls;
            if isls
                head{end+1} = ['Estimand: surface = population average over '...
                    'person residual-SD heterogeneity (gamma location-scale)'];
            else
                head{end+1} = ['Estimand: surface = population average over person '...
                    'dispersion (gamma estimand #2)'];
            end
            if isperpart
                if isls
                    head{end+1} = ['Estimand: per-participant table = conditional '...
                        'on each person''s own log residual SD'];
                else
                    head{end+1} = ['Estimand: per-participant table = conditional on '...
                        'each person''s own log-mean and log-nu'];
                end
            end
        elseif isfield(summ,'marginal') && summ.marginal == 1
            head{end+1} = 'Estimand: population average (lognormal)';
        else
            head{end+1} = 'Estimand: typical person (delta_p = 0)';
        end
        %D-study n' the reliability surface was evaluated at, one line per
        %stratum, from the RESOLVED counts in summ.strata (.obs/.obs_ev/
        %.nocc). Nothing in the export previously recorded these, and the S17
        %curve-description lines appended below (emitted only for the gamma
        %subject-level dynamic difference design, analysis 13 -- the
        %psyrat_provenance_lines gate is an exact match on
        %'ic_diff_dynrel_sserrvar') describe the curve as using the median n',
        %which is only true without an override -- so under an override a
        %reader could not tell what the surface was evaluated at (pre-existing
        %gap adjudicated 2026-08-11). The gamma difference variant computes at
        %PER-EVENT n' (st.obs_ev); report the exact pair when the two differ,
        %since the scalar st.obs is then their mean and names a count neither
        %event has (the same rationale as psyrat_dynrel_varcomp_concat's
        %per-event columns; when the pair is equal the scalar IS both events'
        %count, so the plain line suffices).
        if isfield(summ,'obs_request') && ~isempty(summ.obs_request)
            src = 'user override';
        else
            src = 'default: per-stratum median';
        end
        for s = 1:numel(summ.strata)
            st = summ.strata(s);
            if isfield(st,'obs_ev') && numel(st.obs_ev) == 4
                %DoD variants: four per-cell counts, always itemized (the
                %scalar is their mean and rarely names a count any cell has)
                nstr = sprintf(['cell 1 = %d, cell 2 = %d, cell 3 = %d, '...
                    'cell 4 = %d (%s)'],st.obs_ev(1),st.obs_ev(2),...
                    st.obs_ev(3),st.obs_ev(4),src);
            elseif isfield(st,'obs_ev') && numel(st.obs_ev) == 2 && ...
                    st.obs_ev(1) ~= st.obs_ev(2)
                nstr = sprintf('event 1 = %d, event 2 = %d (%s)',...
                    st.obs_ev(1),st.obs_ev(2),src);
            else
                nstr = sprintf('%d (%s)',st.obs,src);
            end
            if ~isempty(st.nocc)
                nstr = sprintf('%s; n'' (occasions): %d',nstr,st.nocc);
            end
            head{end+1} = sprintf(...
                'Reliability surface n'' (trials), stratum %s: %s',...
                char(string(st.label)),nstr); %#ok<AGROW>
        end
        %run provenance (seed/engine/priors/convergence); fail-open
        head = [head, psyrat_provenance_lines(psyrat_data)'];
        head{end+1} = '';
        head{end+1} = '';

        if strcmpi(ext,'.xlsx')
            writecell(head',outfile);
            writetable(T,outfile,'Range',strcat('A',num2str(numel(head))));
        else
            %CSV: write the header lines, then append the table (its own column
            %names + data) so arbitrary table columns are handled robustly.
            %'WriteMode','append' defaults WriteVariableNames to false, so force
            %it on to keep the column-name header row (matching the .xlsx branch
            %and the headless psyrat_write_table writer).
            fid = fopen(outfile,'w');
            for i = 1:numel(head)
                hline = head{i};
                %quote per the CSV convention (B19): a comma-bearing
                %provenance line must stay ONE spreadsheet field, or the
                %header the on-screen note points readers to arrives
                %scattered across columns. Comma-free lines stay
                %byte-identical; the XLSX branch (writecell) already quotes.
                if contains(hline,',') || contains(hline,'"')
                    hline = ['"' strrep(hline,'"','""') '"'];
                end
                fprintf(fid,'%s\n',hline);
            end
            fclose(fid);
            writetable(T,outfile,'WriteMode','append','WriteVariableNames',true);
        end
        fprintf('Saved %s to %s\n',label,outfile);
    end

end

function [data,cols,widths] = local_ssrel_display(T,ndim)
%Assemble the per-participant table for the uitable view. Generalizability
%(phi_rho,s) and dependability (phi_delta,s) point estimates with credible
%bounds, per subject, at each subject's own standardized dimension value. The
%difference variants carry per-event trial counts (trls1/trls2); the
%non-difference subject-level variants (analyses 26/27) carry a single trial
%count (trls), plus an occasion count (occs) for the two-facet variant.
strat = cellstr(string(T.stratum));
idc = cellstr(string(T.id));
r = @(x) num2cell(round(x,3));
vn = T.Properties.VariableNames;

%leading standardized-dimension column(s)
if ndim == 2
    zdata = [r(T.z1) r(T.z2)];
    zcols = {'z1','z2'};
else
    zdata = r(T.z1);
    zcols = {'z1'};
end

%trial (and occasion) count column(s): per-CELL (four) for the DoD variants -
%tested FIRST, because a DoD table also carries trls1/trls2 and the two-column
%branch would silently show half its cells; per-event for the two-event
%difference variants; a single trial count for the non-difference one-facet
%variant; the per-occasion trial count plus the occasion count for the
%non-difference two-facet variant.
if ismember('trls4',vn)
    ndata = [num2cell(T.trls1) num2cell(T.trls2) ...
        num2cell(T.trls3) num2cell(T.trls4)];
    ncols = {'Trials 1','Trials 2','Trials 3','Trials 4'};
    if ismember('occs1',vn)
        %all four per-cell occasion counts - one column would silently show
        %cell 1's count for all cells, the same failure the trls4-first
        %ordering above exists to prevent
        ndata = [ndata num2cell(T.occs1) num2cell(T.occs2) ...
            num2cell(T.occs3) num2cell(T.occs4)];
        ncols = [ncols {'Occ 1','Occ 2','Occ 3','Occ 4'}];
    end
elseif ismember('trls1',vn)
    ndata = [num2cell(T.trls1) num2cell(T.trls2)];
    ncols = {'Trials 1','Trials 2'};
elseif ismember('occs',vn)
    ndata = [num2cell(T.trls) num2cell(T.occs)];
    ncols = {'Trials/occ','Occasions'};
else
    ndata = num2cell(T.trls);
    ncols = {'Trials'};
end

cdata = [r(T.gen_pt) r(T.gen_ll) r(T.gen_ul) r(T.dep_pt) r(T.dep_ll) r(T.dep_ul)];
ccols = {'G(z)','G low','G high','D(z)','D low','D high'};

data = [strat idc zdata ndata cdata];
cols = [{'Stratum','ID'} zcols ncols ccols];
widths = num2cell([90 70 repmat(65,1,numel(cols)-2)]);
end

function v = local_parse_pos(str,default)
%Parse a positive scalar from an edit field; return default on empty/invalid.
%isreal guards complex literals like '3-4i': str2double parses them, isfinite
%is true for them, and x > 0 compares only the real part, so the accept
%condition passed. A complex n' then cleared every downstream check and errored
%in quantile (via prctile) instead of falling back to the default here.
v = default;
if isempty(str); return; end
x = str2double(str);
if isreal(x) && isfinite(x) && x > 0; v = x; end
end

function psyrat_svb(varargin)
%back button. takes user back to psyrat_start (same local-callback idiom as
%the sibling viewers -- each viewer file carries its own copy)

%close gui
close(varargin{3});

%go back to psyrat_start
psyrat_start;

end

function psyrat_view_loadnewfile(varargin)
%give user the chance to load a new file

%close gui
close(varargin{3});
psyrat_startview;

end
