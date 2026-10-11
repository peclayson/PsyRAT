function psyrat_startview_splits(varargin)
%Viewer for nonparallel data-splits results (Rocha Tables 4/5, observed design;
%analysis cases 23 = ic_splits single occasion, 24 = trt_splits multi occasion).
%Each score is the mean of n_i items in a split, and unequal n_i make the splits
%nonparallel. The default summary uses the observed score; numeric split/occasion
%counts request a reference-design projection. There is no D-study projection
%over hypothetical item counts and no cut-score
%search (split means lose the item-level partition; see psyrat_splits_summary and
%manual sections 13.5 and 13.7). The viewer therefore has no
%projection plot and no trial/split-cutoff table; it shows the per-stratum
%dependability/generalizability coefficients, SEMs, and variance components, and
%exports them.
%
%Case23 retains its existing bound interpretation. Case24 reports D conditional
%on the specified score and fitted restricted model, and the conventional G
%partition without a blanket bound or estimation-bias guarantee. Its default
%score equally weights available occasion means; numeric counts are projections.
%
%psyrat_startview_splits('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data)
%
%Required Input
% psyrat_data - PsyRAT data structure array (REL.analysis = 'ic_splits'
%  or 'trt_splits')
% psyrat_prefs - PsyRAT preferences structure array
%
%Output
% No variables are returned. Tables are created on demand and can be saved from
%  the viewer.

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

%single-occasion (ic_splits) has no occasion facet, so the occasion-n' selector
%is inert there; the multi-occasion (trt_splits) design exposes it.
issingle = strcmp(REL.analysis,'ic_splits');

%analysis-unit labels (split, not trial) for the banner and the table headers.
L = psyrat_unitlabels(REL.splits);
if ~issingle
    L.banner = 'Nonparallel splits: available occasions / count projection';
end

%close the setup gui if open
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end

figwidth = 560;
if ~issingle; figwidth = 700; end
rowspace = 35;
figheight = 430;
if ~issingle
    figheight = figheight + rowspace; %room for the occasion-n' row
end
row = figheight - rowspace*1.4;
lcol = 30;
rcol = (figwidth/9)*5;

psyrat_gui = psyrat_newfigure('unit','pix',...
    'position',[400 400 figwidth figheight],...
    'menub','no','name','Data-Splits Reliability: View Results',...
    'numbertitle','off','resize','off');

%TAG REQUIRED (B48): psyrat_start's teardown matches on a 'psyrat_' Tag
%prefix, so an untagged viewer would survive it. NOT 'psyrat_gui' -- that
%name means the setup screen, which this function looks up and closes above.
psyrat_gui.Tag = 'psyrat_gui_splitsview';

uicontrol(psyrat_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String',['Dataset:  ' REL.filename],...
    'Position',[0 row figwidth 25]);
row = row - rowspace*.8;

%analysis-unit banner (split data, observed design)
uicontrol(psyrat_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','center','String',L.banner,...
    'Position',[0 row figwidth 25]);
row = row - rowspace*.9;

% Case24 interpretations are conditional on its stated score and model.
if issingle
    biasline = 'Dependability is exact; generalizability is a downward-biased lower bound.';
    splitprompt = 'Splits n'' (blank = observed design):';
else
    biasline = 'D is conditional on score/model; G uses the conventional relative partition.';
    splitprompt = 'Splits n'' projection (blank = observed):';
end
uicontrol(psyrat_gui,'Style','text','fontsize',fsize-1,...
    'HorizontalAlignment','center','ForegroundColor',[0.45 0.20 0.20],...
    'String',biasline,...
    'Position',[10 row figwidth-20 25]);
row = row - rowspace*1.0;

%Interval width; method follows the saved engine.
uicontrol(psyrat_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left','String',[psyrat_interval_label(REL) ' (0-1):'],...
    'Position',[lcol row 220 25]);
hci = uicontrol(psyrat_gui,'Style','edit','fontsize',fsize,...
    'String','0.95','Position',[rcol row 200 28]);
row = row - rowspace;

%Blank retains each case's observed score; numeric case24 counts project to a
%balanced design. Case23 retains its corrected effective-count convention.
uicontrol(psyrat_gui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String',splitprompt,...
    'Position',[lcol row 320 25]);
hobs = uicontrol(psyrat_gui,'Style','edit','fontsize',fsize,...
    'String','','Position',[rcol row 200 28]);

%occasions n' (multi-occasion design only; blank = observed occasions)
hnocc = [];
if ~issingle
    row = row - rowspace;
    uicontrol(psyrat_gui,'Style','text','fontsize',fsize,...
        'HorizontalAlignment','left',...
        'String','Occasions n'' projection (blank = observed):',...
        'Position',[lcol row 320 25]);
    hnocc = uicontrol(psyrat_gui,'Style','edit','fontsize',fsize,...
        'String','','Position',[rcol row 200 28]);
end
row = row - rowspace*1.4;

%buttons
bw = (figwidth - 2*lcol - 20)/2;
uicontrol(psyrat_gui,'Style','push','fontsize',fsize,...
    'String','Show coefficients',...
    'Position',[lcol row bw 32],'Callback',@showCoef);
uicontrol(psyrat_gui,'Style','push','fontsize',fsize,...
    'String','Show variance components',...
    'Position',[lcol+bw+20 row bw 32],'Callback',@showVarcomp);
row = row - rowspace*1.1;
saveCoefLabel = 'Save coefficients (CSV)';
saveVarLabel = 'Save variance table (CSV)';
if ~issingle
    saveCoefLabel = 'Save coefficients';
    saveVarLabel = 'Save variance table';
end
uicontrol(psyrat_gui,'Style','push','fontsize',fsize,...
    'String',saveCoefLabel,...
    'Position',[lcol row bw 32],'Callback',@saveCoef);
uicontrol(psyrat_gui,'Style','push','fontsize',fsize,...
    'String',saveVarLabel,...
    'Position',[lcol+bw+20 row bw 32],'Callback',@saveVarcomp);
row = row - rowspace*1.1;
uicontrol(psyrat_gui,'Style','push','fontsize',fsize,...
    'String','Done',...
    'Position',[lcol+(bw+20)/2 row bw 32],'Callback',@(~,~) close(psyrat_gui));

% ---- nested callbacks (share REL/fsize/L) ---------------------------------

%cache the computed summary so repeated clicks (coefficients / variance / save)
%do not recompute the per-stratum coefficients; recompute only when the
%CI / n_s / n_o inputs change.
cachedSumm = [];
cachedKey = '';

    function summ = buildSummary()
        ciperc = local_parse_pos(get(hci,'String'),0.95);
        if ciperc <= 0 || ciperc >= 1; ciperc = 0.95; end
        obsstr = strtrim(get(hobs,'String'));
        if isempty(obsstr); obs = []; else; obs = local_parse_pos(obsstr,[]); end

        %occasion n' override (multi-occasion only; [] -> observed occasions).
        nocc = [];
        if ~issingle
            noccstr = strtrim(get(hnocc,'String'));
            if ~isempty(noccstr); nocc = local_parse_pos(noccstr,[]); end
        end

        %Same defect the dynamic-reliability viewer's buildSummary carries (see
        %psyrat_startview_dynrel), and it applies to BOTH counts here.
        %num2str at DEFAULT precision keeps only about six significant digits,
        %so two different requests could share a key while the summary resolved
        %them to different counts: it prints 20.499999 and 20.5 both as '20.5',
        %but psyrat_splits_summary's max(1,round(...)) turns those into n' = 20
        %vs 21 (and likewise 3.499999 vs 3.5 into n'_o = 3 vs 4). The second
        %save then reused the first summary, and its table and export reported
        %the stale design. 17 digits round-trips any double, so distinct
        %requests now get distinct keys.
        if isempty(obs); obsKey = 'auto'; else; obsKey = num2str(obs,17); end
        if isempty(nocc); noccKey = 'auto'; else; noccKey = num2str(nocc,17); end
        key = sprintf('%.6f|%s|%s',ciperc,obsKey,noccKey);
        if strcmp(key,cachedKey) && ~isempty(cachedSumm)
            summ = cachedSumm;
            return;
        end

        summ = psyrat_splits_summary(REL,'CI',ciperc,'obs',obs,'nocc',nocc);
        cachedSumm = summ;
        cachedKey = key;
    end

    function [data,cols,widths] = coefDisplay(summ)
        %Flatten the per-stratum, per-coefficient summary into uitable rows via
        %the shared builder, which is also used by the headless report
        %(psyrat_report) so the on-screen table and both exports match. The
        %enclosing-scope unit label L names the per-person count column.
        [data,cols,widths] = psyrat_splits_coef_rows(summ, L);
    end

    function showCoef(~,~)
        summ = buildSummary();
        [data,cols,widths] = coefDisplay(summ);
        height = 360;
        tablepos = [0.02 0.12 0.96 0.86];
        notepos = [0.02 0.01 0.96 0.09];
        note = summ.gbias_note;
        if ~issingle
            height = 520;
            tablepos = [0.02 0.34 0.96 0.57];
            notepos = [0.02 0.02 0.96 0.28];
            if strcmp(summ.strata(1).design.mode,'balanced_projection')
                scoreline = 'Balanced projection: explicit counts; common reference support; observed missingness is not projected.';
                modeline = 'Design mode: balanced reference projection';
            else
                scoreline = 'Observed score: equally weighted available occasion means; missing occasions receive no score weight.';
                modeline = 'Design mode: observed available occasions';
            end
            note = {scoreline;...
                'Participant universe/error variances are averaged before coefficients and SEMs are calculated.';...
                'D is exact conditional on this score and the fitted shared-ordinal split model with independent row residuals.';...
                'G uses the conventional relative partition; no universal lower-bound, bias or coverage guarantee is asserted.';...
                'Observed CE/CS may have person-specific fixed-facet support; missingness is not modeled.';...
                'Scalar counts are references. Scroll right for actual count ranges, design mode and full weight metadata.'};
        end
        f = figure('Name','Data-splits reliability: coefficients',...
            'NumberTitle','off','Tag','psyrat_output',...
            'Position',[200 200 1080 height]);
        uitable(f,'Data',data,'ColumnName',cols,...
            'Units','normalized','Position',tablepos,...
            'ColumnWidth',widths);
        uicontrol(f,'Style','text','HorizontalAlignment','left',...
            'Units','normalized','Position',notepos,...
            'String',note,'FontSize',fsize-1);
        if ~issingle
            uicontrol(f,'Style','text','HorizontalAlignment','left',...
                'Units','normalized','Position',[0.02 0.93 0.96 0.05],...
                'String',modeline,'FontSize',fsize,'FontWeight','bold');
        end
    end

    function T = concatVarcomp(summ)
        %Stack the per-stratum variance-component tables via the shared builder
        %(also used by the headless report), prepending the stratum label and its
        %observed design (n_s splits, n_o occasions, harmonic-mean items-per-split
        %n_i) so the export records the design each component was reported at.
        T = psyrat_splits_varcomp_concat(summ);
    end

    function showVarcomp(~,~)
        summ = buildSummary();
        T = concatVarcomp(summ);
        f = figure('Name','Data-splits reliability: variance components',...
            'NumberTitle','off','Tag','psyrat_output',...
            'Position',[200 200 920 360]);
        uitable(f,'Data',[cellstr(string(T.stratum)) ...
            cellstr(string(T.component)) cellstr(string(T.symbol)) ...
            num2cell(round(T.estimate,4)) num2cell(round(T.ci_lower,4)) ...
            num2cell(round(T.ci_upper,4))],...
            'ColumnName',{'Stratum','Component','Symbol','Estimate',...
            'CI lower','CI upper'},...
            'Units','normalized','Position',[0.02 0.02 0.96 0.96],...
            'ColumnWidth',{120 320 110 80 80 80});
    end

    function saveCoef(~,~)
        summ = buildSummary();
        [data,cols,~] = coefDisplay(summ);
        T = cell2table(data,'VariableNames',matlab.lang.makeValidName(cols));
        exportSplitsTable(T,summ,'splits_coefficients.xlsx',...
            'reliability coefficients','Save reliability coefficients');
    end

    function saveVarcomp(~,~)
        summ = buildSummary();
        T = concatVarcomp(summ);
        exportSplitsTable(T,summ,'splits_variance_components.xlsx',...
            'variance components','Save variance components');
    end

    function exportSplitsTable(T,summ,defaultfn,label,dlgtitle)
        %Write a data-splits table to XLSX or CSV with a metadata header,
        %mirroring the Mac-safe writecell/writetable pattern used by the dynamic-
        %reliability viewer (psyrat_startview_dynrel) and the standard viewers.
        [fn,pth] = uiputfile(...
            {'*.xlsx','Excel file (*.xlsx)';'*.csv','CSV file (*.csv)'},...
            dlgtitle,defaultfn);
        if isequal(fn,0); return; end
        outfile = fullfile(pth,fn);
        [~,~,ext] = fileparts(outfile);

        %metadata header (one row-cell per line); writecell(head',...) writes each
        %entry down its own spreadsheet row. Two trailing blanks for spacing.
        head = {};
        head{end+1} = sprintf('PsyRAT data-splits reliability: %s',label);
        head{end+1} = sprintf('Generated: %s',datestr(clock));
        if isfield(psyrat_data,'ver')
            head{end+1} = sprintf('PsyRAT v%s',psyrat_data.ver);
        end
        head{end+1} = sprintf('Dataset: %s',REL.filename);
        head{end+1} = L.banner;
        if isfield(REL,'nchains') && isfield(REL,'niter')
            head{end+1} = psyrat_iterline(REL);
        end
        head{end+1} = sprintf('%s: %.3g',psyrat_interval_label(REL),summ.ci);
        head{end+1} = summ.gbias_note;
        if isfield(summ,'design_note')
            head{end+1} = summ.design_note;
            for ds = 1:numel(summ.strata)
                head{end+1} = sprintf('Design %s: %s',summ.strata(ds).label,...
                    jsonencode(summ.strata(ds).design));
            end
        end
        %run provenance (seed/engine/priors/convergence); fail-open
        head = [head, psyrat_provenance_lines(psyrat_data)'];
        head{end+1} = '';
        head{end+1} = '';

        if strcmpi(ext,'.xlsx')
            %a fresh, checked workbook replaces the chosen file (G110); a
            %refused save raises, so the "Saved" line below is not printed
            psyrat_replace_workbook(outfile,head',T,numel(head));
        else
            %CSV: write the header lines, then append the table (its own column
            %names + data) so arbitrary table columns are handled robustly.
            %'WriteMode','append' defaults WriteVariableNames to false, so force
            %it on to keep the column-name header row (matching the .xlsx branch
            %and the headless psyrat_write_table writer).
            fid = fopen(outfile,'w');
            %quote per the CSV convention (B19): a comma-bearing header line
            %must stay ONE spreadsheet field; comma-free lines stay
            %byte-identical. The XLSX branch (writecell) already quotes.
            for i = 1:numel(head)
                fprintf(fid,'%s\n',psyrat_quote_csv_header(head{i}));
            end
            fclose(fid);
            writetable(T,outfile,'WriteMode','append','WriteVariableNames',true);
        end
        fprintf('Saved %s to %s\n',label,outfile);
    end

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
