function overalltable = psyrat_trt_reloverallt(varargin)
%Display a table with the overall reliability information for data after
% applying the cutoff
%
%psyrat_trt_reloverallt('psyrat_data',psyrat_data,'gui',1);
%
%
%Inputs
% psyrat_data - PsyRAT Toolbox data structure array.
% gui - 0 for off, 1 for on
%
%Outputs
% overalltable - table displaying reliability information
% a gui will also be shown if desired

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
            'See help psyrat_dep for more information about inputs'));
    end
    
    %check if psyrat_data was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('psyrat_data',varargin),1);
    if ~isempty(ind)
        psyrat_data = varargin{ind+1};
    else
        error('varargin:psyrat_data',... %Error code and associated error
            strcat('WARNING: psyrat_data not specified \n\n',...
            'Please input psyrat_data (PsyRAT Toolbox data structure array).\n',...
            'See help psyrat_depvtrialsplot for more information \n'));
    end
    
    %check if gui was specified
    %If it is not found, set display error.
    ind = find(strcmpi('gui',varargin),1);
    if ~isempty(ind)
        gui = varargin{ind+1};
    else
        error('varargin:gui',... %Error code and associated error
            strcat('WARNING: gui not specified \n\n',...
            'Please input gui specifying whether to display a gui.\n',...
            '0 for off, 1 for on\n',...
            'See help psyrat_depvtrialsplot for more information \n'));
    end
    
end

%check whether any groups exist
if strcmpi(psyrat_data.rel.groups,'none')
    ngroups = 1;
    %gnames = cellstr(psyrat_data.rel.groups);
    gnames ={''};
else
    ngroups = length(psyrat_data.rel.groups);
    gnames = psyrat_data.rel.groups(:);
end

%check whether any events exist
if strcmpi(psyrat_data.rel.events,'none')
    nevents = 1;
    %enames = cellstr(psyrat_data.rel.events);
    enames = {''};
else
    nevents = length(psyrat_data.rel.events);
    enames = psyrat_data.rel.events(:);
end

%figure out whether groups or events need to be considered
%1 - no groups or event types to consider
%2 - possible multiple groups but no event types to consider
%3 - possible event types but no groups to consider
%4 - possible groups and event types to consider

if ngroups == 1 && nevents == 1
    analysis = 1;
elseif ngroups > 1 && nevents == 1
    analysis = 2;
elseif ngroups == 1 && nevents > 1
    analysis = 3;
elseif ngroups > 1 && nevents > 1
    analysis = 4;
end

%create placeholders for displaying data in tables in guis
label = {};
overallrel = {};
mintrl = {};
maxtrl = {};
meantrl = {};
medtrl = {};
stdtrl = {};
goodn = {};
badn = {};

%put data together to display in tables
for gloc=1:ngroups
    for eloc=1:nevents
        
        %label for group and/or event
        switch analysis
            case 1
                label{end+1} = 'Measurement';
            case 2
                label{end+1} = gnames{gloc};
            case 3
                label{end+1} = enames{eloc};
            case 4
                label{end+1} = [gnames{gloc} ' - ' enames{eloc}];
        end
        
        %create a string with the reliability point estimate and credible
        %interval for overall data
        overallrel{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            psyrat_data.relsummary.group(gloc).event(eloc).rel.m,...
            psyrat_data.relsummary.group(gloc).event(eloc).rel.ll,...
            psyrat_data.relsummary.group(gloc).event(eloc).rel.ul);
        
        %put together trial summary information
        mintrl{end+1} = psyrat_data.relsummary.group(gloc).event(eloc).trlinfo.min;
        maxtrl{end+1} = psyrat_data.relsummary.group(gloc).event(eloc).trlinfo.max;
        meantrl{end+1} = psyrat_data.relsummary.group(gloc).event(eloc).trlinfo.mean;
        medtrl{end+1} = psyrat_data.relsummary.group(gloc).event(eloc).trlinfo.med;
        stdtrl{end+1} = psyrat_data.relsummary.group(gloc).event(eloc).trlinfo.std;
        
        %pull good and bad ns
        goodn{end+1} = psyrat_data.relsummary.group(gloc).event(eloc).goodn;
        badn{end+1} = length(psyrat_data.relsummary.group(gloc).badids);
        
    end
end


%create table to describe the data including all trials
overalltable = table(label',goodn',badn',overallrel',meantrl',...
    medtrl',stdtrl',mintrl',maxtrl');

switch psyrat_data.relsummary.gcoeff_name
    case 'dep'
        rel_name = 'Dependability';
    case 'gen'
        rel_name = 'Generalizability';
end

%analysis-unit label fragments (trial vs split); a parallel-splits run reports
%splits, not trials. psyrat_unitlabels defaults to "trial" when no splits flag is
%present, so non-splits output is byte-identical.
usplit = 1;
if isfield(psyrat_data.rel,'splits'); usplit = psyrat_data.rel.splits; end
L = psyrat_unitlabels(usplit);

overalltable.Properties.VariableNames = {'Label', ...
    'n_Included','n_Excluded', ...
    rel_name, ['Mean_Num_' L.Units], ['Med_Num_' L.Units],...
    ['Std_Num_' L.Units],['Min_Num_' L.Units],...
    ['Max_Num_' L.Units]};

%display gui if desired
if gui == 1
    
    %define parameters for figure size
    figwidth = 815;
    figheight = 500;
    
    %define space between rows and first row location
    rowspace = 25;
    row = figheight - rowspace*2;
    
    %occasion-estimand label distinguishing single-occasion (n'_o = 1) from a
    %multi-occasion composite (n'_o = k) test-retest score (audit sec.14)
    occlabel = psyrat_trt_occlabel(psyrat_data);

    %create a gui for displaying the overall trial information
    psyrat_overall= psyrat_newfigure('unit','pix',...
        'position',[1150 150 figwidth figheight],...
        'menub','no',...
        'Tag','psyrat_output',...
        'name',[rel_name ' Analyses Including All ' L.Units ' - ' occlabel],...
        'numbertitle','off',...
        'resize','off');

    %Print the occasion-estimand label above the section title
    uicontrol(psyrat_overall,'Style','text','fontsize',11,...
        'HorizontalAlignment','center',...
        'String',occlabel,...
        'Position',[0 row+22 figwidth 22]);

    %Print the name of the loaded dataset
    uicontrol(psyrat_overall,'Style','text','fontsize',16,...
        'HorizontalAlignment','center',...
        'String',sprintf('Overall %s',rel_name),...
        'Position',[0 row figwidth 25]);
    
    %Start a table
    t = uitable('Parent',psyrat_overall,'Position',...
        [25 100 figwidth-50 figheight-175],...
        'Data',table2cell(overalltable));
    set(t,'ColumnName',{'Label' 'n Included' 'n Excluded' ...
        rel_name ['Mean ' L.NumSymbol] ['Med ' L.NumSymbol]...
        ['Std Dev of ' L.Units] ['Min ' L.NumSymbol] ['Max ' L.NumSymbol]});
    set(t,'ColumnWidth',{'auto' 'auto' 'auto' 110 'auto' 'auto' 'auto' 'auto'});
    set(t,'RowName',[]);
    
    %Create a save button that will take save the table
    uicontrol(psyrat_overall,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','Save Table',...
        'Position', [.5*figwidth/8 25 figwidth/4 50],...
        'Callback',{@psyrat_saveoveralltable,psyrat_data,overalltable});
    
    %Create button that will save good/bad ids
    uicontrol(psyrat_overall,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','Save IDs',...
        'Position', [3*figwidth/8 25 figwidth/4 50],...
        'Callback',{@psyrat_saveids,psyrat_data});
    
    str = sprintf(['Use this button to estimate the reliability coefficient\n'...
        'for a given ' L.numberof '. This is helpful when adequate reliability\n'...
        'is not reached but you want to estimate obtained reliability']);
    
    
    if any(strcmp(psyrat_data.relsummary.reltype_name,{'trt','ic_trt'}))
        %Create button that estimates a new reliability coefficient for the
        %current TRT reliability configuration.
        uicontrol(psyrat_overall,'Style','push','fontsize',14,...
            'HorizontalAlignment','center',...
            'String',{'Estimate New';'Reliability Coefficient'},...
            'Tooltip', str,...
            'Position', [5.5*figwidth/8 25 figwidth/4 50],...
            'Callback',{@psyrat_newrel,psyrat_data});
        
    end
end

end



function psyrat_saveoveralltable(varargin)
%if the button to save the overall trial information table was pressed

%parse inputs
psyrat_data = varargin{3};
overalltable = varargin{4};

%analysis-unit label fragments (trial vs split) for the CSV header; defaults to
%"trial" so non-splits output is byte-identical.
usplit = 1;
if isfield(psyrat_data.rel,'splits'); usplit = psyrat_data.rel.splits; end
L = psyrat_unitlabels(usplit);

%ask the user where the file should be saved. XLSX now writes cross-platform
%via writecell/writetable, so both formats are offered on every platform
%(including macOS).
[savename, savepath] = uiputfile(...
    {'*.xlsx','Excel File (.xlsx)';'*.csv',...
    'Comma-Separated Value File (.csv)'},...
    'Where would you like to save table?');

[~,~,ext] = fileparts(fullfile(savepath,savename));

switch psyrat_data.relsummary.gcoeff_name
    case 'dep'
        rel_name = 'Dependability';
    case 'gen'
        rel_name = 'Generalizability';
end

%occasion-estimand label distinguishing single-occasion (n'_o = 1) from a
%multi-occasion composite (n'_o = k) test-retest score (audit sec.14)
occlabel = psyrat_trt_occlabel(psyrat_data);

%save either an excel or csv file
if strcmp(ext,'.xlsx')

    %print header information about the dataset
    %Use the resolved coefficient name: rel_name is already 'Generalizability' on a
    %gcoeff=2 run, and hard-coding "Dependability" here made line 1 of the export
    %contradict both the cutoff line and the coefficient column inside the same file.
    filehead = {[rel_name ' Table Generated on']; datestr(clock);''};
    filehead{end+1} = sprintf('PsyRAT Toolbox v%s',psyrat_data.ver);
    filehead{end+1} = occlabel;
    filehead{end+1} = '';
    filehead{end+1} = sprintf('Dataset: %s',psyrat_data.rel.filename);
    filehead{end+1} = sprintf('%s Cutoff: %0.4f',...
        rel_name,...
        psyrat_data.relsummary.relcutoff);
    filehead{end+1} = sprintf('Cutoff Threshold used the %s',...
        psyrat_data.relsummary.meascutoff);
    filehead{end+1} = psyrat_iterline(psyrat_data.rel);
    %run provenance (seed/engine/priors/convergence); fail-open
    filehead = [filehead; psyrat_provenance_lines(psyrat_data)];
    filehead{end+1}='';
    filehead{end+1}='';

    %write table
    writecell(filehead,fullfile(savepath,savename));
    writetable(overalltable,fullfile(savepath,savename),...
        'Range',strcat('A',num2str(length(filehead))));
    
elseif strcmp(ext,'.csv')
    
    %print header information about dataset
    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s',[rel_name ' Table Generated on ']);
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'PsyRAT Toolbox v%s\n',psyrat_data.ver);
    fprintf(fid,'%s\n',psyrat_quote_csv_header(occlabel));
    fprintf(fid,' \n');
    fprintf(fid,'%s\n',psyrat_quote_csv_header(...
        sprintf('Dataset: %s',psyrat_data.rel.filename)));
    fprintf(fid,'%s Cutoff: %0.4f\n',...
        rel_name,...
        psyrat_data.relsummary.relcutoff);
    fprintf(fid,'%s\n',psyrat_quote_csv_header(...
        sprintf('Cutoff Threshold used the %s',...
        psyrat_data.relsummary.meascutoff)));
    fprintf(fid,'%s\n',psyrat_quote_csv_header(...
        psyrat_iterline(psyrat_data.rel)));
    psyrat_fprintf_provenance(fid, psyrat_data);
    fprintf(fid,' \n');

    %'%s' (was '%s%s%s'): the bracket concatenation builds ONE char array,
    %so two of the three conversions had no data and printed nothing --
    %harmless by MATLAB's argument-exhaustion, but the depoverallt twin
    %uses the correct single conversion (adversarial wave, 2026-08-29)
    fprintf(fid,'%s', ['Label,N Included,N Excluded,'...
        rel_name...
        ',Mean Num of ' L.Units ',Med Num of ' L.Units ','...
        'Std Dev Num of ' L.Units ',Min Num of ' L.Units...
        ',Max Num of ' L.Units]);
    fprintf(fid,' \n');
    
    %write the table information.
    %G47 (latent hardening): unlike the depoverallt twin of this loop, no
    %reachable route feeds this writer a '---' cell today -- trt_diff
    %overall tables go through psyrat_depoverallt (psyrat_relfigures /
    %psyrat_report both route them there), and this builder has no dash
    %block -- but the loop was byte-for-byte the same trap: a char cell
    %pushed through %d prints each character's numeric code and recycles
    %the format, garbling one row across multiple lines. Same fix as
    %there: render every cell to a string first (numeric cells through the
    %exact former numeric conversions, so all current rows stay
    %byte-identical), then join and write the row through ONE %s conversion
    %-- a single pre-joined line cannot lose fields to fprintf's
    %argument-flattening, which an empty char cell would otherwise trigger.
    %Column 6 (median trial count) prints with %g, not %d: a half-integer
    %median (an even trial count) rendered as %e under %d (the one-facet
    %writer's P4 observation, 2026-09-04, applies here identically). %g
    %prints an integer median exactly as %d did.
    numspec = {'%s','%d','%d','%s','%0.2f','%g','%0.2f','%d','%d'};
    for i = 1:height(overalltable)
        rowstr = cell(1,9);
        for j = 1:9
            val = overalltable{i,j};
            if iscell(val); val = val{1}; end
            if ischar(val)
                %Label and coefficient strings, plus any '---' placeholder
                rowstr{j} = val;
            else
                rowstr{j} = sprintf(numspec{j},val);
            end
        end
        fprintf(fid,'%s\n',strjoin(rowstr,','));
    end
    
    fclose(fid);
    
end

end

function psyrat_saveids(varargin)
%if the user pressed the button to save which ids were considered good and
%which were considered bad (based on whether data met the cutoff thresholds

psyrat_data = varargin{3};

%ask the user where the file should be saved. XLSX now writes cross-platform
%via writecell, so both formats are offered on every platform (including macOS).
[savename, savepath] = uiputfile(...
    {'*.xlsx','Excel File (.xlsx)';'*.csv',...
    'Comma-Separated Value File (.csv)'},...
    'Where would you like to save the data?');

[~,~,ext] = fileparts(fullfile(savepath,savename));

switch psyrat_data.relsummary.gcoeff_name
    case 'dep'
        rel_name = 'Dependability';
    case 'gen'
        rel_name = 'Generalizability';
end

%save the information in either an excel or csv format
if strcmp(ext,'.xlsx')
    
    datap{1,1} = 'Data Generated on';
    datap{end+1,1} = datestr(clock);
    datap{end+1,1} = sprintf('PsyRAT Toolbox v%s',psyrat_data.ver);
    datap{end+1,1} = '';
    datap{end+1,1} = sprintf('Dataset: %s',psyrat_data.rel.filename);
    datap{end+1,1} = sprintf('%s Cutoff: %0.4f',...
        rel_name,...
        psyrat_data.relsummary.relcutoff);
    datap{end+1,1} = sprintf('Cutoff Threshold used the %s',...
        psyrat_data.relsummary.meascutoff);
    datap{end+1,1} = psyrat_iterline(psyrat_data.rel);
    datap{end+1,1}='';
    datap{end+1,1}='';
    datap{end+1,1} = 'Good IDs'; datap{end,2} = 'Bad IDs';
    
    gids = [];
    bids = [];
    srow = length(datap);
    
    for j=1:length(psyrat_data.relsummary.group)
        gids = [gids;psyrat_data.relsummary.group(j).goodids(:)];
        bids = [bids;psyrat_data.relsummary.group(j).badids(:)];
    end
    
    for i = 1:length(gids)
        datap{i+srow,1}=char(gids(i));
    end
    
    for i = 1:length(bids)
        datap{i+srow,2}=char(bids(i));
    end
    
    writecell(datap,fullfile(savepath,savename));
    
elseif strcmp(ext,'.csv')
    
    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s\n','Data Generated on');
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'PsyRAT Toolbox v%s\n',psyrat_data.ver);
    fprintf(fid,' \n');
    fprintf(fid,'%s\n',psyrat_quote_csv_header(...
        sprintf('Dataset: %s',psyrat_data.rel.filename)));
    fprintf(fid,'%s Cutoff: %0.4f\n',...
        rel_name,...
        psyrat_data.relsummary.relcutoff);
    fprintf(fid,'%s\n',psyrat_quote_csv_header(...
        sprintf('Cutoff Threshold used the %s',...
        psyrat_data.relsummary.meascutoff)));
    fprintf(fid,'%s\n',psyrat_quote_csv_header(...
        psyrat_iterline(psyrat_data.rel)));
    fprintf(fid,' \n');
    fprintf(fid,'%s\n','Good IDs,Bad IDs');
    
    gids = [];
    bids = [];
    
    for j=1:length(psyrat_data.relsummary.group)
        gids = [gids;psyrat_data.relsummary.group(j).goodids(:)];
        bids = [bids;psyrat_data.relsummary.group(j).badids(:)];
    end
    
    maxlength = max([length(gids) length(bids)]);
    minlength = min([length(gids) length(bids)]);
    
    if minlength == length(gids)
        whichlonger = 1;
    elseif maxlength == length(gids)
        whichlonger = 2;
    end
    
    for i = 1:maxlength
        if i <= minlength
            fprintf(fid,'%s,%s\n',gids{i},bids{i});
        elseif i > minlength && whichlonger == 1
            fprintf(fid,',%s\n',bids{i});
        elseif i > minlength && whichlonger == 2
            fprintf(fid,'%s,\n',gids{i});
        end
    end
    
    fclose(fid);
    
end


end



function psyrat_newrel(varargin)
%if the user pressed the button to calculate a new test-retest reliability
%coefficient, this gui will be pulled up

psyrat_data = varargin{3};

%analysis-unit label fragments (trial vs split) for the recalculation prompts;
%defaults to "trial" so non-splits output is byte-identical.
usplit = 1;
if isfield(psyrat_data.rel,'splits'); usplit = psyrat_data.rel.splits; end
L = psyrat_unitlabels(usplit);

%define parameters for figure position
figwidth = 550;
figheight = 250;

%define space between rows and first row location
rowspace = 40;
row = figheight - rowspace*2;

%define locations of column 1 and 2
lcol = 30;
rcol = (figwidth/8)*5;

%specify font size
fsize = 14;

%create the gui
psyrat_newrelgui = psyrat_newfigure('unit','pix',...
    'position',[400 400 figwidth figheight],...
    'menub','no',...
    'Tag','psyrat_output',...
    'name','Specify Inputs',...
    'numbertitle','off',...
    'resize','off');

%Print the name of the loaded dataset
uicontrol(psyrat_newrelgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String',['Dataset:  ' psyrat_data.rel.filename],...
    'Tooltip','Dataset that was used',...
    'Position',[0 row figwidth 25]);

%next row
row = row - (rowspace*.45);

%Print the name of the measurement analyzed
uicontrol(psyrat_newrelgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String',['Measurement:  ' psyrat_data.proc.measheader],...
    'Tooltip','Dataset that was used',...
    'Position',[0 row figwidth 25]);

%next row
row = row - (rowspace*1.4);

str = sprintf(['Number of ' L.units ' to use for recalculating reliability\n',...
    'This should be the mean or median ' L.numberof ' retained for averaging']);

%Print the text for reliability cutoff with a box for the user to specify
%the input
uicontrol(psyrat_newrelgui,...
    'Style','text',...
    'fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String',[L.NumberOf ':'],...
    'Tooltip',str,...
    'Position', [lcol row figwidth/4 25]);

inputs.trls = uicontrol(psyrat_newrelgui,...
    'Style','edit',...
    'fontsize',fsize,...
    'String','',...
    'Position', [rcol+5 row+6 figwidth/4 25]);

%next row
row = row - rowspace*2;

str = sprintf(['Use this button to estimate the new test-retest reliability\n'...
    'information for the specified ' L.numberof '.']);

%Create button that will save good/bad ids
uicontrol(psyrat_newrelgui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String',{'Estimate New';'Reliability Coefficient'},...
    'Tooltip', str,...
    'Position', [.5*figwidth/8 row figwidth/1.11 50],...
    'Callback',{@psyrat_shownewrel,psyrat_data,inputs.trls});

end

function psyrat_shownewrel(varargin)
%estimate the new reliability estimates and show it

psyrat_data = varargin{3};

%analysis-unit label fragments (trial vs split) for the recalculated table;
%defaults to "trial" so non-splits output is byte-identical.
usplit = 1;
if isfield(psyrat_data.rel,'splits'); usplit = psyrat_data.rel.splits; end
L = psyrat_unitlabels(usplit);

%COMPLEX-LITERAL GUARD, and the field's first validation of any kind.
%
%ntrls is n', the size of the composite the D-study projects to, and it reaches
%psyrat_rel_trt as 'obs', where it DIVIDES the error term in all six
%gcoeff x reltype formulas. Its admissible domain is therefore the positive
%integers: a count of measurement replications entering a mean cannot be zero,
%negative, or fractional. Before this guard nothing checked it and this file has
%no try/catch, so each bad entry produced a plausible-looking answer instead of
%a message -- blank gave "NaN CI [NaN NaN]" under a title reading "for NaN
%Trials"; a negative value pushed the coefficient outside [0,1]; and a complex
%literal crashed later, in quantile().
%
%Zero is worth stating precisely, because it does NOT behave uniformly. For
%reltype 1 and 3 only the error term carries obs, so obs = 0 makes it Inf and
%the coefficient prints as exactly "0.00 CI [0.00 0.00]" -- a publishable-
%looking number. For reltype 2 (coefficient of stability) obs also divides the
%UNIVERSE score (psyrat_rel_trt.m:271, :310), so it is Inf/(Inf+Inf) = NaN
%instead. Both are wrong; only the first is quietly plausible.
%
%~isreal LEADS because isnan is false for a complex and the ordering test reads
%the real part alone. mod() also rejects +Inf (mod(Inf,1) is NaN), and -Inf is
%caught earlier by ntrls < 1, so neither is covered separately.
%
%Floor is 1, not the <= 1 used by the vs-trials PLOT fields in
%psyrat_startview_sing/trt. Those reject 1 because their value becomes a plot
%RANGE endpoint (x = trials(1):trials(2)), where a single point is a degenerate
%line -- not because the colon expression is invalid, since 1:1 is perfectly
%legal. Here obs is a scalar and a single-trial composite is a legitimate, if
%uninformative, D-study target. Note those siblings also carry no integrality
%test, so this guard is their shape plus a mod() term, not a copy.
%
%The field starts blank and has no defensible default, so it is REQUIRED --
%hence errordlg + refocus + return, the house contract for required fields,
%rather than the silent fallback used by optional overrides such as
%psyrat_get_nocc.
ntrls = str2double(varargin{4}.String);
if ~isreal(ntrls) || isnan(ntrls) || ntrls < 1 || mod(ntrls,1) ~= 0
    %B21 stage 1 (2026-08-16): modal + purpose title. The single-argument
    %form was non-modal and untitled, opening at a fixed position the
    %Specify Inputs figure covers at default placement -- the user saw
    %NOTHING on a rejected entry (observed live twice). Modal raises it,
    %and the same-title replacement keeps repeated rejections from
    %stacking into an ignoring-looking pile (the B22 half).
    errordlg(['Please specify a whole number of ' L.units ' (1 or more) ' ...
        'to estimate a new reliability coefficient'],...
        'Invalid trial count','modal');
    if ishghandle(varargin{4})
        uicontrol(varargin{4});
    end
    return;
end

%use the SAME effective number of occasions (n'_o) that the overall
%coefficient used so the recalculated estimate stays on the selected
%single- vs multi-occasion estimand (audit sec.14). Older result structs
%without the setting fall back to the single-occasion default of 1.
if isfield(psyrat_data.relsummary,'nocc') && ...
        ~isempty(psyrat_data.relsummary.nocc)
    nocc = psyrat_data.relsummary.nocc;
else
    nocc = 1;
end

%place psyrat_data.data in REL to work with
REL = psyrat_data.rel;

%check whether any groups exist
if strcmpi(REL.groups,'none')
    ngroups = 1;
    %gnames = cellstr(REL.groups);
    gnames ={''};
else
    ngroups = length(REL.groups);
    gnames = REL.groups(:);
end

%check whether any events exist
if strcmpi(REL.events,'none')
    nevents = 1;
    %enames = cellstr(REL.events);
    enames = {''};
else
    nevents = length(REL.events);
    enames = REL.events(:);
end

%figure out whether groups or events need to be considered
%1 - no groups or event types to consider
%2 - possible multiple groups but no event types to consider
%3 - possible event types but no groups to consider
%4 - possible groups and event types to consider


if ngroups == 1 && nevents == 1
    analysis = 1;
elseif ngroups > 1 && nevents == 1
    analysis = 2;
elseif ngroups == 1 && nevents > 1
    analysis = 3;
elseif ngroups > 1 && nevents > 1
    analysis = 4;
end

%extract information from REL and store in data for crunching
switch analysis
    case 1 %1 - no groups or event types to consider
        
        gloc = 1;
        eloc = 1;
        
        data.g(gloc).e(eloc).label = REL.out.labels{gloc};
        data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,gloc);
        data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,gloc);
        data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,gloc);
        data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,gloc);
        data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,gloc);
        data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,gloc);
        data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,gloc);
        data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,gloc);
        data.g(gloc).e(eloc).elabel = cellstr('none');
        data.g(gloc).glabel = gnames(gloc);
        
    case 2 %2 - possible multiple groups but no event types to consider
        
        eloc = 1;
        
        for gloc=1:length(REL.out.labels)
            
            data.g(gloc).e(eloc).label = REL.out.labels{gloc};
            data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,gloc);
            data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,gloc);
            data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,gloc);
            data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,gloc);
            data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,gloc);
            data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,gloc);
            data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,gloc);
            data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,gloc);
            data.g(gloc).e(eloc).elabel = cellstr('none');
            data.g(gloc).glabel = gnames(gloc);
            
        end
        
    case 3 %3 - possible event types but no groups to consider
        
        gloc = 1;
        
        for eloc=1:length(REL.out.labels)
            
            data.g(gloc).e(eloc).label = REL.out.labels{eloc};
            data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,eloc);
            data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,eloc);
            data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,eloc);
            data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,eloc);
            data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,eloc);
            data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,eloc);
            data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,eloc);
            data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,eloc);
            data.g(gloc).e(eloc).elabel = enames(eloc);
            data.g(gloc).glabel = gnames(gloc);
            
        end
        
    case 4 %4 - possible groups and event types to consider
        for ii=1:length(REL.out.labels)
            
            %use the underscores that were added in psyrat_computevarcomp to
            %differentiate where the group and event the data are for
            lblstr = strsplit(REL.out.labels{ii},'_;_');
            
            eloc = find(ismember(enames,lblstr(2)));
            gloc = find(ismember(gnames,lblstr(1)));
            
            data.g(gloc).e(eloc).label = REL.out.labels{ii};
            data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,ii);
            data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,ii);
            data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,ii);
            data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,ii);
            data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,ii);
            data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,ii);
            data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,ii);
            data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,ii);
            data.g(gloc).e(eloc).elabel = enames(eloc);
            data.g(gloc).glabel = gnames(gloc);
            
        end
end %switch analysis

%compute reliability data for each group and event
switch analysis
    case 1 %no groups or event types to consider
        
        %the same generic structure is used for relsummary, so the event
        %and group locations will both be 1 for the data
        eloc = 1;
        gloc = 1;
        
        %compute reliabiltiy
        [llrel,mrel,ulrel] = psyrat_rel_trt(...
            'gcoeff',psyrat_data.relsummary.gcoeff,...
            'reltype',psyrat_data.relsummary.reltype,...
            'bp',data.g(gloc).e(eloc).sig_id.raw,...
            'bo',data.g(gloc).e(eloc).sig_occ.raw,...
            'bt',data.g(gloc).e(eloc).sig_trl.raw,...
            'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
            'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
            'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
            'err',data.g(gloc).e(eloc).sig_err.raw,...
            'obs', ntrls,'nocc',nocc,'CI',psyrat_data.relsummary.ciperc);
        
        newrelsummary.group(gloc).event(eloc).rel.m = mrel;
        newrelsummary.group(gloc).event(eloc).rel.ll = llrel;
        newrelsummary.group(gloc).event(eloc).rel.ul = ulrel;
        
    case 2 %possible multiple groups but no event types to consider
        
        
        %since the same generic structure is used for relsummary, the event
        %location will be defined as 1.
        eloc = 1;
        
        for gloc=1:ngroups %loop through each group
            
            %compute reliabiltiy
            [llrel,mrel,ulrel] = psyrat_rel_trt(...
                'gcoeff',psyrat_data.relsummary.gcoeff,...
                'reltype',psyrat_data.relsummary.reltype,...
                'bp',data.g(gloc).e(eloc).sig_id.raw,...
                'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                'err',data.g(gloc).e(eloc).sig_err.raw,...
                'obs',ntrls,'nocc',nocc,'CI',psyrat_data.relsummary.ciperc);
            
            newrelsummary.group(gloc).event(eloc).rel.m = mrel;
            newrelsummary.group(gloc).event(eloc).rel.ll = llrel;
            newrelsummary.group(gloc).event(eloc).rel.ul = ulrel;
            
        end
        
    case 3 %possible event types but no groups to consider
        
        %since the relsummary structure is generic for any number of groups
        %or events, the data will count as being from 1 group.
        gloc = 1;
        
        %cylce through each event
        for eloc=1:nevents
            
            %compute reliabiltiy
            [llrel,mrel,ulrel] = psyrat_rel_trt(...
                'gcoeff',psyrat_data.relsummary.gcoeff,...
                'reltype',psyrat_data.relsummary.reltype,...
                'bp',data.g(gloc).e(eloc).sig_id.raw,...
                'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                'err',data.g(gloc).e(eloc).sig_err.raw,...
                'obs',ntrls,'nocc',nocc,'CI',psyrat_data.relsummary.ciperc);
            
            newrelsummary.group(gloc).event(eloc).rel.m = mrel;
            newrelsummary.group(gloc).event(eloc).rel.ll = llrel;
            newrelsummary.group(gloc).event(eloc).rel.ul = ulrel;
            
        end
        
    case 4 %groups and event types to consider
        
        %cycle through each event
        for eloc=1:nevents
            
            %cycle through each group
            for gloc=1:ngroups
                
                %compute reliability
                [llrel,mrel,ulrel] = psyrat_rel_trt(...
                    'gcoeff',psyrat_data.relsummary.gcoeff,...
                    'reltype',psyrat_data.relsummary.reltype,...
                    'bp',data.g(gloc).e(eloc).sig_id.raw,...
                    'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                    'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                    'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                    'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                    'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                    'err',data.g(gloc).e(eloc).sig_err.raw,...
                    'obs',ntrls,'nocc',nocc,'CI',psyrat_data.relsummary.ciperc);
                
                newrelsummary.group(gloc).event(eloc).rel.m = mrel;
                newrelsummary.group(gloc).event(eloc).rel.ll = llrel;
                newrelsummary.group(gloc).event(eloc).rel.ul = ulrel;
                
            end
        end
        
        
end %switch analysis


%create placeholders for displaying data in tables in guis
label = {};
overallrel = {};
new_trls = {};

%put data together to display in tables
for gloc=1:ngroups
    for eloc=1:nevents
        
        %label for group and/or event
        switch analysis
            case 1
                label{end+1} = 'Measurement';
            case 2
                label{end+1} = gnames{gloc};
            case 3
                label{end+1} = enames{eloc};
            case 4
                label{end+1} = [gnames{gloc} ' - ' enames{eloc}];
        end
        
        %create a string with the reliability point estimate and credible
        %interval for overall data
        overallrel{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            newrelsummary.group(gloc).event(eloc).rel.m,...
            newrelsummary.group(gloc).event(eloc).rel.ll,...
            newrelsummary.group(gloc).event(eloc).rel.ul);
        
        new_trls{end+1} = ntrls;
        
    end
end


%create table to describe the data including all trials
overalltable = table(label',overallrel',new_trls');

switch psyrat_data.relsummary.gcoeff_name
    case 'dep'
        rel_name = 'Dependability';
    case 'gen'
        rel_name = 'Generalizability';
end

overalltable.Properties.VariableNames = {'Label', ...
    rel_name, ['Num_' L.Units]};

%define parameters for figure size
figwidth = 415;
figheight = 350;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

%create a gui for displaying the overall trial information
psyrat_overall= psyrat_newfigure('unit','pix',...
    'position',[1150 150 figwidth figheight],...
    'menub','no',...
    'Tag','psyrat_output',...
    'name',sprintf('%s Analyses for %d %s',...
    rel_name,ntrls,L.Units),...
    'numbertitle','off',...
    'resize','off');

%Print the name of the loaded dataset
uicontrol(psyrat_overall,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',sprintf('Overall %s',rel_name),...
    'Position',[0 row figwidth 25]);

%Start a table
t = uitable('Parent',psyrat_overall,'Position',...
    [25 100 figwidth-50 figheight-175],...
    'Data',table2cell(overalltable));
set(t,'ColumnName',{'Label' ...
    rel_name L.NumSymbol});
set(t,'ColumnWidth',{151 110 100});
set(t,'RowName',[]);


end
