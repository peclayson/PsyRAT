function overalltable = psyrat_depoverallt(varargin)
%Display a table with the overall dependability information for data after
%applying the cutoff
%
%psyrat_depoverallt('psyrat_data',psyrat_data,'gui',1);
%
%
%Inputs
% psyrat_data - PsyRAT Toolbox data structure array. 
% gui - 0 for off, 1 for on
%
%Outputs
% overalltable - table displaying dependability information
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
overalldep = {};
mintrl = {};
maxtrl = {};
meantrl = {};
medtrl = {};
stdtrl = {};
goodn = {};
badn = {};

if isfield(psyrat_data,'relsummary') && isfield(psyrat_data.relsummary,'gcoeff_name') && ...
        strcmpi(psyrat_data.relsummary.gcoeff_name,'gen')
    rel_name = 'Generalizability';
else
    rel_name = 'Dependability';
end

if ~any(strcmp(psyrat_data.rel.analysis,{'ic_sserrvar','ic_diff_sserrvar'}))
    %put data together to display in tables
    for gloc=1:ngroups
        for eloc=1:nevents
            
            %label for group and/or event
            switch analysis
                case 1
                    label{end+1} = 'Measurement'; %#ok<*AGROW>
                case 2
                    label{end+1} = gnames{gloc};
                case 3
                    label{end+1} = enames{eloc};
                case 4
                    label{end+1} = [gnames{gloc} ' - ' enames{eloc}];
            end
            
            %create a string with the dependability point estimate and credible
            %interval for overall data
            overalldep{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                psyrat_data.relsummary.group(gloc).event(eloc).dep.m,...
                psyrat_data.relsummary.group(gloc).event(eloc).dep.ll,...
                psyrat_data.relsummary.group(gloc).event(eloc).dep.ul);
            
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
        
        if any(strcmp(psyrat_data.rel.analysis,{'ic_diff','trt_diff'}))

        switch analysis
            case 3
                label{end+1} = 'diff score';
            case 4
                label{end+1} = [gnames{gloc} ' - diff score'];
        end
        
        mintrl{end+1} = '---';
        maxtrl{end+1} = '---';
        meantrl{end+1} = '---';
        medtrl{end+1} = '---';
        stdtrl{end+1} = '---';
        
        %pull good and bad ns
        goodn{end+1} = psyrat_data.relsummary.group(gloc).event(eloc).goodn;
        badn{end+1} = '---';
        
        %create a string with the dependability point estimate and credible
        %interval for cutoff data
        overalldep{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            psyrat_data.relsummary.group(gloc).diffscore.pt,...
            psyrat_data.relsummary.group(gloc).diffscore.ll,...
            psyrat_data.relsummary.group(gloc).diffscore.ul);
    end
    end
    
    
    
elseif any(strcmp(psyrat_data.rel.analysis,{'ic_sserrvar','ic_diff_sserrvar'}))
    for gloc=1:ngroups
        include_event_tables = true;
        if strcmp(psyrat_data.rel.analysis,'ic_diff_sserrvar')
            include_event_tables = isfield(psyrat_data.relsummary.group(gloc),'event') && ...
                isfield(psyrat_data.relsummary.group(gloc).event,'ssrel_table');
        end

        if include_event_tables
            for eloc=1:nevents

                idtable_raw = psyrat_data.relsummary.group(gloc).event(eloc).ssrel_table;

                goodid_rows = ismember(idtable_raw.id,...
                    psyrat_data.relsummary.group(gloc).goodids);

                idtable = idtable_raw(goodid_rows,:);

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

                if isempty(idtable)
                    %G44. A group can retain nobody -- the G39 strict
                    %intersection folds a dead event through the whole group's
                    %retained set -- and mean/std over zero rows would render
                    %' NaN SD: NaN' with empty min/max cells here. Dash the
                    %summary cells out instead, reusing the '---' MARKER the
                    %diff rows in the non-sserr branch already print. WHICH
                    %cells dash differs between the two by design, not by
                    %accident: the diff rows dash badn and keep their
                    %coefficient cell (it is a real draw-based value there),
                    %whereas here the coefficient summary is a mean over zero
                    %retained rows (undefined, dashed) and the counts below
                    %(0 included / N excluded) are real, informative values
                    %that stay numeric and complement.
                    overalldep{end+1} = '---';
                    mintrl{end+1} = '---';
                    maxtrl{end+1} = '---';
                    meantrl{end+1} = '---';
                    medtrl{end+1} = '---';
                    stdtrl{end+1} = '---';
                else
                    %create a string with the dependability point estimate and
                    %credible interval for overall data
                    overalldep{end+1} = sprintf(' %0.2f SD: %0.2f',...
                        mean(idtable.dep_pt),...
                        std(idtable.dep_pt));

                    %put together trial summary information
                    mintrl{end+1} = min(idtable.trls);
                    maxtrl{end+1} = max(idtable.trls);
                    meantrl{end+1} = mean(idtable.trls);
                    medtrl{end+1} = median(idtable.trls);
                    stdtrl{end+1} = std(idtable.trls);
                end

                %pull good and bad ns. goodn is the GROUP-carrier count -- the
                %ids the group retained that appear in this event's table --
                %NOT the sum of this event's own inclusion indicator (G53). On
                %ic_sserrvar the two are provably equal (group goodids are the
                %per-event include set, or their strict intersection, so every
                %retained id is event-included), but on ic_diff_sserrvar the
                %group carrier derives from the DIFFERENCE-score inclusion
                %while each event table is marked independently against the
                %same cutoff; a diff-included but event-excluded subject then
                %fell out of BOTH counts and the rendered row summed short of
                %this event's table. Counting the carrier makes n Included +
                %n Excluded complement to height(idtable_raw) as an ismember
                %identity (the good and bad row masks are the same test
                %negated) -- the event table's roster, which equals the
                %participant roster only when every subject has rows in this
                %event. goodn also matches what Save IDs writes for the good
                %list (per group the carrier count IS numel(goodids) here),
                %and counts exactly the rows the coefficient summary above
                %averages over. badn is NOT the badids export: on the diff
                %route badids come from the diff table, which omits a subject
                %present in only one event, while badn counts every
                %non-carrier row in this event's table.
                goodn{end+1} = height(idtable);

                badid_rows = ~ismember(idtable_raw.id,...
                    psyrat_data.relsummary.group(gloc).goodids);

                badn{end+1} = height(idtable_raw(badid_rows,:));

            end
        end

        if strcmp(psyrat_data.rel.analysis,'ic_diff_sserrvar') && ...
                isfield(psyrat_data.relsummary.group(gloc),'diffscore')

            switch analysis
                case 3
                    label{end+1} = 'diff score';
                case 4
                    label{end+1} = [gnames{gloc} ' - diff score'];
            end

            idtable_raw = psyrat_data.relsummary.group(gloc).diffscore.ssrel_table;
            idtable = idtable_raw(idtable_raw.ind2include,:);

            if isempty(idtable)
                %G57 (same semantics as the G44 arm above). A group can
                %retain nobody on the DIFFERENCE-score inclusion too, and
                %this row's coefficient cell is a mean over the per-subject
                %dep_pt values of the retained rows -- undefined over zero
                %rows, so it is dashed along with the five trial cells. The
                %counts below stay numeric (0 included / N excluded): they
                %come from the raw indicators and remain real, complementing
                %values. This branch used to fall back to the FULL table
                %(idtable = idtable_raw), which printed a coefficient and
                %trial statistics computed over every EXCLUDED subject next
                %to n Included = 0 -- or ' NaN SD: NaN' when the diff draws
                %were all NaN. Unlike the non-sserr diff rows, whose
                %coefficient is a draw-based value that stays numeric by
                %design, this one has no defined value at zero retention.
                overalldep{end+1} = '---';
                mintrl{end+1} = '---';
                maxtrl{end+1} = '---';
                meantrl{end+1} = '---';
                medtrl{end+1} = '---';
                stdtrl{end+1} = '---';
            else
                overalldep{end+1} = sprintf(' %0.2f SD: %0.2f',...
                    mean(idtable.dep_pt),...
                    std(idtable.dep_pt));

                mintrl{end+1} = min(idtable.trls);
                maxtrl{end+1} = max(idtable.trls);
                meantrl{end+1} = mean(idtable.trls);
                medtrl{end+1} = median(idtable.trls);
                stdtrl{end+1} = std(idtable.trls);
            end
            goodn{end+1} = sum(idtable_raw.ind2include);
            badn{end+1} = sum(idtable_raw.ind2exclude);
        end
    end
    
    
end
%analysis-unit label fragments (trial vs split). A parallel-splits run averages
%items within splits, so the per-cell counts are splits, not trials;
%psyrat_unitlabels defaults to "trial" when no splits flag is present, so
%non-splits output is byte-identical.
usplit = 1;
if isfield(psyrat_data.rel,'splits'); usplit = psyrat_data.rel.splits; end
L = psyrat_unitlabels(usplit);

%create table to describe the data including all trials
overalltable = table(label',goodn',badn',overalldep',meantrl',...
    medtrl',stdtrl',mintrl',maxtrl');

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
    
    %create a gui for displaying the overall trial information
    psyrat_overall= psyrat_newfigure('unit','pix',...
        'position',[1150 150 figwidth figheight],...
        'menub','no',...
        'Tag','psyrat_output',...
        'name',[rel_name ' Analyses'],...
        'numbertitle','off',...
        'resize','off');
    
    %Print the name of the loaded dataset
    uicontrol(psyrat_overall,'Style','text','fontsize',16,...
        'HorizontalAlignment','center',...
        'String',['Overall ' rel_name],...
        'Position',[0 row figwidth 25]);
    
    %Start a table
    t = uitable('Parent',psyrat_overall,'Position',...
        [25 100 figwidth-50 figheight-175],...
        'Data',table2cell(overalltable));
    
    if ~any(strcmp(psyrat_data.rel.analysis,{'ic_sserrvar','ic_diff_sserrvar'}))
        set(t,'ColumnName',{'Label' 'n Included' 'n Excluded' ...
            rel_name ['Mean ' L.NumSymbol] ['Med ' L.NumSymbol]...
            ['Std Dev of ' L.Units] ['Min ' L.NumSymbol] ['Max ' L.NumSymbol]});
    elseif any(strcmp(psyrat_data.rel.analysis,{'ic_sserrvar','ic_diff_sserrvar'}))
        set(t,'ColumnName',{'Label' 'n Included' 'n Excluded' ...
            ['Subject Level ' rel_name] ['Mean ' L.NumSymbol] ['Med ' L.NumSymbol]...
            ['Std Dev of ' L.Units] ['Min ' L.NumSymbol] ['Max ' L.NumSymbol]});
    end
    
    set(t,'ColumnWidth',{'auto' 'auto' 'auto' 110 'auto' 'auto' 'auto' 'auto'});
    set(t,'RowName',[]);
    
    %Create a save button that will take save the table
    uicontrol(psyrat_overall,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','Save Table',...
        'Position', [figwidth/8 25 figwidth/4 50],...
        'Callback',{@psyrat_saveoveralltable,psyrat_data,overalltable});
    
    %Create button that will save good/bad ids
    uicontrol(psyrat_overall,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','Save IDs',...
        'Position', [5*figwidth/8 25 figwidth/4 50],...
        'Callback',{@psyrat_saveids,psyrat_data});
    
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

if isfield(psyrat_data,'relsummary') && isfield(psyrat_data.relsummary,'gcoeff_name') && ...
        strcmpi(psyrat_data.relsummary.gcoeff_name,'gen')
    rel_name = 'Generalizability';
else
    rel_name = 'Dependability';
end

%ask the user where the file should be saved. XLSX now writes cross-platform
%via writecell/writetable, so both formats are offered on every platform
%(including macOS).
[savename, savepath] = uiputfile(...
    {'*.xlsx','Excel File (.xlsx)';'*.csv',...
    'Comma-Separated Value File (.csv)'},...
    'Where would you like to save table?');

[~,~,ext] = fileparts(fullfile(savepath,savename));

%save either an excel or csv file
if strcmp(ext,'.xlsx')
    
    %print header information about the dataset
    filehead = {[rel_name ' Table Generated on']; datestr(clock);''}; 
    filehead{end+1} = sprintf('PsyRAT Toolbox v%s',psyrat_data.ver);
    filehead{end+1} = '';
    filehead{end+1} = sprintf('Dataset: %s',psyrat_data.rel.filename);
    %B36: %0.4f, NOT %0.2f. The cutoff is a recorded analytic INPUT, not an
    %estimate. It is consumed at full precision -- psyrat_relsummary's
    %find(llrel >= depcutoff,1) turns it into the minimum trial count that
    %defines the retained sample -- and depcheck accepts any real in (0,1), so
    %a 3-4 dp cutoff is reachable (a sensitivity grid stepping .025 gives .725,
    %.775). At %0.2f a run at 0.9999 wrote "Cutoff: 1.00", a threshold that
    %cannot regenerate the trial cutoff printed elsewhere in the same file.
    %
    %Estimates in the table BODY stay at %0.2f (see the coefficient/CI formats
    %below and in the sibling exporters) -- that is the reporting convention and
    %is unchanged. This is the input/estimate distinction, not "2 dp is wrong".
    %
    %%0.4f rather than %g because the toolbox already writes its OTHER
    %user-typed threshold that way: psyrat_criterionfigures.m's
    %'Criterion Cutoff: %0.4f', pinned as 0.5000 in TestExportCallbacks. %g
    %would render the default as 0.8 and break that suite's sibling pin.
    %All 13 cutoff sites across the 5 exporters moved together; the check that
    %the sweep is complete is that `grep -rn "Cutoff: %0.2f" subroutines/`
    %returns nothing.
    filehead{end+1} = sprintf('%s Cutoff: %0.4f',rel_name,...
        psyrat_data.relsummary.depcutoff);
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
    fprintf(fid,' \n');
    fprintf(fid,'%s\n',psyrat_quote_csv_header(...
        sprintf('Dataset: %s',psyrat_data.rel.filename)));
    fprintf(fid,'%s Cutoff: %0.4f\n',rel_name,...
        psyrat_data.relsummary.depcutoff);
    fprintf(fid,'%s\n',psyrat_quote_csv_header(...
        sprintf('Cutoff Threshold used the %s',...
        psyrat_data.relsummary.meascutoff)));
    fprintf(fid,'%s\n',psyrat_quote_csv_header(...
        psyrat_iterline(psyrat_data.rel)));
    psyrat_fprintf_provenance(fid, psyrat_data);
    fprintf(fid,' \n');

    fprintf(fid,'%s', ['Label,N Included,N Excluded,'...
        rel_name ',Mean Num of ' L.Units ',Med Num of ' L.Units ','...
        'Std Dev Num of ' L.Units ',Min Num of ' L.Units...
        ',Max Num of ' L.Units]);
    fprintf(fid,' \n');
    
    %write the table information.
    %G47: the diff rows (ic_diff/trt_diff) carry '---' in the badn and
    %trial-summary cells, and a zeroed sserr group's row (G44) dashes its
    %COEFFICIENT and trial-summary cells while keeping both counts numeric
    %(the two dash patterns differ by design; see the builder's G44 comment
    %above). The old fixed formatspec pushed any such char cell through %d,
    %which prints each character's NUMERIC CODE (45 for '-') and recycles
    %the format, so one on-screen diff row came out as three garbled CSV
    %lines (ledger G47, verified by reproduction). Render every cell to a
    %string FIRST -- numeric cells through the exact former numeric
    %conversions, so all-numeric rows stay byte-identical -- then join and
    %write the whole row through ONE %s conversion. The strjoin matters:
    %passing nine strings to nine conversions would re-open the trap for an
    %EMPTY cell, because fprintf flattens its arguments and an empty char
    %contributes zero elements, truncating the row mid-line (adversarial
    %wave, 2026-08-29); a single pre-joined line cannot lose fields.
    %Column 6 (median trial count) prints with %g, not %d: the median of an
    %even-count trial vector is x.5, and sprintf renders a non-integer under
    %%d as %e (a P4 observation, 2026-09-04: '2.450000e+01' in a saved CSV).
    %%g prints an integer median exactly as %d did, so integer rows are
    %byte-identical.
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

%Name the coefficient that actually produced this include/exclude partition. The
%threshold below is applied to whichever coefficient the run selected, so labeling it
%"Dependability" on a generalizability run mis-attributes the screening criterion - and
%the partition really does change with the selection, so the two are not interchangeable.
%Derived the same way as in psyrat_depoverallt and psyrat_saveoveralltable above.
if isfield(psyrat_data,'relsummary') && isfield(psyrat_data.relsummary,'gcoeff_name') && ...
        strcmpi(psyrat_data.relsummary.gcoeff_name,'gen')
    rel_name = 'Generalizability';
else
    rel_name = 'Dependability';
end

%ask the user where the file should be saved. XLSX now writes cross-platform
%via writecell, so both formats are offered on every platform (including macOS).
[savename, savepath] = uiputfile(...
    {'*.xlsx','Excel File (.xlsx)';'*.csv',...
    'Comma-Separated Value File (.csv)'},...
    'Where would you like to save the data?');

[~,~,ext] = fileparts(fullfile(savepath,savename));

%save the information in either an excel or csv format
if strcmp(ext,'.xlsx')
    
    datap{1,1} = 'Data Generated on';
    datap{end+1,1} = datestr(clock); 
    datap{end+1,1} = sprintf('PsyRAT Toolbox v%s',psyrat_data.ver);
    datap{end+1,1} = '';
    datap{end+1,1} = sprintf('Dataset: %s',psyrat_data.rel.filename);
    datap{end+1,1} = sprintf('%s Cutoff: %0.4f',rel_name,...
        psyrat_data.relsummary.depcutoff);
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
    fprintf(fid,'%s Cutoff: %0.4f\n',rel_name,...
        psyrat_data.relsummary.depcutoff);
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
