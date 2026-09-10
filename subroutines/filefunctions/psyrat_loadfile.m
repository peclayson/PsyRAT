function dataout = psyrat_loadfile(varargin)
%Loads and prepares the data file for dependability analyses
%
%psyrat_loadfile('file',filename)
%
%
%Required Inputs:
% psyrat_prefs - preferences from PsyRAT Toolbox
% psyrat_data - data from PsyRAT Toolbox
%       OR
% file - location of file to be loaded and prepared for dependability
%  analyses
%
%Optional Inputs:
% idcol - column label for the participant id variable (default: 'id')
% meascol - column label for the measurement to be analyzed (default:
%  'measurement'). The selected column is renamed to 'meas' internally and
%  is the canonical analyzed column. Scripted callers are encouraged to pass
%  meascol explicitly to point at their score column rather than rely on the
%  default.
% groupcol - column label for the group variable. If no label is provided 
%  it is assumed that there is only one group in the data file.
% whichgroups - cell array of labels for the group types found in the 
%  group column to use
% eventcol - column label for the event variable. If no label is provided
%  it is assumed there is only one event type in the data file.
% whichevents - cell array of labels for the event types found in the 
%  event column to use
% timecol - column label for the occasion variable. If no label is provided
%  it is assumed there is only one occasion in the data file.
% whichtimes - cell array of labels for the occasions found in the
%  occasion column to use
% weightcol - column label for the split-weight variable (n_i, the number of
%  items/trials averaged into each split mean). Only used by the splits
%  analysis; carried as a numeric, positive-integer count.
% dim1col - column label for the first continuous between-person covariate used
%  by the dynamic/conditional-reliability analysis (default: none).
% dim2col - column label for the second continuous between-person covariate used
%  by the dynamic/conditional-reliability analysis (default: none).
% dataraw - raw data table outputted from psyrat_startproc (so Matlab doesn't
%  have to re-load entire table)
%
%Output:
% dataout - matlab table with prepared data for reliability analysis.
%  Depending on specifications, the table will have 2 to 8 columns.
%  id: Subject ID (string variable)
%  meas: Measurement
%  group: Group (only when specified, string variable)
%  event: Event Type (only when specified, string variable)
%  time: Occasion (only when specified, string variable)
%  weight: Split weight n_i (only when specified, numeric item count)
%  dim1: First continuous between-person covariate (dynamic reliability, only
%   when specified)
%  dim2: Second continuous between-person covariate (dynamic reliability, only
%   when specified)
%
%Input table requirements:
% The raw input is a long-format table with one row per single-trial score.
%  It MUST contain an id column and a measurement column (located via idcol
%  and meascol; a missing header raises an error). The measurement column
%  MUST be numeric; non-numeric data raise an error during loading. The
%  group, event, and time columns are optional. The splits analysis
%  additionally requires a numeric weight column (n_i; located via
%  weightcol).

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

%if psyrat_prefs and psyrat_data were defined, pull inputs
if ~isempty(psyrat_prefs) && ~isempty(psyrat_data)
    file = psyrat_data.raw.filename;
    idcolname = psyrat_data.proc.idheader;
    meascolname = psyrat_data.proc.measheader;
    eventcolname = psyrat_data.proc.eventheader;
    whichevents = psyrat_data.proc.whichevents;
    groupcolname = psyrat_data.proc.groupheader;
    whichgroups = psyrat_data.proc.whichgroups;
    timecolname = psyrat_data.proc.timeheader;
    whichtimes = psyrat_data.proc.whichtimes;
    %dimension columns are only used by the dynamic-reliability analysis; guard
    %so non-dynrel runs (which lack these proc fields) are unaffected.
    if isfield(psyrat_data.proc,'dim1header')
        dim1colname = psyrat_data.proc.dim1header;
    else
        dim1colname = '';
    end
    if isfield(psyrat_data.proc,'dim2header')
        dim2colname = psyrat_data.proc.dim2header;
    else
        dim2colname = '';
    end
    %split-weight (n_i) column is only used by the splits analysis; guard so
    %non-split runs (which lack this proc field) are unaffected.
    if isfield(psyrat_data.proc,'weightheader')
        weightcolname = psyrat_data.proc.weightheader;
    else
        weightcolname = '';
    end
    dataraw = psyrat_data.raw.data;
end
    

%somersault through varargin inputs to check for which inputs were
%defined and store those values. 
if ~isempty(varargin) && (isempty(psyrat_prefs) || isempty(psyrat_data))
    
    %the optional inputs check assumes that there was an even number of 
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)  
        error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',... 
        'Make sure each variable input is paired with a value \n',...
        'See help psyrat_loadfile for more information on inputs'));
    end
    
    %check if a location for the file to be loaded was specified. 
    %If it is not found, set display error.
    ind = find(strcmp('file',varargin),1);
    if ~isempty(ind)
        file = cell2mat(varargin(ind+1)); 
    else 
        error('varargin:nofile',... %Error code and associated error
        strcat('Error: File location not specified \n\n',... 
        'Please input the full path specifying the file to be loaded \n'));
    end
   
    %check if id column was specified 
    %If it is not found, use default 'id'
    ind = find(strcmp('idcol',varargin),1);
    if ~isempty(ind)
        idcolname = cell2mat(varargin(ind+1)); 
    else 
        idcolname = 'id';
    end
    
    %check if group column was specified 
    %If it is not found, assume only one group is present in file
    ind = find(strcmp('groupcol',varargin),1);
    if ~isempty(ind)
        groupcolname = cell2mat(varargin(ind+1)); 
    else 
        groupcolname = '';
    end
    
    %check if whichgroups was specified 
    %If it is not found, assume only one group is present in file
    ind = find(strcmp('whichgroups',varargin),1);
    if ~isempty(ind)
        whichgroups = varargin{ind+1};
    else 
        whichgroups = '';
    end
    
    %check if measurement column was specified 
    %If it is not found, assume default 'measurement'
    ind = find(strcmp('meascol',varargin),1);
    if ~isempty(ind)
        meascolname = cell2mat(varargin(ind+1)); 
    else 
        meascolname = 'measurement';
    end
    
    %check if event type column was specified 
    %If it is not found, assume only one event type is present in file
    ind = find(strcmp('eventcol',varargin),1);
    if ~isempty(ind)
        eventcolname = cell2mat(varargin(ind+1)); 
    else 
        eventcolname = '';
    end
    
    %check if whichgroups was specified 
    %If it is not found, assume only one group is present in file
    ind = find(strcmp('whichevents',varargin),1);
    if ~isempty(ind)
        whichevents = varargin{ind+1};
    else 
        whichevents = '';
    end
    
    %check if occasion type column was specified
    %If it is not found, assume only one occasion is present in file
    ind = find(strcmp('timecol',varargin),1);
    if ~isempty(ind)
        timecolname = cell2mat(varargin(ind+1));
    else 
        timecolname = '';
    end
    
    %check if whichtimes was specified 
    %If it is not found, assume only one occasion is present in file
    ind = find(strcmp('whichtimes',varargin),1);
    if ~isempty(ind)
        whichtimes = varargin{ind+1};
    else 
        whichtimes = '';
    end
    
    %check if dimension columns were specified (dynamic-reliability analysis).
    %These are continuous, between-person covariates; no which* filter applies.
    ind = find(strcmp('dim1col',varargin),1);
    if ~isempty(ind)
        dim1colname = cell2mat(varargin(ind+1));
    else
        dim1colname = '';
    end

    ind = find(strcmp('dim2col',varargin),1);
    if ~isempty(ind)
        dim2colname = cell2mat(varargin(ind+1));
    else
        dim2colname = '';
    end

    %check if the split-weight (n_i) column was specified (splits analysis).
    %It carries the number of items/trials averaged into each split mean.
    ind = find(strcmp('weightcol',varargin),1);
    if ~isempty(ind)
        weightcolname = cell2mat(varargin(ind+1));
    else
        weightcolname = '';
    end

    %check if dataraw was specified
    %If it is not found, create an empty variable
    ind = find(strcmp('dataraw',varargin),1);
    if ~isempty(ind)
        dataraw = varargin{ind+1};
    else
        dataraw = '';
    end

    %Mirror the GUI workflow's psyrat_data.proc.* selections so the optional
    %which* row-subsetting below applies in the direct-argument workflow too
    %(B6). The filter blocks read psyrat_data.proc.*; in this path psyrat_data
    %is otherwise empty, so without this the parsed which* locals were ignored
    %and programmatic group/event/occasion subsetting silently did nothing.
    %Only non-empty selections are stored, so unspecified filters stay inactive
    %(local_has_proc_field requires the field to be present and non-empty).
    psyrat_data = struct('proc', struct());
    if ~isempty(whichgroups)
        psyrat_data.proc.whichgroups = whichgroups;
    end
    if ~isempty(whichevents)
        psyrat_data.proc.whichevents = whichevents;
    end
    if ~isempty(whichtimes)
        psyrat_data.proc.whichtimes = whichtimes;
    end

elseif isempty(varargin)
    
    error('varargin:incomplete',... %Error code and associated error
    strcat('Error: Optional inputs are incomplete \n\n',... 
    'Make sure each variable input is paired with a value \n',...
    'See help psyrat_loadfile for more information on inputs'));
    
end %if ~isempty(varargin)

%load file if it has not been done already
if isempty(dataraw)
    dataraw = psyrat_readtable('file',file);
end

%grab the filename to store
[~,filename] = fileparts(file); 

%make sure all of the necessary headers are present in the file then load
%the data into a table to be outputted for analysis
colnames = dataraw.Properties.VariableNames;
dataout = table;
dataout.Properties.Description = filename;

%ensure that the id column is present
if ~sum(strcmpi(colnames,idcolname)) 
    if ~exist('headerror','var')
        headerror{1} = 'Subject ID';
    else
        headerror{end+1} = 'Subject ID';
    end
elseif sum(strcmpi(colnames,idcolname)) 
    dataout.id = dataraw{:,strcmpi(colnames,idcolname)};
    colnames = dataraw.Properties.VariableNames;
end

%ensure that the measurement column is present
if ~sum(strcmpi(colnames,meascolname)) 
    if ~exist('headerror','var')
        headerror{1} = 'Measurement';
    else
        headerror{end+1} = 'Measurement';
    end
elseif sum(strcmpi(colnames,meascolname)) 
    dataout.meas = dataraw{:,strcmpi(colnames,meascolname)};
end

%ensure that the group column is present if the group header should be
%there
if ~sum(strcmpi(colnames,groupcolname)) && ~isempty(groupcolname)
    if ~exist('headerror','var')
        headerror{1} = 'Group';
    else
        headerror{end+1} = 'Group';
    end
elseif sum(strcmpi(colnames,groupcolname)) && ~isempty(groupcolname)
    dataout.group = dataraw{:,strcmpi(colnames,groupcolname)};
end

%ensure that the event column is present if the event header should be
%there
if ~sum(strcmpi(colnames,eventcolname)) && ~isempty(eventcolname)
    if ~exist('headerror','var')
        headerror{1} = 'Event';
    else
        headerror{end+1} = 'Event';
    end
elseif sum(strcmpi(colnames,eventcolname)) && ~isempty(eventcolname)
    dataout.event = dataraw{:,strcmpi(colnames,eventcolname)};
end

%ensure that the time column is present if the time header should be
%there
if ~sum(strcmpi(colnames,timecolname)) && ~isempty(timecolname)
    if ~exist('headerror','var')
        headerror{1} = 'Occasion';
    else
        headerror{end+1} = 'Occasion';
    end
elseif sum(strcmpi(colnames,timecolname)) && ~isempty(timecolname)
    dataout.time = dataraw{:,strcmpi(colnames,timecolname)};
end

%ensure the dimension columns are present if requested (dynamic-reliability
%analysis). They are carried as numeric, between-person covariates; the
%between-subject standardization happens later in psyrat_computevarcomp. Use
%find (not sum) so a header that matches MORE than one source column is caught
%as ambiguous rather than silently extracted as an N-by-K matrix, which would
%slip past the numeric guard below and corrupt the later standardization.
dim1idx = find(strcmpi(colnames,dim1colname));
if isempty(dim1idx) && ~isempty(dim1colname)
    if ~exist('headerror','var')
        headerror{1} = 'Dimension 1';
    else
        headerror{end+1} = 'Dimension 1';
    end
elseif numel(dim1idx) > 1 && ~isempty(dim1colname)
    error('dim1:ambiguousheader',... %Error code and associated error
        strcat('Error: The Dimension 1 header matches more than one column\n',...
        'The header specified was',[' ' dim1colname],'\n',...
        'Please give each input a unique column header\n',...
        'See help psyrat_loadfile for more information\n'));
elseif numel(dim1idx) == 1
    dataout.dim1 = dataraw{:,dim1idx};
end

dim2idx = find(strcmpi(colnames,dim2colname));
if isempty(dim2idx) && ~isempty(dim2colname)
    if ~exist('headerror','var')
        headerror{1} = 'Dimension 2';
    else
        headerror{end+1} = 'Dimension 2';
    end
elseif numel(dim2idx) > 1 && ~isempty(dim2colname)
    error('dim2:ambiguousheader',... %Error code and associated error
        strcat('Error: The Dimension 2 header matches more than one column\n',...
        'The header specified was',[' ' dim2colname],'\n',...
        'Please give each input a unique column header\n',...
        'See help psyrat_loadfile for more information\n'));
elseif numel(dim2idx) == 1
    dataout.dim2 = dataraw{:,dim2idx};
end

%ensure the split-weight column (n_i, items per split) is present if requested
%(splits analysis). It is carried as a numeric count and validated below. Use
%find (not sum) so a header matching MORE than one source column is caught as
%ambiguous rather than silently extracted as an N-by-K matrix.
weightidx = find(strcmpi(colnames,weightcolname));
if isempty(weightidx) && ~isempty(weightcolname)
    if ~exist('headerror','var')
        headerror{1} = 'Split Weight (n_i)';
    else
        headerror{end+1} = 'Split Weight (n_i)';
    end
elseif numel(weightidx) > 1 && ~isempty(weightcolname)
    error('weight:ambiguousheader',... %Error code and associated error
        strcat('Error: The Split Weight (n_i) header matches more than one column\n',...
        'The header specified was',[' ' weightcolname],'\n',...
        'Please give each input a unique column header\n',...
        'See help psyrat_loadfile for more information\n'));
elseif numel(weightidx) == 1
    dataout.weight = dataraw{:,weightidx};
end

%error catch in case headers for the columns needed are not specified. Let
%the user know which columns were problematic
if exist('headerror','var')
    error('varargin:colheaders',... %Error code and associated error
    strcat('Error: Column headers not properly specified \n\n',... 
    'Please specify the headers for\n',char(strjoin(headerror,', ')),'\n',...
    'See help psyrat_loadfile \n'));
end

%verify that the data in the measurement column is numeric
if ~isnumeric(dataout.meas(:))
    error('meas:notnumeric',... %Error code and associated error
    strcat('Error: The data in the measurement column is not numeric\n',...
    'Was the wrong column specified as mesurement?\n\n',...
    'The column specified was',[' ' meascolname],'\n',...
    'Please specify a column with numeric data (ERP measurements)\n',...
    'See help psyrat_loadfile for more information\n'));
end

%verify the dimension columns are numeric (continuous between-person
%covariates), mirroring the measurement check so a mis-selected text column is
%reported clearly at load time rather than as a late standardization error.
if ismember('dim1',dataout.Properties.VariableNames) && ~isnumeric(dataout.dim1(:))
    error('dim1:notnumeric',... %Error code and associated error
        strcat('Error: The data in the Dimension 1 column is not numeric\n',...
        'The column specified was',[' ' dim1colname],'\n',...
        'Dimension predictors must be continuous numeric covariates\n',...
        'See help psyrat_loadfile for more information\n'));
end
if ismember('dim2',dataout.Properties.VariableNames) && ~isnumeric(dataout.dim2(:))
    error('dim2:notnumeric',... %Error code and associated error
        strcat('Error: The data in the Dimension 2 column is not numeric\n',...
        'The column specified was',[' ' dim2colname],'\n',...
        'Dimension predictors must be continuous numeric covariates\n',...
        'See help psyrat_loadfile for more information\n'));
end

%verify the split-weight column is a valid item count: numeric, finite, and a
%positive whole number (n_i is the number of items/trials averaged into each
%split mean). This is required because the nonparallel split model uses n_i as
%a precision weight and the split formulas divide item-level error by n_i.
if ismember('weight',dataout.Properties.VariableNames)
    w = dataout.weight(:);
    if ~isnumeric(w) || any(~isfinite(w)) || any(w <= 0) || any(mod(w,1) ~= 0)
        error('weight:invalid',... %Error code and associated error
            strcat('Error: The Split Weight (n_i) column must contain positive\n',...
            'whole numbers (the number of items/trials in each split)\n',...
            'The column specified was',[' ' weightcolname],'\n',...
            'See help psyrat_loadfile for more information\n'));
    end
end

%normalize label columns that are neither cellstr nor numeric -- string,
%categorical, datetime, duration, char (finding B12). Those types otherwise
%survive every layer and fail deep in the estimation core after sampling has
%already completed.
%
%THIS MUST RUN BEFORE THE which* FILTERS BELOW, and the reason is measured,
%not assumed. Those filters compare the column against the user's selection
%with strcmpi, and strcmpi returns a scalar false for a datetime or
%categorical column rather than matching element-wise. Placed after them, a
%dated occasion column plus an occasion selection dies with a misleading
%'times:occasionmismatch' instead of loading.
%
%Numeric columns are deliberately left for the SECOND call at the very end of
%this function: the filters below compare a numeric column with == against a
%numeric selection, so converting numerics HERE would break a path that works
%today. That second call passes convertnumeric = true and converts them once
%the filters have run.
dataout = psyrat_normalize_label_columns(dataout);

%if groups, event and/or time columns are specified, check whether only
%certain groups or event should be processed. If there are specific groups or
%events to process, make sure those groups and events exist in the data (in
%cases which function was not called by gui).
if ~isempty(groupcolname) && local_has_proc_field(psyrat_data,'whichgroups')
    for ii = 1:length(psyrat_data.proc.whichgroups)
        if iscell(psyrat_data.proc.whichgroups)
            if ~isnumeric(psyrat_data.proc.whichgroups{ii})
                if ~any(strcmpi(dataout.group,psyrat_data.proc.whichgroups{ii}))
                    error('groups:groupmismatch',... %Error code and associated error
                        strcat('Error: Specified groups to process do not exist in data\n',...
                        'See help psyrat_loadfile for more information\n'));
                end
            else
                if ~any(find(dataout.group==psyrat_data.proc.whichgroups{ii}))
                    error('groups:groupmismatch',... %Error code and associated error
                        strcat('Error: Specified groups to process do not exist in data\n',...
                        'See help psyrat_loadfile for more information\n'));
                end
            end
        else
            if ~isnumeric(psyrat_data.proc.whichgroups(ii))
                if ~any(strcmpi(dataout.group,psyrat_data.proc.whichgroups(ii)))
                    error('groups:groupmismatch',... %Error code and associated error
                        strcat('Error: Specified groups to process do not exist in data\n',...
                        'See help psyrat_loadfile for more information\n'));
                end
            else
                if ~any(find(dataout.group==psyrat_data.proc.whichgroups(ii)))
                    error('groups:groupmismatch',... %Error code and associated error
                        strcat('Error: Specified groups to process do not exist in data\n',...
                        'See help psyrat_loadfile for more information\n'));
                end
            end
            
        end
    end
    
    try
        if iscell(psyrat_data.proc.whichgroups)
            if ~isnumeric(psyrat_data.proc.whichgroups{:})
                dataout = dataout(ismember(...
                    dataout.group,psyrat_data.proc.whichgroups),:);
            else
                try
                    dataout = dataout(ismember(...
                        dataout.group,psyrat_data.proc.whichgroups{:}),:);
                catch
                    dataout = dataout(ismember(...
                        dataout.group,[psyrat_data.proc.whichgroups{:}]),:);
                end
            end
        else
            if ~isnumeric(psyrat_data.proc.whichgroups(:))
                dataout = dataout(ismember(...
                    dataout.group,psyrat_data.proc.whichgroups),:);
            else
                dataout = dataout(ismember(...
                    dataout.group,psyrat_data.proc.whichgroups(:)),:);
            end
        end
    catch
        if iscell(psyrat_data.proc.whichgroups)
            if ~isnumeric(psyrat_data.proc.whichgroups{1})
                dataout = dataout(ismember(...
                    dataout.group,psyrat_data.proc.whichgroups),:);
            else
                try
                    dataout = dataout(ismember(...
                        dataout.group,psyrat_data.proc.whichgroups{:}),:);
                catch
                    dataout = dataout(ismember(...
                        dataout.group,[psyrat_data.proc.whichgroups{:}]),:);
                end
            end
        else
            if ~isnumeric(psyrat_data.proc.whichgroups(1))
                dataout = dataout(ismember(...
                    dataout.group,psyrat_data.proc.whichgroups),:);
            else
                dataout = dataout(ismember(...
                    dataout.group,psyrat_data.proc.whichgroups(:)),:);
            end
        end
        
    end
end

if ~isempty(eventcolname) && local_has_proc_field(psyrat_data,'whichevents')
    for ii = 1:length(psyrat_data.proc.whichevents)
        if iscell(psyrat_data.proc.whichevents)
            if ~isnumeric(psyrat_data.proc.whichevents{ii})
                if ~any(strcmpi(dataout.event,psyrat_data.proc.whichevents{ii}))
                    error('events:eventmismatch',... %Error code and associated error
                        strcat('Error: Specified events to process do not exist in data\n',...
                        'See help psyrat_loadfile for more information\n'));
                end
            else
                if ~any(find(dataout.event==psyrat_data.proc.whichevents{ii}))
                    error('events:eventmismatch',... %Error code and associated error
                        strcat('Error: Specified events to process do not exist in data\n',...
                        'See help psyrat_loadfile for more information\n'));
                end
            end
        else
            if ~isnumeric(psyrat_data.proc.whichevents(ii))
                if ~any(strcmpi(dataout.event,psyrat_data.proc.whichevents(ii)))
                    error('events:eventmismatch',... %Error code and associated error
                        strcat('Error: Specified events to process do not exist in data\n',...
                        'See help psyrat_loadfile for more information\n'));
                end
            else
                if ~any(find(dataout.event==psyrat_data.proc.whichevents(ii)))
                    error('events:eventmismatch',... %Error code and associated error
                        strcat('Error: Specified events to process do not exist in data\n',...
                        'See help psyrat_loadfile for more information\n'));
                end
            end
        end
    end
    
    
    
    try
        if iscell(psyrat_data.proc.whichevents)
            if ~isnumeric(psyrat_data.proc.whichevents{:})
                dataout = dataout(ismember(...
                    dataout.event,psyrat_data.proc.whichevents),:);
            else
                try
                    dataout = dataout(ismember(...
                        dataout.event,psyrat_data.proc.whichevents{:}),:);
                catch
                    dataout = dataout(ismember(...
                        dataout.event,[psyrat_data.proc.whichevents{:}]),:);
                end
            end
        else
            if ~isnumeric(psyrat_data.proc.whichevents(:))
                dataout = dataout(ismember(...
                    dataout.event,psyrat_data.proc.whichevents),:);
            else
                dataout = dataout(ismember(...
                    dataout.event,psyrat_data.proc.whichevents(:)),:);
            end
        end
        
    catch
        if iscell(psyrat_data.proc.whichevents)
            if ~isnumeric(psyrat_data.proc.whichevents{1})
                dataout = dataout(ismember(...
                    dataout.event,psyrat_data.proc.whichevents),:);
            else
                try
                    dataout = dataout(ismember(...
                        dataout.event,psyrat_data.proc.whichevents{:}),:);
                catch
                    dataout = dataout(ismember(...
                        dataout.event,[psyrat_data.proc.whichevents{:}]),:);
                end
            end
        else
            if ~isnumeric(psyrat_data.proc.whichevents(1))
                dataout = dataout(ismember(...
                    dataout.event,psyrat_data.proc.whichevents),:);
            else
                dataout = dataout(ismember(...
                    dataout.event,psyrat_data.proc.whichevents(:)),:);
            end
        end
    end

end

if ~isempty(timecolname) && local_has_proc_field(psyrat_data,'whichtimes')
    for ii = 1:length(psyrat_data.proc.whichtimes)
        if iscell(psyrat_data.proc.whichtimes)
            if ~isnumeric(psyrat_data.proc.whichtimes{ii})
                if ~any(strcmpi(dataout.time,psyrat_data.proc.whichtimes{ii}))
                    error('times:occasionmismatch',... %Error code and associated error
                        strcat('Error: Specified occasions to process do not exist in data\n',...
                        'See help psyrat_loadfile for more information\n'));
                end
            else
                if ~any(find(dataout.time==psyrat_data.proc.whichtimes{ii}))
                    error('times:occasionmismatch',... %Error code and associated error
                        strcat('Error: Specified occasions to process do not exist in data\n',...
                        'See help psyrat_loadfile for more information\n'));
                end
            end
        else
            if ~isnumeric(psyrat_data.proc.whichtimes(ii))
                if ~any(strcmpi(dataout.time,psyrat_data.proc.whichtimes(ii)))
                    error('times:occasionmismatch',... %Error code and associated error
                        strcat('Error: Specified occasions to process do not exist in data\n',...
                        'See help psyrat_loadfile for more information\n'));
                end
            else
                if ~any(find(dataout.time==psyrat_data.proc.whichtimes(ii)))
                    error('times:occasionmismatch',... %Error code and associated error
                        strcat('Error: Specified occasions to process do not exist in data\n',...
                        'See help psyrat_loadfile for more information\n'));
                end
            end
        end
    end
    
    
    
    try
        if iscell(psyrat_data.proc.whichtimes)
            if ~isnumeric(psyrat_data.proc.whichtimes{:})
                dataout = dataout(ismember(...
                    dataout.time,psyrat_data.proc.whichtimes),:);
            else
                try
                    dataout = dataout(ismember(...
                        dataout.time,psyrat_data.proc.whichtimes{:}),:);
                catch
                    dataout = dataout(ismember(...
                        dataout.time,[psyrat_data.proc.whichtimes{:}]),:);
                end
            end
        else
            if ~isnumeric(psyrat_data.proc.whichtimes(:))
                dataout = dataout(ismember(...
                    dataout.time,psyrat_data.proc.whichtimes),:);
            else
                dataout = dataout(ismember(...
                    dataout.time,psyrat_data.proc.whichtimes(:)),:);
            end
        end
        
    catch
        if iscell(psyrat_data.proc.whichtimes)
            if ~isnumeric(psyrat_data.proc.whichtimes{1})
                dataout = dataout(ismember(...
                    dataout.time,psyrat_data.proc.whichtimes),:);
            else
                try
                    dataout = dataout(ismember(...
                        dataout.time,psyrat_data.proc.whichtimes{:}),:);
                catch
                    dataout = dataout(ismember(...
                        dataout.time,[psyrat_data.proc.whichtimes{:}]),:);
                end
            end
        else
            if ~isnumeric(psyrat_data.proc.whichtimes(1))
                dataout = dataout(ismember(...
                    dataout.time,psyrat_data.proc.whichtimes),:);
            else
                dataout = dataout(ismember(...
                    dataout.time,psyrat_data.proc.whichtimes(:)),:);
            end
        end
    end

end

%verify that each participant has at least one measurement per event type
%(this works because it's already been verified that there are no empty
%cells)
if ~isempty(eventcolname)
    datacheck = varfun(@mean,dataout,'InputVariables','meas',...
       'GroupingVariables',{'id','event'});
    eventcount = varfun(@numel,datacheck,'InputVariables','event',...
        'GroupingVariables','id');
    if length(unique(eventcount.GroupCount)) > 1
        eventcount.GroupCount = [];
        eventcount.Properties.VariableNames{2} = 'Number_of_Events';
        display(eventcount);
        error('meas:mismatchedevents',... %Error code and associated error
        strcat('Error: Participants differ in the number of distinct events\n',...
        '(this check compares the COUNT of events per participant; it does\n',...
        'not verify that every participant has every event)\n',...
        'The ids and number of events found for each participant are printed',...
        '\nabove to help in finding the problem\n'));
    end
end

%verify that each participant has at least one measurement per occasion and
%per event
%(this works because it's already been verified that there are no empty
%cells)
if ~isempty(timecolname) && ~isempty(eventcolname)
    datacheck = varfun(@mean,dataout,'InputVariables','meas',...
       'GroupingVariables',{'id','event','time'});
    timecount = varfun(@numel,datacheck,'InputVariables','time',...
        'GroupingVariables','id');
    if length(unique(timecount.GroupCount)) > 1
        timecount.GroupCount = [];
        timecount.Properties.VariableNames{2} = 'Number_of_Occasions';
        display(timecount);
        error('meas:mismatchedoccasions',... %Error code and associated error
        strcat('Error: Participants differ in the number of distinct occasion\n',...
        'cells (this check compares COUNTS per participant; it does not\n',...
        'verify that every occasion or event-by-occasion cell is present\n',...
        'for every participant)\n',...
        'The ids and number of occasions found for each participant are printed',...
        '\nabove to help in finding the problem\n'));
    end
elseif ~isempty(timecolname) && isempty(eventcolname)
        datacheck = varfun(@mean,dataout,'InputVariables','meas',...
       'GroupingVariables',{'id','time'});
    timecount = varfun(@numel,datacheck,'InputVariables','time',...
        'GroupingVariables','id');
    if length(unique(timecount.GroupCount)) > 1
        timecount.GroupCount = [];
        timecount.Properties.VariableNames{2} = 'Number_of_Occasions';
        display(timecount);
        error('meas:mismatchedoccasions',... %Error code and associated error
        strcat('Error: Participants differ in the number of distinct occasion\n',...
        'cells (this check compares COUNTS per participant; it does not\n',...
        'verify that every occasion or event-by-occasion cell is present\n',...
        'for every participant)\n',...
        'The ids and number of occasions found for each participant are printed',...
        '\nabove to help in finding the problem\n'));
    end
end

%turn all of the ids, groups, events, and occasions into strings. PsyRAT
%toolbox will use this later.
%
%SECOND call to psyrat_normalize_label_columns in this function, and the two
%differ only in the numeric opt-in. The first (above the which* filters)
%takes the default convertnumeric = false, because those filters compare a
%numeric column with == against a numeric selection. This one passes TRUE,
%because by here the filters have run and numeric labels are ready to become
%text. Column presence is handled inside -- group, event and time are skipped
%when the table does not carry them, which is the same condition the
%~isempty(groupcolname) guards used to express.
%
%PLACEMENT IS DELIBERATE: this stays at the END of the function, BELOW the
%varfun completeness checks above, which group by id/event/time. Those checks
%run on the still-numeric columns today, and hoisting the conversion above
%them would switch them to text grouping for no reason.
%
%WHAT THIS REPLACED, AND WHY (finding B12's loader half, ruled 2026-08-12).
%Four hand-rolled cellstr(num2str(...)) blocks, one per label column. The
%ruling was about the SPELLING: num2str right-align PADS a column vector to a
%common width, so [1;2;10] shipped as {' 1';' 2';'10'} with leading spaces,
%while every other label producer in the toolbox emits unpadded
%cellstr(string(...)). Measured consequences of the padding: it breaks four
%cross-source string joins against those unpadded labels (the dodmap checks in
%psyrat_computevarcompwarp and psyrat_computevarcomp, psyrat_report's
%ismember(mapnames,events), and psyrat_startview_dodiff's, which fails
%SILENTLY and falls back to dodmap = 1:4), and it puts a leading space into
%unquoted CSV export fields.
%
%THE HONEST COST, which the ruling accepted: level ORDER moves for numeric
%labels, because sort gives {' 1',' 2','10'} padded but {'1','10','2'}
%unpadded. That was only safe to take AFTER B23. Before it, the
%difference-score cases read event order from a SORTED unique while the Stan
%model matrix was keyed to a 'stable' one, so order carried a VALUE and this
%change would have silently reshuffled which numeric-label datasets were
%accidentally correct. B23 aligned REL.events to REL.diff_names, so order is
%presentational there now.
%
%TWO GUARDS COME WITH THE SWAP -- approved riders, and both user-visible
%changes to previously-accepted input. A NaN numeric label used to become a
%level literally named 'NaN' (num2str renders it, and preflight's missing-cell
%guard cannot see a literal 'NaN'); it now raises labels:missingvalues. Two
%distinct numeric labels rendering to the same text used to merge silently
%into one level; that now raises labels:numericlabelcollision.
dataout = psyrat_normalize_label_columns(dataout,true);

end

function tf = local_has_proc_field(psyrat_data,fieldname)
%Return true only when psyrat_data carries a usable psyrat_data.proc.(fieldname).
%In the GUI workflow psyrat_data is a struct whose proc field holds the
%group/event/occasion selections. In the direct-argument workflow
%psyrat_findprefsdata returns '' (empty), so struct-based row filtering must
%be skipped rather than dereferenced (which would error). Column population
%upstream is unaffected; only the optional which* row subsetting is gated.
tf = isstruct(psyrat_data) && isfield(psyrat_data,'proc') && ...
    isstruct(psyrat_data.proc) && isfield(psyrat_data.proc,fieldname) && ...
    ~isempty(psyrat_data.proc.(fieldname));
end
