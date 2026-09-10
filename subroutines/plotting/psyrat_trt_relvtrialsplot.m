function relplot = psyrat_trt_relvtrialsplot(varargin)
%Plot the reliability as a function of the number of trials included in
%an average (stratified by group and condition, if applicable)
%
%relplot = psyrat_trt_relvtrialsplot('psyrat_data',psyrat_data,'trials',[1 50],...
%   'relline',plotrelline,'relcutoff',relcutoff);
%
%
%Inputs
% psyrat_data - PsyRAT Toolbox data structure array. Variance components should
%  be included.
% trials - range of trials to plot
% relcutoff - reliability threshold for deeming data reliable (used for
%  reference line in plot)
%
%Optional Input
% relline - reliability estimate to plot (default: 2)
%  1 - lower limit
%  2 - point estimate
%  3 - upper limit
% CI - confidence interval width. Decimal from 0 to 1. (default: .95)
%
%
%Outputs
% relplot - figure handle for the plot
% figure that displays the relationship between dependability and the
%  number of trials retained for averaging stratified by group and
%  condition

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
        'See help psyrat_trt_relvtrialsplot for more information about inputs'));
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
            'See help psyrat_trt_relvtrialsplot for more information \n'));
    end
    
    %check if trials was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('trials',varargin),1);
    if ~isempty(ind)
        trials = varargin{ind+1}; 
    else 
        error('varargin:trials',... %Error code and associated error
            strcat('WARNING: Range of trials to plot not specified \n\n',... 
            'Please input trials: [min max]\n',...
            'See help psyrat_trt_relvtrialsplot for more information\n'));
    end
    
    %check if relcutoff was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('relcutoff',varargin),1);
    if ~isempty(ind)
        relcutoff = varargin{ind+1}; 
    else 
        error('varargin:relcutoff',... %Error code and associated error
            strcat('WARNING: Reliability threshold not specified \n\n',... 
            'Please input relcutoff\n',...
            'See help psyrat_trt_relvtrialsplot for more information\n'));
    end
    
    %check if relline was specified. 
    %If it is not found, set as default: 2.
    ind = find(strcmpi('relline',varargin),1);
    if ~isempty(ind)
        relline = varargin{ind+1}; 
        %make sure relline is 1, 2, or 3
        if ~any(relline==[1 2 3])
            error('varargin:relline',... %Error code and associated error
                strcat('WARNING: Reliability estimate to plot not specified \n\n',... 
                'Please enter a valid input for relline\n',...
                '1 - lower limit of credible interval\n',...
                '2 - point estimate of credible interval\n',...
                '3 - upper limit of credible interval\n',...
                'See help psyrat_trt_relvtrialsplot for more information \n'));
        end
    else 
        relline = 2; %default is point estimate
    end
    
    %check if CI was specified. 
    %If it is not found, set default as 95%
    ind = find(strcmpi('CI',varargin),1);
    if ~isempty(ind)
        ciperc = varargin{ind+1};
        if ciperc > 1 || ciperc < 0
            error('varargin:ci',... %Error code and associated error
                strcat('WARNING: Size of credible interval should ',...
                'be a value between 0 and 1\n',...
                'A value of ',sprintf(' %2.2f',ciperc),...
                ' is invalid\n',...
                'See help psyrat_trt_relvtrialsplot for more information \n'));
        end
    else 
        ciperc = .95; %default: 95%
    end
    
end

%make sure the user understood the the CI input is width not edges, if
%inputted incorrectly, provide a warning, but don't change
if ciperc < .5 
    str = sprintf(' %2.0f%%',100*ciperc);
    warning('ci:width',... %Warning code and associated warning
        strcat('WARNING: Size of credible interval is small \n\n',...
        'User specified a credible interval width of ',str,'%\n',...
        'If this was intended, ignore this warning.\n'));
end

%analysis-unit label fragments (trial vs split) for the figure title; defaults
%to "trial" so non-splits output is byte-identical.
usplit = 1;
if isfield(psyrat_data.rel,'splits'); usplit = psyrat_data.rel.splits; end
L = psyrat_unitlabels(usplit);

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

%RC-10. Occasion estimand (n'_o) for the curve. Every other consumer of the
%test-retest coefficient reads relsummary.nocc behind this guard
%(psyrat_depvtrialsplot, psyrat_trt_reloverallt, psyrat_variancet), and
%psyrat_relsummary passes it into every psyrat_rel_trt call that builds the
%table. This plot did not, so it silently used n'_o = 1 while the table beside
%it used the k the user selected: on one run the table read 0.495 and the
%figure 0.197 under an identical "Dependability" header. The figure draws a
%cutoff reference line and the user reads off the trial count where the curve
%crosses it, which is the same question the table's cutoff column answers, so
%the two artifacts were giving different answers to "how many trials do I
%need?". Under the shipped default (noccmode = 1) effnocc is 1 and this is a
%strict no-op.
if isfield(psyrat_data.relsummary,'nocc') && ~isempty(psyrat_data.relsummary.nocc)
    nocc = psyrat_data.relsummary.nocc;
else
    nocc = 1;
end

%create an x-axis for the number of observations
x = trials(1):trials(2);

%create an empty array for storing information into. RC-30. One row per
%requested trial count, i.e. numel(x) = trials(2)-trials(1)+1 rows, NOT
%trials(2) rows. The rows are indexed relative to the start of the range
%(row 1 is trials(1)), matching the convention psyrat_rel_trt already
%documents for its return vectors at psyrat_rel_trt.m:425-428 ("index
%relative to the start of the requested range so the outputs align to
%obs(1):obs(2)", then idx = ii - obs(1) + 1).
%Indexing this buffer by the ABSOLUTE trial number instead grew it to
%trials(2) rows while x had only numel(x) elements, so plot(x,plotrel)
%threw MATLAB:samelen for any trials(1) > 1 - a documented input that did
%not work. Under every production caller (trials(1) == 1) the two indexings
%are identical, so this is a strict no-op there.
%Kept as a deliberate twin of psyrat_depvtrialsplot; change both together.
plotrel = zeros(numel(x),0);

%see how many subplots are needed
if nevents > 2
    xplots = ceil(sqrt(nevents));
    yplots = ceil(sqrt(nevents));
elseif nevents == 2
    xplots = 1;
    yplots = 2;
elseif nevents == 1
    xplots = 1;
    yplots = 1;
end

%create the figure
relplot = figure;
relplot.Tag = 'psyrat_output';
set(gcf,'NumberTitle','Off');
relplot.Position = [125 630 900 450];
fsize = 16;

%extract the data and create the subplots for depplot
for eloc=1:nevents
    for gloc=1:ngroups
        switch relline
            case 1 %lower limit
                [plotrel(1:numel(x),gloc),~,~] = psyrat_rel_trt(...
                    'gcoeff',psyrat_data.relsummary.gcoeff,...
                    'reltype',psyrat_data.relsummary.reltype,...
                    'bp',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_id.raw,...
                    'bo',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_occ.raw,...
                    'bt',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_trl.raw,...
                    'txp',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_trlxid.raw,...
                    'oxp',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_occxid.raw,...
                    'txo',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_trlxocc.raw,...
                    'err',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_err.raw,...
                    'obs',[trials(1) trials(2)],...
                    'nocc',nocc,...
                    'CI',ciperc);
                plottitle = 'Lower Limit of 95% Credible Interval';
            case 2 %point estimate
                [~,plotrel(1:numel(x),gloc),~] = psyrat_rel_trt(...
                    'gcoeff',psyrat_data.relsummary.gcoeff,...
                    'reltype',psyrat_data.relsummary.reltype,...
                    'bp',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_id.raw,...
                    'bo',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_occ.raw,...
                    'bt',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_trl.raw,...
                    'txp',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_trlxid.raw,...
                    'oxp',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_occxid.raw,...
                    'txo',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_trlxocc.raw,...
                    'err',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_err.raw,...
                    'obs',[trials(1) trials(2)],...
                    'nocc',nocc,...
                    'CI',ciperc);
                plottitle = 'Point Estimate';
            case 3 %upper limit
                [~,~,plotrel(1:numel(x),gloc)] = psyrat_rel_trt(...
                    'gcoeff',psyrat_data.relsummary.gcoeff,...
                    'reltype',psyrat_data.relsummary.reltype,...
                    'bp',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_id.raw,...
                    'bo',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_occ.raw,...
                    'bt',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_trl.raw,...
                    'txp',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_trlxid.raw,...
                    'oxp',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_occxid.raw,...
                    'txo',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_trlxocc.raw,...
                    'err',psyrat_data.relsummary.data.g(gloc).e(eloc).sig_err.raw,...
                    'obs',[trials(1) trials(2)],...
                    'nocc',nocc,...
                    'CI',ciperc);
                plottitle = 'Upper Limit of 95% Credible Interval';
        end
    end
    
    switch psyrat_data.relsummary.reltype_name
        case 'ic'
            pref1 = 'Coefficient of Equivalence, ';
        case 'trt'
            pref1 = 'Coefficient of Stability, ';
        case 'ic_trt'
            pref1 = 'Coefficient of Equivalence and Stability, ';
    end
    
    switch psyrat_data.relsummary.gcoeff_name
        case 'dep'
            pref = ['Dependability v ' L.NumberOf ': '];
        case 'gen'
            pref = ['Generalizability v ' L.NumberOf ': '];
    end
    
    %RC-10. Name the occasion estimand too. The single- and multi-occasion
    %curves are different quantities and the axis label says only
    %"Dependability", so without this the window gives the reader no way to
    %tell which n'_o the curve was drawn at.
    set(gcf,'Name',[pref1 pref plottitle ...
        ' (' psyrat_trt_occlabel(psyrat_data) ')']);
    subplot(yplots,xplots,eloc); 
    h = plot(x,plotrel);
    
    %RC-30. The y-axis carries a reliability coefficient, so its upper limit
    %is 1. It previously read trials(1) - the MINIMUM trial count - which is
    %correct only by coincidence, because every production caller passes
    %trials(1) = 1.
    %
    %The two RC-30 defects were STACKED, and this one was never observable:
    %plot(x,plotrel) above throws MATLAB:samelen whenever trials(1) > 1, which
    %is exactly the condition under which this y-limit is wrong, so execution
    %never reached this line with a bad bound. Fixing the buffer indexing alone
    %would have turned an error into a silently nonsense y-scale (0 to 5 for
    %trials = [5 20]). That is why both lines had to move together.
    axis([0 trials(2) 0 1]);
    set(gca,'fontsize',16);
    
    if ~strcmpi(enames{eloc},'none')
        %Interpreter none here and on the legend below: an underscored event
        %or group name (inc_cor) must print literally, where the default TeX
        %interpreter would subscript the character after the underscore.
        title(enames{eloc},'FontSize',20,'Interpreter','none');
    end
    
    switch psyrat_data.relsummary.gcoeff_name
        case 'dep'
            ylabel('Dependability','FontSize',fsize);
        case 'gen'
            ylabel('Generalizability','FontSize',fsize);
    end
    
    xlabel('Number of Observations','FontSize',fsize);
    psyrat_addhline(relcutoff,'Color','b','LineStyle',':');
    if analysis ~= 1 && analysis ~= 3
        leg = legend(gnames{:},'Location','southeast');
        set(leg,'FontSize',fsize,'Interpreter','none');
        
        legend('boxoff');
    end

end

end
