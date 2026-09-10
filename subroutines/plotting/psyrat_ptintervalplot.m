function ptintplot = psyrat_ptintervalplot(varargin)
%Plot the a point estimate with its associated confidence interval
%
%psyrat_ptintervalplot('psyrat_data',psyrat_data,'stat','bet')
%
%
%Inputs
% psyrat_data - PsyRAT Toolbox data structure array. Variance components should
%  be included.
% stat - three options:
%  ICC - to plot ICC estimates
%  Bet - to plot between-subject standard deviations
%  Wit - to plot within-subject standard deviations
%
%Optional input
% CI - confidence interval width. Decimal from 0 to 1. (default: .95)
%
%Outputs
% ptintplot - figure handle for the plot
% figure that displays a point estimate and its credible interval

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
    
    %check if stat was specified.
    %If it is not found, set as default: 2.
    ind = find(strcmpi('stat',varargin),1);
    if ~isempty(ind)
        stat = varargin{ind+1};
        %make sure depline is 1, 2, or 3
        if ~strcmpi('ICC',stat) && ~strcmpi('Bet',stat) &&...
                ~strcmpi('Wit',stat)
            error('varargin:stat',... %Error code and associated error
                strcat('WARNING: stat not properly specified \n\n',...
                'Please input the stat for the plot\n',...
                'ICC - to plot ICC estimates\n',...
                'Bet - to plot between-subject standard deviations\n',...
                'Wit - to plot within-subject standard deviations\n',...
                'See help psyrat_ptintervalplot for more information \n'));
        end
    else
        error('varargin:stat',... %Error code and associated error
            strcat('WARNING: stat not specified \n\n',...
            'Please input the stat for the plot\n',...
            'ICC - to plot ICC estimates\n',...
            'Bet - to plot between-subject standard deviations\n',...
            'Wit - to plot within-subject standard deviations\n',...
            'See help psyrat_ptintervalplot for more information \n'));
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
                'See help psyrat_depvtrialsplot for more information \n'));
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

ptintplot = figure;
set(gcf,'Tag','psyrat_output');
ptintplot.Position = [125 65 900 450];
fsize = 16;
set(gcf,'NumberTitle','Off');
ptest = zeros(nevents,ngroups);
llest = zeros(nevents,ngroups);
ulest = zeros(nevents,ngroups);
offsetm = zeros(nevents,ngroups);

if ~strcmp(psyrat_data.rel.analysis,'ic_sserrvar')
    
    %grab icc information for each event and group
    for gloc=1:ngroups
        for eloc=1:nevents
            %        ptest(eloc,gloc) = relsummary.group(gloc).event(eloc).icc.m;
            %        llest(eloc,gloc) = relsummary.group(gloc).event(eloc).icc.m...
            %            - relsummary.group(gloc).event(eloc).icc.ll;
            %        ulest(eloc,gloc) = relsummary.group(gloc).event(eloc).icc.ul...
            %            - relsummary.group(gloc).event(eloc).icc.m;
            
            %pull stat-specific information
            switch lower(stat)
                case 'icc'
                    llest(eloc,gloc) = ...
                        psyrat_data.relsummary.group(gloc).event(eloc).icc.m -...
                        psyrat_data.relsummary.group(gloc).event(eloc).icc.ll;
                    ptest(eloc,gloc) = psyrat_data.relsummary.group(gloc).event(eloc).icc.m;
                    ulest(eloc,gloc) = ...
                        psyrat_data.relsummary.group(gloc).event(eloc).icc.ul -...
                        psyrat_data.relsummary.group(gloc).event(eloc).icc.m;
                case 'bet'
                    llest(eloc,gloc) = ...
                        psyrat_data.relsummary.group(gloc).event(eloc).betsd.m -...
                        psyrat_data.relsummary.group(gloc).event(eloc).betsd.ll;
                    ptest(eloc,gloc) = psyrat_data.relsummary.group(gloc).event(eloc).betsd.m;
                    ulest(eloc,gloc) = ...
                        psyrat_data.relsummary.group(gloc).event(eloc).betsd.ul -...
                        psyrat_data.relsummary.group(gloc).event(eloc).betsd.m;
                case 'wit'
                    llest(eloc,gloc) = ...
                        psyrat_data.relsummary.group(gloc).event(eloc).witsd.m -...
                        psyrat_data.relsummary.group(gloc).event(eloc).witsd.ll;
                    ptest(eloc,gloc) = psyrat_data.relsummary.group(gloc).event(eloc).witsd.m;
                    ulest(eloc,gloc) = ...
                        psyrat_data.relsummary.group(gloc).event(eloc).witsd.ul -...
                        psyrat_data.relsummary.group(gloc).event(eloc).witsd.m;
            end
            
            %figure out spacing for plot
            if gloc < median(1:ngroups)
                offsetm(eloc,gloc) = eloc - (.4/ngroups);
            elseif gloc > median(1:ngroups)
                offsetm(eloc,gloc) = eloc + (.4/ngroups);
            elseif gloc == median(1:ngroups)
                offsetm(eloc,gloc) = eloc;
            end
            
        end
    end
    
else
    %grab icc information for each event and group
    ciedge = (1-ciperc)/2;
    for gloc=1:ngroups
        for eloc=1:nevents
           
            %pull stat-specific information
            
            switch lower(stat)
                case 'bet'
                    %gro_sds already holds sigma_p as a STANDARD DEVIATION, so
                    %summarize the draws directly. Taking sqrt() here would
                    %report sigma_p^(1/2), which is not a physical quantity.
                    %Matches psyrat_variancet's sserr branch and the 'wit' case
                    %below, both of which summarize on the SD scale.
                    bp_sd = cell2mat(psyrat_data.relsummary.data.g(gloc).e(eloc).gro_sds(:,1));

                    ptest(eloc,gloc) = mean(bp_sd);
                    llest(eloc,gloc) = mean(bp_sd) - ...
                        quantile(bp_sd,ciedge);
                    ulest(eloc,gloc) = quantile(bp_sd,1-ciedge) - ...
                        mean(bp_sd);

                case 'wit'
                    
                    pop_sdlog = psyrat_data.relsummary.data.g(gloc).e(eloc).pop_sdlog;
                    
                    ptest(eloc,gloc) = exp(mean(pop_sdlog));
                    llest(eloc,gloc) = exp(mean(pop_sdlog)) - ...
                        exp(quantile(pop_sdlog,ciedge));
                    ulest(eloc,gloc) = exp(quantile(pop_sdlog,1-ciedge)) - ...
                        exp(mean(pop_sdlog));
            end
            
            %figure out spacing for plot
            if gloc < median(1:ngroups)
                offsetm(eloc,gloc) = eloc - (.4/ngroups);
            elseif gloc > median(1:ngroups)
                offsetm(eloc,gloc) = eloc + (.4/ngroups);
            elseif gloc == median(1:ngroups)
                offsetm(eloc,gloc) = eloc;
            end
            
        end
    end
end

%find the dimensions
[e,g] = size(ptest);

%plot
if ~(g > 1 && e == 1)
    ptintplot = errorbar(offsetm,ptest,llest,ulest,'Marker','.',...
        'MarkerSize',15,'LineWidth',1);
elseif g > 1 && e == 1
    hold on
    for i = 1:g
        ptintplot = errorbar(offsetm(i),ptest(i),llest(i),ulest(i),...
            'Marker','.','MarkerSize',15,'LineWidth',1,...
            'DisplayName',gnames{i});
        legend('-DynamicLegend');
    end
    hold off
end

%add plot-specific title and y-labels
switch lower(stat)
    case 'icc'
        set(gcf,'Name','ICC Estimates');
        ylabel('Intraclass Correlation Coefficient','FontSize',fsize);
    case 'bet'
        set(gcf,'Name','Between-Person Standard Deviations');
        ylabel('Between-Person Standard Deviations','FontSize',fsize);
    case 'wit'
        set(gcf,'Name','Within-Person Standard Deviations');
        ylabel('Within-Person Standard Deviations','FontSize',fsize);
end

%remove extra lines
if nevents > 1
    for i = 1:length(ptintplot)
        ptintplot(i).LineStyle = 'none';
    end
end

%fix axes
ptintplot(1).Parent.XLim = [0 nevents+1+.25];
if nevents == 1
    ptintplot(1).Parent.YLim = [0 max(ulest+ptest)+.05];
elseif nevents > 1
    ptintplot(1).Parent.YLim = [0 max(max(ulest+ptest))+.05];
    xlabel('Event','FontSize',fsize);
end

ptintplot(1).Parent.XTick = 1:nevents;
ptintplot(1).Parent.XTickLabel = enames;
%Interpreter none on the tick labels and (below) the legend: an underscored
%event or group name (inc_cor) must print literally, where the default TeX
%interpreter would subscript the character after the underscore.
ptintplot(1).Parent.TickLabelInterpreter = 'none';
ptintplot(1).Parent.FontSize = fsize;



%add names to legend
if ~(g > 1 && e == 1)
    pl = legend(ptintplot);
    pl.String = gnames;
end

%both legend branches above leave exactly one legend on the axes; the axes
%Legend property reaches it without rebuilding it
ptintplot(1).Parent.Legend.Interpreter = 'none';
legend('boxoff');

%change axis location and rotate plot
ptintplot(1).Parent.YAxisLocation = 'right';
camroll(-90);


end
