function psyrat_ssrelplot(varargin)
%Plot subject-level reliability stratified by group/event, if applicable
%
%psyrat_ssrelplot('psyrat_data',psyrat_data,'stat','dep');
%
%
%Inputs
% psyrat_data - PsyRAT Toolbox data structure array. Variance components should
%  be included.
% stat - 'dep': dependability, 'icc': intraclass correlation coefficient
%
%Optional Input
% CI - confidence interval width. Decimal from 0 to 1. (default: .95)
%
%
%Outputs
% figure that shows subject-level reliability in ascending order

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
            'See help psyrat_ssrelplot for more information about inputs'));
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
            'See help psyrat_ssrelplot for more information \n'));
    end
    
    
    %check whether 'stat' is defined
    %If it is not found, set display error.
    ind = find(strcmpi('stat',varargin),1);
    if ~isempty(ind)
        stat = varargin{ind+1};
    else
        error('varargin:stat_missing',... %Error code and associated error
            strcat('WARNING: stat not specified \n\n',...
            'Please input coefficient type ''dep'' or ''icc''.\n',...
            'See help psyrat_ssrelplot for more information \n'));
    end
    
end


%make sure the user understood the the CI input is width not edges, if
%inputted incorrectly, provide a warning, but don't change
ciperc = psyrat_data.relsummary.ciperc;
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

%The gate tests group(1) while the reads below index group(gloc), and that is
%safe rather than lucky: a MATLAB struct array carries ONE field set shared by
%every element, so isfield(group(1),'diffscore') and isfield(group(gloc),...)
%cannot disagree for any gloc. Element 1 is a valid proxy by construction, not
%by convention, so there is no per-group divergence to guard against here
%(RC-41). The gate is also computed once, outside both gloc loops, so it could
%not take gloc without being moved into them.
%
%What a struct array does NOT promise is a non-empty VALUE at every element. No
%current writer can leave one empty: the four assignments to
%group(gloc).diffscore -- two in psyrat_relsummary's own sing_diff branch (the
%diffscore and its trlcutoff row) and one each in the helpers
%psyrat_relsummary_sing_diff_sserr and psyrat_relsummary_trt_diff -- all run
%unconditionally inside their gloc loops with no continue or early return. So
%this is not guarded, only recorded. Cited by function rather than by line
%number: the line numbers moved twice while this comment was being written, so
%grep "group(gloc).diffscore =" instead. The sibling that DOES guard the empty
%value, and the form to copy if that ever changes, is
%psyrat_admissibility_collect (its isfield-or-isempty continue).
plot_diffscore = strcmpi(psyrat_data.rel.analysis,'ic_diff_sserrvar') && ...
    isfield(psyrat_data,'relsummary') && isfield(psyrat_data.relsummary,'group') && ...
    isfield(psyrat_data.relsummary.group(1),'diffscore');
if plot_diffscore
    nevents = 1;
    enames = {'diff score'};
end

%figure out whether groups or events need to be considered
%1 - no groups or event types to consider
%2 - possible multiple groups but no event types to consider
%3 - possible event types but no groups to consider
%4 - possible groups and event types to consider

%see how many subplots are needed
if (nevents * ngroups) > 2
    xplots = ceil(sqrt((nevents * ngroups)));
    yplots = ceil(sqrt((nevents * ngroups)));
elseif (nevents * ngroups) == 2
    xplots = 1;
    yplots = 2;
elseif (nevents * ngroups) == 1
    xplots = 1;
    yplots = 1;
end

%create the figure
depplot = figure;
depplot.Tag = 'psyrat_output';
set(gcf,'NumberTitle','Off');
depplot.Position = [125 630 900 450];
fsize = 16;
plot_count = 1;

%The dep_* columns hold the selected subject-level coefficient: the absolute
%dependability phi_s (gcoeff=1) or the relative generalizability G_s
%(gcoeff=2). Label and the group reference line follow that choice; pop_errvar
%is the relative residual, so the absolute reference adds the trial main effect
%sigma_i^2 (trl_var).
if isfield(psyrat_data.relsummary,'gcoeff_name') && ...
        strcmpi(psyrat_data.relsummary.gcoeff_name,'gen')
    coeff_label = 'Generalizability';
    gcoeff_abs = false;
else
    coeff_label = 'Dependability';
    gcoeff_abs = true;
end

if strcmp(stat,'dep')
    plottitle = ['Subject-Level ' coeff_label ' Estimates'];
    
    %extract the data and create the subplots for depplot
    for eloc=1:nevents
        %ngroups comes from rel.groups above and the summary is grown to match
        %it; the invariant is stated at psyrat_relsummary's publish site (RC-44).
        for gloc=1:ngroups
            if plot_diffscore
                gtable = psyrat_data.relsummary.group(gloc).diffscore.ssrel_table;
            else
                gtable = psyrat_data.relsummary.group(gloc).event(eloc).ssrel_table;
            end
            gtable = sortrows(gtable,'dep_pt');
            
            plotssrel.subj = 1:height(gtable);
            plotssrel.dep_pt = gtable.dep_pt;
            plotssrel.dep_ll = gtable.dep_pt - gtable.dep_ll;
            plotssrel.dep_ul = gtable.dep_ul - gtable.dep_pt;
            
            %Two-facet (test-retest) tables carry the group reference line as a
            %column, because the two-facet partition needs bo/bt/txp/oxp/txo and
            %none of those reach this function. The bp_var/pop_errvar expression
            %below is the ONE-FACET reference line: correct for psyrat_ssrel
            %tables, but a different quantity from the two-facet dep_* bars it
            %would otherwise be plotted against. Both branches use the same
            %mean(trls) divisor, so this changes the partition, not the n'.
            if any(strcmp(gtable.Properties.VariableNames,'gro_rel'))
                gro_est = gtable.gro_rel(1);
            else
                gro_err = gtable.pop_errvar(1);
                if gcoeff_abs && ...
                        any(strcmp(gtable.Properties.VariableNames,'trl_var'))
                    gro_err = gro_err + gtable.trl_var(1);
                end
                gro_est = gtable.bp_var(1) / ...
                    (gtable.bp_var(1) + (gro_err/mean(gtable.trls)));
            end

            miss = 0;
            miss_ids = [];

            if any(gtable.dep_ul < gro_est) ||...
                    any(gtable.dep_ll > gro_est)
                
                miss = 1;
                miss_ids1 = find(gtable.dep_ul < gro_est);
                miss_ids2 = find(gtable.dep_ll > gro_est);
                miss_ids = [miss_ids1; miss_ids2];
                
                good_ids = 1:height(gtable);
                good_ids(miss_ids) = [];
            end
            
            
            set(gcf,'Name',plottitle);
            subplot(yplots,xplots,plot_count);
            if miss
                hold on
                errorbar(plotssrel.subj(miss_ids),...
                    plotssrel.dep_pt(miss_ids),...
                    plotssrel.dep_ll(miss_ids),...
                    plotssrel.dep_ul(miss_ids),...
                    'Marker','.',...
                    'MarkerSize',25,...
                    'LineWidth',1.5,...
                    'LineStyle','none',...
                    'MarkerEdgeColor','r',...
                    'MarkerFaceColor','r',...
                    'Color','r');
                
                errorbar(plotssrel.subj(good_ids),...
                    plotssrel.dep_pt(good_ids),...
                    plotssrel.dep_ll(good_ids),...
                    plotssrel.dep_ul(good_ids),...
                    'Marker','.',...
                    'MarkerSize',25,...
                    'LineWidth',1.5,...
                    'LineStyle','none',...
                    'MarkerEdgeColor','b',...
                    'MarkerFaceColor','b',...
                    'Color','b');
                
                hold off
            else
                errorbar(plotssrel.subj,...
                    plotssrel.dep_pt,...
                    plotssrel.dep_ll,...
                    plotssrel.dep_ul,...
                    'Marker','.',...
                    'MarkerSize',25,...
                    'LineWidth',1.5,...
                    'LineStyle','none');
                
            end
            axis([0 max(plotssrel.subj)+1 0 1]);
            set(gca,'fontsize',16);
            set(gca,'xticklabel',[])
            
            if plot_diffscore && ngroups > 1
                %Interpreter none on every named title in this block: an
                %underscored group or event name (inc_cor) must print
                %literally, where TeX would subscript the character after it.
                title([gnames{gloc} ': diff score'],'FontSize',20,'Interpreter','none');
            elseif plot_diffscore
                title('diff score','FontSize',20);
            elseif (~strcmpi(enames{eloc},'none') &&...
                    ~strcmpi(gnames{gloc},'none')) &&...
                    nevents > 1 && ngroups > 1
                title([gnames{gloc} ': ' enames{eloc}],'FontSize',20,'Interpreter','none');
            elseif (isempty(gnames{gloc}) || strcmpi(gnames{gloc},'none')) ...
                    && nevents > 1
                %events-only run: the group facet is absent (gnames is {''}
                %above, or 'none' if a caller passes the loader's spelling),
                %so the panel is identified by its event name alone.
                title(enames{eloc},'FontSize',20,'Interpreter','none');
            elseif (isempty(enames{eloc}) || strcmpi(enames{eloc},'none')) ...
                    && ngroups > 1
                %groups-only run: the event facet is absent, so the panel is
                %identified by its group name alone. (Before the pre-beta fix
                %both branches tested the name they were about to PRINT
                %against 'none', which the loader never sets, so neither could
                %fire and every panel on an events-only or groups-only run was
                %left untitled.)
                title(gnames{gloc},'FontSize',20,'Interpreter','none');
            else
                title('');
            end
            
            ylabel(coeff_label,'FontSize',fsize);
            xlabel(['Participants Ordered by ' coeff_label ' Estimate'],...
                'FontSize',fsize);
            psyrat_addhline(gro_est,...
                'Color','k',...
                'LineStyle','-',...
                'LineWidth',1.5);
            
            plot_count = plot_count + 1;
        end
    end
elseif strcmp(stat,'icc')
    plottitle = 'Subject-Level ICCs';
    
    %extract the data and create the subplots for depplot
    for eloc=1:nevents
        for gloc=1:ngroups
            if plot_diffscore
                gtable = psyrat_data.relsummary.group(gloc).diffscore.ssrel_table;
            else
                gtable = psyrat_data.relsummary.group(gloc).event(eloc).ssrel_table;
            end
            gtable = sortrows(gtable,'icc_pt');
            
            plotssrel.subj = 1:height(gtable);
            plotssrel.icc_pt = gtable.icc_pt;
            plotssrel.icc_ll = gtable.icc_pt - gtable.icc_ll;
            plotssrel.icc_ul = gtable.icc_ul - gtable.icc_pt;
            
            %As in the dependability panel above: two-facet tables carry the
            %group ICC as a column (the single-observation coefficient under the
            %same partition as the icc_* bars), and the undivided
            %bp_var/pop_errvar expression is the one-facet fallback.
            if any(strcmp(gtable.Properties.VariableNames,'gro_icc'))
                gro_est = gtable.gro_icc(1);
            else
                gro_err = gtable.pop_errvar(1);
                if gcoeff_abs && ...
                        any(strcmp(gtable.Properties.VariableNames,'trl_var'))
                    gro_err = gro_err + gtable.trl_var(1);
                end
                gro_est = gtable.bp_var(1) / ...
                    (gtable.bp_var(1) + gro_err);
            end

            miss = 0;
            miss_ids = [];

            if any(gtable.icc_ul < gro_est) ||...
                    any(gtable.icc_ll > gro_est)
                
                miss = 1;
                miss_ids1 = find(gtable.icc_ul < gro_est);
                miss_ids2 = find(gtable.icc_ll > gro_est);
                miss_ids = [miss_ids1; miss_ids2];
                
                good_ids = 1:height(gtable);
                good_ids(miss_ids) = [];
            end
            
            
            set(gcf,'Name',plottitle);
            subplot(yplots,xplots,plot_count);
            if miss
                hold on
                errorbar(plotssrel.subj(miss_ids),...
                    plotssrel.icc_pt(miss_ids),...
                    plotssrel.icc_ll(miss_ids),...
                    plotssrel.icc_ul(miss_ids),...
                    'Marker','.',...
                    'MarkerSize',25,...
                    'LineWidth',1.5,...
                    'LineStyle','none',...
                    'MarkerEdgeColor','r',...
                    'MarkerFaceColor','r',...
                    'Color','r');
                
                errorbar(plotssrel.subj(good_ids),...
                    plotssrel.icc_pt(good_ids),...
                    plotssrel.icc_ll(good_ids),...
                    plotssrel.icc_ul(good_ids),...
                    'Marker','.',...
                    'MarkerSize',25,...
                    'LineWidth',1.5,...
                    'LineStyle','none',...
                    'MarkerEdgeColor','b',...
                    'MarkerFaceColor','b',...
                    'Color','b');
                
                hold off
            else
                errorbar(plotssrel.subj,...
                    plotssrel.icc_pt,...
                    plotssrel.icc_ll,...
                    plotssrel.icc_ul,...
                    'Marker','.',...
                    'MarkerSize',25,...
                    'LineWidth',1.5,...
                    'LineStyle','none');
                
            end
            axis([0 max(plotssrel.subj)+1 0 1]);
            set(gca,'fontsize',16);
            set(gca,'xticklabel',[])
            
            if plot_diffscore && ngroups > 1
                %Interpreter none on every named title in this block: an
                %underscored group or event name (inc_cor) must print
                %literally, where TeX would subscript the character after it.
                title([gnames{gloc} ': diff score'],'FontSize',20,'Interpreter','none');
            elseif plot_diffscore
                title('diff score','FontSize',20);
            elseif (~strcmpi(enames{eloc},'none') &&...
                    ~strcmpi(gnames{gloc},'none')) &&...
                    nevents > 1 && ngroups > 1
                title([gnames{gloc} ': ' enames{eloc}],'FontSize',20,'Interpreter','none');
            elseif (isempty(gnames{gloc}) || strcmpi(gnames{gloc},'none')) ...
                    && nevents > 1
                %events-only run: the group facet is absent (gnames is {''}
                %above, or 'none' if a caller passes the loader's spelling),
                %so the panel is identified by its event name alone.
                title(enames{eloc},'FontSize',20,'Interpreter','none');
            elseif (isempty(enames{eloc}) || strcmpi(enames{eloc},'none')) ...
                    && ngroups > 1
                %groups-only run: the event facet is absent, so the panel is
                %identified by its group name alone. (Before the pre-beta fix
                %both branches tested the name they were about to PRINT
                %against 'none', which the loader never sets, so neither could
                %fire and every panel on an events-only or groups-only run was
                %left untitled.)
                title(gnames{gloc},'FontSize',20,'Interpreter','none');
            else
                title('');
            end
            
            ylabel('ICC','FontSize',fsize);
            xlabel('Participants Ordered by ICC',...
                'FontSize',fsize);
            psyrat_addhline(gro_est,...
                'Color','k',...
                'LineStyle','-',...
                'LineWidth',1.5);
            
            plot_count = plot_count + 1;
        end
    end
end

end
