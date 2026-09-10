function depplot = psyrat_depvtrialsplot(varargin)
%Plot the dependability as a function of the number of trials included in
%an average (stratified by group and condition, if applicable)
%
%depplot = psyrat_depvtrialsplot('psyrat_data',psyrat_data,'trials',[1 50],...
%   'depline',plotdepline,'depcutoff',depcutoff);
%
%
%Inputs
% psyrat_data - PsyRAT Toolbox data structure array. Variance components should
%  be included.
% trials - range of trials to plot
% depcutoff - dependality threshold for deeming data reliable (used for
%  reference line in plot)
%
%Optional Input
% depline - dependability estimate to plot (default: 2)
%  1 - lower limit
%  2 - point estimate
%  3 - upper limit
% gcoeff - one-facet coefficient to use
%  1 - dependability/reliability (absolute)
%  2 - generalizability (relative)
% CI - confidence interval width. Decimal from 0 to 1. (default: .95)
%
%
%Outputs
% depplot - figure handle for the plot
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
    
    %check if trials was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('trials',varargin),1);
    if ~isempty(ind)
        trials = varargin{ind+1}; 
    else 
        error('varargin:trials',... %Error code and associated error
            strcat('WARNING: Range of trials to plot not specified \n\n',... 
            'Please input trials: [min max]\n',...
            'See help psyrat_depvtrialsplot for more information\n'));
    end
    
    %check if depcutoff was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('depcutoff',varargin),1);
    if ~isempty(ind)
        depcutoff = varargin{ind+1}; 
    else 
        error('varargin:depcutoff',... %Error code and associated error
            strcat('WARNING: Dependanility threshold not specified \n\n',... 
            'Please input depcutoff\n',...
            'See help psyrat_depvtrialsplot for more information\n'));
    end
    
    %check if depline was specified. 
    %If it is not found, set as default: 2.
    ind = find(strcmpi('depline',varargin),1);
    if ~isempty(ind)
        depline = varargin{ind+1}; 
        %make sure depline is 1, 2, or 3
        if ~any(depline==[1 2 3])
            error('varargin:depline',... %Error code and associated error
                strcat('WARNING: Dependability estimate to plot not specified \n\n',... 
                'Please enter a valid input for depline\n',...
                '1 - lower limit of credible interval\n',...
                '2 - point estimate of credible interval\n',...
                '3 - upper limit of credible interval\n',...
                'See help psyrat_depvtrialsplot for more information \n'));
        end
    else 
        depline = 2; %default is point estimate
    end
    
    %check if gcoeff was specified.
    %If it is not found, set as default: 1 (dependability).
    ind = find(strcmpi('gcoeff',varargin),1);
    if ~isempty(ind)
        gcoeff = varargin{ind+1};
        if ~any(gcoeff == [1 2])
            error('varargin:gcoeff',... %Error code and associated error
                strcat('WARNING: gcoeff is invalid. Valid values are 1 or 2\n'));
        end
    elseif isfield(psyrat_data,'relsummary') && isfield(psyrat_data.relsummary,'gcoeff')
        gcoeff = psyrat_data.relsummary.gcoeff;
    else
        gcoeff = 1;
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

%create an x-axis for the number of observations
x = trials(1):trials(2);

%create an empty array for storing information into. RC-30. One row per
%requested trial count, i.e. numel(x) = trials(2)-trials(1)+1 rows, NOT
%trials(2) rows. The rows are indexed relative to the start of the range
%(row 1 is trials(1)), matching the convention psyrat_rel_sing already
%documents for its return vectors at psyrat_rel_sing.m:167-169 ("index
%relative to the range start so outputs align to obs(1):obs(2)", then
%idx = ii - obs(1) + 1); the trt_diff branch below reads psyrat_rel_trt,
%which uses the same convention at psyrat_rel_trt.m:425-428. Note both cites
%are in files this change does not touch, so neither can be shifted by it.
%Indexing this buffer by the ABSOLUTE trial number
%instead grew it to trials(2) rows while x had only numel(x) elements, so
%plot(x,plotrel) threw MATLAB:samelen for any trials(1) > 1 - a documented
%input that did not work. Under every production caller (trials(1) == 1) the
%two indexings are identical, so this is a strict no-op there.
%Kept as a deliberate twin of psyrat_trt_relvtrialsplot; change both together.
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
depplot = figure;
depplot.Tag = 'psyrat_output';
set(gcf,'NumberTitle','Off');
depplot.Position = [125 630 900 450];
fsize = 16;


%extract the data and create the subplots for depplot
for eloc=1:nevents
    for gloc=1:ngroups
        
        %compute the PER-EVENT reliability curve over the trial range. The
        %one-facet designs use psyrat_rel_sing; trt_diff data use psyrat_rel_trt
        %with the per-condition components (the diagonals of the cross-condition
        %covariances stored as SD draws) at the selected occasion estimand.
        %
        %RC-31. Every branch here is PER-CONDITION, on every analysis this
        %function draws, and this comment used to call the trt_diff branch "the
        %two-facet difference design", which reads as though the curve were a
        %difference score. It is not: psyrat_rel_trt is the plain two-facet
        %kernel and the components are one condition's. That reading is what made
        %it look reasonable to feed this figure diffgcoeff. It follows gcoeff.
        %The genuine difference-score curve is psyrat_trt_diffvtrialsplot, drawn
        %beside this one by psyrat_relfigures.
        if strcmp(psyrat_data.rel.analysis,'trt_diff')
            d = psyrat_data.relsummary.data.g(gloc).e(eloc);
            reltype = psyrat_data.relsummary.reltype;
            if isfield(psyrat_data.relsummary,'nocc') && ~isempty(psyrat_data.relsummary.nocc)
                nocc = psyrat_data.relsummary.nocc;
            else
                nocc = 1;
            end
            [llv,mv,ulv] = psyrat_rel_trt('gcoeff',gcoeff,'reltype',reltype,...
                'bp',d.sig_id.raw,'bo',d.sig_occ.raw,'bt',d.sig_trl.raw,...
                'txp',d.sig_trlxid.raw,'oxp',d.sig_occxid.raw,...
                'txo',d.sig_trlxocc.raw,'err',d.sig_err.raw,...
                'obs',[trials(1) trials(2)],'nocc',nocc,'CI',ciperc);
        else
            d = psyrat_data.relsummary.data.g(gloc).e(eloc);
            if strcmp(psyrat_data.rel.analysis,'ic')
                bp = d.sig_u.raw;
                wp = d.sig_e.raw;
            elseif any(strcmp(psyrat_data.rel.analysis,{'ic_diff','ic_diff_sserrvar'}))
                bp = cell2mat(d.sd_id.raw);
                wp = exp(cell2mat(d.b_sigma.raw));
            end
            %wp above is the person x trial,e residual only. psyrat_rel_sing
            %expects the TOTAL within-person SD, with the trial main effect
            %passed separately so it can be subtracted back out for the relative
            %coefficient. Widening wp and passing 'i' must happen together:
            %either alone gives a curve that is wrong in a different way, and
            %without both this plot draws the relative coefficient under a
            %"Dependability" axis (RC-02).
            sigi = local_sigtrl(d, wp);
            [llv,mv,ulv] = psyrat_rel_sing('gcoeff',gcoeff,'metric','global',...
                'bp',bp,'wp',sqrt(wp.^2 + sigi.^2),'i',sigi,...
                'obs',[trials(1) trials(2)],'CI',ciperc);
        end

        switch depline
            case 1 %lower limit
                plotrel(1:numel(x),gloc) = llv;
                plottitle = 'Lower Limit of 95% Credible Interval';
            case 2 %point estimate
                plotrel(1:numel(x),gloc) = mv;
                plottitle = 'Point Estimate';
            case 3 %upper limit
                plotrel(1:numel(x),gloc) = ulv;
                plottitle = 'Upper Limit of 95% Credible Interval';
        end
    end
    
    if gcoeff == 1
        pref = ['Dependability v ' L.NumberOf ': '];
        ylbl = 'Dependability';
    else
        pref = ['Generalizability v ' L.NumberOf ': '];
        ylbl = 'Generalizability';
    end
    set(gcf,'Name',[pref plottitle]);
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
    
    ylabel(ylbl,'FontSize',fsize);
    xlabel('Number of Observations','FontSize',fsize);
    psyrat_addhline(depcutoff,'Color','b','LineStyle',':');
    %only build a named legend when there are multiple groups to label.
    %For single-condition analyses (one group) there are no series names,
    %so calling legend('boxoff') unconditionally would auto-create MATLAB's
    %default "data1"/"data2" legend. Guard it here so those plots show no
    %legend rather than a meaningless one.
    if analysis ~= 1 && analysis ~= 3
        leg = legend(gnames{:},'Location','southeast');
        set(leg,'FontSize',fsize,'Interpreter','none');
        legend('boxoff');
    end

end


end

function sigi = local_sigtrl(d, wp_resid)
%Trial main effect (sigma_i) draws for one one-facet stratum, returned as an SD
%vector shaped like the residual SD it will be combined with.
%
%The component is stored under a different field name per analysis, because each
%branch of psyrat_relsummary builds its own data struct: the one-facet model
%stores an SD directly (sig_trl); the single-session difference stores a per-
%event SD (sd_trl); the subject-level difference stores the trial covariance and
%the per-event value is its diagonal, a VARIANCE, hence the sqrt. Checked in that
%order because the difference structs carry trl_varcov as well, and sd_trl is the
%unambiguous per-event SD.
%
%Legacy result structs produced before the trial main effect was modeled carry
%none of these; those return zeros, which reduces the coefficient to the prior
%behavior (sigma_i = 0, dependability == generalizability) rather than erroring.
%
%Kept local rather than shared with the identical helper in psyrat_report because
%the two files are otherwise independent and the toolbox convention is private
%locals per surface (compare local_total_within in psyrat_relsummary). Noted in
%the call-site census so a future consolidation is a decision, not a discovery.
sigi = zeros(size(wp_resid));
if isfield(d,'sig_trl') && isfield(d.sig_trl,'raw') && ~isempty(d.sig_trl.raw)
    sigi = d.sig_trl.raw;
elseif isfield(d,'sd_trl') && isfield(d.sd_trl,'raw') && ~isempty(d.sd_trl.raw)
    sigi = cell2mat(d.sd_trl.raw);
elseif isfield(d,'trl_varcov') && isfield(d.trl_varcov,'raw') && ...
        ~isempty(d.trl_varcov.raw)
    sigi = sqrt(max(cell2mat(d.trl_varcov.raw),0));
end
sigi = sigi(:);
end
