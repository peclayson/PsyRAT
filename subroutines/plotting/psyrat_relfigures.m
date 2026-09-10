function psyrat_relfigures(varargin)
%Creates various figures and tables for dependability data.
%
%psyrat_relfigures('psyrat_data',psyrat_data,'analysis','sing')
%
%
%Required Inputs:
% psyrat_data - PsyRAT Toolbox data structure array containing outputs from
%  CmdStan
% analysis - 'sing' for single session data, 'trt' for data with multiple
%  occasions
%
% Option 1:
%  psyrat_prefs - can contain all of the preferences for plotting figures and
%   tables
%
% Option 2 (for analysis = 'sing'):
%  depcutoff - dependability level to use for cutoff when deciding the
%   minimum number of trials needed to achieve this specified level of
%   dependability
%  plotdep - plot the table displaying dependability v number of trials
%   included in average
%  ploticc - plot the intraclass correlation coefficients for data. This ICC
%   can be interpreted as the proportion of total variance that is between
%   persons
%  showinct - display table with information for each event/group
%   combination (cutoffs and dependability)
%  showoverallt - display table with information for data including all
%   trials for those participants that meet cutoff threshhold
%  showstddevt - display table with information about the sources of
%   variance (between person v within person)
%  plotbetstddev - plot the between-person standard deviations stratified
%   by group and event
%  plotdepline - indicate whether to plot 1-lower limit of credible
%   interval, 2-point estimate, 3-upper limit of credible interval for the
%   plot of dependability v number of trials (default: 2)
%  plotntrials - indicate the number of trials to plot (x-axis) in the
%   dependability v number of trials plot (default: 50)
%  meascutoff - which estimate to use to define cutoff for number of trials
%   1 - lower limit of credible interval, 2 - point estimate, 3 - upper
%   limit of credible interval (default: 2)
%  depcentmeas - which measure of central tendency to use to estimate the
%   overall dependability, 1 - mean, 2 - median (default: 1)
%
% Option 3 (for analysis = 'ic_sing_sserr')
%  plotdep - plots the subject-level dependability estimates with credible
%   intervals in ascending order
%  ploticc- plots the subject-level ICCs with credible intervals in
%   ascending order
%
% Option 4 (for analysis = 'trt'):
%  relcutoff - reliability level to use for cutoff when deciding the
%   minimum number of trials needed to achieve this specified level of
%   reliability
%  plotrel - plot the table displaying reliability v number of trials
%   included in average
%  ploticc - plot the intraclass correlation coefficients for data. This ICC
%   can be interpreted as the proportion of total variance that is between
%   persons
%  showinct - display table with information for each event/group
%   combination (cutoffs and reliability)
%  showoverallt - display table with information for data including all
%   trials for those participants that meet cutoff threshhold
%  showstddevt - display table with information about the sources of
%   variance (between person v within person)
%  plotbetstddev - plot the between-person standard deviations stratified
%   by group and event
%  plotrelline - indicate whether to plot 1-lower limit of credible
%   interval, 2-point estimate, 3-upper limit of credible interval for the
%   plot of reliability v number of trials (default: 2)
%  plotntrials - indicate the number of trials to plot (x-axis) in the
%   reliability v number of trials plot (default: 50)
%  meascutoff - which estimate to use to define cutoff for number of trials
%   1 - lower limit of credible interval, 2 - point estimate, 3 - upper
%   limit of credible interval (default: 2)
%  relcentmeas - which measure of central tendency to use to estimate the
%   overall reliability, 1 - mean, 2 - median (default: 1)
%  gcoeff - g-theory coefficient to calculate 1 - dependability, 2 -
%   generalizability
%  reltype - reliability coefficient to plot/calculate - 1 - coefficent of
%   equivalence (internal consistency), 2 - coefficient of stability (test
%   retest reliability), 3 - coefficient of equivalence and stability
%
%
%Output:
% Various figures may be plotted. Tables will also have buttons for
%  saving the table to a file. Additionally the table showing the cutoffs
%  for numbers of trials to achieve a given level of reliability will
%  have a button to save the ids of participants to include and exclude.

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

%analysis-unit label fragments (trial vs split) for the cutoff error dialogs;
%defaults to "trial" so non-splits output is byte-identical.
usplit = 1;
if isfield(psyrat_data,'rel') && isfield(psyrat_data.rel,'splits')
    usplit = psyrat_data.rel.splits;
end
L = psyrat_unitlabels(usplit);

%check whether psyrat_data was provided
%If it is not found, set display error.
if isempty(psyrat_data)
    error('varargin:nofile',... %Error code and associated error
        strcat('WARNING: psyrat_data not specified \n\n',...
        'Please input the psyrat_data to be loaded \n'));
end

%somersault through inputs to find analysis type
if ~isempty(varargin)
    %the optional inputs check assumes that there was an even number of
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)
        error('varargin:incomplete',... %Error code and associated error
            strcat('WARNING: Inputs are incomplete \n\n',...
            'Make sure each variable input is paired with a value \n',...
            'See help psyrat_relfigures for more information on optional inputs'));
    end
    
    %check if the dependability cutoff was specified
    ind = find(strcmp('analysis',varargin),1);
    if ~isempty(ind)
        analysis = varargin{ind+1};
    else
        error('varargin:noanalysistype',... %Error code and associated error
            strcat('WARNING: analysis not specified \n\n',...
            'Please input the analysis to be run:\n',...
            ' ''sing'', ''sing_sserr'', or ''trt'' \n'));
    end
end

%somersault through inputs to find preferences for running data, if
%psyrat_prefs was not provided as input
if ~isempty(varargin) && isempty(psyrat_prefs)
    %check if the dependability cutoff was specified
    ind = find(strcmp('depcutoff',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            depcutoff = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            depcutoff = varargin{ind+1};
        end
    else
        depcutoff = .80; %default level is .80
    end
    
    %check if plotdep is provided
    ind = find(strcmp('plotdep',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            pdep = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            pdep = varargin{ind+1};
        end
    else
        pdep = 1; %default is 1
    end
    
    %check if picc is provided
    ind = find(strcmp('ploticc',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            picc = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            picc = varargin{ind+1};
        end
    else
        picc = 1; %default is 1
    end
    
    %check if showinct is provided
    ind = find(strcmp('showinct',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            showinct = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            showinct = varargin{ind+1};
        end
    else
        showinct = 1; %default is 1
    end
    
    %check if showoverallt is provided
    ind = find(strcmp('showoverallt',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            showoverallt = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            showoverallt = varargin{ind+1};
        end
    else
        showoverallt = 1; %default is 1
    end
    
    %check if plotntrials is provided
    ind = find(strcmp('plotntrials',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            plotntrials = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            plotntrials = varargin{ind+1};
        end
    else
        plotntrials = 50; %default is 50
    end
    
    %check if showstddevt is provided
    ind = find(strcmp('showstddevt',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            showstddevt = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            showstddevt = varargin{ind+1};
        end
    else
        showstddevt = 1; %default is 1
    end
    
    %check if plotbetstddev is provided
    ind = find(strcmp('plotbetstddev',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            plotbetstddev = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            plotbetstddev = varargin{ind+1};
        end
    else
        plotbetstddev = 1; %default is 1
    end
    
    %check if plotwitstddev is provided
    ind = find(strcmp('plotwitstddev',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            plotwitstddev = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            plotwitstddev = varargin{ind+1};
        end
    else
        plotwitstddev = 0; %default is 0
    end
    
    %check if plotdepline is provided
    ind = find(strcmp('plotdepline',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            plotdepline = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            plotdepline = varargin{ind+1};
        end
    else
        plotdepline = 2; %default is 2 (point estimate), matching psyrat_defaults
    end
    
    %check if meascutoff is provided
    ind = find(strcmp('meascutoff',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            meascutoff = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            meascutoff = varargin{ind+1};
        end
    else
        meascutoff = 2; %default is 2 (point estimate of credible interval)
    end
    
    %check if depcentmeas is provided
    ind = find(strcmp('depcentmeas',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            depcentmeas = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            depcentmeas = varargin{ind+1};
        end
    else
        depcentmeas = 1; %default is 1 (mean)
    end

    %check if diffgcoeff is provided (single-session difference scores)
    ind = find(strcmp('diffgcoeff',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            diffgcoeff = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            diffgcoeff = varargin{ind+1};
        end
    else
        diffgcoeff = 1; %default is dependability
    end
    
    %check if gcoeff is provided
    ind = find(strcmp('gcoeff',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            gcoeff = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            gcoeff = varargin{ind+1};
        end
    else
        gcoeff = 1; %default is 1 (dependability)
    end
    
    %check if reltype is provided
    ind = find(strcmp('reltype',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            reltype = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            reltype = varargin{ind+1};
        end
    else
        reltype = 1; %default is 1 (equivalence)
    end
    
    %check if relvalue is provided
    ind = find(strcmp('relvalue',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            relcutoff = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            relcutoff = varargin{ind+1};
        end
    else
        relcutoff = .8; %default is .8
    end
    
    %check if plotrel is provided
    ind = find(strcmp('plotrel',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            plotrel = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            plotrel = varargin{ind+1};
        end
    else
        plotrel = 1; %default is 1
    end
    
    %check if plotrelline is provided
    ind = find(strcmp('plotrelline',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            plotrelline = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            plotrelline = varargin{ind+1};
        end
    else
        plotrelline = 2; %default is 2 (point estimate), matching psyrat_defaults
    end
    
    %check if relcentmeas is provided
    ind = find(strcmp('relcentmeas',varargin),1);
    if ~isempty(ind)
        if iscell(varargin{ind+1})
            relcentmeas = cell2mat(varargin{ind+1});
        elseif isnumeric(varargin{ind+1})
            relcentmeas = varargin{ind+1};
        end
    else
        relcentmeas = 1; %default is 1 (mean)
    end
    
elseif isempty(varargin)
    
    error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',...
        'Make sure each variable input is paired with a value \n',...
        'See help psyrat_relfigures for more information on inputs'));
    
end %if ~isempty(varargin)

if ~isempty(psyrat_prefs)
    switch analysis
        case {'sing','sing_sserr','sing_diff','sing_diff_sserr'}
            depcutoff = psyrat_prefs.view.depvalue;
            pdep = psyrat_prefs.view.plotdep;
            picc = psyrat_prefs.view.ploticc;
            showinct = psyrat_prefs.view.inctrltable;
            if isfield(psyrat_prefs.view,'tablessrel')
                showssrel = psyrat_prefs.view.tablessrel;
            else
                showssrel = 1;
            end
            showoverallt = psyrat_prefs.view.overalltable;
            plotntrials = psyrat_prefs.view.ntrials;
            showstddevt = psyrat_prefs.view.showstddevt;
            plotbetstddev = psyrat_prefs.view.showstddevf;
            plotwitstddev = 0;
            plotdepline = psyrat_prefs.view.plotdepline;
            meascutoff = psyrat_prefs.view.meascutoff;
            depcentmeas = psyrat_prefs.view.depcentmeas;
            if isfield(psyrat_prefs.view,'gcoeff')
                gcoeff = psyrat_prefs.view.gcoeff;
            else
                gcoeff = 1;
            end
            if isfield(psyrat_prefs.view,'diffgcoeff')
                diffgcoeff = psyrat_prefs.view.diffgcoeff;
            else
                diffgcoeff = 1;
            end
            
        case {'trt','trt_diff'}
            gcoeff = psyrat_prefs.view.gcoeff;
            reltype = psyrat_prefs.view.reltype;
            %test-retest occasion estimand: single-occasion (n'_o = 1, default)
            %vs multi-occasion composite (n'_o = k). Back-compat for older prefs.
            if isfield(psyrat_prefs.view,'noccmode') && ...
                    ~isempty(psyrat_prefs.view.noccmode)
                noccmode = psyrat_prefs.view.noccmode;
            else
                noccmode = 1;
            end
            if isfield(psyrat_prefs.view,'nocc')
                noccarg = psyrat_prefs.view.nocc;
            else
                noccarg = [];
            end
            %two-facet difference scores also need the difference-score
            %coefficient choice (dependability vs generalizability).
            if isfield(psyrat_prefs.view,'diffgcoeff')
                diffgcoeff = psyrat_prefs.view.diffgcoeff;
            else
                diffgcoeff = 1;
            end
            relcutoff = psyrat_prefs.view.relvalue;
            plotrel = psyrat_prefs.view.plotrel;
            plotrelline = psyrat_prefs.view.plotrelline;
            relcentmeas = psyrat_prefs.view.relcentmeas;
            picc = psyrat_prefs.view.ploticc;
            showinct = psyrat_prefs.view.inctrltable;
            showoverallt = psyrat_prefs.view.overalltable;
            plotntrials = psyrat_prefs.view.ntrials;
            showstddevt = psyrat_prefs.view.showstddevt;
            plotbetstddev = psyrat_prefs.view.showstddevf;
            plotwitstddev = 0;
            meascutoff = psyrat_prefs.view.meascutoff;
            if isfield(psyrat_prefs.view,'plotssrel')
                plotssrel = psyrat_prefs.view.plotssrel;
            else
                plotssrel = 1;
            end
            if isfield(psyrat_prefs.view,'plotssicc')
                plotssicc = psyrat_prefs.view.plotssicc;
            else
                plotssicc = 1;
            end
            if isfield(psyrat_prefs.view,'tablessrel')
                showssrel = psyrat_prefs.view.tablessrel;
            else
                showssrel = 1;
            end
            
    end
end


switch analysis
    case {'sing','sing_sserr','sing_diff','sing_diff_sserr'}
        %calculate reliabitliy information to be used for plotting and tables
        [psyrat_data, relerr] = psyrat_relsummary('psyrat_data',psyrat_data,...
            'analysis',analysis,...
            'depcutoff',depcutoff,...
            'meascutoff',meascutoff,...
            'depcentmeas',depcentmeas,...
            'gcoeff',gcoeff,...
            'diffgcoeff',diffgcoeff);
        
    case 'trt'
        %thread the test-retest occasion estimand through to relsummary so the
        %coefficient and its SEM share the same n'_o. Only pass an explicit k
        %(nocc) for the multi-occasion composite; otherwise relsummary defaults
        %k to the number of observed occasions.
        trtextra = {'noccmode',noccmode};
        if noccmode == 2 && ~isempty(noccarg)
            trtextra = [trtextra, {'nocc',noccarg}];
        end

        %calculate reliabitliy information to be used for plotting and tables.
        %
        %The subject-level test-retest design (trt_sserrvar, analysis 25) has its
        %own relsummary branch: it carries the location-scale person block
        %(gro_sds/pop_sdlog/ind_sdlog) and never writes sig_id/sig_err, which the
        %'trt' branch dereferences. Routing it to 'trt' therefore died with
        %"Unrecognized field name sig_id" (psyrat_relsummary.m:3203) -- RC-24.
        %
        %Only the relsummary CALL changes; `analysis` stays 'trt' so the
        %sserrvar-aware gates further down (the subject-level plots, the
        %psyrat_variancet table, and the duplicate-suppression check) keep
        %matching, exactly as they were written to.
        %
        %The two branches also take different argument NAMES for the same two
        %quantities: 'trt' takes relcutoff/relcentmeas, 'trt_sserr' takes
        %depcutoff/depcentmeas (psyrat_relsummary.m:261-271). The values are the
        %same GUI inputs.
        if isfield(psyrat_data.rel,'analysis') && ...
                strcmpi(psyrat_data.rel.analysis,'trt_sserrvar')
            [psyrat_data, relerr] = psyrat_relsummary('psyrat_data',psyrat_data,...
                'analysis','trt_sserr',...
                'gcoeff',gcoeff,...
                'reltype',reltype,...
                'depcutoff',relcutoff,...
                'meascutoff',meascutoff,...
                'depcentmeas',relcentmeas,...
                trtextra{:});
        else
            [psyrat_data, relerr] = psyrat_relsummary('psyrat_data',psyrat_data,...
                'analysis','trt',...
                'gcoeff',gcoeff,...
                'reltype',reltype,...
                'relcutoff',relcutoff,...
                'meascutoff',meascutoff,...
                'relcentmeas',relcentmeas,...
                trtextra{:});
        end

    case 'trt_diff'
        %thread the test-retest occasion estimand through to relsummary so the
        %coefficient and its SEM share the same n'_o. Only pass an explicit k
        %(nocc) for the multi-occasion composite.
        trtextra = {'noccmode',noccmode};
        if noccmode == 2 && ~isempty(noccarg)
            trtextra = [trtextra, {'nocc',noccarg}];
        end

        %two-facet difference-score reliability: trt occasion estimand +
        %difference-score coefficient choice (diffgcoeff).
        [psyrat_data, relerr] = psyrat_relsummary('psyrat_data',psyrat_data,...
            'analysis','trt_diff',...
            'gcoeff',gcoeff,...
            'diffgcoeff',diffgcoeff,...
            'reltype',reltype,...
            'relcutoff',relcutoff,...
            'meascutoff',meascutoff,...
            'relcentmeas',relcentmeas,...
            trtextra{:});

end

%if no good data were found then abort so you user can specify a different
%reliability threshold
if relerr.nogooddata == 1
    return;
end

%check which figures and tables have been specified to be shown
%prepare the figures and display them

%RC-24. The subject-level test-retest design (trt_sserrvar, analysis 25) uses
%relsummary's 'trt_sserr' branch, which populates only the per-participant table
%(ssrel_table + the good/bad id and trial-count bookkeeping). It does NOT produce
%the group-level quantities the rest of the test-retest family reads: event.rel,
%event.icc, event.betsd or event.trlcutoff. So the five group-level outputs below
%are unavailable for this design and are gated off here; each one was verified to
%ERROR (not to return a wrong number) without the gate. psyrat_startview_trt
%disables the matching controls so they are not offered.
%
%This list also named relsummary.reltype until RC-26, which stores it -- the
%branch always computed its coefficient FROM a reltype. The gate below keys off
%rel.analysis, never off the presence of that field, so storing it re-enables
%nothing here; only the justification above needed correcting.
%
%The three pre-existing inline tests further down (the subject-level plots, the
%psyrat_variancet table, and the duplicate-suppression check) express this same
%condition in the opposite direction -- they ENABLE subject-level output.
trt_ss = strcmp(analysis,'trt') && isfield(psyrat_data.rel,'analysis') && ...
    strcmpi(psyrat_data.rel.analysis,'trt_sserrvar');

%Plot that shows the relationship between the number of trials retained
%for averaging and dependability
switch analysis
    case {'sing','sing_diff'}
        if pdep == 1
            psyrat_depvtrialsplot('psyrat_data',psyrat_data,...
                'trials',[1 plotntrials],...
                'depline',plotdepline,...
                'gcoeff',gcoeff,...
                'depcutoff',depcutoff);
        end
    case {'sing_sserr','sing_diff_sserr'}
        if pdep == 1
            psyrat_ssrelplot('psyrat_data',psyrat_data,...
                'stat','dep',...
                'depline',plotdepline);
        end
    case 'trt'
        if plotrel == 1 && ~trt_ss
            psyrat_trt_relvtrialsplot('psyrat_data',psyrat_data,...
                'trials',[1 plotntrials],...
                'relline',plotrelline,...
                'relcutoff',relcutoff);
        end
    case 'trt_diff'
        if plotrel == 1
            %RC-31. A 'trt_diff' run reports TWO independently selected error
            %types, and each popup now drives the series it names.
            %
            %PER-CONDITION curve. psyrat_depvtrialsplot's trt_diff branch calls
            %psyrat_rel_trt -- the plain two-facet kernel, not a difference
            %kernel -- on the per-condition components, one subplot per event.
            %It is a per-condition coefficient, so it follows gcoeff, matching
            %the contract the toolbox states at psyrat_trt_coeflabel.m:94-99 and
            %implements in every table. It was fed diffgcoeff here, so the
            %difference-score selector silently moved a per-condition figure
            %while the per-condition rows in the table beside it stayed at
            %gcoeff. Divergence needed no deliberate action: prefs.view.gcoeff
            %is written on every trt view and persists, so a user whose last
            %plain-trt run was at Generalizability opened a trt_diff result with
            %the two popups already disagreeing.
            psyrat_depvtrialsplot('psyrat_data',psyrat_data,...
                'trials',[1 plotntrials],...
                'depline',plotrelline,...
                'gcoeff',gcoeff,...
                'depcutoff',relcutoff);

            %DIFFERENCE-SCORE curve, which is what the difference-score selector
            %was moving the figure above in place of. Genuinely a difference:
            %psyrat_diffrel_trt over the full cross-condition covariance blocks,
            %one line per group (a difference has no per-condition
            %stratification, so it cannot share the layout above). Drawn under
            %the same checkbox rather than a new one, which keeps the viewer
            %layout RC-25 just refit untouched.
            psyrat_trt_diffvtrialsplot('psyrat_data',psyrat_data,...
                'trials',[1 plotntrials],...
                'relline',plotrelline,...
                'diffgcoeff',diffgcoeff,...
                'relcutoff',relcutoff);
        end
end

%Plot that compares the intraclass correlation coefficients for each
%group and/or condition
if picc == 1
    if strcmp(analysis,'sing_sserr') || strcmp(analysis,'sing_diff_sserr')
        psyrat_ssrelplot('psyrat_data',psyrat_data,...
                'stat','icc',...
                'depline',plotdepline);
    elseif ~trt_ss
        psyrat_ptintervalplot('psyrat_data',psyrat_data,'stat','icc');
    end
end

%TRT subject-level plots
if strcmp(analysis,'trt') && isfield(psyrat_data.rel,'analysis') && ...
        strcmpi(psyrat_data.rel.analysis,'trt_sserrvar')
    if exist('plotssrel','var') && plotssrel == 1
        psyrat_ssrelplot('psyrat_data',psyrat_data,'stat','dep');
    end
    if exist('plotssicc','var') && plotssicc == 1
        psyrat_ssrelplot('psyrat_data',psyrat_data,'stat','icc');
    end
end

%plot that shows between-person standard deviation
if plotbetstddev == 1 && ~trt_ss
    psyrat_ptintervalplot('psyrat_data',psyrat_data,'stat','bet');
end

%plot that shows within-person standard deviation
if plotwitstddev == 1
    psyrat_ptintervalplot('psyrat_data',psyrat_data,'stat','wit');
end

%table displaying trial cutoff information
switch analysis
    case {'sing','sing_diff'}
        if showinct == 1
            psyrat_depcutofft('psyrat_data',psyrat_data,'gui',1);
        end
    case 'sing_sserr'
        if exist('showssrel','var') && showssrel == 1
            psyrat_variancet('psyrat_data',psyrat_data,'gui',1);
        end
    case 'sing_diff_sserr'
        if showinct == 1
            psyrat_depcutofft('psyrat_data',psyrat_data,'gui',1);
        end
        if exist('showssrel','var') && showssrel == 1
            psyrat_variancet('psyrat_data',psyrat_data,'gui',1);
        end
    case 'trt'
        if showinct == 1 && ~trt_ss
            psyrat_trt_relcutofft('psyrat_data',psyrat_data,'gui',1);
        end
        if isfield(psyrat_data.rel,'analysis') && ...
                strcmpi(psyrat_data.rel.analysis,'trt_sserrvar') && ...
                exist('showssrel','var') && showssrel == 1
            psyrat_variancet('psyrat_data',psyrat_data,'gui',1);
        end
    case 'trt_diff'
        if showinct == 1
            psyrat_depcutofft('psyrat_data',psyrat_data,'gui',1);
        end
end

%table displaying overall dependability information
switch analysis
    case {'sing','sing_sserr','sing_diff','sing_diff_sserr'}
        if showoverallt == 1
            psyrat_depoverallt('psyrat_data',psyrat_data,'gui',1);
        end
    case 'trt'
        if showoverallt == 1 && ~trt_ss
            psyrat_trt_reloverallt('psyrat_data',psyrat_data,'gui',1);
        end
    case 'trt_diff'
        if showoverallt == 1
            psyrat_depoverallt('psyrat_data',psyrat_data,'gui',1);
        end
end

%table including ICCs and between- and within-person standard deviations.
%psyrat_variancet renders the subject-level reliability table for
%subject-level (sserrvar) analyses and the sources-of-variance table
%otherwise. For sserrvar analyses that subject-level table is already shown
%above, owned by the "Table: Subject-Level Reliability" (showssrel) check, so
%skip the generic call here to avoid opening a byte-identical duplicate.
ssrelshown = any(strcmp(analysis,{'sing_sserr','sing_diff_sserr'})) || ...
    (strcmp(analysis,'trt') && isfield(psyrat_data.rel,'analysis') && ...
    strcmpi(psyrat_data.rel.analysis,'trt_sserrvar'));
if showstddevt == 1 && ~ssrelshown
    psyrat_variancet('psyrat_data',psyrat_data,'gui',1);
end

%show an error if the dependability threshold was never met for the cutoff.
%Titled per purpose but deliberately still NON-modal (2026-08-17 notification
%ruling); raised last so the later-drawn output cannot bury them at creation
%-- see the design comment in psyrat_criterionfigures, which mirrors these
%two branches.
if relerr.trlcutoff == 1
    errorstr = {};
    errorstr{end+1} = [L.Unit ' cutoffs for adequate dependability could not be calculated'];
    errorstr{end+1} = ['Data are too variable or there are not enough ' L.units];

    errordlg(errorstr,'Cutoff not calculable');
end

%show an error if none of the data had enough trials to meet the cutoff
if relerr.trlmax == 1
    errorstr = {};
    errorstr{end+1} = ['Not enough ' L.units ' are present in the current data'];
    errorstr{end+1} = 'Cutoffs represent an extrapolation beyond the data';

    errordlg(errorstr,'Extrapolation beyond data');
end

%Show the admissibility note when a quantity that must be a variance came out
%negative, so a GUI user who never reads the console is not left looking at a
%coefficient that is not interpretable as a reliability. Nothing was clipped to
%produce the numbers already drawn above; this explains them. Empty and silent
%for a clean run.
if isfield(psyrat_data,'relsummary') && ...
        isfield(psyrat_data.relsummary,'admissibility_note') && ...
        ~isempty(psyrat_data.relsummary.admissibility_note)

    warndlg(psyrat_data.relsummary.admissibility_note,...
        'Coefficient is not interpretable as a reliability');
end

end
