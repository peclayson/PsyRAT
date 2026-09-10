function vartable = psyrat_variancet(varargin)
%Table of ICCs, SEMs, and between- and within-person standard deviations
%
%psyrat_variancet('psyrat_data',psyrat_data,'gui',1);
%
%
%Inputs
% psyrat_data - PsyRAT Toolbox data structure array.
% gui - 0 for off, 1 for on
%
%Outputs
% vartable - table displaying information for ICCs, SEMs, and between-/
%  within-person standard deviations
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
            ['Inputs are incomplete. Provide name/value pairs, for example ' ...
            '''psyrat_variancet(''psyrat_data'',psyrat_data,''gui'',0)''.']);
    end

    %check if psyrat_data was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('psyrat_data',varargin),1);
    if ~isempty(ind)
        psyrat_data = varargin{ind+1};
    else
        error('varargin:psyrat_data',... %Error code and associated error
            ['psyrat_data was not provided. Pass the PsyRAT results struct, ' ...
            'for example ''psyrat_variancet(''psyrat_data'',psyrat_data,''gui'',0)''.']);
    end

    %check if gui was specified
    %If it is not found, set display error.
    ind = find(strcmpi('gui',varargin),1);
    if ~isempty(ind)
        gui = varargin{ind+1};
    else
        error('varargin:gui',... %Error code and associated error
            ['gui was not provided. Set ''gui'' to 0 (no figure) or 1 ' ...
            '(show figure) when calling psyrat_variancet.']);
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
icc = {};
betsd = {};
trlsd = {};
witsd = {};
sem = {};

if ~any(strcmp(psyrat_data.rel.analysis,{'ic_sserrvar','trt_sserrvar','ic_diff_sserrvar'}))
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

            %create a string with the icc point estimate and credible interval
            icc{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                psyrat_data.relsummary.group(gloc).event(eloc).icc.m,...
                psyrat_data.relsummary.group(gloc).event(eloc).icc.ll,...
                psyrat_data.relsummary.group(gloc).event(eloc).icc.ul);

            %create a string with the between-person standard devation point
            %estimate and credible interval
            betsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                psyrat_data.relsummary.group(gloc).event(eloc).betsd.m,...
                psyrat_data.relsummary.group(gloc).event(eloc).betsd.ll,...
                psyrat_data.relsummary.group(gloc).event(eloc).betsd.ul);

            %create a string with the trial main-effect SD (sigma_i); present
            %for the one-facet designs, '---' otherwise
            if isfield(psyrat_data.relsummary.group(gloc).event(eloc),'trlsd')
                trlsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                    psyrat_data.relsummary.group(gloc).event(eloc).trlsd.m,...
                    psyrat_data.relsummary.group(gloc).event(eloc).trlsd.ll,...
                    psyrat_data.relsummary.group(gloc).event(eloc).trlsd.ul);
            else
                trlsd{end+1} = '---';
            end

            %create a string with the within-person (person x trial) standard
            %devation point estimate and credible interval
            witsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                psyrat_data.relsummary.group(gloc).event(eloc).witsd.m,...
                psyrat_data.relsummary.group(gloc).event(eloc).witsd.ll,...
                psyrat_data.relsummary.group(gloc).event(eloc).witsd.ul);

            sem{end+1} = psyrat_sem_string(psyrat_data,gloc,eloc,...
                psyrat_data.relsummary.group(gloc).event(eloc).witsd.m,...
                psyrat_data.relsummary.group(gloc).event(eloc).witsd.ll,...
                psyrat_data.relsummary.group(gloc).event(eloc).witsd.ul);

        end

    end

elseif any(strcmp(psyrat_data.rel.analysis,{'ic_sserrvar','trt_sserrvar','ic_diff_sserrvar'}))

    %put data together to display in tables
    for gloc=1:ngroups
        for eloc=1:nevents

            ciedge = (1-psyrat_data.relsummary.ciperc)./2;  %follow the user's CI, not a fixed 95%
            d = psyrat_data.relsummary.data.g(gloc).e(eloc);
            ev = psyrat_data.relsummary.group(gloc).event(eloc);

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

            if isstruct(d) && isfield(d,'gro_sds') && isfield(d,'pop_sdlog')
                %gro_sds holds sigma_p as a STANDARD DEVIATION (declared as an
                %SD in the sserr Stan models and squared on entry by
                %psyrat_ssrel). The ICC below is a ratio of VARIANCES, so square
                %it here, exactly as the sd_id/b_sigma branch does at :241.
                bp_sd = cell2mat(d.gro_sds(:,1));
                bp_var = bp_sd .^ 2;
                pop_sdlog = d.pop_sdlog;

                %RC-06 group twin (owner ruling 2026-07-30). This branch serves
                %BOTH sserr designs, and 'trt_sserrvar' is two-facet: its
                %relsummary branch stores sig_occ, sig_trlxid, sig_occxid and
                %sig_trlxocc alongside reltype, so the one-facet partition below
                %dropped the entire occasion facet and never read reltype. The
                %printed group ICC therefore did not move at all across CE, CS
                %and CES, while the SAME run's caterpillar plot drew the correct
                %two-facet number (gro_icc, psyrat_ssrel_trt.m:132-144) as its
                %group reference line. One run, two different group ICCs.
                %
                %This was code-versus-code, not code-versus-paper: the plain
                %'trt' sibling already fills this column with
                %psyrat_rel_trt(obs=1,nocc=1) (psyrat_relsummary.m:3602-3615),
                %and RC-06's own fix pinned this analysis's SUBJECT-level ICC to
                %the same construction. The field mapping below is copied from
                %local_ssrel_trt_call (psyrat_relsummary.m:4778-4793), which is
                %what feeds gro_icc, so the two are the same quantity by
                %construction rather than by coincidence. psyrat_rel_trt
                %summarizes with mean + the same quantile edges used here
                %(psyrat_rel_trt.m:333-335), so only the PARTITION changes.
                %
                %'ic_sserrvar' is genuinely one-facet and keeps the expression
                %below unchanged. Any struct missing a crossed component, a
                %reltype or a gcoeff falls through to it as well, which is the
                %same fail-open convention the SEM helper uses.
                twofacet = strcmpi(psyrat_data.rel.analysis,'trt_sserrvar') && ...
                    all(isfield(d,{'sig_occ','sig_trl','sig_trlxid',...
                    'sig_occxid','sig_trlxocc'})) && ...
                    isfield(psyrat_data.relsummary,'gcoeff') && ...
                    isfield(psyrat_data.relsummary,'reltype');

                if twofacet
                    [icc_ll,icc_pt,icc_ul] = psyrat_rel_trt(...
                        'gcoeff',psyrat_data.relsummary.gcoeff,...
                        'reltype',psyrat_data.relsummary.reltype,...
                        'bp',bp_sd,...
                        'bo',d.sig_occ,...
                        'bt',d.sig_trl,...
                        'txp',d.sig_trlxid,...
                        'oxp',d.sig_occxid,...
                        'txo',d.sig_trlxocc,...
                        'err',exp(pop_sdlog),...
                        'obs',1,...
                        'nocc',1,...
                        'CI',psyrat_data.relsummary.ciperc);
                else
                    %With the trial main effect now estimated on the mean,
                    %exp(pop_sdlog)^2 is the RELATIVE residual sigma_pi,e^2. The
                    %absolute (dependability) ICC keeps the trial main effect
                    %sigma_i^2; the relative (generalizability) ICC excludes it
                    %(Rocha Table 2). sig_trl is present for the one-facet sserr
                    %designs; older structs without it fall back to the relative
                    %residual alone.
                    icc_err = exp(pop_sdlog).^2;
                    if isfield(psyrat_data.relsummary,'gcoeff') && ...
                            psyrat_data.relsummary.gcoeff == 1 && ...
                            isfield(d,'sig_trl') && ~isempty(d.sig_trl)
                        icc_err = icc_err + d.sig_trl(:).^2;
                    end
                    icc_pt = mean(bp_var ./ (bp_var + icc_err));
                    icc_ll = quantile(bp_var ./ (bp_var + icc_err),ciedge);
                    icc_ul = quantile(bp_var ./ (bp_var + icc_err),1-ciedge);
                end

                %create a string with the icc point estimate and credible interval
                icc{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                    icc_pt, icc_ll, icc_ul);

                %create a string with the between-person standard devation point
                %estimate and credible interval. Summarize the SD draws
                %directly (posterior mean of sigma_p), matching the sibling
                %branch at :267-270 and psyrat_relsummary; taking sqrt(mean(bp_var))
                %would instead report the posterior root-mean-square.
                betsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                    mean(bp_sd),...
                    quantile(bp_sd,ciedge),...
                    quantile(bp_sd,1-ciedge));

                %create a string with the within-person standard devation point
                %estimate and credible interval
                wit_m = exp(mean(pop_sdlog));
                wit_ll = exp(quantile(pop_sdlog,ciedge));
                wit_ul = exp(quantile(pop_sdlog,1-ciedge));
                witsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                    wit_m,wit_ll,wit_ul);
            elseif isstruct(d) && isfield(d,'sd_id') && isfield(d,'b_sigma')
                bp_sd = cell2mat(d.sd_id.raw);
                pop_sdlog = cell2mat(d.b_sigma.raw);
                bp_var = bp_sd.^2;

                %RC-20. This is the 'ic_diff_sserrvar' group row: its branch of
                %psyrat_relsummary stores sd_id/b_sigma/trl_varcov rather than
                %the gro_sds/pop_sdlog pair the branch above expects, so the
                %trial main effect is reached through trl_varcov (stored as a
                %VARIANCE, hence the sqrt inside the resolver). It was
                %previously not reached at all here -- this ICC had no gcoeff
                %branch, so flipping the coefficient selector changed only the
                %column header, not the number under it. Mirrors :209-213:
                %sigma_i^2 joins the error for the absolute (dependability) ICC
                %and is excluded for the relative (generalizability) one (Rocha
                %Table 2). Structs carrying none of the trial fields resolve to
                %sigma_i = 0, which is byte-identical to the prior behavior.
                icc_err = exp(pop_sdlog).^2;
                if isfield(psyrat_data.relsummary,'gcoeff') && ...
                        psyrat_data.relsummary.gcoeff == 1
                    sigi = psyrat_sem_sigtrl(d,pop_sdlog);
                    icc_err = icc_err + sigi.^2;
                end

                icc{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                    mean(bp_var ./ (bp_var + icc_err)),...
                    quantile(bp_var ./ (bp_var + icc_err),ciedge),...
                    quantile(bp_var ./ (bp_var + icc_err),1-ciedge));

                betsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                    mean(bp_sd),...
                    quantile(bp_sd,ciedge),...
                    quantile(bp_sd,1-ciedge));

                wit_m = mean(exp(pop_sdlog));
                wit_ll = quantile(exp(pop_sdlog),ciedge);
                wit_ul = quantile(exp(pop_sdlog),1-ciedge);
                witsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                    wit_m,wit_ll,wit_ul);
            else
                %fallback for structures where summary statistics were already computed
                icc{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                    ev.icc.m, ev.icc.ll, ev.icc.ul);
                betsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                    ev.betsd.m, ev.betsd.ll, ev.betsd.ul);
                wit_m = ev.witsd.m;
                wit_ll = ev.witsd.ll;
                wit_ul = ev.witsd.ul;
                witsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                    wit_m,wit_ll,wit_ul);
            end

            sem{end+1} = psyrat_sem_string(psyrat_data,gloc,eloc,...
                wit_m,wit_ll,wit_ul);

        end

        if strcmp(psyrat_data.rel.analysis,'ic_diff_sserrvar') && ...
                isfield(psyrat_data.relsummary.group(gloc),'diffscore')
            diffscore = psyrat_data.relsummary.group(gloc).diffscore;

            switch analysis
                case 3
                    label{end+1} = 'diff score';
                case 4
                    label{end+1} = [gnames{gloc} ' - diff score'];
            end

            idvar = cell2mat(psyrat_data.rel.out.id_varcov{:,gloc});
            errvar = cell2mat(psyrat_data.rel.out.err_varcov{:,gloc});
            diff_bp = idvar(:,1,1) + idvar(:,2,2) - (2 .* idvar(:,2,1));
            diff_wp = errvar(:,1,1) + errvar(:,2,2) - (2 .* errvar(:,2,1));
            diff_bp = max(diff_bp,0);
            diff_wp = max(diff_wp,0);
            ciedge = (1-psyrat_data.relsummary.ciperc)./2;  %follow the user's CI, not a fixed 95%

            icc{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                diffscore.icc_pt,diffscore.icc_ll,diffscore.icc_ul);
            betsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                mean(sqrt(diff_bp)),...
                quantile(sqrt(diff_bp),ciedge),...
                quantile(sqrt(diff_bp),1-ciedge));
            witsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                mean(sqrt(diff_wp)),...
                quantile(sqrt(diff_wp),ciedge),...
                quantile(sqrt(diff_wp),1-ciedge));

            dtab = diffscore.ssrel_table;
            sem{end+1} = sprintf(' %0.2f SD: %0.2f',...
                mean(dtab.sem_pt),std(dtab.sem_pt));
        end
    end


end

%analysis-unit label fragments (trial vs split). For a parallel-splits run the
%measurement facet is the split, so the trial-main-effect column is relabeled to
%"split"; psyrat_unitlabels defaults to "trial" when no splits flag is present,
%so non-splits output is byte-identical.
usplit = 1;
if isfield(psyrat_data.rel,'splits'); usplit = psyrat_data.rel.splits; end
L = psyrat_unitlabels(usplit);

%ICC and SEM are estimand-dependent: their VALUE changes with gcoeff (the ICC
%branch above adds sigma_i^2 for the absolute coefficient) while the header
%historically did not, so an archived export could not be reconstructed. Suffix
%both with the decision type, reusing the _Dep/_Gen tokens psyrat_dod_observedt
%already emits. Fail-open: when the decision type cannot be resolved (legacy
%structs predating gcoeff_name) the unsuffixed names are kept, so old output is
%byte-identical.
esuf = psyrat_estimand_suffix(psyrat_data);

%create table for displaying between-person, trial, and within-person standard
%deviation information. The one-facet ('ic') analysis adds a trial main-effect
%(sigma_i) column; the other analyses keep the original five-column layout.
if strcmp(psyrat_data.rel.analysis,'ic')
    vartable = table(label',betsd',trlsd',witsd',sem',icc');
    vartable.Properties.VariableNames = {'Label',...
        'Between_StdDev',[L.Unit '_StdDev'],'Within_StdDev',...
        ['SEM' esuf],['ICC' esuf]};
else
    vartable = table(label',betsd',witsd',sem',icc');
    if strcmp(psyrat_data.rel.analysis,'trt')
        vartable.Properties.VariableNames = {'Label',...
            'Between_StdDev','Within_StdDev',['SEM' esuf],['ICC' esuf]};
    elseif any(strcmp(psyrat_data.rel.analysis,{'ic_diff','trt_diff'}))
        vartable.Properties.VariableNames = {'Label',...
            'Between_StdDev','Within_StdDev',['SEM' esuf],['ICC' esuf]};
    elseif any(strcmp(psyrat_data.rel.analysis,{'ic_sserrvar','trt_sserrvar','ic_diff_sserrvar'}))
        vartable.Properties.VariableNames = {'Label',...
            'Group_Between_StdDev','Group_Within_StdDev',...
            ['Group_SEM' esuf],['Group_ICC' esuf]};
    end
end

%only render the standard-deviation figure when a GUI is requested. A headless
%call (gui,0) returns the vartable assembled above without opening a figure,
%mirroring the gui gate in psyrat_depcutofft/psyrat_depoverallt so psyrat_report
%can build the sserr group-level table render-free (C3 step 1). vartable is
%unchanged, so gui==1 output is byte-identical.
if gui == 1

    %define parameters for figure size
    figwidth = 800;
    figheight = 400;

    %define space between rows and first row location
    rowspace = 25;
    row = figheight - rowspace*2;
    %the stated width follows the user's CI setting; every quantity in this
    %table (SDs, SEMs, ICCs) now shares that coverage
    name = ['Point and ' sprintf('%.0f',psyrat_data.relsummary.ciperc*100) ...
        '% Interval Estimates for the Between- '...
        'and Within-Person Standard Deviations, SEMs, and ICCs'];

    %for test-retest the SEM is tied to the selected occasion estimand, so label
    %the table distinctly: single-occasion (n'_o = 1) vs multi-occasion composite
    %(n'_o = k) (audit sec.14). The SEM here uses the SAME n'_o as the reported
    %coefficient.
    if any(strcmp(psyrat_data.rel.analysis,{'trt','trt_diff'}))
        name = [name ' - ' psyrat_trt_occlabel(psyrat_data)];
    elseif strcmp(psyrat_data.rel.analysis,'trt_sserrvar')
        %RC-06 group twin. Until this row's SEM became two-facet it contained no
        %n'_o term at all, so the table was invariant to the occasion selector
        %and needed no label. It is not invariant now: on the standard fixture
        %the SEM moves 0.88 -> 0.66 between n'_o = 1 and n'_o = 2 while the ICC
        %beside it does not move at all, because the ICC is pinned at
        %obs = 1, n'_o = 1 by the same ruling.
        %
        %That is why this cannot reuse the {'trt','trt_diff'} clause above.
        %psyrat_trt_occlabel names the estimand of "the reported coefficient",
        %which would attribute n'_o = k to a row whose coefficient is at 1. Name
        %the two separately instead. (Plain 'trt' carries the same imprecision --
        %its variance-table ICC is also psyrat_rel_trt(obs=1,nocc=1) while its
        %SEM uses relsummary.nocc -- but changing its title would move output
        %this ruling did not scope, so it is left alone.)
        if isfield(psyrat_data.relsummary,'nocc') && ...
                ~isempty(psyrat_data.relsummary.nocc)
            semnocc = psyrat_data.relsummary.nocc;
        else
            semnocc = 1;
        end
        name = [name sprintf(' - SEM at n''_o = %d, ICC at n''_o = 1', semnocc)];
    end

    %create gui for standard-deviation table
    var_gui= psyrat_newfigure('unit','pix',...
        'position',[1250 600 figwidth figheight],...
        'menub','no',...
        'Tag','psyrat_output',...
        'name',name,...
        'numbertitle','off',...
        'resize','off');

    if ~any(strcmp(psyrat_data.rel.analysis,{'ic_sserrvar','trt_sserrvar','ic_diff_sserrvar'}))
        %Print the name of the loaded dataset
        uicontrol(var_gui,'Style','text','fontsize',16,...
            'HorizontalAlignment','center',...
            'String',...
            'Group-Level Between- and Within-Person Standard Deviations, SEMs, and ICCs',...
            'Position',[0 row figwidth 25]);
    end

    %for test-retest, the reported SEM is the single-occasion (n'_o = 1) or
    %multi-occasion composite (n'_o = k) SEM that matches the reported
    %coefficient; label it so the two estimands are not confused (audit sec.14)
    if strcmp(psyrat_data.rel.analysis,'trt')
        uicontrol(var_gui,'Style','text','fontsize',11,...
            'HorizontalAlignment','center',...
            'String',psyrat_trt_occlabel(psyrat_data),...
            'Position',[0 row+22 figwidth 20]);
    end

    %Start a table
    t = uitable('Parent',var_gui,'Position',...
        [25 100 figwidth-50 figheight-175],...
        'Data',table2cell(vartable));
    %drive the on-screen headers from the table's own variable names so the
    %GUI, the CSV and the XLSX cannot drift apart (they were three independent
    %hand-maintained lists, and the sserr layout had already diverged).
    set(t,'ColumnName',psyrat_variancet_displaynames( ...
        vartable.Properties.VariableNames));
    if strcmp(psyrat_data.rel.analysis,'ic')
        set(t,'ColumnWidth',{200 130 130 130 130 130});
    else
        set(t,'ColumnWidth',{200 130 130 130 130});
    end
    set(t,'RowName',[]);

    if ~any(strcmp(psyrat_data.rel.analysis,{'ic_sserrvar','trt_sserrvar','ic_diff_sserrvar'}))

        %Create a save button that will save the table
        uicontrol(var_gui,'Style','push','fontsize',14,...
            'HorizontalAlignment','center',...
            'String','Save Table',...
            'Position', [figwidth/8 25 figwidth/4 50],...
            'Callback',{@psyrat_savevartable,psyrat_data,vartable});

    else

        %Create a save button that will save the table
        uicontrol(var_gui,'Style','push','fontsize',14,...
            'HorizontalAlignment','center',...
            'String',{'Save Group-Level';'Estimates'},...
            'Position', [figwidth/8 25 figwidth/4 50],...
            'Callback',{@psyrat_savevartable,psyrat_data,vartable});

    end

    if strcmp(psyrat_data.rel.analysis,'trt')

        %Create a save button that will view variance components
        uicontrol(var_gui,'Style','push','fontsize',14,...
            'HorizontalAlignment','center',...
            'String',{'View All';'Standard Deviations'},...
            'Position', [(figwidth/8)*5 25 figwidth/4 50],...
            'Callback',{@psyrat_viewvarcomp,psyrat_data});

    elseif any(strcmp(psyrat_data.rel.analysis,{'ic_sserrvar','trt_sserrvar','ic_diff_sserrvar'}))

        %Create a save button that will view variance components
        uicontrol(var_gui,'Style','push','fontsize',14,...
            'HorizontalAlignment','center',...
            'String',{'Save Subject';'Coefficients'},...
            'Position', [(figwidth/8)*5 25 figwidth/4 50],...
            'Callback',{@psyrat_savesscoeffs,psyrat_data});

    elseif strcmp(psyrat_data.rel.analysis,'ic_diff')

        %Create a save button that will view variance components
        uicontrol(var_gui,'Style','push','fontsize',14,...
            'HorizontalAlignment','center',...
            'String',{'View All';'Standard Deviations'},...
            'Position', [(figwidth/8)*5 25 figwidth/4 50],...
            'Callback',{@psyrat_viewdiffcomp,psyrat_data});
    end

end


end



function psyrat_savevartable(varargin)
%if the user pressed the button to save the table with the information
%about between- and within-person standard deviations

%parse inputs
psyrat_data = varargin{3};
vartable = varargin{4};

%ask the user where the file should be saved. XLSX now writes cross-platform
%via writecell/writetable, so both formats are offered on every platform
%(including macOS).
[savename, savepath] = uiputfile(...
    {'*.xlsx','Excel File (.xlsx)';'*.csv',...
    'Comma-Separated Value File (.csv)'},...
    'Where would you like to save table?');

[~,~,ext] = fileparts(fullfile(savepath,savename));

%for test-retest, the SEM is tied to the selected occasion estimand, so
%record the single- vs multi-occasion label in the exported header
%(audit sec.14)
if strcmp(psyrat_data.rel.analysis,'trt')
    occlabel = psyrat_trt_occlabel(psyrat_data);
else
    occlabel = '';
end

%save as either excel or csv file
if strcmp(ext,'.xlsx')

    filehead = {'Table Generated on'; datestr(clock);''};
    filehead{end+1} = sprintf('PsyRAT Toolbox v%s',psyrat_data.ver);
    if ~isempty(occlabel)
        filehead{end+1} = occlabel;
    end
    filehead{end+1} = '';
    filehead{end+1} = sprintf('Dataset: %s',psyrat_data.rel.filename);
    filehead{end+1} = psyrat_iterline(psyrat_data.rel);
    %run provenance (seed/engine/priors/convergence); fail-open
    filehead = [filehead; psyrat_provenance_lines(psyrat_data)];
    filehead{end+1}='';
    filehead{end+1}='';

    %mirror the CSV display-name conversion so XLSX headers match the GUI/CSV
    %("Between_StdDev" -> "Between Std Dev"); see CSV branch below (F-CM-M4xlsxhdr)
    xlsxtable = vartable;
    xlsxtable.Properties.VariableNames = psyrat_variancet_displaynames( ...
        xlsxtable.Properties.VariableNames);

    writecell(filehead,fullfile(savepath,savename));
    writetable(xlsxtable,fullfile(savepath,savename),...
        'Range',strcat('A',num2str(length(filehead))));

elseif strcmp(ext,'.csv')

    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s\n','Table Generated on');
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'PsyRAT Toolbox v%s\n',psyrat_data.ver);
    if ~isempty(occlabel)
        fprintf(fid,'%s\n',psyrat_quote_csv_header(occlabel));
    end
    fprintf(fid,' \n');
    fprintf(fid,'%s\n',psyrat_quote_csv_header(...
        sprintf('Dataset: %s',psyrat_data.rel.filename)));
    fprintf(fid,'%s\n',psyrat_quote_csv_header(...
        psyrat_iterline(psyrat_data.rel)));
    psyrat_fprintf_provenance(fid, psyrat_data);
    fprintf(fid,' \n');

    %build the header and rows from the table itself so the optional trial
    %main-effect (sigma_i) column is exported when present (one-facet designs).
    %Match the on-screen uitable labels: the variable names compress to
    %"...StdDev", whereas the GUI shows "... Std Dev" (F-M4hdr).
    vnames = psyrat_variancet_displaynames(vartable.Properties.VariableNames);
    fprintf(fid,'%s',strjoin(vnames,','));
    fprintf(fid,' \n');

    ncol = width(vartable);
    for i = 1:height(vartable)
        rowvals = cell(1,ncol);
        for c = 1:ncol
            rowvals{c} = char(vartable{i,c});
        end
        fprintf(fid,'%s\n',strjoin(rowvals,','));
    end

    fclose(fid);

end


end



function psyrat_savesscoeffs(varargin)
%if the user pressed the button to save the table with the information
%about between- and within-person standard deviations

%parse inputs
psyrat_data = varargin{3};

[savename, savepath] = uiputfile(...
    {'*.csv',...
    'Comma-Separated Value File (.csv)'},...
    'Where would you like to save table?');

%bail out if the user cancelled the save dialog
if isequal(savename,0) || isequal(savepath,0)
    return;
end

%this helper writes a bespoke multi-section CSV (a per-subject loop, not a plain
%table), so it is CSV-only; normalize the extension so a typed/forced .xlsx name
%cannot produce invalid workbook bytes (F-CM-M4ssxlsx).
[~,~,ext] = fileparts(savename);
if ~strcmpi(ext,'.csv')
    savename = [savename '.csv'];
end

%save as csv file
fid = fopen(fullfile(savepath,savename),'w');
fprintf(fid,'%s\n','Table Generated on');
fprintf(fid,'%s\n',datestr(clock));
fprintf(fid,'PsyRAT Toolbox v%s\n',psyrat_data.ver);
fprintf(fid,' \n');
fprintf(fid,'%s\n',psyrat_quote_csv_header(...
    sprintf('Dataset: %s',psyrat_data.rel.filename)));
fprintf(fid,'%s\n',psyrat_quote_csv_header(...
    psyrat_iterline(psyrat_data.rel)));
psyrat_fprintf_provenance(fid, psyrat_data);
fprintf(fid,'Reliability Cutoff: %0.4f\n',psyrat_data.relsummary.depcutoff);
fprintf(fid,'%s\n','Subject-Specific Estimates');
fprintf(fid,' \n');
fprintf(fid,' \n');

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


if ngroups == 1 && nevents == 1
    analysis = 1;
elseif ngroups > 1 && nevents == 1
    analysis = 2;
elseif ngroups == 1 && nevents > 1
    analysis = 3;
elseif ngroups > 1 && nevents > 1
    analysis = 4;
end

for gloc=1:ngroups
    for eloc=1:nevents
        %label for group and/or event
        switch analysis
            case 1
                label = 'Measurement';
            case 2
                label = gnames{gloc};
            case 3
                label = enames{eloc};
            case 4
                label = [gnames{gloc} ' - ' enames{eloc}];
        end

        t_out = psyrat_data.relsummary.group(gloc).event(eloc).ssrel_table;

        %quote per the CSV convention (B19): the section label embeds
        %runtime group/event names, which can carry commas
        fprintf(fid,'%s\n',psyrat_quote_csv_header(label));

        %the dep_* columns hold the selected coefficient (absolute phi_s or
        %relative G_s); label the export header to match the decision type.
        if isfield(psyrat_data.relsummary,'gcoeff_name') && ...
                strcmpi(psyrat_data.relsummary.gcoeff_name,'gen')
            coeff_label = 'Generalizability';
        else
            coeff_label = 'Dependability';
        end

        fprintf(fid,'%s\n', strcat('ID,Participant Meets Reliability Cutoff',...
            ',Trial Counts,',coeff_label,',',coeff_label,' Lower Limit',...
            ',',coeff_label,' Upper Limit,ICC,ICC Lower Limit,ICC Upper Limit',...
            ',SEM,SEM Lower Limit,SEM Upper Limit,Error Variance'));

        id = t_out.id;
        include = t_out.ind2include;
        trls = t_out.trls;
        dep_pt = t_out.dep_pt;
        dep_ll = t_out.dep_ll;
        dep_ul = t_out.dep_ul;
        icc_pt = t_out.icc_pt;
        icc_ll = t_out.icc_ll;
        icc_ul = t_out.icc_ul;
        sem_pt = t_out.sem_pt;
        sem_ll = t_out.sem_ll;
        sem_ul = t_out.sem_ul;
        ss_errvar = t_out.ss_errvar;


        for ii = 1:height(t_out)
            formatspec = '%s,%d,%d,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n';
            fprintf(fid,formatspec,...
                char(id{ii}),...
                include(ii),...
                trls(ii),...
                dep_pt(ii),...
                dep_ll(ii),...
                dep_ul(ii),...
                icc_pt(ii),...
                icc_ll(ii),...
                icc_ul(ii),...
                sem_pt(ii),...
                sem_ll(ii),...
                sem_ul(ii),...
                ss_errvar(ii));
        end

    end
end

fclose(fid);




end


function psyrat_viewvarcomp(varargin)
%present another gui to view the variance components for trt data

%parse inputs
psyrat_data = varargin{3};

%analysis-unit label fragments (trial vs split) for the crossed-design
%variance-component labels; defaults to "trial" so non-splits output is
%byte-identical.
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


%create placeholders for displaying data in tables in guis
label = {};
bp = {};
bo = {};
bt = {};
txp = {};
oxp = {};
txo = {};
err = {};

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

        %create a string with the icc point estimate and credible interval
        bp{end+1} = sprintf('%0.2f',...
            mean(psyrat_data.relsummary.data.g(gloc).e(eloc).sig_id.raw));

        bo{end+1} = sprintf('%0.2f',...
            mean(psyrat_data.relsummary.data.g(gloc).e(eloc).sig_occ.raw));

        bt{end+1} = sprintf('%0.2f',...
            mean(psyrat_data.relsummary.data.g(gloc).e(eloc).sig_trl.raw));

        txp{end+1} = sprintf('%0.2f',...
            mean(psyrat_data.relsummary.data.g(gloc).e(eloc).sig_trlxid.raw));

        oxp{end+1} = sprintf('%0.2f',...
            mean(psyrat_data.relsummary.data.g(gloc).e(eloc).sig_occxid.raw));

        txo{end+1} = sprintf('%0.2f',...
            mean(psyrat_data.relsummary.data.g(gloc).e(eloc).sig_trlxocc.raw));

        err{end+1} = sprintf('%0.2f',...
            mean(psyrat_data.relsummary.data.g(gloc).e(eloc).sig_err.raw));

    end
end


%create table for displaying standard deviations
vartable = table(label',bp',bo',bt',oxp',txp',txo',err');

vartable.Properties.VariableNames = {'Label',...
    'bp','bo','bt','oxp','txp','txo','err'};

%define parameters for figure size
figwidth = 1000;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;
name = 'Point Estimates of All Standard Deviation Components';

%create gui for standard-deviation table
varall_gui= psyrat_newfigure('unit','pix',...
    'position',[1250 600 figwidth figheight],...
    'menub','no',...
    'Tag','psyrat_output',...
    'name',name,...
    'numbertitle','off',...
    'resize','off');

%Print the name of the loaded dataset
uicontrol(varall_gui,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',...
    'Estimated Standard Deviation Components',...
    'Position',[0 row figwidth 25]);

%Start a table
t = uitable('Parent',varall_gui,'Position',...
    [25 100 figwidth-50 figheight-175],...
    'Data',table2cell(vartable));
set(t,'ColumnName',{'Label'...
    'Between Person' 'Between Session' ['Between ' L.Unit] ...
    'Person x Session' ['Person x ' L.Unit] ['Session x ' L.Unit]...
    'Within Person (Error)'});
set(t,'ColumnWidth',{200 130 130 130 130 130 130 130});
set(t,'RowName',[]);

%Create a save button that will save the table
uicontrol(varall_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Save Table',...
    'Position', [figwidth/8 25 figwidth/4 50],...
    'Callback',{@psyrat_saveallvartable,psyrat_data,vartable});

end

function psyrat_saveallvartable(varargin)
%if the user pressed the button to save all of the standard deviation
%components

%parse inputs
psyrat_data = varargin{3};
vartable = varargin{4};

%analysis-unit label fragments (trial vs split) for the crossed-design
%variance-component headers; defaults to "trial" so non-splits output is
%byte-identical.
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

%for test-retest, the SEM is tied to the selected occasion estimand, so
%record the single- vs multi-occasion label in the exported header
%(audit sec.14)
if strcmp(psyrat_data.rel.analysis,'trt')
    occlabel = psyrat_trt_occlabel(psyrat_data);
else
    occlabel = '';
end

%save as either excel or csv file
if strcmp(ext,'.xlsx')

    filehead = {'Table Generated on'; datestr(clock);''};
    filehead{end+1} = sprintf('PsyRAT Toolbox v%s',psyrat_data.ver);
    if ~isempty(occlabel)
        filehead{end+1} = occlabel;
    end
    filehead{end+1} = '';
    filehead{end+1} = sprintf('Dataset: %s',psyrat_data.rel.filename);
    filehead{end+1} = psyrat_iterline(psyrat_data.rel);
    %run provenance (seed/engine/priors/convergence); fail-open
    filehead = [filehead; psyrat_provenance_lines(psyrat_data)];
    filehead{end+1}='';
    filehead{end+1}='';

    %the table stores compact component codes (bp/bo/bt/...) as variable names;
    %map them to the on-screen ColumnName labels so the XLSX header matches the
    %GUI/CSV instead of writing the raw codes (F-CM-M4xlsxhdr)
    xlsxtable = vartable;
    xlsxtable.Properties.VariableNames = {'Label',...
        'Between Person','Between Session',['Between ' L.Unit],...
        'Person x Session',['Person x ' L.Unit],['Session x ' L.Unit],...
        'Within Person (Error)'};

    writecell(filehead,fullfile(savepath,savename));
    writetable(xlsxtable,fullfile(savepath,savename),...
        'Range',strcat('A',num2str(length(filehead))));

elseif strcmp(ext,'.csv')

    fid = fopen(fullfile(savepath,savename),'w');
    fprintf(fid,'%s\n','Table Generated on');
    fprintf(fid,'%s\n',datestr(clock));
    fprintf(fid,'PsyRAT Toolbox v%s\n',psyrat_data.ver);
    if ~isempty(occlabel)
        fprintf(fid,'%s\n',psyrat_quote_csv_header(occlabel));
    end
    fprintf(fid,' \n');
    fprintf(fid,'%s\n',psyrat_quote_csv_header(...
        sprintf('Dataset: %s',psyrat_data.rel.filename)));
    fprintf(fid,'%s\n',psyrat_quote_csv_header(...
        psyrat_iterline(psyrat_data.rel)));
    psyrat_fprintf_provenance(fid, psyrat_data);
    fprintf(fid,' \n');

    fprintf(fid,'%s', ['Label,Between Person, Between Session,'...
        'Between ' L.Unit ',Person x Session,Person x ' L.Unit...
        ',Session x ' L.Unit ',Within Person (Error)']);
    fprintf(fid,' \n');

    for i = 1:height(vartable)
        formatspec = '%s,%s,%s,%s,%s,%s,%s,%s\n';
        fprintf(fid,formatspec,char(vartable{i,1}),...
            char(vartable{i,2}), char(vartable{i,3}),...
            char(vartable{i,4}), char(vartable{i,5}),...
            char(vartable{i,6}), char(vartable{i,7}),...
            char(vartable{i,8}));
    end

    fclose(fid);

end


end


function psyrat_viewdiffcomp(varargin)
%present another gui to view the variance components for diff score data

%parse inputs
psyrat_data = varargin{3};

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
icc = {};
betsd = {};
witsd = {};
bettrlsd = {};
bp_cov = {};
bt_cov = {};
wp_cov = {};
icc = {};


%put data together to display in tables
for gloc=1:ngroups
    for eloc=1:nevents

        %label for group and/or event
        switch analysis
            case 3
                label{end+1} = enames{eloc};
            case 4
                label{end+1} = [gnames{gloc} ' - ' enames{eloc}];
        end

        %create a string with the icc point estimate and credible interval
        icc{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            psyrat_data.relsummary.group(gloc).event(eloc).icc.m,...
            psyrat_data.relsummary.group(gloc).event(eloc).icc.ll,...
            psyrat_data.relsummary.group(gloc).event(eloc).icc.ul);

        %create a string with the between-person standard devation point
        %estimate and credible interval
        betsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            psyrat_data.relsummary.group(gloc).event(eloc).betsd.m,...
            psyrat_data.relsummary.group(gloc).event(eloc).betsd.ll,...
            psyrat_data.relsummary.group(gloc).event(eloc).betsd.ul);

        %create a string with the within-person standard devation point
        %estimate and credible interval
        witsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
            psyrat_data.relsummary.group(gloc).event(eloc).witsd.m,...
            psyrat_data.relsummary.group(gloc).event(eloc).witsd.ll,...
            psyrat_data.relsummary.group(gloc).event(eloc).witsd.ul);

        if eloc == 1
            bettrlsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                sqrt(psyrat_data.relsummary.group(gloc).diffscore.bt1_var_pt),...
                sqrt(psyrat_data.relsummary.group(gloc).diffscore.bt1_var_ll),...
                sqrt(psyrat_data.relsummary.group(gloc).diffscore.bt1_var_ul));
        elseif eloc == 2
            bettrlsd{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
                sqrt(psyrat_data.relsummary.group(gloc).diffscore.bt2_var_pt),...
                sqrt(psyrat_data.relsummary.group(gloc).diffscore.bt2_var_ll),...
                sqrt(psyrat_data.relsummary.group(gloc).diffscore.bt2_var_ul));
        end

        bp_cov{end+1} = '---';
        bt_cov{end+1} = '---';
        wp_cov{end+1} = '---';
    end

    switch analysis
        case 3
            label{end+1} = 'diff score';
        case 4
            label{end+1} = [gnames{gloc} ' - diff score'];
    end

    betsd{end+1} = '---';
    witsd{end+1} = '---';
    bettrlsd{end+1} = '---';

    %create a string with the dependability point estimate and credible
    %interval for cutoff data
    icc{end+1} = sprintf(' %0.2f CI [%0.2f %0.2f]',...
        psyrat_data.relsummary.group(gloc).diffscore.icc_pt,...
        psyrat_data.relsummary.group(gloc).diffscore.icc_ll,...
        psyrat_data.relsummary.group(gloc).diffscore.icc_ul);

    bp_cov{end+1} = sprintf('%0.2f CI [%0.2f %0.2f]',...
        sqrt(psyrat_data.relsummary.group(gloc).diffscore.bp_cov_pt),...
        sqrt(psyrat_data.relsummary.group(gloc).diffscore.bp_cov_ll),...
        sqrt(psyrat_data.relsummary.group(gloc).diffscore.bp_cov_ul));

    bt_cov{end+1} = sprintf('%0.2f CI [%0.2f %0.2f]',...
        sqrt(psyrat_data.relsummary.group(gloc).diffscore.bt_cov_pt),...
        sqrt(psyrat_data.relsummary.group(gloc).diffscore.bt_cov_ll),...
        sqrt(psyrat_data.relsummary.group(gloc).diffscore.bt_cov_ul));

    wp_cov{end+1} = sprintf('%0.2f',...
        psyrat_data.relsummary.group(gloc).diffscore.wp_cov);

end

%create table for displaying standard deviations
vartable = table(label',betsd',bp_cov',bettrlsd',bt_cov',witsd',wp_cov',icc');

%Every ICC in this difference-score table is a single-observation quantity: the
%per-event rows come from relsummary's event ICC (one trial) and the diff-score
%row from psyrat_diffrel's n'=1 denominators (denom_rel_err/denom_abs_err,
%psyrat_diffrel.m:251-252, explained at its :235-239), which for the
%difference score is a DERIVED quantity, not a Rocha Table 6 entry. The header
%says so; the values are untouched. Keep in step with the CSV header emitted by
%psyrat_savealldifftable below.
vartnames = {'Label',...
    'Between-Person Std Dev','BP Covariance Std Dev',...
    'Between-Trial Std Dev','BT Covariance Std Dev',...
    'Error Std Dev','Error Covariance Std Dev',...
    'ICC (n''=1)'};

vartable.Properties.VariableNames = vartnames;

%define parameters for figure size
figwidth = 1230;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;
name = 'All Standard Deviation Components';

%create gui for standard-deviation table
varall_gui= psyrat_newfigure('unit','pix',...
    'position',[1250 600 figwidth figheight],...
    'menub','no',...
    'Tag','psyrat_output',...
    'name',name,...
    'numbertitle','off',...
    'resize','off');

%Print the name of the loaded dataset
uicontrol(varall_gui,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',...
    'Estimated Standard Deviation Components',...
    'Position',[0 row figwidth 25]);

%Start a table
t = uitable('Parent',varall_gui,'Position',...
    [25 100 figwidth-50 figheight-175],...
    'Data',table2cell(vartable));
set(t,'ColumnName',vartnames);
set(t,'ColumnWidth',{200 140 140 140 140 140 140 140});
set(t,'RowName',[]);

%Create a save button that will save the table
uicontrol(varall_gui,'Style','push','fontsize',14,...
    'HorizontalAlignment','center',...
    'String','Save Table',...
    'Position', [figwidth/8 25 figwidth/4 50],...
    'Callback',{@psyrat_savealldifftable,psyrat_data,vartable});

end


function psyrat_savealldifftable(varargin)
%if the user pressed the button to save all of the standard deviation
%components

%parse inputs
psyrat_data = varargin{3};
vartable = varargin{4};

%ask the user where the file should be saved
[savename, savepath] = uiputfile(...
    {'*.csv',...
    'Comma-Separated Value File (.csv)'},...
    'Where would you like to save table?');

%bail out if the user cancelled the save dialog
if isequal(savename,0) || isequal(savepath,0)
    return;
end

%this helper writes a bespoke CSV; normalize the extension so a typed/forced
%.xlsx name cannot produce invalid workbook bytes (F-CM-M4ssxlsx class)
[~,~,ext] = fileparts(savename);
if ~strcmpi(ext,'.csv')
    savename = [savename '.csv'];
end

%save as csv file
fid = fopen(fullfile(savepath,savename),'w');
fprintf(fid,'%s\n','Table Generated on');
fprintf(fid,'%s\n',datestr(clock));
fprintf(fid,'PsyRAT Toolbox v%s\n',psyrat_data.ver);
fprintf(fid,' \n');
fprintf(fid,'%s\n',psyrat_quote_csv_header(...
    sprintf('Dataset: %s',psyrat_data.rel.filename)));
fprintf(fid,'%s\n',psyrat_quote_csv_header(...
    psyrat_iterline(psyrat_data.rel)));
psyrat_fprintf_provenance(fid, psyrat_data);
fprintf(fid,' \n');

%hand-maintained copy of psyrat_viewdiffcomp's vartnames; keep the two in step
fprintf(fid,'%s', strcat('Label,Between-Person Std Dev,BP Covariance Std Dev,',...
    'Between-Trial Std Dev,BT Covariance Std Dev,Error Std Dev,',...
    'Error Covariance Std Dev,ICC (n''=1)'));
fprintf(fid,' \n');

for i = 1:height(vartable)
    formatspec = '%s,%s,%s,%s,%s,%s,%s,%s\n';
    fprintf(fid,formatspec,char(vartable{i,1}),...
        char(vartable{i,2}), char(vartable{i,3}),...
        char(vartable{i,4}), char(vartable{i,5}),...
        char(vartable{i,6}), char(vartable{i,7}),...
        char(vartable{i,8}));
end

fclose(fid);


end
function sem_str = psyrat_sem_string(psyrat_data,gloc,eloc,wit_m,wit_ll,wit_ul)
%Format SEM with the same point + CI style used for other coefficients.

[sem_ll,sem_pt,sem_ul] = psyrat_sem_stats(psyrat_data,gloc,eloc);

if any(isnan([sem_ll sem_pt sem_ul]))
    obs = psyrat_sem_obs(psyrat_data,gloc,eloc);
    if isnan(obs) || obs <= 0
        sem_str = '---';
        return;
    end

    sem_ll = wit_ll ./ sqrt(obs);
    sem_pt = wit_m ./ sqrt(obs);
    sem_ul = wit_ul ./ sqrt(obs);
end

sem_str = sprintf(' %0.2f CI [%0.2f %0.2f]',sem_pt,sem_ll,sem_ul);

end

function [sem_ll,sem_pt,sem_ul] = psyrat_sem_stats(psyrat_data,gloc,eloc)
%Estimate SEM draws using analysis-specific error terms.

sem_ll = NaN;
sem_pt = NaN;
sem_ul = NaN;

ciperc = psyrat_data.relsummary.ciperc;
ciedge = (1-ciperc)./2;
obs = psyrat_sem_obs(psyrat_data,gloc,eloc);

if isnan(obs) || obs <= 0
    return;
end

switch lower(psyrat_data.rel.analysis)
    case 'ic'
        if ~isfield(psyrat_data.relsummary,'data')
            return;
        end
        d = psyrat_data.relsummary.data.g(gloc).e(eloc);
        wp2 = d.sig_e.raw.^2;
        %Absolute (dependability) SEM includes the trial main effect sigma_i;
        %relative (generalizability) SEM uses only the person x trial residual
        %(Rocha Table 2). sig_trl is present for the one-facet designs; older
        %result structs without it fall back to the residual term alone.
        if psyrat_data.relsummary.gcoeff == 1 && ...
                isfield(d,'sig_trl') && ~isempty(d.sig_trl)
            wp2 = wp2 + d.sig_trl.raw.^2;
        end
        sem_draw = sqrt(wp2 ./ obs);

    case 'trt_sserrvar'
        %RC-06 group twin (owner ruling 2026-07-30). 'trt_sserrvar' used to
        %share the ONE-FACET branch below. That was deliberate under RC-20, but
        %explicitly conditional: the group ICC printed beside this SEM was
        %itself a one-facet partition, and pairing a two-facet SEM with it would
        %have put two different estimands in one row. The ICC is now two-facet
        %(see the display branch), so the condition is discharged and the SEM
        %moves with it, in the SAME change as RC-20 required.
        %
        %The n' convention is the plain-'trt' one, not the ICC's: obs comes from
        %psyrat_sem_obs (the trial-count central tendency) and n'_o from
        %relsummary.nocc, so the SEM is the standard error of the score the user
        %actually has -- a mean over several trials -- rather than a single-trial
        %number no other analysis reports. SEM and ICC therefore share the
        %variance PARTITION but not n', which is exactly the arrangement the
        %plain-'trt' row already ships. Under the default noccmode = 1 the n'_o
        %halves agree anyway.
        %
        %Fail open on the same conditions the one-facet branch does, plus the
        %crossed components and reltype: a struct missing any of them returns
        %NaN and psyrat_sem_string falls back rather than erroring and taking
        %the whole variance table with it.
        if ~isfield(psyrat_data.relsummary,'data')
            return;
        end
        d = psyrat_data.relsummary.data.g(gloc).e(eloc);
        if ~isstruct(d) || ~isfield(d,'pop_sdlog') || isempty(d.pop_sdlog)
            return;
        end
        if ~all(isfield(d,{'sig_occ','sig_trl','sig_trlxid','sig_occxid',...
                'sig_trlxocc'})) || ...
                ~isfield(psyrat_data.relsummary,'gcoeff') || ...
                ~isfield(psyrat_data.relsummary,'reltype')
            return;
        end

        %Same n'_o resolution as the {'trt','trt_diff'} case; legacy structs
        %without relsummary.nocc fall back to the single-occasion default.
        if isfield(psyrat_data.relsummary,'nocc') && ...
                ~isempty(psyrat_data.relsummary.nocc)
            nocc = psyrat_data.relsummary.nocc;
        else
            nocc = 1;
        end

        %The residual here is the POPULATION one, exp(pop_sdlog), matching the
        %group ICC and gro_icc. The per-subject residual belongs to the
        %subject-level table, not to this group row. Fields are bare vectors on
        %this design rather than .raw structs, which is the only reason the
        %argument list differs from the {'trt','trt_diff'} case.
        err_term = local_trt_err_term( ...
            psyrat_data.relsummary.gcoeff, psyrat_data.relsummary.reltype, ...
            d.sig_trlxid(:).^2, d.sig_occxid(:).^2, d.sig_trl(:).^2, ...
            d.sig_occ(:).^2, d.sig_trlxocc(:).^2, exp(d.pop_sdlog(:)).^2, ...
            obs, nocc);
        if isempty(err_term)
            return;
        end
        sem_draw = sqrt(err_term);

    case 'ic_sserrvar'
        if ~isfield(psyrat_data.relsummary,'data')
            return;
        end
        d = psyrat_data.relsummary.data.g(gloc).e(eloc);
        %exp(pop_sdlog)^2 is the RELATIVE residual sigma_pi,e^2 once the trial
        %main effect is modeled. The absolute (dependability) SEM keeps the
        %trial main effect sigma_i^2; the relative (generalizability) SEM uses
        %only the residual (Rocha Table 2). sig_trl is present for the one-facet
        %sserr designs; older structs without it fall back to the residual.
        %
        %RC-20 originally routed 'trt_sserrvar' through this ONE-FACET branch as
        %well, so its SEM would match the one-facet group ICC printed beside it.
        %That pairing was explicitly conditional on the ICC staying one-facet.
        %RC-06's group twin was ruled 2026-07-30 and both halves moved together,
        %so 'trt_sserrvar' now has its own case above and this branch serves
        %only the genuinely one-facet 'ic_sserrvar'. Do not merge them back.
        %
        %Fail open on pop_sdlog for the same reason the difference case below
        %does: the display branch at :192 already anticipates an sserr struct
        %WITHOUT gro_sds/pop_sdlog (that is how ic_diff_sserrvar arrives), so
        %the shape is reachable by this file's own reasoning. Returning here
        %restores the psyrat_sem_string fallback for such a struct instead of
        %erroring and taking the whole variance table with it.
        if ~isstruct(d) || ~isfield(d,'pop_sdlog') || isempty(d.pop_sdlog)
            return;
        end
        wp2 = exp(d.pop_sdlog).^2;
        if isfield(psyrat_data.relsummary,'gcoeff') && ...
                psyrat_data.relsummary.gcoeff == 1 && ...
                isfield(d,'sig_trl') && ~isempty(d.sig_trl)
            wp2 = wp2 + d.sig_trl(:).^2;
        end
        sem_draw = sqrt(wp2 ./ obs);

    case {'trt','trt_diff'}
        %Per-condition test-retest SEM (two-facet). For trt_diff the per-event
        %components are the diagonals of the cross-condition covariances, stored
        %as SD draws in relsummary.data.g.e by psyrat_relsummary_trt_diff.
        d = psyrat_data.relsummary.data.g(gloc).e(eloc);
        txp = d.sig_trlxid.raw.^2;
        oxp = d.sig_occxid.raw.^2;
        bt = d.sig_trl.raw.^2;
        bo = d.sig_occ.raw.^2;
        txo = d.sig_trlxocc.raw.^2;
        err = d.sig_err.raw.^2;

        %The SEM must be drawn from the SAME error variance as the reported
        %coefficient, so it uses the SAME effective number of occasions (n'_o)
        %that psyrat_relsummary threaded into psyrat_rel_trt (scientific audit
        %finding F2 / sec.14). The default estimand is a single-occasion score
        %(n'_o = 1); a multi-occasion composite score uses n'_o = k. Single-
        %and multi-occasion coefficients are never paired with each other's
        %SEM. Older result structs without relsummary.nocc fall back to the
        %single-occasion default of 1 (byte-identical to the prior behavior).
        if isfield(psyrat_data.relsummary,'nocc') && ...
                ~isempty(psyrat_data.relsummary.nocc)
            nocc = psyrat_data.relsummary.nocc;
        else
            nocc = 1;
        end

        err_term = local_trt_err_term( ...
            psyrat_data.relsummary.gcoeff, psyrat_data.relsummary.reltype, ...
            txp, oxp, bt, bo, txo, err, obs, nocc);
        if isempty(err_term)
            return;
        end

        sem_draw = sqrt(err_term);

    case {'ic_diff','ic_diff_sserrvar'}
        if ~isfield(psyrat_data.relsummary,'data')
            return;
        end
        d = psyrat_data.relsummary.data.g(gloc).e(eloc);
        %RC-20. Per-event rows of the two difference-score designs. Both store
        %the person x trial,e residual as a LOG standard deviation in b_sigma
        %(psyrat_relsummary's 'sing_diff' and 'sing_diff_sserr' branches), so
        %exp(b_sigma)^2 is the relative residual sigma_pi,e^2 -- the same
        %quantity the sserr branch above reaches through pop_sdlog. The
        %absolute (dependability) SEM adds the trial main effect sigma_i^2; the
        %relative (generalizability) SEM does not (Rocha Table 2).
        %
        %This is the error term psyrat_rel_sing is already handed for the
        %per-event ICC printed in the SAME row, so the SEM and that ICC share
        %an estimand by construction rather than by coincidence.
        %
        %Reading relsummary.gcoeff is right for BOTH analyses because both
        %difference branches assign it FROM diffgcoeff -- but for opposite
        %reasons, and only one of them is benign. Do not collapse the two:
        %Since RC-32 (closed 2026-07-29) BOTH difference analyses build the
        %"Difference-Score Coefficient" popup, via psyrat_startview_sing's
        %isdiffcoeff flag and its `elseif isdiffcoeff` branch, and that popup
        %writes prefs.view.diffgcoeff -- the value psyrat_relfigures passes and
        %the summary turns into relsummary.gcoeff. So reading gcoeff here
        %follows the one control the user can actually see, on both screens.
        %
        %Cited by symbol rather than line: this comment previously named line
        %ranges in psyrat_startview_sing.m, and RC-32's own +90 lines inverted
        %them so each pointed at the opposite branch from the one described.
        %
        %HISTORY, because it explains why this branch reads gcoeff rather than
        %diffgcoeff and must keep doing so: before RC-32 the ic_diff_sserrvar
        %screen matched `sserrvar` first and built a "G-Theory Coefficient"
        %control writing prefs.view.gcoeff, which this route never reads -- the
        %visible control was inert and the coefficient came from a preference
        %that screen never showed. This SEM followed the coefficient the output
        %was ACTUALLY computed with, which was the only self-consistent choice
        %then and is now simply the correct one.
        %
        %Before this case existed both analyses fell through to the psyrat_
        %sem_string fallback (witsd ./ sqrt(obs)): a residual-only SEM with no
        %gcoeff branch, printed under an _Dep/_Gen header that asserted an
        %estimand it did not honor.
        %
        %Fail open on the field itself, not just on relsummary.data. Unlike the
        %cases above -- whose analyses guarantee the draws they read -- a struct
        %can carry one of these two analysis labels while its data.g.e was built
        %by a different branch (the viewer preserves rel.analysis as the
        %estimation identity, and a result can be relabeled without being
        %re-summarized). Returning here restores the old fallback for such a
        %struct; erroring on the missing field would take the whole variance
        %table down over one column.
        if ~isstruct(d) || ~isfield(d,'b_sigma') || ...
                ~isfield(d.b_sigma,'raw') || isempty(d.b_sigma.raw)
            return;
        end
        resid = exp(cell2mat(d.b_sigma.raw));
        wp2 = resid.^2;
        if isfield(psyrat_data.relsummary,'gcoeff') && ...
                psyrat_data.relsummary.gcoeff == 1
            sigi = psyrat_sem_sigtrl(d,resid);
            wp2 = wp2 + sigi.^2;
        end
        sem_draw = sqrt(wp2 ./ obs);

    otherwise
        return;
end

sem_ll = quantile(sem_draw,ciedge);
sem_pt = mean(sem_draw);
sem_ul = quantile(sem_draw,1-ciedge);

end

function err_term = local_trt_err_term(gcoeff,reltype,txp,oxp,bt,bo,txo,err,obs,nocc)
%Two-facet (test-retest) SEM error variance for one gcoeff x reltype cell.
%
%Extracted so the two designs that need it -- {'trt','trt_diff'} and, since
%RC-06's group twin was ruled, 'trt_sserrvar' -- cannot drift apart. This is a
%deliberate SHARED helper rather than one of the toolbox's deliberate twins:
%the six branches are one estimand table, and a duplicated copy would be a
%silent scientific divergence the moment either side was edited.
%
%The SEM must be drawn from the SAME error variance as the reported
%coefficient, so it uses the SAME effective number of occasions (n'_o) that
%psyrat_relsummary threaded into psyrat_rel_trt (scientific audit finding F2 /
%sec.14). The default estimand is a single-occasion score (n'_o = 1); a
%multi-occasion composite score uses n'_o = k. Single- and multi-occasion
%coefficients are never paired with each other's SEM.
%
%Inputs are VARIANCE draws (already squared by the caller), because the two
%callers reach their components through different field shapes: .raw structs
%for trt/trt_diff, bare vectors for trt_sserrvar. Squaring at the call site
%keeps that difference out of here.
%
%Returns [] for an unrecognized gcoeff/reltype combination, which the callers
%turn into the NaN that makes psyrat_sem_string fall back.
err_term = [];

if gcoeff == 1 && reltype == 1
    err_term = (txp ./ obs) + (err ./ (obs*nocc)) +...
        (bt ./ obs) + (txo ./ (obs*nocc));
elseif gcoeff == 1 && reltype == 2
    err_term = (oxp ./ nocc) + (err ./ (obs*nocc)) +...
        (bo ./ nocc) + (txo ./ (obs*nocc));
elseif gcoeff == 1 && reltype == 3
    err_term = (txp ./ obs) + (oxp ./ nocc) +...
        (err ./ (obs*nocc)) + (bt ./ obs) +...
        (bo ./ nocc) + (txo ./ (obs*nocc));
elseif gcoeff == 2 && reltype == 1
    err_term = (txp ./ obs) + (err ./ (obs*nocc));
elseif gcoeff == 2 && reltype == 2
    err_term = (oxp ./ nocc) + (err ./ (obs*nocc));
elseif gcoeff == 2 && reltype == 3
    err_term = (txp ./ obs) + (oxp ./ nocc) + (err ./ (obs*nocc));
end

end


function obs = psyrat_sem_obs(psyrat_data,gloc,eloc)
%Best available trials/observations count for SEM estimation.

obs = NaN;

if any(strcmpi(psyrat_data.rel.analysis,{'trt','trt_diff'})) &&...
        isfield(psyrat_data.relsummary.group(gloc).event(eloc),'rel') &&...
        isfield(psyrat_data.relsummary.group(gloc).event(eloc).rel,'meas')
    obs = psyrat_data.relsummary.group(gloc).event(eloc).rel.meas;
elseif isfield(psyrat_data.relsummary.group(gloc).event(eloc),'dep') &&...
        isfield(psyrat_data.relsummary.group(gloc).event(eloc).dep,'meas')
    obs = psyrat_data.relsummary.group(gloc).event(eloc).dep.meas;
elseif isfield(psyrat_data.relsummary.group(gloc).event(eloc),'trlinfo') &&...
        isfield(psyrat_data.relsummary.group(gloc).event(eloc).trlinfo,'mean')
    obs = psyrat_data.relsummary.group(gloc).event(eloc).trlinfo.mean;
end

if ~isnan(obs) && obs > 0
    obs = max(obs,1);
else
    obs = NaN;
end

end


function sigi = psyrat_sem_sigtrl(d,wp_resid)
%Trial main effect (sigma_i) draws for one stratum, returned as an SD vector
%shaped like the residual SD it will be combined with.
%
%The component is stored under a DIFFERENT field name per analysis, because
%each branch of psyrat_relsummary builds its own data struct: the one-facet
%'sing' branch stores an SD directly (sig_trl.raw); the single-session
%difference stores a per-event SD (sd_trl.raw); the subject-level difference
%stores the trial covariance, whose diagonal is the per-event value as a
%VARIANCE, hence the sqrt. Checked in that order because the difference structs
%carry trl_varcov as well and sd_trl is the unambiguous per-event SD.
%
%The sserr branches are absent from that list on purpose: they store sig_trl as
%a BARE vector rather than a .raw struct, and their two callers above read it
%inline. isfield on a non-struct is false, so such a struct falls through to
%the zeros default here -- correct only because no sserr caller uses this
%resolver for its own sig_trl.
%
%Legacy result structs produced before the trial main effect was modeled carry
%none of these; those return zeros, which reduces the SEM/ICC to the prior
%behavior (sigma_i = 0, dependability == generalizability) rather than erroring.
%
%TWIN, kept deliberately: local_sigtrl in psyrat_report.m resolves the same
%three conventions for the headless reliability-vs-n curve. The two are not
%consolidated into a shared function -- that would edit a working file for
%elegance rather than correctness, the same call the toolbox already made for
%the three harmonic-mean implementations. Change both together.
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


function suffix = psyrat_estimand_suffix(psyrat_data)
%Resolve the decision-type suffix for the estimand-dependent columns (ICC, SEM).
%
%'_Dep' = dependability (absolute error), '_Gen' = generalizability (relative
%error) -- the same tokens psyrat_dod_observedt already emits, and valid MATLAB
%identifiers so no VariableNamingRule='preserve' is needed.
%
%Difference-score analyses carry their own selector (diffgcoeff), so prefer it
%where present, then fall back to the numeric selector for structs saved before
%gcoeff_name existed. Returns '' only when neither is available -- defensive,
%and in practice unreachable for the analyses that reach this table, since
%psyrat_sem_stats already requires relsummary.gcoeff.

suffix = '';
if ~isfield(psyrat_data,'relsummary'); return; end
rs = psyrat_data.relsummary;

%difference-score views are driven by diffgcoeff; everything else by gcoeff
isdiff = isfield(psyrat_data,'rel') && isfield(psyrat_data.rel,'analysis') && ...
    any(strcmp(psyrat_data.rel.analysis,{'ic_diff','trt_diff','ic_diff_sserrvar'}));

name = '';
if isdiff && isfield(rs,'diffgcoeff_name') && ~isempty(rs.diffgcoeff_name)
    name = rs.diffgcoeff_name;
elseif isfield(rs,'gcoeff_name') && ~isempty(rs.gcoeff_name)
    name = rs.gcoeff_name;
end

if ~isempty(name)
    switch lower(name)
        case 'dep'; suffix = '_Dep';
        case 'gen'; suffix = '_Gen';
    end
    return
end

%fall back to the numeric selector when only it was stored
num = [];
if isdiff && isfield(rs,'diffgcoeff') && ~isempty(rs.diffgcoeff)
    num = rs.diffgcoeff;
elseif isfield(rs,'gcoeff') && ~isempty(rs.gcoeff)
    num = rs.gcoeff;
end
if isequal(num,1)
    suffix = '_Dep';
elseif isequal(num,2)
    suffix = '_Gen';
end

end


function names = psyrat_variancet_displaynames(names)
%Convert table variable names to the human-readable headers shown in the GUI
%and written to CSV/XLSX. Single source of truth for all three surfaces, which
%were previously three independent hand-maintained lists.
%
%Order matters: the _Dep/_Gen expansion must run BEFORE underscores become
%spaces, or the tokens are no longer matchable.

names = strrep(names,'_Dep',' (dependability)');
names = strrep(names,'_Gen',' (generalizability)');
names = strrep(names,'_',' ');
names = strrep(names,'StdDev','Std Dev');

end
