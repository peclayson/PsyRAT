function psyrat_criterionfigures(varargin)
%Generate criterion-score-specific reliability plot/table outputs.
%
%psyrat_criterionfigures('psyrat_prefs',psyrat_prefs,...
%    'psyrat_data',psyrat_data,'analysis','sing')

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

if isempty(psyrat_data)
    error('varargin:psyrat_data',... %Error code and associated error
        strcat('WARNING: psyrat_data not specified \n\n',...
        'Please input psyrat_data (PsyRAT Toolbox data structure array).\n'));
end

if isempty(psyrat_prefs) || ~isstruct(psyrat_prefs)
    psyrat_prefs = psyrat_defaults;
end
if ~isfield(psyrat_prefs,'view')
    psyrat_prefs.view = struct;
end

analysis = 'sing';
ind = find(strcmpi('analysis',varargin),1);
if ~isempty(ind)
    analysis = lower(char(string(varargin{ind+1})));
end
if ~any(strcmp(analysis,{'sing','trt'}))
    error('varargin:analysis',... %Error code and associated error
        strcat('WARNING: analysis is invalid. Valid values are ''sing'' or ''trt''\n'));
end

%G42. Scripted-caller guard for the five GUI-adjacent designs the filing
%named -- NOT a completeness claim over every label the toolbox can stamp.
%The GUI never sends these five here: psyrat_startview_sing gates its
%criterion button on ~sserrvar && ~isdiff && hasEstimatedMu, and
%psyrat_startview_trt hides the button for trt_sserrvar/trt_diff (G38 --
%a two-label BLACKLIST, not a plain-'trt' whitelist; other trt-family
%labels reach that viewer only because psyrat_startview's router keeps
%them out, and this guard backstops the router). But a
%direct call used to fail wrongly on each: trt_sserrvar estimates a
%per-participant residual instead of the group-level components
%psyrat_cutscore_trt consumes (raw MATLAB:nonExistentField), the trt_diff
%D-study stores no mu (same, pinned in TestTrtDiffCurveWiring), and the
%one-facet difference shapes reached the RC-05 sd_id/b_sigma stand-in,
%whose cut-score is not an estimand. Checked on rel.analysis rather than
%the route argument so a mislabeled route cannot slip a blocked design
%past it.
if isfield(psyrat_data,'rel') && isfield(psyrat_data.rel,'analysis') && ...
        any(strcmp(psyrat_data.rel.analysis, ...
        {'ic_sserrvar','ic_diff','ic_diff_sserrvar','trt_sserrvar','trt_diff'}))
    error('criterion:unsupportedDesign',... %Error code and associated error
        strcat('WARNING: Criterion outputs are unavailable for this design (%s)\n\n',...
        'Criterion cut-score outputs are supported for one-facet (''ic'') and\n',...
        'test-retest (''trt'') results only. Subject-level and difference-score\n',...
        'result structures do not carry the components the cut-score kernels\n',...
        'consume.\n'), psyrat_data.rel.analysis);
end

%G50/G55. The categorical refusal for the OTHER unsupported families: DoD
%(ic_dodiff), the sixteen-label dynrel family, and the splits family
%(ic_splits and, since G55, trt_splits). The DoD/dynrel/ic_splits structs
%do not carry the one-facet / test-retest components this screen's
%rebuild and cut-score kernels consume, so a scripted caller used to die
%on a raw error deep inside the rebuild (MATLAB:nonExistentField at the
%first missing component on the 'sing' route -- REL.out.sig_e for
%splits/DoD, REL.groups for the dynrel fixtures -- or
%MATLAB:cellRefFromNonCell on the 'trt' route's char labels). Refused
%HERE, by family, for a reason the adversarial wave forced (2026-08-29):
%an earlier draft translated the rebuild's nonExistentField errors into
%this refusal, which mislabeled SUPPORTED designs on a mismatched route
%(plain trt data on the default 'sing' route dies on sig_e too) and
%masked genuine field problems. A family check cannot capture those:
%mismatched supported designs and legacy label-free structs (family
%'sing' by default) fall through and fail with the raw, factually true
%error, exactly as before G50.
%
%trt_splits is refused since G55 (owner ruling 2026-08-30, reversing the
%G50-era constraint that left it unrefused). Its failure mode is the
%opposite of its family-mates: the case-24 out carries the FULL trt
%component set (mu plus eight sigma components) and the trt-route rebuild
%COMPLETED -- but the numbers it produced are not the splits estimand.
%The trt route passes the per-split-item sig_err unscaled as the
%per-observation residual, dropping the sig_posxid term and the
%items-per-split weighting psyrat_splits_summary folds into its error
%(err = sqrt(sig_posxid^2 + sig_err^2/nbar_i)), and runs a D-study
%projection on a design case 24 documents as observed-design-only. (The
%sig_trl -> bt "trial slot" binding itself is NOT one of the defects:
%psyrat_splits_summary makes the same split-facet-onto-trial-slot mapping
%on purpose, and the unit labels already said "split" via rel.splits.)
%Completing without error was never validation (the pre-G55 guard comment
%said so explicitly); now it does not complete.
%
%Scope note: this refusal is keyed on the NONPARALLEL labels. PARALLEL
%splits runs are stamped 'ic'/'trt' with rel.splits = 2 and legitimately
%reach this screen -- equal n_i makes the unweighted kernels the right
%ones there, and the unit labels switch to "split" (see psyrat_unitlabels
%and the usplit reads below).
if isfield(psyrat_data,'rel') && isfield(psyrat_data.rel,'analysis')
    fam = psyrat_analysis_family(psyrat_data.rel.analysis);
    if any(strcmp(fam, {'dod','dynrel','splits'}))
        error('criterion:unsupportedDesign',... %Error code and associated error
            strcat('WARNING: Criterion outputs are unavailable for this design (%s)\n\n',...
            'Criterion cut-score outputs are supported for one-facet (''ic'') and\n',...
            'test-retest (''trt'') results only. Difference-of-differences,\n',...
            'dynamic-reliability, and nonparallel one-facet splits result\n',...
            'structures do not carry the components the cut-score kernels\n',...
            'consume, and the nonparallel test-retest splits rebuild would\n',...
            'ignore the items-per-split weighting the splits estimand\n',...
            'requires.\n'), psyrat_data.rel.analysis);
    end
end

settings = local_settings(psyrat_prefs,analysis);

%relerr stays empty when psyrat_data already carries a relsummary and no
%rebuild runs, which makes the two error branches at the end of this function
%inert on that path (they are guarded by isfield).
relerr = struct;

if ~isfield(psyrat_data,'relsummary') || ~isfield(psyrat_data.relsummary,'data')
    [psyrat_data, relerr] = local_build_relsummary(psyrat_data, psyrat_prefs, analysis, settings);
    if isfield(relerr,'nogooddata') && relerr.nogooddata == 1
        return;
    end
end

if ~isfield(psyrat_data.relsummary,'ciperc')
    psyrat_data.relsummary.ciperc = 0.95;
end

if strcmp(analysis,'sing')
    psyrat_data.relsummary.gcoeff = settings.gcoeff;
    if settings.gcoeff == 1
        psyrat_data.relsummary.gcoeff_name = 'dep';
    else
        psyrat_data.relsummary.gcoeff_name = 'gen';
    end
elseif strcmp(analysis,'trt')
    psyrat_data.relsummary.gcoeff = settings.gcoeff;
    psyrat_data.relsummary.reltype = settings.reltype;
    if settings.gcoeff == 1
        psyrat_data.relsummary.gcoeff_name = 'dep';
    else
        psyrat_data.relsummary.gcoeff_name = 'gen';
    end
end

if settings.plotcriterion == 1
    local_plotcriterion(psyrat_data, analysis, settings);
end

if settings.tablecriterion == 1
    local_tablecriterion(psyrat_data, analysis, settings);
end

%RC-33. Mirror the two post-rebuild error branches psyrat_relfigures raises
%(:745-751 and :754-760). This function rebuilds the whole D-study itself, so it
%recomputes trlcutoff, trlmax and the eventgoodids/eventbadids partition, but it
%previously reported none of it -- the criterion screens failed more quietly
%than the view screen the user came from.
%
%These sit AFTER the plot and table are drawn, not next to the nogooddata guard
%above. errordlg is non-modal, so a dialog raised before the output figures
%would be buried underneath them, which defeats the point; psyrat_relfigures
%raises both of its dialogs last for the same reason. Neither branch returns:
%the criterion output is still drawn, exactly as on the view screen.
%
%Titled per purpose but deliberately still NON-modal (2026-08-17 notification
%ruling): raising them last keeps the later-drawn output from burying them at
%creation, and a modal notice would lock the user out of inspecting the very
%output the caveat describes. Whether they also rise above the already-open
%windows at real geometries is the live click-through's to verify (B21: a
%non-modal dialog does not raise above an existing figure).
%
%Not gated on analysis, so the one-facet route gets them too. That is
%deliberate -- one-facet data has no occasion facet and so cannot diverge from
%its parent, but psyrat_relfigures already raises these dialogs for it, and a
%one-facet trial-cutoff failure is worth reporting on either screen.
%
%Reachability, re-measured twice on 2026-08-29 (G42, then corrected by the
%adversarial wave). The trlcutoff branch fires on live routes (pinned via
%the trt rebuild in TestSummaryAndPresentationFunctions). The trlmax branch
%is ALSO still live through this function: when a found cutoff excludes
%every participant, the relsummary routes fall back to the bad-id roster
%for the trial table, so trlinfo.max is real and trlcutoff exceeds it (the
%trt producer is pinned by
%testTrtExtrapolatedCutoffWithNoRetainedIdsStillWritesGoodN; the plain-sing
%route has the same fallback shape). The "taken over the good ids, so
%tautological" argument holds only while somebody is retained. What G42
%closed is ONE producer: the 'sing_diff' derived-count comparison, reached
%only by ic_diff data, whose cut-score came from the RC-05 non-estimand
%stand-in -- that closure is pinned in
%testCriterionRouteRefusesTheVirtualDiffRel.
%
%analysis-unit label fragments (trial vs split); defaults to "trial" so
%non-splits output is byte-identical.
usplit = 1;
if isfield(psyrat_data,'rel') && isfield(psyrat_data.rel,'splits')
    usplit = psyrat_data.rel.splits;
end
L = psyrat_unitlabels(usplit);

%show an error if the dependability threshold was never met for the cutoff
if isfield(relerr,'trlcutoff') && relerr.trlcutoff == 1
    errorstr = {};
    errorstr{end+1} = [L.Unit ' cutoffs for adequate dependability could not be calculated'];
    errorstr{end+1} = ['Data are too variable or there are not enough ' L.units];
    errorstr = [errorstr local_divergencenote(analysis, settings, L)];

    errordlg(errorstr,'Cutoff not calculable');
end

%show an error if none of the data had enough trials to meet the cutoff
if isfield(relerr,'trlmax') && relerr.trlmax == 1
    errorstr = {};
    errorstr{end+1} = ['Not enough ' L.units ' are present in the current data'];
    errorstr{end+1} = 'Cutoffs represent an extrapolation beyond the data';
    errorstr = [errorstr local_divergencenote(analysis, settings, L)];

    errordlg(errorstr,'Extrapolation beyond data');
end

end

function settings = local_settings(psyrat_prefs,analysis)
settings = struct;

if ~isfield(psyrat_prefs.view,'criterioncutoff')
    psyrat_prefs.view.criterioncutoff = 0;
end
if ~isfield(psyrat_prefs.view,'plotcriterion')
    psyrat_prefs.view.plotcriterion = 1;
end
if ~isfield(psyrat_prefs.view,'tablecriterion')
    psyrat_prefs.view.tablecriterion = 1;
end
if ~isfield(psyrat_prefs.view,'reltype')
    psyrat_prefs.view.reltype = 1;
end
if ~isfield(psyrat_prefs.view,'plotdepline')
    psyrat_prefs.view.plotdepline = 2;
end
if ~isfield(psyrat_prefs.view,'plotrelline')
    psyrat_prefs.view.plotrelline = 2;
end
if ~isfield(psyrat_prefs.view,'ntrials')
    psyrat_prefs.view.ntrials = 50;
end
if ~isfield(psyrat_prefs.view,'depcentmeas')
    psyrat_prefs.view.depcentmeas = 1;
end
if ~isfield(psyrat_prefs.view,'relcentmeas')
    psyrat_prefs.view.relcentmeas = 1;
end
%RC-33. Read ONLY so the error dialogs can say why this screen may disagree
%with the parent view screen. It is deliberately not passed to
%psyrat_relsummary or to either psyrat_cutscore_trt call -- see the block in
%local_build_relsummary, which is still true as written.
if ~isfield(psyrat_prefs.view,'noccmode') || isempty(psyrat_prefs.view.noccmode)
    psyrat_prefs.view.noccmode = 1;
end

settings.cutoff = psyrat_prefs.view.criterioncutoff;
settings.plotcriterion = psyrat_prefs.view.plotcriterion;
settings.tablecriterion = psyrat_prefs.view.tablecriterion;
settings.gcoeff = 1; %criterion outputs are always dependability (absolute)
settings.reltype = psyrat_prefs.view.reltype;
settings.noccmode = psyrat_prefs.view.noccmode;
settings.ntrials = max(2, round(psyrat_prefs.view.ntrials));

if strcmp(analysis,'sing')
    settings.line = psyrat_prefs.view.plotdepline;
    settings.centmeas = psyrat_prefs.view.depcentmeas;
else
    settings.line = psyrat_prefs.view.plotrelline;
    settings.centmeas = psyrat_prefs.view.relcentmeas;
end

if isnan(settings.cutoff)
    error('varargin:criterioncutoff',... %Error code and associated error
        strcat('WARNING: criterion cutoff must be numeric\n'));
end

end

function [psyrat_data, relerr] = local_build_relsummary(psyrat_data, psyrat_prefs, analysis, settings)
%This runs BEFORE the caller's relsummary.ciperc guard -- it is what builds
%relsummary in the first place -- so the stored width cannot simply be read
%here. Honor one when psyrat_data already carries it (a scripted run that set
%'CI' and is now being re-summarized at the criterion estimand), and fall back
%to the documented .95 otherwise. There is no GUI control for CI on this path.
ciperc = 0.95;
if isfield(psyrat_data,'relsummary') && isfield(psyrat_data.relsummary,'ciperc')
    ciperc = psyrat_data.relsummary.ciperc;
end
if strcmp(analysis,'sing')
    rel_analysis = 'sing';
    if isfield(psyrat_data,'rel') && isfield(psyrat_data.rel,'analysis')
        switch lower(psyrat_data.rel.analysis)
            case 'ic'
                rel_analysis = 'sing';
            case 'ic_diff'
                rel_analysis = 'sing_diff';
            case 'ic_sserrvar'
                rel_analysis = 'sing_sserr';
            otherwise
                rel_analysis = 'sing';
        end
    end

    if ~isfield(psyrat_prefs.view,'depvalue')
        psyrat_prefs.view.depvalue = .8;
    end
    if ~isfield(psyrat_prefs.view,'meascutoff')
        psyrat_prefs.view.meascutoff = 2;
    end
    if ~isfield(psyrat_prefs.view,'depcentmeas')
        psyrat_prefs.view.depcentmeas = 1;
    end
    if ~isfield(psyrat_prefs.view,'diffgcoeff')
        psyrat_prefs.view.diffgcoeff = 1;
    end

    [psyrat_data, relerr] = psyrat_relsummary('psyrat_data',psyrat_data,...
        'analysis',rel_analysis,...
        'depcutoff',psyrat_prefs.view.depvalue,...
        'meascutoff',psyrat_prefs.view.meascutoff,...
        'depcentmeas',psyrat_prefs.view.depcentmeas,...
        'gcoeff',settings.gcoeff,...
        'diffgcoeff',psyrat_prefs.view.diffgcoeff,...
        'CI',ciperc);
else
    if ~isfield(psyrat_prefs.view,'relvalue')
        psyrat_prefs.view.relvalue = .8;
    end
    if ~isfield(psyrat_prefs.view,'meascutoff')
        psyrat_prefs.view.meascutoff = 2;
    end
    if ~isfield(psyrat_prefs.view,'relcentmeas')
        psyrat_prefs.view.relcentmeas = 1;
    end

    %Neither branch below passes noccmode/nocc, so local_parse_nocc defaults to
    %noccmode = 1 and local_resolve_effnocc returns n'_o = 1. That is deliberate
    %and matches the 'nocc',1 the two psyrat_cutscore_trt calls use, so the trial
    %cutoffs, the participant inclusion and the cut-scores on this screen are all
    %computed on one estimand (RC-29, owner ruling 2026-07-30; rationale in
    %local_occnote). Do NOT thread psyrat_prefs.view.nocc in here without also
    %changing both cut-score calls -- half of it produces a mixed estimand.
    %
    %Known consequence of that pin: under noccmode = 2 this rebuild can retain a
    %different participant set than the parent view screen did, because it
    %recomputes the trial cutoffs and the eventgoodids/eventbadids partition at
    %n'_o = 1. It used to do so silently -- unlike psyrat_relfigures, this
    %function ignored both relerr.trlcutoff and relerr.trlmax. RC-33 mirrors
    %those two branches at the end of the main function and adds a line naming
    %the divergence (local_divergencenote). The numbers are unchanged; only the
    %failure modes are now reported.
    %
    %two-facet difference-score data route through the 'trt_diff' branch with the
    %difference-score coefficient choice (diffgcoeff, default dependability).
    if isfield(psyrat_data,'rel') && isfield(psyrat_data.rel,'analysis') && ...
            strcmp(psyrat_data.rel.analysis,'trt_diff')
        if ~isfield(psyrat_prefs.view,'diffgcoeff')
            psyrat_prefs.view.diffgcoeff = 1;
        end
        [psyrat_data, relerr] = psyrat_relsummary('psyrat_data',psyrat_data,...
            'analysis','trt_diff',...
            'gcoeff',settings.gcoeff,...
            'diffgcoeff',psyrat_prefs.view.diffgcoeff,...
            'reltype',settings.reltype,...
            'relcutoff',psyrat_prefs.view.relvalue,...
            'meascutoff',psyrat_prefs.view.meascutoff,...
            'relcentmeas',psyrat_prefs.view.relcentmeas,...
            'CI',ciperc);
    else
        [psyrat_data, relerr] = psyrat_relsummary('psyrat_data',psyrat_data,...
            'analysis','trt',...
            'gcoeff',settings.gcoeff,...
            'reltype',settings.reltype,...
            'relcutoff',psyrat_prefs.view.relvalue,...
            'meascutoff',psyrat_prefs.view.meascutoff,...
            'relcentmeas',psyrat_prefs.view.relcentmeas,...
            'CI',ciperc);
    end
end

end

function local_plotcriterion(psyrat_data, analysis, settings)
%analysis-unit label fragments (trial vs split); defaults to "trial" so
%non-splits output is byte-identical.
usplit = 1;
if isfield(psyrat_data,'rel') && isfield(psyrat_data.rel,'splits')
    usplit = psyrat_data.rel.splits;
end
L = psyrat_unitlabels(usplit);
[ngroups, nevents, gnames, enames, modelcase] = local_groupinfo(psyrat_data);

x = 1:settings.ntrials;

if nevents > 2
    xplots = ceil(sqrt(nevents));
    yplots = ceil(sqrt(nevents));
elseif nevents == 2
    xplots = 1;
    yplots = 2;
else
    xplots = 1;
    yplots = 1;
end

fig = figure;
fig.Tag = 'psyrat_output';
set(fig,'NumberTitle','Off');
fig.Position = [125 630 940 470];

rel_name = local_relname(psyrat_data.relsummary);
pref = ['Criterion ' rel_name ' v ' L.NumberOf];
ylbl = rel_name;

for eloc = 1:nevents
    plotrel = nan(1+settings.ntrials-1, ngroups);
    for gloc = 1:ngroups
        d = psyrat_data.relsummary.data.g(gloc).e(eloc);
        if strcmp(analysis,'sing')
            [ll,pt,ul] = local_cutscore_sing(d, settings.cutoff, [1 settings.ntrials], ...
                psyrat_data.relsummary.ciperc);
        else
            %Cut-scores are always absolute (dependability); force gcoeff = 1
            %regardless of settings.gcoeff. psyrat_cutscore_trt rejects gcoeff = 2.
            %'nocc',1 below is pinned for a DIFFERENT reason -- a ruling, not a
            %validity constraint. See local_occnote, which states it on screen.
            [ll,pt,ul] = psyrat_cutscore_trt( ...
                'gcoeff',1,...
                'reltype',settings.reltype,...
                'bp',d.sig_id.raw,...
                'bo',d.sig_occ.raw,...
                'bt',d.sig_trl.raw,...
                'txp',d.sig_trlxid.raw,...
                'oxp',d.sig_occxid.raw,...
                'txo',d.sig_trlxocc.raw,...
                'err',d.sig_err.raw,...
                'mu',d.mu.raw,...
                'cut',settings.cutoff,...
                'obs',[1 settings.ntrials],...
                'nocc',1,...
                'CI',psyrat_data.relsummary.ciperc);
        end

        switch settings.line
            case 1
                plotrel(:,gloc) = ll(:);
                plottitle = 'Lower Limit';
            case 2
                plotrel(:,gloc) = pt(:);
                plottitle = 'Point Estimate';
            case 3
                plotrel(:,gloc) = ul(:);
                plottitle = 'Upper Limit';
            otherwise
                plotrel(:,gloc) = pt(:);
                plottitle = 'Point Estimate';
        end
    end

    subplot(yplots,xplots,eloc);
    h = plot(x,plotrel,'LineWidth',1.4);
    set(gca,'FontSize',14);
    axis([1 settings.ntrials 0 1]);
    xlabel('Number of Observations','FontSize',14);
    ylabel(ylbl,'FontSize',14);
    if modelcase ~= 1 && modelcase ~= 3
        leg = legend(h,gnames{:},'Location','southeast');
        set(leg,'FontSize',12,'Interpreter','none');
        legend('boxoff');
    end
    if ~strcmpi(enames{eloc},'none')
        %Interpreter none here and on the legend above: an underscored event
        %or group name (inc_cor) must print literally, where the default TeX
        %interpreter would subscript the character after the underscore.
        title(enames{eloc},'FontSize',18,'Interpreter','none');
    else
        title('Measurement','FontSize',18);
    end
end

%The occasion note is appended only when there is one (two-facet data), so the
%one-facet window title is unchanged to the character.
occnote = local_occnote(analysis);
if isempty(occnote)
    fig.Name = sprintf('%s (%s, Cutoff = %0.3f)', pref, plottitle, settings.cutoff);
else
    fig.Name = sprintf('%s (%s, Cutoff = %0.3f, %s)', ...
        pref, plottitle, settings.cutoff, occnote);
end

end

function local_tablecriterion(psyrat_data, analysis, settings)
%analysis-unit label fragments (trial vs split); defaults to "trial" so
%non-splits output is byte-identical.
usplit = 1;
if isfield(psyrat_data,'rel') && isfield(psyrat_data.rel,'splits')
    usplit = psyrat_data.rel.splits;
end
L = psyrat_unitlabels(usplit);
[ngroups, nevents, gnames, enames, modelcase] = local_groupinfo(psyrat_data);
rel_name = local_relname(psyrat_data.relsummary);

labels = {};
cutoffs = [];
trlused = [];
goodns = [];
muvals = [];
ptvals = [];
llvals = [];
ulvals = [];

for gloc = 1:ngroups
    for eloc = 1:nevents
        d = psyrat_data.relsummary.data.g(gloc).e(eloc);
        ev = psyrat_data.relsummary.group(gloc).event(eloc);

        switch modelcase
            case 1
                labels{end+1,1} = 'Measurement'; %#ok<AGROW>
            case 2
                labels{end+1,1} = gnames{gloc}; %#ok<AGROW>
            case 3
                labels{end+1,1} = enames{eloc}; %#ok<AGROW>
            otherwise
                labels{end+1,1} = [gnames{gloc} ' - ' enames{eloc}]; %#ok<AGROW>
        end

        obs = local_obscount(ev, settings.centmeas);
        trlused(end+1,1) = obs; %#ok<AGROW>
        cutoffs(end+1,1) = settings.cutoff; %#ok<AGROW>
        goodns(end+1,1) = local_goodn(ev); %#ok<AGROW>
        muvals(end+1,1) = local_muvalue(d); %#ok<AGROW>

        if strcmp(analysis,'sing')
            [ll,pt,ul] = local_cutscore_sing(d, settings.cutoff, obs, ...
                psyrat_data.relsummary.ciperc);
        else
            %Cut-scores are always absolute (dependability); force gcoeff = 1
            %regardless of settings.gcoeff. psyrat_cutscore_trt rejects gcoeff = 2.
            %'nocc',1 below is pinned for a DIFFERENT reason -- a ruling, not a
            %validity constraint. See local_occnote, which states it on screen.
            [ll,pt,ul] = psyrat_cutscore_trt( ...
                'gcoeff',1,...
                'reltype',settings.reltype,...
                'bp',d.sig_id.raw,...
                'bo',d.sig_occ.raw,...
                'bt',d.sig_trl.raw,...
                'txp',d.sig_trlxid.raw,...
                'oxp',d.sig_occxid.raw,...
                'txo',d.sig_trlxocc.raw,...
                'err',d.sig_err.raw,...
                'mu',d.mu.raw,...
                'cut',settings.cutoff,...
                'obs',obs,...
                'nocc',1,...
                'CI',psyrat_data.relsummary.ciperc);
        end

        ptvals(end+1,1) = pt; %#ok<AGROW>
        llvals(end+1,1) = ll; %#ok<AGROW>
        ulvals(end+1,1) = ul; %#ok<AGROW>
    end
end

%RC-36. n_Included is the retained participant count for this group x event.
%
%It exists because these screens rebuild the whole D-study at n'_o = 1 (see
%local_occnote) and therefore recompute the eventgoodids/eventbadids partition.
%Under noccmode = 2 the parent view screen selected a k-occasion composite and
%can have retained a DIFFERENT set of participants. RC-33 reports the two cases
%where that rebuild fails; the ordinary case -- both screens find a cutoff and
%simply find different ones -- raises no flag, and before this column there was
%no number on either criterion surface a user could compare against the screen
%they came from.
%
%UNCONDITIONAL, deliberately, and NOT gated the way local_divergencenote is.
%local_save_table writes this table to xlsx/csv, so gating the column on
%analysis/noccmode would make the EXPORTED schema vary by route, which is worse
%for downstream scripts than one stable schema. The cost is explicit and
%accepted: this changes the exported column set on every route.
%
%Named and positioned to match the parent screens exactly -- 'n_Included' in
%column 2, as in psyrat_depoverallt and psyrat_trt_reloverallt -- because
%comparing the two Ns is the whole point. It is NOT named off psyrat_unitlabels:
%goodn counts participants, whereas L.Units is trials/splits and is already
%spoken for by the [L.Units '_Used'] column beside it.
tbl = table(labels, goodns, cutoffs, trlused, muvals, ptvals, llvals, ulvals, ...
    'VariableNames', ...
    {'Label','n_Included','Criterion_Cutoff',[L.Units '_Used'],'Mu',...
    [rel_name '_PointEstimate'], [rel_name '_LowerCI'], [rel_name '_UpperCI']});

%widened from 700 when n_Included was added. The ColumnWidth list below already
%summed to more than the table's own width, so this table scrolled horizontally
%even before an eighth column existed.
%
%The arithmetic: sum(ColumnWidth) = 800, and the table Position below leaves
%2*25 px of margin, so figwidth must be at least 850. **850 is not enough.** At
%exactly 850 the table is 800 px for 800 px of columns and the uitable's own
%grid/border eats into the last column -- the final digit of the Upper CI value
%renders clipped. That is invisible to a `sum(ColumnWidth) <= table width`
%check, which passes at equality; it was caught by looking at the window. 870
%leaves 20 px of slack and the last column renders whole.
%
%Every control in this figure is positioned off figwidth, so this one number is
%the whole layout change. If a column is ever added or widened, re-derive this
%AND look at the result -- there is no slack to borrow from.
figwidth = 870;
figheight = 430;

%The occasion note is appended only when there is one (two-facet data), so the
%one-facet window title is unchanged to the character.
occnote = local_occnote(analysis);
figname = ['Overall Criterion ' rel_name ' Summary'];
if ~isempty(occnote)
    figname = [figname ' (' occnote ')'];
end

fig = psyrat_newfigure('unit','pix',...
    'position',[1100 200 figwidth figheight],...
    'menub','no',...
    'Tag','psyrat_output',...
    'name',figname,...
    'numbertitle','off',...
    'resize','off');

uicontrol(fig,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String',sprintf('Criterion %s Summary (Cutoff = %0.3f)', rel_name, settings.cutoff),...
    'Position',[0 figheight-55 figwidth 24]);

%Occasion estimand, on its own line rather than appended to the heading above:
%the heading is already ~47 characters at fontsize 16 in a 700 px control, and
%appending would risk clipping the cutoff value. This sits in the 20 px gap
%between the heading (bottom edge figheight-55 = 375) and the table (top edge
%95 + figheight-170 = 355), so nothing else moves.
if ~isempty(occnote)
    uicontrol(fig,'Style','text','fontsize',11,...
        'HorizontalAlignment','center',...
        'String',['Cut-scores computed at ' occnote],...
        'Position',[0 figheight-75 figwidth 16]);
end

t = uitable('Parent',fig,'Position',...
    [25 95 figwidth-50 figheight-170],...
    'Data',table2cell(tbl));
set(t,'ColumnName',{'Label' 'n Included' 'Criterion Cutoff' [L.Units ' Used']...
    'Mu' ...
    [rel_name ' Point Estimate'] [rel_name ' Lower CI'] [rel_name ' Upper CI']});
set(t,'ColumnWidth',{160 75 105 85 90 95 95 95});
set(t,'RowName',[]);

uicontrol(fig,'Style','push','fontsize',13,...
    'HorizontalAlignment','center',...
    'String','Save Table',...
    'Position', [figwidth/3 24 figwidth/3 50],...
    'Callback',{@local_save_table, psyrat_data, tbl, settings.cutoff});

end

function [ll, pt, ul] = local_cutscore_sing(d, cutoff, obs, ciperc)
%Cut-scores are always absolute (dependability), so this always calls
%psyrat_cutscore_sing with gcoeff = 1 regardless of settings.gcoeff.
if isfield(d,'sig_u') && isfield(d,'sig_e') && isfield(d,'mu')
    bp = d.sig_u.raw;
    mu = d.mu.raw;
    %Total within-person variance for the one-facet cut-score (Rocha 2026
    %Table 2): wp^2 = sigma_pi,e^2 + sigma_i^2. Building wp from sig_e and
    %the trial main effect sig_trl keeps sigma_i^2 in the absolute
    %(dependability) error. This mirrors the relsummary wiring
    %(psyrat_relsummary.m:458). Legacy result structs produced before the
    %one-facet trial main effect was modeled carry no sig_trl; default it to
    %residual-only so the wiring reduces to the prior behavior (sigma_i = 0).
    if isfield(d,'sig_trl') && isfield(d.sig_trl,'raw') && ~isempty(d.sig_trl.raw)
        wp = sqrt(d.sig_e.raw.^2 + d.sig_trl.raw.^2);
    else
        wp = d.sig_e.raw;
    end
elseif isfield(d,'sd_id') && isfield(d,'b_sigma')
    %NOT REACHED BY ANY SHIPPED ROUTE (RC-05), and it must not become reachable
    %as written. mu below is mean(bp): the posterior mean of a between-person
    %STANDARD DEVIATION, standing in for the universe-score location in
    %offsq = (mu - cut)^2. An SD is not a location, so the resulting criterion
    %coefficient is not an estimand. It is left in place rather than fixed
    %because deriving a real mu for a difference-score location-scale model is
    %its own piece of work (compare RC-18, the live mu defect on the DoD route).
    %
    %Why nothing reaches it. This sd_id/b_sigma struct shape is produced by
    %exactly two builders, both difference-score: case 'sing_diff' (ic_diff) and
    %psyrat_relsummary_sing_diff_sserr (ic_diff_sserrvar), both in
    %psyrat_relsummary. psyrat_startview_sing gates the criterion button with
    %
    %  showCriterionButton = ~sserrvar && ~isdiff && hasEstimatedMu
    %
    %and the two producers are suppressed by DIFFERENT terms: ic_diff by
    %~isdiff, which is true for ic_diff ONLY, and ic_diff_sserrvar by ~sserrvar,
    %which covers both sserrvar analyses. Neither term alone covers both
    %producers. Anyone simplifying that gate on the assumption that ~isdiff
    %catches every difference analysis would expose this line for
    %ic_diff_sserrvar. (hasEstimatedMu is RC-18's mu_source guard, added later.)
    %Scripted callers used to bypass the gate entirely; since G42 the guard at
    %the top of this function blocks ic_diff/ic_diff_sserrvar (and the other
    %unsupported designs) by rel.analysis, so reaching this branch now takes a
    %doctored struct that carries the sd_id/b_sigma shape without one of the
    %blocked labels. See RC-05 in PSYRAT_AUDIT_FINDINGS.md.
    bp = cell2mat(d.sd_id.raw);
    wp = exp(cell2mat(d.b_sigma.raw));
    mu = mean(bp) * ones(size(bp));
    %The subject-specific error-variance (location-scale) one-facet model
    %does not separate a trial main effect, so there is no sigma_i^2 to add
    %(consistent with the sd_id/b_sigma path in psyrat_relsummary.m).
else
    error('varargin:analysis',... %Error code and associated error
        strcat('WARNING: Criterion outputs are unavailable for this one-facet data structure.\n'));
end

[ll, pt, ul] = psyrat_cutscore_sing( ...
    'gcoeff',1,...
    'bp',bp,...
    'wp',wp,...
    'mu',mu,...
    'cut',cutoff,...
    'obs',obs,...
    'CI',ciperc);

end

function obs = local_obscount(ev, centmeas)
obs = 1;
if isfield(ev,'trlinfo')
    if centmeas == 2 && isfield(ev.trlinfo,'med')
        obs = ev.trlinfo.med;
    elseif isfield(ev.trlinfo,'mean')
        obs = ev.trlinfo.mean;
    elseif isfield(ev.trlinfo,'med')
        obs = ev.trlinfo.med;
    end
end
obs = max(1, round(obs));
end

function n = local_goodn(ev)
%Retained participant count for one group x event, for the RC-36 column.
%
%goodn is set on every route psyrat_relsummary can take (subject-level and
%difference paths included), and both criterion paths read a relsummary that
%function built -- either the local_build_relsummary rebuild or an already-built
%one carried in on psyrat_data. So the isfield guard should never fire.
%
%It is here anyway because this file guards defensively everywhere else
%(local_obscount on trlinfo, local_muvalue on mu), and because the honest
%fallback for a missing count is NaN, not a plausible-looking number. A NaN in
%one cell tells the user the count is unavailable; a fabricated N would be read
%as a real retained sample. Do NOT "fix" this to 0 or to the group size.
n = NaN;
if isfield(ev,'goodn') && ~isempty(ev.goodn)
    n = double(ev.goodn);
end
end

function muval = local_muvalue(d)
%Extract a scalar average score (mu) for display in summary tables.
muval = NaN;
if isfield(d,'mu') && isfield(d.mu,'raw') && ~isempty(d.mu.raw)
    muval = mean(d.mu.raw(:));
elseif isfield(d,'sd_id') && isfield(d,'b_sigma') && ~isempty(d.sd_id.raw)
    muval = mean(cell2mat(d.sd_id.raw));
end
if isnan(muval)
    muval = 0;
end
end

function [ngroups, nevents, gnames, enames, modelcase] = local_groupinfo(psyrat_data)
if strcmpi(psyrat_data.rel.groups,'none')
    ngroups = 1;
    gnames = {''};
else
    ngroups = length(psyrat_data.rel.groups);
    gnames = psyrat_data.rel.groups(:);
end

if strcmpi(psyrat_data.rel.events,'none')
    nevents = 1;
    enames = {'none'};
else
    nevents = length(psyrat_data.rel.events);
    enames = psyrat_data.rel.events(:);
end

if ngroups == 1 && nevents == 1
    modelcase = 1;
elseif ngroups > 1 && nevents == 1
    modelcase = 2;
elseif ngroups == 1 && nevents > 1
    modelcase = 3;
else
    modelcase = 4;
end
end

function rel_name = local_relname(~)
%criterion outputs always use dependability (absolute decisions)
rel_name = 'Dependability';
end

function note = local_occnote(analysis)
%Occasion-estimand statement for the criterion screens.
%
%These surfaces deliberately pin n'_o = 1: psyrat_cutscore_trt is called with
%'nocc',1 on both the plot and the table path, and the relsummary those paths
%rebuild is built without noccmode/nocc, so local_resolve_effnocc returns 1.
%
%This is NOT the same kind of pin as gcoeff = 1, and conflating the two is the
%mistake this comment exists to prevent. A relative-error cut-score is not a
%valid quantity at all, which is why psyrat_cutscore_trt raises on gcoeff = 2.
%A multi-occasion cut-score IS valid. So n'_o = 1 is an estimand CHOICE, not a
%validity constraint: a criterion-referenced classification is a decision made
%on a single administration, so these screens report the dependability of that
%decision rather than of the composite the user selected for the norm-referenced
%coefficient (RC-27/RC-29, owner ruling 2026-07-30).
%
%Because it is a choice rather than a constraint it has to be visible. Under
%noccmode = 2 the parent view screen reports a k-occasion composite while these
%screens report a single administration, and without this note nothing on screen
%says so. The saved file already carries the same statement in its header.
%
%One-facet ('sing') data has no occasion facet, so it contributes no note and
%those titles stay byte-identical.
note = '';
if strcmp(analysis,'trt')
    note = 'single-occasion, n''_o = 1';
end
end

function lines = local_divergencenote(analysis, settings, L)
%RC-33. Extra error-dialog line naming the one reason these screens can report
%a different retained sample than the view screen the user came from.
%
%local_occnote says n'_o = 1 is pinned here by ruling. What that ruling does not
%cover is that the rebuild in local_build_relsummary re-runs the WHOLE D-study
%at that pin, so the trial cutoffs and the eventgoodids/eventbadids partition
%are recomputed -- and trlinfo from that partition is what local_obscount feeds
%back as obs. Under noccmode = 2 the parent screen selected a k-occasion
%composite and can therefore have retained a different set of participants.
%
%Gated exactly like local_occnote, and for the same reason. Under noccmode = 1
%the parent is at n'_o = 1 as well, so there is nothing to diverge and the line
%would be misleading; one-facet data has no occasion facet at all. The parent's
%k is deliberately not interpolated -- it is on the screen the user just left,
%and a fixed string is one fewer branch to get wrong.
lines = {};
if strcmp(analysis,'trt') && settings.noccmode == 2
    lines{end+1} = ['These screens compute at n''_o = 1, so the ' L.unit ...
        ' cutoffs and the'];
    lines{end+1} = 'retained participants can differ from the previous screen';
end
end

function local_save_table(varargin)
psyrat_data = varargin{3};
tbl = varargin{4};
cutoff = varargin{5};

%XLSX now writes cross-platform via writecell/writetable, so both formats are
%offered on every platform (including macOS).
[savename, savepath] = uiputfile(...
    {'*.xlsx','Excel File (.xlsx)';'*.csv','Comma-Separated Value File (.csv)'},...
    'Where would you like to save table?');

if isequal(savename,0) || isequal(savepath,0)
    return;
end

[~,~,ext] = fileparts(fullfile(savepath,savename));

%Build the header ONCE and use it for both formats. It used to live inside the
%.xlsx branch, which left the .csv branch writing a bare table with no header at
%all -- so a saved criterion table recorded no seed, engine, priors or
%convergence verdict and could not be interpreted on its own (RC-28).
filehead = {'Table Generated on'; datestr(clock);''};
filehead{end+1} = sprintf('PsyRAT Toolbox v%s',psyrat_data.ver);
filehead{end+1} = sprintf('Dataset: %s',psyrat_data.rel.filename);
filehead{end+1} = sprintf('Criterion Cutoff: %0.4f', cutoff);
%run provenance (seed/engine/priors/convergence); fail-open
filehead = [filehead; psyrat_provenance_lines(psyrat_data)];
%Scope of THIS export, stated after the provenance block so it qualifies the
%"Coefficient:" line that block emits on the two-facet route. That line reports
%the user's view-time selections, but the criterion table does not honor all of
%them: local_settings pins gcoeff = 1 (:123) and both cut-score kernels reject
%gcoeff = 2 (psyrat_cutscore_sing/-_trt), and local_tablecriterion hardcodes
%'nocc',1 (:376). The gap that makes this line necessary is the difference-score
%clause: on the trt_diff route psyrat_trt_coeflabel appends "difference score:
%Generalizability (relative error)" whenever the user's diffgcoeff differs from
%gcoeff, which it does by construction here, and no number in this file is a
%relative-error coefficient.
filehead{end+1} = ['Scope: criterion cut-scores are always absolute-error '...
    '(dependability) at n''_o = 1, regardless of any coefficient selection '...
    'reported above.'];
filehead{end+1} = '';

if strcmpi(ext,'.xlsx')
    writecell(filehead,fullfile(savepath,savename));
    writetable(tbl,fullfile(savepath,savename),...
        'Range',strcat('A',num2str(length(filehead)+1)));
else
    if ~strcmpi(ext,'.csv')
        savename = [savename '.csv'];
    end
    %CSV: write the header lines, then append the table (its own column names
    %plus data) so every column is handled without a hand-rolled formatspec --
    %which is also why adding the RC-36 n_Included column needed no change on
    %either export branch. 'WriteMode','append' defaults WriteVariableNames to
    %false, so force it on to keep the column-name row, matching the .xlsx
    %branch. Same idiom as psyrat_startview_splits/-_dynrel.
    fid = fopen(fullfile(savepath,savename),'w');
    %quote per the CSV convention (B19): a comma-bearing header line must
    %stay ONE spreadsheet field; comma-free lines stay byte-identical
    for i = 1:numel(filehead)
        fprintf(fid,'%s\n',psyrat_quote_csv_header(filehead{i}));
    end
    fclose(fid);
    writetable(tbl,fullfile(savepath,savename),...
        'WriteMode','append','WriteVariableNames',true);
end

end
