function diffplot = psyrat_trt_diffvtrialsplot(varargin)
%Plot the two-facet DIFFERENCE-SCORE reliability as a function of the number of
%trials retained per condition (one line per group)
%
%diffplot = psyrat_trt_diffvtrialsplot('psyrat_data',psyrat_data,...
%   'trials',[1 50],'relline',plotrelline,'relcutoff',relcutoff);
%
%
%This is the difference-score companion to psyrat_depvtrialsplot's per-condition
%curve, and the two are NOT interchangeable:
%
% - psyrat_depvtrialsplot draws the reliability of EACH CONDITION separately,
%   calling the plain two-facet kernel psyrat_rel_trt on the per-condition
%   components (the diagonals of the cross-condition covariances). It is
%   stratified by condition, so it draws one subplot per event, and it follows
%   the per-condition error type, relsummary.gcoeff.
% - this function draws the reliability of the DIFFERENCE between the two
%   conditions, calling psyrat_diffrel_trt with the full cross-condition
%   covariance blocks. A difference has no per-condition stratification, so it
%   is a single axes with one line per group, and it follows the separately
%   selected difference-score error type, relsummary.diffgcoeff.
%
%Both are drawn together by psyrat_relfigures for a 'trt_diff' analysis, so each
%of the viewer's two coefficient popups drives the series it names. Before RC-31
%only the per-condition curve existed and it was fed diffgcoeff, so the
%difference-score selector silently moved a per-condition figure.
%
%ESTIMAND: the curve is drawn at EQUAL trial counts in both conditions,
%obs = [n n]. psyrat_diffrel_trt takes a scalar pair [nX nY] and has no range
%mode, so a curve has to choose a path through that plane; equal counts is the
%only one on which the x-axis remains a plain trial count. This matches the
%existing precedent in psyrat_rel_diffdynrel_trt (:178), which projects its
%difference surface the same way. The CONSEQUENCE, stated on the figure because
%it is not otherwise visible: the headline difference score in the tables beside
%this plot is computed at the OBSERVED per-condition counts
%(psyrat_relsummary.m, 'obs',[trls1 trls2]), so when those two counts differ the
%headline value does not lie exactly on this curve.
%
%
%Inputs
% psyrat_data - PsyRAT Toolbox data structure array. Variance components should
%  be included, and psyrat_relsummary must already have been run (this function
%  reads relsummary.group(g).diffcomp).
% trials - range of trials to plot: [min max]
% relcutoff - reliability threshold for deeming data reliable (used for
%  reference line in plot)
%
%Optional Input
% relline - reliability estimate to plot (default: 2)
%  1 - lower limit
%  2 - point estimate
%  3 - upper limit
% diffgcoeff - difference-score coefficient to use. Defaults to
%  relsummary.diffgcoeff, which is what the viewer's difference-score popup
%  writes.
%  1 - dependability/reliability (absolute)
%  2 - generalizability (relative)
% CI - confidence interval width. Decimal from 0 to 1. (default: .95)
%
%
%Outputs
% diffplot - figure handle for the plot
% figure that displays the relationship between the difference-score
%  reliability and the number of trials retained per condition, stratified by
%  group

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
        'See help psyrat_trt_diffvtrialsplot for more information about inputs'));
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
            'See help psyrat_trt_diffvtrialsplot for more information \n'));
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
            'See help psyrat_trt_diffvtrialsplot for more information\n'));
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
            'See help psyrat_trt_diffvtrialsplot for more information\n'));
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
                'See help psyrat_trt_diffvtrialsplot for more information \n'));
        end
    else
        relline = 2; %default is point estimate
    end

    %check if diffgcoeff was specified.
    %If it is not found, fall back to the value relsummary recorded (which is
    %what the viewer's difference-score popup writes), then to dependability.
    ind = find(strcmpi('diffgcoeff',varargin),1);
    if ~isempty(ind)
        diffgcoeff = varargin{ind+1};
        if ~any(diffgcoeff == [1 2])
            error('varargin:diffgcoeff',... %Error code and associated error
                strcat('WARNING: diffgcoeff is invalid. Valid values are 1 or 2\n'));
        end
    elseif isfield(psyrat_data,'relsummary') && ...
            isfield(psyrat_data.relsummary,'diffgcoeff') && ...
            ~isempty(psyrat_data.relsummary.diffgcoeff)
        diffgcoeff = psyrat_data.relsummary.diffgcoeff;
    else
        diffgcoeff = 1;
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
                'See help psyrat_trt_diffvtrialsplot for more information \n'));
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

rs = psyrat_data.relsummary;

%check whether any groups exist
if strcmpi(psyrat_data.rel.groups,'none')
    ngroups = 1;
    gnames = {''};
else
    ngroups = length(psyrat_data.rel.groups);
    gnames = psyrat_data.rel.groups(:);
end

%condition names, used only to label which contrast is plotted. A two-facet
%difference is defined over exactly two events (psyrat_relsummary errors
%otherwise), so this is a label, not a stratification.
if strcmpi(psyrat_data.rel.events,'none')
    enames = {'',''};
else
    enames = cellstr(string(psyrat_data.rel.events(:)));
end

%The difference-score component blocks are published by psyrat_relsummary
%(relsummary.group(g).diffcomp) rather than re-derived here, so this curve and
%the headline difference score in the tables are built from identical
%components and can differ only in the trial counts projected. A struct without
%them was not produced by the 'trt_diff' branch, so say so plainly instead of
%failing on a missing field several frames down.
if ~isfield(rs,'group') || ~isfield(rs.group,'diffcomp')
    error('psyrat_trt_diffvtrialsplot:nodiffcomp', ...
        ['relsummary.group.diffcomp is not populated, so the difference-score ' ...
        'components are unavailable. How to fix: run psyrat_relsummary with ' ...
        '''analysis'',''trt_diff'' on these data first (the results viewer does ' ...
        'this), then call this plot.']);
end

%difference-score error type, and the estimate keyword psyrat_diffrel_trt takes
if diffgcoeff == 2
    diff_est = 'gen';
    ylbl = 'Generalizability';
else
    diff_est = 'dep';
    ylbl = 'Dependability';
end

%Occasion estimand (n'_o), read behind the same guard every other consumer of
%the test-retest coefficient uses (psyrat_depvtrialsplot,
%psyrat_trt_relvtrialsplot, psyrat_trt_reloverallt). Under the shipped default
%(noccmode = 1) effnocc is 1.
if isfield(rs,'nocc') && ~isempty(rs.nocc)
    nocc = rs.nocc;
else
    nocc = 1;
end

%create an x-axis for the number of observations
x = trials(1):trials(2);

%One row per requested trial count, indexed RELATIVE to the start of the range
%(row 1 is trials(1)), which is the RC-30 convention its two sibling curve
%plots now document and psyrat_rel_sing/psyrat_rel_trt already use for their
%return vectors. Indexing by the absolute trial number instead would grow the
%buffer to trials(2) rows while x has only numel(x), which is exactly the
%MATLAB:samelen defect RC-30 fixed in the siblings.
%
%Fully preallocated, unlike the siblings' zeros(numel(x),0): they receive a
%whole column vector back from one range-mode calc call, while the kernel here
%has no range mode and is called once per trial count, so the buffer is filled
%element by element and would otherwise be regrown on every iteration.
plotrel = zeros(numel(x),ngroups);

%create the figure
diffplot = figure;
diffplot.Tag = 'psyrat_output';
set(gcf,'NumberTitle','Off');
diffplot.Position = [125 630 900 450];
fsize = 16;

%psyrat_diffrel_trt emits one admissibility warning per CALL and this loop calls
%it numel(x) times per group, so an unguarded loop would print one block per
%trial count. Suppress, aggregate, and report once -- the same idiom, and the
%same warning identifier, as psyrat_rel_diffdynrel_trt.m:152-202, which loops the
%identical kernel over a grid, and the one psyrat_admissibility_warn's own help
%prescribes for loop callers.
%
%The guard matters MORE here than on an unequal-count caller, not less. At
%obs = [n n] the harmonic mean is n, so every D-study-adjusted block reduces
%exactly to raw/n (and likewise raw/k for the occasion blocks), and the
%projection ALONE therefore cannot produce a negative on this curve.
%
%What does NOT follow -- and what this comment asserted until RC-35 -- is that
%admissibility "does not vary with n," so that the guard either never fires or
%fires at every single point. That is false. The RECORDED quantities are SUMS of
%blocks carrying THREE different divisors (n, the occasion count k, and their
%product n*k), so no single factor cancels and the sign can still flip along the
%curve: at bpi_raw = +18, bpo_raw = -1, k = 2, reltype 3 the relative error
%variance is positive for n <= 20 and negative for n >= 50, from one fixed
%non-PSD block (RC-34). A curve can fire on PART of its range.
%
%The guard and the aggregation below are unchanged by that correction -- an
%unguarded loop still emits up to numel(x) near-identical blocks, which is
%exactly what aggregation exists for. What the correction forbids is
%"simplifying" this to a single evaluation on the reasoning that one point
%stands for all of them.
%
%The reduction to raw/n also means a PSD input cannot trip this on this path: a
%raw variance of a difference is non-negative whenever the cross-condition
%block is a real covariance matrix. What reaches it is a block that is NOT positive
%semi-definite, which the gamma/scaled-chi-square family can produce because its
%diagonals are clamped while its cross-covariances are not (RC-23). Suppressing
%without aggregating would drop exactly that signal.
adm_ws = warning('off','psyrat:negerrvar');
adm_all = [];

for gloc = 1:ngroups

    %the eight blocks exactly as psyrat_relsummary assembled them for the
    %headline difference score
    dc = rs.group(gloc).diffcomp;
    blockargs = {'bp',dc.bp,'bpi',dc.bpi,'bpo',dc.bpo,...
        'bt',dc.bt,'bo',dc.bo,'boi',dc.boi,...
        'er_var',dc.er_var,'er_cov',dc.er_cov};

    for ii = 1:numel(x)
        n = x(ii);

        %obs = [n n]: equal trial counts in both conditions (see the ESTIMAND
        %note in the header). nocc is applied to both conditions for the same
        %reason.
        ds = psyrat_diffrel_trt(blockargs{:},...
            'obs',[n n],'nocc',[nocc nocc],...
            'reltype',rs.reltype,'est',diff_est,'CI',ciperc);

        switch relline
            case 1 %lower limit
                plotrel(ii,gloc) = ds.ll;
            case 2 %point estimate
                plotrel(ii,gloc) = ds.pt;
            case 3 %upper limit
                plotrel(ii,gloc) = ds.ul;
        end

        if isfield(ds,'admissibility')
            adm_all = [adm_all, ds.admissibility]; %#ok<AGROW>
        end
    end
end

warning(adm_ws);
psyrat_admissibility_warn(adm_all);

switch relline
    case 1
        plottitle = 'Lower Limit of 95% Credible Interval';
    case 2
        plottitle = 'Point Estimate';
    case 3
        plottitle = 'Upper Limit of 95% Credible Interval';
end

switch rs.reltype
    case 1
        pref1 = 'Coefficient of Equivalence, ';
    case 2
        pref1 = 'Coefficient of Stability, ';
    case 3
        pref1 = 'Coefficient of Equivalence and Stability, ';
    otherwise
        pref1 = '';
end

%Name the occasion estimand in the window title for the same reason the plain
%test-retest curve does (RC-10): the axis label says only "Dependability", so
%without it the reader cannot tell which n'_o the curve was drawn at.
set(gcf,'Name',['Difference Score, ' pref1 ylbl ' v ' L.NumberOf ': ' ...
    plottitle ' (' psyrat_trt_occlabel(psyrat_data) ')']);

plot(x,plotrel);
axis([0 trials(2) 0 1]);
set(gca,'fontsize',16);

%The axes title carries the two things a reader cannot recover from the axes:
%WHICH contrast is plotted, and that the curve is drawn at equal trial counts
%per condition, so the headline difference score computed at the observed
%(possibly unequal) counts need not lie on it.
if isempty(enames{1}) || isempty(enames{2})
    contraststr = 'Difference Score';
else
    contraststr = ['Difference Score: ' enames{1} ' - ' enames{2}];
end
%Interpreter none here and on the legend below: the contrast string carries
%both event names, and an underscored name (inc_cor) must print literally,
%where the default TeX interpreter would subscript the character after it.
title({contraststr, ...
    ['(equal ' L.units ' per condition)']},'FontSize',16,'Interpreter','none');

ylabel(ylbl,'FontSize',fsize);
xlabel('Number of Observations','FontSize',fsize);
psyrat_addhline(relcutoff,'Color','b','LineStyle',':');

%only build a named legend when there are multiple groups to label; otherwise
%legend('boxoff') would auto-create MATLAB's default "data1" legend for a single
%unnamed series. Same guard as the two sibling curve plots, minus their event
%condition, which does not apply to a per-group-only figure.
if ngroups > 1
    leg = legend(gnames{:},'Location','southeast');
    set(leg,'FontSize',fsize,'Interpreter','none');
    legend('boxoff');
end

end
