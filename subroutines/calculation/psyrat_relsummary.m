function [psyrat_data,relerr] = psyrat_relsummary(varargin)
%Summarize reliabilty of psyrat_data
%
%psyrat_relsummary('psyrat_data',psyrat_data,'analysis','sing',...
%   'depcutoff',depcutoff,'meascutoff',meascutoff,...
%   'depcentmeas',depcentmeas)
%
%
%Inputs
% psyrat_data - PsyRAT Toolbox data structure array. Variance components should
%  be included
% analysis - dispatch key selecting the summary path. Accepted values:
%   'sing' single-occasion; 'sing_sserr' single-occasion subject-level error;
%   'sing_diff' single-occasion difference score; 'sing_diff_sserr'
%   single-occasion subject-level difference; 'trt' test-retest (multi-occasion);
%   'trt_sserr' test-retest subject-level error; 'trt_diff' test-retest
%   difference score. The splits, dynamic-reliability, and
%   difference-of-differences families do not use this function at all:
%   psyrat_report reads their results directly (manual chapter 14). A key
%   outside the accepted set, including a stored rel.analysis value such as
%   'ic', is rejected with error id varargin:analysis.
%
% EITHER (for analysis type 'sing', 'sing_sserr', 'sing_diff')
%  depcutoff - dependability threshold for considering data reliable
%  meascutoff - which estimate of dependability to use for cutoff (1 - lower
%   limit of credible interval, 2 - point estimate, 3 - upper limit of
%   credible interval)
%  depcentmeas - which measure of central tendency to use for estimating
%   overall score dependability after applying trial cutoffs (not used in
%   'sing_sserr')
%  diffgcoeff - (optional for 'sing_diff') coefficient for difference-score
%   reliability: 1 - dependability (absolute), 2 - generalizability
%   (relative). Default is 1.
%
% OR (for analysis type 'trt')
%  gcoeff - g-theory coefficient to calculate 1 - dependability, 2 -
%   generalizability
%  reltype - reliability coefficient to calculate 1 - equivalence, 2 -
%   stability, 3 - equivalence and stability
%  noccmode - (optional) test-retest occasion estimand. 1 - single-occasion
%   score (n'_o = 1), the default: reliability of a score from n'_i trials at
%   one occasion, generalizing over the occasion facet. 2 - multi-occasion
%   composite score (n'_o = k): reliability of a score averaged across k
%   occasions. The same n'_o is applied to the coefficient and its SEM.
%  nocc - (optional) composite size k for noccmode == 2. Positive integer.
%   Defaults to the number of observed occasions (numel(unique(time))). Use an
%   explicit k for a D-study projection. Ignored for the single-occasion score.
%  relcutoff - reliability threshold for considering data reliable
%  meascutoff - which estimate of dependability to use for cutoff (1 - lower
%   limit of credible interval, 2 - point estimate, 3 - upper limit of
%   credible interval)
%  relcentmeas - which measure of central tendency to use for estimating
%   overall score reliability after applying trial cutoffs
%
%Optional Input
% CI - credible interval width. Decimal from 0 to 1. (default: .95)
% interactive - true/false. Whether the threshold-failure notice may be
%  raised as a modal dialog. Default: psyrat_data.proc.interactive when the
%  struct carries that field (psyrat_run writes false for headless runs),
%  otherwise true, which is the GUI's behavior. A .psyrat saved by the GUI
%  carries no such field, so pass false from a script that summarizes
%  GUI-produced files; a batch loop then never blocks on a dialog, and the
%  same condition is still reported through relerr.nogooddata.
%
%
%Outputs
% psyrat_data - rel field will be added to psyrat_data that includes a summary of
%  reliability information for the data
% relerr - structure array with information about reliabitliy analyses.
%  Warnings will be provided to the user if they apply to the analyzed
%  dataset
%

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

%CONVENTION -- THE GOOD/BAD ID CARRIER IS ALWAYS A CELL ARRAY (findings B29
%and B34, both ruled 2026-08-19).
%
%relsummary.group(gloc).goodids and its per-event twin
%...event(eloc).eventgoodids carry THREE states, all of them cell arrays:
%  1. an EMPTY CELL      -- nothing was retained. This covers BOTH "no trial
%                           cutoff could be found" and "a cutoff was found and
%                           it excluded everyone" (RC-38); the carrier does not
%                           distinguish them and does not need to (see below).
%  2. a non-empty cellstr of retained ids
%  3. state 2 where one of those ids is literally the string 'none'
%
%Arm 1 (treat as a real retained set) is correct for states 2 and 3 only. Every
%such test is spelled
%
%    if iscell(<carrier>) && ~isempty(<carrier>)
%
%B34 REMOVED THE CHAR SENTINEL. This carrier used to hold the CHAR 'none' for
%"nothing retained", which made a fourth state and leaked out of this file. The
%consumers that splice the carrier into the Good IDs export
%(psyrat_depoverallt, psyrat_trt_reloverallt) do no type test at all, so the
%char produced THREE different user-visible failures depending on how many
%groups were affected: a silently exported phantom participant named "none"
%when only SOME groups had nothing retained (MATLAB wraps the 4x1 char as one
%cell element), four junk one-letter rows per group when ALL of them did, and
%an outright crash in the CSV branch (gids{i} on a char). It also crashed this
%file, because table('none') errors -- MATLAB reads a char row vector as
%variable NAMES. Producers now assign cell(0,1), so none of that can arise.
%
%WHICH REASON nothing was retained is NOT carried by this field, and must not
%be recovered from it. "No cutoff was found" is recorded by relerr.trlcutoff
%and by the per-cell event(eloc).trlcutoff = -1; "the cutoff extrapolated
%beyond the data" is recorded by relerr.trlmax. Those flags, not the carrier,
%are what raise the user-facing dialogs, and they remain distinct.
%
%B29 is the companion ruling and is why the predicate tests TYPE AND EMPTINESS
%rather than a value: the arm tests were once spelled
%"if ~strcmp(<carrier>,'none')", and with a participant named 'none' among the
%retained ids strcmp returns a logical ARRAY, MATLAB's "if" requires every
%element true, so the test was false and the bad-ids arm was taken -- a silent
%wrong result dressed as a legitimate rejection. Pinned by
%tests/TestRelsummaryNoneIdSentinel.m, which holds the whole-file count at
%23 goodids + 0 eventgoodids sites (re-derived for G39; the census comment in
%that test walks the arithmetic from B29's original 11 + 11).
%
%INTERSECTION SEMANTICS (G39 ruling, 2026-08-29): group-level inclusion is
%STRICT -- a participant is retained for a group only by meeting the cutoff in
%EVERY event, and an event that retained nobody empties the group's retained
%set rather than dropping out of the intersection. The pre-G39 skip made the
%retained sample depend on which events happened to die (retention could
%INCREASE at a stricter threshold), and it disagreed with the trt_diff branch,
%which always intersected strictly. Under strict semantics the per-event rows
%stay informative (a dead event's row carries the -1 sentinels), goodn is the
%group-level count on every row, and n Included + n Excluded complement -- over
%the participants present in every event's data. An id with zero rows in some
%event lands in neither carrier only when that event is the first: the fold
%is seeded from event 1's roster, so an id present in event 1 but absent from a
%later event is folded into badids, while an id absent from event 1 never enters
%the fold at all (audit 2026-09-05; the asymmetry is filed, not fixed here). One
%carve-out to the "group-level goodn" clause (adversarial wave, 2026-08-29):
%the sing_diff_sserr helper's per-event fields carry an EVENT-level goodn
%(the sum of that event's inclusion indicator); no shipped consumer reads
%event.goodn on that route. An earlier revision of this paragraph leaned on
%"psyrat_depoverallt's sserr branch derives its printed counts from the
%stored tables" as the reason the rendered rows were safe -- that assurance
%was FALSE on ic_diff_sserrvar (G53): the branch restricted the stored event
%tables to the GROUP carrier, and on the difference route the group carrier
%derives from the DIFFERENCE-score inclusion while each event table is
%marked independently against the same cutoff, with no enforced ordering
%between the two quantities, so a diff-included but event-excluded subject
%fell out of both printed counts and the row summed short of that event's
%table. G53 fixed this in the RENDERER: psyrat_depoverallt's sserr branch
%now prints the group-carrier count, so its rendered rows complement as an
%ismember identity (the good and bad row masks are the same test negated)
%-- the same group-carrier semantics G39/G45 gave the STORED field on the
%intersection routes, applied at render time; the stored event.goodn on
%this route deliberately keeps the event-level count (the carve-out
%above). The rendered complement target is that event's TABLE roster,
%which equals the participant roster only when every subject has rows in
%both events. The committed fixtures and the manual data are
%complete-by-design, so that distinction is theoretical there, but do not
%lean on the complement for a file with missing cells.
%
%A bare "~isempty(x)" is also correct now that the char is gone, and the two
%places that already used it were left as they were rather than churned.
%Prefer the "iscell(x) && ~isempty(x)" spelling for any NEW arm test: it
%states the invariant rather than relying on it.
%
%DO NOT cite line numbers for the arm-test sites in this comment. Editing this
%file has already invalidated its own citations twice -- re-derive them with
%grep on the predicate above.

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
            'See help psyrat_relsummary for more information \n'));
    end
    
    %check if the dependability cutoff was specified
    ind = find(strcmp('analysis',varargin),1);
    if ~isempty(ind)
        relanalysis = varargin{ind+1};
    else
        error('varargin:noanalysistype',... %Error code and associated error
            strcat('WARNING: analysis not specified \n\n',...
            'Please input the analysis to be run: ''sing'' or ''trt'' \n'));
    end
    %the key is a dispatch key, not the stored rel.analysis string. Neither
    %switch below has an otherwise arm, so an unknown key used to fall through
    %both, leave relerr unassigned and fail with MATLAB's generic
    %output-argument error, or return silently on a one-output call (audit
    %2026-09-05). Reject it here, naming the mapping.
    acceptedkeys = {'sing','sing_sserr','sing_diff','sing_diff_sserr',...
        'trt','trt_sserr','trt_diff'};
    if ~((ischar(relanalysis) || (isstring(relanalysis) && isscalar(relanalysis))) && ...
            any(strcmp(char(relanalysis),acceptedkeys)))
        if ischar(relanalysis) || isstring(relanalysis)
            keytext = char(relanalysis);
        else
            keytext = class(relanalysis);
        end
        error('varargin:analysis',... %Error code and associated error
            strcat('WARNING: unrecognized analysis key ''',keytext,''' \n\n',...
            'analysis is a dispatch key, not the stored rel.analysis value:\n',...
            '  rel.analysis ic / ic_sserrvar / ic_diff / ic_diff_sserrvar ->\n',...
            '    ''sing'' / ''sing_sserr'' / ''sing_diff'' / ''sing_diff_sserr''\n',...
            '  rel.analysis trt / trt_sserrvar / trt_diff ->\n',...
            '    ''trt'' / ''trt_sserr'' / ''trt_diff''\n',...
            'Splits, dynamic-reliability and difference-of-differences results\n',...
            'do not use psyrat_relsummary; pass them to psyrat_report directly.\n',...
            'See help psyrat_relsummary for more information \n'));
    end
    relanalysis = char(relanalysis);
end

%Determine whether this run may show interactive GUI dialogs. Resolution
%order (G3, 2026-09-05): an explicit 'interactive' name-value argument wins;
%otherwise the stored flag psyrat_data.proc.interactive (the processing
%wrapper psyrat_run writes it false for headless/batch runs); otherwise true,
%so the GUI's existing behavior is unchanged. The explicit option exists
%because a .psyrat saved by the GUI carries no proc.interactive field: before
%it, a scripted loop over GUI-produced files resolved to interactive and could
%stall on the threshold-failure modal below, once per dataset, with no way to
%say otherwise. Headless callers must never block on a modal. The
%threshold-failure notices below are gated on this flag, but the control flow
%(relerr.nogooddata = 1; return) stays unconditional, so every caller still
%detects the failure whether or not a dialog was shown.
interactive = [];
ind = find(strcmpi('interactive',varargin),1);
if ~isempty(ind)
    optval = varargin{ind+1};
    %a logical or numeric scalar; ~isreal, NaN and non-scalars are rejected
    %the way the sibling option guards in this file reject them
    if ~(isscalar(optval) && (islogical(optval) || ...
            (isnumeric(optval) && isreal(optval) && ~isnan(optval))))
        error('varargin:interactive',... %Error code and associated error
            strcat('WARNING: interactive should be a logical scalar \n\n',...
            'Pass true to allow the threshold-failure dialog, or false to\n',...
            'suppress it on headless/batch runs.\n',...
            'See help psyrat_relsummary for more information \n'));
    end
    interactive = logical(optval);
end
if isempty(interactive)
    if isfield(psyrat_data,'proc') && isfield(psyrat_data.proc,'interactive') && ...
            ~isempty(psyrat_data.proc.interactive)
        interactive = logical(psyrat_data.proc.interactive);
    else
        interactive = true;
    end
end

switch relanalysis
    case {'sing','sing_sserr','sing_diff','sing_diff_sserr'}
        %check if depcutoff was specified.
        %If it is not found, set display error.
        ind = find(strcmpi('depcutoff',varargin),1);
        if ~isempty(ind)
            depcutoff = varargin{ind+1};
        else
            error('varargin:depvalue',... %Error code and associated error
                strcat('WARNING: depvalue not specified \n\n',...
                'Please input the dependability threshold for trial cutoffs.\n',...
                'See help psyrat_relsummary for more information \n'));
        end
        
        %check if meascutoff was specified.
        %If it is not found, set display error.
        ind = find(strcmpi('meascutoff',varargin),1);
        if ~isempty(ind)
            meascutoff = varargin{ind+1};
            if length(meascutoff) ~= 1 || ~any(meascutoff == [1 2 3])
                %G58: an out-of-range value fell through the consuming
                %switches silently and died far downstream as an undefined
                %variable; fail loudly at the parse block instead. The
                %length test runs first: a vector can slip past the
                %elementwise any() membership test, and an empty errors
                %inside it (MEASURED R2025b: [] == [1 2 3] throws
                %sizeDimensionsMustMatch), so both short-circuit into this
                %named error. Scalar idiom per local_parse_nocc.
                error('varargin:meascutoff',... %Error code and associated error
                    strcat('WARNING: meascutoff should be 1 (lower limit),',...
                    ' 2 (point estimate), or 3 (upper limit)\n',...
                    'See help psyrat_relsummary for more information \n'));
            end
        else
            error('varargin:meascutoff',... %Error code and associated error
                strcat('WARNING: Which dependability estimate not specified \n\n',...
                'Please enter a valid input for meascutoff\n',...
                '1 - lower limit of credible interval\n',...
                '2 - point estimate of credible interval\n',...
                '3 - upper limit of credible interval\n',...
                'See help psyrat_relsummary for more information \n'));
        end
        
        %check if depcentmeas was specified.
        %If it is not found, set display error.
        ind = find(strcmpi('depcentmeas',varargin),1);
        if ~isempty(ind)
            depcentmeas = varargin{ind+1};
            if length(depcentmeas) ~= 1 || ~any(depcentmeas == [1 2])
                %G59: same class as the G58 meascutoff guard above (the
                %length test carries the same rationale).
                error('varargin:depcentmeas',... %Error code and associated error
                    strcat('WARNING: depcentmeas should be 1 (mean) or',...
                    ' 2 (median)\n',...
                    'See help psyrat_relsummary for more information \n'));
            end
        else
            error('varargin:depcentmeas',... %Error code and associated error
                strcat('WARNING: Which central tendency to use for ',...
                'estimating overall score dependability\n',...
                'Please ented a valid input for depcentmeas\n',...
                '1 - mean\n',...
                '2 - median\n',...
                'See help psyrat_relsummary for more information \n'));
        end
        
        %check if gcoeff was specified.
        %If it is not found, set as default: dependability.
        ind = find(strcmpi('gcoeff',varargin),1);
        if ~isempty(ind)
            gcoeff = varargin{ind+1};
            if ~any(gcoeff == [1 2])
                error('varargin:gcoeff',... %Error code and associated error
                    strcat('WARNING: gcoeff not specified \n\n',...
                    'Please input gcoeff as:\n',...
                    '1 - dependability (absolute)\n',...
                    '2 - generalizability (relative)\n',...
                    'See help psyrat_relsummary for more information \n'));
            end
        else
            gcoeff = 1;
        end

        %check if diffgcoeff was specified (single-session difference-score
        %analysis only). If missing, default to dependability.
        ind = find(strcmpi('diffgcoeff',varargin),1);
        if ~isempty(ind)
            diffgcoeff = varargin{ind+1};
            if ~any(diffgcoeff == [1 2])
                error('varargin:diffgcoeff',... %Error code and associated error
                    strcat('WARNING: diffgcoeff not specified \n\n',...
                    'Please input diffgcoeff as:\n',...
                    '1 - dependability (absolute)\n',...
                    '2 - generalizability (relative)\n',...
                    'See help psyrat_relsummary for more information \n'));
            end
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
                    'See help psyrat_depvtrialsplot for more information \n'));
            end
        else
            ciperc = .95; %default: 95%
        end
    case 'trt_sserr'
        %subject-level (per-participant) test-retest reliability (trt_sserrvar,
        %analysis 25). A hybrid of the subject-level (sserr) participant-inclusion
        %contract (depcutoff/meascutoff) and the test-retest coefficient selection
        %(gcoeff/reltype/nocc). See the compute branch below.

        %check if gcoeff was specified. If not, default to dependability.
        ind = find(strcmpi('gcoeff',varargin),1);
        if ~isempty(ind)
            gcoeff = varargin{ind+1};
            if ~any(gcoeff == [1 2])
                error('varargin:gcoeff',... %Error code and associated error
                    strcat('WARNING: gcoeff is invalid \n\n',...
                    'Please input gcoeff as:\n',...
                    '1 - dependability (absolute)\n',...
                    '2 - generalizability (relative)\n',...
                    'See help psyrat_relsummary for more information \n'));
            end
        else
            gcoeff = 1;
        end

        %check if reltype was specified (required: test-retest coefficient type).
        ind = find(strcmpi('reltype',varargin),1);
        if ~isempty(ind)
            reltype = varargin{ind+1};
        else
            error('varargin:reltype',... %Error code and associated error
                strcat('WARNING: reltype not specified \n\n',...
                'Please input the reliability coefficient to calculate.\n',...
                'See help psyrat_relsummary for more information \n'));
        end

        %parse the test-retest occasion estimand (noccmode + nocc), shared with
        %the 'trt'/'trt_diff' branches so the estimand cannot drift.
        [noccmode, noccarg] = local_parse_nocc(varargin{:});

        %check if depcutoff was specified.
        ind = find(strcmpi('depcutoff',varargin),1);
        if ~isempty(ind)
            depcutoff = varargin{ind+1};
        else
            error('varargin:depvalue',... %Error code and associated error
                strcat('WARNING: depvalue not specified \n\n',...
                'Please input the reliability threshold for participant ',...
                'inclusion.\n',...
                'See help psyrat_relsummary for more information \n'));
        end

        %check if meascutoff was specified.
        ind = find(strcmpi('meascutoff',varargin),1);
        if ~isempty(ind)
            meascutoff = varargin{ind+1};
            if length(meascutoff) ~= 1 || ~any(meascutoff == [1 2 3])
                %G58: an out-of-range value fell through the consuming
                %switches silently and died far downstream as an undefined
                %variable; fail loudly at the parse block instead. The
                %length test runs first: a vector can slip past the
                %elementwise any() membership test, and an empty errors
                %inside it (MEASURED R2025b: [] == [1 2 3] throws
                %sizeDimensionsMustMatch), so both short-circuit into this
                %named error. Scalar idiom per local_parse_nocc.
                error('varargin:meascutoff',... %Error code and associated error
                    strcat('WARNING: meascutoff should be 1 (lower limit),',...
                    ' 2 (point estimate), or 3 (upper limit)\n',...
                    'See help psyrat_relsummary for more information \n'));
            end
        else
            error('varargin:meascutoff',... %Error code and associated error
                strcat('WARNING: Which estimate not specified \n\n',...
                'Please enter a valid input for meascutoff\n',...
                '1 - lower limit of credible interval\n',...
                '2 - point estimate of credible interval\n',...
                '3 - upper limit of credible interval\n',...
                'See help psyrat_relsummary for more information \n'));
        end

        %check if CI was specified. If not, default to 95%.
        ind = find(strcmpi('CI',varargin),1);
        if ~isempty(ind)
            ciperc = varargin{ind+1};
            if ciperc > 1 || ciperc < 0
                error('varargin:ci',... %Error code and associated error
                    strcat('WARNING: Size of credible interval should ',...
                    'be a value between 0 and 1\n',...
                    'A value of ',sprintf(' %2.2f',ciperc),...
                    ' is invalid\n',...
                    'See help psyrat_relsummary for more information \n'));
            end
        else
            ciperc = .95; %default: 95%
        end

    case 'trt'
        %check if gcoeff was specified.
        %If it is not found, set display error.
        ind = find(strcmpi('gcoeff',varargin),1);
        if ~isempty(ind)
            gcoeff = varargin{ind+1};
        else
            error('varargin:gcoeff',... %Error code and associated error
                strcat('WARNING: gcoeff not specified \n\n',...
                'Please input the g-theory coefficient to calculate.\n',...
                'See help psyrat_relsummary for more information \n'));
        end
        
        %check if reltype was specified.
        %If it is not found, set display error.
        ind = find(strcmpi('reltype',varargin),1);
        if ~isempty(ind)
            reltype = varargin{ind+1};
        else
            error('varargin:reltype',... %Error code and associated error
                strcat('WARNING: reltype not specified \n\n',...
                'Please input the reliability coefficient to calculate.\n',...
                'See help psyrat_relsummary for more information \n'));
        end

        %parse the test-retest occasion estimand (noccmode + nocc), shared with
        %the 'trt_diff' branch so the two cannot drift. noccmode 1 = single-
        %occasion (n'_o = 1); 2 = multi-occasion composite (n'_o = k, with
        %nocc = k; [] defaults k to the observed occasions, resolved once REL is
        %known). The same n'_o is threaded into the coefficient and its SEM (F2).
        [noccmode, noccarg] = local_parse_nocc(varargin{:});

        %check if relcutoff was specified.
        %If it is not found, set display error.
        ind = find(strcmpi('relcutoff',varargin),1);
        if ~isempty(ind)
            relcutoff = varargin{ind+1};
        else
            error('varargin:relcutoff',... %Error code and associated error
                strcat('WARNING: relcutoff not specified \n\n',...
                'Please input the reliability threshold for trial cutoffs.\n',...
                'See help psyrat_relsummary for more information \n'));
        end
        
        %check if meascutoff was specified.
        %If it is not found, set display error.
        ind = find(strcmpi('meascutoff',varargin),1);
        if ~isempty(ind)
            meascutoff = varargin{ind+1};
            if length(meascutoff) ~= 1 || ~any(meascutoff == [1 2 3])
                %G58: an out-of-range value fell through the consuming
                %switches silently and died far downstream as an undefined
                %variable; fail loudly at the parse block instead. The
                %length test runs first: a vector can slip past the
                %elementwise any() membership test, and an empty errors
                %inside it (MEASURED R2025b: [] == [1 2 3] throws
                %sizeDimensionsMustMatch), so both short-circuit into this
                %named error. Scalar idiom per local_parse_nocc.
                error('varargin:meascutoff',... %Error code and associated error
                    strcat('WARNING: meascutoff should be 1 (lower limit),',...
                    ' 2 (point estimate), or 3 (upper limit)\n',...
                    'See help psyrat_relsummary for more information \n'));
            end
        else
            error('varargin:meascutoff',... %Error code and associated error
                strcat('WARNING: Which estimate not specified \n\n',...
                'Please enter a valid input for meascutoff\n',...
                '1 - lower limit of credible interval\n',...
                '2 - point estimate of credible interval\n',...
                '3 - upper limit of credible interval\n',...
                'See help psyrat_relsummary for more information \n'));
        end
        
        %check if relcentmeas was specified.
        %If it is not found, set display error.
        ind = find(strcmpi('relcentmeas',varargin),1);
        if ~isempty(ind)
            relcentmeas = varargin{ind+1};
            if length(relcentmeas) ~= 1 || ~any(relcentmeas == [1 2])
                %G59: same class as the G58 meascutoff guard above (the
                %length test carries the same rationale).
                error('varargin:relcentmeas',... %Error code and associated error
                    strcat('WARNING: relcentmeas should be 1 (mean) or',...
                    ' 2 (median)\n',...
                    'See help psyrat_relsummary for more information \n'));
            end
        else
            error('varargin:relcentmeas',... %Error code and associated error
                strcat('WARNING: Which central tendency to use for ',...
                'estimating overall score reliability\n',...
                'Please ented a valid input for relcentmeas\n',...
                '1 - mean\n',...
                '2 - median\n',...
                'See help psyrat_relsummary for more information \n'));
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

    case 'trt_diff'
        %Two-facet (test-retest) difference-score analysis. Gathers the
        %test-retest arguments (reltype, gcoeff, the occasion estimand, relcutoff,
        %relcentmeas) PLUS the difference-score coefficient choice (diffgcoeff),
        %so the per-condition reliability uses the trt machinery and the
        %difference score uses psyrat_diffrel_trt.

        ind = find(strcmpi('gcoeff',varargin),1);
        if ~isempty(ind)
            gcoeff = varargin{ind+1};
        else
            error('varargin:gcoeff',... %Error code and associated error
                strcat('WARNING: gcoeff not specified \n\n',...
                'Please input the g-theory coefficient to calculate.\n',...
                'See help psyrat_relsummary for more information \n'));
        end

        ind = find(strcmpi('diffgcoeff',varargin),1);
        if ~isempty(ind)
            diffgcoeff = varargin{ind+1};
            if ~any(diffgcoeff == [1 2])
                error('varargin:diffgcoeff',... %Error code and associated error
                    strcat('WARNING: diffgcoeff not specified \n\n',...
                    'Please input diffgcoeff as:\n',...
                    '1 - dependability (absolute)\n',...
                    '2 - generalizability (relative)\n',...
                    'See help psyrat_relsummary for more information \n'));
            end
        else
            diffgcoeff = 1;
        end

        ind = find(strcmpi('reltype',varargin),1);
        if ~isempty(ind)
            reltype = varargin{ind+1};
        else
            error('varargin:reltype',... %Error code and associated error
                strcat('WARNING: reltype not specified \n\n',...
                'Please input the reliability coefficient to calculate.\n',...
                'See help psyrat_relsummary for more information \n'));
        end

        %occasion estimand (n'_o) parse, shared with the 'trt' branch via
        %local_parse_nocc so the two branches accept identical inputs.
        [noccmode, noccarg] = local_parse_nocc(varargin{:});

        ind = find(strcmpi('relcutoff',varargin),1);
        if ~isempty(ind)
            relcutoff = varargin{ind+1};
        else
            error('varargin:relcutoff',... %Error code and associated error
                strcat('WARNING: relcutoff not specified \n\n',...
                'Please input the reliability threshold for trial cutoffs.\n',...
                'See help psyrat_relsummary for more information \n'));
        end

        ind = find(strcmpi('meascutoff',varargin),1);
        if ~isempty(ind)
            meascutoff = varargin{ind+1};
            if length(meascutoff) ~= 1 || ~any(meascutoff == [1 2 3])
                %G58: an out-of-range value fell through the consuming
                %switches silently and died far downstream as an undefined
                %variable; fail loudly at the parse block instead. The
                %length test runs first: a vector can slip past the
                %elementwise any() membership test, and an empty errors
                %inside it (MEASURED R2025b: [] == [1 2 3] throws
                %sizeDimensionsMustMatch), so both short-circuit into this
                %named error. Scalar idiom per local_parse_nocc.
                error('varargin:meascutoff',... %Error code and associated error
                    strcat('WARNING: meascutoff should be 1 (lower limit),',...
                    ' 2 (point estimate), or 3 (upper limit)\n',...
                    'See help psyrat_relsummary for more information \n'));
            end
        else
            error('varargin:meascutoff',... %Error code and associated error
                strcat('WARNING: Which estimate not specified \n\n',...
                'Please enter a valid input for meascutoff\n',...
                '1 - lower limit of credible interval\n',...
                '2 - point estimate of credible interval\n',...
                '3 - upper limit of credible interval\n',...
                'See help psyrat_relsummary for more information \n'));
        end

        ind = find(strcmpi('relcentmeas',varargin),1);
        if ~isempty(ind)
            relcentmeas = varargin{ind+1};
            if length(relcentmeas) ~= 1 || ~any(relcentmeas == [1 2])
                %G59: same class as the G58 meascutoff guard above (the
                %length test carries the same rationale).
                error('varargin:relcentmeas',... %Error code and associated error
                    strcat('WARNING: relcentmeas should be 1 (mean) or',...
                    ' 2 (median)\n',...
                    'See help psyrat_relsummary for more information \n'));
            end
        else
            error('varargin:relcentmeas',... %Error code and associated error
                strcat('WARNING: Which central tendency to use for ',...
                'estimating overall score reliability\n',...
                'Please ented a valid input for relcentmeas\n',...
                '1 - mean\n',...
                '2 - median\n',...
                'See help psyrat_relsummary for more information \n'));
        end

        ind = find(strcmpi('CI',varargin),1);
        if ~isempty(ind)
            ciperc = varargin{ind+1};
            if ciperc > 1 || ciperc < 0
                error('varargin:ci',... %Error code and associated error
                    strcat('WARNING: Size of credible interval should ',...
                    'be a value between 0 and 1\n',...
                    'A value of ',sprintf(' %2.2f',ciperc),...
                    ' is invalid\n',...
                    'See help psyrat_relsummary for more information \n'));
            end
        else
            ciperc = .95; %default: 95%
        end
end

switch relanalysis
    case 'sing'
        %place psyrat_data.data in REL to work with
        REL = psyrat_data.rel;
        
        %create a data structure for storing outputs
        data = struct;
        
        %create variable to specify whether there are bad/unreliable data
        %default state: 0, will be changed to 1 if there is a problem
        relerr = struct();
        relerr.trlcutoff = 0;
        relerr.trlmax = 0;
        relerr.nogooddata = 0;
        
        %check whether any groups exist
        if strcmpi(REL.groups,'none')
            ngroups = 1;
            %gnames = cellstr(REL.groups);
            gnames ={''};
        else
            ngroups = length(REL.groups);
            gnames = REL.groups(:);
        end
        
        %check whether any events exist
        if strcmpi(REL.events,'none')
            nevents = 1;
            %enames = cellstr(REL.events);
            enames = {''};
        else
            nevents = length(REL.events);
            enames = REL.events(:);
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
        
        %extract information from REL and store in data for crunching.
        %Backward-compat sig_trl default (shared helper); the one-facet model's
        %sig_e fixes the shape so wp = sig_e and the trial term is 0 for legacy
        %structs.
        REL.out = local_default_sigtrl(REL.out, REL.out.sig_e);
        switch analysis
            case 1 %1 - no groups or event types to consider

                gloc = 1;
                eloc = 1;

                data.g(gloc).e(eloc).label = REL.out.labels(gloc);
                data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,gloc);
                data.g(gloc).e(eloc).sig_u.raw = REL.out.sig_u(:,gloc);
                data.g(gloc).e(eloc).sig_e.raw = REL.out.sig_e(:,gloc);
                data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,gloc);
                data.g(gloc).e(eloc).elabel = cellstr('none');
                data.g(gloc).glabel = gnames(gloc);
                
            case 2 %2 - possible multiple groups but no event types to consider
                
                eloc = 1;
                
                for gloc=1:length(REL.out.labels)

                    data.g(gloc).e(eloc).label = REL.out.labels(gloc);
                    data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,gloc);
                    data.g(gloc).e(eloc).sig_u.raw = REL.out.sig_u(:,gloc);
                    data.g(gloc).e(eloc).sig_e.raw = REL.out.sig_e(:,gloc);
                    data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,gloc);
                    data.g(gloc).e(eloc).elabel = cellstr('none');
                    data.g(gloc).glabel = gnames(gloc);
                    
                end
                
            case 3 %3 - possible event types but no groups to consider
                
                gloc = 1;
                
                for eloc=1:length(REL.out.labels)
                    
                    data.g(gloc).e(eloc).label = REL.out.labels(eloc);
                    data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,eloc);
                    data.g(gloc).e(eloc).sig_u.raw = REL.out.sig_u(:,eloc);
                    data.g(gloc).e(eloc).sig_e.raw = REL.out.sig_e(:,eloc);
                    data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,eloc);
                    data.g(gloc).e(eloc).elabel = enames(eloc);
                    data.g(gloc).glabel = gnames(gloc);
                    
                end
                
            case 4 %4 - possible groups and event types to consider
                for i=1:length(REL.out.labels)
                    
                    %use the underscores that were added in psyrat_computevarcomp to
                    %differentiate where the group and event the data are for
                    lblstr = strsplit(REL.out.labels{i},'_;_');
                    
                    eloc = find(ismember(enames,lblstr(1)));
                    gloc = find(ismember(gnames,lblstr(2)));
                    
                    data.g(gloc).e(eloc).label = REL.out.labels(i);
                    data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,i);
                    data.g(gloc).e(eloc).sig_u.raw = REL.out.sig_u(:,i);
                    data.g(gloc).e(eloc).sig_e.raw = REL.out.sig_e(:,i);
                    data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,i);
                    data.g(gloc).e(eloc).elabel = enames(eloc);
                    data.g(gloc).glabel = gnames(gloc);
                    
                end
        end %switch analysis
        
        %store the cutoff in relsummary to pass to other functions more easily
        relsummary.depcutoff = depcutoff;
        relsummary.ciperc = ciperc;

        %Credible-interval edge shared by every variance-component quantile
        %below. The reliability coefficients already honor the user's CI
        %(threaded as 'CI',ciperc); the between-/within-person SDs used to
        %read literal .025/.975, so a run with a non-default CI produced one
        %table whose columns carried two different coverages.
        %
        %(1-ciperc)./2 is the edge convention the rest of the toolbox already
        %uses -- psyrat_rel_sing.m:109, psyrat_cutscore_sing.m:114 and
        %PsyRATAccuracyOracle all spell it this way -- so the component SDs are
        %joining it rather than departing from it.
        %
        %DISCLOSED CONSEQUENCE: on the .95 default this is 6 ULP above the
        %literal .025 it replaces, so a component SD bound can move by ~1 ULP
        %against previously saved output (measured: betsd unchanged, witsd 1
        %ULP on the test fixture). Reported to two decimals, so nothing visible
        %changes; but it is NOT bit-identical, and TestCiCoverageConsistency
        %pins the size of the drift rather than pretending there is none.
        ciedge = (1-ciperc)./2;
        relsummary.gcoeff = gcoeff;
        if gcoeff == 1
            relsummary.gcoeff_name = 'dep';
        else
            relsummary.gcoeff_name = 'gen';
        end
        
        %store which measure was used to specify cutoff
        switch meascutoff
            case 1
                relsummary.meascutoff = strcat('Lower Limit of ',...
                    sprintf(' %2.0f',ciperc*100),'% Credible Interval');
            case 2
                relsummary.meascutoff = 'Point Estimate';
            case 3
                relsummary.meascutoff = strcat('Upper Limit of ',...
                    sprintf(' %2.0f',ciperc*100),'% Credible Interval');
        end
        
        %compute dependability data for each group and event
        switch analysis
            case 1 %no groups or event types to consider
                
                %the same generic structure is used for relsummary, so the event
                %and group locations will both be 1 for the data
                eloc = 1;
                gloc = 1;
                
                %grab the group names (if present)
                relsummary.group(gloc).name = gnames{gloc};
                
                %set the event name as measure
                relsummary.group(gloc).event(eloc).name = 'measure';
                
                try
                    %create empty arrays for storing dependability information
                    trltable = varfun(@length,REL.data{1},'GroupingVariables',{'id'});
                catch
                    %create empty arrays for storing dependability information
                    trltable = varfun(@length,REL.data,'GroupingVariables',{'id'});
                end
                
                %compute dependability
                ntrials = max(trltable.GroupCount(:)) + 1000;
                
                [llrel,mrel,ulrel] = psyrat_rel_sing('gcoeff',gcoeff,'metric','global',...
                    'bp',data.g(gloc).e(eloc).sig_u.raw,...
                    'wp',local_total_within(data.g(gloc).e(eloc).sig_e.raw, data.g(gloc).e(eloc).sig_trl.raw),'i',data.g(gloc).e(eloc).sig_trl.raw,...
                    'obs',[1 ntrials],'CI',ciperc);
                
                %find the number of trials to reach cutoff
                switch meascutoff
                    case 1
                        trlcutoff = find(llrel >= depcutoff, 1);
                    case 2
                        trlcutoff = find(mrel >= depcutoff, 1);
                    case 3
                        trlcutoff = find(ulrel >= depcutoff, 1);
                end
                
                %see whether the trial cutoff was found. If not store all values as
                %-1. If it is found, store the dependability information about the
                %cutoffs.
                if isempty(trlcutoff)
                    
                    relerr.trlcutoff = 1;
                    trlcutoff = -1;
                    relsummary.group(gloc).event(eloc).trlcutoff = -1;
                    relsummary.group(gloc).event(eloc).rel.m = -1;
                    relsummary.group(gloc).event(eloc).rel.ll = -1;
                    relsummary.group(gloc).event(eloc).rel.ul = -1;
                    
                elseif ~isempty(trlcutoff)
                    
                    relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                    relsummary.group(gloc).event(eloc).rel.m = mrel(trlcutoff);
                    relsummary.group(gloc).event(eloc).rel.ll = llrel(trlcutoff);
                    relsummary.group(gloc).event(eloc).rel.ul = ulrel(trlcutoff);
                    
                end
                
                %find the participants without enough trials based on the cutoffs
                %the guard tests the -1 sentinel VALUE, not isempty: the
                %storage guard above (the if/elseif that records the cutoff,
                %two guard clauses up) reassigned trlcutoff to -1 when the
                %cutoff was unreachable, so an isempty test here can never
                %enter this arm and the fall-through marked every participant
                %good against a GroupCount >= -1 comparison (G54). An
                %unreachable cutoff places every participant in Bad IDs,
                %matching cases 2-4 in both families, the trt_diff helper's
                %sentinel arm, and the manual's ch. 9 promise. (trt_sserr is
                %NOT in that list: it has no trial-cutoff search and no
                %sentinel; its route marks inclusion per subject.)
                if trlcutoff == -1 %if the cutoff was not found in the data

                    %Get information for the participants. REL.data arrives as
                    %a bare table on most estimation cases and as a per-event
                    %cell on psyrat_computevarcomp's event-array cases; a cell
                    %reaching analysis 1 holds exactly one table (multi-table
                    %cells force nevents > 1 and route to analyses 3/4)
                    datatrls = REL.data;
                    if iscell(datatrls)
                        datatrls = REL.data{1};
                    end

                    trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});
                    
                    %define all of the data as bad, because the cutoff wasn't even
                    %found
                    relsummary.group(gloc).event(eloc).eventgoodids = cell(0,1);
                    relsummary.group(gloc).event(eloc).eventbadids =...
                        trltable.id;
                    
                    relsummary.group(gloc).goodids = cell(0,1);
                    relsummary.group(gloc).badids = ...
                        relsummary.group(gloc).event(eloc).eventbadids;
                    
                    %calculate the dependability estimates for the overall data
                    datatable = REL.data;
                    if iscell(datatable)
                        datatable = REL.data{1};
                    end

                    badids = table(relsummary.group(gloc).badids);
                    
                    %check whether there are not enough good data after applying
                    %reliability threshold and extrapolating
                    if isempty(badids)
                        dlg = {'Data do not reach reliability threshold after extrapolation';...
                            'Set a lower reliability threshold'};
                        local_relsummary_notice(interactive, dlg);
                        relerr.nogooddata = 1;
                        return;
                    end
                    
                    trlcdata = innerjoin(datatable, badids,...
                        'LeftKeys', 'id', 'RightKeys', 'Var1',...
                        'LeftVariables', {'id' 'meas'});
                    
                    trltable = varfun(@length,trlcdata,...
                        'GroupingVariables',{'id'});
                    
                    trlmean = mean(trltable.GroupCount);
                    trlmed = median(trltable.GroupCount);
                    
                    %calculate dependability using either the mean or median (based
                    %on user input)
                    switch depcentmeas
                        case 1
                            depcent = trlmean;
                        case 2
                            depcent = trlmed;
                    end
                    
                    
                    [relsummary.group(gloc).event(eloc).dep.ll,...
                        relsummary.group(gloc).event(eloc).dep.m,...
                        relsummary.group(gloc).event(eloc).dep.ul] =...
                        psyrat_rel_sing('gcoeff',gcoeff,'metric','global','bp',data.g(gloc).e(eloc).sig_u.raw,...
                        'wp',local_total_within(data.g(gloc).e(eloc).sig_e.raw, data.g(gloc).e(eloc).sig_trl.raw),'i',data.g(gloc).e(eloc).sig_trl.raw,...
                        'obs',depcent,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).dep.meas = depcent;
                    
                    relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                    relsummary.group(gloc).event(eloc).goodn = 0;
                    
                else

                    %Get trial information for those that meet cutoff
                    datatrls = REL.data;

                    try
                        trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});
                    catch
                        trltable = varfun(@length,datatrls{:},'GroupingVariables',{'id'});
                    end
                    
                    ind2include = trltable.GroupCount >= trlcutoff;
                    ind2exclude = trltable.GroupCount < trlcutoff;
                    
                    %store the good ids and the bad ids (ids that don't meet the
                    %cutoff)
                    relsummary.group(gloc).event(eloc).eventgoodids =...
                        trltable.id(ind2include);
                    relsummary.group(gloc).event(eloc).eventbadids =...
                        trltable.id(ind2exclude);
                    
                    relsummary.group(gloc).goodids = ...
                        relsummary.group(gloc).event(eloc).eventgoodids;
                    relsummary.group(gloc).badids = ...
                        relsummary.group(gloc).event(eloc).eventbadids;
                    
                    datatable = REL.data;
                    
                    if iscell(datatable)
                        datatable = REL.data{1};
                    end
                    
                    goodids = table(relsummary.group(gloc).goodids);
                    
                    trlcdata = innerjoin(datatable, goodids,...
                        'LeftKeys', 'id', 'RightKeys', 'Var1',...
                        'LeftVariables', {'id' 'meas'});
                    
                    trltable = varfun(@length,trlcdata,...
                        'GroupingVariables',{'id'});
                    
                    trlmean = mean(trltable.GroupCount);
                    trlmed = median(trltable.GroupCount);
                    
                    %calculate dependability using either the mean or median (based
                    %on user input)
                    switch depcentmeas
                        case 1
                            depcent = trlmean;
                        case 2
                            depcent = trlmed;
                    end
                    
                    [relsummary.group(gloc).event(eloc).dep.ll,...
                        relsummary.group(gloc).event(eloc).dep.m,...
                        relsummary.group(gloc).event(eloc).dep.ul] =...
                        psyrat_rel_sing('gcoeff',gcoeff,'metric','global','bp',data.g(gloc).e(eloc).sig_u.raw,...
                        'wp',local_total_within(data.g(gloc).e(eloc).sig_e.raw, data.g(gloc).e(eloc).sig_trl.raw),'i',data.g(gloc).e(eloc).sig_trl.raw,...
                        'obs',depcent,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).dep.meas = depcent;
                    
                    relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                    relsummary.group(gloc).event(eloc).goodn = height(goodids);
                    
                    %just in case a cutoff was extrapolated, the dependability
                    %information will be overwritten
                    if  relsummary.group(gloc).event(eloc).trlcutoff >...
                            relsummary.group(gloc).event(eloc).trlinfo.max
                        
                        relerr.trlmax = 1;
                        relsummary.group(gloc).event(eloc).rel.m = -1;
                        relsummary.group(gloc).event(eloc).rel.ll = -1;
                        relsummary.group(gloc).event(eloc).rel.ul = -1;
                        
                    end
                    
                    
                end
                
                %calculate iccs, between-subject standard deviations, and
                %within-subject standard deviations
                [relsummary.group(gloc).event(eloc).icc.ll,...
                    relsummary.group(gloc).event(eloc).icc.m,...
                    relsummary.group(gloc).event(eloc).icc.ul] = ...
                    psyrat_rel_sing('gcoeff',gcoeff,'metric','icc','bp',data.g(gloc).e(eloc).sig_u.raw,...
                    'wp',local_total_within(data.g(gloc).e(eloc).sig_e.raw, data.g(gloc).e(eloc).sig_trl.raw),'i',data.g(gloc).e(eloc).sig_trl.raw,'CI',ciperc);
                
                relsummary.group(gloc).event(eloc).betsd.m = ...
                    mean(data.g(gloc).e(eloc).sig_u.raw);
                relsummary.group(gloc).event(eloc).betsd.ll = ...
                    quantile(data.g(gloc).e(eloc).sig_u.raw,ciedge);
                relsummary.group(gloc).event(eloc).betsd.ul = ...
                    quantile(data.g(gloc).e(eloc).sig_u.raw,1-ciedge);
                
                relsummary.group(gloc).event(eloc).witsd.m = ...
                    mean(data.g(gloc).e(eloc).sig_e.raw);
                relsummary.group(gloc).event(eloc).witsd.ll = ...
                    quantile(data.g(gloc).e(eloc).sig_e.raw,ciedge);
                relsummary.group(gloc).event(eloc).witsd.ul = ...
                    quantile(data.g(gloc).e(eloc).sig_e.raw,1-ciedge);
                
                
            case 2 %possible multiple groups but no event types to consider
                
                %since the same generic structure is used for relsummary, the event
                %location will be defined as 1.
                eloc = 1;
                
                for gloc=1:ngroups %loop through each group
                    
                    %store the name of the group
                    if eloc == 1
                        relsummary.group(gloc).name = gnames{gloc};
                    end
                    
                    %store the name of the event
                    relsummary.group(gloc).event(eloc).name = enames{eloc};
                    
                    try
                        %create empty arrays for storing dependability information
                        trltable = varfun(@length,REL.data{1},...
                            'GroupingVariables',{'id'});
                    catch
                        %create empty arrays for storing dependability information
                        trltable = varfun(@length,REL.data,...
                            'GroupingVariables',{'id'});
                    end
                    
                    ntrials = max(trltable.GroupCount(:)) + 1000;
                    
                    %compute dependability
                    
                    
                    [llrel,mrel,ulrel] = psyrat_rel_sing('gcoeff',gcoeff,'metric','global',...
                        'bp',data.g(gloc).e(eloc).sig_u.raw,...
                        'wp',local_total_within(data.g(gloc).e(eloc).sig_e.raw, data.g(gloc).e(eloc).sig_trl.raw),'i',data.g(gloc).e(eloc).sig_trl.raw,...
                        'obs',[1 ntrials],'CI',ciperc);
                    
                    %find the number of trials to reach cutoff
                    switch meascutoff
                        case 1
                            trlcutoff = find(llrel >= depcutoff, 1);
                        case 2
                            trlcutoff = find(mrel >= depcutoff, 1);
                        case 3
                            trlcutoff = find(ulrel >= depcutoff, 1);
                    end
                    
                    %see whether the trial cutoff was found. If not store all
                    %values as -1. If it is found, store the dependability
                    %information about the cutoffs.
                    
                    if isempty(trlcutoff)
                        
                        relerr.trlcutoff = 1;
                        trlcutoff = -1;
                        relsummary.group(gloc).event(eloc).trlcutoff = -1;
                        relsummary.group(gloc).event(eloc).rel.m = -1;
                        relsummary.group(gloc).event(eloc).rel.ll = -1;
                        relsummary.group(gloc).event(eloc).rel.ul = -1;
                        
                    elseif ~isempty(trlcutoff)
                        
                        %store information about cutoffs
                        relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                        relsummary.group(gloc).event(eloc).rel.m = mrel(trlcutoff);
                        relsummary.group(gloc).event(eloc).rel.ll = llrel(trlcutoff);
                        relsummary.group(gloc).event(eloc).rel.ul = ulrel(trlcutoff);
                        
                    end
                    
                    %if the trial cutoff was not found in the specified data
                    if trlcutoff == -1
                        
                        %store all the ids as bad
                        try
                            datatrls = REL.data{1};
                        catch
                            datatrls = REL.data;
                        end
                        ind = strcmp(datatrls.group,gnames{gloc});
                        datatrls = datatrls(ind,:);
                        
                        trltable = varfun(@length,datatrls,...
                            'GroupingVariables',{'id'});
                        
                        relsummary.group(gloc).event(eloc).eventgoodids = cell(0,1);
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            trltable.id;
                        
                    else
                        
                        %store the ids for participants with enough trials as good,
                        %and the ids for participants with too few trials as bad
                        try
                            datatrls = REL.data{1};
                        catch
                            datatrls = REL.data;
                        end
                        
                        ind = strcmp(datatrls.group,gnames{gloc});
                        datatrls = datatrls(ind,:);
                        
                        trltable = varfun(@length,datatrls,...
                            'GroupingVariables',{'id'});
                        
                        ind2include = trltable.GroupCount >= trlcutoff;
                        ind2exclude = trltable.GroupCount < trlcutoff;
                        
                        relsummary.group(gloc).event(eloc).eventgoodids =...
                            trltable.id(ind2include);
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            trltable.id(ind2exclude);
                        
                    end
                    
                    
                end
                
                %since there is only 1 event, those participants with good data for
                %the event will be store as participants having good data for the
                %whole group (this step makes more sense when there is more than
                %one event per group)
                
                for gloc=1:ngroups
                    
                    relsummary.group(gloc).goodids = ...
                        relsummary.group(gloc).event(eloc).eventgoodids;
                    relsummary.group(gloc).badids = ...
                        relsummary.group(gloc).event(eloc).eventbadids;
                    
                end
                
                %calculate the dependability for the overall data for each group
                for gloc=1:ngroups
                    
                    try
                        datatable = REL.data{1};
                    catch
                        datatable = REL.data;
                    end
                    ind = strcmp(datatable.group,gnames{gloc});
                    datasubset = datatable(ind,:);
                    
                    %only factor in the trial counts from those with good data
                    if ~isempty(relsummary.group(gloc).goodids)
                        goodids = table(relsummary.group(gloc).goodids);
                    elseif isempty(relsummary.group(gloc).goodids)
                        goodids = table(relsummary.group(gloc).badids);
                    end
                    
                    %check whether there are not enough good data after applying
                    %reliability threshold and extrapolating
                    if isempty(goodids)
                        dlg = {'Data do not reach reliability threshold after extrapolation';...
                            'Set a lower reliability threshold'};
                        local_relsummary_notice(interactive, dlg);
                        relerr.nogooddata = 1;
                        return;
                    end
                    
                    trlcdata = innerjoin(datasubset, goodids,...
                        'LeftKeys', 'id', 'RightKeys', 'Var1',...
                        'LeftVariables', {'id' 'meas'});
                    
                    trltable = varfun(@length,trlcdata,...
                        'GroupingVariables',{'id'});
                    
                    trlmean = mean(trltable.GroupCount);
                    trlmed = median(trltable.GroupCount);
                    
                    %calculate dependability using either the mean or median (based
                    %on user input)
                    switch depcentmeas
                        case 1
                            depcent = trlmean;
                        case 2
                            depcent = trlmed;
                    end
                    
                    [relsummary.group(gloc).event(eloc).dep.ll,...
                        relsummary.group(gloc).event(eloc).dep.m,...
                        relsummary.group(gloc).event(eloc).dep.ul] =...
                        psyrat_rel_sing('gcoeff',gcoeff,'metric','global','bp',data.g(gloc).e(eloc).sig_u.raw,...
                        'wp',local_total_within(data.g(gloc).e(eloc).sig_e.raw, data.g(gloc).e(eloc).sig_trl.raw),'i',data.g(gloc).e(eloc).sig_trl.raw,...
                        'obs',depcent,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).dep.meas = depcent;
                    relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                    
                    %just in case a cutoff was extrapolated, the dependability
                    %information will be overwritten
                    if  relsummary.group(gloc).event(eloc).trlcutoff >...
                            relsummary.group(gloc).event(eloc).trlinfo.max
                        
                        relerr.trlmax = 1;
                        relsummary.group(gloc).event(eloc).rel.m = -1;
                        relsummary.group(gloc).event(eloc).rel.ll = -1;
                        relsummary.group(gloc).event(eloc).rel.ul = -1;
                        
                    end
                    
                    if ~isempty(relsummary.group(gloc).goodids)
                        relsummary.group(gloc).event(eloc).goodn = height(goodids);
                    elseif isempty(relsummary.group(gloc).goodids)
                        relsummary.group(gloc).event(eloc).goodn = 0;
                    end
                    
                    %calculate iccs, between-subject standard deviations, and
                    %within-subject standard deviations
                    
                    [relsummary.group(gloc).event(eloc).icc.ll,...
                        relsummary.group(gloc).event(eloc).icc.m,...
                        relsummary.group(gloc).event(eloc).icc.ul] = ...
                        psyrat_rel_sing('gcoeff',gcoeff,'metric','icc','bp',data.g(gloc).e(eloc).sig_u.raw,...
                        'wp',local_total_within(data.g(gloc).e(eloc).sig_e.raw, data.g(gloc).e(eloc).sig_trl.raw),'i',data.g(gloc).e(eloc).sig_trl.raw,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).betsd.m = ...
                        mean(data.g(gloc).e(eloc).sig_u.raw);
                    relsummary.group(gloc).event(eloc).betsd.ll = ...
                        quantile(data.g(gloc).e(eloc).sig_u.raw,ciedge);
                    relsummary.group(gloc).event(eloc).betsd.ul = ...
                        quantile(data.g(gloc).e(eloc).sig_u.raw,1-ciedge);
                    
                    relsummary.group(gloc).event(eloc).witsd.m = ...
                        mean(data.g(gloc).e(eloc).sig_e.raw);
                    relsummary.group(gloc).event(eloc).witsd.ll = ...
                        quantile(data.g(gloc).e(eloc).sig_e.raw,ciedge);
                    relsummary.group(gloc).event(eloc).witsd.ul = ...
                        quantile(data.g(gloc).e(eloc).sig_e.raw,1-ciedge);
                    
                    
                end
                
                
            case 3 %possible event types but no groups to consider
                
                %since the relsummary structure is generic for any number of groups
                %or events, the data will count as being from 1 group.
                gloc = 1;
                
                %cylce through each event
                for eloc=1:nevents
                    
                    %get the group name
                    if eloc == 1
                        relsummary.group(gloc).name = gnames{gloc};
                    end
                    
                    %get the event name
                    relsummary.group(gloc).event(eloc).name = enames{eloc};
                    
                    %create empty arrays for storing dependability information
                    trltable = varfun(@length,REL.data{1},...
                        'GroupingVariables',{'id'});
                    
                    %compute dependability
                    ntrials = max(trltable.GroupCount(:)) + 1000;
                    
                    
                    
                    [llrel,mrel,ulrel] = psyrat_rel_sing('gcoeff',gcoeff,'metric','global',...
                        'bp',data.g(gloc).e(eloc).sig_u.raw,...
                        'wp',local_total_within(data.g(gloc).e(eloc).sig_e.raw, data.g(gloc).e(eloc).sig_trl.raw),'i',data.g(gloc).e(eloc).sig_trl.raw,...
                        'obs',[1 ntrials],'CI',ciperc);
                    
                    %find the number of trials to reach cutoffl
                    switch meascutoff
                        case 1
                            trlcutoff = find(llrel >= depcutoff, 1);
                        case 2
                            trlcutoff = find(mrel >= depcutoff, 1);
                        case 3
                            trlcutoff = find(ulrel >= depcutoff, 1);
                    end
                    
                    
                    if isempty(trlcutoff) %if a cutoff wasn't found
                        
                        relerr.trlcutoff = 1;
                        trlcutoff = -1;
                        relsummary.group(gloc).event(eloc).trlcutoff = -1;
                        relsummary.group(gloc).event(eloc).rel.m = -1;
                        relsummary.group(gloc).event(eloc).rel.ll = -1;
                        relsummary.group(gloc).event(eloc).rel.ul = -1;
                        
                        datatrls = REL.data{eloc};
                        
                        trltable = varfun(@length,datatrls,...
                            'GroupingVariables',{'id'});
                        
                        %store all of the participant ids as bad, because none
                        %reached the cutoff
                        relsummary.group(gloc).event(eloc).eventgoodids =...
                            cell(0,1);
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            trltable.id;
                        
                    elseif ~isempty(trlcutoff) %if a cutoff was found
                        
                        %store information about cutoffs
                        relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                        relsummary.group(gloc).event(eloc).rel.m = mrel(trlcutoff);
                        relsummary.group(gloc).event(eloc).rel.ll = llrel(trlcutoff);
                        relsummary.group(gloc).event(eloc).rel.ul = ulrel(trlcutoff);
                        
                        datatrls = REL.data{eloc};
                        
                        trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});
                        
                        ind2include = trltable.GroupCount >= trlcutoff;
                        ind2exclude = trltable.GroupCount < trlcutoff;
                        
                        %store which participants had enough trials to meet cutoff
                        relsummary.group(gloc).event(eloc).eventgoodids =...
                            trltable.id(ind2include);
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            trltable.id(ind2exclude);
                        
                    end
                    
                end
                
                %cycle through the events and store information about which
                %participants have good data across all event types.
                tempids = {};
                badids = [];
                for eloc=1:nevents
                    %STRICT all-events rule (G39 ruling 2026-08-29): an event that
                    %retained nobody folds in as an empty set and empties the running
                    %intersection instead of being skipped; its eventbadids (the full
                    %id list for a dead event, by producer convention) fold in too, so
                    %the group good/bad ids complement.
                    tempids{end+1} = relsummary.group(gloc).event(eloc).eventgoodids;
                    if eloc == 1
                        badids = relsummary.group(gloc).event(eloc).eventbadids;
                    elseif eloc > 1
                        [~,ind]=setdiff(tempids{1},tempids{end});
                        new = tempids{1};
                        badids = vertcat(badids,...
                            relsummary.group(gloc).event(eloc).eventbadids,...
                            new(ind));
                        new(ind) = [];
                        tempids{1} = new;
                    end
                end
                
                %defensive: with the unconditional fold above, tempids is empty only when nevents == 0; keep the empty-cell carrier convention
                if isempty(tempids)
                    tempids{1} = cell(0,1);
                end
                
                %store information about participants with good/bad ids
                relsummary.group(gloc).goodids = tempids{1};
                relsummary.group(gloc).badids = unique(badids);
                %STRICT-fold consequence (G39): a group can retain nobody while
                %its per-event rows stay informative, so warn on the console and
                %continue. Deliberately NOT local_relsummary_notice: its modal
                %contract requires every caller to set relerr.nogooddata and
                %return, and this is not an abort.
                if ~(iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids))
                    fprintf(['\nNo participants met the reliability threshold ',...
                        'for every event type in a group; that group''s ',...
                        'n Included is 0.\n']);
                end
                
                %cycle through each event
                for eloc=1:nevents
                    
                    %get trial information to use for dependability estimates. only
                    %participants with enough data will contribute toward trial
                    %counts
                    datatable = REL.data{eloc};
                    
                    if iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids)
                        goodids = table(relsummary.group(gloc).goodids);
                    else
                        goodids = table(relsummary.group(gloc).badids);
                    end
                    
                    %check whether there are not enough good data after applying
                    %reliability threshold and extrapolating
                    if isempty(goodids)
                        dlg = {'Data do not reach reliability threshold after extrapolation';...
                            'Set a lower reliability threshold'};
                        local_relsummary_notice(interactive, dlg);
                        relerr.nogooddata = 1;
                        return;
                    end
                    
                    trlcdata = innerjoin(datatable, goodids,...
                        'LeftKeys', 'id', 'RightKeys', 'Var1',...
                        'LeftVariables', {'id' 'meas'});
                    
                    trltable = varfun(@length,trlcdata,...
                        'GroupingVariables',{'id'});
                    
                    trlmean = mean(trltable.GroupCount);
                    trlmed = median(trltable.GroupCount);
                    
                    %calculate dependability using either the mean or median (based
                    %on user input)
                    switch depcentmeas
                        case 1
                            depcent = trlmean;
                        case 2
                            depcent = trlmed;
                    end
                    
                    [relsummary.group(gloc).event(eloc).dep.ll,...
                        relsummary.group(gloc).event(eloc).dep.m,...
                        relsummary.group(gloc).event(eloc).dep.ul] =...
                        psyrat_rel_sing('gcoeff',gcoeff,'metric','global','bp',data.g(gloc).e(eloc).sig_u.raw,...
                        'wp',local_total_within(data.g(gloc).e(eloc).sig_e.raw, data.g(gloc).e(eloc).sig_trl.raw),'i',data.g(gloc).e(eloc).sig_trl.raw,...
                        'obs',depcent,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).dep.meas = depcent;
                    
                    relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                    
                    %just in case a cutoff was extrapolated, the dependability
                    %information will be overwritten
                    if  relsummary.group(gloc).event(eloc).trlcutoff >...
                            relsummary.group(gloc).event(eloc).trlinfo.max
                        
                        relerr.trlmax = 1;
                        relsummary.group(gloc).event(eloc).rel.m = -1;
                        relsummary.group(gloc).event(eloc).rel.ll = -1;
                        relsummary.group(gloc).event(eloc).rel.ul = -1;
                        
                    end
                    
                    if iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids)
                        relsummary.group(gloc).event(eloc).goodn = height(goodids);
                    else
                        relsummary.group(gloc).event(eloc).goodn = 0;
                    end
                    
                    %calculate iccs, between-subject standard deviations, and
                    %within-subject standard deviations
                    [relsummary.group(gloc).event(eloc).icc.ll,...
                        relsummary.group(gloc).event(eloc).icc.m,...
                        relsummary.group(gloc).event(eloc).icc.ul] = ...
                        psyrat_rel_sing('gcoeff',gcoeff,'metric','icc','bp',data.g(gloc).e(eloc).sig_u.raw,...
                        'wp',local_total_within(data.g(gloc).e(eloc).sig_e.raw, data.g(gloc).e(eloc).sig_trl.raw),'i',data.g(gloc).e(eloc).sig_trl.raw,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).betsd.m = ...
                        mean(data.g(gloc).e(eloc).sig_u.raw);
                    relsummary.group(gloc).event(eloc).betsd.ll = ...
                        quantile(data.g(gloc).e(eloc).sig_u.raw,ciedge);
                    relsummary.group(gloc).event(eloc).betsd.ul = ...
                        quantile(data.g(gloc).e(eloc).sig_u.raw,1-ciedge);
                    
                    relsummary.group(gloc).event(eloc).witsd.m = ...
                        mean(data.g(gloc).e(eloc).sig_e.raw);
                    relsummary.group(gloc).event(eloc).witsd.ll = ...
                        quantile(data.g(gloc).e(eloc).sig_e.raw,ciedge);
                    relsummary.group(gloc).event(eloc).witsd.ul = ...
                        quantile(data.g(gloc).e(eloc).sig_e.raw,1-ciedge);
                    
                end
                
            case 4 %groups and event types to consider
                
                %cycle through each event
                for eloc=1:nevents
                    
                    %cycle through each group
                    for gloc=1:ngroups
                        
                        %store the group name
                        if eloc == 1
                            relsummary.group(gloc).name = gnames{gloc};
                        end
                        
                        %store the event name
                        relsummary.group(gloc).event(eloc).name = enames{eloc};
                        
                        %create empty arrays for storing dependability information
                        trltable = varfun(@length,REL.data{1},...
                            'GroupingVariables',{'id'});
                        
                        %compute dependability
                        ntrials = max(trltable.GroupCount(:)) + 1000;
                        
                        
                        
                        [llrel,mrel,ulrel] = psyrat_rel_sing('gcoeff',gcoeff,'metric','global',...
                            'bp',data.g(gloc).e(eloc).sig_u.raw,...
                            'wp',local_total_within(data.g(gloc).e(eloc).sig_e.raw, data.g(gloc).e(eloc).sig_trl.raw),'i',data.g(gloc).e(eloc).sig_trl.raw,...
                            'obs',[1 ntrials],'CI',ciperc);
                        
                        %find the number of trials to reach cutoff
                        switch meascutoff
                            case 1
                                trlcutoff = find(llrel >= depcutoff, 1);
                            case 2
                                trlcutoff = find(mrel >= depcutoff, 1);
                            case 3
                                trlcutoff = find(ulrel >= depcutoff, 1);
                        end
                        
                        
                        %if a cutoff was not found
                        if isempty(trlcutoff)
                            
                            %store all the participant ids as having bad data
                            relerr.trlcutoff = 1;
                            trlcutoff = -1;
                            relsummary.group(gloc).event(eloc).trlcutoff = -1;
                            relsummary.group(gloc).event(eloc).rel.m = -1;
                            relsummary.group(gloc).event(eloc).rel.ll = -1;
                            relsummary.group(gloc).event(eloc).rel.ul = -1;
                            
                            datatrls = REL.data{eloc};
                            ind = strcmp(datatrls.group,gnames{gloc});
                            datatrls = datatrls(ind,:);
                            
                            trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});
                            
                            relsummary.group(gloc).event(eloc).eventgoodids =...
                                cell(0,1);
                            relsummary.group(gloc).event(eloc).eventbadids =...
                                trltable.id;
                            
                        elseif ~isempty(trlcutoff) %if a cutoff was found
                            
                            
                            relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                            relsummary.group(gloc).event(eloc).rel.m = mrel(trlcutoff);
                            relsummary.group(gloc).event(eloc).rel.ll = llrel(trlcutoff);
                            relsummary.group(gloc).event(eloc).rel.ul = ulrel(trlcutoff);
                            
                            datatrls = REL.data{eloc};
                            ind = strcmp(datatrls.group,gnames{gloc});
                            datatrls = datatrls(ind,:);
                            
                            %find ids with enough trials based on cutoff
                            trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});
                            
                            ind2include = trltable.GroupCount >= trlcutoff;
                            ind2exclude = trltable.GroupCount < trlcutoff;
                            
                            relsummary.group(gloc).event(eloc).eventgoodids =...
                                trltable.id(ind2include);
                            relsummary.group(gloc).event(eloc).eventbadids =...
                                trltable.id(ind2exclude);
                            
                        end
                        
                    end
                end
                
                %find the ids that have enough trials for each event type
                
                for gloc=1:ngroups
                    tempids = {};
                    badids = [];
                    
                    for eloc=1:nevents
                        %STRICT all-events rule (G39 ruling 2026-08-29): an event that
                        %retained nobody folds in as an empty set and empties the running
                        %intersection instead of being skipped; its eventbadids (the full
                        %id list for a dead event, by producer convention) fold in too, so
                        %the group good/bad ids complement.
                        tempids{end+1} = relsummary.group(gloc).event(eloc).eventgoodids;
                        if eloc == 1
                            badids = relsummary.group(gloc).event(eloc).eventbadids;
                        elseif eloc > 1
                            [~,ind]=setdiff(tempids{1},tempids{end});
                            new = tempids{1};
                            badids = vertcat(badids,...
                                relsummary.group(gloc).event(eloc).eventbadids,...
                                new(ind));
                            new(ind) = [];
                            tempids{1} = new;
                        end
                    end
                    
                    %defensive: with the unconditional fold above, tempids is empty only when nevents == 0; keep the empty-cell carrier convention
                    if isempty(tempids)
                        tempids{1} = cell(0,1);
                    end
                    
                    relsummary.group(gloc).goodids = tempids{1};
                    relsummary.group(gloc).badids = unique(badids);
                    %STRICT-fold consequence (G39): a group can retain nobody while
                    %its per-event rows stay informative, so warn on the console and
                    %continue. Deliberately NOT local_relsummary_notice: its modal
                    %contract requires every caller to set relerr.nogooddata and
                    %return, and this is not an abort.
                    if ~(iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids))
                        fprintf(['\nNo participants met the reliability threshold ',...
                            'for every event type in a group; that group''s ',...
                            'n Included is 0.\n']);
                    end
                    
                end
                
                %calculate dependability estimates for the overall data for each
                %group and event type
                for eloc=1:nevents
                    
                    for gloc=1:ngroups
                        
                        %pull the data to calculate the trial information from
                        %participants with good data
                        datatable = REL.data{eloc};
                        ind = strcmp(datatable.group,gnames{gloc});
                        datasubset = datatable(ind,:);
                        
                        if iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids)
                            goodids = table(relsummary.group(gloc).goodids);
                        else
                            goodids = table(relsummary.group(gloc).badids);
                        end
                        
                        %check whether there are not enough good data after applying
                        %reliability threshold and extrapolating
                        if isempty(goodids)
                            dlg = {'Data do not reach reliability threshold after extrapolation';...
                                'Set a lower reliability threshold'};
                            local_relsummary_notice(interactive, dlg);
                            relerr.nogooddata = 1;
                            return;
                        end
                        
                        trlcdata = innerjoin(datasubset, goodids,...
                            'LeftKeys', 'id', 'RightKeys', 'Var1',...
                            'LeftVariables', {'id' 'meas'});
                        
                        trltable = varfun(@length,trlcdata,...
                            'GroupingVariables',{'id'});
                        
                        trlmean = mean(trltable.GroupCount);
                        trlmed = median(trltable.GroupCount);
                        
                        %calculate dependability using either the mean or median
                        %(based on user input)
                        switch depcentmeas
                            case 1
                                depcent = trlmean;
                            case 2
                                depcent = trlmed;
                        end
                        
                        [relsummary.group(gloc).event(eloc).dep.ll,...
                            relsummary.group(gloc).event(eloc).dep.m,...
                            relsummary.group(gloc).event(eloc).dep.ul] =...
                            psyrat_rel_sing('gcoeff',gcoeff,'metric','global','bp',data.g(gloc).e(eloc).sig_u.raw,...
                            'wp',local_total_within(data.g(gloc).e(eloc).sig_e.raw, data.g(gloc).e(eloc).sig_trl.raw),'i',data.g(gloc).e(eloc).sig_trl.raw,...
                            'obs',depcent,'CI',ciperc);
                        
                        relsummary.group(gloc).event(eloc).dep.meas = depcent;
                        relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                        relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                        relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                        relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                        relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                        
                        %check if trial cutoff for an event exceeded the total
                        %number of trials for any subjects
                        if  relsummary.group(gloc).event(eloc).trlcutoff >...
                                relsummary.group(gloc).event(eloc).trlinfo.max
                            
                            relerr.trlmax = 1;
                            relsummary.group(gloc).event(eloc).rel.m = -1;
                            relsummary.group(gloc).event(eloc).rel.ll = -1;
                            relsummary.group(gloc).event(eloc).rel.ul = -1;
                            
                        end
                        
                        %goodn is the GROUP-level retained count (G39
                        %ruling 2026-08-29): keyed on the group carrier,
                        %exactly as the one-group branches key it, so a
                        %dead event's row and its siblings print the same
                        %complementing count. The event-carrier keying
                        %this replaces printed 0 on a dead event's row
                        %while its siblings printed the group figure, and
                        %under the strict fold it would have counted the
                        %badids fallback table instead.
                        if iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids)
                            relsummary.group(gloc).event(eloc).goodn = height(goodids);
                        else
                            relsummary.group(gloc).event(eloc).goodn = 0;
                        end
                        
                        %calculate iccs, between-subject standard deviations, and
                        %within-subject standard deviations
                        [relsummary.group(gloc).event(eloc).icc.ll,...
                            relsummary.group(gloc).event(eloc).icc.m,...
                            relsummary.group(gloc).event(eloc).icc.ul] = ...
                            psyrat_rel_sing('gcoeff',gcoeff,'metric','icc','bp',data.g(gloc).e(eloc).sig_u.raw,...
                            'wp',local_total_within(data.g(gloc).e(eloc).sig_e.raw, data.g(gloc).e(eloc).sig_trl.raw),'i',data.g(gloc).e(eloc).sig_trl.raw,'CI',ciperc);
                        
                        relsummary.group(gloc).event(eloc).betsd.m = ...
                            mean(data.g(gloc).e(eloc).sig_u.raw);
                        relsummary.group(gloc).event(eloc).betsd.ll = ...
                            quantile(data.g(gloc).e(eloc).sig_u.raw,ciedge);
                        relsummary.group(gloc).event(eloc).betsd.ul = ...
                            quantile(data.g(gloc).e(eloc).sig_u.raw,1-ciedge);
                        
                        relsummary.group(gloc).event(eloc).witsd.m = ...
                            mean(data.g(gloc).e(eloc).sig_e.raw);
                        relsummary.group(gloc).event(eloc).witsd.ll = ...
                            quantile(data.g(gloc).e(eloc).sig_e.raw,ciedge);
                        relsummary.group(gloc).event(eloc).witsd.ul = ...
                            quantile(data.g(gloc).e(eloc).sig_e.raw,1-ciedge);
                        
                    end
                end
                
                
        end %switch analysis

        %Report the trial main-effect SD (sigma_i) as its own component for the
        %one-facet designs. betsd (between-person) and witsd (the person x trial
        %residual) are set per case above; trlsd is added here uniformly for every
        %populated group/event cell that carries a sig_trl draw (the residual
        %witsd is reported separately rather than folded into a combined within SD).
        for psyrat_gg = 1:numel(relsummary.group)
            for psyrat_ee = 1:numel(relsummary.group(psyrat_gg).event)
                if numel(data.g) >= psyrat_gg ...
                        && numel(data.g(psyrat_gg).e) >= psyrat_ee ...
                        && isfield(data.g(psyrat_gg).e(psyrat_ee),'sig_trl') ...
                        && ~isempty(data.g(psyrat_gg).e(psyrat_ee).sig_trl)
                    relsummary.group(psyrat_gg).event(psyrat_ee).trlsd.m = ...
                        mean(data.g(psyrat_gg).e(psyrat_ee).sig_trl.raw);
                    relsummary.group(psyrat_gg).event(psyrat_ee).trlsd.ll = ...
                        quantile(data.g(psyrat_gg).e(psyrat_ee).sig_trl.raw,ciedge);
                    relsummary.group(psyrat_gg).event(psyrat_ee).trlsd.ul = ...
                        quantile(data.g(psyrat_gg).e(psyrat_ee).sig_trl.raw,1-ciedge);
                end
            end
        end

    case 'sing_diff'
        
        %place psyrat_data.data in REL to work with
        REL = psyrat_data.rel;
        
        %create a data structure for storing outputs
        data = struct;
        
        %create variable to specify whether there are bad/unreliable data
        %default state: 0, will be changed to 1 if there is a problem
        relerr = struct();
        relerr.trlcutoff = 0;
        relerr.trlmax = 0;
        relerr.nogooddata = 0;
        
        %check whether any groups exist
        if strcmpi(REL.groups,'none')
            ngroups = 1;
            %gnames = cellstr(REL.groups);
            gnames ={''};
        else
            ngroups = length(REL.groups);
            gnames = REL.groups(:);
        end
        
        %The data have have two events
        nevents = 2;
        enames = REL.events(:);
        
        
        %figure out whether groups or events need to be considered
        %1 - no groups or event types to consider
        %2 - possible multiple groups but no event types to consider
        %3 - possible event types but no groups to consider
        %4 - possible groups and event types to consider
        
        
        if ngroups == 1
            analysis = 3;
        elseif ngroups > 1
            analysis = 4;
        end
        
        %extract information from REL and store in data for crunching
        switch analysis
            case 3 %3 - possible event types but no groups to consider
                
                gloc = 1;
                
                data.g(gloc).e(1).b.raw = REL.out.b{:,1}(:,1);
                data.g(gloc).e(1).b_sigma.raw = REL.out.b_sigma{:,1}(:,1);
                data.g(gloc).e(1).sd_id.raw = REL.out.sd_id{:,1}(:,1);
                data.g(gloc).e(1).sd_trl.raw = REL.out.sd_trl{:,1}(:,1);
                data.g(gloc).e(1).id_varcov.raw = REL.out.id_varcov{:,1}(:,1);
                data.g(gloc).e(1).trl_varcov.raw = REL.out.trl_varcov{:,1}(:,1);
                data.g(gloc).e(1).elabel = enames(1);
                
                data.g(gloc).e(2).b.raw = REL.out.b{:,1}(:,2);
                data.g(gloc).e(2).b_sigma.raw = REL.out.b_sigma{:,1}(:,2);
                data.g(gloc).e(2).sd_id.raw = REL.out.sd_id{:,1}(:,2);
                data.g(gloc).e(2).sd_trl.raw = REL.out.sd_trl{:,1}(:,2);
                data.g(gloc).e(2).id_varcov.raw = REL.out.id_varcov{:,1}(:,2);
                data.g(gloc).e(2).trl_varcov.raw = REL.out.trl_varcov{:,1}(:,2);
                data.g(gloc).e(2).elabel = enames(2);
                data.g(gloc).glabel = gnames(gloc);
                
                
            case 4 %4 - possible groups and event types to consider
                for gloc=1:length(REL.out.glabels)
                    
                    data.g(gloc).e(1).b.raw = REL.out.b{:,gloc}(:,1);
                    data.g(gloc).e(1).b_sigma.raw = REL.out.b_sigma{:,gloc}(:,1);
                    data.g(gloc).e(1).sd_id.raw = REL.out.sd_id{:,gloc}(:,1);
                    data.g(gloc).e(1).sd_trl.raw = REL.out.sd_trl{:,gloc}(:,1);
                    data.g(gloc).e(1).id_varcov.raw = REL.out.id_varcov{:,gloc}(:,1);
                    data.g(gloc).e(1).trl_varcov.raw = REL.out.trl_varcov{:,gloc}(:,1);
                    data.g(gloc).e(1).elabel = enames(1);
                    
                    data.g(gloc).e(2).b.raw = REL.out.b{:,gloc}(:,2);
                    data.g(gloc).e(2).b_sigma.raw = REL.out.b_sigma{:,gloc}(:,2);
                    data.g(gloc).e(2).sd_id.raw = REL.out.sd_id{:,gloc}(:,2);
                    data.g(gloc).e(2).sd_trl.raw = REL.out.sd_trl{:,gloc}(:,2);
                    data.g(gloc).e(2).id_varcov.raw = REL.out.id_varcov{:,gloc}(:,2);
                    data.g(gloc).e(2).trl_varcov.raw = REL.out.trl_varcov{:,gloc}(:,2);
                    data.g(gloc).e(2).elabel = enames(2);
                    data.g(gloc).glabel = gnames(gloc);
                    
                end
        end %switch analysis
        
        %store the cutoff in relsummary to pass to other functions more easily
        relsummary.depcutoff = depcutoff;
        %The per-event coefficients below follow the user's difference-score
        %choice (diffgcoeff) rather than the population gcoeff pref. The viewer
        %for this analysis builds only a "Difference-Score Coefficient" popup
        %(psyrat_startview_sing) and no G-Theory Coefficient popup, so gcoeff
        %arrives here from a stored preference the user cannot see or change on
        %this screen -- and it can carry over from a previously viewed one-facet
        %run, applying a coefficient name to per-event numbers the user never
        %chose. Routing it to diffgcoeff makes the one visible control govern the
        %output it appears to govern, extending to the group-level path the
        %convention already established for the subject-level difference path
        %(see the same assignment in the subject-level helper below).
        gcoeff = diffgcoeff;
        relsummary.gcoeff = gcoeff;
        if gcoeff == 1
            relsummary.gcoeff_name = 'dep';
        else
            relsummary.gcoeff_name = 'gen';
        end
        relsummary.ciperc = ciperc;

        %edge for this arm's variance-component quantiles; see the 'sing' arm
        %above for why the component SDs follow ciperc rather than a fixed 95%
        ciedge = (1-ciperc)./2;
        relsummary.diffgcoeff = diffgcoeff;
        if diffgcoeff == 1
            relsummary.diffgcoeff_name = 'dep';
        else
            relsummary.diffgcoeff_name = 'gen';
        end
        
        %store which measure was used to specify cutoff
        switch meascutoff
            case 1
                relsummary.meascutoff = strcat('Lower Limit of ',...
                    sprintf(' %2.0f',ciperc*100),'% Credible Interval');
            case 2
                relsummary.meascutoff = 'Point Estimate';
            case 3
                relsummary.meascutoff = strcat('Upper Limit of ',...
                    sprintf(' %2.0f',ciperc*100),'% Credible Interval');
        end
        
        %compute dependability data for each group and event
        switch analysis
            
            case 3 %possible event types but no groups to consider
                
                %since the relsummary structure is generic for any number of groups
                %or events, the data will count as being from 1 group.
                gloc = 1;
                
                %cylce through each event
                for eloc=1:nevents
                    
                    %get the group name
                    if eloc == 1
                        relsummary.group(gloc).name = gnames{gloc};
                    end
                    
                    %get the event name
                    relsummary.group(gloc).event(eloc).name = enames{eloc};
                    
                    %create empty arrays for storing dependability information
                    trltable = varfun(@length,REL.data,...
                        'GroupingVariables',{'id','event'});
                    
                    trltable = trltable(strcmp(trltable.event,enames{eloc}),:);
                    
                    %compute dependability
                    ntrials = max(trltable.GroupCount(:)) + 1000;
                    

                    [llrel,mrel,ulrel] = psyrat_rel_sing('gcoeff',gcoeff,'metric','global',...
                        'bp',cell2mat(data.g(gloc).e(eloc).sd_id.raw),...
                        'wp',local_total_within(exp(cell2mat(data.g(gloc).e(eloc).b_sigma.raw)), cell2mat(data.g(gloc).e(eloc).sd_trl.raw)),'i',cell2mat(data.g(gloc).e(eloc).sd_trl.raw),...
                        'obs',[1 ntrials],...
                        'CI',ciperc);
                    
                    %find the number of trials to reach cutoffl
                    switch meascutoff
                        case 1
                            trlcutoff = find(llrel >= depcutoff, 1);
                        case 2
                            trlcutoff = find(mrel >= depcutoff, 1);
                        case 3
                            trlcutoff = find(ulrel >= depcutoff, 1);
                    end
                    
                    
                    if isempty(trlcutoff) %if a cutoff wasn't found
                        
                        relerr.trlcutoff = 1;
                        trlcutoff = -1;
                        relsummary.group(gloc).event(eloc).trlcutoff = -1;
                        relsummary.group(gloc).event(eloc).rel.m = -1;
                        relsummary.group(gloc).event(eloc).rel.ll = -1;
                        relsummary.group(gloc).event(eloc).rel.ul = -1;
                        
                        datatrls = REL.data;
                        
                        trltable = varfun(@length,datatrls,...
                            'GroupingVariables',{'id','event'});
                        
                        trltable = trltable(strcmp(trltable.event,enames{eloc}),:);
                        
                        %store all of the participant ids as bad, because none
                        %reached the cutoff
                        relsummary.group(gloc).event(eloc).eventgoodids =...
                            cell(0,1);
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            trltable.id;
                        
                    elseif ~isempty(trlcutoff) %if a cutoff was found
                        
                        %store information about cutoffs
                        relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                        relsummary.group(gloc).event(eloc).rel.m = mrel(trlcutoff);
                        relsummary.group(gloc).event(eloc).rel.ll = llrel(trlcutoff);
                        relsummary.group(gloc).event(eloc).rel.ul = ulrel(trlcutoff);
                        
                        datatrls = REL.data;
                        
                        trltable = varfun(@length,datatrls,...
                            'GroupingVariables',{'id','event'});
                        
                        trltable = trltable(strcmp(trltable.event,enames{eloc}),:);
                        
                        ind2include = trltable.GroupCount >= trlcutoff;
                        ind2exclude = trltable.GroupCount < trlcutoff;
                        
                        %store which participants had enough trials to meet cutoff
                        relsummary.group(gloc).event(eloc).eventgoodids =...
                            trltable.id(ind2include);
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            trltable.id(ind2exclude);
                        
                    end
                    
                end
                
                %cycle through the events and store information about which
                %participants have good data across all event types.
                tempids = {};
                badids = [];
                for eloc=1:nevents
                    %STRICT all-events rule (G39 ruling 2026-08-29): an event that
                    %retained nobody folds in as an empty set and empties the running
                    %intersection instead of being skipped; its eventbadids (the full
                    %id list for a dead event, by producer convention) fold in too, so
                    %the group good/bad ids complement.
                    tempids{end+1} = relsummary.group(gloc).event(eloc).eventgoodids;
                    if eloc == 1
                        badids = relsummary.group(gloc).event(eloc).eventbadids;
                    elseif eloc > 1
                        [~,ind]=setdiff(tempids{1},tempids{end});
                        new = tempids{1};
                        badids = vertcat(badids,...
                            relsummary.group(gloc).event(eloc).eventbadids,...
                            new(ind));
                        new(ind) = [];
                        tempids{1} = new;
                    end
                end
                
                %defensive: with the unconditional fold above, tempids is empty only when nevents == 0; keep the empty-cell carrier convention
                if isempty(tempids)
                    tempids{1} = cell(0,1);
                end
                
                %store information about participants with good/bad ids
                relsummary.group(gloc).goodids = tempids{1};
                relsummary.group(gloc).badids = unique(badids);
                %STRICT-fold consequence (G39): a group can retain nobody while
                %its per-event rows stay informative, so warn on the console and
                %continue. Deliberately NOT local_relsummary_notice: its modal
                %contract requires every caller to set relerr.nogooddata and
                %return, and this is not an abort.
                if ~(iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids))
                    fprintf(['\nNo participants met the reliability threshold ',...
                        'for every event type in a group; that group''s ',...
                        'n Included is 0.\n']);
                end
                
                %cycle through each event
                for eloc=1:nevents
                    
                    %get trial information to use for dependability estimates. only
                    %participants with enough data will contribute toward trial
                    %counts
                    datatable = REL.data;
                    datatable = datatable(...
                        strcmp(datatable.event,enames{eloc}),:);
                    
                    if iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids)
                        goodids = table(relsummary.group(gloc).goodids);
                    else
                        goodids = table(relsummary.group(gloc).badids);
                    end
                    
                    %check whether there are not enough good data after applying
                    %reliability threshold and extrapolating
                    if isempty(goodids)
                        dlg = {'Data do not reach reliability threshold after extrapolation';...
                            'Set a lower reliability threshold'};
                        local_relsummary_notice(interactive, dlg);
                        relerr.nogooddata = 1;
                        return;
                    end
                    
                    trlcdata = innerjoin(datatable, goodids,...
                        'LeftKeys', 'id', 'RightKeys', 'Var1',...
                        'LeftVariables', {'id' 'meas'});
                    
                    trltable = varfun(@length,trlcdata,...
                        'GroupingVariables',{'id'});
                    
                    trlmean = mean(trltable.GroupCount);
                    trlmed = median(trltable.GroupCount);
                    
                    %calculate dependability using either the mean or median (based
                    %on user input)
                    switch depcentmeas
                        case 1
                            depcent = trlmean;
                        case 2
                            depcent = trlmed;
                    end
                    
                    [relsummary.group(gloc).event(eloc).dep.ll,...
                        relsummary.group(gloc).event(eloc).dep.m,...
                        relsummary.group(gloc).event(eloc).dep.ul] =...
                        psyrat_rel_sing('gcoeff',gcoeff,'metric','global','bp',cell2mat(data.g(gloc).e(eloc).sd_id.raw),...
                        'wp',local_total_within(exp(cell2mat(data.g(gloc).e(eloc).b_sigma.raw)), cell2mat(data.g(gloc).e(eloc).sd_trl.raw)),'i',cell2mat(data.g(gloc).e(eloc).sd_trl.raw),...
                        'obs',depcent,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).dep.meas = depcent;
                    
                    relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                    
                    %just in case a cutoff was extrapolated, the dependability
                    %information will be overwritten
                    if  relsummary.group(gloc).event(eloc).trlcutoff >...
                            relsummary.group(gloc).event(eloc).trlinfo.max
                        
                        relerr.trlmax = 1;
                        relsummary.group(gloc).event(eloc).rel.m = -1;
                        relsummary.group(gloc).event(eloc).rel.ll = -1;
                        relsummary.group(gloc).event(eloc).rel.ul = -1;
                        
                    end
                    
                    if iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids)
                        relsummary.group(gloc).event(eloc).goodn = height(goodids);
                    else
                        relsummary.group(gloc).event(eloc).goodn = 0;
                    end
                    
                    %calculate iccs, between-subject standard deviations, and
                    %within-subject standard deviations
                    [relsummary.group(gloc).event(eloc).icc.ll,...
                        relsummary.group(gloc).event(eloc).icc.m,...
                        relsummary.group(gloc).event(eloc).icc.ul] = ...
                        psyrat_rel_sing('gcoeff',gcoeff,'metric','icc','bp',cell2mat(data.g(gloc).e(eloc).sd_id.raw),...
                        'wp',local_total_within(exp(cell2mat(data.g(gloc).e(eloc).b_sigma.raw)), cell2mat(data.g(gloc).e(eloc).sd_trl.raw)),'i',cell2mat(data.g(gloc).e(eloc).sd_trl.raw),...
                        'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).betsd.m = ...
                        mean(cell2mat(data.g(gloc).e(eloc).sd_id.raw));
                    relsummary.group(gloc).event(eloc).betsd.ll = ...
                        quantile(cell2mat(data.g(gloc).e(eloc).sd_id.raw),ciedge);
                    relsummary.group(gloc).event(eloc).betsd.ul = ...
                        quantile(cell2mat(data.g(gloc).e(eloc).sd_id.raw),1-ciedge);
                    
                    relsummary.group(gloc).event(eloc).witsd.m = ...
                        mean(exp(cell2mat(data.g(gloc).e(eloc).b_sigma.raw)));
                    relsummary.group(gloc).event(eloc).witsd.ll = ...
                        quantile(exp(cell2mat(data.g(gloc).e(eloc).b_sigma.raw)),...
                        ciedge);
                    relsummary.group(gloc).event(eloc).witsd.ul = ...
                        quantile(exp(cell2mat(data.g(gloc).e(eloc).b_sigma.raw)),...
                        1-ciedge);
                    
                end
               
                
            case 4 %groups and event types to consider
                
                %cycle through each event
                for eloc=1:nevents
                    
                    %cycle through each group
                    for gloc=1:ngroups
                        
                        %store the group name
                        if eloc == 1
                            relsummary.group(gloc).name = gnames{gloc};
                        end
                        
                        %store the event name
                        relsummary.group(gloc).event(eloc).name = enames{eloc};
                        
                        %create empty arrays for storing dependability information
                        trltable = varfun(@length,REL.data,...
                            'GroupingVariables',{'group','id','event'});
                        
                        trltable = trltable(strcmp(trltable.event,enames{eloc})...
                            & strcmp(trltable.group,gnames{gloc}),:);
                        
                        %compute dependability
                        ntrials = max(trltable.GroupCount(:)) + 1000;
                        
                        
                        
                        [llrel,mrel,ulrel] = psyrat_rel_sing('gcoeff',gcoeff,'metric','global',...
                            'bp',cell2mat(data.g(gloc).e(eloc).sd_id.raw),...
                            'wp',local_total_within(exp(cell2mat(data.g(gloc).e(eloc).b_sigma.raw)), cell2mat(data.g(gloc).e(eloc).sd_trl.raw)),'i',cell2mat(data.g(gloc).e(eloc).sd_trl.raw),...
                            'obs',[1 ntrials],...
                            'CI',ciperc);
                        
                        %find the number of trials to reach cutoffl
                        switch meascutoff
                            case 1
                                trlcutoff = find(llrel >= depcutoff, 1);
                            case 2
                                trlcutoff = find(mrel >= depcutoff, 1);
                            case 3
                                trlcutoff = find(ulrel >= depcutoff, 1);
                        end
                        
                        %if a cutoff was not found
                        if isempty(trlcutoff)
                            
                            %store all the participant ids as having bad data
                            relerr.trlcutoff = 1;
                            trlcutoff = -1;
                            relsummary.group(gloc).event(eloc).trlcutoff = -1;
                            relsummary.group(gloc).event(eloc).rel.m = -1;
                            relsummary.group(gloc).event(eloc).rel.ll = -1;
                            relsummary.group(gloc).event(eloc).rel.ul = -1;
                            
                            trltable = varfun(@length,REL.data,...
                                'GroupingVariables',{'group','id','event'});
                            
                            trltable = trltable(strcmp(trltable.event,enames{eloc})...
                                & strcmp(trltable.group,gnames{gloc}),:);
                            
                            relsummary.group(gloc).event(eloc).eventgoodids =...
                                cell(0,1);
                            relsummary.group(gloc).event(eloc).eventbadids =...
                                trltable.id;
                            
                        elseif ~isempty(trlcutoff) %if a cutoff was found
                            
                            
                            relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                            relsummary.group(gloc).event(eloc).rel.m = mrel(trlcutoff);
                            relsummary.group(gloc).event(eloc).rel.ll = llrel(trlcutoff);
                            relsummary.group(gloc).event(eloc).rel.ul = ulrel(trlcutoff);
                            
                            trltable = varfun(@length,REL.data,...
                                'GroupingVariables',{'group','id','event'});
                            
                            trltable = trltable(strcmp(trltable.event,enames{eloc})...
                                & strcmp(trltable.group,gnames{gloc}),:);
                            
                            ind2include = trltable.GroupCount >= trlcutoff;
                            ind2exclude = trltable.GroupCount < trlcutoff;
                            
                            relsummary.group(gloc).event(eloc).eventgoodids =...
                                trltable.id(ind2include);
                            relsummary.group(gloc).event(eloc).eventbadids =...
                                trltable.id(ind2exclude);
                            
                        end
                        
                    end
                end
                
                %find the ids that have enough trials for each event type
                
                for gloc=1:ngroups
                    tempids = {};
                    badids = [];
                    
                    for eloc=1:nevents
                        %STRICT all-events rule (G39 ruling 2026-08-29): an event that
                        %retained nobody folds in as an empty set and empties the running
                        %intersection instead of being skipped; its eventbadids (the full
                        %id list for a dead event, by producer convention) fold in too, so
                        %the group good/bad ids complement.
                        tempids{end+1} = relsummary.group(gloc).event(eloc).eventgoodids;
                        if eloc == 1
                            badids = relsummary.group(gloc).event(eloc).eventbadids;
                        elseif eloc > 1
                            [~,ind]=setdiff(tempids{1},tempids{end});
                            new = tempids{1};
                            badids = vertcat(badids,...
                                relsummary.group(gloc).event(eloc).eventbadids,...
                                new(ind));
                            new(ind) = [];
                            tempids{1} = new;
                        end
                    end
                    
                    %defensive: with the unconditional fold above, tempids is empty only when nevents == 0; keep the empty-cell carrier convention
                    if isempty(tempids)
                        tempids{1} = cell(0,1);
                    end
                    
                    relsummary.group(gloc).goodids = tempids{1};
                    relsummary.group(gloc).badids = unique(badids);
                    %STRICT-fold consequence (G39): a group can retain nobody while
                    %its per-event rows stay informative, so warn on the console and
                    %continue. Deliberately NOT local_relsummary_notice: its modal
                    %contract requires every caller to set relerr.nogooddata and
                    %return, and this is not an abort.
                    if ~(iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids))
                        fprintf(['\nNo participants met the reliability threshold ',...
                            'for every event type in a group; that group''s ',...
                            'n Included is 0.\n']);
                    end
                    
                end
                
                %calculate dependability estimates for the overall data for each
                %group and event type
                for eloc=1:nevents
                    
                    for gloc=1:ngroups
                        
                        %pull the data to calculate the trial information from
                        %participants with good data
                        datatable = REL.data;
                        ind = strcmp(datatable.group,gnames{gloc}) & ...
                           strcmp(datatable.event,enames{eloc});
                       
                        datasubset = datatable(ind,:);
                        
                        if iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids)
                            goodids = table(relsummary.group(gloc).goodids);
                        else
                            goodids = table(relsummary.group(gloc).badids);
                        end
                        
                        %check whether there are not enough good data after applying
                        %reliability threshold and extrapolating
                        if isempty(goodids)
                            dlg = {'Data do not reach reliability threshold after extrapolation';...
                                'Set a lower reliability threshold'};
                            local_relsummary_notice(interactive, dlg);
                            relerr.nogooddata = 1;
                            return;
                        end
                        
                        trlcdata = innerjoin(datasubset, goodids,...
                            'LeftKeys', 'id', 'RightKeys', 'Var1',...
                            'LeftVariables', {'id' 'meas'});
                        
                        trltable = varfun(@length,trlcdata,...
                            'GroupingVariables',{'id'});
                        
                        trlmean = mean(trltable.GroupCount);
                        trlmed = median(trltable.GroupCount);
                        
                        %calculate dependability using either the mean or median
                        %(based on user input)
                        switch depcentmeas
                            case 1
                                depcent = trlmean;
                            case 2
                                depcent = trlmed;
                        end
                        
                        [relsummary.group(gloc).event(eloc).dep.ll,...
                            relsummary.group(gloc).event(eloc).dep.m,...
                            relsummary.group(gloc).event(eloc).dep.ul] =...
                            psyrat_rel_sing('gcoeff',gcoeff,'metric','global','bp',cell2mat(data.g(gloc).e(eloc).sd_id.raw),...
                            'wp',local_total_within(exp(cell2mat(data.g(gloc).e(eloc).b_sigma.raw)), cell2mat(data.g(gloc).e(eloc).sd_trl.raw)),'i',cell2mat(data.g(gloc).e(eloc).sd_trl.raw),...
                            'obs',depcent,'CI',ciperc);
                        
                        relsummary.group(gloc).event(eloc).dep.meas = depcent;
                        
                        relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                        relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                        relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                        relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                        relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                        
                        %check if trial cutoff for an event exceeded the total
                        %number of trials for any subjects
                        if  relsummary.group(gloc).event(eloc).trlcutoff >...
                                relsummary.group(gloc).event(eloc).trlinfo.max
                            
                            relerr.trlmax = 1;
                            relsummary.group(gloc).event(eloc).rel.m = -1;
                            relsummary.group(gloc).event(eloc).rel.ll = -1;
                            relsummary.group(gloc).event(eloc).rel.ul = -1;
                            
                        end
                        
                        %goodn is the GROUP-level retained count (G39
                        %ruling 2026-08-29): keyed on the group carrier,
                        %exactly as the one-group branches key it, so a
                        %dead event's row and its siblings print the same
                        %complementing count. The event-carrier keying
                        %this replaces printed 0 on a dead event's row
                        %while its siblings printed the group figure, and
                        %under the strict fold it would have counted the
                        %badids fallback table instead.
                        if iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids)
                            relsummary.group(gloc).event(eloc).goodn = height(goodids);
                        else
                            relsummary.group(gloc).event(eloc).goodn = 0;
                        end
                        
                        %calculate iccs, between-subject standard deviations, and
                        %within-subject standard deviations
                        [relsummary.group(gloc).event(eloc).icc.ll,...
                            relsummary.group(gloc).event(eloc).icc.m,...
                            relsummary.group(gloc).event(eloc).icc.ul] = ...
                            psyrat_rel_sing('gcoeff',gcoeff,'metric','icc','bp',cell2mat(data.g(gloc).e(eloc).sd_id.raw),...
                            'wp',local_total_within(exp(cell2mat(data.g(gloc).e(eloc).b_sigma.raw)), cell2mat(data.g(gloc).e(eloc).sd_trl.raw)),'i',cell2mat(data.g(gloc).e(eloc).sd_trl.raw),...
                            'CI',ciperc);
                        
                        relsummary.group(gloc).event(eloc).betsd.m = ...
                            mean(cell2mat(data.g(gloc).e(eloc).sd_id.raw));
                        relsummary.group(gloc).event(eloc).betsd.ll = ...
                            quantile(cell2mat(data.g(gloc).e(eloc).sd_id.raw),ciedge);
                        relsummary.group(gloc).event(eloc).betsd.ul = ...
                            quantile(cell2mat(data.g(gloc).e(eloc).sd_id.raw),1-ciedge);
                        
                        relsummary.group(gloc).event(eloc).witsd.m = ...
                            mean(exp(cell2mat(data.g(gloc).e(eloc).b_sigma.raw)));
                        relsummary.group(gloc).event(eloc).witsd.ll = ...
                            quantile(exp(cell2mat(data.g(gloc).e(eloc).b_sigma.raw)),...
                            ciedge);
                        relsummary.group(gloc).event(eloc).witsd.ul = ...
                            quantile(exp(cell2mat(data.g(gloc).e(eloc).b_sigma.raw)),...
                            1-ciedge);
                        
                    end
                end
                
                
        end %switch analysis
        
        if ngroups == 0
            ngroups = 1;
        end

        if diffgcoeff == 1
            diff_est = 'dep';
        else
            diff_est = 'gen';
        end
        
        for gloc=1:ngroups
            
            switch depcentmeas
                case 1
                    trls1 = relsummary.group(gloc).event(1).trlinfo.mean;
                    trls2 = relsummary.group(gloc).event(2).trlinfo.mean;
                case 2
                    trls1 = relsummary.group(gloc).event(1).trlinfo.med;
                    trls2 = relsummary.group(gloc).event(2).trlinfo.med;
            end
            
            wp_cov = [];
            if isfield(REL,'diffwpcov') && REL.diffwpcov == 1
                wp_cov = zeros(size(cell2mat(REL.out.b_sigma{:,gloc}),1),1);
            elseif isfield(REL.out,'err_varcov') && ~isempty(REL.out.err_varcov) &&...
                    size(REL.out.err_varcov,2) >= gloc &&...
                    ~isempty(REL.out.err_varcov{1,gloc})
                err_varcov = cell2mat(REL.out.err_varcov{:,gloc});
                wp_cov = err_varcov(:,2,1);
            elseif isfield(REL.out,'rescor') && ~isempty(REL.out.rescor) &&...
                    size(REL.out.rescor,2) >= gloc &&...
                    ~isempty(REL.out.rescor{1,gloc})
                %fallback: reconstruct covariance from residual correlation
                %and event-specific residual SDs
                rescor = cell2mat(REL.out.rescor{:,gloc});
                bsig = cell2mat(REL.out.b_sigma{:,gloc});
                wp_cov = rescor(:) .* exp(bsig(:,1)) .* exp(bsig(:,2));
            end
            
            if isempty(wp_cov)
                %legacy fallback: residual covariance not estimated/stored
                bsig = cell2mat(REL.out.b_sigma{:,gloc});
                wp_cov = zeros(size(bsig,1),1);
            end
            
            diffscore = psyrat_diffrel(...
                'bp',cell2mat(REL.out.id_varcov{:,gloc}),...
                'bt',cell2mat(REL.out.trl_varcov{:,gloc}),...
                'er_var',cell2mat(REL.out.b_sigma{:,gloc}),...
                'wp_cov',wp_cov,...
                'obs',[trls1 trls2],...
                'est',diff_est,...
                'CI',relsummary.ciperc);
            
            relsummary.group(gloc).diffscore = diffscore;

            %Cutoff row for the difference score. This used to evaluate the
            %difference at the two INDEPENDENTLY derived per-event cutoffs --
            %each chosen so that its own event separately cleared the
            %threshold -- and then suppress the trial count in the table. A
            %user reading that row believed they were seeing "trials until the
            %difference score is reliable," which is not what it answered.
            %Derive a single common n' from the difference score itself
            %instead, and report it. Holding the two counts equal also removes
            %the out-of-range case from this row: at a common n' the harmonic,
            %geometric and arithmetic means coincide, so the error variance
            %cannot go negative for positive semi-definite input.
            %
            %Search range mirrors the per-event convention above: the largest
            %observed trial count, plus room to extrapolate beyond it.
            ntrlsearch = max([relsummary.group(gloc).event(1).trlinfo.max ...
                relsummary.group(gloc).event(2).trlinfo.max]) + 1000;

            [diffcutoff,diffcutscore] = psyrat_diffrel_trlcutoff(...
                'bp',cell2mat(REL.out.id_varcov{:,gloc}),...
                'bt',cell2mat(REL.out.trl_varcov{:,gloc}),...
                'er_var',cell2mat(REL.out.b_sigma{:,gloc}),...
                'wp_cov',wp_cov,...
                'est',diff_est,...
                'CI',relsummary.ciperc,...
                'depcutoff',depcutoff,...
                'meascutoff',meascutoff,...
                'ntrials',ntrlsearch);

            %no common trial count reaches the threshold: report the observed
            %trial difference score with the -1 "not reached" sentinel the
            %per-event rows already use. This mirrors the test-retest sibling
            %below, which falls back to the observed-trial difference score
            %when a per-event cutoff is unavailable, and it also keeps a -1
            %from ever reaching psyrat_harmmean, which rejects non-positive
            %trial counts.
            %
            %relerr.trlcutoff is raised alongside it, because every other
            %producer of the -1 sentinel in this file does so, and that flag
            %is what tells the user the cutoffs could not be calculated. The
            %row keeps the observed-trial estimate rather than being blanked
            %to -1 the way a per-event row is; that asymmetry follows the
            %test-retest sibling and is deliberate.
            if diffcutoff < 0
                diffcutscore = diffscore;
                if isfield(diffcutscore,'admissibility')
                    %G60 wave: this arm ALIASES the parent evaluation, so
                    %carrying its admissibility record forward would make
                    %psyrat_admissibility_collect gather the IDENTICAL
                    %record twice (its second gather exists for the found-
                    %cutoff arm's independent psyrat_diffrel_trlcutoff
                    %evaluation). Same defect class as the sing_diff_sserr
                    %self-copy fixed in the same batch.
                    diffcutscore = rmfield(diffcutscore,'admissibility');
                end
                relerr.trlcutoff = 1;

            elseif diffcutoff > min([relsummary.group(gloc).event(1).trlinfo.max ...
                    relsummary.group(gloc).event(2).trlinfo.max])

                %the derived count is a COMMON n', so it needs that many trials
                %in BOTH conditions. It is therefore an extrapolation as soon
                %as it passes the SMALLER of the two observed maxima, not the
                %larger: past that point no participant supplied n' trials in
                %both conditions. Comparing against the larger maximum would
                %leave the whole band between the two unflagged, which is most
                %of the realistic range for the unequal-trial designs this row
                %exists for.
                %
                %A per-event row compares its own cutoff against its own event's
                %maximum, so applying that convention to a count that governs
                %both events gives exactly this min(). Without the flag, this
                %would be the only row in the table that can print an
                %unreachable count while looking reachable.
                relerr.trlmax = 1;
            end

            %the derived count travels with the difference score so the cutoff
            %table can print it instead of suppressing it
            diffcutscore.ntrials = diffcutoff;

            relsummary.group(gloc).diffscore.trlcutoff = diffcutscore;

        end
        
        %here can loop through and add in difference score reliability
        %estimation
        
    case 'sing_diff_sserr'

        [relsummary,data,relerr] = psyrat_relsummary_sing_diff_sserr(...
            psyrat_data,depcutoff,meascutoff,depcentmeas,diffgcoeff,ciperc);

    case 'sing_sserr'
        
        %place psyrat_data.data in REL to work with
        REL = psyrat_data.rel;
        
        %create a data structure for storing outputs
        data = struct;
        
        %create variable to specify whether there are bad/unreliable data
        %default state: 0, will be changed to 1 if there is a problem
        relerr = struct();
        relerr.trlcutoff = 0;
        relerr.trlmax = 0;
        relerr.nogooddata = 0;
        
        %check whether any groups exist
        if strcmpi(REL.groups,'none')
            ngroups = 1;
            %gnames = cellstr(REL.groups);
            gnames ={''};
        else
            ngroups = length(REL.groups);
            gnames = REL.groups(:);
        end
        
        %check whether any events exist
        if strcmpi(REL.events,'none')
            nevents = 1;
            %enames = cellstr(REL.events);
            enames = {''};
        else
            nevents = length(REL.events);
            enames = REL.events(:);
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

        %sserr backward-compat sig_trl default (shared helper); the subject-level
        %model's pop_sdlog fixes the shape so the absolute and relative
        %coefficients coincide for legacy structs (mirrors the 'sing' default).
        REL.out = local_default_sigtrl(REL.out, REL.out.pop_sdlog);

        %extract information from REL and store in data for crunching
        switch analysis

            case {1,2} %2 - possible multiple groups but no event types to consider

                eloc = 1;
                
                for gloc=1:length(REL.out.labels)
                    
                    data.g(gloc).e(eloc).label = REL.out.labels(gloc);
                    data.g(gloc).e(eloc).mu = REL.out.mu(:,gloc);
                    data.g(gloc).e(eloc).ind_bs = REL.out.ind_bs{:,gloc};
                    data.g(gloc).e(eloc).gro_sds = REL.out.gro_sds{:,gloc};
                    data.g(gloc).e(eloc).pop_sdlog = REL.out.pop_sdlog(:,gloc);
                    data.g(gloc).e(eloc).ind_sdlog = REL.out.ind_sdlog{:,gloc};
                    data.g(gloc).e(eloc).sig_trl = REL.out.sig_trl(:,gloc);

                    data.g(gloc).e(eloc).elabel = cellstr('none');
                    data.g(gloc).e(eloc).id_match = REL.out.id_matches{:,gloc};
                    data.g(gloc).glabel = gnames(gloc);
                    
                end
                
            case 3 %3 - possible event types but no groups to consider
                
                gloc = 1;
                
                for eloc=1:length(REL.out.labels)
                    
                    data.g(gloc).e(eloc).label = REL.out.labels(eloc);
                    data.g(gloc).e(eloc).mu = REL.out.mu(:,eloc);
                    data.g(gloc).e(eloc).ind_bs = REL.out.ind_bs{:,eloc};
                    data.g(gloc).e(eloc).gro_sds = REL.out.gro_sds{:,eloc};
                    data.g(gloc).e(eloc).pop_sdlog = REL.out.pop_sdlog(:,eloc);
                    data.g(gloc).e(eloc).ind_sdlog = REL.out.ind_sdlog{:,eloc};
                    data.g(gloc).e(eloc).sig_trl = REL.out.sig_trl(:,eloc);

                    data.g(gloc).e(eloc).elabel = enames(eloc);
                    data.g(gloc).e(eloc).id_match = REL.out.id_matches{:,eloc};
                    data.g(gloc).glabel = cellstr('none');
                    
                end
                
            case 4 %4 - possible groups and event types to consider
                for i=1:length(REL.out.labels)
                    
                    %use the underscores that were added in psyrat_computevarcomp to
                    %differentiate where the group and event the data are for
                    lblstr = strsplit(REL.out.labels{i},'_;_');
                    
                    eloc = find(ismember(enames,lblstr(2)));
                    gloc = find(ismember(gnames,lblstr(1)));
                    
                    data.g(gloc).e(eloc).label = REL.out.labels{i};
                    data.g(gloc).e(eloc).mu = REL.out.mu(:,i);
                    data.g(gloc).e(eloc).ind_bs = REL.out.ind_bs{:,i};
                    data.g(gloc).e(eloc).gro_sds = REL.out.gro_sds{:,i};
                    data.g(gloc).e(eloc).pop_sdlog = REL.out.pop_sdlog(:,i);
                    data.g(gloc).e(eloc).ind_sdlog = REL.out.ind_sdlog{:,i};
                    data.g(gloc).e(eloc).sig_trl = REL.out.sig_trl(:,i);

                    data.g(gloc).e(eloc).elabel = enames(eloc);
                    data.g(gloc).e(eloc).id_match = REL.out.id_matches{:,i};
                    data.g(gloc).glabel = gnames(gloc);
                    
                end
        end %switch analysis
        
        %store the cutoff in relsummary to pass to other functions more easily
        relsummary.depcutoff = depcutoff;
        relsummary.ciperc = ciperc;

        %edge for this arm's variance-component quantiles; see the 'sing' arm
        %above for why the component SDs follow ciperc rather than a fixed 95%
        ciedge = (1-ciperc)./2;
        relsummary.gcoeff = gcoeff;
        if gcoeff == 1
            relsummary.gcoeff_name = 'dep';
        else
            relsummary.gcoeff_name = 'gen';
        end
        
        %store which measure was used to specify cutoff
        switch meascutoff
            case 1
                relsummary.meascutoff = strcat('Lower Limit of ',...
                    sprintf(' %2.0f',ciperc*100),'% Credible Interval');
            case 2
                relsummary.meascutoff = 'Point Estimate';
            case 3
                relsummary.meascutoff = strcat('Upper Limit of ',...
                    sprintf(' %2.0f',ciperc*100),'% Credible Interval');
        end
        
        %compute dependability data for each group and event
        switch analysis
            case 1 %no groups or event types to consider
                
                %the same generic structure is used for relsummary, so the event
                %and group locations will both be 1 for the data
                eloc = 1;
                gloc = 1;
                
                %grab the group names (if present)
                relsummary.group(gloc).name = gnames{gloc};
                
                %set the event name as measure
                relsummary.group(gloc).event(eloc).name = 'measure';
                
                try
                    %create empty arrays for storing dependability information
                    trltable = varfun(@length,REL.data{1},'GroupingVariables',{'id'});
                catch
                    %create empty arrays for storing dependability information
                    trltable = varfun(@length,REL.data,'GroupingVariables',{'id'});
                end
                
                %compute dependability
                idtable = innerjoin(data.g(gloc).e(eloc).id_match,...
                    trltable(:,1:2));
                idtable.Properties.VariableNames{3} = 'trls';
                
                ssrel_table = psyrat_ssrel(...
                    'bp',data.g(gloc).e(eloc).gro_sds(:,1),...
                    'wp_pop',data.g(gloc).e(eloc).pop_sdlog,...
                    'wp_ss',data.g(gloc).e(eloc).ind_sdlog,...
                    'i',data.g(gloc).e(eloc).sig_trl,...
                    'gcoeff',gcoeff,...
                    'idtable',idtable,...
                    'CI',ciperc);
                
                
                
                %determine which participants had adequate reliability.
                %G49: inclusion and exclusion used to come from independent
                %>=/< comparisons, which a NaN dependability draw fails BOTH
                %of, so a partial-NaN subject landed in NEITHER carrier and
                %the good/bad complement broke. The shared marker derives
                %exclusion as ~inclusion, so NaN reads as excluded.
                ssrel_table = psyrat_relsummary_mark_inclusion(...
                    ssrel_table,depcutoff,meascutoff);
                
                %store summary table
                relsummary.group(gloc).event(eloc).ssrel_table = ssrel_table;
                
                if all(ssrel_table.ind2exclude)
                    %nobody met the cutoff -- including the all-NaN
                    %failed-fit chunk, which G49's complement (exclusion =
                    %~inclusion) folds into this arm -- so the whole roster
                    %is bad
                    relsummary.group(gloc).event(eloc).eventgoodids = cell(0,1);
                    relsummary.group(gloc).event(eloc).eventbadids =...
                        ssrel_table.id(:);
                    
                    relsummary.group(gloc).goodids = cell(0,1);
                    relsummary.group(gloc).badids = ...
                        ssrel_table.id(:);
                    
                    %get a little summary information about the number of
                    %trials
                    relsummary.group(gloc).event(eloc).trlinfo.min = min(ssrel_table.trls);
                    relsummary.group(gloc).event(eloc).trlinfo.max = max(ssrel_table.trls);
                    relsummary.group(gloc).event(eloc).trlinfo.mean = mean(ssrel_table.trls);
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(ssrel_table.trls);
                    relsummary.group(gloc).event(eloc).trlinfo.med = median(ssrel_table.trls);
                    relsummary.group(gloc).event(eloc).goodn = 0;
                    
                elseif any(ssrel_table.ind2include)

                    %partition by the indicator: retained subjects to the
                    %good carriers, everyone else -- including any
                    %partial-NaN subject, whom G49's complement marks
                    %excluded -- to the bad carriers
                    relsummary.group(gloc).event(eloc).eventgoodids =...
                        ssrel_table.id(ssrel_table.ind2include);
                    relsummary.group(gloc).event(eloc).eventbadids =...
                        ssrel_table.id(ssrel_table.ind2exclude);
                    
                    relsummary.group(gloc).goodids = ...
                        ssrel_table.id(ssrel_table.ind2include);
                    relsummary.group(gloc).badids = ...
                        ssrel_table.id(ssrel_table.ind2exclude);
                    
                    %get a little summary information about the number of
                    %trials
                    relsummary.group(gloc).event(eloc).trlinfo.min = ...
                        min(ssrel_table.trls(ssrel_table.ind2include));
                    relsummary.group(gloc).event(eloc).trlinfo.max = ...
                        max(ssrel_table.trls(ssrel_table.ind2include));
                    relsummary.group(gloc).event(eloc).trlinfo.mean = ...
                        mean(ssrel_table.trls(ssrel_table.ind2include));
                    relsummary.group(gloc).event(eloc).trlinfo.std = ...
                        std(ssrel_table.trls(ssrel_table.ind2include));
                    relsummary.group(gloc).event(eloc).trlinfo.med = ...
                        median(ssrel_table.trls(ssrel_table.ind2include));
                    relsummary.group(gloc).event(eloc).goodn = ...
                        height(ssrel_table(ssrel_table.ind2include,:));
                    
                else
                    %G46, now a backstop. Since G49 derives ind2exclude as
                    %~ind2include, an all-NaN chunk (nobody includable) makes
                    %all(ind2exclude) true and takes the first arm, whose
                    %assignments are identical -- this arm is unreachable
                    %unless the indicator convention diverges again. Kept as
                    %the documented dead-event producer for that state:
                    %canonical pair, nobody retained, full roster bad.
                    relsummary.group(gloc).event(eloc).eventgoodids = cell(0,1);
                    relsummary.group(gloc).event(eloc).eventbadids =...
                        ssrel_table.id(:);

                    relsummary.group(gloc).goodids = cell(0,1);
                    relsummary.group(gloc).badids = ...
                        ssrel_table.id(:);

                    %get a little summary information about the number of
                    %trials
                    relsummary.group(gloc).event(eloc).trlinfo.min = min(ssrel_table.trls);
                    relsummary.group(gloc).event(eloc).trlinfo.max = max(ssrel_table.trls);
                    relsummary.group(gloc).event(eloc).trlinfo.mean = mean(ssrel_table.trls);
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(ssrel_table.trls);
                    relsummary.group(gloc).event(eloc).trlinfo.med = median(ssrel_table.trls);
                    relsummary.group(gloc).event(eloc).goodn = 0;
                end
                
                
            case 2 %possible multiple groups but no event types to consider
                
                %since the same generic structure is used for relsummary, the event
                %location will be defined as 1.
                eloc = 1;
                
                for gloc=1:ngroups %loop through each group
                    
                    %store the name of the group
                    if eloc == 1
                        relsummary.group(gloc).name = gnames{gloc};
                    end
                    
                    %store the name of the event
                    relsummary.group(gloc).event(eloc).name = enames{eloc};
                    
                    
                    try
                        %create empty arrays for storing dependability information
                        trltable = varfun(@length,REL.data{1},'GroupingVariables',{'group','id'});
                    catch
                        %create empty arrays for storing dependability information
                        trltable = varfun(@length,REL.data,'GroupingVariables',{'group','id'});
                    end

                    %restrict to the current group; this case only runs when
                    %ngroups>1 (analysis==2), so the 'group' column always exists
                    trltable = trltable(strcmp(trltable.group,gnames{gloc}),:);

                    %compute dependability
                    idtable = innerjoin(data.g(gloc).e(eloc).id_match,...
                        trltable(:,{'id','GroupCount'}));
                    idtable.Properties.VariableNames{3} = 'trls';

                    ssrel_table = psyrat_ssrel(...
                        'bp',data.g(gloc).e(eloc).gro_sds(:,1),...
                        'wp_pop',data.g(gloc).e(eloc).pop_sdlog,...
                        'wp_ss',data.g(gloc).e(eloc).ind_sdlog,...
                        'i',data.g(gloc).e(eloc).sig_trl,...
                        'gcoeff',gcoeff,...
                        'idtable',idtable,...
                        'CI',ciperc);
                    
                    
                    
                    %determine which participants had adequate reliability.
                    %G49: inclusion and exclusion used to come from independent
                    %>=/< comparisons, which a NaN dependability draw fails BOTH
                    %of, so a partial-NaN subject landed in NEITHER carrier and
                    %the good/bad complement broke. The shared marker derives
                    %exclusion as ~inclusion, so NaN reads as excluded.
                    ssrel_table = psyrat_relsummary_mark_inclusion(...
                        ssrel_table,depcutoff,meascutoff);
                    
                    %store summary table
                    relsummary.group(gloc).event(eloc).ssrel_table = ssrel_table;
                    
                    if all(ssrel_table.ind2exclude)
                        %nobody met the cutoff -- including the all-NaN
                        %failed-fit chunk, which G49's complement (exclusion
                        %= ~inclusion) folds into this arm -- so the whole
                        %roster is bad
                        relsummary.group(gloc).event(eloc).eventgoodids = cell(0,1);
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            ssrel_table.id(:);
                        
                        relsummary.group(gloc).goodids = cell(0,1);
                        relsummary.group(gloc).badids = ...
                            ssrel_table.id(:);
                        
                        %get a little summary information about the number of
                        %trials
                        relsummary.group(gloc).event(eloc).trlinfo.min = min(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.max = max(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.mean = mean(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.std = std(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.med = median(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).goodn = 0;
                        
                    elseif any(ssrel_table.ind2include)

                        %partition by the indicator: retained subjects to
                        %the good carriers, everyone else -- including any
                        %partial-NaN subject, whom G49's complement marks
                        %excluded -- to the bad carriers
                        relsummary.group(gloc).event(eloc).eventgoodids =...
                            ssrel_table.id(ssrel_table.ind2include);
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            ssrel_table.id(ssrel_table.ind2exclude);
                        
                        relsummary.group(gloc).goodids = ...
                            ssrel_table.id(ssrel_table.ind2include);
                        relsummary.group(gloc).badids = ...
                            ssrel_table.id(ssrel_table.ind2exclude);
                        
                        %get a little summary information about the number of
                        %trials
                        relsummary.group(gloc).event(eloc).trlinfo.min = ...
                            min(ssrel_table.trls(ssrel_table.ind2include));
                        relsummary.group(gloc).event(eloc).trlinfo.max = ...
                            max(ssrel_table.trls(ssrel_table.ind2include));
                        relsummary.group(gloc).event(eloc).trlinfo.mean = ...
                            mean(ssrel_table.trls(ssrel_table.ind2include));
                        relsummary.group(gloc).event(eloc).trlinfo.std = ...
                            std(ssrel_table.trls(ssrel_table.ind2include));
                        relsummary.group(gloc).event(eloc).trlinfo.med = ...
                            median(ssrel_table.trls(ssrel_table.ind2include));
                        relsummary.group(gloc).event(eloc).goodn = ...
                            height(ssrel_table(ssrel_table.ind2include,:));
                        
                    else
                        %G46, now a backstop. Since G49 derives ind2exclude as
                        %~ind2include, an all-NaN chunk (nobody includable) makes
                        %all(ind2exclude) true and takes the first arm, whose
                        %assignments are identical -- this arm is unreachable
                        %unless the indicator convention diverges again. Kept as
                        %the documented dead-event producer for that state:
                        %canonical pair, nobody retained, full roster bad.
                        relsummary.group(gloc).event(eloc).eventgoodids = cell(0,1);
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            ssrel_table.id(:);

                        relsummary.group(gloc).goodids = cell(0,1);
                        relsummary.group(gloc).badids = ...
                            ssrel_table.id(:);

                        %get a little summary information about the number of
                        %trials
                        relsummary.group(gloc).event(eloc).trlinfo.min = min(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.max = max(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.mean = mean(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.std = std(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.med = median(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).goodn = 0;
                    end
                end
                
                
            case 3 %possible event types but no groups to consider
                
                %since the relsummary structure is generic for any number of groups
                %or events, the data will count as being from 1 group.
                gloc = 1;
                
                %cycle through each event
                for eloc=1:nevents
                    
                    %get the group name
                    if eloc == 1
                        relsummary.group(gloc).name = gnames{gloc};
                    end
                    
                    %get the event name
                    relsummary.group(gloc).event(eloc).name = enames{eloc};
                    
                    try
                        %create empty arrays for storing dependability information
                        trltable = varfun(@length,REL.data{1},'GroupingVariables',{'id','event'});
                    catch
                        %create empty arrays for storing dependability information
                        trltable = varfun(@length,REL.data,'GroupingVariables',{'id','event'});
                    end
                    
                    trltable = trltable(strcmp(trltable.event,enames{eloc}),:);
                    
                    %compute dependability
                    idtable = innerjoin(data.g(gloc).e(eloc).id_match,...
                        trltable(:,[1 3]));
                    idtable.Properties.VariableNames{3} = 'trls';

                    ssrel_table = psyrat_ssrel(...
                        'bp',data.g(gloc).e(eloc).gro_sds(:,1),...
                        'wp_pop',data.g(gloc).e(eloc).pop_sdlog,...
                        'wp_ss',data.g(gloc).e(eloc).ind_sdlog,...
                        'i',data.g(gloc).e(eloc).sig_trl,...
                        'gcoeff',gcoeff,...
                        'idtable',idtable,...
                        'CI',ciperc);
                    
                    %determine which participants had adequate reliability.
                    %G49: inclusion and exclusion used to come from independent
                    %>=/< comparisons, which a NaN dependability draw fails BOTH
                    %of, so a partial-NaN subject landed in NEITHER carrier and
                    %the good/bad complement broke. The shared marker derives
                    %exclusion as ~inclusion, so NaN reads as excluded.
                    ssrel_table = psyrat_relsummary_mark_inclusion(...
                        ssrel_table,depcutoff,meascutoff);
                    
                    %store summary table
                    relsummary.group(gloc).event(eloc).ssrel_table = ssrel_table;
                    
                    
                    if all(ssrel_table.ind2exclude)
                        %nobody met the cutoff -- including the all-NaN
                        %failed-fit chunk, which G49's complement (exclusion
                        %= ~inclusion) folds into this arm -- so the whole
                        %roster is bad
                        relsummary.group(gloc).event(eloc).eventgoodids = cell(0,1);
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            ssrel_table.id(:);
                        
                        %get a little summary information about the number of
                        %trials
                        relsummary.group(gloc).event(eloc).trlinfo.min = min(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.max = max(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.mean = mean(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.std = std(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.med = median(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).goodn = 0;
                        
                    elseif any(ssrel_table.ind2include)

                        %partition by the indicator: retained subjects to
                        %the good carriers, everyone else -- including any
                        %partial-NaN subject, whom G49's complement marks
                        %excluded -- to the bad carriers
                        relsummary.group(gloc).event(eloc).eventgoodids =...
                            ssrel_table.id(ssrel_table.ind2include);
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            ssrel_table.id(ssrel_table.ind2exclude);
                        
                        %get a little summary information about the number of
                        %trials
                        relsummary.group(gloc).event(eloc).trlinfo.min = ...
                            min(ssrel_table.trls(ssrel_table.ind2include));
                        relsummary.group(gloc).event(eloc).trlinfo.max = ...
                            max(ssrel_table.trls(ssrel_table.ind2include));
                        relsummary.group(gloc).event(eloc).trlinfo.mean = ...
                            mean(ssrel_table.trls(ssrel_table.ind2include));
                        relsummary.group(gloc).event(eloc).trlinfo.std = ...
                            std(ssrel_table.trls(ssrel_table.ind2include));
                        relsummary.group(gloc).event(eloc).trlinfo.med = ...
                            median(ssrel_table.trls(ssrel_table.ind2include));
                        relsummary.group(gloc).event(eloc).goodn = ...
                            height(ssrel_table(ssrel_table.ind2include,:));
                        
                    else
                        %G46, now a backstop. Since G49 derives ind2exclude as
                        %~ind2include, an all-NaN chunk (nobody includable) makes
                        %all(ind2exclude) true and takes the first arm, whose
                        %assignments are identical -- this arm is unreachable
                        %unless the indicator convention diverges again. Kept as
                        %the documented dead-event producer for that state:
                        %canonical pair, nobody retained, full roster bad.
                        relsummary.group(gloc).event(eloc).eventgoodids = cell(0,1);
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            ssrel_table.id(:);

                        %get a little summary information about the number of
                        %trials
                        relsummary.group(gloc).event(eloc).trlinfo.min = min(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.max = max(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.mean = mean(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.std = std(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).trlinfo.med = median(ssrel_table.trls);
                        relsummary.group(gloc).event(eloc).goodn = 0;
                    end
                    
                end
                
                %cycle through the events and store information about which
                %participants have good data across all event types.
                tempids = {};
                badids = [];
                for eloc=1:nevents
                    %STRICT all-events rule (G39 ruling 2026-08-29): an event that
                    %retained nobody folds in as an empty set and empties the running
                    %intersection instead of being skipped; its eventbadids (the full
                    %id list for a dead event, by producer convention) fold in too, so
                    %the group good/bad ids complement.
                    tempids{end+1} = relsummary.group(gloc).event(eloc).eventgoodids;
                    if eloc == 1
                        badids = relsummary.group(gloc).event(eloc).eventbadids;
                    elseif eloc > 1
                        [~,ind]=setdiff(tempids{1},tempids{end});
                        new = tempids{1};
                        badids = vertcat(badids,...
                            relsummary.group(gloc).event(eloc).eventbadids,...
                            new(ind));
                        new(ind) = [];
                        tempids{1} = new;
                    end
                end
                
                %defensive: with the unconditional fold above, tempids is empty only when nevents == 0; keep the empty-cell carrier convention
                if isempty(tempids)
                    tempids{1} = cell(0,1);
                end
                
                %store information about participants with good/bad ids
                relsummary.group(gloc).goodids = tempids{1};
                relsummary.group(gloc).badids = unique(badids);
                %STRICT-fold consequence (G39): a group can retain nobody while
                %its per-event rows stay informative, so warn on the console and
                %continue. Deliberately NOT local_relsummary_notice: its modal
                %contract requires every caller to set relerr.nogooddata and
                %return, and this is not an abort.
                if ~(iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids))
                    fprintf(['\nNo participants met the reliability threshold ',...
                        'for every event type in a group; that group''s ',...
                        'n Included is 0.\n']);
                end
                
                
            case 4 %groups and event types to consider
                
                %cycle through each event
                for eloc=1:nevents
                    
                    %cycle through each group
                    for gloc=1:ngroups
                        
                        %get the group name
                        if eloc == 1
                            relsummary.group(gloc).name = gnames{gloc};
                        end
                        
                        %get the event name
                        relsummary.group(gloc).event(eloc).name = enames{eloc};
                        
                        try
                            %create empty arrays for storing dependability information
                            trltable = varfun(@length,REL.data{1},'GroupingVariables',{'group','id','event'});
                        catch
                            %create empty arrays for storing dependability information
                            trltable = varfun(@length,REL.data,'GroupingVariables',{'group','id','event'});
                        end

                        trltable = trltable(strcmp(trltable.event,enames{eloc})...
                            & strcmp(trltable.group,gnames{gloc}),:);

                        %compute dependability
                        idtable = innerjoin(data.g(gloc).e(eloc).id_match,...
                            trltable(:,{'id','GroupCount'}));
                        idtable.Properties.VariableNames{3} = 'trls';
                        
                        ssrel_table = psyrat_ssrel(...
                            'bp',data.g(gloc).e(eloc).gro_sds(:,1),...
                            'wp_pop',data.g(gloc).e(eloc).pop_sdlog,...
                            'wp_ss',data.g(gloc).e(eloc).ind_sdlog,...
                            'i',data.g(gloc).e(eloc).sig_trl,...
                            'gcoeff',gcoeff,...
                            'idtable',idtable,...
                            'CI',ciperc);
                        
                        %determine which participants had adequate reliability.
                        %G49: inclusion and exclusion used to come from independent
                        %>=/< comparisons, which a NaN dependability draw fails BOTH
                        %of, so a partial-NaN subject landed in NEITHER carrier and
                        %the good/bad complement broke. The shared marker derives
                        %exclusion as ~inclusion, so NaN reads as excluded.
                        ssrel_table = psyrat_relsummary_mark_inclusion(...
                            ssrel_table,depcutoff,meascutoff);
                        
                        %store summary table
                        relsummary.group(gloc).event(eloc).ssrel_table = ssrel_table;
                        
                        
                        if all(ssrel_table.ind2exclude)
                            %nobody met the cutoff -- including the all-NaN
                            %failed-fit chunk, which G49's complement
                            %(exclusion = ~inclusion) folds into this arm --
                            %so the whole roster is bad
                            relsummary.group(gloc).event(eloc).eventgoodids = cell(0,1);
                            relsummary.group(gloc).event(eloc).eventbadids =...
                                ssrel_table.id(:);
                            
                            %get a little summary information about the number of
                            %trials
                            relsummary.group(gloc).event(eloc).trlinfo.min = min(ssrel_table.trls);
                            relsummary.group(gloc).event(eloc).trlinfo.max = max(ssrel_table.trls);
                            relsummary.group(gloc).event(eloc).trlinfo.mean = mean(ssrel_table.trls);
                            relsummary.group(gloc).event(eloc).trlinfo.std = std(ssrel_table.trls);
                            relsummary.group(gloc).event(eloc).trlinfo.med = median(ssrel_table.trls);
                            relsummary.group(gloc).event(eloc).goodn = 0;
                            
                        elseif any(ssrel_table.ind2include)

                            %partition by the indicator: retained subjects
                            %to the good carriers, everyone else --
                            %including any partial-NaN subject, whom G49's
                            %complement marks excluded -- to the bad
                            %carriers
                            relsummary.group(gloc).event(eloc).eventgoodids =...
                                ssrel_table.id(ssrel_table.ind2include);
                            relsummary.group(gloc).event(eloc).eventbadids =...
                                ssrel_table.id(ssrel_table.ind2exclude);
                            
                            %get a little summary information about the number of
                            %trials
                            relsummary.group(gloc).event(eloc).trlinfo.min = ...
                                min(ssrel_table.trls(ssrel_table.ind2include));
                            relsummary.group(gloc).event(eloc).trlinfo.max = ...
                                max(ssrel_table.trls(ssrel_table.ind2include));
                            relsummary.group(gloc).event(eloc).trlinfo.mean = ...
                                mean(ssrel_table.trls(ssrel_table.ind2include));
                            relsummary.group(gloc).event(eloc).trlinfo.std = ...
                                std(ssrel_table.trls(ssrel_table.ind2include));
                            relsummary.group(gloc).event(eloc).trlinfo.med = ...
                                median(ssrel_table.trls(ssrel_table.ind2include));
                            relsummary.group(gloc).event(eloc).goodn = ...
                                height(ssrel_table(ssrel_table.ind2include,:));
                            
                        else
                            %G46, now a backstop. Since G49 derives ind2exclude as
                            %~ind2include, an all-NaN chunk (nobody includable) makes
                            %all(ind2exclude) true and takes the first arm, whose
                            %assignments are identical -- this arm is unreachable
                            %unless the indicator convention diverges again. Kept as
                            %the documented dead-event producer for that state:
                            %canonical pair, nobody retained, full roster bad.
                            relsummary.group(gloc).event(eloc).eventgoodids = cell(0,1);
                            relsummary.group(gloc).event(eloc).eventbadids =...
                                ssrel_table.id(:);

                            %get a little summary information about the number of
                            %trials
                            relsummary.group(gloc).event(eloc).trlinfo.min = min(ssrel_table.trls);
                            relsummary.group(gloc).event(eloc).trlinfo.max = max(ssrel_table.trls);
                            relsummary.group(gloc).event(eloc).trlinfo.mean = mean(ssrel_table.trls);
                            relsummary.group(gloc).event(eloc).trlinfo.std = std(ssrel_table.trls);
                            relsummary.group(gloc).event(eloc).trlinfo.med = median(ssrel_table.trls);
                            relsummary.group(gloc).event(eloc).goodn = 0;
                        end
                        
                    end
                end
                
                %find the ids that have enough trials for each event type
                
                for gloc=1:ngroups
                    tempids = {};
                    badids = [];
                    
                    for eloc=1:nevents
                        %STRICT all-events rule (G39 ruling 2026-08-29): an event that
                        %retained nobody folds in as an empty set and empties the running
                        %intersection instead of being skipped; its eventbadids (the full
                        %id list for a dead event, by producer convention) fold in too, so
                        %the group good/bad ids complement.
                        tempids{end+1} = relsummary.group(gloc).event(eloc).eventgoodids;
                        if eloc == 1
                            badids = relsummary.group(gloc).event(eloc).eventbadids;
                        elseif eloc > 1
                            [~,ind]=setdiff(tempids{1},tempids{end});
                            new = tempids{1};
                            badids = vertcat(badids,...
                                relsummary.group(gloc).event(eloc).eventbadids,...
                                new(ind));
                            new(ind) = [];
                            tempids{1} = new;
                        end
                    end
                    
                    %defensive: with the unconditional fold above, tempids is empty only when nevents == 0; keep the empty-cell carrier convention
                    if isempty(tempids)
                        tempids{1} = cell(0,1);
                    end
                    
                    relsummary.group(gloc).goodids = tempids{1};
                    relsummary.group(gloc).badids = unique(badids);
                    %STRICT-fold consequence (G39): a group can retain nobody while
                    %its per-event rows stay informative, so warn on the console and
                    %continue. Deliberately NOT local_relsummary_notice: its modal
                    %contract requires every caller to set relerr.nogooddata and
                    %return, and this is not an abort.
                    if ~(iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids))
                        fprintf(['\nNo participants met the reliability threshold ',...
                            'for every event type in a group; that group''s ',...
                            'n Included is 0.\n']);
                    end
                    
                end
                
                
        end %switch analysis
        
    case 'trt_sserr'
        %Subject-level (per-participant) test-retest reliability (trt_sserrvar,
        %analysis 25): the test-retest analog of the single-occasion 'sing_sserr'
        %case. Each participant gets their own residual (from the location-scale
        %person block, exp(pop_sdlog + ind_sdlog)); the five crossed facet SDs are
        %population. Per group and event, build an idtable (with per-participant
        %trial counts) and call psyrat_ssrel_trt for the per-subject coefficient.
        REL = psyrat_data.rel;
        data = struct;
        relerr = struct();
        relerr.trlcutoff = 0;
        relerr.trlmax = 0;
        relerr.nogooddata = 0;

        %check whether any groups exist
        if strcmpi(REL.groups,'none')
            ngroups = 1;
            gnames = {''};
        else
            ngroups = length(REL.groups);
            gnames = REL.groups(:);
        end

        %check whether any events exist
        if strcmpi(REL.events,'none')
            nevents = 1;
            enames = {''};
        else
            nevents = length(REL.events);
            enames = REL.events(:);
        end

        %map each REL.out column (a group/event data chunk) to its (gloc,eloc)
        %and unpack the person block plus the five crossed facet SDs. The person
        %intercept (gro_sds/ind_bs) and log-residual (pop_sdlog/ind_sdlog) are the
        %sserr person block; the facet SDs are the crossed test-retest components.
        for icol = 1:length(REL.out.labels)
            if ngroups > 1 && nevents > 1
                lblstr = strsplit(REL.out.labels{icol},'_;_');
                gloc = find(ismember(gnames,lblstr(1)));
                eloc = find(ismember(enames,lblstr(2)));
            elseif ngroups > 1
                gloc = icol; eloc = 1;
            elseif nevents > 1
                gloc = 1; eloc = icol;
            else
                gloc = 1; eloc = 1;
            end

            data.g(gloc).e(eloc).label = REL.out.labels{icol};
            data.g(gloc).e(eloc).gro_sds = REL.out.gro_sds{:,icol};
            data.g(gloc).e(eloc).pop_sdlog = REL.out.pop_sdlog(:,icol);
            data.g(gloc).e(eloc).ind_sdlog = REL.out.ind_sdlog{:,icol};
            data.g(gloc).e(eloc).sig_occ = REL.out.sig_occ(:,icol);
            data.g(gloc).e(eloc).sig_trl = REL.out.sig_trl(:,icol);
            data.g(gloc).e(eloc).sig_trlxid = REL.out.sig_trlxid(:,icol);
            data.g(gloc).e(eloc).sig_occxid = REL.out.sig_occxid(:,icol);
            data.g(gloc).e(eloc).sig_trlxocc = REL.out.sig_trlxocc(:,icol);
            data.g(gloc).e(eloc).id_match = REL.out.id_matches{:,icol};
        end

        %resolve the test-retest occasion estimand (n'_o) now REL is known, then
        %store the run-level metadata (mirrors the 'trt' and 'sing_sserr' cases).
        [effnocc, nocc_name, ~] = local_resolve_effnocc(noccmode, noccarg, REL);
        relsummary.noccmode = noccmode;
        relsummary.nocc = effnocc;
        relsummary.nocc_name = nocc_name;

        relsummary.depcutoff = depcutoff;
        relsummary.ciperc = ciperc;

        %edge for this arm's variance-component quantiles; see the 'sing' arm
        %above for why the component SDs follow ciperc rather than a fixed 95%
        ciedge = (1-ciperc)./2;
        relsummary.gcoeff = gcoeff;
        if gcoeff == 1
            relsummary.gcoeff_name = 'dep';
        else
            relsummary.gcoeff_name = 'gen';
        end

        %RC-26. The per-participant coefficient below is computed FROM reltype
        %(see the local_ssrel_trt_call further down), so a run is one of CE/CS/
        %CES and its exported tables must be able to say which. reltype was
        %previously consumed here but never stored, so psyrat_trt_coeflabel --
        %which gates on the presence of relsummary.reltype -- emitted no
        %coefficient line for this design, and its exports named n'_o and the
        %error type but not the estimand. Stored with reltype_name for
        %consistency with the 'trt' and 'trt_diff' branches; note that the
        %NUMERIC field is the one label code reads.
        relsummary.reltype = reltype;
        if reltype == 1
            relsummary.reltype_name = 'ic';
        elseif reltype == 2
            relsummary.reltype_name = 'trt';
        elseif reltype == 3
            relsummary.reltype_name = 'ic_trt';
        end
        switch meascutoff
            case 1
                relsummary.meascutoff = strcat('Lower Limit of ',...
                    sprintf(' %2.0f',ciperc*100),'% Credible Interval');
            case 2
                relsummary.meascutoff = 'Point Estimate';
            case 3
                relsummary.meascutoff = strcat('Upper Limit of ',...
                    sprintf(' %2.0f',ciperc*100),'% Credible Interval');
        end

        %compute the per-subject test-retest coefficient for each group and event
        for gloc = 1:ngroups
            relsummary.group(gloc).name = gnames{gloc};
            for eloc = 1:nevents
                if nevents == 1
                    relsummary.group(gloc).event(eloc).name = 'measure';
                else
                    relsummary.group(gloc).event(eloc).name = enames{eloc};
                end

                %per-participant trial counts. The test-retest coefficient uses
                %trials-per-OCCASION (psyrat_ssrel_trt divides by obs and by
                %obs*nocc), so count rows per (id, occasion) and take each
                %participant's mean across occasions - the same trials-per-occasion
                %convention (GroupingVariables include 'time') as the group-level
                %'trt' branch, not the single-occasion total-rows count. This loop
                %runs even when ngroups==1 (no real group facet, placeholder gloc=1),
                %so 'group' is included only when it actually exists as a column.
                grpvars = {'id','time'};
                if ngroups > 1
                    grpvars = [{'group'}, grpvars];
                end
                if nevents > 1
                    grpvars = [grpvars, {'event'}];
                end
                try
                    trltable = varfun(@length,REL.data{1},...
                        'GroupingVariables',grpvars);
                catch
                    trltable = varfun(@length,REL.data,...
                        'GroupingVariables',grpvars);
                end
                if ngroups > 1
                    trltable = trltable(strcmp(trltable.group,gnames{gloc}),:);
                end
                if nevents > 1
                    trltable = trltable(strcmp(trltable.event,enames{eloc}),:);
                end
                %mean trials-per-occasion for each participant
                trlbyid = varfun(@mean,trltable,'GroupingVariables','id',...
                    'InputVariables','GroupCount');
                trlbyid.trls = trlbyid.mean_GroupCount;
                idtable = innerjoin(data.g(gloc).e(eloc).id_match,...
                    trlbyid(:,{'id','trls'}));

                ssrel_table = local_ssrel_trt_call(data.g(gloc).e(eloc),...
                    gcoeff,reltype,idtable,ciperc,effnocc);
                relsummary = local_store_ss_chunk(relsummary,ssrel_table,...
                    gloc,eloc,meascutoff,depcutoff);
            end
        end

        %group-level good/bad ids: participants meeting the cutoff across ALL
        %events in the group (mirrors the sing_sserr cross-event logic).
        for gloc = 1:ngroups
            tempids = {};
            badids = [];
            for eloc = 1:nevents
                %STRICT all-events rule (G39 ruling 2026-08-29): an event that
                %retained nobody folds in as an empty set and empties the running
                %intersection instead of being skipped; its eventbadids (the full
                %id list for a dead event, by producer convention) fold in too, so
                %the group good/bad ids complement.
                eg = relsummary.group(gloc).event(eloc).eventgoodids;
                tempids{end+1} = eg; %#ok<AGROW>
                if numel(tempids) == 1
                    badids = relsummary.group(gloc).event(eloc).eventbadids;
                else
                    [~,ind] = setdiff(tempids{1},tempids{end});
                    new = tempids{1};
                    badids = vertcat(badids,...
                        relsummary.group(gloc).event(eloc).eventbadids,...
                        new(ind));
                    new(ind) = [];
                    tempids{1} = new;
                end
            end
            if isempty(tempids)
                tempids{1} = cell(0,1);
            end
            relsummary.group(gloc).goodids = tempids{1};
            relsummary.group(gloc).badids = unique(badids);
            %STRICT-fold consequence (G39): a group can retain nobody while
            %its per-event rows stay informative, so warn on the console and
            %continue. Deliberately NOT local_relsummary_notice: its modal
            %contract requires every caller to set relerr.nogooddata and
            %return, and this is not an abort.
            if ~(iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids))
                fprintf(['\nNo participants met the reliability threshold ',...
                    'for every event type in a group; that group''s ',...
                    'n Included is 0.\n']);
            end
        end

    case 'trt'
        %place psyrat_data.data in REL to work with
        REL = psyrat_data.rel;

        %create a data structure for storing outputs
        data = struct;

        %create variable to specify whether there are bad/unreliable data
        %default state: 0, will be changed to 1 if there is a problem
        relerr = struct();
        relerr.trlcutoff = 0;
        relerr.trlmax = 0;
        relerr.nogooddata = 0;

        %check whether any groups exist
        if strcmpi(REL.groups,'none')
            ngroups = 1;
            %gnames = cellstr(REL.groups);
            gnames ={''};
        else
            ngroups = length(REL.groups);
            gnames = REL.groups(:);
        end

        %check whether any events exist
        if strcmpi(REL.events,'none')
            nevents = 1;
            %enames = cellstr(REL.events);
            enames = {''};
        else
            nevents = length(REL.events);
            enames = REL.events(:);
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

        %extract information from REL and store in data for crunching
        switch analysis
            case 1 %1 - no groups or event types to consider

                gloc = 1;
                eloc = 1;

                data.g(gloc).e(eloc).label = REL.out.labels{gloc};
                data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,gloc);
                data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,gloc);
                data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,gloc);
                data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,gloc);
                data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,gloc);
                data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,gloc);
                data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,gloc);
                data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,gloc);
                data.g(gloc).e(eloc).elabel = cellstr('none');
                data.g(gloc).glabel = gnames(gloc);
                
            case 2 %2 - possible multiple groups but no event types to consider
                
                eloc = 1;
                
                for gloc=1:length(REL.out.labels)
                    
                    data.g(gloc).e(eloc).label = REL.out.labels{gloc};
                    data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,gloc);
                    data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,gloc);
                    data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,gloc);
                    data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,gloc);
                    data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,gloc);
                    data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,gloc);
                    data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,gloc);
                    data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,gloc);
                    data.g(gloc).e(eloc).elabel = cellstr('none');
                    data.g(gloc).glabel = gnames(gloc);
                    
                end
                
            case 3 %3 - possible event types but no groups to consider
                
                gloc = 1;
                
                for eloc=1:length(REL.out.labels)
                    
                    data.g(gloc).e(eloc).label = REL.out.labels{eloc};
                    data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,eloc);
                    data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,eloc);
                    data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,eloc);
                    data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,eloc);
                    data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,eloc);
                    data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,eloc);
                    data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,eloc);
                    data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,eloc);
                    data.g(gloc).e(eloc).elabel = enames(eloc);
                    data.g(gloc).glabel = gnames(gloc);
                    
                end
                
            case 4 %4 - possible groups and event types to consider
                for ii=1:length(REL.out.labels)
                    
                    %use the underscores that were added in psyrat_computevarcomp to
                    %differentiate where the group and event the data are for
                    lblstr = strsplit(REL.out.labels{ii},'_;_');
                    
                    eloc = find(ismember(enames,lblstr(2)));
                    gloc = find(ismember(gnames,lblstr(1)));
                    
                    data.g(gloc).e(eloc).label = REL.out.labels{ii};
                    data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,ii);
                    data.g(gloc).e(eloc).sig_id.raw = REL.out.sig_id(:,ii);
                    data.g(gloc).e(eloc).sig_occ.raw = REL.out.sig_occ(:,ii);
                    data.g(gloc).e(eloc).sig_trl.raw = REL.out.sig_trl(:,ii);
                    data.g(gloc).e(eloc).sig_trlxid.raw = REL.out.sig_trlxid(:,ii);
                    data.g(gloc).e(eloc).sig_occxid.raw = REL.out.sig_occxid(:,ii);
                    data.g(gloc).e(eloc).sig_trlxocc.raw = REL.out.sig_trlxocc(:,ii);
                    data.g(gloc).e(eloc).sig_err.raw = REL.out.sig_err(:,ii);
                    data.g(gloc).e(eloc).elabel = enames(eloc);
                    data.g(gloc).glabel = gnames(gloc);
                    
                end
        end %switch analysis
        
        %store reliability info in relsummary to pass to other functions more easily
        relsummary.relcutoff = relcutoff;
        relsummary.reltype = reltype;
        relsummary.ciperc = ciperc;

        %edge for this arm's variance-component quantiles; see the 'sing' arm
        %above for why the component SDs follow ciperc rather than a fixed 95%
        ciedge = (1-ciperc)./2;
        
        if reltype == 1
            relsummary.reltype_name = 'ic';
        elseif reltype == 2
            relsummary.reltype_name = 'trt';
        elseif reltype == 3
            relsummary.reltype_name = 'ic_trt';
        end
        
        relsummary.gcoeff = gcoeff;

        if gcoeff == 1
            relsummary.gcoeff_name = 'dep';
        elseif gcoeff == 2
            relsummary.gcoeff_name = 'gen';
        end

        %resolve the effective number of occasions (n'_o) for the D-study
        %scaling of the test-retest coefficient and its SEM. The default
        %estimand is a single-occasion score (n'_o = 1), generalizing over the
        %occasion facet. A multi-occasion composite score averages over k
        %occasions; k defaults to the number of observed occasions but may be
        %set explicitly for a D-study projection. The same n'_o is threaded
        %into the coefficient (psyrat_rel_trt) and its SEM (psyrat_variancet),
        %per the author estimand convention (audit sec.14, finding F2).
        [effnocc, nocc_name, nobservedocc] = local_resolve_effnocc(noccmode, noccarg, REL);
        relsummary.noccmode = noccmode;
        relsummary.nocc = effnocc;
        relsummary.nobservedocc = nobservedocc;
        relsummary.nocc_name = nocc_name;


        %store which measure was used to specify cutoff
        switch meascutoff
            case 1
                relsummary.meascutoff = strcat('Lower Limit of ',...
                    sprintf(' %2.0f',ciperc*100),'% Credible Interval');
            case 2
                relsummary.meascutoff = 'Point Estimate';
            case 3
                relsummary.meascutoff = strcat('Upper Limit of ',...
                    sprintf(' %2.0f',ciperc*100),'% Credible Interval');
        end
        
        
        %compute reliability data for each group and event
        switch analysis
            case 1 %no groups or event types to consider
                
                %the same generic structure is used for relsummary, so the event
                %and group locations will both be 1 for the data
                eloc = 1;
                gloc = 1;
                
                %grab the group names (if present)
                relsummary.group(gloc).name = gnames{gloc};
                
                %set the event name as measure
                relsummary.group(gloc).event(eloc).name = 'measure';
                
                try
                    %create empty arrays for storing reliability information
                    trltable = varfun(@length,REL.data{1},'GroupingVariables',{'id','time'});
                catch
                    %create empty arrays for storing reliability information
                    trltable = varfun(@length,REL.data,'GroupingVariables',{'id','time'});
                end
                
                %compute reliabiltiy
                ntrials = max(trltable.GroupCount(:)) + 1000;
                
                
                [llrel,mrel,ulrel] = psyrat_rel_trt(...
                    'gcoeff',gcoeff,...
                    'reltype',reltype,...
                    'bp',data.g(gloc).e(eloc).sig_id.raw,...
                    'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                    'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                    'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                    'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                    'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                    'err',data.g(gloc).e(eloc).sig_err.raw,...
                    'obs',[1 ntrials],'nocc',effnocc,'CI',ciperc);
                
                %find the number of trials to reach cutoff
                switch meascutoff
                    case 1
                        trlcutoff = find(llrel >= relcutoff, 1);
                    case 2
                        trlcutoff = find(mrel >= relcutoff, 1);
                    case 3
                        trlcutoff = find(ulrel >= relcutoff, 1);
                end
                
                
                %see whether the trial cutoff was found. If not store all values as
                %-1. If it is found, store the reliability information about the
                %cutoffs.
                if isempty(trlcutoff)
                    
                    relerr.trlcutoff = 1;
                    trlcutoff = -1;
                    relsummary.group(gloc).event(eloc).trlcutoff = -1;
                    relsummary.group(gloc).event(eloc).relcutoff.m = -1;
                    relsummary.group(gloc).event(eloc).relcutoff.ll = -1;
                    relsummary.group(gloc).event(eloc).relcutoff.ul = -1;
                    
                elseif ~isempty(trlcutoff)
                    
                    %store information about cutoffs
                    
                    relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                    relsummary.group(gloc).event(eloc).relcutoff.m = mrel(trlcutoff);
                    relsummary.group(gloc).event(eloc).relcutoff.ll = llrel(trlcutoff);
                    relsummary.group(gloc).event(eloc).relcutoff.ul = ulrel(trlcutoff);
                    
                end
                
                %find the participants without enough trials based on the cutoffs
                %the guard tests the -1 sentinel VALUE, not isempty: the
                %storage guard above (the if/elseif that records the cutoff,
                %two guard clauses up) reassigned trlcutoff to -1 when the
                %cutoff was unreachable, so an isempty test here can never
                %enter this arm and the fall-through marked every participant
                %good against a GroupCount >= -1 comparison (G54). An
                %unreachable cutoff places every participant in Bad IDs,
                %matching cases 2-4 in both families, the trt_diff helper's
                %sentinel arm, and the manual's ch. 9 promise. (trt_sserr is
                %NOT in that list: it has no trial-cutoff search and no
                %sentinel; its route marks inclusion per subject.)
                if trlcutoff == -1 %if the cutoff was not found in the data

                    %Get information for the participants. REL.data arrives as
                    %a bare table on most estimation cases and as a per-event
                    %cell on psyrat_computevarcomp's event-array cases; a cell
                    %reaching analysis 1 holds exactly one table. The former
                    %try/catch sat around a plain assignment, which cannot
                    %throw, so it could not tell the shapes apart
                    datatrls = REL.data;
                    if iscell(datatrls)
                        datatrls = REL.data{1};
                    end

                    trltable = varfun(@length,datatrls,'GroupingVariables',{'id','time'});
                    
                    %define all of the data as bad, because the cutoff wasn't even
                    %found
                    relsummary.group(gloc).event(eloc).eventgoodids = cell(0,1);
                    relsummary.group(gloc).event(eloc).eventbadids =...
                        unique(trltable.id);
                    
                    relsummary.group(gloc).goodids = cell(0,1);
                    relsummary.group(gloc).badids = ...
                        unique(relsummary.group(gloc).event(eloc).eventbadids);
                    
                    %calculate the reliability estimates for the overall data
                    %over the (universal) bad set; the arm reads datatrls
                    %normalized above, so no separate datatable copy is needed
                    badids = table(unique(relsummary.group(gloc).badids));
                    
                    %check whether there are not enough good data after applying
                    %reliability threshold and extrapolating
                    if isempty(badids)
                        dlg = {'Data do not reach reliability threshold after extrapolation';...
                            'Set a lower reliability threshold'};
                        local_relsummary_notice(interactive, dlg);
                        relerr.nogooddata = 1;
                        return;
                    end
                    
                    trlcdata = varfun(@length,datatrls,...
                        'GroupingVariables',{'id' 'time'});
                    
                    trltable = innerjoin(trlcdata, badids,...
                        'LeftKeys', {'id'}, 'RightKeys', 'Var1',...
                        'LeftVariables', {'id' 'time' 'GroupCount'});
                    
                    trlmean = mean(trltable.GroupCount);
                    trlmed = median(trltable.GroupCount);
                    
                    %calculate reliability using either the mean or median (based
                    %on user input)
                    switch relcentmeas
                        case 1
                            relcent = trlmean;
                        case 2
                            relcent = trlmed;
                    end
                    
                    
                    [relsummary.group(gloc).event(eloc).rel.ll,...
                        relsummary.group(gloc).event(eloc).rel.m,...
                        relsummary.group(gloc).event(eloc).rel.ul] = ...
                        psyrat_rel_trt(...
                        'gcoeff',gcoeff,...
                        'reltype',reltype,...
                        'bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                        'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                        'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                        'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                        'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                        'err',data.g(gloc).e(eloc).sig_err.raw,...
                        'obs',relcent,'nocc',effnocc,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).rel.meas = relcent;
                    
                    relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                    relsummary.group(gloc).event(eloc).goodn = 0;
                    
                else

                    %Get trial information for those that meet cutoff
                    datatrls = REL.data;

                    try
                        trltable = varfun(@length,datatrls,'GroupingVariables',{'id','time'});
                    catch
                        trltable = varfun(@length,datatrls{:},'GroupingVariables',{'id','time'});
                    end
                    
                    ind2include = trltable.GroupCount >= trlcutoff;
                    ind2exclude = trltable.GroupCount < trlcutoff;
                    
                    prelim_bad = unique(trltable.id(ind2exclude));
                    
                    for jpbad = 1:length(prelim_bad)
                        ind_id = find(strcmp(trltable.id,prelim_bad(jpbad)));
                        ind2include(ind_id) = 0;
                        ind2exclude(ind_id) = 1;
                    end
                    
                    %store the good ids and the bad ids (ids that don't meet the
                    %cutoff)
                    relsummary.group(gloc).event(eloc).eventgoodids =...
                        unique(trltable.id(ind2include));
                    relsummary.group(gloc).event(eloc).eventbadids =...
                        unique(trltable.id(ind2exclude));
                    
                    relsummary.group(gloc).goodids = ...
                        unique(relsummary.group(gloc).event(eloc).eventgoodids);
                    relsummary.group(gloc).badids = ...
                        unique(relsummary.group(gloc).event(eloc).eventbadids);
                    
                    datatable = REL.data;
                    
                    if iscell(datatable)
                        datatable = REL.data{1};
                    end
                    
                    goodids = table(unique(relsummary.group(gloc).goodids));
                    
                    trlcdata = varfun(@length,datatable,...
                        'GroupingVariables',{'id' 'time'});
                    
                    trltable = innerjoin(trlcdata, goodids,...
                        'LeftKeys', {'id'}, 'RightKeys', 'Var1',...
                        'LeftVariables', {'id' 'time' 'GroupCount'});
                    
                    trlmean = mean(trltable.GroupCount);
                    trlmed = median(trltable.GroupCount);
                    
                    %calculate reliability using either the mean or median (based
                    %on user input)
                    switch relcentmeas
                        case 1
                            relcent = trlmean;
                        case 2
                            relcent = trlmed;
                    end
                    
                    [relsummary.group(gloc).event(eloc).rel.ll,...
                        relsummary.group(gloc).event(eloc).rel.m,...
                        relsummary.group(gloc).event(eloc).rel.ul] = ...
                        psyrat_rel_trt(...
                        'gcoeff',gcoeff,...
                        'reltype',reltype,...
                        'bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                        'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                        'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                        'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                        'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                        'err',data.g(gloc).e(eloc).sig_err.raw,...
                        'obs',relcent,'nocc',effnocc,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).rel.meas = relcent;
                    
                    relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                    relsummary.group(gloc).event(eloc).goodn = height(goodids);
                    
                    %just in case a cutoff was extrapolated, the reliability
                    %information will be overwritten
                    if  relsummary.group(gloc).event(eloc).trlcutoff >...
                            relsummary.group(gloc).event(eloc).trlinfo.max
                        
                        relerr.trlmax = 1;
                        relsummary.group(gloc).event(eloc).relcutoff.m = -1;
                        relsummary.group(gloc).event(eloc).relcutoff.ll = -1;
                        relsummary.group(gloc).event(eloc).relcutoff.ul = -1;
                        
                    end
                    
                end

                %calculate iccs, between-subject standard deviations, and
                %within-subject standard deviations. This block sits OUTSIDE the
                %cutoff if/else on purpose (G54): the posterior summaries it
                %reads do not depend on the retained set, and the sentinel arm
                %must also produce these fields, matching case 2 and the sing
                %case-1 sibling, both of which compute them on a common path.
                [relsummary.group(gloc).event(eloc).icc.ll,...
                    relsummary.group(gloc).event(eloc).icc.m,...
                    relsummary.group(gloc).event(eloc).icc.ul] = ...
                    psyrat_rel_trt(...
                    'gcoeff',gcoeff,...
                    'reltype',reltype,...
                    'bp',data.g(gloc).e(eloc).sig_id.raw,...
                    'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                    'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                    'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                    'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                    'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                    'err',data.g(gloc).e(eloc).sig_err.raw,...
                    'obs',1,'nocc',1,'CI',ciperc);

                relsummary.group(gloc).event(eloc).betsd.m = ...
                    mean(data.g(gloc).e(eloc).sig_id.raw);
                relsummary.group(gloc).event(eloc).betsd.ll = ...
                    quantile(data.g(gloc).e(eloc).sig_id.raw,ciedge);
                relsummary.group(gloc).event(eloc).betsd.ul = ...
                    quantile(data.g(gloc).e(eloc).sig_id.raw,1-ciedge);

                relsummary.group(gloc).event(eloc).witsd.m = ...
                    mean(data.g(gloc).e(eloc).sig_err.raw);
                relsummary.group(gloc).event(eloc).witsd.ll = ...
                    quantile(data.g(gloc).e(eloc).sig_err.raw,ciedge);
                relsummary.group(gloc).event(eloc).witsd.ul = ...
                    quantile(data.g(gloc).e(eloc).sig_err.raw,1-ciedge);

            case 2 %possible multiple groups but no event types to consider
                
                
                %since the same generic structure is used for relsummary, the event
                %location will be defined as 1.
                eloc = 1;
                
                for gloc=1:ngroups %loop through each group
                    
                    %store the name of the group
                    if eloc == 1
                        relsummary.group(gloc).name = gnames{gloc};
                    end
                    
                    %store the name of the event
                    relsummary.group(gloc).event(eloc).name = enames{eloc};
                    
                    try
                        %create empty arrays for storing reliability information
                        trltable = varfun(@length,REL.data{1},...
                            'GroupingVariables',{'id','time'});
                    catch
                        %create empty arrays for storing reliability information
                        trltable = varfun(@length,REL.data,...
                            'GroupingVariables',{'id','time'});
                    end
                    
                    
                    %compute reliabiltiy
                    ntrials = max(trltable.GroupCount(:)) + 1000;
                    
                    
                    
                    [llrel,mrel,ulrel] = psyrat_rel_trt(...
                        'gcoeff',gcoeff,...
                        'reltype',reltype,...
                        'bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                        'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                        'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                        'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                        'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                        'err',data.g(gloc).e(eloc).sig_err.raw,...
                        'obs',[1 ntrials],'nocc',effnocc,'CI',ciperc);
                    
                    %find the number of trials to reach cutoff
                    switch meascutoff
                        case 1
                            trlcutoff = find(llrel >= relcutoff, 1);
                        case 2
                            trlcutoff = find(mrel >= relcutoff, 1);
                        case 3
                            trlcutoff = find(ulrel >= relcutoff, 1);
                    end
                    
                    
                    
                    %see whether the trial cutoff was found. If not store all
                    %values as -1. If it is found, store the reliability
                    %information about the cutoffs.
                    
                    if isempty(trlcutoff)
                        
                        relerr.trlcutoff = 1;
                        trlcutoff = -1;
                        relsummary.group(gloc).event(eloc).trlcutoff = -1;
                        relsummary.group(gloc).event(eloc).relcutoff.m = -1;
                        relsummary.group(gloc).event(eloc).relcutoff.ll = -1;
                        relsummary.group(gloc).event(eloc).relcutoff.ul = -1;
                        
                    elseif ~isempty(trlcutoff)
                        
                        %store information about cutoffs
                        
                        relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                        relsummary.group(gloc).event(eloc).relcutoff.m = mrel(trlcutoff);
                        relsummary.group(gloc).event(eloc).relcutoff.ll = llrel(trlcutoff);
                        relsummary.group(gloc).event(eloc).relcutoff.ul = ulrel(trlcutoff);
                        
                    end
                    
                    %if the trial cutoff was not found in the specified data
                    if trlcutoff == -1
                        
                        %store all the ids as bad
                        try
                            datatrls = REL.data{1};
                        catch
                            datatrls = REL.data;
                        end
                        
                        ind = strcmp(datatrls.group,gnames{gloc});
                        datatrls = datatrls(ind,:);
                        
                        trltable = varfun(@length,datatrls,...
                            'GroupingVariables',{'id' 'time'});
                        
                        relsummary.group(gloc).event(eloc).eventgoodids = cell(0,1);
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            unique(trltable.id);
                        
                    else
                        
                        %store the ids for participants with enough trials as good,
                        %and the ids for participants with too few trials as bad
                        try
                            datatrls = REL.data{1};
                        catch
                            datatrls = REL.data;
                        end
                        
                        ind = strcmp(datatrls.group,gnames{gloc});
                        datatrls = datatrls(ind,:);
                        
                        trltable = varfun(@length,datatrls,...
                            'GroupingVariables',{'id' 'time'});
                        
                        ind2include = trltable.GroupCount >= trlcutoff;
                        ind2exclude = trltable.GroupCount < trlcutoff;
                        
                        prelim_bad = unique(trltable.id(ind2exclude));
                        
                        for jpbad = 1:length(prelim_bad)
                            ind_id = find(strcmp(trltable.id,prelim_bad(jpbad)));
                            ind2include(ind_id) = 0;
                            ind2exclude(ind_id) = 1;
                        end
                        
                        relsummary.group(gloc).event(eloc).eventgoodids =...
                            unique(trltable.id(ind2include));
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            unique(trltable.id(ind2exclude));
                        
                    end
                    
                    
                end
                
                %since there is only 1 event, those participants with good data for
                %the event will be store as participants having good data for the
                %whole group (this step makes more sense when there is more than
                %one event per group)
                
                for gloc=1:ngroups
                    
                    relsummary.group(gloc).goodids = ...
                        relsummary.group(gloc).event(eloc).eventgoodids;
                    relsummary.group(gloc).badids = ...
                        relsummary.group(gloc).event(eloc).eventbadids;
                    
                end
                
                %calculate the reliability for the overall data for each group
                for gloc=1:ngroups
                    
                    try
                        datatable = REL.data{1};
                    catch
                        datatable = REL.data;
                    end
                    
                    ind = strcmp(datatable.group,gnames{gloc});
                    datasubset = datatable(ind,:);
                    
                    %only factor in the trial counts from those with good data
                    
                    %goodids is 'none' (a char) when no trial cutoff was found,
                    %but it is an EMPTY CELL when a cutoff was found and then
                    %excluded everyone (an extrapolated cutoff above every
                    %participant's trial count). BOTH must reach arm 2, and an
                    %"elseif" second arm risked leaving goodids unassigned and
                    %crashing the isempty() check below; a plain else is
                    %exhaustive and matches every sibling guard in this file
                    %(e.g. :4050). The arm test keys on TYPE and EMPTINESS
                    %(B29), never on the sentinel's VALUE, so arm 1 is taken
                    %only for a genuine non-empty id cell -- including one that
                    %holds a participant literally named 'none'.
                    if iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids)
                        goodids = table(unique(relsummary.group(gloc).goodids));
                    else
                        goodids = table(unique(relsummary.group(gloc).badids));
                    end

                    %check whether there are not enough good data after applying
                    %reliability threshold and extrapolating
                    %
                    %This guard cannot fire AT THIS SITE, and that is deliberate
                    %-- leave it alone (RC-40). Both arms above yield a non-empty
                    %table here: arm 2 is reached exactly when every participant
                    %was excluded, which is exactly when badids holds all of
                    %them. So isempty() is false whenever the group has anyone in
                    %it at all.
                    %
                    %The deadness is load-bearing, not an oversight waiting to be
                    %corrected. Firing would return from the whole function --
                    %every group, every event -- because one group was bad, and
                    %would replace a run that completes and reports the problem
                    %through relerr.trlmax ("Cutoffs represent an extrapolation
                    %beyond the data") with a silent abort that draws nothing.
                    %ONE committed test depends on it staying dead, not two:
                    %testTrtExtrapolatedCutoffWithNoRetainedIdsStillWritesGoodN
                    %reaches this exact state -- profiler hits = 2 on the guard
                    %test below, once per group, via arm 2 above -- and asserts
                    %goodn and trlmax are written afterward. The other test this
                    %comment used to name, testCriterionScreensReportTrialCutoff
                    %Failure, executes NONE of the ten guard sites (measured
                    %2026-08-18: hits = 0 on all ten across its three calls), so
                    %it constrains nothing here. Its fixture comment still says
                    %otherwise and is wrong about the mechanism.
                    %
                    %DO NOT generalize this to the other guard sites. The
                    %discriminator is whether badids is guaranteed non-empty when
                    %arm 2 is taken, and it is NOT everywhere. The sites fed by
                    %the cross-event aggregation loops start badids = [] (e.g.
                    %:1204) and populate it only inside the eventgoodids arm
                    %test, so when NO event finds a cutoff
                    %badids stays empty, table([]) is empty, and the guard FIRES
                    %normally. Proven by profiler coverage at :1483 --
                    %makeICRelDataMultiGroupEvent at depcutoff 0.9999 executes it
                    %and raises this dialog. Those sites are working as intended
                    %and need nothing.
                    %
                    %What RC-40 actually fixed is a third case: the two
                    %difference-score helpers at the end of this file set
                    %nogooddata with no dialog and no return, so the run carried
                    %on to psyrat_relfigures, aborted there, and told the user
                    %nothing. RC-40 made them raise this notice for themselves;
                    %G45 and then G48 (owner rulings 2026-08-29) replaced both
                    %helpers' flag-and-notice with the strict-fold per-group
                    %console line, so NEITHER helper sets the flag any more and
                    %the ten guarded sites like this one are the only setters.
                    if isempty(goodids)
                        dlg = {'Data do not reach reliability threshold after extrapolation';...
                            'Set a lower reliability threshold'};
                        local_relsummary_notice(interactive, dlg);
                        relerr.nogooddata = 1;
                        return;
                    end
                    
                    trlcdata = varfun(@length,datasubset,...
                        'GroupingVariables',{'id' 'time'});
                    
                    trltable = innerjoin(trlcdata, goodids,...
                        'LeftKeys', {'id'}, 'RightKeys', 'Var1',...
                        'LeftVariables', {'id' 'time' 'GroupCount'});
                    
                    trlmean = mean(trltable.GroupCount);
                    trlmed = median(trltable.GroupCount);
                    
                    %calculate reliability using either the mean or median (based
                    %on user input)
                    switch relcentmeas
                        case 1
                            relcent = trlmean;
                        case 2
                            relcent = trlmed;
                    end
                    
                    [relsummary.group(gloc).event(eloc).rel.ll,...
                        relsummary.group(gloc).event(eloc).rel.m,...
                        relsummary.group(gloc).event(eloc).rel.ul] = ...
                        psyrat_rel_trt(...
                        'gcoeff',gcoeff,...
                        'reltype',reltype,...
                        'bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                        'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                        'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                        'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                        'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                        'err',data.g(gloc).e(eloc).sig_err.raw,...
                        'obs',relcent,'nocc',effnocc,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).rel.meas = relcent;
                    relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                    
                    %just in case a cutoff was extrapolated, the reliability
                    %information will be overwritten
                    if  relsummary.group(gloc).event(eloc).trlcutoff >...
                            relsummary.group(gloc).event(eloc).trlinfo.max
                        
                        relerr.trlmax = 1;
                        relsummary.group(gloc).event(eloc).relcutoff.m = -1;
                        relsummary.group(gloc).event(eloc).relcutoff.ll = -1;
                        relsummary.group(gloc).event(eloc).relcutoff.ul = -1;
                        
                    end
                    
                    %same non-exhaustive guard as the goodids assignment above.
                    %Note this site was MASKED until that one was fixed: the
                    %isempty(goodids) check above threw first, so an unwritten
                    %goodn only becomes reachable once the first guard is
                    %corrected. Fixing one without the other moves the error
                    %rather than removing it. Zero retained participants is the
                    %correct value on the empty-cell path.
                    if iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids)
                        relsummary.group(gloc).event(eloc).goodn = height(goodids);
                    else
                        relsummary.group(gloc).event(eloc).goodn = 0;
                    end
                    
                    %calculate iccs, between-subject standard deviations, and
                    %within-subject standard deviations
                    
                    [relsummary.group(gloc).event(eloc).icc.ll,...
                        relsummary.group(gloc).event(eloc).icc.m,...
                        relsummary.group(gloc).event(eloc).icc.ul] = ...
                        psyrat_rel_trt(...
                        'gcoeff',gcoeff,...
                        'reltype',reltype,...
                        'bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                        'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                        'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                        'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                        'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                        'err',data.g(gloc).e(eloc).sig_err.raw,...
                        'obs',1,'nocc',1,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).betsd.m = ...
                        mean(data.g(gloc).e(eloc).sig_id.raw);
                    relsummary.group(gloc).event(eloc).betsd.ll = ...
                        quantile(data.g(gloc).e(eloc).sig_id.raw,ciedge);
                    relsummary.group(gloc).event(eloc).betsd.ul = ...
                        quantile(data.g(gloc).e(eloc).sig_id.raw,1-ciedge);
                    
                    relsummary.group(gloc).event(eloc).witsd.m = ...
                        mean(data.g(gloc).e(eloc).sig_err.raw);
                    relsummary.group(gloc).event(eloc).witsd.ll = ...
                        quantile(data.g(gloc).e(eloc).sig_err.raw,ciedge);
                    relsummary.group(gloc).event(eloc).witsd.ul = ...
                        quantile(data.g(gloc).e(eloc).sig_err.raw,1-ciedge);
                    
                    
                end
                
            case 3 %possible event types but no groups to consider
                
                %since the relsummary structure is generic for any number of groups
                %or events, the data will count as being from 1 group.
                gloc = 1;
                
                %cylce through each event
                for eloc=1:nevents
                    
                    %get the group name
                    if eloc == 1
                        relsummary.group(gloc).name = gnames{gloc};
                    end
                    
                    %get the event name
                    relsummary.group(gloc).event(eloc).name = enames{eloc};
                    
                    
                    %create empty arrays for storing reliability information
                    trltable = varfun(@length,...
                        REL.data(strcmp(REL.data.event,enames{eloc}),:),...
                        'GroupingVariables',{'id','time'});
                    
                    
                    %compute reliabiltiy
                    ntrials = max(trltable.GroupCount(:)) + 1000;
                    
                    
                    [llrel,mrel,ulrel] = psyrat_rel_trt(...
                        'gcoeff',gcoeff,...
                        'reltype',reltype,...
                        'bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                        'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                        'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                        'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                        'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                        'err',data.g(gloc).e(eloc).sig_err.raw,...
                        'obs',[1 ntrials],'nocc',effnocc,'CI',ciperc);
                    
                    %find the number of trials to reach cutoff
                    switch meascutoff
                        case 1
                            trlcutoff = find(llrel >= relcutoff, 1);
                        case 2
                            trlcutoff = find(mrel >= relcutoff, 1);
                        case 3
                            trlcutoff = find(ulrel >= relcutoff, 1);
                    end
                    
                    
                    if isempty(trlcutoff) %if a cutoff wasn't found
                        
                        relerr.trlcutoff = 1;
                        trlcutoff = -1;
                        relsummary.group(gloc).event(eloc).trlcutoff = -1;
                        relsummary.group(gloc).event(eloc).relcutoff.m = -1;
                        relsummary.group(gloc).event(eloc).relcutoff.ll = -1;
                        relsummary.group(gloc).event(eloc).relcutoff.ul = -1;
                        
                        %store all of the participant ids as bad, because none
                        %reached the cutoff
                        relsummary.group(gloc).event(eloc).eventgoodids =...
                            cell(0,1);
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            unique(trltable.id);
                        
                    elseif ~isempty(trlcutoff) %if a cutoff was found
                        
                        %store information about cutoffs
                        
                        relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                        relsummary.group(gloc).event(eloc).relcutoff.m = mrel(trlcutoff);
                        relsummary.group(gloc).event(eloc).relcutoff.ll = llrel(trlcutoff);
                        relsummary.group(gloc).event(eloc).relcutoff.ul = ulrel(trlcutoff);
                        
                        datatrls = REL.data(strcmp(REL.data.event,enames{eloc}),:);
                        
                        trltable = varfun(@length,datatrls,'GroupingVariables',{'id','time'});
                        
                        ind2include = trltable.GroupCount >= trlcutoff;
                        ind2exclude = trltable.GroupCount < trlcutoff;
                        
                        prelim_bad = unique(trltable.id(ind2exclude));
                        
                        for jpbad = 1:length(prelim_bad)
                            ind_id = find(strcmp(trltable.id,prelim_bad(jpbad)));
                            ind2include(ind_id) = 0;
                            ind2exclude(ind_id) = 1;
                        end
                        
                        %store which participants had enough trials to meet cutoff
                        relsummary.group(gloc).event(eloc).eventgoodids =...
                            unique(trltable.id(ind2include));
                        relsummary.group(gloc).event(eloc).eventbadids =...
                            unique(trltable.id(ind2exclude));
                        
                    end
                    
                end
                
                %cycle through the events and store information about which
                %participants have good data across all event types.
                tempids = {};
                badids = [];
                for eloc=1:nevents
                    %STRICT all-events rule (G39 ruling 2026-08-29): an event that
                    %retained nobody folds in as an empty set and empties the running
                    %intersection instead of being skipped; its eventbadids (the full
                    %id list for a dead event, by producer convention) fold in too, so
                    %the group good/bad ids complement.
                    tempids{end+1} = relsummary.group(gloc).event(eloc).eventgoodids;
                    if eloc == 1
                        badids = relsummary.group(gloc).event(eloc).eventbadids;
                    elseif eloc > 1
                        [~,ind]=setdiff(tempids{1},tempids{end});
                        new = tempids{1};
                        badids = vertcat(badids,...
                            relsummary.group(gloc).event(eloc).eventbadids,...
                            new(ind));
                        new(ind) = [];
                        tempids{1} = new;
                    end
                end
                
                %defensive: with the unconditional fold above, tempids is empty only when nevents == 0; keep the empty-cell carrier convention
                if isempty(tempids)
                    tempids{1} = cell(0,1);
                end
                
                %store information about participants with good/bad ids
                relsummary.group(gloc).goodids = tempids{1};
                relsummary.group(gloc).badids = unique(badids);
                %STRICT-fold consequence (G39): a group can retain nobody while
                %its per-event rows stay informative, so warn on the console and
                %continue. Deliberately NOT local_relsummary_notice: its modal
                %contract requires every caller to set relerr.nogooddata and
                %return, and this is not an abort.
                if ~(iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids))
                    fprintf(['\nNo participants met the reliability threshold ',...
                        'for every event type in a group; that group''s ',...
                        'n Included is 0.\n']);
                end
                
                %cycle through each event
                for eloc=1:nevents
                    
                    %get trial information to use for reliability estimates. only
                    %participants with enough data will contribute toward trial
                    %counts
                    datatable = REL.data(strcmp(REL.data.event,enames{eloc}),:);
                    
                    if iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids)
                        goodids = table(unique(relsummary.group(gloc).goodids));
                    else
                        goodids = table(unique(relsummary.group(gloc).badids));
                    end
                    
                    %check whether there are not enough good data after applying
                    %reliability threshold and extrapolating
                    if isempty(goodids)
                        dlg = {'Data do not reach reliability threshold after extrapolation';...
                            'Set a lower reliability threshold'};
                        local_relsummary_notice(interactive, dlg);
                        relerr.nogooddata = 1;
                        return;
                    end
                    
                    trlcdata = varfun(@length,datatable,...
                        'GroupingVariables',{'id' 'time'});
                    
                    trltable = innerjoin(trlcdata, goodids,...
                        'LeftKeys', {'id'}, 'RightKeys', 'Var1',...
                        'LeftVariables', {'id' 'time' 'GroupCount'});
                    
                    trlmean = mean(trltable.GroupCount);
                    trlmed = median(trltable.GroupCount);
                    
                    %calculate reliability using either the mean or median (based
                    %on user input)
                    switch relcentmeas
                        case 1
                            relcent = trlmean;
                        case 2
                            relcent = trlmed;
                    end
                    
                    [relsummary.group(gloc).event(eloc).rel.ll,...
                        relsummary.group(gloc).event(eloc).rel.m,...
                        relsummary.group(gloc).event(eloc).rel.ul] = ...
                        psyrat_rel_trt(...
                        'gcoeff',gcoeff,...
                        'reltype',reltype,...
                        'bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                        'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                        'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                        'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                        'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                        'err',data.g(gloc).e(eloc).sig_err.raw,...
                        'obs',relcent,'nocc',effnocc,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).rel.meas = relcent;
                    
                    relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                    relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                    relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                    
                    %just in case a cutoff was extrapolated, the reliability
                    %information will be overwritten
                    if  relsummary.group(gloc).event(eloc).trlcutoff >...
                            relsummary.group(gloc).event(eloc).trlinfo.max
                        
                        relerr.trlmax = 1;
                        relsummary.group(gloc).event(eloc).relcutoff.m = -1;
                        relsummary.group(gloc).event(eloc).relcutoff.ll = -1;
                        relsummary.group(gloc).event(eloc).relcutoff.ul = -1;
                        
                    end
                    
                    if iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids)
                        relsummary.group(gloc).event(eloc).goodn = height(goodids);
                    else
                        relsummary.group(gloc).event(eloc).goodn = 0;
                    end
                    
                    %calculate iccs, between-subject standard deviations, and
                    %within-subject standard deviations
                    [relsummary.group(gloc).event(eloc).icc.ll,...
                        relsummary.group(gloc).event(eloc).icc.m,...
                        relsummary.group(gloc).event(eloc).icc.ul] = ...
                        psyrat_rel_trt(...
                        'gcoeff',gcoeff,...
                        'reltype',reltype,...
                        'bp',data.g(gloc).e(eloc).sig_id.raw,...
                        'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                        'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                        'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                        'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                        'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                        'err',data.g(gloc).e(eloc).sig_err.raw,...
                        'obs',1,'nocc',1,'CI',ciperc);
                    
                    relsummary.group(gloc).event(eloc).betsd.m = ...
                        mean(data.g(gloc).e(eloc).sig_id.raw);
                    relsummary.group(gloc).event(eloc).betsd.ll = ...
                        quantile(data.g(gloc).e(eloc).sig_id.raw,ciedge);
                    relsummary.group(gloc).event(eloc).betsd.ul = ...
                        quantile(data.g(gloc).e(eloc).sig_id.raw,1-ciedge);
                    
                    relsummary.group(gloc).event(eloc).witsd.m = ...
                        mean(data.g(gloc).e(eloc).sig_err.raw);
                    relsummary.group(gloc).event(eloc).witsd.ll = ...
                        quantile(data.g(gloc).e(eloc).sig_err.raw,ciedge);
                    relsummary.group(gloc).event(eloc).witsd.ul = ...
                        quantile(data.g(gloc).e(eloc).sig_err.raw,1-ciedge);
                    
                end
                
            case 4 %groups and event types to consider
                
                %cycle through each event
                for eloc=1:nevents
                    
                    %cycle through each group
                    for gloc=1:ngroups
                        
                        %store the group name
                        if eloc == 1
                            relsummary.group(gloc).name = gnames{gloc};
                        end
                        
                        %store the event name
                        relsummary.group(gloc).event(eloc).name = enames{eloc};
                        
                        %create empty arrays for storing reliability information
                        trltable = varfun(@length,REL.data(all(...
                            strcmp(REL.data.group,gnames{gloc}) &...
                            strcmp(REL.data.event,enames{eloc})...
                            ,2),:),...
                            'GroupingVariables',{'id','time'});
                        
                        %compute reliability
                        ntrials = max(trltable.GroupCount(:)) + 1000;
                        
                        
                        
                        
                        [llrel,mrel,ulrel] = psyrat_rel_trt(...
                            'gcoeff',gcoeff,...
                            'reltype',reltype,...
                            'bp',data.g(gloc).e(eloc).sig_id.raw,...
                            'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                            'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                            'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                            'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                            'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                            'err',data.g(gloc).e(eloc).sig_err.raw,...
                            'obs',[1 ntrials],'nocc',effnocc,'CI',ciperc);
                        
                        %find the number of trials to reach cutoff
                        switch meascutoff
                            case 1
                                trlcutoff = find(llrel >= relcutoff, 1);
                            case 2
                                trlcutoff = find(mrel >= relcutoff, 1);
                            case 3
                                trlcutoff = find(ulrel >= relcutoff, 1);
                        end
                        
                        
                        
                        %if a cutoff was not found
                        if isempty(trlcutoff)
                            
                            %store all the participant ids as having bad data
                            relerr.trlcutoff = 1;
                            trlcutoff = -1;
                            relsummary.group(gloc).event(eloc).trlcutoff = -1;
                            relsummary.group(gloc).event(eloc).relcutoff.m = -1;
                            relsummary.group(gloc).event(eloc).relcutoff.ll = -1;
                            relsummary.group(gloc).event(eloc).relcutoff.ul = -1;
                            
                            datatrls = REL.data(all(...
                                strcmp(REL.data.group,gnames{gloc}) &...
                                strcmp(REL.data.event,enames{eloc})...
                                ,2),:);
                            
                            ind = strcmp(datatrls.group,gnames{gloc});
                            datatrls = datatrls(ind,:);
                            
                            trltable = varfun(@length,datatrls,'GroupingVariables',{'id','time'});
                            
                            relsummary.group(gloc).event(eloc).eventgoodids =...
                                cell(0,1);
                            relsummary.group(gloc).event(eloc).eventbadids =...
                                unique(trltable.id);
                            
                        elseif ~isempty(trlcutoff) %if a cutoff was found
                            
                            %store information about cutoffs
                            
                            relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                            relsummary.group(gloc).event(eloc).relcutoff.m = mrel(trlcutoff);
                            relsummary.group(gloc).event(eloc).relcutoff.ll = llrel(trlcutoff);
                            relsummary.group(gloc).event(eloc).relcutoff.ul = ulrel(trlcutoff);
                            
                            datatrls = REL.data(all(...
                                strcmp(REL.data.group,gnames{gloc}) &...
                                strcmp(REL.data.event,enames{eloc})...
                                ,2),:);
                            
                            ind = strcmp(datatrls.group,gnames{gloc});
                            datatrls = datatrls(ind,:);
                            
                            %find ids with enough trials based on cutoff
                            trltable = varfun(@length,datatrls,'GroupingVariables',{'id','time'});
                            
                            ind2include = trltable.GroupCount >= trlcutoff;
                            ind2exclude = trltable.GroupCount < trlcutoff;
                            
                            prelim_bad = unique(trltable.id(ind2exclude));
                            
                            for jpbad = 1:length(prelim_bad)
                                ind_id = find(strcmp(trltable.id,prelim_bad(jpbad)));
                                ind2include(ind_id) = 0;
                                ind2exclude(ind_id) = 1;
                            end
                            
                            relsummary.group(gloc).event(eloc).eventgoodids =...
                                unique(trltable.id(ind2include));
                            relsummary.group(gloc).event(eloc).eventbadids =...
                                unique(trltable.id(ind2exclude));
                            
                        end
                        
                    end
                end
                
                %find the ids that have enough trials for each event type
                
                for gloc=1:ngroups
                    tempids = {};
                    badids = [];
                    
                    for eloc=1:nevents
                        %STRICT all-events rule (G39 ruling 2026-08-29): an event that
                        %retained nobody folds in as an empty set and empties the running
                        %intersection instead of being skipped; its eventbadids (the full
                        %id list for a dead event, by producer convention) fold in too, so
                        %the group good/bad ids complement.
                        tempids{end+1} = relsummary.group(gloc).event(eloc).eventgoodids;
                        if eloc == 1
                            badids = relsummary.group(gloc).event(eloc).eventbadids;
                        elseif eloc > 1
                            [~,ind]=setdiff(tempids{1},tempids{end});
                            new = tempids{1};
                            badids = vertcat(badids,...
                                relsummary.group(gloc).event(eloc).eventbadids,...
                                new(ind));
                            new(ind) = [];
                            tempids{1} = new;
                        end
                    end
                    
                    %defensive: with the unconditional fold above, tempids is empty only when nevents == 0; keep the empty-cell carrier convention
                    if isempty(tempids)
                        tempids{1} = cell(0,1);
                    end
                    
                    relsummary.group(gloc).goodids = tempids{1};
                    relsummary.group(gloc).badids = unique(badids);
                    %STRICT-fold consequence (G39): a group can retain nobody while
                    %its per-event rows stay informative, so warn on the console and
                    %continue. Deliberately NOT local_relsummary_notice: its modal
                    %contract requires every caller to set relerr.nogooddata and
                    %return, and this is not an abort.
                    if ~(iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids))
                        fprintf(['\nNo participants met the reliability threshold ',...
                            'for every event type in a group; that group''s ',...
                            'n Included is 0.\n']);
                    end
                    
                end
                
                %calculate reliability estimates for the overall data for each
                %group and event type
                for eloc=1:nevents
                    
                    for gloc=1:ngroups
                        
                        %pull the data to calculate the trial information from
                        %participants with good data
                        datatable = REL.data(all(...
                            strcmp(REL.data.group,gnames{gloc}) &...
                            strcmp(REL.data.event,enames{eloc})...
                            ,2),:);
                        
                        ind = strcmp(datatable.group,gnames{gloc});
                        datasubset = datatable(ind,:);
                        
                        if iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids)
                            goodids = table(unique(relsummary.group(gloc).goodids));
                        else
                            goodids = table(unique(relsummary.group(gloc).badids));
                        end
                        
                        %check whether there are not enough good data after applying
                        %reliability threshold and extrapolating
                        if isempty(goodids)
                            dlg = {'Data do not reach reliability threshold after extrapolation';...
                                'Set a lower reliability threshold'};
                            local_relsummary_notice(interactive, dlg);
                            relerr.nogooddata = 1;
                            return;
                        end
                        
                        trlcdata = varfun(@length,datasubset,...
                            'GroupingVariables',{'id' 'time'});
                        
                        trltable = innerjoin(trlcdata, goodids,...
                            'LeftKeys', {'id'}, 'RightKeys', 'Var1',...
                            'LeftVariables', {'id' 'time' 'GroupCount'});
                        
                        trlmean = mean(trltable.GroupCount);
                        trlmed = median(trltable.GroupCount);
                        
                        %calculate reliability using either the mean or median
                        %(based on user input)
                        switch relcentmeas
                            case 1
                                relcent = trlmean;
                            case 2
                                relcent = trlmed;
                        end
                        
                        [relsummary.group(gloc).event(eloc).rel.ll,...
                            relsummary.group(gloc).event(eloc).rel.m,...
                            relsummary.group(gloc).event(eloc).rel.ul] = ...
                            psyrat_rel_trt(...
                            'gcoeff',gcoeff,...
                            'reltype',reltype,...
                            'bp',data.g(gloc).e(eloc).sig_id.raw,...
                            'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                            'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                            'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                            'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                            'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                            'err',data.g(gloc).e(eloc).sig_err.raw,...
                            'obs',relcent,'nocc',effnocc,'CI',ciperc);
                        
                        relsummary.group(gloc).event(eloc).rel.meas = relcent;
                        relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                        relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                        relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                        relsummary.group(gloc).event(eloc).trlinfo.std = std(trltable.GroupCount);
                        relsummary.group(gloc).event(eloc).trlinfo.med = trlmed;
                        
                        %check if trial cutoff for an event exceeded the total
                        %number of trials for any subjects
                        if  relsummary.group(gloc).event(eloc).trlcutoff >...
                                relsummary.group(gloc).event(eloc).trlinfo.max
                            
                            relerr.trlmax = 1;
                            relsummary.group(gloc).event(eloc).relcutoff.m = -1;
                            relsummary.group(gloc).event(eloc).relcutoff.ll = -1;
                            relsummary.group(gloc).event(eloc).relcutoff.ul = -1;
                            
                        end
                        
                        %goodn is the GROUP-level retained count (G39
                        %ruling 2026-08-29): keyed on the group carrier,
                        %exactly as the one-group branches key it, so a
                        %dead event's row and its siblings print the same
                        %complementing count. The event-carrier keying
                        %this replaces printed 0 on a dead event's row
                        %while its siblings printed the group figure, and
                        %under the strict fold it would have counted the
                        %badids fallback table instead.
                        if iscell(relsummary.group(gloc).goodids) && ~isempty(relsummary.group(gloc).goodids)
                            relsummary.group(gloc).event(eloc).goodn = height(goodids);
                        else
                            relsummary.group(gloc).event(eloc).goodn = 0;
                        end
                        
                        %calculate iccs, between-subject standard deviations, and
                        %within-subject standard deviations
                        [relsummary.group(gloc).event(eloc).icc.ll,...
                            relsummary.group(gloc).event(eloc).icc.m,...
                            relsummary.group(gloc).event(eloc).icc.ul] = ...
                            psyrat_rel_trt(...
                            'gcoeff',gcoeff,...
                            'reltype',reltype,...
                            'bp',data.g(gloc).e(eloc).sig_id.raw,...
                            'bo',data.g(gloc).e(eloc).sig_occ.raw,...
                            'bt',data.g(gloc).e(eloc).sig_trl.raw,...
                            'txp',data.g(gloc).e(eloc).sig_trlxid.raw,...
                            'oxp',data.g(gloc).e(eloc).sig_occxid.raw,...
                            'txo',data.g(gloc).e(eloc).sig_trlxocc.raw,...
                            'err',data.g(gloc).e(eloc).sig_err.raw,...
                            'obs',1,'nocc',1,'CI',ciperc);
                        
                        relsummary.group(gloc).event(eloc).betsd.m = ...
                            mean(data.g(gloc).e(eloc).sig_id.raw);
                        relsummary.group(gloc).event(eloc).betsd.ll = ...
                            quantile(data.g(gloc).e(eloc).sig_id.raw,ciedge);
                        relsummary.group(gloc).event(eloc).betsd.ul = ...
                            quantile(data.g(gloc).e(eloc).sig_id.raw,1-ciedge);
                        
                        relsummary.group(gloc).event(eloc).witsd.m = ...
                            mean(data.g(gloc).e(eloc).sig_err.raw);
                        relsummary.group(gloc).event(eloc).witsd.ll = ...
                            quantile(data.g(gloc).e(eloc).sig_err.raw,ciedge);
                        relsummary.group(gloc).event(eloc).witsd.ul = ...
                            quantile(data.g(gloc).e(eloc).sig_err.raw,1-ciedge);
                        
                    end
                end
                
                
        end %switch analysis

    case 'trt_diff'
        %Two-facet (test-retest) difference-score reliability. Delegated to a
        %helper, mirroring 'sing_diff_sserr'. Computes the per-condition
        %test-retest reliability (trial cutoffs at the occasion estimand) and the
        %difference-score reliability via psyrat_diffrel_trt (CE/CS/CES x
        %{gen,dep} + ICC), applying the single- vs multi-occasion n'_o convention.
        [relsummary,data,relerr] = psyrat_relsummary_trt_diff(...
            psyrat_data,relcutoff,meascutoff,relcentmeas,gcoeff,diffgcoeff,...
            reltype,noccmode,noccarg,ciperc);

end

%publish contract (RC-44): relsummary.group, and the relsummary.data.g stored on
%the next line, must each carry exactly one element per REL.groups entry, in
%REL.groups order. Consumers split on how they size their group loop, and that
%split is why the contract matters. Most size from REL.groups and index the
%published array with no bounds check (23 such loops across 13 files, among them
%psyrat_variancet (5), psyrat_report (3) and psyrat_ssrelplot (2)), so a short
%publish throws there. The count was 26 when this comment was written; a re-count
%on 2026-07-31 found psyrat_trt_reloverallt contributes ONE such loop (:120), not
%four -- its other three group loops index the LOCAL newrelsummary rebuild
%(:756, :811, :843), not the published array. Do not re-derive this from a bare
%grep for "for gloc": the discriminator is which struct the loop body indexes.
%A few size from the published array instead and are deliberately
%fail-open, so they would drop a group in silence: psyrat_admissibility_collect
%(:57), reached from the last statement of this function, and
%psyrat_startview_criterion (:347), which sits inside a try/catch. A short
%publish would therefore crash some surfaces and quietly omit a group from
%others. It holds at HEAD structurally rather than by enforcement: the group
%loops here are bounded by ngroups (REL.groups under the same
%strcmpi(...,'none') test the consumers use) or by the matching REL.out label
%array, this file has no continue statements, and the ten relerr.nogooddata
%early returns all sit above this line. Because that is an observation about the
%current branches and not a guard, the contract is also pinned by
%testRelsummaryPublishesOneGroupPerRelGroup.
%store relsummary information
psyrat_data.relsummary = relsummary;
psyrat_data.relsummary.data = data;

%Collect any admissibility records the coefficient kernels attached to the
%difference scores in this summary into one note. Empty when nothing is wrong,
%so a clean run is byte-identical to before this diagnostic existed. The note
%is rebuilt every time the summary is, which is the right lifetime: it
%describes the numbers currently on screen, and a coefficient is not persisted
%in the .psyrat file either. Consumers surface it beside the numbers -- see the
%gbias_note precedent in psyrat_splits_summary.
psyrat_data.relsummary.admissibility_note = ...
    psyrat_admissibility_collect(psyrat_data.relsummary);

end

function [relsummary,data,relerr] = psyrat_relsummary_sing_diff_sserr(...
    psyrat_data,depcutoff,meascutoff,depcentmeas,diffgcoeff,ciperc)
REL = psyrat_data.rel;
data = struct;
relerr = struct('trlcutoff',0,'trlmax',0,'nogooddata',0);

if strcmpi(REL.groups,'none')
    ngroups = 1;
    gnames = {''};
else
    ngroups = length(REL.groups);
    gnames = cellstr(string(REL.groups(:)));
end

enames = cellstr(string(REL.events(:)));
if numel(enames) ~= 2
    error('relsummary:diffEvents', ...
        'Subject-level difference-score summaries require exactly two events.');
end

relsummary.depcutoff = depcutoff;
relsummary.ciperc = ciperc;
% The per-event subject-level tables follow the user's difference-score
% coefficient choice (diffgcoeff), so gcoeff/gcoeff_name track it rather than
% being hardcoded to dependability. The headline difference-score coefficient
% and inclusion decisions run on the separate diff_est path below and are
% unaffected.
relsummary.gcoeff = diffgcoeff;
relsummary.diffgcoeff = diffgcoeff;
if diffgcoeff == 1
    diff_est = 'dep';
    relsummary.gcoeff_name = 'dep';
    relsummary.diffgcoeff_name = 'dep';
else
    diff_est = 'gen';
    relsummary.gcoeff_name = 'gen';
    relsummary.diffgcoeff_name = 'gen';
end

switch meascutoff
    case 1
        relsummary.meascutoff = strcat('Lower Limit of ',...
            sprintf(' %2.0f',ciperc*100),'% Credible Interval');
    case 2
        relsummary.meascutoff = 'Point Estimate';
    case 3
        relsummary.meascutoff = strcat('Upper Limit of ',...
            sprintf(' %2.0f',ciperc*100),'% Credible Interval');
end

for gloc = 1:ngroups
    relsummary.group(gloc).name = gnames{gloc};

    idvar = psyrat_relsummary_array(REL.out.id_varcov{:,gloc});
    trlvar = psyrat_relsummary_array(REL.out.trl_varcov{:,gloc});
    bsigma = psyrat_relsummary_array(REL.out.b_sigma{:,gloc});
    er_ss1 = psyrat_relsummary_array(REL.out.er_var_ss1{:,gloc});
    er_ss2 = psyrat_relsummary_array(REL.out.er_var_ss2{:,gloc});
    wp_cov_ss = psyrat_relsummary_array(REL.out.wp_cov_ss{:,gloc});
    id_match = REL.out.id_matches{:,gloc};
    groupdata = psyrat_relsummary_groupdata(REL,gnames,gloc);

    err_varcov = [];
    if isfield(REL.out,'err_varcov') && size(REL.out.err_varcov,2) >= gloc && ...
            ~isempty(REL.out.err_varcov{:,gloc})
        err_varcov = psyrat_relsummary_array(REL.out.err_varcov{:,gloc});
    end
    if isempty(err_varcov)
        wp_cov = zeros(size(bsigma,1),1);
    else
        wp_cov = squeeze(err_varcov(:,2,1));
    end

    event_tables = cell(1,2);
    for eloc = 1:2
        if eloc == 1
            er_ss = er_ss1;
        else
            er_ss = er_ss2;
        end

        data.g(gloc).e(eloc).sd_id.raw = num2cell(sqrt(squeeze(idvar(:,eloc,eloc))));
        data.g(gloc).e(eloc).b_sigma.raw = num2cell(bsigma(:,eloc));
        data.g(gloc).e(eloc).id_varcov.raw = num2cell(squeeze(idvar(:,eloc,eloc)));
        data.g(gloc).e(eloc).trl_varcov.raw = num2cell(squeeze(trlvar(:,eloc,eloc)));
        data.g(gloc).e(eloc).er_var_ss.raw = er_ss;
        data.g(gloc).e(eloc).elabel = enames(eloc);
        data.g(gloc).e(eloc).id_match = id_match;
        data.g(gloc).glabel = gnames(gloc);

        idtable = psyrat_relsummary_event_idtable(groupdata,id_match,enames{eloc});
        wp_ss = er_ss - repmat(bsigma(:,eloc),1,size(er_ss,2));
        ssrel_table = psyrat_ssrel(...
            'bp',data.g(gloc).e(eloc).sd_id.raw,...
            'wp_pop',bsigma(:,eloc),...
            'wp_ss',num2cell(wp_ss),...
            'i',sqrt(squeeze(trlvar(:,eloc,eloc))),...
            'gcoeff',diffgcoeff,...
            'idtable',idtable,...
            'CI',ciperc);
        ssrel_table = psyrat_relsummary_mark_inclusion(ssrel_table,depcutoff,meascutoff);
        event_tables{eloc} = ssrel_table;

        ev = struct();
        ev.name = enames{eloc};
        ev.ssrel_table = ssrel_table;
        ev = psyrat_relsummary_event_fields(...
            ev,ssrel_table,depcutoff,meascutoff,...
            data.g(gloc).e(eloc).sd_id.raw,bsigma(:,eloc),ciperc,...
            sqrt(squeeze(trlvar(:,eloc,eloc))),diffgcoeff);
        relsummary.group(gloc).event(eloc) = ev;
    end

    diff_idtable = psyrat_relsummary_diff_idtable(groupdata,id_match,enames);
    diff_table = psyrat_ssrel_diff(...
        'bp',idvar,...
        'bt',trlvar,...
        'er_var',bsigma,...
        'wp_cov',wp_cov,...
        'idtable',diff_idtable,...
        'CI',ciperc,...
        'est',diff_est,...
        'er_var_ss1',er_ss1,...
        'er_var_ss2',er_ss2,...
        'wp_cov_ss',wp_cov_ss);
    diff_table = psyrat_relsummary_mark_inclusion(diff_table,depcutoff,meascutoff);

    relsummary.group(gloc).goodids = diff_table.id(diff_table.ind2include);
    relsummary.group(gloc).badids = diff_table.id(diff_table.ind2exclude);
    if isempty(relsummary.group(gloc).goodids)
        relsummary.group(gloc).goodids = cell(0,1);
        %G48 (owner ruling 2026-08-29: FULL alignment with the G39/G45
        %strict-fold conventions). This used to set relerr.nogooddata = 1,
        %which raised the modal threshold dialog after the loop and made
        %every caller abort with no output -- and because the flag is
        %run-global, ONE empty group blocked the output for every group.
        %Warn on the console per group and continue instead, exactly like
        %the strict-fold intersection loops and the G45 trt_diff helper:
        %the overall table renders this group as a complementing
        %n Included 0 row (psyrat_depoverallt's G44 dash-out arm) and the
        %per-event rows stay informative. badids already carries the full
        %roster here -- the marker derives exclusion as ~inclusion -- so
        %the complement needs no separate assignment. Deliberately NOT
        %local_relsummary_notice: its modal contract requires
        %abort-and-return, and this is not an abort.
        fprintf(['\nNo participants met the reliability threshold on ',...
            'the difference score for a group; that group''s ',...
            'n Included is 0.\n']);
    end

    included = diff_table.ind2include;
    if any(included)
        switch depcentmeas
            case 1
                obs = [mean(diff_table.trls1(included)) mean(diff_table.trls2(included))];
            case 2
                obs = [median(diff_table.trls1(included)) median(diff_table.trls2(included))];
        end
    else
        %G60 (owner ruling 2026-08-30: convention stands, documented in the
        %ledger). With nobody diff-included, the design point falls back to
        %the FULL table's central trial counts -- the same convention the
        %ic_diff sibling applies (its goodids substitution feeds the same
        %fallback through its depcentmeas switch) and trt_diff never departs
        %from. The coefficient itself is draw-based and stays defined (the
        %G44 kept-coefficient rule); only its design point degrades to the
        %full sample, which manual ch. 9's zero-retention paragraph now
        %states. The switch below is the G60 fix proper: this arm used to
        %hardcode mean, silently ignoring a median request, which was the
        %one real divergence from the sibling convention.
        switch depcentmeas
            case 1
                obs = [mean(diff_table.trls1) mean(diff_table.trls2)];
            case 2
                obs = [median(diff_table.trls1) median(diff_table.trls2)];
        end
    end

    diffscore = psyrat_diffrel(...
        'bp',idvar,...
        'bt',trlvar,...
        'er_var',bsigma,...
        'wp_cov',wp_cov,...
        'obs',obs,...
        'est',diff_est,...
        'CI',ciperc);
    diffscore.ssrel_table = diff_table;
    %The cutoff table's diff row on this route deliberately shows the
    %OBSERVED-design coefficient (this analysis runs no trial-cutoff search;
    %psyrat_depcutofft's placeholder contract prints '---' for the absent
    %ntrials). Store only the three fields that row reads. This used to be a
    %full self-copy (diffscore.trlcutoff = diffscore), which dragged the
    %admissibility record along and made psyrat_admissibility_collect gather
    %the IDENTICAL record twice -- its second gather exists for routes whose
    %cutoff row is an independently evaluated difference score, and its
    %local_gather skips a struct with no admissibility field, which is
    %exactly what carrying pt/ll/ul alone guarantees (G60 wave).
    diffscore.trlcutoff = struct('pt',diffscore.pt,'ll',diffscore.ll,...
        'ul',diffscore.ul);
    diffscore.trlinfo = struct('min',[min(diff_table.trls1) min(diff_table.trls2)],...
        'max',[max(diff_table.trls1) max(diff_table.trls2)],...
        'mean',[mean(diff_table.trls1) mean(diff_table.trls2)],...
        'med',[median(diff_table.trls1) median(diff_table.trls2)],...
        'std',[std(diff_table.trls1) std(diff_table.trls2)]);
    relsummary.group(gloc).diffscore = diffscore;
end

%RC-40 accounting, re-derived at G48 (owner ruling 2026-08-29). NO site in
%this file sets nogooddata outside the ten inline guards any more -- the
%accounting is ten setters total, the ten inline guards. It was eleven until
%G48: this helper set the flag inside its loop when a group's
%difference-score inclusion came up empty and raised the one-per-run modal
%that used to sit at the end of this function, so every caller aborted with
%no output (one empty group blocked every group). G45 had removed the
%twelfth (the trt_diff helper) the same way. The empty state now warns on
%the console per group inside the loop and the summary returns in full; the
%G44 dash-out arm in psyrat_depoverallt renders the zeroed group's row. A
%reviewer who greps isempty( will still find the guard feeding the
%per-group notice above; it belongs to this helper, not to the ten.
%
%The ten inline guards raise this same dialog and then return, and they do NOT
%share one reachability story -- an earlier draft of this comment claimed all
%ten were unreachable, and that universal is false. Adjudicated 2026-08-18
%against profiler coverage rather than a reading of the guards: SIX fire, FOUR
%do not, and the four do not all fail for the same reason. THE PARAGRAPHS UP TO
%THE SUPERSEDED MARKER ARE A HISTORICAL RECORD: their line numbers are as of
%the 2026-08-18 commit (G39's +~280 lines shifted everything), and the
%mechanism they describe includes conditions later removed (the CHAR 'none'
%carrier went under B34; the conditional fold went under G39). Re-derive with
%grep before trusting any of it.
%
%  LIVE (6, pre-G39): :1250, :1482, :1866, :2106, :4120, :4402. Each was fed by
%  a cross-event aggregation loop that did "badids = [];"
%  (:1204/:1426/:1818/:2048/:4074/:4342) and appended ONLY inside
%  the eventgoodids arm test. When no event found a cutoff that
%  conditional never opened, so tempids stayed empty, goodids became the empty
%  carrier (historically the CHAR 'none'; cell(0,1) since B34), and badids
%  stayed unique([]) = []. The arm test was then
%  false, arm 2 builds table([]) -- a 0x1 table, isempty true -- and the guard
%  FIRES. Strongest evidence is at :1482, executed by the committed fixture
%  makeICRelDataMultiGroupEvent at depcutoff 0.9999. The other five were driven
%  under synthesized REL structs with a positive control at a reachable cutoff
%  (0.30) that flipped every indicator, so they are SHOWN EXECUTABLE, not shown
%  reachable from psyrat_start -- no GUI or preflight path was exercised. :4120
%  is the thinnest of the six; treat it as the one to re-measure first.
%
%  DEAD, STRUCTURALLY (2): :728 and :3454. Their isempty(badids) test is not
%  "evaluated and always false" -- it is NEVER EVALUATED. Both sit inside a
%  SECOND "if isempty(trlcutoff)" (:701, :3419) that cannot be entered,
%  because the first arm three lines above reassigns trlcutoff = -1 (:685,
%  :3401) and nothing between re-empties the local. After that if/elseif
%  closes, trlcutoff is either -1 or a non-empty find() index, so control
%  always falls to the "elseif ~isempty(trlcutoff)" sibling (:769, :3503).
%  Measured hits = 0 on both.
%
%  DEAD, BY ARM (2): :1030 and :3862. These are genuinely evaluated, but arm 2
%  is reached only when every participant was excluded, which is exactly when
%  badids holds all of them. The long comment above the :3862 guard works this
%  through. Measured at :3859: evaluated twice, notice raised zero times.
%
%  SUPERSEDED 2026-08-29 (G39 strict fold): the LIVE (6) mechanism above no
%  longer exists. The cross-event loops now append EVERY event's ids -- a dead
%  event folds in its full eventbadids -- so when nothing clears the cutoff
%  the group badids hold every participant, arm 2 builds a NON-empty fallback
%  table, and those six guards are expected DEAD for nevents >= 1 (kept as
%  last-resort aborts). The user-facing report for zero retention is now the
%  per-group console notice beside each intersection loop plus the n Included
%  0 rows; a dead event or dead group no longer aborts the summary. The
%  2026-08-18 profiler record above is kept as the pre-G39 baseline; re-derive
%  before trusting either.
%
%  SUPERSEDED (G54): the DEAD, STRUCTURALLY (2) status above no longer holds.
%  The record correctly derived that the second "if isempty(trlcutoff)" could
%  never be entered, but stopped short of the consequence: the sibling arm the
%  control fell to compared GroupCount >= -1, so an unreachable cutoff marked
%  EVERY participant good in sing case 1 and trt case 1 -- the opposite of
%  cases 2-4 in both families and the trt_diff helper's sentinel arm
%  (trt_sserr has no cutoff search and no sentinel, so it does not belong
%  in this comparison). G54 respelled both
%  second guards on the sentinel VALUE (trlcutoff == -1), reviving the
%  everyone-bad arms. The two isempty(badids) guards inside them are therefore
%  no longer never-evaluated: they now run whenever the sentinel path runs,
%  and they join the never-firing backstop class -- on that path badids holds
%  the entire roster, so the test is false whenever any data exist. The
%  never-evaluated record above is the pre-G54 baseline.
%
%Exactly ONE committed test is load-bearing on any of this, not two:
%testTrtExtrapolatedCutoffWithNoRetainedIdsStillWritesGoodN reaches the trt
%analysis-2 fallback-arm guard and asserts goodn and trlmax are written
%afterward, so that site must stay dead. Line numbers for it are NOT
%restated here on purpose -- earlier revisions of this paragraph carried
%coordinates that went stale twice (adversarial wave, 2026-08-29); the
%current, anchored record lives IN the comment directly above that guard
%(grep for "ONE committed test depends on it staying dead"). The second
%test the older comment named, testCriterionScreensReportTrialCutoffFailure,
%executes NONE of the ten guard sites (measured 2026-08-18: hits = 0 across
%its three psyrat_relsummary calls). Do not repeat the "two tests" claim.
%
%No nogooddata block here any more (G48, owner ruling 2026-08-29). This
%helper was the last outside-the-guards setter and used to raise the
%one-per-run modal notice here (until then the run aborted at
%psyrat_relfigures' nogooddata gate with no figures and no tables); the
%empty state now warns on the console per group inside the loop and the
%summary is returned in full, aligning this branch with the strict-fold
%intersection loops and the G45 trt_diff helper. relerr.nogooddata stays 0
%on every path out of this function.
end

function arr = psyrat_relsummary_array(v)
if iscell(v)
    arr = cell2mat(v);
else
    arr = v;
end
end

function groupdata = psyrat_relsummary_groupdata(REL,gnames,gloc)
groupdata = REL.data;
if any(strcmpi(groupdata.Properties.VariableNames,'group')) && ...
        ~isempty(gnames{gloc}) && ~strcmpi(gnames{gloc},'none')
    groupdata = groupdata(strcmp(string(groupdata.group),string(gnames{gloc})),:);
end
end

function idtable = psyrat_relsummary_event_idtable(groupdata,id_match,eventname)
trltable = varfun(@length,groupdata,'GroupingVariables',{'id','event'});
trltable = trltable(strcmp(string(trltable.event),string(eventname)),:);
idtable = innerjoin(id_match,trltable(:,[1 3]));
idtable.Properties.VariableNames{3} = 'trls';
end

function idtable = psyrat_relsummary_diff_idtable(groupdata,id_match,enames)
trltable = varfun(@length,groupdata,'GroupingVariables',{'id','event'});
tab1 = trltable(strcmp(string(trltable.event),string(enames{1})),[1 3]);
tab2 = trltable(strcmp(string(trltable.event),string(enames{2})),[1 3]);
tab1.Properties.VariableNames = {'id','trls1'};
tab2.Properties.VariableNames = {'id','trls2'};
idtable = innerjoin(tab1,tab2,'Keys','id');
idtable = innerjoin(id_match,idtable,'Keys','id');
idtable.trls = min(idtable.trls1,idtable.trls2);
end

function out = psyrat_relsummary_mark_inclusion(out,depcutoff,meascutoff)
switch meascutoff
    case 1
        out.ind2include = out.dep_ll >= depcutoff;
    case 2
        out.ind2include = out.dep_pt >= depcutoff;
    case 3
        out.ind2include = out.dep_ul >= depcutoff;
end
out.ind2exclude = ~out.ind2include;
end

function out = local_default_sigtrl(out, shapeRef)
%Backward compatibility: result structs produced before the one-facet trial
%main effect (sigma_i) was modeled carry no sig_trl. Default it to zeros of the
%given shape so the coefficient wiring reduces to the prior behavior (trial
%term = 0, generalizability == dependability) instead of erroring. shapeRef is
%the model-appropriate component whose shape sig_trl must match: sig_e for the
%one-facet 'sing' model, pop_sdlog for the subject-level 'sserr' model.
%
%sig_trl_source records WHICH of those two cases produced the value, because a
%backfilled zero and a genuinely estimated near-zero are indistinguishable once
%they reach the coefficient. When the trial term is zero, dependability equals
%generalizability exactly, so an absolute-error label on this path is unearned;
%the flag is what lets the provenance line say so (RC-15). It is only set here
%when it is not already carried by the producer.
if ~isfield(out,'sig_trl') || isempty(out.sig_trl)
    out.sig_trl = zeros(size(shapeRef));
    if ~isfield(out,'sig_trl_source') || isempty(out.sig_trl_source)
        out.sig_trl_source = 'legacy_default';
    end
elseif ~isfield(out,'sig_trl_source') || isempty(out.sig_trl_source)
    out.sig_trl_source = 'estimated';
end
end

function [noccmode, noccarg] = local_parse_nocc(varargin)
%Shared parse of the test-retest occasion estimand inputs, used identically by
%the 'trt' and 'trt_diff' branches so the two cannot silently diverge.
%noccmode: 1 = single-occasion score (n'_o = 1); 2 = multi-occasion composite
%(n'_o = k). noccarg: the composite size k (used only when noccmode == 2; [] =
%default to the number of observed occasions, resolved later once REL is known).
ind = find(strcmpi('noccmode',varargin),1);
if ~isempty(ind)
    noccmode = varargin{ind+1};
    if ~any(noccmode == [1 2])
        error('varargin:noccmode',... %Error code and associated error
            strcat('WARNING: noccmode not specified \n\n',...
            'Please input noccmode as:\n',...
            '1 - single-occasion score (n''_o = 1)\n',...
            '2 - multi-occasion composite score (n''_o = k)\n',...
            'See help psyrat_relsummary for more information \n'));
    end
else
    noccmode = 1;
end

ind = find(strcmpi('nocc',varargin),1);
if ~isempty(ind)
    noccarg = varargin{ind+1};
    if length(noccarg) ~= 1 || noccarg <= 0 || mod(noccarg,1) ~= 0
        error('varargin:nocc',... %Error code and associated error
            strcat('WARNING: nocc should be a positive integer\n',...
            'Please input nocc (number of occasions in the composite).\n',...
            'See help psyrat_relsummary for more information \n'));
    end
else
    noccarg = [];
end
end

function [effnocc, nocc_name, nobservedocc] = local_resolve_effnocc(noccmode, noccarg, REL)
%Shared resolution of the effective occasion estimand n'_o, used identically by
%the 'trt' and 'trt_diff' branches (audit sec.14, finding F2). Single-occasion
%-> 1; multi-occasion composite -> the requested k, or the number of observed
%occasions when k was omitted. The caller stores the values into its own
%relsummary struct.
nobservedocc = 1;
if isfield(REL,'time') && ~isempty(REL.time)
    nobservedocc = max(1,numel(REL.time));
end
if noccmode == 2
    if isempty(noccarg)
        effnocc = nobservedocc;
    else
        effnocc = noccarg;
    end
    nocc_name = 'multi';
else
    effnocc = 1;
    nocc_name = 'single';
end
end

function ssrel_table = local_ssrel_trt_call(d, gcoeff, reltype, idtable, ciperc, effnocc)
%Build one subject-level test-retest ssrel_table (trt_sserrvar, analysis 25) for
%a single group x event chunk. The person block (gro_sds/pop_sdlog/ind_sdlog) is
%stored cell-wrapped by the estimator, so convert to numeric before calling
%psyrat_ssrel_trt (which, unlike psyrat_ssrel, expects numeric draws). Only the
%residual is per-subject: err_ss = exp(pop_sdlog + ind_sdlog) is the per-subject
%residual SD, while the crossed facet SDs are population, so txp_ss/oxp_ss simply
%broadcast the population person x trial / person x occasion SDs across subjects
%(psyrat_ssrel_trt's per-subject path requires all three *_ss inputs).
grosds = cell2mat(d.gro_sds);            % draws x 2 (person intercept SD, log-resid SD)
indsdlog = cell2mat(d.ind_sdlog);        % draws x NSUB (per-subject log-residual)
err_ss = exp(d.pop_sdlog + indsdlog);    % draws x NSUB (per-subject residual SD)
nsub_e = size(err_ss,2);
ssrel_table = psyrat_ssrel_trt(...
    'gcoeff',gcoeff,...
    'reltype',reltype,...
    'bp',grosds(:,1),...
    'bo',d.sig_occ,...
    'bt',d.sig_trl,...
    'txp',d.sig_trlxid,...
    'oxp',d.sig_occxid,...
    'txo',d.sig_trlxocc,...
    'err',exp(d.pop_sdlog),...
    'txp_ss',repmat(d.sig_trlxid,1,nsub_e),...
    'oxp_ss',repmat(d.sig_occxid,1,nsub_e),...
    'err_ss',err_ss,...
    'idtable',idtable,...
    'CI',ciperc,...
    'nocc',effnocc);
end

function relsummary = local_store_ss_chunk(relsummary, ssrel_table, gloc, eloc, meascutoff, depcutoff)
%Apply the participant-inclusion cutoff to a subject-level ssrel_table and store
%it plus the per-chunk good/bad ids and trial-count summary (mirrors the storage
%block of the single-occasion sserr case). Shared by the trt_sserr group/event
%loop above. The group-level goodids (intersection across events) is computed by
%the caller after the loop.
%G49: inclusion and exclusion used to come from independent >=/< comparisons,
%which a NaN dependability draw fails BOTH of, so a partial-NaN subject landed
%in NEITHER carrier and the good/bad complement broke. The shared marker
%derives exclusion as ~inclusion, so NaN reads as excluded.
ssrel_table = psyrat_relsummary_mark_inclusion(ssrel_table,depcutoff,meascutoff);

relsummary.group(gloc).event(eloc).ssrel_table = ssrel_table;

if all(ssrel_table.ind2exclude)
    %no participant met the cutoff -- including the all-NaN failed-fit
    %chunk, which G49's complement folds into this arm; summarize trials
    %over everyone
    relsummary.group(gloc).event(eloc).eventgoodids = cell(0,1);
    relsummary.group(gloc).event(eloc).eventbadids = ssrel_table.id(:);
    relsummary.group(gloc).event(eloc).trlinfo = ...
        local_ss_trlinfo(ssrel_table.trls);
    relsummary.group(gloc).event(eloc).goodn = 0;
elseif any(ssrel_table.ind2include)
    relsummary.group(gloc).event(eloc).eventgoodids = ...
        ssrel_table.id(ssrel_table.ind2include);
    relsummary.group(gloc).event(eloc).eventbadids = ...
        ssrel_table.id(ssrel_table.ind2exclude);
    relsummary.group(gloc).event(eloc).trlinfo = ...
        local_ss_trlinfo(ssrel_table.trls(ssrel_table.ind2include));
    relsummary.group(gloc).event(eloc).goodn = ...
        height(ssrel_table(ssrel_table.ind2include,:));
else
    %G46, now a backstop. Since G49 derives ind2exclude as ~ind2include, an
    %all-NaN chunk (nobody includable) makes all(ind2exclude) true and takes
    %the first arm, whose assignments are identical -- this arm is
    %unreachable unless the indicator convention diverges again. Kept as the
    %documented dead-event producer for that state: canonical pair, nobody
    %retained, full roster bad.
    relsummary.group(gloc).event(eloc).eventgoodids = cell(0,1);
    relsummary.group(gloc).event(eloc).eventbadids = ssrel_table.id(:);
    relsummary.group(gloc).event(eloc).trlinfo = ...
        local_ss_trlinfo(ssrel_table.trls);
    relsummary.group(gloc).event(eloc).goodn = 0;
end
end

function ti = local_ss_trlinfo(trls)
%Trial-count summary block shared by the subject-level storage helper.
ti.min = min(trls);
ti.max = max(trls);
ti.mean = mean(trls);
ti.std = std(trls);
ti.med = median(trls);
end

function wp = local_total_within(sig_e_raw, sig_trl_raw)
%Total within-person SD for the one-facet model: combine the person x trial,e
%residual (sig_e) and the trial main effect (sig_trl) into
%wp = sqrt(sig_e^2 + sig_trl^2). psyrat_rel_sing squares this back to
%wp^2 = sigma_pi,e^2 + sigma_i^2 and subtracts i^2 = sigma_i^2 to recover the
%relative residual sigma_pi,e^2, so generalizability excludes the trial main
%effect while dependability keeps it via wp^2 (Rocha Table 2 / S10/F1). Defining
%this once means the within-person error formula lives in a single place instead
%of being inlined identically at every psyrat_rel_sing call.
%
%CONTRACT, because the two arguments are easy to transpose and the sibling
%subject-level function uses the OPPOSITE convention: psyrat_rel_sing wants the
%TOTAL here and subtracts i^2, whereas psyrat_ssrel wants the RELATIVE residual
%and ADDS i^2. Passing 'i' to psyrat_rel_sing WITHOUT also widening wp to the
%total is therefore not a partial fix but a regression: dependability is left
%unchanged and generalizability has sigma_i^2 subtracted from a quantity that
%never contained it. Always change the pair together (RC-02/RC-19).
%
%Used by both the one-facet 'sing' branch (sig_e / sig_trl draws) and the
%single-session difference 'sing_diff' per-event branch (exp(b_sigma) / sd_trl
%draws); the two differ only in how the caller names and unpacks the draws.
wp = sqrt(sig_e_raw.^2 + sig_trl_raw.^2);
end

function ev = psyrat_relsummary_event_fields(ev,ssrel_table,depcutoff,meascutoff,bp,wp_pop,ciperc,sig_trl,gcoeff)
ev.trlcutoff = min(ssrel_table.trls(ssrel_table.ind2include));
if isempty(ev.trlcutoff)
    ev.trlcutoff = min(ssrel_table.trls);
end
ev.eventgoodids = ssrel_table.id(ssrel_table.ind2include);
ev.eventbadids = ssrel_table.id(ssrel_table.ind2exclude);
if isempty(ev.eventgoodids)
    ev.eventgoodids = cell(0,1);
end
ev.goodn = sum(ssrel_table.ind2include);
if any(ssrel_table.ind2include)
    trls = ssrel_table.trls(ssrel_table.ind2include);
else
    trls = ssrel_table.trls;
end
ev.trlinfo = struct('min',min(trls),'max',max(trls),'mean',mean(trls),...
    'med',median(trls),'std',std(trls));

% Total within-person SD includes the trial main effect (sigma_i), so the
% absolute-error (dependability, gcoeff=1) coefficient carries it;
% psyrat_rel_sing recovers the relative-error (generalizability, gcoeff=2)
% quantity as wp^2 - i^2. Mirrors the standalone one-facet sserr path.
wp_total = sqrt(exp(wp_pop).^2 + sig_trl(:).^2);
[ev.dep.ll,ev.dep.m,ev.dep.ul] = psyrat_rel_sing(...
    'gcoeff',gcoeff,'metric','global','bp',cell2mat(bp),...
    'wp',wp_total,'i',sig_trl,'obs',ev.trlinfo.mean,'CI',ciperc);
ev.dep.meas = ev.trlinfo.mean;
ev.rel = ev.dep;
[ev.icc.ll,ev.icc.m,ev.icc.ul] = psyrat_rel_sing(...
    'gcoeff',gcoeff,'metric','icc','bp',cell2mat(bp),...
    'wp',wp_total,'i',sig_trl,'CI',ciperc);

%component SDs share the caller's CI rather than a fixed 95%; see the main
%function for why the .95 default still reproduces the old .025/.975 output
ciedge = (1-ciperc)./2;
ev.betsd.m = mean(cell2mat(bp));
ev.betsd.ll = quantile(cell2mat(bp),ciedge);
ev.betsd.ul = quantile(cell2mat(bp),1-ciedge);
ev.witsd.m = mean(exp(wp_pop));
ev.witsd.ll = quantile(exp(wp_pop),ciedge);
ev.witsd.ul = quantile(exp(wp_pop),1-ciedge);
end

function [relsummary,data,relerr] = psyrat_relsummary_trt_diff(...
    psyrat_data,relcutoff,meascutoff,relcentmeas,gcoeff,diffgcoeff,...
    reltype,noccmode,noccarg,ciperc)
%Two-facet (test-retest) difference-score relsummary (analysis 'trt_diff').
%For each group it computes (a) per-condition test-retest reliability over
%trials (to locate a trial cutoff at the chosen occasion estimand) using
%psyrat_rel_trt with the per-condition components taken from the diagonals of
%the cross-condition 2x2 covariances, and (b) the difference-score reliability
%(CE/CS/CES x {gen,dep} + ICC) via psyrat_diffrel_trt using the full
%cross-condition covariances and the per-condition residual (+ residual
%covariance when concurrent). The single- vs multi-occasion n'_o convention is
%applied consistently (audit sec.14): n'_o defaults to 1 (single-occasion) and
%is threaded into both the per-condition coefficient and the difference score.
REL = psyrat_data.rel;
data = struct;
relerr = struct('trlcutoff',0,'trlmax',0,'nogooddata',0);

if strcmpi(REL.groups,'none')
    ngroups = 1;
    gnames = {''};
else
    ngroups = length(REL.groups);
    gnames = cellstr(string(REL.groups(:)));
end

enames = cellstr(string(REL.events(:)));
if numel(enames) ~= 2
    error('relsummary:diffEvents', ...
        'Two-facet difference-score summaries require exactly two events.');
end

%resolve the occasion estimand n'_o (shared helper; identical to the 'trt' branch)
[effnocc, nocc_name, nobservedocc] = local_resolve_effnocc(noccmode, noccarg, REL);

relsummary = struct();
relsummary.relcutoff = relcutoff;
%alias depcutoff to relcutoff: the difference-score tables (sing-family) read
%relsummary.depcutoff, while the test-retest path uses relcutoff.
relsummary.depcutoff = relcutoff;
relsummary.reltype = reltype;
if reltype == 1
    relsummary.reltype_name = 'ic';
elseif reltype == 2
    relsummary.reltype_name = 'trt';
else
    relsummary.reltype_name = 'ic_trt';
end
relsummary.ciperc = ciperc;

%component SDs share the caller's CI rather than a fixed 95%; see the main
%function for why the .95 default still reproduces the old .025/.975 output
ciedge = (1-ciperc)./2;
relsummary.gcoeff = gcoeff;
if gcoeff == 1
    relsummary.gcoeff_name = 'dep';
else
    relsummary.gcoeff_name = 'gen';
end
relsummary.diffgcoeff = diffgcoeff;
if diffgcoeff == 1
    diff_est = 'dep';
    relsummary.diffgcoeff_name = 'dep';
else
    diff_est = 'gen';
    relsummary.diffgcoeff_name = 'gen';
end
relsummary.noccmode = noccmode;
relsummary.nocc = effnocc;
relsummary.nobservedocc = nobservedocc;
relsummary.nocc_name = nocc_name;
switch meascutoff
    case 1
        relsummary.meascutoff = strcat('Lower Limit of ',...
            sprintf(' %2.0f',ciperc*100),'% Credible Interval');
    case 2
        relsummary.meascutoff = 'Point Estimate';
    case 3
        relsummary.meascutoff = strcat('Upper Limit of ',...
            sprintf(' %2.0f',ciperc*100),'% Credible Interval');
end

for gloc = 1:ngroups
    relsummary.group(gloc).name = gnames{gloc};

    %per-group cross-condition covariance draws + per-condition residual log-SD
    idvar  = cell2mat(REL.out.id_varcov{:,gloc});   % [draws x 2 x 2]
    tidvar = cell2mat(REL.out.tid_varcov{:,gloc});
    oidvar = cell2mat(REL.out.oid_varcov{:,gloc});
    trlvar = cell2mat(REL.out.trl_varcov{:,gloc});
    occvar = cell2mat(REL.out.occ_varcov{:,gloc});
    tovar  = cell2mat(REL.out.to_varcov{:,gloc});
    bsigma = cell2mat(REL.out.b_sigma{:,gloc});     % [draws x 2]

    %residual covariance: estimated (concurrent) or 0 (non-concurrent)
    er_cov = psyrat_trtdiff_ercov(REL,gloc,bsigma);

    %restrict the trial-count data to this group
    groupdata = REL.data;
    if any(strcmpi(groupdata.Properties.VariableNames,'group')) && ...
            ~isempty(gnames{gloc}) && ~strcmpi(gnames{gloc},'none')
        groupdata = groupdata(strcmp(string(groupdata.group),string(gnames{gloc})),:);
    end

    percell_by_event = cell(1,2);
    uid_by_event = cell(1,2);
    for eloc = 1:2
        c = eloc;
        sig_id  = sqrt(squeeze(idvar(:,c,c)));
        sig_tid = sqrt(squeeze(tidvar(:,c,c)));
        sig_oid = sqrt(squeeze(oidvar(:,c,c)));
        sig_trl = sqrt(squeeze(trlvar(:,c,c)));
        sig_occ = sqrt(squeeze(occvar(:,c,c)));
        sig_to  = sqrt(squeeze(tovar(:,c,c)));
        sig_err = exp(bsigma(:,c));

        %store raw component SD draws (for the variance-components viewer)
        data.g(gloc).e(eloc).sig_id.raw = sig_id;
        data.g(gloc).e(eloc).sig_occ.raw = sig_occ;
        data.g(gloc).e(eloc).sig_trl.raw = sig_trl;
        data.g(gloc).e(eloc).sig_trlxid.raw = sig_tid;
        data.g(gloc).e(eloc).sig_occxid.raw = sig_oid;
        data.g(gloc).e(eloc).sig_trlxocc.raw = sig_to;
        data.g(gloc).e(eloc).sig_err.raw = sig_err;
        data.g(gloc).e(eloc).b_sigma.raw = bsigma(:,c);
        data.g(gloc).e(eloc).elabel = enames(eloc);
        data.g(gloc).glabel = gnames(gloc);

        %per-occasion trial counts per participant (trials within id x time)
        eventdata = groupdata(strcmp(string(groupdata.event),string(enames{eloc})),:);
        trltable = varfun(@length,eventdata,'GroupingVariables',{'id','time'});
        uid = unique(trltable.id,'stable');
        percell = zeros(numel(uid),1);
        for k = 1:numel(uid)
            percell(k) = mean(trltable.GroupCount(ismember(trltable.id,uid(k))));
        end
        percell_by_event{eloc} = percell;
        uid_by_event{eloc} = uid;

        ev = struct();
        ev.name = enames{eloc};

        %per-condition reliability-vs-trials curve to locate the cutoff.
        %percell is a per-participant MEAN trial count across occasions, so it is
        %fractional when occasions are unbalanced (the norm after QC). ceil keeps
        %the obs range bound an integer (every other relsummary branch uses an
        %integer max(GroupCount)); a fractional bound makes psyrat_rel_trt build
        %zeros(1,obs(2)-obs(1)+1) with a non-integer size.
        ntrials = ceil(max(percell)) + 1000;
        [llrel,mrel,ulrel] = psyrat_rel_trt('gcoeff',gcoeff,'reltype',reltype,...
            'bp',sig_id,'bo',sig_occ,'bt',sig_trl,'txp',sig_tid,'oxp',sig_oid,...
            'txo',sig_to,'err',sig_err,'obs',[1 ntrials],'nocc',effnocc,'CI',ciperc);
        switch meascutoff
            case 1
                trlcutoff = find(llrel >= relcutoff,1);
            case 2
                trlcutoff = find(mrel >= relcutoff,1);
            case 3
                trlcutoff = find(ulrel >= relcutoff,1);
        end

        if isempty(trlcutoff)
            relerr.trlcutoff = 1;
            ev.trlcutoff = -1;
            ev.relcutoff = struct('m',-1,'ll',-1,'ul',-1);
            ev.eventgoodids = cell(0,1);
            ev.eventbadids = uid;
        else
            ev.trlcutoff = trlcutoff;
            ev.relcutoff = struct('m',mrel(trlcutoff),'ll',llrel(trlcutoff),...
                'ul',ulrel(trlcutoff));
            inc = percell >= trlcutoff;
            ev.eventgoodids = uid(inc);
            ev.eventbadids = uid(~inc);
            if isempty(ev.eventgoodids)
                ev.eventgoodids = cell(0,1);
            end
        end

        %trial info + the central-tendency reliability estimate
        ev.trlinfo = struct('min',min(percell),'max',max(percell),...
            'mean',mean(percell),'med',median(percell),'std',std(percell));
        switch relcentmeas
            case 1
                relcent = ev.trlinfo.mean;
            case 2
                relcent = ev.trlinfo.med;
        end
        [ev.rel.ll,ev.rel.m,ev.rel.ul] = psyrat_rel_trt('gcoeff',gcoeff,...
            'reltype',reltype,'bp',sig_id,'bo',sig_occ,'bt',sig_trl,...
            'txp',sig_tid,'oxp',sig_oid,'txo',sig_to,'err',sig_err,...
            'obs',relcent,'nocc',effnocc,'CI',ciperc);
        ev.rel.meas = relcent;
        %alias .dep to .rel: the difference-score tables (psyrat_depoverallt)
        %read .dep while the cutoff table (psyrat_depcutofft) reads .rel.
        ev.dep = ev.rel;

        %ICC = single-observation reliability (one trial at one occasion), so it
        %must NOT scale with the multi-occasion composite n'_o. Pin nocc=1 to
        %match the plain 'trt' branch and stay composite-invariant (otherwise a
        %noccmode=2 run divides the occasion variance by k and the same
        %condition's ICC differs between a trt_diff run and a trt run).
        [ev.icc.ll,ev.icc.m,ev.icc.ul] = psyrat_rel_trt('gcoeff',gcoeff,...
            'reltype',reltype,'bp',sig_id,'bo',sig_occ,'bt',sig_trl,...
            'txp',sig_tid,'oxp',sig_oid,'txo',sig_to,'err',sig_err,...
            'obs',1,'nocc',1,'CI',ciperc);

        ev.betsd = struct('m',mean(sig_id),'ll',quantile(sig_id,ciedge),...
            'ul',quantile(sig_id,1-ciedge));
        ev.witsd = struct('m',mean(sig_err),'ll',quantile(sig_err,ciedge),...
            'ul',quantile(sig_err,1-ciedge));

        relsummary.group(gloc).event(eloc) = ev;
    end

    %good/bad ids: a participant is good if it cleared the cutoff in BOTH events
    g1 = relsummary.group(gloc).event(1).eventgoodids;
    g2 = relsummary.group(gloc).event(2).eventgoodids;
    %both carriers are cell arrays (empty when nothing was retained), so
    %intersect is always defined and returns an empty cell in that case, which
    %the isempty(goodids) test below already handles.
    goodids = intersect(g1,g2);
    allids = unique(groupdata.id);
    if isempty(goodids)
        relsummary.group(gloc).goodids = cell(0,1);
        relsummary.group(gloc).badids = allids;
        %G45 (owner ruling 2026-08-29: FULL alignment with the G39
        %strict-fold conventions). This used to set relerr.nogooddata = 1,
        %which raised the modal threshold dialog after the loop and made
        %every caller abort with no output -- and because the flag is
        %run-global, ONE empty group blocked the output for every group.
        %Warn on the console per group and continue instead, exactly like
        %the nine strict-fold intersection loops: the tables render this
        %group as complementing n Included 0 rows and the per-event rows
        %stay informative. Deliberately NOT local_relsummary_notice -- its
        %modal contract requires abort-and-return, and this is not an
        %abort.
        fprintf(['\nNo participants met the reliability threshold in ',...
            'both events for a group; that group''s n Included is 0.\n']);
    else
        relsummary.group(gloc).goodids = goodids;
        relsummary.group(gloc).badids = setdiff(allids,goodids);
    end
    for eloc = 1:2
        %G45 rekey (same ruling): goodn is the GROUP-level retained count
        %(the both-events intersection), not the per-event count -- the same
        %rekey G39 applied at the three case-4 sites. psyrat_depoverallt
        %pairs goodn with the group-level badids count, so a per-event goodn
        %made the printed rows fail to complement under partial retention
        %(they printed retained-in-event alongside excluded-from-
        %intersection). The per-event carriers (eventgoodids/eventbadids)
        %keep the informative per-event partition.
        relsummary.group(gloc).event(eloc).goodn = ...
            numel(relsummary.group(gloc).goodids);
    end

    %difference-score reliability using the observed (central-tendency) trials
    switch relcentmeas
        case 1
            trls1 = relsummary.group(gloc).event(1).trlinfo.mean;
            trls2 = relsummary.group(gloc).event(2).trlinfo.mean;
        case 2
            trls1 = relsummary.group(gloc).event(1).trlinfo.med;
            trls2 = relsummary.group(gloc).event(2).trlinfo.med;
    end
    diffscore = psyrat_diffrel_trt('bp',idvar,'bpi',tidvar,'bpo',oidvar,...
        'bt',trlvar,'bo',occvar,'boi',tovar,'er_var',bsigma,'er_cov',er_cov,...
        'obs',[trls1 trls2],'nocc',[effnocc effnocc],'reltype',reltype,...
        'est',diff_est,'CI',ciperc);
    diffscore.trlinfo = struct(...
        'mean',[relsummary.group(gloc).event(1).trlinfo.mean ...
                relsummary.group(gloc).event(2).trlinfo.mean],...
        'med',[relsummary.group(gloc).event(1).trlinfo.med ...
               relsummary.group(gloc).event(2).trlinfo.med]);

    %difference-score reliability at the per-condition trial cutoffs
    tc1 = relsummary.group(gloc).event(1).trlcutoff;
    tc2 = relsummary.group(gloc).event(2).trlcutoff;
    if tc1 > 0 && tc2 > 0
        diffcut = psyrat_diffrel_trt('bp',idvar,'bpi',tidvar,'bpo',oidvar,...
            'bt',trlvar,'bo',occvar,'boi',tovar,'er_var',bsigma,'er_cov',er_cov,...
            'obs',[tc1 tc2],'nocc',[effnocc effnocc],'reltype',reltype,...
            'est',diff_est,'CI',ciperc);
    else
        diffcut = diffscore;
        if isfield(diffcut,'admissibility')
            %G60 wave: the fallback ALIASES the parent evaluation; strip the
            %admissibility record so the collector does not count it twice
            %(the found-cutoff arm's psyrat_diffrel_trt record is
            %independent and keeps its own). Same class as the
            %sing_diff_sserr self-copy fixed in the same batch.
            diffcut = rmfield(diffcut,'admissibility');
        end
    end
    diffscore.trlcutoff = diffcut;

    relsummary.group(gloc).diffscore = diffscore;

    %RC-31. Publish the difference-score component blocks so the D-study curve
    %can be drawn from them. These are the SAME variables just passed to
    %psyrat_diffrel_trt above, stored rather than re-derived, which is the point:
    %the curve and the headline coefficient are then guaranteed to differ only in
    %the trial counts projected, never in the components. Re-reading REL.out at
    %plot time would work for the six covariance blocks but NOT for er_cov, which
    %comes from the local psyrat_trtdiff_ercov (concurrent vs non-concurrent vs
    %reconstructed-from-rescor); copying that resolution into a plot function
    %would make a fourth undeclared twin of the kind CLAUDE.md warns about.
    %
    %Memory-only: relsummary is rebuilt at view time and is never written into
    %the .psyrat file (psyrat_startproc saves before any summary exists), so this
    %does not grow saved results. Roughly 1 MB per group at typical draw counts.
    relsummary.group(gloc).diffcomp = struct(...
        'bp',idvar,'bpi',tidvar,'bpo',oidvar,...
        'bt',trlvar,'bo',occvar,'boi',tovar,...
        'er_var',bsigma,'er_cov',er_cov);
end

%No nogooddata block here any more (G45, owner ruling 2026-08-29). This helper
%used to be the second of the two sites that set nogooddata outside the ten
%inline guards (a participant must clear the cutoff in BOTH events to count as
%good) and to raise the one-per-run modal notice on it; the empty-intersection
%state now warns on the console per group inside the loop and the summary is
%returned in full, aligning this branch with the nine strict-fold intersection
%loops. No outside-the-guards setter remains anywhere in this file: G48
%(owner ruling 2026-08-29) removed the last one, the sing_diff_sserr
%helper's -- see the re-derived accounting comment there (RC-40) before
%restating any reachability count.
end

function er_cov = psyrat_trtdiff_ercov(REL,gloc,bsigma)
%Resolve the per-draw residual covariance sigma_XY(poi,e) for group gloc:
%estimated from the concurrent (rescor) model, reconstructed from a residual
%correlation if only that is stored, or 0 for the non-concurrent model.
er_cov = [];
if isfield(REL,'diffrescor') && REL.diffrescor == 1
    er_cov = zeros(size(bsigma,1),1);
elseif isfield(REL.out,'err_varcov') && ~isempty(REL.out.err_varcov) && ...
        size(REL.out.err_varcov,2) >= gloc && ~isempty(REL.out.err_varcov{1,gloc})
    errvc = cell2mat(REL.out.err_varcov{:,gloc});
    er_cov = errvc(:,2,1);
elseif isfield(REL.out,'rescor') && ~isempty(REL.out.rescor) && ...
        size(REL.out.rescor,2) >= gloc && ~isempty(REL.out.rescor{1,gloc})
    rescor = cell2mat(REL.out.rescor{:,gloc});
    er_cov = rescor(:) .* exp(bsigma(:,1)) .* exp(bsigma(:,2));
end
if isempty(er_cov)
    er_cov = zeros(size(bsigma,1),1);
end
end

%local_relsummary_interactive was deleted here (G48 adversarial wave,
%2026-08-29): its only caller was the sing_diff_sserr helper's modal block,
%which G45's sibling had already removed from trt_diff and G48 removed from
%sing_diff_sserr. The main function resolves the interactive flag inline
%near its top and passes it to the ten guarded notice sites directly.

function local_relsummary_notice(interactive, dlg)
%Show the "data do not meet reliability threshold" notice only on interactive
%(GUI) runs. Headless/batch callers (psyrat_data.proc.interactive == false) must
%not block on a dialog; they detect the same condition through the unconditional
%relerr.nogooddata flag the caller sets immediately after this call returns.
%MODAL (2026-08-18 ruling; it was non-modal under the 2026-08-17 notification
%ruling, pending this call). Two facts make modal the right answer, both
%verified in this file rather than assumed:
%(1) This is the terminal event of a run that drew NOTHING, which puts it with
%the operation-failure reporters rather than with the post-output notices --
%the reason those stay non-modal, that a modal would lock the user out of the
%very output the notice describes, has no analogue when there is no output.
%(2) Nothing is created after this dialog, so the B20 burial mechanism cannot
%apply. EVERY site that raises this notice sets relerr.nogooddata and returns
%immediately, and both production callers abort before drawing anything:
%psyrat_relfigures returns on nogooddata ahead of every plot/table block, and
%psyrat_criterionfigures returns ahead of local_plotcriterion. Modal also
%answers the B21 risk that a non-modal box does not raise above an ALREADY-open
%viewer when the user re-runs with a stricter cutoff and would otherwise see
%nothing happen at all.
%
%NO uiwait here, deliberately, and the asymmetry with psyrat_startproc's B28
%advisory is principled rather than an oversight. There, execution continues on
%to build new windows, so only uiwait can guarantee acknowledgment before they
%exist. Here the caller returns immediately and draws nothing, so the dialog
%stays visible on its own and modal is sufficient -- and uiwait would block a
%CALCULATION function, hanging the two committed behavioral tests that drive
%this path live (TestSummaryAndPresentationFunctions). errordlg returns
%immediately whatever its CreateMode, so the headless contract below and the
%caller's control flow are both untouched; measured green.
%
%DO NOT restate a reachability count in this comment. An earlier draft asserted
%"the ten unconditional sites are RC-40-dead" and both halves were wrong: this
%file calls them the ten GUARDED sites, and they are not uniformly dead. The
%contradiction that used to be flagged here was adjudicated 2026-08-18 -- four
%are dead, six are live, and the split is worked through once, above
%psyrat_relsummary_sing_diff_sserr. That is the only place it is stated; keep
%it that way rather than paraphrasing a count here. The modality ruling below
%never depended on it either way, because every one of those sites returns
%before drawing regardless of whether it is reachable.
if interactive
    errordlg(dlg, 'Data do not meet reliability threshold','modal');
end
end
