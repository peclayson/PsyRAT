function R = psyrat_gamma_twofacet_restriction(design)
%PSYRAT_GAMMA_TWOFACET_RESTRICTION Disclosure text for the RESTRICTED gamma
% two-facet designs (analyses 14, 25 and 27 under the gamma family,
% location-scale and log-nu alike), one source for every surface.
%
%   R = psyrat_gamma_twofacet_restriction()
%   R = psyrat_gamma_twofacet_restriction(design)
%
% design - optional; the analysis number (14, 25 or 27) or its rel.analysis
%          label ('ic_dynrel_trt', 'trt_sserrvar', 'ic_dynrel_sserrvar_trt').
%          When given, R.dialog names that design; otherwise it names all three.
%
% Output fields:
%   R.designs     - the three rel.analysis labels the label covers. Exact
%                   match is load-bearing at every gate: the labels nest
%                   ('ic_dynrel' prefixes 'ic_dynrel_trt'; 'dynrel_trt' sits
%                   inside the difference variant 'ic_diff_dynrel_trt'), so a
%                   prefix or contains test would fire the caveat on a design
%                   it does not describe (see local_gamma_is_analysis in
%                   psyrat_provenance_lines).
%   R.label       - 'RESTRICTED', the maturity label word (CHANGELOG design
%                   maturity table; manual Section 5.4).
%   R.header      - the one-line export-header disclosure.
%   R.screen      - the short on-screen note for the viewers.
%   R.figure      - the short note on saved figures (the dynamic surface of
%                   14/27 and the subject plots of 25), which travel without
%                   the export header.
%   R.dialog      - cellstr for the Analyze-time notice (warndlg lines).
%   R.dialogtitle - the notice's window title.
%
% Why one function: the same finding is disclosed at selection (psyrat_exec's
% notice before the save dialog), on view (psyrat_startview_dynrel for 14/27,
% psyrat_startview_trt for 25), on saved figures (psyrat_restricted_figure_note)
% and in every export header
% (psyrat_provenance_lines). Separate copies of a 100-word sentence would drift;
% the gamma target captions already follow this rule (psyrat_gamma_target_lines).
% Callers apply their own gate (family = gamma AND one of R.designs); this
% function only owns the words.
%
% The finding (owner rulings 2-4, 2026-09-21; LS27): with few occasion levels
% (three or four in each reviewed capture of these designs), the population
% variance components and the pooled residual reported as posterior MEANS are
% not recovered at level. In the five re-evaluated scope captures the
% posterior mean of the occasion component was 2.4 to 12.1 times its planted
% value, of the highest-order remainder 2.2 to 7.2 times and of the pooled
% residual 1.1 to 1.9 times, with several 95% intervals excluding the truth;
% these ratios are on the variance scale, even where a viewer shows SDs. The
% retained analysis-27 location-scale capture, reviewed as a control, read
% 1.42, 1.66 and 1.02 for the same components, which is not a new pass. The
% two prospective analysis-27 fits put them at 60.1 and 2.45, 31.6 and 3.35,
% and 1.83 and 2.95, the second fit non-converged by the toolbox's own
% criterion. The pattern is consistent with raw-scale components that are
% exponentials of squared log-scale SDs whose posteriors are right-skewed with
% few occasions, so that their means are tail-dominated. The coefficients,
% means of bounded per-draw ratios, are far less affected: G and D within
% about .03 and .06 of the truth for analyses 14 and 25. For analysis 27 both
% G and D are biased downward in the converged fits (medians across
% participants of the posterior-mean coefficient minus its planted value:
% G about -.05 under location-scale at four occasions and -.10 under log-nu at
% three; D -.13 in the converged prospective fit, -.10 under log-nu and -.02
% in the retained capture; the non-converged second prospective fit gave
% D +.08), with participant ordering kept (r >= .86). Wherever production was compared with
% an independent enumeration it agreed to 2e-12 or better: implementation
% consistency, not recovery, calibration, identifiability or sampling
% validation. The five scope captures, the retained control and the first
% prospective fit converged. The one-facet gamma dynamic designs (11, 26) do
% not show the pattern (ruling 4) and are outside this label.
%
% See also PSYRAT_PROVENANCE_LINES, PSYRAT_STARTVIEW_DYNREL,
% PSYRAT_STARTVIEW_TRT, PSYRAT_STARTPROC.

% Copyright (C) 2016-2026 Peter E. Clayson
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

R = struct();
R.designs = {'ic_dynrel_trt','ic_dynrel_sserrvar_trt','trt_sserrvar'};
R.label = 'RESTRICTED';

%the three names as the manual and the analysis reference use them
names = struct('ic_dynrel_trt','analysis 14 (dynamic test-retest reliability)', ...
    'trt_sserrvar','analysis 25 (subject-level test-retest reliability)', ...
    'ic_dynrel_sserrvar_trt','analysis 27 (subject-level dynamic test-retest reliability)');

%resolve the optional design argument to one label, or '' for "all three"
which = '';
if nargin >= 1 && ~isempty(design)
    if isnumeric(design)
        switch design(1)
            case 14; which = 'ic_dynrel_trt';
            case 25; which = 'trt_sserrvar';
            case 27; which = 'ic_dynrel_sserrvar_trt';
        end
    elseif ischar(design) || isstring(design)
        d = char(string(design));
        if any(strcmp(d, R.designs)); which = d; end
    end
end

%---- export header (one line; CSV writers quote it) --------------------------
R.header = ['Design maturity: RESTRICTED (three gamma two-facet designs, analyses ' ...
    '14, 25 and 27, location-scale and log-nu alike; owner rulings of 2026-09-21). ' ...
    'With few occasion levels (three or four in each reviewed capture of these ' ...
    'designs) the population variance components and the pooled residual ' ...
    'reported here as posterior means are not recovered at level. In the five ' ...
    're-evaluated scope captures the posterior mean of the occasion component ' ...
    'was 2.4 to 12.1 times its planted value, of the highest-order remainder 2.2 ' ...
    'to 7.2 times and of the pooled residual 1.1 to 1.9 times; these ratios are ' ...
    'on the variance scale, even where a viewer shows standard deviations. ' ...
    'Several 95% intervals excluded the truth. The retained analysis-27 ' ...
    'location-scale capture, reviewed as a control, read 1.42, 1.66 and 1.02 ' ...
    'for the same three components, which is not a new pass. Two prospective ' ...
    'analysis-27 fits put them at 60.1 and 2.45, 31.6 and 3.35, and 1.83 and ' ...
    '2.95 (the second fit did not converge by the toolbox''s own criterion). ' ...
    'The generalizability and dependability coefficients stay within about .03 ' ...
    'and .06 of the truth for analyses 14 and 25. For analysis 27, G and D ' ...
    'are both biased downward in the converged fits; as medians across ' ...
    'participants of the posterior-mean coefficient minus its planted value, ' ...
    'G is about -.05 under location-scale at four occasions and -.10 under ' ...
    'log-nu at three, and D is -.13 in the converged prospective fit, -.10 ' ...
    'under log-nu and -.02 in the retained capture, while the non-converged ' ...
    'second prospective fit gave D +.08; participant ordering is kept ' ...
    '(r >= .86). The pattern is consistent with raw-scale ' ...
    'components that are exponentials of squared log-scale SDs whose ' ...
    'posteriors are right-skewed with few occasions, which makes their ' ...
    'posterior means tail-dominated. Production matched an independent ' ...
    'enumeration to 2e-12 or better wherever the two were compared; that shows ' ...
    'implementation consistency, not recovery, calibration, identifiability or ' ...
    'sampling validation. Report the coefficients with this caveat; reading the ' ...
    'component estimates against their intervals does not remedy the weak ' ...
    'evidence.'];

%---- viewer note (wrapped on screen by the viewers) ---------------------------
R.screen = ['Restricted design (gamma, two-facet): with few occasions the variance ' ...
    'components and pooled residual shown as posterior means are not recovered ' ...
    'at level in the in-repo checks; see the export header. G and D within about ' ...
    '.03/.06 for analyses 14 and 25; G and D biased low for analysis 27.'];

%---- saved-figure note (surface and subject plots; psyrat_restricted_figure_note)
%A saved PNG, PDF or FIG travels without the export header, so this note points
%at the manual instead.
R.figure = ['Restricted design (gamma two-facet: analyses 14, 25 and 27): variance ' ...
    'components and pooled residual reported as posterior means are not recovered ' ...
    'at level with few occasions; G and D within about .03/.06 for analyses 14 and ' ...
    '25; G and D biased low for analysis 27. See manual Section 5.4.'];

%---- Analyze-time notice (acknowledged before the save dialog) ----------------
if isempty(which)
    first = ['The selected settings fit a gamma two-facet design (analysis 14, ' ...
        '25 or 27).'];
else
    first = ['The selected settings fit ' names.(which) ' under the gamma family.'];
end
R.dialogtitle = 'Restricted design: gamma two-facet';
R.dialog = {first; ...
    ['This design is RESTRICTED: with few occasions, the variance components and ' ...
    'the pooled residual reported as posterior means are not recovered at level ' ...
    'in the in-repo checks (posterior mean over planted value, on the variance ' ...
    'scale: in five re-evaluated scope captures, occasion component 2.4 to 12.1 ' ...
    'times, remainder 2.2 to 7.2 times, pooled residual 1.1 to 1.9 times; in two ' ...
    'prospective analysis-27 fits, occasion component 60.1 and 2.45 times, ' ...
    'pooled residual 1.83 and 2.95 times, and the second fit did not converge).']; ...
    ['G and D stay within about .03 and .06 for analyses 14 and 25; for analysis ' ...
    '27 both G and D are biased low, while participant ordering is kept.']; ...
    ['The same disclosure appears in the results viewer and in every exported ' ...
    'table header. Manual Section 5.4 explains it.']; ...
    ['How to proceed: continue if the coefficients are what you need, or Cancel ' ...
    'at the next (save) dialog to change the settings.']};

end
