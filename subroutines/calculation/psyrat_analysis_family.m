function fam = psyrat_analysis_family(analysis)
%Classify a PsyRAT analysis string into its reliability "family".
%
% fam = psyrat_analysis_family(analysis)
%
%PsyRAT records the design/estimand of a run in REL.analysis (e.g. 'ic',
%'trt_diff', 'ic_dynrel_trt', 'trt_splits'). Several sites need to route on the
%broad family that string belongs to rather than the exact variant: the GUI
%results router (which viewer to open), the dynamic-reliability and data-splits
%summary guards (which results they accept), and the headless report dispatch
%(which table builder to call). Each of those used to carry its own inline
%enumeration of the variant strings, so adding a new variant meant editing every
%site or risking the GUI and the headless report disagreeing. This function is
%the single source of truth for that mapping.
%
%Input
% analysis - the REL.analysis string (char or string scalar). Missing or empty
%            is treated as 'ic', the implicit one-facet default carried by
%            pre-analysis-field result structs (so it maps to 'sing'), matching
%            the GUI router's historical "no analysis field -> single-session
%            viewer" behavior.
%
%Output
% fam - one of:
%   'sing'   - one-facet (single session): 'ic', 'ic_sserrvar', 'ic_diff',
%              'ic_diff_sserrvar' (and the empty/missing default)
%   'trt'    - test-retest: 'trt', 'trt_sserrvar', 'trt_diff'
%   'dod'    - difference-of-differences: 'ic_dodiff'
%   'splits' - nonparallel data splits (observed design): 'ic_splits',
%              'trt_splits'
%   'dynrel' - dynamic/conditional reliability (Rast & Clayson) and its
%              difference/two-facet/rescor/rho variants, plus the
%              person-specific dynamic nonconcurrent difference-of-differences
%              (the 16-string set)
%   ''       - an unrecognized non-empty string (callers keep their own
%              fallback: the router opens no viewer, the summary guards error,
%              the report's slipped-variant guard reports it).
%
%Example
% fam = psyrat_analysis_family(REL.analysis);
% switch fam
%     case 'dynrel'; ...
% end

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

%missing/empty -> the legacy one-facet default ('ic'); accept char or string
if nargin < 1 || isempty(analysis)
    analysis = 'ic';
end
analysis = char(analysis);

switch analysis
    case {'ic','ic_sserrvar','ic_diff','ic_diff_sserrvar'}
        fam = 'sing';
    case {'trt','trt_sserrvar','trt_diff'}
        fam = 'trt';
    case 'ic_dodiff'
        fam = 'dod';
    case {'ic_splits','trt_splits'}
        fam = 'splits';
    case {'ic_dynrel','ic_dynrel_sserrvar',...
            'ic_diff_dynrel','ic_diff_dynrel_sserrvar',...
            'ic_dynrel_trt','ic_dynrel_sserrvar_trt','ic_diff_dynrel_trt',...
            'ic_diff_dynrel_sserrvar_trt',...
            'ic_diff_dynrel_rescor','ic_diff_dynrel_trt_rescor',...
            'ic_diff_dynrel_sserrvar_rescor',...
            'ic_diff_dynrel_sserrvar_trt_rescor',...
            'ic_diff_dynrel_sserrvar_rho','ic_diff_dynrel_sserrvar_trt_rho',...
            'ic_dodiff_dynrel_sserrvar','ic_dodiff_dynrel_sserrvar_trt'}
        %the dynamic DoD strings route as dynrel (they surface through
        %psyrat_dynrel_summary / psyrat_startview_dynrel), NOT as 'dod'
        %(the static viewer's machinery assumes the case-9 conventions)
        fam = 'dynrel';
    otherwise
        fam = '';
end

end
