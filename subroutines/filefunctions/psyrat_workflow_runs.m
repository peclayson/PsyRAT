function runs = psyrat_workflow_runs(relAnalysis)
%Build workflow sweep run configurations for runExhaustiveWorkflowBugSweep.

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

runs = {};

switch lower(relAnalysis)
    case 'ic'
        runs{end+1} = struct('analysis', 'sing', 'depcutoff', 0.5); %#ok<AGROW>
    case 'ic_sserrvar'
        runs{end+1} = struct('analysis', 'sing_sserr', 'depcutoff', 0.5); %#ok<AGROW>
    case 'ic_diff'
        runs{end+1} = struct('analysis', 'sing_diff', 'depcutoff', 0.5, ...
            'diffgcoeff', [1 2]); %#ok<AGROW>
    case 'ic_diff_sserrvar'
        %Only the subject-level builder. psyrat_startview_sing maps
        %rel.analysis to a builder label 1:1 (ic_diff_sservar ->
        %sing_diff_sserr, ic_diff -> sing_diff), and psyrat_report and
        %psyrat_run never produce any other pairing, so sweeping this
        %estimation identity through the 'sing_diff' builder exercised a
        %configuration no product surface can create. That pairing produced a
        %relsummary with no diffscore.ssrel_table while the output functions
        %still dispatched on rel.analysis and took the subject-level branch,
        %which is the mislabeled-viewer failure state that
        %TestSummaryAndPresentationFunctions deliberately pins as a LOUD crash.
        %The 'sing_diff' builder keeps its coverage under case 'ic_diff' above.
        runs{end+1} = struct('analysis', 'sing_diff_sserr', 'depcutoff', 0.5, ...
            'diffgcoeff', [1 2]); %#ok<AGROW>
    case 'trt'
        runs{end+1} = struct('analysis', 'trt', 'relcutoff', 0.5); %#ok<AGROW>
    case 'trt_sserrvar'
        runs{end+1} = struct('analysis', 'trt_sserr', 'depcutoff', 0.5); %#ok<AGROW>
    case 'trt_diff'
        runs{end+1} = struct('analysis', 'trt_diff', 'relcutoff', 0.5, ...
            'diffgcoeff', [1 2]); %#ok<AGROW>
    otherwise
        runs{end+1} = struct('analysis', 'unknown', 'depcutoff', 0.5); %#ok<AGROW>
end

end
