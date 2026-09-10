function labels = psyrat_unitlabels(splits)
%Return the analysis-unit label fragments for plots, tables, and viewers.
%
%labels = psyrat_unitlabels(splits)
%
%PsyRAT was built for single-trial data, where the measurement facet is the
%trial. The splits feature lets the unit of analysis be a data split (a chunk
%of trials averaged into one score). Output must make this distinction clear,
%so every label that used to say "trial" is sourced from this helper instead
%of being hard-coded, and switches to "split" when split data are analyzed.
%
%Input
% splits - analysis-unit code carried on the data/REL struct:
%          1 = single-trial data (default), 2 = parallel splits,
%          3 = nonparallel splits. A logical/empty value is also accepted
%          (true -> split, false/empty -> single-trial).
%
%Output
% labels - struct of label fragments:
%   .issplit     - logical, true for split data (splits == 2 or 3)
%   .Unit        - 'Trial'  | 'Split'
%   .Units       - 'Trials' | 'Splits'
%   .unit        - 'trial'  | 'split'
%   .units       - 'trials' | 'splits'
%   .NumberOf    - 'Number of Trials' | 'Number of Splits'
%   .numberof    - 'number of trials' | 'number of splits'
%   .Cutoff      - 'Trial Cutoff' | 'Split Cutoff'
%   .CutoffVar   - 'Trial_Cutoff' | 'Split_Cutoff' (table variable name)
%   .NumSymbol   - '# of Trials'  | '# of Splits'
%   .banner      - one-line analysis-unit banner for figures/tables
%
%Example
% L = psyrat_unitlabels(REL.splits);
% xlabel(L.NumberOf);

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

%default to single-trial when nothing was passed (backward compatible with
%every pre-splits caller/struct that has no splits field)
if nargin < 1 || isempty(splits)
    issplit = false;
elseif islogical(splits)
    issplit = any(splits);
else
    %numeric code: 1 = single-trial, 2/3 = split data
    issplit = any(splits > 1);
end

labels = struct();
labels.issplit = issplit;

if ~issplit
    labels.Unit      = 'Trial';
    labels.Units     = 'Trials';
    labels.unit      = 'trial';
    labels.units     = 'trials';
    labels.NumberOf  = 'Number of Trials';
    labels.numberof  = 'number of trials';
    labels.Cutoff    = 'Trial Cutoff';
    labels.CutoffVar = 'Trial_Cutoff';
    labels.NumSymbol = '# of Trials';
    labels.banner    = 'Analysis unit: single trials';
else
    labels.Unit      = 'Split';
    labels.Units     = 'Splits';
    labels.unit      = 'split';
    labels.units     = 'splits';
    labels.NumberOf  = 'Number of Splits';
    labels.numberof  = 'number of splits';
    labels.Cutoff    = 'Split Cutoff';
    labels.CutoffVar = 'Split_Cutoff';
    labels.NumSymbol = '# of Splits';
    labels.banner    = 'Analysis unit: SPLITS (each score is the mean of n_i items)';
end

end
