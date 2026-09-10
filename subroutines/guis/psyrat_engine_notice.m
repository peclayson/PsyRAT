function msg = psyrat_engine_notice(rel)
%Human-readable notice flagging that a result was estimated with a native
%(non-CmdStan) engine, for surfacing in the results viewer.
%
%msg = psyrat_engine_notice(rel)
%
%Input:
% rel - the reliability results struct (psyrat_data.rel). Only rel.engine is
%  read: the provenance label stamped by the estimation core
%  (psyrat_fit_engine): 'cmdstan' | 'hmc' | 'fitlme'.
%
%Output:
% msg - a notice string to show the user (e.g. via warndlg or a viewer label),
%  or '' when no notice is needed: a CmdStan run, or a legacy results file with
%  no engine field (which predates the native engines and is therefore CmdStan).
%
%The CmdStan path is the default and produces no notice, so existing CmdStan
%workflows are visually unchanged. Native runs are flagged as not bit-identical
%to CmdStan (HMC) or as approximate frequentist intervals (fitlme).

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

msg = '';
if ~isstruct(rel) || ~isfield(rel,'engine') || isempty(rel.engine)
    return;   % legacy file or unrecorded engine -> treat as CmdStan, no notice
end

switch lower(char(rel.engine))
    case 'cmdstan'
        msg = '';
    case 'hmc'
        msg = ['These results were estimated with the native HMC engine ' ...
            '(hmcSampler), not CmdStan. The draws are genuine Bayesian posterior ' ...
            'samples, but they are NOT bit-identical to a CmdStan run. Re-run ' ...
            'with CmdStan for a definitive analysis.'];
    case 'fitlme'
        msg = ['These results were estimated with the native fitlme engine ' ...
            '(REML + parametric bootstrap), not CmdStan. Its intervals are ' ...
            'APPROXIMATE (frequentist), NOT Bayesian credible intervals, and may ' ...
            'differ from a CmdStan run. Re-run with CmdStan for a definitive ' ...
            'analysis.'];
    otherwise
        msg = sprintf(['These results were estimated with the native ''%s'' ' ...
            'engine, not CmdStan; treat them as approximate and re-run with ' ...
            'CmdStan for a definitive analysis.'], char(rel.engine));
end

end
