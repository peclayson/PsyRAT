function label = psyrat_interval_label(rel,style)
%Name the interval actually produced by the saved estimation engine.
%Legacy results without an engine retain their Bayesian label.
% Copyright (C) 2016-2026 Peter E. Clayson
% SPDX-License-Identifier: GPL-3.0-or-later
if nargin<2; style='sentence'; end
label='Credible interval';
if isstruct(rel) && isfield(rel,'engine') && strcmpi(rel.engine,'fitlme')
    label='Bootstrap interval';
end
if strcmp(style,'lower'); label=lower(label);
elseif strcmp(style,'title'); label=strrep(label,'interval','Interval'); end
end
