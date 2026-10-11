function [rows,info] = psyrat_parse_stansummary_text(output,keepfun)
%Degraded fallback: retain format/precision and use named ESS columns if known.
% Copyright (C) 2016-2026 Peter E. Clayson
% SPDX-License-Identifier: GPL-3.0-or-later
if ischar(output); output=regexp(output,'\r\n|\n','split'); end
info=psyrat_diagnostic_metadata({});
info.format='unknown_positional_text'; info.precision_significant_digits=2;
info.degradation_reason='Printed diagnostic summary at two significant figures.';
info.row_columns={'name','n_eff','r_hat','ess_bulk','ess_tail','mcse_mean'};
hdr={};
for k=1:numel(output)
    f=regexp(strtrim(output{k}),'\S+','match');
    if any(strcmp(f,'R_hat')) && any(strcmp(f,'Mean'))
        hdr=f; break;
    end
end
bulk=find(strcmp(hdr,'ESS_bulk'),1); tail=find(strcmp(hdr,'ESS_tail'),1);
legacy=find(strcmp(hdr,'N_Eff'),1); rh=find(strcmp(hdr,'R_hat'),1);
mc=find(strcmp(hdr,'MCSE'),1);
if ~isempty(bulk)&&~isempty(tail); info.format='modern_text';
elseif ~isempty(legacy); info.format='legacy_text'; end
rows=cell(0,6);
for k=1:numel(output)
    f=regexp(strtrim(output{k}),'\S+','match');
    if isempty(f)||~keepfun(f{1}); continue; end
    b=NaN; t=NaN; m=NaN; e=0; r=Inf;
    malformed=any(strcmp(rows(:,1),f{1})) || ...
        ~all(cellfun(@local_numeric_token,f(2:end)));
    if ~isempty(hdr) && numel(f)==numel(hdr)+1 && ~isempty(rh)
        r=str2double(f{rh+1});
        if ~isempty(mc); m=str2double(f{mc+1}); end
        if strcmp(info.format,'modern_text')
            b=str2double(f{bulk+1}); t=str2double(f{tail+1});
            e=min(b,t); if isnan(b)||isnan(t); e=NaN; end
        elseif strcmp(info.format,'legacy_text'); e=str2double(f{legacy+1}); end
    elseif isempty(hdr) && numel(f)>=3
        %Unknown older text: semantics cannot be recovered from position alone.
        e=str2double(f{end-2}); r=str2double(f{end});
    else
        malformed=true;
    end
    if malformed
        % Invalid strings must not become a source-proven constant's NaN/NaN
        % signature. Keep the named row, with an unmistakable failing result.
        e=0; r=Inf; b=NaN; t=NaN; m=NaN;
        info.degradation_reason=['Printed diagnostic summary at two significant figures. ' ...
            'Malformed numeric fields, row width or duplicate monitored rows.'];
    end
    rows(end+1,:)={f{1},e,r,b,t,m}; %#ok<AGROW>
end
end

function valid=local_numeric_token(token)
value=str2double(token);
valid=isreal(value) && (~isnan(value) || any(strcmpi(token,{'nan','-nan','+nan'})));
end
