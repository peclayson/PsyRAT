function [rows,info] = psyrat_parse_stansummary_csv(csvpath, keepfun)
%Read the CSV that "stansummary --csv_filename" writes and return, for the
%parameters the caller monitors, the convergence diagnostics located by column
%NAME at the precision the file carries.
%
% rows = psyrat_parse_stansummary_csv(csvpath, keepfun)
%
%Inputs
% csvpath - the CSV written by stansummary: one header row, double-quoted
%   parameter names (a matrix element such as "cor_id[1,2]" carries a comma
%   inside its quotes), literal nan for a constant parameter, and trailing
%   '#' comment lines.
% keepfun - function handle; keepfun(name) is true for a parameter row to keep.
%   psyrat_storeconv passes the same prefix test its text parser applies, so
%   the monitored set is unchanged by the route.
%
%Output
% rows - N x 6 cell array {name, n_eff, r_hat, ess_bulk, ess_tail, mcse} in
%   file order, one row per kept parameter. n_eff is min(ess_bulk, ess_tail),
%   the value psyrat_checkconv gates on (owner ruling 2026-09-10, work queue
%   Q4); a row whose ESS is nan on either side keeps nan, so a constant
%   parameter (nan in both diagnostics) is still recognized as such downstream.
%
%Layouts. CmdStan 2.33 and later write
%   name,Mean,MCSE,StdDev,MAD,5%,50%,95%,ESS_bulk,ESS_tail,ESS_bulk/s,R_hat
%and earlier versions
%   name,Mean,MCSE,StdDev,5%,50%,95%,N_Eff,N_Eff/s,R_hat
%where N_Eff supplies only n_eff; unmeasured bulk/tail ESS remain NaN. Every column is found by its header name:
%the positional text parser this route replaces read the third-from-last token
%of the printed table as n_eff, which is N_Eff on one layout and ESS_tail on
%the other, at the two significant figures print(fit) requests (post-beta
%review PBA-01, 2026-09-10).
%
%Errors
% PsyRAT:stansummaryCsvMissing - the file does not exist.
% PsyRAT:stansummaryCsvHeader  - the header names no name/R_hat column, or
%   neither ESS_bulk/ESS_tail nor N_Eff.

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

if ~(ischar(csvpath) || (isstring(csvpath) && isscalar(csvpath))) || ...
        exist(char(csvpath),'file') ~= 2
    error('PsyRAT:stansummaryCsvMissing', ...
        'The stansummary CSV ''%s'' does not exist.', char(csvpath));
end

txt = fileread(char(csvpath));
lines = regexp(txt,'\r\n|\n|\r','split');
lines = strtrim(lines);
lines = lines(~cellfun('isempty',lines));
lines = lines(~startsWith(lines,'#'));   % the trailing legend lines
if isempty(lines)
    error('PsyRAT:stansummaryCsvHeader', ...
        'The stansummary CSV ''%s'' has no header row.', char(csvpath));
end

%fields are either a double-quoted string (which may contain commas) or a run
%of characters up to the next comma; stansummary never writes an empty field
hdr = local_fields(lines{1});
col = @(name) find(strcmp(hdr,name),1);
iname = col('name'); irhat = col('R_hat'); imcse = col('MCSE');
ibulk = col('ESS_bulk'); itail = col('ESS_tail'); ineff = col('N_Eff');
if isempty(iname) || isempty(irhat)
    error('PsyRAT:stansummaryCsvHeader', ...
        ['The stansummary CSV ''%s'' has no name or R_hat column ' ...
        '(header: %s).'], char(csvpath), strjoin(hdr,','));
end
legacy = isempty(ibulk) || isempty(itail);
info=psyrat_diagnostic_metadata({});
info.format='modern_csv'; info.degraded=false; info.degradation_reason='';
info.precision_significant_digits=NaN; % source tokens; caller records requested precision
info.row_columns={'name','n_eff','r_hat','ess_bulk','ess_tail','mcse_mean'};
if legacy
    if isempty(ineff)
        error('PsyRAT:stansummaryCsvHeader', ...
            ['The stansummary CSV ''%s'' has neither ESS_bulk/ESS_tail nor ' ...
            'N_Eff columns (header: %s).'], char(csvpath), strjoin(hdr,','));
    end
    info.format='legacy_csv'; info.degraded=true;
    info.degradation_reason='Legacy N_Eff only; bulk/tail ESS are not measured.';
    ibulk=ineff; itail=ineff;
end

rows = cell(0,6);
for i = 2:numel(lines)
    f = local_fields(lines{i});
    if numel(f)<iname; continue; end
    name = f{iname};
    if ~keepfun(name)
        continue;
    end
    if numel(f)~=numel(hdr)
        error('PsyRAT:stansummaryCsvRow','Malformed monitored row: %s.',name);
    end
    numericcols=unique([ibulk itail irhat imcse]);
    for j=numericcols
        value=str2double(f{j});
        if ~isreal(value) || (isnan(value) && ~any(strcmpi(f{j},{'nan','-nan','+nan'})))
            error('PsyRAT:stansummaryCsvRow','Invalid numeric token for %s: %s.',name,f{j});
        end
    end
    if any(strcmp(rows(:,1),name))
        error('PsyRAT:stansummaryCsvRow','Duplicate monitored row: %s.',name);
    end
    bulk = str2double(f{ibulk});
    tail = str2double(f{itail});
    rhat = str2double(f{irhat});
    mcse = NaN;
    if ~isempty(imcse)
        mcse = str2double(f{imcse});
    end
    %MATLAB's two-input min ignores NaN, so spell the constant case out: nan on
    %either side stays nan (a constant parameter prints nan on both)
    if isnan(bulk) || isnan(tail)
        neff = NaN;
    else
        neff = min(bulk,tail);
    end
    if legacy; bulk=NaN; tail=NaN; end
    rows(end+1,:) = {name, neff, rhat, bulk, tail, mcse}; %#ok<AGROW>
end
end

function f = local_fields(line)
%Split one CSV line into fields: a double-quoted string (quotes removed, commas
%inside kept) or a run of non-comma characters.
f = regexp(line,'"[^"]*"|[^,]+','match');
f = strtrim(f);
for k = 1:numel(f)
    if numel(f{k}) >= 2 && f{k}(1) == '"' && f{k}(end) == '"'
        f{k} = f{k}(2:end-1);
    end
end
end
