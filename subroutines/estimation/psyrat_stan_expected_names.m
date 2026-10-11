function [names,status] = psyrat_stan_expected_names(files,keepfun)
%Read the independent monitored-output inventory from every chain CSV header.
% Copyright (C) 2016-2026 Peter E. Clayson
% SPDX-License-Identifier: GPL-3.0-or-later
names={}; status='unavailable';
if ischar(files)||isstring(files); files=cellstr(files); end
if ~iscell(files)||isempty(files); return; end
for k=1:numel(files)
    fid=fopen(files{k},'r'); if fid<0; status='unavailable'; return; end
    closer=onCleanup(@()fclose(fid));
    hdr=fgetl(fid);
    while ischar(hdr)&&(isempty(strtrim(hdr))||startsWith(strtrim(hdr),'#')); hdr=fgetl(fid); end
    if ~ischar(hdr); status='unavailable'; return; end
    labels=strsplit(strtrim(hdr),',');
    for j=1:numel(labels)
        label=strrep(strtrim(labels{j}),'"','');
        tok=regexp(label,'^([A-Za-z_]\w*)((?:\.\d+)+)$','tokens','once');
        if ~isempty(tok); label=[tok{1} '[' strrep(tok{2}(2:end),'.',',') ']']; end
        labels{j}=label;
    end
    labels=labels(cellfun(keepfun,labels));
    if isempty(labels); status='unavailable'; return; end
    if k>1 && ~isequal(labels,names); status='chain_header_mismatch'; names=unique([names labels],'stable'); return; end
    names=labels;
    clear closer;
end
status='chain_csv_headers';
end
