function contract = psyrat_constant_contract(source,data,names)
%Conservative structural constants from the actual emitted Stan model/data.
%This is deliberately not a general Stan interpreter. Unsupported syntax is
%unassessed; an observed NaN or a correlation-like name grants no exemption.
% Copyright (C) 2016-2026 Peter E. Clayson
% SPDX-License-Identifier: GPL-3.0-or-later
contract=struct('status','unavailable','names',{{}},'values',[],'reasons',{{}});
if iscell(source); source=strjoin(source,newline); end
if isstring(source); source=char(join(source,newline)); end
if ~ischar(source) || isempty(source); return; end
[source,valid]=local_code(source);
if ~valid; return; end
[blocks,valid]=local_blocks(source);
if ~valid; return; end
top=''; topStatements={}; topBlocks={};
for block={'parameters','transformed_parameters','generated_quantities'}
    statements=blocks.(block{1});
    for k=1:numel(statements)
        if ~contains(statements{k},'{')
            top=[top statements{k} newline]; topStatements{end+1}=statements{k}; %#ok<AGROW>
            topBlocks{end+1}=block{1}; %#ok<AGROW>
        end
    end
end
if isempty(top); return; end
contract.status='available_bounded_rules';
kind={}; base={}; dim=[]; ready=[];
for k=1:numel(topStatements)
    t=regexp(topStatements{k},'^(corr_matrix|cholesky_factor_corr)\s*\[([^\]]+)\]\s*([A-Za-z_]\w*)\s*[;=]','tokens','once');
    if isempty(t); continue; end
    n=local_dimension(t{2},data);
    if isnan(n); continue; end
    kind{end+1}=t{1}; base{end+1}=t{3}; dim(end+1)=n; ready(end+1)=Inf; %#ok<AGROW>
    writes=local_writes(source,t{3});
    if strcmp(topBlocks{k},'parameters') && writes==0
        ready(end)=0;
    elseif writes==1
        % A final TP/GQ constraint alone cannot prove an earlier copy. With
        % one complete unconditional assignment, the copied and final parent
        % are identical only after that assignment has executed.
        assignment=['^(?:(?:corr_matrix|cholesky_factor_corr)\s*\[[^\]]+\]\s*)?' t{3} '\s*=(?!=)'];
        for j=k:numel(topStatements)
            if ~isempty(regexp(topStatements{j},assignment,'once')); ready(end)=j; break; end
        end
    end
end
%Current plain-matrix correlation aliases are square principal blocks only.
a={}; aIndex=[];
for k=1:numel(topStatements)
    match=regexp(topStatements{k},['^(?:matrix\s*\[[^\]]+\]\s*)?' ...
        '([A-Za-z_]\w*)\s*=\s*block\(\s*([A-Za-z_]\w*)\s*,\s*1\s*,\s*1\s*,\s*([^,]+),\s*([^\)]+)\)\s*;$'],'tokens','once');
    if ~isempty(match); a{end+1}=match; aIndex(end+1)=k; end %#ok<AGROW>
end
for k=1:numel(a)
    p=find(strcmp(base,a{k}{2}) & strcmp(kind,'corr_matrix'),1);
    n=local_dimension(strtrim(a{k}{3}),data); n2=local_dimension(strtrim(a{k}{4}),data);
    decl=['\<matrix\s*\[\s*' regexptranslate('escape',strtrim(a{k}{3})) '\s*,\s*' regexptranslate('escape',strtrim(a{k}{4})) '\s*\]\s*' a{k}{1} '\s*[;=]'];
    if ~isempty(p) && ready(p)<aIndex(k) && n==n2 && n<=dim(p) && ~isempty(regexp(top,decl,'once')) && local_writes(source,a{k}{1})==1
        kind{end+1}='corr_matrix'; base{end+1}=a{k}{1}; dim(end+1)=n; ready(end+1)=aIndex(k); %#ok<AGROW>
    end
end
%Case8's reviewed zero graph: a disabled data flag fixes transformed rho and
%its residual covariance products. Prior-only rescor_mu/sd_rescor remain free.
d=regexprep(blocks.data,'\s+','');
tp=regexprep(blocks.transformed_parameters,'\s+','');
gq=regexprep(blocks.generated_quantities,'\s+','');
rhoLoop=['for(sin1:NSUB){if(estimate_rescor==1){' ...
    'rho[s]=tanh(rescor_mu+sd_rescor*z_rescor[s]);}else{rho[s]=0;}}'];
covLoop=['for(sin1:NSUB){' ...
    'er_var_ss[s,1]=b_sigma[1]+r_sigma[s,1];' ...
    'er_var_ss[s,2]=b_sigma[2]+r_sigma[s,2];' ...
    'wp_cov_ss[s]=rho[s]*exp(er_var_ss[s,1])*exp(er_var_ss[s,2]);' ...
    'err_varcov[1,1]+=exp(2*er_var_ss[s,1])/NSUB;' ...
    'err_varcov[2,2]+=exp(2*er_var_ss[s,2])/NSUB;' ...
    'err_varcov[1,2]+=wp_cov_ss[s]/NSUB;' ...
    'err_varcov[2,1]+=wp_cov_ss[s]/NSUB;}'];
% Match complete top-level statements in their actual Stan blocks. Snippets
% in strings, a decoy identifier/function or another control flow are no proof.
zerosAllowed=isstruct(data) && isfield(data,'estimate_rescor') && ...
    isequal(double(data.estimate_rescor),0) && isfield(data,'NSUB') && ...
    ~isnan(local_dimension('NSUB',data)) && ...
    all(ismember({'int<lower=1>NSUB;','int<lower=0,upper=1>estimate_rescor;'},d)) && ...
    all(ismember({'vector[NSUB]rho;',rhoLoop},tp)) && ...
    all(ismember({'matrix[2,2]err_varcov;','matrix[NSUB,2]er_var_ss;', ...
    'vector[NSUB]wp_cov_ss;','vector[1]rescor;', ...
    'rescor[1]=mean(rho);','err_varcov=rep_matrix(0,2,2);',covLoop},gq)) && ...
    local_writes(source,'rho')==2 && local_writes(source,'wp_cov_ss')==1 && ...
    local_writes(source,'rescor')==1 && local_writes(source,'err_varcov')==5;
for k=1:numel(names)
    tok=regexp(names{k},'^([A-Za-z_]\w*)\[(\d+),(\d+)\]$','tokens','once');
    value=NaN; why='';
    if ~isempty(tok)
        p=find(strcmp(base,tok{1}),1); i=str2double(tok{2}); j=str2double(tok{3});
        if ~isempty(p) && all([i j]>=1) && all([i j]<=dim(p))
            if strcmp(kind{p},'corr_matrix') && i==j
                value=1; why='unit diagonal of declared correlation matrix/principal block';
            elseif strcmp(kind{p},'cholesky_factor_corr') && (i<j || (i==1 && j==1))
                value=double(i==j); why='structural entry of declared Cholesky correlation factor';
            end
        end
    end
    if zerosAllowed
        one=regexp(names{k},'^(rho|wp_cov_ss|rescor)\[(\d+)\]$','tokens','once');
        if ~isempty(one)
            n=double(data.NSUB); if strcmp(one{1},'rescor'); n=1; end
            idx=str2double(one{2});
            if isscalar(n) && idx>=1 && idx<=n; value=0; why='case8 estimate_rescor=0 transformed residual coupling'; end
        elseif any(strcmp(names{k},{'err_varcov[1,2]','err_varcov[2,1]'}))
            value=0; why='case8 estimate_rescor=0 residual covariance';
        end
    end
    if ~isnan(value)
        contract.names{end+1}=names{k}; contract.values(end+1)=value; contract.reasons{end+1}=why;
    end
end
end
function n=local_dimension(expr,data)
n=str2double(strtrim(expr));
if isnan(n)
    parts=strsplit(regexprep(expr,'\s+',''),'*'); mult=1; name='';
    if numel(parts)==1; name=parts{1};
    elseif numel(parts)==2; mult=str2double(parts{1}); name=parts{2}; end
    if ~isempty(regexp(name,'^[A-Za-z_]\w*$','once')) && ...
            isfinite(mult) && mult>=1 && mult==fix(mult) && ...
            isstruct(data) && isfield(data,name) && isnumeric(data.(name)) && isscalar(data.(name))
        n=mult*double(data.(name));
    end
end
if ~isscalar(n) || ~isreal(n) || ~isfinite(n) || n<1 || n~=fix(n); n=NaN; end
end
function [code,valid]=local_code(source)
% Mask comments and double-quoted strings without joining adjacent tokens.
% A quote inside a comment and a comment delimiter inside a string are inert.
code=source; state=0; k=1;
while k<=numel(source)
    ch=source(k); next=' '; if k<numel(source); next=source(k+1); end
    if state==0
        if ch=='"'; state=1; code(k)=' ';
        elseif ch=='/' && next=='/'; state=2; code(k:k+1)='  '; k=k+1;
        elseif ch=='/' && next=='*'; state=3; code(k:k+1)='  '; k=k+1;
        end
    elseif state==1
        % Stan treats backslashes literally; quotes always close strings.
        if ch==newline || ch==char(13); valid=false; return; end
        code(k)=' ';
        if ch=='"'; state=0; end
    elseif state==2
        if ch==newline || ch==char(13); state=0; else; code(k)=' '; end
    else
        if ch=='*' && next=='/'; code(k:k+1)='  '; k=k+1; state=0;
        elseif ch~=newline && ch~=char(13); code(k)=' '; end
    end
    k=k+1;
end
valid=state==0 || state==2; % An unfinished string/block comment grants no proof.
end

function [blocks,valid]=local_blocks(source)
blocks=struct('data',{{}},'parameters',{{}},'transformed_parameters',{{}}, ...
    'generated_quantities',{{}}); valid=false; seen={}; pos=1;
while pos<=numel(source)
    if isempty(strtrim(source(pos:end))); valid=true; return; end
    opening=find(source(pos:end)=='{',1);
    if isempty(opening); return; end
    opening=opening+pos-1; header=regexprep(strtrim(source(pos:opening-1)),'\s+','_');
    if ~ismember(header,{'functions','data','transformed_data','parameters', ...
            'transformed_parameters','model','generated_quantities'}) || ismember(header,seen); return; end
    seen{end+1}=header; depth=1; closing=opening+1; %#ok<AGROW>
    while closing<=numel(source) && depth>0
        if source(closing)=='{'; depth=depth+1;
        elseif source(closing)=='}'; depth=depth-1; end
        closing=closing+1;
    end
    if depth~=0; return; end
    if isfield(blocks,header)
        body=source(opening+1:closing-2); statements={}; start=1; depth=0;
        for k=1:numel(body)
            if body(k)=='{'; depth=depth+1;
            elseif body(k)=='}'; depth=depth-1; end
            if depth==0 && (body(k)==';' || body(k)=='}')
                statements{end+1}=strtrim(body(start:k)); start=k+1; %#ok<AGROW>
            end
        end
        if ~isempty(strtrim(body(start:end))); return; end
        blocks.(header)=statements;
    end
    pos=closing;
end
valid=true;
end

function n=local_writes(source,name)
%Reject a recognized alias/zero graph with any extra assignment. This bounded
%contract intentionally refuses unsupported syntax rather than interpreting it.
ends=regexp(source,['\<' name '\>'],'end'); n=0;
for k=1:numel(ends)
    pos=ends(k)+1;
    while pos<=numel(source) && isspace(source(pos)); pos=pos+1; end
    % Stan permits nested index expressions and chained indexing. Consume
    % balanced brackets before testing the assignment operator; a flat regex
    % can miss a stochastic write such as alias[index[1],1] = z.
    while pos<=numel(source) && source(pos)=='['
        depth=1; pos=pos+1;
        while pos<=numel(source) && depth>0
            if source(pos)=='['; depth=depth+1;
            elseif source(pos)==']'; depth=depth-1; end
            pos=pos+1;
        end
        if depth>0; n=Inf; return; end % unsupported/malformed source grants no proof
        while pos<=numel(source) && isspace(source(pos)); pos=pos+1; end
    end
    if ~isempty(regexp(source(pos:end),'^(?:(?:\.[*/]|[+*/-])?=|<-)(?!=)','once'))
        n=n+1;
    end
end
end
