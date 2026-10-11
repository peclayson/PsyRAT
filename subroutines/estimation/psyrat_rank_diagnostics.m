function d = psyrat_rank_diagnostics(x)
%Rank/fold split Rhat, bulk/tail ESS and mean MCSE (posterior 1.6 conventions).
%Independent numerical references and finite-sample conventions:
%documentation/diagnostic_output_contract.md and tests/fixtures/rank_diagnostics/.
%Odd chains omit only their middle draw for splitting, never from saved draws.
%Tail ESS diagnoses 5/95% quantiles, not arbitrary interval-endpoint precision.
% Copyright (C) 2016-2026 Peter E. Clayson
% SPDX-License-Identifier: GPL-3.0-or-later

d=struct('rhat',Inf,'ess_bulk',0,'ess_tail',0,'mcse_mean',NaN, ...
    'rhat_rank',Inf,'rhat_fold',Inf,'status','invalid_or_frozen');
[n,m]=size(x); h=floor(n/2);
if ~isnumeric(x) || ~isreal(x) || h<3 || isempty(x) || any(~isfinite(x(:)))
    return;
end
%Original stochastic chains that never move remain a hard failure. Constant
%tail-indicator chains are different: ESS must use between-chain variation.
if any(max(x,[],1)==min(x,[],1)); return; end
split=@(v)[v(1:h,:),v(n-h+1:n,:)];
z=local_rank(split(x));
zfold=local_rank(split(abs(x-median(x(:)))));
d.rhat_rank=local_rhat(z); d.rhat_fold=local_rhat(zfold);
d.rhat=max(d.rhat_rank,d.rhat_fold);
d.ess_bulk=local_ess(z);
q=local_quantile7(x(:),[.05 .95]);
tails=[local_ess(split(double(x<=q(1)))),local_ess(split(double(x<=q(2))))];
if any(isnan(tails)); d.ess_tail=NaN; else; d.ess_tail=min(tails); end
d.mcse_mean=std(x(:),0)/sqrt(local_ess(split(x)));
d.status='assessed';
if ~isfinite(d.ess_tail); d.status='tail_unassessed'; end
end

function z=local_rank(x)
r=tiedrank(x(:)); z=reshape(sqrt(2)*erfinv(2*(r-3/8)/(numel(x)+1/4)-1),size(x));
end
function r=local_rhat(x)
n=size(x,1); w=mean(var(x,0,1));
if w<=0; r=Inf; else; r=sqrt(((n-1)/n*w+var(mean(x,1),0))/w); end
end
function e=local_ess(x)
[n,m]=size(x);
if n<3 || any(~isfinite(x(:))) || max(x(:))==min(x(:)); e=NaN; return; end
ac=zeros(n,m);
for c=1:m
    z=x(:,c)-mean(x(:,c)); f=fft(z,2^nextpow2(2*n));
    a=real(ifft(f.*conj(f))); ac(:,c)=a(1:n)/n;
end
acm=mean(ac,2); w=acm(1)*n/(n-1);
vp=w*(n-1)/n+var(mean(x,1),0);
if vp<=0; e=NaN; return; end
%Geyer initial-positive sequence followed by initial-monotone pairs. Indices
%below are lags (zero based), converted to MATLAB subscripts only at access.
rho=zeros(n,1); t=0; even=1; rho(1)=1;
odd=1-(w-acm(2))/vp; rho(2)=odd;
while t<n-5 && ~isnan(even+odd) && even+odd>0
    t=t+2; even=1-(w-acm(t+1))/vp; odd=1-(w-acm(t+2))/vp;
    if even+odd>=0; rho(t+1)=even; rho(t+2)=odd; end
end
last=t;
if even>0; rho(last+1)=even; end
t=0;
while t<=last-4
    t=t+2;
    if rho(t+1)+rho(t+2)>rho(t-1)+rho(t)
        rho(t+1)=(rho(t-1)+rho(t))/2; rho(t+2)=rho(t+1);
    end
end
%posterior's R index 1:0 retains element 1 when no lag-pair iteration runs.
%Keep its documented finite short-chain convention; this is not adequacy.
if last==0; paired=rho(1); else; paired=sum(rho(1:last)); end
tau=-1+2*paired+rho(last+1);
e=n*m/max(tau,1/log10(n*m));
end
function q=local_quantile7(x,p)
x=sort(x); h=1+(numel(x)-1)*p; lo=floor(h); hi=ceil(h);
q=(1-(h-lo)).*x(lo)' +(h-lo).*x(hi)';
end
