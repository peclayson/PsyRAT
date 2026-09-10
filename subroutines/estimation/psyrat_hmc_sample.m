function [draws,conv] = psyrat_hmc_sample(logpdfgrad,start,nchains,nwarmup, ...
    nsampling,seed,monitor_idx,monitor_names)
%Shared Hamiltonian Monte Carlo sampling harness for the native HMC engine.
%
%[draws,conv] = psyrat_hmc_sample(logpdfgrad,start,nchains,nwarmup,nsampling,...
%    seed,monitor_idx,monitor_names)
%
%This is the reusable Phase-2 core: every native-HMC model supplies a
%log-posterior (with gradient) over its unconstrained parameter vector, and this
%function runs the chains, derives per-chain seeds from the master seed, and
%returns the pooled draws plus a convergence cell in the same layout
%psyrat_storeconv / psyrat_checkconv consume. The per-model code is responsible
%for (a) building logpdfgrad, (b) choosing a start point, and (c) transforming
%the returned draws into the constrained variance components REL.out expects.
%
%Inputs:
% logpdfgrad - function handle theta -> [logp, grad] over the D-dim unconstrained
%  parameter vector (column). grad is the D-dim gradient (column). The model
%  must include any Jacobian terms for constrained parameters in logp/grad.
% start - D x 1 initial point (shared across chains; per-chain RNG jitters them).
% nchains, nwarmup, nsampling - chain count, warmup (burn-in) and post-warmup
%  draws per chain.
% seed - master RNG seed; chain c uses seed+c-1 so runs are reproducible.
% monitor_idx - indices (into theta) of the parameters to report convergence for
%  (typically the variance components, not the thousands of random effects).
% monitor_names - cellstr of names for those indices (for the conv cell labels).
%
%Outputs:
% draws - (nchains*nsampling) x D matrix of pooled post-warmup draws.
% conv - {'name','n_eff','r_hat'; name_i, neff_i, rhat_i; ...} for the monitored
%  parameters, computed with split-R-hat and an autocovariance effective sample
%  size (Stan conventions). psyrat_checkconv flags r_hat>=1.1 or n_eff<10*nchains.

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

D = numel(start);
percnt = zeros(nsampling,D,nchains);   % per-chain draws for R-hat
for c = 1:nchains
    rng(double(seed)+c-1,'twister');
    %a small random jitter of the start point decorrelates the chains so
    %split-R-hat is meaningful; the master+offset seed keeps it reproducible
    st = start(:) + 0.01*randn(D,1);
    smp = hmcSampler(logpdfgrad,st,'CheckGradient',false);
    smp = tuneSampler(smp,'VerbosityLevel',0, ...
        'NumStepSizeTuningIterations',max(50,round(nwarmup/2)));
    chain = drawSamples(smp,'NumSamples',nsampling,'Burnin',nwarmup, ...
        'VerbosityLevel',0);
    percnt(:,:,c) = chain;
end

%pooled draws (nchains*nsampling) x D
draws = reshape(permute(percnt,[1 3 2]),nsampling*nchains,D);

%convergence for the monitored parameters only
hdr = {'name','n_eff','r_hat'};
conv = hdr;
for k = 1:numel(monitor_idx)
    j = monitor_idx(k);
    x = squeeze(percnt(:,j,:));         % nsampling x nchains
    rhat = local_split_rhat(x);
    neff = local_ess(x);
    conv(end+1,1:3) = {monitor_names{k},neff,rhat}; %#ok<AGROW>
end
end

% ------------------------------------------------------------------------------
function rhat = local_split_rhat(x)
%Stan split-R-hat: split each chain in half, then between/within variance ratio.
[n,m] = size(x);
half = floor(n/2);
if half < 2
    rhat = Inf; return;
end
%2*m split chains of length `half`
s = [x(1:half,:), x(half+1:2*half,:)];   % half x (2m)
M = size(s,2);
cmeans = mean(s,1);
cvars  = var(s,0,1);
W = mean(cvars);
B = half*var(cmeans,0);                  % between-chain (scaled)
if W <= 0
    rhat = 1; return;
end
varplus = (half-1)/half*W + B/half;
rhat = sqrt(varplus/W);
end

% ------------------------------------------------------------------------------
function ess = local_ess(x)
%Effective sample size from the combined-chain autocorrelations (Stan-style:
%sum truncated when consecutive autocorrelation-pair sums go negative).
[n,m] = size(x);
N = n*m;
%within-chain mean-centered, then average autocorrelation across chains
acov = zeros(n,1);
for c = 1:m
    xc = x(:,c) - mean(x(:,c));
    a = local_autocov(xc);
    acov = acov + a;
end
acov = acov/m;
if acov(1) <= 0
    ess = N; return;
end
rho = acov/acov(1);
%Geyer initial positive sequence: sum pairs until a pair-sum <= 0
tau = 1;
t = 1;
while t+2 <= n
    pair = rho(t+1) + rho(t+2);
    if pair <= 0
        break;
    end
    tau = tau + 2*pair;
    t = t + 2;
end
ess = N/max(tau,1);
end

% ------------------------------------------------------------------------------
function a = local_autocov(xc)
%Biased autocovariance via FFT, lags 0..n-1.
n = numel(xc);
f = fft(xc,2^nextpow2(2*n));
ac = ifft(f.*conj(f));
a = real(ac(1:n))/n;
end
