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
%  parameters: rank-normalized/folded split Rhat, with n_eff=min(bulk,tail).
%  Data columns 4/5/6 retain bulk ESS, tail ESS and mean MCSE separately.
%  Header cell 4 remains the divergence slot; acceptance ratios retain their
%  existing positions, followed by a schema/provenance struct. Fixed-array
%  values are independently validated against R posterior 1.6.0 conventions.
%  Frozen/nonfinite/too-short chains and invalid acceptance ratios fail.
%  Acceptance below .1 still fails the existing owner-approved floor. The
%  existing Rhat<1.1 / ESS>=10*nchains gate is an early-workflow check only.
%
%DIAGNOSTICS SEAM (tests only): [~,conv] = psyrat_hmc_sample('diagnostics',x,names,accratios)
%  builds the convergence cell from a fixed draw array x (nsampling x nchains x
%  D), a cellstr of D names and a 1 x nchains vector of acceptance ratios
%  through the same local kernels the sampling path uses, so the kernels can
%  be pinned on known arrays (they are local functions and unreachable
%  otherwise). No sampling happens and draws is returned empty.

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

%diagnostics seam (see the header): the positional inputs are reused as
%x, names, accratios
if ischar(logpdfgrad) && strcmp(logpdfgrad,'diagnostics')
    draws = [];
    conv = local_convcell(start,nchains,nwarmup);
    return;
end

D = numel(start);
percnt = zeros(nsampling,D,nchains);   % per-chain draws for R-hat
accratios = zeros(1,nchains);
for c = 1:nchains
    rng(double(seed)+c-1,'twister');
    %a small random jitter of the start point decorrelates the chains so
    %split-R-hat is meaningful; the master+offset seed keeps it reproducible
    st = start(:) + 0.01*randn(D,1);
    smp = hmcSampler(logpdfgrad,st,'CheckGradient',false);
    smp = tuneSampler(smp,'VerbosityLevel',0, ...
        'NumStepSizeTuningIterations',max(50,round(nwarmup/2)));
    %the third output is the fraction of accepted proposals for this chain;
    %requesting it does not touch the random stream, so draws are unchanged
    [chain,~,accratios(c)] = drawSamples(smp,'NumSamples',nsampling, ...
        'Burnin',nwarmup,'VerbosityLevel',0);
    percnt(:,:,c) = chain;
end

%pooled draws (nchains*nsampling) x D
draws = reshape(permute(percnt,[1 3 2]),nsampling*nchains,D);

%convergence for the monitored parameters only: nsampling x nchains x Dmon
x = permute(percnt(:,monitor_idx,:),[1 3 2]);
conv = local_convcell(x,monitor_names,accratios);
end

% ------------------------------------------------------------------------------
function conv = local_convcell(x,names,accratios)
%Build the convergence cell from a draw array x (n x m x D), the D parameter
%names and the m acceptance ratios. Shared by the sampling path and the
%diagnostics seam so that what the tests pin is what production runs.
[~,m,D] = size(x);
if iscell(names) == false
    names = cellstr(names);
end
accratios = double(accratios(:)');
if numel(accratios) ~= m
    error('PsyRAT:hmcAcceptanceShape', ...
        'Expected one acceptance ratio per chain (%d), got %d.',m,numel(accratios));
end
conv = {'name','n_eff','r_hat'};
conv{1,4} = [];                       % divergence slot, stamped by psyrat_storeconv
conv{1,5} = 'accept_ratio';
conv(1,6:5+m) = num2cell(accratios);
%acceptance floor: a chain that rejected more than 90% of its proposals is
%stuck or nearly so, and its draws cannot be trusted whatever the R-hat says
badaccept = ~isfinite(accratios) | imag(accratios)~=0 | real(accratios)<0.1 | real(accratios)>1;
floorfail = any(badaccept);
if floorfail
    bad = find(badaccept);
    warning('PsyRAT:hmcAcceptanceFloor',['Native HMC chain(s) %s accepted ' ...
        'fewer than 10%% of proposals or reported invalid acceptance (%s); the fit is marked not ' ...
        'converged. How to fix: rerun with more warmup, or use CmdStan.'], ...
        mat2str(bad),mat2str(accratios(bad),3));
end
meta=psyrat_diagnostic_metadata({});
meta.format='native_rank_fold'; meta.precision_significant_digits=16;
meta.degraded=false; meta.degradation_reason='';
meta.method='rank/fold split Rhat; bulk/tail ESS; posterior 1.6.0 conventions';
meta.contract_status='no structural constants among unconstrained monitors';
meta.expected_names=names(:)';
meta.row_columns={'name','n_eff','r_hat','ess_bulk','ess_tail','mcse_mean'};
meta.acceptance_ratios=accratios;
meta.monitor_scope='selected unconstrained parameters; returned quantities appended by native engine';
conv{1,6+m}=meta;
for k=1:D
    d=psyrat_rank_diagnostics(x(:,:,k));
    neff=min(d.ess_bulk,d.ess_tail);
    if isnan(d.ess_bulk) || isnan(d.ess_tail); neff=NaN; end
    if floorfail; d.rhat=Inf; end
    conv(end+1,1:6)={names{k},neff,d.rhat,d.ess_bulk,d.ess_tail,d.mcse_mean}; %#ok<AGROW>
end
end
