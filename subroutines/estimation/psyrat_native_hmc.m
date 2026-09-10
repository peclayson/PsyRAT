function fit = psyrat_native_hmc(design,data,seed,nchains,nwarm,nsamp,opts)
%Native Hamiltonian Monte Carlo estimation engine (MATLAB hmcSampler), used as
%a faithful CmdStan fallback for the location-scale and multivariate designs
%that fitlme cannot fit.
%
%Each implemented design supplies, as a local function below, the log-posterior
%(with gradient) over its UNCONSTRAINED parameter vector built to match the
%corresponding generated Stan model term-for-term (same likelihood, same priors
%on the constrained parameters, plus the change-of-variables Jacobians for the
%constrained -> unconstrained map). The shared harness psyrat_hmc_sample runs the
%chains, and the draws are transformed back to the CmdStan-shaped parameter draw
%arrays the per-case extraction in psyrat_computevarcomp consumes, wrapped in a
%PsyRATNativeFit so .block()/.extract()/convstats() are unchanged.
%
%STATUS: Phase 2, incremental. A design is routed here ONLY if
%psyrat_native_supported lists 'hmc' for its key; any design not yet implemented
%hits the guard at the bottom and errors clearly (it is therefore unreachable in
%normal use until its log-posterior lands).
%
%Implemented (both are the one-facet location-scale model: person location +
%log-residual via a 2-D Cholesky, plus a crossed trial main effect on the mean):
% 'sserr'  - subject-level error variance (case 6). The base location-scale
%            model with NO dimension predictors (KDIM = 0). Matches
%            psyrat_build_stan_sserr.
% 'dynrel' - one-facet dynamic/conditional reliability (Rast & Clayson, case 11):
%            sserr EXTENDED with fixed dimension slopes on the mean and the
%            log-residual (KDIM = 1 or 3). Matches psyrat_build_stan_dynrel.
% Both share one local builder (local_hmc_locscale); 'sserr' is the KDIM = 0 case.
% 'diff'   - two-event difference-score location-scale model, no residual
%            covariance (case 7, diffrescor = 1). A univariate Gaussian over the
%            stacked event data with an event-indicator design (per-event mean b
%            and log-residual b_sigma) and TWO correlated 2-D random effects
%            (person and trial, each a 2x2 Cholesky). Matches psyrat_build_stan_diff.
% 'diff_rescor' - two-event difference-score model WITH residual covariance
%            (case 7, diffrescor = 2). A BIVARIATE Gaussian over the paired event
%            responses (meas1,meas2) with constant per-event residual SDs
%            exp(b_sigma) and a residual correlation (Lrescor), plus the same two
%            correlated 2-D person/trial random effects on the two event means.
%            Matches psyrat_build_stan_diff_rescor.
% 'diff_dynrel' - group-level dynamic/conditional reliability of a two-event
%            difference score (Rast & Clayson, case 12). The no-rescor difference
%            model ('diff') EXTENDED with event-specific fixed dimension slopes on
%            BOTH the mean (Xdim*b_dim) and the per-event log-residual
%            (Xdim*b_sigma_dim), exactly as 'dynrel' extends 'sserr'. Matches
%            psyrat_build_stan_diffdynrel.
% 'diff_dynrel_sserr' - SUBJECT-LEVEL dynamic difference score (case 13). The
%            group-level dynamic difference ('diff_dynrel') with the person
%            random-effect block WIDENED from 2-D to 4-D (a joint correlated block
%            over {location event1, location event2, scale event1, scale event2}),
%            so each participant has their own per-event residual via a 4x4 LKJ
%            Cholesky (psyrat_corr_chol_lkj). The trial block stays 2-D
%            (location-only). Matches psyrat_build_stan_diffdynrel_sserr.
% 'dynrel_trt' - two-facet (trial + occasion) dynamic reliability, NON-difference
%            (Rast & Clayson, case 14). The one-facet dynrel locscale model with
%            the single trial main effect GENERALIZED to FIVE scalar crossed-facet
%            mean random effects (occasion, trial, trial x person, occasion x
%            person, trial x occasion), each cauchy-scaled (priors.trt.*). The
%            residual scale carries only the person ind_sd. Matches
%            psyrat_build_stan_dynrel_trt.
% 'dod'    - task-free four-event difference-of-differences (case 9). A 4-cell
%            location-scale model: per-cell mean (b_cell) and log-residual
%            (b_sigma_cell), with person AND trial random effects on BOTH the
%            location (cells 1-4) and the scale (cells 1-4), i.e. an 8-D joint
%            block per grouping factor via an 8x8 LKJ Cholesky
%            (psyrat_corr_chol_lkj). The DoD contrast is applied downstream from
%            the per-cell var-covariances. Prior scales come from priors.dod
%            (attached via native_extra; the defaults reproduce the former
%            fixed brms constants). Matches psyrat_build_stan_dodiff.
% 'diff_dynrel_trt' - group-level dynamic difference, TWO-FACET (trial + occasion)
%            (Rast & Clayson, case 15). The one-facet group-level dynamic
%            difference ('diff_dynrel') generalized from TWO crossed 2-D blocks
%            (person, trial) to SIX (id, trl, occ, tid=trial x id, oid=occ x id,
%            to=trial x occ), each a 2x2 per-event Cholesky, plus the same
%            event-by-dimension slopes. Group-level (fixed-per-event) residual.
%            Matches psyrat_build_stan_diffdynrel_trt.
% 'diff_dynrel_sserr_trt' - SUBJECT-LEVEL two-facet dynamic difference (case 16).
%            The group-level two-facet dynamic difference ('diff_dynrel_trt') with
%            the PERSON (id) factor WIDENED from 2-D to a 4-D joint block (location
%            + log-residual per event via a 4x4 LKJ Cholesky, as case 13), so each
%            participant has their own per-event residual; the OTHER FIVE crossed
%            factors (trl/occ/tid/oid/to) stay 2-D location-only (as case 15).
%            Matches psyrat_build_stan_diffdynrel_sserr_trt.
% 'diff_dynrel_rescor' - CONCURRENT one-facet group-level dynamic difference (case
%            17). The bivariate-residual difference model ('diff_rescor', case 7*)
%            EXTENDED with per-event fixed dimension slopes on the mean (b_dim1/
%            b_dim2) AND on the per-event log-residual SD (b_sigma_dim1/
%            b_sigma_dim2), exactly as 'diff_dynrel' extends 'diff'. The per-event
%            residual SDs (and hence the bivariate residual Cholesky) are per-row
%            (dimension-dependent). Matches psyrat_build_stan_diffdynrel_rescor.
% 'diff_dynrel_trt_rescor' - CONCURRENT two-facet group-level dynamic difference
%            (case 18). The concurrent one-facet model ('diff_dynrel_rescor', case
%            17) with the person+trial 2-D blocks generalized to the SIX crossed
%            cross-condition 2-D factors (id/trl/occ/tid/oid/to, as case 15) in the
%            concurrent (paired) form; same per-row bivariate residual (Lrescor) +
%            per-event dimension slopes. Matches psyrat_build_stan_diffdynrel_trt_rescor.
% 'diff_dynrel_sserr_rescor' - CONCURRENT one-facet SUBJECT-LEVEL dynamic difference
%            (case 19). The concurrent one-facet group model ('diff_dynrel_rescor',
%            case 17) with the person block WIDENED from 2-D to a 4-D joint
%            location+scale block (a 4x4 LKJ Cholesky, as case 13), so each
%            participant has their own per-event residual SDs (the per-subject scale
%            REs enter the per-row sigma); the residual CORRELATION (Lrescor) stays
%            population. Matches psyrat_build_stan_diffdynrel_sserr_rescor.
% 'diff_dynrel_sserr_trt_rescor' - CONCURRENT two-facet SUBJECT-LEVEL dynamic
%            difference (case 20). The concurrent one-facet subject-level model
%            ('diff_dynrel_sserr_rescor', case 19) with the single trial 2-D block
%            generalized to the FIVE crossed cross-condition 2-D factors (trl/occ/
%            tid/oid/to, as case 16); the person (id) factor stays the 4-D joint
%            location+scale block. Same per-row bivariate residual (population
%            Lrescor) + per-event dimension slopes. Matches
%            psyrat_build_stan_diffdynrel_sserr_trt_rescor.
% 'diff_dynrel_sserr_rho' - CONCURRENT one-facet subject-level dynamic difference
%            with a PER-SUBJECT residual correlation (case 21). The concurrent
%            one-facet subject-level model ('diff_dynrel_sserr_rescor', case 19)
%            with the single POPULATION residual correlation replaced by a
%            HIERARCHICAL per-participant correlation rho[s] = tanh(rescor_mu +
%            sd_rescor * z_rescor[s]) (the case-8 rho hierarchy), so each
%            participant has their own residual SDs AND their own residual
%            correlation. Matches psyrat_build_stan_diffdynrel_sserr_rho.
% 'diff_dynrel_sserr_trt_rho' - CONCURRENT two-facet subject-level dynamic
%            difference with a PER-SUBJECT residual correlation (case 22). The
%            concurrent two-facet subject-level model ('diff_dynrel_sserr_trt_rescor',
%            case 20) with the single population residual correlation replaced by
%            the case-21 per-participant rho[s] hierarchy. Matches
%            psyrat_build_stan_diffdynrel_sserr_trt_rho.
% 'diff_sserr' - STATIC subject-level two-event difference (case 8). NOT dynamic
%            (no dimension slopes). Three SEPARATE 2-D blocks: person location
%            (L_id), trial (L_trl), and a per-subject SCALE block (L_id_sigma)
%            feeding the per-event residual log-SD; PLUS a per-subject residual
%            correlation rho[s]=tanh(rescor_mu+sd_rescor*z_rescor[s]) (toggled by
%            estimate_rescor) over PAIRED observations (bivariate) and a univariate
%            normal over UNPAIRED single-event observations. Matches
%            psyrat_build_stan_diff_sserr.
%
%Inputs mirror psyrat_run_native's hand-off:
% design - struct with .key (design id) and any extra descriptor fields the
%  model needs beyond the Stan data struct. For 'dynrel', design.priors carries
%  the configurable scale priors (priors.dynrel: .b, .b_sigma, .sig_trl).
% data   - the SAME Stan data struct the CmdStan path builds (NOBS, NSUB, NTRL,
%  KDIM, id, trl, Xdim_vec, meas).
% seed   - master RNG seed (per-chain seeds are seed+c-1 inside the harness).
% nchains, nwarm, nsamp - chains, warmup, post-warmup draws per chain.
% opts   - native-engine knobs; opts.hmc_gradient is 'auto' (default), 'ad'
%  (automatic differentiation via dlarray; requires Deep Learning Toolbox), or
%  'numer' (central finite differences; slow, no extra toolbox).
%
%See also: psyrat_run_native, psyrat_native_supported, PsyRATNativeFit,
% psyrat_hmc_sample, hmcSampler

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

key = '';
if isstruct(design) && isfield(design,'key')
    key = char(design.key);
end

switch key
    case {'dynrel','sserr'}
        fit = local_hmc_locscale(key,design,data,seed,nchains,nwarm,nsamp,opts);
    case 'trt'
        % plain test-retest crossed variance-components (case 5): a faithful
        % Bayesian alternative to the fitlme 'trt' fallback (preferred over it
        % by the auto-cascade), matching psyrat_build_stan_trt.
        fit = local_hmc_trt(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'diff'
        fit = local_hmc_diff(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'diff_rescor'
        fit = local_hmc_diff_rescor(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'diff_dynrel'
        fit = local_hmc_diffdynrel(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'diff_dynrel_sserr'
        fit = local_hmc_diffdynrel_sserr(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'dynrel_trt'
        fit = local_hmc_dynrel_trt(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'dod'
        fit = local_hmc_dodiff(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'diff_trt'
        % static two-facet difference (case 10, non-concurrent): the KDIM = 0
        % reduction of the dynamic two-facet difference (case 15).
        fit = local_hmc_diffdynrel_trt(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'diff_trt_rescor'
        % static concurrent two-facet difference (case 10, diffrescor = 2): the
        % KDIM = 0 reduction of the dynamic concurrent two-facet difference (case 18).
        fit = local_hmc_diffdynrel_trt_rescor(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'diff_dynrel_trt'
        fit = local_hmc_diffdynrel_trt(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'diff_dynrel_sserr_trt'
        fit = local_hmc_diffdynrel_sserr_trt(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'diff_dynrel_rescor'
        fit = local_hmc_diffdynrel_rescor(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'diff_dynrel_trt_rescor'
        fit = local_hmc_diffdynrel_trt_rescor(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'diff_dynrel_sserr_rescor'
        fit = local_hmc_diffdynrel_sserr_rescor(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'diff_dynrel_sserr_trt_rescor'
        fit = local_hmc_diffdynrel_sserr_trt_rescor(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'diff_dynrel_sserr_rho'
        fit = local_hmc_diffdynrel_sserr_rho(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'diff_dynrel_sserr_trt_rho'
        fit = local_hmc_diffdynrel_sserr_trt_rho(design,data,seed,nchains,nwarm,nsamp,opts);
    case 'diff_sserr'
        fit = local_hmc_diff_sserr(design,data,seed,nchains,nwarm,nsamp,opts);
    otherwise
        error('PsyRAT:hmcNotImplemented',['The native HMC engine is not yet ' ...
            'implemented for design ''%s''. This is an internal guard: the ' ...
            'engine registry (psyrat_native_supported) should not route this ' ...
            'design to the HMC engine until its log-posterior is built.'],key);
end
end

% ===========================================================================
% Designs: sserr (case 6) and dynrel (case 11) - one-facet location-scale
% ===========================================================================
function fit = local_hmc_locscale(key,design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the one-facet location-scale model. With
%dimension predictors (KDIM > 0) this is dynrel (psyrat_build_stan_dynrel); with
%none (KDIM = 0) it is the base subject-level error-variance model
%(psyrat_build_stan_sserr). The two differ only by the Xdim*b / Xdim*b_sigma
%terms, which vanish when KDIM = 0, so a single log-posterior covers both.

% --- assemble the model-data struct the log-posterior closes over ----------
m = struct();
m.NOBS = double(data.NOBS);
m.NSUB = double(data.NSUB);
m.NTRL = double(data.NTRL);
m.meas = double(data.meas(:));
m.id   = double(data.id(:));      % sequential 1..NSUB
m.trl  = double(data.trl(:));     % sequential 1..NTRL
% dimension design: dynrel passes Xdim_vec/KDIM; sserr (KDIM = 0) passes neither,
% so the design is an NOBS x 0 matrix and the slope terms drop out. The flatten
% is column-major (Xdim_vec = Xdim(:)); reshape reads column-major to match.
if isfield(data,'KDIM') && ~isempty(data.KDIM)
    m.KDIM = double(data.KDIM);
else
    m.KDIM = 0;
end
if m.KDIM > 0
    m.Xdim = reshape(double(data.Xdim_vec(:)),m.NOBS,m.KDIM);
else
    m.Xdim = zeros(m.NOBS,0);
end

% configurable scale priors (defaults mirror psyrat_default_priors.dynrel if the
% caller did not attach them; the fixed student_t constants are baked into the
% log-posterior to match the Stan builder).
P = struct('b',10,'b_sigma',1,'sig_trl',10);
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors)
    pr = design.priors;
    if isfield(pr,'b');       P.b = double(pr.b);             end
    if isfield(pr,'b_sigma'); P.b_sigma = double(pr.b_sigma); end
    if isfield(pr,'sig_trl'); P.sig_trl = double(pr.sig_trl); end
end
m.P = P;

% --- unconstrained parameter index map -------------------------------------
m = local_dynrel_index(m);

% --- gradient method: AD (default when available) or numeric fallback -------
gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
% Probe the Deep Learning Toolbox AD stack. NB: dlgradient is a method only
% callable inside dlfeval, so exist('dlgradient','file') is ALWAYS 0 even when
% the toolbox is installed; probe dlfeval (a real file, the AD entry point)
% instead, matching psyrat_estimation_engines_available's avail.ad.
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: ' ...
                'install the Deep Learning Toolbox or use hmc_gradient ' ...
                '''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end

% the log-posterior over the unconstrained vector, returning [logp, grad]
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_dynrel_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_dynrel_logp(t,m), th);
end

% --- start point + a cheap gradient self-check -----------------------------
start = local_dynrel_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

% --- run the shared sampling harness ---------------------------------------
[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

% --- transform draws back to the CmdStan-shaped parameter arrays ------------
drawstruct = local_dynrel_draws(draws,m);

meta = struct('design',key,'gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_dynrel_index(m)
%Lay out theta = [Intercept; Intercept_sigma; b(KDIM); b_sigma(KDIM);
%  log_gro_sds(2); z_rho; log_sig_trl; gro_effs_stndzd(2*NSUB, column-major);
%  trl_raw(NTRL)].
K = m.KDIM; NSUB = m.NSUB; NTRL = m.NTRL;
i = 0;
m.ix.Intercept = i+1;            i = i+1;
m.ix.Isig      = i+1;            i = i+1;   % Intercept_sigma
m.ix.b         = i+(1:K);        i = i+K;
m.ix.b_sigma   = i+(1:K);        i = i+K;
m.ix.lgro      = i+(1:2);        i = i+2;   % log(gro_sds) = log([sigma_p,sigma_dp])
m.ix.zrho      = i+1;            i = i+1;   % atanh(rho)
m.ix.lstrl     = i+1;            i = i+1;   % log(sig_trl)
m.ix.zsub      = i+(1:2*NSUB);   i = i+2*NSUB; % person raw z's (reshape 2xNSUB)
m.ix.ztrl      = i+(1:NTRL);     i = i+NTRL;
m.D = i;
% monitor the population-level parameters (not the per-subject/per-trial raw
% effects) for convergence reporting.
m.mon_idx = [m.ix.Intercept, m.ix.Isig, m.ix.b, m.ix.b_sigma, ...
    m.ix.lgro, m.ix.zrho, m.ix.lstrl];
names = [{'Intercept','Intercept_sigma'}, ...
    local_numbered('b',K), local_numbered('b_sigma',K), ...
    {'log_sigma_p','log_sigma_dp','z_rho','log_sig_trl'}];
m.mon_names = names;
end

function c = local_numbered(base,k)
%{'b1','b2',...} helper for the convergence labels.
c = cell(1,k);
for j = 1:k
    c{j} = sprintf('%s%d',base,j);
end
end

function start = local_dynrel_start(m)
%Sensible start: intercepts at the data scale, raw effects at 0, log-scales at
%rough SD guesses, z_rho at 0 (rho = 0).
start = zeros(m.D,1);
smeas = std(m.meas);
start(m.ix.Intercept) = mean(m.meas);
start(m.ix.Isig)      = log(max(smeas/2,1e-3));
start(m.ix.lgro(1))   = log(max(smeas/2,1e-3));   % between-person mean SD
start(m.ix.lgro(2))   = log(0.3);                 % between-person log-resid SD
start(m.ix.zrho)      = 0;
start(m.ix.lstrl)     = log(max(smeas/2,1e-3));
% b, b_sigma, gro_effs_stndzd, trl_raw start at 0
end

% ---------------------------------------------------------------------------
function lp = local_dynrel_logp(theta,m)
%Unconstrained log-posterior (up to an additive constant). Works on either a
%plain double theta (numeric-gradient path) or a dlarray theta (AD path), since
%it uses only elementwise/linear-algebra ops dlarray supports. Matches
%psyrat_build_stan_dynrel: a Gaussian location-scale likelihood with the same
%priors on the constrained parameters plus the Jacobians for the exp() and tanh()
%maps. Parameter-independent normalizing constants are dropped (irrelevant to HMC
%and to the constrained-space posterior).
ix = m.ix;
Intercept = theta(ix.Intercept);
Isig      = theta(ix.Isig);
b         = theta(ix.b);              % KDIM x 1
b_sigma   = theta(ix.b_sigma);        % KDIM x 1
gro       = exp(theta(ix.lgro));      % 2 x 1 (sigma_p, sigma_dp), >0
rho       = tanh(theta(ix.zrho));     % in (-1,1)
sig_trl   = exp(theta(ix.lstrl));     % >0
zsub      = reshape(theta(ix.zsub),2,m.NSUB);
z1        = zsub(1,:).';              % NSUB x 1 (person mean raw)
z2        = zsub(2,:).';              % NSUB x 1 (person log-resid raw)
ztrl      = theta(ix.ztrl);           % NTRL x 1

% non-centered person effects via the 2-D Cholesky diag(gro)*[1 0; rho sqrt(1-rho^2)]
ind_bs = gro(1).*z1;                                   % NSUB x 1
ind_sd = gro(2).*(rho.*z1 + sqrt(1-rho.^2).*z2);       % NSUB x 1
trl_b  = sig_trl.*ztrl;                                % NTRL x 1

% per-observation mean and log-residual (dimension-conditioned)
mu     = Intercept + m.Xdim*b + ind_bs(m.id) + trl_b(m.trl);
logsig = Isig + m.Xdim*b_sigma + ind_sd(m.id);
sig    = exp(logsig);

% Gaussian observation log-likelihood (keep the param-dependent -sum(logsig))
ll = -sum(logsig) - sum((m.meas - mu).^2 ./ (2*sig.^2));

% priors on the constrained parameters (kernels; fixed student_t constants
% identical to the Stan builder) + change-of-variables Jacobians
P = m.P;
lp = ll ...
    + local_student_t(Intercept,3,15.7,12.4) ...        % Intercept
    + local_student_t(Isig,3,0,2.5) ...                 % Intercept_sigma
    + local_student_t(gro(1),3,0,12.4) + theta(ix.lgro(1)) ... % half-t + Jac
    + local_student_t(gro(2),3,0,12.4) + theta(ix.lgro(2)) ...
    + local_cauchy(sig_trl,0,P.sig_trl) + theta(ix.lstrl) ...  % half-cauchy + Jac
    + sum(local_normal(b,0,P.b)) ...                    % mean dimension slopes
    + sum(local_normal(b_sigma,0,P.b_sigma)) ...        % scale dimension slopes
    + log(1 - rho.^2) ...                               % LKJ(1) uniform + tanh Jac
    - 0.5*sum(z1.^2) - 0.5*sum(z2.^2) ...               % std_normal person raw
    - 0.5*sum(ztrl.^2);                                 % std_normal trial raw
end

% ---- prior log-density kernels (drop parameter-independent constants) ------
function v = local_normal(x,mu,s)
v = -0.5*((x-mu)./s).^2;
end
function v = local_student_t(x,nu,mu,s)
v = -((nu+1)/2).*log(1 + ((x-mu)./s).^2 ./ nu);
end
function v = local_cauchy(x,mu,s)
v = -log(1 + ((x-mu)./s).^2);
end

% ---------------------------------------------------------------------------
function out = local_dynrel_draws(draws,m)
%Transform the pooled unconstrained draws into the CmdStan-shaped parameter draw
%arrays case 11 extracts (Intercept, Intercept_sigma, ind_bs, ind_sd, gro_sds,
%sig_trl, b, b_sigma, chol_corrmat).
ix = m.ix;
ND = size(draws,1);
out = struct();
out.Intercept       = draws(:,ix.Intercept);          % ND x 1
out.Intercept_sigma = draws(:,ix.Isig);               % ND x 1
out.b               = draws(:,ix.b);                  % ND x KDIM
out.b_sigma         = draws(:,ix.b_sigma);            % ND x KDIM
gro                 = exp(draws(:,ix.lgro));          % ND x 2
out.gro_sds         = gro;
out.sig_trl         = exp(draws(:,ix.lstrl));         % ND x 1
rho                 = tanh(draws(:,ix.zrho));         % ND x 1

% chol_corrmat as CmdStan returns a matrix[2,2]: ND x 2 x 2, lower-triangular
% [1 0; rho sqrt(1-rho^2)] (downstream reads rho at (:,2,1)).
chol = zeros(ND,2,2);
chol(:,1,1) = 1;
chol(:,2,1) = rho;
chol(:,2,2) = sqrt(1 - rho.^2);
out.chol_corrmat = chol;

% transformed person effects ind_bs / ind_sd (ND x NSUB), reconstructed from the
% raw z draws and the per-draw scales/correlation (matches the Stan transformed
% parameters block).
Z  = draws(:,ix.zsub);          % ND x 2NSUB, columns [z1_1 z2_1 z1_2 z2_2 ...]
Z1 = Z(:,1:2:end);              % ND x NSUB (person mean raw)
Z2 = Z(:,2:2:end);              % ND x NSUB (person log-resid raw)
out.ind_bs = gro(:,1).*Z1;
out.ind_sd = gro(:,2).*(rho.*Z1 + sqrt(1 - rho.^2).*Z2);
end

% ===========================================================================
% Design: trt (case 5) - plain test-retest crossed variance components
% ===========================================================================
function fit = local_hmc_trt(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the plain test-retest variance-components
%model, matching psyrat_build_stan_trt: a plain-Gaussian crossed random-effects
%ANOVA (Persons x Trials x Occasions) with a fixed intercept, six crossed scalar
%random-effect components (id, occ, trl, trlxid, occxid, trlxocc) in the
%non-centered form (term = sig_k * raw_k, raw_k ~ normal(0,1)), and a single
%constant residual SD sig_err. The seven SDs carry half-Cauchy priors
%(priors.trt.*); the intercept is improper-flat (as in the Stan model). This is
%the faithful Bayesian alternative the auto-cascade prefers over the fitlme 'trt'
%fallback (which is REML + a parametric bootstrap and is weak when sig_occ is
%poorly identified at small occasion counts).

m = struct();
m.NOBS = double(data.NOBS);
m.NSUB = double(data.NSUB);
m.NOCC = double(data.NOCC);
m.NTRL = double(data.NTRL);
m.NTID = double(data.NTID);
m.NOID = double(data.NOID);
m.NTO  = double(data.NTO);
m.meas    = double(data.meas(:));
m.id      = double(data.id(:));        % 1..NSUB
m.occ     = double(data.occ(:));       % 1..NOCC
m.trl     = double(data.trl(:));       % 1..NTRL
m.trlxid  = double(data.trlxid(:));    % 1..NTID
m.occxid  = double(data.occxid(:));    % 1..NOID
m.trlxocc = double(data.trlxocc(:));   % 1..NTO

% half-Cauchy SD prior scales, ORDERED [id occ trl err trlxid occxid trlxocc] to
% match the lsig index block. Defaults mirror psyrat_default_priors.trt; the
% caller attaches priors.trt via native_extra.
P = [10 5 10 20 5 5 1];
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors) && ...
        isfield(design.priors,'trt') && isstruct(design.priors.trt)
    tr = design.priors.trt;
    flds = {'sig_id','sig_occ','sig_trl','sig_err','sig_trlxid','sig_occxid','sig_trlxocc'};
    for k = 1:7
        if isfield(tr,flds{k}); P(k) = double(tr.(flds{k})); end
    end
end
m.Pscale = P(:);   % 7 x 1, aligned with sig = exp(theta(ix.lsig))

m = local_trt_index(m);

gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_trt_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_trt_logp(t,m), th);
end

start = local_trt_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_trt_draws(draws,m);
% report the actual design key ('trt') rather than a hardcoded label.
designname = 'trt';
if isstruct(design) && isfield(design,'key') && ~isempty(design.key)
    designname = char(design.key);
end
meta = struct('design',designname,'gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_trt_index(m)
%theta = [pop_int; log_sig(7) [id occ trl err trlxid occxid trlxocc];
%         id_raw(NSUB); occ_raw(NOCC); trl_raw(NTRL);
%         trlxid_raw(NTID); occxid_raw(NOID); trlxocc_raw(NTO)].
i = 0;
m.ix.pop_int     = i+1;                 i = i+1;
m.ix.lsig        = i+(1:7);             i = i+7;
m.ix.id_raw      = i+(1:m.NSUB);        i = i+m.NSUB;
m.ix.occ_raw     = i+(1:m.NOCC);        i = i+m.NOCC;
m.ix.trl_raw     = i+(1:m.NTRL);        i = i+m.NTRL;
m.ix.trlxid_raw  = i+(1:m.NTID);        i = i+m.NTID;
m.ix.occxid_raw  = i+(1:m.NOID);        i = i+m.NOID;
m.ix.trlxocc_raw = i+(1:m.NTO);         i = i+m.NTO;
m.D = i;
% monitor the population intercept + the seven log-SDs for convergence.
m.mon_idx = [m.ix.pop_int, m.ix.lsig];
m.mon_names = {'pop_int','log_sig_id','log_sig_occ','log_sig_trl', ...
    'log_sig_err','log_sig_trlxid','log_sig_occxid','log_sig_trlxocc'};
end

function start = local_trt_start(m)
%Intercept at the grand mean; log-SDs at rough guesses (residual a touch larger
%than the crossed components); raw effects at 0.
start = zeros(m.D,1);
smeas = std(m.meas);
start(m.ix.pop_int) = mean(m.meas);
start(m.ix.lsig)    = log(max(smeas/3,1e-3));     % all seven SDs
start(m.ix.lsig(4)) = log(max(smeas/2,1e-3));     % residual sig_err
end

function lp = local_trt_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching psyrat_build_stan_trt:
%a plain-Gaussian crossed random-effects ANOVA with a flat intercept, half-Cauchy
%SD priors, and non-centered (sig_k * raw_k) random effects.
ix = m.ix;
pop_int = theta(ix.pop_int);
sig     = exp(theta(ix.lsig));   % 7: [id occ trl err trlxid occxid trlxocc]
id_raw      = theta(ix.id_raw);
occ_raw     = theta(ix.occ_raw);
trl_raw     = theta(ix.trl_raw);
trlxid_raw  = theta(ix.trlxid_raw);
occxid_raw  = theta(ix.occxid_raw);
trlxocc_raw = theta(ix.trlxocc_raw);

% non-centered crossed terms, gathered to the per-observation rows
id_terms      = sig(1).*id_raw;
occ_terms     = sig(2).*occ_raw;
trl_terms     = sig(3).*trl_raw;
trlxid_terms  = sig(5).*trlxid_raw;
occxid_terms  = sig(6).*occxid_raw;
trlxocc_terms = sig(7).*trlxocc_raw;
y_hat = pop_int + id_terms(m.id) + occ_terms(m.occ) + trl_terms(m.trl) ...
    + trlxid_terms(m.trlxid) + occxid_terms(m.occxid) + trlxocc_terms(m.trlxocc);

sig_err = sig(4);
ll = -m.NOBS*log(sig_err) - sum((m.meas - y_hat).^2) ./ (2*sig_err.^2);

% standard-normal raw effects
lp_raw = -0.5*( sum(id_raw.^2) + sum(occ_raw.^2) + sum(trl_raw.^2) ...
    + sum(trlxid_raw.^2) + sum(occxid_raw.^2) + sum(trlxocc_raw.^2) );

% half-Cauchy SD priors (student_t df=1) + the exp() Jacobian (sum of log_sig).
% The intercept is improper-flat, so it contributes nothing.
lp_sd = sum(local_student_t(sig,1,0,m.Pscale)) + sum(theta(ix.lsig));

lp = ll + lp_raw + lp_sd;
end

function out = local_trt_draws(draws,m)
%Transform draws to the CmdStan-shaped scalars case 5 extracts: pop_int and the
%seven component SDs (sig_id/occ/trl/err/trlxid/occxid/trlxocc).
ix = m.ix;
out = struct();
out.pop_int     = draws(:,ix.pop_int);
sig             = exp(draws(:,ix.lsig));   % ND x 7
out.sig_id      = sig(:,1);
out.sig_occ     = sig(:,2);
out.sig_trl     = sig(:,3);
out.sig_err     = sig(:,4);
out.sig_trlxid  = sig(:,5);
out.sig_occxid  = sig(:,6);
out.sig_trlxocc = sig(:,7);
end

% ===========================================================================
% Design: diff (case 7) - two-event difference score, no residual covariance
% ===========================================================================
function fit = local_hmc_diff(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the two-event difference-score location-
%scale model WITHOUT residual covariance, matching psyrat_build_stan_diff: a
%univariate Gaussian over the stacked event data with an event-indicator design
%(per-event mean b and log-residual b_sigma) and two correlated 2-D random
%effects (person r_1, trial r_2), each a non-centered 2x2 Cholesky.

m = struct();
m.NOBS  = double(data.NOBS);
m.NSUB  = double(data.NSUB);
m.NTRL  = double(data.NTRL);
m.meas  = double(data.meas(:));
m.meas1 = double(data.meas1(:));   % event-1 indicator (1/0)
m.meas2 = double(data.meas2(:));   % event-2 indicator (1/0)
m.id    = double(data.ID(:));      % sequential 1..NSUB
m.trl   = double(data.TRL(:));     % sequential 1..NTRL

% configurable scale priors (priors.diff): b ~ normal(0,b);
% b_sigma ~ student_t(3,0,b_sigma); sd_id/sd_trl ~ half-student_t(3,0,.). The
% defaults match psyrat_default_priors.diff so a default run is correct even if
% the caller did not attach priors.
P = struct('b',5,'b_sigma',10,'sd_id',10,'sd_trl',10);
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors)
    pr = design.priors;
    fn = {'b','b_sigma','sd_id','sd_trl'};
    for k = 1:numel(fn)
        if isfield(pr,fn{k}); P.(fn{k}) = double(pr.(fn{k})); end
    end
end
m.P = P;

m = local_diff_index(m);

% gradient method (same AD/numer selection as the locscale builder)
gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
% Probe the Deep Learning Toolbox AD stack. NB: dlgradient is a method only
% callable inside dlfeval, so exist('dlgradient','file') is ALWAYS 0 even when
% the toolbox is installed; probe dlfeval (a real file, the AD entry point)
% instead, matching psyrat_estimation_engines_available's avail.ad.
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_diff_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_diff_logp(t,m), th);
end

start = local_diff_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_diff_draws(draws,m);
meta = struct('design','diff','gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_diff_index(m)
%theta = [b(2); b_sigma(2); log_sd_id(2); z_rho_id; log_sd_trl(2); z_rho_trl;
%         z_id(2*NSUB, column-major); z_trl(2*NTRL, column-major)].
NSUB = m.NSUB; NTRL = m.NTRL; i = 0;
m.ix.b        = i+(1:2);      i = i+2;
m.ix.b_sigma  = i+(1:2);      i = i+2;
m.ix.lsd_id   = i+(1:2);      i = i+2;   % log(sd_id)
m.ix.zrho_id  = i+1;          i = i+1;   % atanh(rho_id)
m.ix.lsd_trl  = i+(1:2);      i = i+2;   % log(sd_trl)
m.ix.zrho_trl = i+1;          i = i+1;   % atanh(rho_trl)
m.ix.z_id     = i+(1:2*NSUB); i = i+2*NSUB;
m.ix.z_trl    = i+(1:2*NTRL); i = i+2*NTRL;
m.D = i;
m.mon_idx = [m.ix.b, m.ix.b_sigma, m.ix.lsd_id, m.ix.zrho_id, ...
    m.ix.lsd_trl, m.ix.zrho_trl];
m.mon_names = {'b1','b2','b_sigma1','b_sigma2','log_sd_id1','log_sd_id2', ...
    'z_rho_id','log_sd_trl1','log_sd_trl2','z_rho_trl'};
end

function start = local_diff_start(m)
%Per-event means at the event-conditional data means; log-SDs at rough guesses.
start = zeros(m.D,1);
m1 = m.meas1 > 0; m2 = m.meas2 > 0;
start(m.ix.b(1)) = mean(m.meas(m1));
start(m.ix.b(2)) = mean(m.meas(m2));
s1 = std(m.meas(m1)); s2 = std(m.meas(m2));
start(m.ix.b_sigma(1)) = log(max(s1/2,1e-3));
start(m.ix.b_sigma(2)) = log(max(s2/2,1e-3));
start(m.ix.lsd_id)  = log(max((s1+s2)/4,1e-3));
start(m.ix.lsd_trl) = log(max((s1+s2)/4,1e-3));
% z_rho_id, z_rho_trl, z_id, z_trl start at 0
end

function lp = local_diff_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching psyrat_build_stan_diff.
ix = m.ix;
b       = theta(ix.b);              % 2 (per-event mean)
b_sigma = theta(ix.b_sigma);        % 2 (per-event log-residual SD)
sd_id   = exp(theta(ix.lsd_id));    % 2 (person event-SDs), >0
rho_id  = tanh(theta(ix.zrho_id));
sd_trl  = exp(theta(ix.lsd_trl));   % 2 (trial event-SDs), >0
rho_trl = tanh(theta(ix.zrho_trl));
z_id  = reshape(theta(ix.z_id),2,m.NSUB);
z_trl = reshape(theta(ix.z_trl),2,m.NTRL);

% non-centered person/trial event-effects via 2-D Cholesky diag(sd)*[1 0; rho .]
r1_1 = sd_id(1).*z_id(1,:).';                                      % NSUB x1 (event1)
r1_2 = sd_id(2).*(rho_id.*z_id(1,:).' + sqrt(1-rho_id.^2).*z_id(2,:).');
r2_1 = sd_trl(1).*z_trl(1,:).';                                    % NTRL x1
r2_2 = sd_trl(2).*(rho_trl.*z_trl(1,:).' + sqrt(1-rho_trl.^2).*z_trl(2,:).');

m1 = m.meas1; m2 = m.meas2;
% event-indicator design: event-1 rows pick b(1)/b_sigma(1)/r*_1, event-2 b(2)/..
mu = m1.*b(1) + m2.*b(2) ...
   + r1_1(m.id).*m1 + r1_2(m.id).*m2 ...
   + r2_1(m.trl).*m1 + r2_2(m.trl).*m2;
logsig = m1.*b_sigma(1) + m2.*b_sigma(2);
sig = exp(logsig);
ll = -sum(logsig) - sum((m.meas - mu).^2 ./ (2*sig.^2));

P = m.P;
lp = ll ...
    + sum(local_normal(b,0,P.b)) ...                          % b ~ normal(0,5)
    + sum(local_student_t(b_sigma,3,0,P.b_sigma)) ...         % b_sigma ~ t(3,0,10)
    + local_student_t(sd_id(1),3,0,P.sd_id)   + theta(ix.lsd_id(1)) ...  % half-t + Jac
    + local_student_t(sd_id(2),3,0,P.sd_id)   + theta(ix.lsd_id(2)) ...
    + local_student_t(sd_trl(1),3,0,P.sd_trl) + theta(ix.lsd_trl(1)) ...
    + local_student_t(sd_trl(2),3,0,P.sd_trl) + theta(ix.lsd_trl(2)) ...
    + 2*log(1 - rho_id.^2) ...                                % LKJ(2) + tanh Jac
    + 2*log(1 - rho_trl.^2) ...
    - 0.5*sum(z_id(:).^2) - 0.5*sum(z_trl(:).^2);             % std_normal raw
end

function out = local_diff_draws(draws,m)
%Transform draws to the CmdStan-shaped arrays case 7 (no-rescor) extracts:
%b, b_sigma (ND x2), sd_id, sd_trl (ND x2), id_varcov, trl_varcov (ND x2x2).
ix = m.ix;
out = struct();
out.b       = draws(:,ix.b);          % ND x2
out.b_sigma = draws(:,ix.b_sigma);    % ND x2
sd_id  = exp(draws(:,ix.lsd_id));     % ND x2
sd_trl = exp(draws(:,ix.lsd_trl));    % ND x2
out.sd_id  = sd_id;
out.sd_trl = sd_trl;
rho_id  = tanh(draws(:,ix.zrho_id));
rho_trl = tanh(draws(:,ix.zrho_trl));
% id_varcov / trl_varcov = quad_form_diag(corr, sd): [[s1^2, r*s1*s2],[.,s2^2]]
out.id_varcov  = local_varcov2(sd_id,rho_id);
out.trl_varcov = local_varcov2(sd_trl,rho_trl);
end

% ===========================================================================
% Design: diff_rescor (case 7, diffrescor = 2) - two-event difference WITH
% residual covariance
% ===========================================================================
function fit = local_hmc_diff_rescor(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the two-event difference-score model WITH
%residual covariance, matching psyrat_build_stan_diff_rescor: a BIVARIATE
%Gaussian over the paired event responses (meas1,meas2) with constant per-event
%residual SDs exp(b_sigma) and a residual correlation (Lrescor), plus the same
%two correlated 2-D person/trial random effects (r_1, r_2) on the two event
%means as the no-rescor sibling. The data are PAIRED here (one row per
%person-trial pair, meas1/meas2 the two event responses), not stacked.

m = struct();
m.NOBS  = double(data.NOBS);
m.NSUB  = double(data.NSUB);
m.NTRL  = double(data.NTRL);
m.meas1 = double(data.meas1(:));   % event-1 response (paired)
m.meas2 = double(data.meas2(:));   % event-2 response (paired)
m.id    = double(data.ID(:));      % sequential 1..NSUB
m.trl   = double(data.TRL(:));     % sequential 1..NTRL

% configurable scale priors (priors.diff): b ~ normal(0,b);
% b_sigma ~ student_t(3,0,b_sigma); sd_id/sd_trl ~ half-student_t(3,0,.). The
% residual correlation (Lrescor) is LKJ(2) with no configurable scale, matching
% the hardcoded lkj_corr_cholesky_lpdf(Lrescor | 2) in the Stan builder. The
% defaults match psyrat_default_priors.diff so a default run is correct even if
% the caller did not attach priors.
P = struct('b',5,'b_sigma',10,'sd_id',10,'sd_trl',10);
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors)
    pr = design.priors;
    fn = {'b','b_sigma','sd_id','sd_trl'};
    for k = 1:numel(fn)
        if isfield(pr,fn{k}); P.(fn{k}) = double(pr.(fn{k})); end
    end
end
m.P = P;

m = local_diff_rescor_index(m);

% gradient method (same AD/numer selection as the other builders)
gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
% Probe the Deep Learning Toolbox AD stack. NB: dlgradient is a method only
% callable inside dlfeval, so exist('dlgradient','file') is ALWAYS 0 even when
% the toolbox is installed; probe dlfeval (a real file, the AD entry point)
% instead, matching psyrat_estimation_engines_available's avail.ad.
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_diff_rescor_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_diff_rescor_logp(t,m), th);
end

start = local_diff_rescor_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_diff_rescor_draws(draws,m);
meta = struct('design','diff_rescor','gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_diff_rescor_index(m)
%theta = [b(2); b_sigma(2); log_sd_id(2); z_rho_id; log_sd_trl(2); z_rho_trl;
%         z_rho_res; z_id(2*NSUB, column-major); z_trl(2*NTRL, column-major)].
%Identical to local_diff_index plus the residual-correlation parameter z_rho_res.
NSUB = m.NSUB; NTRL = m.NTRL; i = 0;
m.ix.b        = i+(1:2);      i = i+2;
m.ix.b_sigma  = i+(1:2);      i = i+2;
m.ix.lsd_id   = i+(1:2);      i = i+2;   % log(sd_id)
m.ix.zrho_id  = i+1;          i = i+1;   % atanh(rho_id)
m.ix.lsd_trl  = i+(1:2);      i = i+2;   % log(sd_trl)
m.ix.zrho_trl = i+1;          i = i+1;   % atanh(rho_trl)
m.ix.zrho_res = i+1;          i = i+1;   % atanh(rho_res), residual correlation
m.ix.z_id     = i+(1:2*NSUB); i = i+2*NSUB;
m.ix.z_trl    = i+(1:2*NTRL); i = i+2*NTRL;
m.D = i;
m.mon_idx = [m.ix.b, m.ix.b_sigma, m.ix.lsd_id, m.ix.zrho_id, ...
    m.ix.lsd_trl, m.ix.zrho_trl, m.ix.zrho_res];
m.mon_names = {'b1','b2','b_sigma1','b_sigma2','log_sd_id1','log_sd_id2', ...
    'z_rho_id','log_sd_trl1','log_sd_trl2','z_rho_trl','z_rho_res'};
end

function start = local_diff_rescor_start(m)
%Per-event means at the per-event data means; log-SDs at rough guesses splitting
%each event's marginal variance across the residual and the two random effects.
start = zeros(m.D,1);
start(m.ix.b(1)) = mean(m.meas1);
start(m.ix.b(2)) = mean(m.meas2);
s1 = std(m.meas1); s2 = std(m.meas2);
start(m.ix.b_sigma(1)) = log(max(s1/2,1e-3));
start(m.ix.b_sigma(2)) = log(max(s2/2,1e-3));
start(m.ix.lsd_id)  = log(max((s1+s2)/4,1e-3));
start(m.ix.lsd_trl) = log(max((s1+s2)/4,1e-3));
% z_rho_id, z_rho_trl, z_rho_res, z_id, z_trl start at 0
end

function lp = local_diff_rescor_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching
%psyrat_build_stan_diff_rescor. Same person/trial random-effect structure and
%priors as local_diff_logp, but the likelihood is a BIVARIATE normal over the
%paired (meas1,meas2) with constant per-event residual SDs exp(b_sigma) and a
%residual correlation rho_res (Lrescor).
ix = m.ix;
b       = theta(ix.b);              % 2 (per-event mean)
b_sigma = theta(ix.b_sigma);        % 2 (per-event log residual SD)
sd_id   = exp(theta(ix.lsd_id));    % 2 (person event-SDs), >0
rho_id  = tanh(theta(ix.zrho_id));
sd_trl  = exp(theta(ix.lsd_trl));   % 2 (trial event-SDs), >0
rho_trl = tanh(theta(ix.zrho_trl));
rho_res = tanh(theta(ix.zrho_res)); % residual correlation between the events
z_id  = reshape(theta(ix.z_id),2,m.NSUB);
z_trl = reshape(theta(ix.z_trl),2,m.NTRL);

% non-centered person/trial event-effects via 2-D Cholesky diag(sd)*[1 0; rho .]
r1_1 = sd_id(1).*z_id(1,:).';                                      % NSUB x1 (event1)
r1_2 = sd_id(2).*(rho_id.*z_id(1,:).' + sqrt(1-rho_id.^2).*z_id(2,:).');
r2_1 = sd_trl(1).*z_trl(1,:).';                                    % NTRL x1
r2_2 = sd_trl(2).*(rho_trl.*z_trl(1,:).' + sqrt(1-rho_trl.^2).*z_trl(2,:).');

% paired event means: mu1 loads the event-1 effects, mu2 the event-2 effects
mu1 = b(1) + r1_1(m.id) + r2_1(m.trl);
mu2 = b(2) + r1_2(m.id) + r2_2(m.trl);
d1 = m.meas1 - mu1;
d2 = m.meas2 - mu2;

% bivariate residual with LSigma = diag(exp(b_sigma)) * Lrescor (lower
% triangular): a = sigma1, off = sigma2*rho_res, c = sigma2*sqrt(1-rho_res^2).
% multi_normal_cholesky_lpdf = -log(det LSigma) - 1/2||LSigma\(Y-Mu)||^2 (the
% -K/2*log(2*pi) constant is dropped). det LSigma = a*c is constant across rows.
sig1 = exp(b_sigma(1)); sig2 = exp(b_sigma(2));
a   = sig1;
off = sig2.*rho_res;
c   = sig2.*sqrt(1-rho_res.^2);
w1 = d1./a;
w2 = (d2 - (off./a).*d1)./c;
% log det LSigma = b_sigma1 + b_sigma2 + 1/2*log(1-rho_res^2); same for every row
logdet = b_sigma(1) + b_sigma(2) + 0.5*log(1-rho_res.^2);
ll = -m.NOBS*logdet - 0.5*sum(w1.^2 + w2.^2);

P = m.P;
lp = ll ...
    + sum(local_normal(b,0,P.b)) ...                          % b ~ normal(0,5)
    + sum(local_student_t(b_sigma,3,0,P.b_sigma)) ...         % b_sigma ~ t(3,0,10)
    + local_student_t(sd_id(1),3,0,P.sd_id)   + theta(ix.lsd_id(1)) ...  % half-t + Jac
    + local_student_t(sd_id(2),3,0,P.sd_id)   + theta(ix.lsd_id(2)) ...
    + local_student_t(sd_trl(1),3,0,P.sd_trl) + theta(ix.lsd_trl(1)) ...
    + local_student_t(sd_trl(2),3,0,P.sd_trl) + theta(ix.lsd_trl(2)) ...
    + 2*log(1 - rho_id.^2) ...                                % L_id  LKJ(2) + tanh Jac
    + 2*log(1 - rho_trl.^2) ...                               % L_trl LKJ(2) + tanh Jac
    + 2*log(1 - rho_res.^2) ...                               % Lrescor LKJ(2) + tanh Jac
    - 0.5*sum(z_id(:).^2) - 0.5*sum(z_trl(:).^2);             % std_normal raw
end

function out = local_diff_rescor_draws(draws,m)
%Transform draws to the CmdStan-shaped arrays case 7 (diffrescor=2) extracts:
%b, b_sigma (ND x2), sd_id, sd_trl (ND x2), id_varcov, trl_varcov (ND x2x2),
%rescor (ND x1), err_varcov (ND x2x2).
ix = m.ix;
out = struct();
out.b       = draws(:,ix.b);          % ND x2
out.b_sigma = draws(:,ix.b_sigma);    % ND x2
sd_id  = exp(draws(:,ix.lsd_id));     % ND x2
sd_trl = exp(draws(:,ix.lsd_trl));    % ND x2
out.sd_id  = sd_id;
out.sd_trl = sd_trl;
rho_id  = tanh(draws(:,ix.zrho_id));
rho_trl = tanh(draws(:,ix.zrho_trl));
rho_res = tanh(draws(:,ix.zrho_res)); % ND x1
% id_varcov / trl_varcov = quad_form_diag(corr, sd)
out.id_varcov  = local_varcov2(sd_id,rho_id);
out.trl_varcov = local_varcov2(sd_trl,rho_trl);
% rescor[1] = cor_res[1,2] = rho_res; err_varcov = quad_form_diag(cor_res, sigma)
out.rescor = rho_res;                 % ND x1 (matches Stan vector[1])
sigma = exp(draws(:,ix.b_sigma));     % ND x2 residual SDs
out.err_varcov = local_varcov2(sigma,rho_res);
end

% ===========================================================================
% Design: diff_dynrel (case 12) - group-level dynamic difference score
% ===========================================================================
function fit = local_hmc_diffdynrel(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the group-level dynamic/conditional
%two-event difference-score model, matching psyrat_build_stan_diffdynrel: the
%no-rescor difference model (local_hmc_diff) EXTENDED with event-specific fixed
%dimension slopes on BOTH the mean (Xdim*b_dim) and the per-event log-residual
%(Xdim*b_sigma_dim). Same two correlated 2-D person/trial random effects as case
%7; the slopes are the only addition (exactly as dynrel extends sserr).

m = struct();
m.NOBS  = double(data.NOBS);
m.NSUB  = double(data.NSUB);
m.NTRL  = double(data.NTRL);
m.meas  = double(data.meas(:));
m.meas1 = double(data.meas1(:));   % event-1 indicator (1/0)
m.meas2 = double(data.meas2(:));   % event-2 indicator (1/0)
m.id    = double(data.ID(:));      % sequential 1..NSUB
m.trl   = double(data.TRL(:));     % sequential 1..NTRL
% event-by-dimension design (KDIM = 2 for one dimension, 6 for two); flattened
% column-major as Xdim_vec = Xdim(:), so reshape reads column-major to match.
m.KDIM  = double(data.KDIM);
m.Xdim  = reshape(double(data.Xdim_vec(:)),m.NOBS,m.KDIM);

% configurable scale priors: the base difference block (b/b_sigma/sd_id/sd_trl)
% uses priors.diff; the dimension slopes use priors.dynrel (b for b_dim,
% b_sigma for b_sigma_dim), matching the Stan builder. Defaults mirror
% psyrat_default_priors so a default run is correct even without attached priors.
P = struct('b',5,'b_sigma',10,'sd_id',10,'sd_trl',10);
Pdim = struct('b',10,'b_sigma',1);
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors)
    pr = design.priors;
    if isfield(pr,'diff') && isstruct(pr.diff)
        fn = {'b','b_sigma','sd_id','sd_trl'};
        for k = 1:numel(fn)
            if isfield(pr.diff,fn{k}); P.(fn{k}) = double(pr.diff.(fn{k})); end
        end
    end
    if isfield(pr,'dynrel') && isstruct(pr.dynrel)
        if isfield(pr.dynrel,'b');       Pdim.b = double(pr.dynrel.b);             end
        if isfield(pr.dynrel,'b_sigma'); Pdim.b_sigma = double(pr.dynrel.b_sigma); end
    end
end
m.P = P;
m.Pdim = Pdim;

m = local_diffdynrel_index(m);

% gradient method (same AD/numer selection as the other builders)
gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
% Probe the Deep Learning Toolbox AD stack. NB: dlgradient is a method only
% callable inside dlfeval, so exist('dlgradient','file') is ALWAYS 0 even when
% the toolbox is installed; probe dlfeval (a real file, the AD entry point)
% instead, matching psyrat_estimation_engines_available's avail.ad.
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_diffdynrel_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_diffdynrel_logp(t,m), th);
end

start = local_diffdynrel_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_diffdynrel_draws(draws,m);
meta = struct('design','diff_dynrel','gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_diffdynrel_index(m)
%theta = [b(2); b_sigma(2); b_dim(KDIM); b_sigma_dim(KDIM); log_sd_id(2);
%         z_rho_id; log_sd_trl(2); z_rho_trl; z_id(2*NSUB); z_trl(2*NTRL)].
%Identical to local_diff_index plus the event-by-dimension slope blocks.
NSUB = m.NSUB; NTRL = m.NTRL; K = m.KDIM; i = 0;
m.ix.b         = i+(1:2);      i = i+2;
m.ix.b_sigma   = i+(1:2);      i = i+2;
m.ix.b_dim     = i+(1:K);      i = i+K;   % event-specific mean dim slopes
m.ix.bs_dim    = i+(1:K);      i = i+K;   % event-specific log-resid dim slopes
m.ix.lsd_id    = i+(1:2);      i = i+2;   % log(sd_id)
m.ix.zrho_id   = i+1;          i = i+1;   % atanh(rho_id)
m.ix.lsd_trl   = i+(1:2);      i = i+2;   % log(sd_trl)
m.ix.zrho_trl  = i+1;          i = i+1;   % atanh(rho_trl)
m.ix.z_id      = i+(1:2*NSUB); i = i+2*NSUB;
m.ix.z_trl     = i+(1:2*NTRL); i = i+2*NTRL;
m.D = i;
m.mon_idx = [m.ix.b, m.ix.b_sigma, m.ix.b_dim, m.ix.bs_dim, m.ix.lsd_id, ...
    m.ix.zrho_id, m.ix.lsd_trl, m.ix.zrho_trl];
m.mon_names = [{'b1','b2','b_sigma1','b_sigma2'}, ...
    local_numbered('b_dim',K), local_numbered('b_sigma_dim',K), ...
    {'log_sd_id1','log_sd_id2','z_rho_id','log_sd_trl1','log_sd_trl2','z_rho_trl'}];
end

function start = local_diffdynrel_start(m)
%Per-event means at the event-conditional data means; log-SDs at rough guesses;
%dimension slopes start at 0 (same as local_diff_start plus the zero slopes).
start = zeros(m.D,1);
m1 = m.meas1 > 0; m2 = m.meas2 > 0;
start(m.ix.b(1)) = mean(m.meas(m1));
start(m.ix.b(2)) = mean(m.meas(m2));
s1 = std(m.meas(m1)); s2 = std(m.meas(m2));
start(m.ix.b_sigma(1)) = log(max(s1/2,1e-3));
start(m.ix.b_sigma(2)) = log(max(s2/2,1e-3));
start(m.ix.lsd_id)  = log(max((s1+s2)/4,1e-3));
start(m.ix.lsd_trl) = log(max((s1+s2)/4,1e-3));
% b_dim, b_sigma_dim, z_rho_*, z_id, z_trl start at 0
end

function lp = local_diffdynrel_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching
%psyrat_build_stan_diffdynrel: local_diff_logp with the dimension-conditioned
%terms Xdim*b_dim on the mean and Xdim*b_sigma_dim on the log-residual, plus the
%normal slope priors.
ix = m.ix;
b       = theta(ix.b);              % 2 (per-event mean)
b_sigma = theta(ix.b_sigma);        % 2 (per-event log-residual SD)
b_dim   = theta(ix.b_dim);          % KDIM (mean dim slopes)
bs_dim  = theta(ix.bs_dim);         % KDIM (log-resid dim slopes)
sd_id   = exp(theta(ix.lsd_id));    % 2 (person event-SDs), >0
rho_id  = tanh(theta(ix.zrho_id));
sd_trl  = exp(theta(ix.lsd_trl));   % 2 (trial event-SDs), >0
rho_trl = tanh(theta(ix.zrho_trl));
z_id  = reshape(theta(ix.z_id),2,m.NSUB);
z_trl = reshape(theta(ix.z_trl),2,m.NTRL);

% non-centered person/trial event-effects via 2-D Cholesky diag(sd)*[1 0; rho .]
r1_1 = sd_id(1).*z_id(1,:).';
r1_2 = sd_id(2).*(rho_id.*z_id(1,:).' + sqrt(1-rho_id.^2).*z_id(2,:).');
r2_1 = sd_trl(1).*z_trl(1,:).';
r2_2 = sd_trl(2).*(rho_trl.*z_trl(1,:).' + sqrt(1-rho_trl.^2).*z_trl(2,:).');

m1 = m.meas1; m2 = m.meas2;
% event-indicator design + the event-by-dimension fixed slopes on both submodels
mu = m1.*b(1) + m2.*b(2) + m.Xdim*b_dim ...
   + r1_1(m.id).*m1 + r1_2(m.id).*m2 ...
   + r2_1(m.trl).*m1 + r2_2(m.trl).*m2;
logsig = m1.*b_sigma(1) + m2.*b_sigma(2) + m.Xdim*bs_dim;
sig = exp(logsig);
ll = -sum(logsig) - sum((m.meas - mu).^2 ./ (2*sig.^2));

P = m.P; Pdim = m.Pdim;
lp = ll ...
    + sum(local_normal(b,0,P.b)) ...                          % b ~ normal(0,5)
    + sum(local_student_t(b_sigma,3,0,P.b_sigma)) ...         % b_sigma ~ t(3,0,10)
    + sum(local_normal(b_dim,0,Pdim.b)) ...                   % b_dim ~ normal(0,.)
    + sum(local_normal(bs_dim,0,Pdim.b_sigma)) ...            % b_sigma_dim ~ normal(0,.)
    + local_student_t(sd_id(1),3,0,P.sd_id)   + theta(ix.lsd_id(1)) ...  % half-t + Jac
    + local_student_t(sd_id(2),3,0,P.sd_id)   + theta(ix.lsd_id(2)) ...
    + local_student_t(sd_trl(1),3,0,P.sd_trl) + theta(ix.lsd_trl(1)) ...
    + local_student_t(sd_trl(2),3,0,P.sd_trl) + theta(ix.lsd_trl(2)) ...
    + 2*log(1 - rho_id.^2) ...                                % L_id  LKJ(2) + tanh Jac
    + 2*log(1 - rho_trl.^2) ...                               % L_trl LKJ(2) + tanh Jac
    - 0.5*sum(z_id(:).^2) - 0.5*sum(z_trl(:).^2);             % std_normal raw
end

function out = local_diffdynrel_draws(draws,m)
%Transform draws to the CmdStan-shaped arrays case 12 extracts: b, b_sigma (ND x2),
%b_dim, b_sigma_dim (ND x KDIM), sd_id, sd_trl (ND x2), id_varcov, trl_varcov
%(ND x2x2).
ix = m.ix;
out = struct();
out.b           = draws(:,ix.b);          % ND x2
out.b_sigma     = draws(:,ix.b_sigma);    % ND x2
out.b_dim       = draws(:,ix.b_dim);      % ND x KDIM
out.b_sigma_dim = draws(:,ix.bs_dim);     % ND x KDIM
sd_id  = exp(draws(:,ix.lsd_id));         % ND x2
sd_trl = exp(draws(:,ix.lsd_trl));        % ND x2
out.sd_id  = sd_id;
out.sd_trl = sd_trl;
rho_id  = tanh(draws(:,ix.zrho_id));
rho_trl = tanh(draws(:,ix.zrho_trl));
out.id_varcov  = local_varcov2(sd_id,rho_id);
out.trl_varcov = local_varcov2(sd_trl,rho_trl);
end

% ===========================================================================
% Design: diff_dynrel_sserr (case 13) - subject-level dynamic difference score
% ===========================================================================
function fit = local_hmc_diffdynrel_sserr(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the SUBJECT-LEVEL dynamic difference-score
%model, matching psyrat_build_stan_diffdynrel_sserr: the group-level dynamic
%difference (local_hmc_diffdynrel) with the person random-effect block WIDENED
%from 2-D to a 4-D joint correlated block over {location event1, location event2,
%scale event1, scale event2}. The scale coefficients (3,4) load on the per-event
%log-residual so each participant has their own residual; the 4x4 person
%correlation uses the general LKJ Cholesky transform (psyrat_corr_chol_lkj). The
%trial block stays 2-D (location-only), as in case 12.

m = struct();
m.NOBS  = double(data.NOBS);
m.NSUB  = double(data.NSUB);
m.NTRL  = double(data.NTRL);
m.meas  = double(data.meas(:));
m.meas1 = double(data.meas1(:));   % event-1 indicator (1/0)
m.meas2 = double(data.meas2(:));   % event-2 indicator (1/0)
m.id    = double(data.ID(:));      % sequential 1..NSUB
m.trl   = double(data.TRL(:));     % sequential 1..NTRL
m.KDIM  = double(data.KDIM);
m.Xdim  = reshape(double(data.Xdim_vec(:)),m.NOBS,m.KDIM);

% configurable scale priors: base difference block (b/b_sigma/sd_trl, the person
% LOCATION SDs sd_id[1:2]) from priors.diff; the person SCALE SDs sd_id[3:4] use
% priors.diff.sd_id_sigma; the dimension slopes use priors.dynrel. Defaults mirror
% psyrat_default_priors.
P = struct('b',5,'b_sigma',10,'sd_id',10,'sd_trl',10,'sd_id_sigma',2.5);
Pdim = struct('b',10,'b_sigma',1);
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors)
    pr = design.priors;
    if isfield(pr,'diff') && isstruct(pr.diff)
        fn = {'b','b_sigma','sd_id','sd_trl','sd_id_sigma'};
        for k = 1:numel(fn)
            if isfield(pr.diff,fn{k}); P.(fn{k}) = double(pr.diff.(fn{k})); end
        end
    end
    if isfield(pr,'dynrel') && isstruct(pr.dynrel)
        if isfield(pr.dynrel,'b');       Pdim.b = double(pr.dynrel.b);             end
        if isfield(pr.dynrel,'b_sigma'); Pdim.b_sigma = double(pr.dynrel.b_sigma); end
    end
end
m.P = P;
m.Pdim = Pdim;

m = local_diffdynrel_sserr_index(m);

gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_diffdynrel_sserr_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_diffdynrel_sserr_logp(t,m), th);
end

start = local_diffdynrel_sserr_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_diffdynrel_sserr_draws(draws,m);
meta = struct('design','diff_dynrel_sserr','gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_diffdynrel_sserr_index(m)
%theta = [b(2); b_sigma(2); b_dim(K); b_sigma_dim(K);
%         log_sd_id(4); yLid(6);          % 4-D joint person block (4 SDs, 6 CPCs)
%         log_sd_trl(2); z_rho_trl;        % 2-D trial block
%         z_id(4*NSUB); z_trl(2*NTRL)].
%sd_id[1:2] are the person LOCATION SDs, sd_id[3:4] the person SCALE (log-residual)
%SDs; yLid are the 4*3/2 = 6 unconstrained CPCs of the 4x4 person correlation.
NSUB = m.NSUB; NTRL = m.NTRL; K = m.KDIM; i = 0;
m.ix.b         = i+(1:2);      i = i+2;
m.ix.b_sigma   = i+(1:2);      i = i+2;
m.ix.b_dim     = i+(1:K);      i = i+K;
m.ix.bs_dim    = i+(1:K);      i = i+K;
m.ix.lsd_id    = i+(1:4);      i = i+4;   % log(sd_id), 4 entries
m.ix.yLid      = i+(1:6);      i = i+6;   % 4x4 LKJ CPCs (row-major lower tri)
m.ix.lsd_trl   = i+(1:2);      i = i+2;
m.ix.zrho_trl  = i+1;          i = i+1;
m.ix.z_id      = i+(1:4*NSUB); i = i+4*NSUB;
m.ix.z_trl     = i+(1:2*NTRL); i = i+2*NTRL;
m.D = i;
m.mon_idx = [m.ix.b, m.ix.b_sigma, m.ix.b_dim, m.ix.bs_dim, m.ix.lsd_id, ...
    m.ix.yLid, m.ix.lsd_trl, m.ix.zrho_trl];
m.mon_names = [{'b1','b2','b_sigma1','b_sigma2'}, ...
    local_numbered('b_dim',K), local_numbered('b_sigma_dim',K), ...
    {'log_sd_id1','log_sd_id2','log_sd_id3','log_sd_id4'}, ...
    local_numbered('yLid',6), {'log_sd_trl1','log_sd_trl2','z_rho_trl'}];
end

function start = local_diffdynrel_sserr_start(m)
%Per-event means at the event-conditional data means; person LOCATION SDs at a
%fraction of the marginal SD, person SCALE SDs small; CPCs and raw effects at 0.
start = zeros(m.D,1);
m1 = m.meas1 > 0; m2 = m.meas2 > 0;
start(m.ix.b(1)) = mean(m.meas(m1));
start(m.ix.b(2)) = mean(m.meas(m2));
s1 = std(m.meas(m1)); s2 = std(m.meas(m2));
start(m.ix.b_sigma(1)) = log(max(s1/2,1e-3));
start(m.ix.b_sigma(2)) = log(max(s2/2,1e-3));
start(m.ix.lsd_id(1:2)) = log(max((s1+s2)/4,1e-3));  % person LOCATION SDs
start(m.ix.lsd_id(3:4)) = log(0.3);                  % person SCALE (log-resid) SDs
start(m.ix.lsd_trl)     = log(max((s1+s2)/4,1e-3));
% yLid (-> identity correlation), z_rho_trl, b_dim, b_sigma_dim, z_id, z_trl: 0
end

function lp = local_diffdynrel_sserr_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching
%psyrat_build_stan_diffdynrel_sserr: local_diffdynrel_logp with the person block
%widened to a 4-D joint correlated block (location + log-residual per event), the
%scale coefficients (3,4) loaded on the per-event log-residual, and the 4x4 LKJ
%person correlation via psyrat_corr_chol_lkj.
ix = m.ix;
b       = theta(ix.b);
b_sigma = theta(ix.b_sigma);
b_dim   = theta(ix.b_dim);
bs_dim  = theta(ix.bs_dim);
sd_id   = exp(theta(ix.lsd_id));    % 4 (1-2 location, 3-4 scale), >0
[L_id, lp_Lid] = psyrat_corr_chol_lkj(theta(ix.yLid), 4, 2);  % 4x4 + LKJ(2) logdens
sd_trl  = exp(theta(ix.lsd_trl));   % 2, >0
rho_trl = tanh(theta(ix.zrho_trl));
z_id  = reshape(theta(ix.z_id),4,m.NSUB);
z_trl = reshape(theta(ix.z_trl),2,m.NTRL);

% person 4-D effects via the non-centered joint Cholesky:
% r_1 = (diag(sd_id) * L_id * z_id)'  -> NSUB x 4. Columns: 1-2 location event1/2,
% 3-4 scale (log-residual) event1/2.
r1 = (sd_id .* (L_id * z_id)).';      % NSUB x 4
r1_1  = r1(:,1); r1_2  = r1(:,2);     % location REs
r1_s1 = r1(:,3); r1_s2 = r1(:,4);     % scale (log-residual) REs
% trial 2-D location effects via the 2x2 Cholesky diag(sd_trl)*[1 0; rho .]
r2_1 = sd_trl(1).*z_trl(1,:).';
r2_2 = sd_trl(2).*(rho_trl.*z_trl(1,:).' + sqrt(1-rho_trl.^2).*z_trl(2,:).');

m1 = m.meas1; m2 = m.meas2;
mu = m1.*b(1) + m2.*b(2) + m.Xdim*b_dim ...
   + r1_1(m.id).*m1 + r1_2(m.id).*m2 ...
   + r2_1(m.trl).*m1 + r2_2(m.trl).*m2;
% per-event log-residual: population + dim slopes + per-subject scale REs
logsig = m1.*b_sigma(1) + m2.*b_sigma(2) + m.Xdim*bs_dim ...
   + r1_s1(m.id).*m1 + r1_s2(m.id).*m2;
sig = exp(logsig);
ll = -sum(logsig) - sum((m.meas - mu).^2 ./ (2*sig.^2));

P = m.P; Pdim = m.Pdim;
lp = ll + lp_Lid ...                                          % L_id LKJ(2) + Jac
    + sum(local_normal(b,0,P.b)) ...                          % b ~ normal(0,5)
    + sum(local_student_t(b_sigma,3,0,P.b_sigma)) ...         % b_sigma ~ t(3,0,10)
    + sum(local_normal(b_dim,0,Pdim.b)) ...                   % b_dim
    + sum(local_normal(bs_dim,0,Pdim.b_sigma)) ...            % b_sigma_dim
    + local_student_t(sd_id(1),3,0,P.sd_id)       + theta(ix.lsd_id(1)) ...  % location SDs
    + local_student_t(sd_id(2),3,0,P.sd_id)       + theta(ix.lsd_id(2)) ...
    + local_student_t(sd_id(3),3,0,P.sd_id_sigma) + theta(ix.lsd_id(3)) ...  % scale SDs
    + local_student_t(sd_id(4),3,0,P.sd_id_sigma) + theta(ix.lsd_id(4)) ...
    + local_student_t(sd_trl(1),3,0,P.sd_trl)     + theta(ix.lsd_trl(1)) ...
    + local_student_t(sd_trl(2),3,0,P.sd_trl)     + theta(ix.lsd_trl(2)) ...
    + 2*log(1 - rho_trl.^2) ...                               % L_trl LKJ(2) + tanh Jac
    - 0.5*sum(z_id(:).^2) - 0.5*sum(z_trl(:).^2);             % std_normal raw
end

function out = local_diffdynrel_sserr_draws(draws,m)
%Transform draws to the CmdStan-shaped arrays case 13 extracts: b, b_sigma (ND x2),
%b_dim, b_sigma_dim (ND x KDIM), sd_id (ND x4), sd_trl (ND x2), L_id (ND x4x4),
%id_varcov, trl_varcov (ND x2x2), er_var_ss (ND x NSUB x 2). The per-draw 4x4 L_id
%and the per-subject scale REs are reconstructed from the raw draws.
ix = m.ix;
ND = size(draws,1); NSUB = m.NSUB;
out = struct();
out.b           = draws(:,ix.b);
out.b_sigma     = draws(:,ix.b_sigma);
out.b_dim       = draws(:,ix.b_dim);
out.b_sigma_dim = draws(:,ix.bs_dim);
sd_id  = exp(draws(:,ix.lsd_id));     % ND x4
sd_trl = exp(draws(:,ix.lsd_trl));    % ND x2
out.sd_id  = sd_id;
out.sd_trl = sd_trl;
rho_trl = tanh(draws(:,ix.zrho_trl));
out.trl_varcov = local_varcov2(sd_trl,rho_trl);

L_id = zeros(ND,4,4);
id_varcov = zeros(ND,2,2);
er = zeros(ND,NSUB,2);
for d = 1:ND
    Ld = psyrat_corr_chol_lkj(draws(d,ix.yLid).',4,2);   % 4x4 correlation Cholesky
    L_id(d,:,:) = Ld;
    cor = Ld*Ld.';                       % 4x4 correlation matrix
    sdd = sd_id(d,:);                    % 1x4
    % id_varcov = quad_form_diag(cor,sd_id)[1:2,1:2] (location-location block)
    id_varcov(d,1,1) = sdd(1)^2;
    id_varcov(d,2,2) = sdd(2)^2;
    id_varcov(d,1,2) = sdd(1)*sdd(2)*cor(1,2);
    id_varcov(d,2,1) = id_varcov(d,1,2);
    % per-subject scale REs from r_1 = (diag(sd_id)*L_id*z_id)' columns 3,4
    z4 = reshape(draws(d,ix.z_id),4,NSUB);
    r1 = (sdd(:) .* (Ld * z4)).';        % NSUB x 4
    er(d,:,1) = draws(d,ix.b_sigma(1)) + r1(:,3).';   % er_var_ss event1
    er(d,:,2) = draws(d,ix.b_sigma(2)) + r1(:,4).';   % er_var_ss event2
end
out.L_id = L_id;
out.id_varcov = id_varcov;
out.er_var_ss = er;
end

% ===========================================================================
% Design: dynrel_trt (case 14) - two-facet (trial + occasion) dynamic
% reliability, non-difference
% ===========================================================================
function fit = local_hmc_dynrel_trt(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the two-facet dynamic reliability model,
%matching psyrat_build_stan_dynrel_trt: the one-facet dynrel locscale model
%(local_hmc_locscale: 2-D person location-scale Cholesky + Xdim mean/scale
%dimension slopes) with the single trial main effect GENERALIZED to FIVE scalar
%crossed-facet mean random effects (occ, trl, trlxid, occxid, trlxocc), each
%cauchy-scaled. The residual scale carries only the person ind_sd (no facet scale
%REs), so sigma stays a single per-person sigma_pto,e^2(z).

m = struct();
m.NOBS = double(data.NOBS);
m.NSUB = double(data.NSUB);
m.meas = double(data.meas(:));
m.id   = double(data.id(:));      % sequential 1..NSUB
m.KDIM = double(data.KDIM);
m.Xdim = reshape(double(data.Xdim_vec(:)),m.NOBS,m.KDIM);
% the five crossed mean facets: per-observation index + level count + the
% priors.trt prior scale (cauchy). Order is fixed (it sets the parameter order and
% the draw-field mapping below).
m.facets = struct( ...
    'name',  {'occ','trl','trlxid','occxid','trlxocc'}, ...
    'idx',   {double(data.occ(:)),double(data.trl(:)),double(data.trlxid(:)), ...
              double(data.occxid(:)),double(data.trlxocc(:))}, ...
    'N',     {double(data.NOCC),double(data.NTRL),double(data.NTID), ...
              double(data.NOID),double(data.NTO)}, ...
    'prior', {5,10,5,5,1});   % defaults: priors.trt.{sig_occ,sig_trl,sig_trlxid,sig_occxid,sig_trlxocc}

% configurable priors: dimension slopes from priors.dynrel (b, b_sigma); the five
% facet cauchy scales from priors.trt (sig_occ/sig_trl/sig_trlxid/sig_occxid/
% sig_trlxocc). The person block + intercepts use the fixed case-11 student_t
% constants baked into the log-posterior.
Pdim = struct('b',10,'b_sigma',1);
trtmap = {'sig_occ','sig_trl','sig_trlxid','sig_occxid','sig_trlxocc'};
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors)
    pr = design.priors;
    if isfield(pr,'dynrel') && isstruct(pr.dynrel)
        if isfield(pr.dynrel,'b');       Pdim.b = double(pr.dynrel.b);             end
        if isfield(pr.dynrel,'b_sigma'); Pdim.b_sigma = double(pr.dynrel.b_sigma); end
    end
    if isfield(pr,'trt') && isstruct(pr.trt)
        for k = 1:5
            if isfield(pr.trt,trtmap{k}); m.facets(k).prior = double(pr.trt.(trtmap{k})); end
        end
    end
end
m.Pdim = Pdim;

m = local_dynrel_trt_index(m);

gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_dynrel_trt_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_dynrel_trt_logp(t,m), th);
end

start = local_dynrel_trt_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_dynrel_trt_draws(draws,m);
meta = struct('design','dynrel_trt','gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_dynrel_trt_index(m)
%theta = [Intercept; Intercept_sigma; b(K); b_sigma(K); log_gro(2); z_rho;
%         log_sig(5);                       % the five facet log-SDs (facet order)
%         gro_effs(2*NSUB);                 % person raw z's (reshape 2 x NSUB)
%         occ_raw; trl_raw; trlxid_raw; occxid_raw; trlxocc_raw].
K = m.KDIM; NSUB = m.NSUB; i = 0;
m.ix.Intercept = i+1;          i = i+1;
m.ix.Isig      = i+1;          i = i+1;
m.ix.b         = i+(1:K);      i = i+K;
m.ix.b_sigma   = i+(1:K);      i = i+K;
m.ix.lgro      = i+(1:2);      i = i+2;
m.ix.zrho      = i+1;          i = i+1;
m.ix.lsig      = i+(1:5);      i = i+5;   % facet log-SDs
m.ix.zsub      = i+(1:2*NSUB); i = i+2*NSUB;
m.ix.fraw = cell(1,5);
for k = 1:5
    m.ix.fraw{k} = i+(1:m.facets(k).N); i = i+m.facets(k).N;
end
m.D = i;
m.mon_idx = [m.ix.Intercept, m.ix.Isig, m.ix.b, m.ix.b_sigma, m.ix.lgro, ...
    m.ix.zrho, m.ix.lsig];
m.mon_names = [{'Intercept','Intercept_sigma'}, ...
    local_numbered('b',K), local_numbered('b_sigma',K), ...
    {'log_sigma_p','log_sigma_dp','z_rho', ...
     'log_sig_occ','log_sig_trl','log_sig_trlxid','log_sig_occxid','log_sig_trlxocc'}];
end

function start = local_dynrel_trt_start(m)
%Intercepts at the data scale, log-scales at rough SD guesses, raw effects + rho
%at 0 (the same shape as the one-facet dynrel start, plus the five facet SDs).
start = zeros(m.D,1);
smeas = std(m.meas);
start(m.ix.Intercept) = mean(m.meas);
start(m.ix.Isig)      = log(max(smeas/2,1e-3));
start(m.ix.lgro(1))   = log(max(smeas/2,1e-3));
start(m.ix.lgro(2))   = log(0.3);
start(m.ix.lsig)      = log(max(smeas/3,1e-3));   % five facet SDs
% b, b_sigma, z_rho, gro_effs, facet raws start at 0
end

function lp = local_dynrel_trt_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching
%psyrat_build_stan_dynrel_trt: the one-facet dynrel locscale likelihood with five
%scalar crossed-facet mean random effects added to mu.
ix = m.ix;
Intercept = theta(ix.Intercept);
Isig      = theta(ix.Isig);
b         = theta(ix.b);
b_sigma   = theta(ix.b_sigma);
gro       = exp(theta(ix.lgro));      % (sigma_p, sigma_dp), >0
rho       = tanh(theta(ix.zrho));
sigf      = exp(theta(ix.lsig));      % 5 facet SDs, >0
zsub      = reshape(theta(ix.zsub),2,m.NSUB);
z1 = zsub(1,:).'; z2 = zsub(2,:).';

% person location + log-residual effects via the 2-D Cholesky (same as case 11)
ind_bs = gro(1).*z1;
ind_sd = gro(2).*(rho.*z1 + sqrt(1-rho.^2).*z2);

% mean: intercept + dimension slopes + person location + the five facet mains
mu = Intercept + m.Xdim*b + ind_bs(m.id);
for k = 1:5
    fk = m.facets(k);
    term = sigf(k).*theta(ix.fraw{k});     % non-centered facet term = sig * raw
    mu = mu + term(fk.idx);
end
% log-residual: intercept + dimension slopes + person scale (no facet scale REs)
logsig = Isig + m.Xdim*b_sigma + ind_sd(m.id);
sig = exp(logsig);
ll = -sum(logsig) - sum((m.meas - mu).^2 ./ (2*sig.^2));

Pdim = m.Pdim;
lp = ll ...
    + local_student_t(Intercept,3,15.7,12.4) ...        % Intercept
    + local_student_t(Isig,3,0,2.5) ...                 % Intercept_sigma
    + local_student_t(gro(1),3,0,12.4) + theta(ix.lgro(1)) ... % half-t + Jac
    + local_student_t(gro(2),3,0,12.4) + theta(ix.lgro(2)) ...
    + log(1 - rho.^2) ...                               % LKJ(1) uniform + tanh Jac
    + sum(local_normal(b,0,Pdim.b)) ...                 % mean dimension slopes
    + sum(local_normal(b_sigma,0,Pdim.b_sigma)) ...     % scale dimension slopes
    - 0.5*sum(z1.^2) - 0.5*sum(z2.^2);                  % std_normal person raw
% the five facet half-cauchy SDs (+ log Jac) and their std_normal raw effects
for k = 1:5
    lp = lp + local_cauchy(sigf(k),0,m.facets(k).prior) + theta(ix.lsig(k)) ...
            - 0.5*sum(theta(ix.fraw{k}).^2);
end
end

function out = local_dynrel_trt_draws(draws,m)
%Transform draws to the CmdStan-shaped arrays case 14 extracts: Intercept,
%Intercept_sigma (ND x1), gro_sds (ND x2), b, b_sigma (ND x KDIM), chol_corrmat
%(ND x2x2), and the five facet SDs sig_occ/sig_trl/sig_trlxid/sig_occxid/
%sig_trlxocc (ND x1 each).
ix = m.ix;
ND = size(draws,1);
out = struct();
out.Intercept       = draws(:,ix.Intercept);
out.Intercept_sigma = draws(:,ix.Isig);
out.gro_sds         = exp(draws(:,ix.lgro));      % ND x2
out.b               = draws(:,ix.b);
out.b_sigma         = draws(:,ix.b_sigma);
rho = tanh(draws(:,ix.zrho));
chol = zeros(ND,2,2);
chol(:,1,1) = 1; chol(:,2,1) = rho; chol(:,2,2) = sqrt(1 - rho.^2);
out.chol_corrmat = chol;
sigf = exp(draws(:,ix.lsig));                     % ND x5 (facet order)
out.sig_occ     = sigf(:,1);
out.sig_trl     = sigf(:,2);
out.sig_trlxid  = sigf(:,3);
out.sig_occxid  = sigf(:,4);
out.sig_trlxocc = sigf(:,5);
end

% ===========================================================================
% Design: dod (case 9) - task-free four-event difference-of-differences
% ===========================================================================
function fit = local_hmc_dodiff(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the four-event difference-of-differences
%model, matching psyrat_build_stan_dodiff: a 4-cell location-scale model with
%per-cell mean (b_cell) and log-residual (b_sigma_cell), and person AND trial
%random effects on BOTH the location (cells 1-4) and the scale (cells 1-4) - i.e.
%an 8-D joint correlated block per grouping factor via an 8x8 LKJ Cholesky. The
%DoD contrast is applied downstream from the per-cell var-covariances, so the
%engine only fits the 4-cell model. Prior SCALES come from priors.dod (b_cell ~
%normal(0,b_cell); b_sigma_cell ~ student_t(10,0,b_sigma_cell); sd_id/sd_trial
%~ student_t(10,0,sd_id/sd_trial); L ~ LKJ(lkj)); the defaults reproduce the
%brms constants this fit hard-coded before 2026-09-05, and the student-t
%degrees of freedom (10) stay baked in to match the Stan builder.

m = struct();
m.N      = double(data.N);
m.C      = double(data.C);             % = 4
m.NB     = 2*m.C;                       % 8-D block dimension (location + scale)
m.Y      = double(data.Y(:));
m.cell   = double(data.cell(:));       % 1..C per observation
m.N_id   = double(data.N_id);
m.id     = double(data.id(:));         % sequential 1..N_id
m.N_trl  = double(data.N_trial);
m.trl    = double(data.trial(:));      % sequential 1..N_trial
% precompute the column-major linear indices for the per-observation gather into
% the (N_id x 8) / (N_trl x 8) random-effect matrices: location uses cell, scale
% uses cell + C.
m.lin_id_loc  = m.id  + (m.cell-1).*m.N_id;
m.lin_id_sca  = m.id  + (m.cell+m.C-1).*m.N_id;
m.lin_trl_loc = m.trl + (m.cell-1).*m.N_trl;
m.lin_trl_sca = m.trl + (m.cell+m.C-1).*m.N_trl;

m = local_dod_index(m);

% configurable prior scales (priors.dod, exposed 2026-09-05). Defaults mirror
% psyrat_default_priors.dod so a default run is correct even if the caller did
% not attach priors; the caller (psyrat_computevarcomp case 9) attaches
% priors.dod via native_extra. The student-t degrees of freedom (10) are
% baked into local_dod_logp to match the Stan builder.
P = struct('b_cell',10,'b_sigma_cell',2,'sd_id',2,'sd_trial',2,'lkj',2);
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors) && ...
        isfield(design.priors,'dod') && isstruct(design.priors.dod)
    dp = design.priors.dod;
    flds = fieldnames(P);
    for k = 1:numel(flds)
        if isfield(dp,flds{k}); P.(flds{k}) = double(dp.(flds{k})); end
    end
end
m.pri = P;

gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_dod_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_dod_logp(t,m), th);
end

start = local_dod_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_dod_draws(draws,m);
meta = struct('design','dod','gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_dod_index(m)
%theta = [b_cell(C); b_sigma_cell(C);
%         log_sd_id(2C); yLid(2C*(2C-1)/2);      % 8-D person block (8 SDs + 28 CPCs)
%         log_sd_trl(2C); yLtrl(2C*(2C-1)/2);    % 8-D trial block
%         z_id(2C*N_id); z_trl(2C*N_trl)].
C = m.C; NB = m.NB; NC = NB*(NB-1)/2;   % NC = 28 CPCs for the 8x8 LKJ
i = 0;
m.ix.b_cell   = i+(1:C);  i = i+C;
m.ix.bs_cell  = i+(1:C);  i = i+C;
m.ix.lsd_id   = i+(1:NB); i = i+NB;
m.ix.yLid     = i+(1:NC); i = i+NC;
m.ix.lsd_trl  = i+(1:NB); i = i+NB;
m.ix.yLtrl    = i+(1:NC); i = i+NC;
m.ix.z_id     = i+(1:NB*m.N_id);  i = i+NB*m.N_id;
m.ix.z_trl    = i+(1:NB*m.N_trl); i = i+NB*m.N_trl;
m.D = i;
m.mon_idx = [m.ix.b_cell, m.ix.bs_cell, m.ix.lsd_id, m.ix.lsd_trl];
m.mon_names = [local_numbered('b_cell',C), local_numbered('b_sigma_cell',C), ...
    local_numbered('log_sd_id',NB), local_numbered('log_sd_trl',NB)];
end

function start = local_dod_start(m)
%Per-cell means/log-SDs at the cell-conditional data moments; log-SDs at rough
%guesses; CPCs and raw effects at 0.
start = zeros(m.D,1);
for c = 1:m.C
    yc = m.Y(m.cell == c);
    if isempty(yc); yc = m.Y; end
    start(m.ix.b_cell(c))  = mean(yc);
    start(m.ix.bs_cell(c)) = log(max(std(yc)/2,1e-3));
end
smeas = std(m.Y);
start(m.ix.lsd_id)  = log(max(smeas/3,1e-3));   % location + scale block SDs
start(m.ix.lsd_trl) = log(max(smeas/3,1e-3));
% yLid, yLtrl (-> identity correlation), z_id, z_trl start at 0
end

function lp = local_dod_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching psyrat_build_stan_dodiff.
ix = m.ix; C = m.C;
b_cell  = theta(ix.b_cell);
bs_cell = theta(ix.bs_cell);
sd_id   = exp(theta(ix.lsd_id));    % 2C, >0
sd_trl  = exp(theta(ix.lsd_trl));   % 2C, >0
[L_id, lp_Lid] = psyrat_corr_chol_lkj(theta(ix.yLid), m.NB, m.pri.lkj);   % 8x8 + LKJ(lkj)
[L_trl,lp_Ltrl]= psyrat_corr_chol_lkj(theta(ix.yLtrl),m.NB, m.pri.lkj);
z_id  = reshape(theta(ix.z_id), m.NB, m.N_id);
z_trl = reshape(theta(ix.z_trl),m.NB, m.N_trl);

% non-centered 8-D effects: r = (diag(sd) * L * z)'  -> (N_id|N_trl) x 8
r_id  = (sd_id  .* (L_id  * z_id )).';   % N_id  x 8
r_trl = (sd_trl .* (L_trl * z_trl)).';   % N_trl x 8

% per-observation gather: location uses cell (cols 1..C), scale uses cell+C
mu = b_cell(m.cell) + r_id(m.lin_id_loc) + r_trl(m.lin_trl_loc);
logsig = bs_cell(m.cell) + r_id(m.lin_id_sca) + r_trl(m.lin_trl_sca);
sig = exp(logsig);
ll = -sum(logsig) - sum((m.Y - mu).^2 ./ (2*sig.^2));

% priors.dod scales (m.pri, defaults = the former brms constants): note the
% student_t nu = 10 here, not 3; it is structural and stays fixed
lp = ll + lp_Lid + lp_Ltrl ...
    + sum(local_normal(b_cell,0,m.pri.b_cell)) ...                 % b_cell ~ normal(0,b_cell)
    + sum(local_student_t(bs_cell,10,0,m.pri.b_sigma_cell)) ...    % b_sigma_cell ~ t(10,0,.)
    + sum(local_student_t(sd_id,10,0,m.pri.sd_id))  + sum(theta(ix.lsd_id)) ...  % half-t + Jac (8)
    + sum(local_student_t(sd_trl,10,0,m.pri.sd_trial)) + sum(theta(ix.lsd_trl)) ...
    - 0.5*sum(z_id(:).^2) - 0.5*sum(z_trl(:).^2);              % std_normal raw
end

function out = local_dod_draws(draws,m)
%Transform draws to the CmdStan-shaped arrays case 9 extracts: b_sigma_cell
%(ND x C), err_var_Cell (ND x C = exp(2*b_sigma_cell)), id_varcov / trl_varcov
%(ND x C x C, the location-cell var-covariances), and cor_ID_Cell /
%cor_TrialNumber_Cell (ND x C x C, the location-cell correlation blocks).
ix = m.ix; C = m.C; NB = m.NB;
ND = size(draws,1);
out = struct();
%b_cell must be packed too: case 9's extraction gained it on 2026-07-27 (the
%DoD -> single-event criterion view needs the universe-score mean), and the
%engine contract is that every parameter the case extracts is populated. The
%draws already carry it (ix.b_cell is monitored); it was only never surfaced.
out.b_cell = draws(:,ix.b_cell);      % ND x C
bs = draws(:,ix.bs_cell);             % ND x C
out.b_sigma_cell = bs;
out.err_var_Cell = exp(2*bs);         % ND x C
sd_id  = exp(draws(:,ix.lsd_id));     % ND x 2C
sd_trl = exp(draws(:,ix.lsd_trl));    % ND x 2C

id_varcov  = zeros(ND,C,C);
trl_varcov = zeros(ND,C,C);
cor_ID = zeros(ND,C,C);
cor_TR = zeros(ND,C,C);
loc = 1:C;                            % location-cell coordinates (first C of 2C)
for d = 1:ND
    Lid = psyrat_corr_chol_lkj(draws(d,ix.yLid).', NB, 2);   % 8x8
    cId = Lid*Lid.';                  % 8x8 correlation
    cIdL = cId(loc,loc);              % C x C location-location correlation block
    cor_ID(d,:,:) = cIdL;
    sdL = sd_id(d,loc);               % 1 x C
    id_varcov(d,:,:) = (sdL.'*sdL) .* cIdL;   % quad_form_diag(corr, sd)
    Ltr = psyrat_corr_chol_lkj(draws(d,ix.yLtrl).', NB, 2);
    cTr = Ltr*Ltr.'; cTrL = cTr(loc,loc);
    cor_TR(d,:,:) = cTrL;
    sdT = sd_trl(d,loc);
    trl_varcov(d,:,:) = (sdT.'*sdT) .* cTrL;
end
out.id_varcov = id_varcov;
out.trl_varcov = trl_varcov;
out.cor_ID_Cell = cor_ID;
out.cor_TrialNumber_Cell = cor_TR;
end

% ===========================================================================
% Design: diff_dynrel_trt (case 15) - group-level two-facet dynamic difference
% ===========================================================================
function fit = local_hmc_diffdynrel_trt(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the group-level two-facet (trial +
%occasion) dynamic difference-score model, matching
%psyrat_build_stan_diffdynrel_trt: the one-facet group-level dynamic difference
%(local_hmc_diffdynrel) generalized from TWO crossed 2-D blocks (person, trial)
%to SIX (id, trl, occ, tid, oid, to), each a 2x2 per-event Cholesky on the mean,
%plus the same event-by-dimension slopes. The residual is group-level
%(fixed-per-event): sig_err = exp(X*b_sigma + Xdim*b_sigma_dim).

m = struct();
m.NOBS  = double(data.NOBS);
m.meas  = double(data.meas(:));
m.meas1 = double(data.meas1(:));
m.meas2 = double(data.meas2(:));
% dimension design: the dynamic difference (case 15) passes Xdim_vec/KDIM; the
% STATIC two-facet difference (case 10, diff_trt) has no dimension slopes and
% passes neither, so KDIM = 0 and the Xdim*b_dim / Xdim*b_sigma_dim terms drop
% out, reducing this log-posterior to psyrat_build_stan_diff_trt exactly. The
% flatten is column-major (Xdim_vec = Xdim(:)); reshape reads column-major.
if isfield(data,'KDIM') && ~isempty(data.KDIM)
    m.KDIM = double(data.KDIM);
    m.Xdim = reshape(double(data.Xdim_vec(:)),m.NOBS,m.KDIM);
else
    m.KDIM = 0;
    m.Xdim = zeros(m.NOBS,0);
end
% the six crossed cross-condition 2x2 factors (fixed order = the parameter and
% draw-field order). Each: per-observation index, level count, half-student_t
% prior scale (priors.diff_trt.<prior>), and the draw field <name>_varcov.
m.factors = struct( ...
    'name',  {'id','trl','occ','tid','oid','to'}, ...
    'idx',   {double(data.ID(:)),double(data.TRL(:)),double(data.OCC(:)), ...
              double(data.TID(:)),double(data.OID(:)),double(data.TO(:))}, ...
    'N',     {double(data.NSUB),double(data.NTRL),double(data.NOCC), ...
              double(data.NTID),double(data.NOID),double(data.NTO)}, ...
    'prior', {10,10,5,5,5,1}, ...   % defaults: priors.diff_trt.{sd_id,sd_trl,sd_occ,sd_tid,sd_oid,sd_to}
    'field', {'id_varcov','trl_varcov','occ_varcov', ...
              'tid_varcov','oid_varcov','to_varcov'});

% configurable priors: b/b_sigma from priors.diff_trt; the six factor SD scales
% from priors.diff_trt.sd_*; the dimension slopes from priors.dynrel.
Pb = 5; Pbs = 10; Pdim = struct('b',10,'b_sigma',1);
priormap = {'sd_id','sd_trl','sd_occ','sd_tid','sd_oid','sd_to'};
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors)
    pr = design.priors;
    if isfield(pr,'diff_trt') && isstruct(pr.diff_trt)
        if isfield(pr.diff_trt,'b');       Pb  = double(pr.diff_trt.b);       end
        if isfield(pr.diff_trt,'b_sigma'); Pbs = double(pr.diff_trt.b_sigma); end
        for k = 1:6
            if isfield(pr.diff_trt,priormap{k})
                m.factors(k).prior = double(pr.diff_trt.(priormap{k}));
            end
        end
    end
    if isfield(pr,'dynrel') && isstruct(pr.dynrel)
        if isfield(pr.dynrel,'b');       Pdim.b = double(pr.dynrel.b);             end
        if isfield(pr.dynrel,'b_sigma'); Pdim.b_sigma = double(pr.dynrel.b_sigma); end
    end
end
m.Pb = Pb; m.Pbs = Pbs; m.Pdim = Pdim;

m = local_diffdynrel_trt_index(m);

gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_diffdynrel_trt_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_diffdynrel_trt_logp(t,m), th);
end

start = local_diffdynrel_trt_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_diffdynrel_trt_draws(draws,m);
% report the actual design key (case 15 'diff_dynrel_trt' or, at KDIM = 0, the
% static case-10 'diff_trt') rather than a hardcoded label.
designname = 'diff_dynrel_trt';
if isstruct(design) && isfield(design,'key') && ~isempty(design.key)
    designname = char(design.key);
end
meta = struct('design',designname,'gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_diffdynrel_trt_index(m)
%theta = [b(2); b_sigma(2); b_dim(K); b_sigma_dim(K);
%         per factor f (in order): log_sd_f(2), z_rho_f(1);
%         per factor f: z_f(2*N_f)].
K = m.KDIM; i = 0;
m.ix.b        = i+(1:2);  i = i+2;
m.ix.b_sigma  = i+(1:2);  i = i+2;
m.ix.b_dim    = i+(1:K);  i = i+K;
m.ix.bs_dim   = i+(1:K);  i = i+K;
m.ix.lsd  = cell(1,6);
m.ix.zrho = cell(1,6);
for k = 1:6
    m.ix.lsd{k}  = i+(1:2); i = i+2;
    m.ix.zrho{k} = i+1;     i = i+1;
end
m.ix.z = cell(1,6);
for k = 1:6
    m.ix.z{k} = i+(1:2*m.factors(k).N); i = i+2*m.factors(k).N;
end
m.D = i;
% monitor the population params + each factor's log-SDs and correlation
mon = [m.ix.b, m.ix.b_sigma, m.ix.b_dim, m.ix.bs_dim];
names = [{'b1','b2','b_sigma1','b_sigma2'}, ...
    local_numbered('b_dim',K), local_numbered('b_sigma_dim',K)];
for k = 1:6
    mon = [mon, m.ix.lsd{k}, m.ix.zrho{k}]; %#ok<AGROW>
    names = [names, {sprintf('log_sd_%s1',m.factors(k).name), ...
        sprintf('log_sd_%s2',m.factors(k).name), ...
        sprintf('z_rho_%s',m.factors(k).name)}]; %#ok<AGROW>
end
m.mon_idx = mon;
m.mon_names = names;
end

function start = local_diffdynrel_trt_start(m)
%Per-event means at the event-conditional data means; all factor log-SDs at a
%rough small guess; correlations and raw effects at 0.
start = zeros(m.D,1);
m1 = m.meas1 > 0; m2 = m.meas2 > 0;
start(m.ix.b(1)) = mean(m.meas(m1));
start(m.ix.b(2)) = mean(m.meas(m2));
s1 = std(m.meas(m1)); s2 = std(m.meas(m2));
start(m.ix.b_sigma(1)) = log(max(s1/2,1e-3));
start(m.ix.b_sigma(2)) = log(max(s2/2,1e-3));
% spread the marginal SD across the six crossed factors
g = log(max((s1+s2)/8,1e-3));
for k = 1:6
    start(m.ix.lsd{k}) = g;
end
% b_dim, b_sigma_dim, z_rho_*, z_* start at 0
end

function lp = local_diffdynrel_trt_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching
%psyrat_build_stan_diffdynrel_trt: the event-indicator mean with six crossed 2x2
%per-event random-effect blocks + event-by-dimension slopes, group-level residual.
ix = m.ix;
b       = theta(ix.b);
b_sigma = theta(ix.b_sigma);
b_dim   = theta(ix.b_dim);
bs_dim  = theta(ix.bs_dim);
m1 = m.meas1; m2 = m.meas2;

mu = m1.*b(1) + m2.*b(2) + m.Xdim*b_dim;
lp_re = 0;
for k = 1:6
    f = m.factors(k);
    sd  = exp(theta(ix.lsd{k}));      % 2
    rho = tanh(theta(ix.zrho{k}));
    z   = reshape(theta(ix.z{k}),2,f.N);
    r1 = sd(1).*z(1,:).';
    r2 = sd(2).*(rho.*z(1,:).' + sqrt(1-rho.^2).*z(2,:).');
    mu = mu + r1(f.idx).*m1 + r2(f.idx).*m2;
    lp_re = lp_re ...
        + local_student_t(sd(1),3,0,f.prior) + theta(ix.lsd{k}(1)) ...   % half-t + Jac
        + local_student_t(sd(2),3,0,f.prior) + theta(ix.lsd{k}(2)) ...
        + 2*log(1 - rho.^2) ...                                          % LKJ(2) + tanh Jac
        - 0.5*sum(z(:).^2);                                              % std_normal raw
end
logsig = m1.*b_sigma(1) + m2.*b_sigma(2) + m.Xdim*bs_dim;
sig = exp(logsig);
ll = -sum(logsig) - sum((m.meas - mu).^2 ./ (2*sig.^2));

Pdim = m.Pdim;
lp = ll + lp_re ...
    + sum(local_normal(b,0,m.Pb)) ...                          % b ~ normal(0,.)
    + sum(local_student_t(b_sigma,3,0,m.Pbs)) ...              % b_sigma ~ t(3,0,.)
    + sum(local_normal(b_dim,0,Pdim.b)) ...                    % b_dim
    + sum(local_normal(bs_dim,0,Pdim.b_sigma));               % b_sigma_dim
end

function out = local_diffdynrel_trt_draws(draws,m)
%Transform draws to the CmdStan-shaped arrays case 15 extracts: b, b_sigma (ND x2),
%b_dim, b_sigma_dim (ND x KDIM), and the six crossed 2x2 cross-condition
%var-covariances (<name>_varcov, ND x2x2).
ix = m.ix;
out = struct();
out.b           = draws(:,ix.b);
out.b_sigma     = draws(:,ix.b_sigma);
out.b_dim       = draws(:,ix.b_dim);
out.b_sigma_dim = draws(:,ix.bs_dim);
for k = 1:6
    sd  = exp(draws(:,ix.lsd{k}));     % ND x2
    rho = tanh(draws(:,ix.zrho{k}));   % ND x1
    out.(m.factors(k).field) = local_varcov2(sd,rho);
end
end

% ===========================================================================
% Design: diff_dynrel_sserr_trt (case 16) - subject-level two-facet dynamic
% difference
% ===========================================================================
function fit = local_hmc_diffdynrel_sserr_trt(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the SUBJECT-LEVEL two-facet dynamic
%difference-score model, matching psyrat_build_stan_diffdynrel_sserr_trt: the
%group-level two-facet dynamic difference (local_hmc_diffdynrel_trt) with the
%PERSON (id) factor widened from a 2-D location block to a 4-D joint location +
%log-residual block (a 4x4 LKJ Cholesky via psyrat_corr_chol_lkj, exactly as case
%13), so each participant has their own per-event residual. The OTHER FIVE crossed
%factors (trl/occ/tid/oid/to) stay 2-D location-only (as case 15).

m = struct();
m.NOBS  = double(data.NOBS);
m.NSUB  = double(data.NSUB);
m.meas  = double(data.meas(:));
m.meas1 = double(data.meas1(:));
m.meas2 = double(data.meas2(:));
m.id    = double(data.ID(:));
m.KDIM  = double(data.KDIM);
m.Xdim  = reshape(double(data.Xdim_vec(:)),m.NOBS,m.KDIM);
% the FIVE crossed location-only 2x2 factors (the id factor is the 4-D block,
% handled separately). Fixed order = the parameter and draw-field order.
m.factors = struct( ...
    'name',  {'trl','occ','tid','oid','to'}, ...
    'idx',   {double(data.TRL(:)),double(data.OCC(:)),double(data.TID(:)), ...
              double(data.OID(:)),double(data.TO(:))}, ...
    'N',     {double(data.NTRL),double(data.NOCC),double(data.NTID), ...
              double(data.NOID),double(data.NTO)}, ...
    'prior', {10,5,5,5,1}, ...     % priors.diff_trt.{sd_trl,sd_occ,sd_tid,sd_oid,sd_to}
    'field', {'trl_varcov','occ_varcov','tid_varcov','oid_varcov','to_varcov'});

% configurable priors: b/b_sigma from priors.diff_trt; the five crossed factor SD
% scales from priors.diff_trt.sd_*; the PERSON block location/scale SDs from
% priors.diff (sd_id / sd_id_sigma); the dimension slopes from priors.dynrel.
Pb = 5; Pbs = 10; Psd_id = 10; Psd_id_sigma = 2.5; Pdim = struct('b',10,'b_sigma',1);
priormap = {'sd_trl','sd_occ','sd_tid','sd_oid','sd_to'};
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors)
    pr = design.priors;
    if isfield(pr,'diff_trt') && isstruct(pr.diff_trt)
        if isfield(pr.diff_trt,'b');       Pb  = double(pr.diff_trt.b);       end
        if isfield(pr.diff_trt,'b_sigma'); Pbs = double(pr.diff_trt.b_sigma); end
        for k = 1:5
            if isfield(pr.diff_trt,priormap{k})
                m.factors(k).prior = double(pr.diff_trt.(priormap{k}));
            end
        end
    end
    if isfield(pr,'diff') && isstruct(pr.diff)
        if isfield(pr.diff,'sd_id');       Psd_id = double(pr.diff.sd_id);             end
        if isfield(pr.diff,'sd_id_sigma'); Psd_id_sigma = double(pr.diff.sd_id_sigma); end
    end
    if isfield(pr,'dynrel') && isstruct(pr.dynrel)
        if isfield(pr.dynrel,'b');       Pdim.b = double(pr.dynrel.b);             end
        if isfield(pr.dynrel,'b_sigma'); Pdim.b_sigma = double(pr.dynrel.b_sigma); end
    end
end
m.Pb = Pb; m.Pbs = Pbs; m.Psd_id = Psd_id; m.Psd_id_sigma = Psd_id_sigma; m.Pdim = Pdim;

m = local_diffdynrel_sserr_trt_index(m);

gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_diffdynrel_sserr_trt_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_diffdynrel_sserr_trt_logp(t,m), th);
end

start = local_diffdynrel_sserr_trt_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_diffdynrel_sserr_trt_draws(draws,m);
meta = struct('design','diff_dynrel_sserr_trt','gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_diffdynrel_sserr_trt_index(m)
%theta = [b(2); b_sigma(2); b_dim(K); b_sigma_dim(K);
%         log_sd_id(4); yLid(6);                 % 4-D person block (location+scale)
%         per crossed factor f (5): log_sd_f(2), z_rho_f(1);
%         z_id(4*NSUB);                          % person raw
%         per crossed factor f: z_f(2*N_f)].
K = m.KDIM; NSUB = m.NSUB; i = 0;
m.ix.b        = i+(1:2);  i = i+2;
m.ix.b_sigma  = i+(1:2);  i = i+2;
m.ix.b_dim    = i+(1:K);  i = i+K;
m.ix.bs_dim   = i+(1:K);  i = i+K;
m.ix.lsd_id   = i+(1:4);  i = i+4;
m.ix.yLid     = i+(1:6);  i = i+6;
m.ix.lsd  = cell(1,5);
m.ix.zrho = cell(1,5);
for k = 1:5
    m.ix.lsd{k}  = i+(1:2); i = i+2;
    m.ix.zrho{k} = i+1;     i = i+1;
end
m.ix.z_id = i+(1:4*NSUB); i = i+4*NSUB;
m.ix.z = cell(1,5);
for k = 1:5
    m.ix.z{k} = i+(1:2*m.factors(k).N); i = i+2*m.factors(k).N;
end
m.D = i;
mon = [m.ix.b, m.ix.b_sigma, m.ix.b_dim, m.ix.bs_dim, m.ix.lsd_id];
names = [{'b1','b2','b_sigma1','b_sigma2'}, ...
    local_numbered('b_dim',K), local_numbered('b_sigma_dim',K), ...
    {'log_sd_id1','log_sd_id2','log_sd_id3','log_sd_id4'}];
for k = 1:5
    mon = [mon, m.ix.lsd{k}, m.ix.zrho{k}]; %#ok<AGROW>
    names = [names, {sprintf('log_sd_%s1',m.factors(k).name), ...
        sprintf('log_sd_%s2',m.factors(k).name), ...
        sprintf('z_rho_%s',m.factors(k).name)}]; %#ok<AGROW>
end
m.mon_idx = mon; m.mon_names = names;
end

function start = local_diffdynrel_sserr_trt_start(m)
%Per-event means at the event-conditional data means; person location SDs at a
%fraction of the marginal SD, person scale SDs small; crossed-factor SDs small;
%correlations/CPCs/raw effects at 0.
start = zeros(m.D,1);
m1 = m.meas1 > 0; m2 = m.meas2 > 0;
start(m.ix.b(1)) = mean(m.meas(m1));
start(m.ix.b(2)) = mean(m.meas(m2));
s1 = std(m.meas(m1)); s2 = std(m.meas(m2));
start(m.ix.b_sigma(1)) = log(max(s1/2,1e-3));
start(m.ix.b_sigma(2)) = log(max(s2/2,1e-3));
start(m.ix.lsd_id(1:2)) = log(max((s1+s2)/6,1e-3));   % person LOCATION SDs
start(m.ix.lsd_id(3:4)) = log(0.3);                   % person SCALE SDs
g = log(max((s1+s2)/8,1e-3));
for k = 1:5
    start(m.ix.lsd{k}) = g;
end
end

function lp = local_diffdynrel_sserr_trt_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching
%psyrat_build_stan_diffdynrel_sserr_trt: the 4-D person block (location + scale)
%plus the five crossed 2x2 location factors + event-by-dimension slopes.
ix = m.ix;
b       = theta(ix.b);
b_sigma = theta(ix.b_sigma);
b_dim   = theta(ix.b_dim);
bs_dim  = theta(ix.bs_dim);
m1 = m.meas1; m2 = m.meas2;

% 4-D person block (location events 1-2, scale events 3-4)
sd_id = exp(theta(ix.lsd_id));
[L_id, lp_Lid] = psyrat_corr_chol_lkj(theta(ix.yLid), 4, 2);
z_id = reshape(theta(ix.z_id),4,m.NSUB);
r1 = (sd_id .* (L_id * z_id)).';      % NSUB x 4
r1_1 = r1(:,1); r1_2 = r1(:,2);       % location REs
r1_s1 = r1(:,3); r1_s2 = r1(:,4);     % scale (log-residual) REs

mu = m1.*b(1) + m2.*b(2) + m.Xdim*b_dim ...
   + r1_1(m.id).*m1 + r1_2(m.id).*m2;
logsig = m1.*b_sigma(1) + m2.*b_sigma(2) + m.Xdim*bs_dim ...
   + r1_s1(m.id).*m1 + r1_s2(m.id).*m2;

% five crossed 2x2 location-only factors
lp_re = 0;
for k = 1:5
    f = m.factors(k);
    sd  = exp(theta(ix.lsd{k}));
    rho = tanh(theta(ix.zrho{k}));
    z   = reshape(theta(ix.z{k}),2,f.N);
    rA = sd(1).*z(1,:).';
    rB = sd(2).*(rho.*z(1,:).' + sqrt(1-rho.^2).*z(2,:).');
    mu = mu + rA(f.idx).*m1 + rB(f.idx).*m2;
    lp_re = lp_re ...
        + local_student_t(sd(1),3,0,f.prior) + theta(ix.lsd{k}(1)) ...
        + local_student_t(sd(2),3,0,f.prior) + theta(ix.lsd{k}(2)) ...
        + 2*log(1 - rho.^2) ...
        - 0.5*sum(z(:).^2);
end

sig = exp(logsig);
ll = -sum(logsig) - sum((m.meas - mu).^2 ./ (2*sig.^2));

Pdim = m.Pdim;
lp = ll + lp_re + lp_Lid ...                                  % L_id LKJ(2) + Jac
    + sum(local_normal(b,0,m.Pb)) ...                         % b ~ normal(0,.)
    + sum(local_student_t(b_sigma,3,0,m.Pbs)) ...             % b_sigma ~ t(3,0,.)
    + sum(local_normal(b_dim,0,Pdim.b)) ...                   % b_dim
    + sum(local_normal(bs_dim,0,Pdim.b_sigma)) ...            % b_sigma_dim
    + local_student_t(sd_id(1),3,0,m.Psd_id)       + theta(ix.lsd_id(1)) ...  % person location SDs
    + local_student_t(sd_id(2),3,0,m.Psd_id)       + theta(ix.lsd_id(2)) ...
    + local_student_t(sd_id(3),3,0,m.Psd_id_sigma) + theta(ix.lsd_id(3)) ...  % person scale SDs
    + local_student_t(sd_id(4),3,0,m.Psd_id_sigma) + theta(ix.lsd_id(4)) ...
    - 0.5*sum(z_id(:).^2);                                    % std_normal person raw
end

function out = local_diffdynrel_sserr_trt_draws(draws,m)
%Transform draws to the CmdStan-shaped arrays case 16 extracts: b, b_sigma (ND x2),
%b_dim, b_sigma_dim (ND x KDIM), sd_id (ND x4), L_id (ND x4x4), id_varcov (ND x2x2,
%location block), the five crossed 2x2 var-covariances (trl/occ/tid/oid/to), and
%er_var_ss (ND x NSUB x 2).
ix = m.ix;
ND = size(draws,1); NSUB = m.NSUB;
out = struct();
out.b           = draws(:,ix.b);
out.b_sigma     = draws(:,ix.b_sigma);
out.b_dim       = draws(:,ix.b_dim);
out.b_sigma_dim = draws(:,ix.bs_dim);
sd_id = exp(draws(:,ix.lsd_id));      % ND x4
out.sd_id = sd_id;
% per-draw 4x4 person Cholesky -> L_id, id_varcov (location block), er_var_ss
L_id = zeros(ND,4,4);
id_varcov = zeros(ND,2,2);
er = zeros(ND,NSUB,2);
for d = 1:ND
    Ld = psyrat_corr_chol_lkj(draws(d,ix.yLid).',4,2);
    L_id(d,:,:) = Ld;
    cor = Ld*Ld.';
    sdd = sd_id(d,:);
    id_varcov(d,1,1) = sdd(1)^2;
    id_varcov(d,2,2) = sdd(2)^2;
    id_varcov(d,1,2) = sdd(1)*sdd(2)*cor(1,2);
    id_varcov(d,2,1) = id_varcov(d,1,2);
    z4 = reshape(draws(d,ix.z_id),4,NSUB);
    r1 = (sdd(:) .* (Ld * z4)).';     % NSUB x 4
    er(d,:,1) = draws(d,ix.b_sigma(1)) + r1(:,3).';
    er(d,:,2) = draws(d,ix.b_sigma(2)) + r1(:,4).';
end
out.L_id = L_id;
out.id_varcov = id_varcov;
out.er_var_ss = er;
% five crossed 2x2 var-covariances
for k = 1:5
    sd  = exp(draws(:,ix.lsd{k}));
    rho = tanh(draws(:,ix.zrho{k}));
    out.(m.factors(k).field) = local_varcov2(sd,rho);
end
end

% ===========================================================================
% Design: diff_dynrel_rescor (case 17) - concurrent one-facet group-level
% dynamic difference (bivariate residual)
% ===========================================================================
function fit = local_hmc_diffdynrel_rescor(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the CONCURRENT one-facet group-level
%dynamic difference model, matching psyrat_build_stan_diffdynrel_rescor: the
%bivariate-residual difference model (local_hmc_diff_rescor, case 7*) EXTENDED
%with per-event fixed dimension slopes on the mean (b_dim1/b_dim2) AND on the
%per-event log-residual SD (b_sigma_dim1/b_sigma_dim2). The per-event residual SDs
%are per-row (dimension-dependent), so the bivariate residual Cholesky LSigma_n =
%diag([sigma1_n, sigma2_n]) * Lrescor varies by row.

m = struct();
m.NOBS  = double(data.NOBS);
m.NSUB  = double(data.NSUB);
m.NTRL  = double(data.NTRL);
m.meas1 = double(data.meas1(:));   % event-1 response (paired)
m.meas2 = double(data.meas2(:));   % event-2 response (paired)
m.id    = double(data.ID(:));
m.trl   = double(data.TRL(:));
m.KDIM  = double(data.KDIM);
m.Xdim  = reshape(double(data.Xdim_vec(:)),m.NOBS,m.KDIM);  % between-person dim design

% configurable scale priors: base diff block from priors.diff; per-event dimension
% slopes from priors.dynrel (b for the mean slopes, b_sigma for the log-resid
% slopes). Lrescor is LKJ(2), no configurable scale.
P = struct('b',5,'b_sigma',10,'sd_id',10,'sd_trl',10);
Pdim = struct('b',10,'b_sigma',1);
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors)
    pr = design.priors;
    if isfield(pr,'diff') && isstruct(pr.diff)
        fn = {'b','b_sigma','sd_id','sd_trl'};
        for k = 1:numel(fn)
            if isfield(pr.diff,fn{k}); P.(fn{k}) = double(pr.diff.(fn{k})); end
        end
    end
    if isfield(pr,'dynrel') && isstruct(pr.dynrel)
        if isfield(pr.dynrel,'b');       Pdim.b = double(pr.dynrel.b);             end
        if isfield(pr.dynrel,'b_sigma'); Pdim.b_sigma = double(pr.dynrel.b_sigma); end
    end
end
m.P = P; m.Pdim = Pdim;

m = local_diffdynrel_rescor_index(m);

gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_diffdynrel_rescor_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_diffdynrel_rescor_logp(t,m), th);
end

start = local_diffdynrel_rescor_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_diffdynrel_rescor_draws(draws,m);
meta = struct('design','diff_dynrel_rescor','gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_diffdynrel_rescor_index(m)
%theta = [b(2); b_sigma(2); b_dim1(K); b_dim2(K); b_sigma_dim1(K); b_sigma_dim2(K);
%         log_sd_id(2); z_rho_id; log_sd_trl(2); z_rho_trl; z_rho_res;
%         z_id(2*NSUB); z_trl(2*NTRL)]. Case 7* parameters plus the four per-event
%dimension-slope vectors.
K = m.KDIM; NSUB = m.NSUB; NTRL = m.NTRL; i = 0;
m.ix.b        = i+(1:2);  i = i+2;
m.ix.b_sigma  = i+(1:2);  i = i+2;
m.ix.b_dim1   = i+(1:K);  i = i+K;
m.ix.b_dim2   = i+(1:K);  i = i+K;
m.ix.bs_dim1  = i+(1:K);  i = i+K;
m.ix.bs_dim2  = i+(1:K);  i = i+K;
m.ix.lsd_id   = i+(1:2);  i = i+2;
m.ix.zrho_id  = i+1;      i = i+1;
m.ix.lsd_trl  = i+(1:2);  i = i+2;
m.ix.zrho_trl = i+1;      i = i+1;
m.ix.zrho_res = i+1;      i = i+1;
m.ix.z_id     = i+(1:2*NSUB); i = i+2*NSUB;
m.ix.z_trl    = i+(1:2*NTRL); i = i+2*NTRL;
m.D = i;
m.mon_idx = [m.ix.b, m.ix.b_sigma, m.ix.b_dim1, m.ix.b_dim2, ...
    m.ix.bs_dim1, m.ix.bs_dim2, m.ix.lsd_id, m.ix.zrho_id, ...
    m.ix.lsd_trl, m.ix.zrho_trl, m.ix.zrho_res];
m.mon_names = [{'b1','b2','b_sigma1','b_sigma2'}, ...
    local_numbered('b_dim1',K), local_numbered('b_dim2',K), ...
    local_numbered('b_sigma_dim1',K), local_numbered('b_sigma_dim2',K), ...
    {'log_sd_id1','log_sd_id2','z_rho_id','log_sd_trl1','log_sd_trl2', ...
     'z_rho_trl','z_rho_res'}];
end

function start = local_diffdynrel_rescor_start(m)
%Per-event means at the per-event data means; log-SDs at rough guesses; dimension
%slopes/CPCs/raw effects at 0.
start = zeros(m.D,1);
start(m.ix.b(1)) = mean(m.meas1);
start(m.ix.b(2)) = mean(m.meas2);
s1 = std(m.meas1); s2 = std(m.meas2);
start(m.ix.b_sigma(1)) = log(max(s1/2,1e-3));
start(m.ix.b_sigma(2)) = log(max(s2/2,1e-3));
start(m.ix.lsd_id)  = log(max((s1+s2)/4,1e-3));
start(m.ix.lsd_trl) = log(max((s1+s2)/4,1e-3));
end

function lp = local_diffdynrel_rescor_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching
%psyrat_build_stan_diffdynrel_rescor: local_diff_rescor_logp with per-event
%dimension slopes on the mean and the per-row log-residual SDs.
ix = m.ix;
b       = theta(ix.b);
b_sigma = theta(ix.b_sigma);
b_dim1  = theta(ix.b_dim1);  b_dim2  = theta(ix.b_dim2);
bs_dim1 = theta(ix.bs_dim1); bs_dim2 = theta(ix.bs_dim2);
sd_id   = exp(theta(ix.lsd_id));
rho_id  = tanh(theta(ix.zrho_id));
sd_trl  = exp(theta(ix.lsd_trl));
rho_trl = tanh(theta(ix.zrho_trl));
rho_res = tanh(theta(ix.zrho_res));
z_id  = reshape(theta(ix.z_id),2,m.NSUB);
z_trl = reshape(theta(ix.z_trl),2,m.NTRL);

% non-centered person/trial event-effects (same as case 7*)
r1_1 = sd_id(1).*z_id(1,:).';
r1_2 = sd_id(2).*(rho_id.*z_id(1,:).' + sqrt(1-rho_id.^2).*z_id(2,:).');
r2_1 = sd_trl(1).*z_trl(1,:).';
r2_2 = sd_trl(2).*(rho_trl.*z_trl(1,:).' + sqrt(1-rho_trl.^2).*z_trl(2,:).');

% dimension-conditioned per-event means + per-row log-residual SDs
mu1 = b(1) + m.Xdim*b_dim1 + r1_1(m.id) + r2_1(m.trl);
mu2 = b(2) + m.Xdim*b_dim2 + r1_2(m.id) + r2_2(m.trl);
logsig1 = b_sigma(1) + m.Xdim*bs_dim1;       % NOBS-vector
logsig2 = b_sigma(2) + m.Xdim*bs_dim2;
sig1 = exp(logsig1); sig2 = exp(logsig2);
d1 = m.meas1 - mu1; d2 = m.meas2 - mu2;

% per-row bivariate residual: LSigma_n = diag([sig1_n,sig2_n]) * Lrescor.
% multi_normal_cholesky: -sum(log det LSigma_n) - 1/2 sum ||LSigma_n\(Y-Mu)||^2
% (drop -K/2*log(2*pi)). log det LSigma_n = logsig1_n + logsig2_n + 1/2 log(1-rho^2).
w1 = d1./sig1;
w2 = (d2 - (sig2.*rho_res./sig1).*d1) ./ (sig2.*sqrt(1-rho_res.^2));
ll = -sum(logsig1) - sum(logsig2) - 0.5*m.NOBS*log(1-rho_res.^2) ...
   - 0.5*sum(w1.^2 + w2.^2);

P = m.P; Pdim = m.Pdim;
lp = ll ...
    + sum(local_normal(b,0,P.b)) ...                          % b ~ normal(0,5)
    + sum(local_student_t(b_sigma,3,0,P.b_sigma)) ...         % b_sigma ~ t(3,0,10)
    + sum(local_normal(b_dim1,0,Pdim.b)) + sum(local_normal(b_dim2,0,Pdim.b)) ...
    + sum(local_normal(bs_dim1,0,Pdim.b_sigma)) + sum(local_normal(bs_dim2,0,Pdim.b_sigma)) ...
    + local_student_t(sd_id(1),3,0,P.sd_id)   + theta(ix.lsd_id(1)) ...
    + local_student_t(sd_id(2),3,0,P.sd_id)   + theta(ix.lsd_id(2)) ...
    + local_student_t(sd_trl(1),3,0,P.sd_trl) + theta(ix.lsd_trl(1)) ...
    + local_student_t(sd_trl(2),3,0,P.sd_trl) + theta(ix.lsd_trl(2)) ...
    + 2*log(1 - rho_id.^2) ...                                % L_id  LKJ(2) + tanh Jac
    + 2*log(1 - rho_trl.^2) ...                               % L_trl LKJ(2) + tanh Jac
    + 2*log(1 - rho_res.^2) ...                               % Lrescor LKJ(2) + tanh Jac
    - 0.5*sum(z_id(:).^2) - 0.5*sum(z_trl(:).^2);             % std_normal raw
end

function out = local_diffdynrel_rescor_draws(draws,m)
%Transform draws to the CmdStan-shaped arrays case 17 extracts: b, b_sigma (ND x2),
%b_dim1/b_dim2/b_sigma_dim1/b_sigma_dim2 (ND x KDIM), id_varcov/trl_varcov (ND x2x2),
%rescor (ND x1), err_varcov (ND x2x2, at z = 0).
ix = m.ix;
out = struct();
out.b       = draws(:,ix.b);
out.b_sigma = draws(:,ix.b_sigma);
out.b_dim1  = draws(:,ix.b_dim1);
out.b_dim2  = draws(:,ix.b_dim2);
out.b_sigma_dim1 = draws(:,ix.bs_dim1);
out.b_sigma_dim2 = draws(:,ix.bs_dim2);
sd_id  = exp(draws(:,ix.lsd_id));
sd_trl = exp(draws(:,ix.lsd_trl));
rho_id  = tanh(draws(:,ix.zrho_id));
rho_trl = tanh(draws(:,ix.zrho_trl));
rho_res = tanh(draws(:,ix.zrho_res));
out.id_varcov  = local_varcov2(sd_id,rho_id);
out.trl_varcov = local_varcov2(sd_trl,rho_trl);
out.rescor = rho_res;                          % ND x1
sigma = exp(draws(:,ix.b_sigma));              % residual SDs at z = 0
out.err_varcov = local_varcov2(sigma,rho_res); % at z = 0
end

% ===========================================================================
% Design: diff_dynrel_trt_rescor (case 18) - concurrent two-facet group-level
% dynamic difference (bivariate residual)
% ===========================================================================
function fit = local_hmc_diffdynrel_trt_rescor(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the CONCURRENT two-facet group-level dynamic
%difference model, matching psyrat_build_stan_diffdynrel_trt_rescor: the concurrent
%one-facet model (local_hmc_diffdynrel_rescor, case 17) with the person+trial 2-D
%blocks generalized to the SIX crossed cross-condition 2-D factors (id/trl/occ/
%tid/oid/to, as case 15) in the concurrent (paired) form (each factor's event-1
%effect loads mu1, its event-2 effect loads mu2); same per-row bivariate residual
%(Lrescor) + per-event dimension slopes as case 17.

m = struct();
m.NOBS  = double(data.NOBS);
m.meas1 = double(data.meas1(:));
m.meas2 = double(data.meas2(:));
% dimension design: the dynamic concurrent difference (case 18) passes
% Xdim_vec/KDIM; the STATIC concurrent two-facet difference (case 10,
% diff_trt_rescor) has no dimension slopes and passes neither, so KDIM = 0 and
% the per-event slope terms drop out, reducing this log-posterior to
% psyrat_build_stan_diff_trt_rescor exactly (the per-row residual collapses to a
% constant bivariate Cholesky). Column-major flatten, column-major reshape.
if isfield(data,'KDIM') && ~isempty(data.KDIM)
    m.KDIM = double(data.KDIM);
    m.Xdim = reshape(double(data.Xdim_vec(:)),m.NOBS,m.KDIM);
else
    m.KDIM = 0;
    m.Xdim = zeros(m.NOBS,0);
end
m.factors = struct( ...
    'name',  {'id','trl','occ','tid','oid','to'}, ...
    'idx',   {double(data.ID(:)),double(data.TRL(:)),double(data.OCC(:)), ...
              double(data.TID(:)),double(data.OID(:)),double(data.TO(:))}, ...
    'N',     {double(data.NSUB),double(data.NTRL),double(data.NOCC), ...
              double(data.NTID),double(data.NOID),double(data.NTO)}, ...
    'prior', {10,10,5,5,5,1}, ...
    'field', {'id_varcov','trl_varcov','occ_varcov', ...
              'tid_varcov','oid_varcov','to_varcov'});

Pb = 5; Pbs = 10; Pdim = struct('b',10,'b_sigma',1);
priormap = {'sd_id','sd_trl','sd_occ','sd_tid','sd_oid','sd_to'};
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors)
    pr = design.priors;
    if isfield(pr,'diff_trt') && isstruct(pr.diff_trt)
        if isfield(pr.diff_trt,'b');       Pb  = double(pr.diff_trt.b);       end
        if isfield(pr.diff_trt,'b_sigma'); Pbs = double(pr.diff_trt.b_sigma); end
        for k = 1:6
            if isfield(pr.diff_trt,priormap{k})
                m.factors(k).prior = double(pr.diff_trt.(priormap{k}));
            end
        end
    end
    if isfield(pr,'dynrel') && isstruct(pr.dynrel)
        if isfield(pr.dynrel,'b');       Pdim.b = double(pr.dynrel.b);             end
        if isfield(pr.dynrel,'b_sigma'); Pdim.b_sigma = double(pr.dynrel.b_sigma); end
    end
end
m.Pb = Pb; m.Pbs = Pbs; m.Pdim = Pdim;

m = local_diffdynrel_trt_rescor_index(m);

gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_diffdynrel_trt_rescor_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_diffdynrel_trt_rescor_logp(t,m), th);
end

start = local_diffdynrel_trt_rescor_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_diffdynrel_trt_rescor_draws(draws,m);
% report the actual design key (case 18 'diff_dynrel_trt_rescor' or, at KDIM = 0,
% the static case-10 'diff_trt_rescor') rather than a hardcoded label.
designname = 'diff_dynrel_trt_rescor';
if isstruct(design) && isfield(design,'key') && ~isempty(design.key)
    designname = char(design.key);
end
meta = struct('design',designname,'gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_diffdynrel_trt_rescor_index(m)
%theta = [b(2); b_sigma(2); b_dim1(K); b_dim2(K); b_sigma_dim1(K); b_sigma_dim2(K);
%         per factor f (6): log_sd_f(2), z_rho_f(1);
%         z_rho_res;
%         per factor f: z_f(2*N_f)].
K = m.KDIM; i = 0;
m.ix.b        = i+(1:2);  i = i+2;
m.ix.b_sigma  = i+(1:2);  i = i+2;
m.ix.b_dim1   = i+(1:K);  i = i+K;
m.ix.b_dim2   = i+(1:K);  i = i+K;
m.ix.bs_dim1  = i+(1:K);  i = i+K;
m.ix.bs_dim2  = i+(1:K);  i = i+K;
m.ix.lsd  = cell(1,6);
m.ix.zrho = cell(1,6);
for k = 1:6
    m.ix.lsd{k}  = i+(1:2); i = i+2;
    m.ix.zrho{k} = i+1;     i = i+1;
end
m.ix.zrho_res = i+1; i = i+1;
m.ix.z = cell(1,6);
for k = 1:6
    m.ix.z{k} = i+(1:2*m.factors(k).N); i = i+2*m.factors(k).N;
end
m.D = i;
mon = [m.ix.b, m.ix.b_sigma, m.ix.b_dim1, m.ix.b_dim2, m.ix.bs_dim1, ...
    m.ix.bs_dim2, m.ix.zrho_res];
names = [{'b1','b2','b_sigma1','b_sigma2'}, ...
    local_numbered('b_dim1',K), local_numbered('b_dim2',K), ...
    local_numbered('b_sigma_dim1',K), local_numbered('b_sigma_dim2',K), {'z_rho_res'}];
for k = 1:6
    mon = [mon, m.ix.lsd{k}, m.ix.zrho{k}]; %#ok<AGROW>
    names = [names, {sprintf('log_sd_%s1',m.factors(k).name), ...
        sprintf('log_sd_%s2',m.factors(k).name), ...
        sprintf('z_rho_%s',m.factors(k).name)}]; %#ok<AGROW>
end
m.mon_idx = mon; m.mon_names = names;
end

function start = local_diffdynrel_trt_rescor_start(m)
start = zeros(m.D,1);
start(m.ix.b(1)) = mean(m.meas1);
start(m.ix.b(2)) = mean(m.meas2);
s1 = std(m.meas1); s2 = std(m.meas2);
start(m.ix.b_sigma(1)) = log(max(s1/2,1e-3));
start(m.ix.b_sigma(2)) = log(max(s2/2,1e-3));
g = log(max((s1+s2)/8,1e-3));
for k = 1:6
    start(m.ix.lsd{k}) = g;
end
end

function lp = local_diffdynrel_trt_rescor_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching
%psyrat_build_stan_diffdynrel_trt_rescor: the six crossed 2x2 factors in concurrent
%(paired) form + per-event dimension slopes + per-row bivariate residual (Lrescor).
ix = m.ix;
b       = theta(ix.b);
b_sigma = theta(ix.b_sigma);
b_dim1  = theta(ix.b_dim1);  b_dim2  = theta(ix.b_dim2);
bs_dim1 = theta(ix.bs_dim1); bs_dim2 = theta(ix.bs_dim2);
rho_res = tanh(theta(ix.zrho_res));

mu1 = b(1) + m.Xdim*b_dim1;
mu2 = b(2) + m.Xdim*b_dim2;
lp_re = 0;
for k = 1:6
    f = m.factors(k);
    sd  = exp(theta(ix.lsd{k}));
    rho = tanh(theta(ix.zrho{k}));
    z   = reshape(theta(ix.z{k}),2,f.N);
    rA = sd(1).*z(1,:).';                                         % event-1 effect
    rB = sd(2).*(rho.*z(1,:).' + sqrt(1-rho.^2).*z(2,:).');       % event-2 effect
    mu1 = mu1 + rA(f.idx);
    mu2 = mu2 + rB(f.idx);
    lp_re = lp_re ...
        + local_student_t(sd(1),3,0,f.prior) + theta(ix.lsd{k}(1)) ...
        + local_student_t(sd(2),3,0,f.prior) + theta(ix.lsd{k}(2)) ...
        + 2*log(1 - rho.^2) ...
        - 0.5*sum(z(:).^2);
end

logsig1 = b_sigma(1) + m.Xdim*bs_dim1;
logsig2 = b_sigma(2) + m.Xdim*bs_dim2;
sig1 = exp(logsig1); sig2 = exp(logsig2);
d1 = m.meas1 - mu1; d2 = m.meas2 - mu2;
% per-row bivariate residual LSigma_n = diag([sig1_n,sig2_n]) * Lrescor (case 17)
w1 = d1./sig1;
w2 = (d2 - (sig2.*rho_res./sig1).*d1) ./ (sig2.*sqrt(1-rho_res.^2));
ll = -sum(logsig1) - sum(logsig2) - 0.5*m.NOBS*log(1-rho_res.^2) ...
   - 0.5*sum(w1.^2 + w2.^2);

Pdim = m.Pdim;
lp = ll + lp_re ...
    + sum(local_normal(b,0,m.Pb)) ...
    + sum(local_student_t(b_sigma,3,0,m.Pbs)) ...
    + sum(local_normal(b_dim1,0,Pdim.b)) + sum(local_normal(b_dim2,0,Pdim.b)) ...
    + sum(local_normal(bs_dim1,0,Pdim.b_sigma)) + sum(local_normal(bs_dim2,0,Pdim.b_sigma)) ...
    + 2*log(1 - rho_res.^2);                                      % Lrescor LKJ(2) + tanh Jac
end

function out = local_diffdynrel_trt_rescor_draws(draws,m)
%Transform draws to the CmdStan-shaped arrays case 18 extracts: b, b_sigma (ND x2),
%b_dim1/b_dim2/b_sigma_dim1/b_sigma_dim2 (ND x KDIM), the six crossed 2x2
%var-covariances, rescor (ND x1), err_varcov (ND x2x2 at z = 0).
ix = m.ix;
out = struct();
out.b           = draws(:,ix.b);
out.b_sigma     = draws(:,ix.b_sigma);
out.b_dim1      = draws(:,ix.b_dim1);
out.b_dim2      = draws(:,ix.b_dim2);
out.b_sigma_dim1 = draws(:,ix.bs_dim1);
out.b_sigma_dim2 = draws(:,ix.bs_dim2);
for k = 1:6
    sd  = exp(draws(:,ix.lsd{k}));
    rho = tanh(draws(:,ix.zrho{k}));
    out.(m.factors(k).field) = local_varcov2(sd,rho);
end
rho_res = tanh(draws(:,ix.zrho_res));
out.rescor = rho_res;
sigma = exp(draws(:,ix.b_sigma));
out.err_varcov = local_varcov2(sigma,rho_res);   % at z = 0
end

% ===========================================================================
% Design: diff_dynrel_sserr_rescor (case 19) - concurrent one-facet subject-level
% dynamic difference (bivariate residual + 4-D person block)
% ===========================================================================
function fit = local_hmc_diffdynrel_sserr_rescor(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the CONCURRENT one-facet SUBJECT-LEVEL
%dynamic difference model, matching psyrat_build_stan_diffdynrel_sserr_rescor: the
%concurrent one-facet group model (local_hmc_diffdynrel_rescor, case 17) with the
%person block widened from 2-D to a 4-D joint location+scale block (a 4x4 LKJ
%Cholesky, as case 13). The person scale REs (coords 3-4) enter the per-row
%residual SDs sigma1_n/sigma2_n = exp(b_sigma[e] + Xdim*b_sigma_dim_e + r_1_s_e[id]);
%the residual correlation (Lrescor) stays population.

m = struct();
m.NOBS  = double(data.NOBS);
m.NSUB  = double(data.NSUB);
m.NTRL  = double(data.NTRL);
m.meas1 = double(data.meas1(:));
m.meas2 = double(data.meas2(:));
m.id    = double(data.ID(:));
m.trl   = double(data.TRL(:));
m.KDIM  = double(data.KDIM);
m.Xdim  = reshape(double(data.Xdim_vec(:)),m.NOBS,m.KDIM);

% configurable priors: base diff block from priors.diff (b/b_sigma/sd_trl, person
% LOCATION sd_id, person SCALE sd_id_sigma); dimension slopes from priors.dynrel.
P = struct('b',5,'b_sigma',10,'sd_id',10,'sd_trl',10,'sd_id_sigma',2.5);
Pdim = struct('b',10,'b_sigma',1);
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors)
    pr = design.priors;
    if isfield(pr,'diff') && isstruct(pr.diff)
        fn = {'b','b_sigma','sd_id','sd_trl','sd_id_sigma'};
        for k = 1:numel(fn)
            if isfield(pr.diff,fn{k}); P.(fn{k}) = double(pr.diff.(fn{k})); end
        end
    end
    if isfield(pr,'dynrel') && isstruct(pr.dynrel)
        if isfield(pr.dynrel,'b');       Pdim.b = double(pr.dynrel.b);             end
        if isfield(pr.dynrel,'b_sigma'); Pdim.b_sigma = double(pr.dynrel.b_sigma); end
    end
end
m.P = P; m.Pdim = Pdim;

m = local_diffdynrel_sserr_rescor_index(m);

gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_diffdynrel_sserr_rescor_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_diffdynrel_sserr_rescor_logp(t,m), th);
end

start = local_diffdynrel_sserr_rescor_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_diffdynrel_sserr_rescor_draws(draws,m);
meta = struct('design','diff_dynrel_sserr_rescor','gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_diffdynrel_sserr_rescor_index(m)
%theta = [b(2); b_sigma(2); b_dim1(K); b_dim2(K); b_sigma_dim1(K); b_sigma_dim2(K);
%         log_sd_id(4); yLid(6);          % 4-D person block (location+scale)
%         log_sd_trl(2); z_rho_trl;        % trial 2-D location block
%         z_rho_res;                       % population Lrescor
%         z_id(4*NSUB); z_trl(2*NTRL)].
K = m.KDIM; NSUB = m.NSUB; NTRL = m.NTRL; i = 0;
m.ix.b        = i+(1:2);  i = i+2;
m.ix.b_sigma  = i+(1:2);  i = i+2;
m.ix.b_dim1   = i+(1:K);  i = i+K;
m.ix.b_dim2   = i+(1:K);  i = i+K;
m.ix.bs_dim1  = i+(1:K);  i = i+K;
m.ix.bs_dim2  = i+(1:K);  i = i+K;
m.ix.lsd_id   = i+(1:4);  i = i+4;
m.ix.yLid     = i+(1:6);  i = i+6;
m.ix.lsd_trl  = i+(1:2);  i = i+2;
m.ix.zrho_trl = i+1;      i = i+1;
m.ix.zrho_res = i+1;      i = i+1;
m.ix.z_id     = i+(1:4*NSUB); i = i+4*NSUB;
m.ix.z_trl    = i+(1:2*NTRL); i = i+2*NTRL;
m.D = i;
m.mon_idx = [m.ix.b, m.ix.b_sigma, m.ix.b_dim1, m.ix.b_dim2, ...
    m.ix.bs_dim1, m.ix.bs_dim2, m.ix.lsd_id, m.ix.lsd_trl, ...
    m.ix.zrho_trl, m.ix.zrho_res];
m.mon_names = [{'b1','b2','b_sigma1','b_sigma2'}, ...
    local_numbered('b_dim1',K), local_numbered('b_dim2',K), ...
    local_numbered('b_sigma_dim1',K), local_numbered('b_sigma_dim2',K), ...
    {'log_sd_id1','log_sd_id2','log_sd_id3','log_sd_id4', ...
     'log_sd_trl1','log_sd_trl2','z_rho_trl','z_rho_res'}];
end

function start = local_diffdynrel_sserr_rescor_start(m)
start = zeros(m.D,1);
start(m.ix.b(1)) = mean(m.meas1);
start(m.ix.b(2)) = mean(m.meas2);
s1 = std(m.meas1); s2 = std(m.meas2);
start(m.ix.b_sigma(1)) = log(max(s1/2,1e-3));
start(m.ix.b_sigma(2)) = log(max(s2/2,1e-3));
start(m.ix.lsd_id(1:2)) = log(max((s1+s2)/6,1e-3));  % person LOCATION SDs
start(m.ix.lsd_id(3:4)) = log(0.3);                  % person SCALE SDs
start(m.ix.lsd_trl)     = log(max((s1+s2)/4,1e-3));
end

function lp = local_diffdynrel_sserr_rescor_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching
%psyrat_build_stan_diffdynrel_sserr_rescor: case-17 bivariate residual + dimension
%slopes with the person block widened to the case-13 4-D block (scale REs feed the
%per-row residual SDs).
ix = m.ix;
b       = theta(ix.b);
b_sigma = theta(ix.b_sigma);
b_dim1  = theta(ix.b_dim1);  b_dim2  = theta(ix.b_dim2);
bs_dim1 = theta(ix.bs_dim1); bs_dim2 = theta(ix.bs_dim2);
% 4-D person block (location 1-2, scale 3-4)
sd_id = exp(theta(ix.lsd_id));
[L_id, lp_Lid] = psyrat_corr_chol_lkj(theta(ix.yLid), 4, 2);
z_id = reshape(theta(ix.z_id),4,m.NSUB);
r1 = (sd_id .* (L_id * z_id)).';      % NSUB x 4
r1_1 = r1(:,1); r1_2 = r1(:,2); r1_s1 = r1(:,3); r1_s2 = r1(:,4);
% trial 2-D location block
sd_trl  = exp(theta(ix.lsd_trl));
rho_trl = tanh(theta(ix.zrho_trl));
z_trl = reshape(theta(ix.z_trl),2,m.NTRL);
r2_1 = sd_trl(1).*z_trl(1,:).';
r2_2 = sd_trl(2).*(rho_trl.*z_trl(1,:).' + sqrt(1-rho_trl.^2).*z_trl(2,:).');
rho_res = tanh(theta(ix.zrho_res));

mu1 = b(1) + m.Xdim*b_dim1 + r1_1(m.id) + r2_1(m.trl);
mu2 = b(2) + m.Xdim*b_dim2 + r1_2(m.id) + r2_2(m.trl);
% per-row residual log-SDs: population + dim slopes + per-subject scale RE
logsig1 = b_sigma(1) + m.Xdim*bs_dim1 + r1_s1(m.id);
logsig2 = b_sigma(2) + m.Xdim*bs_dim2 + r1_s2(m.id);
sig1 = exp(logsig1); sig2 = exp(logsig2);
d1 = m.meas1 - mu1; d2 = m.meas2 - mu2;
% per-row bivariate residual (case 17)
w1 = d1./sig1;
w2 = (d2 - (sig2.*rho_res./sig1).*d1) ./ (sig2.*sqrt(1-rho_res.^2));
ll = -sum(logsig1) - sum(logsig2) - 0.5*m.NOBS*log(1-rho_res.^2) ...
   - 0.5*sum(w1.^2 + w2.^2);

P = m.P; Pdim = m.Pdim;
lp = ll + lp_Lid ...                                          % L_id LKJ(2) + Jac
    + sum(local_normal(b,0,P.b)) ...
    + sum(local_student_t(b_sigma,3,0,P.b_sigma)) ...
    + sum(local_normal(b_dim1,0,Pdim.b)) + sum(local_normal(b_dim2,0,Pdim.b)) ...
    + sum(local_normal(bs_dim1,0,Pdim.b_sigma)) + sum(local_normal(bs_dim2,0,Pdim.b_sigma)) ...
    + local_student_t(sd_id(1),3,0,P.sd_id)       + theta(ix.lsd_id(1)) ...  % person location SDs
    + local_student_t(sd_id(2),3,0,P.sd_id)       + theta(ix.lsd_id(2)) ...
    + local_student_t(sd_id(3),3,0,P.sd_id_sigma) + theta(ix.lsd_id(3)) ...  % person scale SDs
    + local_student_t(sd_id(4),3,0,P.sd_id_sigma) + theta(ix.lsd_id(4)) ...
    + local_student_t(sd_trl(1),3,0,P.sd_trl) + theta(ix.lsd_trl(1)) ...
    + local_student_t(sd_trl(2),3,0,P.sd_trl) + theta(ix.lsd_trl(2)) ...
    + 2*log(1 - rho_trl.^2) ...                               % L_trl LKJ(2) + tanh Jac
    + 2*log(1 - rho_res.^2) ...                               % Lrescor LKJ(2) + tanh Jac
    - 0.5*sum(z_id(:).^2) - 0.5*sum(z_trl(:).^2);             % std_normal raw
end

function out = local_diffdynrel_sserr_rescor_draws(draws,m)
%Transform draws to the CmdStan-shaped arrays case 19 extracts: b, b_sigma (ND x2),
%b_dim1/b_dim2/b_sigma_dim1/b_sigma_dim2 (ND x KDIM), sd_id (ND x4), L_id (ND x4x4),
%id_varcov (ND x2x2 location block), trl_varcov (ND x2x2), rescor (ND x1),
%err_varcov (ND x2x2 at z=0), er_var_ss (ND x NSUB x 2).
ix = m.ix;
ND = size(draws,1); NSUB = m.NSUB;
out = struct();
out.b           = draws(:,ix.b);
out.b_sigma     = draws(:,ix.b_sigma);
out.b_dim1      = draws(:,ix.b_dim1);
out.b_dim2      = draws(:,ix.b_dim2);
out.b_sigma_dim1 = draws(:,ix.bs_dim1);
out.b_sigma_dim2 = draws(:,ix.bs_dim2);
sd_id = exp(draws(:,ix.lsd_id));      % ND x4
out.sd_id = sd_id;
sd_trl = exp(draws(:,ix.lsd_trl));
rho_trl = tanh(draws(:,ix.zrho_trl));
out.trl_varcov = local_varcov2(sd_trl,rho_trl);
rho_res = tanh(draws(:,ix.zrho_res));
out.rescor = rho_res;
sigma = exp(draws(:,ix.b_sigma));
out.err_varcov = local_varcov2(sigma,rho_res);   % at z = 0 (scale RE = 0)
L_id = zeros(ND,4,4);
id_varcov = zeros(ND,2,2);
er = zeros(ND,NSUB,2);
for d = 1:ND
    Ld = psyrat_corr_chol_lkj(draws(d,ix.yLid).',4,2);
    L_id(d,:,:) = Ld;
    cor = Ld*Ld.';
    sdd = sd_id(d,:);
    id_varcov(d,1,1) = sdd(1)^2;
    id_varcov(d,2,2) = sdd(2)^2;
    id_varcov(d,1,2) = sdd(1)*sdd(2)*cor(1,2);
    id_varcov(d,2,1) = id_varcov(d,1,2);
    z4 = reshape(draws(d,ix.z_id),4,NSUB);
    r1 = (sdd(:) .* (Ld * z4)).';
    er(d,:,1) = draws(d,ix.b_sigma(1)) + r1(:,3).';
    er(d,:,2) = draws(d,ix.b_sigma(2)) + r1(:,4).';
end
out.L_id = L_id;
out.id_varcov = id_varcov;
out.er_var_ss = er;
end

% ===========================================================================
% Design: diff_dynrel_sserr_trt_rescor (case 20) - concurrent two-facet
% subject-level dynamic difference (bivariate residual + 4-D person block)
% ===========================================================================
function fit = local_hmc_diffdynrel_sserr_trt_rescor(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the CONCURRENT two-facet SUBJECT-LEVEL
%dynamic difference model, matching psyrat_build_stan_diffdynrel_sserr_trt_rescor:
%the concurrent one-facet subject-level model (local_hmc_diffdynrel_sserr_rescor,
%case 19) with the single trial 2-D block generalized to the FIVE crossed
%cross-condition 2-D factors (trl/occ/tid/oid/to, as case 16); the person (id)
%factor stays the 4-D joint location+scale block. Same per-row bivariate residual
%(population Lrescor) + per-event dimension slopes.

m = struct();
m.NOBS  = double(data.NOBS);
m.NSUB  = double(data.NSUB);
m.meas1 = double(data.meas1(:));
m.meas2 = double(data.meas2(:));
m.id    = double(data.ID(:));
m.KDIM  = double(data.KDIM);
m.Xdim  = reshape(double(data.Xdim_vec(:)),m.NOBS,m.KDIM);
% the FIVE crossed location-only 2x2 factors (the id factor is the 4-D block)
m.factors = struct( ...
    'name',  {'trl','occ','tid','oid','to'}, ...
    'idx',   {double(data.TRL(:)),double(data.OCC(:)),double(data.TID(:)), ...
              double(data.OID(:)),double(data.TO(:))}, ...
    'N',     {double(data.NTRL),double(data.NOCC),double(data.NTID), ...
              double(data.NOID),double(data.NTO)}, ...
    'prior', {10,5,5,5,1}, ...
    'field', {'trl_varcov','occ_varcov','tid_varcov','oid_varcov','to_varcov'});

Pb = 5; Pbs = 10; Psd_id = 10; Psd_id_sigma = 2.5; Pdim = struct('b',10,'b_sigma',1);
priormap = {'sd_trl','sd_occ','sd_tid','sd_oid','sd_to'};
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors)
    pr = design.priors;
    if isfield(pr,'diff_trt') && isstruct(pr.diff_trt)
        if isfield(pr.diff_trt,'b');       Pb  = double(pr.diff_trt.b);       end
        if isfield(pr.diff_trt,'b_sigma'); Pbs = double(pr.diff_trt.b_sigma); end
        for k = 1:5
            if isfield(pr.diff_trt,priormap{k})
                m.factors(k).prior = double(pr.diff_trt.(priormap{k}));
            end
        end
    end
    if isfield(pr,'diff') && isstruct(pr.diff)
        if isfield(pr.diff,'sd_id');       Psd_id = double(pr.diff.sd_id);             end
        if isfield(pr.diff,'sd_id_sigma'); Psd_id_sigma = double(pr.diff.sd_id_sigma); end
    end
    if isfield(pr,'dynrel') && isstruct(pr.dynrel)
        if isfield(pr.dynrel,'b');       Pdim.b = double(pr.dynrel.b);             end
        if isfield(pr.dynrel,'b_sigma'); Pdim.b_sigma = double(pr.dynrel.b_sigma); end
    end
end
m.Pb = Pb; m.Pbs = Pbs; m.Psd_id = Psd_id; m.Psd_id_sigma = Psd_id_sigma; m.Pdim = Pdim;

m = local_diffdynrel_sserr_trt_rescor_index(m);

gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_diffdynrel_sserr_trt_rescor_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_diffdynrel_sserr_trt_rescor_logp(t,m), th);
end

start = local_diffdynrel_sserr_trt_rescor_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_diffdynrel_sserr_trt_rescor_draws(draws,m);
meta = struct('design','diff_dynrel_sserr_trt_rescor','gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_diffdynrel_sserr_trt_rescor_index(m)
%theta = [b(2); b_sigma(2); b_dim1(K); b_dim2(K); b_sigma_dim1(K); b_sigma_dim2(K);
%         log_sd_id(4); yLid(6);          % 4-D person block
%         per crossed factor f (5): log_sd_f(2), z_rho_f(1);
%         z_rho_res;
%         z_id(4*NSUB);
%         per crossed factor f: z_f(2*N_f)].
K = m.KDIM; NSUB = m.NSUB; i = 0;
m.ix.b        = i+(1:2);  i = i+2;
m.ix.b_sigma  = i+(1:2);  i = i+2;
m.ix.b_dim1   = i+(1:K);  i = i+K;
m.ix.b_dim2   = i+(1:K);  i = i+K;
m.ix.bs_dim1  = i+(1:K);  i = i+K;
m.ix.bs_dim2  = i+(1:K);  i = i+K;
m.ix.lsd_id   = i+(1:4);  i = i+4;
m.ix.yLid     = i+(1:6);  i = i+6;
m.ix.lsd  = cell(1,5); m.ix.zrho = cell(1,5);
for k = 1:5
    m.ix.lsd{k}  = i+(1:2); i = i+2;
    m.ix.zrho{k} = i+1;     i = i+1;
end
m.ix.zrho_res = i+1; i = i+1;
m.ix.z_id = i+(1:4*NSUB); i = i+4*NSUB;
m.ix.z = cell(1,5);
for k = 1:5
    m.ix.z{k} = i+(1:2*m.factors(k).N); i = i+2*m.factors(k).N;
end
m.D = i;
mon = [m.ix.b, m.ix.b_sigma, m.ix.b_dim1, m.ix.b_dim2, m.ix.bs_dim1, ...
    m.ix.bs_dim2, m.ix.lsd_id, m.ix.zrho_res];
names = [{'b1','b2','b_sigma1','b_sigma2'}, ...
    local_numbered('b_dim1',K), local_numbered('b_dim2',K), ...
    local_numbered('b_sigma_dim1',K), local_numbered('b_sigma_dim2',K), ...
    {'log_sd_id1','log_sd_id2','log_sd_id3','log_sd_id4','z_rho_res'}];
for k = 1:5
    mon = [mon, m.ix.lsd{k}, m.ix.zrho{k}]; %#ok<AGROW>
    names = [names, {sprintf('log_sd_%s1',m.factors(k).name), ...
        sprintf('log_sd_%s2',m.factors(k).name), ...
        sprintf('z_rho_%s',m.factors(k).name)}]; %#ok<AGROW>
end
m.mon_idx = mon; m.mon_names = names;
end

function start = local_diffdynrel_sserr_trt_rescor_start(m)
start = zeros(m.D,1);
start(m.ix.b(1)) = mean(m.meas1);
start(m.ix.b(2)) = mean(m.meas2);
s1 = std(m.meas1); s2 = std(m.meas2);
start(m.ix.b_sigma(1)) = log(max(s1/2,1e-3));
start(m.ix.b_sigma(2)) = log(max(s2/2,1e-3));
start(m.ix.lsd_id(1:2)) = log(max((s1+s2)/6,1e-3));
start(m.ix.lsd_id(3:4)) = log(0.3);
g = log(max((s1+s2)/8,1e-3));
for k = 1:5
    start(m.ix.lsd{k}) = g;
end
end

function lp = local_diffdynrel_sserr_trt_rescor_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching
%psyrat_build_stan_diffdynrel_sserr_trt_rescor: the 4-D person block (location+
%scale) + five crossed 2x2 location factors (concurrent) + per-event dimension
%slopes + per-row bivariate residual (the person scale REs feed the per-row sigma).
ix = m.ix;
b       = theta(ix.b);
b_sigma = theta(ix.b_sigma);
b_dim1  = theta(ix.b_dim1);  b_dim2  = theta(ix.b_dim2);
bs_dim1 = theta(ix.bs_dim1); bs_dim2 = theta(ix.bs_dim2);
% 4-D person block
sd_id = exp(theta(ix.lsd_id));
[L_id, lp_Lid] = psyrat_corr_chol_lkj(theta(ix.yLid), 4, 2);
z_id = reshape(theta(ix.z_id),4,m.NSUB);
r1 = (sd_id .* (L_id * z_id)).';
r1_1 = r1(:,1); r1_2 = r1(:,2); r1_s1 = r1(:,3); r1_s2 = r1(:,4);
rho_res = tanh(theta(ix.zrho_res));

mu1 = b(1) + m.Xdim*b_dim1 + r1_1(m.id);
mu2 = b(2) + m.Xdim*b_dim2 + r1_2(m.id);
lp_re = 0;
for k = 1:5
    f = m.factors(k);
    sd  = exp(theta(ix.lsd{k}));
    rho = tanh(theta(ix.zrho{k}));
    z   = reshape(theta(ix.z{k}),2,f.N);
    rA = sd(1).*z(1,:).';
    rB = sd(2).*(rho.*z(1,:).' + sqrt(1-rho.^2).*z(2,:).');
    mu1 = mu1 + rA(f.idx);
    mu2 = mu2 + rB(f.idx);
    lp_re = lp_re ...
        + local_student_t(sd(1),3,0,f.prior) + theta(ix.lsd{k}(1)) ...
        + local_student_t(sd(2),3,0,f.prior) + theta(ix.lsd{k}(2)) ...
        + 2*log(1 - rho.^2) ...
        - 0.5*sum(z(:).^2);
end
% per-row residual log-SDs (population + dim slopes + per-subject scale RE)
logsig1 = b_sigma(1) + m.Xdim*bs_dim1 + r1_s1(m.id);
logsig2 = b_sigma(2) + m.Xdim*bs_dim2 + r1_s2(m.id);
sig1 = exp(logsig1); sig2 = exp(logsig2);
d1 = m.meas1 - mu1; d2 = m.meas2 - mu2;
w1 = d1./sig1;
w2 = (d2 - (sig2.*rho_res./sig1).*d1) ./ (sig2.*sqrt(1-rho_res.^2));
ll = -sum(logsig1) - sum(logsig2) - 0.5*m.NOBS*log(1-rho_res.^2) ...
   - 0.5*sum(w1.^2 + w2.^2);

Pdim = m.Pdim;
lp = ll + lp_re + lp_Lid ...
    + sum(local_normal(b,0,m.Pb)) ...
    + sum(local_student_t(b_sigma,3,0,m.Pbs)) ...
    + sum(local_normal(b_dim1,0,Pdim.b)) + sum(local_normal(b_dim2,0,Pdim.b)) ...
    + sum(local_normal(bs_dim1,0,Pdim.b_sigma)) + sum(local_normal(bs_dim2,0,Pdim.b_sigma)) ...
    + local_student_t(sd_id(1),3,0,m.Psd_id)       + theta(ix.lsd_id(1)) ...
    + local_student_t(sd_id(2),3,0,m.Psd_id)       + theta(ix.lsd_id(2)) ...
    + local_student_t(sd_id(3),3,0,m.Psd_id_sigma) + theta(ix.lsd_id(3)) ...
    + local_student_t(sd_id(4),3,0,m.Psd_id_sigma) + theta(ix.lsd_id(4)) ...
    + 2*log(1 - rho_res.^2) ...                                  % Lrescor LKJ(2) + tanh Jac
    - 0.5*sum(z_id(:).^2);                                       % std_normal person raw
end

function out = local_diffdynrel_sserr_trt_rescor_draws(draws,m)
%Transform draws to the CmdStan-shaped arrays case 20 extracts: b, b_sigma (ND x2),
%b_dim1/b_dim2/b_sigma_dim1/b_sigma_dim2 (ND x KDIM), sd_id (ND x4), L_id (ND x4x4),
%id_varcov (ND x2x2 location block), the five crossed 2x2 var-covariances, rescor
%(ND x1), err_varcov (ND x2x2 at z=0), er_var_ss (ND x NSUB x 2).
ix = m.ix;
ND = size(draws,1); NSUB = m.NSUB;
out = struct();
out.b           = draws(:,ix.b);
out.b_sigma     = draws(:,ix.b_sigma);
out.b_dim1      = draws(:,ix.b_dim1);
out.b_dim2      = draws(:,ix.b_dim2);
out.b_sigma_dim1 = draws(:,ix.bs_dim1);
out.b_sigma_dim2 = draws(:,ix.bs_dim2);
sd_id = exp(draws(:,ix.lsd_id));
out.sd_id = sd_id;
rho_res = tanh(draws(:,ix.zrho_res));
out.rescor = rho_res;
sigma = exp(draws(:,ix.b_sigma));
out.err_varcov = local_varcov2(sigma,rho_res);
for k = 1:5
    sd  = exp(draws(:,ix.lsd{k}));
    rho = tanh(draws(:,ix.zrho{k}));
    out.(m.factors(k).field) = local_varcov2(sd,rho);
end
L_id = zeros(ND,4,4); id_varcov = zeros(ND,2,2); er = zeros(ND,NSUB,2);
for d = 1:ND
    Ld = psyrat_corr_chol_lkj(draws(d,ix.yLid).',4,2);
    L_id(d,:,:) = Ld;
    cor = Ld*Ld.';
    sdd = sd_id(d,:);
    id_varcov(d,1,1) = sdd(1)^2;
    id_varcov(d,2,2) = sdd(2)^2;
    id_varcov(d,1,2) = sdd(1)*sdd(2)*cor(1,2);
    id_varcov(d,2,1) = id_varcov(d,1,2);
    z4 = reshape(draws(d,ix.z_id),4,NSUB);
    r1 = (sdd(:) .* (Ld * z4)).';
    er(d,:,1) = draws(d,ix.b_sigma(1)) + r1(:,3).';
    er(d,:,2) = draws(d,ix.b_sigma(2)) + r1(:,4).';
end
out.L_id = L_id;
out.id_varcov = id_varcov;
out.er_var_ss = er;
end

% ===========================================================================
% Design: diff_dynrel_sserr_rho (case 21) - concurrent one-facet subject-level
% dynamic difference with a PER-SUBJECT residual correlation
% ===========================================================================
function fit = local_hmc_diffdynrel_sserr_rho(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the concurrent one-facet subject-level
%dynamic difference model with a PER-SUBJECT residual correlation, matching
%psyrat_build_stan_diffdynrel_sserr_rho: the case-19 model
%(local_hmc_diffdynrel_sserr_rescor) with the single population residual
%correlation replaced by a hierarchical per-participant correlation
%rho[s] = tanh(rescor_mu + sd_rescor * z_rescor[s]) (the case-8 rho hierarchy), so
%the per-row bivariate residual uses rho[ID[n]].

m = struct();
m.NOBS  = double(data.NOBS);
m.NSUB  = double(data.NSUB);
m.NTRL  = double(data.NTRL);
m.meas1 = double(data.meas1(:));
m.meas2 = double(data.meas2(:));
m.id    = double(data.ID(:));
m.trl   = double(data.TRL(:));
m.KDIM  = double(data.KDIM);
m.Xdim  = reshape(double(data.Xdim_vec(:)),m.NOBS,m.KDIM);

% configurable priors: base diff block from priors.diff (incl. rescor_mu /
% sd_rescor for the per-subject correlation hierarchy); slopes from priors.dynrel.
P = struct('b',5,'b_sigma',10,'sd_id',10,'sd_trl',10,'sd_id_sigma',2.5, ...
    'rescor_mu',1,'sd_rescor',1);
Pdim = struct('b',10,'b_sigma',1);
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors)
    pr = design.priors;
    if isfield(pr,'diff') && isstruct(pr.diff)
        fn = {'b','b_sigma','sd_id','sd_trl','sd_id_sigma','rescor_mu','sd_rescor'};
        for k = 1:numel(fn)
            if isfield(pr.diff,fn{k}); P.(fn{k}) = double(pr.diff.(fn{k})); end
        end
    end
    if isfield(pr,'dynrel') && isstruct(pr.dynrel)
        if isfield(pr.dynrel,'b');       Pdim.b = double(pr.dynrel.b);             end
        if isfield(pr.dynrel,'b_sigma'); Pdim.b_sigma = double(pr.dynrel.b_sigma); end
    end
end
m.P = P; m.Pdim = Pdim;

m = local_diffdynrel_sserr_rho_index(m);

gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_diffdynrel_sserr_rho_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_diffdynrel_sserr_rho_logp(t,m), th);
end

start = local_diffdynrel_sserr_rho_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_diffdynrel_sserr_rho_draws(draws,m);
meta = struct('design','diff_dynrel_sserr_rho','gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_diffdynrel_sserr_rho_index(m)
%theta = [b(2); b_sigma(2); b_dim1(K); b_dim2(K); b_sigma_dim1(K); b_sigma_dim2(K);
%         log_sd_id(4); yLid(6); log_sd_trl(2); z_rho_trl;
%         rescor_mu; log_sd_rescor; z_rescor(NSUB);   % per-subject rho hierarchy
%         z_id(4*NSUB); z_trl(2*NTRL)].
K = m.KDIM; NSUB = m.NSUB; NTRL = m.NTRL; i = 0;
m.ix.b        = i+(1:2);  i = i+2;
m.ix.b_sigma  = i+(1:2);  i = i+2;
m.ix.b_dim1   = i+(1:K);  i = i+K;
m.ix.b_dim2   = i+(1:K);  i = i+K;
m.ix.bs_dim1  = i+(1:K);  i = i+K;
m.ix.bs_dim2  = i+(1:K);  i = i+K;
m.ix.lsd_id   = i+(1:4);  i = i+4;
m.ix.yLid     = i+(1:6);  i = i+6;
m.ix.lsd_trl  = i+(1:2);  i = i+2;
m.ix.zrho_trl = i+1;      i = i+1;
m.ix.rescor_mu  = i+1;        i = i+1;
m.ix.lsd_rescor = i+1;        i = i+1;
m.ix.z_rescor   = i+(1:NSUB); i = i+NSUB;
m.ix.z_id     = i+(1:4*NSUB); i = i+4*NSUB;
m.ix.z_trl    = i+(1:2*NTRL); i = i+2*NTRL;
m.D = i;
m.mon_idx = [m.ix.b, m.ix.b_sigma, m.ix.b_dim1, m.ix.b_dim2, ...
    m.ix.bs_dim1, m.ix.bs_dim2, m.ix.lsd_id, m.ix.lsd_trl, ...
    m.ix.zrho_trl, m.ix.rescor_mu, m.ix.lsd_rescor];
m.mon_names = [{'b1','b2','b_sigma1','b_sigma2'}, ...
    local_numbered('b_dim1',K), local_numbered('b_dim2',K), ...
    local_numbered('b_sigma_dim1',K), local_numbered('b_sigma_dim2',K), ...
    {'log_sd_id1','log_sd_id2','log_sd_id3','log_sd_id4', ...
     'log_sd_trl1','log_sd_trl2','z_rho_trl','rescor_mu','log_sd_rescor'}];
end

function start = local_diffdynrel_sserr_rho_start(m)
start = zeros(m.D,1);
start(m.ix.b(1)) = mean(m.meas1);
start(m.ix.b(2)) = mean(m.meas2);
s1 = std(m.meas1); s2 = std(m.meas2);
start(m.ix.b_sigma(1)) = log(max(s1/2,1e-3));
start(m.ix.b_sigma(2)) = log(max(s2/2,1e-3));
start(m.ix.lsd_id(1:2)) = log(max((s1+s2)/6,1e-3));
start(m.ix.lsd_id(3:4)) = log(0.3);
start(m.ix.lsd_trl)     = log(max((s1+s2)/4,1e-3));
start(m.ix.lsd_rescor)  = log(0.3);
end

function lp = local_diffdynrel_sserr_rho_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching
%psyrat_build_stan_diffdynrel_sserr_rho: case-19 model with the population residual
%correlation replaced by a per-subject rho[s] = tanh(rescor_mu + sd_rescor*z_rescor[s]).
ix = m.ix;
b       = theta(ix.b);
b_sigma = theta(ix.b_sigma);
b_dim1  = theta(ix.b_dim1);  b_dim2  = theta(ix.b_dim2);
bs_dim1 = theta(ix.bs_dim1); bs_dim2 = theta(ix.bs_dim2);
sd_id = exp(theta(ix.lsd_id));
[L_id, lp_Lid] = psyrat_corr_chol_lkj(theta(ix.yLid), 4, 2);
z_id = reshape(theta(ix.z_id),4,m.NSUB);
r1 = (sd_id .* (L_id * z_id)).';
r1_1 = r1(:,1); r1_2 = r1(:,2); r1_s1 = r1(:,3); r1_s2 = r1(:,4);
sd_trl  = exp(theta(ix.lsd_trl));
rho_trl = tanh(theta(ix.zrho_trl));
z_trl = reshape(theta(ix.z_trl),2,m.NTRL);
r2_1 = sd_trl(1).*z_trl(1,:).';
r2_2 = sd_trl(2).*(rho_trl.*z_trl(1,:).' + sqrt(1-rho_trl.^2).*z_trl(2,:).');
% per-subject residual correlation hierarchy
rescor_mu = theta(ix.rescor_mu);
sd_rescor = exp(theta(ix.lsd_rescor));
z_rescor  = theta(ix.z_rescor);
rho_s = tanh(rescor_mu + sd_rescor.*z_rescor);   % NSUB x 1
rho_n = rho_s(m.id);                             % per observation

mu1 = b(1) + m.Xdim*b_dim1 + r1_1(m.id) + r2_1(m.trl);
mu2 = b(2) + m.Xdim*b_dim2 + r1_2(m.id) + r2_2(m.trl);
logsig1 = b_sigma(1) + m.Xdim*bs_dim1 + r1_s1(m.id);
logsig2 = b_sigma(2) + m.Xdim*bs_dim2 + r1_s2(m.id);
sig1 = exp(logsig1); sig2 = exp(logsig2);
d1 = m.meas1 - mu1; d2 = m.meas2 - mu2;
% per-row bivariate residual with the per-subject rho_n (logdet term per row)
w1 = d1./sig1;
w2 = (d2 - (sig2.*rho_n./sig1).*d1) ./ (sig2.*sqrt(1-rho_n.^2));
ll = -sum(logsig1) - sum(logsig2) - 0.5*sum(log(1-rho_n.^2)) ...
   - 0.5*sum(w1.^2 + w2.^2);

P = m.P; Pdim = m.Pdim;
lp = ll + lp_Lid ...
    + sum(local_normal(b,0,P.b)) ...
    + sum(local_student_t(b_sigma,3,0,P.b_sigma)) ...
    + sum(local_normal(b_dim1,0,Pdim.b)) + sum(local_normal(b_dim2,0,Pdim.b)) ...
    + sum(local_normal(bs_dim1,0,Pdim.b_sigma)) + sum(local_normal(bs_dim2,0,Pdim.b_sigma)) ...
    + local_student_t(sd_id(1),3,0,P.sd_id)       + theta(ix.lsd_id(1)) ...
    + local_student_t(sd_id(2),3,0,P.sd_id)       + theta(ix.lsd_id(2)) ...
    + local_student_t(sd_id(3),3,0,P.sd_id_sigma) + theta(ix.lsd_id(3)) ...
    + local_student_t(sd_id(4),3,0,P.sd_id_sigma) + theta(ix.lsd_id(4)) ...
    + local_student_t(sd_trl(1),3,0,P.sd_trl) + theta(ix.lsd_trl(1)) ...
    + local_student_t(sd_trl(2),3,0,P.sd_trl) + theta(ix.lsd_trl(2)) ...
    + 2*log(1 - rho_trl.^2) ...                                  % L_trl LKJ(2) + tanh Jac
    + local_normal(rescor_mu,0,P.rescor_mu) ...                  % rescor_mu ~ normal(0,.)
    + local_student_t(sd_rescor,3,0,P.sd_rescor) + theta(ix.lsd_rescor) ... % half-t + Jac
    - 0.5*sum(z_rescor.^2) ...                                   % std_normal z_rescor
    - 0.5*sum(z_id(:).^2) - 0.5*sum(z_trl(:).^2);                % std_normal raw
end

function out = local_diffdynrel_sserr_rho_draws(draws,m)
%Transform draws to the CmdStan-shaped arrays case 21 extracts: b, b_sigma (ND x2),
%b_dim1/b_dim2/b_sigma_dim1/b_sigma_dim2 (ND x KDIM), sd_id (ND x4), L_id (ND x4x4),
%id_varcov/trl_varcov (ND x2x2), rho (ND x NSUB, the per-subject correlations),
%rescor_mu (ND x1), sd_rescor (ND x1), rescor (ND x1 = mean(rho)), err_varcov
%(ND x2x2, typical-person at z=0 using mean(rho)), er_var_ss (ND x NSUB x 2).
ix = m.ix;
ND = size(draws,1); NSUB = m.NSUB;
out = struct();
out.b           = draws(:,ix.b);
out.b_sigma     = draws(:,ix.b_sigma);
out.b_dim1      = draws(:,ix.b_dim1);
out.b_dim2      = draws(:,ix.b_dim2);
out.b_sigma_dim1 = draws(:,ix.bs_dim1);
out.b_sigma_dim2 = draws(:,ix.bs_dim2);
sd_id = exp(draws(:,ix.lsd_id));
out.sd_id = sd_id;
sd_trl = exp(draws(:,ix.lsd_trl));
rho_trl = tanh(draws(:,ix.zrho_trl));
out.trl_varcov = local_varcov2(sd_trl,rho_trl);
out.rescor_mu = draws(:,ix.rescor_mu);
out.sd_rescor = exp(draws(:,ix.lsd_rescor));
% per-subject correlations rho(d,s) = tanh(rescor_mu + sd_rescor*z_rescor[s])
rho = tanh(draws(:,ix.rescor_mu) + out.sd_rescor .* draws(:,ix.z_rescor));  % ND x NSUB
out.rho = rho;
mean_rho = mean(rho,2);                          % ND x1
out.rescor = mean_rho;
sigma = exp(draws(:,ix.b_sigma));
out.err_varcov = local_varcov2(sigma,mean_rho);  % typical person at z=0, mean(rho)
% person 4x4 -> L_id, id_varcov (location block), er_var_ss
L_id = zeros(ND,4,4); id_varcov = zeros(ND,2,2); er = zeros(ND,NSUB,2);
for d = 1:ND
    Ld = psyrat_corr_chol_lkj(draws(d,ix.yLid).',4,2);
    L_id(d,:,:) = Ld;
    cor = Ld*Ld.';
    sdd = sd_id(d,:);
    id_varcov(d,1,1) = sdd(1)^2;
    id_varcov(d,2,2) = sdd(2)^2;
    id_varcov(d,1,2) = sdd(1)*sdd(2)*cor(1,2);
    id_varcov(d,2,1) = id_varcov(d,1,2);
    z4 = reshape(draws(d,ix.z_id),4,NSUB);
    r1 = (sdd(:) .* (Ld * z4)).';
    er(d,:,1) = draws(d,ix.b_sigma(1)) + r1(:,3).';
    er(d,:,2) = draws(d,ix.b_sigma(2)) + r1(:,4).';
end
out.L_id = L_id;
out.id_varcov = id_varcov;
out.er_var_ss = er;
end

% ===========================================================================
% Design: diff_dynrel_sserr_trt_rho (case 22) - concurrent two-facet
% subject-level dynamic difference with a PER-SUBJECT residual correlation
% ===========================================================================
function fit = local_hmc_diffdynrel_sserr_trt_rho(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the concurrent two-facet subject-level
%dynamic difference model with a PER-SUBJECT residual correlation, matching
%psyrat_build_stan_diffdynrel_sserr_trt_rho: the case-20 model
%(local_hmc_diffdynrel_sserr_trt_rescor) with the single population residual
%correlation replaced by the case-21 hierarchical per-participant correlation
%rho[s] = tanh(rescor_mu + sd_rescor * z_rescor[s]).

m = struct();
m.NOBS  = double(data.NOBS);
m.NSUB  = double(data.NSUB);
m.meas1 = double(data.meas1(:));
m.meas2 = double(data.meas2(:));
m.id    = double(data.ID(:));
m.KDIM  = double(data.KDIM);
m.Xdim  = reshape(double(data.Xdim_vec(:)),m.NOBS,m.KDIM);
m.factors = struct( ...
    'name',  {'trl','occ','tid','oid','to'}, ...
    'idx',   {double(data.TRL(:)),double(data.OCC(:)),double(data.TID(:)), ...
              double(data.OID(:)),double(data.TO(:))}, ...
    'N',     {double(data.NTRL),double(data.NOCC),double(data.NTID), ...
              double(data.NOID),double(data.NTO)}, ...
    'prior', {10,5,5,5,1}, ...
    'field', {'trl_varcov','occ_varcov','tid_varcov','oid_varcov','to_varcov'});

Pb = 5; Pbs = 10; Psd_id = 10; Psd_id_sigma = 2.5; Prescor_mu = 1; Psd_rescor = 1;
Pdim = struct('b',10,'b_sigma',1);
priormap = {'sd_trl','sd_occ','sd_tid','sd_oid','sd_to'};
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors)
    pr = design.priors;
    if isfield(pr,'diff_trt') && isstruct(pr.diff_trt)
        if isfield(pr.diff_trt,'b');       Pb  = double(pr.diff_trt.b);       end
        if isfield(pr.diff_trt,'b_sigma'); Pbs = double(pr.diff_trt.b_sigma); end
        for k = 1:5
            if isfield(pr.diff_trt,priormap{k})
                m.factors(k).prior = double(pr.diff_trt.(priormap{k}));
            end
        end
    end
    if isfield(pr,'diff') && isstruct(pr.diff)
        if isfield(pr.diff,'sd_id');       Psd_id = double(pr.diff.sd_id);             end
        if isfield(pr.diff,'sd_id_sigma'); Psd_id_sigma = double(pr.diff.sd_id_sigma); end
        if isfield(pr.diff,'rescor_mu');   Prescor_mu = double(pr.diff.rescor_mu);     end
        if isfield(pr.diff,'sd_rescor');   Psd_rescor = double(pr.diff.sd_rescor);     end
    end
    if isfield(pr,'dynrel') && isstruct(pr.dynrel)
        if isfield(pr.dynrel,'b');       Pdim.b = double(pr.dynrel.b);             end
        if isfield(pr.dynrel,'b_sigma'); Pdim.b_sigma = double(pr.dynrel.b_sigma); end
    end
end
m.Pb = Pb; m.Pbs = Pbs; m.Psd_id = Psd_id; m.Psd_id_sigma = Psd_id_sigma;
m.Prescor_mu = Prescor_mu; m.Psd_rescor = Psd_rescor; m.Pdim = Pdim;

m = local_diffdynrel_sserr_trt_rho_index(m);

gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_diffdynrel_sserr_trt_rho_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_diffdynrel_sserr_trt_rho_logp(t,m), th);
end

start = local_diffdynrel_sserr_trt_rho_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_diffdynrel_sserr_trt_rho_draws(draws,m);
meta = struct('design','diff_dynrel_sserr_trt_rho','gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_diffdynrel_sserr_trt_rho_index(m)
%theta = case 20's layout but with z_rho_res replaced by the case-21 per-subject
%rho hierarchy (rescor_mu; log_sd_rescor; z_rescor(NSUB)).
K = m.KDIM; NSUB = m.NSUB; i = 0;
m.ix.b        = i+(1:2);  i = i+2;
m.ix.b_sigma  = i+(1:2);  i = i+2;
m.ix.b_dim1   = i+(1:K);  i = i+K;
m.ix.b_dim2   = i+(1:K);  i = i+K;
m.ix.bs_dim1  = i+(1:K);  i = i+K;
m.ix.bs_dim2  = i+(1:K);  i = i+K;
m.ix.lsd_id   = i+(1:4);  i = i+4;
m.ix.yLid     = i+(1:6);  i = i+6;
m.ix.lsd  = cell(1,5); m.ix.zrho = cell(1,5);
for k = 1:5
    m.ix.lsd{k}  = i+(1:2); i = i+2;
    m.ix.zrho{k} = i+1;     i = i+1;
end
m.ix.rescor_mu  = i+1;        i = i+1;
m.ix.lsd_rescor = i+1;        i = i+1;
m.ix.z_rescor   = i+(1:NSUB); i = i+NSUB;
m.ix.z_id = i+(1:4*NSUB); i = i+4*NSUB;
m.ix.z = cell(1,5);
for k = 1:5
    m.ix.z{k} = i+(1:2*m.factors(k).N); i = i+2*m.factors(k).N;
end
m.D = i;
mon = [m.ix.b, m.ix.b_sigma, m.ix.b_dim1, m.ix.b_dim2, m.ix.bs_dim1, ...
    m.ix.bs_dim2, m.ix.lsd_id, m.ix.rescor_mu, m.ix.lsd_rescor];
names = [{'b1','b2','b_sigma1','b_sigma2'}, ...
    local_numbered('b_dim1',K), local_numbered('b_dim2',K), ...
    local_numbered('b_sigma_dim1',K), local_numbered('b_sigma_dim2',K), ...
    {'log_sd_id1','log_sd_id2','log_sd_id3','log_sd_id4','rescor_mu','log_sd_rescor'}];
for k = 1:5
    mon = [mon, m.ix.lsd{k}, m.ix.zrho{k}]; %#ok<AGROW>
    names = [names, {sprintf('log_sd_%s1',m.factors(k).name), ...
        sprintf('log_sd_%s2',m.factors(k).name), ...
        sprintf('z_rho_%s',m.factors(k).name)}]; %#ok<AGROW>
end
m.mon_idx = mon; m.mon_names = names;
end

function start = local_diffdynrel_sserr_trt_rho_start(m)
start = zeros(m.D,1);
start(m.ix.b(1)) = mean(m.meas1);
start(m.ix.b(2)) = mean(m.meas2);
s1 = std(m.meas1); s2 = std(m.meas2);
start(m.ix.b_sigma(1)) = log(max(s1/2,1e-3));
start(m.ix.b_sigma(2)) = log(max(s2/2,1e-3));
start(m.ix.lsd_id(1:2)) = log(max((s1+s2)/6,1e-3));
start(m.ix.lsd_id(3:4)) = log(0.3);
start(m.ix.lsd_rescor)  = log(0.3);
g = log(max((s1+s2)/8,1e-3));
for k = 1:5
    start(m.ix.lsd{k}) = g;
end
end

function lp = local_diffdynrel_sserr_trt_rho_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching
%psyrat_build_stan_diffdynrel_sserr_trt_rho: case-20 model with the population
%residual correlation replaced by the per-subject rho[s].
ix = m.ix;
b       = theta(ix.b);
b_sigma = theta(ix.b_sigma);
b_dim1  = theta(ix.b_dim1);  b_dim2  = theta(ix.b_dim2);
bs_dim1 = theta(ix.bs_dim1); bs_dim2 = theta(ix.bs_dim2);
sd_id = exp(theta(ix.lsd_id));
[L_id, lp_Lid] = psyrat_corr_chol_lkj(theta(ix.yLid), 4, 2);
z_id = reshape(theta(ix.z_id),4,m.NSUB);
r1 = (sd_id .* (L_id * z_id)).';
r1_1 = r1(:,1); r1_2 = r1(:,2); r1_s1 = r1(:,3); r1_s2 = r1(:,4);
% per-subject residual correlation hierarchy
rescor_mu = theta(ix.rescor_mu);
sd_rescor = exp(theta(ix.lsd_rescor));
z_rescor  = theta(ix.z_rescor);
rho_s = tanh(rescor_mu + sd_rescor.*z_rescor);
rho_n = rho_s(m.id);

mu1 = b(1) + m.Xdim*b_dim1 + r1_1(m.id);
mu2 = b(2) + m.Xdim*b_dim2 + r1_2(m.id);
lp_re = 0;
for k = 1:5
    f = m.factors(k);
    sd  = exp(theta(ix.lsd{k}));
    rho = tanh(theta(ix.zrho{k}));
    z   = reshape(theta(ix.z{k}),2,f.N);
    rA = sd(1).*z(1,:).';
    rB = sd(2).*(rho.*z(1,:).' + sqrt(1-rho.^2).*z(2,:).');
    mu1 = mu1 + rA(f.idx);
    mu2 = mu2 + rB(f.idx);
    lp_re = lp_re ...
        + local_student_t(sd(1),3,0,f.prior) + theta(ix.lsd{k}(1)) ...
        + local_student_t(sd(2),3,0,f.prior) + theta(ix.lsd{k}(2)) ...
        + 2*log(1 - rho.^2) ...
        - 0.5*sum(z(:).^2);
end
logsig1 = b_sigma(1) + m.Xdim*bs_dim1 + r1_s1(m.id);
logsig2 = b_sigma(2) + m.Xdim*bs_dim2 + r1_s2(m.id);
sig1 = exp(logsig1); sig2 = exp(logsig2);
d1 = m.meas1 - mu1; d2 = m.meas2 - mu2;
w1 = d1./sig1;
w2 = (d2 - (sig2.*rho_n./sig1).*d1) ./ (sig2.*sqrt(1-rho_n.^2));
ll = -sum(logsig1) - sum(logsig2) - 0.5*sum(log(1-rho_n.^2)) ...
   - 0.5*sum(w1.^2 + w2.^2);

lp = ll + lp_re + lp_Lid ...
    + sum(local_normal(b,0,m.Pb)) ...
    + sum(local_student_t(b_sigma,3,0,m.Pbs)) ...
    + sum(local_normal(b_dim1,0,m.Pdim.b)) + sum(local_normal(b_dim2,0,m.Pdim.b)) ...
    + sum(local_normal(bs_dim1,0,m.Pdim.b_sigma)) + sum(local_normal(bs_dim2,0,m.Pdim.b_sigma)) ...
    + local_student_t(sd_id(1),3,0,m.Psd_id)       + theta(ix.lsd_id(1)) ...
    + local_student_t(sd_id(2),3,0,m.Psd_id)       + theta(ix.lsd_id(2)) ...
    + local_student_t(sd_id(3),3,0,m.Psd_id_sigma) + theta(ix.lsd_id(3)) ...
    + local_student_t(sd_id(4),3,0,m.Psd_id_sigma) + theta(ix.lsd_id(4)) ...
    + local_normal(rescor_mu,0,m.Prescor_mu) ...                 % rescor_mu ~ normal(0,.)
    + local_student_t(sd_rescor,3,0,m.Psd_rescor) + theta(ix.lsd_rescor) ... % half-t + Jac
    - 0.5*sum(z_rescor.^2) ...                                   % std_normal z_rescor
    - 0.5*sum(z_id(:).^2);                                       % std_normal person raw
end

function out = local_diffdynrel_sserr_trt_rho_draws(draws,m)
%Transform draws to the CmdStan-shaped arrays case 22 extracts: case-20 draws plus
%the case-21 per-subject rho hierarchy (rho, rescor_mu, sd_rescor, rescor, and the
%mean(rho) err_varcov).
ix = m.ix;
ND = size(draws,1); NSUB = m.NSUB;
out = struct();
out.b           = draws(:,ix.b);
out.b_sigma     = draws(:,ix.b_sigma);
out.b_dim1      = draws(:,ix.b_dim1);
out.b_dim2      = draws(:,ix.b_dim2);
out.b_sigma_dim1 = draws(:,ix.bs_dim1);
out.b_sigma_dim2 = draws(:,ix.bs_dim2);
sd_id = exp(draws(:,ix.lsd_id));
out.sd_id = sd_id;
out.rescor_mu = draws(:,ix.rescor_mu);
out.sd_rescor = exp(draws(:,ix.lsd_rescor));
rho = tanh(draws(:,ix.rescor_mu) + out.sd_rescor .* draws(:,ix.z_rescor));   % ND x NSUB
out.rho = rho;
mean_rho = mean(rho,2);
out.rescor = mean_rho;
sigma = exp(draws(:,ix.b_sigma));
out.err_varcov = local_varcov2(sigma,mean_rho);
for k = 1:5
    sd  = exp(draws(:,ix.lsd{k}));
    rh  = tanh(draws(:,ix.zrho{k}));
    out.(m.factors(k).field) = local_varcov2(sd,rh);
end
L_id = zeros(ND,4,4); id_varcov = zeros(ND,2,2); er = zeros(ND,NSUB,2);
for d = 1:ND
    Ld = psyrat_corr_chol_lkj(draws(d,ix.yLid).',4,2);
    L_id(d,:,:) = Ld;
    cor = Ld*Ld.';
    sdd = sd_id(d,:);
    id_varcov(d,1,1) = sdd(1)^2;
    id_varcov(d,2,2) = sdd(2)^2;
    id_varcov(d,1,2) = sdd(1)*sdd(2)*cor(1,2);
    id_varcov(d,2,1) = id_varcov(d,1,2);
    z4 = reshape(draws(d,ix.z_id),4,NSUB);
    r1 = (sdd(:) .* (Ld * z4)).';
    er(d,:,1) = draws(d,ix.b_sigma(1)) + r1(:,3).';
    er(d,:,2) = draws(d,ix.b_sigma(2)) + r1(:,4).';
end
out.L_id = L_id;
out.id_varcov = id_varcov;
out.er_var_ss = er;
end

% ===========================================================================
% Design: diff_sserr (case 8) - static subject-level two-event difference
% ===========================================================================
function fit = local_hmc_diff_sserr(design,data,seed,nchains,nwarm,nsamp,opts)
%Build and run the native HMC fit for the static subject-level two-event
%difference model, matching psyrat_build_stan_diff_sserr: three separate 2-D
%blocks - person location (L_id -> r_id), trial (L_trl -> r_trl), and a per-subject
%SCALE block (L_id_sigma -> r_sigma) feeding the per-event residual log-SD - plus a
%per-subject residual correlation rho[s] = tanh(rescor_mu + sd_rescor*z_rescor[s])
%(toggled by estimate_rescor). The likelihood is a bivariate multi_normal_cholesky
%over the PAIRED observations and a univariate normal over the UNPAIRED single-event
%observations. NOT dynamic (no dimension slopes).

m = struct();
m.NSUB = double(data.NSUB);
m.NTRL = double(data.NTRL);
m.NPAIR = double(data.NPAIR);
m.NUNP  = double(data.NUNP);
m.estimate_rescor = double(data.estimate_rescor) == 1;
% the data are PADDED to NPAIR_SAFE / NUNP_SAFE; use the first NPAIR / NUNP rows.
np = m.NPAIR; nu = m.NUNP;
if np > 0
    m.meas1 = double(data.meas1(1:np)); m.meas2 = double(data.meas2(1:np));
    m.id_pair = double(data.ID_PAIR(1:np)); m.trl_pair = double(data.TRL_PAIR(1:np));
else
    m.meas1 = []; m.meas2 = []; m.id_pair = []; m.trl_pair = [];
end
if nu > 0
    m.meas_unp = double(data.meas_unp(1:nu));
    m.event_unp = double(data.EVENT_UNP(1:nu));
    m.id_unp = double(data.ID_UNP(1:nu));
    m.trl_unp = double(data.TRL_UNP(1:nu));
    m.e1 = double(m.event_unp == 1); m.e2 = double(m.event_unp == 2);
else
    m.meas_unp = []; m.event_unp = []; m.id_unp = []; m.trl_unp = [];
    m.e1 = []; m.e2 = [];
end

% configurable priors (all from priors.diff)
P = struct('b',5,'b_sigma',10,'sd_id',10,'sd_trl',10,'sd_id_sigma',2.5, ...
    'rescor_mu',1,'sd_rescor',1);
if isstruct(design) && isfield(design,'priors') && isstruct(design.priors)
    pr = design.priors;
    fn = fieldnames(P);
    for k = 1:numel(fn)
        if isfield(pr,fn{k}); P.(fn{k}) = double(pr.(fn{k})); end
    end
end
m.P = P;

m = local_diff_sserr_index(m);

gradmethod = 'auto';
if isstruct(opts) && isfield(opts,'hmc_gradient') && ~isempty(opts.hmc_gradient)
    gradmethod = lower(char(opts.hmc_gradient));
end
have_ad = (exist('dlarray','file') > 0) && (exist('dlfeval','file') > 0);
switch gradmethod
    case 'ad'
        if ~have_ad
            error('PsyRAT:hmcNoAD',['The native HMC engine was asked for ' ...
                'automatic-differentiation gradients (''ad''), but the Deep ' ...
                'Learning Toolbox (dlarray) is not available. How to fix: use ' ...
                'hmc_gradient ''numer'' (slower) or the fitlme/CmdStan engines.']);
        end
        use_ad = true;
    case 'numer'
        use_ad = false;
    otherwise % 'auto'
        use_ad = have_ad;
end
if use_ad
    logpfun = @(th) local_ad_logpgrad(@(t) local_diff_sserr_logp(t,m), th);
else
    logpfun = @(th) local_numer_logpgrad(@(t) local_diff_sserr_logp(t,m), th);
end

start = local_diff_sserr_start(m);
local_hmc_gradcheck(logpfun,start,m,use_ad);

[draws,conv] = psyrat_hmc_sample(logpfun,start,nchains,nwarm,nsamp,seed, ...
    m.mon_idx,m.mon_names);

drawstruct = local_diff_sserr_draws(draws,m);
meta = struct('design','diff_sserr','gradient',ternary(use_ad,'ad','numer'), ...
    'nchains',nchains,'nwarmup',nwarm,'nsampling',nsamp);
fit = PsyRATNativeFit(drawstruct,conv,'hmc',meta);
end

% ---------------------------------------------------------------------------
function m = local_diff_sserr_index(m)
%theta = [b(2); b_sigma(2);
%         log_sd_id(2); z_rho_id;              % person location 2-D block
%         log_sd_trl(2); z_rho_trl;            % trial 2-D block
%         log_sd_id_sigma(2); z_rho_id_sigma;  % person scale 2-D block
%         rescor_mu; log_sd_rescor; z_rescor(NSUB);  % per-subject rho hierarchy
%         z_id(2*NSUB); z_trl(2*NTRL); z_id_sigma(2*NSUB)].
NSUB = m.NSUB; NTRL = m.NTRL; i = 0;
m.ix.b         = i+(1:2);  i = i+2;
m.ix.b_sigma   = i+(1:2);  i = i+2;
m.ix.lsd_id    = i+(1:2);  i = i+2;
m.ix.zrho_id   = i+1;      i = i+1;
m.ix.lsd_trl   = i+(1:2);  i = i+2;
m.ix.zrho_trl  = i+1;      i = i+1;
m.ix.lsd_ids   = i+(1:2);  i = i+2;   % log(sd_id_sigma)
m.ix.zrho_ids  = i+1;      i = i+1;   % atanh(rho_id_sigma)
m.ix.rescor_mu  = i+1;        i = i+1;
m.ix.lsd_rescor = i+1;        i = i+1;
m.ix.z_rescor   = i+(1:NSUB); i = i+NSUB;
m.ix.z_id     = i+(1:2*NSUB); i = i+2*NSUB;
m.ix.z_trl    = i+(1:2*NTRL); i = i+2*NTRL;
m.ix.z_ids    = i+(1:2*NSUB); i = i+2*NSUB;
m.D = i;
m.mon_idx = [m.ix.b, m.ix.b_sigma, m.ix.lsd_id, m.ix.zrho_id, ...
    m.ix.lsd_trl, m.ix.zrho_trl, m.ix.lsd_ids, m.ix.zrho_ids, ...
    m.ix.rescor_mu, m.ix.lsd_rescor];
m.mon_names = {'b1','b2','b_sigma1','b_sigma2','log_sd_id1','log_sd_id2', ...
    'z_rho_id','log_sd_trl1','log_sd_trl2','z_rho_trl', ...
    'log_sd_id_sigma1','log_sd_id_sigma2','z_rho_id_sigma', ...
    'rescor_mu','log_sd_rescor'};
end

function start = local_diff_sserr_start(m)
%Per-event means at the event-conditional data means; log-SDs at rough guesses.
start = zeros(m.D,1);
allm = [m.meas1; m.meas2; m.meas_unp];
if isempty(allm); allm = 0; end
start(m.ix.b(1)) = mean([m.meas1; m.meas_unp(m.e1>0)]);
start(m.ix.b(2)) = mean([m.meas2; m.meas_unp(m.e2>0)]);
if any(isnan(start(m.ix.b))); start(m.ix.b) = mean(allm); end
s = std(allm); if ~isfinite(s) || s <= 0; s = 1; end
start(m.ix.b_sigma) = log(max(s/2,1e-3));
start(m.ix.lsd_id)  = log(max(s/3,1e-3));
start(m.ix.lsd_trl) = log(max(s/3,1e-3));
start(m.ix.lsd_ids) = log(0.3);
start(m.ix.lsd_rescor) = log(0.3);
end

function lp = local_diff_sserr_logp(theta,m)
%Unconstrained log-posterior (up to a constant) matching psyrat_build_stan_diff_sserr.
ix = m.ix;
b       = theta(ix.b);
b_sigma = theta(ix.b_sigma);
% three 2-D non-centered Cholesky blocks (location, trial, scale)
sd_id  = exp(theta(ix.lsd_id));   rho_id  = tanh(theta(ix.zrho_id));
sd_trl = exp(theta(ix.lsd_trl));  rho_trl = tanh(theta(ix.zrho_trl));
sd_ids = exp(theta(ix.lsd_ids));  rho_ids = tanh(theta(ix.zrho_ids));
z_id  = reshape(theta(ix.z_id),2,m.NSUB);
z_trl = reshape(theta(ix.z_trl),2,m.NTRL);
z_ids = reshape(theta(ix.z_ids),2,m.NSUB);
r_id_1  = sd_id(1).*z_id(1,:).';
r_id_2  = sd_id(2).*(rho_id.*z_id(1,:).' + sqrt(1-rho_id.^2).*z_id(2,:).');
r_trl_1 = sd_trl(1).*z_trl(1,:).';
r_trl_2 = sd_trl(2).*(rho_trl.*z_trl(1,:).' + sqrt(1-rho_trl.^2).*z_trl(2,:).');
r_sig_1 = sd_ids(1).*z_ids(1,:).';
r_sig_2 = sd_ids(2).*(rho_ids.*z_ids(1,:).' + sqrt(1-rho_ids.^2).*z_ids(2,:).');
% per-subject residual correlation hierarchy (0 when estimate_rescor is off)
rescor_mu = theta(ix.rescor_mu);
sd_rescor = exp(theta(ix.lsd_rescor));
z_rescor  = theta(ix.z_rescor);
if m.estimate_rescor
    rho_s = tanh(rescor_mu + sd_rescor.*z_rescor);   % NSUB x 1
else
    rho_s = 0.*z_rescor;
end

ll = 0;
% PAIRED observations: bivariate multi_normal_cholesky (per-row sigma + rho[sid])
if m.NPAIR > 0
    id = m.id_pair; trl = m.trl_pair;
    mu1 = b(1) + r_id_1(id) + r_trl_1(trl);
    mu2 = b(2) + r_id_2(id) + r_trl_2(trl);
    logsig1 = b_sigma(1) + r_sig_1(id);
    logsig2 = b_sigma(2) + r_sig_2(id);
    sig1 = exp(logsig1); sig2 = exp(logsig2);
    rho_n = rho_s(id);
    d1 = m.meas1 - mu1; d2 = m.meas2 - mu2;
    w1 = d1./sig1;
    w2 = (d2 - (sig2.*rho_n./sig1).*d1) ./ (sig2.*sqrt(1-rho_n.^2));
    ll = ll - sum(logsig1) - sum(logsig2) - 0.5*sum(log(1-rho_n.^2)) ...
       - 0.5*sum(w1.^2 + w2.^2);
end
% UNPAIRED single-event observations: univariate normal (event-selected effects)
if m.NUNP > 0
    e1 = m.e1; e2 = m.e2; id = m.id_unp; trl = m.trl_unp;
    mu  = e1.*(b(1) + r_id_1(id) + r_trl_1(trl)) ...
        + e2.*(b(2) + r_id_2(id) + r_trl_2(trl));
    logsig = e1.*(b_sigma(1) + r_sig_1(id)) + e2.*(b_sigma(2) + r_sig_2(id));
    sig = exp(logsig);
    ll = ll - sum(logsig) - sum((m.meas_unp - mu).^2 ./ (2*sig.^2));
end

P = m.P;
lp = ll ...
    + sum(local_normal(b,0,P.b)) ...
    + sum(local_student_t(b_sigma,3,0,P.b_sigma)) ...
    + local_student_t(sd_id(1),3,0,P.sd_id)   + theta(ix.lsd_id(1)) ...
    + local_student_t(sd_id(2),3,0,P.sd_id)   + theta(ix.lsd_id(2)) ...
    + local_student_t(sd_trl(1),3,0,P.sd_trl) + theta(ix.lsd_trl(1)) ...
    + local_student_t(sd_trl(2),3,0,P.sd_trl) + theta(ix.lsd_trl(2)) ...
    + local_student_t(sd_ids(1),3,0,P.sd_id_sigma) + theta(ix.lsd_ids(1)) ...
    + local_student_t(sd_ids(2),3,0,P.sd_id_sigma) + theta(ix.lsd_ids(2)) ...
    + 2*log(1 - rho_id.^2) ...                                % L_id  LKJ(2) + tanh Jac
    + 2*log(1 - rho_trl.^2) ...                               % L_trl LKJ(2) + tanh Jac
    + 2*log(1 - rho_ids.^2) ...                               % L_id_sigma LKJ(2) + tanh Jac
    + local_normal(rescor_mu,0,P.rescor_mu) ...               % rescor_mu ~ normal(0,.)
    + local_student_t(sd_rescor,3,0,P.sd_rescor) + theta(ix.lsd_rescor) ... % half-t + Jac
    - 0.5*sum(z_rescor.^2) ...                                % std_normal z_rescor
    - 0.5*sum(z_id(:).^2) - 0.5*sum(z_trl(:).^2) - 0.5*sum(z_ids(:).^2);   % std_normal raw
end

function out = local_diff_sserr_draws(draws,m)
%Transform draws to the CmdStan-shaped arrays case 8 extracts: b, b_sigma (ND x2),
%sd_id, sd_trl (ND x2), id_varcov, trl_varcov (ND x2x2), rescor (ND x1), err_varcov
%(ND x2x2, averaged over subjects), er_var_ss (ND x NSUB x2), wp_cov_ss (ND x NSUB).
ix = m.ix;
ND = size(draws,1); NSUB = m.NSUB;
out = struct();
out.b       = draws(:,ix.b);
out.b_sigma = draws(:,ix.b_sigma);
sd_id  = exp(draws(:,ix.lsd_id));   sd_trl = exp(draws(:,ix.lsd_trl));
out.sd_id  = sd_id;
out.sd_trl = sd_trl;
rho_id  = tanh(draws(:,ix.zrho_id));
rho_trl = tanh(draws(:,ix.zrho_trl));
out.id_varcov  = local_varcov2(sd_id,rho_id);
out.trl_varcov = local_varcov2(sd_trl,rho_trl);
% per-subject scale REs r_sigma -> er_var_ss, and rho_s -> wp_cov_ss / rescor /
% the subject-averaged err_varcov.
sd_ids  = exp(draws(:,ix.lsd_ids));
rho_ids = tanh(draws(:,ix.zrho_ids));
sd_rescor = exp(draws(:,ix.lsd_rescor));
er = zeros(ND,NSUB,2);
wp = zeros(ND,NSUB);
err_varcov = zeros(ND,2,2);
rescor = zeros(ND,1);
estrho = m.estimate_rescor;
for d = 1:ND
    Z = draws(d,ix.z_ids); Z1 = Z(1:2:end).'; Z2 = Z(2:2:end).';   % NSUB raw (col-major)
    rsig1 = sd_ids(d,1).*Z1;
    rsig2 = sd_ids(d,2).*(rho_ids(d).*Z1 + sqrt(1-rho_ids(d).^2).*Z2);
    e1 = draws(d,ix.b_sigma(1)) + rsig1;     % NSUB x1 (er_var_ss event1)
    e2 = draws(d,ix.b_sigma(2)) + rsig2;
    er(d,:,1) = e1.'; er(d,:,2) = e2.';
    if estrho
        rho_s = tanh(draws(d,ix.rescor_mu) + sd_rescor(d).*draws(d,ix.z_rescor).');
    else
        rho_s = zeros(NSUB,1);
    end
    wp(d,:) = (rho_s .* exp(e1) .* exp(e2)).';
    rescor(d) = mean(rho_s);
    err_varcov(d,1,1) = mean(exp(2*e1));
    err_varcov(d,2,2) = mean(exp(2*e2));
    err_varcov(d,1,2) = mean(wp(d,:));
    err_varcov(d,2,1) = err_varcov(d,1,2);
end
out.er_var_ss = er;
out.wp_cov_ss = wp;
out.err_varcov = err_varcov;
out.rescor = rescor;
end

function V = local_varcov2(sd,rho)
%ND x 2 x 2 covariance per draw from per-draw SDs (ND x2) and correlation (ND x1).
ND = size(sd,1);
V = zeros(ND,2,2);
V(:,1,1) = sd(:,1).^2;
V(:,2,2) = sd(:,2).^2;
V(:,1,2) = rho .* sd(:,1) .* sd(:,2);
V(:,2,1) = V(:,1,2);
end

% ===========================================================================
% Gradient plumbing
% ===========================================================================
function [lp,g] = local_ad_logpgrad(fad,theta)
%Wrap an AD log-posterior fad(theta_dl)->scalar dlarray into the [logp,grad]
%signature hmcSampler / psyrat_hmc_sample expect, returning plain doubles.
theta = theta(:);
[lp_dl,g_dl] = dlfeval(@local_ad_helper,fad,dlarray(theta));
lp = double(extractdata(lp_dl));
g  = double(extractdata(g_dl));
g  = g(:);
end

function [lp,g] = local_ad_helper(fad,theta_dl)
%Traced inside dlfeval: scalar log-posterior and its reverse-mode gradient.
lp = fad(theta_dl);
g  = dlgradient(lp,theta_dl);
end

function [lp,g] = local_numer_logpgrad(fdbl,theta)
%Central finite-difference gradient fallback (no Deep Learning Toolbox). O(D)
%log-posterior evaluations per gradient, so this is slow for large designs; the
%AD path is strongly preferred and used by default when available.
theta = theta(:);
lp = fdbl(theta);
D = numel(theta);
g = zeros(D,1);
h = 1e-5;
%Perturb theta in place one coordinate at a time and restore it, rather than
%allocating a fresh D-length basis vector per coordinate (which is O(D^2)
%allocation per gradient). The perturbed values are identical to theta+/-e, so
%the finite-difference result is bit-for-bit unchanged.
for i = 1:D
    ti = theta(i);
    theta(i) = ti + h; fp = fdbl(theta);
    theta(i) = ti - h; fm = fdbl(theta);
    theta(i) = ti;
    g(i) = (fp - fm) / (2*h);
end
end

function local_hmc_gradcheck(logpfun,start,m,use_ad)
%One-time, cheap sanity check that the supplied gradient agrees with central
%finite differences on a random subset of coordinates. Warns (does not error) on
%a mismatch so a real run is never silently driven by a wrong gradient. Skipped
%for the numeric path (its gradient IS finite differences).
if ~use_ad
    return;
end
[lp0,g] = logpfun(start);
if ~isfinite(lp0)
    warning('PsyRAT:hmcStartNonfinite',['The native HMC start point has a ' ...
        'non-finite log-posterior (%.3g); the sampler may fail to initialize.'],lp0);
end
%Pick the check coordinates from a fixed seed for reproducibility, but save and
%restore the global RNG stream so this one-time diagnostic does not perturb the
%state the sampler draws from (the sampler reseeds per chain regardless, so this
%keeps the gradcheck a true no-op on the estimation output).
rng_state = rng;
rng(12345,'twister');
k = min(20,m.D);
idx = randperm(m.D,k);
rng(rng_state);
h = 1e-5; gfd = zeros(k,1);
for a = 1:k
    e = zeros(m.D,1); e(idx(a)) = h;
    gfd(a) = (local_first_output(logpfun,start+e) - ...
        local_first_output(logpfun,start-e)) / (2*h);
end
relerr = norm(g(idx) - gfd) / max(1,norm(gfd));
if relerr > 1e-3
    warning('PsyRAT:hmcGradMismatch',['The native HMC automatic-differentiation ' ...
        'gradient disagreed with finite differences (rel. err %.2e on %d ' ...
        'coordinates). Estimation will proceed but should be checked.'],relerr,k);
end
end

function lp = local_first_output(fun,th)
[lp,~] = fun(th);
end

function v = ternary(cond,a,b)
if cond; v = a; else; v = b; end
end
