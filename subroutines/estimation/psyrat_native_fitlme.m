function fit = psyrat_native_fitlme(design,data,seed,nchains,opts)
%Native frequentist estimation engine: REML variance components via fitlme,
%with parametric-bootstrap draws, used as a CmdStan fallback for the simple
%plain-Gaussian, univariate designs (one-facet; later test-retest and splits).
%
%fit = psyrat_native_fitlme(design,data,seed,nchains,opts)
%
%Why a bootstrap, not a point estimate:
% The downstream calculation layer (psyrat_rel_sing, ...) consumes full draw
% VECTORS of each variance component (in SD units) and computes reliability per
% draw, then takes quantiles for the interval. A single REML point estimate per
% component cannot reproduce the interval, because reliability is a nonlinear
% function rel = bp/(bp + wp/n) and the interval must propagate the JOINT
% distribution of the components. A parametric bootstrap (simulate from the
% fitted model, refit, collect the refit components) yields aligned draw vectors
% that do propagate that joint uncertainty. These are frequentist sampling
% draws, not a Bayesian posterior, which is why the cascade/forced paths warn
% that the resulting intervals are approximate.
%
%Inputs:
% design - struct with .key selecting the model ('sing' = one-facet, analysis
%  case 1; more keys added as fitters land).
% data - the Stan data struct the matching case built (e.g. for 'sing':
%  .id (1..NSUB), .trl (1..NTRL), .meas). The fitter reads the same canonical
%  inputs CmdStan receives.
% seed - run seed. Native reproducibility uses the MATLAB RNG seeded here, so
%  the same seed + data give the same bootstrap draws.
% nchains - chains setting (carried only to size the convergence sentinel so
%  psyrat_checkconv's n_eff threshold, 10*nchains, is met on success).
% opts - struct of native knobs; uses opts.bootstrap_reps (default 1000).
%
%Output:
% fit - a PsyRATNativeFit with .draws holding SD draw vectors named exactly as
%  the matching case extracts them (e.g. mu, sig_u, sig_trl, sig_e), and a
%  convergence cell psyrat_storeconv/psyrat_checkconv can read.
%
%See also: psyrat_run_native, PsyRATNativeFit, psyrat_rel_sing

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

key = '';
if isstruct(design) && isfield(design,'key')
    key = char(design.key);
end

%number of parametric-bootstrap resamples
B = 1000;
if isstruct(opts) && isfield(opts,'bootstrap_reps') && ~isempty(opts.bootstrap_reps)
    B = round(opts.bootstrap_reps);
end
if B < 2
    error('PsyRAT:bootstrapReps',['bootstrap_reps must be >= 2 (got %g). ' ...
        'How to fix: pass a larger ''bootstrap_reps''.'],B);
end

[oldmsg,oldid] = lastwarn;
warncleanup = onCleanup(@()lastwarn(oldmsg,oldid)); %#ok<NASGU>

%reproducibility: the native engine uses the MATLAB RNG seeded with the run seed
rng(double(seed),'twister');

switch key
    case 'sing'
        fit = local_fit_sing(data,B,nchains);
    case 'sing_facet'
        fit = local_fit_sing_facet(design,data,B,nchains);
    case 'trt'
        fit = local_fit_trt(data,B,nchains);
    case 'splits'
        fit = local_fit_splits(data,B,nchains);
    case 'splits_trt'
        fit = local_fit_splits_trt(data,B,nchains);
    otherwise
        error('PsyRAT:fitlmeUnsupportedDesign',['The native fitlme engine ' ...
            'does not implement design ''%s''.'],key);
end
fit.meta.seed = double(seed);
fit.meta.design = key;
end

% ------------------------------------------------------------------------------
function fit = local_fit_trt(data,B,nchains)
%Test-retest crossed model (case 5): person x occasion x trial with all three
%two-way interactions. Components map to pop_int (intercept), sig_id, sig_occ,
%sig_trl (main-effect SDs), sig_trlxid, sig_occxid, sig_trlxocc (two-way
%interaction SDs), and sig_err (residual SD) - the same names the case-5
%extraction pulls. The CmdStan model uses the pre-computed interaction-level
%indices (trlxid/occxid/trlxocc) as flat random effects; this fitter does the
%same so the variance decomposition matches.
%
%DEFENSE-IN-DEPTH: the splits test-retest (case 24) has a near-identical call
%block but adds a per-observation 'weight' field (residual scaled by
%1/sqrt(weight)); routing it here would silently drop the weighting. Refuse a
%data struct that carries 'weight' so a mis-wiring fails loudly instead.
%
%PERFORMANCE: the interaction factors (especially trlxid) are high-cardinality,
%so each crossed-RE REML fit is comparatively slow and the parametric bootstrap
%multiplies that cost. For large test-retest datasets prefer a smaller
%bootstrap_reps, or CmdStan/HMC if available.
if isfield(data,'weight')
    error('PsyRAT:fitlmeTrtWeighted',['The native ''trt'' fitlme fitter was ' ...
        'given a weighted (splits) data struct. The splits test-retest design ' ...
        'needs a weighted fitter; this is a wiring error.']);
end
reqfields = {'id','occ','trl','trlxid','occxid','trlxocc','meas'};
missing = reqfields(~isfield(data,reqfields));
if ~isempty(missing)
    error('PsyRAT:fitlmeTrtFields',['The native ''trt'' fitlme fitter is ' ...
        'missing required data fields: %s.'],strjoin(missing,', '));
end

id   = double(data.id(:));
occ  = double(data.occ(:));
trl  = double(data.trl(:));
txi  = double(data.trlxid(:));
oxi  = double(data.occxid(:));
txo  = double(data.trlxocc(:));
meas = double(data.meas(:));
NOBS = numel(meas);

%reindex every grouping factor to 1..K for fast simulation
[~,~,idI]  = unique(id);
[~,~,occI] = unique(occ);
[~,~,trlI] = unique(trl);
[~,~,txiI] = unique(txi);
[~,~,oxiI] = unique(oxi);
[~,~,txoI] = unique(txo);

%Few-occasions caveat: the occasion variance (sig_occ) is estimated from only
%NOCC levels. With the handful of occasions typical of test-retest designs it is
%weakly identified (e.g. 1 df at NOCC=2). REML collapses such a component toward
%its boundary, whereas the Bayesian CmdStan model regularizes it via its prior -
%so the occasion variance, and the STABILITY/test-retest reliability that uses
%it, can differ MATERIALLY between this fallback and a CmdStan run. Warn so the
%result is read as a rough approximation, not a CmdStan-equivalent estimate.
NOCC = max(occI);
if NOCC < 5
    warning('PsyRAT:fitlmeTrtFewOccasions',['The native fitlme test-retest ' ...
        'engine estimated the occasion variance from only %d occasion level(s). ' ...
        'With few occasions this component is weakly identified and the ' ...
        'stability/test-retest reliability can differ MATERIALLY from the ' ...
        'Bayesian CmdStan estimate. Treat the result as a rough approximation; ' ...
        'final-inference readiness is unassessed. A CmdStan/HMC run also ' ...
        'requires diagnostic and uncertainty validation.'],NOCC);
end

formula = ['meas ~ 1 + (1|id) + (1|occ) + (1|trl) + (1|trlxid) + ' ...
    '(1|occxid) + (1|trlxocc)'];
Tbl = table(meas,categorical(idI),categorical(occI),categorical(trlI), ...
    categorical(txiI),categorical(oxiI),categorical(txoI), ...
    'VariableNames',{'meas','id','occ','trl','trlxid','occxid','trlxocc'});

%REML is requested EXPLICITLY at every fitlme call in this file (original fits
%and bootstrap refits alike): fitlme's default FitMethod is ML, and until
%2026-09-11 (post-beta review PBA-04, work queue Q1) every call here omitted the
%option, so the engine fitted by ML under its REML label.
lastwarn('');   % clear so we capture only this fit's warnings (REML convergence)
try
    lme = fitlme(Tbl,formula,'FitMethod','REML');
catch err
    error('PsyRAT:fitlmeFailed',['fitlme failed to fit the test-retest ' ...
        'model: %s'],err.message);
end
[fitmsg,fitwarn] = lastwarn;
hat = local_extract_trt(lme);
original = psyrat_fitlme_validity(hat,fitwarn,fitmsg);
converged = original.valid;

%grouping index vectors and the matching point-estimate SDs, in formula order
idx = {idI,occI,trlI,txiI,oxiI,txoI};
sds = [hat.sig_id,hat.sig_occ,hat.sig_trl,hat.sig_trlxid,hat.sig_occxid, ...
    hat.sig_trlxocc];

%parametric bootstrap: simulate y from the fitted crossed model, refit, collect
P = zeros(B,1); sid = zeros(B,1); socc = zeros(B,1); strl = zeros(B,1);
stxi = zeros(B,1); soxi = zeros(B,1); stxo = zeros(B,1); serr = zeros(B,1);
ok = true(B,1); failids = {};
records = repmat(psyrat_fitlme_validity(struct(),'',''),B,1);
for b = 1:B
    y = hat.pop_int + hat.sig_err*randn(NOBS,1);
    for g = 1:6
        eff = sds(g)*randn(max(idx{g}),1);
        y = y + eff(idx{g});
    end
    Tb = Tbl;
    Tb.meas = y;
    try
        lastwarn('');
        lb = fitlme(Tb,formula,'FitMethod','REML');
        [wmsg,wid] = lastwarn;
        hb = local_extract_trt(lb);
        records(b) = psyrat_fitlme_validity(hb,wid,wmsg);
        if ~records(b).valid
            ok(b) = false; failids{end+1} = records(b).error_id; %#ok<AGROW>
            continue;
        end
        P(b)=hb.pop_int; sid(b)=hb.sig_id; socc(b)=hb.sig_occ; strl(b)=hb.sig_trl;
        stxi(b)=hb.sig_trlxid; soxi(b)=hb.sig_occxid; stxo(b)=hb.sig_trlxocc; serr(b)=hb.sig_err;
    catch err
        %a failed refit is DROPPED, not replaced by the point estimate (Q3)
        ok(b) = false; failids{end+1} = err.identifier; %#ok<AGROW>
        [wmsg,wid] = lastwarn;
        records(b) = psyrat_fitlme_validity(struct(),wid,wmsg);
        records(b).valid = false;
        records(b).error_id = err.identifier; records(b).error_message = err.message;
    end
end
local_fitlme_bootcheck(ok,B,'test-retest');
[P,sid,socc,strl,stxi,soxi,stxo,serr] = local_keepok(ok,P,sid,socc,strl,stxi,soxi,stxo,serr);

draws = struct('pop_int',P,'sig_id',sid,'sig_occ',socc,'sig_trl',strl, ...
    'sig_trlxid',stxi,'sig_occxid',soxi,'sig_trlxocc',stxo,'sig_err',serr);
conv = local_fitlme_conv(converged,B,nchains);
meta = struct('method','fitlme_reml_parametric_bootstrap','formula',formula, ...
    'bootstrap_reps',B,'point',hat,'bootstrap_ok',nnz(ok), ...
    'bootstrap_failed',nnz(~ok),'bootstrap_fail_ids',{unique(failids)}, ...
    'original',original,'bootstrap_success_mask',ok,'bootstrap_draw_indices',find(ok), ...
    'bootstrap_records',records);
fit = PsyRATNativeFit(draws,conv,'fitlme',meta);
end

% ------------------------------------------------------------------------------
function hat = local_extract_trt(lme)
%SD-scale point estimates from a fitted test-retest model. RE covariances come
%back in formula order: id, occ, trl, trlxid, occxid, trlxocc.
[psi,mse] = covarianceParameters(lme);
if numel(psi) ~= 6
    error('PsyRAT:fitlmeComponents',['Expected six random-effect grouping ' ...
        'terms in the test-retest model, found %d.'],numel(psi));
end
fe = fixedEffects(lme);
hat = struct( ...
    'pop_int',     fe(1), ...
    'sig_id',      sqrt(psi{1}(1,1)), ...
    'sig_occ',     sqrt(psi{2}(1,1)), ...
    'sig_trl',     sqrt(psi{3}(1,1)), ...
    'sig_trlxid',  sqrt(psi{4}(1,1)), ...
    'sig_occxid',  sqrt(psi{5}(1,1)), ...
    'sig_trlxocc', sqrt(psi{6}(1,1)), ...
    'sig_err',     sqrt(mse));
end

% ------------------------------------------------------------------------------
function fit = local_fit_splits_trt(data,B,nchains)
%Splits test-retest (case 24, p x o x (i:s) observed design): the test-retest
%crossed model PLUS a person x occasion x split interaction (sig_posxid), with
%each row a split mean of n_i items (residual variance ~ sig_err^2/n_i). Emits
%pop_int + sig_id/sig_occ/sig_trl/sig_trlxid/sig_occxid/sig_trlxocc/sig_posxid/
%sig_err, matching the case-24 extraction. posxid = unique (id,occ,trl).
%
%CAVEATS: this is the heaviest native fitter (seven crossed REs, one with one row
%per level) and inherits BOTH the few-occasions identifiability of test-retest
%(sig_occ) and the weighted near-confounding of splits (sig_posxid vs sig_err).
%Treat as a rough approximation; prefer CmdStan/HMC. Requires a 'weight' field.
if ~isfield(data,'weight')
    error('PsyRAT:fitlmeSplitsTrtWeight',['The native ''splits_trt'' fitter ' ...
        'requires a ''weight'' (items-per-split n_i) field.']);
end
id  = double(data.id(:));  occ = double(data.occ(:));  trl = double(data.trl(:));
txi = double(data.trlxid(:)); oxi = double(data.occxid(:)); txo = double(data.trlxocc(:));
meas = double(data.meas(:)); w = double(data.weight(:));
NOBS = numel(meas);

[~,~,idI]  = unique(id);   [~,~,occI] = unique(occ);  [~,~,trlI] = unique(trl);
[~,~,txiI] = unique(txi);  [~,~,oxiI] = unique(oxi);  [~,~,txoI] = unique(txo);
[~,~,posI] = unique([idI occI trlI],'rows');   % person x occasion x split

NOCC = max(occI);
if NOCC < 5
    warning('PsyRAT:fitlmeTrtFewOccasions',['The native splits test-retest ' ...
        'engine estimated the occasion variance from only %d occasion level(s); ' ...
        'with few occasions the test-retest reliability can differ MATERIALLY ' ...
        'from the Bayesian CmdStan estimate.'],NOCC);
end
warning('PsyRAT:fitlmeSplitsIdentifiability',['The native splits test-retest ' ...
    'engine separates the person x occasion x split variance from the per-item ' ...
    'residual only through the items-per-split variation (a near-confounded ' ...
    'decomposition), so the observed-design reliability is approximate.']);

formula = ['meas ~ 1 + (1|id) + (1|occ) + (1|trl) + (1|trlxid) + (1|occxid) + ' ...
    '(1|trlxocc) + (1|posxid)'];
Tbl = table(meas,categorical(idI),categorical(occI),categorical(trlI), ...
    categorical(txiI),categorical(oxiI),categorical(txoI),categorical(posI), ...
    'VariableNames',{'meas','id','occ','trl','trlxid','occxid','trlxocc','posxid'});

lastwarn('');   % clear so we capture only this fit's warnings (REML convergence)
try
    lme = fitlme(Tbl,formula,'Weights',w,'FitMethod','REML');
catch err
    error('PsyRAT:fitlmeFailed',['fitlme failed to fit the splits test-retest ' ...
        'model: %s'],err.message);
end
[fitmsg,fitwarn] = lastwarn;
hat = local_extract_splits_trt(lme);
original = psyrat_fitlme_validity(hat,fitwarn,fitmsg);
converged = original.valid;

idx = {idI,occI,trlI,txiI,oxiI,txoI,posI};
sds = [hat.sig_id,hat.sig_occ,hat.sig_trl,hat.sig_trlxid,hat.sig_occxid, ...
    hat.sig_trlxocc,hat.sig_posxid];

P=zeros(B,1); sid=zeros(B,1); socc=zeros(B,1); strl=zeros(B,1);
stxi=zeros(B,1); soxi=zeros(B,1); stxo=zeros(B,1); spos=zeros(B,1); serr=zeros(B,1);
ok = true(B,1); failids = {};
records = repmat(psyrat_fitlme_validity(struct(),'',''),B,1);
for b = 1:B
    y = hat.pop_int + (hat.sig_err ./ sqrt(w)) .* randn(NOBS,1);
    for g = 1:7
        eff = sds(g)*randn(max(idx{g}),1);
        y = y + eff(idx{g});
    end
    Tb = Tbl; Tb.meas = y;
    try
        lastwarn('');
        lb = fitlme(Tb,formula,'Weights',w,'FitMethod','REML');
        [wmsg,wid] = lastwarn;
        hb = local_extract_splits_trt(lb);
        records(b) = psyrat_fitlme_validity(hb,wid,wmsg);
        if ~records(b).valid
            ok(b) = false; failids{end+1} = records(b).error_id; %#ok<AGROW>
            continue;
        end
        P(b)=hb.pop_int; sid(b)=hb.sig_id; socc(b)=hb.sig_occ; strl(b)=hb.sig_trl;
        stxi(b)=hb.sig_trlxid; soxi(b)=hb.sig_occxid; stxo(b)=hb.sig_trlxocc;
        spos(b)=hb.sig_posxid; serr(b)=hb.sig_err;
    catch err
        ok(b) = false; failids{end+1} = err.identifier; %#ok<AGROW>
        [wmsg,wid] = lastwarn;
        records(b) = psyrat_fitlme_validity(struct(),wid,wmsg);
        records(b).valid = false;
        records(b).error_id = err.identifier; records(b).error_message = err.message;
    end
end
local_fitlme_bootcheck(ok,B,'splits test-retest');
[P,sid,socc,strl,stxi,soxi,stxo,spos,serr] = local_keepok(ok,P,sid,socc,strl,stxi,soxi,stxo,spos,serr);

draws = struct('pop_int',P,'sig_id',sid,'sig_occ',socc,'sig_trl',strl, ...
    'sig_trlxid',stxi,'sig_occxid',soxi,'sig_trlxocc',stxo,'sig_posxid',spos, ...
    'sig_err',serr);
conv = local_fitlme_conv(converged,B,nchains);
meta = struct('method','fitlme_reml_parametric_bootstrap','formula',formula, ...
    'bootstrap_reps',B,'point',hat,'bootstrap_ok',nnz(ok), ...
    'bootstrap_failed',nnz(~ok),'bootstrap_fail_ids',{unique(failids)}, ...
    'original',original,'bootstrap_success_mask',ok,'bootstrap_draw_indices',find(ok), ...
    'bootstrap_records',records);
fit = PsyRATNativeFit(draws,conv,'fitlme',meta);
end

% ------------------------------------------------------------------------------
function hat = local_extract_splits_trt(lme)
%SD-scale estimates for the splits test-retest model. RE order: id, occ, trl,
%trlxid, occxid, trlxocc, posxid.
[psi,mse] = covarianceParameters(lme);
if numel(psi) ~= 7
    error('PsyRAT:fitlmeComponents',['Expected seven random-effect grouping ' ...
        'terms in the splits test-retest model, found %d.'],numel(psi));
end
fe = fixedEffects(lme);
hat = struct('pop_int',fe(1), ...
    'sig_id',sqrt(psi{1}(1,1)),'sig_occ',sqrt(psi{2}(1,1)),'sig_trl',sqrt(psi{3}(1,1)), ...
    'sig_trlxid',sqrt(psi{4}(1,1)),'sig_occxid',sqrt(psi{5}(1,1)), ...
    'sig_trlxocc',sqrt(psi{6}(1,1)),'sig_posxid',sqrt(psi{7}(1,1)), ...
    'sig_err',sqrt(mse));
end

% ------------------------------------------------------------------------------
function fit = local_fit_splits(data,B,nchains)
%Nonparallel splits, single occasion (case 23, p x (i:s) observed design). Each
%row is a split whose score is the split MEAN of n_i items, so its residual
%variance scales as sig_err^2 / n_i. Components map to mu, sig_u (person),
%sig_trl (split main effect), sig_splitxid (person x split interaction), sig_err
%(per-item residual). The weighted residual is reproduced with fitlme's
%'Weights' = n_i (which gives residual variance sigma^2 / n_i).
%
%IDENTIFIABILITY CAVEAT: sig_splitxid (constant variance) and sig_err (variance
%~ 1/n_i) are separated ONLY by the variation in n_i, and each splitxid level has
%one row - a near-confounded structure. REML trades the two off (sig_splitxid
%tends to be underestimated and sig_err inflated), so the observed-design
%reliability can differ from the Bayesian CmdStan estimate. Warned below.
if ~isfield(data,'weight')
    error('PsyRAT:fitlmeSplitsWeight',['The native ''splits'' fitlme fitter ' ...
        'requires a ''weight'' (items-per-split n_i) field.']);
end
id   = double(data.id(:));
trl  = double(data.trl(:));
meas = double(data.meas(:));
w    = double(data.weight(:));
NOBS = numel(meas);

[~,~,idI]  = unique(id);
[~,~,trlI] = unique(trl);
[~,~,sxiI] = unique([idI trlI],'rows');   % person x split interaction levels
NSUB = max(idI); NTRL = max(trlI); NSXI = max(sxiI);

warning('PsyRAT:fitlmeSplitsIdentifiability',['The native fitlme splits engine ' ...
    'separates the person x split variance (sig_splitxid) from the per-item ' ...
    'residual (sig_err) only through the variation in items-per-split; this is ' ...
    'a weakly identified (near-confounded) decomposition, so the observed-design ' ...
    'reliability is approximate and can differ from a CmdStan run.']);

formula = 'meas ~ 1 + (1|id) + (1|trl) + (1|splitxid)';
Tbl = table(meas,categorical(idI),categorical(trlI),categorical(sxiI), ...
    'VariableNames',{'meas','id','trl','splitxid'});

lastwarn('');   % clear so we capture only this fit's warnings (REML convergence)
try
    lme = fitlme(Tbl,formula,'Weights',w,'FitMethod','REML');
catch err
    error('PsyRAT:fitlmeFailed',['fitlme failed to fit the splits model: %s'], ...
        err.message);
end
[fitmsg,fitwarn] = lastwarn;
hat = local_extract_splits(lme);
original = psyrat_fitlme_validity(hat,fitwarn,fitmsg);
converged = original.valid;

%bootstrap: simulate with the weighted (split-mean) residual var = sig_err^2/n_i
mu = zeros(B,1); su = zeros(B,1); st = zeros(B,1); ssx = zeros(B,1); se = zeros(B,1);
ok = true(B,1); failids = {};
records = repmat(psyrat_fitlme_validity(struct(),'',''),B,1);
for b = 1:B
    ub  = hat.sig_u        * randn(NSUB,1);
    tb  = hat.sig_trl      * randn(NTRL,1);
    sxb = hat.sig_splitxid * randn(NSXI,1);
    eb  = (hat.sig_err ./ sqrt(w)) .* randn(NOBS,1);
    y = hat.mu + ub(idI) + tb(trlI) + sxb(sxiI) + eb;
    Tb = Tbl;
    Tb.meas = y;
    try
        lastwarn('');
        lb = fitlme(Tb,formula,'Weights',w,'FitMethod','REML');
        [wmsg,wid] = lastwarn;
        hb = local_extract_splits(lb);
        records(b) = psyrat_fitlme_validity(hb,wid,wmsg);
        if ~records(b).valid
            ok(b) = false; failids{end+1} = records(b).error_id; %#ok<AGROW>
            continue;
        end
        mu(b)=hb.mu; su(b)=hb.sig_u; st(b)=hb.sig_trl; ssx(b)=hb.sig_splitxid; se(b)=hb.sig_err;
    catch err
        ok(b) = false; failids{end+1} = err.identifier; %#ok<AGROW>
        [wmsg,wid] = lastwarn;
        records(b) = psyrat_fitlme_validity(struct(),wid,wmsg);
        records(b).valid = false;
        records(b).error_id = err.identifier; records(b).error_message = err.message;
    end
end
local_fitlme_bootcheck(ok,B,'splits');
[mu,su,st,ssx,se] = local_keepok(ok,mu,su,st,ssx,se);

draws = struct('mu',mu,'sig_u',su,'sig_trl',st,'sig_splitxid',ssx,'sig_err',se);
conv = local_fitlme_conv(converged,B,nchains);
meta = struct('method','fitlme_reml_parametric_bootstrap','formula',formula, ...
    'bootstrap_reps',B,'point',hat,'bootstrap_ok',nnz(ok), ...
    'bootstrap_failed',nnz(~ok),'bootstrap_fail_ids',{unique(failids)}, ...
    'original',original,'bootstrap_success_mask',ok,'bootstrap_draw_indices',find(ok), ...
    'bootstrap_records',records);
fit = PsyRATNativeFit(draws,conv,'fitlme',meta);
end

% ------------------------------------------------------------------------------
function hat = local_extract_splits(lme)
%SD-scale point estimates from a fitted splits model. RE covariances in formula
%order: id, trl (split main effect), splitxid (person x split).
[psi,mse] = covarianceParameters(lme);
if numel(psi) ~= 3
    error('PsyRAT:fitlmeComponents',['Expected three random-effect grouping ' ...
        'terms in the splits model, found %d.'],numel(psi));
end
fe = fixedEffects(lme);
hat = struct('mu',fe(1),'sig_u',sqrt(psi{1}(1,1)),'sig_trl',sqrt(psi{2}(1,1)), ...
    'sig_splitxid',sqrt(psi{3}(1,1)),'sig_err',sqrt(mse));
end

% ------------------------------------------------------------------------------
function fit = local_fit_sing(data,B,nchains)
%One-facet internal consistency (case 1): meas ~ 1 + (1|id) + (1|trl). Draws are
%named mu/sig_u/sig_trl/sig_e, matching the case-1 extraction.
core = local_fit_sing_core(data.id,data.trl,data.meas,B);
local_fitlme_bootcheck(core.ok,B,'one-facet');
draws = local_maskdraws(core.draws,core.ok);
conv = local_fitlme_conv(core.converged,B,nchains);
meta = struct('method','fitlme_reml_parametric_bootstrap','formula',core.formula, ...
    'bootstrap_reps',B,'point',core.point,'bootstrap_ok',nnz(core.ok), ...
    'bootstrap_failed',nnz(~core.ok),'bootstrap_fail_ids',{unique(core.failids)}, ...
    'original',core.original,'bootstrap_success_mask',core.ok, ...
    'bootstrap_draw_indices',find(core.ok),'bootstrap_records',core.records);
fit = PsyRATNativeFit(draws,conv,'fitlme',meta);
end

% ------------------------------------------------------------------------------
function fit = local_fit_sing_facet(design,data,B,nchains)
%Group/event one-facet internal consistency (cases 2 and 3). The CmdStan path
%fits one multi-stratum model with per-stratum parameters named with a suffix
%(G for groups, E for events) and a 1-based index, e.g. sig_u_G1, mu_E2. Here
%each stratum is an independent one-facet fit; draws are named identically so
%the per-stratum case extraction (fit.extract('pars','sig_u_G1') etc.) is
%unchanged. The Stan data struct stores each stratum's columns as id_<S>i,
%trl_<S>i, meas_<S>i.
if ~isfield(design,'suffix') || ~isfield(design,'nstrata')
    error('PsyRAT:fitlmeFacetSpec',['The sing_facet design requires .suffix ' ...
        '(''G'' or ''E'') and .nstrata in the native spec.']);
end
S = design.suffix;
n = design.nstrata;

draws = struct();
points = cell(n,1);
allconverged = true;
%Each stratum is fitted and simulated independently. Apply only its own mask.
%Component vectors within that fit retain a common replicate index.
fits = struct([]);
for i = 1:n
    id   = data.(sprintf('id_%s%d',S,i));
    trl  = data.(sprintf('trl_%s%d',S,i));
    meas = data.(sprintf('meas_%s%d',S,i));
    core = local_fit_sing_core(id,trl,meas,B);
    local_fitlme_bootcheck(core.ok,B,sprintf('%s%d one-facet',S,i));
    kept = local_maskdraws(core.draws,core.ok);
    names = fieldnames(kept);
    for k = 1:numel(names)
        draws.(sprintf('%s_%s%d',names{k},S,i)) = kept.(names{k});
    end
    points{i} = core.point;
    allconverged = allconverged && core.converged;
    record = struct('label',sprintf('%s%d',S,i),'formula',core.formula, ...
        'point',core.point,'original',core.original,'bootstrap_reps',B, ...
        'bootstrap_ok',nnz(core.ok),'bootstrap_failed',nnz(~core.ok), ...
        'bootstrap_success_mask',core.ok,'bootstrap_draw_indices',find(core.ok), ...
        'bootstrap_records',core.records,'bootstrap_fail_ids',{unique(core.failids)});
    if isempty(fits); fits = record; else; fits(end+1) = record; end %#ok<AGROW>
end

conv = local_fitlme_conv(allconverged,B,nchains);
meta = struct('method','fitlme_reml_parametric_bootstrap', ...
    'suffix',S,'nstrata',n,'bootstrap_reps',B,'points',{points}, ...
    'bootstrap_ok',[fits.bootstrap_ok],'bootstrap_failed',[fits.bootstrap_failed], ...
    'fits',fits);
fit = PsyRATNativeFit(draws,conv,'fitlme',meta);
end

% ------------------------------------------------------------------------------
function core = local_fit_sing_core(id,trl,meas,B)
%Fit one one-facet model (meas ~ 1 + (1|id) + (1|trl)) by REML and draw B
%parametric-bootstrap replicates. Components map to mu (intercept), sig_u
%(person SD), sig_trl (trial main-effect SD = sigma_i), sig_e (residual SD).
%Returns a struct with .draws (Bx1 SD vectors), .point (REML estimates),
%.converged, and .formula. Shared by the one-facet and group/event fitters.
id   = double(id(:));
trl  = double(trl(:));
meas = double(meas(:));
NOBS = numel(meas);

%canonical 1..K integer indices for fast bootstrap simulation (the cases pass
%sequential ids, but reindex defensively)
[~,~,idIdx]  = unique(id);
[~,~,trlIdx] = unique(trl);
NSUB = max(idIdx);
NTRL = max(trlIdx);

formula = 'meas ~ 1 + (1|id) + (1|trl)';
Tbl = table(meas,categorical(idIdx),categorical(trlIdx), ...
    'VariableNames',{'meas','id','trl'});

%fit the observed data once for the point estimates that drive the bootstrap
lastwarn('');   % clear so we capture only this fit's warnings (REML convergence)
try
    lme = fitlme(Tbl,formula,'FitMethod','REML');
catch err
    error('PsyRAT:fitlmeFailed',['fitlme failed to fit the one-facet model: ' ...
        '%s'],err.message);
end
[fitmsg,fitwarn] = lastwarn;
hat = local_extract_sing(lme);
original = psyrat_fitlme_validity(hat,fitwarn,fitmsg);
converged = original.valid;

%parametric bootstrap: simulate from the fitted model, refit, collect draws
mu = zeros(B,1); su = zeros(B,1); st = zeros(B,1); se = zeros(B,1);
ok = true(B,1); failids = {};
records = repmat(psyrat_fitlme_validity(struct(),'',''),B,1);
for b = 1:B
    ub = hat.sig_u   * randn(NSUB,1);
    tb = hat.sig_trl * randn(NTRL,1);
    eb = hat.sig_e   * randn(NOBS,1);
    yb = hat.mu + ub(idIdx) + tb(trlIdx) + eb;
    Tb = Tbl;
    Tb.meas = yb;
    try
        lastwarn('');
        lb = fitlme(Tb,formula,'FitMethod','REML');
        [wmsg,wid] = lastwarn;
        hb = local_extract_sing(lb);
        records(b) = psyrat_fitlme_validity(hb,wid,wmsg);
        if ~records(b).valid
            ok(b) = false; failids{end+1} = records(b).error_id; %#ok<AGROW>
            continue;
        end
        mu(b) = hb.mu; su(b) = hb.sig_u; st(b) = hb.sig_trl; se(b) = hb.sig_e;
    catch err
        %a failed refit is DROPPED, not replaced by the point estimate (Q3,
        %2026-09-11): the replacement was an artificial point mass that
        %narrowed the interval (post-beta review PBA-05). The callers apply
        %the mask and enforce the adequacy floor (local_fitlme_bootcheck).
        ok(b) = false; failids{end+1} = err.identifier; %#ok<AGROW>
        [wmsg,wid] = lastwarn;
        records(b) = psyrat_fitlme_validity(struct(),wid,wmsg);
        records(b).valid = false;
        records(b).error_id = err.identifier; records(b).error_message = err.message;
    end
end
%The caller enforces the floor on each independent fit before dropping rows.
core = struct('draws',struct('mu',mu,'sig_u',su,'sig_trl',st,'sig_e',se), ...
    'point',hat,'converged',converged,'formula',formula,'ok',ok,'failids',{failids}, ...
    'original',original,'records',records);
end

% ------------------------------------------------------------------------------
function hat = local_extract_sing(lme)
%Pull SD-scale point estimates from a fitted one-facet LinearMixedModel.
%covarianceParameters returns the random-effect covariances in formula order, so
%for the fixed formula 'meas ~ 1 + (1|id) + (1|trl)' authored in local_fit_sing,
%psi{1} is the person (id) variance and psi{2} is the trial (trl) variance. The
%order is deterministic because this function controls the formula; guard the
%count so any future formula change fails loudly rather than silently mislabeling
%components.
[psi,mse] = covarianceParameters(lme);
if numel(psi) ~= 2
    error('PsyRAT:fitlmeComponents',['Expected exactly two random-effect ' ...
        'grouping terms (id, trl) in the one-facet model, found %d.'],numel(psi));
end
sig_u   = sqrt(psi{1}(1,1));   % person SD (sigma_p)
sig_trl = sqrt(psi{2}(1,1));   % trial main-effect SD (sigma_i)
sig_e   = sqrt(mse);           % residual SD (sigma_e)
fe = fixedEffects(lme);
hat = struct('mu',fe(1),'sig_u',sig_u,'sig_trl',sig_trl,'sig_e',sig_e);
end

% ------------------------------------------------------------------------------
function conv = local_fitlme_conv(converged,B,nchains)
%Build the convergence cell in psyrat_storeconv layout. fitlme is not an MCMC
%method, so there is no genuine R-hat/n_eff; instead encode REML convergence so
%psyrat_checkconv's thresholds (r_hat < 1.1 and n_eff >= 10*nchains) resolve to
%the correct converged/not decision. A converged fit reports r_hat = 1 and a
%large n_eff; a failed extraction reports r_hat = Inf / n_eff = 0.
hdr = {'name','n_eff','r_hat'};
if converged
    neff = max(B,10*nchains + 1);
    conv = [hdr; {'fitlme_reml',neff,1.0}];
else
    conv = [hdr; {'fitlme_reml',0,Inf}];
end
end

% ------------------------------------------------------------------------------
function local_fitlme_bootcheck(ok,B,label)
%Enforce the parametric-bootstrap adequacy floor and disclose every failed refit
%(Q3, 2026-09-11; post-beta review PBA-05). A failed refit is dropped, never
%replaced by the point estimate: the replacement put an artificial point mass
%at the point fit and narrowed the interval. Below the floor the run is refused
%with an error rather than a "not converged" verdict, because the automatic
%reruns of the warp driver double iterations that this engine ignores and would
%repeat the same failures. The floor is the owner's ruling: at least 90% of the
%requested replicates must succeed, and at least 200 when 200 or more were
%requested (smaller requests keep only the 90% rule, so a small test-size
%bootstrap stays usable).
nok = nnz(ok);
nfail = B - nok;
nmin = ceil(0.9*B);
if B >= 200
    nmin = max(nmin,200);
end
if nok < nmin
    error('PsyRAT:fitlmeBootstrapFloor',['%d of %d parametric-bootstrap refits ' ...
        'of the native fitlme %s model failed; at least %d successful refits ' ...
        'are required (90%% of the requested count, and 200 when 200 or more are ' ...
        'requested), so no interval was formed from the %d surviving replicates. ' ...
        'How to fix: use CmdStan or the native HMC engine, request more ' ...
        'replicates, or simplify the design.'],nfail,B,label,nmin,nok);
end
if nfail > 0
    warning('PsyRAT:fitlmeBootstrapFailures', ...
        ['%d of %d (%.0f%%) parametric-bootstrap refits of the native fitlme ' ...
        '%s model failed and were dropped; the interval uses the %d successful ' ...
        'replicates. Failed refits can signal a weakly identified model; treat ' ...
        'the interval as approximate and prefer CmdStan/HMC.'], ...
        nfail,B,100*nfail/B,label,nok);
end
end

% ------------------------------------------------------------------------------
function varargout = local_keepok(ok,varargin)
%Keep the rows of every draw vector whose replicate succeeded (same mask for
%all of them, so the vectors stay aligned for the calculation layer).
varargout = cell(1,numel(varargin));
for k = 1:numel(varargin)
    v = varargin{k};
    varargout{k} = v(ok);
end
end

% ------------------------------------------------------------------------------
function d = local_maskdraws(d,ok)
%Apply one success mask to every field of a draws struct.
f = fieldnames(d);
for k = 1:numel(f)
    v = d.(f{k});
    d.(f{k}) = v(ok);
end
end
