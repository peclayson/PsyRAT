function psyrat_data = psyrat_run(varargin)
%Run a complete PsyRAT reliability analysis from the command line (no GUI).
%
%psyrat_run is the scripted entry point that mirrors the GUI "Process New
%Data" workflow. It builds the same psyrat_prefs/psyrat_data scaffold the GUI
%builds, then calls the same shared estimation path
%(psyrat_loadfile -> psyrat_computevarcompwarp -> psyrat_computevarcomp). Given
%the same inputs and seed, a psyrat_run analysis produces the same variance
%components, convergence diagnostics, and saved .psyrat file as the GUI.
%
%Example (one-facet dependability from a long-format CSV):
%  psyrat_data = psyrat_run( ...
%      'file','/data/study/ern_scores.csv', ...
%      'idcol','id','meascol','meas','eventcol','event', ...
%      'savepath','/data/study','savename','ern_study.psyrat', ...
%      'chains',4,'warmup',5000,'sampling',5000,'seed',12345);
%
%Example (in-memory table, test-retest with an occasion facet):
%  tbl = readtable('scores.csv');
%  psyrat_data = psyrat_run('data',tbl,'meascol','meas', ...
%      'timecol','session','savepath',pwd,'savename','rt.psyrat');
%
%==================== Required inputs (one data source) ====================
% file - full path to a long-format data file (one row per single-trial
%  score) that readtable can parse (.csv, .txt, .xlsx, ...).
%      OR
% data - a long-format MATLAB table already in memory (use instead of 'file').
%
% savepath - existing directory where the .psyrat output is written. Must
%  contain no whitespace (CmdStan restriction).
% savename - output filename (e.g. 'study.psyrat'). A '.psyrat' extension is
%  added if omitted. Must contain no whitespace.
%
%==================== Column-mapping inputs (optional) =====================
% idcol    - participant id column name (default 'id').
% meascol  - numeric score column to analyze (default 'meas'). The selected
%  column is renamed to 'meas' internally.
% eventcol - event/condition column name (default '' = single event).
% groupcol - group column name (default '' = single group).
% timecol  - occasion column name for test-retest designs (default '' = one
%  occasion).
% weightcol- items-per-split (n_i) column name; required for splits analyses
%  (default '').
% dim1col, dim2col - continuous between-person covariate column(s) for
%  dynamic/conditional reliability (default '').
% whichevents/whichgroups/whichtimes - cell arrays restricting which labels in
%  the corresponding column are analyzed (default {} = use all).
%
%==================== Estimation settings (optional) ======================
% chains   - number of CmdStan chains (default 4). Must be an integer of 3 or
%  higher, the same floor the GUI enforces: R-hat is a ratio of between-chain
%  to within-chain variance, so it is undefined for one chain and unstable for
%  two, and a run with fewer than three chains would report a convergence
%  verdict that carries no information.
% warmup   - warmup iterations per chain (default 5000). Integer, 0 or higher.
% sampling - post-warmup sampling iterations per chain (default 5000). Integer,
%  1 or higher.
% seed     - CmdStan RNG seed for reproducibility (default 12345). Integer, 0
%  or higher.
% adapt_delta - NUTS target acceptance rate handed to CmdStan (delta=). Empty
%  [] (the default) keeps each model's own tuned value, so default runs are
%  unchanged. A supplied value must lie strictly between 0 and 1 (e.g. 0.95);
%  raising it usually removes divergent transitions at the cost of speed.
%  REPRODUCIBILITY: this alters step-size adaptation during warmup, so it
%  changes the draws even at a fixed seed. It is therefore recorded in the
%  run-config sidecar and restored by psyrat_run('config',...).
% max_treedepth - NUTS maximum tree depth handed to CmdStan (max_depth=). Empty
%  [] (the default) keeps CmdStan's default of 10. A supplied value must be a
%  positive integer (e.g. 12). Same reproducibility note as adapt_delta.
% verbose  - 1 = quiet, 2 = print CmdStan iterations (default 2, matching the
%  GUI default in psyrat_defaults).
% withinchain - 1 = off, 2 = within-chain parallelization (default 1).
% priors   - prior structure (default psyrat_default_priors via
%  psyrat_defaults). See documentation/priors_and_sensitivity.md.
% engine   - estimation engine. Accepts 'auto'/'cmdstan'/'fitlme'/'hmc' or the
%  numeric codes 0/1/2/3 (default: from psyrat_defaults = 0, auto). 'auto'
%  tries CmdStan, then the native HMC engine, then the native fitlme engine,
%  using only engines installed on this machine, and warns whenever it falls
%  back off CmdStan. A forced value is an explicit opt-in and errors if that
%  engine is unavailable or cannot handle the requested design. The native
%  fitlme engine supports only the simple one-facet/test-retest/splits designs
%  and reports approximate (non-Bayesian) intervals.
% bootstrap_reps - parametric-bootstrap resamples for the native fitlme engine
%  (default 1000). Ignored by the CmdStan and HMC engines.
% hmc_gradient - gradient method for the native HMC engine: 'auto' (default),
%  'ad' (automatic differentiation; needs the Deep Learning Toolbox), or
%  'numer' (numerical). Ignored by the CmdStan and fitlme engines.
% family   - observation-likelihood family: 'gaussian' (default, normal model)
%  or 'gamma' (scaled chi-square / Gamma log-link, for strictly positive,
%  mean-linked scores such as time-frequency power). The engine validates the
%  value and rejects an unsupported design.
%
%==================== Design / analysis flags (optional) ==================
% ssubjrel - subject-level reliability: 1 = off (default), 2 = on. ('sserrvar'
%  is accepted as a legacy alias.)
% diffest  - difference scores: 1 = off (default), 2 = 2-event difference,
%  3 = 4-event difference-of-differences.
% dodmap   - 1x4 cell of event labels in ERP_1..ERP_4 order; required when
%  diffest = 3.
% diffrescor / diffwpcov - the SAME setting under two names: whether the
%  within-person residual covariance between co-occurring difference-score
%  events is estimated. 1 = off, covariance fixed to 0 (default); 2 = on
%  ("correlated residuals"). Only relevant when diffest = 2; a diffest = 3
%  (difference-of-differences) run rejects a concurrent request outright
%  rather than silently ignoring it. The GUI exposes a
%  single control ("Estimate residual covariance for co-occurring events") and
%  writes both names from it, and psyrat_computevarcomp treats them as
%  synonyms, so pass EITHER one and the other follows. Passing both with
%  contradictory values is an error: the two values select structurally
%  different models, so there is no safe way to guess which was meant.
% dynrel   - dynamic/conditional reliability: 1 = off (default), 2 = on.
% splits   - reliability from splits: 1 = single-trial (default), 2 = parallel,
%  3 = nonparallel. Requires weightcol.
% ssrescor - residual-correlation parameterization for the subject-level
%  concurrent dynamic difference design: 1 = population (default), 2 =
%  per-subject.
% dispersion - copula residual-dispersion parameterization for the GROUP-LEVEL
%  concurrent gamma difference designs (family='gamma', diffest=2, diffrescor=2
%  -- equivalently diffwpcov=2, see above -- sserrvar=1): 1 = per-person log-nu
%  (experimental), 2 = global fixed dispersion
%  (owner reference; group-level only, not for subject-level reliability). Omit for
%  the design default, which is GLOBAL (2) for that design and per-person (1)
%  elsewhere.
% gammascale - scale parameterization for the Gamma dispersion submodel
%  (family='gamma'). Supported for the DYNAMIC designs -- one-facet
%  non-difference (analysis 11, or 26 with sserrvar=2), two-facet non-difference
%  (14, or 27 with sserrvar=2; both need a time/occasion column), and the
%  one-facet dynamic DIFFERENCE at the group (12) or subject (13, sserrvar=2)
%  level -- and for the STATIC subject-level designs (6, or 25 with a
%  time/occasion column; both need sserrvar=2) and the non-concurrent static
%  DIFFERENCE designs (7, or 8 with sserrvar=2; two-facet 10 with a
%  time/occasion column; four-event difference-of-differences 9) -- and for the
%  subject-level dynamic CONCURRENT difference designs, one-facet (19) and
%  trial + occasion two-facet (20), which use the modular cut-copula workflow
%  and are location-scale ONLY, and the person-specific dynamic NONCONCURRENT
%  difference-of-differences (28, or 29 with a time/occasion column), also
%  location-scale only. Not supported for analyses 1-5, which are excluded on
%  the merits rather than as unbuilt work,
%  nor for the STATIC concurrent difference designs (7/8/10 with diffrescor=2,
%  equivalently diffwpcov=2), whose copula model has no observed-scale residual
%  covariance defined for the direct residual-SD parameterization. NOTE this
%  clause used to read "nor for any CONCURRENT variant", which stopped being
%  true when 19/20 landed: they are concurrent AND location-scale only. The
%  authoritative list is the accept-list guard in psyrat_computevarcomp -- grep
%  `if gammascale == 2 && ~any(analysis ==` -- and the capability overview is
%  documentation/model_availability_matrix.md.
%  1 = log-nu dispersion, 2 = log residual SD (location-scale).
%  DEFAULT, since 2026-08-07 (rulings 33-36): omitting gammascale selects 2 on
%  every SUPPORTED design listed above, and 1 everywhere else -- analyses 1-5,
%  which have no location-scale variant, and the STATIC concurrent variants,
%  which reject 2. The default is therefore per-design, not global. This CHANGED: it
%  was 1 everywhere until 2026-08-07, so a script that relied on the old default
%  and did not name gammascale explicitly now fits a different estimand. Pass
%  gammascale=1 to restore the previous behavior -- EXCEPT on the
%  location-scale-only designs, which reject an explicit 1: analyses 13, 19,
%  20, 28 and 29 (the guards in psyrat_computevarcomp; 13 has no derived
%  log-nu estimand). Runs made from
%  2026-08-03 onward record the RESOLVED value in the saved result and the
%  run-config sidecar, so replaying one of those reproduces its own run. A sidecar
%  written BEFORE that date carries no gammascale at all and replays under the
%  current default -- pass gammascale=1 explicitly when replaying an older config.
%  Under 1 the conditional variance is 2*mu^2/nu, which makes
%  the model a constant coefficient-of-variation one: CV = sqrt(2/nu) (= 1/3 at
%  the default nu = 18), so every log-nu term is a statement about RELATIVE
%  dispersion. WHICH term carries it depends on the design, and the parameter
%  NAMES are reused across roles -- check the design before reading a posterior:
%  non-difference DYNAMIC (11/26, 14/27) declares b_nu / b_sigma as the dimension
%  SLOPES; the DIFFERENCE designs (7, 8, 9, 10, 12, 13) declare them as per-event
%  (per-cell, length 4, for 9) log-scale INTERCEPTS, with 12/13's dimension slopes
%  carried by the separate b_nu_dim / b_sigma_dim; and the static subject-level
%  designs (6, 25) declare neither, carrying dispersion in Intercept_nu (or
%  Intercept_sigma) plus the person effects ind_nu (or ind_sd). What blends is
%  correspondingly the dimension slope, the event contrast, or the person
%  contrast -- one parameterization issue, three predictors. Because
%  log sigma = log mu + 0.5*log(2) - 0.5*log(nu), a log-nu term blends amplitude
%  with absolute residual SD and so cannot be read as a statement about precision
%  alone. The identity is b_nu = 2*b - 2*b_sigma for a slope, and in CONTRAST
%  form for the intercept designs (the 0.5*log(2) cancels only under
%  differencing). Under 1 a POSITIVE b_nu slope means DECREASING relative error
%  (the implied log-CV slope is -b_nu/2). Under 2 the conditional variance is
%  exactly sigma^2 and the corresponding log-sigma term is a PURE absolute
%  residual-SD quantity. Option 1 is a valid
%  model, not an error -- it simply does not support an absolute-precision
%  reading; 2 is preferred whenever the output is described as absolute
%  precision, individual error variance, or Gaussian-comparable dependability,
%  which is why it is now the default where it is available. The two are
%  DIFFERENT ESTIMANDS and their results are not comparable.
%
%==================== Execution control (optional) ========================
% config   - path to a *_runconfig.json sidecar written by a previous run
%  (see psyrat_write_runconfig). Any option the caller does not explicitly pass
%  is filled in from the recorded configuration; explicit arguments always win,
%  so a replay can be redirected to new data or a new output file. A replay
%  warns (it does not stop) when this machine cannot reproduce the recorded run
%  bit for bit: when the detected CmdStan version differs from the recorded
%  one, and when the recorded run was produced by an engine the restored
%  preference will not necessarily select again (the sidecar records both the
%  requested 'engine' and the 'engine_used'; only the request is restored).
% interactive - false (default) runs fully headless: on non-convergence the
%  model auto-reruns with more iterations up to a cap and then errors, with no
%  GUI prompts. Set true to allow the interactive rerun dialog. It does NOT
%  enable trace-plot review: this function pins psyrat_prefs.proc.traceplots
%  to 1 (see below) and exposes no option to change it, and every trace-plot
%  branch in psyrat_computevarcompwarp requires 2, so trace plots are
%  unreachable from psyrat_run on either setting.
%
%==================== Output =============================================
% psyrat_data - the toolbox data structure with estimated variance components
%  and convergence diagnostics in psyrat_data.rel. The same structure is saved
%  to fullfile(savepath, savename). Open it later with the GUI "View Results"
%  button or pass psyrat_data.rel to the calculation/table functions.

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

%-------------------------------------------------------------------------
% 1. Parse inputs. inputParser gives named options, defaults, and clear
%    errors for a user-facing command-line entry point. Heavy data-level
%    validation is intentionally deferred to the shared validator
%    (psyrat_preflight_validate, invoked inside psyrat_computevarcompwarp) so
%    the CLI and GUI enforce exactly the same rules.
%-------------------------------------------------------------------------
p = inputParser;
p.FunctionName = 'psyrat_run';

% data source (exactly one of file/data must be supplied; checked below)
addParameter(p,'file','',@(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p,'data',[],@(x) istable(x));

% replay a previously written run-config sidecar (provides settings only; the
% data source and output location are still supplied/overridable as normal
% arguments). See psyrat_write_runconfig.
addParameter(p,'config','',@local_ischarlike);

% output location
addParameter(p,'savepath','',@(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p,'savename','',@(x) ischar(x) || (isstring(x) && isscalar(x)));

% column mapping
addParameter(p,'idcol','id',@local_ischarlike);
addParameter(p,'meascol','meas',@local_ischarlike);
addParameter(p,'eventcol','',@local_ischarlike);
addParameter(p,'groupcol','',@local_ischarlike);
addParameter(p,'timecol','',@local_ischarlike);
addParameter(p,'weightcol','',@local_ischarlike);
addParameter(p,'dim1col','',@local_ischarlike);
addParameter(p,'dim2col','',@local_ischarlike);
addParameter(p,'whichevents',{},@(x) iscell(x) || isempty(x));
addParameter(p,'whichgroups',{},@(x) iscell(x) || isempty(x));
addParameter(p,'whichtimes',{},@(x) iscell(x) || isempty(x));

% estimation settings
addParameter(p,'chains',4,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'warmup',5000,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'sampling',5000,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'seed',12345,@(x) isnumeric(x) && isscalar(x));
% optional CmdStan sampler overrides. [] means "use the default" (each model's
% tuned adapt_delta; CmdStan's max_depth = 10), which is what keeps default runs
% and the golden Stan snapshots byte-identical. Bounds are enforced after the
% run-config merge below, with the same limits the GUI applies.
addParameter(p,'adapt_delta',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p,'max_treedepth',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p,'verbose',2,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'withinchain',1,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'priors',[]);
% estimation engine ([] = use the psyrat_defaults value; accepts a numeric code
% 0/1/2/3 or a name 'auto'/'cmdstan'/'fitlme'/'hmc', normalized below)
addParameter(p,'engine',[],@(x) (isnumeric(x) && isscalar(x)) || local_ischarlike(x));
addParameter(p,'bootstrap_reps',[],@(x) isnumeric(x) && isscalar(x));
addParameter(p,'hmc_gradient','',@local_ischarlike);

% design / analysis flags
addParameter(p,'ssubjrel',[],@(x) isnumeric(x) && isscalar(x));
addParameter(p,'sserrvar',[],@(x) isnumeric(x) && isscalar(x)); % legacy alias
addParameter(p,'diffest',1,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'dodmap',{},@(x) iscell(x));
addParameter(p,'diffrescor',1,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'diffwpcov',1,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'dynrel',1,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'splits',1,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'ssrescor',1,@(x) isnumeric(x) && isscalar(x));
% copula residual-dispersion parameterization for the GROUP-LEVEL concurrent gamma
% difference designs (family='gamma', diffest=2, diffrescor=2, sserrvar=1): 1 =
% per-person log-nu (experimental); 2 = global fixed dispersion (owner reference,
% group-level only). Default [] = unset, so the engine applies its context-aware
% default (GLOBAL for that design; per-person elsewhere). The engine validates that 2
% is used only in that context.
addParameter(p,'dispersion',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x)));
% scale parameterization for the Gamma family: 1 = log-nu dispersion; 2 = log
% residual SD (location-scale). Default [] = unset, so the engine resolves it by
% design (location-scale where the design supports it, log-nu otherwise; see the
% accept lists in psyrat_computevarcomp - this comment twice went stale tracking
% the analysis numbers in prose, so it no longer names them).
addParameter(p,'gammascale',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x)));
% observation-likelihood family ('gaussian' default / 'gamma'); the engine
% validates the value and rejects an unsupported design, so accept any charlike
addParameter(p,'family','gaussian',@local_ischarlike);

% execution control
addParameter(p,'interactive',false,@(x) islogical(x) || (isnumeric(x) && isscalar(x)));

parse(p,varargin{:});
opt = p.Results;

% Which options the caller passed EXPLICITLY. Captured before the run-config
% merge below, because the merge rewrites opt but not this list: several checks
% further down have to distinguish "the caller asked for 1" from "1 is just the
% default", and a value's identity alone cannot tell them apart.
usingDefaults = p.UsingDefaults;

% If a run-config sidecar was supplied, fill in any option the caller did NOT
% explicitly set from the recorded configuration. Explicit arguments always
% win, so a replay can be redirected to new data or a new output file.
% cfgfilled lists the options the sidecar actually supplied (as opposed to the
% ones left at their defaults), which the residual-covariance reconciliation
% below needs for the same reason usingDefaults is needed.
cfgfilled = {};
if ~isempty(char(opt.config))
    [opt,cfgfilled] = local_apply_config(opt,char(opt.config),usingDefaults);
end

%-------------------------------------------------------------------------
% 1b. Validate the estimation settings. This runs AFTER the run-config merge
%     so a replayed sidecar is held to exactly the same rules as an explicit
%     argument. The bounds mirror the GUI's preference-window validation
%     (psyrat_startproc, psyrat_prefs_save), so a setting the GUI refuses
%     cannot be smuggled in through the command line instead.
%-------------------------------------------------------------------------
local_check_integer(opt.chains,'chains',3, ...
    'pass ''chains'' as an integer of 3 or higher (default 4)');
local_check_integer(opt.warmup,'warmup',0, ...
    'pass ''warmup'' as an integer of 0 or higher (default 5000)');
local_check_integer(opt.sampling,'sampling',1, ...
    'pass ''sampling'' as an integer of 1 or higher (default 5000)');
local_check_integer(opt.seed,'seed',0, ...
    'pass ''seed'' as an integer of 0 or higher (default 12345)');

% Optional sampler overrides: empty is the documented default and is left
% alone; a supplied value is checked against the GUI's limits.
%
% COMPLEX-LITERAL GUARD. isreal() sits inside the shape clause so it
% short-circuits everything after it. Without it, isfinite() is TRUE for a
% complex and the ordering comparisons read only the REAL part, so 0.5+0.5i
% cleared the adapt_delta check here AND the engine's own guard in
% psyrat_computevarcomp -- an unpatched engine run was measured completing a
% CmdStan fit instead of erroring (2026-08-12); what CmdStan then made of the
% token delta=0.5+0.5i was NOT established. max_treedepth did fail, but inside
% mod() with MATLAB:math:mustBeReal rather than the message below.
% local_check_integer (further down this file) already guards chains / warmup /
% sampling / seed this way, so these two were the odd pair out. Note the
% addParameter validators near the top of this file are weaker still
% (isnumeric && isscalar only) for those four; local_check_integer, not the
% parser, is what enforces realness there.
if ~isempty(opt.adapt_delta)
    ad = opt.adapt_delta;
    if ~(isnumeric(ad) && isscalar(ad) && isreal(ad)) || ~isfinite(ad) || ...
            ad <= 0 || ad >= 1
        error('psyrat_run:adapt_delta', ...
            ['Invalid ''adapt_delta'' value. It must lie strictly between 0 ' ...
            'and 1. How to fix: omit ''adapt_delta'' to keep each model''s ' ...
            'tuned default, or pass a value such as 0.95.']);
    end
end
if ~isempty(opt.max_treedepth)
    mt = opt.max_treedepth;
    if ~(isnumeric(mt) && isscalar(mt) && isreal(mt)) || ~isfinite(mt) || ...
            mt <= 0 || mod(mt,1) ~= 0
        error('psyrat_run:max_treedepth', ...
            ['Invalid ''max_treedepth'' value. It must be a positive integer. ' ...
            'How to fix: omit ''max_treedepth'' to keep CmdStan''s default of ' ...
            '10, or pass a value such as 12.']);
    end
end

%-------------------------------------------------------------------------
% 2. Resolve and validate the data source.
%-------------------------------------------------------------------------
hasFile = ~isempty(char(opt.file));
hasData = ~isempty(opt.data);
if hasFile == hasData
    error('psyrat_run:datasource', ...
        ['Provide exactly one data source. ' ...
        'How to fix: pass either ''file'',''/path/to/scores.csv'' or ' ...
        '''data'',aMatlabTable (not both, not neither).']);
end

if hasData
    rawtable = opt.data;
    rawfilename = 'psyrat_run_input';
    rawsourcepath = '';
else
    rawfile = char(opt.file);
    if exist(rawfile,'file') ~= 2
        error('psyrat_run:filenotfound', ...
            ['Input file not found: ''%s''. ' ...
            'How to fix: pass ''file'' with the full path to an existing ' ...
            'long-format data file.'],rawfile);
    end
    %read through the toolbox reader, as the GUI does (psyrat_startproc), so
    %the scripted route gets the same supported-extension check and delimiter
    %cycling chapter 6 documents; a bare readtable here had neither (audit
    %2026-09-05). For a plain comma-separated file the two are identical.
    rawtable = psyrat_readtable('file',rawfile);
    [~,fn,fe] = fileparts(rawfile);
    rawfilename = [fn fe];
    rawsourcepath = rawfile;
end

%-------------------------------------------------------------------------
% 3. Resolve and normalize the output location. Whitespace/existence are also
%    checked by the shared validator and the wrapper; a friendly up-front
%    error here saves the caller a deep stack trace.
%-------------------------------------------------------------------------
savepath = char(opt.savepath);
savename = char(opt.savename);
if isempty(savepath)
    error('psyrat_run:savepath', ...
        ['No ''savepath'' was provided. ' ...
        'How to fix: pass ''savepath'' with an existing directory (no spaces).']);
end
if isempty(savename)
    error('psyrat_run:savename', ...
        ['No ''savename'' was provided. ' ...
        'How to fix: pass ''savename'' such as ''study.psyrat'' (no spaces).']);
end
% add the .psyrat extension if the caller omitted any extension
[~,~,se] = fileparts(savename);
if isempty(se)
    savename = [savename '.psyrat'];
end

%-------------------------------------------------------------------------
% 4. Build psyrat_prefs from the toolbox defaults, then override the
%    processing fields the caller set. Starting from psyrat_defaults keeps the
%    CLI defaults identical to the GUI defaults (priors, traceplots, etc.).
%-------------------------------------------------------------------------
psyrat_prefs = psyrat_defaults();

psyrat_prefs.proc.nchains   = opt.chains;
psyrat_prefs.proc.nwarmup   = opt.warmup;
psyrat_prefs.proc.nsampling = opt.sampling;
psyrat_prefs.proc.niter     = opt.warmup + opt.sampling;
psyrat_prefs.proc.verbose   = opt.verbose;
psyrat_prefs.proc.traceplots = 1;            % trace-plot review is GUI-only
psyrat_prefs.proc.seed      = opt.seed;
% optional CmdStan sampler overrides; [] (the default) leaves estimation exactly
% as it was before these were exposed here.
psyrat_prefs.proc.adapt_delta   = opt.adapt_delta;
psyrat_prefs.proc.max_treedepth = opt.max_treedepth;
psyrat_prefs.proc.withinchain = opt.withinchain;
psyrat_prefs.proc.diffest   = opt.diffest;
% Residual covariance for co-occurring difference-score events. 'diffrescor'
% and 'diffwpcov' are one setting under two names, and BOTH must be written
% consistently here: psyrat_computevarcompwarp re-derives diffrescor from
% diffwpcov, so setting only diffrescor would be discarded downstream and the
% run would silently fit the uncorrelated-residual model. See the resolver.
[diffrescor,diffwpcov] = local_reconcile_rescor(opt,usingDefaults,cfgfilled);
psyrat_prefs.proc.diffrescor = diffrescor;
psyrat_prefs.proc.diffwpcov = diffwpcov;
% Reconcile the dynamic-reliability flag with the dimension column(s), matching
% the GUI coupling (psyrat_startproc.m:2176-2194): supplying a Dimension-1
% column turns dynamic reliability on, a Dimension-2 column requires Dimension-1,
% and dynrel=2 is meaningless without a dimension column. Without this a
% 'dim1col' passed without 'dynrel',2 would be silently ignored and a standard
% one-facet analysis would run instead.
if ~isempty(char(opt.dim1col))
    psyrat_prefs.proc.dynrel = 2;    % a Dimension-1 column forces dynamic reliability
elseif ~isempty(char(opt.dim2col))
    error('psyrat_run:dim2WithoutDim1', ...
        ['Dimension 2 (''%s'') was supplied without Dimension 1. ' ...
        'How to fix: pass ''dim1col'' with a Dimension-1 column, or drop ' ...
        '''dim2col''.'],char(opt.dim2col));
elseif opt.dynrel == 2
    error('psyrat_run:dynrelNeedsDim1', ...
        ['Dynamic reliability (''dynrel'',2) requires a between-person ' ...
        'dimension column, but none was supplied. ' ...
        'How to fix: pass ''dim1col'' with a continuous covariate column, or ' ...
        'set ''dynrel'',1.']);
else
    psyrat_prefs.proc.dynrel = opt.dynrel;
end
psyrat_prefs.proc.splits    = opt.splits;
psyrat_prefs.proc.ssrescor  = opt.ssrescor;
psyrat_prefs.proc.dispersion = opt.dispersion;
psyrat_prefs.proc.gammascale = opt.gammascale;
% observation-likelihood family, normalized to a lowercase char (the warp reads
% proc.family and forwards it to the engine + preflight)
psyrat_prefs.proc.family    = lower(char(opt.family));
if ~isempty(opt.priors)
    psyrat_prefs.proc.priors = opt.priors;
end

% estimation engine and its native-engine sub-settings. Leave the
% psyrat_defaults values in place unless the caller set them (so the default
% stays 0 = auto, with CmdStan preferred whenever present).
if ~isempty(opt.engine)
    psyrat_prefs.proc.engine = local_normalize_engine(opt.engine);
end
if ~isempty(opt.bootstrap_reps)
    psyrat_prefs.proc.bootstrap_reps = opt.bootstrap_reps;
end
if ~isempty(char(opt.hmc_gradient))
    psyrat_prefs.proc.hmc_gradient = char(opt.hmc_gradient);
end

% subject-level reliability: prefer 'ssubjrel', fall back to the legacy
% 'sserrvar' alias, then to the default (off). Keep both fields synchronized so
% the wrapper's legacy-reconciliation branch is a no-op.
if ~isempty(opt.ssubjrel)
    ssubjrel = opt.ssubjrel;
elseif ~isempty(opt.sserrvar)
    ssubjrel = opt.sserrvar;
else
    ssubjrel = 1;
end
psyrat_prefs.proc.ssubjrel = ssubjrel;
psyrat_prefs.proc.sserrvar = ssubjrel;

% difference-of-differences event mapping (only when diffest == 3)
if opt.diffest == 3 && ~isempty(opt.dodmap)
    psyrat_prefs.proc.dodmap = opt.dodmap;
    psyrat_prefs.proc.dodset = 1;
    psyrat_prefs.proc.dodmapcol = [];
end

%-------------------------------------------------------------------------
% 5. Build psyrat_data with the same fields the GUI populates, so the shared
%    GUI-mode psyrat_loadfile path applies the same column selection and
%    which* row-filtering. Optional facets are set to '' (= use all), matching
%    the GUI convention.
%-------------------------------------------------------------------------
psyrat_data = struct();
psyrat_data.raw.filename = rawfilename;
psyrat_data.raw.filepath = savepath;
psyrat_data.raw.data     = rawtable;
%full path of the input file, for the sidecar's file hash. raw.filename is the
%bare name the viewers display; hashing from it only worked when the file
%happened to sit in the current folder or on the MATLAB path, and could hash a
%same-named file elsewhere (audit 2026-09-05). Empty for an in-memory table.
psyrat_data.raw.sourcepath = rawsourcepath;

% Version stamp. The GUI path sets this from psyrat_prefs.ver in
% psyrat_startproc; the table export callbacks read psyrat_data.ver
% unguarded, so the scripted path must set it too or every Save Table /
% Save IDs press on a CLI-produced result errors.
psyrat_data.ver = psyrat_defineversion;

psyrat_data.proc.idheader    = char(opt.idcol);
psyrat_data.proc.measheader  = char(opt.meascol);
psyrat_data.proc.measheaders = {char(opt.meascol)};
psyrat_data.proc.eventheader = char(opt.eventcol);
psyrat_data.proc.groupheader = char(opt.groupcol);
psyrat_data.proc.timeheader  = char(opt.timecol);

% which* selections are read unconditionally by GUI-mode psyrat_loadfile, so
% they must exist; '' means "use every label found in the column".
psyrat_data.proc.whichevents = local_which(opt.whichevents);
psyrat_data.proc.whichgroups = local_which(opt.whichgroups);
psyrat_data.proc.whichtimes  = local_which(opt.whichtimes);

% dimension and split-weight headers are optional (isfield-guarded in
% psyrat_loadfile); only set them when the caller named a column.
if ~isempty(char(opt.dim1col))
    psyrat_data.proc.dim1header = char(opt.dim1col);
end
if ~isempty(char(opt.dim2col))
    psyrat_data.proc.dim2header = char(opt.dim2col);
end
if ~isempty(char(opt.weightcol))
    psyrat_data.proc.weightheader = char(opt.weightcol);
end

psyrat_data.proc.savepath = savepath;
psyrat_data.proc.savename = savename;

% interactive == false marks this as a headless/batch run: the wrapper skips
% all GUI prompts, auto-reruns to a cap on non-convergence, and errors (rather
% than silently saving a non-converged result) if it still does not converge.
psyrat_data.proc.batchmode   = ~logical(opt.interactive);
psyrat_data.proc.interactive = logical(opt.interactive);

%-------------------------------------------------------------------------
% 6. Prepare the analysis-ready table via the same loader the GUI uses, then
%    hand off to the shared estimation wrapper. The wrapper runs the shared
%    preflight validator, executes CmdStan, checks convergence, and saves the
%    .psyrat output to fullfile(savepath, savename).
%-------------------------------------------------------------------------
psyrat_data.proc.data = psyrat_loadfile('psyrat_prefs',psyrat_prefs, ...
    'psyrat_data',psyrat_data);

psyrat_data = psyrat_computevarcompwarp('psyrat_prefs',psyrat_prefs, ...
    'psyrat_data',psyrat_data);

if ~isfield(psyrat_data,'rel')
    error('psyrat_run:missingReliabilityOutput', ...
        ['No reliability output was generated for measurement ''%s''. ' ...
        'How to fix: review the CmdStan output above for convergence or ' ...
        'model errors.'],char(opt.meascol));
end

end

%-------------------------------------------------------------------------
% Local helpers
%-------------------------------------------------------------------------
function tf = local_ischarlike(x)
%Accept char row vectors and scalar strings (including '') for column names.
tf = ischar(x) || (isstring(x) && isscalar(x));
end

function w = local_which(x)
%Normalize a which* selection: empty -> '' (use all labels), else a cell array
%passed straight through to psyrat_loadfile.
if isempty(x)
    w = '';
else
    w = x;
end
end

function e = local_normalize_engine(x)
%Map an engine option to its numeric code (0=auto,1=cmdstan,2=fitlme,3=hmc).
%Accepts either the numeric code or a friendly name so the CLI can be called as
%'engine','cmdstan' or 'engine',1.
if isnumeric(x)
    e = x;
    if ~ismember(e,[0 1 2 3])
        error('psyrat_run:engine', ...
            ['Invalid engine code %g. Valid values: 0 (auto), 1 (cmdstan), ' ...
            '2 (fitlme), 3 (hmc).'],e);
    end
    return;
end
s = lower(strtrim(char(x)));
switch s
    case {'auto','0'}
        e = 0;
    case {'cmdstan','stan','1'}
        e = 1;
    case {'fitlme','lme','2'}
        e = 2;
    case {'hmc','hmcsampler','3'}
        e = 3;
    otherwise
        error('psyrat_run:engine', ...
            ['Invalid engine ''%s''. Valid values: auto, cmdstan, fitlme, hmc ' ...
            '(or the numeric codes 0/1/2/3).'],s);
end
end

function [opt,filled] = local_apply_config(opt,configfile,usingDefaults)
%Merge a run-config sidecar into the parsed options. Only options the caller
%left at their default (i.e. listed in usingDefaults) are filled from the
%config, so explicitly passed arguments always take precedence.
%
%Second output: the names of the options this function actually took from the
%config. Callers that must distinguish "the sidecar recorded this value" from
%"nobody set it, so it is the parser default" need that list, because the value
%alone cannot tell the two apart.
filled = {};
if exist(configfile,'file') ~= 2
    error('psyrat_run:confignotfound', ...
        ['Run-config file not found: ''%s''. ' ...
        'How to fix: pass ''config'' with the path to a *_runconfig.json ' ...
        'written by a previous run.'],configfile);
end
txt = fileread(configfile);
cfg = jsondecode(txt);

% Map each psyrat_run option name to its location in the sidecar
% (section, field). Options absent from the config are left untouched.
map = {
    'idcol','columns','idheader';
    'meascol','columns','measheader';
    'eventcol','columns','eventheader';
    'groupcol','columns','groupheader';
    'timecol','columns','timeheader';
    'weightcol','columns','weightheader';
    'dim1col','columns','dim1header';
    'dim2col','columns','dim2header';
    'whichevents','columns','whichevents';
    'whichgroups','columns','whichgroups';
    'whichtimes','columns','whichtimes';
    'chains','estimation','nchains';
    'warmup','estimation','nwarmup';
    'sampling','estimation','nsampling';
    'seed','estimation','seed';
    'adapt_delta','estimation','adapt_delta';
    'max_treedepth','estimation','max_treedepth';
    'verbose','estimation','verbose';
    'withinchain','estimation','withinchain';
    'engine','estimation','engine';
    'bootstrap_reps','estimation','bootstrap_reps';
    'hmc_gradient','estimation','hmc_gradient';
    'ssubjrel','design','ssubjrel';
    'diffest','design','diffest';
    'diffrescor','design','diffrescor';
    'diffwpcov','design','diffwpcov';
    'dynrel','design','dynrel';
    'splits','design','splits';
    'ssrescor','design','ssrescor';
    'dispersion','design','dispersion';
    'gammascale','design','gammascale';
    'family','design','family';
    'dodmap','design','dodmap';
    'priors','','priors';
    'savepath','output','savepath';
    'savename','output','savename';
    };

for k = 1:size(map,1)
    optname  = map{k,1};
    section  = map{k,2};
    cfgfield = map{k,3};

    % skip anything the caller set explicitly
    if ~ismember(optname,usingDefaults)
        continue;
    end

    % locate the value in the (optionally nested) config structure
    if isempty(section)
        if ~isfield(cfg,cfgfield)
            continue;
        end
        val = cfg.(cfgfield);
    else
        if ~isfield(cfg,section) || ~isfield(cfg.(section),cfgfield)
            continue;
        end
        val = cfg.(section).(cfgfield);
    end

    % normalize empty DoD mapping back to a cell array (jsondecode yields [])
    if strcmp(optname,'dodmap') && isempty(val)
        val = {};
    end

    opt.(optname) = val;
    filled{end+1} = optname; %#ok<AGROW>
end

% Warn about environment differences that break bit-reproducibility of the
% replay. Done here, after the merge, so the checks see the engine preference
% the replay will actually use.
local_check_replay_environment(cfg,opt);
end

function local_check_replay_environment(cfg,opt)
%Compare the machine replaying a run-config against the machine that wrote it.
%
%Both checks warn rather than error. Replaying a recorded configuration on a
%different machine is a legitimate thing to do (often the caller only wants the
%settings), so the right behavior is to proceed and say what changed. What is
%NOT acceptable is silence: psyrat_run advertises "same inputs and seed ->
%same variance components", and both differences below break that promise
%invisibly.
%
% 1. CmdStan version. The sidecar records the version DETECTED when the run was
%    written. NUTS draws are not guaranteed identical across CmdStan versions
%    at a fixed seed, so a version change can move every estimate slightly.
% 2. Engine actually used. The replay restores the engine PREFERENCE
%    (cfg.estimation.engine), never the engine that produced the recorded run
%    (cfg.estimation.engine_used). Under the default 'auto' preference the two
%    differ whenever the original run fell back off CmdStan, and the cascade is
%    re-derived from THIS machine's installed engines -- so an auto replay can
%    quietly switch between the Bayesian CmdStan fit and an approximate native
%    one.

%--- 1. CmdStan version --------------------------------------------------
recorded = '';
if isfield(cfg,'cmdstan_version') && ~isempty(cfg.cmdstan_version) && ...
        (ischar(cfg.cmdstan_version) || isstring(cfg.cmdstan_version))
    recorded = char(cfg.cmdstan_version);
end
if ~isempty(recorded)
    %detection fails open ('' = no readable CmdStan here); never let a probe of
    %the local install stop a replay from starting.
    try
        detected = char(psyrat_detect_cmdstan_version());
    catch
        detected = '';
    end
    if isempty(detected)
        warning('psyrat_run:cmdstanVersionMismatch', ...
            ['The recorded run used CmdStan %s, but no CmdStan version could ' ...
            'be detected on this machine. Estimates from this replay are not ' ...
            'guaranteed to match the recorded ones, even at the same seed. ' ...
            'How to fix: install CmdStan %s to reproduce the recorded draws, ' ...
            'or treat this run as a new analysis rather than a replication.'], ...
            recorded,recorded);
    elseif ~strcmp(detected,recorded)
        warning('psyrat_run:cmdstanVersionMismatch', ...
            ['The recorded run used CmdStan %s; this machine has CmdStan %s. ' ...
            'CmdStan draws are not guaranteed to be identical across versions ' ...
            'at the same seed, so this replay may not reproduce the recorded ' ...
            'estimates exactly. How to fix: install CmdStan %s to reproduce ' ...
            'the recorded draws, or report the version actually used.'], ...
            recorded,detected,recorded);
    end
end

%--- 2. Engine actually used ---------------------------------------------
if ~isfield(cfg,'estimation') || ~isstruct(cfg.estimation) || ...
        ~isfield(cfg.estimation,'engine_used') || isempty(cfg.estimation.engine_used)
    return;
end
used = lower(strtrim(char(cfg.estimation.engine_used)));
switch used
    case 'cmdstan'
        usedcode = 1;
    case 'fitlme'
        usedcode = 2;
    case 'hmc'
        usedcode = 3;
    otherwise
        return;  %unrecognized name: say nothing rather than guess
end

%the engine preference this replay will run under ([] = psyrat_defaults, auto)
try
    if isempty(opt.engine)
        req = 0;
    else
        req = local_normalize_engine(opt.engine);
    end
catch
    return;  %an invalid engine value is reported by the normalizer later
end

if req == 0
    %auto: CmdStan when present, so anything else means the recorded run fell back
    if usedcode ~= 1
        warning('psyrat_run:engineUsedNotRestored', ...
            ['The recorded run was produced by the ''%s'' engine, but the ' ...
            'run-config restores the engine PREFERENCE (auto), not the engine ' ...
            'that ran. This replay re-derives the cascade from the engines ' ...
            'installed here, so it may use a different one. How to fix: pass ' ...
            '''engine'',''%s'' to force the recorded engine.'],used,used);
    end
elseif usedcode ~= req
    warning('psyrat_run:engineUsedNotRestored', ...
        ['The recorded run was produced by the ''%s'' engine, but this replay ' ...
        'is forced to engine %d. Results will not match the recorded ones. ' ...
        'How to fix: pass ''engine'',''%s'' to force the recorded engine, or ' ...
        'accept the difference deliberately.'],used,req,used);
end
end

function local_check_integer(v,name,minval,howtofix)
%Shared bound check for the integer estimation settings, so chains/warmup/
%sampling/seed all report the same shape of message.
ok = isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && ...
    mod(v,1) == 0 && v >= minval;
if ok
    return;
end
if isnumeric(v) && isscalar(v)
    shown = num2str(v);
else
    shown = ['a ' class(v)];
end
error(['psyrat_run:' name], ...
    ['Invalid ''%s'' value (%s). It must be an integer of %g or higher. ' ...
    'How to fix: %s.'],name,shown,minval,howtofix);
end

function [diffrescor,diffwpcov] = local_reconcile_rescor(opt,usingDefaults,cfgfilled)
%Resolve the two names for ONE setting: whether the within-person residual
%covariance between co-occurring difference-score events is estimated.
%
%'diffrescor' and 'diffwpcov' are not two independent knobs. The GUI shows a
%single control ("Estimate residual covariance for co-occurring events") and
%writes both fields from it, and psyrat_computevarcomp's error messages name
%them together as "''diffrescor''=2 / ''diffwpcov''=2".
%
%That engine reconciliation is NOT symmetric, which is why the resolution below
%happens here rather than being left to it: psyrat_computevarcomp resolves a
%LONE diffrescor=2 by promoting diffwpcov to 2, but when BOTH names are supplied
%and disagree it follows diffwpcov and clamps diffrescor with no warning. The
%same precedence is stated at psyrat_startproc's ssrescor gate. Reconciling here
%means the engine never sees a disagreeing pair from this entry point.
%
%The reason this resolver exists: psyrat_computevarcompwarp re-derives
%diffrescor FROM diffwpcov, and psyrat_run always writes proc.diffwpcov. So
%before this, psyrat_run(...,'diffrescor',2) without 'diffwpcov',2 was reset to
%1 on the way to the engine and the run fit the UNCORRELATED-residual model
%with no runtime signal -- every difference-score coefficient, SEM and
%cut-score computed under the wrong residual covariance structure. Writing both
%fields consistently here makes that clamp a no-op.
%
%Precedence follows the documented run-config rule: an explicit argument beats
%a value restored from a sidecar. Two explicit arguments that disagree are an
%error, not a coin flip.

rescor = opt.diffrescor;
wpcov  = opt.diffwpcov;
local_check_onoff('diffrescor',rescor);
local_check_onoff('diffwpcov',wpcov);

%where each value came from: an explicit argument, a run-config sidecar, or
%neither (the parser default). The value alone cannot say.
argRescor = ~ismember('diffrescor',usingDefaults);
argWpcov  = ~ismember('diffwpcov',usingDefaults);
setRescor = argRescor || ismember('diffrescor',cfgfilled);
setWpcov  = argWpcov  || ismember('diffwpcov',cfgfilled);

if setRescor && setWpcov && rescor ~= wpcov
    if argRescor && argWpcov
        error('psyrat_run:diffrescorConflict', ...
            ['Contradictory residual-covariance settings: ''diffrescor'',%d ' ...
            'and ''diffwpcov'',%d. These are one setting under two names ' ...
            '(1 = residual covariance fixed to 0, 2 = estimated), and the two ' ...
            'values fit structurally different models, so PsyRAT will not ' ...
            'pick one for you. How to fix: pass only ''diffrescor'',%d (the ' ...
            'other name follows automatically), or set both to the same ' ...
            'value.'],rescor,wpcov,rescor);
    elseif argRescor
        warning('psyrat_run:diffrescorConfigOverride', ...
            ['The run-config recorded diffwpcov=%d, but ''diffrescor'',%d was ' ...
            'passed explicitly. These are one setting under two names, so the ' ...
            'explicit argument wins and diffwpcov is set to %d.'], ...
            wpcov,rescor,rescor);
        wpcov = rescor;
    elseif argWpcov
        warning('psyrat_run:diffrescorConfigOverride', ...
            ['The run-config recorded diffrescor=%d, but ''diffwpcov'',%d was ' ...
            'passed explicitly. These are one setting under two names, so the ' ...
            'explicit argument wins and diffrescor is set to %d.'], ...
            rescor,wpcov,wpcov);
        rescor = wpcov;
    else
        %both came from the sidecar and they disagree: the recorded
        %configuration is internally inconsistent (a hand-edited or corrupted
        %file), so which model was fit cannot be recovered from it.
        error('psyrat_run:diffrescorConfigConflict', ...
            ['The run-config records contradictory residual-covariance ' ...
            'settings (diffrescor=%d, diffwpcov=%d). These are one setting ' ...
            'under two names, so a sidecar written by PsyRAT always records ' ...
            'them equal; this one cannot say which model was fit. How to fix: ' ...
            'correct the sidecar, or pass ''diffrescor'' (or ''diffwpcov'') ' ...
            'explicitly to override it.'],rescor,wpcov);
    end
elseif setRescor && ~setWpcov
    wpcov = rescor;      %the fix: a lone diffrescor now reaches estimation
elseif setWpcov && ~setRescor
    rescor = wpcov;
end

diffrescor = rescor;
diffwpcov  = wpcov;
end

function local_check_onoff(name,v)
%Both names of the residual-covariance setting are 1 (off) or 2 (on); reject
%anything else here so a typo cannot fall through to the engine as a silent 1.
if ~(isnumeric(v) && isscalar(v) && any(v == [1 2]))
    if isnumeric(v) && isscalar(v)
        shown = num2str(v);
    else
        shown = ['a ' class(v)];
    end
    error(['psyrat_run:' name], ...
        ['Invalid ''%s'' value (%s). It must be 1 (residual covariance fixed ' ...
        'to 0) or 2 (residual covariance estimated). How to fix: pass ' ...
        '''%s'',1 or ''%s'',2, or omit it for the default (1).'], ...
        name,shown,name,name);
end
end
