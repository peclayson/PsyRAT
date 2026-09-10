function sidecarfile = psyrat_write_runconfig(psyrat_data,psyrat_prefs)
%Write a human-readable run-configuration sidecar next to a saved .psyrat file.
%
%The sidecar is a JSON file recording the seed, estimation settings, design
%flags, column mapping, viewer/decision settings, and priors used for a run. It
%makes any analysis (GUI or command line) reproducible: the values can be read
%by eye for a methods section, and psyrat_run('config',sidecarfile,...) can
%replay the settings on the same (or new) data. The viewer block is record-only
%(psyrat_run exposes no options for it).
%
%psyrat_write_runconfig(psyrat_data,psyrat_prefs)
%
%Input
% psyrat_data  - toolbox data structure after estimation (uses .proc and,
%   when present, .rel.analysis).
% psyrat_prefs - toolbox preferences (uses .proc for settings and priors, and
%   .view for the viewer/decision block; a prefs struct without .view falls back
%   to the psyrat_defaults values).
%
%Output
% sidecarfile - full path to the JSON sidecar that was written. The file is
%   named <savename-without-extension>_runconfig.json and placed in savepath.

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

proc = psyrat_data.proc;

% Resolve the sidecar path from the .psyrat output location.
[~,base,~] = fileparts(proc.savename);
sidecarfile = fullfile(proc.savepath,[base '_runconfig.json']);

cfg = struct();
% Schema version. 1.1.0 added estimation.adapt_delta, estimation.max_treedepth,
% view.reltype and view.diffgcoeff. The bump is a MINOR one because the change
% is purely additive: the reader (psyrat_run local_apply_config) walks a fixed
% option map and skips any section or field it does not recognize, so a 1.0.0
% sidecar still loads and replays unchanged, and nothing in the toolbox branches
% on this string.
% 1.2.0 added source_revision; additive, same reader-skips-unknown contract.
cfg.schema = 'psyrat_runconfig/1.2.0';
cfg.created_utc = local_utc_now();
cfg.psyrat_version = local_safe_version();

% Source revision (best effort): git short hash of the toolbox checkout with
% '-dirty' appended when the tree has local changes; '' when the toolbox is
% not a git checkout (every released copy) or git is unavailable. The version
% string alone is shared by every commit between releases; this field is what
% lets a sidecar identify the exact code state.
cfg.source_revision = local_safe_source_revision();
cfg.matlab_version = version();

% CmdStan version actually detected in this environment (best effort; '' when
% CmdStan is not installed or its version cannot be read). Recorded for every
% run regardless of engine, so a native-engine run still documents whether a
% CmdStan was available. This is the DETECTED version, not the recommended
% version pinned in psyrat_dependentsversions.
cfg.cmdstan_version = local_safe_cmdstan_version();

if isfield(psyrat_data,'rel') && isfield(psyrat_data.rel,'analysis')
    cfg.analysis = psyrat_data.rel.analysis;
else
    cfg.analysis = '';
end

% Output location (lets psyrat_run('config',...) and the reader find the
% companion .psyrat result).
cfg.output = struct('savepath',proc.savepath,'savename',proc.savename);

% Column mapping (the user's source column names). Optional facets are stored
% as '' when absent, matching the psyrat_run/loadfile convention.
cfg.columns = struct();
cfg.columns.idheader     = local_field(proc,'idheader','');
cfg.columns.measheader   = local_field(proc,'measheader','');
cfg.columns.eventheader  = local_field(proc,'eventheader','');
cfg.columns.groupheader  = local_field(proc,'groupheader','');
cfg.columns.timeheader   = local_field(proc,'timeheader','');
cfg.columns.weightheader = local_field(proc,'weightheader','');
cfg.columns.dim1header   = local_field(proc,'dim1header','');
cfg.columns.dim2header   = local_field(proc,'dim2header','');
cfg.columns.whichevents  = local_field(proc,'whichevents','');
cfg.columns.whichgroups  = local_field(proc,'whichgroups','');
cfg.columns.whichtimes   = local_field(proc,'whichtimes','');

% Input-data provenance: a content hash of the loaded long-format table plus,
% when the run loaded from a file on disk, the source file name and a hash of
% that file. The table hash (psyrat_hash_struct) is deterministic and present
% even for in-memory CLI runs (psyrat_run('data',table,...)); the file hash
% (psyrat_hash_file) ties the run to the exact source file when one exists.
% nrows is the number of single-trial rows that entered estimation. Every piece
% is best-effort: a missing/unhashable piece is recorded as '' (or [] for nrows).
cfg.data = local_data_provenance(psyrat_data);

% Estimation settings (the reproducibility-critical block).
pp = psyrat_prefs.proc;
cfg.estimation = struct();
cfg.estimation.seed        = local_field(pp,'seed',12345);
cfg.estimation.nchains     = local_field(pp,'nchains',[]);
cfg.estimation.nwarmup     = local_field(pp,'nwarmup',[]);
cfg.estimation.nsampling   = local_field(pp,'nsampling',[]);
cfg.estimation.niter       = local_field(pp,'niter',[]);
cfg.estimation.withinchain = local_field(pp,'withinchain',1);
% The RESOLVED within-chain thread count (1 when threading was off, since
% psyrat_computevarcomp initializes REL.threads_per_chain = 1 and overwrites it
% only when threading is enabled; empty only for a result saved before the
% field existed). withinchain above is the user's 1/2 request; the threads per chain
% the core derived from the machine's core count set the reduce_sum grainsize
% and therefore the partition of the log-density sum, so a replay on another
% machine cannot reproduce a threaded run without it (audit 2026-09-05).
cfg.estimation.threads_per_chain = [];
if isfield(psyrat_data,'rel') && isstruct(psyrat_data.rel) && ...
        isfield(psyrat_data.rel,'threads_per_chain') && ...
        ~isempty(psyrat_data.rel.threads_per_chain)
    cfg.estimation.threads_per_chain = psyrat_data.rel.threads_per_chain;
end
cfg.estimation.verbose     = local_field(pp,'verbose',[]);

% Optional CmdStan sampler overrides. These belong in the reproducibility block,
% not with tuning trivia: both change NUTS step-size adaptation during warmup,
% so a run that used them does NOT produce the same draws at the same seed as
% one that did not. A sidecar that omitted them let a replay resample different
% draws while reporting the same seed. Empty ([]) is the documented default
% (each model's tuned adapt_delta; CmdStan's max_depth = 10) and is recorded as
% empty, so a default run's sidecar replays as a default run.
cfg.estimation.adapt_delta   = local_field(pp,'adapt_delta',[]);
cfg.estimation.max_treedepth = local_field(pp,'max_treedepth',[]);

% Estimation engine: the requested preference (0 auto / 1 cmdstan / 2 fitlme /
% 3 hmc) plus the native-engine knobs, so a replay (psyrat_run 'config')
% reproduces the same engine choice. engine_used records which engine actually
% produced the run, which can differ from the preference when the auto cascade
% fell back off CmdStan.
cfg.estimation.engine         = local_field(pp,'engine',0);
cfg.estimation.bootstrap_reps = local_field(pp,'bootstrap_reps',[]);
cfg.estimation.hmc_gradient   = local_field(pp,'hmc_gradient','auto');
if isfield(psyrat_data,'rel') && isstruct(psyrat_data.rel) && ...
        isfield(psyrat_data.rel,'engine') && ~isempty(psyrat_data.rel.engine)
    cfg.estimation.engine_used = psyrat_data.rel.engine;
else
    cfg.estimation.engine_used = 'cmdstan';
end

% Design / analysis flags.
cfg.design = struct();
cfg.design.ssubjrel   = local_field(pp,'ssubjrel',1);
cfg.design.diffest    = local_field(pp,'diffest',1);
cfg.design.diffrescor = local_field(pp,'diffrescor',1);
cfg.design.diffwpcov  = local_field(pp,'diffwpcov',1);
cfg.design.dynrel     = local_field(pp,'dynrel',1);
cfg.design.splits     = local_field(pp,'splits',1);
cfg.design.ssrescor   = local_field(pp,'ssrescor',1);
% Copula residual-dispersion parameterization for the group-level concurrent gamma
% difference (1 = per-person log-nu, 2 = global fixed dispersion, [] = unset ->
% engine context-default); recorded as-is so a replay reproduces the same estimand
% (a run that relied on the default records [] and re-derives the default on replay,
% rather than being pinned to a concrete value).
cfg.design.dispersion = local_field(pp,'dispersion',[]);
% Gamma scale parameterization for the one-facet dynamic designs (1 = log-nu
% dispersion, 2 = log residual SD, [] = unset -> the engine's CONTEXT-AWARE default,
% which since rulings 33-36 is location-scale on the designs that support it and
% log-nu elsewhere; it was a flat 1 before 2026-08-07). Recorded for
% the same reason as family above, and with more force: 1 and 2 are DIFFERENT
% ESTIMANDS whose variance components are not comparable, so a replay that lost this
% would silently produce different numbers from the run it claims to reproduce. The
% warp pins the RESOLVED value here after estimation.
cfg.design.gammascale = local_field(pp,'gammascale',[]);
% Observation-likelihood family ('gaussian' default / 'gamma'); recorded so a
% replay reproduces the same estimand (gamma uses a log-linked mean, so a run
% replayed as Gaussian would give different variance components).
cfg.design.family     = local_field(pp,'family','gaussian');
cfg.design.dodmap     = local_field(pp,'dodmap',{});

% Viewer/decision settings. These live in psyrat_prefs.VIEW (not .proc) and they
% pick the estimand the output is reported under, so a sidecar without them
% documents the fit but not the numbers a reader was looking at:
%  gcoeff    - 1 dependability (absolute error) / 2 generalizability (relative)
%  reltype   - test-retest coefficient type: 1 equivalence (CE) / 2 stability
%              (CS) / 3 both (CES). This selects WHICH variance terms enter the
%              coefficient, and the three are not near-neighbors -- on one
%              dataset the same components gave 0.40 / 0.35 / 0.24.
%  diffgcoeff- the difference-score error type (1 absolute / 2 relative), which
%              REPLACES gcoeff on the difference-score routes: psyrat_relsummary
%              assigns "gcoeff = diffgcoeff" on the group-level route and passes
%              'gcoeff',diffgcoeff into psyrat_ssrel on the subject-level one
%  noccmode  - test-retest occasion convention (1 = default n'_o = 1)
%  nocc      - explicit occasion n' when noccmode overrides the default ([] = unset)
%  depvalue  - reliability cutoff on the ONE-FACET viewer
%              (psyrat_startview_sing:199): the threshold the trial-count
%              cutoff and the include/exclude split that follows are read off
%  relvalue  - the same cutoff on the TEST-RETEST viewer
%              (psyrat_startview_trt:154)
%
% Neither cutoff names a coefficient of its own. They are paired by VIEWER, not
% by estimand: each is compared against whichever coefficient that screen's
% selector picked. On the plain routes that is gcoeff, which defaults to 1, so
% by default BOTH are cutoffs on dependability. On the difference-score routes
% it is diffgcoeff instead: psyrat_relsummary assigns "gcoeff = diffgcoeff" on
% the group-level route and passes 'gcoeff',diffgcoeff into psyrat_ssrel on the
% subject-level route, then applies depcutoff to what came back (cited by
% symbol, not line -- that file moves). Read depvalue/relvalue together with the
% recorded coefficient selectors; neither is a dependability or a
% generalizability threshold on its own.
%
% reltype and diffgcoeff are recorded below. Their earlier omission
% (PSYRAT_AUDIT_FINDINGS.md section F, CC2) was a MACHINE-REPLAY gap only, not a
% reporting one: an exported table already names the estimand a human is looking
% at -- psyrat_provenance_lines emits a "Coefficient:" header from
% psyrat_trt_coeflabel, which reads reltype and adds the difference-score error
% type when diffgcoeff differs from gcoeff, and psyrat_variancet suffixes the
% difference-score columns from relsummary.diffgcoeff_name. What a program
% reading only this sidecar could not previously recover was which coefficient
% those numbers were.
%
% Recorded as-is, with the psyrat_defaults values as fallbacks. Two honest
% caveats: (1) this sidecar is written at estimation time (see
% psyrat_computevarcompwarp), so these are the viewer settings IN EFFECT for the
% run, not a capture of a later viewer session -- a user who changes gcoeff in
% the viewer after the fact is not re-recorded here; (2) psyrat_run has no
% options for them, so they are record-only and the config replay in
% local_apply_config ignores them by construction (it walks a fixed option map).
% Additive: the reader skips any section/field it does not map, so a sidecar
% written before this block still loads and replays unchanged.
if isfield(psyrat_prefs,'view') && isstruct(psyrat_prefs.view)
    pv = psyrat_prefs.view;
else
    pv = struct();
end
cfg.view = struct();
cfg.view.gcoeff     = local_field(pv,'gcoeff',1);
cfg.view.reltype    = local_field(pv,'reltype',1);
cfg.view.diffgcoeff = local_field(pv,'diffgcoeff',1);
cfg.view.noccmode   = local_field(pv,'noccmode',1);
cfg.view.nocc       = local_field(pv,'nocc',[]);
cfg.view.depvalue   = local_field(pv,'depvalue',.8);
cfg.view.relvalue   = local_field(pv,'relvalue',.8);

% Priors (whole structure, so the exact prior configuration is recorded).
cfg.priors = local_field(pp,'priors',struct());

% Encode (pretty-print when the MATLAB release supports it) and write.
try
    txt = jsonencode(cfg,'PrettyPrint',true);
catch
    txt = jsonencode(cfg);
end

fid = fopen(sidecarfile,'w');
if fid == -1
    error('psyrat_write_runconfig:fopen', ...
        'Could not open run-config sidecar for writing: %s',sidecarfile);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid,'%s',txt);

end

%-------------------------------------------------------------------------
% Local helpers
%-------------------------------------------------------------------------
function v = local_field(s,name,default)
%Return s.(name) if present and non-empty, else the default.
if isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function s = local_utc_now()
%ISO-8601 UTC timestamp, with a best-effort fallback on older releases.
try
    s = char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
catch
    s = datestr(now,'yyyy-mm-ddTHH:MM:SS');
end
end

function v = local_safe_version()
%PsyRAT version string, tolerant of psyrat_defineversion being unavailable.
try
    v = psyrat_defineversion();
catch
    v = '';
end
end

function rev = local_safe_source_revision()
%Git short hash (+ '-dirty') of the toolbox checkout; '' outside a checkout.
rev = '';
try
    tooldir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    [st, out] = system(sprintf('git -C "%s" rev-parse --short HEAD', tooldir));
    if st ~= 0
        return;
    end
    rev = strtrim(out);
    [st2, out2] = system(sprintf('git -C "%s" status --porcelain', tooldir));
    if st2 == 0 && ~isempty(strtrim(out2))
        rev = [rev '-dirty'];
    end
catch
    rev = '';
end
end

function v = local_safe_cmdstan_version()
%Detected CmdStan version string, or '' when CmdStan is absent/unreadable.
%psyrat_detect_cmdstan_version already fails open (returns '' rather than
%erroring); the extra try/catch guards against an unexpected throw so writing
%the sidecar can never break a finished run.
try
    v = psyrat_detect_cmdstan_version();
    if isempty(v)
        v = '';
    end
catch
    v = '';
end
end

function d = local_data_provenance(psyrat_data)
%Best-effort provenance for the input data: source file name + file hash (when
%the run loaded from a file that still exists), a deterministic content hash of
%the loaded table, and the row count. Any piece that cannot be produced is left
%as its empty default so the sidecar is always written.
d = struct('source_file','','source_path','','file_sha256','','table_sha256','','nrows',[]);

if ~isfield(psyrat_data,'raw') || ~isstruct(psyrat_data.raw)
    return;
end
raw = psyrat_data.raw;

% Source file name and file hash (file hash only when the path exists on disk;
% in-memory CLI runs have no source file). raw.filename is the bare name the
% screens display; both producers (psyrat_run, psyrat_startproc) also record
% raw.sourcepath, the full path of the input file, and the hash is taken from
% that. Hashing from the bare name only worked when the file sat in the
% current folder or on the MATLAB path, and would hash a same-named file found
% there instead (audit 2026-09-05). A struct without sourcepath keeps the old
% best-effort behavior.
if isfield(raw,'filename') && ~isempty(raw.filename) && ...
        (ischar(raw.filename) || isstring(raw.filename))
    d.source_file = char(raw.filename);
    hashsrc = d.source_file;
    if isfield(raw,'sourcepath') && ~isempty(raw.sourcepath) && ...
            (ischar(raw.sourcepath) || isstring(raw.sourcepath))
        d.source_path = char(raw.sourcepath);
        hashsrc = d.source_path;
    end
    if exist(hashsrc,'file') == 2
        try
            d.file_sha256 = psyrat_hash_file(hashsrc);
        catch
            d.file_sha256 = '';
        end
    end
end

% Loaded long-format table: content hash + row count. Robust for in-memory runs.
if isfield(raw,'data') && ~isempty(raw.data)
    try
        d.table_sha256 = psyrat_hash_struct(raw.data);
    catch
        d.table_sha256 = '';
    end
    try
        d.nrows = height(raw.data);
    catch
        try
            d.nrows = size(raw.data,1);
        catch
            d.nrows = [];
        end
    end
end
end
