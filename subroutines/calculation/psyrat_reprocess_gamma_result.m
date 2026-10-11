function [updated, record] = psyrat_reprocess_gamma_result(source, varargin)
%PSYRAT_REPROCESS_GAMMA_RESULT Explicit, versioned reprocessing of a SAVED gamma
% reliability result to the gamma_target_version 2 residual contract
% (Q8/PBA-08 and Q9/GUI-1 owner rulings, 2026-09-13), with the original
% preserved and a migration record returned and stored.
%
%   [updated, record] = psyrat_reprocess_gamma_result(source, ...
%       'output', freshPath, 'surface_target', 'population_marginal'|'typical_person')
%
% source
%   - a full psyrat_data struct (fields rel/proc/...), or
%   - a bare REL struct (fields analysis/family/out), or
%   - a .psyrat or .mat file path holding the variable psyrat_data (the file
%     contract psyrat_startview loads with load(file,'-mat')).
%   Summary-only exports (.xlsx/.xls/.csv/.txt/.tsv) are REFUSED: they carry
%   no posterior draws, and no reliability can be re-derived from them.
% 'output'  (optional for struct input, REQUIRED for file input)
%   A fresh .psyrat/.mat path. The source file is never overwritten and an
%   existing destination is never replaced. Every other variable in the source
%   file is carried into the destination unchanged. Struct input without
%   'output' writes nothing; the result is only returned. A bare REL source
%   cannot be written (no wrapper to preserve) and refuses 'output'.
% 'surface_target'  (dynamic non-difference gamma routes only)
%   'population_marginal' (default) or 'typical_person'. If the source already
%   carries an explicit surface_target and none is supplied, it is kept.
%   The historical summary option marginal = 0 (and the ignored marginal = 1)
%   both denote the OLD population surface; that flag is never translated
%   into a target. Static analysis 25 has no surface target ('not_applicable').
%
% updated  - the same kind of object as source, with REL.out.gamma_target
%            stamped and, for analysis 25, the residual encoding reprocessed.
%            record is also stored at REL.out.gamma_target.migration_history(end).
% record   - migration record: version_from/to, status, source, output, time,
%            analysis, parameterization, method, strata, fields read/written/
%            derived/absent, per-stratum transforms (compatibility statuses per
%            pair with counts and margins), count_source, id_mapping,
%            invalidation, limitations, notes.
%
% VERSIONS. The only stamp version this tool accepts is EXACTLY 2. An absent
% stamp is the historical version-1 (observation-only) encoding, the only
% legacy form ever produced. Any other explicit version is refused; no other
% legacy contract is invented.
%
% WHAT IS REPROCESSED
% Static analysis 25 (trt_sserrvar), version 1 -> 2. Version 1 pop_sdlog /
% ind_sdlog encode only the observation variances (C_typ, C_p). Version 2
% encodes E = R + C with R the population highest-order (person x occasion x
% trial) mean-surface remainder. R is reconstructed PER STRATUM AND PER DRAW
% from the retained observed-scale quantities (no refit, no person effects),
% in the normalized form of the independent reference (the frozen 2C
% static-25 migration reference, a development record):
%   a = (sd_P/mu)^2, b = (sd_O/mu)^2, c = (sd_I/mu)^2,
%   h_PO = (sd_PO/mu)^2 (h_PI, h_OI likewise),
%   q_P = a, q_O = b, q_I = c,
%   q_PO = (h_PO - a b) / ((1+a)(1+b))   (q_PI, q_OI likewise),
%   v_A = log1p(q_A),  alpha = log(mu) - sum(v)/2,
%   R = psyrat_gamma_mean_remainder_trt(alpha, sqrt(v))   (support-union sum)
% with sd_P = gro_sds(:,1), sd_O = sig_occ, sd_I = sig_trl, sd_PO = sig_occxid,
% sd_PI = sig_trlxid, sd_OI = sig_trlxocc (the extractor's slot assignment).
% Then E_typ = R + exp(2 pop_v1), E_p = R + exp(2 (pop_v1 + ind_v1)),
% pop_v2 = .5 log E_typ, ind_v2 = .5 log E_p - pop_v2. Column order of the
% [ndraws x NSUB] blocks (id2) and the id_matches row order are preserved.
%
% COMPONENT COMPATIBILITY (checked BEFORE any negative inverse is touched).
% Each pair deficit h - a b is classified against the declared per-draw
% operation-norm bound of the frozen reference
% (science/frozen_inverse_compatibility.json, "declared_error_model"):
%   u = 2^-53, gamma_n = n u / (1 - n u),
%   condition = 1 + 2|log mu| + sum |v_raw|   (v_raw: unclamped inverse),
%   eta = 2 gamma_32 condition,
%   bound = eta S + |b| eta (1+a) + |a| eta (1+b) + eta^2 (1+a)(1+b)
%           + gamma_3 (|h| + |a b|),   S = (1+a+b+h) + (1+a) + (1+b) + 1.
%   deficit >  bound : compatible_positive   (q_pair retained, including tiny positives)
%   deficit < -bound : materially_incompatible -> REFUSED, naming the component,
%                      stratum, draw, deficit and bound
%   otherwise        : roundoff_unresolved; a negative deficit is canonicalized
%                      to zero and RECORDED as unresolved legacy precision (the
%                      reconstruction is not claimed exact for that draw)
% The model assumes finite normal binary64 arithmetic with at most one ulp per
% elementary operation; it is a conservative proposal grounded in the legacy
% producer, not a certified MATLAB/libm bound. Overflow/underflow (non-finite
% ratios, B, R or E) is refused as an unsupported numeric regime; nothing is
% fabricated. When REL.out.gamma_disp carries the raw log-scale pieces of the
% LAST stratum (alpha, s_o, s_t, s_po, s_pt, s_ot), the normalized components
% they imply are compared with the stored ones under twice the same bound and
% a discrepancy is refused.
%
% Already version 2 (static): every stratum is validated (alignment, IDs,
% finite positive C_typ/C_p, finite nonnegative R, finite positive E, finite
% encodings that decode to finite positive values) and the encoding identity
% |log(R + C) - 2 (pop + ind)| <= gamma_8 (1 + |log(R + C)| + 2|pop| + 2|ind|)
% is checked elementwise (operation-count bound for the producer's sum, log,
% scaling and subtraction and the validator's sum, log and re-summation; the
% pop/ind magnitudes enter because their sum cancels; no absolute floor).
% Numerics are returned byte-equal; nothing is pooled twice.
%
% NUMERIC REGIME. The error model assumes finite NORMAL binary64 arithmetic.
% Every stored input, its square, the normalized ratios, the prefactor
% B = mu^2, the reconstructed q and v, the decoded C, the residual E and the
% encodings must be finite normal (|x| in [realmin, realmax]); a component
% may be exactly zero only where it is structurally absent, i.e. where the
% SOURCE value is zero or a negative pair deficit inside the bound was
% canonicalized. A strictly positive source that maps to zero, a subnormal or
% Inf under any transform (square, normalization by mu, pair quotient, log1p,
% the remainder) is positive-to-zero underflow / overflow and is refused
% (unsupported_numeric_regime), never carried as a structural zero. Every positive
% remainder term (a product of q over a factor subset with support union
% {P,O,I}) is pre-checked in log space, including its left-to-right partial
% products and its product with B, so a tiny-but-normal q whose product would
% underflow inside the remainder helper is refused (unsupported_numeric_regime)
% instead of yielding a truncated R. A complete zero-facet R = 0 stays valid.
%
% STORAGE CONTRACT. Per-stratum draw blocks that the summary/table/provenance
% consumers read through cell2mat (gro_sds, b, b_nu/b_sigma, chol_corrmat,
% ind_bs, ind_nu, ind_sdlog) must be stored exactly as the estimator stores
% them (num2cell blocks of scalar draws). Numeric inner blocks are refused
% (unsupported_storage): a migrated result that the consumers cannot open is
% not a valid migration, and nothing is normalized silently.
%
% PROVENANCE. An ABSENT gamma_target field is version 1. A present but empty
% or malformed gamma_target is refused (partial provenance), and a gamma_resid
% block without a stamp is refused as inconsistent provenance (that schema
% exists only under version 2), never reinterpreted as observation-only.
%
% Dynamic non-difference gamma routes (ic_dynrel 11, ic_dynrel_sserrvar 26,
% ic_dynrel_trt 14, ic_dynrel_sserrvar_trt 27): raw draws are validated in
% full (finite population draws and slopes, nonnegative SDs, per-stratum
% labels and counts, ndim/dim_names/dimz/slope-basis consistency, complete
% draw-aligned Cholesky correlation draws for BOTH parameterizations and BOTH
% surface choices, own effects and ssinfo id/id2/z/count mapping for subject
% routes) and a versioned gamma_target is stamped. The raw draws are not
% changed and no refit is needed, but the stamp decides what the next summary
% derives from them: a stamped analysis-27 result is summarized under the
% version-2 remainder scheme, so its per-participant coefficients (and their
% comparison with any cutoff; no inclusion flag is stored) can change, and an
% explicit surface_target, on a first migration or a restamp, can change the
% dynamic surface of any of these designs.
%
% WHAT IS NOT DONE. No person effects are reconstructed (u_p / v_p are not
% identifiable from C_p). Raw scale draws that version 1 did not retain
% (log-nu Intercept_nu / ind_nu under analysis 25) are recorded as absent,
% never manufactured; under location-scale they are exact (pop_v1 =
% Intercept_sigma, ind_v1 = ind_sd). Gaussian, difference, copula, DoD, splits
% and the static one-facet gamma design are refused. A cached
% psyrat_data.relsummary derived under version 1 is removed from the NEW
% result (recorded). Nothing migrates automatically on open.
%
% See also PSYRAT_GAMMA_MEAN_REMAINDER_TRT, PSYRAT_RELSUMMARY, PSYRAT_DYNREL_SUMMARY.

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

ID = 'psyrat_reprocess_gamma_result';
TARGET_VERSION = 2;

% ------------------------------------------------------------------ options
outputPath = '';
surfaceTargetIn = '';
k = 1;
while k <= numel(varargin)
    if ~ischar(varargin{k}) && ~isstring(varargin{k})
        error([ID ':bad_option'], 'Options must be name/value pairs; argument %d is not a name.', k);
    end
    name = lower(char(varargin{k}));
    if k == numel(varargin)
        error([ID ':bad_option'], 'Option ''%s'' has no value.', name);
    end
    val = varargin{k+1};
    switch name
        case 'output'
            outputPath = char(val);
        case 'surface_target'
            surfaceTargetIn = lower(char(val));
            if ~ismember(surfaceTargetIn, {'population_marginal', 'typical_person'})
                error([ID ':bad_surface_target'], ...
                    'surface_target must be ''population_marginal'' or ''typical_person''; got ''%s''.', surfaceTargetIn);
            end
        otherwise
            error([ID ':bad_option'], 'Unknown option ''%s''. Accepted: ''output'', ''surface_target''.', name);
    end
    k = k + 2;
end

% ------------------------------------------------------------------ source
src = struct('kind', '', 'path', '', 'description', '', 'file_vars', {{}});
fileVars = struct();
if ischar(source) || isstring(source)
    src.kind = 'file';
    src.path = char(source);
    [~, ~, ext] = fileparts(src.path);
    if ismember(lower(ext), {'.xlsx', '.xls', '.csv', '.txt', '.tsv'})
        error([ID ':summary_only_input'], ...
            ['''%s'' is a summary/spreadsheet export. It contains no posterior draws, so no ', ...
            'reliability can be reprocessed from it. How to fix: reprocess the saved .psyrat ', ...
            'result the export came from, or refit.'], src.path);
    end
    if ~ismember(lower(ext), {'.psyrat', '.mat'})
        error([ID ':bad_source_file'], 'Source file must be .psyrat or .mat; got ''%s''.', src.path);
    end
    if ~isfile(src.path)
        error([ID ':source_not_found'], 'Source file ''%s'' does not exist.', src.path);
    end
    if isempty(outputPath)
        error([ID ':output_required'], ...
            ['File input requires an explicit ''output'' path. The source file is never ', ...
            'modified in place. How to fix: pass ''output'', ''<fresh>.psyrat''.']);
    end
    fileVars = load(src.path, '-mat');
    src.file_vars = fieldnames(fileVars)';
    if isfield(fileVars, 'psyrat_data') && isstruct(fileVars.psyrat_data)
        wrapper = fileVars.psyrat_data;
        src.description = sprintf('file %s (variable psyrat_data)', src.path);
    else
        error([ID ':no_psyrat_data'], ...
            ['Source file ''%s'' holds no ''psyrat_data'' variable (variables: %s). Only files ', ...
            'written by PsyRAT (load(file,''-mat'') contract) are reprocessed.'], src.path, strjoin(src.file_vars, ', '));
    end
    if ~isfield(wrapper, 'rel') || ~isstruct(wrapper.rel)
        error([ID ':missing_field'], 'psyrat_data in ''%s'' has no ''rel'' struct; nothing to reprocess.', src.path);
    end
    rel = wrapper.rel;
elseif isstruct(source) && isfield(source, 'rel') && isstruct(source.rel)
    src.kind = 'psyrat_data';
    wrapper = source;
    rel = source.rel;
    src.description = 'psyrat_data struct';
elseif isstruct(source) && isfield(source, 'out') && isfield(source, 'analysis')
    src.kind = 'rel';
    wrapper = struct();
    rel = source;
    src.description = 'REL struct';
    if ~isempty(outputPath)
        error([ID ':output_not_allowed'], ...
            ['A bare REL struct has no psyrat_data wrapper to preserve, so it cannot be ', ...
            'written. How to fix: pass the full psyrat_data struct or the .psyrat file.']);
    end
else
    error([ID ':bad_source'], ...
        'source must be a psyrat_data struct, a REL struct (with .out and .analysis), or a .psyrat/.mat path.');
end

% ------------------------------------------------------------------ output rules
if ~isempty(outputPath)
    [~, ~, oext] = fileparts(outputPath);
    if ~ismember(lower(oext), {'.psyrat', '.mat'})
        error([ID ':bad_output'], 'output must end in .psyrat or .mat; got ''%s''.', outputPath);
    end
    if strcmp(src.kind, 'file') && strcmp(local_fullpath(outputPath), local_fullpath(src.path))
        error([ID ':output_is_source'], ...
            'output must differ from the source file (''%s''); in-place migration is refused.', src.path);
    end
    if exist(outputPath, 'file') || exist(outputPath, 'dir')
        error([ID ':output_exists'], ...
            ['Refusing to write: ''%s'' already exists. Existing files are never overwritten. ', ...
            'How to fix: choose a fresh destination.'], outputPath);
    end
end

% ------------------------------------------------------------------ route
if ~isfield(rel, 'family') || ~strcmp(local_str(rel, 'family'), 'gamma')
    error([ID ':unsupported_family'], ...
        ['Only family = ''gamma'' results are reprocessed (source family: %s). Gaussian results ', ...
        'are unaffected by the gamma residual contract.'], local_str(rel, 'family'));
end
if ~isfield(rel, 'out') || ~isstruct(rel.out)
    error([ID ':missing_field'], 'REL.out is missing; nothing to reprocess.');
end
analysis = local_str(rel, 'analysis');
switch analysis
    case 'trt_sserrvar',           route = 'static25';
    case 'ic_dynrel',              route = 'dyn11';
    case 'ic_dynrel_trt',          route = 'dyn14';
    case 'ic_dynrel_sserrvar',     route = 'dyn26';
    case 'ic_dynrel_sserrvar_trt', route = 'dyn27';
    otherwise
        error([ID ':unsupported_analysis'], ...
            ['Analysis ''%s'' is not covered by the gamma residual-target contract (covered: ', ...
            'trt_sserrvar, ic_dynrel, ic_dynrel_trt, ic_dynrel_sserrvar, ic_dynrel_sserrvar_trt). ', ...
            'Difference, copula, DoD, splits and the static one-facet designs are not reprocessed.'], analysis);
end

% ------------------------------------------------------------------ stamp version (exact)
% Only an ABSENT gamma_target field is the historical version-1 encoding. A
% field that is present but empty or malformed is not "unstamped": it is a
% partial/unknown provenance and is refused rather than reinterpreted (an
% empty marker on an already-migrated result would otherwise add R again).
hasStamp = isfield(rel.out, 'gamma_target');
if hasStamp
    gt0 = rel.out.gamma_target;
    if isempty(gt0) || ~isstruct(gt0) || ~isscalar(gt0)
        error([ID ':unsupported_target_version'], ...
            ['REL.out.gamma_target is present but empty or malformed (%s). Only an absent field ', ...
            'denotes the historical version-1 encoding; a present marker must be a version-%d stamp. ', ...
            'Refusing to interpret partial provenance.'], class(gt0), TARGET_VERSION);
    end
    if ~isfield(gt0, 'gamma_target_version')
        error([ID ':unsupported_target_version'], ...
            'REL.out.gamma_target exists but carries no gamma_target_version; refusing to interpret it.');
    end
    ver = gt0.gamma_target_version;
    if ~isnumeric(ver) || ~isscalar(ver) || ~isfinite(ver) || ver ~= TARGET_VERSION
        error([ID ':unsupported_target_version'], ...
            ['REL.out.gamma_target.gamma_target_version = %s is not supported. Only version %d ', ...
            'is recognized; an ABSENT stamp is the historical observation-only encoding. No other ', ...
            'version was ever produced, so none is interpreted.'], mat2str(ver), TARGET_VERSION);
    end
    versionFrom = TARGET_VERSION;
else
    % the gamma_resid raw-piece schema exists only under version 2: its presence
    % without a stamp is inconsistent provenance (the composite encoding must
    % not be reinterpreted as observation-only C and pooled again)
    if isfield(rel.out, 'gamma_resid')
        error([ID ':inconsistent_provenance'], ...
            ['REL.out.gamma_resid is present but REL.out.gamma_target is absent. The raw residual-piece ', ...
            'schema was introduced only with gamma_target_version %d, so this result''s encoding cannot ', ...
            'be identified as version 1; refusing rather than pooling R again.'], TARGET_VERSION);
    end
    versionFrom = 1;
end

% ------------------------------------------------------------------ parameterization (all provenance must agree)
parameterization = local_parameterization(rel, hasStamp, ID);

% ------------------------------------------------------------------ record skeleton
record = struct();
record.tool = ID;
record.tool_version = local_toolbox_version();
record.time = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss'));
record.source = src.description;
record.source_kind = src.kind;
record.output = outputPath;
record.analysis = analysis;
record.route = route;
record.parameterization = parameterization;
record.version_from = versionFrom;
record.version_to = TARGET_VERSION;
record.status = '';
record.method = '';
record.surface_target = '';
record.count_source = '';
record.fields_read = {};
record.fields_written = {};
record.fields_derived = {};
record.fields_absent = {};
record.transforms = struct([]);
record.id_mapping = {};
record.invalidation = {};
record.limitations = {};
record.notes = {};

% historical 'marginal' flags (0 or 1) are NOT a target; note and ignore
for holder = {rel, rel.out}
    if isfield(holder{1}, 'marginal')
        record.notes{end+1} = sprintf(['historical marginal=%s flag found and ignored: both 0 and 1 ', ...
            'denote the old population surface, not a target choice'], mat2str(holder{1}.marginal));
    end
end

% ------------------------------------------------------------------ dispatch
newrel = rel;
switch route
    case 'static25'
        [newrel, record] = local_static25(newrel, record, hasStamp, ID, TARGET_VERSION);
    otherwise
        [newrel, record] = local_dynamic(newrel, record, hasStamp, route, surfaceTargetIn, ID, TARGET_VERSION);
end

% ------------------------------------------------------------------ cached summaries
if strcmp(record.status, 'migrated') && strcmp(route, 'static25') && isfield(wrapper, 'relsummary')
    wrapper = rmfield(wrapper, 'relsummary');
    record.invalidation{end+1} = 'psyrat_data.relsummary removed (derived under version 1; recompute with psyrat_relsummary)';
end

% ------------------------------------------------------------------ store record
if ~strcmp(record.status, 'already_current')
    gt = newrel.out.gamma_target;
    if isfield(gt, 'migration_history') && isstruct(gt.migration_history) && ~isempty(gt.migration_history)
        hist = gt.migration_history;
        hist(end+1) = local_conform(record, hist(1)); %#ok<AGROW>
        gt.migration_history = hist;
    else
        gt.migration_history = record;
    end
    newrel.out.gamma_target = gt;
end

% ------------------------------------------------------------------ assemble / write
switch src.kind
    case 'rel'
        updated = newrel;
    otherwise
        wrapper.rel = newrel;
        updated = wrapper;
end
if ~isempty(outputPath)
    if strcmp(src.kind, 'file')
        fileVars.psyrat_data = updated;     % every other variable carried unchanged
        save(outputPath, '-struct', 'fileVars');
    else
        psyrat_data = updated; %#ok<NASGU>
        save(outputPath, 'psyrat_data');
    end
end
end

% =========================================================================
%                              STATIC ANALYSIS 25
% =========================================================================
function [rel, record] = local_static25(rel, record, hasStamp, ID, TARGET_VERSION)
out = rel.out;
req = {'mu', 'gro_sds', 'pop_sdlog', 'ind_sdlog', 'sig_occ', 'sig_trl', ...
    'sig_trlxid', 'sig_occxid', 'sig_trlxocc', 'id_matches', 'labels'};
local_require(out, req, ID, 'analysis 25');
record.fields_read = req;

[ndraws, nstrata] = local_mu_shape(out, ID);
local_check_labels(out, nstrata, ID);
for f = {'pop_sdlog', 'sig_occ', 'sig_trl', 'sig_trlxid', 'sig_occxid', 'sig_trlxocc'}
    local_check_matrix(out, f{1}, ndraws, nstrata, ID);
end
for f = {'gro_sds', 'ind_sdlog', 'id_matches'}
    local_check_cellcount(out, f{1}, nstrata, ID);
end
hasIndBs = isfield(out, 'ind_bs') && iscell(out.ind_bs) && numel(out.ind_bs) == nstrata;
if ~hasIndBs
    record.fields_absent{end+1} = 'ind_bs (person log-mean effects; not required for reprocessing)';
end

% ---- common per-stratum alignment / finiteness / ID validation (ALWAYS, before any early return) ----
S = cell(1, nstrata);
for s = 1:nstrata
    st = struct();
    st.mu  = out.mu(:, s);
    gro = local_cellmat(out.gro_sds{s}, sprintf('gro_sds{%d}', s), ID);
    if size(gro, 1) ~= ndraws || size(gro, 2) < 2
        error([ID ':misaligned_draws'], ...
            'gro_sds{%d} is %s; expected [%d x 2] (observed person SD, scale heterogeneity SD).', ...
            s, mat2str(size(gro)), ndraws);
    end
    local_finite(gro, sprintf('gro_sds{%d}', s), ID, 'nonneg');
    st.sdP  = gro(:, 1);
    st.sdO  = out.sig_occ(:, s);
    st.sdI  = out.sig_trl(:, s);
    st.sdPO = out.sig_occxid(:, s);     % person x occasion (extractor slot)
    st.sdPI = out.sig_trlxid(:, s);     % person x trial
    st.sdOI = out.sig_trlxocc(:, s);
    local_finite(st.mu, sprintf('mu(:,%d)', s), ID, 'pos');
    local_normal(st.mu, sprintf('mu(:,%d)', s), ID, false);
    for f = {'sdP', 'sdO', 'sdI', 'sdPO', 'sdPI', 'sdOI'}
        local_finite(st.(f{1}), sprintf('%s stratum %d', f{1}, s), ID, 'nonneg');
        local_normal(st.(f{1}), sprintf('%s stratum %d', f{1}, s), ID, true);
    end
    st.pop = out.pop_sdlog(:, s);
    local_finite(st.pop, sprintf('pop_sdlog(:,%d)', s), ID, 'any');
    st.ind = local_cellmat(out.ind_sdlog{s}, sprintf('ind_sdlog{%d}', s), ID);
    if size(st.ind, 1) ~= ndraws
        error([ID ':misaligned_draws'], 'ind_sdlog{%d} has %d draw rows; REL.out.mu has %d.', s, size(st.ind, 1), ndraws);
    end
    local_finite(st.ind, sprintf('ind_sdlog{%d}', s), ID, 'any');
    st.nsub = size(st.ind, 2);
    idt = out.id_matches{s};
    if ~istable(idt) || ~all(ismember({'id', 'id2'}, idt.Properties.VariableNames))
        error([ID ':invalid_id_matches'], 'id_matches{%d} must be a table with columns id and id2.', s);
    end
    if height(idt) ~= st.nsub
        error([ID ':misaligned_subjects'], ...
            'ind_sdlog{%d} has %d subject columns but id_matches{%d} has %d rows.', s, st.nsub, s, height(idt));
    end
    local_check_ids(idt.id, st.nsub, sprintf('id_matches{%d}.id', s), ID);
    local_check_id2(idt.id2, st.nsub, sprintf('id_matches{%d}.id2', s), ID);
    if hasIndBs
        ib = local_cellmat(out.ind_bs{s}, sprintf('ind_bs{%d}', s), ID);
        if ~isequal(size(ib), [ndraws st.nsub])
            error([ID ':misaligned_subjects'], 'ind_bs{%d} is %s; expected [%d x %d].', s, mat2str(size(ib)), ndraws, st.nsub);
        end
        local_finite(ib, sprintf('ind_bs{%d}', s), ID, 'any');
    end
    st.idt = idt;
    S{s} = st;
end

% ---- already version 2: validate every identifier and every stratum, never pool again ----
if hasStamp
    local_validate_stamp_identifiers(out.gamma_target, local_static25_stamp(record.parameterization, 1), ...
        {'gamma_target_version', 'surface_target', 'reference_target', 'subject_target', ...
        'remainder_scheme', 'parameterization', 'count_source', 'surface_count_source', 'analysis'}, ID);
    local_validate_static25_v2(out, S, ndraws, nstrata, record.parameterization, ID);
    record.status = 'already_current';
    record.method = 'validated existing version-2 encoding; no numeric change';
    record.count_source = out.gamma_target.count_source;
    record.surface_target = 'not_applicable';
    record.notes{end+1} = sprintf(['source already carries gamma_target_version %d with valid, aligned ', ...
        'gamma_resid per stratum; numerics returned unchanged'], TARGET_VERSION);
    return;
end

% ---- optional last-stratum raw pieces for the cross-check ----
gd = struct();
if isfield(out, 'gamma_disp') && isstruct(out.gamma_disp), gd = out.gamma_disp; end
canCross = all(isfield(gd, {'alpha', 's_o', 's_t', 's_po', 's_pt', 's_ot'}));

transforms = struct([]);
newpop = out.pop_sdlog;
newind = out.ind_sdlog;
gresid = cell(1, nstrata);
idmap = cell(1, nstrata);
for s = 1:nstrata
    st = S{s};
    inv = local_static25_inverse(st, s, ID);

    crossChecked = false;
    if canCross && s == nstrata
        local_crosscheck_gamma_disp(gd, st, inv, s, ndraws, ID);
        crossChecked = true;
    end

    R = psyrat_gamma_mean_remainder_trt(inv.alpha, sqrt(inv.v(:,1)), sqrt(inv.v(:,2)), sqrt(inv.v(:,3)), ...
        sqrt(inv.v(:,4)), sqrt(inv.v(:,5)), sqrt(inv.v(:,6)));
    if any(~isfinite(R) | R < 0)
        error([ID ':unsupported_numeric_regime'], ...
            'Stratum %d: the reconstructed remainder R is not finite/nonnegative (overflow or underflow); not reprocessed.', s);
    end
    % R must be strictly positive normal wherever a structurally nonzero term
    % exists, and exactly zero only where none does (complete zero facet)
    local_positive_transform(R, inv.hasLive, sprintf('remainder R stratum %d', s), ID);
    Ctyp = exp(2 .* st.pop);                        % version-1 encoding: C_typ only
    Cp   = exp(2 .* (st.pop + st.ind));             % version-1 encoding: C_p only
    local_normal(Ctyp, sprintf('decoded C_typ stratum %d', s), ID, false);
    local_normal(Cp, sprintf('decoded C_p stratum %d', s), ID, false);
    Etyp = R + Ctyp;
    Ep   = R + Cp;
    local_normal(Etyp, sprintf('E_typ stratum %d', s), ID, false);
    local_normal(Ep, sprintf('E_p stratum %d', s), ID, false);
    pop2 = 0.5 .* log(Etyp);
    ind2 = 0.5 .* log(Ep) - pop2;
    local_finite(pop2, sprintf('pop_sdlog_v2 stratum %d', s), ID, 'any');
    local_finite(ind2, sprintf('ind_sdlog_v2 stratum %d', s), ID, 'any');
    newpop(:, s) = pop2;
    newind{s} = num2cell(ind2);                     % same storage convention, same column order

    gr = struct('mean_remainder', R(:), 'obs_var_typ', Ctyp(:), 'obs_var_ss', Cp, ...
        'scale0', [], 'ind_scale', [], 'parameterization', record.parameterization);
    if strcmp(record.parameterization, 'ls')
        gr.scale0 = st.pop(:);                      % version-1 LS stored Intercept_sigma / ind_sd directly
        gr.ind_scale = st.ind;
        scaleNote = 'scale0/ind_scale exact from version-1 pop_sdlog/ind_sdlog (Intercept_sigma, ind_sd)';
    else
        scaleNote = 'scale0/ind_scale (Intercept_nu, ind_nu) absent: version 1 retained only C_typ/C_p; not manufactured';
    end
    gresid{s} = gr;

    t = struct('stratum', s, 'label', local_label(out.labels, s), 'ndraws', ndraws, 'nsub', st.nsub, ...
        'compatibility', inv.compat, 'crosschecked_gamma_disp', crossChecked, ...
        'reconstruction_exact_within_model', inv.compat.n_unresolved == 0, ...
        'R_median', median(R), 'R_min', min(R), 'R_max', max(R), ...
        'Etyp_median_v1', median(Ctyp), 'Etyp_median_v2', median(Etyp), ...
        'scale_provenance', scaleNote, ...
        'transform', 'pop_v2 = .5 log(R + exp(2 pop_v1)); ind_v2 = .5 log(R + exp(2 (pop_v1 + ind_v1))) - pop_v2');
    if isempty(transforms), transforms = t; else, transforms(end+1) = t; end %#ok<AGROW>
    idmap{s} = st.idt;
end

out.pop_sdlog = newpop;
out.ind_sdlog = newind;
out.gamma_resid = gresid;
if isfield(gd, 'sigma_pot_e2_pop')
    out.gamma_disp.sigma_pot_e2_pop_v1_observation_only = gd.sigma_pot_e2_pop;
    out.gamma_disp.sigma_pot_e2_pop = gresid{nstrata}.mean_remainder + gresid{nstrata}.obs_var_typ;
    record.fields_written{end+1} = 'gamma_disp.sigma_pot_e2_pop (now E_typ of the last stratum; v1 value kept as sigma_pot_e2_pop_v1_observation_only)';
end
out.gamma_target = local_static25_stamp(record.parameterization, 1);
rel.out = out;

record.status = 'migrated';
record.method = ['normalized per-draw inverse (a,b,c,h from (sd/mu)^2), pair compatibility against the ', ...
    'declared operation-norm bound (frozen_inverse_compatibility.json), R via psyrat_gamma_mean_remainder_trt, ', ...
    'one-time variance-domain addition, log-SD re-encoding'];
record.count_source = 'mean_trials_per_occupied_occasion';
record.surface_target = 'not_applicable';
record.fields_written = [record.fields_written, {'pop_sdlog', 'ind_sdlog', 'gamma_resid', 'gamma_target'}];
record.fields_derived = {'gamma_resid{s}.mean_remainder (R from mu and the six observed-scale SDs)', ...
    'gamma_resid{s}.obs_var_typ = exp(2 pop_v1)', 'gamma_resid{s}.obs_var_ss = exp(2 (pop_v1 + ind_v1))'};
if strcmp(record.parameterization, 'lognu')
    record.fields_absent{end+1} = 'gamma_resid{s}.scale0 / ind_scale (Intercept_nu, ind_nu not retained by version 1)';
end
record.transforms = transforms;
record.id_mapping = idmap;
record.limitations = {['R reconstructed from posterior draws of observed-scale components under a declared ', ...
    'finite-normal one-ulp error model (not a certified MATLAB/libm bound); draws whose pair deficit falls ', ...
    'inside the bound are recorded as roundoff-unresolved, not exact'], ...
    'person effects u_p/v_p are not reconstructed (unidentifiable from C_p alone); consumers must not require them', ...
    'the six observed components and their SD encodings are unchanged; only the residual encoding changes', ...
    'subject counts remain the source''s mean trials per occupied occasion; no surface override enters the subject table'};
end

% -------------------------------------------------------------------------
function gt = local_static25_stamp(param, fromVersion)
gt = struct( ...
    'gamma_target_version', 2, ...
    'surface_target', 'not_applicable', ...
    'reference_target', 'typical_person_zero_person_effects', ...
    'subject_target', 'own_observation_variance_with_population_mean_components_at_own_z', ...
    'remainder_scheme', 'population_highest_order_plus_person_observation_v2', ...
    'parameterization', param, ...
    'count_source', 'mean_trials_per_occupied_occasion', ...
    'surface_count_source', 'not_applicable', ...
    'analysis', 'trt_sserrvar', ...
    'migrated_from_version', fromVersion);
end

% -------------------------------------------------------------------------
function inv = local_static25_inverse(st, s, ID)
% Normalized per-draw inverse with component compatibility classified BEFORE
% any negative pair inverse is canonicalized (see header, "COMPONENT
% COMPATIBILITY"). Returns v [ndraws x 6] (log variances, order P O I PO PI OI),
% alpha, and the compatibility record.
u = 2^-53;
gam = @(n) n * u / (1 - n * u);
mu = st.mu;
% Declared numeric regime (finite NORMAL binary64, or exact zero for a
% component that is structurally absent): every stored input, its square,
% the normalized ratios and the prefactor B = mu^2 must lie in
% [realmin, realmax] in magnitude. Subnormal or overflowed intermediates
% violate the error model's assumptions and are refused, never carried.
local_normal(mu, sprintf('mu stratum %d', s), ID, false);
local_normal(mu.^2, sprintf('mu^2 (B) stratum %d', s), ID, false);
sds = {st.sdP, st.sdO, st.sdI, st.sdPO, st.sdPI, st.sdOI};
sdn = {'sd_P', 'sd_O', 'sd_I', 'sd_PO', 'sd_PI', 'sd_OI'};
% A zero is admissible ONLY where the source itself is zero (structurally
% absent component). A positive source that maps to zero, a subnormal or Inf
% under a transform (square, normalization, pair quotient, log1p) is a
% positive-to-zero/overflow underflow and is refused; it must never be carried
% forward as a structural zero (the remainder preflight would then wrongly
% classify products involving it as structurally absent).
for i = 1:6
    local_normal(sds{i}, sprintf('%s stratum %d', sdn{i}, s), ID, true);
    srcPos = sds{i} > 0;
    local_positive_transform(sds{i}.^2, srcPos, sprintf('%s^2 stratum %d', sdn{i}, s), ID);
    local_positive_transform((sds{i} ./ mu).^2, srcPos, sprintf('(%s/mu)^2 stratum %d', sdn{i}, s), ID);
end
a = (st.sdP ./ mu).^2;   b = (st.sdO ./ mu).^2;   c = (st.sdI ./ mu).^2;
h = [(st.sdPO ./ mu).^2, (st.sdPI ./ mu).^2, (st.sdOI ./ mu).^2];   % PO PI OI
qPositive = [st.sdP > 0, st.sdO > 0, st.sdI > 0, false(size(a, 1), 3)];   % filled per pair below
vP = log1p(a); vO = log1p(b); vI = log1p(c);
% unclamped pair inverse, used only for the condition factor of the error model
vraw = [vP vO vI, log1p(a + b + h(:,1)) - vP - vO, log1p(a + c + h(:,2)) - vP - vI, log1p(b + c + h(:,3)) - vO - vI];
if any(~isfinite(vraw(:)))
    error([ID ':unsupported_numeric_regime'], 'Stratum %d: unclamped inverse log variances are not finite; not reprocessed.', s);
end
condition = 1 + 2 .* abs(log(mu)) + sum(abs(vraw), 2);
eta = 2 * gam(32) .* condition;
if any(eta >= 1)
    error([ID ':unsupported_numeric_regime'], 'Stratum %d: error-model factor eta >= 1 (extreme magnitudes); not reprocessed.', s);
end
pairs = {'PO', 'PI', 'OI'};
left  = {a, a, b};
right = {b, c, c};
q = [a b c zeros(size(a, 1), 3)];
compat = struct('n_compatible_positive', 0, 'n_unresolved', 0, 'n_unresolved_negative_canonicalized', 0, ...
    'max_unresolved_negative_margin', 0, 'bound_min', Inf, 'bound_max', 0, ...
    'unresolved', struct('pair', {}, 'draw', {}, 'deficit', {}, 'bound', {}), ...
    'model', 'eta*S+|b|eta(1+a)+|a|eta(1+b)+eta^2(1+a)(1+b)+gamma3(|h|+|ab|); eta=2 gamma32 (1+2|log mu|+sum|v_raw|); u=2^-53');
for p = 1:3
    A = left{p}; Bb = right{p}; H = h(:, p);
    da = eta .* (1 + A); db = eta .* (1 + Bb);
    Ssum = (1 + A + Bb + H) + (1 + A) + (1 + Bb) + 1;
    bound = eta .* Ssum + abs(Bb) .* da + abs(A) .* db + da .* db + gam(3) .* (abs(H) + abs(A .* Bb));
    deficit = H - A .* Bb;
    bad = deficit < -bound;
    if any(bad)
        d = find(bad, 1);
        error([ID ':invalid_inverse'], ...
            ['Stratum %d, draw %d, component %s: normalized deficit h - a*b = %.6g is below -bound = %.6g ', ...
            '(materially incompatible with the gamma log-link model, e.g. a missing or corrupted ', ...
            'interaction component). Not reprocessed.'], s, d, pairs{p}, deficit(d), -bound(d));
    end
    ok = deficit > bound;
    unres = ~ok;
    compat.n_compatible_positive = compat.n_compatible_positive + nnz(ok);
    compat.n_unresolved = compat.n_unresolved + nnz(unres);
    negU = unres & deficit < 0;
    compat.n_unresolved_negative_canonicalized = compat.n_unresolved_negative_canonicalized + nnz(negU);
    if any(negU)
        compat.max_unresolved_negative_margin = max(compat.max_unresolved_negative_margin, max(-deficit(negU) ./ bound(negU)));
    end
    compat.bound_min = min(compat.bound_min, min(bound));
    compat.bound_max = max(compat.bound_max, max(bound));
    for d = find(unres)'
        compat.unresolved(end+1) = struct('pair', pairs{p}, 'draw', d, 'deficit', deficit(d), 'bound', bound(d)); %#ok<AGROW>
    end
    qd = deficit;
    qd(negU) = 0;                                   % canonicalize negatives INSIDE the bound only
    qp = qd ./ ((1 + A) .* (1 + Bb));
    % a strictly positive deficit must stay a strictly positive NORMAL pair q;
    % only a zero deficit (source zero or canonicalized negative) may give zero
    local_positive_transform(qp, qd > 0, sprintf('pair quotient q_%s stratum %d', pairs{p}, s), ID);
    q(:, 3 + p) = qp;
    qPositive(:, 3 + p) = qd > 0;
end
v = log1p(q);
alpha = log(mu) - sum(v, 2) ./ 2;
if any(~isfinite(v(:))) || any(~isfinite(alpha))
    error([ID ':unsupported_numeric_regime'], 'Stratum %d: reconstructed log variances/alpha are not finite; not reprocessed.', s);
end
local_positive_transform(q, qPositive, sprintf('normalized q stratum %d', s), ID);
local_positive_transform(v, qPositive, sprintf('reconstructed log variances stratum %d', s), ID);
hasLive = local_range_check_remainder(q, mu.^2, s, ID);
inv = struct('v', v, 'q', q, 'alpha', alpha, 'eta', eta, 'compat', compat, 'hasLive', hasLive);
end

% -------------------------------------------------------------------------
function local_positive_transform(x, srcPositive, what, ID)
% Transform-output guard with a source-positive mask: where the source was
% positive, the transformed value must be finite, strictly positive and NORMAL
% (>= realmin, <= realmax); where the source was zero, it must be exactly 0.
% A positive source underflowing to zero/subnormal or overflowing is refused.
if ~isnumeric(x) || ~isreal(x) || ~isequal(size(x), size(srcPositive))
    error([ID ':unsupported_numeric_regime'], '%s: transform output is not real numeric of the source shape.', what);
end
xp = x(srcPositive);
if any(~isfinite(xp))
    error([ID ':unsupported_numeric_regime'], '%s: a positive source overflowed or is not finite under the transform.', what);
end
if any(xp <= 0) || any(xp < realmin)
    error([ID ':unsupported_numeric_regime'], ...
        ['%s: a strictly positive source maps to zero or a subnormal (< realmin = %.3g) under the ', ...
        'transform. Positive-to-zero underflow is outside the declared finite-normal regime and ', ...
        'cannot be carried as a structural zero; not reprocessed.'], what, realmin);
end
if any(x(~srcPositive) ~= 0)
    error([ID ':unsupported_numeric_regime'], '%s: a zero source produced a nonzero transform.', what);
end
end

% -------------------------------------------------------------------------
function hasLive = local_range_check_remainder(q, B, s, ID)
% The remainder R = B * sum over the 63 factor subsets whose support union is
% {P,O,I} of prod q_A is evaluated by psyrat_gamma_mean_remainder_trt as
% plain products. A product of tiny NORMAL q values can underflow to zero (or
% a subnormal) before the large prefactor B restores it, silently dropping a
% positive term even though the true R is a normal number. That is outside
% the declared finite-normal regime, so it is refused here, per draw, by
% evaluating every structurally nonzero subset product in log space and
% requiring it and B times it to lie within [realmin, realmax]. Subsets with a
% zero factor are structurally zero (an absent facet) and stay exempt, so a
% complete zero-facet R = 0 remains valid.
% Returns hasLive [ndraws x 1]: true where at least one structurally nonzero
% term exists (so R must come out strictly positive), false where R = 0 is the
% legitimate structural value. Callers have already refused any positive
% source that underflowed, so a zero q here is a true structural zero.
lo = log(realmin); hi = log(realmax);
ndraw = size(q, 1);
hasLive = false(ndraw, 1);
supp = logical([1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1]);
logB = log(B);
for m = 1:63
    inT = logical(bitget(m, 1:6));
    if ~all(any(supp(inT, :), 1)), continue; end          % union must be {P,O,I}
    members = find(inT);
    qm = q(:, members);
    live = all(qm > 0, 2);                               % structurally nonzero products only
    hasLive = hasLive | live;
    if ~any(live), continue; end
    % prefix products in the helper's left-to-right column order, so no
    % intermediate partial product may leave the range either
    cs = cumsum(log(qm(live, :)), 2);
    badRow = any(cs < lo | cs > hi, 2);
    if any(badRow)
        d = find(live); d = d(find(badRow, 1));
        error([ID ':unsupported_numeric_regime'], ...
            ['Stratum %d, draw %d: a positive term of the remainder (product of q over factors [%s]) has ', ...
            'log magnitude %.4g, outside the finite-normal range [%.4g, %.4g]; the product evaluation ', ...
            'would underflow/overflow it. Refusing rather than reporting a truncated R.'], ...
            s, d, num2str(members), cs(find(badRow, 1), end), lo, hi);
    end
    ls = cs(:, end);
    lsB = ls + logB(live);
    if any(lsB < lo) || any(lsB > hi - log(63))
        d = find(live); d = d(find(lsB < lo | lsB > hi - log(63), 1));
        error([ID ':unsupported_numeric_regime'], ...
            ['Stratum %d, draw %d: B times a remainder term (factors [%s]) leaves the finite-normal range; ', ...
            'refusing rather than reporting a truncated R.'], s, d, num2str(members));
    end
end
end

% -------------------------------------------------------------------------
function local_normal(x, what, ID, allowZero)
% finite NORMAL binary64 (|x| in [realmin, realmax]); zero allowed only when
% the quantity may be structurally absent
if ~isnumeric(x) || ~isreal(x) || any(~isfinite(x(:)))
    error([ID ':unsupported_numeric_regime'], '%s is not finite real numeric.', what);
end
ax = abs(x(:));
sub = ax > 0 & ax < realmin;
if any(sub)
    error([ID ':unsupported_numeric_regime'], ...
        '%s contains subnormal values (|x| < realmin = %.3g); outside the declared finite-normal regime.', what, realmin);
end
if ~allowZero && any(ax == 0)
    error([ID ':unsupported_numeric_regime'], '%s must be nonzero (finite normal).', what);
end
end

% -------------------------------------------------------------------------
function local_crosscheck_gamma_disp(gd, st, inv, s, ndraws, ID)
% Compare the normalized components implied by the LAST stratum's raw log-scale
% pieces with the stored ones, elementwise, under twice the one-sided source
% bound (a discrepancy larger than both sides' allowances is a contradiction).
raw = {gd.alpha(:), gd.s_o(:), gd.s_t(:), gd.s_po(:), gd.s_pt(:), gd.s_ot(:)};
if any(cellfun(@numel, raw) ~= ndraws)
    error([ID ':inverse_crosscheck_failed'], ...
        'Stratum %d: gamma_disp raw pieces are not draw-aligned with REL.out.mu (%d draws).', s, ndraws);
end
if any(~isfinite([raw{:}]), 'all')
    error([ID ':inverse_crosscheck_failed'], 'Stratum %d: gamma_disp raw pieces are not finite.', s);
end
u = 2^-53; gam = @(n) n * u / (1 - n * u);
a = inv.q(:, 1);
bR = expm1(raw{2}.^2); cR = expm1(raw{3}.^2);
qR = [expm1(raw{4}.^2), expm1(raw{5}.^2), expm1(raw{6}.^2)];
hStored = [(st.sdPO ./ st.mu).^2, (st.sdPI ./ st.mu).^2, (st.sdOI ./ st.mu).^2];
hRaw = [a .* bR + qR(:,1) .* (1 + a) .* (1 + bR), a .* cR + qR(:,2) .* (1 + a) .* (1 + cR), ...
    bR .* cR + qR(:,3) .* (1 + bR) .* (1 + cR)];
left = {a, a, bR}; right = {bR, cR, cR};
names = {'PO', 'PI', 'OI'};
for p = 1:3
    A = left{p}; Bb = right{p}; H = hStored(:, p);
    da = inv.eta .* (1 + A); db = inv.eta .* (1 + Bb);
    Ssum = (1 + A + Bb + H) + (1 + A) + (1 + Bb) + 1;
    bound = inv.eta .* Ssum + abs(Bb) .* da + abs(A) .* db + da .* db + gam(3) .* (abs(H) + abs(A .* Bb));
    dev = abs(H - hRaw(:, p));
    if any(dev > 2 .* bound)
        d = find(dev > 2 .* bound, 1);
        error([ID ':inverse_crosscheck_failed'], ...
            ['Stratum %d, draw %d, component %s: the normalized component implied by gamma_disp raw pieces ', ...
            '(%.6g) differs from the stored one (%.6g) by %.3g, beyond the two-sided source bound %.3g. ', ...
            'Stored components and raw provenance contradict; not reprocessed.'], s, d, names{p}, hRaw(d, p), H(d), dev(d), 2 * bound(d));
    end
end
% main effects and the mean: |b - b_raw| <= 2 eta (1+b); |log mu - (alpha_raw + V_raw/2)| <= 2 eta (1 + |log mu|)
for p = 1:2
    stored = inv.q(:, 1 + p); implied = {bR, cR};
    dev = abs(stored - implied{p});
    if any(dev > 2 .* inv.eta .* (1 + stored))
        d = find(dev > 2 .* inv.eta .* (1 + stored), 1);
        error([ID ':inverse_crosscheck_failed'], ...
            'Stratum %d, draw %d: main-effect component %s implied by gamma_disp differs from the stored one by %.3g.', ...
            s, d, names{p}(2), dev(d));
    end
end
Vraw = inv.v(:, 1) + raw{2}.^2 + raw{3}.^2 + raw{4}.^2 + raw{5}.^2 + raw{6}.^2;
devMu = abs(log(st.mu) - (raw{1} + Vraw ./ 2));
if any(devMu > 2 .* inv.eta .* (1 + abs(log(st.mu))))
    d = find(devMu > 2 .* inv.eta .* (1 + abs(log(st.mu))), 1);
    error([ID ':inverse_crosscheck_failed'], ...
        'Stratum %d, draw %d: log mu differs from gamma_disp alpha + V/2 by %.3g.', s, d, devMu(d));
end
end

% -------------------------------------------------------------------------
function local_validate_static25_v2(out, S, ndraws, nstrata, param, ID)
% A version-2 claim must be backed by aligned, finite, consistent gamma_resid.
if ~isfield(out, 'gamma_resid') || ~iscell(out.gamma_resid) || numel(out.gamma_resid) ~= nstrata
    error([ID ':invalid_v2_claim'], ...
        ['gamma_target claims version 2 but gamma_resid is missing or has %d strata (expected %d). ', ...
        'The claim is not backed by raw residual pieces; refusing.'], local_count(out, 'gamma_resid'), nstrata);
end
need = {'mean_remainder', 'obs_var_typ', 'obs_var_ss', 'scale0', 'ind_scale', 'parameterization'};
% Log-domain encoding identity bound (operation count, no absolute floor).
% Producer: E = R + C (1 op), log E (elementary, <= 1 ulp = 2u relative),
% .5*log E (exact), ind = .5 log E_p - pop (1 op, cancels against pop).
% Validator: E = R + C (1 op), log E (elementary), pop + ind (1 op, cancels),
% 2*(.) (exact). Each operation contributes at most u times the magnitude of
% its operands, so the accumulated absolute error is bounded by
%   gamma_8 * (1 + |log E| + 2|pop| + 2|ind|)
% with gamma_8 = 8u/(1-8u): the 1 covers the relative error of the sums seen
% through log, |log E| the elementary-function ulps on both sides, and the
% pop/ind magnitudes the two cancelling sums (a large |pop| with an almost
% equal |ind| loses digits in both the producer's subtraction and the
% validator's re-summation; INVERSE_ROUNDOFF_REVIEW.md requires those
% magnitudes in the bound). The 100% mismatch E = 1e-12 vs R + C = 2e-12
% (|log 2| = .69) stays refused because no magnitude term there exceeds ~30.
u = 2^-53; gam8 = 8 * u / (1 - 8 * u);
for s = 1:nstrata
    g = out.gamma_resid{s};
    st = S{s};
    if ~isstruct(g) || ~all(isfield(g, need))
        error([ID ':invalid_v2_claim'], 'gamma_resid{%d} lacks the raw-piece schema fields (%s).', s, strjoin(need, ', '));
    end
    if ~strcmp(local_str(g, 'parameterization'), param)
        error([ID ':contradictory_parameterization'], ...
            'gamma_resid{%d}.parameterization = ''%s'' contradicts the resolved parameterization ''%s''.', ...
            s, local_str(g, 'parameterization'), param);
    end
    R = g.mean_remainder(:); Ct = g.obs_var_typ(:); Cs = g.obs_var_ss;
    if numel(R) ~= ndraws || numel(Ct) ~= ndraws || ~isequal(size(Cs), [ndraws st.nsub])
        error([ID ':invalid_v2_claim'], ...
            'gamma_resid{%d} draw/subject dimensions do not match REL.out.mu / ind_sdlog{%d}.', s, s);
    end
    local_finite(R, sprintf('gamma_resid{%d}.mean_remainder', s), ID, 'nonneg');
    local_finite(Ct, sprintf('gamma_resid{%d}.obs_var_typ', s), ID, 'pos');
    local_finite(Cs, sprintf('gamma_resid{%d}.obs_var_ss', s), ID, 'pos');
    local_normal(R, sprintf('gamma_resid{%d}.mean_remainder', s), ID, true);
    local_normal(Ct, sprintf('gamma_resid{%d}.obs_var_typ', s), ID, false);
    local_normal(Cs, sprintf('gamma_resid{%d}.obs_var_ss', s), ID, false);
    Et = R + Ct; Ep = R + Cs;
    local_normal(Et, sprintf('gamma_resid{%d} E_typ', s), ID, false);
    local_normal(Ep, sprintf('gamma_resid{%d} E_p', s), ID, false);
    decT = exp(2 .* st.pop); decP = exp(2 .* (st.pop + st.ind));
    local_normal(decT, sprintf('decoded pop_sdlog(:,%d)', s), ID, false);
    local_normal(decP, sprintf('decoded pop+ind stratum %d', s), ID, false);
    % log-domain identity, elementwise, cancellation-aware (see above)
    lt = log(Et); lp = log(Ep);
    devT = abs(lt - 2 .* st.pop); devP = abs(lp - 2 .* (st.pop + st.ind));
    bT = gam8 .* (1 + abs(lt) + 2 .* abs(st.pop));
    bP = gam8 .* (1 + abs(lp) + 2 .* abs(st.pop) + 2 .* abs(st.ind));
    if any(devT > bT)
        d = find(devT > bT, 1);
        error([ID ':invalid_v2_claim'], ...
            ['gamma_resid{%d}: stored pop_sdlog does not encode R + C_typ at draw %d ', ...
            '(|log(R+C) - 2 pop| = %.3g > %.3g). The version-2 claim is inconsistent; refusing.'], s, d, devT(d), bT(d));
    end
    if any(devP(:) > bP(:))
        [d, c] = find(devP > bP, 1);
        error([ID ':invalid_v2_claim'], ...
            ['gamma_resid{%d}: stored ind_sdlog does not encode R + C_p at draw %d, column %d ', ...
            '(|log(R+C) - 2 (pop+ind)| = %.3g > %.3g). The version-2 claim is inconsistent; refusing.'], ...
            s, d, c, devP(d, c), bP(d, c));
    end
    if ~isempty(g.scale0) && (numel(g.scale0) ~= ndraws || any(~isfinite(g.scale0(:))))
        error([ID ':invalid_v2_claim'], 'gamma_resid{%d}.scale0 is present but not a finite [%d x 1] block.', s, ndraws);
    end
    if ~isempty(g.ind_scale) && (~isequal(size(g.ind_scale), [ndraws st.nsub]) || any(~isfinite(g.ind_scale(:))))
        error([ID ':invalid_v2_claim'], 'gamma_resid{%d}.ind_scale is present but not a finite [%d x %d] block.', s, ndraws, st.nsub);
    end
end
end

% =========================================================================
%                          DYNAMIC 11 / 14 / 26 / 27
% =========================================================================
function [rel, record] = local_dynamic(rel, record, hasStamp, route, surfaceTargetIn, ID, TARGET_VERSION)
out = rel.out;
param = record.parameterization;
twoFacet = ismember(route, {'dyn14', 'dyn27'});
subject  = ismember(route, {'dyn26', 'dyn27'});

req = {'mu', 'gro_sds', 'b', 'dimz', 'labels', 'ntrials', 'chol_corrmat'};
if twoFacet
    req = [req, {'sig_occ', 'sig_trl', 'sig_trlxid', 'sig_occxid', 'sig_trlxocc', 'nocc'}];
else
    req = [req, {'sig_trl'}];
end
if strcmp(param, 'lognu'), req = [req, {'pop_lognu', 'b_nu'}]; else, req = [req, {'pop_sdlog', 'b_sigma'}]; end
if subject
    req = [req, {'ssinfo'}];
    if strcmp(param, 'lognu'), req = [req, {'ind_bs', 'ind_nu'}]; else, req = [req, {'ind_sdlog'}]; end
end
local_require(out, req, ID, sprintf('analysis %s', record.analysis));
record.fields_read = req;

% dimensions
if ~isfield(rel, 'ndim') || ~isnumeric(rel.ndim) || ~isscalar(rel.ndim) || ~ismember(rel.ndim, [1 2])
    error([ID ':invalid_dimensions'], 'REL.ndim must be exactly 1 or 2; got %s.', local_str(rel, 'ndim'));
end
ndim = rel.ndim;
if ~isfield(rel, 'dim_names') || ~iscell(rel.dim_names) || numel(rel.dim_names) ~= ndim ...
        || ~all(cellfun(@(x) (ischar(x) || isstring(x)) && ~isempty(char(x)), rel.dim_names))
    error([ID ':invalid_dimensions'], ...
        'REL.dim_names must be a cell of %d nonempty name(s) matching REL.ndim; got %s.', ndim, local_str(rel, 'dim_names'));
end
if ndim == 1, nbasis = 1; else, nbasis = 3; end   % z1 | z1, z2, z1*z2

[ndraws, nstrata] = local_mu_shape(out, ID);
local_check_labels(out, nstrata, ID);
local_check_counts(out, 'ntrials', nstrata, ID);
if twoFacet, local_check_counts(out, 'nocc', nstrata, ID); end
popScale = 'pop_lognu'; slopeScale = 'b_nu'; ownScale = 'ind_nu';
if strcmp(param, 'ls'), popScale = 'pop_sdlog'; slopeScale = 'b_sigma'; ownScale = 'ind_sdlog'; end
local_check_matrix(out, popScale, ndraws, nstrata, ID);
local_finite(out.(popScale), popScale, ID, 'any');
facets = {'sig_trl'};
if twoFacet, facets = [facets, {'sig_occ', 'sig_trlxid', 'sig_occxid', 'sig_trlxocc'}]; end
for f = facets
    local_check_matrix(out, f{1}, ndraws, nstrata, ID);
    local_finite(out.(f{1}), f{1}, ID, 'nonneg');
end
cellFields = [{'gro_sds', 'b', slopeScale, 'dimz', 'chol_corrmat'}, req(ismember(req, {'ssinfo', 'ind_bs', 'ind_nu', 'ind_sdlog'}))];
for f = cellFields
    local_check_cellcount(out, f{1}, nstrata, ID);
end

% resolve the surface target: explicit argument > existing explicit stamp > default
if ~isempty(surfaceTargetIn)
    surfaceTarget = surfaceTargetIn;
elseif hasStamp && isfield(out.gamma_target, 'surface_target') ...
        && ismember(local_str(out.gamma_target, 'surface_target'), {'population_marginal', 'typical_person'})
    surfaceTarget = out.gamma_target.surface_target;
else
    surfaceTarget = 'population_marginal';
end

transforms = struct([]);
idmap = cell(1, nstrata);
for s = 1:nstrata
    local_finite(out.mu(:, s), sprintf('mu(:,%d)', s), ID, 'any');
    gro = local_cellmat(out.gro_sds{s}, sprintf('gro_sds{%d}', s), ID);
    if size(gro, 1) ~= ndraws || size(gro, 2) < 2
        error([ID ':misaligned_draws'], ...
            'gro_sds{%d} is %s; expected [%d x 2] (person log-mean SD, scale heterogeneity SD).', s, mat2str(size(gro)), ndraws);
    end
    local_finite(gro, sprintf('gro_sds{%d}', s), ID, 'nonneg');
    bm = local_cellmat(out.b{s}, sprintf('b{%d}', s), ID);
    bs = local_cellmat(out.(slopeScale){s}, sprintf('%s{%d}', slopeScale, s), ID);
    if size(bm, 1) ~= ndraws || size(bm, 2) ~= nbasis
        error([ID ':misaligned_draws'], 'b{%d} is %s; expected [%d x %d] for ndim = %d.', s, mat2str(size(bm)), ndraws, nbasis, ndim);
    end
    if size(bs, 1) ~= ndraws || size(bs, 2) ~= nbasis
        error([ID ':misaligned_draws'], '%s{%d} is %s; expected [%d x %d] for ndim = %d.', slopeScale, s, mat2str(size(bs)), ndraws, nbasis, ndim);
    end
    local_finite(bm, sprintf('b{%d}', s), ID, 'any');
    local_finite(bs, sprintf('%s{%d}', slopeScale, s), ID, 'any');
    dz = out.dimz{s};
    if ~isnumeric(dz) || ~isreal(dz) || isempty(dz) || size(dz, 2) ~= ndim || any(~isfinite(dz(:)))
        error([ID ':invalid_dimensions'], 'dimz{%d} must be a finite real numeric [n x %d] standardized-dimension matrix.', s, ndim);
    end
    % complete draw-aligned Cholesky correlation draws, required for BOTH
    % parameterizations and BOTH surface choices (the summary, variance table
    % and provenance consumers read them unconditionally; a reusable result
    % must be able to switch to the population-marginal primary)
    cholDev = local_check_chol(out.chol_corrmat{s}, ndraws, s, ID);
    nsub = NaN;
    if subject
        si = out.ssinfo{s};
        if ~istable(si) || ~all(ismember({'id', 'id2', 'z1', 'trls'}, si.Properties.VariableNames))
            error([ID ':invalid_ssinfo'], 'ssinfo{%d} must be a table with columns id, id2, z1, trls.', s);
        end
        if ndim == 2 && ~ismember('z2', si.Properties.VariableNames)
            error([ID ':invalid_ssinfo'], 'ssinfo{%d} lacks z2 for a two-dimension result.', s);
        end
        nsub = height(si);
        local_check_ids(si.id, nsub, sprintf('ssinfo{%d}.id', s), ID);
        local_check_id2(si.id2, nsub, sprintf('ssinfo{%d}.id2', s), ID);
        zcols = {'z1'}; if ndim == 2, zcols = {'z1', 'z2'}; end
        for zc = zcols
            local_check_column(si.(zc{1}), nsub, sprintf('ssinfo{%d}.%s', s, zc{1}), ID, 'any');
        end
        local_check_column(si.trls, nsub, sprintf('ssinfo{%d}.trls', s), ID, 'pos');
        ownFields = {ownScale};
        if strcmp(param, 'lognu'), ownFields = {'ind_bs', 'ind_nu'}; end
        for f = ownFields
            m = local_cellmat(out.(f{1}){s}, sprintf('%s{%d}', f{1}, s), ID);
            if ~isequal(size(m), [ndraws nsub])
                error([ID ':misaligned_subjects'], ...
                    '%s{%d} is %s; expected [%d x %d] (draws x subjects, columns indexed by id2).', ...
                    f{1}, s, mat2str(size(m)), ndraws, nsub);
            end
            local_finite(m, sprintf('%s{%d}', f{1}, s), ID, 'any');
        end
        idmap{s} = si(:, {'id', 'id2'});
    end
    t = struct('stratum', s, 'label', local_label(out.labels, s), 'ndraws', ndraws, 'nsub', nsub, ...
        'chol_max_deviation', cholDev, 'transform', 'none (raw draws unchanged; target stamped)');
    if isempty(transforms), transforms = t; else, transforms(end+1) = t; end %#ok<AGROW>
end

expected = local_dynamic_stamp(route, record.analysis, param, surfaceTarget, record.version_from);

if hasStamp
    gt = out.gamma_target;
    % every scientific identifier must agree with the resolved contract; only
    % surface_target may differ, and only when explicitly requested
    fixedFields = {'gamma_target_version', 'reference_target', 'subject_target', 'remainder_scheme', ...
        'parameterization', 'count_source', 'surface_count_source', 'analysis'};
    local_validate_stamp_identifiers(gt, expected, fixedFields, ID);
    stampTarget = local_str(gt, 'surface_target');
    if ~ismember(stampTarget, {'population_marginal', 'typical_person'})
        error([ID ':contradictory_stamp'], ...
            'gamma_target.surface_target = ''%s'' is not a recognized dynamic surface target.', stampTarget);
    end
    if strcmp(stampTarget, surfaceTarget)
        record.status = 'already_current';
        record.method = 'validated existing version-2 stamp and raw draws; no change';
        record.count_source = expected.count_source;
        record.surface_target = surfaceTarget;
        record.notes{end+1} = 'source already carries a valid version-2 stamp with this target; numerics and stamp returned unchanged';
        return;
    end
    if isempty(surfaceTargetIn)
        error([ID ':contradictory_stamp'], 'internal: surface target differs without an explicit request.');
    end
    record.status = 'target_restamped';
    record.method = 'explicit surface_target change on a validated version-2 result; raw draws unchanged; the dynamic surface can change at the next summary';
    record.notes{end+1} = sprintf('surface_target changed %s -> %s (no refit; raw draws unchanged)', stampTarget, surfaceTarget);
    if isfield(gt, 'migration_history'), expected.migration_history = gt.migration_history; end
else
    record.status = 'migrated';
    record.method = ['validated raw draws of a historical (unstamped) result; version-2 target stamped; ', ...
        'raw draws unchanged, no refit needed; the next summary follows the stamp, so analysis-27 ', ...
        'per-participant coefficients (and their comparison with any cutoff; no inclusion flag is stored) ', ...
        'can change, and a surface_target other than the historical population_marginal can change ', ...
        'the dynamic surface'];
end
out.gamma_target = expected;
rel.out = out;

record.surface_target = surfaceTarget;
record.count_source = expected.count_source;
record.fields_written = {'gamma_target'};
record.transforms = transforms;
record.id_mapping = idmap;
record.limitations = {'raw draws are stamped, not changed; the resolved target takes effect only through a summary that honors gamma_target', ...
    'an unstamped legacy analysis-27 result is interpreted as legacy_observation_only_v1 by the resolver until reprocessed'};
if strcmp(route, 'dyn26')
    record.notes{end+1} = 'analysis 26 keeps its previously approved S17/S20 remainder (population mean interaction pooled per person)';
end
end

% -------------------------------------------------------------------------
function gt = local_dynamic_stamp(route, analysis, param, surfaceTarget, fromVersion)
switch route
    case 'dyn11'
        scheme = 'population_mean_remainder_plus_selected_surface_observation_v1';
        countSource = 'not_applicable'; subjTarget = 'not_applicable';
    case 'dyn14'
        scheme = 'population_mean_remainder_plus_selected_surface_observation_v1';
        countSource = 'not_applicable'; subjTarget = 'not_applicable';
    case 'dyn26'
        scheme = 'population_mean_remainder_plus_person_observation_v1';
        countSource = 'median_trials_per_occupied_occasion';
        subjTarget = 'own_observation_variance_with_population_mean_components_at_own_z';
    case 'dyn27'
        scheme = 'population_highest_order_plus_person_observation_v2';
        countSource = 'median_trials_per_occupied_occasion';
        subjTarget = 'own_observation_variance_with_population_mean_components_at_own_z';
end
gt = struct( ...
    'gamma_target_version', 2, ...
    'surface_target', surfaceTarget, ...
    'reference_target', 'not_applicable', ...
    'subject_target', subjTarget, ...
    'remainder_scheme', scheme, ...
    'parameterization', param, ...
    'count_source', countSource, ...
    'surface_count_source', 'default_stratum_median', ...
    'analysis', analysis, ...
    'migrated_from_version', fromVersion);
end

% -------------------------------------------------------------------------
function local_validate_stamp_identifiers(gt, expected, fields, ID)
for f = fields
    if ~isfield(gt, f{1})
        error([ID ':contradictory_stamp'], ...
            'gamma_target lacks the required identifier ''%s'' (expected %s); the version-2 claim is incomplete.', ...
            f{1}, local_str(expected, f{1}));
    end
    got = gt.(f{1}); want = expected.(f{1});
    if isnumeric(want)
        same = isnumeric(got) && isscalar(got) && isequal(got, want);
    else
        same = (ischar(got) || isstring(got)) && strcmp(char(got), want);
    end
    if ~same
        error([ID ':contradictory_stamp'], ...
            'gamma_target.%s = %s contradicts the contract for this result (%s); refusing.', ...
            f{1}, local_str(gt, f{1}), local_str(expected, f{1}));
    end
end
end

% -------------------------------------------------------------------------
function dev = local_check_chol(x, ndraws, s, ID)
% Complete [ndraws x 2 x 2] Cholesky factor of a 2x2 correlation matrix per
% draw, as a num2cell block or numeric array. Rows must be unit-norm lower
% triangular with unit first diagonal; rho = +-1 (L22 = 0) is a valid limit.
% Tolerance 1e-5 reflects CmdStan's default six-significant-digit CSV output
% of cholesky_factor_corr draws; it is a structural check, not a numeric floor.
tol = 1e-5;
what = sprintf('chol_corrmat{%d}', s);
% storage contract (see local_cellmat): the summary reads this block through
% cell2mat, so only the estimator's num2cell [ndraws x 2 x 2] form is supported
if ~iscell(x)
    error([ID ':unsupported_storage'], ...
        ['%s is stored as %s; the toolbox storage contract is a num2cell [%d x 2 x 2] block of ', ...
        'scalar draws. Refusing: the summary reads it through cell2mat.'], what, class(x), ndraws);
end
if ~isequal(size(x), [ndraws 2 2]) || ~all(cellfun(@(v) isnumeric(v) && isreal(v) && isscalar(v), x(:)))
    error([ID ':invalid_correlation_draws'], ...
        '%s must be a complete [%d x 2 x 2] num2cell block of real scalar draws; got a %s cell.', what, ndraws, mat2str(size(x)));
end
L = zeros(ndraws, 2, 2);
for i = 1:2
    for j = 1:2
        L(:, i, j) = cell2mat(x(:, i, j));
    end
end
if ~isreal(L) || any(~isfinite(L(:)))
    error([ID ':invalid_correlation_draws'], '%s contains non-finite or complex draws.', what);
end
L11 = L(:,1,1); L12 = L(:,1,2); L21 = L(:,2,1); L22 = L(:,2,2);
dev = max([max(abs(L11 - 1)), max(abs(L12)), max(max(abs(L21) - 1, 0)), max(max(-L22, 0)), max(abs(L21.^2 + L22.^2 - 1))]);
if dev > tol
    error([ID ':invalid_correlation_draws'], ...
        ['%s is not a valid Cholesky correlation factor per draw (max structural deviation %.3g > %g): ', ...
        'need L11 = 1, L12 = 0, |L21| <= 1, L22 >= 0, L21^2 + L22^2 = 1.'], what, dev, tol);
end
end

% =========================================================================
%                          PARAMETERIZATION PROVENANCE
% =========================================================================
function param = local_parameterization(rel, hasStamp, ID)
% Collect every available identifier; all must agree. No default is assumed.
evidence = struct('source', {}, 'param', {});
if isfield(rel, 'gammascale')
    g = rel.gammascale;
    if ~isnumeric(g) || ~isscalar(g) || ~isfinite(g) || ~ismember(g, [1 2])
        error([ID ':unknown_parameterization'], ...
            'REL.gammascale = %s is not a recognized scale parameterization (1 = log-nu, 2 = location-scale).', mat2str(g));
    end
    if g == 2, p = 'ls'; else, p = 'lognu'; end
    evidence(end+1) = struct('source', 'REL.gammascale', 'param', p); %#ok<AGROW>
end
out = rel.out;
if isfield(out, 'gamma_disp') && isstruct(out.gamma_disp)
    gd = out.gamma_disp;
    isLs = (isfield(gd, 'is_ls') && isequal(gd.is_ls, true)) || isfield(gd, 'log_sigma') || isfield(gd, 's_sd') || isfield(gd, 'rho_ps');
    isNu = isfield(gd, 'nu') || isfield(gd, 's_v') || isfield(gd, 'V_w') || isfield(gd, 'rho_pv');
    if isLs && isNu
        error([ID ':contradictory_parameterization'], 'REL.out.gamma_disp carries both location-scale and log-nu provenance fields.');
    elseif isLs
        evidence(end+1) = struct('source', 'gamma_disp', 'param', 'ls'); %#ok<AGROW>
    elseif isNu
        evidence(end+1) = struct('source', 'gamma_disp', 'param', 'lognu'); %#ok<AGROW>
    end
end
hasNuRaw = isfield(out, 'pop_lognu') || isfield(out, 'ind_nu') || isfield(out, 'b_nu');
hasLsRaw = isfield(out, 'b_sigma');
if hasNuRaw && hasLsRaw
    error([ID ':contradictory_parameterization'], 'REL.out carries both log-nu (pop_lognu/b_nu/ind_nu) and location-scale (b_sigma) raw fields.');
elseif hasNuRaw
    evidence(end+1) = struct('source', 'raw fields', 'param', 'lognu'); %#ok<AGROW>
elseif hasLsRaw
    evidence(end+1) = struct('source', 'raw fields', 'param', 'ls'); %#ok<AGROW>
end
if hasStamp && isfield(out.gamma_target, 'parameterization')
    p = local_str(out.gamma_target, 'parameterization');
    if ~ismember(p, {'lognu', 'ls'})
        error([ID ':contradictory_stamp'], 'gamma_target.parameterization = ''%s'' is not recognized.', p);
    end
    evidence(end+1) = struct('source', 'gamma_target', 'param', p); %#ok<AGROW>
end
if isfield(out, 'gamma_resid') && iscell(out.gamma_resid)
    for s = 1:numel(out.gamma_resid)
        g = out.gamma_resid{s};
        if isstruct(g) && isfield(g, 'parameterization')
            evidence(end+1) = struct('source', sprintf('gamma_resid{%d}', s), 'param', local_str(g, 'parameterization')); %#ok<AGROW>
        end
    end
end
if isempty(evidence)
    error([ID ':ambiguous_parameterization'], ...
        ['Cannot determine the gamma scale parameterization: REL.gammascale is absent and no provenance ', ...
        '(gamma_disp, raw scale fields, stamp) disambiguates log-nu from location-scale. A default is not ', ...
        'assumed; how to fix: set REL.gammascale (1 = log-nu, 2 = location-scale) from the run''s provenance, or refit.']);
end
params = {evidence.param};
if numel(unique(params)) > 1
    parts = arrayfun(@(e) sprintf('%s=%s', e.source, e.param), evidence, 'UniformOutput', false);
    error([ID ':contradictory_parameterization'], ...
        'Scale parameterization provenance is contradictory (%s); refusing rather than restamping.', strjoin(parts, ', '));
end
param = params{1};
end

% =========================================================================
%                              SHARED CHECKS
% =========================================================================
function [ndraws, nstrata] = local_mu_shape(out, ID)
mu = out.mu;
if ~isnumeric(mu) || ~isreal(mu) || isempty(mu) || ndims(mu) ~= 2
    error([ID ':misaligned_draws'], 'REL.out.mu must be a nonempty real numeric [draws x strata] matrix.');
end
[ndraws, nstrata] = size(mu);
end

function local_check_labels(out, nstrata, ID)
lab = out.labels;
if ~iscell(lab) || numel(lab) ~= nstrata
    error([ID ':invalid_strata'], ...
        'REL.out.labels must hold exactly one label per stratum (%d); found %s of size %s.', nstrata, class(lab), mat2str(size(lab)));
end
for s = 1:nstrata
    if ~(ischar(lab{s}) || isstring(lab{s})) || isempty(char(lab{s}))
        error([ID ':invalid_strata'], 'REL.out.labels{%d} must be a nonempty text label.', s);
    end
end
end

function local_check_counts(out, name, nstrata, ID)
x = out.(name);
if ~isnumeric(x) || ~isreal(x) || numel(x) ~= nstrata || any(~isfinite(x(:)) | x(:) <= 0)
    error([ID ':invalid_strata'], ...
        'REL.out.%s must hold exactly one finite positive count per stratum (%d); got %s.', name, nstrata, mat2str(x));
end
end

function local_require(out, names, ID, what)
missing = names(~isfield(out, names));
if ~isempty(missing)
    error([ID ':missing_field'], ...
        'REL.out is missing required field(s) for %s: %s. The saved draws are insufficient; refit is required.', ...
        what, strjoin(missing, ', '));
end
for i = 1:numel(names)
    if isempty(out.(names{i}))
        error([ID ':missing_field'], 'REL.out.%s is empty (%s).', names{i}, what);
    end
end
end

function local_check_matrix(out, name, ndraws, nstrata, ID)
x = out.(name);
if ~isnumeric(x) || ~isreal(x) || ~isequal(size(x), [ndraws nstrata])
    error([ID ':misaligned_draws'], ...
        'REL.out.%s is %s; expected a real numeric [%d draws x %d strata] matrix aligned with REL.out.mu.', ...
        name, mat2str(size(x)), ndraws, nstrata);
end
end

function local_check_cellcount(out, name, nstrata, ID)
x = out.(name);
if ~iscell(x) || numel(x) ~= nstrata
    error([ID ':misaligned_draws'], ...
        'REL.out.%s must be a cell with one entry per stratum (%d); found %s of size %s.', name, nstrata, class(x), mat2str(size(x)));
end
end

function local_finite(x, what, ID, kind)
% kind: 'any' (finite real), 'nonneg' (finite real >= 0), 'pos' (finite real > 0)
if ~isnumeric(x) || ~isreal(x) || any(~isfinite(x(:)))
    error([ID ':invalid_source_values'], '%s must be finite real numeric (NaN/Inf/complex found).', what);
end
switch kind
    case 'nonneg'
        if any(x(:) < 0)
            error([ID ':invalid_source_values'], '%s must be nonnegative (a negative value was found).', what);
        end
    case 'pos'
        if any(x(:) <= 0)
            error([ID ':invalid_source_values'], '%s must be strictly positive.', what);
        end
end
end

function local_check_ids(ids, nsub, what, ID)
% subject labels: exactly one per row ([nsub x 1]), a supported type (real
% numeric column, cellstr column, string column or categorical column),
% non-missing and unique within the stratum. The same label may recur in
% another stratum; the key is stratum plus label.
if ~isequal(size(ids), [nsub 1])
    error([ID ':invalid_ids'], '%s must be one ID per row ([%d x 1]); got %s.', what, nsub, mat2str(size(ids)));
end
if isnumeric(ids)
    if ~isreal(ids) || any(~isfinite(ids))
        error([ID ':invalid_ids'], '%s contains missing/non-finite subject IDs.', what);
    end
    key = double(ids);
elseif iscellstr(ids) || isstring(ids) || iscategorical(ids) %#ok<ISCLSTR>
    key = string(ids);
    if any(ismissing(key)) || any(strlength(key) == 0)
        error([ID ':invalid_ids'], '%s contains missing/empty subject IDs.', what);
    end
else
    error([ID ':invalid_ids'], ...
        '%s has an unsupported ID type (%s); supported: numeric, cellstr, string or categorical column.', what, class(ids));
end
if numel(unique(key)) ~= numel(key)
    error([ID ':invalid_ids'], '%s contains duplicate subject IDs within the stratum.', what);
end
end

function local_check_id2(id2, nsub, what, ID)
% draw-column index: a [nsub x 1] real integer column forming a complete
% permutation of 1:nsub (never flattened from a wider block)
if ~isnumeric(id2) || ~isreal(id2) || ~isequal(size(id2), [nsub 1]) || any(~isfinite(id2)) ...
        || any(id2 ~= round(id2)) || ~isequal(sort(double(id2))', 1:nsub)
    error([ID ':invalid_id_matches'], ...
        '%s must be a [%d x 1] complete permutation of 1:%d (draw-column index, one per subject).', what, nsub, nsub);
end
end

function local_check_column(x, nsub, what, ID, kind)
% ssinfo numeric columns: exactly [nsub x 1], finite real; kind 'pos' also > 0
if ~isnumeric(x) || ~isreal(x) || ~isequal(size(x), [nsub 1]) || any(~isfinite(x))
    error([ID ':invalid_ssinfo'], '%s must be a finite real [%d x 1] column; got %s %s.', what, nsub, class(x), mat2str(size(x)));
end
if strcmp(kind, 'pos') && any(x <= 0)
    error([ID ':invalid_ssinfo'], '%s must be strictly positive counts.', what);
end
end

function m = local_cellmat(x, what, ID)
% Storage contract: every per-stratum draw block the summary/table/provenance
% consumers read through cell2mat (gro_sds, b, b_nu/b_sigma, chol_corrmat,
% ind_bs, ind_nu, ind_sdlog) must be stored EXACTLY as the estimator stores it,
% a num2cell block of real scalar draws. A numeric inner block is refused as
% unsupported storage rather than silently accepted: the migrated result
% would then fail to open in the consumers, and nothing is normalized here.
if ~iscell(x)
    error([ID ':unsupported_storage'], ...
        ['%s is stored as %s; the toolbox storage contract for this field is a num2cell block of ', ...
        'scalar draws (as psyrat_computevarcomp writes it). Refusing: consumers read this field ', ...
        'through cell2mat and would fail on the migrated result.'], what, class(x));
end
if isempty(x) || ~all(cellfun(@(v) isnumeric(v) && isreal(v) && isscalar(v), x(:)))
    error([ID ':unsupported_storage'], '%s must be a nonempty num2cell block of real scalar draws.', what);
end
m = cell2mat(x);
if ~isnumeric(m) || isempty(m) || ~isreal(m)
    error([ID ':unsupported_storage'], '%s did not concatenate into a real numeric block.', what);
end
end

function s = local_str(x, f)
if isstruct(x) && isfield(x, f)
    v = x.(f);
    if ischar(v), s = v;
    elseif isstring(v) && isscalar(v), s = char(v);
    elseif isnumeric(v) || islogical(v), s = mat2str(v);
    elseif iscell(v), s = sprintf('cell(%s)', mat2str(size(v)));
    else, s = class(v);
    end
else
    s = '(absent)';
end
end

function n = local_count(out, f)
if isfield(out, f), n = numel(out.(f)); else, n = 0; end
end

function lab = local_label(labels, s)
lab = '';
if iscell(labels) && numel(labels) >= s
    lab = char(string(labels{s}));
end
end

function v = local_toolbox_version()
v = 'unavailable';
if exist('psyrat_defineversion', 'file')
    try
        v = psyrat_defineversion();
    catch
    end
end
end

function p = local_fullpath(p)
p = char(p);
if isempty(regexp(p, '^(/|[A-Za-z]:[\\/]|\\\\)', 'once'))
    p = fullfile(pwd, p);
end
end

function r = local_conform(r, template)
tf = fieldnames(template);
for i = 1:numel(tf)
    if ~isfield(r, tf{i}), r.(tf{i}) = []; end
end
rf = fieldnames(r);
extra = rf(~ismember(rf, tf));
if ~isempty(extra), r = rmfield(r, extra); end
r = orderfields(r, template);
end
