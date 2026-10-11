function v = psyrat_gamma_target_validate(REL)
%PSYRAT_GAMMA_TARGET_VALIDATE Strict, non-throwing validation of a gamma
% result's REL.gammascale and REL.out.gamma_target marker against the schema
% its design implies (Q8/Q9, gamma_target_version 2). Shared by the summary
% resolver (psyrat_gamma_target_resolve, which turns a failure into the
% precise error) and by psyrat_provenance_lines (which must never print
% POOLED/v2 or LEGACY from a marker it cannot interpret).
%
%   v = psyrat_gamma_target_validate(REL)
%
% v.status   'absent'  - no gamma_target field at all: a historical, unstamped
%                        result (the ONLY case that means legacy)
%            'valid'   - a 1x1 struct marker with gamma_target_version 2, every
%                        contract field present, and every identifier coherent
%                        with this result's design and parameterization
%            'invalid' - anything else (present but empty, non-struct,
%                        non-scalar, unsupported version, missing field,
%                        contradictory value, or an uninterpretable gammascale)
% v.errid    error identifier the resolver raises for 'invalid'
%            ('dynrel:gammatargetmalformed' | 'dynrel:gammatargetversion' |
%             'dynrel:gammatargetincoherent' | 'dynrel:badgammascale'); ''
% v.reason   one-sentence explanation ('' unless invalid)
% v.stamp    the marker struct when valid ([] otherwise)
% v.parameterization   'lognu' | 'ls' resolved from REL.gammascale (a MISSING
%            gammascale is the documented historical fallback to log-nu; a
%            PRESENT gammascale must be a finite real scalar 1 or 2)
% v.design   struct of the identifiers the design implies (expected values),
%            or [] when the analysis is not a gamma design this schema covers
% v.covered  true when REL is a gamma result of analysis trt_sserrvar (static
%            25) or ic_dynrel / ic_dynrel_sserrvar / ic_dynrel_trt /
%            ic_dynrel_sserrvar_trt (dynamic 11/26/14/27)
%
% SCHEMA (all values exact):
%   static 25:  surface_target 'not_applicable'; reference_target
%               'typical_person_zero_person_effects'; subject_target 'own_...';
%               remainder_scheme 'population_highest_order_plus_person_observation_v2';
%               count_source 'mean_trials_per_occupied_occasion';
%               surface_count_source 'not_applicable'; analysis 'trt_sserrvar'
%   dynamic:    surface_target 'population_marginal' | 'typical_person';
%               reference_target 'not_applicable'; surface_count_source
%               'default_stratum_median'; subject_target 'own_...' (26/27) |
%               'not_applicable' (11/14); remainder_scheme v2 (27) |
%               'population_mean_remainder_plus_person_observation_v1' (26) |
%               'population_mean_remainder_plus_selected_surface_observation_v1'
%               (11/14); count_source 'median_trials_per_occupied_occasion'
%               (26/27) | 'not_applicable' (11/14); analysis = REL.analysis
%   parameterization must equal the value resolved from REL.gammascale.
%
% See also PSYRAT_GAMMA_TARGET_RESOLVE, PSYRAT_PROVENANCE_LINES.

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

v = struct('status','absent','errid','','reason','','stamp',[],...
    'parameterization','lognu','design',[],'covered',false);

analysis = '';
if isstruct(REL) && isfield(REL,'analysis') && ~isempty(REL.analysis)
    analysis = char(string(REL.analysis));
end
isgamma = isstruct(REL) && isfield(REL,'family') && ...
    strcmp(char(string(REL.family)),'gamma');

%---- parameterization from gammascale: missing -> log-nu (historical fallback);
%---- present -> must be exactly a finite real scalar 1 or 2
if isstruct(REL) && isfield(REL,'gammascale')
    gs = REL.gammascale;
    if ~(isnumeric(gs) && isscalar(gs) && isreal(gs) && isfinite(gs) && any(gs == [1 2]))
        v.status = 'invalid';
        v.errid = 'dynrel:badgammascale';
        v.reason = sprintf(['REL.gammascale is present but is not exactly 1 (log-nu) ' ...
            'or 2 (location-scale): %s. A present gammascale must be a finite ' ...
            'real scalar 1 or 2; only an ABSENT gammascale falls back to log-nu.'], ...
            local_describe(gs));
        return;
    end
    if gs == 2
        v.parameterization = 'ls';
    end
end

%---- design schema ---------------------------------------------------------
own = 'own_observation_variance_with_population_mean_components_at_own_z';
v2  = 'population_highest_order_plus_person_observation_v2';
switch analysis
    case 'trt_sserrvar'
        d = struct('surface_target', {{'not_applicable'}}, ...
            'reference_target', 'typical_person_zero_person_effects', ...
            'subject_target', own, 'remainder_scheme', v2, ...
            'count_source', 'mean_trials_per_occupied_occasion', ...
            'surface_count_source', 'not_applicable', 'analysis', analysis);
    case 'ic_dynrel_sserrvar_trt'
        d = local_dynamic_design(analysis, own, v2, 'median_trials_per_occupied_occasion');
    case 'ic_dynrel_sserrvar'
        d = local_dynamic_design(analysis, own, ...
            'population_mean_remainder_plus_person_observation_v1', ...
            'median_trials_per_occupied_occasion');
    case {'ic_dynrel','ic_dynrel_trt'}
        d = local_dynamic_design(analysis, 'not_applicable', ...
            'population_mean_remainder_plus_selected_surface_observation_v1', ...
            'not_applicable');
    otherwise
        d = [];
end
v.covered = isgamma && ~isempty(d);
v.design = d;
if ~v.covered
    return;   % status 'absent': nothing to validate for this design
end

%---- marker: absent vs present ----------------------------------------------
if ~(isfield(REL,'out') && isstruct(REL.out) && isfield(REL.out,'gamma_target'))
    return;   % 'absent' - historical unstamped result
end
st = REL.out.gamma_target;
if ~(isstruct(st) && isscalar(st))
    v.status = 'invalid';
    v.errid = 'dynrel:gammatargetmalformed';
    v.reason = sprintf(['REL.out.gamma_target is present but is not a 1x1 struct ' ...
        '(%s). A present but empty or malformed marker is not evidence of a ' ...
        'historical result and is not interpreted.'], local_describe(st));
    return;
end
if ~isfield(st,'gamma_target_version')
    v.status = 'invalid';
    v.errid = 'dynrel:gammatargetversion';
    v.reason = ['REL.out.gamma_target carries no gamma_target_version, so it ' ...
        'cannot be trusted.'];
    return;
end
ver = st.gamma_target_version;
if ~(isnumeric(ver) && isscalar(ver) && isreal(ver) && isequal(ver,2))
    v.status = 'invalid';
    v.errid = 'dynrel:gammatargetversion';
    v.reason = sprintf(['REL.out.gamma_target.gamma_target_version is %s; this ' ...
        'toolbox understands version 2 only, and a versioned marker of unknown ' ...
        'meaning is not read as an unstamped legacy result.'], local_describe(ver));
    return;
end
required = {'surface_target','reference_target','subject_target',...
    'remainder_scheme','parameterization','count_source',...
    'surface_count_source','analysis'};
for k = 1:numel(required)
    if ~isfield(st,required{k})
        v.status = 'invalid';
        v.errid = 'dynrel:gammatargetincoherent';
        v.reason = sprintf(['REL.out.gamma_target (version 2) lacks the field ' ...
            '''%s''; a partial marker is not interpreted.'], required{k});
        return;
    end
end
checks = {
    'analysis',             d.analysis
    'parameterization',     v.parameterization
    'reference_target',     d.reference_target
    'subject_target',       d.subject_target
    'remainder_scheme',     d.remainder_scheme
    'count_source',         d.count_source
    'surface_count_source', d.surface_count_source
    };
for k = 1:size(checks,1)
    f = checks{k,1};
    stored = local_str(st.(f));
    if ~strcmp(stored, checks{k,2})
        v.status = 'invalid';
        v.errid = 'dynrel:gammatargetincoherent';
        v.reason = sprintf(['REL.out.gamma_target.%s is ''%s'' but this result ' ...
            'implies ''%s''.'], f, stored, checks{k,2});
        return;
    end
end
stored = local_str(st.surface_target);
if ~any(strcmp(stored, d.surface_target))
    v.status = 'invalid';
    v.errid = 'dynrel:gammatargetincoherent';
    v.reason = sprintf(['REL.out.gamma_target.surface_target is ''%s'' but this ' ...
        'result (%s) allows only: %s.'], stored, analysis, strjoin(d.surface_target,', '));
    return;
end
v.status = 'valid';
v.stamp = st;

end

function d = local_dynamic_design(analysis, subject_target, scheme, count_source)
d = struct('surface_target', {{'population_marginal','typical_person'}}, ...
    'reference_target', 'not_applicable', ...
    'subject_target', subject_target, 'remainder_scheme', scheme, ...
    'count_source', count_source, ...
    'surface_count_source', 'default_stratum_median', 'analysis', analysis);
end

function s = local_str(x)
%Text of a stored identifier; non-text values render as a description so the
%reason names what was found instead of erroring inside the validator.
if ischar(x) || (isstring(x) && isscalar(x))
    s = char(string(x));
else
    s = local_describe(x);
end
end

function s = local_describe(x)
if isempty(x)
    s = sprintf('empty %s', class(x));
elseif isnumeric(x) && isscalar(x)
    s = mat2str(x);
elseif ischar(x) || isstring(x)
    s = sprintf('text ''%s''', char(strjoin(string(x),' ')));
else
    s = sprintf('%s of size %s', class(x), mat2str(size(x)));
end
end
