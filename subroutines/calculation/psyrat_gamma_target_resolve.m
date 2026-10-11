function [target, gt] = psyrat_gamma_target_resolve(REL, target_in, scheme_in, obs_in)
%PSYRAT_GAMMA_TARGET_RESOLVE Resolve the explicit Q9 SURFACE target and the Q8
% per-participant REMAINDER scheme for a gamma non-difference dynamic result
% (analyses 11/26/14/27), validating any stored REL.out.gamma_target stamp, and
% build the summary-level gamma_target carrier (the shared metadata contract
% read by the summary, viewer, plot, report, exports and provenance).
%
%   [target, gt] = psyrat_gamma_target_resolve(REL, target_in, scheme_in, obs_in)
%
% Inputs
%   REL       - the result struct (REL.analysis, REL.family, REL.gammascale,
%               REL.out.gamma_target when stamped)
%   target_in - explicit 'surface_target' option ([] = not given)
%   scheme_in - explicit 'remainder_scheme' option ([] = not given)
%   obs_in    - the surface n' request ([] = per-stratum median default)
% Outputs
%   target - 'population_marginal' | 'typical_person' ('' for other designs)
%   gt     - the carrier ([] for other designs):
%     gamma_target_version   2
%     surface_target         'population_marginal' | 'typical_person'
%     surface_target_source  'explicit_option' | 'stored_stamp' | 'historical_default'
%     subject_target         'own_observation_variance_with_population_mean_components_at_own_z' (26/27) | 'not_applicable' (11/14)
%     remainder_scheme       27: 'population_highest_order_plus_person_observation_v2' | 'legacy_observation_only_v1'
%                            26: 'population_mean_remainder_plus_person_observation_v1'
%                            11/14: 'population_mean_remainder_plus_selected_surface_observation_v1'
%     parameterization       'lognu' | 'ls'
%     count_source           'median_trials_per_occupied_occasion' (26/27) | 'not_applicable' (11/14)
%     surface_count_source   'default_stratum_median' | 'user_override'
%     analysis               REL.analysis
%     subject_scheme_source  'stamped_v2' | 'legacy_unstamped' | 'explicit_override' | 'unstamped_s20_convention' | 'not_applicable'
%     stamped                true when a VALID version-2 stamp is present
%
% PRECEDENCE of the surface target: an explicit option > a valid stored
% REL.out.gamma_target.surface_target > the historical population-marginal
% default. Legacy gamma results carry no stamp; their marginal Boolean (0 or 1)
% never selected anything and is NOT reinterpreted as a target.
%
% STORED STAMP VALIDATION (psyrat_gamma_target_validate; fail loudly, never
% silently downgrade to legacy): a present marker must be a 1x1 struct with
% gamma_target_version 2 and every contract field (surface_target,
% reference_target, subject_target, remainder_scheme, parameterization,
% count_source, surface_count_source, analysis), each coherent with THIS
% result's design and with REL.gammascale (dynamic: reference_target
% 'not_applicable', surface_count_source 'default_stratum_median',
% surface_target population_marginal | typical_person). A present but empty
% or malformed marker ('dynrel:gammatargetmalformed'), an unsupported version
% ('dynrel:gammatargetversion'), a partial or contradictory marker
% ('dynrel:gammatargetincoherent') and a present gammascale that is not
% exactly 1 or 2 ('dynrel:badgammascale') are all refused; only the ABSENCE
% of the marker means "historical, unstamped" (and only an absent gammascale
% falls back to log-nu). The summary-level surface_count_source may say
% 'user_override' for the current request; that never bypasses the stored
% validation, which checks the stamp's own 'default_stratum_median'.
%
% RULES. (1) An explicit target or scheme on any other design (Gaussian,
% difference, copula, DoD) is an error, never a silent no-op. (2) An unstamped
% analysis-27 result gets the LEGACY observation-only subject scheme until an
% explicit migration stamps it; requesting v2 on unstamped draws is refused
% (that is the migration's job). (3) The surface target never changes the
% subject scheme, so per-participant coefficients are fixed when only the
% surface target moves. (4) Analysis 26 has exactly one scheme (S20).
%
% See also PSYRAT_DYNREL_SUMMARY, PSYRAT_GAMMA_SURFACE_TARGET,
% PSYRAT_STARTVIEW_DYNREL, PSYRAT_GAMMA_TARGET_LINES.

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

if nargin < 2; target_in = []; end
if nargin < 3; scheme_in = []; end
if nargin < 4; obs_in = []; end

analysis = '';
if isfield(REL,'analysis') && ~isempty(REL.analysis)
    analysis = char(string(REL.analysis));
end
isgamma = isfield(REL,'family') && strcmp(char(string(REL.family)),'gamma');
issubjnondiff     = strcmp(analysis,'ic_dynrel_sserrvar');
issubjnondiff_trt = strcmp(analysis,'ic_dynrel_sserrvar_trt');
isnondiffgamma = isgamma && any(strcmp(analysis,{'ic_dynrel','ic_dynrel_sserrvar',...
    'ic_dynrel_trt','ic_dynrel_sserrvar_trt'}));

if ~isnondiffgamma
    if ~isempty(target_in)
        error('dynrel:surfacetargetunsupported',...
            ['''surface_target'' applies only to the gamma non-difference dynamic ' ...
            'designs (ic_dynrel, ic_dynrel_sserrvar, ic_dynrel_trt, ' ...
            'ic_dynrel_sserrvar_trt); this result is analysis ''%s'' under ' ...
            'family ''%s''. How to fix: omit ''surface_target'' (Gaussian ' ...
            'non-difference runs select their estimand with ''marginal'').'],...
            analysis, local_family_name(REL));
    end
    if ~isempty(scheme_in)
        error('dynrel:remainderschemeunsupported',...
            ['''remainder_scheme'' applies only to the gamma non-difference ' ...
            'dynamic designs; this result is analysis ''%s'' under family ' ...
            '''%s''. How to fix: omit ''remainder_scheme''.'],...
            analysis, local_family_name(REL));
    end
    target = '';
    gt = [];
    return;
end

%---- strict validation of gammascale and any present marker ------------------
%(psyrat_gamma_target_validate holds the schema; this turns a failure into the
%precise error BEFORE any numerical work). A present-but-empty, malformed,
%unsupported-version, partial or contradictory marker is refused; only an
%ABSENT marker means historical/unstamped. A present gammascale must be exactly
%1 or 2; only an absent one falls back to log-nu.
vv = psyrat_gamma_target_validate(REL);
if strcmp(vv.status,'invalid')
    error(vv.errid, ['%s How to fix: re-run the versioned migration on the ' ...
        'saved result (it writes a coherent marker), restore the file, or ' ...
        'update PsyRAT if the marker is newer than this toolbox. The result ' ...
        'is not summarized under any scheme until then.'], vv.reason);
end
param = vv.parameterization;
v2 = 'population_highest_order_plus_person_observation_v2';
legacy = 'legacy_observation_only_v1';
if issubjnondiff_trt
    design_scheme = v2;
    subject_target = 'own_observation_variance_with_population_mean_components_at_own_z';
    count_source = 'median_trials_per_occupied_occasion';
elseif issubjnondiff
    design_scheme = 'population_mean_remainder_plus_person_observation_v1';
    subject_target = 'own_observation_variance_with_population_mean_components_at_own_z';
    count_source = 'median_trials_per_occupied_occasion';
else
    design_scheme = 'population_mean_remainder_plus_selected_surface_observation_v1';
    subject_target = 'not_applicable';
    count_source = 'not_applicable';
end

%---- stored stamp: absent -> unstamped/legacy; valid -> its saved target ----
stamped = strcmp(vv.status,'valid');
stored_target = '';
if stamped
    stored_target = char(string(vv.stamp.surface_target));
end

%---- surface target: explicit > stored > historical default ----------------
if ~isempty(target_in)
    target = psyrat_gamma_surface_target(target_in);   % validates
    target_source = 'explicit_option';
elseif stamped
    target = stored_target;
    target_source = 'stored_stamp';
else
    target = 'population_marginal';
    target_source = 'historical_default';
end

%---- per-participant remainder scheme -------------------------------------
if ~isempty(scheme_in) && ...
        ~(ischar(scheme_in) || (isstring(scheme_in) && isscalar(scheme_in)))
    error('dynrel:badremainderscheme',...
        ['''remainder_scheme'' must be a text scheme name. How to fix: pass ' ...
        'one of the documented scheme strings or omit the option.']);
end
if ~isempty(scheme_in); scheme_in = char(scheme_in); end

if issubjnondiff_trt
    if stamped
        scheme = v2; src = 'stamped_v2';
    else
        scheme = legacy; src = 'legacy_unstamped';
    end
    if ~isempty(scheme_in)
        if ~any(strcmp(scheme_in,{v2,legacy}))
            error('dynrel:badremainderscheme',...
                ['Unknown remainder_scheme ''%s'' for analysis 27. Valid values are ' ...
                '''%s'' and ''%s''. How to fix: pass one of those or omit the option.'],...
                scheme_in, v2, legacy);
        end
        if strcmp(scheme_in,v2) && ~stamped
            error('dynrel:remainderschemeunstamped',...
                ['This analysis-27 result carries no gamma_target stamp (it was ' ...
                'saved before gamma_target_version 2), so its per-participant ' ...
                'residuals cannot be summarized under ''%s'' without an explicit ' ...
                'migration. How to fix: run the versioned migration on the saved ' ...
                'result (it stamps the draws and records the conversion), or omit ' ...
                '''remainder_scheme'' to keep the legacy observation-only scheme ' ...
                'the file was produced under.'], v2);
        end
        scheme = scheme_in; src = 'explicit_override';
    end
elseif issubjnondiff
    scheme = design_scheme;
    if stamped; src = 'stamped_v2'; else; src = 'unstamped_s20_convention'; end
    if ~isempty(scheme_in) && ~strcmp(scheme_in,scheme)
        error('dynrel:badremainderscheme',...
            ['Analysis 26 has a single remainder scheme, ''%s''; ''%s'' is not ' ...
            'valid for it. How to fix: omit ''remainder_scheme''.'],...
            scheme, scheme_in);
    end
else
    scheme = design_scheme;
    src = 'not_applicable';
    if ~isempty(scheme_in) && ~strcmp(scheme_in,scheme)
        error('dynrel:badremainderscheme',...
            ['The group-level designs 11/14 have a single remainder scheme, ' ...
            '''%s''; ''%s'' is not valid for them. How to fix: omit ' ...
            '''remainder_scheme''.'], scheme, scheme_in);
    end
end

if isempty(obs_in)
    surface_count_source = 'default_stratum_median';
else
    surface_count_source = 'user_override';
end

gt = struct( ...
    'gamma_target_version', 2, ...
    'surface_target', target, ...
    'surface_target_source', target_source, ...
    'subject_target', subject_target, ...
    'remainder_scheme', scheme, ...
    'parameterization', param, ...
    'count_source', count_source, ...
    'surface_count_source', surface_count_source, ...
    'analysis', analysis, ...
    'subject_scheme_source', src, ...
    'stamped', stamped);

end

function name = local_family_name(REL)
%Family string for error text; '(unrecorded)' when the result carries none.
if isfield(REL,'family') && ~isempty(REL.family)
    name = char(string(REL.family));
else
    name = '(unrecorded)';
end
end
