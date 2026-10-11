function lines = psyrat_gamma_target_lines(summ, isperpart)
%PSYRAT_GAMMA_TARGET_LINES Estimand caption lines for a gamma non-difference
% dynamic reliability summary (analyses 11/26/14/27), read from the explicit
% summary.gamma_target contract (Q9/GUI-1 and Q8, 2026-09-13).
%
%   lines = psyrat_gamma_target_lines(summ, isperpart)
%
% summ      - psyrat_dynrel_summary output. Returns {} (nothing to say) when
%             summ carries no gamma_target, so callers can append blindly.
% isperpart - true when the caption accompanies the PER-PARTICIPANT table
%             (adds the subject-target and remainder-scheme lines); false for
%             the surface / variance-component outputs.
%
% Output: a column cellstr. Line 1 names the surface target and its residual
% construction; line 2 is the ACTIVE summary contract as literal key = value
% identifiers (gamma_target_version, surface_target, subject_target,
% remainder_scheme, parameterization, count_source, surface_count_source,
% reference_target, analysis, surface_target_source, subject_scheme_source),
% emitted for every covered design and output with not_applicable where a
% design has no subject target; for the subject-level designs line 3 names
% the conditioning set and remainder scheme (including the LEGACY
% observation-only scheme of an unstamped analysis-27 result, which is
% reported as such and never as v2) and line 4 the subject count source. The
% active contract reflects view-time selections (surface target, n'
% override, explicit scheme); the STORED-result contract the file carries is
% printed separately by psyrat_provenance_lines, so the two are
% distinguishable when they differ. One source of truth for the viewer export
% header, the headless report and the figure caption, so they cannot drift.
%
% See also PSYRAT_DYNREL_SUMMARY, PSYRAT_STARTVIEW_DYNREL, PSYRAT_REPORT,
% PSYRAT_DYNRELPLOT.

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

lines = {};
if nargin < 2; isperpart = false; end
if ~isstruct(summ) || ~isfield(summ,'gamma_target') || isempty(summ.gamma_target)
    return;
end
gt = summ.gamma_target;
isls = strcmp(gt.parameterization,'ls');
twofacet = any(strcmp(gt.analysis,{'ic_dynrel_trt','ic_dynrel_sserrvar_trt'}));

%---- surface target ---------------------------------------------------------
if twofacet
    means = ['the six population mean components at z (person, occasion, trial ' ...
        'and the three two-way interactions)'];
    mterm = 'the pooled population highest-order (person x occasion x trial) mean remainder R';
else
    means = 'the population person and trial main-effect components at z';
    mterm = 'the pooled population person x trial mean-surface term';
end
if isls
    cpop = 'exp(2 log sigma(z) + 2 s_sd^2)';
    ctyp = 'exp(2 log sigma(z)) (w_p = 0)';
    over = 'the person log residual-SD distribution';
else
    cpop = '2 exp(2 alpha(z) + 2 V - log nu(z) + s_v^2/2 - 2 cov_pv)';
    ctyp = '2 exp(2 alpha(z) + 2 W - log nu(z)) (u_p = v_p = 0)';
    over = 'the JOINT person location (log-mean) and log-nu distribution';
end
%how the companion relates to the primary differs by parameterization: under
%log-nu BOTH person effects must be zero (setting s_v = 0 alone leaves the
%person log-mean integrated), whereas under location-scale the primary is
%exp(2 log sigma + 2 s_sd^2) and reduces to the companion exactly when
%s_sd = 0 - the zero-scale-heterogeneity equality control.
if isls
    relation = ['The observation-variance terms coincide exactly when the ' ...
        'person log residual-SD SD s_sd is zero; the surfaces therefore ' ...
        'coincide then (the zero-scale-heterogeneity equality control).'];
else
    relation = ['Not the population-marginal primary, and not obtained by ' ...
        'zeroing the person log-nu variance alone: that would still integrate ' ...
        'over the person log-mean. The observation-variance terms coincide ' ...
        'exactly when 4 s_p^2 + s_v^2 - 4 cov_pv = 0; the surfaces ' ...
        'therefore coincide then.'];
end
switch gt.surface_target
    case 'typical_person'
        lines{end+1,1} = ['Estimand: surface = TYPICAL PERSON (companion): ' ...
            means ' plus ' mterm ', with the observation variance at zero ' ...
            'person effects, ' ctyp '. ' relation];
    otherwise
        lines{end+1,1} = ['Estimand: surface = POPULATION-MARGINAL (primary): ' ...
            means ' plus ' mterm ', with the observation variance ' ...
            'marginalized over ' over ' on each draw before the coefficient, ' ...
            cpop '. This is not an average of per-participant coefficients.'];
end
%line 2 (index kept stable: the per-participant note is line 3): the ACTIVE
%summary contract as literal key = value identifiers, emitted for every
%covered design and every output (variance table included), with
%not_applicable where a design has no subject target. This is what was
%summarized NOW (view-time selections applied); the stored-result contract
%the file carries is printed separately by psyrat_provenance_lines.
lines{end+1,1} = sprintf(['Active summary contract: gamma_target_version = %d; ' ...
    'surface_target = %s; subject_target = %s; remainder_scheme = %s; ' ...
    'parameterization = %s; count_source = %s; surface_count_source = %s; ' ...
    'reference_target = %s; analysis = %s; surface_target_source = %s; ' ...
    'subject_scheme_source = %s'], ...
    gt.gamma_target_version, gt.surface_target, gt.subject_target, ...
    gt.remainder_scheme, gt.parameterization, gt.count_source, ...
    gt.surface_count_source, local_opt_field(gt,'reference_target','not_applicable'), ...
    gt.analysis, local_opt_field(gt,'surface_target_source','not_recorded'), ...
    gt.subject_scheme_source);

%---- per-participant target --------------------------------------------------
if isperpart && ~strcmp(gt.subject_target,'not_applicable')
    if isls
        cond = 'their own log residual SD';
    else
        cond = 'their own log-mean and log-nu';
    end
    switch gt.remainder_scheme
        case 'legacy_observation_only_v1'
            if strcmp(gt.subject_scheme_source,'explicit_override')
                why = ['the legacy observation-only scheme was explicitly ' ...
                    'requested for this stamped result'];
            else
                why = ['this saved result carries no gamma_target stamp and has ' ...
                    'not been migrated'];
            end
            sch = ['LEGACY observation-only residual (' why '): each ' ...
                'participant''s error is their own observation variance ALONE; ' ...
                'the population highest-order mean remainder R the surface ' ...
                'carries is OMITTED from it, so these coefficients can be larger ' ...
                'than under scheme v2 (they coincide when R is zero, so a ' ...
                'zero-facet control provides an equality check), and an ' ...
                'equal-count zero-effect ' ...
                'participant need not sit on a reference surface'];
        case 'population_highest_order_plus_person_observation_v2'
            sch = ['residual = pooled population highest-order (person x ' ...
                'occasion x trial) mean remainder at their own z + their own ' ...
                'observation variance (scheme v2); an equal-count participant ' ...
                'with zero person effects sits on the typical-person surface'];
        otherwise
            sch = ['residual = pooled population person x trial mean-surface ' ...
                'term at their own z + their own observation variance (S20 ' ...
                'pooling, scheme v1); an equal-count participant with zero ' ...
                'person effects sits on the typical-person surface'];
    end
    lines{end+1,1} = ['Estimand: per-participant table = MIXED population/person ' ...
        'target: ' means ' are retained as population quantities at each ' ...
        'participant''s own dimension value(s); the participant''s own effects ' ...
        '(' cond ') enter ONLY the observation variance; ' sch '. It is not a ' ...
        'fully person-conditioned decomposition of the mean surface. The ' ...
        'surface target selector does not change these coefficients.'];
    lines{end+1,1} = sprintf(['Per-participant metadata: remainder_scheme = %s ' ...
        '(source: %s); subject n'' = %s; selected occasion count as shown ' ...
        'per stratum'], gt.remainder_scheme, ...
        strrep(gt.subject_scheme_source,'_',' '), strrep(gt.count_source,'_',' '));
end

end

function v = local_opt_field(gt, f, fallback)
%Text of an optional carrier field; the fallback names its absence rather
%than inventing a value.
if isfield(gt,f) && ~isempty(gt.(f))
    v = char(string(gt.(f)));
else
    v = fallback;
end
end
