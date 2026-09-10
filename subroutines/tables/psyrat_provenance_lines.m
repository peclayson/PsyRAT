function lines = psyrat_provenance_lines(psyrat_data)
%Standard run-provenance lines for exported-table headers and the headless
%report writer.
%
%lines = psyrat_provenance_lines(psyrat_data)
%
%Returns a cell column of char rows recording the reproducibility-relevant
%settings for a run: the estimation seed, the engine that produced it, whether
%the priors were the PsyRAT defaults or customized, the convergence verdict,
%the CmdStan divergent-transition count, the credible-interval width the
%summary was built at (when a summary records one), and (test-retest only)
%which of the three coefficient estimands was selected. These are appended to
%the table-specific header block each exporter already builds, so a saved
%table carries enough provenance to be interpreted on its own.
%
%Every line is best-effort: a field that is absent (for example in a legacy
%saved file produced before a field existed) yields a "not recorded" or "not
%assessed" line rather than an error, so adding these lines can never break a
%table export.
%
%Input
% psyrat_data - toolbox data structure after estimation; reads psyrat_data.rel
%   (.seed, .analysis, .engine, .family, .dispersion, .priors, .out.conv.converged,
%   .out.conv.ndivergent) and, for the two-facet designs only,
%   psyrat_data.relsummary (.reltype, .gcoeff, .noccmode, .nocc). The
%   relsummary read is the sole exception to the rel-only rule below: reltype
%   is a view-time selection stored on the summary, so it exists nowhere on
%   rel, and it is gated so no other design is affected.
%
%Output
% lines - N-by-1 cell array of char strings.

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

lines = cell(0,1);
rel = local_rel(psyrat_data);

% Estimation seed.
if isfield(rel,'seed') && ~isempty(rel.seed)
    lines{end+1,1} = sprintf('Seed: %g', rel.seed);
else
    lines{end+1,1} = 'Seed: not recorded';
end

% Analysis design label. Without it, exports of the same dataset under
% different designs are indistinguishable on their face.
if isfield(rel,'analysis') && ~isempty(rel.analysis)
    lines{end+1,1} = sprintf('Analysis: %s', char(string(rel.analysis)));
else
    lines{end+1,1} = 'Analysis: not recorded';
end

% Estimation engine (cmdstan / hmc / fitlme). Defaults to cmdstan when the
% field is absent, matching the runconfig writer's convention for legacy runs.
if isfield(rel,'engine') && ~isempty(rel.engine)
    lines{end+1,1} = sprintf('Estimation engine: %s', char(rel.engine));
    % Within-chain threading, only when it was actually enabled: the threaded
    % models use reduce_sum, whose partitioning the scheduler chooses, so the
    % draws of such a run are not bit-reproducible and the header must say so
    % (audit 2026-09-05). Unthreaded runs print nothing here, so their headers
    % are unchanged.
    if isfield(rel,'withinchain_enabled') && isequal(rel.withinchain_enabled,true) && ...
            isfield(rel,'threads_per_chain') && ~isempty(rel.threads_per_chain)
        lines{end+1,1} = sprintf(['Within-chain threads per chain: %d ', ...
            '(reduce_sum partitioning; draws not bit-reproducible across runs)'], ...
            rel.threads_per_chain);
    end
else
    lines{end+1,1} = 'Estimation engine: cmdstan';
end

% Observation likelihood family. Gaussian is the historical default; a line is
% emitted only for a non-Gaussian family so existing Gaussian output stays
% byte-identical (absence/empty == Gaussian, matching the engine-defaults-to-
% cmdstan convention above). The gamma family stores observed-score-scale
% SD-equivalents in the standard REL.out slots, so without this label a reader
% cannot tell a gamma table from a Gaussian one. Normalize to a trimmed char
% first so a string-typed empty ("") is treated as absent (best-effort contract)
% rather than emitting a blank "Likelihood family:" line.
famstr = '';
if isfield(rel,'family') && ~isempty(rel.family)
    famstr = strtrim(char(string(rel.family)));
end
if ~isempty(famstr) && ~strcmpi(famstr,'gaussian')
    lines{end+1,1} = sprintf('Likelihood family: %s', local_family_label(famstr));
end

% Gamma dispersion-heterogeneity diagnostics. The marginal (estimand #2)
% observation term rescales the per-cell residual by F = exp(0.5*s_v^2 - 2*cov_pv),
% built from the person log-nu SD (s_v) and the person mean<->log-nu covariance.
% Surface s_v, the mean-dispersion correlation, and the implied F so a reader can
% tell whether the correction is data-driven or (with a weakly identified s_v)
% prior-driven. F near 1 means the marginal correction is negligible (dispersion
% roughly homogeneous), i.e. the population-nu plug-in would have been adequate.
if ~isempty(famstr) && strcmpi(famstr,'gamma')
    lines = [lines; local_gamma_dispersion_lines(rel)];
end

% Modular cut-copula disclosure (any gamma design storing REL.out.copula;
% today analyses 19/20). MUST be its own block, not a
% branch of local_gamma_dispersion_lines: that subfunction returns empty for
% these designs (they store no out.gamma/gamma_disp dispersion struct), so its
% is_diff_copula branch can never fire here. The cut-posterior status is the
% defining inferential fact of these designs -- rho_e is estimated conditionally
% on retained margin draws with no feedback, so intervals are not full
% joint-posterior credible intervals -- and the stage's header requires the
% labels to reach every output surface.
if ~isempty(famstr) && strcmpi(famstr,'gamma')
    lines = [lines; local_modular_cut_lines(rel)];
end

% Nonconcurrent-DoD disclosure (analyses 28/29). UNCONDITIONAL on family since
% the Gaussian arm shipped (2026-08-20 owner ruling): the block is
% presence-keyed on out.dod_labels and returns empty for every other run, so
% all other output stays byte-identical - but a family gate here would silently
% strip the counterfactual-estimand disclosure from a Gaussian DoD result, the
% exact class of drift the B15 review caught on the estimand note.
lines = [lines; local_nonconcurrent_dod_lines(rel)];

% One-facet trial main effect (sigma_i). Emitted ONLY when sigma_i is absent or
% identically zero, mirroring the emit-only-when-notable convention used for the
% family line above so existing output stays byte-identical in the normal case.
% It matters because absolute error is what separates dependability from
% generalizability: with sigma_i = 0 the two coefficients are equal to the last
% bit, so a "Dependability" heading on such a table is unearned and a reader
% cannot tell that from the number alone.
lines = [lines; local_sigtrl_lines(rel)];

% Test-retest coefficient estimand: which of CE / CS / CES was selected, the
% error type, and n'_o. Emitted only for designs that record a reltype, so
% one-facet, splits and dynamic output stays byte-identical (same emit-only-
% when-applicable convention as the family and sigma_i lines above). It matters
% because the three coefficients are not near-neighbors -- on one dataset the
% same components gave 0.40 / 0.35 / 0.24 -- and nothing else in an exported
% two-facet table records which one produced the number.
lines = [lines; local_coef_lines(psyrat_data)];

% Credible-interval width (G41). Every interval in a summary-derived table was
% computed at relsummary.ciperc, but no header line recorded it, so a saved
% table's intervals could not be interpreted on their own (and a
% psyrat_report run at a different 'CI' could put two widths in one export
% with nothing in the file saying so -- the report now discloses that case in
% its own header note). Emitted only when a built summary records the width
% (the same emit-only-when-applicable convention as the coefficient line
% above), so exports from structs with no relsummary are unchanged.
if isstruct(psyrat_data) && isfield(psyrat_data,'relsummary') && ...
        isfield(psyrat_data.relsummary,'ciperc') && ...
        isnumeric(psyrat_data.relsummary.ciperc) && ...
        isscalar(psyrat_data.relsummary.ciperc) && ...
        ~isnan(psyrat_data.relsummary.ciperc)
    lines{end+1,1} = sprintf('Credible interval: %g%%', ...
        100 * psyrat_data.relsummary.ciperc);
end

% Priors: default vs custom relative to the toolbox defaults.
lines{end+1,1} = local_priors_line(rel);

% Convergence verdict (0/1 stamped by psyrat_checkconv for all engines).
if local_has(rel,{'out','conv','converged'}) && ~isempty(rel.out.conv.converged)
    if isequal(rel.out.conv.converged,1)
        lines{end+1,1} = 'Converged: Yes';
    else
        lines{end+1,1} = 'Converged: No (see estimation diagnostics)';
    end
else
    lines{end+1,1} = 'Converged: not assessed';
end

% CmdStan divergent transitions (report-only; NaN/absent -> not assessed, e.g.
% native engines that expose no divergence count).
nd = [];
if local_has(rel,{'out','conv','ndivergent'})
    nd = rel.out.conv.ndivergent;
end
if isempty(nd) || (isnumeric(nd) && all(isnan(nd(:))))
    lines{end+1,1} = 'Divergent transitions: not assessed';
else
    lines{end+1,1} = sprintf('Divergent transitions: %g', nd);
end
end

%-------------------------------------------------------------------------
% Local helpers
%-------------------------------------------------------------------------
function rel = local_rel(psyrat_data)
if isstruct(psyrat_data) && isfield(psyrat_data,'rel') && isstruct(psyrat_data.rel)
    rel = psyrat_data.rel;
else
    rel = struct();
end
end

function tf = local_has(s,path)
%True when the nested field path exists on struct s.
tf = true;
for k = 1:numel(path)
    if isstruct(s) && isfield(s,path{k})
        s = s.(path{k});
    else
        tf = false;
        return;
    end
end
end

function lines = local_coef_lines(psyrat_data)
%Coefficient-estimand line for the two-facet (test-retest) designs.
%
%This helper takes psyrat_data rather than rel because reltype/gcoeff/n'_o live
%on relsummary, not rel -- they are view-time selections, not estimation
%settings. psyrat_trt_coeflabel returns '' for any struct that records no
%reltype, which is the gate: every one-facet, splits and dynamic export
%contributes no line and is unchanged.
lines = {};
lbl = psyrat_trt_coeflabel(psyrat_data);
if ~isempty(lbl)
    lines{end+1,1} = sprintf('Coefficient: %s', lbl);
end
end

function s = local_priors_line(rel)
if ~isfield(rel,'priors') || isempty(rel.priors)
    s = 'Priors: not recorded';
    return;
end
try
    if local_priors_are_default(rel.priors)
        s = 'Priors: PsyRAT defaults';
    else
        s = 'Priors: custom (see run-config sidecar for full specification)';
    end
catch
    s = 'Priors: recorded (see run-config sidecar for full specification)';
end
end

function tf = local_priors_are_default(p)
%A stored prior set counts as the defaults when every family it carries equals
%the current default for that family. A family ABSENT from the stored set is at
%its default by construction: before a family is exposed on
%psyrat_default_priors its builder hard-codes the very numbers the default
%block later records (priors.dod, 2026-09-05), so a result saved before that
%exposure must not print as custom when it is reopened. A family the defaults
%do not know, or any field that differs, is custom. Field order is ignored, as
%isequal ignores it, so a jsondecode'd replay compares the same way. A struct
%with no families passes vacuously and prints as the defaults, which is also
%what merging it over the defaults would fit; the engine never stores one
%(REL.priors is the merged full set), so that case is documentary only.
d = psyrat_default_priors;
tf = isstruct(p) && isscalar(p);
if ~tf
    return;
end
fams = fieldnames(p);
for i = 1:numel(fams)
    if ~isfield(d,fams{i}) || ~isequal(p.(fams{i}), d.(fams{i}))
        tf = false;
        return;
    end
end
end

function lines = local_sigtrl_lines(rel)
%Caveat line for the one-facet trial main effect. Scoped to the two one-facet
%non-difference designs ('ic' and its subject-level twin), because those are the
%analyses whose absolute-vs-relative distinction rests on sigma_i; test-retest,
%splits and dynamic designs partition error differently and would be misdescribed
%by this line. The DoD -> single-event route sets rel.analysis = 'ic', so it is
%covered here.
%
%State is read from rel.out directly rather than from the relsummary backfill,
%so the line does not depend on a summary having been built: an explicit
%sig_trl_source wins; otherwise a present, non-degenerate sig_trl means it was
%estimated, and an absent or all-zero one means the coefficient reduces to the
%pre-sigma_i behavior.
lines = {};
if ~isfield(rel,'analysis') || isempty(rel.analysis) || ...
        ~any(strcmpi(char(string(rel.analysis)),{'ic','ic_sserrvar'}))
    return;
end

if local_has(rel,{'out','sig_trl_source'}) && ~isempty(rel.out.sig_trl_source)
    src = strtrim(char(string(rel.out.sig_trl_source)));
elseif local_has(rel,{'out','sig_trl'})
    st = rel.out.sig_trl;
    if iscell(st)
        st = cell2mat(st);
    end
    if ~isempty(st) && isnumeric(st) && any(st(:) ~= 0)
        src = 'estimated';
    else
        src = 'legacy_default';
    end
else
    src = 'legacy_default';
end

if strcmpi(src,'estimated')
    return;   % nothing notable: sigma_i was modeled and enters absolute error
end

lines{end+1,1} = ['Trial main effect (sigma_i): NOT estimated for this result ' ...
    '- assumed 0. Absolute error therefore contains no trial term, so the ' ...
    'dependability and generalizability coefficients are numerically identical ' ...
    'here and the "Dependability" label does not reflect an estimated ' ...
    'absolute-error component.'];
end

function lines = local_gamma_dispersion_lines(rel)
%Emit the gamma dispersion-heterogeneity summary (person log-nu SD s_v, the
%mean<->dispersion correlation, and the implied marginal correction factor F),
%read from whichever slot the design populated: rel.out.gamma (one-facet, test-
%retest), rel.out.gamma_disp (difference / difference-of-differences, per-cell),
%or, for the DYNAMIC/conditional gamma designs (analyses 11/14), the per-stratum
%gro_sds/chol_corrmat cells (there is no single rel.out.gamma slot; s_v/rho_pv are
%pooled across strata for this header summary). Draw vectors are summarized by
%their median. Absent diagnostics -> no lines.
lines = {};
%LOCATION-SCALE gamma dynrel (gammascale = 2) first, as its own branch: it has no
%log-nu parameter at all, so none of the s_v / cov_pv / F machinery below applies
%and the quantities it would report do not exist. It is recognized by the
%Gaussian dynrel's slot names (pop_sdlog + b_sigma) appearing on a GAMMA run -
%safe only because the caller has already gated on family = 'gamma'; the same
%slots on a Gaussian run must never reach here.
%
%THIS BRANCH IS TESTED FIRST AND RETURNS, so it - not the absence of pop_lognu -
%is what keeps the log-nu chain below off a location-scale run. (Corrected
%2026-08-05: this note used to say the log-nu chain "would return [] for such a
%run anyway, since local_dynrel_gamma_disp requires out.pop_lognu". True of the
%field, but irrelevant to the outcome, because the chain is never reached while
%this branch fires. Verified by adding pop_lognu = [] to a real location-scale
%layout: the emitted lines are byte-identical.) The pop_lognu absence is a
%SECOND lock that matters only if this branch's four-field gate ever stops
%firing; a run that lost one of those four AND carried a pop_lognu placeholder
%would report "person log-nu SD, s_v" built from gro_sds column 2, which is s_sd.
lsg = local_dynrel_gamma_ls_disp(rel);
if ~isempty(lsg)
    lines{end+1,1} = ['Gamma scale parameterization: LOCATION-SCALE (log ' ...
        'residual SD). The residual is modeled directly as sigma(z) rather ' ...
        'than through the dispersion nu, so the dimension slope b_sigma is a ' ...
        'pure residual-SD effect. This is a DIFFERENT ESTIMAND from the ' ...
        'log-nu parameterization (gammascale=1); results from the two are not ' ...
        'comparable.'];
    lines{end+1,1} = sprintf(['Gamma dispersion (person log residual-SD SD, ' ...
        's_sd): median %.3f'], median(lsg.s_sd(:)));
    lines{end+1,1} = sprintf(['Gamma mean-scale correlation (rho_p,sigma): ' ...
        'median %.3f - estimated and reported, but it does NOT enter the ' ...
        'observed-scale conversion, because the residual is independent of ' ...
        'the mean under this parameterization'], median(lsg.rho_ps(:)));
    lines{end+1,1} = sprintf(['Marginal residual correction factor ' ...
        'F = exp(2*s_sd^2): median %.3f (F = 1 is the homogeneous-residual ' ...
        'plug-in)'], median(lsg.disp_factor(:)));
    %The counterpart of the is_dynrel line below, and it says the OPPOSITE
    %thing. Under log-nu every observed-scale term including the residual
    %carries exp(2*alpha(z)), so a mean-level slope cancels. Here the residual
    %is exp(2*log sigma(z) + 2*s_sd^2) and carries no alpha, so it does not.
    %
    %CORRECTED TWICE on 2026-08-05 (increment 8), both times by reading the
    %emitted text rather than by any gate.
    %
    %(a) SCOPE. This line used to enumerate "the between-person and trial
    %components", which are the only two the ONE-FACET design (analysis 11/26)
    %has. This branch now also serves the TWO-FACET dynamic designs (14/27),
    %whose signal side additionally carries occasion and the three two-way
    %interactions - six components, not two. Fixed by not enumerating: "every
    %signal component" is true of both, and stays true of case 12 when it lands.
    %
    %(b) PRECISION. It also said "but the residual does not [carry
    %exp(2*alpha(z))]", which overstates it. The residual is recovered by
    %subtraction as B*(...) + e_cond_var with B = exp(2*alpha + V), so part of it
    %DOES scale with amplitude; what has no counterpart on the signal side is the
    %amplitude-free e_cond_var term, and THAT is what breaks the cancellation.
    %The conclusion was right and the stated reason was not - the same defect
    %class increment 7 recorded. The calculator headers
    %(psyrat_rel_dynrel_gamma_ls, psyrat_rel_dynrel_trt_gamma_ls) already gave
    %the algebra correctly; only this user-facing line was loose.
    lines{end+1,1} = ['Gamma dynamic reliability (location-scale): every ' ...
        'signal component carries exp(2*alpha(z)), but the residual also ' ...
        'carries an amplitude-FREE term exp(2*log sigma(z) + 2*s_sd^2) that ' ...
        'the signal side has no counterpart for, so BOTH the mean slope (b) ' ...
        'and the residual-SD slope (b_sigma) move the conditional ' ...
        'coefficient. The mean-slope cancellation that holds under the log-nu ' ...
        'parameterization does not hold here; report both slope posteriors.'];
    %Subject-level pooled-residual disclosure, analysis 26 ONLY (S17-FIX). The
    %two-facet sibling 27 breaks person x trial out as its own component and has
    %no pooling counterpart, so the gate is the exact one-facet label.
    if local_gamma_is_analysis(rel,'ic_dynrel_sserrvar')
        lines{end+1,1} = ['Gamma subject-level residual (location-scale, ' ...
            'POOLED): each participant''s error is the person-by-trial term ' ...
            'of the log-linked mean surface at their own dimension value (a ' ...
            'POPULATION quantity - the same term the population curve ' ...
            'carries) plus their own conditional variance sigma_s(z)^2. The ' ...
            'participant''s OWN amplitude does not enter. The curve and the ' ...
            'per-participant points remain different estimands - the curve ' ...
            'marginalizes the conditional term over the person residual-SD ' ...
            'population and uses a reference trial count (the per-group ' ...
            'median unless overridden) - so points fall on BOTH sides of the ' ...
            'curve; that is not estimation error.'];
    end
    return;
end
%STATIC DIFFERENCE-OF-DIFFERENCES LOCATION-SCALE (analysis 9, gammascale = 2)
%first among the static LS branches. Its own four-column stamp means it cannot be
%confused with either branch below - each of those requires an exact column count
%(two and one respectively) and so already fails closed on a four-cell layout - but
%it is placed first so the reading order matches the fail-closed widths: 4, then 2,
%then 1.
%
%A SEPARATE branch rather than widening the two-event one, deliberately. That
%branch's text is written for exactly two EVENTS ("by event: median %.3f, %.3f");
%analysis 9's four cells are 2 events x 2 conditions combined through the
%c = [1 -1 -1 1] contrast, so calling them events would be wrong, and making the
%existing text variadic would mean editing a working, pinned branch to serve a
%layout it was not written for. One stamp per layout, each failing closed, is the
%convention increments 4 and 5 established.
lsdodg = local_static_gamma_dod_ls_disp(rel);
if ~isempty(lsdodg)
    lines{end+1,1} = ['Gamma scale parameterization: LOCATION-SCALE (log ' ...
        'residual SD). Each cell''s residual is modeled directly as sigma_c ' ...
        'rather than through the dispersion nu, so Var(Y | mu, sigma) = ' ...
        'sigma_c^2 exactly. This is a DIFFERENT ESTIMAND from the default ' ...
        'log-nu parameterization; results from the two are not comparable.'];
    %Per-cell, never pooled, for the same reason the two-event branch reports per
    %event: the cell contrasts are what the design exists to estimate, and the
    %difference-of-differences is built from all four of them at once.
    lines{end+1,1} = sprintf(['Gamma dispersion (person log residual-SD SD, ' ...
        's_sd) by cell: median %.3f, %.3f, %.3f, %.3f'], ...
        median(lsdodg.s_sd(:,1)), median(lsdodg.s_sd(:,2)), ...
        median(lsdodg.s_sd(:,3)), median(lsdodg.s_sd(:,4)));
    lines{end+1,1} = sprintf(['Gamma mean-scale correlation (rho_p,sigma) by ' ...
        'cell: median %.3f, %.3f, %.3f, %.3f - estimated and reported, but it ' ...
        'does NOT enter the observed-scale conversion, because the residual is ' ...
        'independent of the mean under this parameterization'], ...
        median(lsdodg.rho_ps(:,1)), median(lsdodg.rho_ps(:,2)), ...
        median(lsdodg.rho_ps(:,3)), median(lsdodg.rho_ps(:,4)));
    lines{end+1,1} = sprintf(['Marginal residual correction factor ' ...
        'F = exp(2*s_sd^2) by cell: median %.3f, %.3f, %.3f, %.3f (F = 1 is ' ...
        'the homogeneous-residual plug-in)'], ...
        median(lsdodg.disp_factor(:,1)), median(lsdodg.disp_factor(:,2)), ...
        median(lsdodg.disp_factor(:,3)), median(lsdodg.disp_factor(:,4)));
    %The scientific point for a DIFFERENCE-OF-DIFFERENCES design. The two-event
    %argument applies to every one of the six cell pairs, and the DoD contrast
    %combines all four cells, so there is no cell whose amplitude drops out.
    lines{end+1,1} = ['Gamma difference-of-differences (location-scale): ' ...
        'because log sigma_c = log mu_c + 0.5*log(2) - 0.5*log(nu_c), any ' ...
        'log-nu cell contrast equals 2*(b_j - b_k) - 2*(b_sigma,j - ' ...
        'b_sigma,k), so it moves when the cells differ in AMPLITUDE even at ' ...
        'identical residual precision. The difference-of-differences contrast ' ...
        '(default [1 -1 -1 1]) carries a nonzero weight on every cell, so no ' ...
        'cell''s amplitude cancels out of it. Under this parameterization the ' ...
        'b_sigma contrasts are pure absolute residual-SD contrasts.'];
    lines{end+1,1} = ['Gamma difference-of-differences residual ' ...
        '(non-concurrent): the residual covariance between cells is ' ...
        'structurally zero, not merely assumed - correlating the person scale ' ...
        'effects across cells correlates the residual VARIANCES, not the ' ...
        'residuals - so the contrast residual is the weighted sum of the four ' ...
        'per-cell residual variances. A concurrent difference-of-differences ' ...
        'is rejected outright for this design, in either parameterization.'];
    %The increment-6 finding, restated where a user of this design will meet it.
    %It is a property of how the coefficient and the contrast aggregate per-cell
    %error, not of any one run, so it is stated unconditionally rather than
    %gated on a diagnostic.
    lines{end+1,1} = ['Gamma per-cell residual-SD contrasts are noisier than ' ...
        'the coefficient itself: the coefficient depends on the SUM of the ' ...
        'per-cell residuals, where estimation errors partly cancel, while a ' ...
        'contrast depends on their DIFFERENCE, where the same errors compound. ' ...
        'A well-recovered generalizability or dependability estimate is ' ...
        'therefore NOT evidence that the per-cell absolute-precision ' ...
        'comparison is sharp; read the contrast posteriors directly.'];
    return;
end
%STATIC DIFFERENCE LOCATION-SCALE (analyses 7/8, gammascale = 2) next. It must sit
%ABOVE the subject-level LS branch below, because that one would otherwise absorb
%these designs and report them wrongly in three ways: it prints SUBJECT-LEVEL text
%("each participant's residual ... their coefficient") at what is, for analysis 7,
%a GROUP-LEVEL design; it flattens the two events with median(s_sd(:)), reporting
%a pooled median of two different events' scale SDs as though it were one number;
%and it returns, making every difference-score caveat further down unreachable.
%The subject-level branch is additionally tightened to require a single column, so
%it now fails closed rather than silently absorbing a two-event layout.
%
%Keyed on the explicit is_diff stamp (paired with is_ls), not on the column count.
%The stamp is the signature, matching the convention the LS branches already use.
lsdg = local_static_gamma_diff_ls_disp(rel);
if ~isempty(lsdg)
    lines{end+1,1} = ['Gamma scale parameterization: LOCATION-SCALE (log ' ...
        'residual SD). Each event''s residual is modeled directly as sigma_e ' ...
        'rather than through the dispersion nu, so Var(Y | mu, sigma) = ' ...
        'sigma_e^2 exactly. This is a DIFFERENT ESTIMAND from the default ' ...
        'log-nu parameterization; results from the two are not comparable.'];
    %Per-event, never pooled: the two events are the contrast this design exists
    %to estimate, so collapsing them to one median would hide exactly the
    %quantity the parameterization was chosen to separate.
    lines{end+1,1} = sprintf(['Gamma dispersion (person log residual-SD SD, ' ...
        's_sd) by event: median %.3f, %.3f'], ...
        median(lsdg.s_sd(:,1)), median(lsdg.s_sd(:,2)));
    lines{end+1,1} = sprintf(['Gamma mean-scale correlation (rho_p,sigma) by ' ...
        'event: median %.3f, %.3f - estimated and reported, but it does NOT ' ...
        'enter the observed-scale conversion, because the residual is ' ...
        'independent of the mean under this parameterization'], ...
        median(lsdg.rho_ps(:,1)), median(lsdg.rho_ps(:,2)));
    lines{end+1,1} = sprintf(['Marginal residual correction factor ' ...
        'F = exp(2*s_sd^2) by event: median %.3f, %.3f (F = 1 is the ' ...
        'homogeneous-residual plug-in)'], ...
        median(lsdg.disp_factor(:,1)), median(lsdg.disp_factor(:,2)));
    %The scientific point of this parameterization for a DIFFERENCE design, and it
    %is not the same as the subject-level one below. There the blended contrast is
    %a person contrast; here it is the EVENT contrast, which is the construct.
    lines{end+1,1} = ['Gamma difference score (location-scale): because ' ...
        'log sigma_e = log mu_e + 0.5*log(2) - 0.5*log(nu_e), a log-nu event ' ...
        'contrast (b_nu,2 - b_nu,1) equals 2*(b_2 - b_1) - 2*(b_sigma,2 - ' ...
        'b_sigma,1), so it moves when the two events differ in AMPLITUDE even ' ...
        'at identical residual precision and cannot be read as a precision ' ...
        'contrast. Under this parameterization b_sigma,2 - b_sigma,1 is a pure ' ...
        'absolute residual-SD contrast.'];
    lines{end+1,1} = ['Gamma difference residual (non-concurrent): the ' ...
        'residual covariance between the two events is structurally zero, not ' ...
        'merely assumed - correlating the person scale effects across events ' ...
        'correlates the residual VARIANCES, not the residuals - so the ' ...
        'difference residual is sigma^2_Delta = sigma^2_1 + sigma^2_2.'];
    if local_has(rel,{'out','gamma_disp'}) && isstruct(rel.out.gamma_disp) ...
            && isfield(rel.out.gamma_disp,'is_diff_sserr') ...
            && isequal(rel.out.gamma_disp.is_diff_sserr,true)
        lines{end+1,1} = ['Gamma subject-level difference reliability ' ...
            '(location-scale, POOLED): each participant''s per-event residual ' ...
            'is the person-by-trial term of the log-linked mean surface (a ' ...
            'POPULATION quantity - the same term the group residual carries) ' ...
            'plus their own conditional variance read off their own log ' ...
            'residual-SD effect. The participant''s OWN amplitude does NOT ' ...
            'enter their coefficient - unlike under the log-nu ' ...
            'parameterization - so screening participants on subject-level ' ...
            'dependability no longer partly ranks them by amplitude. The ' ...
            'conditional term carries no trial averaging factor, because ' ...
            'sigma has no facet terms; that is the main reason the two ' ...
            'parameterizations give different numbers on the same data.'];
    end
    return;
end
%STATIC SUBJECT-LEVEL LOCATION-SCALE (analyses 6/25, gammascale = 2) next, and it
%must sit above the log-nu chain for the same reason the dynrel LS branch does:
%there is no log-nu parameter, so none of the s_v / cov_pv / F machinery below
%applies. It is recognized by the explicit is_ls stamp the extractor sets rather
%than by slot shapes, because these designs populate pop_sdlog and gro_sds just
%like their log-nu twins and would otherwise be indistinguishable from them.
%
%It CANNOT be confused with the dynrel LS branch above, but NOT because a static
%design never creates out.b_sigma - analyses 7/8 do create that slot (as the
%observed-scale residual encoding, a different quantity from the dynrel dimension
%slope of the same name). The dynrel branch is safe because it requires FOUR
%fields: pop_sdlog, b_sigma, gro_sds AND chol_corrmat. Analyses 6/25 supply the
%first, third and fourth but no b_sigma; analyses 7/8 supply only b_sigma. Do not
%weaken that four-field gate - the dynrel text asserts "BOTH the mean slope and
%the residual-SD slope move the conditional coefficient", which is meaningless for
%any static design.
%
%Without this branch a location-scale case 6/25 run would emit NOTHING: the log-nu
%chain resolves g from fields named s_v/disp_factor, which the LS extractor
%deliberately does not create, so it would fall straight through the empty-g
%return. Silent is better than wrong, but neither is acceptable for an estimand
%this run was explicitly asked for.
lssg = local_static_gamma_ls_disp(rel);
if ~isempty(lssg)
    lines{end+1,1} = ['Gamma scale parameterization: LOCATION-SCALE (log ' ...
        'residual SD). Each participant''s CONDITIONAL (single-trial) ' ...
        'residual is modeled directly as sigma_p rather than through the ' ...
        'dispersion nu, so Var(Y | mu, sigma) = sigma_p^2 exactly. This is a ' ...
        'DIFFERENT ESTIMAND from the default ' ...
        'log-nu parameterization; results from the two are not comparable.'];
    lines{end+1,1} = sprintf(['Gamma dispersion (person log residual-SD SD, ' ...
        's_sd): median %.3f'], median(lssg.s_sd(:)));
    lines{end+1,1} = sprintf(['Gamma mean-scale correlation (rho_p,sigma): ' ...
        'median %.3f - estimated and reported, but it does NOT enter the ' ...
        'observed-scale conversion, because the residual is independent of ' ...
        'the mean under this parameterization'], median(lssg.rho_ps(:)));
    lines{end+1,1} = sprintf(['Marginal residual correction factor ' ...
        'F = exp(2*s_sd^2): median %.3f (F = 1 is the homogeneous-residual ' ...
        'plug-in)'], median(lssg.disp_factor(:)));
    %The counterpart of the log-nu subject-level caveat further down, and it says
    %the OPPOSITE thing. That one warns that the participant's own amplitude does
    %not cancel; here it does not enter at all. The residual sentence is split by
    %analysis: 6 (one-facet) pools the person x trial mean-surface term (S17-FIX),
    %25 (two-facet) reports person x trial as its own component and stays a pure
    %conditional read-out - there is no two-facet pooling counterpart.
    if strcmpi(char(string(rel.analysis)),'ic_sserrvar')
        lines{end+1,1} = ['Gamma subject-level reliability (location-scale, ' ...
            'POOLED): each participant''s residual is the person-by-trial ' ...
            'term of the log-linked mean surface (a POPULATION quantity - the ' ...
            'same term the group residual carries) plus their own conditional ' ...
            'variance read off their own log residual-SD effect. The ' ...
            'participant''s OWN amplitude does NOT enter their coefficient - ' ...
            'unlike under the log-nu parameterization - and the person prior ' ...
            'shrinks toward constant absolute SD rather than constant CV, ' ...
            'matching the Gaussian subject-level model. The conditional term ' ...
            'carries no trial averaging factor, because sigma has no facet ' ...
            'terms; that is the main reason the two parameterizations give ' ...
            'different numbers on the same data.'];
    else
        lines{end+1,1} = ['Gamma subject-level reliability (location-scale): each ' ...
            'participant''s residual is read off their own log residual-SD effect ' ...
            'and never touches the mean, so - unlike the log-nu parameterization - ' ...
            'the participant''s own amplitude does NOT enter their coefficient and ' ...
            'the person prior shrinks toward constant absolute SD rather than ' ...
            'constant CV, matching the Gaussian subject-level model. The residual ' ...
            'also carries no trial/occasion averaging factor, because sigma has no ' ...
            'facet terms; that is the main reason the two parameterizations give ' ...
            'different numbers on the same data.'];
    end
    return;
end
%SUBJECT-LEVEL DYNAMIC DIFFERENCE (analysis 13, location-scale only) first. It is
%the sibling of the case-12 branch immediately below and shares its reasons for
%existing: it must sit above the empty-g return, and it needs its OWN branch rather
%than a widened one, because local_gamma_is_analysis matches EXACTLY and case 13's
%label is 'ic_diff_dynrel_sserrvar', not 'ic_diff_dynrel'. Case 13's gamma init
%block creates sd_id/cor_id/cor_i/er_var_ss1, so the four-field dynrel LS gate
%(pop_sdlog/b_sigma/gro_sds/chol_corrmat) does not fire, and neither does the
%static subject-level difference gate, which keys on is_diff_sserr -- a stamp this
%path does not set. Without this branch a case-13 run emits NO dispersion caveat at
%all, which is the failure mode increment 5 recorded for case 9.
%
%There is no log-nu counterpart: analysis 13 is location-scale only.
if local_gamma_is_analysis(rel,'ic_diff_dynrel_sserrvar')
    lines{end+1,1} = ['Gamma scale parameterization: LOCATION-SCALE (log ' ...
        'residual SD). Each participant''s residual is modeled directly as ' ...
        'sigma_e,s(z) rather than through the dispersion nu, so the CONDITIONAL ' ...
        'variance is Var(Y | mu, sigma) = sigma_e,s(z)^2 exactly. This is a ' ...
        'DIFFERENT ESTIMAND from the log-nu parameterization (gammascale=1); ' ...
        'results from the two are not comparable. Within the Gamma family this design ' ...
        'is available only under this parameterization; the Gaussian family ' ...
        'supports the same design independently.'];
    lines{end+1,1} = ['Gamma subject-level dynamic difference: each ' ...
        'participant''s coefficient is evaluated at THEIR OWN dimension value ' ...
        'and conditioned on THEIR OWN residual-SD effects for the two events, ' ...
        'so two participants at the same dimension value can differ. The ' ...
        'between-person and trial components are POPULATION quantities ' ...
        'evaluated at that participant''s dimension value but NOT at their own ' ...
        'amplitude: a participant does not have their own between-person ' ...
        'variance.'];
    lines{end+1,1} = ['Gamma subject-level residual (location-scale, POOLED): ' ...
        'each participant''s per-event error is the person-by-trial term of ' ...
        'the log-linked mean surface (a population quantity at that ' ...
        'participant''s dimension value) plus their OWN conditional variance ' ...
        'sigma_e,s(z)^2 - restoring the combined residual Rast & Clayson (in ' ...
        'press, p. 4) define at the population level, instantiated per ' ...
        'participant the way the owner''s reference implementation, the ' ...
        'population curve and the concurrent designs all do. The conditional ' ...
        'term carries no ' ...
        'trial-averaging factor - sigma has no trial index to average over - ' ...
        'and the whole pooled quantity is a SINGLE-TRIAL error variance, ' ...
        'divided by each participant''s own per-event trial count.'];
    lines{end+1,1} = ['Gamma subject-level dynamic difference residual ' ...
        '(non-concurrent): the conditional residual covariance between the two ' ...
        'events is taken as zero, justified by the events being measured on ' ...
        'different physical trials - correlating the person scale effects ' ...
        'across events correlates the residual VARIANCES, not the residuals ' ...
        'themselves. The difference residual is therefore sigma^2_Delta,s(z) = ' ...
        'sigma^2_1,s(z) + sigma^2_2,s(z).'];
    lines{end+1,1} = ['Per-participant points and the population curve remain ' ...
        'DIFFERENT estimands, so the vertical gap between a point and the ' ...
        'curve is NOT estimation error and is NOT one-signed. Both carry the ' ...
        'same person-by-trial mean-surface term, but the curve marginalizes ' ...
        'the conditional variance over the population of person residual-SD ' ...
        'effects and divides by a reference per-event n'' (the median unless ' ...
        'overridden), while each participant uses their OWN residual-SD ' ...
        'effects and their OWN per-event trial counts. A participant whose ' ...
        'own residual-SD effect exceeds the population value, or whose trial ' ...
        'counts fall below the reference n'', can and routinely will sit ' ...
        'BELOW the curve. Do not read a point below the curve as a problem, ' ...
        'and do not read the curve as an average of the points.'];
    lines{end+1,1} = ['Gamma coefficient of variation (location-scale): CV is ' ...
        'DERIVED here, CV(Y | mu, sigma) = sigma(z)/mu(z), rather than being a ' ...
        'modeled constant as under log-nu (CV = sqrt(2/nu)). It therefore ' ...
        'varies across the mean surface and across z.'];
    return;
end
%DYNAMIC DIFFERENCE (analysis 12) next, and it MUST sit above the empty-g return
%below. Case 12 populates none of the three slots that return a non-empty g -- its
%gamma init block creates sd_id/cor_id/cor_i, not pop_lognu/gro_sds/chol_corrmat --
%so g is empty, the function returns, and a branch placed further down would be
%dead code. It therefore emitted NO dispersion caveat at all, which is the wrong
%way round: this is the one dynamic design where the exp(2*alpha(z)) cancellation
%demonstrably BREAKS (asymmetric event mean slopes move the surface with
%b_nu_dim = 0; see SCIENTIFIC_FORMULA_TRACEABILITY_TABLE.md, row GAMMA-DIFF-DYNREL),
%while analyses 11/14, where it provably holds, did get one.
%Keyed on an EXACT analysis match, never startsWith/contains: 'ic_diff_dynrel'
%is a strict prefix of ic_diff_dynrel_sserrvar / _trt / _rescor / _rho, so a prefix
%test would fire on cases 13/15/16 and every concurrent variant. Keying on
%REL.out.b_nu_dim instead would be wrong for the same reason in the other
%direction: gamma_chisquare_remaining_models.md specifies b_nu_dim for case 15 too.
if local_gamma_is_analysis(rel,'ic_diff_dynrel') ...
        && local_has(rel,{'out','b_sigma_dim'})
    %LOCATION-SCALE case 12 (gammascale = 2). This MUST be tested before the
    %log-nu branch below: both parameterizations tag REL.analysis =
    %'ic_diff_dynrel', so the exact-match test alone does not separate them and a
    %location-scale run would otherwise be handed the log-nu caveat, which talks
    %about b_nu_dim and mean-coupling - neither of which exists in this model.
    %
    %Keyed on b_sigma_dim, which within the gamma family is created ONLY by the
    %case-12 location-scale init branch. It cannot collide with the Gaussian
    %case-12 slot of the same name: the caller reaches this function only when
    %REL.family is 'gamma'. The four-field dynrel LS gate used elsewhere
    %(pop_sdlog/b_sigma/gro_sds/chol_corrmat) does NOT fire here, because case 12
    %stores sd_id/cor_id instead of gro_sds/chol_corrmat - hence a new branch
    %rather than a widened one.
    lines{end+1,1} = ['Gamma scale parameterization: LOCATION-SCALE (log ' ...
        'residual SD). Each event''s residual is modeled directly as ' ...
        'sigma_e(z) rather than through the dispersion nu, so the CONDITIONAL ' ...
        'variance is Var(Y | mu, sigma) = sigma_e(z)^2 exactly and no longer ' ...
        'moves with the amplitude. This is a DIFFERENT ESTIMAND from the ' ...
        'log-nu parameterization (gammascale=1); results from the two are not ' ...
        'comparable.'];
    lines{end+1,1} = ['Gamma residual component (location-scale): the REPORTED ' ...
        'residual sigma^2_pi,e(z) CAN vary with z even when b_sigma_dim = 0. ' ...
        'It confounds the person x trial interaction of the log-linked mean ' ...
        'surface, which scales as exp(2*alpha(z)), with the conditional ' ...
        'variance, which does not - so a nonzero mean slope b_dim moves it ' ...
        'whenever both the person and trial log-SDs are nonzero. What this ' ...
        'parameterization decouples from the mean is the conditional part: the ' ...
        'piece representing measurement error rather than mean-surface curvature.'];
    lines{end+1,1} = ['Gamma dynamic difference (location-scale): the ' ...
        'conditional coefficient can move through (a) the absolute ' ...
        'residual-SD dimension slopes b_sigma_dim, (b) asymmetry between the ' ...
        'events'' log-mean slopes b_dim, and (c) the resulting change in the ' ...
        'observed-scale signal contrast and cross-event covariance. A ' ...
        'mean-only asymmetry moves the surface even at b_sigma_dim = 0, ' ...
        'exactly as under log-nu. What DOES differ: b_sigma_dim is a pure ' ...
        'residual-SD slope, so a dimension effect carried by it IS ' ...
        'attributable to precision - unlike b_nu_dim, which blends amplitude ' ...
        'with precision. Report b_dim and b_sigma_dim together.'];
    lines{end+1,1} = ['Gamma dynamic difference residual (non-concurrent): ' ...
        'the conditional residual covariance between the two events is ' ...
        'structurally zero at every z - correlating the person scale effects ' ...
        'across events correlates the residual VARIANCES, not the residuals - ' ...
        'so the difference residual is sigma^2_Delta(z) = sigma^2_1(z) + ' ...
        'sigma^2_2(z). The separate model-implied nonlinear interaction ' ...
        'cross-covariance induced by the log link is omitted; it is a ' ...
        'documented approximation, measured at 4e-06 to 6e-05 of the residual ' ...
        'variance it sits beside.'];
    lines{end+1,1} = ['Gamma coefficient of variation (location-scale): CV is ' ...
        'DERIVED here, CV(Y | mu, sigma) = sigma(z)/mu(z), rather than being a ' ...
        'modeled constant as under log-nu (CV = sqrt(2/nu)). It therefore ' ...
        'varies across the mean surface and across z.'];
    return;
end
if local_gamma_is_analysis(rel,'ic_diff_dynrel')
    lines{end+1,1} = ['Gamma dynamic difference: the observed-score residual is ' ...
        'mean-coupled (Var = 2*mu(z)^2/nu(z) per event) and the two events scale ' ...
        'ASYMMETRICALLY, so the exp(2*alpha(z)) cancellation that holds for the ' ...
        'single-event dynamic designs does NOT hold here. The conditional ' ...
        'coefficient can move through (a) the log-nu dimension slope b_nu_dim, ' ...
        '(b) asymmetry between the events'' log-mean slopes b_dim, and (c) the ' ...
        'resulting change in the observed-scale signal contrast and cross-event ' ...
        'covariance. A mean-only asymmetry moves the surface even at b_nu_dim = 0. ' ...
        'Report b_dim and b_nu_dim together; a dimension effect on reliability ' ...
        'here is NOT attributable to precision alone.'];
    lines{end+1,1} = local_gamma_cv_line();
    return;
end
if local_has(rel,{'out','gamma'}) && isstruct(rel.out.gamma) ...
        && isfield(rel.out.gamma,'s_v') && isfield(rel.out.gamma,'disp_factor')
    g = rel.out.gamma;
elseif local_has(rel,{'out','gamma_disp'}) && isstruct(rel.out.gamma_disp) ...
        && isfield(rel.out.gamma_disp,'s_v') && isfield(rel.out.gamma_disp,'disp_factor')
    g = rel.out.gamma_disp;
else
    g = local_dynrel_gamma_disp(rel);   % [] unless the dynrel gamma layout is present
end
if isempty(g)
    return;
end
%GLOBAL fixed dispersion (dispersion=2): the fitted model has ONE nu per event and no
%person log-nu parameter at all, so s_v / cov_pv / disp_factor are STRUCTURAL
%constants the extractor supplies for layout compatibility, not posterior summaries.
%Reporting them through the estimated-dispersion lines below would state "person
%log-nu SD: median 0.000" and "F: median 1.000", which reads as an empirical finding
%of homogeneous dispersion. Say what was actually fit instead.
if isfield(g,'is_global_nu') && g.is_global_nu
    lines{end+1,1} = ['Gamma dispersion: GLOBAL (fixed) nu per event - one dispersion ' ...
        'estimated per event and held constant across persons. There is no person ' ...
        'log-nu parameter in this model, so no s_v / mean-dispersion correlation is ' ...
        'estimated and the marginal residual correction is F = 1 by construction ' ...
        '(the population-nu plug-in), not as an empirical result.'];
else
    sv = median(g.s_v(:));
    F  = median(g.disp_factor(:));
    lines{end+1,1} = sprintf('Gamma dispersion (person log-nu SD, s_v): median %.3f', sv);
    if isfield(g,'rho_pv') && ~isempty(g.rho_pv)
        lines{end+1,1} = sprintf('Gamma mean-dispersion correlation: median %.3f', ...
            median(g.rho_pv(:)));
    end
    lines{end+1,1} = sprintf(['Marginal residual correction factor ' ...
        'F = exp(0.5*s_v^2 - 2*cov): median %.3f (F = 1 is the population-nu plug-in)'], F);
end
%is_dynrel is set only by local_dynrel_gamma_disp below, i.e. the NON-difference
%dynamic designs (analyses 11 and 14), which carry a single event. This line used
%to tell users a dimension effect "blends mean-level and dispersion effects".
%That is wrong for exactly these designs: under the log link every observed-scale
%component carries the same exp(2*alpha(z)) factor (psyrat_gamma_varcomps :96,
%:99, :114, :123), so a mean-level slope scales numerator and denominator alike
%and cancels in the ratio. Checked numerically against the production surface
%calculators: at b_nu = 0 a mean slope leaves the one-facet G/D/ICC invariant to
%~1e-15 and the two-facet CES to ~1e-15, while the residual itself moves by
%nearly a factor of two. Only the DISPERSION slope moves the coefficient here,
%because nu(z) enters the residual alone. Do not restore the old wording; the
%blend claim IS correct for the two-event difference designs below, where the
%events scale asymmetrically and the cancellation breaks.
if isfield(g,'is_dynrel') && g.is_dynrel
    lines{end+1,1} = ['Gamma dynamic reliability: every observed-score component ' ...
        'varies with z and the residual is mean-coupled (Var = 2*mu(z)^2/nu(z)), ' ...
        'but they share a common exp(2*alpha(z)) factor, so a mean-level dimension ' ...
        'effect scales them together and cancels in the coefficient. The conditional ' ...
        'coefficient moves with the RELATIVE-dispersion (log-nu) slope; report the ' ...
        'b_nu / s_v posteriors (the dispersion side is weakly identified at low ' ...
        'trial counts).'];
    %Ruling 8 (2026-08-04): b_nu is a RELATIVE-dispersion slope, and saying only
    %"the dispersion slope moves the coefficient" invites the absolute-precision
    %reading this parameterization does not support. Give the sign and the factor.
    lines{end+1,1} = ['Gamma dynamic reliability, interpreting b_nu: it is a slope ' ...
        'on log-nu, i.e. on RELATIVE dispersion, so the implied log-CV slope is ' ...
        '-b_nu/2 and a POSITIVE b_nu means DECREASING relative error. It is not a ' ...
        'statement about absolute residual SD: the log residual-SD slope is ' ...
        'b_sigma = b - b_nu/2, which includes the log-mean slope b, so -b_nu/2 and ' ...
        'b_sigma differ by exactly b. For an absolute-precision estimand refit with ' ...
        '''gammascale'',2.'];
    %Subject-level pooled-residual disclosure, analysis 26 ONLY (S20, the
    %log-nu counterpart of the S17-FIX line in the location-scale branch). The
    %two-facet sibling 27 breaks person x trial out as its own component and
    %has no pooling counterpart, so the gate is the exact one-facet label.
    if local_gamma_is_analysis(rel,'ic_dynrel_sserrvar')
        lines{end+1,1} = ['Gamma subject-level residual (log-nu, POOLED): each ' ...
            'participant''s error is the person-by-trial term of the log-linked ' ...
            'mean surface at their own dimension value (a POPULATION quantity - ' ...
            'the same term the population curve carries) plus their own ' ...
            'conditional variance 2*mu_s(z)^2/nu_s(z), evaluated at their own ' ...
            'log-mean and log-nu effects. The participant''s own amplitude enters ' ...
            'through the conditional part only; the pooled term does not depend ' ...
            'on the participant. The curve and the per-participant points remain ' ...
            'different estimands, so points fall on both sides of the curve.'];
    end
elseif isfield(g,'is_diff_sserr') && g.is_diff_sserr
    lines{end+1,1} = ['Gamma subject-level difference: the observed-score residual is ' ...
        'mean-coupled (Var = 2*mu^2/nu per event, per subject), so an event/condition ' ...
        'effect on the difference blends mean-level and dispersion (log-nu) effects; ' ...
        'report the s_v posterior (the dispersion side is weakly identified at low ' ...
        'trial counts).'];
    lines{end+1,1} = ['Gamma subject-level difference residual (log-nu, POOLED): ' ...
        'each participant''s per-event error is the per-event person-by-trial ' ...
        'term of the log-linked mean surface (a POPULATION quantity - the same ' ...
        'term the group residual carries) plus their own conditional variance ' ...
        '2*mu_e,s^2/nu_e,s at their own per-event log-mean and log-nu effects ' ...
        '(S20). The difference residual then follows the design''s ' ...
        'residual-covariance treatment: zero cross-covariance when the events ' ...
        'are non-concurrent, the Monte-Carlo copula cross-covariance when they ' ...
        'are concurrent. The pooled term''s cross-event half stays at zero by ' ...
        'convention (S14); the direction of that omission follows the sign of ' ...
        'the omitted cross-covariance (omitting a positive one lowers the ' ...
        'difference coefficient, omitting a negative one raises it).'];
%STATIC designs, added 2026-08-04 (S16). These populate g, so the dispersion
%numbers above were already printed for them, but neither is_dynrel nor
%is_diff_sserr is set, so they fell out of this chain with NO mean-coupling
%caveat -- and 6/25 are where the measured effect is largest.
%DELIBERATELY keyed on the analysis label rather than added as a blanket else:
%analyses 1-4 ('ic') and 5 ('trt') are SINGLE-MEAN designs where the
%exp(2*alpha) cancellation genuinely holds (a stratum coefficient is a function
%of s_p, s_i, nu and n only), so a caveat there would be false and alarming.
%The split below is not cosmetic -- the two families fail differently. Static
%DIFFERENCES break through the event/cell mean CONTRAST; SUBJECT-LEVEL designs
%break through the participant's own log-mean, which the population numerator
%does not carry.
elseif local_gamma_is_analysis(rel,{'ic_diff','trt_diff','ic_dodiff'})
    lines{end+1,1} = ['Gamma difference score: the observed-score residual is ' ...
        'mean-coupled (Var = 2*mu^2/nu per event/cell), so the event with the ' ...
        'larger mean is FORCED to carry proportionally larger absolute error ' ...
        'variance. A shift common to every event cancels out of the coefficient, ' ...
        'but the event CONTRAST does not -- so the reported difference-score ' ...
        'coefficient is partly determined by the condition effect being measured. ' ...
        'The per-event log-nu contrast describes RELATIVE dispersion, not absolute ' ...
        'precision. Report the per-event log-nu posteriors alongside the ' ...
        'coefficient.'];
    lines{end+1,1} = local_gamma_cv_line();
elseif local_gamma_is_analysis(rel,{'ic_sserrvar','trt_sserrvar'})
    lines{end+1,1} = ['Gamma subject-level reliability: each participant''s residual ' ...
        'is mean-coupled and evaluated at that participant''s OWN log-mean ' ...
        '(Var = 2*mu^2/nu), while the reliability numerator is a population ' ...
        'quantity. The grand mean therefore cancels but the participant''s own ' ...
        'amplitude does NOT: at equal log-nu, a higher-amplitude participant ' ...
        'receives a lower coefficient. The person log-nu prior shrinks toward ' ...
        'constant CV, where the Gaussian subject-level model shrinks toward ' ...
        'constant absolute SD. Report the s_v posterior, and treat participant ' ...
        'rankings and any dependability-based screening as partly amplitude-driven.'];
    %Pooled-residual disclosure, analysis 6 ONLY (S20): the two-facet
    %sibling 25 reports person x trial as its own component and pools nothing.
    if local_gamma_is_analysis(rel,'ic_sserrvar')
        lines{end+1,1} = ['Gamma subject-level residual (log-nu, POOLED): each ' ...
            'participant''s error is the person-by-trial term of the log-linked ' ...
            'mean surface (a POPULATION quantity - the same term the group ' ...
            'residual carries) plus their own conditional variance 2*mu_s^2/nu_s ' ...
            'at their own log-mean and log-nu effects. The amplitude statement ' ...
            'above concerns the conditional part; the pooled term does not ' ...
            'depend on the participant.'];
    end
    lines{end+1,1} = local_gamma_cv_line();
end
%The copula caveat is INDEPENDENT of the group-vs-subject-level split, so it must
%be its own test rather than a branch of the chain above: the CONCURRENT
%subject-level difference (analysis 8, diffrescor=2) sets BOTH is_diff_sserr and
%is_diff_copula (psyrat_gamma_extract_diff_copula_sserr inherits gamma_disp from
%the group copula extractor and then adds is_diff_sserr), so as an elseif this
%line was unreachable for the one design whose residual really is MC-projected.
if isfield(g,'is_diff_copula') && g.is_diff_copula
    copula_line = ['Gamma concurrent difference (Gaussian copula): the ' ...
        'observed-score residuals of the two concurrently measured events are coupled ' ...
        'by a Gaussian copula (residual correlation rho_e), and the log link makes ' ...
        'their observed residual cross-covariance nonzero even at rho_e = 0 (a ' ...
        'mean-surface interaction), so it is not rho_e*sigma1*sigma2 and has no closed ' ...
        'form - it is recovered by Monte-Carlo projection; report the residual ' ...
        'correlation (rescor) posterior'];
    %Under GLOBAL nu there is no s_v to report alongside rescor, so do not ask for it.
    if isfield(g,'is_global_nu') && g.is_global_nu
        copula_line = [copula_line '.'];
    else
        copula_line = [copula_line ' alongside s_v (the dispersion side is weakly ' ...
            'identified at low trial counts).'];
    end
    lines{end+1,1} = copula_line;
    %S11 experimental-status label. The two copula parameterizations are NOT
    %interchangeable in how well they sample, and nothing else in an exported table
    %names which one produced the number, so key the label off the RESOLVED dispersion
    %the estimator stamped onto rel (psyrat_computevarcomp REL.dispersion, the same
    %value the run-config sidecar records as cfg.design.dispersion). It MUST stay
    %inside this is_diff_copula branch: REL.dispersion is stamped on EVERY run,
    %defaulting to 1 outside the group-level concurrent gamma difference, so keying on
    %it alone would put an "experimental" label on every Gaussian run and on every
    %gamma design where the option is meaningless. Same emit-only-when-applicable
    %gating as the reltype-gated coefficient line above.
    lines{end+1,1} = local_copula_dispersion_line(rel, g);
end
end

function lines = local_modular_cut_lines(rel)
%Disclosure lines for the MODULAR CUT-COPULA designs (today analyses 19/20,
%family = 'gamma'; keyed on the stage output itself so future designs built on
%the stage inherit the disclosure): what kind of posterior produced the
%numbers, the location-scale estimand statement, and the stage-2 diagnostic
%flag counts. Three lines when a stage output is present, empty otherwise, per
%the emit-only-when-applicable convention every other branch in this file
%follows (existing output stays byte-identical for every other run).
%
%Why this discloses at all: the two-stage CUT posterior fits the gamma margins
%with no copula term, then estimates the copula correlation rho_e conditionally
%on each retained marginal draw. There is no feedback from stage 2 to stage 1,
%so the reported intervals are NOT full joint-posterior credible intervals, and
%rho_e is known to attenuate at moderate designs. A methods section written
%from output that omits this would misdescribe the model.
%
%The framework token is read from the archived stage labels, never restated
%here: psyrat_gamma_cut_copula_stage is the single source of truth for the
%strings, and a change there must flow through rather than silently diverge.
%Best-effort like every other branch: a stored result missing the fields emits
%what it can and never errors.
lines = cell(0,1);
%Keyed on the PRESENCE of out.copula - the field only the cut-copula stage arms
%write (verified: four sites, all in the 19/20 gamma arms) - rather than on
%exact analysis strings. The field is the exact invariant, and string-keying
%had a concrete future failure: a design later built on the same stage (the
%barred 21/22) would carry cut-posterior intervals with the disclosure
%silently absent (2026-08-11 code review). The exact-match doctrine elsewhere
%in this file guards against PREFIX matching the wrong design; a
%presence-of-output test cannot mismatch that way.
if ~local_has(rel,{'out','copula'}) || isempty(rel.out.copula) ...
        || ~iscell(rel.out.copula)
    return;
end
%The labels are run-constants written identically into every stratum's stage
%output; psyrat_copula_labels (shared with psyrat_dynrel_summary, so the two
%carriers cannot drift) takes the first stratum that carries them.
labels = psyrat_copula_labels(rel.out.copula);
if ~isempty(labels)
    framework = char(labels.inference_framework);
else
    %labels absent (pre-release stored result): still disclose the estimator
    %class, which is determined by the stage output's presence alone.
    framework = 'modular cut two-stage (gamma margins + Gaussian copula)';
end
lines{end+1,1} = sprintf(['Estimation framework: %s - a modular CUT '...
    'posterior: gamma margins fitted with no copula feedback, and the copula '...
    'residual correlation rho_e estimated conditionally on the retained '...
    'margin draws. Not a full joint Bayesian posterior; intervals do not '...
    'propagate copula-stage uncertainty into the margins.'], framework);
%The location-scale estimand statement every other gamma LS design (analyses
%12/13) emits; without it a reader comparing a 13 export with a 19 export
%cannot tell the two residual estimands have the same standing (2026-08-11
%code review). These designs are gammascale=2 by construction (the resolution
%forces it and rejects an explicit 1), so no gammascale check is needed here.
lines{end+1,1} = ['Gamma scale parameterization: LOCATION-SCALE (log ' ...
    'residual SD). Each participant''s residual is modeled directly as ' ...
    'sigma_e,s(z) rather than through the dispersion nu, so the CONDITIONAL ' ...
    'variance is Var(Y | mu, sigma) = sigma_e,s(z)^2 exactly. This is a ' ...
    'DIFFERENT ESTIMAND from the log-nu parameterization (gammascale=1); ' ...
    'results from the two are not comparable. Within the Gamma family this ' ...
    'design is available only under this parameterization; the Gaussian ' ...
    'family supports the same design independently.'];
%Aggregate the stage-2 diagnostic flags across strata via the shared counter.
%The affirmative "none flagged" is earned ONLY when every stratum was
%summarizable with both flag fields present: a skipped stratum or a missing
%flag field previously still printed the all-clear, rendering absence of
%evidence as evidence of absence (2026-08-11 code review). invalidity =
%numerical trouble (excessive PIT clipping or copula posterior mass at the rho
%support bound); calibration = descriptive PIT mean/SD departure. Neither
%blocks computation (diagnose-not-hide, per the stage contract).
ndraws = 0; ninv = 0; ncal = 0; counted = false; complete = true;
for c = 1:numel(rel.out.copula)
    f = psyrat_copula_flag_counts(rel.out.copula{c});
    counted = counted || f.recorded;
    complete = complete && f.recorded && f.flags_complete;
    ndraws = ndraws + f.n_draws;
    ninv = ninv + f.n_invalid;
    ncal = ncal + f.n_calibration;
end
if counted
    if ~complete
        lines{end+1,1} = sprintf(['Copula stage-2 diagnostics: not fully '...
            'recorded by this stored result - %d invalid and %d '...
            'calibration-flagged among the %d summarizable retained draws; '...
            'unrecorded strata or flag fields are not counted.'], ...
            ninv, ncal, ndraws);
    elseif ninv == 0 && ncal == 0
        lines{end+1,1} = sprintf(['Copula stage-2 diagnostics: none of %d '...
            'retained draws flagged (numerical invalidity or PIT '...
            'calibration).'], ndraws);
    else
        lines{end+1,1} = sprintf(['Copula stage-2 diagnostics: %d of %d '...
            'retained draws flagged numerically invalid, %d flagged for PIT '...
            'calibration; estimates are reported, not blocked - inspect the '...
            'summary copula_flags.'], ninv, ndraws, ncal);
    end
end
end

function lines = local_nonconcurrent_dod_lines(rel)
%Disclosure lines for the person-specific dynamic NONCONCURRENT
%difference-of-differences designs (analyses 28/29, gamma location-scale OR
%Gaussian): what kind of posterior produced the numbers, the counterfactual
%estimand statement, the scale-parameterization statement (family-specific
%wording), and where projection invalidity is recorded. Empty for every other
%run, per the emit-only-when-applicable convention (existing output stays
%byte-identical).
%
%Why this discloses at all: the estimand fixes all six pairwise
%observation-level residual covariances to exactly zero - the four
%constituents are TREATED as measured on distinct trials whether or not they
%were. A methods section written from output that omits this would present a
%counterfactual sensitivity estimand as the reliability of concurrent
%measurements.
%
%Keyed on the PRESENCE of out.dod_labels - the field only the 28/29 arms
%write - through the shared reader psyrat_dod_labels (also the summary
%carrier's reader, so the two cannot drift). Never keyed on analysis strings,
%for the same reason the modular-cut block is not. NOTE these designs have no
%out.copula, so the modular-cut block correctly stays silent for them; the
%two blocks are mutually exclusive by construction.
lines = cell(0,1);
labels = [];
if local_has(rel,{'out'})
    labels = psyrat_dod_labels(rel.out);
end
isdodstring = isfield(rel,'analysis') && ...
    (strcmp(rel.analysis,'ic_dodiff_dynrel_sserrvar') || ...
    strcmp(rel.analysis,'ic_dodiff_dynrel_sserrvar_trt'));
if isempty(labels)
    if isdodstring
        %DEGRADED MODE (labels absent on a damaged/pre-release stored
        %result): still disclose the estimand, keyed on the two EXACT
        %analysis strings (exact-match per this file's doctrine; the
        %presence-keyed path above remains primary). The copula block has
        %the same belt-and-braces arrangement.
        lines{end+1,1} = ['NONCONCURRENT COUNTERFACTUAL ESTIMAND '...
            '(disclosure labels missing from this stored result): the four '...
            'difference-of-differences constituents are treated as if '...
            'measured on distinct trials; all six pairwise conditional '...
            'residual covariances are fixed to exactly zero by assumption. '...
            'Must not be reported as concurrent reliability.'];
    end
    return;
end
%the family noun in the framework sentence follows the stamp; the quoted
%framework token (written by the engine's case arm per family) is the
%machine-readable authority and the noun must agree with it
famnoun = 'Gaussian';   %absence/empty == Gaussian, this file's doctrine
if isfield(rel,'family') && ...
        strcmpi(strtrim(char(string(rel.family))),'gamma')
    famnoun = 'Gamma';
end
lines{end+1,1} = sprintf(['Estimation framework: %s - the FULL Bayesian '...
    'posterior of four conditionally independent %s margins fitted '...
    'jointly (8-D person block, correlated facet blocks). No dependence '...
    'module exists; this is not a cut posterior.'], ...
    char(labels.inference_framework), famnoun);
%map the cell indices back to the user's event names: every per-cell output
%column says "cell q", and without this line the mapping lives only in the
%user's memory of the configured contrast (pre-merge review, 2026-08-13)
if isfield(rel,'dod_map') && numel(rel.dod_map) == 4
    lines{end+1,1} = sprintf(['DoD cells (constituent order): 1 = %s, '...
        '2 = %s, 3 = %s, 4 = %s; contrast (1 - 2) - (3 - 4).'],...
        char(string(rel.dod_map{1})),char(string(rel.dod_map{2})),...
        char(string(rel.dod_map{3})),char(string(rel.dod_map{4})));
end
lines{end+1,1} = ['NONCONCURRENT COUNTERFACTUAL ESTIMAND: the four '...
    'difference-of-differences constituents are treated as if measured on '...
    'distinct trials, so all six pairwise conditional residual covariances '...
    'are fixed to exactly zero by assumption. Person- and facet-level '...
    'covariances remain estimated (a different level of variation). These '...
    'results must not be reported as the reliability of physically '...
    'concurrent measurements without making this assumption explicit.'];
%The scale-parameterization statement, family-specific wording: the gamma arm
%contrasts against its log-nu sibling (these designs are gammascale=2 by
%construction; the resolution forces it and a dedicated guard rejects an
%explicit 1); the Gaussian family has no such sibling, so its line states the
%parameterization without the contrast.
if strcmp(famnoun,'Gamma')
    lines{end+1,1} = ['Gamma scale parameterization: LOCATION-SCALE (log ' ...
        'residual SD). Each participant''s residual is modeled directly as ' ...
        'sigma_e,s(z) rather than through the dispersion nu, so the CONDITIONAL ' ...
        'variance is Var(Y | mu, sigma) = sigma_e,s(z)^2 exactly. This is a ' ...
        'DIFFERENT ESTIMAND from the log-nu parameterization (gammascale=1); ' ...
        'results from the two are not comparable.'];
else
    lines{end+1,1} = ['Scale parameterization: person-specific log ' ...
        'residual SD. Each participant''s residual is modeled directly as ' ...
        'sigma_e,s(z), so the conditional variance is ' ...
        'Var(Y | mu, sigma) = sigma_e,s(z)^2 exactly.'];
end
%Projection invalidity (a materially negative signed composite, or a
%coefficient outside [0,1]) is a property of the STAGE-3 projection at the
%chosen D-study counts, not of the stored fit, so no count is claimed here -
%unlike the copula stage flags, which are archived per draw at estimation
%time. Point to where the counts live instead.
lines{end+1,1} = ['Signed-projection validity: composites are reported '...
    'unrepaired (record-never-repair); draws whose signed composite is '...
    'materially negative or whose coefficient leaves [0,1] are counted at '...
    'summary time in summary.strata(s).dod_flags and per participant in '...
    'the invalid_frac columns of the per-participant table.'];
end

function s = local_copula_dispersion_line(rel, g)
%Parameterization + experimental-status label for the Gaussian-copula concurrent
%gamma difference designs (audit finding S11). Called only from inside the
%is_diff_copula branch, so the design is already known to be one of that family;
%this function only names WHICH of the two structurally different models ran and
%what is known about its sampling behavior.
%
%Wording constraints, taken from S11 and its 2026-07-31 scope narrowing:
%  - label the PARAMETERIZATION, not the family. "The copula family is unvalidated"
%    is explicitly wrong now: global fixed dispersion is merged and is the
%    group-level default for analyses 7/10.
%  - do NOT say the design is wrong or the estimates invalid. The claim is about
%    sampling geometry and what has been checked here, not about the estimand.
%  - the R-hat ~2.5 / ESS ~2 figure is ONE live case-8 fit, so report it as one fit
%    and point the reader at this run's own diagnostics.
%  - do not let the global-nu evidence read as a convergence guarantee: S11 says its
%    live check "establishes tractability, not robustness" - one simulated dataset,
%    one settings choice, and the two-facet run still produced divergent transitions.
%
%An absent/empty rel.dispersion (a file saved before the field existed) yields the
%"not recorded" line rather than an error, per the best-effort contract in the header.
%A value other than 1 or 2 cannot reach here through the pipeline (psyrat_computevarcomp
%rejects it at input), so it falls into the same not-recorded branch.
d = [];
if isfield(rel,'dispersion') && ~isempty(rel.dispersion) && isnumeric(rel.dispersion)
    d = rel.dispersion(1);
end

if isequal(d, 2)
    s = ['Gamma copula dispersion parameterization: GLOBAL fixed nu per event ' ...
        '(dispersion = 2), the group-level default. There is no person log-nu ' ...
        'parameter, so the per-person funnel geometry is absent. This ' ...
        'parameterization reached max R-hat = 1.000 in the one in-repo check ' ...
        'performed (simulated known-truth data), but that is a single ' ...
        'dataset at a single choice of sampler settings and the two-facet check ' ...
        'still produced divergent transitions, so it establishes tractability, not ' ...
        'robustness, and is not a convergence guarantee for this dataset. Judge ' ...
        'this run on its own convergence verdict and divergent-transition count ' ...
        'and on the per-parameter R-hat / ESS in the estimation diagnostics.'];
elseif isequal(d, 1)
    s = ['Gamma copula dispersion parameterization: PER-PERSON log-nu ' ...
        '(dispersion = 1) - EXPERIMENTAL. The person-level log-nu term gives the ' ...
        'posterior a funnel geometry that is hard for NUTS, and convergence of ' ...
        'this parameterization has NOT been established in this repository: one ' ...
        'in-repo subject-level fit reached R-hat of about 2.5 with ESS of about 2, ' ...
        'and did not improve when the simulated participant and trial counts were ' ...
        'increased. That is one fit rather than a general property of the model, ' ...
        'so judge this run on its own convergence verdict and ' ...
        'divergent-transition count and on the per-parameter R-hat / ESS in the ' ...
        'estimation diagnostics before interpreting these estimates; long warmup ' ...
        'and/or threaded compute may be needed.'];
    if isfield(g,'is_diff_sserr') && g.is_diff_sserr
        s = [s ' Subject-level reliability requires per-person nu (global fixed ' ...
            'dispersion is a group-level estimand), so this design has no ' ...
            'non-experimental alternative parameterization.'];
    end
else
    s = ['Gamma copula dispersion parameterization: not recorded (saved before ' ...
        'this setting was stored), so it cannot be reported here whether the ' ...
        'EXPERIMENTAL per-person log-nu model or the global fixed-nu model ' ...
        'produced this result. Judge this run on its own convergence verdict and ' ...
        'divergent-transition count and on the per-parameter R-hat / ESS in the ' ...
        'estimation diagnostics.'];
end
end

function g = local_dynrel_gamma_disp(rel)
%Build the gamma dispersion summary from the DYNAMIC gamma per-stratum layout
%(rel.out.gro_sds{k} = [s_p s_v] draws; rel.out.chol_corrmat{k} = 2x2 Cholesky;
%rel.out.pop_lognu present). s_v / rho_pv / the implied F are pooled across strata
%for a single header summary. Returns [] when the layout is absent (so the caller
%falls through to no lines for non-dynrel gamma designs).
g = [];
if ~(local_has(rel,{'out','pop_lognu'}) && local_has(rel,{'out','gro_sds'}) ...
        && local_has(rel,{'out','chol_corrmat'}))
    return;
end
gro  = rel.out.gro_sds;
chol = rel.out.chol_corrmat;
if ~iscell(gro) || ~iscell(chol) || isempty(gro)
    return;
end
sv = []; rho = []; cov = [];
for k = 1:numel(gro)
    grok  = cell2mat(gro{k});    % draws x 2 [s_p, s_v]
    cholk = cell2mat(chol{k});   % draws x 2 x 2
    if size(grok,2) < 2
        continue;
    end
    s_p_k = grok(:,1);
    s_v_k = grok(:,2);
    rho_k = cholk(:,2,1);
    sv  = [sv;  s_v_k];             %#ok<AGROW>
    rho = [rho; rho_k];            %#ok<AGROW>
    cov = [cov; rho_k .* s_p_k .* s_v_k]; %#ok<AGROW>
end
if isempty(sv)
    return;
end
g = struct('s_v', sv, 'rho_pv', rho, ...
    'disp_factor', exp(0.5.*sv.^2 - 2.*cov), 'is_dynrel', true);
end

function g = local_dynrel_gamma_ls_disp(rel)
%Build the scale summary from the LOCATION-SCALE dynamic gamma layout
%(gammascale = 2): rel.out.gro_sds{k} = [s_p s_sd] draws, rel.out.chol_corrmat{k}
%= 2x2 Cholesky of the person (mean, log-sigma) correlation, plus pop_sdlog and
%b_sigma. s_sd / rho_p,sigma / the implied F are pooled across strata for a
%single header summary. Returns [] when the layout is absent.
%
%pop_sdlog ALONE is not a sufficient key - the Gaussian dynrel writes it too, as
%does every Gaussian subject-level design. What makes this branch safe is the
%caller's family = 'gamma' gate; b_sigma is required as well so that a future
%gamma design carrying a residual intercept but no dimension slopes cannot be
%mistaken for the dynamic one.
g = [];
if ~(local_has(rel,{'out','pop_sdlog'}) && local_has(rel,{'out','b_sigma'}) ...
        && local_has(rel,{'out','gro_sds'}) && local_has(rel,{'out','chol_corrmat'}))
    return;
end
gro  = rel.out.gro_sds;
chol = rel.out.chol_corrmat;
if ~iscell(gro) || ~iscell(chol) || isempty(gro)
    return;
end
ssd = []; rho = [];
for k = 1:numel(gro)
    grok  = cell2mat(gro{k});    % draws x 2 [s_p, s_sd]
    cholk = cell2mat(chol{k});   % draws x 2 x 2
    if size(grok,2) < 2
        continue;
    end
    ssd = [ssd; grok(:,2)];      %#ok<AGROW>
    rho = [rho; cholk(:,2,1)];   %#ok<AGROW>
end
if isempty(ssd)
    return;
end
%F is the factor by which marginalizing over the person residual-SD population
%inflates the residual relative to the homogeneous plug-in:
%E[sigma_p^2] = exp(2*log sigma) * exp(2*s_sd^2). It plays the same interpretive
%role as the log-nu F but is a different quantity, and is always >= 1 (the log-nu
%F can fall below 1 through the mean-dispersion covariance term, which has no
%counterpart here).
g = struct('s_sd', ssd, 'rho_ps', rho, 'disp_factor', exp(2.*ssd.^2));
end

function label = local_family_label(family)
%Human-readable label for the observation likelihood family. Falls back to the
%raw family name for anything other than the known 'gamma' key, so a future
%family surfaces its name rather than erroring. ASCII "chi-square" keeps the
%line encoding-safe for plain-text CSV headers.
switch lower(family)
    case 'gamma'
        label = ['Gamma / scaled chi-square '...
            '(observed-score-scale variance components)'];
    otherwise
        label = family;
end
end

function tf = local_gamma_is_analysis(rel,names)
%True when rel.analysis EXACTLY matches one of the supplied analysis labels.
%
%Exact match is load-bearing and must not be relaxed to startsWith/contains.
%The analysis labels nest as strict prefixes of one another -- 'ic_diff' prefixes
%'ic_diff_dynrel', 'ic_diff_sserrvar' and 'ic_diff_dynrel_sserrvar_trt_rescor';
%'ic_diff_dynrel' prefixes cases 13/15/16 and every concurrent variant -- so a
%prefix test would fire the wrong caveat on a design it does not describe.
%
%Callers are already gated on family = 'gamma' by psyrat_provenance_lines, so
%this does not re-check the family. Accepts a char label or a cellstr of labels.
tf = false;
if ~local_has(rel,{'analysis'})
    return;
end
a = rel.analysis;
if isstring(a)
    a = char(a);
end
if ~ischar(a) || isempty(a)
    return;
end
if ~iscell(names)
    names = {names};
end
for k = 1:numel(names)
    if strcmpi(a,names{k})
        tf = true;
        return;
    end
end
end

function s = local_gamma_cv_line()
%The constant-CV reading of nu, shared by every log-nu caveat.
%
%Var(Y | mu, nu) = 2*mu^2/nu makes the gamma family a constant
%coefficient-of-variation (proportional-error) model, so nu is the whole error
%model expressed in RELATIVE units rather than a nuisance shape parameter. This
%line exists so that no caveat says "the dispersion contrast" without saying
%which dispersion: users reading "dispersion" as absolute precision is exactly
%the misreading the location-scale option was added to prevent.
%
%Stated as a MODEL PROPERTY, not a defect. The log-nu models are coherent
%constant-CV models; what they do not support is an absolute-precision reading.
s = ['Gamma dispersion, how to read it: Var = 2*mu^2/nu makes this a constant ' ...
    'coefficient-of-variation model, CV = sqrt(2/nu) (= 1/3 at the default ' ...
    'nu = 18). A log-nu contrast is therefore a statement about RELATIVE ' ...
    'dispersion; since log sigma = log mu + 0.5*log(2) - 0.5*log(nu), every ' ...
    'log-nu contrast blends amplitude with absolute residual SD. This is a valid ' ...
    'model, not an error -- it simply does not support an absolute ' ...
    'residual-precision interpretation.'];
end

function g = local_static_gamma_dod_ls_disp(rel)
%Build the dispersion summary for a STATIC four-cell DIFFERENCE-OF-DIFFERENCES
%gamma fit under the location-scale parameterization (analysis 9,
%gammascale = 2). Returns [] for anything else, so the caller falls through to the
%two-event difference branch, then the subject-level LS branch, then the log-nu
%chain.
%
%Keyed on the is_dod stamp PAIRED with is_ls, both set by
%psyrat_gamma_extract_dodiff_ls. Keying on the column count instead would work
%today but would silently mis-route the moment another design gains four scale
%columns; the stamp is the signature, matching the convention the other LS
%branches use.
%
%Returns the per-cell MATRICES unflattened ([draws x 4]), like the two-event
%sibling and unlike the subject-level one. The four cells are what the
%c = [1 -1 -1 1] contrast is built from, so pooling them with (:) would hide
%exactly the quantity this parameterization was chosen to separate.
%
%Deliberately does NOT accept the log-nu field names, for the same fail-closed
%reason as the other two LS branches: an LS extractor that emitted s_v/cov_pv
%"for layout compatibility" would fall through to the log-nu chain, which would
%then print a person log-nu SD and an F = exp(0.5*s_v^2 - 2*cov) correction for a
%model that has no nu at all.
g = [];
if ~local_has(rel,{'out','gamma_disp'})
    return;
end
gd = rel.out.gamma_disp;
if ~isstruct(gd) || ~isfield(gd,'is_ls') || ~isequal(gd.is_ls,true)
    return;
end
if ~isfield(gd,'is_dod') || ~isequal(gd.is_dod,true)
    return;
end
if ~(isfield(gd,'s_sd') && isfield(gd,'rho_ps') && isfield(gd,'disp_factor'))
    return;
end
%Require the four-cell layout the caller indexes by column. A stamped-but-
%malformed struct falls through rather than erroring on gd.s_sd(:,4).
if size(gd.s_sd,2) ~= 4 || size(gd.rho_ps,2) ~= 4 || size(gd.disp_factor,2) ~= 4
    return;
end
g = struct('s_sd', gd.s_sd, 'rho_ps', gd.rho_ps, ...
    'disp_factor', gd.disp_factor);
end

function g = local_static_gamma_diff_ls_disp(rel)
%Build the dispersion summary for a STATIC two-event DIFFERENCE gamma fit under
%the location-scale parameterization (analyses 7 and 8, gammascale = 2). Returns
%[] for anything else, so the caller falls through to the subject-level LS branch
%and then to the log-nu chain.
%
%Keyed on the is_diff stamp PAIRED with is_ls, both set by
%psyrat_gamma_extract_diff_ls. Keying on the column count of s_sd instead would
%work today but would silently mis-route the moment a non-difference design gains
%a second scale column; the stamp is the signature, matching the convention the
%other LS branches use.
%
%Unlike its subject-level sibling this returns the per-event MATRICES unflattened
%([draws x 2]), because the two events are the contrast the design exists to
%estimate. The caller reports each event separately; pooling them with (:) would
%hide exactly the quantity this parameterization was chosen to separate.
%
%Deliberately does NOT accept the log-nu field names, for the same fail-closed
%reason as local_static_gamma_ls_disp: an LS extractor that emitted s_v/cov_pv
%"for layout compatibility" would fall through to the log-nu chain, which would
%then print a person log-nu SD and an F = exp(0.5*s_v^2 - 2*cov) correction for a
%model that has no nu at all.
g = [];
if ~local_has(rel,{'out','gamma_disp'})
    return;
end
gd = rel.out.gamma_disp;
if ~isstruct(gd) || ~isfield(gd,'is_ls') || ~isequal(gd.is_ls,true)
    return;
end
if ~isfield(gd,'is_diff') || ~isequal(gd.is_diff,true)
    return;
end
if ~(isfield(gd,'s_sd') && isfield(gd,'rho_ps') && isfield(gd,'disp_factor'))
    return;
end
%Require the two-event layout the caller indexes by column. A stamped-but-
%malformed struct falls through rather than erroring on gd.s_sd(:,2).
if size(gd.s_sd,2) ~= 2 || size(gd.rho_ps,2) ~= 2 || size(gd.disp_factor,2) ~= 2
    return;
end
g = struct('s_sd', gd.s_sd, 'rho_ps', gd.rho_ps, ...
    'disp_factor', gd.disp_factor);
end

function g = local_static_gamma_ls_disp(rel)
%Build the dispersion summary for a STATIC subject-level gamma fit under the
%location-scale parameterization (analyses 6 and 25, gammascale = 2). Returns []
%for anything else, so the caller falls through to the log-nu chain.
%
%Keyed on the explicit is_ls stamp the LS extractors set on REL.out.gamma_disp,
%NOT on slot shapes. These designs populate pop_sdlog / ind_sdlog / gro_sds
%exactly as their log-nu twins do - that is the whole point, since it lets
%psyrat_ssrel be reused unchanged - so there is no structural signature to key on.
%The stamp is the signature.
%
%Deliberately does NOT accept the log-nu field names. If a future edit makes the
%LS extractor emit s_v/rho_pv "for layout compatibility", this returns [] and the
%run falls through to the log-nu chain, which would then print a person log-nu SD
%and an F = exp(0.5*s_v^2 - 2*cov) correction for a model that has no nu at all.
%Failing closed here is the point.
g = [];
if ~local_has(rel,{'out','gamma_disp'})
    return;
end
gd = rel.out.gamma_disp;
if ~isstruct(gd) || ~isfield(gd,'is_ls') || ~isequal(gd.is_ls,true)
    return;
end
if ~(isfield(gd,'s_sd') && isfield(gd,'rho_ps') && isfield(gd,'disp_factor'))
    return;
end
%SINGLE-COLUMN ONLY. Added with the difference designs (analyses 7/8), whose
%gamma_disp carries one column PER EVENT and which are served by
%local_static_gamma_diff_ls_disp above. Without this check those designs would be
%absorbed here and reported wrongly: median(s_sd(:)) would pool two different
%events' scale SDs into one number, and the caller's subject-level text would be
%printed for the group-level analysis 7. Fail closed rather than summarize a
%layout this branch was not written for.
if size(gd.s_sd,2) ~= 1
    return;
end
g = struct('s_sd', gd.s_sd, 'rho_ps', gd.rho_ps, ...
    'disp_factor', gd.disp_factor);
end
