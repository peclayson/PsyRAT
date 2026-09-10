classdef TestDodPriorsExposed < PsyRATTestBase
    %Audit 2026-09-05 (owner ruling: expose before the beta). The Gaussian
    %four-event difference-of-differences model (analysis 9) was the one
    %design whose priors were hard-coded inside its Stan builder and its native
    %HMC twin, so a 'priors' override silently never reached it. The scales are
    %now priors.dod (b_cell, b_sigma_cell, sd_id, sd_trial, lkj) at the former
    %constants, threaded into psyrat_build_stan_dodiff and, through
    %native_extra, into local_hmc_dodiff. The Stan-text substitution itself is
    %pinned by TestStanModelSnapshots (dod_event at the defaults must stay
    %byte-identical to the pre-exposure golden; dod_event_priors pins a
    %non-default block). This file pins the pieces that test cannot see.
    %
    %RED against the pre-fix code: the defaults test (no dod family), the
    %provenance test (rmfield on a missing family errors, and a stored set
    %lacking the family compared unequal to the defaults), and both source
    %pins (the builder took one argument; the HMC arm read no priors).

    methods (Test)

        function testDefaultsCarryTheDodBlockAtTheFormerConstants(testCase)
            p = psyrat_default_priors();
            testCase.assertTrue(isfield(p, 'dod'), 'priors.dod must exist');
            testCase.verifyEqual(sort(fieldnames(p.dod)), ...
                sort({'b_cell'; 'b_sigma_cell'; 'sd_id'; 'sd_trial'; 'lkj'}));
            %the values the builder hard-coded before the exposure, so a
            %default run emits the same model
            testCase.verifyEqual(p.dod.b_cell, 10);
            testCase.verifyEqual(p.dod.b_sigma_cell, 2);
            testCase.verifyEqual(p.dod.sd_id, 2);
            testCase.verifyEqual(p.dod.sd_trial, 2);
            testCase.verifyEqual(p.dod.lkj, 2);
        end

        function testProvenanceTreatsAMissingFamilyAsTheDefaults(testCase)
            %A .psyrat saved before priors.dod existed stores rel.priors without
            %that family. Its header must still read "PsyRAT defaults": the
            %builder of that day used the very values the default block now
            %records. A changed value, or a family the defaults do not know,
            %must read as custom (the differential control).
            pd = localMinimalRel();
            pd.rel.priors = rmfield(psyrat_default_priors(), 'dod');
            testCase.verifyTrue(any(strcmp(psyrat_provenance_lines(pd), ...
                'Priors: PsyRAT defaults')), ...
                'a stored set lacking a later-exposed family is at its defaults');

            pd.rel.priors = psyrat_default_priors();
            testCase.verifyTrue(any(strcmp(psyrat_provenance_lines(pd), ...
                'Priors: PsyRAT defaults')), 'the full default set reads as defaults');

            pd.rel.priors = psyrat_default_priors();
            pd.rel.priors.dod.b_cell = 11;
            lines = psyrat_provenance_lines(pd);
            testCase.verifyTrue(any(startsWith(lines, 'Priors: custom')), ...
                'a changed dod scale must read as custom');
            testCase.verifyFalse(any(strcmp(lines, 'Priors: PsyRAT defaults')));

            pd.rel.priors = psyrat_default_priors();
            pd.rel.priors.notafamily = struct('x', 1);
            testCase.verifyTrue(any(startsWith(psyrat_provenance_lines(pd), 'Priors: custom')), ...
                'a family the defaults do not know must read as custom');
        end

        function testMergeFillsAPartialDodOverride(testCase)
            %The merge contract the docs promise (partial structs are safe)
            %holds for the new family: an override of one field keeps the
            %other four at their defaults. psyrat_merge_priors is a local
            %function of psyrat_computevarcomp, so the contract is exercised
            %through the emit hook: a partial dod override must produce the
            %substituted line AND the untouched default lines.
            stan = localEmitDodStan(testCase, struct('dod', struct('b_cell', 7)));
            testCase.verifyTrue(contains(stan, 'b_cell ~ normal(0, 7);'), ...
                'the overridden scale must reach the Stan text');
            testCase.verifyTrue(contains(stan, 'b_sigma_cell ~ student_t(10, 0, 2);'), ...
                'an omitted field keeps its default');
            testCase.verifyTrue(contains(stan, 'L_id ~ lkj_corr_cholesky(2);'), ...
                'the omitted LKJ shape keeps its default');
        end

        function testBothEnginesReadTheSameBlock(testCase)
            %Source-level pins for the native twin, which no headless test can
            %run to completion: the case-9 site must attach priors.dod to the
            %native spec, and the HMC log-posterior must read every one of the
            %five scales from the attached block rather than from a literal.
            src = fileread(which('psyrat_computevarcomp'));
            src = regexprep(src, '\.\.\.\s*[\r\n]+\s*', ' ');
            testCase.verifyTrue(contains(src, ...
                "native_extra = struct('priors',struct('dod',priors.dod));"), ...
                'case 9 must attach priors.dod to the native spec');
            testCase.verifyTrue(contains(src, 'psyrat_build_stan_dodiff(priors,true)') && ...
                contains(src, 'psyrat_build_stan_dodiff(priors,false)'), ...
                'both builder calls must pass the priors struct');
            hmc = fileread(which('psyrat_native_hmc'));
            testCase.verifyTrue(contains(hmc, "isfield(design.priors,'dod')"), ...
                'the HMC arm must read design.priors.dod');
            for fld = {'b_cell', 'b_sigma_cell', 'sd_id', 'sd_trial', 'lkj'}
                testCase.verifyTrue(contains(hmc, ['m.pri.' fld{1}]), ...
                    sprintf('the HMC log-posterior must consume m.pri.%s', fld{1}));
            end
        end

    end
end

function pd = localMinimalRel()
%The smallest struct psyrat_provenance_lines reads: every field is best-effort,
%so only .rel with a priors set is needed for the line under test.
rel = struct();
rel.seed = 12345;
rel.analysis = 'ic';
rel.engine = 'cmdstan';
pd = struct('rel', rel);
end

function stanText = localEmitDodStan(testCase, priors)
%Drive psyrat_computevarcomp for the four-event design with the emit hook
%active (the mechanism TestStanModelSnapshots uses) and return the generated
%Stan text. The run aborts before any CmdStan interaction.
events = {'E1', 'E2', 'E3', 'E4'};
ids = {}; evs = {}; meas = [];
for s = 1:4
    for e = 1:4
        for t = 1:6
            ids{end+1, 1} = sprintf('s%d', s); %#ok<AGROW>
            evs{end+1, 1} = events{e}; %#ok<AGROW>
            meas(end+1, 1) = 0.1 * t + 0.01 * e; %#ok<AGROW>
        end
    end
end
tbl = table(ids, meas, evs, 'VariableNames', {'id', 'meas', 'event'});

emitFile = [tempname '.stan'];
setenv('PSYRAT_EMIT_MODEL_FILE', emitFile);
restoreEmit = onCleanup(@() setenv('PSYRAT_EMIT_MODEL_FILE', '')); %#ok<NASGU>
cleanupFile = onCleanup(@() localDeleteIfExists(emitFile)); %#ok<NASGU>

runner = @() psyrat_computevarcomp('data', tbl, 'chains', 1, 'warmup', 1, ...
    'sampling', 1, 'seed', 12345, 'showgui', 0, 'verbose', 1, ...
    'sserrvar', 1, 'diffest', 3, 'diffwpcov', 1, 'diffrescor', 1, ...
    'dodmap', events, 'priors', priors);
testCase.verifyError(runner, 'PsyRAT:emitModelOnly');
testCase.assertTrue(isfile(emitFile), 'the emit hook must write the model file');
stanText = fileread(emitFile);
end

function localDeleteIfExists(f)
if isfile(f)
    delete(f);
end
end
