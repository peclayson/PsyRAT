function psyrat_data = psyrat_computevarcompwarp(varargin)
%Function to parse inputs from the toolbox into psyrat_computevarcomp
%
%psyrat_data = psyrat_computevarcompwarp('psyrat_prefs',psyrat_prefs,'psyrat_data',psyrat_data)
%
%
%Input
% psyrat_prefs - toolbox preferences
% psyrat_data - toolbox data
%
%Output
% psyrat_data - PsyRAT Toolbox dataset with the variance components and
%  convergence checks from Stan in psyrat_data.rel
% A directory will be created where the temporary files will be saved for
%  running the Stan model (stan creates various temp files). After the
%  model is done running (executed with psyrat_computevarcomp), the temporary
%  directory and files will be removed.

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

[psyrat_prefs, psyrat_data] = psyrat_findprefsdata(varargin);

%In case psyrat_computevarcompwarp was executed using the CLI, check to make sure
%the data are set up.

%check for psyrat_prefs
if isempty(psyrat_prefs)
    error('psyrat_prefs:notfound',...
        ['Input ''psyrat_prefs'' was not provided. ' ...
        'How to fix: call psyrat_computevarcompwarp with both ',...
        '''psyrat_prefs'' and ''psyrat_data'' arguments.']);
end

%check for psyrat_data
if isempty(psyrat_data)
    error('psyrat_data:notfound',...
        ['Input ''psyrat_data'' was not provided. ' ...
        'How to fix: call psyrat_computevarcompwarp with both ',...
        '''psyrat_prefs'' and ''psyrat_data'' arguments.']);
end
  
%check that data are loaded
if ~isfield(psyrat_data.proc,'data')
    %if the data have not been loaded where they are expected, try loading
    %the file specified in psyrat_data
    try
        psyrat_data.proc.data = psyrat_loadfile('psyrat_prefs',psyrat_prefs,...
        'psyrat_data',psyrat_data);
    catch
        error('psyrat_data:datanotfound',...
            ['Loaded dataset was not found in psyrat_data.proc.data and ',...
            'automatic load failed. How to fix: ensure input file metadata ',...
            'are set correctly, then run psyrat_loadfile or restart from ',...
            'psyrat_startproc.']);
    end
end

%reconcile the two alias pairs (diffrescor/diffwpcov and ssubjrel/sserrvar)
%in the shared canonicalizer; the policy tables and rationale live there (the
%two inline blocks that used to sit here were extracted verbatim, end states
%pinned by TestCanonicalizeProc's golden table against the pre-change code).
%This call site keeps the warp's historical warn-on-conflict behavior as the
%backstop for callers that build psyrat_prefs by hand: psyrat_run reconciles
%before calling here and the GUI writes both names of each pair from one
%control, so both entry points arrive consistent and this is a no-op for them.
psyrat_prefs.proc = psyrat_canonicalize_proc(psyrat_prefs.proc);
psyrat_prefs.proc = psyrat_sync_samplingprefs(psyrat_prefs.proc);

%check if there is a name of the dataset to be saved
if ~isfield(psyrat_data.proc,'savename')
    error('psyrat_data:savename',...
        strcat('Output filename is missing.\n\n',...
        'Set psyrat_data.proc.savename to a non-empty filename, such as ',...
        '''study.psyrat'', then rerun.'));
else
    %ensure something is in savename
    if isempty(psyrat_data.proc.savename)
        error('psyrat_data:savename',...
            strcat('Output filename is missing.\n\n',...
            'Set psyrat_data.proc.savename to a non-empty filename, such as ',...
            '''study.psyrat'', then rerun.'));
    end
    %check for whitespace for cmdstan
    if any(isspace(psyrat_data.proc.savename))
        error('psyrat_data:savename',...
            strcat('Output filename contains whitespace: ''',psyrat_data.proc.savename,'''\n\n',...
            'CmdStan does not accept whitespace in filenames.\n',...
            'Rename the file to remove spaces and rerun.'));
    end
end

%check if the path for the new dataset exists
if ~isfield(psyrat_data.proc,'savepath')
    error('psyrat_data:savepath',...
        strcat('Output save path is missing.\n\n',...
        'Set psyrat_data.proc.savepath to an existing directory path, then rerun.'));
else
    %check if the path is in savepath
    if isempty(psyrat_data.proc.savepath)
        error('psyrat_data:savepath',...
            strcat('Output save path is empty.\n\n',...
            'Set psyrat_data.proc.savepath to an existing directory path, then rerun.'));
    end
    %check if the path exists
    if ~exist(psyrat_data.proc.savepath,'dir')
        error('psyrat_data:savepath',...
            strcat('Output save path does not exist: ''',psyrat_data.proc.savepath,'''\n\n',...
            'Create this directory or choose a different existing folder, then rerun.'));
    end
    %check for whitespace for cmdstan
    if any(isspace(psyrat_data.proc.savepath))
        error('psyrat_data:savepath',...
            strcat('Output save path contains whitespace: ''',psyrat_data.proc.savepath,'''\n\n',...
            'CmdStan does not accept whitespace in paths.\n',...
            'Use a save path without spaces and rerun.'));
    end
end

%resolve the dynamic/conditional-reliability flag (Rast & Clayson; default off,
%1) before the preflight so the dynamic-reliability relaxations apply in the GUI
%path: the subject-level two-facet difference combinations (sserrvar = 2 with a
%time facet) are valid only under dynrel. The engine's own preflight
%(psyrat_computevarcomp) already passes dynrel; mirror it here. The while loop
%below re-reads dynrel defensively for the engine call.
if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'dynrel') && ...
        ~isempty(psyrat_prefs.proc.dynrel)
    dynrel = psyrat_prefs.proc.dynrel;
else
    dynrel = 1;
end

%resolve the reliability-with-splits flag (data splits; default off, 1) before
%the preflight so the splits checks (n_i column required, n_i-variation gate)
%apply in the GUI path, mirroring dynrel above. The while loop below re-reads it
%defensively for the engine call.
if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'splits') && ...
        ~isempty(psyrat_prefs.proc.splits)
    splits = psyrat_prefs.proc.splits;
else
    splits = 1;
end

%resolve the observation-likelihood family ('gaussian' default / 'gamma') before
%the preflight so the gamma positivity guard applies in the GUI path (a legacy
%prefs struct without the field therefore defaults to Gaussian). The family is
%forwarded to psyrat_computevarcomp for the engine call below.
if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'family') && ...
        ~isempty(psyrat_prefs.proc.family)
    family = lower(char(psyrat_prefs.proc.family));
else
    family = 'gaussian';
end

%run consolidated preflight checks so CLI and GUI workflows share logic
preflight = psyrat_preflight_validate(...
    'datatable', psyrat_data.proc.data, ...
    'sserrvar', psyrat_prefs.proc.ssubjrel, ...
    'diffest', psyrat_prefs.proc.diffest, ...
    'diffwpcov', psyrat_prefs.proc.diffwpcov, ...
    'diffrescor', psyrat_prefs.proc.diffrescor, ...
    'dynrel', dynrel, ...
    'splits', splits, ...
    'family', family, ...
    'savename', psyrat_data.proc.savename, ...
    'savepath', psyrat_data.proc.savepath, ...
    'requireSaveName', true, ...
    'requireSavePath', true);
if ~preflight.ok
    diag = preflight.errors(1);
    error(diag.id,'%s',diag.message);
end

%Non-fatal preflight warnings (few participants, a low median trial count,
%incomplete event x occasion cells; see psyrat_preflight_validate) were
%computed and dropped here before the pre-beta batch, because only errors(1)
%was read (council G8). Print them so the scripted route sees them before the
%fit; the GUI additionally raises them as a modal dialog before calling this
%function (psyrat_startproc). They never stop the run.
for wk = 1:numel(preflight.warnings)
    fprintf('\nPreflight warning (%s): %s\n', ...
        preflight.warnings(wk).id, preflight.warnings(wk).message);
end

%DoD runs require an explicit ERP_1..ERP_4 mapping.
dodmap = {};
if isfield(psyrat_prefs.proc,'diffest') && psyrat_prefs.proc.diffest == 3
    [dodmap,map_valid,map_msg] = psyrat_resolve_dodmap_pref(...
        psyrat_prefs.proc,psyrat_data.proc.data);
    if ~map_valid
        error('psyrat_prefs:dodmap',['DoD contrast mapping is invalid. ' map_msg]);
    end
end

%Change working dir for temporary Stan files
%change to vaoid problems with running multiple instances
tempdir = tempname(psyrat_data.proc.savepath);
[~,tempdir] = fileparts(tempdir);
tempdir = ['psyrat_stan_' tempdir];
mkdir(psyrat_data.proc.savepath,tempdir);

% if ~exist(fullfile(psyrat_data.proc.savepath,'Temp_StanFiles'), 'dir')
%   mkdir(psyrat_data.proc.savepath,'Temp_StanFiles');
% end
origdir = cd(fullfile(psyrat_data.proc.savepath,tempdir));

%Restore the working directory on EVERY exit path, including an error thrown by
%estimation below. The explicit cd() near the end of this function covers the
%normal path only; when psyrat_computevarcomp throws (a rejected option, a
%preflight failure, a CmdStan error) that line never runs and MATLAB is left
%INSIDE this temp directory. The directory is then removed, so the damage lands
%on the NEXT run that gets far enough to restore: it records the dangling path as
%its own origdir and dies in cd() with 'MATLAB:cd:NonExistentFolder' naming a
%directory from an unrelated earlier analysis. onCleanup fires on normal return
%and on error alike.
restoredir = onCleanup(@() psyrat_restore_dir(origdir)); %#ok<NASGU>

%Determine whether this run may show interactive GUI dialogs. Headless/batch
%runs (psyrat_data.proc.interactive == false, set by psyrat_run) must never
%block on a prompt: on non-convergence they auto-rerun with more iterations up
%to a fixed cap, then stop and let the post-loop logic report the
%non-convergence. Runs without the field default to interactive, so the GUI's
%existing rerun/trace-plot behavior is unchanged.
if isfield(psyrat_data.proc,'interactive') && ~isempty(psyrat_data.proc.interactive)
    interactive = logical(psyrat_data.proc.interactive);
else
    interactive = true;
end
maxreruns = 4;       %cap on automatic iteration-doublings for non-interactive runs
rerunattempt = 0;    %count of automatic reruns already performed

%set initial state of rerun to 1
%this will run psyrat_computevarcomp
%if chains properly converged then rerun will be changed to 0 and the while
%loop will be exited
rerun = 1;

%whether chains converged will be checked each time psyrat_computevarcomp is run
while rerun ~= 0
    
    if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'diffrescor')
        diffrescor = psyrat_prefs.proc.diffrescor;
    else
        diffrescor = 1;
    end

    if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'withinchain')
        withinchain = psyrat_prefs.proc.withinchain;
    else
        withinchain = 1;
    end

    if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'seed') && ...
            ~isempty(psyrat_prefs.proc.seed)
        seed = psyrat_prefs.proc.seed;
    else
        seed = 12345;
    end

    %optional CmdStan sampler overrides. Empty ([]) is the intended default
    %(keep each model's tuned adapt_delta / CmdStan's max_depth=10), so a missing
    %field maps to [] and estimation behavior is unchanged.
    if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'adapt_delta')
        adapt_delta = psyrat_prefs.proc.adapt_delta;
    else
        adapt_delta = [];
    end

    if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'max_treedepth')
        max_treedepth = psyrat_prefs.proc.max_treedepth;
    else
        max_treedepth = [];
    end

    %estimation engine and native-engine knobs. Default 0 (auto) for the
    %GUI/CLI path so CmdStan is used whenever present (unchanged behavior) and a
    %native fallback is used only when CmdStan is unavailable. A legacy prefs
    %struct without these fields therefore defaults to auto. See
    %psyrat_resolve_engine.
    if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'engine') && ...
            ~isempty(psyrat_prefs.proc.engine)
        engine = psyrat_prefs.proc.engine;
    else
        engine = 0;
    end
    if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'bootstrap_reps') && ...
            ~isempty(psyrat_prefs.proc.bootstrap_reps)
        bootstrap_reps = psyrat_prefs.proc.bootstrap_reps;
    else
        bootstrap_reps = 1000;
    end
    if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'hmc_gradient') && ...
            ~isempty(psyrat_prefs.proc.hmc_gradient)
        hmc_gradient = psyrat_prefs.proc.hmc_gradient;
    else
        hmc_gradient = 'auto';
    end

    %dynamic/conditional reliability (Rast & Clayson); default off (1).
    if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'dynrel') && ...
            ~isempty(psyrat_prefs.proc.dynrel)
        dynrel = psyrat_prefs.proc.dynrel;
    else
        dynrel = 1;
    end

    %reliability with data splits (1 = off/single-trial, 2 = parallel,
    %3 = nonparallel); default off (1).
    if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'splits') && ...
            ~isempty(psyrat_prefs.proc.splits)
        splits = psyrat_prefs.proc.splits;
    else
        splits = 1;
    end

    %residual-correlation parameterization for the subject-level concurrent
    %dynamic difference variants; default Population (1). Per-subject (2) selects
    %cases 21/22.
    if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'ssrescor') && ...
            ~isempty(psyrat_prefs.proc.ssrescor)
        ssrescor = psyrat_prefs.proc.ssrescor;
    else
        ssrescor = 1;
    end
    %Per-subject (2) is only valid for the subject-level concurrent dynamic
    %difference combination (dynrel=2, diffest=2, sserrvar=2, diffrescor=2); the
    %engine hard-errors otherwise. Clamp to Population (1) for any other run so a
    %stale Per-subject selection (the main-window popup keeps its value when
    %disabled, or no Dimension is selected) cannot trigger that error from the GUI
    %path. The engine combination guard remains the backstop for direct calls.
    if ssrescor == 2 && ~(dynrel == 2 && ...
            isfield(psyrat_prefs.proc,'diffest') && psyrat_prefs.proc.diffest == 2 && ...
            isfield(psyrat_prefs.proc,'ssubjrel') && psyrat_prefs.proc.ssubjrel == 2 && ...
            diffrescor == 2)
        ssrescor = 1;
    end

    %copula residual-dispersion parameterization for the GROUP-LEVEL concurrent gamma
    %difference designs. Forward the user's explicit choice VERBATIM; when unset, omit
    %it so the engine applies its context-aware default (GLOBAL for that design,
    %per-person elsewhere). The engine is the single authority on which designs accept
    %dispersion=2: it hard-errors ('varargin:dispersion') on an explicit 2 outside the
    %group-level design. Do NOT clamp/drop the value here - silently substituting a
    %different parameterization would fit a structurally different model than the user
    %asked for, and the sidecar written below would then record a value the run never
    %used. (No GUI control sets proc.dispersion, so a non-empty value is always a
    %deliberate CLI argument or a replayed config, never a stale popup selection.)
    dispersion_arg = {};
    if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'dispersion') && ...
            ~isempty(psyrat_prefs.proc.dispersion)
        dispersion_arg = {'dispersion', psyrat_prefs.proc.dispersion};
    end

    %gamma scale parameterization (1 = log-nu dispersion, 2 = log residual SD) for the
    %one-facet dynamic designs. Forwarded VERBATIM on the same terms as dispersion
    %above: the engine is the single authority on which designs accept gammascale=2 and
    %hard-errors ('varargin:gammascale') elsewhere, so clamping here would fit a
    %different estimand than the user asked for. Omitting it when unset lets the engine
    %apply its CONTEXT-AWARE default: since rulings 33-36 (2026-08-07) that is
    %location-scale on the designs that support it and log-nu everywhere else, NOT a
    %flat 1 as this comment previously said. Omitting is therefore the correct
    %behavior for an unset value, and the reason the GUI's "Automatic" option stores
    %[] rather than a concrete number - the engine is the only place that knows which
    %designs are movable. Without this forwarding the option would be unreachable from
    %the GUI/CLI route entirely.
    gammascale_arg = {};
    if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'gammascale') && ...
            ~isempty(psyrat_prefs.proc.gammascale)
        gammascale_arg = {'gammascale', psyrat_prefs.proc.gammascale};
    end

    %population-level priors. Forwarded on the same conditional terms as dispersion
    %and gammascale above, with one extra clause: proc.priors is NEVER empty on a
    %normal run (psyrat_defaults seeds it with the full psyrat_default_priors
    %struct), so emptiness alone cannot distinguish "user customized the priors"
    %from "untouched defaults". The isequal test makes a pure default run emit
    %exactly the argument list it emits today - the engine then applies
    %psyrat_default_priors itself - while a customized struct (set via
    %psyrat_run('priors',...) or a replayed sidecar) is merged onto the defaults
    %by psyrat_merge_priors inside the engine, so partial overrides stay safe.
    %isequal ignores struct field order, which keeps a jsondecode'd replay of a
    %default sidecar on the default path. The isstruct clause matters because the
    %engine performs no validation of 'priors' (a non-struct would be silently
    %dropped back to the defaults - the exact defect class this forwarding fixes).
    %Without this forwarding, priors set through psyrat_run or a replayed config
    %were recorded in the sidecar but never reached the engine (finding G33).
    priors_arg = {};
    if isfield(psyrat_prefs,'proc') && isfield(psyrat_prefs.proc,'priors') && ...
            ~isempty(psyrat_prefs.proc.priors) && ...
            isstruct(psyrat_prefs.proc.priors) && ...
            ~isequal(psyrat_prefs.proc.priors, psyrat_default_priors)
        priors_arg = {'priors', psyrat_prefs.proc.priors};
    end

    %pass the data to psyrat_computevarcomp for analysis
    REL = psyrat_computevarcomp('data',psyrat_data.proc.data,...
        'chains',psyrat_prefs.proc.nchains,...
        'warmup',psyrat_prefs.proc.nwarmup,...
        'sampling',psyrat_prefs.proc.nsampling,...
        'verbose',psyrat_prefs.proc.verbose,...
        'sserrvar',psyrat_prefs.proc.ssubjrel,...
        'diffest',psyrat_prefs.proc.diffest,...
        'diffrescor',diffrescor,...
        'withinchain',withinchain,...
        'diffwpcov',psyrat_prefs.proc.diffwpcov,...
        'dodmap',dodmap,...
        'dynrel',dynrel,...
        'ssrescor',ssrescor,...
        dispersion_arg{:},...
        gammascale_arg{:},...
        priors_arg{:},...
        'splits',splits,...
        'family',family,...
        'seed',seed,...
        'adapt_delta',adapt_delta,...
        'max_treedepth',max_treedepth,...
        'engine',engine,...
        'bootstrap_reps',bootstrap_reps,...
        'hmc_gradient',hmc_gradient,...
        'showgui',1);
    
    %check convergence of chains
    RELout = psyrat_checkconv(REL);

    %total the per-fit HMC divergent-transition counts psyrat_storeconv stamped
    %into each convergence table, and report them. This is independent of the
    %converged flag above: divergences are informational and do not trigger a
    %rerun, but they signal possible posterior bias even when R-hat looks fine.
    %NaN (native engines / unavailable / no conv table) prints nothing. The
    %viewer surfaces the same count via psyrat_convergence_notice.
    if isfield(RELout.out.conv,'data')
        RELout.out.conv.ndivergent = psyrat_total_divergences(RELout.out.conv.data);
    else
        RELout.out.conv.ndivergent = NaN;
    end
    if ~isnan(RELout.out.conv.ndivergent) && RELout.out.conv.ndivergent > 0
        fprintf(['\nWarning: %d divergent transition(s) detected. Variance ' ...
            'components and reliability estimates may be biased; consider more ' ...
            'warmup iterations or stronger priors.\n'], ...
            RELout.out.conv.ndivergent);
    end

    if RELout.out.conv.converged == 0

        if interactive
            %interactive run: ask the user whether to rerun with more iterations
            psyrat_reruncheck;
            psyrat_gui = findobj('Tag','psyrat_gui');
            if isempty(psyrat_gui)
                %the dialog is already gone (council G13): read that as
                %"do not rerun" rather than erroring after a finished fit
                rerun = 0;
            else
                rerun = guidata(psyrat_gui);
                close(psyrat_gui);
                delete(psyrat_gui);
            end
        else
            %headless/batch run: never prompt. Automatically rerun (iterations
            %are doubled below) up to maxreruns, then stop and let the post-loop
            %logic report the non-convergence.
            rerunattempt = rerunattempt + 1;
            if rerunattempt > maxreruns
                rerun = 0;
                fprintf(['\nModel did not converge after %d automatic reruns ' ...
                    '(non-interactive run); stopping.\n'],maxreruns);
            else
                rerun = 1;
                fprintf(['\nModel did not converge; non-interactive run, ' ...
                    'automatically rerunning (attempt %d of %d).\n'], ...
                    rerunattempt,maxreruns);
            end
        end

    else

        %if chains converged, do not rerun
        rerun = 0;
        fprintf('\nModel converged\n');

    end
    
    %check if the user wanted to see the trace plots (interactive runs only;
    %the trace-plot prompt is a GUI dialog, so headless runs skip it)
    if interactive && psyrat_prefs.proc.traceplots == 2 && rerun == 0 && ...
            strcmp(RELout.analysis,'ic')
        
        psyrat_checktraceplots(RELout,'askuser',1);
        
        %close gui after pulling input from user
        psyrat_gui = findobj('Tag','psyrat_gui');
        rerun = guidata(psyrat_gui);
        close(psyrat_gui);
        delete(psyrat_gui);
        
        %close window for trace plots
        psyrat_tplots = findobj('Tag','psyrat_tplots');
        if ~isempty(psyrat_tplots)
            close(psyrat_tplots);
            delete(psyrat_tplots);
        end
        
    elseif psyrat_prefs.proc.traceplots == 2 && rerun == 0 && ...
            (strcmp(RELout.analysis,'trt') || strcmp(RELout.analysis,'trt_sserrvar'))
        fprintf('\nViewing trace plots for trt analyses is not currently supported\n');
    elseif psyrat_prefs.proc.traceplots == 2 && rerun == 0 && ...
            strcmp(RELout.analysis,'ic_sserrvar')
        fprintf(strcat('\nViewing trace plots for subject-level',...
            ' reliability analyses is not currently supported\n'));
    elseif psyrat_prefs.proc.traceplots == 2 && rerun == 0 && ...
            (strcmp(RELout.analysis,'ic_diff') || ...
            strcmp(RELout.analysis,'ic_diff_sserrvar') || ...
            strcmp(RELout.analysis,'ic_dodiff'))
        fprintf(strcat('\nViewing trace plots for difference score',...
            ' reliability analyses is not currently supported\n'));
    elseif psyrat_prefs.proc.traceplots == 2 && rerun == 0 && ...
            contains(RELout.analysis,'dynrel')
        %dynamic/conditional reliability family (Rast & Clayson, cases 11-22).
        %Every dynamic-reliability analysis name contains 'dynrel'. Trace plots
        %are not rendered for these (as for the difference/subject-level
        %analyses); print an informative message so the run is never a silent
        %no-op. Convergence was already enforced by the rerun loop above.
        fprintf(strcat('\nTrace plots are not currently rendered for the',...
            ' dynamic/conditional reliability family.\nThe model converged',...
            ' (the run loop reruns with more iterations until it does).\n'));
    elseif psyrat_prefs.proc.traceplots == 2 && rerun == 0
        %any other analysis with trace plots requested: avoid a silent no-op.
        fprintf(strcat('\nTrace plots are not currently supported for this',...
            ' analysis type (',RELout.analysis,').\n'));
    end
    
    %if convergence was not met and the user would like to rerun the model,
    %double the number of iterations (if the doubled number is less than
    %1000, then run 1000 iterations)
    if rerun == 1
        psyrat_prefs.proc.nwarmup = psyrat_prefs.proc.nwarmup * 2;
        psyrat_prefs.proc.nsampling = psyrat_prefs.proc.nsampling * 2;
        psyrat_prefs.proc.niter = psyrat_prefs.proc.nwarmup + ...
            psyrat_prefs.proc.nsampling;
        if psyrat_prefs.proc.niter < 1000
            psyrat_prefs.proc.nwarmup = 500;
            psyrat_prefs.proc.nsampling = 500;
            psyrat_prefs.proc.niter = 1000;
        end
        fprintf('\nIncreasing warmup to %d and sampling to %d (total iterations=%d)\n',...
            psyrat_prefs.proc.nwarmup,psyrat_prefs.proc.nsampling,...
            psyrat_prefs.proc.niter);
    end
end

%sometimes the psyrat_gui doesn't close
psyrat_gui = findobj('Tag','psyrat_gui');
if ~isempty(psyrat_gui)
    close(psyrat_gui);
end
    
%change the working directory back to the original directory 
%(the onCleanup registered next to origdir also covers this, so the call here is
%the normal-path restore and the cleanup is the error-path safety net)
psyrat_restore_dir(origdir);

%attempt to remove the temporary directory and its files
try
    
    rmdir(fullfile(psyrat_data.proc.savepath,tempdir),'s');

catch
    
    fprintf('\n\nTemporary directory could not be removed.');
    fprintf('\nPath: %s\n',psyrat_data.proc.savepath);
end

%if the chains did not converge send the user back to psyrat_startproc_gui
if RELout.out.conv.converged == 0

    if isfield(psyrat_data.proc,'batchmode') && psyrat_data.proc.batchmode
        error('psyrat:modelDidNotConverge',...
            ['Model did not converge for measurement ''%s''.\n\n',...
            'How to fix: increase warmup/sampling iterations and rerun this measurement.'],...
            psyrat_data.proc.measheader);
    else
        fprintf('\nModel did not converge. Saving the current non-converged result because rerun was declined.\n');
    end
    
else
    
    %chains converged. Let the user know.
    %str = 'Models successfully converged with %d chains and %d iterations';
    %fprintf(strcat('\n',str,'\n'),procpref.nchains,procprefs.niter);
    
end

fprintf('\nSaving Processed Data...\n\n');

psyrat_data.rel = RELout;

save(fullfile(psyrat_data.proc.savepath,psyrat_data.proc.savename),'psyrat_data');

%write a human-readable run-configuration sidecar next to the .psyrat output so
%any run (GUI or CLI) records its seed and settings and can be replayed with
%psyrat_run('config',...). Best-effort: a sidecar failure must not discard the
%just-computed results.
%Pin the RESOLVED copula dispersion (the engine applies a context-aware default when
%the user does not pass one) so the sidecar records the estimand that was actually
%fit. Writing the unset [] instead would make a replay re-derive whatever default the
%running version happens to have, which is not the same guarantee.
if isfield(REL,'dispersion') && ~isempty(REL.dispersion)
    psyrat_prefs.proc.dispersion = REL.dispersion;
end
%Pin the RESOLVED gamma scale parameterization for the same reason. This one matters
%more than dispersion: gammascale=1 and 2 are different estimands whose results are not
%comparable, so a sidecar that omitted it would replay into whichever default the
%running version happens to have and silently produce different variance components.
if isfield(REL,'gammascale') && ~isempty(REL.gammascale)
    psyrat_prefs.proc.gammascale = REL.gammascale;
end
%Pin the priors the engine ACTUALLY FIT (REL.priors is the merged result of
%psyrat_merge_priors, or the defaults when nothing was forwarded) so the sidecar
%records what ran rather than the raw preference. Before this pin the sidecar
%could claim a custom prior set on a run that fit the defaults (finding G33).
%isfield-guarded because the headless test engine stubs return a REL without it.
if isfield(REL,'priors') && ~isempty(REL.priors)
    psyrat_prefs.proc.priors = REL.priors;
end
try
    psyrat_write_runconfig(psyrat_data,psyrat_prefs);
catch sidecarME
    fprintf('\nRun-config sidecar could not be written: %s\n',sidecarME.message);
end

end

function psyrat_restore_dir(origdir)
%Return to the directory the run started in. Tolerant on purpose: if origdir has
%been removed while the run was in progress, fall back to a directory that exists
%rather than throwing. This is called from a cleanup path, where an exception
%would either mask the real error or leave the session stranded in the very temp
%directory it was supposed to escape - which is the failure mode a bare cd() here
%produced ('MATLAB:cd:NonExistentFolder' naming an unrelated run's directory).
if exist(origdir,'dir') == 7
    cd(origdir);
else
    fprintf(['\n\nThe directory this run started in is no longer available, so '...
        'PsyRAT\nreturned to your home directory instead. The missing '...
        'directory was:\n  %s\n'],origdir);
    cd(psyrat_fallback_dir());
end
end

function d = psyrat_fallback_dir()
%A directory that is essentially always present, for the fallback above.
d = userpath;
if isempty(d) || exist(d,'dir') ~= 7
    if ispc
        d = getenv('USERPROFILE');
    else
        d = getenv('HOME');
    end
end
if isempty(d) || exist(d,'dir') ~= 7
    d = tempdir;
end
end

function proc = psyrat_sync_samplingprefs(proc)
%Ensure warmup/sampling preferences are available from legacy niter.

if ~isfield(proc,'niter') || isempty(proc.niter)
    proc.niter = 10000;
end

%optional sampler overrides default to empty ([]) when absent from legacy
%preference files; empty preserves current behavior (see psyrat_computevarcomp)
if ~isfield(proc,'adapt_delta')
    proc.adapt_delta = [];
end
if ~isfield(proc,'max_treedepth')
    proc.max_treedepth = [];
end

if ~isfield(proc,'nwarmup') || isempty(proc.nwarmup) || ...
        ~isfield(proc,'nsampling') || isempty(proc.nsampling)
    total_iters = max(round(proc.niter),2);
    proc.nwarmup = max(floor(total_iters/2),1);
    proc.nsampling = total_iters - proc.nwarmup;
end

proc.nwarmup = max(0,round(proc.nwarmup));
proc.nsampling = max(1,round(proc.nsampling));
proc.niter = proc.nwarmup + proc.nsampling;
end

function [dodmap,map_valid,msg] = psyrat_resolve_dodmap_pref(proc,datatable)
%Validate DoD mapping in processing preferences against loaded event labels.

dodmap = {};
map_valid = false;
msg = 'How to fix: set diffest=3 and configure "DoD Contrast" before running analysis.';

if ~istable(datatable) || ...
        ~any(strcmpi(datatable.Properties.VariableNames,'event'))
    msg = ['The processed data table does not include an event column. ',...
        'How to fix: define an event column and select exactly four events.'];
    return;
end

eventnames = cellstr(string(unique(datatable.event(:),'stable')));
if numel(eventnames) ~= 4
    msg = sprintf(['DoD requires exactly 4 selected events, but %d event type(s) were found. ',...
        'How to fix: select exactly four events for this run.'],numel(eventnames));
    return;
end

if ~isfield(proc,'dodset') || proc.dodset ~= 1
    msg = ['DoD contrast has not been configured for this run. ',...
        'How to fix: click "DoD Contrast..." and map ERP_1..ERP_4.'];
    return;
end

if ~isfield(proc,'dodmap') || numel(proc.dodmap) ~= 4
    msg = ['DoD mapping is missing or incomplete. ',...
        'How to fix: configure all four ERP mapping fields in the Configure DoD Contrast window.'];
    return;
end

if isfield(proc,'dodmapcol') && isfield(proc,'inp') && isfield(proc.inp,'event') && ...
        ~isempty(proc.dodmapcol) && ~isempty(proc.inp.event) && ...
        proc.dodmapcol ~= proc.inp.event
    msg = ['The event column changed after DoD contrast configuration. ',...
        'How to fix: click "DoD Contrast..." and save ERP_1..ERP_4 again.'];
    return;
end

dodmap = cell(1,4);
for i = 1:4
    dodmap{i} = char(string(proc.dodmap{i}));
end

if numel(unique(dodmap)) ~= 4
    dodmap = {};
    msg = ['DoD mapping must contain four unique events. ',...
        'How to fix: assign a different event to each ERP slot.'];
    return;
end

if ~all(ismember(dodmap,eventnames))
    dodmap = {};
    msg = ['DoD mapping does not match currently selected events. ',...
        'How to fix: click "DoD Contrast..." and remap ERP_1..ERP_4.'];
    return;
end

map_valid = true;
msg = '';
end
