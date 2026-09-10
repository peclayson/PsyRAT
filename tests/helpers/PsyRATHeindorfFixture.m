classdef PsyRATHeindorfFixture
    % External validation fixture from an independent (under-review) analysis.
    %
    % Source: Heindorf, G., Rocha, H. A., Chen, T., Vispoel, W. P., &
    %   Clayson, P. E. (under review). Reliability of outcome monitoring
    %   indices after effort exertion: A psychometric evaluation of the
    %   effort-doors task for the study of individual differences.
    %
    % The effort-doors task crosses effort (Low/High) with feedback
    % (Loss/Gain), giving four conditions. The published tables report, for
    % the reward positivity (RewP) and feedback-P3 (fb-P3) components:
    %   - between-person, between-trial, and residual (co)variances, and
    %   - the number of trials per condition required to reach reliability
    %     thresholds of 0.70/0.80/0.90/0.95, for each raw condition, for the
    %     two within-effort difference scores, and for the difference-of-
    %     differences (DoD) contrast.
    %
    % These published trial counts are an EXTERNAL reference: they were
    % produced by a separate analysis, not by this toolbox, so reproducing
    % them validates correctness rather than internal self-consistency.
    %
    % Important caveats encoded by the consuming tests:
    %   - The (co)variances above are rounded to two decimals in the source,
    %     so reproduction from these inputs is exact to within ~1 trial for
    %     single conditions.
    %   - The published counts derive from full posterior summaries, whereas
    %     feeding these point estimates as a single degenerate draw computes
    %     reliability at the posterior mean. For difference scores and the
    %     DoD contrast this produces counts a few percent below the published
    %     values (Jensen-type gap), so those tests use a relative tolerance.
    %
    % Condition order used throughout: [LE-Loss, LE-Gain, HE-Loss, HE-Gain]
    % (LE = low effort, HE = high effort).

    methods (Static)
        function c = components()
            thr = [0.70 0.80 0.90 0.95];

            % ----- Reward positivity (RewP) -----
            rewp.name = 'RewP';
            rewp.condition_order = {'LE-Loss','LE-Gain','HE-Loss','HE-Gain'};
            rewp.thresholds = thr;
            rewp.bp = [ ...
                5.30 5.53 4.49 4.99; ...
                5.53 6.90 4.65 5.92; ...
                4.49 4.65 4.88 4.54; ...
                4.99 5.92 4.54 6.16];
            rewp.bt = [ ...
                0.14 0.05  0.06 0.01; ...
                0.05 0.20  0.05 0.001; ...
                0.06 0.05  0.33 0.01; ...
                0.01 0.001 0.01 0.13];
            rewp.resid = [32.61 33.88 29.54 29.98];
            % published single-condition counts: rows = conditions, cols = thr
            rewp.single = [ ...
                15 25 56 118; ...
                12 20 45 94; ...
                15 25 56 117; ...
                12 20 44 93];
            % difference scores (balanced trials per constituent condition)
            rewp.diff_LEgain_LEloss = [140 239 537 1134];
            rewp.diff_HEgain_HEloss = [73 125 280 592];
            % difference-of-differences: (HE-Gain - HE-Loss) - (LE-Gain - LE-Loss)
            rewp.dod = [197 337 758 1599];

            % ----- Feedback-P3 (fb-P3) -----
            fbp3.name = 'fb-P3';
            fbp3.condition_order = {'LE-Loss','LE-Gain','HE-Loss','HE-Gain'};
            fbp3.thresholds = thr;
            fbp3.bp = [ ...
                7.06 6.88 6.87 6.63; ...
                6.88 7.52 7.12 6.98; ...
                6.87 7.12 7.59 6.91; ...
                6.63 6.98 6.91 7.39];
            fbp3.bt = [ ...
                0.60 0.31 0.36 0.40; ...
                0.31 0.70 0.33 0.37; ...
                0.36 0.33 0.53 0.22; ...
                0.40 0.37 0.22 0.70];
            fbp3.resid = [32.41 29.76 32.29 34.28];
            fbp3.single = [ ...
                11 19 43 90; ...
                10 17 37 78; ...
                11 18 39 83; ...
                12 19 43 90];
            fbp3.diff_LEgain_LEloss = [183 314 705 1488];
            fbp3.diff_HEgain_HEloss = [137 235 528 1114];
            fbp3.dod = [170 292 656 1385];

            c = [rewp fbp3];
        end

        function n = singleConditionCount(bp_var, bt_var, resid_var, threshold, ci)
            % Smallest trial count at which the point-estimate absolute
            % dependability reaches the threshold. Mirrors how psyrat_relsummary
            % derives the recommended trial cutoff: dependability across a trial
            % range, first crossing of the threshold. Within-person variance is
            % the trial main effect plus residual (Table 2 absolute error).
            bp = sqrt(bp_var);
            wp = sqrt(bt_var + resid_var);
            nmax = 4000;
            [~, pt, ~] = psyrat_dep('bp', bp, 'wp', wp, 'obs', [1 nmax], 'CI', ci);
            n = find(pt >= threshold, 1);
            if isempty(n)
                n = NaN;
            end
        end

        function n = diffScoreCount(bp2, bt2, erv2, threshold, ci)
            % bp2/bt2 are [1 x 2 x 2] (co)variance draws for the two
            % constituent conditions; erv2 is [1 x 2] log residual SD. Trials
            % are balanced across the two constituents. Non-concurrent
            % conditions, so the residual covariance is left at its default 0.
            fn = @(nn) localDiffRel(bp2, bt2, erv2, nn, ci);
            n = PsyRATHeindorfFixture.crossingBinary(fn, threshold, 4000);
        end

        function n = dodCount(bp4, bt4, erv4, threshold, ci)
            % bp4/bt4 are [1 x 4 x 4]; erv4 is [1 x 4] log residual SD.
            % Balanced trials across the four conditions; default contrast.
            fn = @(nn) localDodRel(bp4, bt4, erv4, nn, ci);
            n = PsyRATHeindorfFixture.crossingBinary(fn, threshold, 4000);
        end

        function n = crossingBinary(relFn, threshold, nmax)
            % Smallest integer n in [1, nmax] with relFn(n) >= threshold.
            % relFn is monotonically increasing in n (more trials reduce
            % measurement error), so a binary search is exact.
            if relFn(1) >= threshold
                n = 1;
                return;
            end
            if relFn(nmax) < threshold
                n = NaN;
                return;
            end
            lo = 1;
            hi = nmax;
            while hi - lo > 1
                mid = floor((lo + hi) / 2);
                if relFn(mid) >= threshold
                    hi = mid;
                else
                    lo = mid;
                end
            end
            n = hi;
        end
    end
end

function r = localDiffRel(bp2, bt2, erv2, nn, ci)
ds = psyrat_diffrel('bp', bp2, 'bt', bt2, 'er_var', erv2, ...
    'obs', [nn nn], 'est', 'dep', 'CI', ci);
r = ds.pt;
end

function r = localDodRel(bp4, bt4, erv4, nn, ci)
dd = psyrat_dodiffrel('bp', bp4, 'bt', bt4, 'er_var', erv4, ...
    'obs', [nn nn nn nn], 'est', 'dep', 'CI', ci);
r = dd.pt;
end
