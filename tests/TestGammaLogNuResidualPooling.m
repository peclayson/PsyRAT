classdef TestGammaLogNuResidualPooling < matlab.unittest.TestCase
    %TESTGAMMALOGNURESIDUALPOOLING The four log-nu subject-level gamma
    %residual sites must pool the population person x trial term (S20).
    %
    % S20 (PSYRAT_AUDIT_FINDINGS.md, owner-adjudicated 2026-08-12, fixed in
    % the pre-beta batch): under gammascale = 1 the subject-level designs 6,
    % 8 (both extractors) and 26 built each participant's residual as a single
    % exponential, E_t[Var(Y | mu, nu_p) | p], and dropped the induced
    % person x trial term of the lognormal mean surface that the group-level
    % residual carries and that S17-FIX pooled under location-scale. Every
    % per-participant coefficient was optimistically biased.
    %
    % Three of the four sites are local functions inside
    % psyrat_computevarcomp.m and cannot be called from a test (the repo's
    % standing limitation for gamma extractors), so this class pins the
    % SOURCE of each: the residual assignment must start from pi_mean and the
    % un-pooled single-exponential form must be gone. The location-scale
    % twin, already pooled by S17-FIX, is the positive control that the
    % predicate can pass. The analysis-26 calculator is a real function and
    % is pinned numerically in TestGammaDynrelSserrReliability; the live
    % recoveries (TestGammaSserrRecovery and siblings) close the loop.

    methods (Test)
        function testAnalysis6ExtractorPoolsPiMean(testCase)
            body = localFunctionBody(testCase, 'psyrat_gamma_extract_sserr');
            localAssertPooled(testCase, body, 'analysis 6, psyrat_gamma_extract_sserr');
        end

        function testAnalysis8ExtractorPoolsPiMean(testCase)
            body = localFunctionBody(testCase, 'psyrat_gamma_extract_diff_sserr');
            localAssertPooled(testCase, body, 'analysis 8, psyrat_gamma_extract_diff_sserr');
        end

        function testAnalysis8CopulaExtractorPoolsPiMean(testCase)
            body = localFunctionBody(testCase, 'psyrat_gamma_extract_diff_copula_sserr');
            localAssertPooled(testCase, body, 'analysis 8 concurrent, psyrat_gamma_extract_diff_copula_sserr');
        end

        function testAnalysis26CalculatorPoolsPiMean(testCase)
            src = localReadSource('subroutines/calculation/psyrat_ssrel_dynrel_gamma.m');
            localAssertPooled(testCase, src, 'analysis 26, psyrat_ssrel_dynrel_gamma');
        end

        function testLocationScaleTwinPassesThePredicatePositiveControl(testCase)
            %S17-FIX pooled this twin; the same predicate must accept it, so a
            %predicate that nothing could satisfy cannot go green vacuously.
            body = localFunctionBody(testCase, 'psyrat_gamma_extract_sserr_ls');
            testCase.verifyTrue(~isempty(regexp(body, ...
                'pi_mean\s*=\s*max\(vc\.var_mu\s*-\s*vc\.sigma_p2\s*-\s*vc\.sigma_i2,\s*0\)', 'once')), ...
                'positive control: the LS twin derives pi_mean by the helper subtraction');
            testCase.verifyTrue(~isempty(regexp(body, 'pi_mean\s*\+\s*exp\(', 'once')), ...
                'positive control: the LS twin adds pi_mean to its residual');
        end
    end
end

function src = localReadSource(relpath)
here = fileparts(mfilename('fullpath'));
src = fileread(fullfile(fileparts(here), relpath));
end

function body = localFunctionBody(testCase, fname)
%The text of local function FNAME in psyrat_computevarcomp.m, from its
%signature to the next top-level 'function' line.
src = localReadSource('subroutines/estimation/psyrat_computevarcomp.m');
starts = regexp(src, ['(?m)^function [^\n]*[ =]' fname '\s*\('], 'start');
testCase.assertEqual(numel(starts), 1, sprintf('expected exactly one definition of %s', fname));
rest = src(starts(1):end);
nxt = regexp(rest(2:end), '(?m)^function ', 'start', 'once');
if isempty(nxt); body = rest; else; body = rest(1:nxt); end
end

function localAssertPooled(testCase, body, label)
code = regexprep(body, '(?m)%[^\n]*', '');   % comments cannot satisfy a pin
code = regexprep(code, '\.\.\.\s*\n\s*', ' '); % join line continuations, so formatting cannot fail a pin
testCase.verifyTrue(~isempty(regexp(code, ...
    'pi_mean(_e)?\s*=\s*max\(vc(_e)?\.var_mu\s*-\s*vc(_e)?\.sigma_p2\s*-\s*vc(_e)?\.sigma_i2,\s*0\)', 'once')), ...
    sprintf('%s: pi_mean must be derived as max(var_mu - sigma_p2 - sigma_i2, 0) from psyrat_gamma_varcomps', label));
testCase.verifyTrue(~isempty(regexp(code, 'sig_pi_e2\s*=\s*pi_mean(_e)?\s*\+\s*2\s*\.\*\s*exp\(', 'once')), ...
    sprintf('%s: the per-participant residual must start from pi_mean', label));
testCase.verifyTrue(isempty(regexp(code, 'sig_pi_e2\s*=\s*2\s*\.\*\s*exp\(', 'once')), ...
    sprintf('%s: the un-pooled single-exponential residual must be gone', label));
end
