classdef PsyRAT3ARecoveryReference
    % Independent planted-truth reference for the bounded 3A STATIC subject-level
    % gamma recovery cells: the one-facet sserr design (analysis 6), the
    % non-concurrent two-event difference (analysis 8) and its concurrent
    % Gaussian-copula twin (analysis 8, diffrescor = 2).
    %
    % WHY A SEPARATE HELPER. Those live recovery tests used to take their planted
    % truth from the production converters (psyrat_gamma_varcomps,
    % psyrat_gamma_crosscov) and left the population person x trial mean term M
    % out of the per-subject residual, while production (since S20) pools it. A
    % truth that calls production cannot detect a production error, and a truth
    % without M is stale against the accepted target. Every quantity here is
    % written from the approved log-nu DGP algebra
    % (the 2C S20 recovery reference and its scientific contract, local evidence)
    % as explicit exp/expm1 products. No production
    % routine, no accuracy oracle, no clipping. Deterministic pins live in
    % tests/TestGamma3ARecoveryReference.m.
    %
    % APPROVED DGP (continuous scaled chi-square, one event):
    %   log mu_pi = alpha + U_p + U_i,   U_p ~ N(0, sP^2),  U_i ~ N(0, sI^2)
    %   nu_p      = exp(beta + V_p)
    %   Y | mu, nu_p ~ Gamma(shape = nu_p/2, scale = 2 mu/nu_p)   E = mu, Var = 2 mu^2/nu_p
    %
    % POPULATION mean components (expanding the two independent centered
    % lognormal factors; B = (E mu)^2, the SQUARED marginal mean of mu at one
    % cell, not the second moment E(mu^2) = exp(2 alpha + 2 sP^2 + 2 sI^2)):
    %   B = exp(2 alpha + sP^2 + sI^2)
    %   P = B expm1(sP^2)                 person universe-score variance (signal)
    %   I = B expm1(sI^2)                 trial main effect (absolute error only)
    %   M = B expm1(sP^2) expm1(sI^2)     person x trial mean interaction
    %
    % OWN (person-conditioned) observation variance, averaged over trials:
    %   C_p = 2 exp(2 (alpha + u_p) + 2 sI^2 - beta - v_p)
    %
    % POOLED own residual, added exactly once:   E_p = M + C_p.
    % M is a POPULATION quantity at the population intercept: it is NOT
    % recomputed at alpha + u_p, and the exp(2 sI^2) inside C_p is a lognormal
    % second-moment factor that does not contain M.
    %
    % Coefficients of the n-trial mean (relative / absolute error):
    %   G_p = P / (P + E_p / n),        D_p = P / (P + E_p / n + I / n)
    %
    % TWO EVENTS (difference score), m_e = exp(alpha_e + (sP_e^2 + sI_e^2)/2):
    %   CovP12 = m_1 m_2 expm1(cP12),   cP12 = rho_mean sP_1 sP_2
    %   signal = P_1 + P_2 - 2 CovP12   (the person covariance is SIGNAL, not error)
    %   non-concurrent error  = E_1/n_1 + E_2/n_2          (residual covariance zero)
    %   concurrent, equal n   = (E_1 + E_2 - 2 K)/n,       K = E_i Cov(Y1, Y2 | p, i)
    % The population mean-interaction cross-covariance is
    %   m_1 m_2 expm1(cP12) expm1(cI12),
    % zero in the recovery fixtures ONLY because their trial effects are
    % independent across events (cI12 = 0). It is nonzero only when BOTH cP12
    % and cI12 are nonzero. It is exposed so a test can show that the fixture
    % zero follows from cI12 = 0 and is not a general rule.
    %
    % MONTE-CARLO UNCERTAINTY of a simulated K (sample covariance of N iid pairs):
    % influence-function variance mcse^2 = var((x - xbar)(y - ybar)) / N, and the
    % first-order propagation dG/dK = 2 S / [n (S + err)^2].
    %
    % INPUT CONTRACT. Log-scale SDs (sP, sI) must be finite and nonnegative;
    % variance-type inputs (M, C, E, P, I, S, err) finite and nonnegative; trial
    % counts finite and strictly positive. Every public method checks its inputs
    % against that contract and its result for finiteness, and errors otherwise.
    % Nothing is clipped: an invalid reference must fail, not become valid.

    methods (Static)

        function c = oneFacetMeanComponents(alpha, sP, sI)
            % Population observed-scale mean components of one event (see header;
            % B = (E mu)^2 is the squared marginal mean, not the second moment).
            % expm1(0) is exactly 0, so a zero factor SD gives an exactly zero
            % component without any clipping. Scalars, or equal-size arrays.
            % Negative SDs are refused (a squared negative SD would silently
            % evaluate as a valid variance).
            PsyRAT3ARecoveryReference.requireFinite('oneFacetMeanComponents', alpha, sP, sI);
            PsyRAT3ARecoveryReference.requireNonnegative('oneFacetMeanComponents', sP, sI);
            vP = sP .^ 2;
            vI = sI .^ 2;
            B  = exp(2 .* alpha + vP + vI);
            c  = struct();
            c.B = B;
            c.P = B .* expm1(vP);
            c.I = B .* expm1(vI);
            c.M = B .* expm1(vP) .* expm1(vI);
            PsyRAT3ARecoveryReference.requireFinite('oneFacetMeanComponents', c.B, c.P, c.I, c.M);
        end

        function C = ownObservationVariance(alpha, u, sI, beta, v)
            % Person p's OWN trial-averaged observation variance,
            %   C_p = 2 exp(2 (alpha + u_p) + 2 sI^2 - beta - v_p)
            %       = E_i[ 2 mu_pi^2 / nu_p | u_p, v_p ],
            % exact because nu_p is constant over trials (ID-only nu). u and v are
            % per-person arrays of equal size; alpha, sI and beta are the
            % population constants of the event. A negative sI is refused.
            PsyRAT3ARecoveryReference.requireFinite('ownObservationVariance', alpha, u, sI, beta, v);
            PsyRAT3ARecoveryReference.requireNonnegative('ownObservationVariance', sI);
            C = 2 .* exp(2 .* (alpha + u) + 2 .* sI .^ 2 - beta - v);
            PsyRAT3ARecoveryReference.requireFinite('ownObservationVariance', C);
        end

        function E = pooledOwnResidual(M, C)
            % Pooled own residual E_p = M + C_p: the POPULATION mean interaction
            % (one scalar per event, at the population intercept) added exactly
            % once to each person's own observation variance. Both must be finite
            % and nonnegative variances.
            PsyRAT3ARecoveryReference.requireFinite('pooledOwnResidual', M, C);
            PsyRAT3ARecoveryReference.requireNonnegative('pooledOwnResidual', M, C);
            E = M + C;
            PsyRAT3ARecoveryReference.requireFinite('pooledOwnResidual', E);
        end

        function [G, D] = oneFacetCoefficients(P, I, E, n)
            % Per-person coefficients of the n-trial mean:
            %   G = P / (P + E/n)          relative error  (psyrat_ssrel gcoeff 2)
            %   D = P / (P + E/n + I/n)    absolute error  (psyrat_ssrel gcoeff 1)
            % A zero P gives exactly 0; a zero denominator is refused, not repaired.
            % n must be a finite, strictly positive count.
            PsyRAT3ARecoveryReference.requireFinite('oneFacetCoefficients', P, I, E);
            PsyRAT3ARecoveryReference.requireNonnegative('oneFacetCoefficients', P, I, E);
            PsyRAT3ARecoveryReference.requireCount('oneFacetCoefficients', n);
            G = P ./ (P + E ./ n);
            D = P ./ (P + E ./ n + I ./ n);
            PsyRAT3ARecoveryReference.requireFinite('oneFacetCoefficients', G, D);
        end

        function m = eventMean(alpha, sP, sI)
            % Marginal expected score of one event, m = exp(alpha + (sP^2 + sI^2)/2).
            % Negative SDs are refused.
            PsyRAT3ARecoveryReference.requireFinite('eventMean', alpha, sP, sI);
            PsyRAT3ARecoveryReference.requireNonnegative('eventMean', sP, sI);
            m = exp(alpha + (sP .^ 2 + sI .^ 2) ./ 2);
            PsyRAT3ARecoveryReference.requireFinite('eventMean', m);
        end

        function c = personCrossCovariance(m1, m2, cP12)
            % Observed-scale cross-event PERSON covariance (the difference-signal
            % off-diagonal): CovP12 = m1 m2 expm1(cP12), with cP12 = rho_mean sP1 sP2
            % the log-scale person covariance. Exactly zero at cP12 = 0; takes the
            % sign of cP12.
            PsyRAT3ARecoveryReference.requireFinite('personCrossCovariance', m1, m2, cP12);
            c = m1 .* m2 .* expm1(cP12);
            PsyRAT3ARecoveryReference.requireFinite('personCrossCovariance', c);
        end

        function c = trialMeanCrossCovariance(m1, m2, cP12, cI12)
            % Population cross-event covariance of the person x trial mean
            % interactions, m1 m2 expm1(cP12) expm1(cI12). It vanishes when the
            % trial effects are independent across events (cI12 = 0), which is the
            % ONLY reason the recovery fixtures carry no mean-interaction
            % off-diagonal; it is nonzero only when BOTH cP12 and cI12 are nonzero
            % (correlated trial effects alone, with cP12 = 0, still give zero).
            PsyRAT3ARecoveryReference.requireFinite('trialMeanCrossCovariance', m1, m2, cP12, cI12);
            c = m1 .* m2 .* expm1(cP12) .* expm1(cI12);
            PsyRAT3ARecoveryReference.requireFinite('trialMeanCrossCovariance', c);
        end

        function err = differenceError(E1, n1, E2, n2, K)
            % Per-person relative error of the difference of two n-trial means:
            %   non-concurrent (K = 0):     E1/n1 + E2/n2      (residual covariance zero)
            %   concurrent, equal counts:   (E1 + E2 - 2 K)/n
            % K = E_i Cov(Y1, Y2 | p, i) is the person's concurrent observation
            % cross-covariance. The concurrent form is defined here for EQUAL
            % counts only (the fixture design) and is deliberately not generalized
            % to unequal counts. A negative error is refused, never clipped.
            if nargin < 5 || isempty(K)
                K = 0;
            end
            PsyRAT3ARecoveryReference.requireFinite('differenceError', E1, E2, K);
            PsyRAT3ARecoveryReference.requireNonnegative('differenceError', E1, E2);
            PsyRAT3ARecoveryReference.requireCount('differenceError', n1, n2);
            if any(K(:) ~= 0) && ~isequal(n1, n2)
                error('PsyRAT3ARecoveryReference:unequalCounts', ...
                    ['The concurrent difference error is defined for equal trial ', ...
                    'counts only in this reference (n1 = %g, n2 = %g).'], n1, n2);
            end
            err = E1 ./ n1 + E2 ./ n2 - 2 .* K ./ n1;
            PsyRAT3ARecoveryReference.requireFinite('differenceError', err);
            PsyRAT3ARecoveryReference.requireNonnegative('differenceError', err);
        end

        function G = differenceCoefficient(S, err)
            % G = S / (S + err), S = P1 + P2 - 2 CovP12 the difference signal.
            PsyRAT3ARecoveryReference.requireFinite('differenceCoefficient', S, err);
            PsyRAT3ARecoveryReference.requireNonnegative('differenceCoefficient', S, err);
            G = S ./ (S + err);
            PsyRAT3ARecoveryReference.requireFinite('differenceCoefficient', G);
        end

        function [mcse, ifvar, K] = covarianceMcse(x, y)
            % Monte-Carlo standard error of the sample covariance of N iid pairs,
            % from its influence function IF(x, y) = (x - mu_x)(y - mu_y) - sigma_xy:
            %   p_i   = (x_i - xbar)(y_i - ybar)
            %   K     = sum(p_i)/(N - 1)        (the off-diagonal cov(x, y) returns)
            %   ifvar = var(p)                   ((N - 1)-normalized)
            %   mcse  = sqrt(ifvar / N)
            % Deterministic in the draws: no random numbers are consumed.
            x = x(:);
            y = y(:);
            N = numel(x);
            if numel(y) ~= N || N < 3
                error('PsyRAT3ARecoveryReference:badMcseInput', ...
                    'covarianceMcse needs two equal-length vectors of at least 3 pairs.');
            end
            p = (x - mean(x)) .* (y - mean(y));
            K = sum(p) ./ (N - 1);
            ifvar = var(p);
            mcse = sqrt(ifvar ./ N);
            PsyRAT3ARecoveryReference.requireFinite('covarianceMcse', mcse, ifvar, K);
        end

        function d = dGdK(S, err, n)
            % First-order sensitivity of the concurrent difference coefficient to
            % the simulated cross-covariance K at equal counts n:
            %   G = S / (S + err),  err = (E1 + E2 - 2 K)/n
            %   dG/dK = 2 S / [n (S + err)^2]
            % Multiply by the K Monte-Carlo standard error to propagate it to G.
            % S and err must be finite nonnegative; n a finite positive count.
            PsyRAT3ARecoveryReference.requireFinite('dGdK', S, err);
            PsyRAT3ARecoveryReference.requireNonnegative('dGdK', S, err);
            PsyRAT3ARecoveryReference.requireCount('dGdK', n);
            d = 2 .* S ./ (n .* (S + err) .^ 2);
            PsyRAT3ARecoveryReference.requireFinite('dGdK', d);
        end

    end

    methods (Static, Access = private)

        function requireFinite(where, varargin)
            % Every reference value must be a finite real number. Nothing is
            % repaired: a non-finite value is an error in the fixture or the
            % algebra and must surface as one.
            for k = 1:numel(varargin)
                x = varargin{k};
                if ~isnumeric(x) || ~isreal(x) || ~all(isfinite(x(:)))
                    error('PsyRAT3ARecoveryReference:nonFinite', ...
                        '%s produced a non-finite or non-real reference value (argument %d).', ...
                        where, k);
                end
            end
        end

        function requireNonnegative(where, varargin)
            % Variance-type quantities and SDs are refused when negative rather
            % than clipped to zero (the production converters clip; this
            % reference must not, or it could not detect an invalid fixture).
            for k = 1:numel(varargin)
                x = varargin{k};
                if any(x(:) < 0)
                    error('PsyRAT3ARecoveryReference:negativeValue', ...
                        ['%s received or produced a negative variance-type quantity ', ...
                        '(argument %d); it is refused rather than clipped.'], where, k);
                end
            end
        end

        function requireCount(where, varargin)
            % Trial counts must be finite real numbers strictly greater than
            % zero; a zero, negative, NaN or Inf count would otherwise produce a
            % silently invalid (Inf/NaN/negative) projection.
            for k = 1:numel(varargin)
                n = varargin{k};
                if ~isnumeric(n) || ~isreal(n) || isempty(n) || ~all(isfinite(n(:))) || any(n(:) <= 0)
                    error('PsyRAT3ARecoveryReference:badCount', ...
                        '%s requires finite, strictly positive trial counts (argument %d).', where, k);
                end
            end
        end

    end
end
