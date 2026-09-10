function [nodes, weights] = psyrat_gauss_hermite(n)
%PSYRAT_GAUSS_HERMITE Gauss-Hermite quadrature nodes/weights for a STANDARD NORMAL
% weight, computed by the Golub-Welsch algorithm.
%
%   [nodes, weights] = psyrat_gauss_hermite(n)
%
% Returns n nodes x_k and weights w_k (both n-by-1, nodes ascending) such that
% for a standard normal Z,
%
%   E[f(Z)] = integral f(z) phi(z) dz  ~=  sum_k w_k * f(x_k),
%
% exact for polynomials f of degree <= 2n-1. The weights sum to 1, so they are
% probability weights rather than the classical physicists' weights (which sum
% to sqrt(pi)); this is the "probabilists'" normalization and is what the
% expectation form above requires.
%
% ALGORITHM (Golub-Welsch). The physicists' Hermite polynomials (weight
% exp(-x^2)) satisfy a three-term recurrence whose Jacobi matrix J is symmetric
% tridiagonal with zero diagonal and off-diagonal sqrt(k/2), k = 1..n-1. Its
% eigenvalues are the physicists' nodes and the squared first components of its
% normalized eigenvectors are the weights, normalized to sum to 1. Converting
% the weight exp(-x^2) to the standard normal density rescales the nodes by
% sqrt(2) and leaves the (already normalized) weights unchanged.
%
% This is a direct port of gauss_hermite_normal() in the owner's reference bundle
% (person_specific_dynamic_concurrent_modular/R/reliability.R), including its
% node ordering, so quadrature results are comparable term by term.
%
% Results are cached across calls: the nodes depend only on n, and the copula
% covariance converter calls this once per quadrature evaluation, which happens
% once per person per posterior draw.
%
% BASE MATLAB ONLY (eig). No Statistics Toolbox dependency; see
% tests/TestBaseMatlabOnly.m for why that matters.
%
% See also PSYRAT_GAMMA_COPULA_RESCOV.

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

% The node set is a pure function of n, so cache it. containers.Map keyed by the
% node count; base MATLAB, and cheap relative to the eigendecomposition.
persistent cache
if isempty(cache)
    cache = containers.Map('KeyType', 'double', 'ValueType', 'any');
end

if ~isnumeric(n) || ~isscalar(n) || ~isfinite(n) || n ~= fix(n) || n < 2
    error('psyrat_gauss_hermite:nodes', ...
        ['Gauss-Hermite quadrature requires an integer node count of at '...
        'least two. How to fix: pass n >= 2.']);
end

if isKey(cache, n)
    cached = cache(n);
    nodes = cached.nodes;
    weights = cached.weights;
    return
end

% Symmetric tridiagonal Jacobi matrix for the physicists' Hermite recurrence.
% Zero diagonal (the recurrence has no shift term); off-diagonal sqrt(k/2).
offdiag = sqrt((1:n-1)/2);
jacobi = diag(offdiag, 1) + diag(offdiag, -1);

% Golub-Welsch: eigenvalues are the nodes, squared first eigenvector components
% are the weights. MATLAB dispatches to the symmetric solver for a real
% symmetric matrix, so the eigenvalues are real and the eigenvectors orthonormal.
[vectors, values] = eig(jacobi);
[lambda, order] = sort(diag(values));

% sqrt(2) converts physicists' nodes (weight exp(-x^2)) to standard-normal nodes.
nodes = sqrt(2) * lambda(:);
% Already sum to 1 because the eigenvectors are orthonormal; no mu_0 factor is
% applied, which is exactly what the probability normalization requires.
weights = (vectors(1, order).^2)';

cache(n) = struct('nodes', nodes, 'weights', weights);

end
