function comp = psyrat_dod_composite(varargin)
%PSYRAT_DOD_COMPOSITE Signed four-cell composite variance w'*Sigma*w for the
% person-specific dynamic difference-of-differences designs (analyses 28/29).
%
%   comp = psyrat_dod_composite('varcov',V,'cvec',w,'rule','unscaled')
%   comp = psyrat_dod_composite('varcov',V,'cvec',w,'rule','harmonic', ...
%       'n',[n1 n2 n3 n4])
%
% Required Inputs:
%  varcov - component covariance matrix draws (ndraws x 4 x 4). Diagonal =
%           per-cell variances; off-diagonals = SIGNED cross-cell covariances.
%  rule   - 'unscaled'  the plain quadratic form w'*Sigma*w
%           'harmonic'  the D-study scaled form: variance terms divided by
%                       their own cell's replication count, covariance terms by
%                       the PAIRWISE HARMONIC MEAN h(n_a,n_b) = 2/(1/n_a+1/n_b)
%
% Optional Inputs:
%  cvec - signed contrast vector (default [1 -1 -1 1]) corresponding to
%         (ERP_1 - ERP_2) - (ERP_3 - ERP_4)
%  n    - per-cell replication counts [n1 n2 n3 n4]; REQUIRED for 'harmonic'
%  tol  - relative tolerance for the tiny-negative floor (default 1e-10)
%
% Output struct comp:
%  value            - the composite variance draws (ndraws x 1). Tiny
%                     floating-point negatives (within tol of zero at the
%                     draw's own magnitude scale) are set to exactly zero;
%                     MATERIAL negatives are returned UNCHANGED and flagged.
%                     Nothing else is modified.
%  invalid_negative - 1 where the draw is materially negative (ndraws x 1)
%  scale            - the magnitude scale used for the tolerance (ndraws x 1)
%
% WHY THIS KERNEL EXISTS RATHER THAN AN EDIT TO PSYRAT_DODIFFREL. That kernel
% serves the STATIC DoD design (analysis 9): its residual channel is
% diagonal-only (er_var * (c.^2./obs)), and its admissibility commentary is
% load-bearing precisely BECAUSE that channel carries no signed off-diagonals.
% The dynamic nonconcurrent designs need a residual channel with signed
% off-diagonals (the mean-surface person-by-trial cross terms survive even
% though the six OBSERVATION-level residual covariances are fixed to zero), so
% the shipped kernel cannot be reused without rewriting its audited
% reachability arguments. This kernel is pinned AGAINST it instead: on the
% shared subspace (harmonic rule, a covariance matrix with the static kernel's
% structure) the two must agree to 1e-12, and TestGammaDodDynrelLsConversion
% asserts exactly that. Same construction as PSYRAT_GAMMA_CROSSCOV_RESCOR
% versus PSYRAT_GAMMA_CROSSCOV: never move shipped numbers.
%
% THE FLAGGING CONVENTION (reference: reliability.R, composite_component and
% scaled_composite_component; provenance in tests/baselines/nonconcurrent_dod).
% The scale is the sum of absolute variance terms (scaled by n under
% 'harmonic') plus twice the sum of absolute UNSCALED covariances - the
% reference's own choice, transcribed rather than improved, so the frozen
% fixtures pin it. A value below -tol*max(1,scale) is a genuinely negative
% composite: it is reported as-is with invalid_negative = 1, mirroring
% PSYRAT_ADMISSIBILITY's record-never-repair rule.
%
% COUNT RULE SCOPE (owner ruling 2026-08-12): the harmonic rule is the only
% scaled rule implemented. The reference bundle additionally defines an
% exact-shared-count sensitivity (Cov * n_shared/(n_a*n_b)); it is frozen in
% the fixtures but deliberately NOT implemented here. A later increment can
% add rule 'exactshared' without touching the two rules above.
%
% See also PSYRAT_DODIFFREL, PSYRAT_GAMMA_DOD_COMPONENTS_LS,
% PSYRAT_HARMMEAN.

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

if mod(length(varargin),2)
    error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',...
        'Make sure each variable input is paired with a value \n',...
        'See help psyrat_dod_composite for more information about inputs'));
end

varcov = local_req(varargin,'varcov','Component covariance draws (varcov) were not specified.');
rule = local_req(varargin,'rule','Scaling rule (rule) was not specified.');

ind = find(strcmpi('cvec',varargin),1);
if ~isempty(ind)
    cvec = varargin{ind+1};
else
    cvec = [1 -1 -1 1];
end

ind = find(strcmpi('tol',varargin),1);
if ~isempty(ind)
    tol = varargin{ind+1};
else
    tol = 1e-10;
end

if ~any(strcmpi(rule,{'unscaled','harmonic'}))
    error('varargin:rule',... %Error code and associated error
        'rule must be ''unscaled'' or ''harmonic''.');
end

if ndims(varcov) ~= 3 || size(varcov,2) ~= 4 || size(varcov,3) ~= 4
    error('varargin:varcov',... %Error code and associated error
        'varcov must be ndraws x 4 x 4 covariance draws.');
end

if numel(cvec) ~= 4
    error('varargin:cvec',... %Error code and associated error
        'cvec must contain four signed contrast weights.');
end
cvec = cvec(:);

if ~isscalar(tol) || ~isreal(tol) || isnan(tol) || tol < 0
    error('varargin:tol',... %Error code and associated error
        'tol must be a real nonnegative scalar.');
end

if strcmpi(rule,'harmonic')
    ind = find(strcmpi('n',varargin),1);
    if isempty(ind)
        error('varargin:n',... %Error code and associated error
            'Per-cell replication counts (n) are required for rule=''harmonic''.');
    end
    n = varargin{ind+1};
    if numel(n) ~= 4
        error('varargin:n',... %Error code and associated error
            'n must include four per-cell replication counts [n1 n2 n3 n4].');
    end
    n = n(:)';
    if any(~isreal(n)) || any(isnan(n)) || any(n <= 0)
        error('varargin:n',... %Error code and associated error
            'All n values must be real, numeric and > 0.');
    end
else
    %the unscaled rule divides every term by 1
    n = ones(1,4);
end

ndraws = size(varcov,1);

%Variance terms carry their own cell's count; each covariance pair carries the
%pairwise harmonic mean (which degenerates to the shared count when the two
%counts are equal, so the same expression serves both same-count and
%different-count pairs). Under 'unscaled' every divisor is 1 and this is the
%plain quadratic form.
vars = zeros(ndraws,4);
for q = 1:4
    vars(:,q) = varcov(:,q,q);
end

value = (vars ./ n) * cvec.^2;
abscov = zeros(ndraws,1);
for a = 1:3
    for b = (a+1):4
        nharm = 2 / ((1/n(a)) + (1/n(b)));
        value = value + 2 .* cvec(a) .* cvec(b) .* varcov(:,a,b) ./ nharm;
        abscov = abscov + abs(varcov(:,a,b));
    end
end

%The reference's magnitude scale: scaled absolute variance terms plus twice the
%absolute UNSCALED covariances (reliability.R, scaled_composite_component).
scale = sum(abs(vars ./ n),2) + 2 .* abscov;

if any(~isfinite(value))
    error('psyrat_dod_composite:nonfinite',... %Error code and associated error
        ['A composite variance draw is non-finite. How to fix: check the '...
        'component covariance draws for Inf/NaN before calling.']);
end

%Tiny-negative floor at the draw's own magnitude; material negatives are
%FLAGGED and returned unchanged (record, never repair).
threshold = tol .* max(1, abs(scale));
tiny = value < 0 & value >= -threshold;
value(tiny) = 0;
invalid_negative = double(value < -threshold);

comp = struct;
comp.value = value;
comp.invalid_negative = invalid_negative;
comp.scale = scale;

end

function val = local_req(args,name,msg)
ind = find(strcmpi(name,args),1);
if isempty(ind)
    error(['varargin:' lower(name)],... %Error code and associated error
        msg);
end
val = args{ind+1};
end
