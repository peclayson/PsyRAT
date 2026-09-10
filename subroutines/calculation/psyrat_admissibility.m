function adm = psyrat_admissibility(adm,label,values,fcn,role)
%Record whether a quantity that must be a variance actually is one
%
%adm = psyrat_admissibility(adm,label,values,fcn)
%adm = psyrat_admissibility(adm,label,values,fcn,role)
%
%Inputs
% adm - diagnostic struct to accumulate into; pass [] to start a new one
% label - what the quantity is, in words, as it should appear to a user
%  (e.g. 'relative error variance', 'universe-score variance'). Records
%  sharing a label are ROLLED UP by psyrat_admissibility_note, which sums
%  their draw counts, so a label must name one quantity and one only.
% values - vector of per-draw values for that quantity
% fcn - name of the calling function, recorded so a message can say where
%  the quantity came from
%
%Optional Inputs
% role - where the quantity sits in the coefficient: 'numerator' (a
%  universe-score variance) or 'denominator' (an error variance). Default
%  'numerator', which is the wording every record got before this argument
%  existed. It changes nothing about what is counted; it decides only how the
%  exactly-zero case is described, and that case means OPPOSITE things in the
%  two positions. See "Zero is recorded separately from negative" below.
%
%Outputs
% adm - the diagnostic struct with one more entry appended
%
%What this does and, more importantly, what it does NOT do
% It RECORDS. It does not clip, clamp, floor, or otherwise modify the values,
% and it does not substitute NaN for them. The reliability coefficient is
% reported exactly as the formula produces it, and anything inadmissible is
% surfaced to the user with an explanation so they can judge it themselves.
% Silently repairing an inadmissible quantity produces a number that looks
% ordinary and is not, which is worse than the raw arithmetic: clipping a
% negative error variance to zero, for instance, drives the coefficient to
% exactly 1 when the universe term is also negative.
%
%Why a variance can stop being one
% Difference-score components divide each variance by its own n' and subtract
% a cross-condition covariance scaled by the HARMONIC mean of the two counts
% (Rocha et al., 2026, Table 6). A raw variance of a difference is non-negative
% whenever the cross-condition block is a real covariance matrix, so there are
% two ways to get a negative one, and which applies is decided by whether the
% two counts are equal.
%
% UNEQUAL counts: the harmonic mean falls below the geometric mean, so that
% subtraction can exceed the variances it is taken from. The estimated
% components are not at fault here; it is the PROJECTION to n' observations per
% condition that leaves the space of variances.
%
% EQUAL counts: the harmonic, geometric and arithmetic means coincide and every
% projected block reduces to the raw block divided by its own count, so if every
% raw block were a valid covariance structure each term would be non-negative
% and nothing here could go negative. The projection ALONE therefore cannot be
% the cause -- at least one raw block reaching this recorder is not positive
% semi-definite. That is not hypothetical: the gamma / scaled-chi-square blocks
% are assembled in MATLAB with their variance diagonals clamped at zero while
% the matching cross-covariances are not (RC-23), which is precisely how a block
% stops being a covariance matrix. Equal-count callers are ordinary, not exotic
% -- the difference-score trials curves evaluate at obs = [n n] by construction,
% and the headline difference score does whenever the two conditions have equal
% observed trial counts.
%
% What equal counts do NOT buy you is independence from the counts. Do not
% conclude that admissibility at obs = [n n] is a property of the raw blocks
% alone: the quantities recorded here are SUMS of blocks carrying DIFFERENT
% divisors (n, the occasion count k, and their product n*k), so no single factor
% cancels and the sign can still flip with n. Demonstrated on the production
% kernel at bpi_raw = +18, bpo_raw = -1, k = 2, reltype 3: the relative error
% variance is positive for n <= 20 and negative for n >= 50, from one fixed
% non-PSD block. A warning that disappears when the count changes has not been
% fixed, and a curve over n can fire on part of its range rather than all of it.
%
%Zero is recorded separately from negative, and deliberately
% A universe-score variance of exactly zero is a legitimate, reportable
% result: two conditions that are perfectly correlated with equal variance
% produce a difference score with no between-person variance, and a
% reliability of zero is then the correct answer rather than a failure. A
% NEGATIVE value is arithmetic that has left the space of variances. The two
% deserve different messages, so they are counted separately here.
%
% Exactly zero means the OPPOSITE thing in the denominator, which is why the
% role argument exists. An error variance of exactly zero removes that term
% from the denominator and drives the coefficient toward ONE, not toward zero,
% so describing it with the numerator's wording tells the user the reverse of
% what happened. The mechanism is the same -- error components that are
% perfectly correlated with equal variance cancel in the difference -- but the
% consequence is not, and this recorder is the only place that knows which
% position the quantity was taken from.

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

%A mis-typed role would silently produce user-facing prose that says the
%reverse of what happened, which is the exact failure this argument exists to
%prevent, so it is validated rather than accepted and stored.
if nargin < 5 || isempty(role)
    role = 'numerator';
end

if ~ischar(role) && ~isstring(role)
    error('psyrat_admissibility:role',... %Error code and associated error
        'role must be the text ''numerator'' or ''denominator''.');
end

role = lower(char(role));

if ~any(strcmp(role,{'numerator','denominator'}))
    error('psyrat_admissibility:role',... %Error code and associated error
        'role must be ''numerator'' or ''denominator'', not ''%s''.',role);
end

if isempty(adm) || ~isstruct(adm)
    adm = struct('entries',struct('label',{},'fcn',{},'role',{},'ndraw',{},...
        'nneg',{},'nzero',{},'minval',{}));
end

values = values(:);

%NaN draws are counted separately from negative ones and are NOT treated as
%inadmissible here: a NaN means the quantity was never computed, which is a
%different problem from one that was computed and came out impossible.
finitevals = values(~isnan(values));

entry = struct();
entry.label = label;
entry.fcn = fcn;
entry.role = role;
entry.ndraw = numel(values);
entry.nneg = sum(finitevals < 0);
entry.nzero = sum(finitevals == 0);

if isempty(finitevals)
    entry.minval = NaN;
else
    entry.minval = min(finitevals);
end

adm.entries(end+1) = entry;

end
