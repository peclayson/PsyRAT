function proc = psyrat_canonicalize_proc(proc, varargin)
%Canonicalize the two alias pairs in a proc preferences struct
%
%proc = psyrat_canonicalize_proc(proc)
%proc = psyrat_canonicalize_proc(proc,'quiet',true)
%proc = psyrat_canonicalize_proc(proc,'quiet',true,'pairs','ssubjrel')
%
%This is the single shared reconciler for the two settings that live in the
%proc struct under two names each (council code-architecture finding #5):
%
%  diffrescor / diffwpcov - ONE setting under two names (residual covariance
%   for co-occurring difference-score events: 1 = fixed to 0, 2 = estimated).
%   diffwpcov (the user-facing GUI name) is decisive on a disagreeing
%   in-range pair, the same direction as the engine's own normalization in
%   psyrat_computevarcomp, which is the policy of record (B33 ruling).
%  ssubjrel / sserrvar - canonical name + legacy alias for the subject-level
%   reliability toggle (1 = off, 2 = on). ssubjrel wins on a disagreeing
%   pair; sserrvar is always rewritten as its mirror.
%
%The blocks this function extracted lived inline in psyrat_computevarcompwarp
%(which warned on conflicts) and as five inline blocks across four
%psyrat_startproc callbacks (all silent); the end states reproduce those blocks
%exactly, pinned
%row-for-row by TestCanonicalizeProc's golden table against the pre-change
%code. Only warning EMISSION differs by call site, via 'quiet'.
%
%Input
% proc - the proc sub-struct of psyrat_prefs (fields other than the four
%  alias fields pass through untouched)
%
%Optional name/value inputs
% quiet - true: suppress the conflict warnings (the GUI call sites' historical
%  behavior); false (default): warn on conflicts (the warp's behavior)
% pairs - 'both' (default): reconcile both pairs; 'ssubjrel': reconcile only
%  the ssubjrel/sserrvar pair. The restriction is load-bearing, not cosmetic:
%  the psyrat_prefs_back and psyrat_exec call sites historically reconciled
%  only that pair, and clamping the diff pair there would change what
%  psyrat_preflight_validate sees for hand-built diffrescor=2/diffwpcov=1
%  structs (its localDiffCovRequested treats the names symmetrically).
%
%Output
% proc - the same struct with diffrescor/diffwpcov and ssubjrel/sserrvar
%  canonicalized (per 'pairs') and mirrored
%

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

%lightweight option parse. Misuse fails loudly: a typo'd option silently
%reconciling the diff pair at a call site that must not is exactly the
%hazard class this function exists to remove. All production call sites are
%enumerated (the warp plus four psyrat_startproc callbacks), so this guard
%is a developer seam, not a user surface.
quiet = false;
pairs = 'both';
if mod(numel(varargin), 2)
    error('psyrat_canonicalize_proc:badoption', ...
        ['Optional inputs must be name/value pairs. How to fix: call as ', ...
        'psyrat_canonicalize_proc(proc,''quiet'',true) or with ', ...
        '''pairs'',''ssubjrel''.']);
end
for k = 1:2:numel(varargin)
    switch varargin{k}
        case 'quiet'
            quiet = logical(varargin{k+1});
        case 'pairs'
            pairs = varargin{k+1};
        otherwise
            error('psyrat_canonicalize_proc:badoption', ...
                ['Unknown option ''%s''. How to fix: only ''quiet'' and ', ...
                '''pairs'' are accepted.'], char(string(varargin{k})));
    end
end
if ~any(strcmp(pairs, {'both', 'ssubjrel'}))
    error('psyrat_canonicalize_proc:badoption', ...
        ['Option ''pairs'' must be ''both'' or ''ssubjrel'', not ''%s''. ', ...
        'How to fix: omit it to reconcile both pairs.'], char(string(pairs)));
end

if strcmp(pairs, 'both')
    %diffwpcov and diffrescor are ONE setting under two names (residual
    %covariance for co-occurring difference-score events: 1 = fixed to 0,
    %2 = estimated). The GUI exposes a single control and writes both fields
    %from it, and psyrat_run reconciles them before its warp call, so both
    %entry points arrive consistent and the clamp below is a no-op for them.
    %It stays as the backstop for callers that build psyrat_prefs by hand --
    %and it ANNOUNCES itself (unless quiet) when it actually changes
    %diffrescor. Swapping the residual-covariance structure silently would
    %refit a structurally different model than the caller asked for, and
    %every difference-score coefficient, SEM and cut-score from that run
    %would be computed under the wrong structure with no runtime signal.
    %Values outside {1,2} are left untouched, exactly as the pre-extraction
    %warp block left them (a deliberate, pinned non-fix: the engine now
    %range-rejects both names at its own boundary, which is where the
    %protection belongs because it also covers direct engine callers).
    %NOTE the order is load-bearing: diffwpcov's missingness is captured
    %BEFORE diffrescor is backfilled, so a lone diffrescor=2 still promotes.
    diffwpcov_missing = ~isfield(proc,'diffwpcov') || isempty(proc.diffwpcov);
    if ~isfield(proc,'diffrescor') || isempty(proc.diffrescor)
        proc.diffrescor = 1;
    end
    if diffwpcov_missing
        if proc.diffrescor == 2
            proc.diffwpcov = 2;
        else
            proc.diffwpcov = 1;
        end
    end
    if any(proc.diffwpcov == [1 2])
        if proc.diffrescor ~= proc.diffwpcov && ~quiet
            warning('psyrat_prefs:diffrescorClamped',...
                ['Found conflicting residual-covariance settings ' ...
                '(diffrescor=%d, diffwpcov=%d). These are one setting under two ' ...
                'names, so PsyRAT is following diffwpcov and estimating with ' ...
                'diffrescor=%d. How to fix: set one of the two and leave the ' ...
                'other unset, or set both to the same value.'],...
                proc.diffrescor,proc.diffwpcov,...
                proc.diffwpcov);
        end
        proc.diffrescor = proc.diffwpcov;
    end
end

%ssubjrel is the canonical public name for the subject-level reliability
%toggle; sserrvar is the legacy alias kept for backwards compatibility
%(and baked into persisted REL fields and analysis labels, which is why it
%cannot simply be retired). The canonical name wins on a disagreement --
%the OPPOSITE direction from the diff pair above, on purpose: there the
%newer user-facing name is decisive, here the canonical one is.
hasLegacy = isfield(proc,'sserrvar') && ~isempty(proc.sserrvar);
hasCanonical = isfield(proc,'ssubjrel') && ~isempty(proc.ssubjrel);

if ~hasCanonical && hasLegacy
    proc.ssubjrel = proc.sserrvar;
elseif ~hasCanonical && ~hasLegacy
    proc.ssubjrel = 1;
elseif hasCanonical && hasLegacy && (proc.ssubjrel ~= proc.sserrvar) && ~quiet
    warning('psyrat_prefs:sserrvarMismatch',...
        ['Found conflicting subject-level reliability settings ',...
        '(ssubjrel=%d, sserrvar=%d). PsyRAT will use ssubjrel. ',...
        'How to fix: set one value and keep both fields synchronized.'],...
        proc.ssubjrel,proc.sserrvar);
end
proc.sserrvar = proc.ssubjrel;

end
