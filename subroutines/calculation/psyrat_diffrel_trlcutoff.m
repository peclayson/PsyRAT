function [trlcutoff,diffscore] = psyrat_diffrel_trlcutoff(varargin)
%Find the number of trials at which a difference score reaches a reliability
%cutoff, using a single trial count common to both measures
%
%[trlcutoff,diffscore] = psyrat_diffrel_trlcutoff('bp',id_varcov,...
%  'bt',trl_varcov,'er_var',b_sigma,'wp_cov',wp_cov,'est','dep','CI',.95,...
%  'depcutoff',.70,'meascutoff',1,'ntrials',1100)
%
%Inputs
% bp - between-person (co)variance estimates from CmdStan (id_varcov)
% bt - between-trial (co)variance estimates from CmdStan (trl_varcov)
% er_var - log residual SD draws from CmdStan (b_sigma); exponentiated and
%  squared internally
% est - type of estimate: 'dep' - dependability, 'gen' - generalizability
% CI - size of the credible interval in decimal format: .95 = 95%
% depcutoff - reliability threshold the difference score has to reach
% meascutoff - which part of the credible interval has to reach the
%  threshold: 1 - lower limit, 2 - point estimate, 3 - upper limit
% ntrials - largest trial count to search
%
%Optional Inputs
% wp_cov - within-person (residual) covariance between the two measures,
%  passed through to psyrat_diffrel. Default: 0 (non-concurrent events).
%
%Outputs
% trlcutoff - smallest common number of trials per measure at which the
%  difference score reaches depcutoff, or -1 if no count up to ntrials does
% diffscore - the psyrat_diffrel output at that trial count, or [] when no
%  cutoff was found. Callers are expected to fall back to the difference
%  score at the OBSERVED trial counts in that case.
%
%Why a single common trial count
% The per-event trial cutoffs are derived independently, each chosen so that
% its own event separately clears the threshold. Evaluating the difference
% score at that pair answers a question nobody asked, because neither count
% was chosen with the difference score in mind. This routine instead asks the
% question a user reading a difference-score cutoff row believes is being
% answered: how many trials per condition does the DIFFERENCE score need?
%
% Holding the two counts equal also removes ONE of the two ways the coefficient
% can leave [0,1] -- not both. The difference-score error variance scales its
% covariance term by the harmonic mean of the two counts (Rocha et al., 2026,
% Table 6), which falls below their geometric mean whenever the counts differ,
% so the error variance could go negative. psyrat_diffrel warns and RECORDS
% when that happens; its recorder's help is explicit that "It RECORDS. It does
% not clip, clamp, floor" (see psyrat_admissibility), so the out-of-range value
% reaches the caller with its diagnostic attached rather than being repaired.
% At a common n' THAT cause does not arise at all: the harmonic, geometric and
% arithmetic means coincide, the error variance reduces to the
% single-observation variance of the difference divided by n', and it is
% non-negative for any positive semi-definite input.
%
% What a common n' does NOT rule out is the other cause -- a raw cross-condition
% block that is not itself positive semi-definite, which the
% gamma/scaled-chi-square family can produce because its variance diagonals are
% clamped while its cross-covariances are not (RC-23). That one is negative at
% ANY counts. Since this routine evaluates at equal counts BY DESIGN, it is
% precisely the regime where the non-PSD cause is the only one left, so a
% negative here should never be read as a projection artifact. See
% psyrat_admissibility's help for the full two-cause statement. The cutoff row
% is therefore in range whenever the raw blocks are valid covariance
% structures -- whenever a cutoff is found. The not-found fallback reports the
% observed-trial score instead, and that one IS evaluated at unequal counts
% (see psyrat_relsummary).
%
%Note on the search
% Every posterior draw's coefficient is uni / (uni + errvar / n) with both
% terms non-negative, which is monotonically increasing in n. Point estimates
% and credible-interval limits inherit that monotonicity, so the first trial
% count that reaches the threshold is the smallest one that does, and the
% search can stop there. This mirrors how the per-event cutoffs are found in
% psyrat_relsummary (the psyrat_rel_sing call over a range of trial counts,
% followed by find(...,1)).

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

if isempty(varargin) || mod(length(varargin),2)
    error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',...
        'Make sure each variable input is paired with a value \n',...
        'See help psyrat_diffrel_trlcutoff for more information about inputs'));
end

%the (co)variance inputs are handed straight to psyrat_diffrel, which does its
%own checking of them, so they are only pulled out here
bp = local_reqinput(varargin,'bp','Between-person (co)variances');
bt = local_reqinput(varargin,'bt','Between-trial (co)variances');
er_var = local_reqinput(varargin,'er_var','Error variances');
est = local_reqinput(varargin,'est','Reliability estimate type');
ciperc = local_reqinput(varargin,'CI','Size of the credible interval');
depcutoff = local_reqinput(varargin,'depcutoff','Reliability threshold');
meascutoff = local_reqinput(varargin,'meascutoff','Credible-interval measure');
ntrials = local_reqinput(varargin,'ntrials','Largest trial count to search');

%residual covariance is optional and defaults to the non-concurrent case
ind = find(strcmpi('wp_cov',varargin),1);
if ~isempty(ind)
    wp_cov = varargin{ind+1};
else
    wp_cov = [];
end

if length(meascutoff) ~= 1 || ~any(meascutoff == [1 2 3])
    %G58 alignment: this guard is the precedent the psyrat_relsummary parse
    %guards mirror, so it gains the same scalar length test (a vector slips
    %past the elementwise any(); an empty errors inside it) and the same
    %strcat spacing fix (strcat strips a char argument's TRAILING
    %whitespace, MEASURED R2025b, so the separating space lives at the
    %head of the next fragment).
    error('varargin:meascutoff',... %Error code and associated error
        strcat('WARNING: meascutoff should be 1 (lower limit),',...
        ' 2 (point estimate), or 3 (upper limit)\n',...
        'See help psyrat_diffrel_trlcutoff for more information \n'));
end

%round down and bound the search: a fractional or non-positive ntrials would
%otherwise make the loop below silently do nothing
ntrials = floor(ntrials);
if ~isfinite(ntrials) || ntrials < 1
    error('varargin:ntrials',... %Error code and associated error
        strcat('WARNING: ntrials should be a finite value of at least 1\n',...
        'See help psyrat_diffrel_trlcutoff for more information \n'));
end

%default state: no cutoff found. -1 is the same "not reached" sentinel the
%per-event trial cutoffs use in psyrat_relsummary.
trlcutoff = -1;
diffscore = [];

for ntrl = 1:ntrials

    %evaluate the difference score with BOTH measures at the same trial count.
    %psyrat_diffrel is called rather than reimplementing the coefficient here,
    %so this routine cannot drift away from the production formula.
    ds = psyrat_diffrel(...
        'bp',bp,...
        'bt',bt,...
        'er_var',er_var,...
        'wp_cov',wp_cov,...
        'obs',[ntrl ntrl],...
        'est',est,...
        'CI',ciperc);

    %compare whichever part of the credible interval the user asked for
    switch meascutoff
        case 1
            relmeas = ds.ll;
        case 2
            relmeas = ds.pt;
        case 3
            relmeas = ds.ul;
    end

    %monotone in ntrl, so the first hit is the smallest trial count that works
    if relmeas >= depcutoff
        trlcutoff = ntrl;
        diffscore = ds;
        break;
    end

end

end

function val = local_reqinput(args,name,description)
%Pull a required name/value pair out of varargin, erroring by name if it is
%missing. Written once here rather than repeated per input, because this
%routine takes eight required inputs and eight near-identical error blocks
%would obscure the search itself.

ind = find(strcmpi(name,args),1);

if isempty(ind)
    error('varargin:missing',... %Error code and associated error
        strcat('WARNING:',sprintf(' %s',description),' not specified \n\n',...
        'Please input',sprintf(' %s',name),...
        '. See help psyrat_diffrel_trlcutoff for more information \n'));
end

val = args{ind+1};

end
