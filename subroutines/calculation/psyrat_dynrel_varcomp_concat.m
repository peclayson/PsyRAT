function T = psyrat_dynrel_varcomp_concat(summ)
%Stack the per-stratum dynamic-reliability variance-component tables into one.
%
% T = psyrat_dynrel_varcomp_concat(summ)
%
%Single source of truth for the dynamic-reliability variance-component export:
%each stratum's variance-component table, prepended with the stratum label and
%its observed number of trials (n_trials = stratum .obs). Used by BOTH the
%dynrel viewer (psyrat_startview_dynrel) and the headless report (psyrat_report's
%local_dynrel_tables) so the GUI export and the headless export are identical.
%
%Input
% summ - the struct returned by psyrat_dynrel_summary (per-stratum .varcomp
%        table plus .obs trial count).
%
%Output
% T - the stacked table (stratum, n_trials, then the original variance-component
%     columns).

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

parts = cell(numel(summ.strata),1);
for s = 1:numel(summ.strata)
    st = summ.strata(s);
    vc = st.varcomp;
    h = height(vc);
    vc.stratum = repmat(string(st.label),h,1);
    %The gamma difference variant carries PER-EVENT n' (st.obs_ev) and its surface
    %is computed at those exact counts; the scalar st.obs is only their mean and
    %would stamp a count neither event has (the reference ERP design is 49 error
    %vs 341 correct trials). Emit the two exact counts as their own NUMERIC
    %columns there -- mirroring the surface table in psyrat_report -- rather than
    %stringifying n_trials, so no column ever changes type between designs.
    if isfield(st,'obs_ev') && numel(st.obs_ev) == 4
        %DoD variants: four per-cell counts, itemized for the same reason
        vc.n_trials_cell1 = repmat(st.obs_ev(1),h,1);
        vc.n_trials_cell2 = repmat(st.obs_ev(2),h,1);
        vc.n_trials_cell3 = repmat(st.obs_ev(3),h,1);
        vc.n_trials_cell4 = repmat(st.obs_ev(4),h,1);
        lead = {'stratum','n_trials_cell1','n_trials_cell2',...
            'n_trials_cell3','n_trials_cell4'};
    elseif isfield(st,'obs_ev') && numel(st.obs_ev) == 2
        vc.n_trials_event1 = repmat(st.obs_ev(1),h,1);
        vc.n_trials_event2 = repmat(st.obs_ev(2),h,1);
        lead = {'stratum','n_trials_event1','n_trials_event2'};
    else
        vc.n_trials = repmat(st.obs,h,1);
        lead = {'stratum','n_trials'};
    end
    parts{s} = movevars(vc,lead,'Before','component');
end
T = vertcat(parts{:});

end
