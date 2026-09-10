function lbl = psyrat_trt_coeflabel(psyrat_data)
%Coefficient-estimand label for test-retest (TRT) reliability outputs.
%
%lbl = psyrat_trt_coeflabel(psyrat_data)
%
%A test-retest run reports ONE of three user-selected coefficients, and they
%are not near-neighbors: on a single dataset the same variance components gave
%0.40 (equivalence), 0.35 (stability) and 0.24 (equivalence and stability).
%The selection is recorded nowhere else in an exported table, so an archived
%file reading "Dependability 0.24" cannot be reconstructed without this label.
%
%The label names three things that together fix the estimand:
% - the coefficient type (reltype): which facets are generalized over
%     CE  - coefficient of equivalence (generalize over trials)
%     CS  - coefficient of stability (generalize over occasions)
%     CES - coefficient of equivalence and stability (both)
% - the error type (gcoeff): dependability (absolute) vs generalizability
%   (relative), which differ by whether the item/trial main effect is included
% - for difference-score designs, the SECOND error type (diffgcoeff), which the
%   user selects independently and which governs the "diff score" row that sits
%   in the same exported table as the per-condition rows
% - the occasion estimand (n'_o), delegated to psyrat_trt_occlabel
%
%UNLIKE psyrat_trt_occlabel, this function returns an EMPTY label rather than
%falling back to a default. That asymmetry is deliberate. n'_o = 1 was the
%historical behavior, so a struct with no occasion setting can be described
%accurately as single-occasion. There is no comparable historical default for
%reltype -- CE, CS and CES are three different estimands and no one of them is
%what a struct without the field "must have been" -- so guessing would be worse
%than silence. Callers are expected to skip the line when the label is empty.
%
%Note that the coefficient type is read from the NUMERIC relsummary.reltype
%rather than from relsummary.reltype_name. psyrat_criterionfigures overrides
%the numeric field on its local copy without updating the name, so on that path
%the two can disagree and only the numeric field reflects what was computed.
%
%Inputs
% psyrat_data - the toolbox results struct (uses psyrat_data.relsummary)
%
%Outputs
% lbl - a human-readable label for the selected TRT coefficient estimand, or
%   '' when the results struct records no reltype (all one-facet designs, and
%   the subject-level test-retest design, which computes from a reltype but
%   does not store it)

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

%no reltype recorded -> no label (see the header note above)
lbl = '';

if ~isstruct(psyrat_data) || ~isfield(psyrat_data,'relsummary')
    return;
end
rs = psyrat_data.relsummary;
if ~isstruct(rs) || ~isfield(rs,'reltype') || ~isscalar(rs.reltype)
    return;
end

%coefficient type. An unrecognized code falls through with an empty label for
%the same reason an absent one does: better silent than mislabeled.
switch rs.reltype
    case 1
        lbl = 'Coefficient of Equivalence (CE)';
    case 2
        lbl = 'Coefficient of Stability (CS)';
    case 3
        lbl = 'Coefficient of Equivalence and Stability (CES)';
    otherwise
        return;
end

%error type. Absent on a legacy struct, in which case the clause is simply
%omitted rather than guessed.
gname = local_errtype(rs,'gcoeff');
if ~isempty(gname)
    lbl = [lbl ', ' gname];
end

%Difference-score designs carry a SECOND, independently selected error type.
%psyrat_relfigures reads gcoeff and diffgcoeff from separate preferences and
%passes BOTH into psyrat_relsummary, which stores both; psyrat_depoverallt then
%prints the per-condition rows at gcoeff and a "diff score" row at diffgcoeff in
%the SAME table. Naming only gcoeff would put the wrong error type beside half
%of that table, which is precisely the failure this label exists to prevent.
%Emitted only when the two differ, so the common case (both dependability,
%their shared default) keeps a single clause and existing output is unchanged.
dname = local_errtype(rs,'diffgcoeff');
if ~isempty(dname) && isfield(rs,'gcoeff') && ~isempty(rs.gcoeff) && ...
        ~isequal(rs.diffgcoeff, rs.gcoeff)
    lbl = [lbl '; difference score: ' dname];
end

%occasion estimand, from the existing helper rather than a second copy of the
%same logic. It always returns a label, so no guard is needed here.
lbl = [lbl '; ' psyrat_trt_occlabel(psyrat_data)];

end

%-------------------------------------------------------------------------
% Local helpers
%-------------------------------------------------------------------------
function name = local_errtype(rs,fieldname)
%Human-readable error type for a g-coefficient selector field: 1 = absolute
%error (dependability), 2 = relative error (generalizability). Returns '' when
%the field is absent, empty, non-scalar or unrecognized, so the caller omits the
%clause rather than asserting an error type the run did not record.
name = '';
if ~isfield(rs,fieldname) || ~isscalar(rs.(fieldname))
    return;
end
switch rs.(fieldname)
    case 1
        name = 'Dependability (absolute error)';
    case 2
        name = 'Generalizability (relative error)';
end
end
