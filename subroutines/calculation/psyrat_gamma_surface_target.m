function target = psyrat_gamma_surface_target(target_in)
%PSYRAT_GAMMA_SURFACE_TARGET Validate and canonicalize the gamma dynamic
% reliability SURFACE target (Q9/GUI-1, owner rulings 9-12 and 2026-09-13).
%
%   target = psyrat_gamma_surface_target(target_in)
%
% target_in - [] or '' (default -> 'population_marginal'), or one of
%             'population_marginal'  the primary surface: each draw's residual is
%                                    M + C_pop, the population mean remainder plus
%                                    the observation variance MARGINALIZED over the
%                                    person scale population BEFORE the coefficient
%                                    is formed (not an average of person
%                                    coefficients);
%             'typical_person'       the companion: M + C_typ, the observation
%                                    variance at zero person effects (log-nu: both
%                                    u_p = v_p = 0; location-scale: w_p = 0).
% Anything else errors with 'dynrel:badsurfacetarget'; nothing is coerced.
%
% This is the surface target only. It never reaches the per-participant
% calculators, whose subject target is fixed by the design and the stored
% remainder scheme (see psyrat_dynrel_summary).
%
% See also PSYRAT_DYNREL_SUMMARY, PSYRAT_REL_DYNREL_GAMMA,
% PSYRAT_REL_DYNREL_GAMMA_LS, PSYRAT_REL_DYNREL_TRT_GAMMA,
% PSYRAT_REL_DYNREL_TRT_GAMMA_LS.

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

if nargin < 1 || isempty(target_in)
    target = 'population_marginal';
    return;
end
if ~(ischar(target_in) || (isstring(target_in) && isscalar(target_in)))
    error('dynrel:badsurfacetarget',...
        ['The gamma surface target must be the text ''population_marginal'' or ' ...
        '''typical_person''. How to fix: pass one of those two strings as ' ...
        '''surface_target'' (omit it for the population-marginal primary).']);
end
target = char(target_in);
if ~any(strcmp(target,{'population_marginal','typical_person'}))
    error('dynrel:badsurfacetarget',...
        ['Unknown gamma surface target ''%s''. Valid values are ' ...
        '''population_marginal'' (primary; observation variance marginalized ' ...
        'over the person scale population before the coefficient) and ' ...
        '''typical_person'' (companion; zero person effects in the observation ' ...
        'variance). How to fix: pass one of those two strings as ' ...
        '''surface_target''.'], target);
end

end
