function methods = psyrat_scoring_methods()
%Canonical list of supported single-trial scoring methods.
%
% methods = psyrat_scoring_methods()
%
% Single source of truth shared by psyrat_scoring_spec_create (method
% validation) and the import/scoring GUI (method picker), so the two cannot
% drift out of sync. Order is the display order used in the GUI.

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

methods = {'absolute_peak_amplitude', ...
    'local_peak_amplitude', ...
    'time_window_mean_amplitude', ...
    'peak_latency', ...
    'centroid_latency', ...
    'fractional_area_latency', ...
    'integrated_area_amplitude'};
end
