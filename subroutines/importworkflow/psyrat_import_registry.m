function adapters = psyrat_import_registry()
%Return registered ERP import adapters.

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

adapters = struct([]);

adapters(1).name = 'eeglab';
adapters(1).display_name = 'EEGLAB';
adapters(1).detect_fn = @psyrat_adapter_eeglab_detect;
adapters(1).read_fn = @psyrat_adapter_eeglab_read;
adapters(1).version = '1.0.0';

adapters(2).name = 'erplab';
adapters(2).display_name = 'ERPLAB';
adapters(2).detect_fn = @psyrat_adapter_erplab_detect;
adapters(2).read_fn = @psyrat_adapter_erplab_read;
adapters(2).version = '1.0.0';

adapters(3).name = 'ep_toolkit_mat';
adapters(3).display_name = 'EP Toolkit MAT';
adapters(3).detect_fn = @psyrat_adapter_epmat_detect;
adapters(3).read_fn = @psyrat_adapter_epmat_read;
adapters(3).version = '1.0.0';

end
