function hex = psyrat_hash_struct(in)
%Hash a MATLAB variable deterministically using JSON encoding.
%
%Returns a lowercase SHA-256 hex digest.

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

try
    payload = jsonencode(in);
catch
    payload = evalc('disp(in)');
end

bytes = uint8(payload);
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(bytes);
hex = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));

end
