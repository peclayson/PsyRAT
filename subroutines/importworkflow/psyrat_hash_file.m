function [hex, info] = psyrat_hash_file(filepath)
%Return SHA-256 hash plus basic file metadata.

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

if exist(filepath,'file') ~= 2
    error('psyrat_import:fileMissing', ...
        'Source file was not found: %s', filepath);
end

fid = fopen(filepath,'r');
if fid < 0
    error('psyrat_import:fileOpenFailed', ...
        'Could not open source file for hashing: %s', filepath);
end

c = onCleanup(@() fclose(fid));
md = java.security.MessageDigest.getInstance('SHA-256');

while true
    [chunk,count] = fread(fid, 1048576, '*uint8');
    if count <= 0
        break;
    end
    md.update(chunk);
end

hex = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));

d = dir(filepath);
if isempty(d)
    info = struct('bytes', NaN, 'datenum', NaN, 'date', '');
else
    info = struct('bytes', d.bytes, 'datenum', d.datenum, 'date', d.date);
end

end
