function [data, comment, descriptor, dimension, private, samplerate, unit] = loadsdata(filename);
% LOADSDATA Load AG500 'IPSK'-structured data from the given 'mat'-file.  
%           The function guarantees that all of the expected variables are
%           defined, whether they actually exist in the file or not. 
%           If the file defines additional variables, they are omitted. 
%           If the given file does not exist, the function just displays a 
%           warning and proceeds. 
%     
%           [data, comment, descriptor, dimension, private, samplerate, unit]...
%                = loadsdata(filename);
%
%           If called without argument, the function returns a set of empty
%           variables, this feature can be used for a defined initialization.

%---------------------------------------------------------------------
% Copyright © 2005-2007 by Andreas Zierdt (Anderas.Zierdt@phonetik.uni-muenchen.de)
% 
% This program is free software; you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 2 of the License, or (at your option) 
% any later version.
%
% This program is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with 
% this program; if not, write to the Free Software Foundation, Inc., 
% 51 Franklin Street, Fifth Floor, Boston, MA 
%---------------------------------------------------------------------

data = []; comment = ''; descriptor = ''; dimension = struct('descriptor', '', 'unit', '', 'axis', '');
private = struct('Iterations', -1); samplerate = []; unit = '';

if (nargin ~= 1),
	return
end

if (isempty(strfind(filename, '.mat'))) % MATLAB position data
   error('Unsupported file type.');
end   

if exist(filename, 'file')
	FData = load(filename); 
	fn = fieldnames(FData);
	if (~isempty(strmatch('data', fn, 'exact'))) data = double(FData.data); end
	if (~isempty(strmatch('comment', fn, 'exact'))) comment = FData.comment; end
	if (~isempty(strmatch('descriptor', fn, 'exact'))) descriptor = FData.descriptor; end
	if (~isempty(strmatch('dimension', fn, 'exact'))) dimension = FData.dimension; end
	if (~isempty(strmatch('private', fn, 'exact'))) private = FData.private; end
	if (~isempty(strmatch('samplerate', fn, 'exact'))) samplerate = FData.samplerate; end
	if (~isempty(strmatch('unit', fn, 'exact'))) unit = FData.unit; end
else
	disp(['nonexistig file: ''' filename '''']);
end
	
