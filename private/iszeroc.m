function liszc = iszeroc(data);
% ISZEROC True for zero crossing.
%         ISZEROC(X) returns 1's where two elements of X are of different sign for
%         the smaller one. Returns 0's otherwise. 

%---------------------------------------------------------------------
% Copyright © 2005 by Andreas Zierdt (Anderas.Zierdt@phonetik.uni-muenchen.de)
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

shifted_data = data; shifted_data(length(data)+1) = data(1); shifted_data(1) = [];
sign_changed = (sign(data) ~= sign(shifted_data)); sign_changed(end) = 0;
idx = find(sign_changed);					% vector with indices where the sign changes
reord = abs(shifted_data(sign_changed)) < abs(data(sign_changed));	% If the 2'nd point is lesser,
idx = idx + reord;																	% take it.

iszc = zeros(size(sign_changed));		% change indices back to logical vector
iszc(idx) = ones(size(idx));
liszc = logical(iszc);
