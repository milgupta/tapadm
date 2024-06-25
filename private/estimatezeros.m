function [XZeros, XMax] = estimatezeros(x, y, accuracy);
% ESTIMATEZEROS estimates all nulls of the data with the given accuracy.
%           The function uses spline-interpolation to estimate all zeros
%           of y. The Result is returned as cell array, containing a vector
%           with x-values (zeros) for every column of y.

%---------------------------------------------------------------------
% Copyright © 2006 by Andreas Zierdt (Anderas.Zierdt@phonetik.uni-muenchen.de)
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

xi = (x(1) : accuracy : (x(end)-accuracy))';
yi = interp1(x, y, xi, 'spline');

n = find(iszeroc(yi(:,1)));

[NumPoints, NumCoils] = size(yi);

idx1 = (1:(NumPoints-1))'; idx2 = (2:NumPoints)';
S = abs((sign(yi(idx1,:)) - sign(yi(idx2,:)))) /2;
sign_changed = logical(S); 

XZeros = cell(1, NumCoils);
XMax = zeros(1, NumCoils);
for coil = 1:NumCoils
	XZeros{1, coil} = xi(find(sign_changed(:,coil))) + 0.5 *accuracy;
	[max_val, max_idx] = max(abs(yi(:,coil))); XMax(coil) = xi(max_idx);
end

