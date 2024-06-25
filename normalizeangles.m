function  RA = NormalizeAngles(A)
% NORMALIZEANGLES a small function to make sure, that phi and theta, are in the
%            range: -180 < phi <= +180 and -90 < theta <= +90.
%
%            [phi theta] = NormalizeAngles([phi theta])

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

phi_low_idx = find(A(:, 1) <= -180); phi_high_idx = find(A(:, 1) > 180); 
if (~isempty(phi_low_idx))
A(phi_low_idx, 1) = A(phi_low_idx, 1) + 360;
end
if (~isempty(phi_high_idx))
	A(phi_high_idx, 1) = A(phi_high_idx, 1) - 360; 
end	

[x, y, z] = sph2cart(A(:, 1)*pi/180, A(:, 2)*pi/180, 1); O = [x y z];

theta_out_idx = find((A(:, 2) <= -90) | (A(:, 2) > 90)); 
O(theta_out_idx, :) = O(theta_out_idx, :) * -1; 		
[p, t, z] = cart2sph(O(:,1), O(:,2), O(:,3));

RA = [p t] * 180/pi;

