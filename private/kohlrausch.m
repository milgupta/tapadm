function H = kohlrausch(LocalPos);
% KOHLRAUSCH calculates the magnetic field vector H for the given
%           local coordinates. (see trans2local.m)  
%           If size(LocalPos) = Nx3, the result will be a a Nx3 matrix, 
%           containing the local x-,y- and z-component of the magnetic 
%           field-vectors.
%     
%           H = kohlrausch(LocalPos);
%
%           The function uses a fast approximation for the exact formular, 
%           without the need to solve complete elliptic integrals. The 
%           approximation becomes less accurate in the vicinity of the              
%           transmitter coil, but has proved to yield a fair basic 
%           field-model for 3D-EMA.

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

% switch to local polar-coordinates
[phi, theta, r] = Cart2Sph(LocalPos(:,1), LocalPos(:,2), LocalPos(:,3));

Er = angles2vec(phi, theta);			% radial unit vector (in cartesian coordinates)

thetat = abs(theta) - pi/2;
idx = find(theta<0);
thetat(idx) = -thetat(idx);
Et = angles2vec(phi, thetat);			% tangential unit vector
Et(idx, :) = -Et(idx, :);

Fr = 1./(r.*r.*r); 						% f(r) =  1 / r.^3;	

Hr = repmat((2 .* sin(theta)), 1, 3) .* Er;			% radial field component	
Ht = repmat(cos(theta), 1, 3) .* Et;					% tangential

H = repmat(Fr, 1, 3) .* (Hr + Ht);

%------------- subfunctions -----------
function V = angles2vec(phi, theta)
	[N, dummy] = size(phi);
	V = zeros(N, 3);	
	[V(:,1), V(:,2), V(:,3)] = sph2cart(phi, theta, ones(N,1));
