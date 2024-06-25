function CA = calcamps(Coordinates, CalibrationFactors);
% CALCAMPS calculate signal amplitudes for the given coordinates.
%           Coordinates has to be a Nx5xM-matrix, containing for N samples
%           the spatial x-,y- and z-position (in mm) and the orientation 
%           angles phi and theta (in degrees) for M channels.
%           i.e input has scale and orientation (so-called autokal
%           coordinates) of calcpos or tapad (matlab version) output.
%           If Coordinates has more than 5 columns, redundant columns 
%           are ignored. M can be singleton.
%           CalibrationFactors is an optional 1x6 vector. When omitted,
%           the calculated amplitude for a sensor parallel to a 
%           transmitter coil at position [0 0 0] is set to 1.
%
%           CA = calcamps(Coordinates, CalibrationFactors);
%
%           The function returns a Nx6 matrix with amplitudes for every
%           spatial point (row) and transmitter (column). 

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
%munlock;
[NumPoints, NumCoordinates, NumChannels] = size(Coordinates);
if (NumCoordinates < 5)
   error('Matrix must contain at least 5 columns, not %d.', NumCoordinates);
end

persistent CFK;   % compute default calibration factors
if (isempty(CFK)) 
	[LocalPos, LocalOrientVec] = trans2local([0 0 0], [1 0 0]); % Nx3x6
	HK = kohlrausch(LocalPos(:,:,1)); 
	HK = sum(HK(:,:) .* LocalOrientVec(:, :, 1), 2);
  	CFK = 1/HK * [1 -1 -1 1 1 -1];
end   
if (nargin < 2)
  	CalibrationFactors = CFK;
else
  	CalibrationFactors = CalibrationFactors .* CFK;
end   

CA= zeros(NumPoints, 6, NumChannels);	% calculated amplitudes

for Channel=1:NumChannels
	SensorPositions = Coordinates(:, 1:3, Channel);
	phi = Coordinates(:, 4, Channel) * pi/180; theta = Coordinates(:, 5, Channel) * pi/180; 
	[x, y, z] = sph2cart(phi, theta, 1); 
   SensorOrientations = [x y z]; 
   SensorPositions = 1/1000 * rotat(SensorPositions); SensorOrientations = rotat(SensorOrientations);
	[LocalPos, LocalOrientVec] = trans2local(SensorPositions, SensorOrientations); % Nx3x6
	for Coil=1:6
		H = kohlrausch(LocalPos(:,:,Coil));
		CA(:, Coil, Channel) = sum(H(:,:) .* LocalOrientVec(:, :, Coil), 2) * CalibrationFactors(Coil);
	end
end

