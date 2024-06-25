function [result, residual, iter] = calcpos(amps, amp_derivative, startpoint, options);
% CALCPOS calculates sensor positions and orientations for amplitudes.
%           Coordinates has to be a Nx5-matrix, containing for N samples
%           the spatial x-,y- and z-position and the orientation 
%           angles phi and theta (in radians).
%           CalibrationFactors is an optional 1x6 vector.
%
%           [result, residual, iter] = calcpos(amps, startpoint);
%
%           The function returns a Nx7 matrix with amplitudes for every
%           spatial point (row) and transmitter (column), plus 'rms-value'
%           plus exit condition (column 7): 
%            >0 The function converged to a solution x.
%            0  The maximum number of function calls/iterations was exceeded.
%            <0 The function did not converge to a solution.

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

if (nargin < 4)
	options = optimset('lsqnonlin');
	options = optimset(options, 'Display', 'off', 'LargeScale', 'on', 'LevenbergMarquardt', 'off');
	options = optimset(options, 'TolX', 1E-4, 'TolFun', 1E-6, 'MaxFunEvals', 1000);
end
if ((nargin < 3) | isempty(startpoint))
   startpoint = zeros(1, 5);
end
persistent rmsfac;
if (isempty(rmsfac))
	rmsfac = 2500 / sqrt(6);
end

weight = abs(amp_derivative);
[result, resnorm, residual, exitflag, output] = lsqnonlin(@amperr, startpoint, [], [], options, amps, weight);
result = [result (sqrt(resnorm)*rmsfac) exitflag]; 
iter = output.iterations;

%------------- subfunctions -----------
function da = amperr(point, amps, weight)
   a0 = calcamps(point);
   da = (amps - a0) .* weight;

% performs a stereographic projection on the orientation angles   
function oplanepos = angles2plane(oangles)
   oangles = oangles * pi/180;
   [x, y, z] = sph2cart(rem(oangles(:,1), 2*pi), rem(oangles(:,2), 2*pi), 1);
   o = rotat([x y z]);                 % transform orientation to TAPAD coordinates
   d = 1 - o(:,3);
	oplanepos(:, 1) = o(:, 1) .* 2 ./ d;
	oplanepos(:, 2) = o(:, 2) .* 2 ./ d;
   opplanepos = oplanepos * 100;
   
 function oangles = plane2angles(oplanepos)
   opplanepos = oplanepos / 100;
   d = oplanepos(:, 1) .* oplanepos(:, 1) + oplanepos(:, 2) .* oplanepos(:, 2); 
	o(:,3) = d; d = d + 4;
   o(:,1) = oplanepos(:, 1) .* 4 ./ d;
   o(:,2) = oplanepos(:, 2) .* 4 ./ d; d = d ./2;
   o(:,3) = o(:,3) ./ d -1;
   o = rotta(o);                       % transform orientation to Autokal coordinates
   [phi, theta, r] = cart2sph(o(:,1), o(:,2), o(:,3));
   oangles = [phi theta] * 180/pi;
