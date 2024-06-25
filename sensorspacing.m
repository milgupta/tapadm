function [EstPosErr, Result] = SensorSpacing(DataFileName, channels, result, comment, samplerate);
% SENSORSPACING calculates the inter-channel spacing in terms of position and orientation.
%           Plots spatial spacing and orientation gap between the given channels.  
%           If channels is a 1x2-vector the function opens two figures,
%           plotting the distance and the cross correlation of distance and
%           coordinates.
%           If more than two channels are given, i.e. channels is a 1xN
%           vector, no graphical output is generated, but results
%           are printed for all combinations of the given channels. Both
%           graphical and textual output is supressed, when called with one
%           or more output arguments.
%
%           [EstPosErr, Result] = SensorSpacing(DataFileName, channels);
%
%           The 1st return value is mean(spatial_dist_iqr)/2 - an
%           rough estimation of the spatial accuracy.
%
%           The Function returns a complex structure as 2nd result:
%              Channels        Nx2 matrix with the different channel pairs
%              SpatialDist     Median and IQR for each pair, total range 
%              OrientationGap   (the same for orientations)
%              Comment         the comment of the original data
%
%           Instead of loading the data from a file, it can be directly
%           supplied to the function, e.g. to use it as an objective
%           function for optimization: 
%
%           EstPosErr = SensorSpacing('', channels, result, comment, samplerate);
%           
%           SensorSpacing needs the statistic toolbox

%---------------------------------------------------------------------
% Copyright © 2005-2006 by Andreas Zierdt (Anderas.Zierdt@phonetik.uni-muenchen.de)
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

% load data and build time axis
if (nargin < 5)
	[result, comment, descriptor, dimension, private, samplerate] = loadsdata(DataFileName);
end
NPoints = size(result, 1); time = (1 / samplerate) * (1 : NPoints)';	

if (length(channels) > 2)							% compare more channels
	combination = sortrows(combnk(channels, 2)); N = size(combination, 1);
	spatial_dist_med = zeros(1, N);		spatial_dist_iqr = zeros(1, N);
	orientation_dist_med = zeros(1, N);	orientation_dist_iqr = zeros(1, N);
	if (nargout == 0)			% if an output argument is ommitted, print statistics
		disp('chan	PosIQR[mm]	OriIQR[°]	MedianPos[mm]	MedianOri[°]');
	end
	for (i=1:N)
		[spatial_dist, spatial_dist_med(i), spatial_dist_iqr(i), orientation_dist, orientation_dist_med(i),...
				orientation_dist_iqr(i)] = calvValues(result, combination(i, 1), combination(i, 2));		
		if (nargout == 0)			% if an output argument is ommitted, print statistics
			disp(sprintf('%d-%d:\t   %3.2f   \t  %3.2f   \t  %3.2f   \t     %3.2f',...
				[combination(i,:) spatial_dist_iqr(i) orientation_dist_iqr(i) spatial_dist_med(i) orientation_dist_med(i)]))

		end
	end
else												% compare two channels
	combination = channels;
	spatial_dist_med = 0;		spatial_dist_iqr = 0;
	orientation_dist_med = 0;	orientation_dist_iqr = 0;

	[spatial_dist, spatial_dist_med, spatial_dist_iqr, orientation_dist, orientation_dist_med,...
			orientation_dist_iqr] = calvValues(result, combination(1), combination(2));
	
	if (nargout == 0)			% only generate pots if a output argument is ommitted
		doPlot(time, NPoints, DataFileName, combination(1), combination(2), samplerate, result, spatial_dist, spatial_dist_med,...
				 spatial_dist_iqr, orientation_dist, orientation_dist_med, orientation_dist_iqr);
	end
end

PosData = struct('EstErr', mean(spatial_dist_iqr)/2, 'Median', spatial_dist_med, 'IQR', spatial_dist_iqr, 'Range', range(spatial_dist), 'SDev', std(spatial_dist), 'Unit', 'mm');
OrientData = struct('EstErr', mean(orientation_dist_iqr)/2, 'Median', orientation_dist_med, 'IQR', orientation_dist_iqr, 'Range', range(orientation_dist), 'SDev', std(orientation_dist), 'Unit', '°');
Result = struct('Channels', combination', 'SpatialDist', PosData, 'OrientationGap', OrientData, 'Comment', comment);
EstPosErr = Result.SpatialDist.EstErr;

%------------- subfunctions -----------

function [spatial_dist, spatial_dist_med, spatial_dist_iqr, orientation_dist, orientation_dist_med,...
			orientation_dist_iqr] = calvValues(result, c1, c2)
	% subtract the data of the two channels to compare, i.e. calculate the difference
	D = result(:, :, c1) - result(:, :, c2);
	m = nanmedian(D); q = iqr(D); ur= m+q/2; lr= m-q/2;

	% calculate the spatial distance of the two channels and it's Median 
	% and Interquartile Range (both NaN insensitive)
	spatial_dist = sqrt(D(:, 1) .^2 + D(:, 2) .^2 + D(:, 3).^2); 
	spatial_dist_med = nanmedian(spatial_dist); spatial_dist_iqr = iqr(spatial_dist); 

	% for the orientation, calculate the angle between the sensors of both channels
	phi = result(:, 4, c1) * pi/180; theta = result(:, 5, c1) * pi/180; [x1, y1, z1] = sph2cart(phi, theta, 1);
	phi = result(:, 4, c2) * pi/180; theta = result(:, 5, c2) * pi/180; [x2, y2, z2] = sph2cart(phi, theta, 1);
	orientation_dist = acos(x1 .* x2 + y1 .* y2 + z1 .* z2) * 180/pi; 

	orientation_dist_med = nanmedian(orientation_dist); orientation_dist_iqr = iqr(orientation_dist); 
return



function doPlot(time, NPoints, DataFileName, c1, c2, samplerate, result, spatial_dist, spatial_dist_med, spatial_dist_iqr, orientation_dist, orientation_dist_med, orientation_dist_iqr)
	phi = result(:, 4, c2) * pi/180; theta = result(:, 5, c2) * pi/180; 
	ur= spatial_dist_med + spatial_dist_iqr/2; lr= spatial_dist_med - spatial_dist_iqr/2;
	figure(1);
	subplot(2, 1, 1)
	plot(time, spatial_dist, 'b-', get(gca, 'XLim'), repmat(spatial_dist_med, 1, 2), 'r-', get(gca, 'XLim'), [lr lr], 'g:', get(gca, 'XLim'), [ur ur], 'g:'); 
	grid on; xlabel('t [s]'); ylabel('spatial spacing [mm]'); 
	legend('sensor spacing', ['Median ' num2str(spatial_dist_med, '%4.1f') 'mm' ], ['IQR ' num2str(spatial_dist_iqr, '%4.1f') ' mm'], 0);
	t = title([DataFileName '    spatial spacing between channel ' num2str(c1) ' and ' num2str(c2) ' (total range is ' num2str(range(spatial_dist), '%4.1f') ' mm, std dev '...
			num2str(std(spatial_dist), '%4.1f') ' mm)']);  set(t, 'Interpreter', 'none');

	subplot(2, 1, 2)
	ur= orientation_dist_med + orientation_dist_iqr; lr= orientation_dist_med - orientation_dist_iqr;
	plot(time, orientation_dist, 'b-', get(gca, 'XLim'), repmat(orientation_dist_med, 1, 2), 'r-', get(gca, 'XLim'), [lr lr], 'k:', get(gca, 'XLim'), [ur ur], 'k:'); 
	grid on; xlabel('t [s]'); ylabel('orientation gap [°]'); 
	%axis([get(gca, 'XLim') (orientation_dist_med - 2*orientation_dist_iqr) (orientation_dist_med + 2*orientation_dist_iqr)]); 
	legend('sensor spacing', ['Median ' num2str(orientation_dist_med, '%4.1f') '°' ], ['IQR ' num2str(orientation_dist_iqr, '%4.1f') ' °'], 0);
	title(['orientation gap between channel ' num2str(c1) ' and ' num2str(c2) ' (total range is ' num2str(range(orientation_dist), '%4.1f') ' °, std dev '...
			num2str(std(orientation_dist), '%4.1f') ' °)']);
	
	figure(2);
	halftime = -(NPoints-1) : (NPoints-1); halftime = (1 / samplerate) * halftime;
	subplot(2, 1, 1)	
	c = xcorr(spatial_dist, orientation_dist, 'coeff');
	c_phi = xcorr(spatial_dist, phi, 'coeff');
	c_theta = xcorr(spatial_dist, theta, 'coeff');
	c_x = xcorr(spatial_dist, result(:, 1, c2), 'coeff');
	c_y = xcorr(spatial_dist, result(:, 2, c2), 'coeff');
	c_z = xcorr(spatial_dist, result(:, 3, c2), 'coeff');
	plot(halftime, c, halftime, c_phi, halftime, c_theta, halftime, c_x, halftime, c_y, halftime, c_z, halftime, -c, 'k:'); grid on;
	xlabel('lag [s]'); ylabel('correlation [1]'); 
	legend('\Delta Pos \otimes \Delta Orientation', '\Delta Pos \otimes \Phi', '\Delta Pos \otimes \Theta',...
			 '\Delta Pos \otimes X', '\Delta Pos \otimes Y', '\Delta Pos \otimes Z', 1);
	title('cross correlation of sensor spacing and coordinates');
	
	subplot(2, 1, 2)
	c = xcorr(orientation_dist, orientation_dist, 'coeff');
	c_phi = xcorr(orientation_dist, phi, 'coeff');
	c_theta = xcorr(orientation_dist, theta, 'coeff');
	c_x = xcorr(orientation_dist, result(:, 1, c2), 'coeff');
	c_y = xcorr(orientation_dist, result(:, 2, c2), 'coeff');
	c_z = xcorr(orientation_dist, result(:, 3, c2), 'coeff');
	plot(halftime, c, halftime, c_phi, halftime, c_theta, halftime, c_x, halftime, c_y, halftime, c_z, halftime, -c, 'k:'); grid on;
	xlabel('lag [s]'); ylabel('correlation [1]'); 
	legend('\Delta Orient \otimes \Delta Orient', '\Delta Orient \otimes \Phi', '\Delta Orient \otimes \Theta',...
			 '\Delta Orient \otimes X', '\Delta Orient \otimes Y', '\Delta Orient \otimes Z', 1);
	title('cross correlation of orientation gap and coordinates');
return