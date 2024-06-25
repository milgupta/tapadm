function simcal;
% SMARTCAL  Three-dimensional Artikulographic Position and Align Determination.
%           Calculates positions for a set of AG500 Amp-files.
%           BasePath is somethig like 'C:\MyStuf\EMA', i.e. the base 
%           for all data. Hold your BasePath in an global Variable,
%           or use the pwd command.
%           AmpPath and ResPath defines the relative path for input and
%           output data, e.g. AmpPath = 'Data\Amp'; ResPath = 'Data\Pos'; 
%
%           tapad(BasePath, AmpPath, ResPath, TrialIdx, ChanIdx, options);
%
%           If there is already an existing TAPADM result file and TAPADM
%           is used to calculate position data just for some channels, the 
%           function will keep existing position data which is not affected 
%           by this run. Thus, position calculating can be performed in a 
%           sequential order: channel by channel. Please notice, that
%           TAPADM will not preserve other entries in the Result-File! 
%
%           options: -s don't use last result as start point (significantly
%                       increases computation time!) 
%                    -f flip time, i. e. process data onwards from the last point
%                    -d use amplitude derivatives to weight errors
%                    -l use Levenberg-Marquardt instead of Newton method

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
COLORS = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 0 0 0; 1 0 1];
set(0, 'DefaultAxesColorOrder', COLORS);

% parameters of calibration and their start-values
r = 		80;									% disk radius [mm]	
a0 = 		 0 * pi/180;	%65				% disk starting angle (offset)
phi =		45 *pi/180;	 						% sensor orientation (Autokal-system)
theta =	 0 *pi/180;
z0 = 		 0;									% ZA-Offset of the disk


arg = [r 	 a0		phi	 theta 	z0];		% parameter vector for optimisation
lb  = [r-5	 a0-10	phi-10	theta-10	z0-5];		% lower bound for optimisation
ub  = [r+5	 a0+10 	phi+10	theta+10	z0+5];		% upper bound for optimisation

lb = []; ub = [];

da = 3.6 *pi/180;								% desired dalpha 
alphas = 0 : da : 2*pi-da;

Amps = SC_CalcAmps(arg, alphas);
na = NormalizeAngles([arg(3)*180/pi arg(4)*180/pi]);
info = ['Disk radius= ' num2str(arg(1), '%6.2f') ' mm Z-Offset= ' num2str(arg(5), '%6.2f') ' mm alpha0= ' num2str(arg(2)*180/pi, '%6.2f') '°'];
info = [info 'Sensor at disk is ' num2str(na(1), '%6.2f') '° horizontally rotated and ' num2str(na(2), '%6.2f') '° vertically elevated'];
figure; plotAmps(info, alphas, Amps);
figure; drapefig(80); hold on; plotVec(arg, alphas, 6, Amps); hold off; 

figure

Psi = repmat(alphas', 1, 6);
C = abs(Amps) .* (cos(Psi) + i * sin(Psi));
C = [C; C(1, :)];
plot(C); grid on; axis square;
legend('T1', 'T2', 'T3', 'T4', 'T5', 'T6', 3); 
xlabel('cos \alpha  [1]'); ylabel(' i * sin \alpha  [i]'); 
title(info);
hold on;
for coil=1:6
	fill(real(C(:, coil)), imag(C(:, coil)), COLORS(coil, :)); alpha(0.1);
end
	plot(C); 
plot(1.0 * (cos(Psi) + i * sin(Psi)), 'k:');
plot(2.0 * (cos(Psi) + i * sin(Psi)), 'k:');
plot(3.0 * (cos(Psi) + i * sin(Psi)), 'k:');
text(3.2*cos(pi/4), 3.2*sin(pi/4), '\alpha = 45°');
text(2.8*cos(pi/2), 2.8*sin(pi/2), '\alpha = 90°');
text(2.8*cos(3*pi/2), 2.8*sin(3*pi/2), '\alpha = 270°');
hold off;
return

deps = 0.005;
dAdR = (SC_CalcAmps([80+deps 0	45 0 	0], alphas) - SC_CalcAmps([80-deps 0	45 0 	0], alphas)) ./ (2*deps);
dAdA = (SC_CalcAmps([80 0+deps	45 0 	0], alphas) - SC_CalcAmps([80 0-deps	45 0 	0], alphas)) ./ (2*deps);
dAdP = (SC_CalcAmps([80 0	45+deps 0 	0], alphas) - SC_CalcAmps([80 0	45-deps 0 	0], alphas)) ./ (2*deps);
dAdT = (SC_CalcAmps([80 0	45 0+deps 	0], alphas) - SC_CalcAmps([80 0	45 0-deps 	0], alphas)) ./ (2*deps);
dAdZ = (SC_CalcAmps([80 0	45 0 	0+deps], alphas) - SC_CalcAmps([80 0	45 0 	0-deps], alphas)) ./ (2*deps);

figure; plotData(alphas, abs(dAdR)); grid on; title('dAdR');
figure; plotData(alphas, abs(dAdA)); grid on; title('dAdA');
figure; plotData(alphas, abs(dAdP)); grid on; title('dAdP');
figure; plotData(alphas, abs(dAdT)); grid on; title('dAdT');
figure; plotData(alphas, abs(dAdZ)); grid on; title('dAdZ');

return

%------------- subfunctions -----------
%function f= SC_RMSCalFunction(arg, Interp_Alphas, Interp_Amps, nulls)		% arg = [R alpha0 PhiStart ThetaStart L]
function f= SC_RMSCalFunction(arg, nulls)		% arg = [R alpha0 PhiStart ThetaStart L]
	da = 3.6 *pi/180;								% desired resolution dalpha 
	a = 0 : da : 2*pi-da;
	%c_nulls = estimate_nulls(alphas, Amps)

% calculate signal-amplitudes for positions of the 'nulls'-vector (root of the measured data)
	% and build a 12 component error vector with these amps (2 roots per transmitter) 
	%a = nulls; % + arg(2);				% add angle-offset
	CAmps = SC_CalcAmps(arg, a);	
	c_nulls = estimate_nulls(a, CAmps);
	
	f = nulls-c_nulls; %[CAmps(1:2, 1); CAmps(3:4, 2); CAmps(5:6, 3); CAmps(7:8, 4); CAmps(9:10, 5); CAmps(11:12, 6)];
return

function CAmps = SC_CalcAmps(arg, alphas)  % arg = [R alpha0 PhiStart ThetaStart L];
	a0 = arg(2);								% disk starting angle (offset)
	phi 	= arg(3);	 						% sensor orientation (Autokal-system)
	theta = arg(4);

	[x, y, z] = sph2cart(phi, theta, 1); o  = [x y z]; 

	a = alphas + a0;	
	P = SC_CalcPosFromPar(arg, a);		 % generate calibration positions and orientations
	O = SC_CalcOrientFromPar(o, a);

	[phi, theta] = cart2sph(O(:,1), O(:,2), O(:,3));  
	CAmps = calcamps([P phi*180/pi theta*180/pi]);
return

function Positions = SC_CalcPosFromPar(par, alphas)
% SC_CalcPosFromPar(par, dalphas) calculates Autokal-XYZ positions
%           for the given relative rotations of the calibration disk.
% 
%           par = [r a0 phi theta ZOffset];
%           alphas row vector with disk angles relative alpha0 around Za
%           (radian measure)
	R = par(1); alpha0 = par(2); ZOffset = par(5);
	Positions = ZOffset * ones(length(alphas), 3);
	Positions(:,1) = R .* cos(alpha0 + alphas');
	Positions(:,2) = R .* sin(alpha0 + alphas');
return

function Orientations = SC_CalcOrientFromPar(StartOrientation, alphas)
% SC_CalcOrientFromPar(par, dalphas) calculates sensor orientation
%           for the given relative rotations of the calibration disk.
% 
%           StartOrientation = [ox oy oz] (Autokal co-ordinates) 
%           alphas row vector with disk angles relative alpha around Za
%           (radian measure)
	ox = StartOrientation(1); oy = StartOrientation(2); oz = StartOrientation(3); 
	Orientations = oz * ones(length(alphas), 3);
	Orientations(:,1) = ox * cos(alphas') - oy * sin(alphas');
	Orientations(:,2) = ox * sin(alphas') + oy * cos(alphas');
return

function [alphas, Amps] = loaddata(filename) 
	cdata = load(filename);						% Load text-file with AG500 calibration data
	Amps = cdata(:, [2 4 6 8 10 12]);		% get signal amplitudes
	Amps = Amps .* repmat([-1 1 -1 1 -1 1], size(Amps, 1), 1);		% get signal amplitudes
	% get alpha, the rotation angle of the calibration device and convert it to radians, 
	% relative to the first point. That will be the abscissa of the calibration data.
	ENCODER_INCREMENTS_PER_REVOLUTION = 8000;
	alphas = (cdata(:,1) - cdata(1,1)) .* 2*pi/ENCODER_INCREMENTS_PER_REVOLUTION; 
	% map data to one revolution: 0 <= alpha < 2pi 
	idx = find(alphas < 0); alphas(idx) = alphas(idx) + 2*pi;
	[alphas, idx] = sort(alphas); Amps = Amps(idx, :); clear idx; clear cdata;
return

function nulls = estimate_nulls(alphas, Amps)
% interpolate signal amplitudes with high resolution to estimate nulls of the functions
	da = 0.0005 *pi/180;								% desired resolution dalpha 
	Interp_Alphas = 0 : da : 2*pi-da;
	Interp_Amps = interp1(alphas, Amps, Interp_Alphas, 'spline');

	ZIdx = ones(2, 6) * NaN;							% find zero-crossing amplitudes
	for coil=1:6
		n = find(iszeroc(Interp_Amps(:,coil)));
		if (size(n,1) == 2)
			ZIdx(:,coil) = n;
		else
			error(['found ' num2str(size(n,1)) ' nulls, expected 2 for transmitter ' num2str(coil)]);
		end
	end

	nulls = Interp_Alphas(ZIdx(:));				% Estimate alphas where f(alpha)=0 for the measured data
return

function plotAmps(PlotTitle, alphas, Amps);
	% define line styles for plots 
	lsp = ['b:'; 'r:'; 'g:'; 'c:'; 'k:'; 'm:'];
	lsd = ['b.'; 'r.'; 'g.'; 'c.'; 'k.'; 'm.'];
	lsc = ['b-'; 'r-'; 'g-'; 'c-'; 'k-'; 'm-'];
	a = alphas * 180/pi;
	plot(a, Amps(:,1), lsc(1,:), a, Amps(:,2), lsc(2,:), a, Amps(:,3), lsc(3,:),... 
		  a, Amps(:,4), lsc(4,:), a, Amps(:,5), lsc(5,:), a, Amps(:,6), lsc(6,:));
	legend('T1', 'T2', 'T3', 'T4', 'T5', 'T6', 3); 
	xlabel('\alpha [°]'); ylabel('Amp [1]'); grid; li = get(gca, 'YLim'); li = max(abs(li));
	set(gca, 'XLim', [0 360], 'XTick', [0 45 90 135 180 225 270 315 360], 'YLim', [-li +li]); 
	
	title(PlotTitle);
	
function plotVec(arg, alphas, Coil, Amps)
	figure(2);
	P = SC_CalcPosFromPar(arg, alphas);		 % generate calibration positions and orientations
	phi 	= arg(3);	 						% sensor orientation (Autokal-system)
	theta = arg(4);

	[x, y, z] = sph2cart(phi, theta, 1); o  = [x y z]; 
	O = SC_CalcOrientFromPar(o, alphas);

	maxd = 100;
	x = P(:, 1); y = P(:, 2); z = P(:, 3);
	d = sqrt(x .* x + y .* y + z .*z); maxd = ceil(max(d));
	ox = O(:, 1); oy = O(:, 2); oz = O(:, 3);

	SensorPositions = P;   SensorOrientations = O; 
   SensorPositions = 1/1000 * rotat(SensorPositions); SensorOrientations = rotat(SensorOrientations);
	[LocalPos, LocalOrientVec] = trans2local(SensorPositions, SensorOrientations); % Nx3x6
	H = kohlrausch(LocalPos(:,:,Coil));
	whos
	H = rotta(H);
	plot3(x, y, z, 'k-'); disk = fill3(x, y, z, [0.5 0.5 0.5]); alpha(0.5);
	%axis([-maxd maxd -maxd maxd -maxd maxd]);
	axis([-maxd maxd -maxd maxd -maxd maxd]);
	xlabel('Xa [mm]'); ylabel('Ya [mm]'); zlabel('Za [mm]');  grid on;
	title(['Field vectors and resulting amplitude of T' num2str(Coil)]);
	hold on;
	
	quiver3(x, y, z, ox, oy, oz, 0.01*maxd, 'r.'); quiver3(x, y, z, -ox, -oy, -oz, 0.01*maxd, 'b.');
	pidx = find( H(:,3) >= 0 ); nidx = find( H(:,3) < 0 );
	plot3(x(19), y(19), z(19), 'rx');									plot3(x(38), y(38), z(38), 'ro');
	quiver3(x, y, z, H(:,1), H(:,2), H(:,3), 0.01*maxd, 'm.');	quiver3(x, y, z, -H(:,1), -H(:,2), -H(:,3), 0.01*maxd, 'c.');  


	lsc = ['b-'; 'r-'; 'g-'; 'c-'; 'k-'; 'm-'];
	plot3(x, y, Amps(:,Coil)*40, 'g--');  alpha(0.2);
%quiver3(x(pidx), y(pidx), z(pidx), H(pidx,1), H(pidx,2), H(pidx,3), 0.01*maxd, 'm.'); 
%quiver3(x(nidx), y(nidx), z(nidx), H(nidx,1), H(nidx,2), H(nidx,3), 0.01*maxd, 'c.'); 

	view(90, 5);
	hold off;

return

function plotData(alphas, Amps);
	% define line styles for plots 
	lsp = ['b:'; 'r:'; 'g:'; 'c:'; 'k:'; 'm:'];
	lsd = ['b.'; 'r.'; 'g.'; 'c.'; 'k.'; 'm.'];
	lsc = ['b-'; 'r-'; 'g-'; 'c-'; 'k-'; 'm-'];
	a = alphas*180/pi;

	plot(a, Amps(:,1), lsc(1,:), a, Amps(:,2), lsc(2,:), a, Amps(:,3), lsc(3,:),... 
		  a, Amps(:,4), lsc(4,:), a, Amps(:,5), lsc(5,:), a, Amps(:,6), lsc(6,:));
	legend('T1', 'T2', 'T3', 'T4', 'T5', 'T6', 0); 
	xlabel('\alpha [°]'); grid;  set(gca, 'XLim', [0 360], 'XTick', [0 45 90 135 180 225 270 315 360], 'YLim', [0 4]); 
	


	