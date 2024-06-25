function smartcal(CalDataFile, ChanIdx, ResultFile, SmartcalOptions, AlphaPreOffsetDeg);
% SMARTCAL  Calculates calibration factors from a calibration data file.
%           The function reads calibration data from a MAT-file
%           (CalDataFile) and operates on the channels named in ChanIdx. 
%           It then writes the computed calibration data to ResultFile.
%
%           smartcal(CalDataFile, ChanIdx, ResultFile, SmartcalOptions, AlphaPreOffset);
%           
%           AlphaPreOffset defines an alpha0 start value (given in
%           degrees). The Circal supports a pre-offset for alpha which 
%           can be found in the file 'circal.ini' on the IDA-PC. (LogicZero)
%           Multiply this value by (360/800), to get the appropriate AlphaPreOffset.   
%
%           see cal_prepdata.m

% smartcal('E:\2006\cal_ap_ema3\calibration.mat', 2, 'smartcaltest', '-p -h', 0);

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
COLORS = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 0 0 0; 1 0 1];
set(0, 'DefaultAxesColorOrder', COLORS);
disp(['        SMARTCAL ' tapadversion ]);	

if (nargin < 2)
	ChanIdx = 1:12;
end

if (nargin < 3)
	ResultFile = 'SmartCalData.mat';
end

if (nargin < 4)								% evaluate SmartcalOptions
	SmartcalOptions = '';
end

use_history = 	~isempty(findstr('-h', SmartcalOptions));
generate_plots = 	~isempty(findstr('-p', SmartcalOptions));

if ((length(ChanIdx) > 1) & generate_plots)
	warn('Can not plot multi-channel calibration data.');
	generate_plots = logical(0);
end

if (nargin < 5)
	AlphaPreOffset = 345 * pi / 400;		% The offeset (encoder steps) as stated in 'circal.ini' on our IDA-PC
else
	AlphaPreOffset = AlphaPreOffsetDeg * pi / 180;
end

%--------------------------------------

alpha_options = optimset('fminsearch');			% option set for the first pass (alpha-optimization)
if (generate_plots)
	alpha_options = optimset(alpha_options, 'Display', 'iter', 'Diagnostics', 'on');
else
	alpha_options = optimset(alpha_options, 'Display', 'notify');
end
alpha_options = optimset(alpha_options, 'MaxIter', 50, 'MaxFunEvals', 10000, 'TolX', 1E-8, 'TolFun', 5*pi/180);

tricomp_options = optimset('fminsearch');			% option set for the tricomp-optimization
if (generate_plots)
	tricomp_options = optimset(tricomp_options, 'Display', 'iter', 'Diagnostics', 'on');
else
	tricomp_options = optimset(tricomp_options, 'Display', 'notify');
end
%tricomp_options = optimset(tricomp_options, 'MaxIter', 4, 'MaxFunEvals', 10000, 'TolX', 1E-8, 'TolFun', 1E-10);
tricomp_options = optimset(tricomp_options, 'MaxIter', 4, 'MaxFunEvals', 1E6, 'TolX', 1E-9, 'TolFun', 1E-10); % 100 = 1/4h

%-------------------------------------- define parameters of calibration and their  initial values

alpha0 = 0;								 		% disk starting angle (offset), measured counter-clockwise around Z
DiskPosition = [0 0 0];						% spatial position of the calibration disk
DiskAxis = [0 pi/2];							% unit vector, points direction of disk's pivot  
DiskR	= 78;										% radius of the disk

SensorOrientations = repmat([pi/4 0], [1 1 12]); % individual sensor orientations (phi, theta) in rad for 12 channels
SensorPositions = zeros(1, 3, 12); 		% individual sensor positions relative to the 'Magazin' 
CalibrationFactors = ones(1, 6, 12);	% expected signal at the center when sensor is parallel transmitter (TAPAD calculation) 
ReCalibrationFactors = ones(1, 6, 12);	% to perform a supplementary calibration, i.e. 'correct' Carstens Calibration Factors

if (use_history)
	if (exist(ResultFile, 'file'))
		f = load(ResultFile);
		fn = fieldnames(FData);
		if (~isempty(strmatch('alpha0', fn, 'exact'))) alpha0 = double(FData.alpha0); end
		if (~isempty(strmatch('DiskR', fn, 'exact'))) DiskR = double(FData.DiskR); end
		if (~isempty(strmatch('DiskPosition', fn, 'exact'))) DiskPosition = double(FData.DiskPosition); end
		if (~isempty(strmatch('DiskAxis', fn, 'exact'))) DiskAxis = double(FData.DiskAxis); end
		if (~isempty(strmatch('SensorOrientations', fn, 'exact'))) SensorOrientations = double(FData.SensorOrientations); end
		if (~isempty(strmatch('CalibrationFactors', fn, 'exact'))) CalibrationFactors = double(FData.CalibrationFactors); end
		if (~isempty(strmatch('ReCalibrationFactors', fn, 'exact'))) ReCalibrationFactors = double(FData.ReCalibrationFactors); end
		clear f;
	else
	disp(['nonexistig file: ''' ResultFile '''']);
	end
end

%-------------------------------------- Read the AG500 (Circal) data

if (exist(CalDataFile, 'file'))
	calfile = load(CalDataFile); 
	disp(['loaded ' calfile.VERSIONINFO]);
else
	error(['nonexistig file: ''' CalDataFile '''']);
end
CircalCalibrationFactors = ones(1, 6, 12);	% calibration factors from the "Circal"-Program

AlphaOffsets = zeros(4, 6); AlphaWhereNull = zeros(6, 2); AlphaWhereMax = zeros(1, 6);

%-------------------------------------- Iterate channels...

for i_chan = 1:length(ChanIdx)
   chan = ChanIdx(i_chan);
	
	% apply the alpha pre-offset and wrap negative angles.
	CircalCalibrationFactors(1, :, chan) = calfile.CalData(chan).OriginalCalibrationFactor;
	alphas = calfile.CalData(chan).alphas; 	Amps = calfile.CalData(chan).Amps; 
	NumPoints = calfile.CalData(chan).NumPoints;	alpha_resolution = alphas(2)-alphas(1); alphas = alphas'; 
	alphas = alphas - AlphaPreOffset; 
	alphas = [alphas(find(alphas>=0)) alphas(find(alphas<0))+2*pi]; % (-a) = 360° + (-a)

	if (generate_plots)
		% Now find the circal angles where the measured amps are maximal and calculate
		% amplitudes for these angles. The ratio between measured and calculated amps
		% gives a first estimation of the calibration factors. (Measured amplitudes are 
		% in 'ADC-Digits', while calculated amps are monic to one)
		[AlphaWhereNull, AlphaWhereMax] = estimatezeros(alphas, Amps,  0.0005*pi/180);
		CAmps = CalcCalAmps(alphas); 	[CAlphaWhereNull, CAlphaWhereMax] = estimatezeros(alphas, CAmps,  0.0005*pi/180);
		Similarity = zeros(1, 6);
		for (coil=1:6)
			Similarity(coil) = tricomp([AlphaWhereNull{coil}' AlphaWhereMax(coil)], [CAlphaWhereNull{coil}' CAlphaWhereMax(coil)]);
		end

		[amps_max_amp, idx] = max(abs(Amps));		amps_max_alpha = alphas(idx);
		[camps_max_amp, idx] = max(abs(CAmps));	camps_max_alpha = alphas(idx); clear idx;
		CalibrationFactors(1, :, chan) = amps_max_amp ./ camps_max_amp;

		disp(['channel ' num2str(chan) ' initial setting']);	
		disp(['Similarity [%]=                  ' num2str(Similarity*100, '%6.1f ') ]);	
		disp(['CircalCalibrationFactors=   ' num2str(CircalCalibrationFactors(:, :, chan), '%6.0f ')]);	
		disp(['estimated CalibrationFactors= ' num2str(CalibrationFactors(:, :, chan), '%6.0f ') ]);	
		disp(['difference =                  ' num2str(CircalCalibrationFactors(:, :, chan)-CalibrationFactors(:, :, chan), '%6.0f ') ]);	
	
		n = sort([AlphaWhereNull{1} AlphaWhereNull{2} AlphaWhereNull{3} AlphaWhereNull{4} AlphaWhereNull{5} AlphaWhereNull{6}]) * 180/pi; 
		CAmps = CalcCalAmps(alphas, alpha0, DiskAxis, DiskR, DiskPosition, SensorOrientations(:, :, chan), SensorPositions(:, :, chan), CalibrationFactors(:, :, chan)); 	
		DX = CAmps - Amps;
		figure(1);
		subplot(2,1,1);
		plot2SignalSets(alphas, Amps, CAmps , 'measured amp', 'calculated Amps', ['SmartCal ' tapadversion ' ' CalDataFile ' channel ' num2str(chan)], 'initial settimg');
		hold on; plot(n(1), 0, 'bd', n(2), 0, 'bd', n(3), 0, 'rd', n(4), 0, 'rd', n(5), 0, 'gd', n(6), 0, 'gd');  hold off;
		subplot(2,1,2);
		plot2SignalSets(alphas, Amps, DX , 'measured amp', 'dx', '', '', '-d');
		disp('Now performing the optimization of the calibration parameters. First estimate alpha0...');
	end		

	% Perform optimization of the calibration parameters (see bottom of this script)
	
	try											% first estimate alpha0
		[alpha0, fval, exitflag]  = fminsearch(@OFunc_SimpleAlphaParam, alpha0, alpha_options, alphas, Amps); 
	catch
		exitflag = -1;
	end

	arg = [alpha0 DiskR DiskPosition(3) SensorOrientations(:, :, chan)];
	dispArg(arg);
	tic; 
	try
		[param, fval, exitflag]  = fminsearch(@OFunc_TricompClassicParams, arg, tricomp_options, AlphaWhereNull, AlphaWhereMax); 
	catch
		exitflag = -1;
	end
	t = toc;
	if (exitflag >= 0)
		alpha0 = param(1);					 					% disk starting angle (offset), measured counter-clockwise around Z
		DiskR = param(2);
		DiskPosition = [0 0 param(3)];						% spatial position of the calibration disk
		SensorOrientations(:, :, chan) = param(4:5);
		if (exitflag == 0)
			warning('optimization abort: maximum number of function evaluations or iterations was exceeded.');
		end
		if (generate_plots)
			disp(['optimization took ' num2str(round(t/60)) ' min and ' num2str(rem(t, 60)) ' s']);
		
			CAmps = CalcCalAmps(alphas, alpha0, DiskAxis, DiskR, DiskPosition, SensorOrientations(:, :, chan));
 			[CAlphaWhereNull, CAlphaWhereMax] = estimatezeros(alphas, CAmps,  0.0005*pi/180);
			Similarity = zeros(1, 6);
			for (coil=1:6)
				Similarity(coil) = tricomp([AlphaWhereNull{coil}' AlphaWhereMax(coil)], [CAlphaWhereNull{coil}' CAlphaWhereMax(coil)]);
			end
	
			[amps_max_amp, idx] = max(abs(Amps));		amps_max_alpha = alphas(idx);
			[camps_max_amp, idx] = max(abs(CAmps));	camps_max_alpha = alphas(idx); clear idx;
			CalibrationFactors(1, :, chan) = amps_max_amp ./ camps_max_amp;
			dispArg(param);
			disp(['Similarity after optimization[%]=                  ' num2str(Similarity*100, '%6.1f ') ]);	
			disp(['CircalCalibrationFactors=   ' num2str(CircalCalibrationFactors(:, :, chan), '%6.0f ')]);	
			disp(['estimated CalibrationFactors= ' num2str(CalibrationFactors(:, :, chan), '%6.0f ') ]);	
			disp(['difference =                  ' num2str(CircalCalibrationFactors(:, :, chan)-CalibrationFactors(:, :, chan), '%6.0f ') ]);	
			
			n = sort([AlphaWhereNull{1} AlphaWhereNull{2} AlphaWhereNull{3} AlphaWhereNull{4} AlphaWhereNull{5} AlphaWhereNull{6}]) * 180/pi; 
			CAmps = CalcCalAmps(alphas, alpha0, DiskAxis, DiskR, DiskPosition, SensorOrientations(:, :, chan), SensorPositions(:, :, chan), CalibrationFactors(:, :, chan)); 	
			DX = CAmps - Amps;
			figure(2); 
			subplot(2,1,1);
			plot2SignalSets(alphas, Amps, CAmps , 'measured amp', 'calculated Amps', ['SmartCal ' tapadversion ' ' CalDataFile ' channel ' num2str(chan)], 'TriComp algorithm');
			hold on; plot(n(1), 0, 'bx', n(2), 0, 'bx', n(3), 0, 'rx', n(4), 0, 'rx', n(5), 0, 'gx', n(6), 0, 'gx');  hold off;
			subplot(2,1,2);
			plot2SignalSets(alphas, Amps, DX , 'measured amp', 'dx', '', '', '-d');			
		end
	else
		warning('fminsearch found no solution.');
	end

	
	% remove all outlier LCFs and all where measured signals are small 

	CAmps = CalcCalAmps(alphas, alpha0, DiskAxis, DiskR, DiskPosition, SensorOrientations(:, :, chan));
	LCF = Amps ./ (CAmps+eps);
	N = size(LCF,1);
	LCF(abs(LCF - repmat(mean(LCF), N, 1)) > repmat(2*std(LCF), N, 1)) = NaN;
	LCF(abs(Amps) < 1000) = NaN;
	CalibrationFactors(1, :, chan) = nanmedian(LCF);
	if (generate_plots)
		figure(3); plot(LCF);
		disp(['CircalCalibrationFactors=   ' num2str(CircalCalibrationFactors(:, :, chan), '%6.0f ')]);	
		disp(['smartcal CalibrationFactors = ' num2str(CalibrationFactors(:, :, chan), '%6.0f ') ]);	
		disp(['difference =                  ' num2str(CircalCalibrationFactors(:, :, chan)-CalibrationFactors(:, :, chan), '%6.0f ') ]);	
	end
	
end												% of iterate channels 
	
% Hassle with calibration factors:
% Carstens established relative signal amplitudes, where the measured
% signal for a sensor at position [0 0 0] and orientation parallel to the
% magnetic field vector is set to 1. To use a different calibration
% procedure (and/or different field model), we have to reverse that first
% before we can apply different calibration factors. 
% This can be performed by using 'ReCalibrationFactors' as new calibration
% factors, when operating on relative amplitude data.

ReCalibrationFactors = CalibrationFactors ./ CircalCalibrationFactors;	% divide measured Data by

disp(['saving  calibration data to ''' ResultFile '''']);
v_tmp = version; Matlab_Is_Newer_Than_6 = str2num(v_tmp(1)) > 6; clear v_tmp;
if Matlab_Is_Newer_Than_6	
	save(ResultFile, 'alpha0', 'DiskR', 'DiskPosition', 'DiskAxis', 'SensorOrientations', 'SensorPositions', 'CalibrationFactors', 'ReCalibrationFactors', '-v6');
else	
	save(ResultFile, 'alpha0', 'DiskR', 'DiskPosition', 'DiskAxis', 'SensorOrientations', 'SensorPositions', 'CalibrationFactors', 'ReCalibrationFactors');
end

if (generate_plots)
	disp(['alpha0= ' num2str(alpha0*180/pi, '%6.2f') '° DiskRadius= ' num2str(DiskR, '%6.2f') 'mm DiskPosition= ['...
			num2str(DiskPosition, '%6.2f') ']mm elevation of DiskAxis= ' num2str(DiskAxis(2)*180/pi, '%6.2f') '°']);
	disp(['distance 1-2= ' num2str(norm(SensorPositions(:, :, 1) - SensorPositions(:, :, 2)), '%6.2f') ' mm '...
		' distance 2-3= ' num2str(norm(SensorPositions(:, :, 2) - SensorPositions(:, :, 3)), '%6.2f') ' mm '...
		' distance 3-4= ' num2str(norm(SensorPositions(:, :, 3) - SensorPositions(:, :, 4)), '%6.2f') ' mm ']);
end
for (chan=1:12)	
	disp(['channel ' num2str(chan) ': SensorPosition= [' num2str(SensorPositions(:, :, chan), '%6.2f') ']mm SensorOrientation: phi= '...
			num2str(SensorOrientations(:, 1, chan)*180/pi, '%6.2f') '° theta= ' num2str(SensorOrientations(:, 2, chan)*180/pi, '%6.2f')...
			'° CalibrationFactors= ' num2str(CalibrationFactors(:, :, chan), '%6.0f ')]);
end						

disp('Difference to Carstens-CalibrationFactors')
disp(num2str(shiftdim(CircalCalibrationFactors(:,:,1:4)-CalibrationFactors(:,:,1:4))', '%6.2f '))
%keyboard
return 


%------------- subfunctions -----------
	
function d = AngleDiff(a, b)				% returns the difference of two vectors (d = a - b). 
	% The function assumes that a and b are containing angles in [rad] and makes shure that 
	% all componets of the result are in the range (-pi .. pi) 
	
	d =  a - b; 
	idx = (d <= pi); d(idx) = d(idx) + 2*pi;
	idx = (d > pi); d(idx) = d(idx) - 2*pi;
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
	

function dispArg(arg)
	na = NormalizeAngles([arg(4)*180/pi arg(5)*180/pi]);
	disp(['DiskR: ' num2str(arg(2), '%6.2f') 'mm Z: ' num2str(arg(3), '%6.2f/') 'mm']);
	disp(['alpha0 is ' num2str(arg(1)*180/pi, '%6.2f') '° Sensor at disk is ' num2str(na(1), '%6.2f') '° horizontally rotated and '...
			num2str(na(2), '%6.2f') '° vertically elevated']);


%---------------------------------------------------------------------
% C a l i b r a t i o n   P a r a m e t e r s
%---------------------------------------------------------------------
% alpha0						disk starting angle (offset), measured counter-clockwise around Z
% DiskPosition				spatial position of the calibration disk
% DiskAxis 					unit vector, points direction of disk's pivot  
% DiskR						radius of the disk
% SensorOrientations		individual sensor orientations (phi, theta) in rad for 12 channels
% SensorPositions 		individual sensor positions relative to the 'Magazin' 
% CalibrationFactors		expected signal at the center when sensor is parallel transmitter

% calculate signal amplitudes as function of alpha, the rotation angle of the calibration disk
function CAmps = CalcCalAmps(alphas, alpha0, DiskAxis, DiskR, DiskPosition, SensorOrientation, SensorPosition, CalibrationFactors) 
	if (nargin < 7)
		SensorPosition = zeros(1, 3); 		% individual sensor positions relative to the 'Magazin' 
	end
	if (nargin < 6)
		SensorOrientation = [pi/4 0]; 		% individual sensor orientations (phi, theta) in rad for 12 channels
	end
	if (nargin < 5)
		DiskPosition = [0 0 1];					% spatial position of the calibration disk	
	end
	if (nargin < 4)
		DiskR	= 80;									% radius of the disk
	end
	if (nargin < 3)
		DiskAxis = [0 pi/2];						% unit vector, points direction of disk's pivot
	end
	if (nargin < 2)
		alpha0 = 0;				 					% disk starting angle (offset), measured counter-clockwise around Z
	end

	[ox, oy, oz] = sph2cart(SensorOrientation(1), SensorOrientation(2), 1);
	[axo, aoy, aoz] = sph2cart(DiskAxis(1), DiskAxis(2), 1); DiskAxisVec = [axo aoy aoz];

	P = CalcPos4Cal(alphas'+alpha0, DiskAxisVec, DiskR, DiskPosition+SensorPosition);		 % generate calibration positions and orientations
	O = CalcOrient4Cal(alphas+alpha0, DiskAxisVec, [ox oy oz]);

	[phi, theta] = cart2sph(O(:,1), O(:,2), O(:,3));  
	if (nargin >7)
		CAmps = calcamps([P phi*180/pi theta*180/pi], CalibrationFactors);
	else
		CAmps = calcamps([P phi*180/pi theta*180/pi]);
	end
return

% generate position (x/y/z) data as function of alpha, the rotation angle of the calibration disk  
function P = CalcPos4Cal(alphas, DiskAxis, DiskR, PositionOffset)
	NumPoints = size(alphas, 1);	P = zeros(NumPoints, 3);	O = zeros(NumPoints, 3);
	s = sin(alphas'); c = cos(alphas'); t = 1 - c;

	pivot = cross(DiskAxis, [0 0 1]);	x = pivot(1); y = pivot(2); z= pivot(3);
	s = norm(pivot); c = cos(asin(s));	t = 1 - c;
	
	R = [t*x^2 + c			t*x*y + s*z			t*x*z - s*y;...
		  t*x*y - s*z		t*y^2 + c			t*y*z + s*x;...
  		  t*x*z + s*y		t*y*z - s*x			t*z^2 + c];  

	P(:,1) = DiskR .* cos(alphas);
	P(:,2) = DiskR .* sin(alphas);
	
	for (i=1:NumPoints)
		P(i, :) = P(i, :) * R + PositionOffset; 
	end	
return

% generate sensor orientation (phi/theta) data as function of alpha, the rotation angle of the calibration disk  
function O = CalcOrient4Cal(alphas, DiskAxis,  StartOrientation)
	NumPoints = size(alphas, 2);
	pivot = cross(DiskAxis, [0 0 1]);	x = pivot(1); y = pivot(2); z= pivot(3);
	s = norm(pivot); c = cos(asin(s));	t = 1 - c;
	
	R = [t*x^2 + c			t*x*y + s*z			t*x*z - s*y;...
		  t*x*y - s*z		t*y^2 + c			t*y*z + s*x;...
  		  t*x*z + s*y		t*y*z - s*x			t*z^2 + c];  

	ox = StartOrientation(1); oy = StartOrientation(2); oz = StartOrientation(3); 
	O = oz * ones(NumPoints, 3);
	O(:,1) = ox * cos(alphas') - oy * sin(alphas');
	O(:,2) = ox * sin(alphas') + oy * cos(alphas');

	for (i=1:NumPoints)
		O(i, :) = O(i, :) * R; 
	end	
return

%------------------- Objective Functions -----------------------------
% Not all of these objective functions are actually used in this function, but
% they may be a starting point for own experiments. 

% Objective for a single channel
function f = OFunc_SingleChanParams(arg, alphas, Amps, alpha0, DiskAxis, DiskR, DiskPosition)
	CalibrationFactors = arg(1:6);	% calibration parameters
	SensorPosition = arg(7:9);
	SensorOrientation = arg(10:11);
	
	CAmps = CalcCalAmps(alphas, alpha0, DiskAxis, DiskR, DiskPosition, SensorOrientation, SensorPosition, CalibrationFactors);
	D = Amps - CAmps;							% The Amplitude Error, i.e. measured minus calculated amps for all transmitters 
	% For lsqnonlin we need a error vector, so now we have to merge the error-signals of the 6 transmitters. This makes sense anyway,
	% because the sensor position and orientation obviously is all the same for every transmitter signal, just the calibration 
	% factors are different. Several ways are thinkable to compute the error norm, quadratic or not. Since a quadratic norm like                               
	% RMS is very sensitive to outliers, it may distort the whole parameter-set to fit one transmitter which has a large error.
	% Also, RMS assumes errors are random and Gaussian distributed. Since this objective uses only some of the calibration parameters,
	% the errors are most probably not normal distributed. 
	f = median(D, 2);							
return

% Objective for global parameters, works on 4 combined channels
function f = OFunc_MultiChanParams(arg, CalData, channels, CalibrationFactors, SensorPositions, SensorOrientations)
	alpha0 = arg(1);			% calibration parameters
	DiskAxis = arg(2:3);
	DiskR = arg(4);
	DiskPosition = arg(5:7);
	D = zeros(100, 6, length(channels));
	for (i=1: length(channels))
		chan = channels(i);
		alphas = CalData(chan).alphas; 	Amps = CalData(chan).Amps; 
		CAmps = CalcCalAmps(alphas', alpha0, DiskAxis, DiskR, DiskPosition, SensorOrientations(:, :, chan), SensorPositions(:, :, chan), CalibrationFactors(:, :, chan)); 
		D(:, :, i) = (Amps - CAmps);							% The Amplitude Error, i.e. measured minus calculated amps for all transmitters 
	end
%	E = sum(D.^2, 3).^0.5;
	E = median(D, 2);
	f = median(E, 3);							
%	f = sum(E.^2, 3).^0.5; %median(E, 2);							
return

% Objective for all 51 parameters, works on 4 combined channels
function f = OFunc_4ChanAllParams(arg, CalData)
	alpha0 = arg(1);
	DiskAxis = arg(2:3);
	DiskR = arg(4);
	DiskPosition = arg(5:7);
	CalibrationFactors = reshape(arg(8:31), [1 6 4]);
	SensorPositions = reshape(arg(32:43),  [1 3 4]);
	SensorOrientations = reshape(arg(44:51),  [1 2 4]);
	
	D = zeros(100, 6, 4);
	for (chan=1:4)
		alphas = CalData(chan).alphas; 	Amps = CalData(chan).Amps; 
		CAmps = CalcCalAmps(alphas', alpha0, DiskAxis, DiskR, DiskPosition, SensorOrientations(:, :, chan), SensorPositions(:, :, chan), CalibrationFactors(:, :, chan)); 
		D(:, :, chan) = (Amps - CAmps);							% The Amplitude Error, i.e. measured minus calculated amps for all transmitters 
	end
	E = sum(D.^2, 3).^0.5;	% merge channels
%	E = median(D, 3);
	f = median(E, 2);							
%	f = sum(E.^2, 3).^0.5; %median(E, 2);							
return
%arg = [alpha0 DiskAxis DiskR DiskPosition reshape(CalibrationFactors(:,:,1:4), 1, 24) reshape(SensorPositions(:,:,1:4), 1, 12) reshape(SensorOrientations(:,:,1:4), 1, 8)];
%tic; [param, resnorm, residual, exitflag]  = lsqnonlin(@OFunc_4ChanAllParams, arg, [], [], options, calfile.CalData); t = toc
%if (exitflag >= 0)
%	alpha0 = param(1);
%	DiskAxis = param(2:3);
%	DiskR = param(4);
%	DiskPosition = param(5:7);
%	CalibrationFactors(:, :, 1:4) = reshape(param(8:31), [1 6 4]);
%	SensorPositions(:, :, 1:4) = reshape(param(32:43),  [1 3 4]);
%	SensorOrientations(:, :, 1:4) = reshape(param(44:51),  [1 2 4]);
%	if (exitflag == 0)
%		warning('optimization abort: maximum number of function evaluations or iterations was exceeded.');
%	end
%	disp(['optimization took ' num2str(round(t/60)) ' min and ' num2str(rem(t, 60)) ' s residual norm is ' num2str(resnorm)]);
%else
%	warning('Calibration failed due to unsuccessful lsqnonlin-optimization.');
%end


% Objective for a single channel
function f = OFunc_ClassicParams(arg, alphas, Amps)
	DiskAxis = [0 pi/2];	SensorPosition = [0 0 0];

	CalibrationFactors = arg(1:6);	% calibration parameters
	alpha0 = arg(7);					 					% disk starting angle (offset), measured counter-clockwise around Z
	DiskR = arg(8);
	DiskPosition = [0 0 arg(9)];						% spatial position of the calibration disk
	SensorOrientation = arg(10:11);
	

	CAmps = CalcCalAmps(alphas, alpha0, DiskAxis, DiskR, DiskPosition, SensorOrientation, SensorPosition, CalibrationFactors);
	D = Amps - CAmps;							% The Amplitude Error, i.e. measured minus calculated amps for all transmitters 
	% For lsqnonlin we need a error vector, so now we have to merge the error-signals of the 6 transmitters. This makes sense anyway,
	% because the sensor position and orientation obviously is all the same for every transmitter signal, just the calibration 
	% factors are different. Several ways are thinkable to compute the error norm, quadratic or not. Since a quadratic norm like                               
	% RMS is very sensitive to outliers, it may distort the whole parameter-set to fit one transmitter which has a large error.
	% Also, RMS assumes errors are random and Gaussian distributed. Since this objective uses only some of the calibration parameters,
	% the errors are most probably not normal distributed. 
	f = sum(D.^2, 2).^0.5;	% merge channels
%	f = median(D, 2);							
return


% Objective for 'Tricomp'-method, i.e.  channel
% arg = [alpha0 DiskR DiskPosition(3) SensorOrientations(:, :, chan)];
function f = OFunc_TricompClassicParams(arg, AlphaWhereNull, AlphaWhereMax)
	DiskAxis = [0 pi/2];	SensorPosition = [0 0 0];

	alpha0 = arg(1);					 					% disk starting angle (offset), measured counter-clockwise around Z
	DiskR = arg(2);
	DiskPosition = [0 0 arg(3)];						% spatial position of the calibration disk
	SensorOrientation = arg(4:5);

	da = 3.6 *pi/180;								% desired dalpha 
	alphas = 0 : da : 2*pi-da;

	CAmps = CalcCalAmps(alphas, alpha0, DiskAxis, DiskR, DiskPosition, SensorOrientation, SensorPosition);
	[CAlphaWhereNull, CAlphaWhereMax] = estimatezeros(alphas, CAmps,  0.0005*pi/180);
	Strangeness = zeros(1, 6);
	for (coil=1:6)
		Strangeness(coil) = 1 - tricomp([AlphaWhereNull{coil}' AlphaWhereMax(coil)], [CAlphaWhereNull{coil}' CAlphaWhereMax(coil)]);
	end
	f = sum(Strangeness.^2, 2).^0.5;	% merge channels
%	f = median(D, 2);							
return

% Objective for 'Tricomp'-method, i.e.  channel
% arg = [alpha0 DiskR DiskPosition(3) SensorOrientations(:, :, chan)];
function f = OFunc_TricompSensorParams(arg, AlphaWhereNull, AlphaWhereMax, alpha0, DiskAxis, DiskR, DiskPosition)
	SensorPosition = arg(1:3);
	SensorOrientation = arg(4:5);

	da = 3.6 *pi/180;								% desired dalpha 
	alphas = 0 : da : 2*pi-da;

	CAmps = CalcCalAmps(alphas, alpha0, DiskAxis, DiskR, DiskPosition, SensorOrientation, SensorPosition);
	[CAlphaWhereNull, CAlphaWhereMax] = estimatezeros(alphas, CAmps,  0.0005*pi/180);
	Strangeness = zeros(1, 6);
	for (coil=1:6)
		Strangeness(coil) = 1 - tricomp([AlphaWhereNull{coil}' AlphaWhereMax(coil)], [CAlphaWhereNull{coil}' CAlphaWhereMax(coil)]);
	end
	f = sum(Strangeness.^2, 2).^0.5;	% merge channels
return

% Objective for estimation of alpha0, uses only distance of the maxima
function f = OFunc_SimpleAlphaParam(alpha0, alphas, Amps)
	CAmps = CalcCalAmps(alphas, alpha0);
	[amps_max_amp, idx] = max(abs(Amps));		amps_max_alpha = alphas(idx);
	[camps_max_amp, idx] = max(abs(CAmps));	camps_max_alpha = alphas(idx); clear idx;
	f = norm(AngleDiff(amps_max_alpha, camps_max_alpha));
return





