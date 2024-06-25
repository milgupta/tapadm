function [DParams, Residuals] = derivative(data, Start, options);
% DERIVATIVE calculates a approximation of the derivative of data.
%           The function tries to fit data to a taylor series and
%           returns the accordant parameters to approximate the data
%           along with it's 1'st and 2'nd derivative.
%
%           [DParams, Residuals] = derivative(data, Start, options);
%
%            see trig_approx_func.m

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

[L, N, M] = size(data);
DParams = zeros(16, N, M);
Residuals = zeros(1, N, M);
t = (1/200) * (1 : L);

if (nargin < 3)
	options = optimset('lsqcurvefit');
	options = optimset(options, 'Display', 'off');
	options = optimset(options, 'LargeScale', 'on', 'LevenbergMarquardt', 'off');
	options = optimset(options, 'TolX', 1E-4, 'TolFun', 1E-6, 'MaxFunEvals', 1000);
	if (nargin < 2)
		Start = zeros(size(DParams));				% Start values for approximation
		Start = estimate_start_params(t, data);
	end
end	
	
for i=1:N
	for j=1:M
		[DParams(:, i, j), Residuals(1, i, j)] = calc_single_derivative( data(:, i, j)', t,	options,	Start(:, i, j)' );
	end
end

%------------- subfunctions -----------
function Start = estimate_start_params(data_t, data)    	% Start values for approximation
	[L, N, M] = size(data); Start = zeros(16, N, M);
	for (j = 1 : M)
		for (i = 1 : N)
		AmpRange = max(data(:, i, j)) - min(data(:, i, j)); 	
		Start(1, i, j) = (max(data(:, i, j)) - abs(min(data(:, i, j))))/2;
		Start(2, i, j) = 2*pi/data_t(end);		Start(3, i, j) = -AmpRange/2;				Start(4, i, j) = AmpRange / 2;  		
		Start(5, i, j) = Start(2, i, j) * 2;	Start(6, i, j) = Start(3, i, j)* 0.25; Start(7, i, j) = Start(6, i, j);  		
		Start(8, i, j) = Start(2, i, j) * 2.5;	Start(9, i, j) = Start(6, i, j);			Start(10, i, j) = Start(9, i, j);  		
		Start(11, i, j) = Start(2, i, j) * 3;	Start(12, i, j) = Start(6, i, j);		Start(13, i, j) = Start(9, i, j);  		
		Start(14, i, j) = Start(2, i, j)* 3.5;	Start(15, i, j) = Start(6, i, j);		Start(16, i, j) = Start(9, i, j);  		
		end
	end	
return	

function [x, resnorm] = calc_single_derivative(ydata, t, options, x0)
	[x, resnorm] = lsqcurvefit(@trig_approx_func, x0, t, ydata, [], [], options);
