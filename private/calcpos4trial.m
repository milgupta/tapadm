function [Result, Residuals, Iterations] = calcpos4trial(Result, Residuals, Iterations, Amps, FirstPointNo, LastPointNo, Channel, StartPoint, Options,  DispFunc, DispArg);
% CALCPOS4TRIAL calculates sensor positions and orientations from amplitudes.
%           The function computes a whole trial or a trial segment on a
%           given data structure. The first three arguments Result,
%           Residuals and Iterations are returned after calculating values
%           for the segment between FirstPointNo and LastPointNo and the
%           given Channel. If StartPoint is singleton, the segment is processed 
%           from FirstPointNo to LastPointNo (First may be greater than
%           Last) using the predecessor as start value for the calculation.
%
%           If StartPoint is a vector, it's elements are taken as individual 
%           start values. (In the order the segment is processed.) 
%
%           [Result, Residuals, Iterations] = calcpos4trial(Result, Residuals, 
%                                              Iterations, Amps, FirstPointNo, 
%                                              LastPointNo, Channel, StartPoint, Options);
%
%           The function returns the Nx7 matrix 'Result' with amplitudes for every
%           spatial point (row) and transmitter (column), plus 'rms-value'
%           plus exit condition (column 7): 
%            
%   See Also
%            CALCPOS

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
	Constraint = 0.08; %1.5;						% equals vmax=300mm/s and 10° in 50ms

if (nargin < 11)	
	DispArg = '';
end
if (nargin < 10)	
	DispFunc = [];
end

if (size(StartPoint, 1)>1)
	[Result, Residuals, Iterations] = calc_pos_individual(Result, Residuals, Iterations, Amps, FirstPointNo, LastPointNo, Channel, StartPoint, Options, @amperr,0);
	
else	% if StartPoint is singleton
	
	% process the first (and maybe only) data point														
	[Result, Residuals, Iterations] = calc_single_pos(Result, Residuals, Iterations, ... 
																Amps, FirstPointNo, Channel, StartPoint, Options, [], []);
	if (FirstPointNo == LastPointNo)
		return;
	elseif (FirstPointNo < LastPointNo)	% if the part contains more than one point, identify their ordering
		step = 1;
	else
		step = -1;
	end
	
	for (PointNo = FirstPointNo+step:step:LastPointNo)	% process all the rest
		% if last calculation ended successfully, its result can be used as start point
		if (Iterations(PointNo-step, 1, Channel)  > -1)		
			thisStartPoint = Result(PointNo-step, 1:5, Channel);
		else
			thisStartPoint = StartPoint;
		end
			
		[Result, Residuals, Iterations] = calc_single_pos(Result, Residuals, Iterations, ...
																Amps, PointNo, Channel, thisStartPoint, Options, [], []);
															
		if (rem(PointNo, 100) == 0)     
		   pause(1); 									% a short break to catch a possible <Ctrl-C>
			if (~isempty(DispFunc))
 		 	 	feval(DispFunc, PointNo, DispArg);
			end
		end
	end

end

%------------- subfunctions -----------
% calculate trial with individual startpoints
function [Result, Residuals, Iterations] = calc_pos_individual(Result, Residuals, Iterations, ... 
																Amps, FirstPointNo, LastPointNo, Channel, StartPoints, Options, ErrorFunc, ErrorArg)
	if (FirstPointNo < LastPointNo)
		step = 1;
	else
		step = -1;
	end
	Constraint = 0.08; %1.5;						% equals vmax=300mm/s and 10° in 50ms
	for (PointNo = FirstPointNo+step:step:LastPointNo)					% iterate samples
		[Result, Residuals, Iterations] = calc_single_pos(Result, Residuals, Iterations, ...
																Amps, PointNo, Channel, StartPoints(PointNo, :, Channel), Options, ErrorFunc, Constraint);
		if (rem(PointNo, 5) == 0)     
		   pause(1); 									% a short break to catch a possible <Ctrl-C>
		end
	end

% calculate results for a part of a trial
function [Result, Residuals, Iterations] = calc_segment_pos(Result, Residuals, Iterations, ... 
																Amps, FirstPointNo, LastPointNo, Channel, StartPoint, Options, DispFunc, DispArg)
	% process the first (and maybe only) data point														
	[Result, Residuals, Iterations] = calc_single_pos(Result, Residuals, Iterations, ... 
																Amps, FirstPointNo, Channel, StartPoint, Options, [], []);
	if (FirstPointNo == LastPointNo)
		return;
	elseif (FirstPointNo < LastPointNo)	% if the part contains more than one point, identify their ordering
		step = 1;
	else
		step = -1;
	end
	
	for (PointNo = FirstPointNo+step:step:LastPointNo)	% process all the rest
		% if last calculation ended successfully, its result can be used as start point
		if (Iterations(PointNo-step, 1, Channel)  > -1)		
			thisStartPoint = Result(PointNo-step, 1:5, Channel);
		else
			thisStartPoint = StartPoint;
		end
			
		[Result, Residuals, Iterations] = calc_single_pos(Result, Residuals, Iterations, ...
																Amps, PointNo, Channel, thisStartPoint, Options, [], []);
															
		if (rem(PointNo, 10) == 0)     
		   pause(1); 									% a short break to catch a possible <Ctrl-C>
			if (~isempty(DispFunc))
 		 	 	feval(DispFunc, PointNo, DispArg);
			end
		end
	end

% calculate one single result point from/within an existing Amp/Result array
function [Result, Residuals, Iterations] = calc_single_pos(Result, Residuals, Iterations, ...
																Amps, PointNo, Channel, StartPoint, Options, DummyFunc, Constraint)
	persistent lb ub;
	if (isempty(Constraint))
		lb = []; ub = [];
	else			
		lb = StartPoint - Constraint; ub = StartPoint + Constraint;
	end
	[Result(PointNo, 1:5, Channel), resnorm, Residuals(PointNo, :, Channel), Result(PointNo, 7, Channel), output] = ...
												lsqnonlin(@amperr, StartPoint, lb, ub, Options, Amps(PointNo, :, Channel));	
											
	if (Result(PointNo, 7, Channel) >= 0)
		Result(PointNo, 4:5, Channel) = NormalizeAngles(Result(PointNo, 4:5, Channel));
		Result(PointNo, 6, Channel) = sqrt(resnorm) * 2500/sqrt(6);
		Iterations(PointNo, 1, Channel) = output.iterations;	
	else											% if lsqnonlin  did not converge to a solution,
		Result(PointNo, 1:6, Channel) = [NaN NaN NaN NaN NaN NaN];	% set position to missing data and 
		Iterations(PointNo, 1, Channel) = -output.iterations;			% store number of iterations tried as negative value, i.e. Iteration > 0 means failure
	end
		
function da = amperr(point, amp)
   ca = calcamps(point);
   da = (amp - ca);
	