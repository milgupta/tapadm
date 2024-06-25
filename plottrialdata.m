function PlotTrialData(BasePath, AmpPath, ResPath, TrialNo, ChanNo, Options)
% PLOTTRIALDATA Plots all avaible data for one trial.
%           The function opens three windows and plots both signal amplitudes and
%           position data (in a flat and spatial view). A cursor can be
%           synchronously moved through 
%
%           PlotTrialData(BasePath, AmpPath, ResPath, TrialNo, ChanNo, Options)

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

% check function arguments and evaluate Options
errormsg = nargchk(5, 6, nargin);
if (~isempty(errormsg))
	error(errormsg);
end
if (nargin < 6)								
	Options = '';
end

% Load resultdata for the trial

T = loadtrial(fullfile(BasePath, ResPath, trialfile(TrialNo)),... 
              fullfile(BasePath, ResPath, trialfile(TrialNo, 'resi')));

figure(1)										% Amplitude plot
[amp_fh, amp_sfh] = plotamps(fullfile(BasePath, AmpPath, trialfile(TrialNo)), ChanNo, gcf);
set(gca, 'UserData', amp_sfh, 'Tag', 'AmpCursorCallback');

figure(3)										% spatial plot
point_idx = 1:T.NumPoints;		
Result = T.Result; Result(:,7,:) = T.Iterations;
Styles = ['r-'; 'g-'; 'b-'; 'm-'; 'y-'; 'c-'; 'r:'; 'g:'; 'b:'; 'm:'; 'y:'; 'c:'];
MinP = min(min(min(Result(point_idx, 1:3, ChanNo))));
MaxP = max(max(max(Result(point_idx, 1:3, ChanNo))));
MeanP = mean(Result(point_idx, 1:3, ChanNo));
PosRange = (MaxP-MinP);

drapefig(PosRange); fh = gcf;
title(['Trial ' num2str(TrialNo)]); 
plotsc(Result, ChanNo, point_idx, fh, 'g-');
set(gca, 'UserData', Result(:,:, ChanNo), 'Tag', 'AxisWithResultData');

L = calc_3d_cursor_line(Result(1, 1:5, ChanNo), 1); %calculate cursor position for the 1st point
camtarget(mean(L)); 							% set the camera target to this point
line(L(:,1), L(:,2), L(:,3),...			% create the spatial cursor
	'LineWidth', 4, 'Color', [.8 .8 .8], 'EraseMode', 'xor', 'Tag', '3DCursor');
sfh = @set_cursors;							% define a callback function for the cursor

figure(2)										% flat result plot
plotsc(Result, ChanNo, point_idx, gcf, 'flat', sfh);			  
			  

%------------- subfunctions -----------
function L = calc_3d_cursor_line(d, scale)
	[x, y, z] = sph2cart(d(4)*pi/180, d(5)*pi/180, scale); ov = [x y z];
	p1 = d(1:3) - ov;	p2 = d(1:3) + ov;
	L = [p1; p2];
return

% set cursors in all subplots (can be called from outside)  	
function set_cursors(idx, cf)
	da = findobj('Tag', 'AxisWithResultData'); data = get(da, 'UserData');
	L = calc_3d_cursor_line(data(idx, 1:5), 1);
	c = findobj('Tag', '3DCursor');	set(c, 'XData', L(:,1)); set(c, 'YData', L(:,2)); set(c, 'ZData', L(:,3));
	da = findobj('Tag', 'AmpCursorCallback'); amp_sfh = get(da, 'UserData');
	if (~isempty(amp_sfh))						% if defined, call the cursor callback function of the amplitude figure
		feval(amp_sfh, idx, gcf);
	end
return	
