function plot2SignalSets(alphas, Data1, Data2, Label1, Label2, PlotTitle, PlotSubTitle, plot2sOptions);
% PLOT2SIGNALSETS  plots two sets of 6 signals.
%           The function reads calibration data from a MAT-file
%           (CalDataFile) and operates on the channels named in ChanIdx. 
%           It then writes the computed calibration data to ResultFile.
%
%           plot2SignalSets(alphas, Data1, Data2, Label1, Label2, PlotTitle, PlotSubTitle, plot2sOptions);
%           
%           see cal_prepdata.m

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


% compare Measured and calculated amplitudes	
	% define line styles for plots 
	lsp = ['b:'; 'r:'; 'g:'; 'c:'; 'k:'; 'm:'];
	lsd = ['b.'; 'r.'; 'g.'; 'c.'; 'k.'; 'm.'];
	lsc = ['b-'; 'r-'; 'g-'; 'c-'; 'k-'; 'm-'];
	
	if (nargin < 8)								% evaluate SmartcalOptions
		plot2sOptions = '';
	end
	dualScale = 			~isempty(findstr('-d', plot2sOptions));	
	
	a = alphas*180/pi;
	
	if (dualScale)
		plot(a, Data1(:,1), lsd(1,:), a, Data1(:,2), lsd(2,:), a, Data1(:,3), lsd(3,:),... 
			  a, Data1(:,4), lsd(4,:), a, Data1(:,5), lsd(5,:), a, Data1(:,6), lsd(6,:));
		legend([Label1 ' T1'], [Label1 ' T2'], [Label1 ' T3'], [Label1 ' T4'], [Label1 ' T5'], [Label1 ' T6'], 3); 
		xlabel('\alpha [°]'); ylabel(Label1); grid; li = get(gca, 'YLim'); li = max(abs(li));
		set(gca, 'XLim', [0 360], 'XTick', [0 45 90 135 180 225 270 315 360], 'YLim', [-li +li]); 
	
		h1 = gca; h2 = axes('Position', get(h1, 'Position'));
		plot(a, Data2(:,1), lsc(1,:), a, Data2(:,2), lsc(2,:), a, Data2(:,3), lsc(3,:),... 
			  a, Data2(:,4), lsc(4,:), a, Data2(:,5), lsc(5,:), a, Data2(:,6), lsc(6,:));  ylabel([Label2 '(solid lines)']);
		legend([Label2 ' T1'], [Label2 ' T2'], [Label2 ' T3'], [Label2 ' T4'], [Label2 ' T5'], [Label2 ' T6'], 1); 
		set(gca, 'XLim', [0 360], 'XTick', [0 45 90 135 180 225 270 315 360]); li = get(gca, 'YLim'); li = max(abs(li));
		set(h2, 'YAxisLocation', 'right', 'XTickLabel', [], 'XGrid', 'off', 'Color', 'none', 'YColor', 'k', 'YLim', [-li +li]);
	else
		plot(a, Data1(:,1), lsd(1,:), a, Data1(:,2), lsd(2,:), a, Data1(:,3), lsd(3,:),... 
			  a, Data1(:,4), lsd(4,:), a, Data1(:,5), lsd(5,:), a, Data1(:,6), lsd(6,:),...
			  a, Data2(:,1), lsc(1,:), a, Data2(:,2), lsc(2,:), a, Data2(:,3), lsc(3,:),... 
			  a, Data2(:,4), lsc(4,:), a, Data2(:,5), lsc(5,:), a, Data2(:,6), lsc(6,:));  
		legend('measured T1', 'measured T2', 'measured T3', 'measured T4', 'measured T5', 'measured T6',...
			'calculated  T1', 'calculated  T2', 'calculated  T3', 'calculated  T4', 'calculated  T5', 'calculated  T6', 0);
		xlabel('\alpha [°]'); ylabel(Label2); grid; li = get(gca, 'YLim'); li = max(abs(li));
		set(gca, 'XLim', [0 360], 'XTick', [0 45 90 135 180 225 270 315 360], 'YLim', [-li +li]); 
	end
	if (nargin > 4)
		t = title({PlotTitle, PlotSubTitle}); 
	elseif (nargin > 3)
		t = title(PlotTitle); 
	else
		t = title('plot of 2 signal sets'); 
	end
	set(t, 'Interpreter', 'none');
return
	
