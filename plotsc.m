function [fh, sfh] = plotsc(data, channel, point_idx, fh, style, ccfh);
% PLOTSC plots single channel position data.  
%           This function plots the result data, i.e. position and orientation of 
%           a TAPADM run for a single channel. It can be used in two different ways:
%           In the default (3d) mode, it plots the 3D-Trace of a channel.
%           In 'flat' mode, it plots the 7 columns of a result data set in different 
%           subfigures, and allows to browse the data with synchronized cursors.
%            
%            fh = plotsc(data, channel, point_idx, fh, style);
%            fh = plotsc(data, channel, point_idx, fh, 'flat', ccfh);
%
%           To maintain compatibility with TAPAD and CalcPos, the function
%           plots column no 7 as 'iterations', but TAPADM stores there a
%           exit flag and saves the iteration depth in a different file.
%           For a proper display, copy the iterations intow row 7 before calling 
%           this function in flat mode. 
%
%           In flat mode, a cursor-callback function handle (ccfh) can be
%           used, which function is called every time the slider is moved.
%           The function is called with two arguments: 
%           1. the actual point index
%           2. the figure handle
%           This callback can be used to pass the selected cursor position
%           to other functions.
%
%           Again, the function returns as 2nd value the handle of the
%           internal set_cursors function, which can be used to remotely
%           set/move the cursors. (But without moving the slider)


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

if  (nargin < 6)
	ccfh = [];
end
if (nargin < 5)
	style = 'g-';
end
if (nargin < 4)
	fh = [];	
end
if (strcmpi('flat', style))
	flat_plot(data, channel, point_idx, fh, style, ccfh);
	fh = gcf;
	sfh = @set_cursors;
else
	spatial_plot(data, channel, point_idx, fh, style);
	sfh = [];
end	
fh = gcf;
return

%------------- subfunctions -----------
function spatial_plot(data, channel, point_idx, fh, style)
	l = findstr('+q', style); addquiver = ~isempty(l);
	if (~isempty(l))
		style(l:l+1) = [];
	end

	if (isempty(fh))
		fh = figure;	
	else
		figure(fh);
		hold on;
	end	
	if ((nargin < 3) | isempty(point_idx))
		[NumPoints, dummy1, dummy2] = size(data);
		point_idx = 1 : NumPoints;
	end	

	x = data(point_idx, 1, channel); y = data(point_idx, 2, channel); z = data(point_idx, 3, channel);
	d = sqrt(x .* x + y .* y + z .*z); maxd = ceil(max(d));
	[ox, oy, oz] = sph2cart(data(point_idx, 4, channel)' * 180/pi, data(point_idx, 5, channel)' * 180/pi, 1);

	plot3(x, y, z, style);
	if (nargin < 4)
		axis([-maxd maxd -maxd maxd -maxd maxd]);
		xlabel('Xa [mm]'); ylabel('Ya [mm]'); zlabel('Za [mm]');  grid on;
	end
	if (addquiver) 
		hold on;
		quiver3(x, y, z, ox', oy', oz', 0.05*maxd, 'r.'); quiver3(x, y, z, -ox', -oy', -oz', 0.05*maxd, 'b.');
		hold off;
	end
return

% calculate the decimal cap (upper limit) of the data, e.g. coorcap(3412) = 4000
function cap = coorcap(x) 
	d = max(abs(x)); e = 10^floor(log10(d));	cap = ceil(d/e)*e;

% variant for multi-figure plot	
function flat_plot(data, channel, point_idx, fh, style, ccfh)
	if (isempty(fh))						
		fh = gcf;								% open a new figure
	else
		figure(fh);								% or re-use an existing one
	end	

	% calculate
	NumPoints = size(data, 1);						data_idx = 1 : NumPoints;
	IDX_CEIL = decround(data_idx);				
	
	subplot(2, 3, 1)
	h = plot(data(:, 1, channel), data(:, 2, channel), 'k.'); title('Spatial position');
	axis square; grid on; xlabel('X_A [mm]'); ylabel('Y_A [mm]');
	set(h(1), 'MarkerSize', 1);
	line([data(1, 1, channel) data(1, 1, channel)], get(gca, 'YLim'), 'EraseMode', 'xor', 'Tag', 'XCursor_1'); 
	line(get(gca, 'XLim'), [data(1, 2, channel) data(1, 2, channel)], 'EraseMode', 'xor', 'Tag', 'YCursor_1'); 
	set(gca, 'UserData', data(:,:, channel), 'Tag', 'AxisWithData');

	subplot(2, 3, 4)
	h = plot(data(:, 1, channel), data(:, 3, channel), 'k.'); 
	axis square; grid on; xlabel('X_A [mm]'); ylabel('Z_A [mm]');
	set(h(1), 'MarkerSize', 1);
	line([data(1, 1, channel) data(1, 1, channel)], get(gca, 'YLim'), 'EraseMode', 'xor', 'Tag', 'XCursor_4'); 
	line(get(gca, 'XLim'), [data(1, 3, channel) data(1, 3, channel)], 'EraseMode', 'xor', 'Tag', 'YCursor_4'); 
	
	subplot(2, 3, 2)
	h = polar(data(:, 5, channel)*pi/180, data_idx', 'k.'); set(h, 'MarkerSize', 1); title('\Phi');
	line([0 0], get(gca, 'XLim'), 'EraseMode', 'xor', 'Tag', 'XCursor_2'); 
	line(get(gca, 'XLim'), [0 0], 'EraseMode', 'xor', 'Tag', 'YCursor_2'); 
	
	subplot(2, 3, 5)
	h = polar(data(:, 6, channel)*pi/180, data_idx', 'k.');  set(h, 'MarkerSize', 1); title('\Theta');
	line([0 0], get(gca, 'XLim'), 'EraseMode', 'xor', 'Tag', 'XCursor_5'); 
	line(get(gca, 'XLim'), [0 0], 'EraseMode', 'xor', 'Tag', 'YCursor_5'); 
		
	subplot(2, 3, 3)
	h = semilogy(data_idx, data(:, 6, channel), 'r.'); 
	set(h, 'MarkerSize', 1);
	grid on; title('''RMS-Value''')
	axis([0 IDX_CEIL 1E-3 1E3]); 
	set(gca, 'YMinorGrid', 'off')
%	keyboard
	c = line([0 0], get(gca, 'YLim'), 'EraseMode', 'xor', 'Tag', 'XCursor_3'); 

	subplot(2, 3, 6)
	h = plot(data_idx, data(:, 7, channel), 'k');  set(h, 'MarkerSize', 1); grid on;
	axis([0 IDX_CEIL 0 20]); title('Iterations'); xlabel('idx [1]'); 
	c = line([0 0], get(gca, 'YLim'), 'EraseMode', 'xor', 'Tag', 'XCursor_6'); 

	p = get(gcf, 'Position'); w=p(3); h=p(4);
	h = uicontrol('Style', 'slider', 'Position', [w-(w/30) h/50 w/60 h-h/20],...
		'Tag', 'IdxSlider', 'Min', 1, 'Max', NumPoints, 'Value', 1);
	set(h, 'sliderstep', [1/NumPoints 100/NumPoints]); set(h, 'Callback', {@slider_callback, data(:,:, channel), ccfh});
	
	set(gcf, 'ResizeFcn', @resize_callback);
return

% Manage to adjust the slider size and position to the actual window size
function resize_callback(h, d)
	u = findobj('Tag','IdxSlider');
	h = gcbo;	old_units = get(h,'Units');
	set(h,'Units','pixels');	p = get(h, 'Position');
	set(u, 'Position', [p(3)-(p(3)/30) p(4)/50 p(3)/60 p(4)-p(4)/20]);
	set(h, 'Units', old_units);
return

% set cursors in all subplots (can be called from outside)  	
function set_cursors(idx)
	da = findobj('Tag', 'AxisWithData'); data = get(da, 'UserData');
	cx = findobj('Tag', 'XCursor_1'); set(cx, 'XData', [data(idx, 1) data(idx, 1)]); 
	cy = findobj('Tag', 'YCursor_1'); set(cy, 'YData', [data(idx, 2) data(idx, 2)]);
	cx = findobj('Tag', 'XCursor_4'); set(cx, 'XData', [data(idx, 1) data(idx, 1)]); 
	cy = findobj('Tag', 'YCursor_4'); set(cy, 'YData', [data(idx, 3) data(idx, 3)]);
	c = findobj('Tag', 'XCursor_3'); set(c, 'XData', [idx idx]);
	c = findobj('Tag', 'XCursor_6'); set(c, 'XData', [idx idx]);
	x = cos(data(idx, 5)*pi/180)*idx;	y = sin(data(idx, 5)*pi/180)*idx;
	cx = findobj('Tag', 'XCursor_2'); set(cx, 'XData', [x x]); 
	cy = findobj('Tag', 'YCursor_2'); set(cy, 'YData', [y y]);
	x = cos(data(idx, 6)*pi/180)*idx;	y = sin(data(idx, 6)*pi/180)*idx;
	cx = findobj('Tag', 'XCursor_5'); set(cx, 'XData', [x x]); 
	cy = findobj('Tag', 'YCursor_5'); set(cy, 'YData', [y y]);
return	
	
% move cursors in all subplots when the slider has moved 	
function slider_callback(h, d, data, ccfh)
	idx = round(get(h, 'Value'));			% slider position gives the index of the actual selected point
	set_cursors(idx);
	
	if (~isempty(ccfh))						% if defined, call the cursor callback function
		feval(ccfh, idx, gcf);
	end
return	
	