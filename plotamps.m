function [fh, sfh] = plotamps(FileName, channel, fh);
% PLOTAMPS plots amplitude data. 
%           The function opens a new figure (if fh is omitted), or it plots 
%           in the figure given by figure handle fh. It generates six sub-plots
%           and plots the signal amplitude. If possible, it also evaluates the 
%           DerivativParameters of the data and shows the estimated amplitude.
%
%            plotamps(FileName, channel, fh);
%
%           The function returns the current figure and a callback handle,
%           which can be used to move a cursor through the sub-plots.
%
%---------------------------------------------------------------------
% Copyright © 2007 by Andreas Zierdt (Anderas.Zierdt@phonetik.uni-muenchen.de)
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

errormsg = nargchk(2, 3, nargin);
if (nargin < 3)
	fh = [];
end

if (length(channel) ~=1)
	error('channel must be a scalar');
end
[MAmps, amp_comment, amp_descriptor, amp_dimension, amp_private,...
		amp_samplerate, amp_unit] = loadsdata(FileName); NumPoints = size(MAmps, 1); % load the data
data_idx = 1 : NumPoints;	t = (1/amp_samplerate) * data_idx';

% Does the data contains DerivativParameters?
fn = fieldnames(amp_private);
has_derivative = ~isempty(strmatch('DerivativParameters', fn, 'exact')); 
if (has_derivative) DParams = amp_private.DerivativParameters; end

IDX_CEIL = decround(NumPoints, NumPoints/10); % calculate the decimal cap (upper limit)

% if no figure handle is given, create a new figure 
if (isempty(fh))
	fh = figure;	
else
	figure(fh);
end	

a = MAmps(:, 1:6, channel); AMP_CEIL = decround(a(:)); clear a;
for (coil=1:6)
	subplot(3, 2, coil)	
	if (has_derivative)
		[approx_amp, approx_vel] = trig_approx_func(DParams(:, coil, channel)', t);
		h = plot(data_idx, MAmps(:, coil, channel), 'b.', data_idx, approx_amp, 'k-');
		set(h(2), 'Color', 0.5*[1 1 1]);
	else
		h = plot(data_idx, MAmps(:,coil, channel), 'b.');
	end
	set(h(1), 'MarkerSize', 1);	
	axis([0 IDX_CEIL -AMP_CEIL AMP_CEIL]); ylabel(['Amp T' num2str(coil) ' [1]']);
	line([0 0], get(gca, 'YLim'), 'EraseMode', 'xor', 'Tag', ['AmpXCursor_' num2str(coil)]);
end
xlabel('idx [1]');
subplot(3, 2, 1);
t = title([FileName ' channel ' num2str(channel)]); set(t, 'Interpreter', 'none');
if (has_derivative)
	legend('measured Amp', 'approx. Amp', 0);
else
	legend('measured Amp', 0);
end
fh = gcf; sfh = @set_cursors;

return
%------------- subfunctions -----------

% set cursors in all subplots (can be called from outside)  	
function set_cursors(idx, cf)
	for (coil=1:6)
		c = findobj('Tag', ['AmpXCursor_' num2str(coil)]);		set(c, 'XData', [idx idx]);
	end	
return	
