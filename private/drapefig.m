function drapefig(scale);
% DRAPEFIG plots stylized shadows of a head on the 'walls' of the 3-D-plot. 
%            Scale defines the axis extent of the plot.   
%            
%            fh = drapefig(scale);

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

p0 = [0 0 0];									% axis scaling
PRange = scale * ones(3,1);
Range = max(PRange) * 0.6;
XYAxisScale = [(p0(1)-Range) (p0(1)+Range) (p0(2)-Range) (p0(2)+Range)];
YZAxisScale = [(p0(2)-Range) (p0(2)+Range) (p0(3)-Range) (p0(3)+Range)];
XZAxisScale = [(p0(1)-Range) (p0(1)+Range) (p0(3)-Range) (p0(3)+Range)];


circ = 0 : pi/16 : 2*pi;			% Head-shape for sagittal section
rscale = 0.8 * ones(size(circ)); rscale(1) = 1; 
circ = [circ circ(1)]; rscale = [rscale rscale(1)];
rscale = rscale * min([(XZAxisScale(2) - XZAxisScale(1)) (XZAxisScale(4) - XZAxisScale(3))]) * 0.5; 
HeadX = rscale .* cos(circ); HeadZ = rscale .* sin(circ);
XZAxisLen = [(XZAxisScale(2) - XZAxisScale(1)) (XZAxisScale(4) - XZAxisScale(3))]  ;
HeadZ = HeadZ + XZAxisScale(4) - abs(XZAxisLen(1)/2);
HeadX = HeadX + XZAxisScale(2) - abs(XZAxisLen(2)/2);

circ = 0 : pi/16 : 2*pi;			% Head-shape for axial section
rscale = 0.8 * ones(size(circ)); rscale(1) = 1; rscale(end) = 1;
rscale(8) = 0.9; rscale(9) = rscale(8); rscale(10) = rscale(8);
rscale(24) = rscale(8); rscale(25) = rscale(8); rscale(26) = rscale(8);
rscale = rscale * min([(XYAxisScale(2) - XYAxisScale(1)) (XYAxisScale(4) - XYAxisScale(3))]) * 0.5; 
AHeadX = rscale .* cos(circ); AHeadY = rscale .* sin(circ);
XYAxisLen = [(XYAxisScale(2) - XYAxisScale(1)) (XYAxisScale(4) - XYAxisScale(3))]  ;
AHeadY = AHeadY + XYAxisScale(4) - abs(XYAxisLen(1)/2);
AHeadX = AHeadX + XYAxisScale(2) - abs(XYAxisLen(2)/2);

circ = 0 : pi/16 : 2*pi;			% Head-shape for coronal section
rscale = 0.8 * ones(size(circ)); %rscale(1) = 1; rscale(end) = 1;
rscale(16) = 0.9; rscale(17) = rscale(16); rscale(18) = rscale(16);
rscale(32) = rscale(16); rscale(33) = rscale(16); rscale(1) = rscale(16); rscale(2) = rscale(16);
rscale =  rscale * min([(YZAxisScale(2) - YZAxisScale(1)) (YZAxisScale(4) - YZAxisScale(3))]) * 0.5; 
CHeadY = rscale .* cos(circ); CHeadZ = rscale .* sin(circ);
YZAxisLen = [(YZAxisScale(2) - YZAxisScale(1)) (YZAxisScale(4) - YZAxisScale(3))]  ;
CHeadY = CHeadY + YZAxisScale(2) - abs(YZAxisLen(1)/2);
CHeadZ = CHeadZ + YZAxisScale(4) - abs(YZAxisLen(2)/2);

hold off;
plot3(AHeadX, AHeadY, -scale*ones(size(AHeadX)), ':k'); grid on;
axis([-scale scale -scale scale -scale scale]); 
xlabel('X_{A} [mm]'); ylabel('Y_{A} [mm]'); zlabel('Z_{A} [mm]');
hold on;
plot3(scale*ones(size(CHeadY)), CHeadY, CHeadZ, ':k'); 
plot3(HeadX, scale*ones(size(HeadX)), HeadZ, ':k'); 
hold off;
