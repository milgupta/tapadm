function AVec = rotta(TVec);
% RotTA     rotates the Nx3 matrix TVec with N vectors
%           from the TAPAD into the Autokal coordinate system. 
%           The 'Autokal' coordinate system is a rough approximation for the 
%           skull-fixed coordinate system. The transformation from 
%           TAPAD to Autokal related coordinates is performed by a
%           135° rotation around the Z-Axis, followed by a rotation 
%           of approximately 33° around Y.
%           Notice that the function does not rescale the vectors,
%           so to transform position-vectors use:
%
%           Autokal_Positions = 1000 * RotTA(TAPAD_Positions);     
%
%           see RotAT.m
 
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

persistent RotTA;		% Matrix for coordinate transformation TAPAD -> Autokal 
if (isempty(RotTA)) 
   TA = [ -1/sqrt(3)   1/sqrt(3)  -1/sqrt(3);...
          -1/sqrt(2)  -1/sqrt(2)      0;...
          -1/sqrt(6)   1/sqrt(6)   2/sqrt(6)]; 
end
AVec = (TA * TVec')';								
