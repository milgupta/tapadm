function [LocalPos, LocalOrientVec] = trans2local(GlobalPos, GlobalOrientVec);
% TRANS2LOCAL takes a Nx3 matrix with positions and another Nx3 matrix with 
%           orientation vectors and transforms them to local coordinates according to
%           spatial placement and alignment of the six transmitter coils. Thus, the 
%           result will consist of two a Nx3x6 arrays, containing the given coordinates
%           relative to the local coordinate system of each transmitter.
%     
%           [LocalPos, LocalOrientVec] = trans2local(GlobalPos, GlobalOrientVec);
%
%           Notice that the function applies a translation followed by a rotation 
%           on the first argument GlobalPos, but only a rotation on GlobalOrientVec,
%           since orientations are translation-invariant.
%           Also notice that the local z-axis is the symmetry axis of the coil with 
%           its origin at the coil-center. The local coordinate systems is cartesian and 
%           'right handed'. So if you spread the first 3 fingers of your right hand 
%           orthogonal, thumb is x, forefinger y and middle finger z. 
%           (Due to the cylindrical symmetry of the magnetic field, the exact
%           direction of x and y does not matter for field computations.)

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

% Define spatial placement and alignment of the six transmitter coils. 
% The resulting transformation-matrices shall persist between function calls
% to increase speed after the first execution
persistent RK CoilPos InvCoilTrans;		
if (isempty(RK) | isempty(CoilPos) | isempty(InvCoilTrans)) 
	RK		= 0.3375;							% radius of the coil-mounting structure
%	SYSTEM = 1;									% uncomment this line for old AG500 systems
	SYSTEM = -1;								% uncomment this line for new AG500 systems
	CoilPos = [0 0 -1; 0 -1 0; -1 0 0;...	% Transmitter - coil positions (CoilPos*RK)
				  0 0 1;  0 1 0;   1 0 0]; 
	ODS2	= 1/sqrt(2);							% "OneDividedbySqrt2" equals cos(45)
	CoilAlign= [1 0 0; ODS2    ODS2 0;...	% Transmitter - coil alignment
					0 1 0; 0       ODS2 ODS2;...
					0 0 1; SYSTEM*ODS2    0   ODS2];
	CoilTrans= zeros(3,3,6);				% 6 rotary matrices for local coordinate transformation
	InvCoilTrans= zeros(3,3,6);			% each column contains the x-, y- and z-coordinate of the
													% local coordinate basis vector expressed in terms of the
	for coil=1:6								% global (TAPAD) coordinate system!
		CoilTrans(:, 3, coil) = CoilAlign(coil, :)'; % local z = orientation of coil axis	
		% the local y axis is perpendicular to a line connecting the coil with TAPAD-Orign and local z
		CoilTrans(:, 2, coil) = cross(CoilPos(coil, :), CoilAlign(coil, :))';
		CoilTrans(:, 2, coil) = CoilTrans(:, 2, coil) ./ norm(CoilTrans(:, 2, coil));
		% local x is perpendicular to z and y (and equal -CoilAlign for T1,3,5) 
		CoilTrans(:, 1, coil) = cross(CoilTrans(:, 2, coil), CoilTrans(:, 3, coil)); 
		if norm(CoilTrans(:, 1, coil)) > 0;
			CoilTrans(:, 1, coil) = CoilTrans(:, 1, coil) ./ norm(CoilTrans(:, 1, coil));
		end
		% Calc inverse matrix, which describes the global (TAPAD) coordinate basis vectors in terms of
		% each local coordinatesystem
		InvCoilTrans(:,:,coil) = inv(CoilTrans(:,:,coil));
	end
end 	% of initialization

[NumPoints, dummy] = size(GlobalPos);
%LocalPos = zeros(NumPoints, 3, 6);		LocalOrientVec = zeros(NumPoints, 3, 6);

% transform positions
LocalPos = translate(GlobalPos, CoilPos, RK, NumPoints); LocalPos = rotate(LocalPos, InvCoilTrans, NumPoints);

% transform orientations
GlobalOrientVec = cat(3, GlobalOrientVec, GlobalOrientVec, GlobalOrientVec, GlobalOrientVec, GlobalOrientVec, GlobalOrientVec);
LocalOrientVec = rotate(GlobalOrientVec, InvCoilTrans, NumPoints);

%------------- subfunctions -----------
% performs translatory-part of the coordinate-transformation without time-consuming for-loops	
function  T = translate(GlobalPos, CoilPos, RK, NumPoints);
	CPbig = RK*cat(3,	repmat(CoilPos(1, :), NumPoints, 1), repmat(CoilPos(2, :), NumPoints, 1),...
							repmat(CoilPos(3, :), NumPoints, 1), repmat(CoilPos(4, :), NumPoints, 1),...
							repmat(CoilPos(5, :), NumPoints, 1), repmat(CoilPos(6, :), NumPoints, 1));
	T = cat(3, GlobalPos, GlobalPos, GlobalPos, GlobalPos, GlobalPos, GlobalPos) - CPbig;

% performs rotation-part of the coordinate-transformation without time-consuming for-loops	
function  R = rotate(LocalPos, InvCoilTrans, NumPoints);
	LPbig = cat(3,	LocalPos(:, :, 1)', LocalPos(:, :, 2)', LocalPos(:, :, 3)',...
						LocalPos(:, :, 4)', LocalPos(:, :, 5)', LocalPos(:, :, 6)');
	LPbig = repmat(LPbig(:)', 3, 1);
	ICTbig = repmat(InvCoilTrans, 1, NumPoints);
	ICTbig = cat(2,	ICTbig(:, :, 1), ICTbig(:, :, 2), ICTbig(:, :, 3),...
							ICTbig(:, :, 4), ICTbig(:, :, 5), ICTbig(:, :, 6)); 
	R = ICTbig .* LPbig;
	R = reshape(R', 3, NumPoints*3*6);
	R = sum(R);
	R = reshape(R', NumPoints*6, 3);
	np = 1:NumPoints;
	R = cat(3,	R(0*NumPoints+np, :), R(1*NumPoints+np, :), R(2*NumPoints+np, :),...
					R(3*NumPoints+np, :), R(4*NumPoints+np, :), R(5*NumPoints+np, :));

% performs translatory-part of the coordinate-transformation slow but understandable	
function  T = translate_slow(GlobalPos, CoilPos, RK, NumPoints);
	for Coil=1:6		
		LocalPos(:, :, Coil) = GlobalPos -  repmat(RK*CoilPos(Coil, :), NumPoints, 1); 
	end
	T = LocalPos;

% performs rotational-part of the coordinate-transformation slow but understandable	
function  R = rotate_slow(LocalPos, InvCoilTrans, NumPoints);
	for Coil=1:6						
		for i=1:NumPoints
			LocalPos(i, :, Coil) = (InvCoilTrans(:,:,Coil) * LocalPos(i, :, Coil)')'; 
		end
	end
	R = LocalPos;
	