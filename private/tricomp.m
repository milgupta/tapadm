function coverage = tricomp(Tri1, Tri2)
% TRICOMP   Compares two triangles within the unit circle, each given by
%           a vector with 3 angles (in rad). The function returns a
%           value between 0 and 1, which describes, how many of the area  
%           of the 1st triangle is covered by the 2nd. 
%
%           coverage = tricomp(Tri1, Tri2);

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
% Cartesian coordinates of the triangle points, columns containig x and y
T1_Points = [Tri1 Tri1(1)]; T1_Points = [cos(T1_Points); sin(T1_Points)]; 
T2_Points = [Tri2 Tri2(1)]; T2_Points = [cos(T2_Points); sin(T2_Points)]; 

% Given the 3 points of each triangle, we build linear equations
% to describe the edges of the two triangles. (Point and direction vector form)
% straight line i:  T1_Points(:, i) + c * T1_Edges(:, i)
T1_Edges = CalcEdges(T1_Points); T2_Edges = CalcEdges(T2_Points); 

% Calculate the area of the first triangle, which will be the reference for
% the degree of overlap
T1_Area = CalcArea(T1_Edges);

% Calculate all Intersections of the two triangles, i.e. find points where:
% T1_Points(:,i) + c1 * T1_Edges(:,i)     = T2_Points(:,i) + c2 * T2_Edges(:,i)
% c1 * T1_Edges(:,i) - c2 * T2_Edges(:,i) = T2_Points(:,i) - T1_Points(:,i)
% ==> find c with:             E * c      =    p
Intersections = [];
for (i = 1:3)	% i iterates edges of T1
	for (j = 1:3)	% j iterates edges of T2
		E = [T1_Edges(:, i) T2_Edges(:, j)]; p = T2_Points(:, j) - T1_Points(:, i); c = E\p; 
		sol1 = T1_Points(:, i) + c(1) * T1_Edges(:, i);	sol2 = T2_Points(:, j) - c(2) * T2_Edges(:, j);
		if (norm(sol1-sol2) > sqrt(eps))
			[i j]
			[sol1 sol2]
			error('solutions differ');
		end
		if (norm(sol1) <= 1)
			if (isempty(Intersections) | isempty(find((Intersections(1, :) == sol1(1)))))
				Intersections = [Intersections sol1];
			end
		end
	end
end
NI = size(Intersections, 2);

if (NI < 3)
	warning('No overlap, return 0');
	coverage = 0;
elseif (NI == 3)
	I_Edges = CalcEdges([Intersections Intersections(:,1)]); 
	coverage = CalcArea(I_Edges) / T1_Area;
else	
	coverage = PolyArea(Intersections) / T1_Area;
end

%------------- subfunctions -----------
function a = PolyArea(Polygon)
	N = size(Polygon, 2);
	if (N < 3)
		a = 0;
	elseif (N == 3)
		edges = CalcEdges([Polygon Polygon(:,1)]); 
		a = CalcArea(edges);
	else	
		[Poly_A, Poly_B] = DecomposePolygon(Polygon);
		a = PolyArea(Poly_A);
		a = a + PolyArea(Poly_B);
	end
return

function [Poly_A, Poly_B] = DecomposePolygon(Point)
	NI = size(Point, 2);
	Distance = zeros(NI, NI); % compute euclidian distances between all points
	Edge = []; EdgeIdx = [];
	for i=1:NI-1
		for j=i+1:NI
			Distance(i, j) = norm( Point(:, i) - Point(:, j) );
			EdgeIdx = [EdgeIdx; i j]; 
			Edge = [Edge Point(:, i) - Point(:, j)];
		end
	end
	NE = size(Edge, 2);
	
	Affiliation = zeros(NE, NI);
	for (EdgeNr = 1:NE)
		Diagonal = [Point(:, EdgeIdx(EdgeNr, 1)) Point(:, EdgeIdx(EdgeNr, 2))];
		% construct a direction vector on this diagonal line
		g = Diagonal(:, 1) - Diagonal(:, 2); 	
		n = [-g(2) g(1)]; n = n / norm(n);	% a normal-vector, perpendicular to the diagonal line
		% Each intersection which lays not on the main diagonal is assigned to one of the two (sub-)polygons,
		% that lay at right and left hand side of the main diagonal line 	
		for i=1:NI									% for all intersections...
			if ((EdgeIdx(EdgeNr, 1) ~= i) & (EdgeIdx(EdgeNr, 2) ~= i))
				l = (Point(:, i) - Diagonal(:, 1))' * g / (g' * g); 
				d = Point(:, i) - ( Diagonal(:, 1) + l * g);
				Affiliation(EdgeNr, i) = n * d;	% sign marks points 'right' or 'left' of the line, zero if point lays on the edge
			end
		end
	end
	PossibleEdgeIdx = find(sum(Affiliation>0,2).*sum(Affiliation<0,2)); % Edges, with at least one point on both sides

	Balance = zeros(1, length(PossibleEdgeIdx)); 
	for (i=1:length(PossibleEdgeIdx))
		EdgeNr = PossibleEdgeIdx(i);
		a = Affiliation(EdgeNr, :); ap = find(a>0); am = find(a<0); 
		Balance(i) = length(a(ap)) * length(a(am));
	end
	[m, i] = max(Balance); EdgeNr = PossibleEdgeIdx(i);
	a = Affiliation(EdgeNr, :); [ap, Poly_A_Idx] = find(a>0); [am, Poly_B_Idx] = find(a<0); 
	Diagonal = [Point(:, EdgeIdx(EdgeNr, 1)) Point(:, EdgeIdx(EdgeNr, 2))];
	
	Poly_A = zeros(2, length(Poly_A_Idx)+2);
	for i=1:length(Poly_A_Idx)
		Poly_A(:, i+1) = Point(:, Poly_A_Idx(i));
	end	
	if (norm(Diagonal(:, 1) - Poly_A(:, 2)) <= norm(Diagonal(:, 1) - Poly_A(:, end-1)) )
		Poly_A(:, 1) = Diagonal(:, 1);	Poly_A(:, end) = Diagonal(:, 2);
	else
		Poly_A(:, 1) = Diagonal(:, 2); 	Poly_A(:, end) = Diagonal(:, 1);
	end
	
	Poly_B = zeros(2, length(Poly_B_Idx)+2);
	for i=1:length(Poly_B_Idx)
		Poly_B(:, i+1) = Point(:, Poly_B_Idx(i));
	end	
	if (norm(Diagonal(:, 1) - Poly_B(:, 2)) <= norm(Diagonal(:, 1) - Poly_B(:, end-1)) )
		Poly_B(:, 1) = Diagonal(:, 1);	Poly_B(:, end) = Diagonal(:, 2);
	else
		Poly_B(:, 1) = Diagonal(:, 2); 	Poly_B(:, end) = Diagonal(:, 1);
	end	
return

function Edges = CalcEdges(Points)
	Edges = zeros(2, 3);
	for (i = 1:3)
		Edges(:, i) = Points(:, i+1) - Points(:, i);
	end
return

function area = CalcArea(Edges)		% Calculate the area of a triangle
	s = ( norm(Edges(:, 1)) + norm(Edges(:, 2)) + norm(Edges(:, 3)) ) / 2;
	area = sqrt( s * (s-norm(Edges(:, 1))) * (s-norm(Edges(:, 2))) * (s-norm(Edges(:, 3))));
return