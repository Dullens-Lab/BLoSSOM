function [Clustcell, uncorrelated, degM12ARR]...
    = PositionOrientationCorrelation(pos, Orientation,...
    Operators1, Operators2, Rcut, anglethresh)
%NAME: OrientationCorrelation
%Function: Assigns points with orinetations to a grain. It does this by a
%triangulation followed by a misorientation calculation over each
%triangulation edge. Edges with low misorientations form a graph, from
%which the connected regions, the grains, are resolved. 
%INPUTS:
%   pos - x y z Cartesian coordinates as an n by 3 array. 
%   Orientation - quaternion representations of the orientations in an n by
%   4 array.
%   Operators1 - rotational symmetry operators required for calculating the
%   misorientation angles, expressed as quaternions in an n by 4 array.
%   Please see the variable saves. Operators for the symmetry of point1; 
%   Operators2 - likewise for point2.
%   Rcut - the max distance between two connected coordinates.
%   anglethresh - the maximum misorientation angle tolerated for two points
%   to be considered connected. Degrees.
%OUTPUTS: 
%   ClustCell - The same information but given in a cell. Each cell
%   pertains to a new grain. The members of each cell are the IDs of the
%   member points of that grain.
%   Uncorrelated  - the points with no orientation correlation. Don't
%   forget to combine these points with any that failed regstration if you
%   want to look for all of the "unordered points".
%   degM12ARR - misorientation angles between connected points in degrees.
%HISTORY: 
%   Written by Nick Orr 1st April 2020.  
%   BACKGROUND: 
%   The idea for this function was to replace the crystallinity search
%   using bond coherence Steinhardt order parameters. Grains may be
%   identified only by looking for structure, orientation and positional
%   correlations. The missing piece was a robust way to measure orientation
%   correlations. 
%   Nick Orr 7 April 2020. Removed the clustV output. Changed the function
%   name to PositionOrientationCorrelation.
%CODE START: 

%perform triangulation 
DT = delaunayTriangulation(pos(:,1), pos(:,2), pos(:,3)); 
DTedges = edges(DT); %traingulation edges 
pointsL = pos(DTedges(:,1),1:3);
pointsR = pos(DTedges(:,2),1:3);
pointsL2 = pointsL.^2;
pointsR2 = pointsR.^2;
DTdistances = sqrt(sum((pointsL - pointsR).^2,2)); %Distance calc 
figure; histogram(DTdistances, 0:(20*Rcut)/100:(50*Rcut)); xlabel 'Delaunay distance'; ylabel freq %%!!!!!!!!!
connected_pairs = DTedges(DTdistances < Rcut, :); %connect close points
GT = connected_pairs; 
%numel(GT(:,1))/ numel(pos(:,1))% the mean number of connections. ~10
%then plow through these pairs and save the ones with low misorientation.
tic
M12ARR = zeros(numel(GT(:,1)), 1); 
paircoords = zeros(numel(GT(:,1)), 3); 
for a = 1:numel(GT(:,1))
    q1 = Orientation(GT(a,1), :);
    q2 = Orientation(GT(a,2), :);
    q2 = [q2(1), -q2(2:4)]; %take the complex conjugate
    M12ARR(a,1) = MisorientationQA(q1,q2, Operators1, Operators2); %calc mis angle
    paircoords(a,:) = (pos(GT(a,1),1:3) + pos(GT(a,2), 1:3))./2;
end 
toc % shows you how long it took for the misorientatin function to complete 
%crash% DEBUG
degM12ARR = (M12ARR./(2*pi)) * 360;
figure; histogram(degM12ARR, 0:20/100:20); xlabel misorientation; ylabel freq %%!!!!!!!!!!!!!
GTM = GT(degM12ARR < anglethresh,:); % set the angle thresh
Gr = graph(GTM(:,1), GTM(:,2));
bin = conncomp(Gr);
binU = unique(bin); 
%bin only contains correlated members if it has more than one member.
N = zeros(numel(binU),1); 
for a = 1:numel(binU)
   N(a) = numel(find(bin == binU(a)));  
end
binUL = binU(N >= 2); %must have more than 1 member!
ID = 1:numel(pos(:,1));
%ClustV = binUL;
Clustcell = cell(numel(binUL), 1);
Vc = zeros(numel(ID), 1); 
d = 1; 
for a = 1:numel(binUL)
    Clustcell{a} = (ID(bin == binUL(a)))'; %invert to keep Grain_Plot happy 
    V = ID(bin == binUL(a)); 
    Vc(d: numel(V) + d - 1) = V;
    d = d + numel(V); 
end
Vcp = Vc(Vc > 0); 
IDuc = setdiff(ID, Vcp); 
uncorrelated = pos(IDuc, 1:3); 
%ClustV = bin(Vcp); 
%package for output.
end

