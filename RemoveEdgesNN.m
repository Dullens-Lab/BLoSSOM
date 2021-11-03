function [NN_B, bulkid] = RemoveEdgesNN(coords, NN, depth)
%NAME: Remove Edges Nearest Neighbours
%FUNCTION: remove particles that are on the edges of the FOV.
%INPUTS: 
%   coords - xyz coordinates in an n by 3 array. 
%   NN - Nearest neighbour structure, from the output of NNsearch 
%   depth - the depth of the edge cut off. The bigger the depth, the more
%   edge will be removed. 
%OUTPUTS: 
%   NN_B - Nearest neighbour from the bulk as a strcuture with the same fom
%   as the input structure NN.
%   bulkid - the logical index for the bulk coordinates. 
%HISTORY: 
%Written by Nick Orr 14/04/20; 
%
maxx = max(coords(:,1)) - depth; 
maxy = max(coords(:,2)) - depth; 
maxz = max(coords(:,3)) - depth; 
minx = min(coords(:,1)) + depth; 
miny = min(coords(:,2)) + depth; 
minz = min(coords(:,3)) + depth; 
bulkid = coords(:,1) < maxx & coords(:,1) > minx  & ...
    coords(:,2) < maxy & coords(:,2) > miny &...
    coords(:,3) < maxz & coords(:,3) > minz; 
LISTN = NN.LISTN; 
D = NN.D; 
LIST = NN.LIST; 
LISTNB = LISTN(bulkid, :); 
DB = D(bulkid, :); 
LISTB = LIST(bulkid, :); 
NN_B.LISTN = LISTNB; 
NN_B.D = DB; 
NN_B.LIST = LISTB; 

end

