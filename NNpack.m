function [ NNp, pos ] = NNpack( coordinates, NN, minN)
%Name:NNPackage, short for Nearest Neighbour package
%FUNCTION: 
%Takes the position data and the NN ids and packages it into a form where
%the coordinates are given in a structure form. It gives the data in the
%same form as the outpur of qcorrelation. 
%Written by Nick Orr 17 March 2020;
%removes any points with fewer than minN number of neighbours

%INPUTS:
%coordinates, the position data from which NN is created. 
%NN, the nearest neighbour structure from the coordinates. 
%minN, the minimum number of nearest neighbours a cluster must have to be
%packed into the NNp structure. 
%
%OUTPUTS: 
%NNp a long structure with the positions of the centreal particles and its
%neighbours. 
%pos, an array: the positions of the particles with more than or equal to
%the minN number of neighbours. The fourth column is the id from the
%orignal array.



L = numel(NN.LISTN);
NNlongP = struct; 
%NNp = repmat(NNlongP, 1, L); 
LISTN = NN.LISTN;
%D = NN.D; 
LIST = NN.LIST; 
LISTNabove = LISTN(LISTN >= minN);
LISTNaboveID = find(LISTN >= minN);
NNp = repmat(NNlongP, 1, numel(LISTNaboveID));
pos = zeros(numel(LISTNabove), 4); 
for a = 1:numel(LISTNabove)
    NNp(a).pos = coordinates(LISTNaboveID(a),:); 
    pos(a, 1:3) = coordinates(LISTNaboveID(a),:);
    pos(a, 4) = LISTNaboveID(a);
   %NNp(a).Nc = LISTN(a); 
    NNa_ids = LIST(LISTNaboveID(a),:);
    NNa_idsNE0 = NNa_ids(NNa_ids ~= 0);
    NNp(a).Npos = coordinates(NNa_idsNE0,:);
    NNp(a).ID = LISTNaboveID(a); 
end
end 
