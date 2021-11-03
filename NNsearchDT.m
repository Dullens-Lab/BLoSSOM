function [NN] = NNsearchDT(Data, maxN, Rcut)
%NAME: NNsearchDT
%FUNCTION: Search for the nearest neighbours for each particle. The metric
%to find nearest neighbours is: 
% "The up to maxN nearest points connected by a Delaunay triangulation edge 
%within a distance Rcut of each particle"
%INPUTS: 
%   Data - n by 3 array of coordinates, new row new particle. 
%   maxN - The maximum number of neighbours a particle can have. 
%   Rcut - the cut of distance for the neighbour search. 
%OUTPUTS: 
%   NN - strutcure with three fields. 
%       1) NNListN - the number of neughbours for the particle 
%       2) NNList - the list of neighbour particle IDs
%       3) NDlist - the distances between neihgbours
%HISTORY: 
%Written by Nick Orr, July 2021. 
%
DT = delaunayTriangulation(Data(:,1), Data(:,2), Data(:,3)); 
DTedges = edges(DT); %traingulation edges 
pointsL = Data(DTedges(:,1),1:3);
pointsR = Data(DTedges(:,2),1:3);
DTdistances = sqrt(sum((pointsL - pointsR).^2,2)); %Distance calc 
connected_pairs = DTedges(DTdistances < Rcut, :); %connect close points
GT = connected_pairs;
GTD = DTdistances(DTdistances < Rcut, :);
GTD_tot = [GTD; GTD];
NNlist = zeros(numel(Data(:,1)),maxN);
NDlist = zeros(numel(Data(:,1)),maxN);
NNlistN = zeros(numel(Data(:,1)),1);
NN = struct;
GT_tot = [GT;[GT(:,2), GT(:,1)]];
GTU = unique(GT_tot(:,1)); 
GTU = sort(GTU, 'ascend');
%crash
%Idx = 1:numel(Data(:,1));
for a = 1:numel(GTU(:,1))
    NNa = GT_tot(GT_tot(:,1) == GTU(a), 2);
    [GTDs,id2] = sort(GTD_tot(GT_tot(:,1) == GTU(a)), 'ascend');
    if numel(NNa) > maxN
        %GTDs = GTDs(1:maxN);
        id2 = id2(1:maxN);
    end
    NNlistN(GTU(a)) = numel(id2); 
    NNlist(GTU(a), 1:NNlistN(GTU(a))) = NNa(id2);
    NDlist(GTU(a), 1:NNlistN(GTU(a))) = GTDs(id2); 
end 
NN.LISTN = NNlistN; 
NN.D = NDlist; 
NN.LIST = NNlist;  
end

%I dont think the ids are right here. 