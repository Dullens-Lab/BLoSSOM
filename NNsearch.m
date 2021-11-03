function [ NN ] = NNsearch( Data, maxN, Rcut )
%NAME: NNsearch
%FUNCTION: Search for the nearest neighbours for each particle. The metric
%to find nearest neighbours is: 
% "The up to maxN nearest points within a distance Rcut of each particle"
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
%Written by Nick Orr, 2019. 
%
n = numel(Data(:,1)); 
NNlist = zeros(n,maxN);
NDlist = zeros(n,maxN);
NNlistN = zeros(n,1);
NN = struct;
MA = createns(Data, 'NSMethod','kdtree','Distance','euclidean');
[Idx, w] = knnsearch(MA,Data,'K',maxN + 1); 
%crash
for a = 1:n
    adist = w(a,:);
    Dists = adist(adist <= Rcut & ne(adist, 0) );
    NNlistN(a) = numel(Dists); 
    NNlist(a, 1:NNlistN(a)) = Idx(a, adist <= Rcut & ne(adist, 0)); 
    NDlist(a, 1:NNlistN(a)) = Dists; 
end 
NN.LISTN = NNlistN; 
NN.D = NDlist; 
NN.LIST = NNlist;  

end

