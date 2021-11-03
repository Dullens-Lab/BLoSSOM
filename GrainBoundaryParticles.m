function [GBpts] = GrainBoundaryParticles(Apos, GBpos, Rcut)
%NAME: GrainBoundaryParticles
%FUNCTION: Identified grain boundary particles as those within Rcut of its
%closest grain boundary point. 
%INPUTS 
%   Apos: n by 3 array of positions of particles not within any grain
%   GBpos: n by 3 array of grain boundary points
%   Rcut: maximum distance a grain boundary point can be away from the
%   grain boundary.
%OUTPUTS: 
%   GBpts: n by 3 array of grain boundary points
%HISTORY: 
%Written by Nick Orr 13/ 07/ 21
nApos = numel(Apos(:,1)); 
nGBpos = numel(GBpos(:,1)); 

%nApos = 3
%nGBpos = 4

A = 1:nApos;
G = 1:nGBpos;

A1 = repmat(A',nGBpos, 1);
G1 = repmat(G, nApos, 1);
G2 = reshape(G1, numel(G1), 1);

dists = sqrt(sum((Apos(A1, :) - GBpos(G2, :)).^2, 2)); 
GBptsID = A1(dists < Rcut); 
GBptsIDU = unique(GBptsID); 
GBpts = Apos(GBptsIDU, :); 

end

