function [ GC, xyzC ] = Grain_Combine( Grain_s1, Grain_s2, xyzRF1, xyzRF2, s1, s2 )
%stitch the grains of different symmetries together
%INPUTS: 
%   Grains_s1 - the grains of particles of symmetry 1
%   Grains_s2 - the grains of particles of symmetry 2
%   xyzRF1 - the xyz and RF parameters of symmetry 1 n by 6 array, new
%   particle new row.
%   xyzRF2 - the xyz and RF parameters of symmetry 2 n by 6 array, new
%   particle new row.
%   s1 - symmetry 1; (the index of the reference structures)
%   s2 - symmetry 2; (the index of the reference structures)
%OUTPUTS: 
%   GC - the combined grain cell, with indeces updated to reflect the
%   positions of particles in the array xyzC.
%   xyzC - a concatenated array containing xyzRF1 and xyzRF2. 
Ns1 = numel(xyzRF1(:,1));
xyzRFs1 = [xyzRF1, zeros(numel(xyzRF1(:,1)), 1) + s1];
xyzRFs2 = [xyzRF2, zeros(numel(xyzRF2(:,1)), 1) + s2];
for a = 1:numel(Grain_s2)
    Grain_s22{a} = Grain_s2{a} + Ns1; 
end 
%crash
GC = [Grain_s1, Grain_s22];
xyzC = [xyzRFs1; xyzRFs2];
end

