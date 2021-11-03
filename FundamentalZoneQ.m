function [qFZ, RFFZ] = FundamentalZoneQ(qARR,qoperators)
%NAME: FundamentalZoneQ
%FUNCTION: Find quaternions in the fundamental zone of orientation 
%INPUTS: 
%   qARR - n by 4 array containing quaternions
%   qoperators - the symmetry opertaors as an m by 4 array. 
%OUTPUTS: 
%   qFZ - n by 4 arry of quaternions in the fundamental zone of orientation
%   RFFZ - n by 3 array of RF vectors. 
%HISTORY: 
%Written by Nick Orr 03/04/20. The motiviation of this is two fold: 
%1/ to visualise orientations in fundemantal zone (must convert to RF)
%2/ to calculate the average orientation of a grain. Simply taking the
%average orientation of a mix of all of the symmetry equivalent
%orientations is not sufficient. 
qFZ = zeros(numel(qARR(:,1)), 4); 
RFFZ = zeros(numel(qARR(:,1)), 3); 
for a = 1:numel(qARR(:,1))
    qa = qARR(a,:); 
    qS = quatmultiply(qoperators, qa); 
    RFS = qS(:,2:4)./qS(:,1); 
    mod = sum(RFS.^2, 2); 
    minmod = min(mod); 
    qFZ(a,:) = qS(mod == minmod, :); 
    RFFZ(a,:) = RFS(mod == minmod, :); 
end

end

