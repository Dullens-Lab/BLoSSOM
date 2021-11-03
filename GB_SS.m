function [BptsSS] = GB_SS(BptsRF,s1, s2)
%NAME: GB_SS Grain boundary symmetry select (2nd generation)
%FUNCTION: To find the grain boundary misorientation for a particular
%symmetry combination. 
%INPUTS: 
%   BptsRF, columns 1 to 3 the boundary points, 4:6 the disorientation in
%   RF fundamental zone. 7 and 8 contain the two symmetries between which
%   the misorientation was calculated. 
%   s1 s2, the symmetry identifier for the two grain symmetries. 
%OUTPUTS: 
%   BptsSS the boundary points and misorientations for the selected
%   symmetries. 
%HITSORY:
%Written by Nick Orr April 2020
BS1 = BptsRF(:,7); 
BS2 = BptsRF(:,8); 
if s1 == s2
    BptsSS = BptsRF(BS1 == s1 & BS2 == s1,:); 
elseif s1 ~= s2
    BptsSS =  BptsRF(BS1 ~= BS2,:);
end
end

