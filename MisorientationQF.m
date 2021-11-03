function [M12min] = MisorientationQF(q1, q2, qOps1, qOps2, s1, s2)
%NAME: MisorientationQF, Q for quaternion, F for fundamental zone.
%function: 
%Calculates misorientations in the fundamental zone of misorientation as
%described by Heinz and Neumann in their 1991 paper. It Only works for O
%and D6 rotational symmetry.
%INPUT: 
% q1 - the quaternion for orientation 1
% q2 - the CONJUGATE quaternion for orientation 2
% qOps1, qOps2 - the symmetry operators of points 1 and 2 respectively.
% Note the convention of 1 having lower symmetry than 2. 
% s1, s2 - The symmetry code for points 1 and 2 respectively. This is to
% allow for the correct standard sterographic volume to be selected for the
% fundamental zone. NB 1 for O, 2 for D6
%OUTPUT: 
%   M12min - the misorientation between orientations 1 and 2 exfressed as an
%   RF vector within the fundamental zone of misorientation.
%History:
%Written by Nick Orr 1 04 20. 
%hijack of some code from misorientationQA. 
M12_proto = quatmultiply(q1, q2);%q1 * conj(q2);
M12_i = quatmultiply(M12_proto, qOps2);
IDs1 = 1:numel(qOps1(:,1)); 
IDs1r = repmat(IDs1, numel(qOps2(:,1)),1);
IDs1s = reshape(IDs1r, numel(IDs1r), 1);
qOps1REP = qOps1(IDs1s, :); 
M12_iREP = repmat(M12_i, numel(qOps1(:,1)), 1);
M12_S = quatmultiply(qOps1REP, M12_iREP);
%convert from quaternions to RF:
M12RFi = M12_S(:,2:4)./M12_S(:,1);

%standard stereographic volume inequaliteis, H&N 1991.
if s1 == 1 && s2 == s1
    M12RF = [M12RFi; -M12RFi]; % forwards and backwards relation is the same 
    M12 = M12RF(M12RF(:,1) >= M12RF(:,2) & M12RF(:,2) >= M12RF(:,3)...
        & M12RF(:,3) >= 0, :); 
end

if s1 == 2 && s2 == s1
    M12RF = [M12RFi; -M12RFi]; % forwards and backwards relation is the same 
    M12 = M12RF(M12RF(:,2) >= 0 & M12RF(:,2) <= M12RF(:,1) * (1/sqrt(3))...
        & M12RF(:,3) >= 0, :); 
end

if s1 == 2 && s2 == 1
    aa = (sqrt(3) - 1)/(sqrt(3) + 1);
    bb = (3 - sqrt(3))/(sqrt(3) + 1);
    M12RF = -M12RFi; % forwards and backwards relation is not the same 
    %the -ve sign is there to fit the H and N results, somewhere along the
    %way I have reversed the misorientation. It is corrected here. 
    M12 = M12RF(M12RF(:,1) <= (sqrt(2) - 1) & ...
        M12RF(:,2) <= (sqrt(2) - 1) &...
        (2 * sqrt(2) - sqrt(3) - 1)/(sqrt(3)-1) >= M12RF(:,3) &...
        bb >= aa * M12RF(:,1) &...
        bb >= 1 * M12RF(:,2) &...
        bb >= aa * M12RF(:,3) &...
        bb >= 1 * M12RF(:,1) &...
        bb >= aa * M12RF(:,2) &...
        bb >= -aa * M12RF(:,3) &... 
        M12RF(:,2) >= 0 & M12RF(:,3) >= 0, :) ; 
end

M12mod = sqrt(sum(M12.^2,2));
minmod = min(M12mod);
M12mini = M12(M12mod == minmod, :); 
M12min = M12mini(1,:); %incase there are some multiple equivalents, there 
%will be if one of your orientations is equal to your reference. 
%crash 
end

