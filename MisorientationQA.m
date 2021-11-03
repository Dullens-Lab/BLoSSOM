function [M12] = MisorientationQA(q1, q2, qOps1, qOps2)
%NAME: MisorientationQA
%function: 
%Calculates the smallest magnitude symmetry equivalent misorientation using
%quaterions. 
%INPUT: 
% q1 - the quaternion for orientation 1
% q2 - the quaternion for orientation 2
% qOps1- the symmetry operators for grain 1 
% qOps2 - likewise for grain 2
%History:
%Mostly copied from MisorientationQ, this version should be faster.
M12_proto = quatmultiply(q1, q2);%q1 * conj(q2);
M12_i = quatmultiply(M12_proto, qOps2);
IDs1 = 1:numel(qOps1(:,1)); 
IDs1r = repmat(IDs1, numel(qOps2(:,1)),1);
IDs1s = reshape(IDs1r, numel(IDs1r), 1);
qOps1REP = qOps1(IDs1s, :); 
M12_iREP = repmat(M12_i, numel(qOps1(:,1)), 1);
M12_S = quatmultiply(qOps1REP, M12_iREP);
%so if there are multiple with the same misorientation (which there will
%be) you dont need to select anyone in particluar because we only care about
%the angle of misorientation, not the axis.
M12ang = acos(M12_S(:,1)) * 2; 
M12 =  min(M12ang);

end

