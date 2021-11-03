function [ avRF ] = RF_av_plus( RF )
%A function to average multiple Rodrigues Frank orientations. Uses quatrenion
%averaging, see below.
%RF is a concatenated n by 3 array of RF vectors.
%HISTORY: 
%Modified from RF_av 20/03/19.
Q = zeros(numel(RF(:,1)), 4);
for a = 1:numel(RF(:,1))
    rot = RF_to_Rot(RF(a,:));
    Q(a,:) = rotm2quat(rot);
    %crash
end 
avQ = avg_quaternion_markley(Q); 
avrot = quat2rotm(avQ'); 
avRF = Rot_to_AngAx(avrot, 1); 

end


function [Qavg]=avg_quaternion_markley(Q)
% by Tolga Birdal
% Q is an Mx4 matrix of quaternions. Qavg is the average quaternion
% Based on 
% Markley, F. Landis, Yang Cheng, John Lucas Crassidis, and Yaakov Oshman. 
% "Averaging quaternions." Journal of Guidance, Control, and Dynamics 30, 
% no. 4 (2007): 1193-1197.
% Form the symmetric accumulator matrix
A = zeros(4,4);
M = size(Q,1);

for i=1:M
    q = Q(i,:)';
    if(q(1)<0) % handle the antipodal configuration
		q = -q;
	end
    A = q*q'+A; % rank 1 update
end

% scale
A=(1.0/M)*A;

% Get the eigenvector corresponding to largest eigen value
[Qavg, Eval] = eigs(A,1);

end