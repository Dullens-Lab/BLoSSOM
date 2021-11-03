function [ resmat ] = Rot_to_AngAx( rotation_matrix, RF )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
g = rotation_matrix;
phi = acos((g(1,1) + g(2,2) + g(3,3) -1) /2);
if phi == 180
    r1 = sqrt(g(1,1) + 1);
    r2 = sqrt(g(2,2) + 1);
    r3 = sqrt(g(3,3) + 1);
else
r1 = g(2,3) - g(3,2);
r2 = g(3,1) - g(1,3);
r3 = g(1,2) - g(2,1);

mr = sqrt(r1^2 + r2^2 +r3^2);

r1n = r1/mr;
r2n = r2/mr;
r3n = r3/mr;
resmat = [phi, r1n, r2n, r3n];
if RF == 1 
    R1 = tan(phi/2)*r1n;
    R2 = tan(phi/2)*r2n;
    R3 = tan(phi/2)*r3n;
    resmat = [R1, R2, R3];
end 
    
    
end

