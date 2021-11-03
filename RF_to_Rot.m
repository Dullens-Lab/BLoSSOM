function [RotMat] = RF_to_Rot( RFvector )
%NAME:RF_to_Rot
%FUNCTION: Turn a RF vector into a rotation matrix
%History: written by Nick Orr druing pt II 2016-17
if sum(RFvector.^2) ~= 0 
    tanpsiby2 = sqrt(sum(RFvector.^2));
    UV = RFvector./tanpsiby2;
    sqrt(sum(UV.^2));
    psi = atan(tanpsiby2) * 2;
    RotMat = vrrotvec2mat([UV,psi])';
else
    RotMat = eye(3); 
end
%RotMat = [((UV(1)^2) * (1-cos(psi))+cos(psi)), (UV(1) * UV(2) * (1-cos(psi)) - UV(3) * sin(psi)),...
    %(UV(1) * UV(3) * (1-cos(psi)) + UV(2) * sin(psi));...
    %(UV(1) * UV(2) * (1-cos(psi)) + UV(3) * sin(psi)), ((UV(2)^2) * (1 - cos(psi)) + cos(psi)),...
    %(UV(2) * UV(3) * (1-cos(psi)) - UV(1) * sin(psi));...
    %(UV(1) * UV(3) * (1-cos(psi)) - UV(2) * sin(psi)), (UV(2) * UV(3) * (1-cos(psi)) +UV(1) * sin(psi)),...
    %((UV(3)^2) * (1-cos(psi)) + cos(psi))]';
%dont know why this doesnt quite work, there might be a typo in the book
end    