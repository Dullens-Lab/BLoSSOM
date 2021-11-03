function [BtotROI,BtotconROI] = GB_MROI(Btot,BtotCon, minmaxR1, minmaxR2, minmaxR3)
%NAME: GB_MROI, short for grain boundary misorientation region of interest.
%FUNCTION: to select grain boundaries within a misorientation region of
%interest. 
%INPUTS: 
%Btot, xyzR123 position and misorientation in columns 1:3 and 4:6
%respectively, columns 7:8 contain the symmetries between which the
%misorientation was measured. 
%Btotcon, the conectivity between the Btot points
%minmaxR1/2/3, are the bounds for the ROI, each as a 1 by 2 vector. [min,
%max]. 
%HISTORY: 
%Written by Nick Orr May 2020

%for Btot we need a refercne Id vector. 
RefID = 1:numel(Btot(:,1));
RefID = RefID';
ROI = Btot(:,4) > minmaxR1(1) & Btot(:,4) < minmaxR1(2)...
    & Btot(:,5) > minmaxR2(1) & Btot(:,5) < minmaxR2(2)...
    & Btot(:,6) > minmaxR3(1) & Btot(:,6) < minmaxR3(2);
BtotROI = Btot(ROI, :); 
RefIDselect = RefID(ROI, :); 

Btotconselect  = BtotCon; 
renamearray = zeros(numel(Btotconselect(:,1)), 2); 
for a = 1:numel(RefIDselect)
    Refa = RefIDselect(a);
    B1 = Btotconselect(:,1) == Refa;
    B2 = Btotconselect(:,2) == Refa;
    Btotconselect(B1, 1) = a; 
    Btotconselect(B2, 2) = a; 
    renamearray(B1, 1) = 1; 
    renamearray(B2, 2) = 1; 
end 
%crash 
BtotconROI = Btotconselect(renamearray(:,1) == 1 & renamearray(:,2) == 1, :); 

end

%for all of the reference IDs find where they exist in the connection matrix and change 
%their value to where the ID of the reference ID. Delete anything that wasnt changed  