function [ Ar, Br, Fit, minDiffSq ] = register_points_sumoutput_4_loop( A,B, dimension, indeces, symmineqidx)
%NAME: register_points
%FUNCTION: register pointsin A with points in B.
%INPUTS: 
%   A,B Arrays with row vectors. The number of points in A can
%   be smaller than the number of points in B. 
%   dimension -3D points is 3 etc...
%   NB! The centroid of the object is not based on centre of mass of
%   the points so absent points can be accounted for. its is at the
%   orirgin.
%   googenough -An average mean squared displacement fit that is good
%   enough for satisfactory assignment for the user. Only checks once every
%   slection of 3 from A points. 
%OUTPUTS: 
%   Ar,Br -registered a and B points. 
%   Fit, the mean squared distance between the registered points. If Fit is
%   returned as -1, registration failed. Ar and Br will not change A and B.
%   minDiffSq - the sum of the numel(A) shortest distances for the
%   regsitration for each permutation attempted.
%
%DETAILS: 
%The function works by selecting "dimention", d, points from A to
%test. This avoids computing huge numbers of posibilities. Once d points
%are selcted thier optimal rotation to d points in B is calculated for all
%posibilities of d points in B. ICP is performed and the total lsf is
%calculated. lsf is minimised over all d points in B and A. 
%entension into 1D an d 4+D should be straightforward. I cant think of any
%uses currently so I'm not going to code it. 
%History: 
%Written by Nick Orr in Jan 2018.
%Modified from register_points.m by Nick Orr in Jan 2018. 
%Modified ftom register_points_2.m by Nick Orr Jan 2018. This version is
%register_points_2 but with the indeces generation removed. The loop only
%chooses three points from the data. Tests have shown that the cluster is
%either assigned straight away or not at all. 
%Looking to write another one where the indeces for B combinations are
%saved and not generated over and over.
%dbstop if error %debug. 
%v4 - to work with the test version of orientation_4 Nick Orr March 2019. 

%Generation of the indeces
 npointsA = numel(A)/dimension; 
 npointsB = numel(B)/dimension;


indecesA = indeces.indecesA;
%indecesB = indeces.indecesB; 
indecesC = indeces.indecesC; 
indecesD = indeces.indecesD;
indecesB = symmineqidx;

a = 1;  
pointsA =  A(indecesA(a*dimension - (dimension-1): a*dimension), :);
minDiffSq = ones(1,numel(symmineqidx)/ dimension)*-1; 
registration1 = cell(1,nchoosek(npointsA, dimension) * factorial(dimension)); 
for b = 1:numel(symmineqidx)/ dimension
    pointsB = B(indecesB(b*dimension - (dimension-1) : b*dimension), :);
    %need to use kabsch where the centroid is modified. 
    %find rotation. 
    %rotate all of B onto A 
    %find lsf by icp. 
    %store lsf. 
    Rot = Kabsch_centroid(pointsA', pointsB'); 
    RotA = zeros(1, numel(B)); %generating a padded A array.
    for c = 1:numel(A(:,1))
        RotA(c*dimension - (dimension-1) : c*dimension) = Rot*A(c, :)'; %longv      
    end 
    longvectA = zeros(1,npointsB^2*dimension);
    longvectB = longvectA; 
    longvectA = RotA(indecesC(:,1));

    Bv = reshape(B', 1, numel(B)); 
    longvectB = Bv(indecesC(:,2));
    Difflong = (longvectA-longvectB).^2; %vectorising the difference calculation; 
    Diff = reshape(Difflong, dimension, npointsB*npointsA); 
    Diffsquared = sum(Diff);
    [Diff_sorted, Diff_order] = sort(Diffsquared); 
    indecesD_sorted = indecesD(Diff_order(1:npointsA),:); 
    UD_1 = unique(indecesD_sorted(1:npointsA, 1)); 
    UD_2 = unique(indecesD_sorted(1:npointsA, 2)); 
    if numel(UD_1) < npointsA || numel(UD_2) <npointsA %Fail to convincingly reg'
        minDiffSq(b) = -1; %fail detection will be a negative value. 
    else
        minDiffSq(b) = sum(Diff_sorted(1:npointsA)); %sum of the npointsA smallest
        registration1{b} = indecesD_sorted(1:npointsA, :); %The registration indeces
    end 


end 
pminDiffSq = minDiffSq(minDiffSq >= 0); %removing all the failed attampts
pregistration1 = registration1(minDiffSq>= 0); %and un-assigned entries.
sizepminDiffSq = size(pminDiffSq);
if sizepminDiffSq(1) * sizepminDiffSq(2) > 0; 
    mindimselect = min(pminDiffSq);
    registration2i = pregistration1((pminDiffSq == min(pminDiffSq))); 
    registration2 = registration2i{1}; 
    Fit = mindimselect/ npointsA; 
    ArBr = registration2;
    Ar = ArBr(:,1); 
    Br = ArBr(:,2);
else 
    Fit = -1; 
    Ar = 1:npointsA; 
    Br = 1:npointsA; %so dont have to use a padded kabsch calc 
end   
end