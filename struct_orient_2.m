function [ pos, orientation_array_min, fitTmin, structcode, minDiffSqARR2 ] = struct_orient_2( coordinates, NN, refcell, minN, symmineqidxCell)
%NAME: Structure_orient_2
%FUNCTION: 
%Finds the structure and orientation of nearest neighbour clusters.
%INPUTS: 
%   Data - the data with nearest neighbours all of one symmetry. This
%   should be in a structure with the same fields as the output from
%   q4w4sort. needs to be a long structure : 1 bby n where n is the number
%   of nearest neighbour clusters.
%   reference - The reference structures as a cell array.
%   minN - the minimum number of neighbours a partcile will have.
%   symmineqidxCell - The symmety equiv indeces as a cell array.
%OUTPUT: 
%   pos - The positions of all of the coordinates where orientations have
%   been attempted to be assigned. This happens for all NNCs that have 3 or
%   more members. 3D orientation constraint. The 4th column is the ID from
%   the original dataset.
%   orientation_array - an array containg the orientation matrices. Every
%   3 rows the next particle's orientation matrix is found. 
%   symmcode - the symmetry code for the orientation. 
%   fit - the fit of the point set registration, the fit wil be -1 if the
%   registration has failed.
%   minDiffSqARR - the sum of the squared distances between registered
%   points for each registration. 

%HISTORY: 
%orientation functions:
%V_0 written by Nick Orr in January 2018. 
%V_1 modified by Nick Orr to remove the goodenough input. 
%V_2 Added parellisation. Nick Orr April 2018.
%V_4 Tweaked number of combinations of three points to oonly run only
%symmetry inequivalent, Nick Orr March 2019. 
%get the dimension from the reference structures.
%Hijacked orientation_4.m to make the struct_orient code.
%Written by Nick Orr sept 2019. Instead of using Steinhardt bond
%orientation order parameters this function minimises over a set of
%reference structures.
%struct_orient_2. 17 March 2020. Nick Orr. Made it so that structures with
%different numbers of nearest neighbours can be compared.
minDiffSqARRcell = cell(numel(refcell), 1); 
% need to know the number of members of each reference structure. 
Nref = numel(refcell) ;
N = zeros(Nref, 1); 
for a = 1:Nref;
    refa = refcell{a};
    N(a) = numel(refa(:,1)); 
end 
%Now I need to make a list of the NNCs that can be compared for each
%reference structure. i.e. ones that will be compared with reference
%structre 1, 2.. 3 etc ...
%Now data input isnt necessarily going to be from q_corelation, it can just
%be a NNlist no change to long structure. Maybe a use for the ID,
%orientation_4 doesnt seen to use it.

%then you only need to compare Ids that have two or more entries.
if minN < 3
    disp('minN must be >= 3, setting minN = 3')
    minN = 3;  
end
[DataLong, pos] = NNpack(coordinates, NN, minN); %I NEED TO EDIT POS CALC 
DataLISTN = NN.LISTN;
DataLISTN_Min = DataLISTN(DataLISTN >=  minN); %particles with a number of neighbours above the minimum
OrignialIDs_Above = find(DataLISTN >=  minN); % the IDs for the orignial coordinate array.
fitTARR = zeros(numel(DataLISTN_Min), numel(refcell)) - 2; 
orientation_arrayCONC = zeros(numel(DataLISTN_Min)*3, numel(refcell)*3); 

for a = 1:numel(refcell)
    refa = refcell{a}; 
    COMP_IDs = find(DataLISTN_Min <= numel(refa(:,1))); %particles with a number or neighbours above the minimum and below the number in the reference structure. 
    OriginalIDs_AboveAndBelow = OrignialIDs_Above(COMP_IDs);
    DataLongCOMP = DataLong(COMP_IDs);
    symmineqidxa = symmineqidxCell{a};
    [orientation_array, fitT, minDiffSqARR] = orientation_4( DataLongCOMP, refa, minN, symmineqidxa);
    fitTARR(COMP_IDs,a) = fitT; 
    %crash
    %dealing with the ids for the orientation array.
    COMP_IDs3 = 3*COMP_IDs;
    COMP_IDs2 = COMP_IDs3-1;
    COMP_IDs1 = COMP_IDs3-2;
    COMP_IDs3C = [COMP_IDs1, COMP_IDs2, COMP_IDs3];
    COMP_IDs3CV = reshape(COMP_IDs3C', 1, numel(COMP_IDs3C));
    orientation_arrayCONC(COMP_IDs3CV,a*3-2 : a*3) = orientation_array; %Saving the orientation array 
    minDiffSqARRcell{a} = minDiffSqARR;
end
%[~, minfitID] = min(fitTARR'); % gives the id of which symmetry is the best fitting
minfitID = zeros(numel(fitTARR(:,1)), 1);
for a = 1:numel(fitTARR(:,1))
    fitARRa = fitTARR(a,:); 
    fitARRapos = fitARRa(fitARRa >= 0);
    minpos = min(fitARRapos); 
    if numel(minpos) >= 1; 
        fmin = find(fitARRa == minpos); 
        minfitID(a) = fmin(1); %if two mins are the same, select the first one.
    else 
        minfitID(a) = 1; %just a default. 
    end
end
%now I have a list of Ids 
%minfitID = minfitID'; %this gives the symmetry identity
orientation_array_min = zeros(size(orientation_array));
fitTmin = zeros(numel(fitTARR(:,1)),1);
for a = 1:numel(minfitID(:,1))
    minfita = minfitID(a);
    orientation_array_min(a*3 - 2 : a*3,:) = orientation_arrayCONC(a*3 - 2 : a*3, minfita*3 - 2 : minfita*3);
    fitTmin(a) = fitTARR(a,minfita); 
    minDiffSqARR2 = minDiffSqARRcell{minfitID}; 
end 
structcode = minfitID; 
%crash
end

%NEED TO FIGURE OUT the mindiffsrARR stuff 


