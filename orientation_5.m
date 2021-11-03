function [ orientation_array, fitT] = orientation_5( Data, reference, minN, symmineqidx)
%NAME: orientation
%FUNCTION: 
%To calculate the orientation of grains within a ploycrystalline crystal. 
%INPUTS: 
%   Data - the data with nearest neighbours all of one symmetry. This
%   should be in a structure with the same fields as the output from
%   q4w4sort.
%   reference - The reference structures. the index
%   minN - the minimum number of neighbours a partcile will have.
%   goodenough - mean squared displacement good enough for registration.
%OUTPUT: 
%   orientation_array - an array containg the orientation matrices. Every
%   3 rows the next particle's orientation matrix is found. 
%   symmcode - the symmetry code for the orientation. 
%   fit - the fit of the point set registration, the fit wil be -1 if the
%   registration has failed.
%   minDiffSqARR - the sum of the minimum specimen- reference differences 
%   of the registration permutations trialed in the registration stage.
%   Each new row is are the numbers for each cluster
%HISTORY: 
%V_0 written by Nick Orr in January 2018. 
%V_1 modified by Nick Orr to remove the goodenough input. 
%V_2 Added parellisation. Nick Orr April 2018.
%V_4 Tweaked number of combinations of three points to oonly run only
%symmetry inequivalent, Nick Orr March 2019. 
%get the dimension from the reference structures.
%V5 Added parallelisation. Nick Orr April 2019. 

dimension = numel(reference(1,:)); 
maxN = numel(reference(:,1)); 
n = numel(Data); 
orientation_array = zeros(dimension*n, dimension); 
fitT = zeros(1,n); 

indecesT = generate_indeces(minN, maxN, dimension);
U = cell(n, 1); %V5 addition April 2019
parfor a = 1:n %V5 addition April 2019
    cent = Data(a).pos; 
    Neighbours = Data(a).Npos; 
    rcent = repmat(cent, numel(Neighbours(:,1)), 1); 
    rpoints = Neighbours - rcent;
    npointsA = numel(Neighbours(:,1)); 
    index = npointsA - minN + 1; 
    indeces = indecesT(index); 
    [Ar, Br, fit, ~] = register_points_sumoutput_4_loop(rpoints, reference, dimension, indeces, symmineqidx);
    sA = rpoints(Ar, :); 
    sB = reference(Br, :); 
    fitT(a) = fit; 
    U{a} = Kabsch_centroid(sA', sB');
    %minDiffSqARR(a,1:numel(minDiffSq)) = minDiffSq'; 
    if mod(a/1, 100) == 0  
         feedback = strcat('Run_',num2str(a),'_of_', num2str(n), '_particles');
         disp(feedback);
    end 
end 
for a = 1:numel(U)
    orientation_array(a*dimension - dimension +1: a*dimension, :) = U{a}; 
end
end 



