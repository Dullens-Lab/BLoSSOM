function [ idx_symm ] = symmetry_unique_pselect( points, symmetry_op )
%NAME: symmetry_unique_pselect
%FUNCTION: find the symmetry unique selections of D = 3 points from the
%reference structure given in points. D is the dimension of the points. 
%
%INPUTS: 
%   points - n by D array of points. Your reference NNC
%   symmetry_op - the symmetry operators in an (3*n)*3 matrix. Vertically 
%   concatenated rotation matrices. In this version I have hard coded the O
%   and D3 rotational symmetry operators for FCC and HCP NNCs,
%   respectively. You may want to add others if you need to. 
%OUTPUTS: 
%   The ids to compare with a specimen. horizontally concatenated.
%   %vertically? This contains the symmetry inequivalent permutations of
%   D points from the reference. 
%
%NOTES: 
%using symmetry operators to find these equivalent selection is the
%quickest brute force way to do this. 
%number 
%HISTORY:
%modified 11/03/20 Nick Orr 
%changed the generate indeces line so that the function recognises the
%number of points in the input 'points'.
if symmetry_op == 1
    symmetry_op = [1, 0, 0; 0, 1, 0; 0, 0, 1;
        0,0,-1;0,-1,0;-1,0,0;
        0,0,-1;0,1,0;1,0,0;
        -1,0,0;0,1,0;0,0,-1;
        0,0,1;0,1,0;-1,0,0;
        1,0,0;0,0,-1;0,1,0;
        1,0,0;0,-1,0;0,0,-1;
        1,0,0;0,0,1;0,-1,0;
        0,-1,0;1,0,0;0,0,1;
        -1,0,0;0,-1,0;0,0,1;
        0,1,0;-1,0,0;0,0,1;
        0,0,1;1,0,0;0,1,0;
        0,1,0;0,0,1;1,0,0;
        0,0,-1;-1,0,0;0,1,0;
        0,-1,0;0,0,1;-1,0,0;
        0,1,0;0,0,-1;-1,0,0;
        0,0,-1;1,0,0;0,-1,0;
        0,0,1;-1,0,0;0,-1,0;
        0,-1,0;0,0,-1;1,0,0;
        0,1,0;1,0,0;0,0,-1;
        -1,0,0;0,0,1;0,1,0;
        0,0,1;0,-1,0;1,0,0;
        0,-1,0;-1,0,0;0,0,-1;
        -1,0,0;0,0,-1;0,-1,0];
elseif symmetry_op == 2
    rtt = sqrt(3)/2;
    symmetry_op = [1,0,0;0,1,0;0,0,1;    
        -0.5,rtt,0;-rtt,-0.5,0;0,0,1;    
        -0.5,-rtt,0;rtt,-0.5,0;0,0,1;    
        1,0,0;0,-1,0;0,0,-1;    
        -0.5,rtt,0;rtt,0.5,0;0,0,-1;    
        -0.5,-rtt,0;-rtt,0.5,0;0,0,-1];
end
D = numel(points(1, :)); 
symmineq = zeros(0,D); 
%longpoints = reshape(points', 3*12, 1);
np = numel(points(:,1)); 
indeces = generate_indeces(np, np, D); %Nick Orr 11/03/20
indecesP = indeces.indecesB; 
for a = 1:D:numel(indecesP)
    pointsa = points(indecesP(a: (a+(D-1))), :);
    longpointsa = reshape(pointsa', numel(pointsa(:,1))*D, 1);
    if numel(symmineq) == 0; 
        symmineq = longpointsa; 
    else 
        pointsarARR = zeros(numel(symmetry_op(:,1)), D);
        for b = 1:D:numel(symmetry_op(:,1))
            opb = symmetry_op(b:b+(D-1), :);
            for c = 1:D
                pointsar(c, :) = (opb * pointsa(c, :)')';
            end
            pointsarARR(b:b+(D-1), :) = pointsar;       
        end
        longpointsar = reshape(pointsarARR', numel(pointsarARR(:,1))*D, 1);
        %are there any symmetry equivs that we already have in our array? 
        replongpointsar = repmat(longpointsar, numel(symmineq)/D^2, 1);
        repsymmineq = repmat(symmineq, numel(longpointsar)/D^2, 1);
        diff = repsymmineq - replongpointsar; 
        diffr = reshape(diff',D, numel(diff)/D)';
        mods = sum(diffr'.^2)'; 
        modstripl = reshape(mods',D, numel(mods)/D)';
        summodstripl = sum(modstripl,2);
        mods_sort = sort(summodstripl); 
        if mods_sort(1) > 1e-20 %&& mods_sort(2) > 1e-20 && mods_sort(3) > 1e-20
            symmineq = [symmineq; longpointsa]; %a lon vector 123123 etc 
        end
%         if numel(symmineq(:,1)) > 225 
%             crash
%         end
    end
end 
symmineqr_i = reshape(symmineq',D, numel(symmineq)/D)';
idx_symm = xyz_index(symmineqr_i, points);
end 

function [ idx_v ] = xyz_index( fccsub_array, fccxyz )
%NAME: fccxyz_index 
%FUNCTION: find the index of a point from the fccxyz array. 
%INPUT:     
%   fccsub_array - an n by 3 array containing n points of the fccxyz array.
%   
%OUTPUT:    
%   idx_v - long vector where each entry corresponds to the index of the
%   fccsub_array row (point); 
%HISTORY:
%Written by Nick Orr. March 2019. 
D = numel(fccxyz(1, :)); 
fccxyzr = reshape(fccxyz', numel(fccxyz), 1);
%fccsubr = reshape(fccsub_array', numel(fccsub_array), 1);
idx_v = zeros(numel(fccsub_array(:,1)), 1);
for a = 1:numel(fccsub_array(:,1))
    longfccsuba = repmat(fccsub_array(a,:)', numel(fccxyz(:,1)), 1);
    diffs = longfccsuba - fccxyzr; 
    diffsr = reshape(diffs',D, numel(diffs)/D)';
    mods = sum(diffsr.^2,2)'; 
    idx_v(a) = find((mods < 1e-20)); 
end 


end

