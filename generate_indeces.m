function [ indeces ] = generate_indeces( low, high, dimension )
%
%indicesA are selections of 3 points from the number of points in the NNC
%(not used in the current version of the script)
%indicesB generates the permutations of points in a long string. The number
%of permutations of D=dimension points is the number of elements of
%indicesB / D. 
%indicesD is the possible registrations. 
%indicesC is the possible registrations extended to vectorise difference
%calculations between points of D dimension. So that all of the distances
%can be identified for the registration. 
n = high - low + 1; 
indeces = struct;
indeces = repmat(indeces, 1, n); 
npointsB = high; 
for a = 1:n; 
    npointsA = a + low - 1; 
    if dimension == 3; 
        indecesB = zeros(1, factorial(dimension) * nchoosek(npointsB, dimension));
        c = 1:npointsB;
        count = 1; 
        for x = 1:npointsB
            for y = 1:npointsB
                for z = 1:npointsB
                    if ne(x,y) && ne(x,z) && ne(y,z);
                        indecesB(count*3-2 : count*3) = [c(x), c(y), c(z)]; 
                        count = count + 1; 
                    end 
                end 
            end 
        end 
    elseif dimension == 2; 
        indecesB = zeros(1, dimension * nchoosek(npointsB, dimension)); 
        c = 1:npointsB;
        count = 1; 
        for x = 1:npointsB
            for y = 1:npointsB            
                if ne(x,y)
                    indecesB(count*2-1 : count*2) = [c(x), c(y)]; 
                    count = count + 1; 
                end  
            end 
        end 
    end 

    I = 1:npointsA;
    indecesA = nchoosek(I, dimension); 
    indecesA = reshape(indecesA',1, numel(indecesA(:,1)) * numel(indecesA(1,:)));


    indecesC = zeros(npointsA*npointsB*dimension,2); % for the ICP and vectorising diff
    count = 1;
    dim = dimension; %so i can fit subsequent calcs onto one line
    for x = 1:npointsA
        for y = 1:npointsB    
            for z = 0:(dimension-1)
                indecesC(count, :) = [x * dim - (dim - 1) + z, y*dim-(dim-1) + z]; 
                count = count + 1; 
            end
        end 
    end 
    indecesD = zeros(npointsA * npointsB, 2); 
    count = 1; 
    for x = 1:npointsA
        for y = 1:npointsB    
            indecesD(count, :) = [x,y]; 
            count = count + 1; 
        end 
    end
    indeces(a).indecesA = indecesA; 
    indeces(a).indecesB = indecesB; 
    indeces(a).indecesC = indecesC; 
    indeces(a).indecesD = indecesD; 
end

