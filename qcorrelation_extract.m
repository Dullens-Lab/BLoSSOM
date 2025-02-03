function [ XYZC, qcor_xyz_ARR_OUT ] = qcorrelation_extract(dataArr,L,ll,NN,Sthresh, Nthresh) 
%NAME: qcorrelation_order
%FUNCTION: measures the correlation of the bond orientational order of 
%a particle with its neighbours.
%INPUTS:
    %dataArr - 3 by n array of xyz coordiantes as columsn and each row is a
    %new particle
    %L is the size of the box
    %ll is the order (usually 4 or 6)
    %NN is the nearest neighbour structure from NNsearch
    %Sthresh is the correlation threshold, commonly 0.7
    %Nthresh is the number of neighbours which must have a correlation
    %threshold above Sthresh for the particle to be considered crystalline.
    %
%OUTPUTS: 
    %XYZC a structure with 4 fields: pos, NC, Npos and ID. pos contains the
    %position of the crystalline particle. NC gives the number of
    %neighbours with correlated crystalline order. Npos gives the positions
    %of the neighbours with correlated order. ID gives the id of the
    %crystalline particle from its row in dataArr.
    %Obtain the positions of the crystalline points as an array:
    %   pos = zeros(numel(cor), 3); 
    %   for A = 1:numel(cor)
    %        pos(A, :) = cor(A).pos;
    %   end 
%History: 
%April 2019 Nick Orr, changed the input from vector to n by 3 array.
%Verisons of qcorrelation_order not in Orientation Packages will require
%the use of a data vector input.
%
%September 2024 Nick Orr (in Nijmegen). Modified from qcorrelation_order so
%that the values of each d6 (bond coherence) and the interparticle
%coordinate associated with the d6 value is recorded. 



data = reshape(dataArr', (numel(dataArr(:,1))*3), 1);
CCount = 1; 
%crash
[q, qc] = BOO('q', data, L, ll, NN, 1);
np = numel(data)/ 3;
nbs = NN.LIST;
%crash
%qcor_xyz_cell = cell(1,np);
qcor_xyz_ARR = zeros(np*12, 4);
count = 1;

for A = 1:np
    %coordinates of the central particle of the NNC 
    central_xyz = data((3*A) -2:(3*A));
    Anbs = nbs(A,:);
    AnbsI = Anbs(Anbs ~= 0); % the neighbouring particles
    nAnbsI = numel(AnbsI);
    %Aqc = [real(qc{A}), imag(qc{A})];
    Aqc = qc{A};
    
    normA = 0;
        for C = 2:ll+1
            %DOT(C) = real(Bqc(C, :)* conj(Aqc(C, :))); %/(hypot(Bqc(C,1), Bqc(C,2))*hypot(Aqc(C,1), Aqc(C,2)));
            normA = normA + 2 * real(Aqc(C) *conj(Aqc(C)));
        end
        normA = normA + real(Aqc(1)*conj(Aqc(1)));
    normA = sqrt(normA);
    Aqc = Aqc / normA;
    
    DOT = zeros(nAnbsI, 1); 
    for B = 1:nAnbsI
        %nqc = numel(qc{AnbsI(B)});
        %Bqc = [real(qc{AnbsI(B)}), imag(qc{AnbsI(B)})];
        Bqc = qc{AnbsI(B)};
        
    normB = 0;
        for C = 2:ll+1
            %DOT(C) = real(Bqc(C, :)* conj(Aqc(C, :))); %/(hypot(Bqc(C,1), Bqc(C,2))*hypot(Aqc(C,1), Aqc(C,2)));
            normB = normB + 2 * real(Bqc(C) * conj(Bqc(C)));
        end
        normB = normB + real(Bqc(1)*conj(Bqc(1)));
    normB = sqrt(normB);
    Bqc = Bqc / normB;
    
        
        SC = 0;
        for C = 2:ll+1
            %DOT(C) = real(Bqc(C, :)* conj(Aqc(C, :))); %/(hypot(Bqc(C,1), Bqc(C,2))*hypot(Aqc(C,1), Aqc(C,2)));
            SC = SC + 2 * real(Aqc(C) * conj(Bqc(C)));
        end
        SC = SC + real(Aqc(1)*conj(Bqc(1)));
        %sumDOT = sum(DOT(2:nqc)) * 2 + DOT(1); %/nqc
        %DOT(B) = dot(Bqc, conj(Aqc));
        DOT(B) = SC; % These are all of the bond coherences between neighbours.

    end 
   
    %getting the poition of each of the neighbours 
    nbsxyz = zeros(nAnbsI, 3);
    intermediate_xyz = zeros(nAnbsI, 3); %The coordinates of the intermediate points. 
    for n = 1:nAnbsI
        nbxyz = data((3*AnbsI(n)) -2:(3*AnbsI(n)));
        nbsxyz(n,:) = nbxyz'; 
        intermediate_xyz(n,:) = (central_xyz' + nbxyz')/2;
    end
    %so for each of these neighbours I calculate the intermediate position to
    %the central particle. 
    qcor_xyz_ARR(count:1:(count + numel(intermediate_xyz(:,1))-1),:) = [intermediate_xyz, DOT];
    count = count + numel(intermediate_xyz(:,1));
    %qcor_xyz_cell{A} = [intermediate_xyz, DOT]; % there should be as many values as there are neighbours here 
    %crash
    
    %O = ones(nAnbsI);
    %fC1 = find(real(DOT) > Sthresh);
    fC1 = find(DOT >= Sthresh);
    nC1 = numel(fC1);
    if nC1 >= Nthresh
        [~, q_order] = sort(DOT,'descend'); 
        AnbsI = AnbsI(q_order); 
        Npos = zeros(nAnbsI, 3); 
        for NCount = 1:nAnbsI
            Npos(NCount, :) = data((3*AnbsI(NCount) -2 :(3*AnbsI(NCount)))); 
        end 
        XYZC(CCount) = struct('pos', (data((3*A) -2:(3*A)))', 'Nc', nC1, 'Npos', Npos, 'ID', A );
        CCount = CCount + 1; 
    end
    clear DOT
end 
qcor_xyz_ARR_OUT = qcor_xyz_ARR(qcor_xyz_ARR(:,1) ~= 0,:);
end