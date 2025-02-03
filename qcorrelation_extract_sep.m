function [ XYZC, qcor_xyz_ARR_OUT ] = qcorrelation_extract_sep(dataArr, L, ll, NN) 
    %NAME: qcorrelation_extract_sep
%FUNCTION: measures the correlation of the bond orientational order of 
%a particle with its neighbours.
%INPUTS:
    %dataArr - 3 by n array of xyz coordiantes as columsn and each row is a
    %new particle
    %L is the size of the box
    %ll is the order (usually 4 or 6)
    %NN is the nearest neighbour structure from NNsearch
   
%OUTPUTS:
    %XYZC - struct with two fileds: pos
    %                               avg_d6
    
%%%UPDATES:
    % Modified version of Nick's original script
    % Oct 2024 Merin: Get the d6 of each particle by averaging d6
    % over all its neighbours. Output should be xyz coordinates of each
    % particle (pos) and its d6 (avg_d6)


    data = reshape(dataArr', (numel(dataArr(:,1)) * 3), 1);

    np = numel(data) / 3; % Number of particles
    qcor_xyz_ARR = zeros(np, 4); % Initialize 4-column output: [X, Y, Z, avg_d6]
    % Calculate bond orientational order for each particle
    [q, qc] = BOO('q', data, L, ll, NN, 1);
    nbs = NN.LIST; 
    
   
    count = 1;
    
    for A = 1:np
        % Coordinates of the central particle
        central_xyz = data((3*A) - 2 : (3*A))';
        
        % Neighbors of the current particle
        Anbs = nbs(A, :);
        AnbsI = Anbs(Anbs ~= 0); 
        nAnbsI = numel(AnbsI); % Number of neighbors
        
        % Normalized bond orientational order for particle A
        Aqc = qc{A};
        normA = sqrt(sum(2 * real(Aqc(2:ll+1) .* conj(Aqc(2:ll+1)))) + real(Aqc(1) * conj(Aqc(1))));
        Aqc = Aqc / normA;
        
        % Initialize the sum of d6 correlations for particle A
        sum_d6 = 0;
        
        for B = 1:nAnbsI
            % Normalized bond orientational order for neighbor B
            Bqc = qc{AnbsI(B)};
            normB = sqrt(sum(2 * real(Bqc(2:ll+1) .* conj(Bqc(2:ll+1)))) + real(Bqc(1) * conj(Bqc(1))));
            Bqc = Bqc / normB;
            
            % Calculate the bond order correlation (d6) between A and B
            d6 = sum(2 * real(Aqc(2:ll+1) .* conj(Bqc(2:ll+1)))) + real(Aqc(1) * conj(Bqc(1)));
            
            % Add the d6 correlation to the sum for particle A
            sum_d6 = sum_d6 + d6;
        end
        
        % Calculate the average d6 for particle A 
        if nAnbsI > 0
            avg_d6 = sum_d6 / nAnbsI;
        else
            avg_d6 = 0; % If no neighbors, set average to 0
        end
        
        % Store the XYZ coordinates and avg_d6 for particle A
        qcor_xyz_ARR(count, :) = [central_xyz, avg_d6];
        count = count + 1;
    end
    
    
    qcor_xyz_ARR_OUT = qcor_xyz_ARR(qcor_xyz_ARR(:, 1) ~= 0, :);
    
   
    XYZC = struct('pos', qcor_xyz_ARR_OUT(:, 1:3), 'avg_d6', qcor_xyz_ARR_OUT(:, 4));
end
