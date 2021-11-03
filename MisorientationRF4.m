function [minRF] = MisorientationRF4( rotation_matrix1, rotation_matrix2, p_group1, p_group2); 
%NAME: MisorientationRF4
%FUNCTION: 
%finds dosorientation in RF space from rotation matrices of the two
%orientattions.
%p_group1 and p_group2 are the symmetry identifiers 1 = O, 3 = D6, 2 = D3
%minRF n by 3 array
form = 3; %only coding for RF
CRO_Oh = [1, 0, 0; 0, 1, 0; 0, 0, 1;
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

rtt = sqrt(3)/2;
CRO_D3h = [1,0,0;0,1,0;0,0,1;    
    -0.5,rtt,0;-rtt,-0.5,0;0,0,1;    
    -0.5,-rtt,0;rtt,-0.5,0;0,0,1;    
    1,0,0;0,-1,0;0,0,-1;    
    -0.5,rtt,0;rtt,0.5,0;0,0,-1;    
    -0.5,-rtt,0;-rtt,0.5,0;0,0,-1];
CRO_D6h = [1,0,0;0,1,0;0,0,1;    %with more symmetry so hcp has one orientation
    -0.5,rtt,0;-rtt,-0.5,0;0,0,1;    
    -0.5,-rtt,0;rtt,-0.5,0;0,0,1;    
    1,0,0;0,-1,0;0,0,-1;    
    -0.5,rtt,0;rtt,0.5,0;0,0,-1;    
    -0.5,-rtt,0;-rtt,0.5,0;0,0,-1;
    0.5,-rtt,0;rtt,0.5,0;0,0,1;
    -1,0,0;0,-1,0;0,0,1;
    0.5,rtt,0;-rtt,0.5,0;0,0,1;
    0.5,rtt,0;rtt,-0.5,0;0,0,-1;
    -1,0,0;0,1,0;0,0,-1;
    0.5,-rtt,0;-rtt,-0.5,0;0,0,-1];
M = rotation_matrix1 * rotation_matrix2';
if p_group1 == 1
    S1 = CRO_Oh;
    n1 = 24;
elseif p_group1 == 2
    S1 = CRO_D3h;
    n1 = 6; 
elseif p_group1 == 3
    S1 = CRO_D6h;
    n1 = 12; 
end 
if p_group2 == 1 
    S2 = CRO_Oh;
    n2 = 24; 
elseif p_group2 == 2
    S2 = CRO_D3h;
    n2 = 6; 
elseif p_group2 == 3
    S2 = CRO_D6h;
    n2 = 12; 
end 
%D = 1; 
E = 1; 
for A = 1:3:n1*3
    S1A = S1(A:A+2, :); 
    for B = 1:3:n2*3
        S2B = S2(B:B+2, :); 
        if form == 3; 
            RFs(E, :) = Rot_to_AngAx(S1A*M*S2B, 1);
        elseif form == 2; 
            RFs(E, :) = Rot_to_AngAx(S1A*M*S2B*S1A', 0);
        end
        E = E + 1;
    end 
end 
if p_group1 == p_group2
    RFs = [RFs; -RFs]; 
end
%whos RFs
c = 1;
if p_group1 == 1 && p_group2 == 1
    for a = 1:numel(RFs(:,1)); 
        RFa = RFs(a,:); %applying inequality conditions from Heinz and Neumann.
        %if RFa(1) <= (sqrt(2) - 1) || -RFa(1) <= (sqrt(2) - 1)
            %if RFa(2) <= (sqrt(2) - 1) || -RFa(2) <= (sqrt(2) - 1)
                %if RFa(3) <= (sqrt(2) - 1) || -RFa(3) <= (sqrt(2) - 1)
                    %if  RFa(1) + RFa(2) + RFa(3) <= 1 ||...
                        %-RFa(1) - RFa(2) - RFa(3) <= 1 ||...
                        %RFa(1) + RFa(2) - RFa(3) <= 1 ||...
                        %RFa(1) - RFa(2) + RFa(3) <= 1 ||...
                        %RFa(1) - RFa(2) - RFa(3) <= 1 ||...
                        %-RFa(1) + RFa(2) + RFa(3) <= 1 ||...
                        %-RFa(1) + RFa(2) - RFa(3) <= 1 ||...
                        %-RFa(1) - RFa(2) + RFa(3) <= 1 
                        if RFa(1) >= RFa(2) && RFa(2) >= RFa(3) && RFa(3) >= 0
                            RFdis(c, :) = RFa;
                            c = c + 1;
                        end 
                    %end  
                %end
            %end
        %end 
    end
elseif p_group1 == 3 && p_group2 == 3
    for a = 1:numel(RFs(:,1)); 
        RFa = RFs(a,:);
        if RFa(2) >= 0 && RFa(2) <= RFa(1)*(1/sqrt(3)) && RFa(3) >= 0
            RFdis(c, :) = RFa;
            c = c + 1;
        end
    end
    
elseif p_group1 ~= p_group2
aa = (sqrt(3) - 1)/(sqrt(3) + 1);
bb = (3 - sqrt(3))/(sqrt(3) + 1);
    if p_group1 == 3 && p_group2 == 1
        for a = 1:numel(RFs(:,1)); 
            RFa = RFs(a,:);
            if RFa(1) <= (sqrt(2) - 1) && RFa(2) <= (sqrt(2) - 1) &&...
                (2 * sqrt(2) - sqrt(3) - 1)/(sqrt(3)-1) >= RFa(3) &&...
                bb >= aa * RFa(1) &&...
                bb >= 1 * RFa(2) &&...
                bb >= aa * RFa(3) &&...
                bb >= 1 * RFa(1) &&...
                bb >= aa * RFa(2) &&...
                bb >= -aa * RFa(3) &&... 
                RFa(2) >= 0 && RFa(3) >= 0 
                RFdis(c, :) = RFa;
                c = c + 1;
            end
        end
    elseif p_group1 == 1 && p_group2 == 3
        disp('please use convention low symmetry to high symmetry')
    end
end
%RFdis
for a = 1:numel(RFdis(:,1))
        modRF(a) = sqrt(sum(RFdis(a,:).^2));   
end 
    minRF = RFdis(modRF==min(modRF), :);


end

