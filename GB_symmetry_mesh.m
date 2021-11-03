function [ Bxyz_mid, connect_select_mid, C12 ] = GB_symmetry_mesh(GBptsRF, connectivity, select1, select2 )
%NAME: GB_symmetry_mesh
%FUNCTION: decompose total GB meshes into their symmetry pairs. ONLY WORKS
%WITH TWO KINDS OF SYMMETRIES. Any more must be done manually. The function
%also calculates the misorientation associated with each mesh edge.
%INPUTS: 
    %GBptsRF - n by 8 array with x y z R1 R2 R3 s1 s2 columns. 
    %connectivity - m by 4 array with m being the number of connections, i.e
    %mesh edges. The first two columns conatain ids of connected BptsRF
    %points. The second two shows the symmetry of the two grains that were
    %connected accross the boundary. 
%OUTPUTS: 
    %Bxyz_mid, the GBptsRF vertically concatenated with the midpoints of
    %the mesh edges given by the 'connectivity' input.
    %connect_select_mid - the ids of GBptsRF that are connected by a line in
    %the GBmesh
    %C12 - the colour for the line in an n by 3 array. 
%HISTORY: 
    %written Nick Orr, April 2019. 
    %Hijacked from GB_symmetry_select. Averaging of misorientations is to
    %be avoided due to issues with periodic boundaries
s1 = connectivity(:,3); %symmetry one
s2 = connectivity(:,4); %symmetry two 
if select1 ~= select2
    s12 = s1 ~= s2;
else
    s12 = s1 == s2 & s1 == select1;
end 
connect_select12 = connectivity(s12, :);
if numel(connect_select12) > 0  %Nick Orr March08 2021
    
    %connect_select12 = GB_SS(GBptsRF, select1, select2); 

    if select1 == select2 && select1 == 1
        Angmax = tand(0.5 * 62.8); %ref Heinz & Neumann
    elseif select1 ~= select2
        Angmax = tand(56.60/2); %ref Heinz & Neumann
    elseif select1 == select2 && select1 == 2
        Angmax = tand(0.5 * 93.84); %ref Heinz & Neumann
    end

    Mid = zeros(numel(connect_select12(:,1)), 3); %an array containing the midpoints 
    %of the Mesh edges.
    for a = 1:numel(connect_select12(:,1))
        xyz1 = GBptsRF(connect_select12(a,1), 1:3); 
        xyz2 = GBptsRF(connect_select12(a,2), 1:3); 
        Mid(a, 1:3) = (xyz1 + xyz2)/ 2;
    end 
    BptsRF_1 = GBptsRF(connect_select12(:,1), :);
    BptsRF_2 = GBptsRF(connect_select12(:,2), :);
    BptsRF_12 = [BptsRF_1;BptsRF_2]; 
    pre_connect_select_mid = [connect_select12(:,1);connect_select12(:,2)];
    pre_connect_select_mid2 = [[1:numel(connect_select12(:,1))]'; [1:numel(connect_select12(:,1))]'];
    connect_select_mid = [pre_connect_select_mid, pre_connect_select_mid2 + numel(GBptsRF(:,1))];
    Bxyz_mid = [GBptsRF(:,1:3); Mid]; %an array with the Boundary mesh corners and
    %midpoints in the the same array. all of the Mids follow the midpoints.
    M12_t = [BptsRF_12, [1:numel(BptsRF_12(:,1))]']; % adding numbers on the end to be able to recover
    %the orieginal order. The reason all of BptsRF is put into angle colour is
    %thats just how i wrote Angle colour and i dont want to cnage it just for
    %this algorithm. 
    %crash
    [M12_2, colors_total] = Angle_Colour(M12_t, Angmax, 1000, 0, 'xk', 5); 
    [~, so] = sort(M12_2(:,9)); %getting the orignal order back as Angle color removes the ordering
    %BptsRF_3 = BptsRF_2(so,:); 
    C12 = colors_total(so,:);
    %Disorientation = M12_t(:,4:6);
else %Nick Orr March 2021
    disp('!! No boundaries between selected symmetries !!')
    Bxyz_mid = -1;
    connect_select_mid = -1;
    C12 = -1;
end
end

