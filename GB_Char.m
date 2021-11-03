function [BptsRF,con] = GB_Char(G, xyzqs, qOperatorCell, minmaxx, minmaxy, minmaxz)
%NAME: GB_Char
%FUNCTION: Characterise the grain boundaries and create the mesh
%connectivity.
%INPUTS: 
%   G - The grains as a cell. Each cell is a new grain 
%   xyzqs - coordinates of the points columns 1:3, quaternion orientation
%   columns 4:7, symmetry identifier column 8. symmetry identifieir
%   corresponds to the index of qOperatorCell
%   qOperatorCell - a cell containing the symmetry operators expressed as
%   quaternions each of the symmetries in column 8 of xyzqs
%   minmaxx/y/z - the min and maximum bound for the grain boundary. This
%   prevents spurious connections.
%OUTPUTS: 
%   BptsRF - xyz coordinates of the boundary points columns 1:3, RF values
%   for the misorientations columns 4:6. 
%   con - the connectivity between the Bpts. ids in coulmns 1 and 2,
%   symmetries in columns 3 and 4. 
%   
%HISTORY: 
%Written by Nick Orr April 2020 in response to finding bugs in the previous
%versions of this code. 

%create an array of grain ids
G_idp = zeros(numel(xyzqs(:,1)), 1); 
for a = 1:numel(G)
    G_idp(G{a}) = a; 
end 
G_id = G_idp(G_idp ~=  0);    %Removing the grain unassigned rows.
xyzqsG = xyzqs(G_idp ~= 0, :);%Removing the grain unassigned rows.

%Triangulate and find all edges. 
DT = delaunayTriangulation(xyzqsG(:,1), xyzqsG(:,2), xyzqsG(:,3)); 
DTedges = edges(DT); 
DTedgesG_id = G_id(DTedges);
IGedges = DTedges(DTedgesG_id(:,1) ~= DTedgesG_id(:,2),:); %intergrain edges
IGedges_id = DTedgesG_id(DTedgesG_id(:,1) ~= DTedgesG_id(:,2),:);
%crash 
%remove the particles at the edge of  the FOV. 
EDGE_R = xyzqsG(IGedges(:,1), 1) > minmaxx(:,1) & ...
    xyzqsG(IGedges(:,1), 1) < minmaxx(:,2) & ...
    xyzqsG(IGedges(:,2), 1) > minmaxx(:,1) & ...
    xyzqsG(IGedges(:,2), 1) < minmaxx(:,2) & ...
    xyzqsG(IGedges(:,1), 2) > minmaxy(:,1) & ...
    xyzqsG(IGedges(:,1), 2) < minmaxy(:,2) & ...
    xyzqsG(IGedges(:,2), 2) > minmaxy(:,1) & ...
    xyzqsG(IGedges(:,2), 2) < minmaxy(:,2) & ...
    xyzqsG(IGedges(:,1), 3) > minmaxz(:,1) & ...
    xyzqsG(IGedges(:,1), 3) < minmaxz(:,2) & ...
    xyzqsG(IGedges(:,2), 3) > minmaxz(:,1) & ...
    xyzqsG(IGedges(:,2), 3) < minmaxz(:,2);
IGedges_B = IGedges(EDGE_R,:);
IGedges_id_B = IGedges_id(EDGE_R,:);
%finding the boundary points and misorientations
BptsRF = zeros(numel(IGedges_B(:,1)), 10); 
for a = 1:numel(IGedges_B(:,1))
    E1 = IGedges_B(a, 1); 
    E2 = IGedges_B(a, 2); 
    p1 = xyzqsG(E1, 1:3); 
    p2 = xyzqsG(E2, 1:3); 
    q1 = xyzqsG(E1, 4:7);
    q2 = xyzqsG(E2, 4:7);
    s1 = xyzqsG(E1, 8); 
    s2 = xyzqsG(E2, 8); 
    g1g2 = IGedges_id_B(a,:); 
    BptsRF(a,1:3) = (p1 + p2)/2;
    if s1 == s2
        q2c = [q2(1), -q2(2:4)];
        [dis] = MisorientationQF(q1,q2c,...
            qOperatorCell{s1},...
            qOperatorCell{s2}, s1, s2); 
    elseif s1 == 2 && s2 == 1 %must follow low to high symmetry convention
        q2c = [q2(1), -q2(2:4)];
        [dis] = MisorientationQF(q1,q2c,...
            qOperatorCell{s1},...
            qOperatorCell{s2}, s1, s2); 
    elseif s1 == 1 && s2 == 2
        q1c = [q1(1), -q1(2:4)];
        [dis] = MisorientationQF(q2,q1c,...
            qOperatorCell{s2},...
            qOperatorCell{s1}, s2, s1); 
    end
    BptsRF(a, 4:6) = dis;
    BptsRF(a, 7:8) = [s1, s2]; 
    BptsRF(a, 9:10) = g1g2; %the grains to which they belong
end
%crash
%finding the tetrahedra that each edge belong to 
Td = DT.ConnectivityList; %connectivity list for the original DT
Tdid = repmat([1:numel(Td(:,1))]',1, 4); 
Tdidlong = reshape(Tdid', numel(Td), 1); 
Tdlong = reshape(Td', numel(Td), 1); 
Tdshare = cell(numel(IGedges_B(:,1)),1);
TdshareN = zeros(numel(IGedges_B(:,1)), 1); 
for a = 1:numel(IGedges_B(:,1))
    idsa = IGedges_B(a,:); 
    idlong1 = Tdidlong(Tdlong == idsa(1));
    idlong2 = Tdidlong(Tdlong == idsa(2));
    Tdshare{a} = intersect(idlong1, idlong2); %the tetrahedra with that edge
    %i.e. the Td that the boundary point belogs to. 
    %I want the boundary points for each tetrahedra. 
    TdshareN(a) = numel(Tdshare{a}); %the number of tetrahedra that share that edge 
end 
TotalTdshareN = sum(TdshareN); 
Tdsharelong = zeros(TotalTdshareN, 1); 
Tdsharelongid = zeros(TotalTdshareN, 1);
b = 1; 
for a = 1:numel(Tdshare)
    Tdsharea = Tdshare{a}; 
    Tdsharelong(b:b + TdshareN(a)-1) = Tdsharea;
    Tdsharelongid(b:b + TdshareN(a)-1) = a; %the Bpt id
    b = b + TdshareN(a);
end 
TdU = unique(Tdsharelong); %the unique td with bpts in them 
TdUN = zeros(numel(TdU), 1);
TdUcell = cell(numel(TdU),1); 
numgrains = zeros(numel(TdU), 1);
%crash
for a = 1:numel(TdU)
   TdUa = TdU(a); 
   TdUcella = Tdsharelongid(Tdsharelong == TdUa); 
   TdUcell{a} = TdUcella;
   TdUN(a) = numel(TdUcell{a});
   grains = BptsRF(TdUcella,9:10);  
   Ugrains = unique(grains); 
   numgrains(a) = numel(Ugrains); 
end
TdUcell2 = TdUcell(numgrains == 2 & TdUN > 1); 
TdUN2 = TdUN(numgrains == 2 & TdUN > 1);
%TdUNtotal = sum(TdUN); 
totalconnect = numel(find(TdUN2 == 2)) + numel(find(TdUN2 == 3))*3 +...
    numel(find(TdUN2 == 4)) * 5;
%crash
con = zeros(totalconnect,4);
TetEdgeID = [1,2;1,3;1,4;2,3;2,4;3,4];  %used for drawing meshes with
%4 Bpts in a Td w/ two flavours of verticies. 
b = 1; %the counter 
for a = 1:numel(TdUcell2)
    TdUcella = TdUcell2{a}; %the points in a Td
    %if TdUN2(a) == 2      
    %    conBpts = [TdUcella(1), TdUcella(2)]; %invert just so the columns are right 
    %    symmetries = BptsRF(TdUcella(1),7:8);
    %    con(b,:) = [conBpts, symmetries];
    %    b = b +1;
    if TdUN2(a) == 3
        %grains = BptsRF(TdUcella(1:3),9:10);  
        %Due to edges it is also possible that there are three grains
        %connected in this tetreahedron. 
        %Ugrains = unique(grains); 
        %if numel(Ugrains) == 2
            symmetries = BptsRF(TdUcella(1:3),7:8);
            conBpts = [TdUcella(1), TdUcella(2); TdUcella(2),...
                TdUcella(3); TdUcella(3), TdUcella(1)];        
            con(b:b+2,:) = [conBpts, symmetries];
            b = b + 3;
        %end
        %crash
    elseif TdUN2(a) == 4
        points = BptsRF(TdUcella,1:3); 
        dists = sum((points(TetEdgeID(:,1)) -...
            points(TetEdgeID(:,2))).^2, 2);
        %must select the shortest five; 
        [~, distsIDO] = sort(dists);
        shortedges = TetEdgeID(distsIDO(1:5), :);
        conBpts = [TdUcella(shortedges(1,1)), TdUcella(shortedges(1,2));...
            TdUcella(shortedges(2,1)), TdUcella(shortedges(2,2));...
            TdUcella(shortedges(3,1)), TdUcella(shortedges(3,2));...
            TdUcella(shortedges(4,1)), TdUcella(shortedges(4,2));...
            TdUcella(shortedges(5,1)), TdUcella(shortedges(5,2))];
        symmetries = BptsRF(TdUcella,7:8);
        symmetries = [symmetries;symmetries(1,:)]; %they should all have...
        %the same symmetry  
        con(b:b+4,:) = [conBpts, symmetries];
        b = b + 5;
    end
    %conBpts
    
    %xyzqs(IGedges_B(? ,8)
end
%crash
con = con(con(:,1) ~=0,:); %just while im debugging
BptsRF = BptsRF(:,1:8); 
%crash 
end

