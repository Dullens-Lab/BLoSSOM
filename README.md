# BLoSSOM
Boundary, Location, Structure, Surface, Orientation and Misorientation. 

Please see our paper at https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.123605 for more details.

# Boundary, Location, Structure, Surface, Orientation and Misorientation.

## ---------------Tutorial for BLoSSOM 1.0-----------------------------------
Written by: Nicholas Orr, September 2021, Oxford University, UK. 
email: nicholasorrchem@gmail.com
 
## ---------------GETTING STARTED--------------------------------------------
This code accompanies the paper: "Grain boundary characterisation from particle coordinates" by N. H. P. Orr, T. Yanagishima, E. Maire and R. P. A. Dullens. 
 
This package of functions was developed in MatLabR1029b and the robotics toolbox. The functions are presented in a  modular form so that the results from each step of the method are accessible. The following demonstrates the basic function of the grain boundary identification code on an example hard sphere colloidal crystal data set.

Load the example data. You will find the variable `dat` in your workspace.

`load('example_data.mat')`
 
Plot the example data. The function p3d is a fast way of making 3d plots.

`p3d(dat) ;` 
 
The coordinates show the positions of colloidal particles in a crystal at the bottom and with liquid like order at the top. The coordinates have units of micrometres, um.
 
## -------------NEAREST NEIGHBOURS-------------------------------------------
 
Find the radial distribution function, `g(r)`.

```Matlab
[xx, g] = gofr3D_NP(dat, 500, 250,50, 0.1, 5);

plot(g, xx); xlabel('r'); ylabel('g(r)')
```
 
We should see a figure containing the g(r) with a peak at r = 1.8 um, and a first minimum at 2.3 um. The first peak tells us the mean nearest neighbour distance and the first minimum tells us the cut off distance for the nearest neighbour search.
 
Next, we identify the nearest neighbour clusters, NNCs.

```Matlab
maxN = 12; 
%the maximum number of neighbours. For close packed systems (which hard sphere colloidal crystals often are) maxN = 12 due to geometric constraints.

Rcut = 2.3; 
%identify the cut off for neighbour search. NN = NNsearch(dat, maxN, Rcut); %NN contains the neighbours. 
```

The function `NNsearch` finds neighbours as:

"The up to maxN points within a cut off distance"

Alternatively, the function `NNsearchDT` may be used that calculates the neighbours as

"The up to maxN points connected by a Delaunay triangulation edge within a cut off distance".

For this tutorial, we will simply use `NNsearch`.

We may look at how many neighbours each particle has. To do so meaningfully, we must remove particles that are on the edge of the field of view.

```Matlab
depth = Rcut; 
[NN_B, bulkid] = RemoveEdgesNN(dat, NN, depth); 
```

Plot a histogram of the number of nearest neighbours.

```Matlab
figure
histogram(NN_B.LISTN, [0:12])
set(gca,'FontSize',18)
set(gca, 'LineWidth', 1)
set(gcf,'renderer','Painters')
box on
ylabel('N_p')
xlabel('N_n')
``` 
 
We see a large peak at 12 as we expected. In fact, hard sphere colloidal crystals are most often hexagonal close packed, HCP, or face centered cubic, FCC. Therefore, we are going to test against these two structures. There is already functionality for simple cubic, SC, and body centered cubic, BCC. Searching for other structures is possible. You will need to find a reference NNC and the symmetry operators as well as finding the symmetry in-equivalent permutations of points for the point set registration step. The basics of how to do this is in the supplementary information of the paper, and implemented in the function `symmetry_unique_pselect.m`. See the documentation within.

It is important that the number of neighbours is equal to or below the number of points in the reference NNC, otherwise, the code with fail to register the points. The way you find NNCs is crucial to determining structure using this method and forms part of the structure assginment process. Take care!
 
## -----------------POINT SET REGISTRATION-----------------------------------
 
Fistly, we load in the reference NNCs:

```Matlab
load('Reference_Structures.mat', 'fccxyz')
load('Reference_Structures.mat', 'hcpxyz')
```

Secondly, we load in the symmetry in-equivalent permutations for efficient point set registration.

```Matlab
load('idx_symm.mat', 'idx_symm_FCC')
load('idx_symm.mat', 'idx_symm_HCP')
symmIDcell{1} = idx_symm_FCC; %setting up the symmetry selection cell
symmIDcell{2} = idx_symm_HCP;
```

`symmIDcell` contains the symmetry in-equivalent permutations of points. Scaling the reference NNCs to be as large as the specimen NNCs.

```Matlab
d = 1.8; %Nearest neighbour distance, see radial distribution function
refcell{1} = fccxyz*d; 
refcell{2} = hcpxyz*d;
```

IMPORTANT
We have set the reference NNCs so that reference 1 is FCC and reference 2 is HCP. Note that in refcell and symmIDcell the FCC reference NNC is in position 1 and the HCP reference NNC is in position 2.
 
`struct_orient` will find the structure and orientation of each NNC.

```Matlab
tic
[ pos, orientation_array, fitTsqr, structcode] =...
    struct_orient( dat,...
    NN, refcell, 3, symmIDcell, Rcut);
toc
%The output, FitTsqr is the mean square of the differences between the
%registered points at optimum overlap. The following converts to root mean
%square. 
fitT = zeros(numel(fitTsqr),1);
fitT(fitTsqr >=0) = sqrt(fitTsqr(fitTsqr >= 0))/d;
fitT(fitTsqr < 0) = -1; 
```

We can check how well the NNCs match the reference NNCs by plotting histograms of the RMS deviations

```Matlab
figure
hist(fitT, 100)
xlabel('RMS displacement of registered points')
ylabel('number')
set(gca, 'fontsize', 18)
%The peak at -1 shows all of the points for which registration has failed.
%We do not have a structure or orientation measurements for these points. 
%We can also look at the assignments for FCC and HCP individually: 
figure;
subplot(2, 1, 1)
histogram(fitT(structcode == 1 & fitT >= 0),0:0.01:1, 'facecolor','green')
%note structcode == 1 is for FCC
xlabel('RMS displacement of registered points')
ylabel('number')
set(gca, 'fontsize', 18)
subplot(2, 1, 2)
histogram(fitT(structcode == 2 & fitT >= 0),0:0.01:1, 'facecolor','red')
%note structcode == 2 is for HCP
xlabel('RMS displacement of registered points')
ylabel('number')
set(gca, 'fontsize', 18)
```
 
%-------------STRUCTURE AGGSINMENT-----------------------------------------
 
%We can find FCC and HCP particles using the following. 
%note that structcode == 1 and 2 are for FCC and HCP respectively and that
%fitT >= 0 selects for particles with successful NNC registrations.
posFCC = pos(structcode == 1 & fitT >= 0, :); 
posHCP = pos(structcode == 2 & fitT >= 0, :); 
fitTFCC = fitT(structcode == 1 & fitT >= 0, :);
fitTHCP = fitT(structcode == 2 & fitT >= 0, :);
 
%separating the orientation arrays by structure. 
structcode_rep = repmat(structcode, 1, 3);
structcode_long = reshape(structcode_rep', 3* numel(structcode), 1);  
% making the symmcode repeat three times for each entry so the rotation 
%matrix can be recalled.
OFCC = orientation_array(structcode_long == 1, :); %the orientation array 
%for FCC 
OHCP = orientation_array(structcode_long == 2, :);  %the orientation array 
%for HCP
 
%The next code converts the FCC orientation arrays into quaternions. 
b = 1;
qARRFCC = zeros(numel(OFCC(:,1))/3,4); 
for a = 1:3:numel(OFCC(:,1))
    qARRFCC(b,:) = rotm2quat(OFCC(a:a+2, :));
    b = b + 1; 
end
fitfcc = fitT(structcode == 1); %some indexing 
qARRFCC_s = qARRFCC(fitfcc >= 0, :);
 
%orientational fundamental zone
load('QuaternionSymmetyOperator.mat') % Loading symmetry elements. 
%qO, qD3 and qD6 contain the rotational symmetry elements of the O D3 and
%D6 rotational symmetry groups respectively. The symmetry elements
%correspond to the axes of symmetry in the reference NNCs
load('Fundamental_Zones.mat'); %This loads the bounding boxes for 
%fundamental zones
%plotting the FCC NNC orientations within the fundamental zone of 
%orientation. 
[qfccFZ, RFfccFZ] = FundamentalZoneQ(qARRFCC_s, qO); 
figure; hold on 
p3d(RFfccFZ); 
p3d(O_ofz, '-')
title 'Correlated and uncorrelated FCC orientations FZ'
[qhcpFZ, RFhcpFZ] = FundamentalZoneQ(qARRHCP_s, qD6); 
figure; hold on 
p3d(RFhcpFZ); 
p3d(D6_ofz, '-')
title 'Correlated and uncorrelated HCP orientations FZ'
 
%--------------GRAIN SEPARATION--------------------------------------------
 
%Here we perform position and orientation correlations to separate the
%particles into grains.
 
thetacut = 4; % the misorientation cut off in degrees. 
 
%Start with FCC
[G_fcc, UC, MO] = ...
    PositionOrientationCorrelation(posFCC,...
    qARRFCC_s, qO, qO, Rcut, thetacut);  %This function finds grains. 
figure; hold on; 
[ col_GFCC, xyzabg_GFCC ] = Grain_Plot(G_fcc, [posFCC, qARRFCC_s, ...
    fitTFCC], 1, 1, 5); % This function plots the grains.
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
 
%In this example, xyzabg_GFCC contains columns 1-3 for the x y z Cartesian 
%coordinate, column 4 contains the row ID from the original data set,
%columns 5 - 8 contain the quaternion orientation of the particle's NNC,
%and column 9 for the RMS deviation at optimum reference/ specimen NNC
%overlap (see the 2nd input for Grain_Plot). Each row is a new particle in
%an FCC grain. 
 
%Each grain is contained in a cell of G_fcc
%we can plot the 1st grain: 
figure;
Grain_Plot(G_fcc(1), [posFCC, qARRFCC_s, fitTFCC],...
    1, 1, 5);
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
%we can plot the 2nd grain: 
figure;
Grain_Plot(G_fcc(2), [posFCC, qARRFCC_s, fitTFCC],...
    1, 1, 5);
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
 
%now we convert orientation arrays to quaternions and perform orientation 
%correlations for HCP 
b = 1;
qARRHCP = zeros(numel(OHCP(:,1))/3,4); 
for a = 1:3:numel(OHCP(:,1))
    qARRHCP(b,:) = rotm2quat(OHCP(a:a+2, :));
    b = b + 1; 
end
fithcp = fitT(structcode == 2); %some indexing 
qARRHCP_s = qARRHCP(fithcp >= 0, :);
[G_hcp, UChcp, MOhcp] =...
    PositionOrientationCorrelation(posHCP,...
    qARRHCP_s, qD6, qD6, Rcut, thetacut); %Note that we use D6 instead of 
%D3 for HCP because the crystal has D6 symmetry but each NNC only has D3 
%symmetry.
figure; hold on; 
[ col_GHCP, xyzabg_GHCP ] = Grain_Plot(G_hcp, [posHCP, qARRHCP_s, ...
    fitTHCP], 1, 1, 5); 
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
%xyzabg_GHCP contains the same columns as for xyzabg_GFCC. Each row is 
%instead a new particle in an FCC grain. 
 
%Now we have separated the particles into grains and performed orientation
%correlations, we can re plot the orientation fundamental zones, with
%spuriously identified orientations removed. 
[qfccFZ_C, RFfccFZ_C] = FundamentalZoneQ(xyzabg_GFCC(:, 5:8), qO); 
figure; hold on 
p3d(RFfccFZ_C); 
p3d(O_ofz, '-')
title 'Correlated FCC orientations FZ'
xlabel('RF_x'); ylabel('RF_y'); zlabel('RF_z')
%We should now see that only the densely clustered regions of orientation
%survive. Try playing around with thetacut and see how this effects the
%plots. 
%We do the same for HCP now. 
[qhcpFZ_C, RFhcpFZ_C] = FundamentalZoneQ(xyzabg_GHCP(:, 5:8), qD6); 
figure; hold on 
p3d(RFhcpFZ_C); 
p3d(D6_ofz, '-')
title 'Correlated FCC orientations FZ'
xlabel('RF_x'); ylabel('RF_y'); zlabel('RF_z')
 
%--------------AMORPHOUS PARTICLES-----------------------------------------
 
%Any particles that are not assigned to a grain may be considered amorphous. 
[GC, xyzC] = Grain_Combine(G_fcc', G_hcp', [posFCC(:,1:3), qARRFCC_s],...
    [posHCP(:,1:3), qARRHCP_s], 1,2);
 
crystpos = zeros(0,3); 
for a = 1:numel(GC)
    crystpos = [crystpos ;xyzC(GC{a}, 1:3)]; 
end 
Ampos = setdiff(dat,crystpos(:,1:3), 'rows'); 
 
figure; hold on 
p3d(Ampos, 'ok', 'MarkerFaceColor', 'b')
p3d(crystpos, 'ok', 'MarkerFaceColor', 'y'); 
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum'); axis equal
%We see a crystal at the bottom of the field of view with an amorphous
%fluid region above. Particles that are too close to the edge are assigned
%as amorphous. 
 
%--------------GRAIN BOUNDARIES--------------------------------------------
 
%We can find the grain boundaries between all grains in the following way: 
 
%We can look for the grain boundaries between larger grain to simplify the
%picture. 
min_grainsize = 9; %Only consider grains with a size greater than 9.
sizegrain = zeros(1,numel(GC));
for a = 1:numel(GC)
sizegrain(a) = numel(GC{a});
end
Large = GC(sizegrain > min_grainsize);
figure;
Grain_Plot(Large, xyzC, 1, 3, 5); % we can see which grains we are 
%considering
title('Large FCC and HCP grains')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
 
 
L = 1.5; %the number of lattice constants to cut off around the edges. 
%Removes edge effects
LatticeConstant = d;
minmaxx = [(min(xyzC(:,1)) + LatticeConstant*L), (max(xyzC(:,1)) - ...
LatticeConstant*L)];
minmaxy = [(min(xyzC(:,2)) + LatticeConstant*L), (max(xyzC(:,2)) - ...
LatticeConstant*L)];
minmaxz = [(min(xyzC(:,3)) + LatticeConstant*L), (max(xyzC(:,3)) - ...
LatticeConstant*L)]; %edge effects defined by extent of grains not FOV
%GB_NQ
qOperatorCell{1} = qO; 
qOperatorCell{2} = qD6;
[Btot, BtotCon] = GB_Char(Large, xyzC, qOperatorCell, minmaxx, minmaxy,...
    minmaxz);
[BptsRF11, B11, Col11] = GB_symmetry_mesh(Btot, BtotCon, 1, 1); %FCC:FCC
[BptsRF21, B21, Col21] = GB_symmetry_mesh(Btot, BtotCon, 2, 1); %HCP:FCC
[BptsRF22, B22, Col22] = GB_symmetry_mesh(Btot, BtotCon, 2, 2); %HCP:HCP
[D11] = GB_SS(Btot, 1, 1 ); %selecting grain boundaries between FCC
[D21] = GB_SS( Btot, 1, 2 ); %selecting grain boundaries between HCP and FCC
[D22] = GB_SS(Btot, 2, 2 ); %selecting grain boundaries between HCP
 
%Grain boundaries between FCC grains
figure; hold on 
Angle_Colour(D11, tand(62.8/2), 200, 2, '.', 5);
title('FCC:FCC Disorientation Fundamental Zone') %Rodrigues Frank FZ 
xlabel('RF_x'); ylabel('RF_y'); zlabel('RF_z')
p3d(OO_dfz, '-'); axis equal
figure; hold on 
Angle_Colour(D11, tand(62.8/2), 200, 1, '.', 5); axis equal
title('FCC:FCC Grain Boundary Points')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
figure; hold on 
MeshPlot(BptsRF11, B11, Col11); axis equal 
title('FCC:FCC Grain boundary surfaces')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
 
%Grain boundary between FCC and HCP grains
figure; hold on 
Angle_Colour(D21, tand(56.60/2), 200, 2, '.', 5);
p3d(D6O_dfz, '-');  axis equal
title('HCP:FCC Disorientation Fundamental Zone') %Rodrigues Frank FZ
xlabel('RF_x'); ylabel('RF_y'); zlabel('RF_z')
figure; hold on 
Angle_Colour(D21, tand(56.60/2), 200, 1, '.', 5); axis equal
title('HCP:FCC Grain Boundary Points')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
figure; hold on 
MeshPlot(BptsRF21, B21, Col21); axis equal
title('HCP:FCC Grain boundary surfaces')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
 
%Grain boundary between HCP grains
figure; hold on 
Angle_Colour(D22, tand(93.8/2), 200, 2, '.', 5);
p3d(D6D6_dfz, '-');  axis equal
title('HCP:HCP Disorientation Fundamental Zone') %Rodrigues Frank FZ 
xlabel('RF_x'); ylabel('RF_y'); zlabel('RF_z')
figure; hold on 
Angle_Colour(D22, tand(93.8/2), 200, 1, '.', 5); axis equal
title('HCP:HCP Grain Boundary Points')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
figure; hold on 
MeshPlot(BptsRF22, B22, Col22); axis equal
title('HCP:HCP Grain boundary surfaces')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
 
%Above we have found the boundaries considering FCC and HCP grains.
%Stacking faults, where HCP layers can separate FCC grains may prevent twin
%boundaries from being identified. Therefore we can find te grain
%boundaries between FCC grains, ignoring HCP grains: 
 
%Look for the large FCC grains. 
sizegrain = zeros(1,numel(G_fcc));
for a = 1:numel(G_fcc)
sizegrain(a) = numel(G_fcc{a});
end
Large_fcc = G_fcc(sizegrain > min_grainsize);
figure;
[ col_GFCC_L, xyzabg_GFCC_L ] = Grain_Plot(Large_fcc, posFCC, 1, 3, 5); 
title('Large FCC grains')
% we can see which grains we are considering
 
LFCC = [posFCC(:,1:3), qARRFCC_s, ones(numel(posFCC(:,1)),1)]; %This line 
%of code does what Grain_Combine does, by appending the strcuture id to the
%end. However, here we only want FCC grains. 
LatticeConstant = d;
minmaxxLFCC = [(min(LFCC(:,1)) + LatticeConstant*1.5), (max(LFCC(:,1)) - ...
LatticeConstant*1.5)];
minmaxyLFCC = [(min(LFCC(:,2)) + LatticeConstant*1.5), (max(LFCC(:,2)) - ...
LatticeConstant*1.5)];
minmaxzLFCC = [(min(LFCC(:,3)) + LatticeConstant*1.5), (max(LFCC(:,3)) - ...
LatticeConstant*1.5)]; %edge effects defined by extent of grains not FOV
%Above we have dealt with edge effects (if you dont do this grain 
%boundaries) can wrap around the edges of grains near the edge of the FOV. 
 
qOperatorCell{1} = qO; %These don't need to be redefined again. I just 
%include it here for convenience
qOperatorCell{2} = qD6;
[BtotLFCC, BtotConLFCC] = GB_Char(Large_fcc, LFCC, qOperatorCell, ...
    minmaxxLFCC, minmaxyLFCC, minmaxzLFCC);
[BptsRF11LFCC, B11LFCC, Col11LFCC] = GB_symmetry_mesh(BtotLFCC, ...
    BtotConLFCC, 1, 1);
[D11_LFCC] = GB_SS( BtotLFCC, 1, 1 );
 
%Grain boundary between FCC grains, ignoring the HCP grains. 
figure; hold on 
Angle_Colour(D11_LFCC, tand(62.8/2), 200, 2, '.', 5); 
p3d(OO_dfz, '-k'); axis equal
title('FCC:FCC no HCP Disorientation Fundamental Zone') %Rodrigues Frank FZ 
xlabel('RF_x'); ylabel('RF_y'); zlabel('RF_z')
figure; hold on 
Angle_Colour(D11_LFCC, tand(62.8/2), 200, 1, '.', 5); axis equal
title('FCC:FCC no HCP Grain Boundary Points')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
figure; 
hold on; 
MeshPlot(BptsRF11LFCC, B11LFCC, Col11LFCC); axis equal
title('FCC:FCC no HCP Grain boundary surfaces')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
%Notice now how the grain boundary surface joins into one now. 
 
%we can plot the grains with the boundary surfaces too. 
figure
hold on; 
MeshPlot(BptsRF11LFCC, B11LFCC, Col11LFCC); axis equal
Grain_Plot(Large_fcc, posFCC, 1, 1, 5); 
title('FCC:FCC no HCP Grain boundary surfaces and grains')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
 
%We now do the same for HCP, we ignore the FCC grains
 
%Look for the large HCP grains. 
sizegrain = zeros(1,numel(G_hcp));
for a = 1:numel(G_hcp)
sizegrain(a) = numel(G_hcp{a});
end
Large_hcp = G_hcp(sizegrain > min_grainsize);
figure;
[ col_GHCP_L, xyzabg_GHCP_L ] = Grain_Plot(Large_hcp, posHCP, 1, 3, 5); 
% we can see which grains we are considering
L = 1.5; 
LHCP = [posHCP(:,1:3), qARRHCP_s, ones(numel(posHCP(:,1)),1) + 1];
LatticeConstant = d;
minmaxxLHCP = [(min(LHCP(:,1)) + LatticeConstant*L), (max(LHCP(:,1)) - ...
LatticeConstant*L)];
minmaxyLHCP = [(min(LHCP(:,2)) + LatticeConstant*L), (max(LHCP(:,2)) - ...
LatticeConstant*L)];
minmaxzLHCP = [(min(LHCP(:,3)) + LatticeConstant*L), (max(LHCP(:,3)) - ...
LatticeConstant*L)]; %edge effects defined by extent of grains not FOV
%GB_NQ
qOperatorCell{1} = qO; 
qOperatorCell{2} = qD6;
[BtotLHCP, BtotConLHCP] = GB_Char(Large_hcp, LHCP, qOperatorCell, ...
    minmaxxLHCP, minmaxyLHCP, minmaxzLHCP);
[BptsRF22LHCP, B22LHCP, Col22LHCP] = GB_symmetry_mesh(BtotLHCP, ...
    BtotConLHCP, 2, 2);
[D22_LHCP] = GB_SS( BtotLHCP, 2, 2 );
figure; hold on 
Angle_Colour(D22_LHCP, tand(93.8/2), 200, 2, '.', 5);
p3d(D6D6_dfz, '-k'); axis equal
title('HCP:HCP no FCC Disorientation Fundamental Zone') %Rodrigues Frank FZ 
xlabel('RF_x'); ylabel('RF_y'); zlabel('RF_z')
figure; hold on 
Angle_Colour(D22_LHCP, tand(93.8/2), 200, 1, '.', 5); axis equal
title('HCP:HCP no FCC Grain Boundary Points')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
figure; 
hold on; 
MeshPlot(BptsRF22LHCP, B22LHCP, Col22LHCP); axis equal 
title('HCP:HCP no FCC Grain Boundary Surfaces ')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
 
%Plot the boundaries with the grains. 
figure; 
hold on; 
MeshPlot(BptsRF22LHCP, B22LHCP, Col22LHCP); 
Grain_Plot(Large_hcp, posHCP, 1, 1, 5); axis equal
title('HCP:HCP no FCC Grain Boundary Surfaces and Grains')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
 
%-----------GRAIN BOUNDARY PARTICLES---------------------------------------
 
%Find the GB particles for the GB between the large FCC grains. 
GBptcl = GrainBoundaryParticles(Ampos, D11_LFCC(:,1:3), Rcut);
figure; 
p3d(GBptcl, 'ok', 'MarkerFaceColor', 'g'); 
title('Large FCC:FCC GB particles')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
 
%Plot them with the two large FCC grains. 
figure; hold on 
p3d(GBptcl, 'ok', 'MarkerFaceColor', 'g'); 
Grain_Plot(Large_fcc, posFCC, 1, 1, 5); 
title('Large FCC:FCC GB particles with Grains')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
 
%We can find more of the GB particles by including more GBpoints as the
%input for GrainBoundaryParticles: 
GBtot =  [D11_LFCC(:,1:3); D21(:,1:3); D22_LHCP(:,1:3)];
GBptcltot = GrainBoundaryParticles(Ampos, GBtot, Rcut);
figure; 
p3d(GBptcltot, 'ok', 'MarkerFaceColor', 'g'); 
title('GB particles')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
 
figure; 
hold on; 
MeshPlot(BptsRF22LHCP, B22LHCP, Col22LHCP); 
MeshPlot(BptsRF21, B21, Col21); 
MeshPlot(BptsRF11LFCC, B11LFCC, Col11LFCC); 
p3d(GBptcltot, 'ok', 'MarkerFaceColor', 'g'); 
title('GB particles')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
 
GBptcl21 = GrainBoundaryParticles(Ampos, D21(:,1:3), Rcut);
figure; hold on; 
p3d(GBptcl21, 'ok', 'MarkerFaceColor', 'g'); 
title('GB particles HCP:FCC')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
MeshPlot(BptsRF21, B21, Col21); 
%Notice how there are not grain boundary points on the stacking faults! 
 
GBptcl22 = GrainBoundaryParticles(Ampos, D22(:,1:3), Rcut);
figure; hold on; 
p3d(GBptcl22, 'ok', 'MarkerFaceColor', 'g'); 
title('GB particles HCP:FCC')
xlabel('x /\mum'); ylabel('y /\mum'); zlabel('z /\mum')
MeshPlot(BptsRF22, B22, Col22); 
p3d(D22_LHCP); %The grain boundary points
%NB that the grain boundary points extend further to the edges than the
%mesh surface does. We may change the cut off for the edge effects for the
%grain boundary surface calculation if desired. However, removing too
%little edge can cause GB surfaces to wrap around grains in an unphysical
%way. 
 
 
%-----------MISORIENTATION ANGLE HISTOGRAMS--------------------------------
 
%It can be quite useful to see what the distribution of misorientation
%angles looks like. 
%These are probablity distributions with the area under the curves equaling
%1.
load('Random_Misorientation_Distributions.mat'); 
 
n = 70; %number of bins. 
fontsz = 18;
binnumbersHC = 0 + (56.6/(2*n)):56.6/n:56.6 - (56.6/(2*n));
Disorientation_Angle_Histogram(D21(:,4:6), 3, 1, binnumbersHC, ...
    D6O_rd);
hold on
xlim([0, 56.6])
set(gca,'FontSize',fontsz)
box on 
xtickformat('%.1f')
xticks([0, 30, 56.6])
ylabel('P(\theta)')
xlabel('\theta /^o')
set(gcf,'renderer','Painters')
title('HCP:FCC')
 
 
n = 70;  %number of bins. 
binnumbersCC = 0 + (62.8/(2*n)):62.8/n:62.8 - (62.8/(2*n));
Disorientation_Angle_Histogram(D11_LFCC(:,4:6), 3, 1, binnumbersCC, ...
    OO_rd);
xlim([0, 62.8])
set(gca,'FontSize',fontsz)
box on 
xtickformat('%.1f')
xticks([0, 62.8])
ylabel('P(\theta)')
xlabel('\theta /^o')
set(gcf,'renderer','Painters')
title('FCC:FCC no HCP')
 
n = 70;
MAXA = 93.84;
binnumbersHH = 0 + (MAXA/(2*n)):MAXA/n:MAXA - (MAXA/(2*n));
Angrr33_50mil = Disorientation_Angle_Histogram(D22_LHCP(:,4:6), 3, 3,...
    binnumbersHH, D6D6_rd); 
xlim([0, MAXA])
set(gca,'FontSize',fontsz)
box on
xtickformat('%.1f')
ylabel('P(\theta)')
xlabel('\theta /^o')
set(gcf,'renderer','Painters')
title('HCP:HCP no FCC')
 
%-----------USING PyMol TO CREATE PLOTS------------------------------------
 
%This if for plotting things with the python program PyMol. These functions 
%produce .py and .xyz files that render the grains and grain boundaries.
%Please see PyMol's own documentation. 
 
%plot all of the FCC and HCP grains
PymolPlot(col_GHCP, xyzabg_GHCP, 'HCP_G'); 
PymolPlot(col_GFCC, xyzabg_GFCC, 'FCC_G'); 
 
%plot the large FCC and HCP grains
PymolPlot(col_GHCP_L, xyzabg_GHCP_L, 'HCP_G_L'); 
PymolPlot(col_GFCC_L, xyzabg_GFCC_L, 'FCC_G_L'); 
 
%plot the grain boundary surfaces
GBmesh_Pymol(BptsRF21, B21, Col21, 'GBFCCHCP')
GBmesh_Pymol(BptsRF11LFCC, B11LFCC, Col11LFCC, 'FCCFCCnoHCP');
GBmesh_Pymol(BptsRF22LHCP, B22LHCP, Col22LHCP, 'HCPHCPnoFCC'); 
 
%-------------TROUBLE SHOOTING---------------------------------------------
%contact Nick by email at nicholasorrchem@gmail.com with the 
%subject line: BLoSSOM
%
%
%Nicholas Orr, September 2021. 
%--------------------------------------------------------------------------
 

