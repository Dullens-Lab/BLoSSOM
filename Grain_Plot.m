function [ colours, xyzabgS_s ] = Grain_Plot( G, xyzabgS, m, seed, markersize )
%NAME: Grain_Plot
%FUNCTION: Plot grains in different colours
%INPUTS: 
%   G - Cell array containing the indeces of the grain members. 
%   xyzabgS - matrix array contains the coordinates and orientation
%   parameters and symmetry.
%   m - = 0 to not plot anything, only generates the colours for the
%       particles. 
%       = 1 for physical space, 
%       = 2 for orientation space,
%   seed - version of colourmap 3 and 4 to just plot symm 1 or 3 by
%   themselves 
    
%OUTPUTS: 
%   colours - matrix containing the rgb colours for each particle by row
%   in 0- 225.
%   xyzabgS_s - The particles assigned to grains (some particles are not 
%   assigned to grains). Contains the rows of the input xyzabgS that are 
%   assigned to a grain. 
%HISTORY: 
%   Written by Nick Orr 26/09/18
%   overwrote the seed so that you can put in some ids, line 21 to line 25
%   04/09/19
%   Nick Orr 09/04/19 added the axis equal command. 
%   Nick Orr 22/09/21 edited Grain_Plot_Legacy to remove some of the unused
%   plotting modes for different values of the input m. Updated the
%   documentation
Data = xyzabgS;
cmap = hsv(numel(G)); %hsv
if numel(seed) > 1
    cmap_s = cmap(seed, :); 
else
    cmap_s = scramble_rep(cmap, numel(cmap(:,1)), seed);
end
%cmap_s = [cmap_s(1,:); cmap_s(2,:) + [0, 0.25, 0.25]; cmap_s(3:9,:)]; %EDIT THE COLOUR!
colours = zeros(numel(xyzabgS(:,1)), 3) - 1; 
%figure
id_s = zeros(0,1);
if m == 2; 
    hold on
    for b = 1:numel(G)      
        plot3(Data(G{b}, 4),Data(G{b}, 5),...
        Data(G{b}, 6),'o','MarkerfaceColor',cmap_s(b,:), 'color', cmap_s(b, :), 'MarkerSize', markersize);
        id_s = [id_s; G{b}]; 
        colours(G{b}', 1) = cmap_s(b, 1) * 225; 
        colours(G{b}', 2) = cmap_s(b, 2) * 225; 
        colours(G{b}', 3) = cmap_s(b, 3) * 225; %for some reason doing this in one line doesnt work
    end
    axis equal
elseif m == 1; 
    hold on
    for b = 1:numel(G)
        plot3(Data(G{b}, 1),Data(G{b}, 2),...
        Data(G{b}, 3),'o','MarkerfaceColor',cmap_s(b,:), 'color', 'k', 'MarkerSize', markersize); %cmap_s(b, :)
        id_s = [id_s; G{b}];
        colours(G{b}', 1) = cmap_s(b, 1) * 255; 
        colours(G{b}', 2) = cmap_s(b, 2) * 255; 
        colours(G{b}', 3) = cmap_s(b, 3) * 255; %for some reason doing this in one line doesnt work
    end
    axis equal
    elseif m == 0 
    for b = 1:numel(G)
        id_s = [id_s; G{b}];
    end
end 
colours = colours(id_s, :);
xyzabgS_s = xyzabgS(id_s, :); 


end

