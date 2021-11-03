function [  ] = MeshPlot( Bpts, con, col )
%NAME: MeshPlot
%FUNCTION: Create a plot of the grain boundary mesh. 
% NOTE! You must create a figure with hold on for this plotting tool to
% work. 
%INPUTS : 
%   Bpts - The coordinates to be connected with a mesh n by 3 array. 
%   con -  n by 2 array with the ids of Bpts to be connected with a line
%   (Mesh edge)
%   col - the colour of the mesh edge. 
%OUTPUTS:
%   Creates a plot in a figure
%HISTORY: 
%Written by Nick Orr April 2019. 
for a = 1:numel(con(:,1))
        plot3([Bpts(con(a,1), 1), Bpts(con(a,2), 1)], [Bpts(con(a,1), 2), Bpts(con(a,2), 2)], [Bpts(con(a,1), 3), Bpts(con(a,2), 3)], 'Color',col(a,:)); 
end 

end

