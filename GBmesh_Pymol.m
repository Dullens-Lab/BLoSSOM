function [  ] = GBmesh_Pymol( xyz, DEdgeConnect, colors, prefix )
%NAME: GBmesh_Pymol
%FUNCTION: form an xyztinker format file from particle cooridnates and 
%connectivities from delaunay triangulation. S
%INPUTS: 
%xyz - cartesian coordinates, new row new point. 
%DEdgeconnect - Delaunay edge indeces. 
%prefix - the file name prefix. 
%OUTPUTS:
%three files, prefixmesh.py, prefixcolor.py and prefix.xyz

filename = strcat(prefix,'mesh.py');
f = fopen(filename, 'w');
%s = sprintf('%d', numel(xyz(:,1)));
%fprintf(f, '%s\n', s);
for a = 1:numel(DEdgeConnect(:,1))
    s = sprintf('select aa, model %s and idx %d \nselect bb, model %s and idx %d \nbond aa, bb',prefix, DEdgeConnect(a,1),...
        prefix, DEdgeConnect(a,2)); 
    fprintf(f, '%s\n', s);
end
%fprintf(f, 'show wire \nhide sphere');
fclose(f); 

filename2 = strcat(prefix, '.xyz');
f = fopen(filename2, 'w');
s = sprintf('%d', numel(xyz(:,1)));
fprintf(f, '%s\n', s);
for a = 1:numel(xyz(:,1))
    s = sprintf('H %d %d %d', xyz(a,1), xyz(a,2), xyz(a,3)); 
    fprintf(f, '%s\n', s);  
end 
fclose(f); 

% filename3 = strcat(prefix, 'color.py'); 
% f = fopen(filename3, 'w'); 
% for a = 1:numel(DEdgeConnect(:,1))
%     s = sprintf('set_color bondcolor, [%d, %d, %d]', colors(a,1), colors(a,2), colors(a,3));
%     fprintf(f, '%s\n', s);
%     s = sprintf('set_bond stick_color, bondcolor, model %s and idx %d, model %s and idx %d',...
%         prefix, DEdgeConnect(a,1), prefix, DEdgeConnect(a,2));
%     fprintf(f, '%s\n', s);
% end 
% fclose(f); 

filename3 = strcat(prefix, 'color.py'); 
f = fopen(filename3, 'w'); 
for a = 1:numel(DEdgeConnect(:,1))
    s = sprintf('set_color bondcolor%d, [%d, %d, %d]',a, colors(a,1), colors(a,2), colors(a,3));
    fprintf(f, '%s\n', s);
end 
for a = 1:numel(DEdgeConnect(:,1))
    s = sprintf('set_bond stick_color, bondcolor%d, model %s and idx %d, model %s and idx %d',...
        a, prefix, DEdgeConnect(a,1), prefix, DEdgeConnect(a,2));
    fprintf(f, '%s\n', s);
end 
fclose(f); 

% filename4 = strcat(prefix, 'dis.xyz');
% f = fopen(filename4, 'w');
% s = sprintf('%d', numel(xyz(:,1)));
% fprintf(f, '%s\n', s);
% for a = 1:numel(xyz(:,1))
%     s = sprintf('H %d %d %d', xyz(a,4), xyz(a,5), xyz(a,6)); 
%     fprintf(f, '%s\n', s);  
% end 
% fclose(f); 
end

