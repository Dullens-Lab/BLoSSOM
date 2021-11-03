function [  ] = xyzwrite( xyz, prefix )
%NAME: xyzwrite
%FUNCTION: writes an xyz file from the input coordinates. All the
%coordinates are given the same atom type.
%INPUTS: 
%   xyz - n by 3 array of coordinates where each row corresponds to a new 
%       entry.
%   prefix - a string used to name the file
%OUTPUTS: 
%   prefix.xyz file named as specified in the input.
%HISTORY:
%   Written by Nick Orr 16/ 01/ 19. It is a rainy day today!

filename_xyz = strcat(prefix, '.xyz'); 
f = fopen(filename_xyz, 'w'); 
s = sprintf('%d', numel(xyz(:,1)));
fprintf(f, '%s\n', s);
for a = 1:numel(xyz(:,1))
    s = sprintf('H %d %d %d', xyz(a,1), xyz(a,2), xyz(a,3)); 
    fprintf(f, '%s\n', s); 
end 
fclose(f);
end