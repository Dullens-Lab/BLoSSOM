function [ ] = PymolPlot( Colors, xyz, prefix )
%NAME: PymolPlot
%FUNCTION: 
%convert from matlab variables to formats and cmmands that pymol can
%interpret. 
%INPUTS: 
%   Colors, n by 3 array with the colors r g b from 1 to 255.
%   xyz, the coodrinates associated with each color. 
%   prefix, the prefix for the filenames of the xyz and python files to be
%   generated. It is good practice to include the seed used for the
%   generation of the grain colours in the prefix.
%OUTPUTS: 
%   two text files, one in xyz format and one with pymol script format that
%   will color the coordinates. 
[colorU, ~, ids] = unique(Colors, 'rows');
%generating scripts that defines the colours: 
filename_ext = strcat(prefix, '.py'); 
f = fopen(filename_ext, 'w'); 
for a = 1:numel(colorU(:,1))
    s = sprintf('set_color nick%d, [%d, %d, %d]\n', a, colorU(a,1), colorU(a,2), colorU(a,3));
    fprintf(f, '%s', s);
end 
%appending the revious script with one that colours the relevent points.
for a = 1:numel(xyz(:,1))
    s = sprintf('color nick%d, model %s and idx %d\n', ids(a), prefix, a); 
    fprintf(f, '%s', s);
end 
fclose(f);
%generating the xyzfile: 
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

