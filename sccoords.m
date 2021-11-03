function [sc] = sccoords(hi)
%generate sccoordinates, simple cubic coordiantes. 

d1 = 0:hi;
d2 = repmat(d1, hi+1, 1);
d3 = reshape(d2, numel(d2),1);
d4 = d2'; 
d5 = reshape(d4, numel(d3), 1);
d35 = [d3,d5]; 
d35r = repmat(d35, hi+1, 1); 
d6 = repmat(d1, (hi+1)^2, 1);
d7 = reshape(d6, numel(d6), 1); 
d35r7 = [d35r, d7]; 
sc = d35r7;
end

