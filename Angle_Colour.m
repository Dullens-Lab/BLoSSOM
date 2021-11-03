function [ data_out, color_out ] = Angle_Colour( data, max, Ncolour_bins, pmode, tick, markersize)
%NAME: Angle_Colour
%Function: 
%   plot quantities by RGB values related to thier quantity.
%INPUT:
%   data - xyzR1R2R3 in an n by 6 array. columns 1 to 3 are the coordinates
%   in physical space and columns 4 to 6 are the RF values 1 to 3.
%   max - the maximum value which corresponds to the end of the colourmap. 
%   This is the maximum length of Rodrigues Frank vectors for this GB
%   sytmmetry. 
%   Ncolor bins, the maximum number of color bins. 
%   pmode - 1 for real space plot, 2 for orientation space plot, 0 for no
%   plot
%   tick - string input for the plot marker, same as conventional matlab
%   inputs
%   markersize - the size of the marke in pixels
%OUTPUT:
%   data_out - coloured grain boundary points
%   color_out - colours for the grain boundary points
%HISTORY: 
%   Written by Nick Orr 2018
%   Nick Orr sept 2021: added documentation. 

cbin = 1/Ncolour_bins;
%cbiny = rangey/Ncolour_bins;
%cbinz = rangez/Ncolour_bins;
L = sqrt(sum(data(:,4:6).^2,2))/max;
for a = 1:numel(L); 
    C(a) = ceil(L(a)/cbin);%/Ncolour_bins;
end 
UC = unique(C); 
%inc = 1:Ncolor_bins;
col = jet(Ncolour_bins); 
U_col = zeros(numel(UC), 3);
for a = 1:numel(UC)
    U_col(a,:) = col(UC(a), :);
end 
%U_col = jet(numel(UC)); 
if pmode == 1; 
   ii = 1; jj = 2; kk = 3;  
elseif pmode == 2
   ii = 4; jj = 5; kk = 6;
end
data_CO_c = zeros(0, numel(data(1,:)));
color_CO_c = zeros(0, 3);
for a = 1:numel(U_col(:,1)); 
    if pmode ~= 0 
        plot3(data(C == UC(a), ii),...
            data(C == UC(a), jj),...
            data(C == UC(a), kk),...
            tick, 'MarkerFaceColor', U_col(a,:), 'Color', U_col(a,:), 'MarkerSize', markersize)
    end
    data_CO = data(C == UC(a), :);
    data_CO_c = [data_CO_c; data_CO];
    color_CO = U_col(a, :); 
    re_color_CO = repmat(color_CO, numel(data_CO(:,1)), 1); 
    color_CO_c = [color_CO_c; re_color_CO]; 
end 
data_out = data_CO_c;
color_out = color_CO_c;% * 255;
end

