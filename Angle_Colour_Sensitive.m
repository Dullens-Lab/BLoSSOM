function [ data_out, color_out ] = Angle_Colour_Sensitive( data, minmax, Ncolour_bins, pmode, tick, markersize)
%NAME: RGB_Colour
%Function: 
%   plot quantities by RGB values related to thier quantity.
%INPUTS:
%   data - xyzR1R2R3 in an n by 6 array. columns 1 to 3 are the coordinates
%   in physical space and columns 4 to 6 are the RF values 1 to 3.
%   Ncolor bins, the maximum number of color bins. 
%   pmode - 1 for real space plot, 2 for orientation space plot, 0 for no
%   plot
%   tick - string input for the plot marker, same as conventional matlab
%   inputs
%   markersize - the size of the marke in pixels


cbin = 1/Ncolour_bins;
%cbiny = rangey/Ncolour_bins;
%cbinz = rangez/Ncolour_bins;
L = atand(sqrt(sum(data(:,4:6).^2,2))) * 2; %/max;
L(L<minmax(1)) = minmax(1);
L(L>minmax(2)) = minmax(2);
L = L-minmax(1);
L = L/(minmax(2)-minmax(1)); %nick orr july 2021
%Lmin =  tand(2/2)
%L(L>1) =1; 
for a = 1:numel(L); 
    C(a) = ceil(L(a)/cbin);%/Ncolour_bins;
end 
C(C==0) = 1; % maybe a fudge? 
UC = unique(C); 
%inc = 1:Ncolor_bins;
col = parula(Ncolour_bins); %jet
U_col = zeros(numel(UC), 3);
for a = 1:numel(UC)
    U_col(a,:) = col(UC(a), :);
end 
%U_col = jet(numel(UC)); 
if pmode == 1; 
   ii = 1; jj = 2; kk = 3;  
elseif pmode == 2
   ii = 4; jj = 5; kk = 6;
elseif pmode == -1
   ii = 1; jj = 2; kk = 3;  
end
data_CO_c = zeros(0, numel(data(1,:)));
color_CO_c = zeros(0, 3);
for a = 1:numel(U_col(:,1)); 
    if a == 1; 
        hold off %edit july 2021
    else 
        hold on
    end 
    if pmode > 0 
        plot3(data(C == UC(a), ii),...
            data(C == UC(a), jj),...
            data(C == UC(a), kk),...
            tick, 'MarkerFaceColor', U_col(a,:), 'Color', U_col(a,:), 'MarkerSize', markersize)
    elseif pmode  == -1
        hold on 
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

