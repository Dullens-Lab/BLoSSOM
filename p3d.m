function [] = p3d( data, varargin )
%Quick way of plotting points in 3D. 
%INPUTS: data: n by 3 array of coordinates.
%        varagin: the colour and symbol code used for matlab plots.
%OUTPUTS: plots a figure in the current fugure. 
%HISTORY: Written by Nick Orr in 2016, added varargin functionality in Nov
%2017.
if numel(varargin) == 0
    plot3(data(:,1), data(:,2), data(:,3), 'xk') 
else
    plot3(data(:,1), data(:,2), data(:,3), varargin{:})
end 
axis equal 

end

