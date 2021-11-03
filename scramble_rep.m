function [ scrambled ] = scramble_rep( coords, number, seed )
%NAME: scramble_rep
%FUNCTION: to scramble to order of points in an array
%INPUTS: 
%   coords - The coordinates of the points, each new row is a new point,
%   the number of columns is the dimention. 
%   number - The number of points to be returned. If numbe is set to less
%   than the orieginal number of points, random points will be ommitted. 
%   seed - positive inter number to seed the rng.
%OUTPUTS: 
%   scrambled - the scrambles array with "number" of points in the array. 
%HISTROY: 
%Written by Nick Orr in Jan 2018. 
%s = rng(seed); 
npoints = numel(coords(:,1)); 
rng(seed)
random_numbers = rand(1,npoints); 
[~, order_random] = sort(random_numbers); 
scrambled = coords(order_random(1:number), :); 
end

