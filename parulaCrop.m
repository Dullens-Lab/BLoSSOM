function [colours] = parulaCrop(Ncolours,frac,res)
%NAME: parulaCrop
%Function: produce colours from a cropped parula colourmap to remove the
%yellows at the top 
%INPUTS: 
% Ncolours, the number of discrete colours returned in the colourmap
% frac, the fraction of the colourmap you want to use should be a number
% between 0 and 1. I reccomend 0.85
% res, the resolution of the colourmap parula, set to a high value.
%OUTPUTS: 
% colours, the colours in the cropped parula map. 

N = Ncolours; 
coloursL = parula(res/frac);
coloursC = coloursL(1:res,:);
I_1 = 1:999/(N-1):1000;
I_2 = round(I_1 - I_1(1)/2)';
colours = coloursC(I_2,:);

end

