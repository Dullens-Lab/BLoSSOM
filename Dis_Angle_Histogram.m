function [xyDist] = Dis_Angle_Histogram( dis, bins, Angrr, LW)
%NAME: Disorientation_Angle_Histogram
%FUNCTION: Create a histogram of a disorientation angle disribution. It is
%like a radial distribution plot for disorientation space. It is useful for
%comparing distibutions, and highlihting key disorientations.
%INPUTS:
%   dis - n by 3 array containing tthe disorientation vectors in RF
%   parameters. 
%   bins - the bins for the histogram
%   Angrr - the random misorientation distribution, as a n by 2 array of
%   coordinates.
%   LW - the line width for the plot
%OUTPUTS:
%   plots a histogram in a matlab figure.
%   xyDist - the distribution vlues to be plotted
%HISTORY: 
%Written by Nick Orr, 18/01/19. First frost of the winter this morning!
%Nick Orr 10 Dec 2020: Added the functionality to output the plotting 
%variables xyDist and xyRandDist
%Nick Orr 03 11 2021, changed so it works with the new, lighter Angrr
%variable saves. 
Ang = zeros(numel(dis(:,1)), 1);
for a = 1:numel(dis(:,1))
    modulus = sqrt(sum(dis(a,:).^2)); 
    Ang(a) = atand(modulus) * 2;
end
[counts,centers] = hist(Ang, bins, '-k');

        figure
        hold on 
plot(centers, counts/(numel(Ang)*(bins(2)-bins(1))),'-k', 'LineWidth', LW);% (bins(2) - bins(1))/2
plot(Angrr(:,1), Angrr(:,2), '--k', 'LineWidth', LW)
xyDist = [centers;counts/(numel(Ang)*(bins(2)-bins(1)))];
end


