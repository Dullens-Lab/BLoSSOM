function [ Angrr, random_dist, xyDist, xyRandDist ] = Disorientation_Angle_Histogram( dis, s1, s2, bins, Angrr)
%NAME: Disorientation_Angle_Histogram
%FUNCTION: Create a histogram of a disorientation angle disribution. It is
%like a radial distribution plot for disorientation space. It is useful for
%comparing distibutions, and highlihting key disorientations.
%INPUTS:
%   dis - n by 3 array containing tthe disorientation vectors in RF
%   parameters. 
%   s1 - symmetry 1, hard coded 1 is for O, 2 is for D3 and 3 is for D6 point groups.  
%   s2 - symmetry 2
%   bins
%OUTPUTS:
%   plots a histogram in a matlab figure.
%   Angrr - the distribution for random rotations
%   random_dist - the x y histogram coordinates.
%   xyDist - the distribution vlues to be plotted
%   xyRandDist - the random distribution to be plotted
%HISTORY: 
%Written by Nick Orr, 18/01/19. First frost of the winter this morning!
%Nick Orr 10 Dec 2020: Added the functionality to output the plotting 
%variables xyDist and xyRandDist
Ang = zeros(numel(dis(:,1)), 1);
for a = 1:numel(dis(:,1))
    modulus = sqrt(sum(dis(a,:).^2)); 
    Ang(a) = atand(modulus) * 2;
end
N = numel(dis(:,1));
if Angrr == 0;
    for a = 1:50000000 %this brute force calculates the random distribution. 
        %It takes a long time but you can save the variable. 
        rr1 = random_rotation(3);
        rr2 = random_rotation(3);
        M12 = MisorientationRF4(rr1, rr2, s1, s2);
        modulus = sqrt(sum(M12.^2)); 
        Angrr(a) = atand(modulus) * 2;
    end
end
%crash
[counts,centers] = hist(Ang, bins, '-k');
[countsrr,centersrr] = hist(Angrr, bins, 'm');
random_dist = [centersrr; countsrr/numel(Angrr)];
% whos centers;
% whos centersrr;
        figure
        hold on 
plot(centers, counts/(numel(Ang)*(bins(2)-bins(1))),'-k', 'LineWidth', 1);% (bins(2) - bins(1))/2
plot(centersrr, countsrr/(numel(Angrr) *(bins(2)-bins(1))), '--k', 'LineWidth', 1)
% figure
% plot(centers, (counts/N) - (countsrr/numel(Angrr)), '-k', 'LineWidth', 1)
%NEED TO MAKE THE HISTOGRAM BINS THE SAME
xyDist = [centers;counts/(numel(Ang)*(bins(2)-bins(1)))];
xyRandDist = [centersrr; countsrr/(numel(Angrr) *(bins(2)-bins(1)))];
end


