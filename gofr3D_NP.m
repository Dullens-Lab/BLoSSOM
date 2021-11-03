function [g,xx] = gofr3D_NP(dat,Lxy,Lz,res,minR,maxR)
%NAME: gofr3D_NP
%Function: generate a g of r radisal distribution function from a set of
%points.
%INPUTS: 
%   dat - n by 3 array of points. 
%   Lxy - the length of the box containing the data in x and y (must be 
%   square).
%   Lz - the height of the box. 
%   res - resolution of the g of r. i.e. how many bins of your histogram.
%   minR - the minimum distance of the g of r.
%   maxR - the maximum distance of the g of r. 
%OUTPUTS: 
%   g - the g of r 
%   xx - the bin centers. 
%HISTORY: 
%Written by Taiki Yanagishima, inherited to Nick Orr as
%gofr3D_NonPeriodic.m 
%Nick Orr 04/04/20 changed the input to be an n by 3 array. 
%Nick Orr 04/04/20 changed the name of the function to gof3D_NP.
%Nick Orr 04/04/20 Added documentation.


X = reshape(dat',(size(dat,1)*3), 1); %reshape it for legacy reasons
n = max(size(X))/3;
disp(n);
g = zeros(res,1);
edge = maxR;
dD = (maxR - minR)/res;

xx = (linspace(1,res,res)'-0.5)*dD+minR;
count =0; 
TOTN = 0;
for I = 1:n-1
    if X(3*I-2) > edge && X(3*I-2) < Lxy-edge && ...
        X(3*I-1) > edge && X(3*I-1) < Lxy-edge && ...
        X(3*I) > edge && X(3*I) < Lz-edge
        for J = I+1:n
            vec = X(3*I-2:3*I)-X(3*J-2:3*J);  % Get vector of pairs
            d = norm(vec);
            if d < maxR && d > minR
                II = ceil((d - minR)/dD);
                g(II) = g(II) + 2;
            end
        end
        TOTN = TOTN + 1;
    end
    if count > 100
        count = 0;
        I
    end 
    count = count + 1;
end

v = zeros(res,1);
for I = 1:res
v(I) = (I^3 - (I-1)^3) * dD^3;
end
nid = 4/3 * pi * v * (TOTN/((Lxy-2*edge)^2 * (Lz-2*edge)));
g = g/TOTN ./nid;

end