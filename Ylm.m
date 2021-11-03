function OUT = Ylm(l,phi,theta)
%Calculate spherical harmonic values
%degree l
%
%theta = colatitude, polar angle, 0 to pi
%phi = longitude, 0 to 2pi
%
%terms refer to orders -m =< l =< m

%UNNORMALISED legendre
L = legendre(l,cos(theta));
OUT = zeros(l+1,1);

COUNT = 1;
for m = 0:l
    OUT(COUNT) = ((-1)^m) * sqrt((2*l+1)/(4*pi) * factorial(l-m) / (factorial(l+m))) ...
        * exp(1i * m * phi) * L(m+1);
    %OUT(COUNT) = sqrt((2*l+1)/(4*pi) * factorial(l-m) / (factorial(l+m))) ...

    COUNT = COUNT + 1;
end

end