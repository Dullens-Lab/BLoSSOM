function [Qavg]=avg_quaternion_markley(Q)
% by Tolga Birdal
% Q is an Mx4 matrix of quaternions. Qavg is the average quaternion
% Based on 
% Markley, F. Landis, Yang Cheng, John Lucas Crassidis, and Yaakov Oshman. 
% "Averaging quaternions." Journal of Guidance, Control, and Dynamics 30, 
% no. 4 (2007): 1193-1197.
% Form the symmetric accumulator matrix
A = zeros(4,4);
M = size(Q,1);

for i=1:M
    q = Q(i,:)';
    if(q(1)<0) % handle the antipodal configuration
		q = -q;
	end
    A = q*q'+A; % rank 1 update
end

% scale
A=(1.0/M)*A;

% Get the eigenvector corresponding to largest eigen value
[Qavg, Eval] = eigs(A,1);

end

%The MIT License (MIT)

%Copyright (c) 2014 Tolga Birdal

%Permission is hereby granted, free of charge, to any person obtaining a copy
%of this software and associated documentation files (the "Software"), to deal
%in the Software without restriction, including without limitation the rights
%to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all
%copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%SOFTWARE.