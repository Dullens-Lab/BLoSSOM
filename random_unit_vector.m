function an=random_unit_vector(varargin)
% random_unit_vector
% random_unit_vector(n)
% random_unit_vector(m,n)
% random_unit_vector([m n])
% random_unit_vector('double')
% random_unit_vector(n,'double')
% random_unit_vector(m,n,'double')
% random_unit_vector([m n],'double')
% random_unit_vector('single')
% random_unit_vector(n,'single')
% random_unit_vector(m,n,'single')
% random_unit_vector([m n],'single')

% m - dimentionarity 
% n - number of unit vectors

md=3; % default m

isg=false; % if single
switch nargin
    case 0
        % if no inputs
        n=1;
        m=md;
    case 1
        %
        i1=varargin{1};
        if ischar(i1)
            if strcmpi(i1,'single')
                isg=true;
            end
        else
            if length(i1)==1
                m=md;
                n=i1;
            else
                m=i1(1);
                n=i1(2);
            end
        end
    case 2
        i1=varargin{1};
        i2=varargin{2};
        if ischar(i2)
            if length(i1)==1
                m=md;
                n=i1;
            else
                m=i1(1);
                n=i1(2);
            end
            if strcmpi(i2,'single')
                isg=true;
            end
        else
            m=i1;
            n=i2;
        end
    case 3
        m=varargin{1};
        n=varargin{2};
        if strcmpi(varargin{3},'single')
            isg=true;
        end
            
        
end

% simple case of 1d
if m==1
    if isg
        an=single(2*(randn(1,n)>0)-1);
    else
        an=2*(randn(1,n)>0)-1;
    end
    return;
end

if isg
    v=randn(m,n,'single');
else
    v=randn(m,n);
end

% normalize:
if isg
    an=zeros(m,n,'single');
else
    an=zeros(m,n);
end
for nc=1:n
    while 1
        v2=v(:,nc)'*v(:,nc);
        if v2>1e-10 % too small values must be excluded 
            % because it will have discretization errors
            an(:,nc)=v(:,nc)/sqrt(v2);
            break;
        else
            % otherwise repeat random generation
            if isg
                v(:,nc)=randn(m,1,'single');
            else
                v(:,nc)=randn(m,1);
            end
        end
    end
end

%Copyright (c) 2009, Maxim Vedenyov
%All rights reserved.

%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:

%* Redistributions of source code must retain the above copyright notice, this
%  list of conditions and the following disclaimer.

%* Redistributions in binary form must reproduce the above copyright notice,
%  this list of conditions and the following disclaimer in the documentation
%  and/or other materials provided with the distribution
%* Neither the name of http://simulations.narod.ru/ nor the names of its
%  contributors may be used to endorse or promote products derived from this
%  software without specific prior written permission.
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
%FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
%DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
%SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
%CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
