function [out,outCOMP] = BOO(boo,X,L,ll,input,NNmode)
%X is a 3n long matrix
%L is period of box
%
%NNmode = 1,2
%1 - input is NN list from NNlist2
%2 - input is file name of Voronoi face data


switch boo

    case 'q'
        
        n = ms(X)/3;
        out = zeros(n,1);
        outCOMP = cell(n,1);
        
        switch NNmode
            %case 1 - neighbour list based on threshold, max set
            case(1)
                
                %Need a neighbour list from NNlist2
                %Make sure it is in NN with the right struct.
                NN = input;
                
                %Calculate q_6 for each particle
                for I = 1:n
                    qlm = zeros(ll+1,1);
                    for J = 1:NN.LISTN(I)
                        n1 = NN.LIST(I,J);
                        dx = X(3*n1-2:3*n1) - X(3*I-2:3*I);
                        dx = dx - L * round(dx/L);
                        
                        [phi,theta,~] = cart2sph(dx(1),dx(2),dx(3));
                        theta = pi/2 - theta;
                        %phi = phi + pi;
                        qlm = qlm + Ylm(ll,phi,theta);
                    end
                    qlm = qlm / NN.LISTN(I);
                    %qlm is a vector, complex terms in different orders, 0 to m
                    %For Associated legendre, m range from 0 to m, not -m to m
                    
                    outCOMP{I} = qlm;
                    %Sum over orders
                    qlm2 = 0;
                    %Start with terms 1 to l, double to account for -ve m values
                    for J = 2:ll+1
                        qlm2 = qlm2 + 2 * abs(qlm(J))^2;
                    end
                    %Add m = 0 contribution afterwards
                    qlm2 = qlm2 + abs(qlm(1))^2;
                    
                    out(I) = sqrt(qlm2 * 4*pi/(2*ll+1));
                end
                
                
                %Voronoi tesselation
            case(2)
                
                stream = input{1};
                
                P = 1;
                VDATAn = 1;
                
                while P <= n
                    %Skip particle index
                    VDATAn = VDATAn + 1;
                    nn = stream(VDATAn);
                    
                    %Nlist contains indices of neighbours of particle P
                    Nlist = zeros(nn,1);
                    for I = 1:nn
                        VDATAn = VDATAn + 1;
                        Nlist(I) = stream(VDATAn);
                    end
                    
                    %Flist contains
                    Flist = zeros(nn,1);
                    for I = 1:nn
                        VDATAn = VDATAn + 1;
                        Flist(I) = stream(VDATAn);
                    end
                    FA = sum(Flist);
                    
                    qlm = zeros(ll+1,1);
                    for J = 1:nn
                        n1 = Nlist(J);
                        dx = X(3*n1-2:3*n1) - X(3*P-2:3*P);
                        dx = dx - L * round(dx/L);
                        [phi,theta,~] = cart2sph(dx(1),dx(2),dx(3));
                        theta = pi/2 - theta;
                        %phi = phi + pi;
                        qlm = qlm + Flist(J) / FA * Ylm(ll,phi,theta);
                    end
                    %qlm = qlm / nn;
                    %qlm is a vector, complex terms in different orders, 0 to m
                    %For Associated legendre, m range from 0 to m, not -m to m
                    outCOMP{P} = qlm;
                    
                    %Sum over orders
                    qlm2 = 0;
                    %Start with terms 1 to l, double to account for -ve m values
                    for J = 2:ll+1
                        qlm2 = qlm2 + 2 * abs(qlm(J))^2;
                    end
                    %Add m = 0 contribution afterwards
                    qlm2 = qlm2 + abs(qlm(1))^2;
                    out(P) = sqrt(qlm2 * 4*pi/(2*ll+1));
                    
                    %Next particle...
                    P = P + 1;
                    VDATAn = VDATAn + 1;
                end
                
        end

        
        
        
    case 'Q'
        
        n = ms(X)/3;
        out = zeros(n,1);
        outCOMP = cell(n,1);
        
        switch NNmode
            %case 1 - neighbour list based on threshold, max set
            case(1)
                
                %Need a neighbour list from NNlist2
                %Make sure it is in NN with the right struct.
                NN = input;
                
                %Calculate q_6 for each particle
                qlm = zeros(ll+1,n);
                for I = 1:n
                    for J = 1:NN.LISTN(I)
                        n1 = NN.LIST(I,J);
                        dx = X(3*n1-2:3*n1) - X(3*I-2:3*I);
                        dx = dx - L * round(dx/L);
                        
                        [phi,theta,~] = cart2sph(dx(1),dx(2),dx(3));
                        theta = pi/2 - theta;
                        %phi = phi + pi;
                        qlm(:,I) = qlm(:,I) + Ylm(ll,phi,theta);
                    end
                    qlm(:,I) = qlm(:,I) / NN.LISTN(I);
                    %qlm is a vector, complex terms in different orders, 0 to m
                    %For Associated legendre, m range from 0 to m, not -m to m
                end
                
                for I = 1:n
                    qlmav = zeros(ll+1,1);
                    for J = 1:NN.LISTN(I)
                        n1 = NN.LIST(I,J);
                        qlmav = qlmav + qlm(:,n1);
                    end
                    qlmav = qlmav + qlm(:,I);
                    qlmav = qlmav / (NN.LISTN(I)+1);
                    
                    outCOMP{I} = qlmav;
                    
                    %Sum over orders
                    qlm2 = 0;
                    %Start with terms 1 to l, double to account for -ve m values
                    for J = 2:ll+1
                        qlm2 = qlm2 + 2 * abs(qlmav(J))^2;
                    end
                    %Add m = 0 contribution afterwards
                    qlm2 = qlm2 + abs(qlmav(1))^2;
                    out(I) = sqrt(qlm2 * 4*pi/(2*ll+1));
                end
                
                
                %Voronoi tesselation
            case(2)
                
                stream = input{1};
                
                P = 1;
                VDATAn = 1;
                qlm = zeros(ll+1,n);
                nn = zeros(n,1);
                NlistMEM = cell(n,1);
                
                while P <= n
                    %Skip particle index
                    VDATAn = VDATAn + 1;
                    nn(P) = stream(VDATAn);
                    
                    %Nlist contains indices of neighbours of particle P
                    Nlist = zeros(nn(P),1);
                    for I = 1:nn(P)
                        VDATAn = VDATAn + 1;
                        Nlist(I) = stream(VDATAn);
                    end
                    NlistMEM{P} = Nlist;
                    
                    %Flist contains
                    Flist = zeros(nn(P),1);
                    for I = 1:nn(P)
                        VDATAn = VDATAn + 1;
                        Flist(I) = stream(VDATAn);
                    end
                    FA = sum(Flist);
                    
                    for J = 1:nn(P)
                        n1 = Nlist(J);
                        dx = X(3*n1-2:3*n1) - X(3*P-2:3*P);
                        dx = dx - L * round(dx/L);
                        [phi,theta,~] = cart2sph(dx(1),dx(2),dx(3));
                        theta = pi/2 - theta;
                        %phi = phi + pi;
                        qlm(:,P) = qlm(:,P) + Flist(J) / FA * Ylm(ll,phi,theta);
                    end
                    %qlm = qlm / nn;
                    %qlm is a vector, complex terms in different orders, 0 to m
                    %For Associated legendre, m range from 0 to m, not -m to m
                    
                    %Next particle...
                    P = P + 1;
                    VDATAn = VDATAn + 1;
                end
                
                %For Minkowski metric, sum over NN neighbours as judged by
                %Voronoi (?)
                for I = 1:n
                    qlmav = zeros(ll+1,1);
                    for J = 1:nn(I)
                        n1 = NlistMEM{I}(J);
                        qlmav = qlmav + qlm(:,n1);
                    end
                    qlmav = qlmav + qlm(:,I);
                    qlmav = qlmav / (nn(I)+1);
                    
                    outCOMP{I} = qlmav;
                    
                    %Sum over orders
                    qlm2 = 0;
                    %Start with terms 1 to l, double to account for -ve m values
                    for J = 2:ll+1
                        qlm2 = qlm2 + 2 * abs(qlmav(J))^2;
                    end
                    %Add m = 0 contribution afterwards
                    qlm2 = qlm2 + abs(qlmav(1))^2;
                    out(I) = sqrt(qlm2 * 4*pi/(2*ll+1));
                end
                
        end
          
end

end