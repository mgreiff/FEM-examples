function [ Ke ] = plan4gie( ec , t , D , ed , es )
    x = ec(1,:);
    y = ec(2,:);
    
    Ke = zeros(8,8);
    
    p = [-sqrt(1/3),sqrt(1/3)];

    wxi = [1,1,1,1];
    wmu = [1,1,1,1];
    
    xipoints = [p,-p];
    mupoints = [p(1),p(1),p(2),p(2)];
    
    for ii  = 1:4
        xi = xipoints(ii);
        mu = mupoints(ii);
        
        dNxi = [mu-1,-mu+1,mu+1,-mu-1]./4;
        dNmu = [xi-1,-xi-1,xi+1,-xi+1]./4;

        J11 = x*dNxi';
        J12 = y*dNxi';
        J21 = x*dNmu';
        J22 = y*dNmu';

        Jdet = J11*J22-J12*J21;
        
        Jinv = [J22,-J12;-J21,J11]./Jdet;
        
        gradN = zeros(2,4);
        for jj = 1:4
            gradN(:,jj) = Jinv*[dNxi(jj);dNmu(jj)]; 
        end

        Ndx = gradN(1,:);
        Ndy = gradN(2,:);
        
        B0 = [Ndx(1),0     ,Ndx(2),0     ,Ndx(3),0     ,Ndx(4),0     ;
              0     ,Ndy(1),0     ,Ndy(2),0     ,Ndy(3),0     ,Ndy(4);
              Ndy(1),Ndx(1),Ndy(2),Ndx(2),Ndy(3),Ndx(3),Ndy(4),Ndx(4)];

        H = [Ndx(1),0     ,Ndx(2),0     ,Ndx(3),0     ,Ndx(4),0     ;
             Ndy(1),0     ,Ndy(2),0     ,Ndy(3),0     ,Ndy(4),0     ;
             0     ,Ndx(1),0     ,Ndx(2),0     ,Ndx(3),0     ,Ndx(4);
             0     ,Ndy(1),0     ,Ndy(2),0     ,Ndy(3),0     ,Ndy(4)];

        duxdx = Ndx*ed(1:2:7)';
        duxdy = Ndy*ed(1:2:7)';
        duydx = Ndx*ed(2:2:8)';
        duydy = Ndy*ed(2:2:8)';

        A = [duxdx,0    ,duydx,0    ;
             0    ,duxdy,0    ,duydy;
             duxdy,duxdx,duydy,duydx];

        B = B0 + A*H;
        
        S = [es{ii}(1), es{ii}(3);
             es{ii}(3), es{ii}(2)];

        R=[S, eye(2); eye(2), S];

        Ke = Ke + wxi(ii).*wmu(ii).*(B'*D{ii}*B+H'*R*H)*Jdet*t;
    end
end

  