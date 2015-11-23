% Newton rhapson solver for exercise 6.1
run('data.m')

% Defines additional material properties
E = 210*1e9; %[Pa]; ?? D-matrisen st?mmer d? med control.mat
v = 0.3;

D = (E/((1+v)*(1-2*v))) .* [1-v, v, 0;   v, 1-v, 0; 0, 0, (1-2*v)/2];

nelm = length(edof(:,1));
ndof = max(max(edof));

nmax = 10;
imax = 20;
LIMIT = 1e-4;

u = zeros(ndof, 1);
f = zeros(ndof,1);
df = zeros(ndof,1);
%forceIndices=[321;365;367;369;305];
%df(forceIndices) = 1e9;
%bc = bc(1:end-5,:);
%bc(end-4:end,2) = 1000;

t = 0.1; % OBS! arbitrarily set

uval = [];
fval = [];

for n = 1:nmax
    %Applies force
    f = f + df;
    
    for i = 1:imax
        % Creates the global stiffness matrix
        K = zeros(ndof);
        for ii = 1:nelm
            ec = [ex(ii,:)', ey(ii,:)'];
            ed = u(edof(ii,2:7));
            [ ee , ~ ] = plan3gs(ec' , ed);
            es = D * ee; % Compute second order poila-kirchoff stress ??
            Ke = plan3ge( ec' , t , D , ed , es );
            K(edof(ii,2:7),edof(ii,2:7)) = K(edof(ii,2:7),edof(ii,2:7)) + Ke;
        end
        
        fint = zeros(ndof,1);
        for jj = 1:nelm
            ec = [ex(ii,:)', ey(ii,:)']';
            
            [ee, ~] = plan3gs( ec , u(edof(ii,2:7)));
            A = 0.5*abs(ec(1,1)*(ec(2,2)-ec(2,3))+...
                ec(1,2)*(ec(2,3)-ec(2,1))+...
                ec(1,3)*(ec(2,1)-ec(2,2)));
            es = A.*ee;
            finte = plan3gf( ec, t, ed , es );
            fint(edof(ii,2:7)) = fint(edof(ii,2:7)) + finte;
        end

        r = f - fint;
        du = solveq(K , r , bc);
        u = u + du;

        if norm(r) < LIMIT * norm(df)
            break
        end
        if i == imax
            disp(['Warning in iteration n =', num2str(n),'. Reached iteration limit at ||r||_2 = ', num2str(norm(r)), ' without convergence.'])
        end
    end
    uval = [uval, u];
    fval = [fval, f];
end
%% Computes the deformed x and y coordinates
defx = ex;
defy = ey;
for ii = 1:nelm
    defx(ii,:) = ex(ii,:) + u(edof(ii, 2:2:6))';
    defx(ii,:) = ex(ii,:) + u(edof(ii, 3:2:7))';
end
hold on;
 for ii = 1:length(defx(:,1))
    ind = [1,2,3,1];
    for jj = 1:3
        plot([defx(ii,ind(jj)) , defx(ii, ind(jj + 1))],...
             [defy(ii,ind(jj)) , defy(ii, ind(jj + 1))], 'b')
    end
end