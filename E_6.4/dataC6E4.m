load('geom7e1.mat')
% Transforms the problem into four node elements

E = 210*1e9;
v = 0.3;
t = 0.1; % OBS! arbitrarily set
D = (E/((1+v)*(1-2*v))) .* [1-v, v, 0;   v, 1-v, 0; 0, 0, (1-2*v)/2];

% Plots the original geometry
if 1
    hold on;
    figure(1)
    for ii = 1:length(ex(:,1))
        ind = [1,2,3,1];
        for jj = 1:3
            plot([ex(ii,ind(jj)) , ex(ii, ind(jj + 1))],...
                 [ey(ii,ind(jj)) , ey(ii, ind(jj + 1))], 'b')
        end
    end
end

% Transforms geometry to four-node elements
edofHat = [];
exHat = [];
eyHat = [];
nelm = length(edof(:,1))/4;
for n = 1:nelm
    dofs = [];
    for ii  = 4*(n - 1) + (1:4)
        dofs = [dofs , edof(ii,2:7)];
    end

    validDofs= [];
    for dof = unique(dofs)
        if sum(dofs == dof) == 2
            validDofs = [validDofs, dof];
        end
    end

    validCoords = [];
    for ii  = 4*(n - 1) + (1:4)
        for jj = 1:3
            if sum(edof(ii,2*jj) == validDofs) > 0
                validCoords = [validCoords, [ex(ii,jj);ey(ii,jj)]];
            end
        end
    end
    validCoords = unique(validCoords','rows')';
    edofHat = [edofHat;n, validDofs];
    exHat = [exHat;validCoords(1,:)];
    eyHat = [eyHat;validCoords(2,:)];
end
if 1
    hold on;
    figure(1)
    for ii = 1:length(exHat(:,1))
        ind = [1,2,3,4,1];
        for jj = 1:4
            plot([exHat(ii,ind(jj)) , exHat(ii, ind(jj + 1))],...
                 [eyHat(ii,ind(jj)) , eyHat(ii, ind(jj + 1))], 'r--')
        end
    end
end