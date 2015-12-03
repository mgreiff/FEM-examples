load('geom7e1.mat')
% Transforms the problem into four node elements

E = 210*1e9;
v = 0.3;
t = 0.1; % OBS! arbitrarily set
D = (E/((1+v)*(1-2*v))) .* [1-v, v, 0;   v, 1-v, 0; 0, 0, (1-2*v)/2];

nelm = length(edof(:,1));
ndof = max(max(edof));

% Plots the original geometry
if 0
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
edofHat = [];
coord = [];
for ii  = (1:4)
    coord = [coord, [ex(ii,:); ey(ii,:)]];
end
coordHat = unique(coord' , 'rows')';
plot(coordHat(1,:),coordHat(2,:),'*')