clear all;
close all;
load('geom7e1.mat')

% Defines additional material properties
E = 210; %[GPa];
v = 0.3;
% The geometry is defined in [mm]

hold on;
figure(1)

% Displays geometry
for ii = 1:length(ex(:,1))
    ind = [1,2,3,1];
    for jj = 1:3
        plot([ex(ii,ind(jj)) , ex(ii, ind(jj + 1))],...
             [ey(ii,ind(jj)) , ey(ii, ind(jj + 1))], 'b')
    end
end

exhat = reshape(ex, [numel(ex), 1]);
eyhat = reshape(ey, [numel(ey), 1]);
plot(exhat, eyhat, '*r')
hold off;