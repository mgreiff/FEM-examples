clear all; close all;
load('geom7e1.mat')

% Scales geometry to SI-units
ex = ex.*1e-3;
ey = ey.*1e-3;

if 0
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
end