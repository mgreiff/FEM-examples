clear all; close all;
load('geom7e1.mat')

% Scales geometry
ex = ex;
ey = ey;

%% Plots the elements
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

    % Plots the vrtices
    %exhat = reshape(ex, [numel(ex), 1]);
    %eyhat = reshape(ey, [numel(ey), 1]);
    %plot(exhat, eyhat, '*r')
    
    % plots the boundary conditins
    for ii = 1:length(edof(:,1))
        for jj = 1:3
            if find(bc(:,1) == edof(ii,jj*2))
                plot(ex(ii,jj),ey(ii,jj), 'gx')
            end
            if find(bc(:,1) == edof(ii,jj*2+1))
                plot(ex(ii,jj),ey(ii,jj), 'r+')
            end
        end
    end
    hold off;
end