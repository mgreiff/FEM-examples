function eldrawTri2D(ex, ey, edof)
    hold on;
    indices = [1,2,3,1];
    for ii = 1:length(edof(:,1))
        for jj = 1:3
            plot([ex(ii,indices(jj)),ex(ii,indices(jj + 1))],...
                 [ey(ii,indices(jj)),ey(ii,indices(jj + 1))], 'b')
        end
    end
    % Plots vertex where pressure is to be applied
    pressureIndex = 64;
    plot(ex(pressureIndex,1),ey(pressureIndex,1),'r*')
    
    % Plots boundary conditions
    plot(ex(1,1),ey(1,1),'kx')
    plot(ex(1,1),ey(1,1),'g+')
    plot(ex(end,2),ey(end,2),'kx')
    plot(ex(end,2),ey(end,2),'g+')
    hold off;
    
end