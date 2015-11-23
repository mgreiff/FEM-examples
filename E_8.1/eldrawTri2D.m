function eldrawTri2D(ex, ey, edof)
    hold on;
    indices = [1,2,3,1];
    for ii = 1:length(edof(:,1))
        for jj = 1:3
            plot([ex(ii,indices(jj)),ex(ii,indices(jj + 1))],...
                 [ey(ii,indices(jj)),ey(ii,indices(jj + 1))], 'b')
        end
    end
    hold off;
end