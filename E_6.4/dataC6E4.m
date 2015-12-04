load('geom7e4.mat')
if 1
    hold on;
    figure(1)
    for ii = 1:length(ex(:,1))
        ind = [1,2,3,4,1];
        for jj = 1:4
            plot([ex(ii,ind(jj)) , ex(ii, ind(jj + 1))],...
                 [ey(ii,ind(jj)) , ey(ii, ind(jj + 1))], 'r')
        end
    end
end