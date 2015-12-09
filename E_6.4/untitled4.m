ec = [86.7312   99.0405   97.9542   85.3416;
      39.6184   46.5523   58.1482   49.4416];
  
hold on;
x = ec(1,:);
y = ec(2,:);
plot3(x,y,zeros(1,length(x)),'*')
Ne = @(a,b) ((a-x(2))/(x(1) - x(2)))*((b-y(4))/(y(1) - y(4)));
sx = linspace(min(x),max(x),10);
sy = linspace(min(y),max(y),10);
Z = zeros(10,10);
for ii = 1:length(sx)
    for jj = 1:length(sy)
        Z(ii,jj) = Ne(sx(ii),sy(jj));
    end
end
surf(sx,sy,Z)
hold off;