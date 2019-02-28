close all;

%% Geometry and eslasticity
b = 1.0;
a = 0.1;
l0 = sqrt(a.^2 + b.^2);
EA = 1.0;

%% Solver parameters
nmax=40;
imax=8;
eps=10^-5;
deltaf = -1e-4;
uval=[];
fval=[];
fnew=0;
unew=0;

%% Find the force-displacement curve
for n=1:nmax;
    fnew=fnew+deltaf;
    for i=1:imax
        K      = kt_bar(unew, a, b, l0, EA);
        res    = fnew - g_bar(unew, a, b, l0, EA);
        deltau = K\res;
        unew   = unew + deltau;
        if norm(res) < eps * norm(deltaf)
            break
        end
    end
    uval=[uval abs(unew/a)];
    fval=[fval abs(fnew*((l0/a)^3/EA))];
end

%% Visualize results
figure(1);
plot(uval,fval,'bo:')
title('Two bars when increasing force / force-displacement curve', 'Interpreter', 'latex')
xlabel('u', 'Interpreter', 'latex')
ylabel('f(u)', 'Interpreter', 'latex')

%% Save plot to the figures directory
saveas(gcf,'../../figures/example_1_3.eps','epsc')
