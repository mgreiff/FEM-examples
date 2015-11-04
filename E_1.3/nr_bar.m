clear all;
run data.m
nmax = 20;
imax = 50;
LIM = 1e-4;

fval = [];
uval = [];
gval = [];
uold = 0;
f = 0;
for n = 1:nmax
    u = uold;
    f = f + deltaf - k * uold;
    fval = [fval, f];
    i = 1;
    for i =1:imax
        K = kt_bar( u );
        g = g_bar( u );
        r = f - g;
        deltau = K \ r;
        u = u + deltau;
        if norm(r, 2) < LIM * norm(f, 2)
            uold = u;
            uval = [uval, uold];
            gval = [gval, g];
            break
        end
        if i == imax
            disp('Reached iteration limit')
        end
    end
end
figure(1);
plot(uval, fval)
xlabel('u')
ylabel('f(u)')

figure(2);
plot(-uval/a, -(gval/EA)*(l0/a).^3)
xlabel('-u/a')
ylabel('-(g(u)/EA)*(l0/a)^3')
