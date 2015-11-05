clear all;
close all;
run data.m;
nmax=40;
imax=8;
eps=10^-5;
uval=[];
fval=[];
fnew=0;
unew=0;
for n=1:nmax;
    fnew=fnew+deltaf;
    %disp(['loadstep ',num2str(n)])
    for i=1:imax
        K=kt_bar(unew);
        res=fnew-g_bar(unew);
        deltau=K\res;
        unew=unew+deltau;
        if norm(res) < eps * norm(deltaf)
            break
        end
        %disp(['iterationstep ',num2str(norm(res))])
    end
    uval=[uval abs(unew/a)];
    fval=[fval abs(fnew*((l0/a)^3/EA))];
end
plot(uval,fval,'*')
xlabel('u')
ylabel('f(u)')
title('Two bars when increasing force / force-displacement curve')