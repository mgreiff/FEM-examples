load('dataE6.4_01.mat')

disp('plan4gif')
compfe = plan4gif( ec, t, ed , es )';
disp(fe./compfe)

disp('plan4gie')
compKe = plan4gie( ec , t , D , ed , es );
disp(Ke./compKe)

disp('plan4gis')
[compee, compef] = plan4gis( ec , ed );
for ii = 1:4
    disp(ee{ii}./compee{ii})
    disp(ef{ii}./compef{ii})
end
