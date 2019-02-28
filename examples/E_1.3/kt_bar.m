function [ K ] = kt_bar( u )
run data.m
K = (2.*EA/l0).*(a/l0).^2.*(1 + 3*(u/a) + (3/2).*(u/a).^2) + k;
end

