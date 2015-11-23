function [ g ] = g_bar( u )
run data.m
g = 2.*EA.*(a/l0).^3*(u/a + 3/2 .* (u/a).^2 + (1/2).*(u/a).^3) + k*u;
end

