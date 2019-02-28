function [ g ] = g_bar( u, a, b, l0, EA)
    k = 1.5 * (EA / l0).*(a/l0).^2;
    g = 2.*EA.*(a/l0).^3*(u/a + 3/2 .* (u/a).^2 + (1/2).*(u/a).^3) + k*u;
end

