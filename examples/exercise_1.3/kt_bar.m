function [ K ] = kt_bar( u, a, b, l0, EA)
    k = 1.5 * (EA / l0).*(a/l0).^2;
    K = (2.*EA/l0).*(a/l0).^2.*(1 + 3*(u/a) + (3/2).*(u/a).^2) + k;
end

