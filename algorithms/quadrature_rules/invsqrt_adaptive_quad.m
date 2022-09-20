function [c, z, ell] = invsqrt_adaptive_quad(SV,SAV,Sv,tol)


ell1 = 32;
ell2 = 45;

quad_err = inf;
first = true;
while quad_err > tol && ell2 < 1000
    c1 = pi/ell1*ones(1,ell1);
    c2 = pi/ell2*ones(1,ell2);
    z1 = cos((2*(1:ell1)-1)/(2*ell1) * pi);
    z2 = cos((2*(1:ell2)-1)/(2*ell2) * pi);

    if first
        h1 = 0;
        for j = 1:ell1
            h1 = h1 + c1(j)*((SV'*(-(1-z1(j))*SV - (1+z1(j))*SAV))\(SV'*Sv));
        end
        h1 = -2/pi*h1;
        first = false;
    end

    h2 = 0;
    for j = 1:ell2
        h2 = h2 + c2(j)*((SV'*(-(1-z2(j))*SV - (1+z2(j))*SAV))\(SV'*Sv));
    end
    h2 = -2/pi*h2;

    quad_err = norm(h1-h2);
    if quad_err <= tol || floor(sqrt(2)*ell2) >= 1000
        ell = ell2;
        c = c2;
        z = z2;
        disp(['Number of quadrature nodes: ', num2str(ell)])
        break
    else
        ell1 = ell2;
        ell2 = floor(sqrt(2)*ell2);
        h1 = h2;
    end
end