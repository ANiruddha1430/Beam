function f=ld(le,q,xi)
    N1 = (2-3*xi+xi^3)/4;
    N2 = (1-xi -xi^2 +xi^3)/4;
    N3 = (2 + 3*xi -xi^3)/4;
    N4 = (-1 -xi + xi^2 + xi^3)/4;
    Nu = [N1; le*N2/2; N3; le*N4/2];
    f=Nu*q*le/2;