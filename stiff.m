function k=stiff(E,Ie,L)
k=E*Ie/(L^3)*[12, 6*L, -12, 6*L;
              6*L, 4*L^2, -6*L, 2*L^2;
              -12, -6*L, 12, -6*L;
              6*L, 2*L^2, -6*L, 4*L^2];