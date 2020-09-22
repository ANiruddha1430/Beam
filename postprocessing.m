function [xnume, unume] = postprocessing(nele,coord,connect,un,xi)
m=size(xi);
xn=[];
for i=1:nele
    x1=coord(connect(i,2),2);
    x2=coord(connect(i,3),2);
    le=x2-x1;
    u(1:4,1)=un(2*i-1:2*i+2,1);
    for j=1:m(1)
        N1 = (2-3*xi(j,1)+xi(j,1)^3)/4;
    N2 = (1-xi(j,1) -xi(j,1)^2 +xi(j,1)^3)/4;
    N3 = (2 + 3*xi(j,1) -xi(j,1)^3)/4;
    N4 = (-1 -xi(j,1) + xi(j,1)^2 + xi(j,1)^3)/4;
    Nu = [N1, le*N2/2, N3, le*N4/2];
    
        xn(j,i)=((1-(xi(j,1)))/2)*x1+((1+xi(j,1))/2)*x2;
        u1(j,i)=Nu*u;
    end
end
xnume=reshape(xn,1,[]);
xnume=xnume';
unume=reshape(u1,1,[]);
unume=unume';




    