function [K,F,k]=stiff_load(nele,ngauss,coord,connect,xivec,wvec,E,Ie,q_load)

K=zeros(2*(nele+1),2*(nele+1));

for i=1:nele
    L=coord(connect(i,3),2)-coord(connect(i,2),2);
    k(1:4,1:4,i)=stiff(E(i,1),Ie(i,1),L);
    K(2*i-1:2*i+2,2*i-1:2*i+2)=K(2*i-1:2*i+2,2*i-1:2*i+2)+k(1:4,1:4,i); 
end

m=size(q_load);
n=m(1);
F=zeros(2*(nele+2),1);

for i=1:n
    x1=coord(connect(i,2),2);
    x2=coord(connect(i,3),2);
    le=x2-x1;
    f=zeros(4,1);
    
    for j=1:ngauss
        xi=xivec(j);
        w=wvec(j);
        x=(1-xi)*x1/2+(1+xi)*x2/2;
        q=q_load(i,2)+q_load(i,3)*x+q_load(i,4)*x*x;
        f=f+w*ld(le,q,xi);
    end
    
    F((2*i-1):(2*i+2),1)=F((2*i-1):(2*i+2),1)+f(1:4,1);
end

F;

