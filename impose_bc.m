function [K,F] = impose_bc(nele,K,F,BC_data)
m=size(BC_data);
n=m(1);
A=[1:2*(nele+1)];
B=zeros(n,1);
for i=1:n
    if (BC_data(i,2)==1)
        B(i,1)=2*BC_data(i,1)-1;
    else
        B(i,1)=2*BC_data(i,1);
    end
end
for i=1:2*(nele+1)
    for j=1:n
        if BC_data(j,2)==1
            F(i)=F(i)-K(i,2*BC_data(j,1)-1)*BC_data(j,3);
        else
            F(i)=F(i)-K(i,2*BC_data(j,1))*BC_data(j,3);
        end
    end
end
        
c=setdiff(A',B);
K=K(c,c);
F=F(c,1);



