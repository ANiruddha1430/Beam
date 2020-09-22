function un = bc_update(ureduce,BC_data);
m=size(BC_data);
n=size(ureduce);
x=m(1);
y=n(1);
un=zeros((x+y),1);
A=[1:m+n]';
for i=1:m
    if (BC_data(i,2)==1)
        B(i,1)=2*BC_data(i,1)-1;
    else
        B(i,1)=2*BC_data(i,1);
    end
end
c=setdiff(A,B);
for i=1:n
    un(c(i))=ureduce(i);
end
for i=1:m
    un(B(i))=BC_data(i,3);
end
