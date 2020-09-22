function F = point_ld_mom(F,P_load,P_moment)

fnd=size(P_load);
nd=fnd(1);
if nd>=0
    for i=1:nd
        F(2*(P_load(i,1))-1,1)=F(2*(P_load(i,1))-1,1)+P_load(i,2);
    end
end
F;

mnd=size(P_moment);
md=mnd(1);

if md>0
    for i=1:md
       F(2*(P_moment(i,1)),1)=F(2*(P_moment(i,1)),1)+P_moment(i,2);
    end
end
