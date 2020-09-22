% ===================
% Reading Input File:
% ===================
for i=1:2
    
A=sprintf('input_file_beam_prob_c%d',i);
eval(A)

% Global Stiffnes Matrix and Global load vector
% ---------------------------------------------
% Function "stiff_load" calculate Global Stiffness Matrix and 
%     Global load vector due to distributed load
% ----------
% I N P U T
% = = = = = 
% nele   = No. of elements 
% ngauss = No. of gauss points for integration
% coord  = Nodal coordinates    % First Column is Node numbers
%                                 Second Column is Co-ordinate
% connect = Nodal Connectivities    % First Column is element number    
                                    % Second & Third Column are Nodes (in sequence)  
                                    % For that element. 
% xivec = Gauss points
% wvec  = weights
% E = Young's Modulus of the element
% Ie = Area Moment of inertia of the element
% q0 = Maximum distributed load (Triangular) magnitude
% L  = Length
%
% ----------
% O U T P U T
% = = = = =
% K = Global stiffness matrix
% F = Global load vector
%

[K,F,k] = stiff_load(nele,ngauss,coord,connect,xivec,wvec,E,Ie,q_load);

% Point load and Point moment
% ---------------------------
% This function "point_ld_mom" update Global load vector after incorporating point load 
% and point moment data
%
% ----------
% I N P U T
% = = = = = 
% F      = Global load vector before implementing point load and point moment
% P_load = Point load data  % First Column is Node number
                            % Second Column is Point load value
% P_moment = Point moment data       % First Column is Node number
                                     % Second Column is Point moment value
% -------------
% O U T P U T
% = = = = =====
% F = Global load vector after implementing point load and point moment

F = point_ld_mom(F,P_load,P_moment);


% ===========================
% Imposition of B.C.
K_glob = K;
F_glob = F;

% This function "impose_bc" update Global stiffness matrix and Global load vector 
% after incorporating boundary condition data
%
% ----------
% I N P U T
% = = = = = 
% K       = Global stiffness matrix before implementing Boundary condition data
% F       = Global load vector before implementing Boundary condition data
% BC_data = Boundary condition data        % First Column is Node number
                                           % Second Column is the prescribed D.O.F
                                           % Third Column is value of the prescribed D.O.F 
%
% -------------
% O U T P U T
% = = = = =====
% K       = Global stiffness matrix after implementing Boundary condition data
% F       = Global load vector after implementing Boundary condition data

[K,F] = impose_bc(nele,K,F,BC_data);


% Finding Solution
ureduce = inv(K)*F;

% Full Solution vector (Free + Prescribed D. O. F.)
% -------------------------------------------------
% This function "bc_update" update solution vector with values of prescribed DOFs
% 
% I N P U T
% =========
% ureduce = Solution vector just after inversion
%           It contains only free DOFs
% BC_data = Boundary condition data        % First Column is Node number
                                           % Second Column is the prescribed D.O.F
                                           % Third Column is value of the prescribed D.O.F 
% -------------
% O U T P U T
% = = = = =====
% un = Full solution vectors with Free and Prescribed DOF values
un = bc_update(ureduce,BC_data)

% Finding Reaction Force
Freac = K_glob*un;

% Post Processing: FEM displacement
xi = [-1:0.2:1]';          % Distribution of data points

% This function "postprocessing" calculate variable u at diffent distributed points across
% element from nodal values of u
% This function calculate variable u at diffent distributed points across
% element from nodal values of u
% ----------
% I N P U T
% = = = = = 
% nele   = No. of elements
% coord  = Nodal coordinates    % First Column is Node numbers
%                                 Second Column is Co-ordinate
% connect = Nodal Connectivities    % First Column is element number    
                                    % Second & Third Column are Nodes (in sequence)  
                                    % For that element. 
% xi = Points distributed for an element in master domain
% un = Nodal values of u
%
% ----------
% O U T P U T
% = = = = =
% xnume = x coordinates of the distributed points
% unume = values of u at distributed points

[xnume, unume] = postprocessing(nele,coord,connect,un,xi);
rotation = cal_theta(nele,coord,connect,un,xi);
figure(1)

if i==1
    defl=figure(1);
    hold on
    plot(xnume,unume,'-r*')
    rota=figure(2);
    hold on
    plot(xnume,rotation,'-r*')
else
    defl=figure(1);
    plot(xnume,unume,'-bo')
    hold off
    rota=figure(2);
    plot(xnume,rotation,'-bo')
    hold off
end
%writing a file

if i==1
    fid=fopen('For 2 element.txt ','w');
    fprintf(fid,'\n');
    fprintf(fid,'Elemental stiffness matrix\n');
    fprintf(fid,'**************************\n');
    size_k=size(k);
    for l=1:nele
        fprintf(fid,'The Element %d stiffness matrix is\n',l);
        for m=1:4
            for n=1:4
                fprintf(fid,'%12.4e\t',k(m,n,l));
            end
                 fprintf(fid,'\n');
        end
             fprintf(fid,'\n');
    end
    
    fprintf(fid,'Global stiffness matrix\n');
    fprintf(fid,'****************************\n');
    for m=1:2*(nele+1)
        for n=1:2*(nele+1)
            fprintf(fid,'%12.4e\t',K_glob(m,n));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'Global load vector\n');
    fprintf(fid,'**************************\n');
    for m=1:2*(nele+1)
        fprintf(fid,'%12.4e\n',F_glob(m,1));
    end
    fprintf(fid,'\n');
    fprintf(fid,'Global stiffness matrix after imposing boundary conditions\n');
    fprintf(fid,'***********************************************************\n');
    size_K=size(K)
    for m=1:size_K(1,1)
        for n=1:size_K(1,1)
            fprintf(fid,'%12.4e\t',K(m,n));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'Global load matrix after imposing boundary conditions\n');
    fprintf(fid,'***********************************************************\n');
    size_F=size(F)
    for m=1:size_F(1,1)
            fprintf(fid,'%12.4e\n',F(m,1));
    end
    fprintf(fid,'\n');
    fprintf(fid,'u_reduce matrix is\n');
    fprintf(fid,'**********************');
    size_ureduce=size(ureduce);
    for m=1:size_ureduce(1,1)
            fprintf(fid,'%12.4e\n',ureduce(m,1));
    end
    fprintf(fid,'\n');
    fprintf(fid,' un\n');
    fprintf(fid,'**********************');
    size_un=size(un);
    for m=1:size_un(1,1)
            fprintf(fid,'%12.4e\n',un(m,1));
    end
    fprintf(fid,'\n');
    fprintf(fid,' Reaction Force\n');
    fprintf(fid,'**********************');
    size_Freac=size(Freac);
    for m=1:size_Freac(1,1)
            fprintf(fid,'%12.4e\n',Freac(m,1));
    end
    fprintf(fid,'\n');
end
if i==2
    fid=fopen('For 6 element.txt ','w');
    fprintf(fid,'\n');
    fprintf(fid,'Elemental stiffness matrix\n');
    fprintf(fid,'**************************\n');
    size_k=size(k);
    for l=1:nele
        fprintf(fid,'The Element %d stiffness matrix is\n',l);
        for m=1:4
            for n=1:4
                fprintf(fid,'%12.4e\t',k(m,n,l));
            end
                 fprintf(fid,'\n');
        end
             fprintf(fid,'\n');
    end
    
    fprintf(fid,'Global stiffness matrix\n');
    fprintf(fid,'****************************\n');
    for m=1:2*(nele+1)
        for n=1:2*(nele+1)
            fprintf(fid,'%12.4e\t',K_glob(m,n));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'Global load vector\n');
    fprintf(fid,'**************************\n');
    for m=1:2*(nele+1)
        fprintf(fid,'%12.4e\n',F_glob(m,1));
    end
    fprintf(fid,'\n');
    fprintf(fid,'Global stiffness matrix after imposing boundary conditions\n');
    fprintf(fid,'***********************************************************\n');
    size_K=size(K)
    for m=1:size_K(1,1)
        for n=1:size_K(1,1)
            fprintf(fid,'%12.4e\t',K(m,n));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'Global load matrix after imposing boundary conditions\n');
    fprintf(fid,'***********************************************************\n');
    size_F=size(F)
    for m=1:size_F(1,1)
            fprintf(fid,'%12.4e\n',F(m,1));
    end
    fprintf(fid,'\n');
    fprintf(fid,'u_reduce matrix is\n');
    fprintf(fid,'**********************');
    size_ureduce=size(ureduce);
    for m=1:size_ureduce(1,1)
            fprintf(fid,'%12.4e\n',ureduce(m,1));
    end
    fprintf(fid,'\n');
    fprintf(fid,' un\n');
    fprintf(fid,'**********************');
    size_un=size(un);
    for m=1:size_un(1,1)
            fprintf(fid,'%12.4e\n',un(m,1));
    end
    fprintf(fid,'\n');
    fprintf(fid,' Reaction Force\n');
    fprintf(fid,'**********************');
    size_Freac=size(Freac);
    for m=1:size_Freac(1,1)
            fprintf(fid,'%12.4e\n',Freac(m,1));
    end
    fprintf(fid,'\n');
end

end
defl=figure(1);
legend('2-elements','6-element');
xlabel('xnume');
ylabel('unume');
saveas(defl,'deflection','png');


rota=figure(2);
legend('2-elements','6-element');
xlabel('xnume');
ylabel('theta');
saveas(rota,'rotation','png');

