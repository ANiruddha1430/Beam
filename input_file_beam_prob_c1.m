% Material properties and other inputs
%--------------------------------------
q0 = 500;  
L1 = 0.1; L2 = 0.2; 
I1 = 5e-9; I2 = 2e-8; E0 = 200e9;
P0 = 1000; M0 = 5000;

nele  = 2;                       % No. of Elements

% Gauss Points and weights vector          || C H A N G E ||
% --------------------------------
ngauss = 3;                              % No. of Gauss points
xivec = [-0.774597, 0, 0.774597];       % Gauss points
wvec = [5/9, 8/9, 5/9];                 % Weights


 

% Co-ordinates for Nodes              || C H A N G E ||
% -------------------------
coord = [1,   0.0;           % First Column is Node numbers
         2,   L1;           % Second Column is Co-ordinate
         3,   L1+L2];


connect = [1,  1,  2;      % First Column is element number
           2,  2,  3];     % Second & Third Column are Nodes (in sequence)  
                           %       For that element.
           
       
% Material Properties and Area moment of Inertia
% ----------------------------------------------
E = E0*ones(nele,1);
Ie = [I1; I2];

% Boundary Condition Suppressed
% -----------------------------
BC_data = [1, 1, 0;        % First Column is Node number
           1, 2, 0;       % Second Column is the prescribed D.O.F
           3, 1, 0];      % Third Column is value of the prescribed D.O.F
           
       
% Point Load and Point Moment Data       
P_load = [2, P0];       % First Column is Node number
                        % Second Column is Point load value
       
P_moment = [3, M0];     % First Column is Node number
                        % Second Column is Point moment value
       
% Distributed Load data
q_load = [1, q0,   0,  0];  % First Column is element number
                            % For quadratic load  q = a + bx + cx^2
                            % 2nd to 4th Columns are a, b, c 
