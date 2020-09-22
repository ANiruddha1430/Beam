% Material properties and other inputs
%--------------------------------------
q0 = 500;  
L1 = 0.1/3; L2 = .1/3; L3 = .1/3; L4 = 0.2/3; L5 = .2/3; L6 = .2/3;
I1 = 5e-9; I2 = 2e-8; E0 = 200e9;
P0 = 1000; M0 = 5000;

nele  = 6;                       % No. of Elements

% Gauss Points and weights vector          || C H A N G E ||
% --------------------------------
ngauss = 3;                              % No. of Gauss points
xivec = [-0.774597, 0, 0.774597];       % Gauss points
wvec = [5/9, 8/9, 5/9];                 % Weights


 

% Co-ordinates for Nodes              || C H A N G E ||
% -------------------------
coord = [1,   0.0;           % First Column is Node numbers
         2,   L1;           % Second Column is Co-ordinate
         3,   L1+L2;
         4,   L1+L2+L3;
         5,   L1+L2+L3+L4;
         6,   L1+L2+L3+L4+L5;
         7,   L1+L2+L3+L4+L5+L6];


connect = [1,  1,  2;      % First Column is element number
           2,  2,  3;     % Second & Third Column are Nodes (in sequence)  
           3,  3,  4;      %       For that element.
           4,  4,  5;
           5,  5,  6;
           6,  6,  7];
           
       
% Material Properties and Area moment of Inertia
% ----------------------------------------------
E = E0*ones(nele,1);
Ie = [I1; I1; I1; I2; I2; I2];

% Boundary Condition Suppressed
% -----------------------------
BC_data = [1, 1, 0;        % First Column is Node number
           1, 2, 0;       % Second Column is the prescribed D.O.F
           7, 1, 0];      % Third Column is value of the prescribed D.O.F
           
       
% Point Load and Point Moment Data       
P_load = [4, P0];       % First Column is Node number
                        % Second Column is Point load value
       
P_moment = [7, M0];     % First Column is Node number
                        % Second Column is Point moment value
       
% Distributed Load data
q_load = [1, q0,   0,  0;   % First Column is element number
          1, q0,   0,  0];  % For quadratic load  q = a + bx + cx^2
                            % 2nd to 4th Columns are a, b, c 
