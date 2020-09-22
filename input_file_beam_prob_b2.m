% Material properties and other inputs
%--------------------------------------
E0 = 70e9; Ie0 = 8.1e-6; q0 = 1000;  
L1 = 0.8; L2 = 0.6; F0 = 1000; q1 = 444.44; q2 = 111.111;

nele  = 6;                       % No. of Elements

% Gauss Points and weights vector          || C H A N G E ||
% --------------------------------
ngauss = 3;                              % No. of Gauss points
xivec = [-0.774597, 0, 0.774597];       % Gauss points
wvec = [5/9, 8/9, 5/9];                 % Weights


 

% Co-ordinates for Nodes              || C H A N G E ||
% -------------------------
coord = [1,   0.0;           % First Column is Node numbers
         2,   L1/3;           % Second Column is Co-ordinate
         3,   2*L1/3;
         4,   L1;
         5,   L1+L2/3;
         6,   L1+2*L2/3;
         7,   L1+L2];

connect = [1,  1,  2;      % First Column is element number
           2,  2,  3;     % Second & Third Column are Nodes (in sequence)  
           3,  3,  4;                %       For that element.
           4,  4,  5;
           5,  5,  6;
           6,  6,  7];
       
% START CHANGE %         

% Material Properties and Area moment of Inertia
% ----------------------------------------------
E = E0*ones(nele,1);
Ie = Ie0*ones(nele,1);


% Boundary Condition Suppressed
% -----------------------------
BC_data = [1, 1, 0;        % First Column is Node number
           1, 2, 0];       % Second Column is the prescribed D.O.F
                           % Third Column is value of the prescribed D.O.F
           
       
% Point Load and Point Moment Data       
P_load = [7, F0];       % First Column is Node number
                           % Second Column is Point load value
       
P_moment = [];     % First Column is Node number
                           % Second Column is Point moment value
       
% Distributed Load data
q_load = [1, q0,           -2*q0/(L1),  q0/(L1)^2;  % First Column is element number
          2, q0,           -2*q0/(L1),  q0/(L1)^2;  % For quadratic load  q = a + bx + cx^2
          3, q0,           -2*q0/(L1),  q0/(L1)^2]; % 2nd to 4th Columns are a, b, c 

% Analytical Solution
% End Slope:
L = L1 + L2;
M1 = F0*L + q0*L1^2/12; R1 = F0 + q0*L1/3;

AreaM = -M1*L + R1*L^2/2 - 11*q0*L1^3/120 -q0*L1^2*L2/4 - q0*L1*L2^2/6;
theta = AreaM/E(1)/Ie(1);

% End Displacement:
MomentM = -M1*L^2/2 + R1*L^3/6 + (q0*L1^3/2)*(-11*L2/60 - L1/20)-q0*L1^2*L2^2/8 - q0*L1*L2^3/18;
delta = MomentM/E(1)/Ie(1);                                         