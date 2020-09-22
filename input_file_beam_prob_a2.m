% Material properties and other inputs
%--------------------------------------
E0 = 90e6; Ie0 = 8.24e-4; q0 = 300; q1=250; q2 = 200; 
L1 = 0.1; L2 = 0.1; L3 = 0.12; L4 = .12; L5 = .08; L6 = .08; F0 = 5000;

nele  = 6;                       % No. of Elements

% Gauss Points and weights vector          || C H A N G E ||
% --------------------------------
ngauss = 3;                              % No. of Gauss points
xivec = [-0.774597, 0, 0.774597];       % Gauss points
wvec = [5/9, 8/9, 5/9];                 % Weights


 

% Co-ordinates for Nodes              || C H A N G E ||
% -------------------------
coord = [1,   0.0;           % First Column is Node numbers
         2,   L1;            % Second Column is Co-ordinate
         3,   L1+L2;
         4,   L1+L2+L3;
	     5,   L1+L2+L3+L4;
	     6,   L1+L2+L3+L4+L5;
	     7,   L1+L2+L3+L4+L5+L6];

connect = [1,  1,  2;      % First Column is element number
           2,  2,  3;      % Second & Third Column are Nodes (in sequence)  
           3,  3,  4;     %       For that element.
	       4,  4,  5;
	       5,  5,  6;
           6,  6,  7];         

% START CHANGE %         

% Boundary Condition Suppressed
% -----------------------------
BC_data = [1, 1, 0;        % First Column is Node number
           1, 2, 0;        % Second Column is the prescribed D.O.F
           5, 1, 0];       % Third Column is value of the prescribed D.O.F
           
% Material Properties and Area moment of Inertia
% ----------------------------------------------
E = E0*ones(nele,1);
Ie = Ie0*ones(nele,1);
       
       
% Point Load and Point Moment Data       
P_load = [7, F0];       % First Column is Node number
                           % Second Column is Point load value
       
P_moment = [];     % First Column is Node number
                           % Second Column is Point moment value
       
% Distributed Load data
q_load = [1, q0,            -(q0-q1)/L1,  0;   % First Column is element number
          2, q1,            -(q1-q2)/L2,  0;   % For quadratic load  q = a + bx + cx^2                       
          3, q2,                      0,  0;   % 2nd to 4th Columns are a, b, c 
          4, q2,                      0,  0];
