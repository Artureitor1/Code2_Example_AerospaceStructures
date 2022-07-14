%-------------------------------------------------------------------------%
% ASSIGNMENT 02 (B)
%-------------------------------------------------------------------------%
% Date:
% Author/s: Pol Padilla; Arturo Baltanas
%

clear;
close all;

%% INPUT DATA

% Geometric data
H = 0.9;
W = 0.85;
B = 3.2;
d = 0;
D1 =18e-3 ; 
d1 = 7.5e-3;
D2 = 3e-3;

% Mass
M = 150;

% Other
g = 9.81;

%% PREPROCESS

% Nodal coordinates matrix 
%  x(a,j) = coordinate of node a in the dimension j
x = [%     X      Y      Z
         2*W,  -W/2,     0; % (1)
         2*W,   W/2,     0; % (2)
         2*W,     0,     H; % (3)
           0,     0,     H; % (4)
           0,    -B,     H; % (5)
           0,     B,     H; % (6)
           W,     0,     H; % (7)
];

% Nodal connectivities  
%  Tnod(e,a) = global nodal number associated to node a of element e
Tnod = [%     a      b
1 2
1 3
2 3
3 7
4 7
4 5
3 5
3 6
4 6
5 7
6 7
1 4
1 5
1 7
2 4
2 6
2 7
];


% Fix nodes matrix creation
%  fixNod(k,1) = node at which some DOF is prescribed
%  fixNod(k,2) = DOF prescribed (local, 1-2)
%  fixNod(k,3) = prescribed displacement in the corresponding DOF (0 for fixed)
fixNod = [% Node        DOF  Magnitude
          % Write the data here...
          4 1 0
          4 2 0
          4 3 0
          3 2 0
          3 3 0
          6 3 0
];

% Material properties matrix
%  mat(m,1) = Young modulus of material m
%  mat(m,2) = Section area of material m
%  mat(m,3) = Density of material m
%  --more columns can be added for additional material properties--
mat = [% Young M.        Section A.         Density      Inertia
       75000e6,    pi*((D1/2)^2-(d1/2)^2),   3350,  (pi/4)*((D1/2)^4-(d1/2)^4);    % Material (1)
          147000e6,      pi*(D2/2)^2,        950,       (pi/4)*(D2/2)^4;          % Material (2)
];

% Material connectivities
%  Tmat(e) = Row in mat corresponding to the material associated to element e 
Tmat = [% Mat. index
 1;1;1;1;1;1;1;1;1;1;1
 2;2;2;2;2;2
];

%% SOLVER

% Dimensions
n_d = size(x,2);              % Number of dimensions
n_i = n_d;                    % Number of DOFs for each node
n = size(x,1);                % Total number of nodes
n_dof = n_i*n;                % Total number of degrees of freedom
n_el = size(Tnod,1);          % Total number of elements
n_nod = size(Tnod,2);         % Number of nodes for each element
n_el_dof = n_i*n_nod;         % Number of DOFs for each element 


% Computation of the DOFs connectivities
Td = connectDOFs(n_el,n_nod,n_i,Tnod);

% Computation of element stiffness matrices
Kel = computeKelBar(n_d,n_el,n_el_dof,x,Tnod,mat,Tmat);

% Global matrix assembly
KG = assemblyKG(n_el,n_el_dof,n_dof,Td,Kel);

% Global force vector assembly
Fext = computeF(n_i,n_dof,Tmat,mat,x,Tnod,n_el,g,M);

% Apply conditions 
[vL,vR,uR] = applyCond(n_i,n_dof,fixNod);

% System resolution
[u,R] = solveSys(vL,vR,uR,KG,Fext);

% Compute strain and stresses
[eps,sig] = computeStrainStressBar(n_d,n_el,u,Td,x,Tnod,mat,Tmat);

% Check of buckling failure
[FB] = bucklingFailure(mat,Tmat,x,Tnod,n_el, sig, D1, d1, D2);


%% POSTPROCESS

% Plot deformed structure with stress of each bar
scale = 10; % Adjust this parameter for properly visualizing the deformation
plotBarStress3D(x,Tnod,u,sig,scale);