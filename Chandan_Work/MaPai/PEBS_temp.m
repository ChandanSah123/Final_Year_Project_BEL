Y = [ ...
    -1i*17.361     0             0            1i*17.361      0              0              0              0              0;
     0            -1i*16         0            0              0              0           1i*16                0              0;
     0             0            -1i*17.065    0              0              0               0              0              1i*17.065;
     1i*17.361     0             0            3.307 - 1i*39.309  -1.365 + 1i*11.604  -1.942 + 1i*10.511  0  0  0;
     0             0             0           -1.365 + 1i*11.604   2.553 - 1i*17.338     0           -1.188 + 1i*5.975   0  0;
     0             0             0           -1.942 + 1i*10.511    0    3.224 - 1i*15.841   0   0   -1.282 + 1i*5.588;
     0             1i*16         0            0          -1.188 + 1i*5.975   0   2.805 - 1i*35.446  -1.617 + 1i*13.698   0;
     0             0             0            0          0          0   -1.617 + 1i*13.698   2.772 - 1i*23.303  -1.155 + 1i*9.784;
     0             0             1i*17.065    0          0   -1.282 + 1i*5.588   0   -1.155 + 1i*9.784   2.437 - 1i*32.154
];

% Create readable MATLAB table
Y_table = array2table(Y, ...
    'VariableNames', compose('Bus%d', 1:9), ...
    'RowNames', compose('Bus%d', 1:9));

% Export to Excel file
writetable(Y_table, 'Ybus_matrix.xlsx', 'WriteRowNames', true);

%Optional: display to check
disp(Y_table)
save('Y.mat','Y');

%%

% Y previously entered Ybus matrix

% Define ybar (3x3)
ybar = diag([-1i*16.45, -1i*8.35, -1i*5.52]);
load_ybar=complex(zeros(9,9));
load_ybar(5,5)=(1.26-0.504i);
load_ybar(6,6)=(0.877-0.292i);
load_ybar(8,8)=(0.969-0.339i);

% --- Construct the augmented matrix ---
Y_aug = zeros(12,12);

% Map the existing Y into Y_aug (rows 4:12, cols 4:12)
Y(5,5)=2.553-17.338i;
Y(6,6)=3.224-15.841i;
Y(8,8)=2.772-23.303i;
YN = Y;

% Add ybar and -ybar blocks
Y_aug(1:3, 1:3) = ybar;        % top-left block
Y_aug(1:3, 4:6) = -ybar;       % top-middle block (coupling)
Y_aug(4:6, 1:3) = -ybar;       % middle-left block (transpose coupling)

% Keep others as zero per the paper figure
YNadd=zeros(9,9);
YNadd(1:3,1:3)=ybar;
YN1=YN+YNadd;
YN2=YN1+load_ybar;
Y_aug(4:12, 4:12) = YN2;

%%
% ybar (machine terminal admittances) you provided earlier:
ybar = diag([-1i*16.45, -1i*8.35, -1i*5.52]);
% load admittance additions (as you had)
load_ybar = complex(zeros(9,9));
load_ybar(5,5) = (1.26 - 0.504i);
load_ybar(6,6) = (0.877 - 0.292i);
load_ybar(8,8) = (0.968 - 0.339i);

% transformer branch Y for 5-7 (used when removing branch)
Y78 = 1.187 - 5.975i;
B78=0.153i;

% Pre-fault machine terminal voltages (angles) from Table 7.1 / Example 7.6
Vt1 = 1.04*exp(1j*deg2rad(0.0));
Vt2 = 1.025*exp(1j*deg2rad(9.3));
Vt3 = 1.02*exp(1j*deg2rad(4.7));
Vt = [Vt1, Vt2, Vt3];

%%%%% ------------------ BUILD AUGMENTED MATRIX & Yint ------------------ %%%%%
% The augmented Y (12x12) format: top-left machine terminal nodes (3), rest 9 network nodes
YN = Y;             % base network Y (9x9)
% ybar on the network diagonal for terminals
YNadd = zeros(9,9);
YNadd(1:3,1:3) = ybar;
YN1 = YN + YNadd;
YN2 = YN1 + load_ybar;   % final network block used in your prior code

% Build augmented Y (12x12)
Y_aug = complex(zeros(12,12));
Y_aug(1:3,1:3) = ybar;
Y_aug(1:3,4:6) = -ybar;
Y_aug(4:6,1:3) = -ybar;
Y_aug(4:12,4:12) = YN2;

% Reduce to internal Yint (post-fault etc uses same reduction pattern)
YA = ybar;
YB = complex(zeros(3,9));
YB(:,1:3) = -ybar;
YC = complex(zeros(9,3));
YC(1:3,:) = -ybar;
YD = YN2;
Yint_pre = YA - YB * (YD \ YC);   % internal (pre-fault) Y seen by internal EMFs

% Save these for debugging/inspection
Yint_pre_full = Yint_pre;

%%%%% ------------------ CREATE FAULT-ON and POST-FAULT Yint ------------------ %%%%%
% Fault at bus 7 (bus indexing: follow your Y bus indexes so 7 means the same)
% Create YN2F (remove bus 7 from network) - for the fault-on reduction we need the network with bus 7
% ------------------ CORRECT 3-PHASE FAULT MODEL (BUS 7) ------------------
YN2F = YN2;
YN2F(7,:) = [];
YN2F(:,7) = [];

% Build YintF using same block partitioning but adapted sizes:
% Note: for fault-on Yint, the network block dimension shrinks to 8 x 8 (machine block still 3)
YA_F = ybar;
YB_F = complex(zeros(3,8));
YB_F(:,1:3) = -ybar;
YC_F = complex(zeros(8,3));
YC_F(1:3,:) = -ybar;
YD_F = YN2F;
YintF = YA_F - YB_F * (YD_F \ YC_F);


% Post-fault: line 5-7 opened => modify YN2 accordingly
YN2PF = YN2;
% remove branch admittance Y78 from off-diagonals and diagonals
YN2PF(5,7) = 0;
YN2PF(7,5) = 0;
YN2PF(5,5) = YN2PF(5,5) - (Y78+B78);
YN2PF(7,7) = YN2PF(7,7) - (Y78+B78);

% Now form YintPF (post-fault reduced internal)
YAPF = ybar;
YBPF = complex(zeros(3,9));
YBPF(:,1:3) = -ybar;
YCPF = complex(zeros(9,3));
YCPF(1:3,:) = -ybar;
YDPF = YN2PF;
YintPF = YAPF - YBPF * (YDPF \ YCPF);
%disp(YintPF);
%disp(YintF );

G = real(YintPF);
B = imag(YintPF);

%%
close all; clc;
% Inertia constants (book example)
H= [23.64, 6.4, 3.01]; 
M = 2*H;
Pm = [0.716, 1.63, 0.85]; 
%internal emf and rotor angle
E1 = 1.054*exp(1j*deg2rad(2.2670));
E2 = 1.050*exp(1j*deg2rad(19.750));
E3 = 1.017*exp(1j*deg2rad(13.20));
E = [abs(E1), abs(E2), abs(E3)];
angleE=[angle(E1),angle(E2),angle(E3)]; 
%intializing the rotor angles
delpre=angleE;
del0=(delpre(1)*M(1)+delpre(2)*M(2)+delpre(3)*M(3))/sum(M);
del10=delpre(1);
del20=delpre(2);
del30=delpre(3);
theta10=del10-del0;
theta20=del20-del0;
theta30=del30-del0;
theta0=[theta10;theta20;theta30];
w10=0;
w20=0;
w30=0;
%thetas1=-0.1782;
%thetas2=0.5309;
%thetas3=0.2711;


% forming Cij and Dij
C=zeros(3,3);
D=zeros(3,3);
for i = 1:3
     
        for j = 1:3
          C(i,j)=E(i)*E(j)*B(i,j);
          D(i,j)=E(i)*E(j)*G(i,j);
        end
     
end
%%

% Solve fi(theta) = 0 using Newton-Raphson

%=== Symbolic variables ===%
syms theta1 theta2 theta3

%=== Define equations ===%
f1 = Pm(1)-D(1,1)-( C(1,2)*sin(theta1-theta2) + D(1,2)*cos(theta1-theta2) ...
                  + C(1,3)*sin(theta1-theta3) + D(1,3)*cos(theta1-theta3) );

f2 = Pm(2)-D(2,2)-( C(2,1)*sin(theta2-theta1) + D(2,1)*cos(theta2-theta1) ...
                  + C(2,3)*sin(theta2-theta3) + D(2,3)*cos(theta2-theta3) );

f3 = Pm(3)-D(3,3)-( C(3,1)*sin(theta3-theta1) + D(3,1)*cos(theta3-theta1) ...
                  + C(3,2)*sin(theta3-theta2) + D(3,2)*cos(theta3-theta2) );

F = [f1; f2; f3];

%=== Jacobian ===%
J = jacobian(F,[theta1,theta2,theta3]);

%=== Convert symbolic to functions ===%
F_fun = matlabFunction(F,'Vars',{theta1,theta2,theta3});
J_fun = matlabFunction(J,'Vars',{theta1,theta2,theta3});

%=== Initial guess ===%
theta = theta0;
maxIter = 50;
tol = 1e-12;
damping = 0.7;

%=== Newton-Raphson iterations (COI enforced) ===%
for k = 1:maxIter
  

    F_val = F_fun(theta(1),theta(2),theta(3));
    J_val = J_fun(theta(1),theta(2),theta(3));

    % COI constraint: add row to eliminate singularity
    J_mod = J_val;
    F_mod = F_val;
    J_mod(1,:) = [M(1), M(2), M(3)];
    F_mod(1)   = 0;

    delta = -J_mod \ F_mod;
    theta = theta + damping*delta;

    if norm(delta) < tol 
        break;
    end
end

disp('Solved theta angles (radians):')
disp(theta);
disp(YintF);
disp(YintPF);
%% there is certain problem either in code or in parameter value which is making it difficult to converge the sep to the exact book value
%the mismatch is about (-0.1782  0.5309  0.2711)-(-0.2023 0.5982 0.3154)=
%(-0.0241 -0.067 -0.0443)
%calculation of fault on trajectory

