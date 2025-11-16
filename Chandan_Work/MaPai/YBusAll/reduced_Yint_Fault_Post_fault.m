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
GF=real(YintF);
BF=imag(YintF);
disp(YintF);
disp(YintPF);

