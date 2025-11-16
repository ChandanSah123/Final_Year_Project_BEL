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

thetas1=-0.1782;
thetas2=0.5309;
thetas3=0.2711;
theta_s=[thetas1;thetas2;thetas3];

%% ---------------- Step-2: Fault-on trajectory + energy traces (PEBS) -------------
% Requires existing variables in workspace:
%   YintF (3x3)    - reduced internal admittance during fault
%   YintPF (3x3)   - reduced internal admittance after clearing (post-fault)
%   theta0 (3x1)   - initial rotor angles (COI-shifted)
%   M (3x1)        - inertia weights (same used earlier)
%   E (3x1) or (3x1) - internal EMF magnitudes (column)
%   Pm (3x1)       - mechanical powers (column)
%   theta_s (3x1)  - post-fault SEP (COI frame), V_PE(theta_s) will be reference zero
%
% If any of the above are row vectors, convert to columns:
E = E(:);
Pm = Pm(:);
M  = M(:);

% faulted network matrices
Gf = real(YintF);
Bf = imag(YintF);

% post-fault network for V_PE evaluation
Gpost = real(YintPF);
Bpost = imag(YintPF);

% helper: compute Pe vector for any theta using given G,B,E
function Pe = compute_Pe(theta_vec, Evec, Gmat, Bmat)
    mloc = length(theta_vec);
    Pe = zeros(mloc,1);
    for ii = 1:mloc
        s = 0;
        for jj = 1:mloc
            phi = theta_vec(ii) - theta_vec(jj);
            s = s + Evec(ii)*Evec(jj)*( Gmat(ii,jj)*cos(phi) + Bmat(ii,jj)*sin(phi) );
        end
        Pe(ii) = s;
    end
end

% helper: compute V_PE(theta) by straight-line path integral from theta_s to theta
% V_PE(theta) = - \int_0^1 (Pm - Pe(theta_s + s*(theta-theta_s)))' * (theta - theta_s) ds
function Vpe = compute_VPE_numeric(theta, theta_s, Evec, Gmat, Bmat, Pmvec, nsteps)
    if nargin < 7, nsteps = 200; end
    svec = linspace(0,1,nsteps);
    integrand = zeros(size(svec));
    delta = theta - theta_s;
    for k = 1:length(svec)
        s = svec(k);
        thetas = theta_s + s*delta;
        Pe_s = compute_Pe(thetas, Evec, Gmat, Bmat);
        integrand(k) = (Pmvec - Pe_s).' * delta;   % scalar
    end
    % trapezoid integrate: integral_0^1 integrand(s) ds
    I = trapz(svec, integrand);
    Vpe = - I;   % by formula above; Vpe(theta_s) = 0
end

% -- simulate fault-on swing equations (classical model) --
% state x = [theta1;theta2;theta3; omega1;omega2;omega3]
swing_fault = @(t,x) swing_rhs(x, E, Pm, M, Gf, Bf);

% integration window and options
tspan = [0 0.6];         % 0.6s enough to capture example (increase if needed)
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3);

x0 = [theta0(:); zeros(3,1)];  % initial angles and zero speeds
[tf, xf] = ode45(swing_fault, tspan, x0, opts);

% postprocess: compute Pe(t), f(t), KE, V_PE(using post-fault network), total V
N = size(xf,1);
Pe_fault = zeros(N,3);
fvecs    = zeros(N,3);
VKE      = zeros(N,1);
VPE_post = zeros(N,1);
Vtot     = zeros(N,1);
ftheta_dot = zeros(N,1);

for k = 1:N
    thk = xf(k,1:3).';
    wk  = xf(k,4:6).';
    Pe_k = compute_Pe(thk, E, Gf, Bf);      % electrical power during fault
    Pe_fault(k,:) = Pe_k.';
    f_k = Pm - Pe_k;                        % Pm - Pe during fault
    fvecs(k,:) = f_k.';
    VKE(k) = 0.5 * sum( M .* (wk.^2) );     % kinetic energy
    % compute V_PE at this theta using POST-FAULT network (book computes V_PE wrt post-fault SEP)
    VPE_post(k) = compute_VPE_numeric(thk, theta_s, E, Gpost, Bpost, Pm, 300);
    Vtot(k) = VKE(k) + VPE_post(k);
    ftheta_dot(k) = f_k.' * wk;            % f^T * theta_dot  (scalar)
end

% find V_PE_max and its time (book: ~1.0377 at ~0.36s)
[VPE_max, idx_max] = max(VPE_post);
t_VPE_max = tf(idx_max);
fprintf('V_PE_max = %.6f at t = %.6f s (index %d)\n', VPE_max, t_VPE_max, idx_max);

% find first zero crossing of f^T * theta_dot (PEBS crossing)
signs = sign(ftheta_dot);
idx_cross = find(signs(1:end-1).*signs(2:end) <= 0, 1, 'first');
if isempty(idx_cross)
    warning('No zero crossing of f^T*theta_dot found in tspan.');
    t_star = NaN;
else
    % linear interpolation for t*
    t1 = tf(idx_cross); t2 = tf(idx_cross+1);
    y1 = ftheta_dot(idx_cross); y2 = ftheta_dot(idx_cross+1);
    t_star = t1 - y1*(t2-t1)/(y2-y1);
    fprintf('PEBS zero crossing at t* ≈ %.6f s\n', t_star);
end

% find t_cr: first time V_total reaches VPE_max (interpolated)
idx_cr = find(Vtot >= VPE_max, 1, 'first');
if isempty(idx_cr)
    warning('V_total never reached V_PE_max in tspan; increase tspan.');
    t_cr = NaN;
else
    if idx_cr == 1
        t_cr = tf(1);
    else
        % linear interp between tf(idx_cr-1) and tf(idx_cr)
        t1 = tf(idx_cr-1); t2 = tf(idx_cr);
        v1 = Vtot(idx_cr-1); v2 = Vtot(idx_cr);
        t_cr = t1 + (VPE_max - v1) * (t2-t1) / (v2 - v1);
    end
    fprintf('t_cr (first V = V_PE_max) ≈ %.6f s\n', t_cr);
end

% --- plots (book-like) ---
figure('Units','normalized','Position',[0.05 0.05 0.9 0.8]);
subplot(3,1,1);
plot(tf, xf(:,1:3));
grid on; legend('\theta_1','\theta_2','\theta_3'); ylabel('rad');
title('Rotor angles during fault (fault-on)');

subplot(3,1,2);
plot(tf, ftheta_dot);
hold on; plot([tf(1) tf(end)], [0 0], 'k--');
xlabel('time (s)'); ylabel('f^T \thetȧ');
title('f^T(\theta) \cdot \dot\theta'); grid on;
if ~isnan(t_star)
    plot(t_star, 0, 'ro', 'MarkerSize',8, 'LineWidth',1.5);
    text(t_star, 0, sprintf(' t* = %.3fs', t_star), 'VerticalAlignment','bottom');
end

subplot(3,1,3);
plot(tf, VPE_post, 'b-', tf, Vtot, 'r-'); grid on;
legend('V_{PE}(θ) [post-fault ref]','V(θ,ω)=KE+V_{PE}');
xlabel('time (s)');
ylabel('Energy (pu)');
title('Energy traces during fault');
% mark V_PE_max and t_cr
hold on;
plot(t_VPE_max, VPE_max, 'ks', 'MarkerFaceColor','k');
if ~isnan(t_cr)
    plot(t_cr, VPE_max, 'mo', 'MarkerFaceColor','m');
    text(t_cr, VPE_max, sprintf(' t_{cr}=%.3fs', t_cr), 'VerticalAlignment','bottom');
end

% Save results to struct for later use
faultTrajectory.tf = tf;
faultTrajectory.xf = xf;
faultTrajectory.Pe_fault = Pe_fault;
faultTrajectory.fvecs = fvecs;
faultTrajectory.VPE_post = VPE_post;
faultTrajectory.Vtot = Vtot;
faultTrajectory.ftheta_dot = ftheta_dot;
faultTrajectory.t_star = t_star;
faultTrajectory.VPE_max = VPE_max;
faultTrajectory.t_VPE_max = t_VPE_max;
faultTrajectory.t_cr = t_cr;

% --------------------- helper RHS used by ode45 -------------------------
function dx = swing_rhs(x, Evec, Pmvec, Mvec, Gmat, Bmat)
    th = x(1:3);
    w  = x(4:6);
    Pe_local = compute_Pe(th, Evec, Gmat, Bmat);
    dth = w;
    domega = (Pmvec - Pe_local) ./ Mvec;
    dx = [dth; domega];
end
