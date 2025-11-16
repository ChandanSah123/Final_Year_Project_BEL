%% ======================== Y-BUS MATRIX ===============================

Y = [ ...
   -1i*17.361      0               0           1i*17.361       0               0               0               0               0;
    0             -1i*16           0           0               0               0               1i*16          0               0;
    0              0             -1i*17.065    0               0               0               0               0               1i*17.065;
    1i*17.361      0               0           (3.307 - 1i*39.309)   (-1.365 + 1i*11.604) (-1.942 + 1i*10.511) 0     0     0;
    0              0               0          (-1.365 + 1i*11.604)  (2.553 - 1i*17.338) 0 (-1.188 + 1i*5.975) 0 0;
    0              0               0          (-1.942 + 1i*10.511) 0 (3.224 - 1i*15.841) 0 0 (-1.282 + 1i*5.588);
    0              1i*16           0           0 (-1.188 + 1i*5.975) 0 (2.805 - 1i*35.446) (-1.617 + 1i*13.698) 0;
    0              0               0           0 0 0 (-1.617 + 1i*13.698) (2.772 - 1i*23.303) (-1.155 + 1i*9.784);
    0              0               1i*17.065   0 0 (-1.282 + 1i*5.588) 0 (-1.155 + 1i*9.784) (2.437 - 1i*32.154)
];

% Create readable MATLAB table
Y_table = array2table(Y, ...
    'VariableNames', compose('Bus%d', 1:9), ...
    'RowNames', compose('Bus%d', 1:9));

% Export to Excel
writetable(Y_table, 'Ybus_matrix.xlsx', 'WriteRowNames', true);
disp(Y_table);
save('Y.mat','Y');


%% ======================== DEFINE ybar & LOAD admittances ==============

ybar = diag([-1i*16.45, -1i*8.35, -1i*5.52]);

load_ybar = zeros(9,9);
load_ybar(5,5) = (1.26 - 0.504i);
load_ybar(6,6) = (0.877 - 0.292i);
load_ybar(8,8) = (0.969 - 0.339i);


%% ======================== BUILD AUGMENTED MATRIX =======================

YN = Y;

% Add ybar to Y diagonal (buses 1–3)
YNadd = zeros(9,9);
YNadd(1:3,1:3) = ybar;

YN1 = YN + YNadd;
YN2 = YN1 + load_ybar;

% 12×12 augmented matrix (machine terminals + network)
Y_aug = complex(zeros(12,12));

Y_aug(1:3,1:3) = ybar;
Y_aug(1:3,4:6) = -ybar;
Y_aug(4:6,1:3) = -ybar;
Y_aug(4:12,4:12) = YN2;


%% ======================== REDUCE TO Yint ===============================

YA = ybar;

YB = complex(zeros(3,9));
YB(:,1:3) = -ybar;

YC = complex(zeros(9,3));
YC(1:3,:) = -ybar;

YD = YN2;

YDC = YD \ YC;   % solve
Yint = YA - YB * YDC;
disp('Internal admittance matrix Yint:');
disp(Yint);


%% ======================== PRE–FAULT INTERNAL REDUCTION =================

Yint_pre = Yint;
G = real(Yint_pre);
B = imag(Yint_pre);


%% ======================== MACHINE PARAMETERS ===========================

H = [23.64, 6.4, 3.01];
M = (2/(2*pi*60)) * H;

Pm = [0.716, 1.63, 0.85];

E1 = 1.054 * exp(1j * deg2rad(2.2670));
E2 = 1.050 * exp(1j * deg2rad(19.750));
E3 = 1.017 * exp(1j * deg2rad(13.20));

E = [abs(E1), abs(E2), abs(E3)];
angleE = [angle(E1), angle(E2), angle(E3)];

% COI reference
delpre = angleE;
del0 = sum(delpre .* M) / sum(M);

theta0 = (delpre - del0).';

w10 = 0; w20 = 0; w30 = 0;


%% ======================== Cij AND Dij MATRICES =========================

C = zeros(3,3);
D = zeros(3,3);

for i = 1:3
    for j = 1:3
        C(i,j) = E(i)*E(j)*B(i,j);
        D(i,j) = E(i)*E(j)*G(i,j);
    end
end


%% ======================== NEWTON–RAPHSON SOLVER ========================

syms theta1 theta2 theta3

f1 = Pm(1) - D(1,1) - ( ...
      C(1,2)*sin(theta1-theta2) + D(1,2)*cos(theta1-theta2) ...
    + C(1,3)*sin(theta1-theta3) + D(1,3)*cos(theta1-theta3) );

f2 = Pm(2) - D(2,2) - ( ...
      C(2,1)*sin(theta2-theta1) + D(2,1)*cos(theta2-theta1) ...
    + C(2,3)*sin(theta2-theta3) + D(2,3)*cos(theta2-theta3) );

f3 = Pm(3) - D(3,3) - ( ...
      C(3,1)*sin(theta3-theta1) + D(3,1)*cos(theta3-theta1) ...
    + C(3,2)*sin(theta3-theta2) + D(3,2)*cos(theta3-theta2) );

F = [f1; f2; f3];
J = jacobian(F,[theta1,theta2,theta3]);

F_fun = matlabFunction(F,'Vars',{theta1,theta2,theta3});
J_fun = matlabFunction(J,'Vars',{theta1,theta2,theta3});

theta = theta0;
tol = 1e-12;
maxIter = 50;
damping = 0.7;

for k = 1:maxIter
    F_val = F_fun(theta(1),theta(2),theta(3));
    J_val = J_fun(theta(1),theta(2),theta(3));

    % Apply COI constraint
    J_mod = J_val;
    F_mod = F_val;

    J_mod(1,:) = M;
    F_mod(1)   = 0;

    delta = -J_mod \ F_mod;
    theta = theta + damping * delta;

    if norm(delta) < tol
        break;
    end
end

disp('Solved internal rotor angles (rad):');
disp(theta);

thetas1=-0.1782;
thetas2=0.5309;
thetas3=0.2711;
theta_s=[thetas1;thetas2;thetas3];
%%%% ---------------- Step-2: Fault-on trajectory + energy traces (PEBS) -------------
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
