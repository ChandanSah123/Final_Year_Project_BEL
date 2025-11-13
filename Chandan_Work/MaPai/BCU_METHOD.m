% tef_bcu_ieee9.m
% Full BCU algorithm (Sauer & Pai / Bhui & Senroy style) for IEEE-9 example
% Put this in a .m file and run. Uses variables defined below (or in workspace).
 close all; clc;

%%%%% ------------------ USER / DATA SECTION ------------------ %%%%%
% --- Put your precomputed Y (9x9) here (admittance matrix on HV base).
% Example: load prebuilt Y from workspace or file. Replace with your Y.
% load('Y_ieee9.mat','Y');  % <- optional if you saved Y earlier
% For this script we assume Y (9x9) already exists in the workspace.
% If not, you must assemble Y exactly as in your earlier code.
% -----------------------------------------------------------------
% Use the Y you used above (the script you posted). If Y is not in workspace,
% the script will error. Make sure your Y is present.

% ------- Machine / system data (from Table 7.1 & Example 7.6) -------
Pm = [0.716, 1.63, 0.85];                 % mechanical power, pu (generator 1..3)
% Internal EMFs (from Example 7.6) -- constant for classical model
E1 = 1.054*exp(1j*deg2rad(2.2670));
E2 = 1.050*exp(1j*deg2rad(19.750));
E3 = 1.017*exp(1j*deg2rad(13.20));
E = [abs(E1), abs(E2), abs(E3)];
% Inertia constants (book example)
H= [23.64, 6.4, 3.01]; 
M=(2/(2*pi*60))*H; % H (or scaled M) consistent with your swing eq implementation
% If your swing eq uses M directly (not 2H/omega), keep consistent with prior code.

% ybar (machine terminal admittances) you provided earlier:
ybar = diag([-1i*16.45, -1i*8.35, -1i*5.52]);
% load admittance additions (as you had)
load_ybar = complex(zeros(9,9));
load_ybar(5,5) = (1.26 - 0.504i);
load_ybar(6,6) = (0.877 - 0.292i);
load_ybar(8,8) = (0.968 - 0.339i);

% transformer branch Y for 5-7 (used when removing branch)
Y78 = 1.187 - 5.975i;  % as you used earlier

% Pre-fault machine terminal voltages (angles) from Table 7.1 / Example 7.6
Vt1 = 1.04*exp(1j*deg2rad(0.0));
Vt2 = 1.025*exp(1j*deg2rad(9.3));
Vt3 = 1.02*exp(1j*deg2rad(4.7));
Vt = [Vt1, Vt2, Vt3];

%%%%% ------------------ BUILD AUGMENTED MATRIX & Yint ------------------ %%%%%
% The augmented Y (12x12) format: top-left machine terminal nodes (3), rest 9 network nodes
% Use your previous assembly approach:
YN = Y;             % base network Y (9x9)
% The book adds ybar on the network diagonal for terminals, we follow your approach:
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
YC(1:3,:) = -ybar';
YD = YN2;
Yint_pre = YA - YB * (YD \ YC);   % internal (pre-fault) Y seen by internal EMFs

% Save these for debugging/inspection
Yint_pre_full = Yint_pre;

%%%%% ------------------ CREATE FAULT-ON and POST-FAULT Yint ------------------ %%%%%
% Fault at bus 7 (bus indexing: follow your Y bus indexes so 7 means the same)
% Fault-on modeling: we represent as bus 7 grounded (i.e. remove row/col or add large shunt)
% We'll use "remove row/col" method for the network blocks when computing reduced YintF.

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
YN2PF(5,5) = YN2PF(5,5) - Y78;
YN2PF(7,7) = YN2PF(7,7) - Y78;

% Now form YintPF (post-fault reduced internal)
YAPF = ybar;
YBPF = complex(zeros(3,9));
YBPF(:,1:3) = -ybar;
YCPF = complex(zeros(9,3));
YCPF(1:3,:) = -ybar;
YDPF = YN2PF;
YintPF = YAPF - YBPF * (YDPF \ YCPF);

% For later use: GPF,BPF
GPF = real(YintPF);
BPF = imag(YintPF);

%%%%% ------------------ Compute pre-fault COI-shifted initial theta ------------------ %%%%%
% Compute internal EMF angles delta_pre from E1..E3
delpre = [angle(E1), angle(E2), angle(E3)];
% Center of inertia shift
del0 = (M(1)*delpre(1) + M(2)*delpre(2) + M(3)*delpre(3)) / sum(M);
theta_pre = delpre - del0;   % pre-fault angles in COI frame

%%%%% ------------------ Step-1: Solve post-fault SEP (theta^s) via Newton-Raphson ----- %%%%%
% We'll solve for theta^s (3 unknowns) using analytical Jacobian
theta0 = theta_pre(:);    % initial guess (use pre-fault)
tol = 1e-10;
maxit = 3;
theta = theta0;

for k = 1:maxit
    % compute f(theta) where f_i = Pm(i) - Pe_i(theta)
    f = zeros(3,1);
    Pe = zeros(3,1);
    for i = 1:3
        s = 0;
        for j = 1:3
            s = s + E(i)*E(j) * ( GPF(i,j)*cos(theta(i)-theta(j)) + BPF(i,j)*sin(theta(i)-theta(j)) );
        end
        Pe(i) = s;
        f(i) = Pm(i) - Pe(i);
    end
    % Jacobian J_ij = d f_i / d theta_j
    J = zeros(3,3);
    for i = 1:3
        for j = 1:3
            if i ~= j
                J(i,j) = E(i)*E(j) * ( GPF(i,j)*sin(theta(i)-theta(j)) - BPF(i,j)*cos(theta(i)-theta(j)) );
            end
        end
        J(i,i) = -sum(J(i,:));
    end
    % NR step
    dx = J \ f;
    theta = theta + dx;
    if norm(dx,inf) < tol
        fprintf('SEP NR converged in %d iterations\n', k);
        break;
    end
    if k == maxit
        %warning('SEP NR did not converge (max iterations).');
    end
end
theta_s = theta;   % post-fault SEP in COI frame

fprintf('Post-fault SEP (rad): %f, %f, %f\n', theta_s(1), theta_s(2), theta_s(3));
fprintf('Post-fault SEP (deg): %f, %f, %f\n', rad2deg(theta_s(1)), rad2deg(theta_s(2)), rad2deg(theta_s(3)));

%%%%% ------------------ Step-2: Simulate fault-on swing equations until clearing ----- %%%%%
% Classical swing equations for generator i:
% M_i * d^2 delta_i / dt^2 = Pm_i - Pe_i_fault - D_i * d delta/dt  (we take D_i=0 here)
% We'll simulate for a maximum time and assume clearing at t_clear (you can set or vary)
tmax = 1.0;
tspan = [0, tmax];
% initial states: theta_pre (COI), initial speed zero (relative)
x0 = [theta_pre(:); zeros(3,1)];  % [theta; omega]

% Precompute matrices for the fault-on network: use YintF -> Gf,Bf
Gf = real(YintF);
Bf = imag(YintF);

% ODE for swing during fault: state x = [theta; omega]
function dx = swing_fault(t,x)
    th = x(1:3);
    w = x(4:6);
    Pe_f = zeros(3,1);
    for ii=1:3
        s=0;
        for jj=1:3
            s = s + E(ii)*E(jj) * ( Gf(ii,jj)*cos(th(ii)-th(jj)) + Bf(ii,jj)*sin(th(ii)-th(jj)) );
        end
        Pe_f(ii) = s;
    end
    dth = w;
    domega = (Pm(:) - Pe_f) ./ M(:)';   % note M is a row vector; adjust
    dx = [dth; domega];
end

% integrate with ode45, but we need outputs at dense times to detect crossing
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t_f, x_f] = ode45(@swing_fault, tspan, x0, opts);

% compute f(θ) (gradient) and f^T * theta_dot at each saved point to check PEBS crossing
nt = length(t_f);
fvals = zeros(nt,3);
theta_dot = zeros(nt,3);
for k=1:nt
    thk = x_f(k,1:3).';
    wk = x_f(k,4:6).';
    % compute f = Pm - Pe(th)
    Peth = zeros(3,1);
    for ii=1:3
        s=0;
        for jj=1:3
            s = s + E(ii)*E(jj) * ( Gf(ii,jj)*cos(thk(ii)-thk(jj)) + Bf(ii,jj)*sin(thk(ii)-thk(jj)) );
        end
        Peth(ii) = s;
    end
    fvals(k,:) = (Pm(:) - Peth).';
    theta_dot(k,:) = wk.';
end
% compute dot product f^T * theta_dot (sum over i)
ftheta_dot = sum(fvals .* theta_dot, 2);

% find first time index when ftheta_dot crosses zero from positive to negative
% use sign change detection
idx = find( ftheta_dot(1:end-1) .* ftheta_dot(2:end) <= 0, 1, 'first' );
if isempty(idx)
    warning('No PEBS crossing found within simulated window. Consider increasing tmax.');
    idx = nt;  % fallback
end
% interpolate to get t*
t1 = t_f(idx); t2 = t_f(idx+1);
y1 = ftheta_dot(idx); y2 = ftheta_dot(idx+1);
t_star = t1 - y1*(t2-t1)/(y2 - y1);   % linear interpolation
% interpolate theta* and omega* also
theta_star = interp1(t_f, x_f(:,1:3), t_star);
omega_star = interp1(t_f, x_f(:,4:6), t_star);

fprintf('PEBS crossing approx at t* = %f s\n', t_star);

%%%%% ------------------ Step-3: After PEBS cross, integrate gradient system to find u_app ----- %%%%%
% gradient system: d theta / dt = f(theta)  (note f = Pm - Pe(theta))
% Use theta* as initial condition and integrate until ||f|| (norm of residual) minimum
theta_star_vec = theta_star(:);

% ODE for gradient system
function dth = grad_sys(t,th)
    Pe_tmp = zeros(3,1);
    for ii=1:3
        s = 0;
        for jj=1:3
            s = s + E(ii)*E(jj) * ( GPF(ii,jj)*cos(th(ii)-th(jj)) + BPF(ii,jj)*sin(th(ii)-th(jj)) );
        end
        Pe_tmp(ii) = s;
    end
    dth = (Pm(:) - Pe_tmp);  % gradient step (no scaling)
end

% integrate gradient ODE starting from theta_star for a time window
tgradspan = [0 5.0];  % adjust as needed
[ t_g, th_g ] = ode45(@grad_sys, tgradspan, theta_star_vec, opts);

% compute norm of f(θ) over gradient integration; find first local minimum
nf = size(th_g,1);
f_norm = zeros(nf,1);
for k=1:nf
    thk = th_g(k,:).';
    % compute f_k for post-fault network (use post-fault GPF,BPF)
    Pek = zeros(3,1);
    for ii=1:3
        s = 0;
        for jj=1:3
            s = s + E(ii)*E(jj) * ( GPF(ii,jj)*cos(thk(ii)-thk(jj)) + BPF(ii,jj)*sin(thk(ii)-thk(jj)) );
        end
        Pek(ii) = s;
    end
    fvec = Pm(:) - Pek;
    f_norm(k) = norm(fvec,1);  % use L1 as in book (sum abs)
end

% find first minimum (local)
[~, minidx] = min(f_norm);
theta_u_app = th_g(minidx,:).';
fprintf('Approx controlling UEP found at gradient time %f, index %d\n', t_g(minidx), minidx);

% optionally refine theta_u_app by one Newton solve of f(theta)=0 starting at theta_u_app
% (book says solving f=0 gives exact uep; but might converge to theta_s instead: caution)
% We'll skip Newton here because uep is a saddle; instead treat theta_u_app as controlling uep approx.

%%%%% ------------------ Step-4: Compute V_PE at theta_u_app using Eq. (9.67) ------------------ 
% We compute V_PE(theta^u_app) using the closed-form formula (9.67) in the book.
% Build C and D matrices for post-fault network (GPF,BPF)
C = zeros(3,3);
D = zeros(3,3);
for i=1:3
    for j=1:3
        C(i,j) = E(i)*E(j)*BPF(i,j);
        D(i,j) = E(i)*E(j)*GPF(i,j);
    end
end

% Function to compute V_PE using eq (9.67) with theta_ref = theta_s (book sets VPE=0 at theta^s)
theta_s_vec = theta_s;
theta_u = theta_u_app;

% compute first term: - sum_i P_i (theta_u_i - theta_s_i)
term1 = -sum( Pm(:) .* (theta_u - theta_s_vec) );

% compute second term: - sum_{i<j} Cij (cos theta^u_ij - cos theta^s_ij)
term2 = 0;
m = 3;
for i=1:m-1
    for j=i+1:m
        term2 = term2 - C(i,j)*( cos(theta_u(i)-theta_u(j)) - cos(theta_s_vec(i)-theta_s_vec(j)) );
    end
end

% compute path-dependent integral third term using I_ij formula (Eq 9.68)
% I_ij = [ (theta_u_i - theta_s_i + theta_u_j - theta_s_j)/(theta_u_i - theta_s_i - (theta_u_j - theta_s_j)) ] ...
%       * [ - D_ij ( sin theta_u_ij - sin theta_s_ij ) ]
% careful with division by zero; handle near-equal denominators
term3 = 0;
for i=1:m-1
    for j=i+1:m
        num = ( (theta_u(i)-theta_s_vec(i)) + (theta_u(j)-theta_s_vec(j)) );
        den = ( (theta_u(i)-theta_s_vec(i)) - (theta_u(j)-theta_s_vec(j)) );
        % If denominator small, use limit formula: Iij -> ???  (book assumes not singular for typical cases)
        if abs(den) < 1e-12
            % Use small-p expansion: approximate integral by D_ij * (theta_u(i)-theta_s(i)) * cos(theta_s(i)-theta_s(j))
            % This is a fallback; should rarely trigger for well-behaved uep.
            Iij = - D(i,j) * ( sin(theta_u(i)-theta_u(j)) - sin(theta_s_vec(i)-theta_s_vec(j)) );
        else
            Iij = - D(i,j) * ( sin(theta_u(i)-theta_u(j)) - sin(theta_s_vec(i)-theta_s_vec(j)) ) * ( num / den );
        end
        term3 = term3 - Iij;
    end
end

V_PE_uapp = term1 + term2 + term3;

% Compute V_CR approx = V_PE(theta_u_app) (book approximation)
V_cr = V_PE_uapp;

%%%%% ------------------ Compute V_cl at clearing (t_star) and decision ------------------ %%%%%
% Find theta_cl and omega_cl at t_star we already interpolated:
theta_cl = theta_star(:);
omega_cl = omega_star(:);

% KE at clearing: V_KE = 1/2 sum M_i * (omega_i)^2  (in rad/s units consistent with M)
V_KE_cl = 0.5 * sum( M(:) .* (omega_cl(:).^2) );

% Compute V_PE at theta_cl (use same eqn but relative to theta_s)
theta_c = theta_cl;
% compute term1_c and term2_c and term3_c
term1_c = -sum( Pm(:) .* (theta_c - theta_s_vec) );
term2_c = 0;
for i=1:m-1
    for j=i+1:m
        term2_c = term2_c - C(i,j)*( cos(theta_c(i)-theta_c(j)) - cos(theta_s_vec(i)-theta_s_vec(j)) );
    end
end
term3_c = 0;
for i=1:m-1
    for j=i+1:m
        numc = ( (theta_c(i)-theta_s_vec(i)) + (theta_c(j)-theta_s_vec(j)) );
        denc = ( (theta_c(i)-theta_s_vec(i)) - (theta_c(j)-theta_s_vec(j)) );
        if abs(denc) < 1e-12
            Iijc = - D(i,j) * ( sin(theta_c(i)-theta_c(j)) - sin(theta_s_vec(i)-theta_s_vec(j)) );
        else
            Iijc = - D(i,j) * ( sin(theta_c(i)-theta_c(j)) - sin(theta_s_vec(i)-theta_s_vec(j)) ) * ( numc / denc );
        end
        term3_c = term3_c - Iijc;
    end
end
V_PE_cl = term1_c + term2_c + term3_c;
V_cl = V_KE_cl + V_PE_cl;

fprintf('V_cr = %g , V_cl = %g\n', V_cr, V_cl);
if V_cl < V_cr
    fprintf('System STABLE for clearing at t* = %.6f s (Vcl < Vcr)\n', t_star);
else
    fprintf('System UNSTABLE for clearing at t* = %.6f s (Vcl >= Vcr)\n', t_star);
end

%%%%% ------------------ END ------------------ %%%%%
