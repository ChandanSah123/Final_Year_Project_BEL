
%% During Fault
Pm = [0.716, 1.63, 0.85];   
YaugF=Y_aug;
YaugF(10,:)=[];
YaugF(:,10)=[];

YNF=Y;
YNF(:,7)=[];
YNF(7,:)=[];

YN1F=YN1;
YN1F(:,7)=[];
YN1F(7,:)=[];

YN2F=YN2;
YN2F(:,7)=[];
YN2F(7,:)=[];

ybarF=ybar;
YAF=ybar;
YBF=complex(zeros(3,8));
YBF(1:3,1:3)=-ybarF;
YCF=complex(zeros(8,3));
YCF(1:3,1:3)=-ybar;
YDF=YN2F;
YDCF=YDF\YCF;
YintF=YAF-YBF*YDCF;


%% Post Fault

Y78= 1.187-5.975*1i;
YN2PF=YN2;
YN2PF(5,7)=0;
YN2PF(7,5)=0;
YN2PF(5,5)=YN2PF(5,5)-Y78;
YN2PF(7,7)=YN2PF(7,7)-Y78;

ybarPF=ybar;
YAPF=ybar;
YBPF=complex(zeros(3,9));
YBPF(1:3,1:3)=-ybarF;
YCPF=complex(zeros(9,3));
YCPF(1:3,1:3)=-ybar;
YDPF=YN2PF;
YDCPF=YDPF\YCPF;
YintPF=YAPF-YBPF*YDCPF;

%% Post Fault Electrical Power calculation
BPF=imag(YintPF);
GPF=real(YintPF);


E1=1.054*exp(1j*deg2rad(2.2670));
E2=1.050*exp(1j*deg2rad(19.750));
E3=1.017*exp(1j*deg2rad(13.20));

E=[abs(E1),abs(E2),abs(E3)];
delpre=[angle(E1),angle(E2),angle(E3)];
M=[23.64,6.4,3.01];

del0=(M(1)*delpre(1)+M(2)*delpre(2)+M(3)*delpre(3))/(M(1)+M(2)+M(3));

thetaP=[delpre(1)-del0,delpre(2)-del0,delpre(3)-del0];  % PREFAULT THEATA WITH RESPECT TO COI

% let the post fault theta be.....
syms th1 th2 th3

% Cij=EiEjBij and Dij=EiEjGij
C=zeros(3,3);
D=zeros(3,3);
for i=1:3
    for j=1:3
        C(i,j)=E(i)*E(j)*BPF(i,j);
        D(i,j)=E(i)*E(j)*GPF(i,j);
    end
end


   f1= Pm(1) - ( D(1,1) + ( C(1,2)*sin(th1-th2) + D(1,2)*cos(th1-th2) ) + ( C(1,3)*sin(th1-th3) + D(1,3)*cos(th1-th3) ) ); 
   f2=Pm(2) - ( D(2,2) + ( C(2,1)*sin(th2-th1) + D(2,1)*cos(th2-th1) ) + ( C(2,3)*sin(th2-th3) + D(2,3)*cos(th2-th3) ) );
   f3=Pm(3) - ( D(3,3) + ( C(3,1)*sin(th3-th1) + D(3,1)*cos(th3-th1) ) + ( C(3,2)*sin(th3-th2) + D(3,2)*cos(th3-th2) ) );


   %initial guesses delth1 delth2 delth3

  Jf11= diff(f1,th1);
  Jf12=diff(f1,th2);
  Jf13=diff(f1,th3);

  Jf21= diff(f2,th1);
  Jf22=diff(f2,th2);
  Jf23=diff(f2,th3);

  Jf31= diff(f3,th1);
  Jf32=diff(f3,th2);
  Jf33=diff(f3,th3);


  % chatgpt hero
%% Robust Newton-Raphson with COI anchoring, step-clamp and angle-wrapping
% Assumes: Pm, C, D, thetaP, M are defined

% Initial guess (radians) — use prefault COI relative angles
x = double(real(thetaP));        % [th1 th2 th3]
x = mod(x + pi, 2*pi) - pi;      % normalize to (-pi, pi]

tol = 1e-8;
maxIter = 100;
min_alpha = 1e-6;
c_red = 0.5;            % line search shrink
verbose = true;
step_max = 0.2;         % max allowed step in radians (adjustable)

% helper residual function
residual = @(xx) [
    Pm(1) - ( D(1,1) + ( C(1,2)*sin(xx(1)-xx(2)) + D(1,2)*cos(xx(1)-xx(2)) ) + ( C(1,3)*sin(xx(1)-xx(3)) + D(1,3)*cos(xx(1)-xx(3)) ) );
    Pm(2) - ( D(2,2) + ( C(2,1)*sin(xx(2)-xx(1)) + D(2,1)*cos(xx(2)-xx(1)) ) + ( C(2,3)*sin(xx(2)-xx(3)) + D(2,3)*cos(xx(2)-xx(3)) ) );
    Pm(3) - ( D(3,3) + ( C(3,1)*sin(xx(3)-xx(1)) + D(3,1)*cos(xx(3)-xx(1)) ) + ( C(3,2)*sin(xx(3)-xx(2)) + D(3,2)*cos(xx(3)-xx(2)) ) )
];

% COI weight normalization function (subtract weighted average)
coi_subtract = @(xx) ( (M(1)*xx(1)+M(2)*xx(2)+M(3)*xx(3)) / (M(1)+M(2)+M(3)) );

% initial residual
F = residual(x);
res_norm = norm(F);
if verbose, fprintf('Initial residual norm = %.3e\n', res_norm); end

converged = false;
for k = 1:maxIter
    th1 = x(1); th2 = x(2); th3 = x(3);

    % ---- Build Jacobian (analytical) ----
    J = zeros(3,3);
    J(1,1) = -( C(1,2)*cos(th1-th2) - D(1,2)*sin(th1-th2) ) ...
              -( C(1,3)*cos(th1-th3) - D(1,3)*sin(th1-th3) );
    J(1,2) =  ( C(1,2)*cos(th1-th2) - D(1,2)*sin(th1-th2) );
    J(1,3) =  ( C(1,3)*cos(th1-th3) - D(1,3)*sin(th1-th3) );

    J(2,1) =  ( C(2,1)*cos(th2-th1) - D(2,1)*sin(th2-th1) );
    J(2,2) = -( C(2,1)*cos(th2-th1) - D(2,1)*sin(th2-th1) ) ...
              -( C(2,3)*cos(th2-th3) - D(2,3)*sin(th2-th3) );
    J(2,3) =  ( C(2,3)*cos(th2-th3) - D(2,3)*sin(th2-th3) );

    J(3,1) =  ( C(3,1)*cos(th3-th1) - D(3,1)*sin(th3-th1) );
    J(3,2) =  ( C(3,2)*cos(th3-th2) - D(3,2)*sin(th3-th2) );
    J(3,3) = -( C(3,1)*cos(th3-th1) - D(3,1)*sin(th3-th1) ) ...
              -( C(3,2)*cos(th3-th2) - D(3,2)*sin(th3-th2) );

    % ---- Defensive checks ----
    if any(isnan(J(:))) || any(isinf(J(:)))
        warning('Jacobian contains NaN/Inf at iter %d. Aborting.', k); break;
    end

    % Regularize if near-singular
    rcondJ = rcond(J);
    if rcondJ < 1e-12
        eps_reg = 1e-6 * max(1, max(abs(diag(J))));
        J_reg = J + eps_reg * eye(3);
    else
        J_reg = J;
    end

    % ---- Newton step ----
    dx = - (J_reg \ F);   % column vector (3x1)

    % clamp the step to avoid huge jumps
    if norm(dx, inf) > step_max
        dx = dx * (step_max / norm(dx, inf));
        if verbose, fprintf('Step clamped to %.3f rad (iter %d)\n', step_max, k); end
    end

    % ---- Damped update with backtracking ----
    alpha = 1.0;
    x_trial = x + (alpha * dx).';    % dx is column, x row
    % anchor to remove COI rigid rotation BEFORE evaluating residual:
    x_trial = x_trial - coi_subtract(x_trial);
    % wrap angles
    x_trial = mod(x_trial + pi, 2*pi) - pi;

    F_trial = residual(x_trial);
    while norm(F_trial) >= norm(F) && alpha > min_alpha
        alpha = alpha * c_red;
        x_trial = x + (alpha * dx).';
        x_trial = x_trial - coi_subtract(x_trial);
        x_trial = mod(x_trial + pi, 2*pi) - pi;
        F_trial = residual(x_trial);
    end

    % If no improvement, fallback to pseudoinverse step with heavy damping
    if norm(F_trial) >= norm(F)
        dx_alt = - pinv(J_reg) * F;
        dx_alt = dx_alt * min(1, step_max / max(1e-12, norm(dx_alt, inf)));
        alpha = 0.5;
        x_trial = x + (alpha * dx_alt).';
        x_trial = x_trial - coi_subtract(x_trial);
        x_trial = mod(x_trial + pi, 2*pi) - pi;
        F_trial = residual(x_trial);
    end

    % ---- Accept update ----
    x = x_trial;
    F = F_trial;
    res_norm = norm(F);

    if verbose
        fprintf('Iter %2d: alpha=%.4f, rcond(J)=%.2e, res=%.3e, th=[%.6f, %.6f, %.6f]\n', ...
                k, alpha, rcondJ, res_norm, x(1), x(2), x(3));
    end

    % ---- Convergence test ----
    if norm(dx, inf) < tol && res_norm < 1e-6
        converged = true;
        fprintf('Converged in %d iterations. Final residual = %.3e\n', k, res_norm);
        break;
    end
end

if ~converged
    warning('Newton-Raphson did not converge in %d iterations. Final residual = %.3e', maxIter, res_norm);
end

% Final results (angle vector already COI-anchored and wrapped)
theta_sol = x;
fprintf('\nFinal Angles (radians):\n');
fprintf('θ1 = %.6f\nθ2 = %.6f\nθ3 = %.6f\n', theta_sol(1), theta_sol(2), theta_sol(3));
fprintf('\nFinal Angles (degrees):\n');
fprintf('θ1 = %.6f°\nθ2 = %.6f°\nθ3 = %.6f°\n', rad2deg(theta_sol(1)), rad2deg(theta_sol(2)), rad2deg(theta_sol(3)));




