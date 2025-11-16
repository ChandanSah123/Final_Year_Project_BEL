close all; clc;
% Inertia constants (book example)
H= [23.64, 6.4, 3.01]; 
M=(2/(2*pi*60))*H;
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
%% there is certain problem either in code or in parameter value which is making it difficult to converge the sep to the exact book value
%the mismatch is about (-0.1782  0.5309  0.2711)-(-0.2023 0.5982 0.3154)=
%(-0.0241 -0.067 -0.0443)
%%calculation of fault on trajectory  