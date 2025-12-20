%% Settings
t_fault=1; % time that the fault is applied
t_cl=0.168; % clearing time

% Time domain simulation setting (TDS)
t_end=4; % end time for time domain simulation
del_t=0.001; % times tep for tds
del_t_fault=0.001; % times step for fault tds
del_t_BCU=0.01; % times step for BCU tds

num_bus=39;
num_gen=10;
slack_bus=31;
plot_bus=2;
%% Set up variables
sys_case=39;
H = [42.0; 30.3; 35.8; 28.6; 26.0; 34.8; 26.4; 24.3; 34.5; 500.0];
M_gen=H.*(1/(pi*60));
M_T=sum(M_gen);
D_gen = zeros(10,1);
v_gen=[];
xd_p = [0.031; 0.0697; 0.0531; 0.0436; 0.132; 0.050; 0.049; 0.057; 0.057; 0.006];
idx_gen=[30 31 32 33 34 35 36 37 38 39];
idx_load=setdiff(1:num_bus,idx_gen);
idx_delta=1:num_gen; idx_omega=num_gen+1:2*num_gen; num_var=2*num_gen;
%% theta from Ma Pai
M=M_gen;
%E = [abs(E1), abs(E2), abs(E3)];
E=[ 1.0562;
    1.1782;
    1.1507;
    1.1004;
    1.3594;
    1.1806;
    1.1288;
    1.0721;
    1.1382;
    1.0215
  ];

%angleE=[angle(E1),angle(E2),angle(E3)];
%angleE=[-0.0440; 0.3124; 0.3554; 0.3292; 0.653;  0.328;  0.3283; 0.3261; 0.5540;  -0.1722];
angleE=data(5,)
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
%%
% Apply Kron Reduction
Y_pre_kron=Yint_pre;
Y_fault_kron=Yint_fault;
Y_post_kron=Yint_post;
%Pm =[2.50; 5.71; 6.50; 6.32; 5.08; 6.50; 5.60; 5.40; 8.30; 10.00];
Pm=data(3,11:20);
% Calculate initial electrical power using the Y-bus and initial voltages
V_complex = E .* exp(1j * angleE);
I_inj = Y_pre_kron * V_complex;
S_inj = V_complex .* conj(I_inj);
P_electrical_initial = real(S_inj);

% FORCE mechanical power to match this initial state
Pgen_pre = P_electrical_initial; 
Pgen_post = Pgen_pre; % Assume constant mechanical power
Pgen= Pm;
%Sload=zeros(num_bus,1); Sload(PQ.con(:,1))=(PQ.con(:,4)+1i*PQ.con(:,5));
% Pre-contingency Equilibrium
%Pgen_pre=data(1,4:6)';
Eeq_pre=E;
%delta_eq_pre=data(5,1:3)';
delta_eq_pre=delpre;

x_eq_pre=[delta_eq_pre-M_gen'*delta_eq_pre/M_T; zeros(num_gen,1)];

% Post-contingency Equilibrium

%Pgen_post=data(135,4:6)';
Eeq_post=Eeq_pre;
%delta_eq_post=data(130,1:3)';




%Y_pre_kron=Yint_pre;
%Y_fault_kron=YintF;
%Y_post_kron=YintPF;
edge_kron=nchoosek(1:num_gen,2);
E_kron=zeros(size(edge_kron,1),num_gen);
for i=1:size(E_kron,1)
    E_kron(i,edge_kron(i,1))=1;  E_kron(i,edge_kron(i,2))=-1;
    gij_pre_kron(i,1)=real(Y_pre_kron(E_kron(i,:)==1,E_kron(i,:)==-1))*Eeq_pre(E_kron(i,:)==1)*Eeq_pre(E_kron(i,:)==-1);
    gij_fault_kron(i,1)=real(Y_fault_kron(E_kron(i,:)==1,E_kron(i,:)==-1))*Eeq_post(E_kron(i,:)==1)*Eeq_post(E_kron(i,:)==-1);
    gij_post_kron(i,1)=real(Y_post_kron(E_kron(i,:)==1,E_kron(i,:)==-1))*Eeq_post(E_kron(i,:)==1)*Eeq_post(E_kron(i,:)==-1);
    bij_pre_kron(i,1)=imag(Y_pre_kron(E_kron(i,:)==1,E_kron(i,:)==-1))*Eeq_pre(E_kron(i,:)==1)*Eeq_pre(E_kron(i,:)==-1);
    bij_fault_kron(i,1)=imag(Y_fault_kron(E_kron(i,:)==1,E_kron(i,:)==-1))*Eeq_post(E_kron(i,:)==1)*Eeq_post(E_kron(i,:)==-1);
    bij_post_kron(i,1)=imag(Y_post_kron(E_kron(i,:)==1,E_kron(i,:)==-1))*Eeq_post(E_kron(i,:)==1)*Eeq_post(E_kron(i,:)==-1);
end
P_pre_kron=Pgen_pre-real(diag(Y_pre_kron)).*Eeq_pre.^2;
P_fault_kron=Pgen_post-real(diag(Y_fault_kron)).*Eeq_post.^2;
P_post_kron=Pgen_post-real(diag(Y_post_kron)).*Eeq_post.^2;

 %% Model
f_pre_kron= @(x) [x(idx_omega)-1/M_T*M_gen'*x(idx_omega); diag(1./M_gen)*(-D_gen.*x(idx_omega)-E_kron'*diag(bij_pre_kron)*sin(E_kron*x(idx_delta))-abs(E_kron)'*diag(gij_pre_kron)*cos(E_kron*x(idx_delta))+P_pre_kron)-1/M_T*(sum(P_pre_kron)-2*ones(size(E_kron'))*diag(gij_pre_kron)*cos(E_kron*x(idx_delta)))];
J_pre_kron=@(x) [zeros(num_gen) eye(num_gen)-ones(num_gen,1)*M_gen'/M_T; diag(1./M_gen)*(-E_kron'*diag(bij_pre_kron)*diag(cos(E_kron*x(idx_delta)))*E_kron+abs(E_kron)'*diag(gij_pre_kron)*diag(sin(E_kron*x(idx_delta)))*E_kron)-2/M_T*ones(size(E_kron'))*diag(gij_pre_kron)*diag(sin(E_kron*x(idx_delta)))*E_kron -diag(D_gen./M_gen)];
f_fault_kron= @(x) [x(idx_omega)-1/M_T*M_gen'*x(idx_omega); diag(1./M_gen)*(-D_gen.*x(idx_omega)-E_kron'*diag(bij_fault_kron)*sin(E_kron*x(idx_delta))-abs(E_kron)'*diag(gij_fault_kron)*cos(E_kron*x(idx_delta))+P_fault_kron)-1/M_T*(sum(P_pre_kron)-2*ones(size(E_kron'))*diag(gij_fault_kron)*cos(E_kron*x(idx_delta)))];
J_fault_kron=@(x) [zeros(num_gen) eye(num_gen)-ones(num_gen,1)*M_gen'/M_T; diag(1./M_gen)*(-E_kron'*diag(bij_fault_kron)*diag(cos(E_kron*x(idx_delta)))*E_kron+abs(E_kron)'*diag(gij_fault_kron)*diag(sin(E_kron*x(idx_delta)))*E_kron)-2/M_T*ones(size(E_kron'))*diag(gij_fault_kron)*diag(sin(E_kron*x(idx_delta)))*E_kron -diag(D_gen./M_gen)];
f_post_kron= @(x) [x(idx_omega)-1/M_T*M_gen'*x(idx_omega); diag(1./M_gen)*(-D_gen.*x(idx_omega)-E_kron'*diag(bij_post_kron)*sin(E_kron*x(idx_delta))-abs(E_kron)'*diag(gij_post_kron)*cos(E_kron*x(idx_delta))+P_post_kron)-1/M_T*(sum(P_post_kron)-2*ones(size(E_kron'))*diag(gij_post_kron)*cos(E_kron*x(idx_delta)))];
J_post_kron=@(x) [zeros(num_gen) eye(num_gen)-ones(num_gen,1)*M_gen'/M_T; diag(1./M_gen)*(-E_kron'*diag(bij_post_kron)*diag(cos(E_kron*x(idx_delta)))*E_kron+abs(E_kron)'*diag(gij_post_kron)*diag(sin(E_kron*x(idx_delta)))*E_kron)-2/M_T*ones(size(E_kron'))*diag(gij_post_kron)*diag(sin(E_kron*x(idx_delta)))*E_kron -diag(D_gen./M_gen)];

f_post_theta= @(x) -E_kron'*diag(bij_post_kron)*sin(E_kron*x)-abs(E_kron)'*diag(gij_post_kron)*cos(E_kron*x)+P_post_kron-(M_gen/M_T).*(sum(P_post_kron)-2*ones(size(E_kron'))*diag(gij_post_kron)*cos(E_kron*x));
J_post_theta=@(x) -E_kron'*diag(bij_post_kron)*diag(cos(E_kron*x))*E_kron+abs(E_kron)'*diag(gij_post_kron)*diag(sin(E_kron*x))*E_kron-2*diag(M_gen)/M_T*ones(size(E_kron'))*diag(gij_post_kron)*diag(sin(E_kron*x))*E_kron;

%% ---- Solve post-fault SEP using NR ----
theta01 = x_eq_pre(idx_delta);   % initial guess from PSSE
theta_sep = NR(f_post_theta, J_post_theta, theta01);

if any(isinf(theta_sep))
    error('Post-fault SEP did not converge.');
end

x_eq_post = [theta_sep; zeros(num_gen,1)];

%% Direct Simulation
tic; num_sim_idx=ceil((t_end-t_cl)/del_t)+ceil(t_cl/del_t_fault);
x_sim=zeros(num_var,num_sim_idx); t_sim=zeros(1,num_sim_idx);
x_sim(:,1)=x_eq_pre; t_sim(1)=0; t_idx=1;

while t_sim(t_idx)<=t_end
    if t_sim(t_idx)<t_fault
        f=@(x) -x+x_sim(:,t_idx)+del_t/2*(f_pre_kron(x_sim(:,t_idx))+f_pre_kron(x));
        J=@(x) -eye(num_var)+del_t/2*J_pre_kron(x);
        t_sim(t_idx+1)=t_sim(t_idx)+del_t;
    elseif  t_sim(t_idx)<t_fault+t_cl
        f=@(x) -x+x_sim(:,t_idx)+del_t_fault/2*(f_fault_kron(x_sim(:,t_idx))+f_fault_kron(x));
        J=@(x) -eye(num_var)+del_t_fault/2*J_fault_kron(x);
        t_sim(t_idx+1)=t_sim(t_idx)+del_t_fault;
    else
        f=@(x) -x+x_sim(:,t_idx)+del_t/2*(f_post_kron(x_sim(:,t_idx))+f_post_kron(x));
        J=@(x) -eye(num_var)+del_t/2*J_post_kron(x);
        t_sim(t_idx+1)=t_sim(t_idx)+del_t;
    end
    x_sim(:,t_idx+1)=NR(f,J,x_sim(:,t_idx));
    t_idx=t_idx+1;
end
t_sim=t_sim(1:t_idx); x_sim=x_sim(:,1:t_idx);
tcomp_TDS=toc;
result_TDS=(max(x_sim(idx_delta,t_idx-1))-min(x_sim(idx_delta,t_idx-1)))<2*pi;
%% PEBS (Potential Energy Boundary Surface)
tic; x_fault=x_eq_pre;
max_PEBS_sim=10000;
V_ke=zeros(1,max_PEBS_sim); V_p=zeros(1,max_PEBS_sim); V_d=zeros(1,max_PEBS_sim);
V_ke(1)=0.5*sum(M_gen.*x_fault(idx_omega,1).^2);
V_p(1)=-P_post_kron'*(x_fault(idx_delta,1)-x_eq_post(idx_delta))-bij_post_kron'*(cos(E_kron*x_fault(idx_delta,1))-cos(E_kron*x_eq_post(idx_delta)));
V_d(1)=0;
V_pe0=-P_post_kron'*(x_eq_pre(idx_delta)-x_eq_post(idx_delta))-bij_post_kron'*(cos(E_kron*x_eq_pre(idx_delta))-cos(E_kron*x_eq_post(idx_delta)))+0.5*(gij_post_kron.*(cos(E_kron*x_eq_pre(idx_delta))+cos(E_kron*x_eq_post(idx_delta))))'*(abs(E_kron)*x_eq_pre(idx_delta)-abs(E_kron)*x_eq_post(idx_delta));
V_total=V_ke+V_p+V_d;
PEBS=zeros(1,max_PEBS_sim);
PEBS(1)=f_post_theta(x_fault(idx_delta,1))'*(x_fault(idx_delta,1)-x_eq_post(idx_delta));
PEBS_idx=0; t_idx=1;

idx_tcl=round((t_cl)/del_t_fault);
idx_tf=max(round(t_fault/del_t_fault),1);

while t_idx<max_PEBS_sim && (PEBS_idx==0 || t_idx<=idx_tcl)
    f=@(x) -x+x_fault(:,t_idx)+del_t_fault/2*(f_fault_kron(x_fault(:,t_idx))+f_fault_kron(x));
    J=@(x) -eye(num_var)+del_t_fault/2*J_fault_kron(x);
    x_fault(:,t_idx+1)=NR(f,J,x_fault(:,t_idx));
    V_ke(t_idx+1)=0.5*sum(M_gen.*x_fault(idx_omega,t_idx+1).^2);
    V_p(t_idx+1)=-P_post_kron'*(x_fault(idx_delta,t_idx+1)-x_eq_post(idx_delta))-bij_post_kron'*(cos(E_kron*x_fault(idx_delta,t_idx+1))-cos(E_kron*x_eq_post(idx_delta)));
    V_d(t_idx+1)=V_d(t_idx)+0.5*(gij_post_kron.*(cos(E_kron*x_fault(idx_delta,t_idx+1))+cos(E_kron*x_fault(idx_delta,t_idx))))'*(abs(E_kron)*x_fault(idx_delta,t_idx+1)-abs(E_kron)*x_fault(idx_delta,t_idx));
    V_total(t_idx+1)=V_ke(t_idx+1)+V_p(t_idx+1)+V_d(t_idx+1);
    PEBS(t_idx+1)=f_post_theta(x_fault(idx_delta,t_idx+1))'*(x_fault(idx_delta,t_idx+1)-x_eq_post(idx_delta));
    if PEBS(t_idx)<0 && PEBS(t_idx+1)>=0 && PEBS_idx==0; PEBS_idx=t_idx; end
    t_idx=t_idx+1;
end
t_fault_sim=(0:t_idx-1)*del_t_fault; V_ke=V_ke(:,1:t_idx); V_p=V_p(:,1:t_idx); V_d=V_d(:,1:t_idx); V_total=V_total(:,1:t_idx); PEBS=PEBS(:,1:t_idx);

Vcr_PEBS=V_p(PEBS_idx)+V_d(PEBS_idx)-V_pe0; 
idx_tcr_PEBS=find(V_total<Vcr_PEBS); idx_tcr_PEBS=idx_tcr_PEBS(end);
tcr_PEBS=del_t_fault*idx_tcr_PEBS;
result_PEBS=V_total(idx_tcl)<Vcr_PEBS;
tcomp_PEBS=toc;


%% BCU (Boundary Controlling u.e.p. Method)
tic; t_idx=1; BCU_idx=0;
max_BCU_sim=round(5/del_t_BCU);
theta_sim=zeros(num_gen,max_BCU_sim); theta_abs=zeros(1,max_BCU_sim);
theta_sim(:,1)=x_fault(idx_delta,PEBS_idx); theta_abs(1)=sum(abs(f_post_theta(theta_sim(:,1))));

while t_idx<max_BCU_sim && PEBS_idx~=0 && BCU_idx==0
    f=@(x) -x+theta_sim(:,t_idx)+del_t_BCU/2*(f_post_theta(theta_sim(:,t_idx))+f_post_theta(x)); % Trapezoidal rule
    J=@(x) -eye(num_gen)+del_t_BCU/2*J_post_theta(x);
    theta_sim(:,t_idx+1)=NR(f,J,theta_sim(:,t_idx));
    theta_abs(t_idx+1)=sum(abs(f_post_theta(theta_sim(:,t_idx+1))));
    t_idx=t_idx+1;
    if theta_abs(t_idx)>theta_abs(t_idx-1) && theta_abs(max(1,t_idx-2))>theta_abs(t_idx-1); BCU_idx=t_idx-2; end
end
t_bcu_sim=(0:t_idx-1)*del_t_BCU; theta_sim=theta_sim(:,1:t_idx); theta_abs=theta_abs(:,1:t_idx);

theta_u=theta_sim(:,BCU_idx);
x_eq_unstable=[theta_u; zeros(num_gen,1)];
Vcr_BCU=-P_post_kron'*(theta_u-x_eq_post(idx_delta))-bij_post_kron'*(cos(E_kron*theta_u)-cos(E_kron*x_eq_post(idx_delta)))+((abs(E_kron)*(theta_u-x_eq_post(idx_delta)))./(E_kron*(theta_u-x_eq_post(idx_delta))))'*diag(gij_post_kron)*(sin(E_kron*theta_u)-sin(E_kron*x_eq_post(idx_delta)))-V_pe0;
idx_tcr_BCU=find(V_total<Vcr_BCU); idx_tcr_BCU=idx_tcr_BCU(end);
tcr_BCU=del_t_fault*idx_tcr_BCU;
result_BCU=V_total(idx_tcl)<Vcr_BCU;
tcomp_BCU=toc; tcomp_BCU=tcomp_BCU+tcomp_PEBS;

%% Display Plot
figure;
subplot(2,1,1); hold all; grid on; box on;
plot(t_sim,x_sim(idx_delta,:))
xlim([0 t_end]); set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('time (sec)'); ylabel('\delta');
subplot(2,1,2); hold all; grid on; box on;
plot(t_sim,x_sim(idx_omega,:))
xlim([0 t_end]); set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('time (sec)'); ylabel('\omega');

figure;
subplot(3,1,1); hold all; grid on; box on;
plot(t_fault_sim,V_ke,'r--');
plot(t_fault_sim,V_p+V_d,'b--');
plot(t_fault_sim,V_total);
% xlim([t_fault t_fault+t_cl])
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('time (sec)'); ylabel('Energy function');
legend('V_{ke}','V_{pe}','V_{ke}+V_{pe}')

subplot(3,1,2); hold all; grid on; box on;
plot(t_fault_sim,PEBS);
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('time (sec)');  ylabel('PEBS function');

subplot(3,1,3); hold all; grid on; box on;
plot(t_bcu_sim,theta_abs);
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('iteration');  ylabel('|f|');


% Phase Portrait
figure;
plot_axis=[plot_bus plot_bus+num_gen];
plot_rng=[-pi pi -10 10];
x_plot_pre=x_eq_pre;
x_plot_post=x_eq_post;
[x1_plot_ij,x2_plot_ij]=meshgrid(linspace(plot_rng(1),plot_rng(2),30),linspace(plot_rng(3),plot_rng(4),30));
for i=1:size(x1_plot_ij,1)
    for j=1:size(x2_plot_ij,2)
        x_plot_pre(plot_axis)=[x1_plot_ij(i,j); x2_plot_ij(i,j)];
        x_plot_post(plot_axis)=[x1_plot_ij(i,j); x2_plot_ij(i,j)];
        dx_plot_pre=f_fault_kron(x_plot_pre);
        dx_plot_post=f_post_kron(x_plot_post);
        dx1_plot_pre(i,j)=dx_plot_pre(plot_axis(1));
        dx2_plot_pre(i,j)=dx_plot_pre(plot_axis(2));
        dx1_plot_post(i,j)=dx_plot_post(plot_axis(1));
        dx2_plot_post(i,j)=dx_plot_post(plot_axis(2));
        num_prevent_zero=E_kron*(x_plot_post(idx_delta)-x_eq_post(idx_delta));
        num_prevent_zero(num_prevent_zero==0)=inf;
        V_plot(i,j)=0.5*sum(M_gen.*x_plot_post(idx_omega).^2)-P_post_kron'*(x_plot_post(idx_delta)-x_eq_post(idx_delta))-bij_post_kron'*(cos(E_kron*x_plot_post(idx_delta))-cos(E_kron*x_eq_post(idx_delta)))+((abs(E_kron)*(x_plot_post(idx_delta)-x_eq_post(idx_delta)))./(num_prevent_zero))'*diag(gij_post_kron)*(sin(E_kron*x_plot_post(idx_delta))-sin(E_kron*x_eq_post(idx_delta)));
    end
end
subplot(1,2,1); hold all; grid on; box on;
plot(x_sim(plot_axis(1),idx_tf : (idx_tf+idx_tcl)), x_sim(plot_axis(2),idx_tf : (idx_tf+idx_tcl)),'LineWidth',3);
streamslice(x1_plot_ij,x2_plot_ij,dx1_plot_pre,dx2_plot_pre);
axis(plot_rng); 
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('\delta (rad)'); ylabel('\omega (rad/s)'); title('Fault-on Dynamics')
subplot(1,2,2); hold all; grid on; box on;
scatter(theta_u(plot_axis(1)),0,100,'r','filled')
scatter(x_fault(plot_axis(1),PEBS_idx),0,100,'b','filled')
plot(x_sim(plot_axis(1),:),x_sim(plot_axis(2),:),'b','LineWidth',2);
streamslice(x1_plot_ij,x2_plot_ij,dx1_plot_post,dx2_plot_post);
axis(plot_rng); 
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('\delta (rad)'); ylabel('\omega (rad/s)'); title('Post-contingency Dynamics'); legend('u.e.p. (BCU)','u.e.p. (PEBS)')



%% Corrected Hybrid Method (Refining the Result)
tic
disp('--- Starting Hybrid Search (TDS Refinement) ---');

% 1. Set the search window based on your BCU (low) and PEBS (high) results
t_lower = tcr_BCU;       % Safe lower bound (0.163s)
t_upper = tcr_PEBS;      % Unsafe upper bound (0.184s)
precision = 0.001;       % Stop when we are within 1ms accuracy

% 2. Loop until we narrow it down
while (t_upper - t_lower) > precision
    
    % Pick a test clearing time in the middle
    t_test = (t_upper + t_lower) / 2 ;
    
    % --- RUN MINI-TDS FOR t_test ---
    % Initialize simulation variables for this specific run
    num_sim_hybrid = ceil((t_end-t_test)/del_t) + ceil(t_test/del_t_fault);
    x_sim = zeros(num_var, num_sim_hybrid); 
    t_sim = zeros(1, num_sim_hybrid);
    x_sim(:,1) = x_eq_pre; 
    t_sim(1) = 0; 
    t_idx = 1;
    
    % Run the simulation loop
    while t_sim(t_idx) <= t_end
        if t_sim(t_idx) < t_fault
            % Pre-fault
            f=@(x) -x+x_sim(:,t_idx)+del_t/2*(f_pre_kron(x_sim(:,t_idx))+f_pre_kron(x));
            J=@(x) -eye(num_var)+del_t/2*J_pre_kron(x);
            dt_curr = del_t;
        elseif t_sim(t_idx) < t_fault + t_test  % <--- Using t_test here
            % Fault-on
            f=@(x) -x+x_sim(:,t_idx)+del_t_fault/2*(f_fault_kron(x_sim(:,t_idx))+f_fault_kron(x));
            J=@(x) -eye(num_var)+del_t_fault/2*J_fault_kron(x);
            dt_curr = del_t_fault;
        else
            % Post-fault
            f=@(x) -x+x_sim(:,t_idx)+del_t/2*(f_post_kron(x_sim(:,t_idx))+f_post_kron(x));
            J=@(x) -eye(num_var)+del_t/2*J_post_kron(x);
            dt_curr = del_t;
        end
        
        % Solve step
        x_sim(:,t_idx+1) = NR(f,J,x_sim(:,t_idx));
        t_sim(t_idx+1) = t_sim(t_idx) + dt_curr;
        t_idx = t_idx + 1;
        
        % Optimization: Stop early if already unstable (angles blow up)
        curr_spread = max(x_sim(idx_delta,t_idx)) - min(x_sim(idx_delta,t_idx));
        if curr_spread > 2*pi
            break; % Stop simulation, we know it's unstable
        end
    end
    
    % --- CHECK STABILITY ---
    final_spread = max(x_sim(idx_delta,t_idx)) - min(x_sim(idx_delta,t_idx));
    is_stable = final_spread < 2*pi;
    
    % --- UPDATE BOUNDS ---
    if is_stable
        % System survived. Try a longer duration.
        t_lower = t_test;
        % fprintf('t = %.4f : Stable\n', t_test);
    else
        % System crashed. Try a shorter duration.
        t_upper = t_test;
        % fprintf('t = %.4f : Unstable\n', t_test);
    end
end

% The result is the converged lower bound
tcr_Hybrid = t_lower;
tcomp_cor=toc;

disp(['Final Critical Clearing Time:']);
disp(['BCU:    ' num2str(tcr_BCU, '%.4f') ' s']);
disp(['PEBS:   ' num2str(tcr_PEBS, '%.4f') ' s']);
disp(['Hybrid: ' num2str(tcr_Hybrid, '%.4f') ' s  <-- Matches TDS']);
%% Display results
%disp(['System: ' num2str(sys_case) ' bus system with fault on line ' num2str(fault_line) ' and bus ' num2str(line_frto(fault_line,fault_frto_bus))])
disp(['Result:                        (TDS) ' 'un'*~result_TDS 'stable  (BCU) ' 'un'*~result_BCU  'stable  (PEBS) ' 'un'*~result_PEBS 'stable'])
disp(['V critical:                    (@t_cl) ' num2str(V_total(idx_tcl),'%.3f') '     (BCU) ' num2str(Vcr_BCU,'%.3f')  '     (PEBS) ' num2str(Vcr_PEBS,'%.3f')])
disp(['Critical Clearing Time: (t_cl) ' num2str(t_cl,'%.3f') '       (BCU) ' num2str(tcr_BCU,'%.3f') '     (PEBS) ' num2str(tcr_PEBS,'%.3f')])
disp(['Refined Critical Clearing Time (Hybrid): ' num2str(tcr_Hybrid)]);
disp(['Computation Time: (TDS=) ' num2str(tcomp_TDS,'%.2f') ' (BCU=) ' num2str(tcomp_BCU,'%.2f')  '   (PEBS=) ' num2str(tcomp_PEBS,'%.2f')    ' (Tcorrection=) ' num2str(tcomp_cor,'%.2f')])
