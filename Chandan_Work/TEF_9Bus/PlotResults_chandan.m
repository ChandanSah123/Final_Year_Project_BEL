% Code to plot simulation results from IEEE9BusSyste
% Generate new simulation results if they don't exist or if they need to be updated
if ~exist('simlog_IEEE9BusSystem', 'var') || ...
        simlogNeedsUpdate(simlog_IEEE9BusSystem, 'IEEE9BusSystem')
    sim('IEEE9BusSystem')
end

% Reuse figure if it exists, else create new figure
if ~exist('h1_IEEE9BusSystem', 'var') || ...
        ~isgraphics(h1_IEEE9BusSystem, 'figure')
    h1_IEEE9BusSystem = figure('Name', 'IEEE9BusSystem');
end
figure(h1_IEEE9BusSystem)
clf(h1_IEEE9BusSystem)

% Get simulation results
t = simlog_IEEE9BusSystem.Gen1Bus1_Swing.Rotor_velocity.pu_output.series.time;
w_Gen1Bus1 = simlog_IEEE9BusSystem.Gen1Bus1_Swing.Rotor_velocity.pu_output.series.values;
Vt_Gen1Bus1 = simlog_IEEE9BusSystem.Gen1Bus1_Swing.Terminal_voltage.pu_output.series.values;
theta_Gen1Bus1 = simlog_IEEE9BusSystem.Gen1Bus1_Swing.Rotor_electrical_angle.pu_output.series.values;

w_Gen2Bus2 = simlog_IEEE9BusSystem.Gen2Bus2_PV_1_025_pu_163_MW.Rotor_velocity.pu_output.series.values;
Vt_Gen2Bus2 = simlog_IEEE9BusSystem.Gen2Bus2_PV_1_025_pu_163_MW.Terminal_voltage.pu_output.series.values;
theta_Gen2Bus2 = simlog_IEEE9BusSystem.Gen2Bus2_PV_1_025_pu_163_MW.Rotor_electrical_angle.pu_output.series.values;


w_Gen3Bus3 = simlog_IEEE9BusSystem.Gen3Bus3_PV_1_025_pu_85_MW.Rotor_velocity.pu_output.series.values;
Vt_Gen3Bus3 = simlog_IEEE9BusSystem.Gen3Bus3_PV_1_025_pu_85_MW.Terminal_voltage.pu_output.series.values;
theta_Gen3Bus3 = simlog_IEEE9BusSystem.Gen3Bus3_PV_1_025_pu_85_MW.Rotor_electrical_angle.pu_output.series.values;

delta1=theta_Gen1Bus1-theta_Gen1Bus1;
delta2=theta_Gen2Bus2-theta_Gen1Bus1;
delta3=theta_Gen3Bus3-theta_Gen1Bus1;


% Plot results
% Plot W
figure(1)
handles(1)=subplot(2, 2, 1);
plot(t, w_Gen1Bus1, 'LineWidth', 1)
hold on
plot(t, w_Gen2Bus2, 'LineWidth', 1)
hold on
plot(t, w_Gen3Bus3, 'LineWidth', 1)
hold off
grid on
title('Generator rotor speed and terminal voltage')
ylabel('Rotor speed (p.u.)')
legend({'Gen1@Bus1', 'Gen2@Bus2', 'Gen3@Bus3'},'Location','southwest');

%Plot Vt
handles(2)=subplot(2, 2, 2);
plot(t, Vt_Gen1Bus1, 'LineWidth', 1)
hold on
plot(t, Vt_Gen2Bus2, 'LineWidth', 1)
hold on
plot(t, Vt_Gen3Bus3, 'LineWidth', 1)
hold off
grid on
ylabel('Terminal voltage (p.u.)')
xlabel('Time (s)')
legend({'Gen1@Bus1', 'Gen2@Bus2', 'Gen3@Bus3'},'Location','southwest');

%Plot Theta
handles(3)=subplot(2, 2, 3);
plot(t,theta_Gen1Bus1,'LineWidth',1)
hold on
plot(t, theta_Gen2Bus2, 'LineWidth', 1)
hold on
plot(t, theta_Gen3Bus3, 'LineWidth', 1)
hold off

grid on
ylabel('Rotor electrical angle')
xlabel('Time (s)')
legend({'Gen1@Bus1', 'Gen2@Bus2', 'Gen3@Bus3'},'Location','southwest');

% Delta and W plot together
handle(1)=subplot(2,2,4);
plot(t,delta1,'LineWidth',1)
hold on
plot(t,delta2,'LineWidth',1)
hold on
plot(t,delta3,'LineWidth',1)

grid on
ylabel('delta')
xlabel('Time (s)')
linkaxes(handles, 'x')