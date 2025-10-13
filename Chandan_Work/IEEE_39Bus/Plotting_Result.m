% Code to plot simulation results from IEEE39BusSystem
if ~exist('simlog_IEEE9BusSystem', 'var') || ...
        simlogNeedsUpdate(simlog_IEEE9BusSystem, 'IEEE9BusSystem')
 sim('Simulink_model')
end
% Reuse figure if it exists, else create new figure
if ~exist('h1_IEEE39BusSystem', 'var') || ...
        ~isgraphics(h1_IEEE39BusSystem, 'figure')
    h1_IEEE39BusSystem = figure('Name', 'IEEE39BusSystem');
end
figure(h1_IEEE39BusSystem)
clf(h1_IEEE39BusSystem)

plotData(simlog_IEEE39BusSystem)

% Plot rotor speeds and terminal voltages of generators
function plotData(simlog)

% Get simulation results
t = simlog.Generators.Gen1Bus39.Rotor_velocity.pu_output.series.time;

w_Gen10Bus30 = simlog.Generators.Gen10Bus30.Rotor_velocity.pu_output.series.values;
Vt_Gen10Bus30 = simlog.Generators.Gen10Bus30.Terminal_voltage.pu_output.series.values;
theta_Gen10Bus30=simlog.Generators.Gen10Bus30.Rotor_electrical_angle.pu_output.series.values;

w_Gen2Bus31 = simlog.Generators.Gen2Bus31.Rotor_velocity.pu_output.series.values;
Vt_Gen2Bus31 = simlog.Generators.Gen2Bus31.Terminal_voltage.pu_output.series.values;
theta_Gen2Bus31=simlog.Generators.Gen2Bus31.Rotor_electrical_angle.pu_output.series.values;

w_Gen3Bus32 = simlog.Generators.Gen3Bus32.Rotor_velocity.pu_output.series.values;
Vt_Gen3Bus32 = simlog.Generators.Gen3Bus32.Terminal_voltage.pu_output.series.values;
theta_Gen3Bus32=simlog.Generators.Gen3Bus32.Rotor_electrical_angle.pu_output.series.values;

w_Gen4Bus33 = simlog.Generators.Gen4Bus33.Rotor_velocity.pu_output.series.values;
Vt_Gen4Bus33 = simlog.Generators.Gen4Bus33.Terminal_voltage.pu_output.series.values;
theta_Gen4Bus33=simlog.Generators.Gen4Bus33.Rotor_electrical_angle.pu_output.series.values;

w_Gen5Bus34 = simlog.Generators.Gen5Bus34.Rotor_velocity.pu_output.series.values;
Vt_Gen5Bus34 = simlog.Generators.Gen5Bus34.Terminal_voltage.pu_output.series.values;
theta_Gen5Bus34=simlog.Generators.Gen5Bus34.Rotor_electrical_angle.pu_output.series.values;


w_Gen6Bus35 = simlog.Generators.Gen6Bus35.Rotor_velocity.pu_output.series.values;
Vt_Gen6Bus35 = simlog.Generators.Gen6Bus35.Terminal_voltage.pu_output.series.values;
theta_Gen6Bus35=simlog.Generators.Gen6Bus35.Rotor_electrical_angle.pu_output.series.values;


w_Gen7Bus36 = simlog.Generators.Gen7Bus36.Rotor_velocity.pu_output.series.values;
Vt_Gen7Bus36 = simlog.Generators.Gen7Bus36.Terminal_voltage.pu_output.series.values;
theta_Gen7Bus36=simlog.Generators.Gen7Bus36.Rotor_electrical_angle.pu_output.series.values;


w_Gen8Bus37 = simlog.Generators.Gen8Bus37.Rotor_velocity.pu_output.series.values;
Vt_Gen8Bus37 = simlog.Generators.Gen8Bus37.Terminal_voltage.pu_output.series.values;
theta_Gen8Bus37=simlog.Generators.Gen8Bus37.Rotor_electrical_angle.pu_output.series.values;


w_Gen9Bus38 = simlog.Generators.Gen9Bus38.Rotor_velocity.pu_output.series.values;
Vt_Gen9Bus38 = simlog.Generators.Gen9Bus38.Terminal_voltage.pu_output.series.values;
theta_Gen9Bus38=simlog.Generators.Gen9Bus38.Rotor_electrical_angle.pu_output.series.values;


w_Gen1Bus39 = simlog.Generators.Gen1Bus39.Rotor_velocity.pu_output.series.values;
Vt_Gen1Bus39 = simlog.Generators.Gen1Bus39.Terminal_voltage.pu_output.series.values;
theta_Gen1Bus39=simlog.Generators.Gen1Bus39.Rotor_electrical_angle.pu_output.series.values;


%  Delta Calculation
Delta1=theta_Gen1Bus39-(2*pi*60*t);
Delta2=theta_Gen2Bus31-(2*pi*60*t);
Delta3=theta_Gen3Bus32-(2*pi*60*t);
Delta4=theta_Gen4Bus33-(2*pi*60*t);
Delta5=theta_Gen5Bus34-(2*pi*60*t);
Delta6=theta_Gen6Bus35-(2*pi*60*t);
Delta7=theta_Gen7Bus36-(2*pi*60*t);
Delta8=theta_Gen8Bus37-(2*pi*60*t);
Delta9=theta_Gen9Bus38-(2*pi*60*t);
Delta10=theta_Gen10Bus30-(2*pi*60*t);

%Relative Delta Calculation
delta1=Delta1-Delta1;
delta2=Delta2-Delta1;
delta3=Delta3-Delta1;
delta4=Delta4-Delta1;
delta5=Delta5-Delta1;
delta6=Delta6-Delta1;
delta7=Delta7-Delta1;
delta8=Delta8-Delta1;
delta9=Delta9-Delta1;
delta10=Delta10-Delta1;

% Plot results
figure(1); 
% W plot
handles(1) = subplot(3, 1, 1);
plot(t, w_Gen1Bus39, 'LineWidth', 1)
hold on
plot(t, w_Gen2Bus31, 'LineWidth', 1)
hold on
plot(t, w_Gen3Bus32, 'LineWidth', 1)
hold on
plot(t, w_Gen4Bus33, 'LineWidth', 1)
hold on
plot(t, w_Gen5Bus34, 'LineWidth', 1)
hold on
plot(t, w_Gen6Bus35, 'LineWidth', 1)
hold on
plot(t, w_Gen7Bus36, 'LineWidth', 1)
hold on
plot(t, w_Gen8Bus37, 'LineWidth', 1)
hold on
plot(t, w_Gen9Bus38, 'LineWidth', 1)
hold on
plot(t, w_Gen10Bus30, 'LineWidth', 1)
hold off

grid on
title('Generator rotor speed,terminal voltage and rotor electrical angle')
ylabel('Rotor speed (p.u.)')
legend({'Gen1Bus39','Gen2@Bus31','Gen3@Bus32','Gen4@Bus33','Gen5@Bus34','Gen6@Bus35','Gen7@Bus36','Gen8@Bus37','Gen9@Bus38','Gen10Bus30'},'Location','southwest');

%Vt plot
handles(2) = subplot(3, 1, 2);
plot(t, Vt_Gen1Bus39, 'LineWidth', 1)
hold on
plot(t, Vt_Gen2Bus31, 'LineWidth', 1)
hold on
plot(t, Vt_Gen3Bus32, 'LineWidth', 1)
hold on
plot(t, Vt_Gen4Bus33, 'LineWidth', 1)
hold on
plot(t, Vt_Gen5Bus34, 'LineWidth', 1)
hold on
plot(t, Vt_Gen6Bus35, 'LineWidth', 1)
hold on
plot(t, Vt_Gen7Bus36, 'LineWidth', 1)
hold on
plot(t, Vt_Gen8Bus37, 'LineWidth', 1)
hold on
plot(t, Vt_Gen9Bus38, 'LineWidth', 1)
hold on
plot(t, Vt_Gen10Bus30, 'LineWidth', 1)
hold off

grid on
ylabel('Terminal voltage (p.u.)')
xlabel('Time (s)')
legend({'Gen1Bus39','Gen2@Bus31','Gen3@Bus32','Gen4@Bus33','Gen5@Bus34','Gen6@Bus35','Gen7@Bus36','Gen8@Bus37','Gen9@Bus38','Gen10@Bus30'},'Location','southwest');

%theta plot
handles(3)=subplot(3,1,3);
plot(t,theta_Gen1Bus39,'LineWidth',1)
hold on
plot(t,theta_Gen2Bus31,'LineWidth',1)
hold on
plot(t,theta_Gen3Bus32,'LineWidth',1)
hold on
plot(t,theta_Gen4Bus33,'LineWidth',1)
hold on
plot(t,theta_Gen5Bus34,'LineWidth',1)
hold on
plot(t,theta_Gen6Bus35,'LineWidth',1)
hold on
plot(t,theta_Gen7Bus36,'LineWidth',1)
hold on
plot(t,theta_Gen8Bus37,'LineWidth',1)
hold on
plot(t,theta_Gen9Bus38,'LineWidth',1)
hold on
plot(t,theta_Gen10Bus30,'LineWidth',1)
hold off

grid on
ylabel('Rotor Electrical angle')
xlabel('Time (s)')
legend({'Gen1Bus39','Gen2@Bus31','Gen3@Bus32','Gen4@Bus33','Gen5@Bus34','Gen6@Bus35','Gen7@Bus36','Gen8@Bus37','Gen9@Bus38','Gen10@Bus30'},'Location','southwest');

% Delta Plot
figure(2); 
handles(1)=subplot(5,2,1);
plot(t,delta1,'LineWidth',1);
grid on
ylabel('delta1')
xlabel('Time (s)')

handles(2)=subplot(5,2,2);
plot(t,delta2,'LineWidth',1);
grid on
ylabel('Delta2')
xlabel('Time (s)')

handles(3)=subplot(5,2,3);
plot(t,delta3,'LineWidth',1);
grid on
ylabel('Delta3')
xlabel('Time (s)')

handles(4)=subplot(5,2,4);
plot(t,delta4,'LineWidth',1);
grid on
ylabel('Delta4')
xlabel('Time (s)')

handles(5)=subplot(5,2,5);
plot(t,delta5,'LineWidth',1);
grid on
ylabel('Delta5')
xlabel('Time (s)')

handles(6)=subplot(5,2,6);
plot(t,delta6,'LineWidth',1);
grid on
ylabel('Delta6')
xlabel('Time (s)')

handles(7)=subplot(5,2,7);
plot(t,delta7,'LineWidth',1);
grid on
ylabel('Delta7')
xlabel('Time (s)')

handles(8)=subplot(5,2,8);
plot(t,delta8,'LineWidth',1);
grid on
ylabel('Delta8')
xlabel('Time (s)')

handles(9)=subplot(5,2,9);
plot(t,delta9,'LineWidth',1);
grid on
ylabel('Delta9')
xlabel('Time (s)')

handles(10)=subplot(5,2,10);
plot(t,delta10,'LineWidth',1);
grid on
ylabel('Delta10')
xlabel('Time (s)')

% Plot of W of each generator
figure(3); 
handles(1)=subplot(5,2,1);
plot(t, w_Gen1Bus39,'LineWidth',1);
grid on
ylabel(' w1')
xlabel('Time (s)')

handles(2)=subplot(5,2,2);
plot(t, w_Gen2Bus31,'LineWidth',1);
grid on
ylabel(' w2')
xlabel('Time (s)')

handles(3)=subplot(5,2,3);
plot(t, w_Gen3Bus32,'LineWidth',1);
grid on
ylabel(' w3')
xlabel('Time (s)')

handles(4)=subplot(5,2,4);
plot(t, w_Gen4Bus33,'LineWidth',1);
grid on
ylabel(' w4')
xlabel('Time (s)')

handles(5)=subplot(5,2,5);
plot(t, w_Gen5Bus34,'LineWidth',1);
grid on
ylabel(' w5')
xlabel('Time (s)')

handles(6)=subplot(5,2,6);
plot(t, w_Gen6Bus35,'LineWidth',1);
grid on
ylabel(' w5')
xlabel('Time (s)')

handles(7)=subplot(5,2,7);
plot(t, w_Gen7Bus36,'LineWidth',1);
grid on
ylabel(' w7')
xlabel('Time (s)')

handles(8)=subplot(5,2,8);
plot(t, w_Gen8Bus37,'LineWidth',1);
grid on
ylabel(' w8')
xlabel('Time (s)')

handles(9)=subplot(5,2,9);
plot(t, w_Gen9Bus38,'LineWidth',1);
grid on
ylabel(' w9')
xlabel('Time (s)')

handles(10)=subplot(5,2,10);
plot(t, w_Gen10Bus30,'LineWidth',1);
grid on
ylabel(' w10')
xlabel('Time (s)')

% Rotor Angle delta and W on The same plot
%delta plot
figure(4)
subplot(2,1,1)
plot(t,delta1,'LineWidth',1)
hold on
plot(t,delta2,'LineWidth',1)
hold on
plot(t,delta3,'LineWidth',1)
hold on
plot(t,delta4,'LineWidth',1)
hold on
plot(t,delta5,'LineWidth',1)
hold on
plot(t,delta6,'LineWidth',1)
hold on
plot(t,delta7,'LineWidth',1)
hold on
plot(t,delta8,'LineWidth',1)
hold on
plot(t,delta9,'LineWidth',1)
hold on
plot(t,delta10,'LineWidth',1)
hold off

grid on
ylabel('delta(rad)')
xlabel('Time (s)')
legend({'Gen1Bus39','Gen2@Bus31','Gen3@Bus32','Gen4@Bus33','Gen5@Bus34','Gen6@Bus35','Gen7@Bus36','Gen8@Bus37','Gen9@Bus38','Gen10@Bus30'},'Location','southwest');

% W plot
subplot(2,1,2)
plot(t, w_Gen1Bus39,'LineWidth',1);
hold on
plot(t, w_Gen2Bus31,'LineWidth',1);
hold on
plot(t, w_Gen3Bus32,'LineWidth',1);
hold on
plot(t, w_Gen4Bus33,'LineWidth',1);
hold on
plot(t, w_Gen5Bus34,'LineWidth',1);
hold on
plot(t, w_Gen6Bus35,'LineWidth',1);
hold on
plot(t, w_Gen7Bus36,'LineWidth',1);
hold on
plot(t, w_Gen8Bus37,'LineWidth',1);
hold on
plot(t, w_Gen9Bus38,'LineWidth',1);
hold on
plot(t, w_Gen10Bus30,'LineWidth',1);
hold off
grid on
ylabel('omega(w) P.U')
xlabel('Time (s)')
legend({'Gen1Bus39','Gen2@Bus31','Gen3@Bus32','Gen4@Bus33','Gen5@Bus34','Gen6@Bus35','Gen7@Bus36','Gen8@Bus37','Gen9@Bus38','Gen10@Bus30'},'Location','southwest');


linkaxes(handles, 'x')

end