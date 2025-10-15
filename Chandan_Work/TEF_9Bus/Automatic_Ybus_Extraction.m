%% ========================================================================
%  extract_Ybus_for_TEF.m
%
%  Purpose:
%  --------
%  Automatically extracts the Y-bus (admittance matrix) from a
%  Simscape Electrical / SimPowerSystems model (e.g., IEEE 9-bus),
%  and computes the G (conductance) and B (susceptance) matrices.
%
%  Works with models built using "Fundamental Parameterization"
%  synchronous machines (classical model).
%
%  Output:
%    - Ybus, G, B matrices
%    - C and D matrices for TEF: C_ij = Ei*Ej*B_ij, D_ij = Ei*Ej*G_ij
%
%  Reference:
%    P. Bhui and N. Senroy, “Real-Time Prediction and Control of
%    Transient Stability Using TEF,” IEEE TPWRS, 2016.
%
%  ------------------------------------------------------------------------
%% ========================================================================


%% ---- USER INPUTS -------------------------------------------------------

% Name of your Simulink model (without .slx)
modelName = bdroot;   % <-- replace with your model name

% Make sure your model is on MATLAB path
load_system(modelName);

% Optional: ensure model is compiled before calling power_ybus
set_param(modelName,'SimulationCommand','update');

% If using Simscape Electrical (Specialized Power Systems toolbox)
% you can use 'power_ybus' to extract Ybus automatically.
% The model must be loaded in memory for this to work.

fprintf('\nExtracting Y-bus for model: %s\n', modelName);
Ybus_info = power_analyze(modelName,'Y');   % structure with Ybus and bus names
try
    % Newer versions: returns Ybus directly
    Ybus = power_analyze(modelName,'Y');
    busNames = {};   % Not returned by default
    fprintf('Ybus extracted as a direct matrix (no structure fields).\n');
catch
    % Older versions: returns structure
    Ybus_info = power_analyze(modelName,'Y');
    if isstruct(Ybus_info)
        Ybus = Ybus_info.Ybus;
        if isfield(Ybus_info,'BusName')
            busNames = Ybus_info.BusName;
        else
            busNames = {};
        end
        fprintf('Ybus extracted from structure (old-style output).\n');
    else
        error('Unexpected output type from power_analyze.');
    end
end

nbus = length(busNames);

fprintf('Y-bus extracted successfully: %d x %d complex matrix.\n', nbus, nbus);

%% ---- DISPLAY BASIC PROPERTIES -----------------------------------------

fprintf('\nFirst few bus names:\n');
disp(busNames(1:min(5,nbus)));

fprintf('Sample element Y(1,2) = %.5f + j%.5f\n', real(Ybus(1,2)), imag(Ybus(1,2)));

%% ---- COMPUTE REAL AND IMAG PARTS --------------------------------------

G = real(Ybus);
B = imag(Ybus);

fprintf('\nAverage line conductance G magnitude: %.4e\n', mean(abs(G(:))));
fprintf('Average susceptance B magnitude: %.4e\n', mean(abs(B(:))));

%% ---- TEF MATRICES (C and D) -------------------------------------------
% C_ij = Ei * Ej * B_ij
% D_ij = Ei * Ej * G_ij
%
% For demo, assume |E| = 1.0 pu for all generator buses.
% Replace this with measured internal EMFs or PMU voltage magnitudes later.

E = ones(nbus,1);    % placeholder (use PMU/steady-state magnitudes if available)

C = (E * E.') .* B;
D = (E * E.') .* G;

fprintf('\nC and D matrices constructed for TEF framework.\n');

%% ---- SAVE OUTPUTS ------------------------------------------------------

save('Ybus_TEF_data.mat','Ybus','G','B','C','D','busNames','E');

fprintf('\nSaved: Ybus_TEF_data.mat\n');
fprintf('Variables: Ybus, G, B, C, D, busNames, E\n');
fprintf('\nUse these in your TEF computation script.\n');

%% ---- OPTIONAL: Visualization ------------------------------------------
figure;
imagesc(abs(B)); colorbar; title('|B| (Imag part of Ybus)');
xlabel('Bus index'); ylabel('Bus index');

figure;
imagesc(abs(G)); colorbar; title('|G| (Real part of Ybus)');
xlabel('Bus index'); ylabel('Bus index');
