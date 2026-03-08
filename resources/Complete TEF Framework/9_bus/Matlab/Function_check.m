load('Y_all.mat')
xd_prime = [0.0608, 0.1198, 0.1813];
gen_buses = [1, 2, 3];

% Prepare Load Admittance (Matrix Format: [Bus, ComplexValue])
load_adm_matrix = [
    5, (1.26 - 0.504i);
    6, (0.877 - 0.292i);
    8, (0.969 - 0.339i)
];

% --- 2. Call the Function ---
Yint_Post_fault = calculate_kron_reduction(Y_post, xd_prime, gen_buses, load_adm_matrix);

% --- 3. Display Result ---
disp('Reduced Internal Y-Matrix (Post Fault):');
disp(Yint_Post_fault);


% 1. Generator Data
gen_buses = [30, 31, 32, 33, 34, 35, 36, 37, 38, 39];
xd_prime  = [0.025, 0.05, 0.045, 0.035, 0.089, 0.04, 0.044, 0.045, 0.045, 0.004];

% 2. Load Admittance Data
% I converted your Python dictionary to a MATLAB Matrix: [Bus_Index, Complex_Y]
% Note: I excluded the '0' values for cleanliness, but including them is fine too.
load_adm_matrix = [
    3,  3.034 - 0.0226i;
    4,  4.9612 - 1.8257i;
    7,  2.3521 - 0.8451i;
    8,  5.262 - 1.7742i;
    15, 3.1037 - 1.4839i;
    16, 3.0903 - 0.3034i;
    18, 1.4867 - 0.2823i;
    20, 6.392 - 1.0484i;
    21, 2.5737 - 1.0802i;
    23, 2.2673 - 0.775i;
    24, 2.8681 + 0.855i;
    25, 2.0027 - 0.422i;
    26, 1.2557 - 0.1536i;
    27, 2.6095 - 0.7011i;
    28, 1.8681 - 0.2503i;
    29, 2.5719 - 0.244i;
    31, 0.0838 - 0.0419i;   % Load on Generator Bus
    39, 10.4063 - 2.3565i   % Load on Generator Bus
];

% 3. Assume Y_pre is your 39x39 Pre-fault Ybus Matrix
% (You must load this from your data source)
% Y_pre = ... 

% 4. Call the Function
Yint_39 = calculate_kron_reduction(Y_pre, xd_prime, gen_buses, load_adm_matrix);

% 5. Check Result
% The Output Yint_39 will be a 10x10 Matrix
disp('Reduced Yint Dimensions:');
disp(size(Yint_39)); % Should be 10 10