
S_base = 100e6;      % 100 MVA base
V_base = 230e3;% 230 kV base (line-to-line)
Z_base = (V_base^2) / S_base;
f = 60;              % frequency (Hz)
omega = 2*pi*f;
% -------------- 1) Bus and connection definitions ----------------- %
nBus = 9;

% Transmission lines: order must match line_data order (6 entries)
line_conn = [
    4 5;
    4 6;
    5 7;
    6 9;
    7 8;
    8 9;
];

%primary (230 kV bus))
trans_conn = [
    1 4;  % Gen1 (16.5kV) secondary <-> primary bus 4 (230kV)
    2 7;  % Gen2 (18kV)   secondary <-> primary bus 7 (230kV)
    3 9   % Gen3 (13.8kV) secondary <-> primary bus 9 (230kV)
];

% Per-transformer secondary nominal line-to-line voltages (match order above)
V_sec = [16.5e3, 18e3, 13.8e3];  % [V] for buses [1,2,3] respectively
V_pri = 230e3;                    % primary line-to-line (230 kV)

% Initialize Ybus
Ybus = complex(zeros(nBus,nBus));

% -------------- 2) Process transmission lines (convert to pu & stamp) -- %
for k = 1:length(line_data)
    % Extract values (assumed numeric)
    len_km = line_data(k).Length;            % km
    R_ohm_per_km = line_data(k).Resistance;  % ohm/km
    L_h_per_km = line_data(k).Inductance*1e-3;    % H/km
    M_h_per_km = line_data(k).Mutual_inductance*1e-3;  % H/km (mutual)
    Rm_ohm_per_km = line_data(k).Mutual_Resistance; % if any (ohm/km)
    Cg_f_per_km = line_data(k).CapacitanceCg*1e-6; % F/km (to ground)
    CL_f_per_km = line_data(k).CapacitanceCL*1e-6; % F/km (mutual or line-line total)

    % positive-sequence approximation: series reactance uses (L - M)
    % total series per-phase impedance (ohms)
    R_total = (R_ohm_per_km) * len_km;        % ohms
    X_total = omega * (L_h_per_km - M_h_per_km) * len_km; % ohms (positive-seq approx)
    Z_line_ohm = R_total + 1j*X_total;

    % Convert series impedance to per-unit on system base (230kV, S_base)
    Zpu_line = Z_line_ohm / Z_base;   % pu
    Ypu_series = 1 / Zpu_line;        % pu admittance (series)

    % Shunt capacitance: total B (siemens) for the whole line
    C_total = (Cg_f_per_km + CL_f_per_km) * len_km;    % F
    B_total_S = omega * C_total;                       % S (imag admittance magnitude)
    % per-unit shunt admittance (pure imaginary) on system base:
    Ysh_total_pu = 1j * B_total_S * Z_base;  % pu
    Ysh_end_pu = Ysh_total_pu / 2;           % split half at each end (pi-model)

    % stamp into Ybus
    from = line_conn(k,1);
    to   = line_conn(k,2);

    % Off-diagonal (mutual) entry for series branch
    Ybus(from,to) = Ybus(from,to) - Ypu_series;
    Ybus(to,from) = Ybus(from,to);

    % Diagonals add series admittance + half shunt
    Ybus(from,from) = Ybus(from,from) + Ypu_series + Ysh_end_pu;
    Ybus(to,to)     = Ybus(to,to)     + Ypu_series + Ysh_end_pu;
end


for k = 1:length(trans_data)
    sec_bus = trans_conn(k,1);  % secondary bus (generator side)
    pri_bus = trans_conn(k,2);  % primary bus (230 kV side)
    Vsec = V_sec(k);
    Vpri = V_pri;

    % Extract per-winding pu (may be zero or small)
    Rw1 = trans_data(k).PrimaryWindingResistance;   % pu (primary winding)
    Xl1 = trans_data(k).PrimaryLeakageReactance;    % pu (primary)
    Rw2 = trans_data(k).SecondaryWindingResistance; % pu (secondary)
    Xl2 = trans_data(k).SecondaryLeakageReactance;  % pu (secondary)

    % Build complex pu impedances on transformer base (per-winding)
    Zpu_w1 = complex(Rw1, Xl1);   % pu on transformer base (primary winding)
    Zpu_w2 = complex(Rw2, Xl2);   % pu on transformer base (secondary winding)

    % Total series pu on primary side (pu on same base as primary)
    Zpu_total = Zpu_w1 + Zpu_w2 ;

    % Convert series pu into series admittance (pu)
    Yt_pu = 1 / Zpu_total;

    % Stamp transformer as a series branch between pri_bus and sec_bus
    Ybus(pri_bus, pri_bus) = Ybus(pri_bus, pri_bus) + Yt_pu;
    Ybus(sec_bus, sec_bus) = Ybus(sec_bus, sec_bus) + Yt_pu;
    Ybus(pri_bus, sec_bus) = Ybus(pri_bus, sec_bus) - Yt_pu;
    Ybus(sec_bus, pri_bus) = Ybus(sec_bus, pri_bus) - Yt_pu;
end
%disp(Ybus)