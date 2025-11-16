%%................Ybus-formation...................

S_base = 1e+8;       % 100MVA base power
V_base = 345e3;      % 345KV  base voltage
Z_base = (V_base^2)/S_base;

f= 60;
omega = 2*pi*f;

%----------1) Bus and Connection Defination--------------------%
nBus = 39;

% Transmission lines: order must match line_data order (6 entries)
line_conn = [
1    2
10  11
10  13
11  6
13  14
14  15
15  16
16  17
16  24
17  27
18  17
19  16
2  25
21  16
21  22
23  22
23  24
25  26
26  28
26  29
27  26
28  29
3  18
3  2
39  1
4  14
4  3
5  4
6  5
7  6
8  5
8  7
9  39
9  8
];

%primary (345 kV bus))
trans_conn = [
11  12
13  12
2  30
20  19
22  35
25  37
31  6
32  10
33  19
34  20
36  23
38  29
];

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
    pri_bus = trans_conn(k,2);  % primary bus
   

    % Extract per-winding pu (may be zero or small)
    Rw1 = trans_data(k).PrimaryWindingResistance;   % pu (primary winding)
    Xl1 = trans_data(k).PrimaryLeakageReactance;    % pu (primary)
    Rw2 = trans_data(k).SecondaryWindingResistance; % pu (secondary)
    Xl2 = trans_data(k).SecondaryLeakageReactance;  % pu (secondary)
    Vp =trans_data(k).PrimaryRatedVoltageVp;        % KV (primary)
    Vs =trans_data(k).SecondaryRatedVoltageVs;      % KV (secondary)
    Vsec = Vs*1000;   % in Volt
    Vpri = Vp*1000;

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
%disp(Ybus);


