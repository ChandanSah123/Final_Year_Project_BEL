% read_system_data.m
% Minimal loader for IEEE-9 or other systems. Should return a struct 'sys' with network & machine params.
function sys = read_system_data(caseName)
% TODO: Replace with a robust loader: MATPOWER case, PST data, or hand-coded data.
% Required fields in returned sys:
% sys.Ybus - full complex admittance matrix (nb x nb)
% sys.bus  - bus data table (nb x ?)
% sys.line - line data table
% sys.gen  - struct array with fields: bus, M, D, E, Xdprime, Pm, H, inertia
% sys.baseMVA

switch lower(caseName)
    case 'ieee9'
        % Lightweight IEEE-9 classical data (placeholders) - replace with real data
        % For demonstration, provide a minimal structure with 3 gen buses.
        [Ybus, bus, line] = load_ieee9_admittance(); % implement or load MAT file
        gen(1).bus = 1; gen(1).M = 10; gen(1).D = 0; gen(1).E = 1.1; gen(1).Pm = 0.716;
        gen(2).bus = 2; gen(2).M = 8; gen(2).D = 0; gen(2).E = 1.05; gen(2).Pm = 1.63;
        gen(3).bus = 3; gen(3).M = 6; gen(3).D = 0; gen(3).E = 1.01; gen(3).Pm = 0.85;
        sys.Ybus = Ybus;
        sys.bus = bus;
        sys.line = line;
        sys.gen = gen;
        sys.baseMVA = 100;
    otherwise
        error('Unknown caseName. Implement loader or use MATPOWER.');
end
end
