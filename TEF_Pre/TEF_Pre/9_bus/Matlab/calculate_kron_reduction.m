function Yint = calculate_kron_reduction(Ybus, xd_prime, gen_buses, load_admittance_data)
    % 1. Setup Sizes and Constants
    [N, ~] = size(Ybus);      % Number of buses in the network
    K = length(xd_prime);     % Number of generators
    ybar_diag_vals = 1 ./ (1j * xd_prime);
    ybar = diag(ybar_diag_vals);

    YD = Ybus;

    if isa(load_admittance_data, 'containers.Map')
        % If input is a Map (Dictionary)
        load_buses = cell2mat(keys(load_admittance_data));
        for i = 1:length(load_buses)
            bus_idx = load_buses(i);
            val = load_admittance_data(bus_idx);
            YD(bus_idx, bus_idx) = YD(bus_idx, bus_idx) + val;
        end
    else
        % If input is a Matrix [Bus_Idx, Value]
        for i = 1:size(load_admittance_data, 1)
            bus_idx = real(load_admittance_data(i, 1));
            val = load_admittance_data(i, 2);
            YD(bus_idx, bus_idx) = YD(bus_idx, bus_idx) + val;
        end
    end

    % B. Add Generator Internal Admittances to Network Diagonals
    % (Corresponds to: YN1 = YN + YNadd)
    for i = 1:K
        bus_idx = gen_buses(i);
        YD(bus_idx, bus_idx) = YD(bus_idx, bus_idx) + ybar(i, i);
    end

    % 4. Build YA (Internal Nodes Diagonal)
    YA = ybar;

    % 5. Build YB (Coupling Matrix: Generators <-> Network)
    % Size: K x N
    % Logic: Places -ybar at the connection point
    YB = complex(zeros(K, N));
    
    for i = 1:K
        bus_idx = gen_buses(i);
        % Map the internal gen node (row i) to the physical bus (col bus_idx)
        YB(i, bus_idx) = -ybar(i, i);
    end

    % 6. Build YC (Transpose of YB)
    YC = YB.';

    % 7. Perform Kron Reduction
    % Formula: Yint = YA - YB * inv(YD) * YC
    % Using MATLAB backslash operator for efficiency (YD \ YC is inv(YD)*YC)
    
    Yint = YA - YB * (YD \ YC);

end