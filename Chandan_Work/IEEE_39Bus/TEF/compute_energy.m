% compute_energy.m
% Compute total transient energy V = VKE + VPE (structure-preserving individual-machine style)
% theta, omega: vectors in COI or machine angle coordinates (consistent with theta_s)
function V = compute_energy(theta, omega, theta_s, sys)
% VKE
VKE = 0;
for i=1:numel(sys.gen)
    M = sys.gen(i).M;
    VKE = VKE + 0.5 * M * (omega(i))^2;
end
% VPE: implement Bhui / Sauer expression
VPE = compute_VPE(theta, theta_s, sys);
V = VKE + VPE;
end

% compute_VPE: potential energy relative to SEP using machine EMFs and network B/G
function VPE = compute_VPE(theta, theta_s, sys)
n = numel(sys.gen);
% mapping from generator index to bus index
theta_ij = @(i,j) theta(i) - theta(j);
theta_s_ij = @(i,j) theta_s(i) - theta_s(j);
VPE = 0;
% Precompute machine parameters C_ij = E_i*E_j*B_ij, D_ij = E_i*E_j*G_ij
% TODO: Build B_ij and G_ij from sys.Ybus reduced to internal generator nodes
for i=1:n
    Pi = sys.gen(i).Pm; % approx P_m - E^2 Gii included if available
    VPE = VPE + Pi*(theta(i)-theta_s(i));
end
% pairwise coupling (placeholder)
for i=1:n-1
    for j=i+1:n
        % TODO: compute C_ij and D_ij accurately from network and EMFs
        Cij = 0; Dij = 0;
        VPE = VPE + Cij*(cos(theta_ij(i,j)) - cos(theta_s_ij(i,j))) + ...
            Dij*(theta_ij(i,j) - theta_s_ij(i,j) - (sin(theta_ij(i,j)) - sin(theta_s_ij(i,j))));
    end
end
end
