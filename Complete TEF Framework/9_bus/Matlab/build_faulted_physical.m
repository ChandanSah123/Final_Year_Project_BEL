function Yf = build_faulted_physical(Y_full, fault_bus_num)
% BUILD_FAULTED_PHYSICAL Removes the row and column of the faulted bus.
%
% Inputs:
%   Y_full:        NxN complex double matrix (The full Ybus)
%   fault_bus_num: Integer (The bus number to remove)
%
% Output:
%   Yf:            (N-1)x(N-1) complex matrix

    % Create a copy of the input matrix to avoid modifying the original
    Yf = Y_full;

    % 1. Remove the Row
    % MATLAB handles the resizing automatically.
    Yf(fault_bus_num, :) = [];

    % 2. Remove the Column
    % Note: After removing the row, the matrix is (N-1)xN.
    % The column indices remain valid relative to the original width.
    Yf(:, fault_bus_num) = [];
    
end