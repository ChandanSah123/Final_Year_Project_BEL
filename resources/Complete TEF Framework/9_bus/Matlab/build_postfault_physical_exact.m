function Yp = build_postfault_physical_exact(Y_full, from_bus, To_bus, y_line, b_charging)
    % MATLAB arrays are passed by value when modified, so Yp starts as a copy.
    Yp = Y_full;
    idx_from = from_bus;
    idx_To = To_bus;
    
    % 1. Subtract Y_series and B_charging from diagonals
    % (Removing the shunt and series contributions of the line from the nodes)
    y_total = y_line + b_charging;
    
    Yp(idx_from, idx_from) = Yp(idx_from, idx_from) - y_total;
    Yp(idx_To, idx_To) = Yp(idx_To, idx_To) - y_total;
    Yp(idx_from, idx_To) = 0;
    Yp(idx_To, idx_from) = 0;

end