Y = [ ...
    -1i*17.361     0             0            1i*17.361      0              0              0              0              0;
     0            -1i*16         0            0              0              0           1i*16                0              0;
     0             0            -1i*17.065    0              0              0               0              0              1i*17.065;
     1i*17.361     0             0            3.307 - 1i*39.309  -1.365 + 1i*11.604  -1.942 + 1i*10.511  0  0  0;
     0             0             0           -1.365 + 1i*11.604   2.553 - 1i*17.338     0           -1.188 + 1i*5.975   0  0;
     0             0             0           -1.942 + 1i*10.511    0    3.224 - 1i*15.841   0   0   -1.282 + 1i*5.588;
     0             1i*16         0            0          -1.188 + 1i*5.975   0   2.805 - 1i*35.446  -1.617 + 1i*13.698   0;
     0             0             0            0          0          0   -1.617 + 1i*13.698   2.772 - 1i*23.303  -1.155 + 1i*9.784;
     0             0             1i*17.065    0          0   -1.282 + 1i*5.588   0   -1.155 + 1i*9.784   2.437 - 1i*32.154
];

% Create readable MATLAB table
Y_table = array2table(Y, ...
    'VariableNames', compose('Bus%d', 1:9), ...
    'RowNames', compose('Bus%d', 1:9));

% Export to Excel file
writetable(Y_table, 'Ybus_matrix.xlsx', 'WriteRowNames', true);

%Optional: display to check
disp(Y_table)
save('Y.mat','Y');