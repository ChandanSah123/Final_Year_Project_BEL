%% --- 2. Load Simulink Model and Extract Parameters ------------------- %%
try
    modelName = bdroot;
catch
    error("No model is currently open. Please open your IEEE9Bus Simulink model.");
end


% ----- Transmission Lines -----
lineBlocks = find_system(modelName, 'MaskType', 'Transmission Line (Three-Phase)');
line_data = struct();

for i = 1:length(lineBlocks)
    block_path = lineBlocks{i};
    line_data(i).BlockName = get_param(block_path, 'Name');
    line_data(i).Length = str2double(get_param(block_path, 'length')); 
    line_data(i).Resistance = str2double(get_param(block_path, 'r'));
    line_data(i).Inductance = str2double(get_param(block_path, 'l'));
    line_data(i).Mutual_inductance = str2double(get_param(block_path, 'M'));
    line_data(i).CapacitanceCg = str2double(get_param(block_path, 'Cg'));
    line_data(i).CapacitanceCL = str2double(get_param(block_path, 'CL'));
    line_data(i).Mutual_Resistance = str2double(get_param(block_path, 'Rm'));
end

% ----- Transformers -----
transBlocks = find_system(modelName, 'MaskType', 'Two-Winding Transformer (Three-Phase)');
trans_data = struct();

for i = 1:length(transBlocks)
    block_path1 = transBlocks{i};
    trans_data(i).PrimaryWindingResistance = str2double(get_param(block_path1, 'pu_RW1'));
    trans_data(i).SecondaryWindingResistance = str2double(get_param(block_path1, 'pu_RW2'));
    trans_data(i).PrimaryLeakageReactance = eval(get_param(block_path1, 'pu_Xl1'));   % has '/2' form
    trans_data(i).SecondaryLeakageReactance = eval(get_param(block_path1, 'pu_Xl2'));
end

