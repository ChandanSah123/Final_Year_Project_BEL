%%----Load simulink IEEE39bus model and extract paramters----------%%
try
    %open_system('D:\Final year project\Final_Year_Project_BEL[1]\Chandan_Work\IEEE_39Bus\Simulink_model.slx')
    modelName = bdroot;
    
catch
    error("no model is currently working")
end

% Transmission line parameter
lineBlocks=find_system(modelName,'Masktype','Transmission Line (Three-Phase)');
line_data = struct();
for i=1:length(lineBlocks)
    block_path = lineBlocks{i};
    line_data(i).BlockName=get_param(block_path,'Name');
    line_data(i).Length=str2double(get_param(block_path,'length'));
    line_data(i).Resitance=str2double(get_param(block_path,'R'));
    line_data(i).Inductance=str2double(get_param(block_path,'L'));
    line_data(i).Mutual_inductance=str2double(get_param(block_path,'M'));
    line_data(i).CapacitanceCg=str2double(get_param(block_path,'Cg'));
    line_data(i).CapacitanceCl=str2double(get_param(block_path,'Cl'));
    line_data(i).Mutual_resistance=str2double(get_param(block_path,'Rm'));
end

% For Transformer
transBlocks=find_system(modelName,'Masktype','Two-Winding Transformer (Three-Phase)');
trans_data = struct();
for i=1:length(transBlocks)
    block_path1=transBlocks{i};
    trans_data(i).BlockName=get_param(block_path1,'Name');
    trans_data(i).PrimaryWindingResistance = str2double(get_param(block_path1, 'pu_RW1'));
    trans_data(i).SecondaryWindingResistance = str2double(get_param(block_path1, 'pu_RW2'));
    trans_data(i).PrimaryLeakageReactance = eval(get_param(block_path1, 'pu_Xl1'));   % has '/2' form
    trans_data(i).SecondaryLeakageReactance = eval(get_param(block_path1, 'pu_Xl2'));
end

%loading active and reactive power of load
load_block= find_system(modelName, 'MaskType', 'Wye-Connected Load');
load_data=struct();
for i=1:length(load_block)
    block_path2=load_block{i};
    load_data(i).BlockName= get_param(block_path2, 'Name');
    load_data(i).Real_power=str2double(get_param(block_path2,'P'));
    load_data(i).Reactive_power=str2double(get_param(block_path2,'Qpos'));
end

%loading generator resistance, reactances, transient and sub-transeint reactances
generator_block1=find_system(modelName,'MaskType','Synchronous Machine Salient Pole');
generator_block2=find_system(modelName,'MaskType','Synchronous Machine Round Rotor');
generator_data1=struct();
generator_data2=struct();

%for silent rotor type
for i=1:length(generator_block1)
    block_path4=generator_block1{i};
    generator_data1(i).BlockName=get_param(block_path4, 'Name');
    generator_data1(i).stator_resistance=str2double(get_param(block_path4,'Ra'));
    generator_data1(i).leakage_reactance=str2double(get_param(block_path4,'Xl'));
    generator_data1(i).Xd=str2double(get_param(block_path4,'Xd'));
    generator_data1(i).Xq=str2double(get_param(block_path4,'Xq'));
    generator_data1(i).Xdd=str2double(get_param(block_path4,'Xdd'));
    generator_data1(i).Xddd=str2double(get_param(block_path4,'Xddd'));
    generator_data1(i).Xqdd=str2double(get_param(block_path4,'Xqdd'));
end
%for round rotor type
for i=1:length(generator_block2)
    block_path3=generator_block2{i};
    generator_data2(i).BlockName=get_param(block_path3, 'Name');
    generator_data2(i).stator_resistance=str2double(get_param(block_path3,'Ra'));
    generator_data2(i).leakage_reactance=str2double(get_param(block_path3,'Xl'));
    generator_data2(i).Xd=str2double(get_param(block_path3,'Xd'));
    generator_data2(i).Xq=str2double(get_param(block_path3,'Xq'));
    generator_data2(i).Xdd=str2double(get_param(block_path3,'Xdd'));
    generator_data2(i).Xqd=str2double(get_param(block_path3,'Xqd'));
    generator_data2(i).Xddd=str2double(get_param(block_path3,'Xddd'));
    generator_data2(i).Xqdd=str2double(get_param(block_path3,'Xqdd'));

end

%Reading Mechanical Power
mechanical_block=find_system(modelName,'MaskType','Machine Mechanical Power');
mechanical_data=struct();
for i=1:length(mechanical_block)
    block_path5=mechanical_block{i};
    mechanical_data(i).BlockName=get_param(block_path5,'Name');
    mechanical_data(i).ApparantPower=str2double(get_param(block_path5,'SRated'));
    mechanical_data(i).Inetria_H=str2double(get_param(block_path5,'H'));
end
