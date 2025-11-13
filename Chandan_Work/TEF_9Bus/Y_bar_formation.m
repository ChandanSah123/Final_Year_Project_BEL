try
    modelName = bdroot;
catch
    error("No model is currently open. Please open your IEEE9Bus Simulink model.");
end
load("IEEE9BusSystem_loadflow_Bus_result.mat");
V=busbars(:,'Voltage Magnitude, pu');
Vpf=[V{5,1} V{6,1} V{8,1}];
Sbase=100e6;
load_conn=[
     5 5;
     6 6;
     8 8
         ];
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
ybar=(zeros(3,3));
ybar(1,1)=(1/(1i*generator_data1.Xdd));
ybar(2,2)=(1/(1i*generator_data2(1).Xdd));
ybar(3,3)=(1/(1i*generator_data2(2).Xdd));

%loading active and reactive power of load
load_block= find_system(modelName, 'MaskType', 'Wye-Connected Load');
load_data=struct();
for i=1:length(load_block)
    block_path2=load_block{i};
    load_data(i).BlockName= get_param(block_path2, 'Name');
    load_data(i).Real_power=str2double(get_param(block_path2,'P'));
    load_data(i).Reactive_power=str2double(get_param(block_path2,'Qpos'));
end

Yloadbar=zeros(9,9);
%complex load Y bar.
for i=1:length(load_data)
    P=load_data(i).Real_power;
    Q=load_data(i).Reactive_power;
    Sconj=P-1i*Q;
    Sconjpu=Sconj/Sbase;
    Ylpu=Sconjpu/(Vpf(i)*Vpf(i));
    from=load_conn(i,1);
    to=load_conn(i,2);
    Yloadbar(from,to)=Ylpu;
end

%disp(Yloadbar);