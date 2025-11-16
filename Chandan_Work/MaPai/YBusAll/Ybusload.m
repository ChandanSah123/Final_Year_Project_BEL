
load_ybar=complex(zeros(9,9));
load_ybar(5,5)=(1.26-0.504i);
load_ybar(6,6)=(0.877-0.292i);
load_ybar(8,8)=(0.969-0.339i);

Yload=Y+load_ybar;

gen_buses  = [1 2 3];
load_buses = [4 5 6 7 8 9];

% === Step 3: Partition the Y matrix ===
Ygg = Yload(gen_buses, gen_buses);
Ygl = Yload(gen_buses, load_buses);
Ylg = Yload(load_buses, gen_buses);
Yll = Yload(load_buses, load_buses);

% === Step 4: Kron Reduction ===
Y_redL = Ygg - Ygl * (Yll \ Ylg);
disp(Y_redL);