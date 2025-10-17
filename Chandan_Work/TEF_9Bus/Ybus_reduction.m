Y=Ybus;
gen_buses  = [1 2 3];
load_buses = [4 5 6 7 8 9];

% === Step 3: Partition the Y matrix ===
Ygg = Y(gen_buses, gen_buses);
Ygl = Y(gen_buses, load_buses);
Ylg = Y(load_buses, gen_buses);
Yll = Y(load_buses, load_buses);

% === Step 4: Kron Reduction ===
Y_red = Ygg - Ygl * (Yll \ Ylg);
disp(Y_red);

G = real(Y_red);
B = imag(Y_red);
disp(B);
disp(G);