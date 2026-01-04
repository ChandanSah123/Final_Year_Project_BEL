clear;
global Yint_post;
load('data1.mat'); 
load('Y_all.mat');
Y1=Yint_post;
g=3;
E=[1.101 1.095 1.125]';
 %E=[1.054 1.050 1.017];
 C=zeros(g,g);
 D=zeros(g,g);
 for i=1:g
    for j=1:g
         C(i,j)=E(i)*E(j)*imag(Y1(i,j));
        D(i,j)=E(i)*E(j)*real(Y1(i,j)); 
    end
 end
Pi=[-0.642497899332592 1.19790121698815 0.527705756316238]';
th_vec=[-0.821804214185591  1.80954146737929   1.13324638419523]';
ths=[-0.171021845218841  0.518621972296451  0.240457075839241]';
Vcr = Calculate_PE_single_point(th_vec, ths, Pi, C, D, g);
disp(Vcr)