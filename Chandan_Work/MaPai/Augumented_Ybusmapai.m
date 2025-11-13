
% Y previously entered Ybus matrix

% Define ybar (3x3)
ybar = diag([-1i*16.45, -1i*8.35, -1i*5.52]);
load_ybar=complex(zeros(9,9));
load_ybar(5,5)=(1.26-0.504i);
load_ybar(6,6)=(0.877-0.292i);
load_ybar(8,8)=(0.968-0.339i);

% --- Construct the augmented matrix ---
Y_aug = zeros(12,12);

% Map the existing Y into Y_aug (rows 4:12, cols 4:12)
Y(5,5)=2.553-17.338i;
Y(6,6)=3.224-15.841i;
Y(8,8)=2.772-23.303i;
YN = Y;

% Add ybar and -ybar blocks
Y_aug(1:3, 1:3) = ybar;        % top-left block
Y_aug(1:3, 4:6) = -ybar;       % top-middle block (coupling)
Y_aug(4:6, 1:3) = -ybar;       % middle-left block (transpose coupling)

% Keep others as zero per the paper figure
YNadd=zeros(9,9);
YNadd(1:3,1:3)=ybar;
YN1=YN+YNadd;
YN2=YN1+load_ybar;
Y_aug(4:12, 4:12) = YN2;


%% Display result
%disp('Augmented Admittance Matrix (YN1):');
%disp(Y_aug);
%Reducing to form Yint
YA=ybar;
YB=complex(zeros(3,9));
YB(1:3,1:3)=-ybar;
YC=complex(zeros(9,3));
YC(1:3,1:3)=-ybar;
YD=YN2;
YDC=YD\YC;
Yint=YA-YB*YDC;

