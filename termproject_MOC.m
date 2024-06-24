% Specify filename
filename = 'C:\Users\hp\Downloads\samplefile4.txt';

% Read data into a matrix
data = readmatrix(filename)

%-------------------------------------------------------------------------%

% Mtrl. properites
Ef = data(1,1)*10.^-3 % Elastic Modulus of Fiber in GPa
nu_f = data(1,2) % poisson's ratio of Fiber
Em = data(2,1)*10.^-3 % Elastic Modulus of Matrix in GPa
nu_m = data(2,2) % poisson's ratio of Matrix
Gf = Ef/(2*(1+nu_f)) % shear Modulus of Fiber in GPa
Gm = Em/(2*(1+nu_m)) % shear Modulus of Matrix in GPa

% (directions 1 and 2 refer to the principal material directions.)

Xt = data(4,1) % Strength of the lamina in tension in the 1-direction in MPa
Xc = -data(4,2) % Strength of the lamina in compression in the 1-direction in MPa
Yt = data(4,3) % Strength of the lamina in tension in the 2-direction in MPa
Yc = -data(4,4) % Strength of the lamina in compression in the 2-direction in MPa
S = data(4,5) % Strength of the lamina in shear in MPa

vf = data(3,1) % volume fraction of Fiber

%-------------------------------------------------------------------------%

% Elastic Modulus of the Composite in 1-direction (using ROM)
E1 = Ef*vf + Em*(1-vf)

% Elastic Modulus of the Composite in 1-direction (using IROM)
E2 = (vf/Ef + (1-vf)/Em).^-1

% Shear Modulus of the composite
G12 = (vf/Gf + (1-vf)/Gm).^-1

% poissons ratio of the composite
nu12 = nu_f*vf + nu_m*(1-vf)

% Define vf range and number of points
vf_ = 0:0.05:1;

E1_ = Ef*vf_ + Em*(1-vf_);
E2_ = (vf_/Ef + (1-vf_)/Em).^-1;
G12_ = (vf_/Gf + (1-vf_)/Gm).^-1;
nu12_ = nu_f*vf_ + nu_m*(1-vf_);

figure

subplot(2,2,1)
plot(vf_,E1_,'r')
xlabel('vf')
ylabel('E1 (in GPa)')
title("E1 vs vf")

subplot(2,2,2)
plot(vf_,E2_,'b')
xlabel('vf')
ylabel('E2 (in GPa)')
title("E2 vs vf")

subplot(2,2,3)
plot(vf_,G12_,'g')
xlabel('vf')
ylabel('G12 (in GPa)')
title("G12 vs vf")

subplot(2,2,4)
plot(vf_,nu12_,'m')
xlabel('vf')
ylabel('nu12 (in GPa)')
title("nu12 vs vf")

hold off

%-------------------------------------------------------------------------%

% Since thickness is quite small, therefore plane-stress condition is considered.
% Now for plane stress condition stress components sigma3, sigma4, sigma5 = 0.
% In the principal material coordinate system (1-2), the components of normal and shear stress are uncoupled.
% Thus the Q matrix will be having only Q11, Q12, Q22, Q66 elements present.

Q = [E1/(1-(nu12.^2)*E2/E1) nu12.*E2/(1-(nu12.^2).*E2/E1) 0; 
     nu12.*E2/(1-(nu12.^2)*E2/E1) E2/(1-(nu12.^2).*E2/E1) 0;
     0 0 G12]
Q_reduced = [Q(1,1), Q(1,2), Q(2,2), Q(3,3)]

% Now we need to transform Q mstrix for any arbitary direction say X-Y system
% Let's denote it as Qt
theta = [0 0 0 0];
numlayers = length(theta);
base_string = 'Qt_';

for i = 1:numlayers
    theta_rad = deg2rad(theta(i)); % Converts theta to radians
    filename = strcat(base_string, num2str(i)); % Concatenates index number as suffix 
    m = cos(theta_rad); n = sin(theta_rad);
    Qxy = Q_transformed(Q, m, n);
    disp(filename)
    disp(Qxy)
end

% Define theta range to find variation of Qx-y wrt theta
theta_range = -90:30:90 % Adjust the range and step size as needed
num_pts = length(theta_range)

% Initialize empty arrays for elements of Qxy matrix
Qxx = zeros(1, num_pts);
Qxy = zeros(1, num_pts);
Qyy = zeros(1, num_pts);
Qxs = zeros(1, num_pts);
Qys = zeros(1, num_pts);
Qss = zeros(1, num_pts);

for i = 1:num_pts
    theta_rad = deg2rad(theta_range(i)); % Converts theta to radians
    m = cos(theta_rad); n = sin(theta_rad);
    Qxy = Q_transformed(Q, m, n);

    Qxx(i) = Qxy(1,1);
    Qxx(1,i) = Qxx(i);

    Qxy(i) = Qxy(1,2);
    Qxy(1,i) = Qxy(i);

    Qyy(i) = Qxy(2,2);
    Qyy(1,i) = Qyy(i);

    Qxs(i) = Qxy(1,3);
    Qxs(1,i) = Qxs(i);

    Qys(i) = Qxy(2,3);
    Qys(1,i) = Qys(i);

    Qss(i) = Qxy(3,3);
    Qss(1,i) = Qss(i);
end

figure 

plot(theta_range,Qxx,'r')
hold on

plot(theta_range,Qxy,'b')
hold on

plot(theta_range,Qxs,'g')
hold on

plot(theta_range,Qyy,'y')
hold on

plot(theta_range,Qys,'c')
hold on

plot(theta_range,Qss,'m')

title('Q(X-Y) components vs theta')
xlabel('theta in degree')
ylabel('in GPa')
legend('Qxx', 'Qxy', 'Qxs', 'Qyy', 'Qys', 'Qss')
hold off

%-------------------------------------------------------------------------%

% Formation of A, B, D matrices to calculate ply stresses in the Laminate

t_l = data(6,1); % thickness of the Laminate in mm
t = t_l/3; % thickness of each in mm

% We previously defined theta in a form of 1D array as = [0, 15, 45, 90]

% Defining z as peripendicular distance from mid-plane (downward from mid-plane a +ve)
z = [-0.6 -0.3 0.0 0.3 0.6];
length_z = length(z);

% Initializing A, B, D as zero matrix
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

% The Laminate is neither symmetric nor balanced, so B matrix will exist
for i = 1:numlayers
    theta_rad = deg2rad(theta(i)); % Converts theta to radians
    m = cos(theta_rad); n = sin(theta_rad);
    Qxy = Q_transformed(Q, m, n);
    A = A + Qxy*(z(i+1)- z(i));
    B = B + 0.5*Qxy*(z(i+1).^2- z(i).^2);
    D = D + 0.334*Qxy*(z(i+1).^3- z(i).^3);
end
A % in GPa
B % in GPa
D % in GPa

% The Global Laminate stiffness matrix or ABBD
G = [A B; B D] % in GPa

% Defining force vector and moment matrix (unit width)
N = [100, 100, 0]'; % in N/mm 
M = [50, 25, 12]'; % N-m/m

% Evaluting strain vector of mid-plane {ephsilon}
eph = [inv(G)*[N; M]]*10.^-3
eph0 = eph(1:3,1) % mid-plane strain matrix
kap = eph(4:end,1) % mid-plane curvature matrix

sigmaxx = zeros(1,5);
sigmaxy = zeros(1,5);
sigmayy = zeros(1,5);

for i = 1:numlayers
    for j = i:i+1
    theta_rad = deg2rad(theta(i));
    m = cos(theta_rad); n = sin(theta_rad);
    Qxy = Q_transformed(Q, m, n);
    sigmaxx(j) = [Qxy(1,:)*eph0 + Qxy(1,:)*kap*z(j)]*10.^3;
    sigmayy(j) = [Qxy(2,:)*eph0 + Qxy(2,:)*kap*z(j)]*10.^3;
    sigmaxy(j) = [Qxy(3,:)*eph0 + Qxy(3,:)*kap*z(j)]*10.^3;
    end
end

figure

plot(z,sigmaxx,'r')
hold on

plot(z,sigmaxy,'b')
hold on

plot(z,sigmayy,'g')

title('sigma(X-Y) vs z')
xlabel('z in mm')
ylabel('in GPa')
legend('sigmaxx', 'sigmaxy', 'sigmayy')
hold off

sigmaxx_avg = zeros(1,numlayers);
sigmaxy_avg = zeros(1,numlayers);
sigmayy_avg = zeros(1,numlayers);

for i = 1:numlayers
    sigmaxx_avg(i) = (sigmaxx(i)+sigmaxx(i+1))/2;
    sigmaxy_avg(i) = (sigmaxy(i)+sigmaxy(i+1))/2;
    sigmayy_avg(i) = (sigmayy(i)+sigmayy(i+1))/2;
end

sigmaxx_avg % taken average 
sigmaxy_avg % taken average
sigmayy_avg % taken average

sigma1 = zeros(1,numlayers);
sigma2 = zeros(1,numlayers);
sigma12 = zeros(1,numlayers);

for i = 1:4
    theta_rad = deg2rad(theta(i));
    m = cos(theta_rad); n = sin(theta_rad);

    sigma1(i) = [m.^2, n.^2, 2.*m.*n]*[sigmaxx(i) sigmaxy(i) sigmayy(i)]';
    sigma2(i) = [n.^2, m.^2, -2.*m.*n]*[sigmaxx(i) sigmaxy(i) sigmayy(i)]';
    sigma12(i) = [-m.*n, m.*n, (m.^2-n.^2)]*[sigmaxx(i) sigmaxy(i) sigmayy(i)]';
end

sigma1  % Applied stress matrix along fiber direction
sigma2  % Applied stress matrix transverse to fiber direction
sigma12 % Applied shear stress matrix

%-------------------------------------------------------------------------%

% Plane stress condition and Orthotropic lamina is taken into account
% Where F1, F2, F6, F12 are strength factor
% According to Hoffman failure criteria, above mentioned parameters are stated in the function script file.

MOS = zeros(1,numlayers);

for i = 1:numlayers
   sigma1_ = sigma1(i);
   sigma2_ = sigma2(i);
   sigma12_ = sigma12(i);
   safety_factor(sigma1_, sigma2_, sigma12_, Xt, Xc, Yt, Yc, S);
   MOS(i) = ans;
end

MOS
MOS_comp = min(MOS);
disp(['MOS of composite = ' num2str(MOS_comp)])

%-------------------------------------------------------------------------%

% checking condition of each laminas

for i = 1:numlayers
    if (MOS(i)>1)
        disp(['the lamina of orientation ' num2str(theta(i)) ' degree is predicted to be safe'])
    elseif (MOS(i)==1)
        disp(['the lamina of orientation ' num2str(theta(i)) ' degree is prone to fail at the current load level'])
    else
        disp(['the lamina of orientation ' num2str(theta(i)) ' degree has failed under the applied loads.'])
    end
end


%-------------------------------------------------end of program--------------------------------------------------------%



