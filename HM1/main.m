% Fluide and Transfert Mechanical 2
%  by Philippe de Posson   ????-??-??
%     Thanh-Son Tran       8116-12-00
%% =======================================
% Data
% Fluide Air
gamma = 1.4;
R = 287.1;


% Nozzle  
nozzle_de = 0.015;          % [m]
nozzle_he = nozzle_de / 2;  % [m]
nozzle_dt = 0.006;          % [m]
nozzle_ht = nozzle_dt / 2;  % [m]

% Chenal Duct
duct_w = 0.050;     % [m]
duct_h = nozzle_he; % [m]
duct_l = 0.270;     % [m]

% For this problematic
T0 = 300; % [K]
pa = 1;   % [atm] ~= [bar]
%% =======================================
% Isentropique relation
 T0T   = @(x)  (1+((gamma - 1)/2 )*x*x);
 P0P   = @(x) ((1+((gamma - 1)/2 )*x*x).^(gamma/(gamma-1)));
 ro0ro = @(x) ((1+((gamma - 1)/2 )*x*X).^(    1/(gamma-1))); 
 Qm    = @(p0,t0,A_star) (  ((2/(gamma+1)).^((gamma+1)/(2*(gamma-1)))) * sqrt(gamma/R) * p0 / sqrt(T0)*A_star );
% Lambda Colebrook
 lambda_cole = @(lambda_init,Re) ((-3*log10(2.03/Re * 1/sqrt(lambda_init))^(-1)))^2;  
% Sutherland Formula Dry Air
mu_ref = 1.716*10^(-5); % [N.s/m²]
S_ref  = 111.0;  %[K]
T_ref  = 273.15; %[K]
mu_T = @(T) ( ((T/T_ref)^(3/2))*((T_ref + S_ref)/(T + S_ref)) );
% Fanno
 f_M = @(M) ( (1/gamma)*( ((1-(M*M) )/(M*M)) + ((gamma+1)/2)*log(  ((gamma+1)/2*(M*M))/(1+((gamma-1)/2) *(M*M)) ) )); 
%% =======================================
% 3 Questions:
% 1) Determine p0 so
%    a/ sonic at throat and subsonic downstream
%    b/ sonic at throat and supersonic up to a shock located at 70.0 mm
%    c/ sonic at throat and supersonic up to a shock located at 120.0 mm
% 2) Determine the mass flow rate Qm for each case
% 3) Compute the flow condition p(x),p0(x) and M(x).
%% =======================================
% Question 1)
% a/ sonic at throat and subsonic downstream M < 1
% Let chose a M Guess 1 = 0.5 and Lambda Guess 1
T0e = T0; pe = pa; 

Mg1 = 0.5; 
Lambdag1 = 0.03;

Dh_duct = 4*(duct_w * duct_h) / (2*( duct_w + duct_h)) ; % diamètre hydrolique pour un rectangle

lambda = Lambdag1;


FM = lambda/2 * (duct_l)/Dh_duct;



% Test du Fanno (pour être sur d'avoir le même graphe)
x = linspace(0.3,2,100); 
for i = 1 : 100 
    y(i) = f_M(x(i)); 
end
plot(x,y)

%% =======================================
% Question 2)

%% =======================================
% Question 3)



