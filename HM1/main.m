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

%% =======================================
% Isentropique relation
 T0T   = @(x)  (1+((gamma - 1)/2 )*x*x);
 P0P   = @(x) ((1+((gamma - 1)/2 )*x*x).^(gamma/(gamma-1)));
 ro0ro = @(x) ((1+((gamma - 1)/2 )*x*X).^(    1/(gamma-1))); 
 QmP    = @(p0,t0,A_star) (  ((2/(gamma+1)).^((gamma+1)/(2*(gamma-1)))) * sqrt(gamma/R) * p0 / sqrt(T0)*A_star );
 AstarA = @(M) ( (((gamma+1)/2)/( 1+ (gamma-1)/2 * M*M ))^((gamma+1)/(2*(gamma-1)))*M);
 veloC = @(T) sqrt( gamma*R*T);
% Lambda Colebrook
 lambda_cole = @(lambda_init,Re) ((-3*log10(2.03/Re * 1/sqrt(lambda_init))^(-1)))^2;  
% Sutherland Formula Dry Air
mu_ref = 1.716*10^(-5); % [N.s/m²]
S_ref  = 111.0;  %[K]
T_ref  = 273.15; %[K]
mu_T = @(T) ( ((T/T_ref)^(3/2))*((T_ref + S_ref)/(T + S_ref)) );
% Fanno
 f_M = @(M) ( (1/gamma)*( ((1-(M*M) )/(M*M)) + ((gamma+1)/2)*log(  ((gamma+1)/2*(M*M))/(1+((gamma-1)/2) *(M*M)) ) )); 
 fp0pstar0 = @(M) 1/M * (( 1 + (gamma-1)/2 *M*M )/((gamma+1)/2) )^((gamma+1)/(2*(gamma-1)));
 fTTstar = @(M) ( ((gamma+1)/2) / (1+ (gamma-1)/2 *M*M) );
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

% TO DO: il faut encore faire l'itération pour vérifié Mguess et Lambda
% guess
% TO DO: Err.. comment rafiner M_ex ? 
options = optimset('Display','off');


T0e = T0; pe = pa; 

% For this problematic
T0 = 273; % [K]
pa = 1;   % [atm] ~= [bar]
Mg1 = 0.01; 
M_ex = Mg1;


Astar = pi*nozzle_ht*nozzle_ht;
Aexha = pi*nozzle_he*nozzle_he;
Aduct = duct_w * duct_h;
Dh_duct = 4*(duct_w * duct_h) / (2*( duct_w + duct_h)); % diamètre hydrolique pour un rectangle

lambda0 = 0.5;
lambda = lambda0;
error = 10;
while error > 0.0001
fM = lambda/2 * (duct_l)/Dh_duct;
% On va essayer de trouver M_in 
% pour f(M) = f(M1) - f(M2)  
fM1 = fM + f_M(Mg1);
% On établit l'équation à faire entrer dans fsolve 
f_Mnd = @(M) ( (1/gamma)*( ((1-(M*M) )/(M*M)) + ((gamma+1)/2)*log(  ((gamma+1)/2*(M*M))/(1+((gamma-1)/2) *(M*M)) ) )) - fM1 ; 
M_in = fsolve(f_Mnd,0.1,options);

% calcule des températures Te et Ti
Te = T0e/T0T(M_ex);
Ti = fTTstar(M_in) / fTTstar(M_ex) * Te;
% calcule de ro
ro_e = pe/ (Te* R);
ro_i = pe/ (Ti* R);
% calcule du Reynolds
Re_ex = M_ex * veloC(Te) * Dh_duct *ro_e / mu_T(Te);
Re_in = M_in * veloC(Ti) * Dh_duct *ro_i / mu_T(Ti);
Re_duct = (Re_ex + Re_in)/2;
% recalcule du lambda 
lambda_cole = @(lambda_init) ((-3*log10(2.03/Re_duct * 1/sqrt(lambda_init))^(-1)))-1/sqrt(lambda_init);
lambda = fsolve(lambda_cole,0.1,options);
error = abs(lambda - lambda0);
lambda0 = lambda;
end

T0i= T0T(M_in) * Ti;
p0e = P0P(M_ex)*pe;


p0i = fp0pstar0(M_in)/fp0pstar0(M_ex) * p0e;
pi  = p0i/P0P(M_in);

M_in
M_ex
p0



%calcule dans le nozzle
ro_i = pi/(Ti*R);
ro_t = Aexha/Astar * ro_i * M_in * veloC(Ti);

QmP= @(p0) ( ((2/(gamma+1)).^((gamma+1)/(2*(gamma-1)))) * sqrt(gamma/R) * p0 / sqrt(T0e)) - ro_t;
%p0 = fsolve(QmP,0,options);

% Test du Fanno (pour être sur d'avoir le même graphe)
%x = linspace(0.3,2,100); 
%for i = 1 : 100 
%    y(i) = f_M(x(i)); 
%end
%plot(x,y)

%% =======================================
% Question 2)

%% =======================================
% Question 3)



