% Fluide and Transfert Mechanical 2
%  by Philippe de Posson   5706-10-00
%     Thanh-Son Tran       8116-12-00
%% =======================================
% Data
% Fluide Air
gamma = 1.4;
R = 287.1;


% Nozzle (le nozzle est rectangulaire, pas circulaire)
nozzle_de = 0.015;          % [m]
nozzle_he = nozzle_de / 2;  % [m]
nozzle_dt = 0.006;          % [m]
nozzle_ht = nozzle_dt / 2;  % [m]
nozzle_w  = 0.050;          % [m] largeur constante du nozzle
A_t       = nozzle_dt * nozzle_w; %[m�]

% Chenal Duct
duct_w = 0.050;     % [m]
duct_h = nozzle_he; % [m]
duct_d = 0.015;     % [m]
duct_l = 0.270;     % [m]
duct_coeff = (duct_d + duct_w) / duct_w; %[adimensionnel]
Dh_duct = 4*(duct_w * duct_d) / (2*( duct_w + duct_d)); %diam hydrolique pour Reynolds
A_duct = duct_w * duct_d; %[m�]

% For this problematic
T0 = 300; % [K]
pa = 101325;   % [Pa] ~= 1.01325[bar]
%% =======================================
% Isentropique relation
 T0T   = @(x)  (1+((gamma - 1)/2 )*x*x);
 P0P   = @(x) ((1+((gamma - 1)/2 )*x*x).^(gamma/(gamma-1)));
 ro0ro = @(x) ((1+((gamma - 1)/2 )*x*x).^(    1/(gamma-1))); 
 QmP    = @(p0,T0,A_star) (  ((2/(gamma+1)).^((gamma+1)/(2*(gamma-1)))) * sqrt(gamma/R) * p0 / sqrt(T0)*A_star );
 AstarA = @(M) ( (((gamma+1)/2)/( 1+ (gamma-1)/2 * M*M ))^((gamma+1)/(2*(gamma-1)))*M);
 veloC  = @(T) sqrt( gamma*R*T);
% Lambda Colebrook
% Attention formule erronn�e  lambda_cole = @(lambda_init,Re) ((-3*log10(2.03/Re * 1/sqrt(lambda_init))^(-1)))^2;  
% Sutherland Formula Dry Air
mu_ref  = 1.716*10^(-5); % [N.s/m²]
S_ref   = 111.0;  %[K]
T_ref   = 273.15; %[K]
mu_T    = @(T) ( ((T/T_ref)^(3/2))*((T_ref + S_ref)/(T + S_ref)) )*mu_ref;
% Fanno
 f_M        = @(M) ( (1/gamma)*( ((1-(M*M) )/(M*M)) + ((gamma+1)/2)*log(  ((gamma+1)/2*(M*M))/(1+((gamma-1)/2) *(M*M)) ) )); 
 fp0pstar0  = @(M) 1/M * (( 1 + (gamma-1)/2 *M*M )/((gamma+1)/2) )^((gamma+1)/(2*(gamma-1)));
 fTTstar    = @(M) ( ((gamma+1)/2) / (1+ (gamma-1)/2 *M*M) );
% Shock Relations
P0sh2_P0sh1 = @(M) ( ((gamma+1)/2) / (gamma*M^2-(gamma-1)/2))^(1/(gamma-1)) * ( ( (gamma+1)/2*M^2 ) / ( 1+(gamma-1)/2*M^2 ) )^(gamma/(gamma-1));
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
% 
options = optimset('Display','off');
T0e = T0; pe = pa;


% Calcule du M_in
% Le fanno est dict� par le comportement du nozzle
Astar = nozzle_dt*nozzle_w;
Aexha = nozzle_de*nozzle_w;

AstarAfsolve = @(M) ( (((gamma+1)/2)/( 1+ (gamma-1)/2 * M*M ))^((gamma+1)/(2*(gamma-1)))*M) - 0.4;
M_in = fsolve(AstarAfsolve,0.5,options);
% Calcule du M_ex et lambda
lambda0 = 0.009;
lambda = lambda0;
error = 1;
while error > 0.0001
   %fM = lambda/2 * (duct_l)/Dh_duct;
   fM = lambda * duct_coeff * duct_l / duct_d;
   % On va essayer de trouver M_in 
   % pour f(M) = f(M1) - f(M2)  
   fM1 = f_M(M_in) - fM;

   % On établit l'équation à faire entrer dans fsolve 
   
   f_Mnd = @(M) ( (1/gamma)*( ((1-(M*M) )/(M*M)) + ((gamma+1)/2)*log(  ((gamma+1)/2*(M*M))/(1+((gamma-1)/2) *(M*M)) ) )) - fM1 ; 
   M_ex = fsolve(f_Mnd,0.5,options);

   % calcule des températures Te et Ti
   Te = T0e/T0T(M_ex);
   %Ti = fTTstar(M_in) / fTTstar(M_ex) * Te;
   Ti = T0e/T0T(M_in); 

   % calcule de ro
   ro_e = pe/ (Te* R);
   ro_i = pe/ (Ti* R);
 
   % calcule du Reynolds
   Re_ex = M_ex * veloC(Te) * Dh_duct *ro_e / mu_T(Te);
   Re_in = M_in * veloC(Ti) * Dh_duct *ro_i / mu_T(Ti);
   Re_duct = (Re_ex + Re_in)/2;

   % recalcule du lambda 
   
   %lambda_cole = @(lambda_init) ((-3*log10(2.03/Re_duct * 1/sqrt(lambda_init))))-1/sqrt(lambda_init);
   %lambda = fsolve(lambda_cole,lambda0)
   interLambda_old= lambda0;
   errorLambda = 1;
   while errorLambda > 0.001
       interLambda_new = (-3*log10(2.03/Re_duct * 1/sqrt(interLambda_old)))^(-2);
       errorLambda = interLambda_new - interLambda_old; 
       interLambda_old = interLambda_new;
   end
   lambda = interLambda_new;
   error = abs(lambda - lambda0);
   lambda0 = lambda;
end
  
T0i= T0T(M_in) * Ti;
p0e = P0P(M_ex)*pe;

p0i = fp0pstar0(M_in)/fp0pstar0(M_ex) * p0e;
pi  = p0i/P0P(M_in);
p0 = p0i;
bar = 10^5;

%disp('       p0       p0i       p0e       pi        pe       M_in      M_ex     lambda ');
%disp([p0/bar p0i/bar p0e/bar pi/bar pe/bar M_in M_ex lambda ]);
%% =======================================
%x = linspace(0.3,2,100); 
%for i = 1 : 100 
%    y(i) = f_M(x(i)); 
%end
%plot(x,y)

%% =======================================
%x = linspace(0,0.12,1000); 
%for i = 1 : 1000 

 %   y4(i) = height(x(i));
%end
%hold;
%plot(x,y4);
%axis('equal');
%hold;
%% =======================================
% Question 2)
% Sonic at throat, shock at x = 0.07, subsonic until the end.

height007 = 2 * height(0.07);
A007      = height007 * nozzle_w;
A_tA007 = A_t / A007;

A_tA007fsolve = @(M) ( (((gamma+1)/2)/( 1+ (gamma-1)/2 * M*M ))^((gamma+1)/(2*(gamma-1)))*M) - A_tA007;
M_sh12 = fsolve(A_tA007fsolve,1.2,options);

M_sh22 = (1 + (gamma - 1)/2*M_sh12^2 ) / (gamma*M_sh12^2 - (gamma - 1)/2 );
Astar22 = A007 * ( ( ((gamma+1)/2) / (1+(gamma-1)/2 * M_sh22*M_sh22) )^((gamma+1)/(2*(gamma-1))) * M_sh22 );
Astar22Aduct = Astar22 / A_duct;

Astar22Aductfsolve = @(M) ( (((gamma+1)/2)/( 1+ (gamma-1)/2 * M*M ))^((gamma+1)/(2*(gamma-1)))*M) - Astar22Aduct;
M_ind2 = fsolve(Astar22Aductfsolve,0.5,options);

lambda02 = 0.009;
lambda2 = lambda02;
error2 = 1;
while error2 > 0.0001
    
   fM2 = lambda2 * duct_coeff * duct_l / duct_d;

   % pour f(M) = f(M1) - f(M2)  
   fM12 = f_M(M_ind2) - fM2;

   % On établit l'équation à faire entrer dans fsolve 
   
   f_Mnd = @(M) ( (1/gamma)*( ((1-(M*M) )/(M*M)) + ((gamma+1)/2)*log(  ((gamma+1)/2*(M*M))/(1+((gamma-1)/2) *(M*M)) ) )) - fM12 ; 
   M_ex2 = fsolve(f_Mnd,0.5,options);

   % calcule des températures Te et Ti
   Te2 = T0/T0T(M_ex2);
   Ti2 = T0/T0T(M_ind2); 

   % calcule de ro
   ro_e2 = pe/ (Te2* R);
   ro_i2 = pe/ (Ti2* R);
 
   % calcule du Reynolds
   Re_ex2 = M_ex2 * veloC(Te2) * Dh_duct *ro_e2 / mu_T(Te2);
   Re_in2 = M_ind2 * veloC(Ti2) * Dh_duct *ro_i2 / mu_T(Ti2);
   Re_duct2 = (Re_ex2 + Re_in2)/2;

   % recalcule du lambda 
   
   %lambda_cole = @(lambda_init) ((-3*log10(2.03/Re_duct * 1/sqrt(lambda_init))))-1/sqrt(lambda_init);
   %lambda = fsolve(lambda_cole,lambda0)
   interLambda_old2= lambda02;
   errorLambda2 = 1;
   while errorLambda2 > 0.001
       interLambda_new2 = (-3*log10(2.03/Re_duct2 * 1/sqrt(interLambda_old2)))^(-2);
       errorLambda2 = interLambda_new2 - interLambda_old2; 
       interLambda_old2 = interLambda_new2;
   end
   lambda2 = interLambda_new2;
   error2 = abs(lambda2 - lambda02);
   lambda02 = lambda2;
end


p0e2    = P0P(M_ex2)*pe;                               % Pression_0 sortie duct 
p0i2    = fp0pstar0(M_ind2)/fp0pstar0(M_ex2) * p0e2;   % Pression_0 entree duct and after shock
pi2     = p0i2/P0P(M_ind2);                            % Pression entree duct
p0sh12   = p0i2 / P0sh2_P0sh1(M_sh12);
p02     = p0sh12;

%% =======================================
% Question 3)

A_tA_duct = A_t / A_duct;

A_tA_ductfsolve = @(M) ( (((gamma+1)/2)/( 1+ (gamma-1)/2 * M*M ))^((gamma+1)/(2*(gamma-1)))*M) - A_tA_duct;
M_sh13 = fsolve(A_tA_ductfsolve,1.2,options);

M_sh23 = (1 + (gamma - 1)/2*M_sh13^2 ) / (gamma*M_sh13^2 - (gamma - 1)/2 );
Astar23 = A_duct * ( ( ((gamma+1)/2) / (1+(gamma-1)/2 * M_sh23*M_sh23) )^((gamma+1)/(2*(gamma-1))) * M_sh23 );
Astar23Aduct = Astar23 / A_duct;

%Astar23Aductfsolve = @(M) ( (((gamma+1)/2)/( 1+ (gamma-1)/2 * M*M ))^((gamma+1)/(2*(gamma-1)))*M) - Astar23Aduct;
%M_ind3 = fsolve(Astar23Aductfsolve,0.5,options);
M_ind3  = M_sh23;

lambda03 = 0.009;
lambda3 = lambda03;
error3 = 1;
while error3 > 0.0001
    
   fM3 = lambda3 * duct_coeff * duct_l / duct_d;

   % pour f(M) = f(M1) - f(M2)  
   fM13 = f_M(M_ind3) - fM3;

   % On établit l'équation à faire entrer dans fsolve 
   
   f_Mnd = @(M) ( (1/gamma)*( ((1-(M*M) )/(M*M)) + ((gamma+1)/2)*log(  ((gamma+1)/2*(M*M))/(1+((gamma-1)/2) *(M*M)) ) )) - fM13 ; 
   M_ex3 = fsolve(f_Mnd,0.5,options);

   % calcule des températures Te et Ti
   Te3 = T0/T0T(M_ex3);
   Ti3 = T0/T0T(M_ind3); 

   % calcule de ro
   ro_e3 = pe/ (Te3* R);
   ro_i3 = pe/ (Ti3* R);
 
   % calcule du Reynolds
   Re_ex3 = M_ex3  * veloC(Te3) * Dh_duct *ro_e3 / mu_T(Te3);
   Re_in3 = M_ind3 * veloC(Ti3) * Dh_duct *ro_i3 / mu_T(Ti3);
   Re_duct3 = (Re_ex3 + Re_in3)/2;

   % recalcule du lambda 
   
   %lambda_cole = @(lambda_init) ((-3*log10(2.03/Re_duct * 1/sqrt(lambda_init))))-1/sqrt(lambda_init);
   %lambda = fsolve(lambda_cole,lambda0)
   interLambda_old3 = lambda03;
   errorLambda3 = 1;
   while errorLambda3 > 0.001
       interLambda_new3 = (-3*log10(2.03/Re_duct3 * 1/sqrt(interLambda_old3)))^(-2);
       errorLambda3 = interLambda_new3 - interLambda_old3; 
       interLambda_old3 = interLambda_new3;
   end
   lambda3 = interLambda_new3;
   error3 = abs(lambda3 - lambda03);
   lambda03 = lambda3;
end


p0e3    = P0P(M_ex3)*pe;                               % Pression_0 sortie duct 
p0i3    = fp0pstar0(M_ind3)/fp0pstar0(M_ex3) * p0e3;   % Pression_0 entree duct and after shock
pi3     = p0i3/P0P(M_ind3);                            % Pression entree duct
p0sh13   = p0i3 / P0sh2_P0sh1(M_sh13);
p03     = p0sh13;

