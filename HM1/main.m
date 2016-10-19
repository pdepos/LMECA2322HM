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

% Chenal Duct
duct_w = 0.050;     % [m]
duct_h = nozzle_he; % [m]
duct_d = 0.015;     % [m]
duct_l = 0.270;     % [m]
duct_coeff = (duct_d + duct_w) / duct_w; %[adimensionnel] permet de calculer F(M) en channel avec 2*duct_h

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
 veloC = @(T) sqrt( gamma*R*T);
% Lambda Colebrook
 lambda_cole = @(lambda_init,Re) ((-3*log10(2.03/Re * 1/sqrt(lambda_init))^(-1)))^2;  
% Sutherland Formula Dry Air
mu_ref = 1.716*10^(-5); % [N.s/m²]
S_ref  = 111.0;  %[K]
T_ref  = 273.15; %[K]
mu_T = @(T) ( ((T/T_ref)^(3/2))*((T_ref + S_ref)/(T + S_ref)) )*mu_ref;
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
% 
options = optimset('Display','off');
T0e = T0; pe = pa;


% Calcule du M_in
% Le fanno est dict� par le comportement du nozzle
Astar = nozzle_ht*nozzle_w;
Aexha = nozzle_he*nozzle_w;
Aduct = duct_w * duct_h;
Dh_duct = 4*(duct_w * duct_h) / (2*( duct_w + duct_h));

AstarAfsolve = @(M) ( (((gamma+1)/2)/( 1+ (gamma-1)/2 * M*M ))^((gamma+1)/(2*(gamma-1)))*M) - 0.4;
M_in = fsolve(AstarAfsolve,0.5,options);
% Calcule du M_ex et lambda
lambda0 = 0.5;
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

disp('       p0       p0i       p0e       pi        pe       M_in      M_ex     lambda ');
disp([p0/bar p0i/bar p0e/bar pi/bar pe/bar M_in M_ex lambda ]);
%% =======================================
%x = linspace(0.3,2,100); 
%for i = 1 : 100 
%    y(i) = f_M(x(i)); 
%end
%plot(x,y)

%% =======================================
c1 = @(x) (-1)*sqrt((0.2543^2)-((x-0.03)^2))+0.2573;
c2 = @(x)  (1)*sqrt((0.1537^2)-((x-0.12)^2))-0.1462;
c3 = @(x) sqrt(1-x^2)
x = linspace(-0.5,0.5,1000); 
for i = 1 : 1000 
    y1(i) = c1(x(i));
    y2(i) = c2(x(i));
    y3(i) = c3(x(i));
end
hold;
plot(x,y1);
plot(x,y2);
axis('equal');
hold;
%% =======================================
% Question 2)

%% =======================================
% Question 3)
