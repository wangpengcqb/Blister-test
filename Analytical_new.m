syms theta r D sigma0 p C1 K a v S E2D;






% theta = C1*besselj(1,(-sigma0/D)^0.5*r) - p/2/sigma0*r;
% 
% theta_1 = diff(theta,r);
% theta_1a = subs(theta_1,r,a);
% theta_a = subs(theta,r,a);


% equ = sprintf('((-sigma0/D)^(1/2)*besselj(0, a*(-sigma0/D)^(1/2)) - besselj(1, a*(-sigma0/D)^(1/2))/a + (K/D + v/a)*besselj(1, a*(-sigma0/D)^(1/2)))*C1 - p/(2*sigma0) - (a*p*(K/D + v/a))/(2*sigma0) = 0');
% equ = subs(equ);
% C_multi = solve(equ,C1);

% CC =  (p/(2*sigma0) + (a*p*(K/D + v/a))/(2*sigma0))/((-sigma0/D)^(1/2)*besselj(0, a*(-sigma0/D)^(1/2)) - besselj(1, a*(-sigma0/D)^(1/2))/a + (K/D + v/a)*besselj(1, a*(-sigma0/D)^(1/2))); 

% C1 = (a*p*(D + K*a + D*v))/(2*(sigma0*besselj(1, a*(-sigma0/D)^(1/2))*(K*a - D + D*v) + D*a*sigma0*(-sigma0/D)^(1/2)*besselj(0, a*(-sigma0/D)^(1/2))));


theta = C1*besselj(1,(-sigma0/D)^0.5*r) - p/2/sigma0*r;

w = - int(theta,r,r,a);

epsilon_r = sigma0/S/(1+v);
epsilon_z = sigma0/S/(1+v);
e0 = sigma0/S/(1+v);

% u_stretch =  0;%sigma0*(epsilon_r+epsilon_z);
u_stretch =  E2D/2/(1-v^2)*(2*(1+v)*e0^2 + (1+v)*e0*theta^2 + 1/4*theta^4);

theta_1 = diff(theta,r);
u_bend = D/2*(theta_1^2 + theta^2/r^2 + 2*v/r*theta*theta_1);

theta_a = subs(theta,r,a);
u_k = 1/2*K*theta_a^2;

Pi_t = 2*pi*int((u_stretch + u_bend)*r,r,0,a) + 2*pi*a*u_k - 2*pi*p*int(w*r,r,0,a);

G = -1/2/pi/a*diff(Pi_t,a);

D = 1.30e-5;
h = 0.0333e-3;
v = 0.30;
K = 0.32;
E2D = 12*(1-v^2)/h^2*D;
S = E2D/(1-v^2);
E_y = E2D/h;

prestrain = 0.0020;
sigma0 = prestrain*(1+v)*S;

% D = 1.33e-5;
% h = 0.0332e-3;
% v = 0.30;
% K = 0.40;
% E2D = 12*(1-v^2)/h^2*D;
% S = E2D/(1-v^2);
% E_y = E2D/h;
% 
% prestrain = 0.0015;
% sigma0 = prestrain*(1+v)*S;





% pp = [40058.5556,61053.0998,69154.4428,81461.5894,89873.1966,109833.5268,116073.2846,118900.1362,123347.2564];
% aa = [0.391173008e-3/h,0.391372954e-3/h,0.391029699e-3/h,0.393262505e-3/h,0.395843515e-3/h,0.413545881e-3/h,0.428987845e-3/h,0.432423458e-3/h,0.448750042e-3/h];

% pp = [55571.7656,62742.316,63328.3706,69188.9166,81151.3252,90666.094,102421.6598,108178.7844,112556.957,115728.5466];
% aa = [0.401486445e-3/h,0.404923708e-3/h,0.406988046e-3/h,0.40681478e-3/h,0.409219048e-3/h,0.413174458e-3/h,0.422798133e-3/h,0.423311329e-3/h,0.427783236e-3/h,0.429157811e-3/h];

for i=1:length(pp)
    p = pp(i);
    a = aa(i);
    

C1 = (a*p*(D + K*a + D*v))/(2*(sigma0*besselj(1, a*(-sigma0/D)^(1/2))*(K*a - D + D*v) + D*a*sigma0*(-sigma0/D)^(1/2)*besselj(0, a*(-sigma0/D)^(1/2))));

Gn(i) = subs(G);
end

Gn = double(Gn);
Gnr = real(Gn);

bb = (aa - aa(1))*1e3;


% theta_1 = diff(theta,r);
% theta_1a = subs(theta_1,r,a);
% theta_a = subs(theta,r,a);
% 
% aaa = theta_1a + (v/a + K/D)*theta_a; 

% w = - int(theta,r,r,a);
% w = (a^2*p)/(4*sigma0) - (p*r^2)/(4*sigma0) + (a*p*(besselj(0, r*(-sigma0/D)^(1/2)) - besselj(0, a*(-sigma0/D)^(1/2)))*(D + K*a + D*v))/(2*(a*sigma0^2*besselj(0, a*(-sigma0/D)^(1/2)) - sigma0*(-sigma0/D)^(1/2)*besselj(1, a*(-sigma0/D)^(1/2))*(K*a - D + D*v)));
% 
% h = (a^2*p)/(4*sigma0) - (a*p*(besselj(0, a*(-sigma0/D)^(1/2)) - 1)*(D + K*a + D*v))/(2*a*sigma0^2*besselj(0, a*(-sigma0/D)^(1/2)) - 2*sigma0*(-sigma0/D)^(1/2)*besselj(1, a*(-sigma0/D)^(1/2))*(K*a - D + D*v));

% diff(theta,r,2) + 1/r*diff(theta,r) - 1/r^2*theta - sigma0/D*theta - p*r/2/D