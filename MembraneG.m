% nonlinear plate model 

n = 1000;    % nunber of notes

h = 0.0332e-3;


D = 1.33e-5;
v = 0.35;

E2D = 12*(1-v^2)/h^2*D;
E_y = E2D/h;

% E2D = 3.608e9*h;
% v = 0.22;
% D = E2D*h^2/12/(1-v^2);


% pp = [132724.13,129311.2238,115728.5466,112556.957,102421.6598,81151.3252,55571.7656];
% aa = [0.454e-3/h,0.447e-3/h,0.430e-3/h,0.428e-3/h,0.423e-3/h,0.409e-3/h,0.401e-3/h];

% pp = [61053.0998,69154.4428,81461.5894,89873.1966,109833.5268,116073.2846];
% aa = [0.391e-3/h,0.391e-3/h,0.393e-3/h,0.396e-3/h,0.414e-3/h,0.429e-3/h];

% pp = [62742.316,63328.3706,69188.9166,81151.3252,90666.094,102421.6598,108178.7844,112556.957,115728.5466];
% aa = [0.405e-3/h,0.407e-3/h,0.407e-3/h,0.409e-3/h,0.413e-3/h,0.423e-3/h,0.423e-3/h,0.428e-3/h,0.429e-3/h];



K_stiff = 0.40;
% K_stiff*a*h/D



prestrain = 0.0015;

presigma0 = prestrain*E2D/(1-v);
presigma = presigma0*h^2/D;


z = zeros(n-1,1);
delta = zeros(n-1,1);
u = zeros(n-1,1);
f = zeros(n-1,1);
g = zeros(n-1,1);
e = zeros(2*n-2,2*n-2);
b = ones(2*n-2,1);
R = ones(2*n-2,1);
m = zeros(2*n-2,1);
w = zeros(n+1,1);
defl = zeros(21,1);
pressure = zeros(21,1);
epsilon_r = zeros(n+1,1);
epsilon_z = zeros(n+1,1);
z_g = zeros(n-1,1);



for j=1:1:length(pp)
    j
%     a = 0.454e-3/h;
%     p = 132724.13;
    p = pp(j);
    a = aa(j);
    

     
    delta_a = a/1000;
%     K_stiff = KK(j);
    
    b = ones(2*n-2,1); 
    R = ones(2*n-2,1);
    u = zeros(n-1,1);
    q = p*h^3/D;
    
    
    
z = 1./16.*q.*a^3.*(-1+1./n.^2.*(1:n-1).^2)./n.*(1:n-1);  %Using the linear solution to initialize
z = z';
zinitial = z;
% for i=1:1:tinetotal
while(norm(R)>1e-6||norm(b)>1e-6)
%Row 1   
      k = 1;
      f(k) = 1/a^2*(n^2+n^2/2/k)*z(k+1) - 1/a^2*(2*n^2+n^2/k^2)*z(k) - presigma*z(k) - 6*n/a*z(k)*(u(k+1)) - 12*v*n/k/a*u(k)*z(k) - 6*z(k)^3 - 1/2*q*k*a/n;
      g(k) = 1/a^2*(n^2+n^2/2/k)*u(k+1) - 1/a^2*(2*n^2+n^2/k^2)*u(k) + n/k/a*(1-v)/2*z(k)^2 + n/2/a*z(k)*z(k+1);
   
      
      e(k,k) = -1/a^2*(2*n^2+n^2/k^2)- presigma - 6*n/a*(u(k+1)) - 12*v*n/k/a*u(k) - 18*z(k)^2;
      e(k,k+1) = 1/a^2*(n^2+n^2/2/k);
      e(k,k+n-1) = -12*v*n/k/a*z(k);
      e(k,k+n) = -6*n/a*z(k);
      
      
      e(k+n-1,k) = n/k/a*(1-v)*z(k) + n/2/a*(z(k+1));
      e(k+n-1,k+1) = n/2/a*z(k);
      e(k+n-1,k+n-1) = -1/a^2*(2*n^2+n^2/k^2);
      e(k+n-1,k+n) = 1/a^2*(n^2+n^2/2/k);
   
   %Row n-1
      k = n-1;
      f(k) = 1/a^2*(n^2+n^2/2/k)*n/(n+v+K_stiff*a*h/D)*z(k) - 1/a^2*(2*n^2+n^2/k^2)*z(k) + 1/a^2*(n^2-n^2/2/k)*z(k-1) - presigma*z(k) - 6*n/a*z(k)*(-u(k-1)) - 12*v*n/k/a*u(k)*z(k) - 6*z(k)^3 - 1/2*q*k*a/n;
      g(k) = - 1/a^2*(2*n^2+n^2/k^2)*u(k) + 1/a^2*(n^2-n^2/2/k)*u(k-1) + n/k/a*(1-v)/2*z(k)^2 + n/2/a*z(k)*(n/(n+v+K_stiff*a*h/D)*z(k)-z(k-1));
   
      e(k,k-1) = 1/a^2*(n^2-n^2/2/k);
      e(k,k) = 1/a^2*(n^2+n^2/2/k)*n/(n+v+K_stiff*a*h/D) - 1/a^2*(2*n^2+n^2/k^2)- presigma - 6*n/a*(-u(k-1)) - 12*v*n/k/a*u(k) - 18*z(k)^2;
      e(k,k+n-2) = 6*n/a*z(k);
      e(k,k+n-1) = -12*v*n/k/a*z(k);
      
      e(k+n-1,k-1) = -n/2/a*z(k);
      e(k+n-1,k) = n/k/a*(1-v)*z(k) + n/a*n/(n+v+K_stiff*a*h/D)*z(k) + n/2/a*(-z(k-1));
      e(k+n-1,k+n-2) = 1/a^2*(n^2-n^2/2/k);
      e(k+n-1,k+n-1) = -1/a^2*(2*n^2+n^2/k^2);
   

   for k=2:1:(n-2) 
      f(k) = 1/a^2*(n^2+n^2/2/k)*z(k+1) - 1/a^2*(2*n^2+n^2/k^2)*z(k) + 1/a^2*(n^2-n^2/2/k)*z(k-1) - presigma*z(k) - 6*n/a*z(k)*(u(k+1)-u(k-1)) - 12*v*n/k/a*u(k)*z(k) - 6*z(k)^3 - 1/2*q*k*a/n;
      g(k) = 1/a^2*(n^2+n^2/2/k)*u(k+1) - 1/a^2*(2*n^2+n^2/k^2)*u(k) + 1/a^2*(n^2-n^2/2/k)*u(k-1) + n/k/a*(1-v)/2*z(k)^2 + n/2/a*z(k)*(z(k+1)-z(k-1));
   
      e(k,k-1) = 1/a^2*(n^2-n^2/2/k);
      e(k,k) = -1/a^2*(2*n^2+n^2/k^2)- presigma - 6*n/a*(u(k+1)-u(k-1)) - 12*v*n/k/a*u(k) - 18*z(k)^2;
      e(k,k+1) = 1/a^2*(n^2+n^2/2/k);
      e(k,k+n-2) = 6*n/a*z(k);
      e(k,k+n-1) = -12*v*n/k/a*z(k);
      e(k,k+n) = -6*n/a*z(k);
      
      e(k+n-1,k-1) = -n/2/a*z(k);
      e(k+n-1,k) = n/k/a*(1-v)*z(k) + n/2/a*(z(k+1)-z(k-1));
      e(k+n-1,k+1) = n/2/a*z(k);
      e(k+n-1,k+n-2) = 1/a^2*(n^2-n^2/2/k);
      e(k+n-1,k+n-1) = -1/a^2*(2*n^2+n^2/k^2);
      e(k+n-1,k+n) = 1/a^2*(n^2+n^2/2/k);
   end
   
   m(1:n-1) = z(1:n-1);
   m(n:2*n-2) = u(1:n-1);
   
   b(1:n-1) = f(1:n-1);
   b(n:2*n-2) = g(1:n-1);
   
   delta = -e\b;
   
   R = delta./m;
   norm(R)
   norm(b)
   m = m + delta;
   
   z(1:n-1) = m(1:n-1);
   u(1:n-1) = m(n:2*n-2);
   
   w(n) = 0;
   w(n-1) = 0 - 1/2*(n/(n+v+K_stiff*a*h/D)*z(n-1) + z(n-1))*a/n;

   for i = n-2:-1:1
       w(i) = w(i+1) - 1/2*(z(i+1)+z(i))*a/n;
   end
    w0 = w(1) - 1/2*z(1)*a/n;
   
end


epsilon_r(1) = (u(2)-0)*n/2/a + 1/2*z(1)^2;
epsilon_z(1) = u(1)*n/a;
for k =2:n-2
epsilon_r(k) = (u(k+1)-u(k-1))*n/2/a + 1/2*z(k)^2;
epsilon_z(k) = u(k)*n/k/a;
end
epsilon_r(n-1) = (0-u(n-2))*n/2/a + 1/2*z(n-1)^2;
epsilon_z(n-1) = u(n-1)*n/(n-1)/a;

z_g(1) = (z(2)-0)*n/2/a;
for k =2:n-2
z_g(k) = (z(k+1)-z(k-1))*n/2/a;
end
z_g(n-1) = (n/(n+v+K_stiff*a*h/D)*z(n-1)-z(n-2))*n/2/a;

r_bend(1:n-1,1) = 1:n-1;
r_bend = r_bend*a/n;


u_strech = E2D/2/(1-v^2)*(epsilon_r.^2 + 2*v*epsilon_r.*epsilon_z + epsilon_z.^2) + E2D/2/(1-v)*prestrain*(epsilon_r+epsilon_z);

u_bend = D/2*(z_g.^2 + z.^2./r_bend.^2 + 2*v./r_bend.*z.*z_g)/h^2;


volume = 0;
for i=1:n-1
    volume = volume + 2*pi*(w(i)*h)*i*(a*h)^2/n^2;
end
    volume0 = volume;
    
U_temp = 0;
for i=1:n-1
    U_temp = U_temp + 2*pi*(u_strech(i)+u_bend(i))*i*(a*h)^2/n^2;% + 2*pi*u_bend(i)*i*a^2/n^2;
end
    U_temp = U_temp + pi*(u_strech(n-1)+u_bend(n-1))*n*(a*h)^2/n^2;% + pi*u_bend(n-1)*(n-1)*a^2/n^2;
    U_strain0(j) = U_temp;

% u_vdw = -gamma*(3/2*(delta0./(delta0+w)).^3 - 1/2*(delta0./(delta0+w)).^9);
% 
% for i=1:n-1
%     U_temp = U_temp + 2*pi*u_vdw(i)*i*(a*h)^2/n^2;
% end

U_total0(j) = U_temp - p*volume0 + 1/2*K_stiff*(n/(n+v+K_stiff*a*h/D)*z(n-1))^2*2*pi*(a*h);


% epsilon_t_total = epsilon_r + prestrain;
% epsilon_z_total = epsilon_z + prestrain;

% len = 0;
% for i = 1:n
%     len = len + sqrt((w(i)-w(i+1))^2 + (a/n)^2);
% end
% 
% epsilon_0 = (len - a)/a;


wlinear = 1./64.*q.*a^4.*(1-1./n.^2.*(0:n).^2).^2;
wlinear = wlinear';

% u(2:n) = u(1:n-1);
% u(1) = 0;
% u(n+1) = 0;

w(n+1) = 0;
w(n) = 0 - 1/2*(n/(n+v+K_stiff*a*h/D)*z(n-1) + z(n-1))*a/n;

for i = n-1:-1:2
    w(i) = w(i+1) - 1/2*(z(i)+z(i-1))*a/n;
end
w(1) = w(2) - 1/2*z(1)*a/n;

height(j+1) = w(1)*h;

xh = 0:a*h/n:a*h;
defm(j,:) = w*h/1e-6;
xm(j,:) = xh/1e-3;

x = -a*h:a*h/n:a*h;
ddef(n+2:2*n+1) = w(2:n+1)*h;
w = w(end:-1:1);
ddef(1:n+1) = w(1:n+1)*h;

x = x/1e-3;
ddef = ddef/1e-6;



% x1 = 0:a*h/n:a*h;
% x1 = xl*1e3;
% defm(j,:) = ddef;
% xm(j,:) = x;



a = a + delta_a;

    b = ones(2*n-2,1); 
    R = ones(2*n-2,1);
    u = zeros(n-1,1);
    q = p*h^3/D;
    
    
    
z = 1./16.*q.*a^3.*(-1+1./n.^2.*(1:n-1).^2)./n.*(1:n-1);  %Using the linear solution to initialize
z = z';
% for i=1:1:tinetotal
while(norm(R)>1e-6||norm(b)>1e-6)
%Row 1   
      k = 1;
      f(k) = 1/a^2*(n^2+n^2/2/k)*z(k+1) - 1/a^2*(2*n^2+n^2/k^2)*z(k) - presigma*z(k) - 6*n/a*z(k)*(u(k+1)) - 12*v*n/k/a*u(k)*z(k) - 6*z(k)^3 - 1/2*q*k*a/n;
      g(k) = 1/a^2*(n^2+n^2/2/k)*u(k+1) - 1/a^2*(2*n^2+n^2/k^2)*u(k) + n/k/a*(1-v)/2*z(k)^2 + n/2/a*z(k)*z(k+1);
   
      
      e(k,k) = -1/a^2*(2*n^2+n^2/k^2)- presigma - 6*n/a*(u(k+1)) - 12*v*n/k/a*u(k) - 18*z(k)^2;
      e(k,k+1) = 1/a^2*(n^2+n^2/2/k);
      e(k,k+n-1) = -12*v*n/k/a*z(k);
      e(k,k+n) = -6*n/a*z(k);
      
      
      e(k+n-1,k) = n/k/a*(1-v)*z(k) + n/2/a*(z(k+1));
      e(k+n-1,k+1) = n/2/a*z(k);
      e(k+n-1,k+n-1) = -1/a^2*(2*n^2+n^2/k^2);
      e(k+n-1,k+n) = 1/a^2*(n^2+n^2/2/k);
   
   %Row n-1
      k = n-1;
      f(k) = 1/a^2*(n^2+n^2/2/k)*n/(n+v+K_stiff*a*h/D)*z(k) - 1/a^2*(2*n^2+n^2/k^2)*z(k) + 1/a^2*(n^2-n^2/2/k)*z(k-1) - presigma*z(k) - 6*n/a*z(k)*(-u(k-1)) - 12*v*n/k/a*u(k)*z(k) - 6*z(k)^3 - 1/2*q*k*a/n;
      g(k) = - 1/a^2*(2*n^2+n^2/k^2)*u(k) + 1/a^2*(n^2-n^2/2/k)*u(k-1) + n/k/a*(1-v)/2*z(k)^2 + n/2/a*z(k)*(n/(n+v+K_stiff*a*h/D)*z(k)-z(k-1));
   
      e(k,k-1) = 1/a^2*(n^2-n^2/2/k);
      e(k,k) = 1/a^2*(n^2+n^2/2/k)*n/(n+v+K_stiff*a*h/D) - 1/a^2*(2*n^2+n^2/k^2)- presigma - 6*n/a*(-u(k-1)) - 12*v*n/k/a*u(k) - 18*z(k)^2;
      e(k,k+n-2) = 6*n/a*z(k);
      e(k,k+n-1) = -12*v*n/k/a*z(k);
      
      e(k+n-1,k-1) = -n/2/a*z(k);
      e(k+n-1,k) = n/k/a*(1-v)*z(k) + n/a*n/(n+v+K_stiff*a*h/D)*z(k) + n/2/a*(-z(k-1));
      e(k+n-1,k+n-2) = 1/a^2*(n^2-n^2/2/k);
      e(k+n-1,k+n-1) = -1/a^2*(2*n^2+n^2/k^2);
   

   for k=2:1:(n-2) 
      f(k) = 1/a^2*(n^2+n^2/2/k)*z(k+1) - 1/a^2*(2*n^2+n^2/k^2)*z(k) + 1/a^2*(n^2-n^2/2/k)*z(k-1) - presigma*z(k) - 6*n/a*z(k)*(u(k+1)-u(k-1)) - 12*v*n/k/a*u(k)*z(k) - 6*z(k)^3 - 1/2*q*k*a/n;
      g(k) = 1/a^2*(n^2+n^2/2/k)*u(k+1) - 1/a^2*(2*n^2+n^2/k^2)*u(k) + 1/a^2*(n^2-n^2/2/k)*u(k-1) + n/k/a*(1-v)/2*z(k)^2 + n/2/a*z(k)*(z(k+1)-z(k-1));
   
      e(k,k-1) = 1/a^2*(n^2-n^2/2/k);
      e(k,k) = -1/a^2*(2*n^2+n^2/k^2)- presigma - 6*n/a*(u(k+1)-u(k-1)) - 12*v*n/k/a*u(k) - 18*z(k)^2;
      e(k,k+1) = 1/a^2*(n^2+n^2/2/k);
      e(k,k+n-2) = 6*n/a*z(k);
      e(k,k+n-1) = -12*v*n/k/a*z(k);
      e(k,k+n) = -6*n/a*z(k);
      
      e(k+n-1,k-1) = -n/2/a*z(k);
      e(k+n-1,k) = n/k/a*(1-v)*z(k) + n/2/a*(z(k+1)-z(k-1));
      e(k+n-1,k+1) = n/2/a*z(k);
      e(k+n-1,k+n-2) = 1/a^2*(n^2-n^2/2/k);
      e(k+n-1,k+n-1) = -1/a^2*(2*n^2+n^2/k^2);
      e(k+n-1,k+n) = 1/a^2*(n^2+n^2/2/k);
   end
   
   m(1:n-1) = z(1:n-1);
   m(n:2*n-2) = u(1:n-1);
   
   b(1:n-1) = f(1:n-1);
   b(n:2*n-2) = g(1:n-1);
   
   delta = -e\b;
   
   R = delta./m;
   m = m + delta;
   
   z(1:n-1) = m(1:n-1);
   u(1:n-1) = m(n:2*n-2);
   
   w(n) = 0;
   w(n-1) = 0 - 1/2*(n/(n+v+K_stiff*a*h/D)*z(n-1) + z(n-1))*a/n;

   for i = n-2:-1:1
       w(i) = w(i+1) - 1/2*(z(i+1)+z(i))*a/n;
   end
    w0 = w(1) - 1/2*z(1)*a/n;
   
end


epsilon_r(1) = (u(2)-0)*n/2/a + 1/2*z(1)^2;
epsilon_z(1) = u(1)*n/a;
for k =2:n-2
epsilon_r(k) = (u(k+1)-u(k-1))*n/2/a + 1/2*z(k)^2;
epsilon_z(k) = u(k)*n/k/a;
end
epsilon_r(n-1) = (0-u(n-2))*n/2/a + 1/2*z(n-1)^2;
epsilon_z(n-1) = u(n-1)*n/(n-1)/a;

z_g(1) = (z(2)-0)*n/2/a;
for k =2:n-2
z_g(k) = (z(k+1)-z(k-1))*n/2/a;
end
z_g(n-1) = (n/(n+v+K_stiff*a*h/D)*z(n-1)-z(n-2))*n/2/a;  

r_bend(1:n-1,1) = 1:n-1;
r_bend = r_bend*a/n;

u_strech = E2D/2/(1-v^2)*(epsilon_r.^2 + 2*v*epsilon_r.*epsilon_z + epsilon_z.^2) + E2D/2/(1-v)*prestrain*(epsilon_r+epsilon_z);

u_bend = D/2*(z_g.^2 + z.^2./r_bend.^2 + 2*v./r_bend.*z.*z_g)/(h)^2;


volume = 0;
for i=1:n-1
    volume = volume + 2*pi*(w(i)*h)*i*(a*h)^2/n^2;
end
    volume_a_plus = volume;
    
U_temp = 0;
for i=1:n-1
    U_temp = U_temp + 2*pi*(u_strech(i)+u_bend(i))*i*(a*h)^2/n^2;% + 2*pi*u_bend(i)*i*a^2/n^2;
end
    U_temp = U_temp + pi*(u_strech(n-1)+u_bend(n-1))*n*(a*h)^2/n^2;% + pi*u_bend(n-1)*(n-1)*a^2/n^2;
    U_strain_a_plus(j) = U_temp;


U_total_a_plus(j) = U_temp - p*volume_a_plus + 1/2*K_stiff*(n/(n+v+K_stiff*a*h/D)*z(n-1))^2*2*pi*(a*h);


G(j) = -(U_total_a_plus(j) - U_total0(j))/(delta_a*h)/2/pi/((a-delta_a)*h);

end
% x = 0:1/n:1;
% ddef = w*h;
% ddef = ddef/ddef(1);

% x = -a*h:a*h/n:a*h;
% ddef(n+2:2*n+1) = w(2:n+1)*h;
% w = w(end:-1:1);
% ddef(1:n+1) = w(1:n+1)*h;
% 
% x = x/1e-3;
% ddef = ddef/1e-6;

% ddef = ddef*10;
% ddef = ddef';
% 

% plot(x,ddef)

for j=1:1:length(pp)
    
    plot(xx(:,j),yy(:,j),'o');
    hold on;
    plot(xm(j,:),defm(j,:));
    hold on;
    
end

% j=1;
% plot(xm(j,:),defm(j,:),xx(:,j),yy(:,j),'o');

%  plot(xm(j,:),defm(j,:),xx(:,j),yy(:,j),'o');


% plot(x,ddef,xx(:,5),xx(:,6),'o')
% plot(xm(1,:),defm(1,:),xm(2,:),defm(2,:),xm(3,:),defm(3,:),xm(4,:),defm(4,:),xm(5,:),defm(5,:),xm(6,:),defm(6,:),xm(7,:),defm(7,:),xx(:,1),xx(:,2),'o',xx(:,3),xx(:,4),'o',xx(:,5),xx(:,6),'o',xx(:,7),xx(:,8),'o',xx(:,9),xx(:,10),'o',xx(:,11),xx(:,12),'o',xx(:,13),xx(:,14),'o')
