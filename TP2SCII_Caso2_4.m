% Caso de estudio 2 - Péndulo invertido apartado 4
clc; clear ; close all

% Parámetros
m = 0.1;
F = 0.1;
l = 1.6;
g = 9.8;
M = 1.5;  
tf = 30; dt = 1*10^-4; t = 0:dt:(tf-dt); 
n=round(tf/dt);

% Matrices x = [delta ; delta_p ; phi ; phi_p]
A = [0,1,0,0 ; 0,(-F/M),(-m*g/M),0 ; 0,0,0,1 ; 0,(F/(l*M)),(((M+m)*g)/(l*M)),0 ];
B = [0 ; (1/M) ; 0 ; (-1/(l*M))];
C = [1,0,0,0];
D = 0;

% Controlabilidad y observabilidad
Co = ctrb(A, B);
rank(Co); % = 4 por ende es controlable
Ob = obsv(A,C);
rank(Ob); % = 4 por ende es observable

% Matrices del sistema observador
Ao=A';
Bo=C';
Co=B';

% LQR
Q = diag([1 1 1000 10000]);
R = 10;
K = lqr(A,B,Q,R);

Qo = diag([1 1/1 1 1/1]); 
Ro = 1;
Ko = lqr(Ao,Bo,Qo,Ro);

%Ganancia de prealimentacion
G = -inv(C*inv(A-B*K)*B);

% Condidiones iniciales
tf = 30; dt = 1*10^-3; t = 0:dt:(tf-dt); 
n=round(tf/dt);
X=zeros(4,n);
X(1,1)=0;   %delta    inicial
X(2,1)=0;   %delta_p  inicial
X(3,1)=0;   %phi      inicial
X(4,1)=0;   %phi_p    inicial

Xhat=zeros(4,n); 
Xhat(1,1)=0;   %deltahat       inicial
Xhat(2,1)=0;   %deltahat_p     inicial
Xhat(3,1)=0;   %phihat         inicial
Xhat(4,1)=0;   %phihat_p       inicial

Xc=zeros(4,n);
Xc(1,1)=0;   %deltac     inicial
Xc(2,1)=0;   %deltac_p   inicial
Xc(3,1)=0;   %phic       inicial
Xc(4,1)=0;   %phic_p     inicial

% Referencia
ref = -10*ones(1,n);

% Iteración
for i=1:1:n-1
    X_a = X(:,i);         %[delta ; delta_p ; phi ; phi_p ]
    Xhat_a = Xhat(:,i);   %[deltahat deltahat_p phihat phihat_p] 
    Y = X_a*C;
    Yhat = Xhat_a*C;
    err = Y-Yhat;
    %Ua=-K*X_a+Ref(i)*G;             %U del sistema normal
    %Ua=-K*Xhat_a+ref(i)*G;          %U del sistema observado
    Ua = -K(2:4)*Xhat_a(2:4)-K(1)*X_a(1)+ref(i)*G;%U para el sistema mixto
    
    Xp_a = A*X_a+B*Ua;
    X(:,i+1) = X_a+ dt*Xp_a;
    % %Xf= X_a+ dt*Xp_a; % Realizamos la integracion de euler y actualizamos matriz X
    % %X(:,i+1)=Xf;
    %OBSERVADOR
    Xhat_p = err*Ko'+A*Xhat_a+Ua*B;
    Xhat(:,i+1) = Xhat_a+dt*Xhat_p;
    
%     %sistema de COMPARACION sin observador
%     Xc_a = Xc(:,i);       %[deltac ; deltac_p ; phic ; phic_p ]
%     Uc_a = -K*Xc_a+ref(i)*G;
%     Uc = [Uc Uc_a];
%     Xc_p = A*Xc_a+B*Uc_a;
%     Xcf = Xc_a+ dt*Xc_p;
%     Xc(:,i+1) = Xcf;
end

% Condidiones iniciales
x0 = [0 0 0 0];
delta(1)=0;
delta_p(1)=0;
phi(1)=0;
phi_p(1) = 0;
u(1) = 0;
phi_pp = 0;

% Iteración
for i=1:1:n-1
x = [delta(i); delta_p(i); phi(i); phi_p(i)];

u(i+1) = -K*x+ref(i)*G;

delta_pp = 1/(M+m) *(-m*l*phi_pp*cos(phi(i))+m*l*(phi_p(i))^2*sin(phi(i))-F*delta_p(i)+u(i));
phi_pp = (1/l)* (g*sin(phi(i))-delta_pp*cos(phi(i)));

delta_p(i+1) = delta_p(i)+dt*delta_pp;
delta(i+1) = delta(i) + dt*delta_p(i);
phi_p(i+1) = phi_p(i) + dt*phi_pp;
phi(i+1) = phi(i) + dt*phi_p(i);
end


% Gráficas
figure
subplot(2,1,1);
plot(t,ref,'k');
grid on
hold on
plot(t,delta,'g');
hold on
plot(t,X(1,:),'r');
title('posicion');xlabel('tiempo[s]');ylabel('posicion[m]');

legend('Referencia','Sin Observador','Con Observador');

subplot(2,1,2);
%plot(t,Ref,'k');
grid on
hold on
plot(t,phi,'g');
hold on
plot(t,X(3,:),'r');
title('angulo');xlabel('tiempo[s]');ylabel('angulo[rad]');
legend('Sin Observador','Con Observador');

% figure
% 
% plot(t,U,'b');title('accion de control');