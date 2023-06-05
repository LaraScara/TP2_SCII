% Caso de estudio 2 - Péndulo invertido apartado 3
clc;clear all;close all

% Parámetros
m = 0.1;
F = 0.1;
l = 1.6;
g = 9.8;
M = 1.5; 
t_etapa = 1e-4; tF = 30;
t = 0:t_etapa:tF;
tiempo = round(tF/t_etapa);


% Matrices x = [delta ; delta_p ; phi ; phi_p]
A = [0 1 0 0;0 -F/M -m*g/M 0; 0 0 0 1; 0 F/(l*M) g*(m+M)/(l*M) 0];
B = [0; 1/M; 0; -1/(l*M)];
C = [1 0 0 0];
D = 0;

% Controlabilidad y observabilidad
Co = ctrb(A, B);
rank(Co); % = 4 por ende es controlable
Ob = obsv(A,C);
rank(Ob); % = 4 por ende es observable

% LQR
Q = diag([1 1 1000 10000]);
R = 10;
K = lqr(A,B,Q,R);

% Ganancia de prealimentación
G=-inv(C*inv(A-B*K)*B);

% Referencia
ref=0*t - 10;

% Condidiones iniciales
x0 = [0 0 0 0];
delta(1)=0;
delta_p(1)=0;
phi(1)=0;
phi_p(1) = 0;
u(1) = 0;
phi_pp = 0;

for i=1:1:tiempo
x=[delta(i); delta_p(i); phi(i); phi_p(i)];

u(i+1)=-K*x+ref(i)*G;

delta_pp = 1/(M+m) *(-m*l*phi_pp*cos(phi(i))+m*l*(phi_p(i))^2*sin(phi(i))-F*delta_p(i)+u(i));
phi_pp = (1/l)* (g*sin(phi(i))-delta_pp*cos(phi(i)));

delta_p(i+1) = delta_p(i)+t_etapa*delta_pp;
delta(i+1) = delta(i) + t_etapa*delta_p(i);
phi_p(i+1) = phi_p(i) + t_etapa*phi_pp;
phi(i+1) = phi(i) + t_etapa*phi_p(i);
end

% Gráficas
figure
plot(t,u,'r');
title('Acción de control');
xlabel('Tiempo [s]');
legend('u(t)');
grid on;

figure
plot(t,delta,'r');
grid on;
hold on;
plot(t,ref);
title('Posición del carro');
xlabel('Tiempo [s]');
ylabel('Distancia [m]');
legend('?(t)', "referencia");

figure
plot(t,phi,'r');
title('Ángulo del péndulo');
xlabel('Tiempo [s]');
ylabel('Ángulo [rad]');
legend('?(t)');
grid on;


disp("Terminado")