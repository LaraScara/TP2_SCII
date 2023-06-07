% Caso de estudio 2 - Péndulo invertido apartado 5
% masa constante en 1
clc;clear all;close all

% Parámetros
m = 1;
F = 0.1;
l = 1.6;
g = 9.8;
M = 1.5; 
t_etapa = 1e-4; tF = 40;
t = 0:t_etapa:tF;
tiempo = round(tF/t_etapa);

% Matrices expandidas (hay un integrador) x = [delta ; delta_p ; phi ; phi_p]
AA = [0 1 0 0 0;0 -F/M -m*g/M 0 0; 0 0 0 1 0; 0 -F/(l*M) -g*(m+M)/(l*M) 0 0;-1 0 0 0 0];
BB = [0; 1/M; 0; 1/(l*M); 0];
CC = [1 0 0 0 0];
D = 0;

% LQR
QQ = 1*diag([.5    5    5   0.1    0.5]);  
RR = 1;
KK = lqr(AA, BB, QQ, RR);

% Referencia
ent =1*square(2*pi*t/40)+1;

% Condidiones iniciales
x0 = [0 0 pi 0];
delta(1) = 0;
delta_p(1) = 0;
phi(1) = pi;
phi_p(1) = 0;
dseta(1) = 0;
u(1) = 0;
phi_pp = 0;

limite = tiempo/2;

% Iteración
for i=1:1:tiempo
    x = [delta(i); delta_p(i); phi(i); phi_p(i); dseta(i)];

    dseta_p = ent(i) - CC(1:4)*x(1:4);
    dseta(i+1) = dseta(i) + t_etapa*dseta_p;

    u(i+1) = -KK(1:4)*x(1:4)-KK(5)*dseta(i);

    delta_pp = 1/(M+m) *(-m*l*phi_pp*cos(phi(i))+m*l*(phi_p(i))^2*sin(phi(i))-F*delta_p(i)+u(i));
    phi_pp = (1/l)* (g*sin(phi(i))-delta_pp*cos(phi(i)));

    delta_p(i+1) = delta_p(i)+t_etapa*delta_pp;
    delta(i+1) = delta(i) + t_etapa*delta_p(i);
    phi_p(i+1) = phi_p(i) + t_etapa*phi_pp;
    phi(i+1) = phi(i) + t_etapa*phi_p(i);
end

figure
plot(t,ent,'b');
hold on
plot(t,delta,'r');
title('Desplazamiento del carro');
xlabel('Tiempo [s]');
ylabel('Posición [m]');
legend('referencia','?(t)');
grid on;

figure
plot(t,phi,'r');
title('Ángulo del péndulo');
xlabel('Tiempo [s]');
ylabel('Ángulo [rad]');
legend('?(t)');
grid on;