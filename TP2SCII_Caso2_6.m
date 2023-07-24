% Caso de estudio 2 - Péndulo invertido apartado 6
% clc; clear ; close all

% Parámetros
m = 0.1;
F = 0.1;
l = 1.6;
g = 9.8;
M = 1.5;  
tf = 100; dt = 1*10^-3;
t = 0:dt:(tf-dt);
n = round(tf/dt);

% Matrices x = [delta ; delta_p ; phi ; phi_p]
A = [0 1 0 0;0 -F/M -m*g/M 0; 0 0 0 1; 0 -F/(l*M) -g*(m+M)/(l*M) 0];
B = [0; 1/M; 0; -1/(l*M)];
C = [1 0 0 0; 0 0 1 0];
D = 0;

% Matrices ampliadas debido al integrador
An = [A zeros(4,1) ; -(C(1,:)) 0];
Bn = [B ; 0];
Cn = [1 0];

% Matrices del sistema observador
Ao = A';
Bo = C';
Co = B';

% LQR
Q = diag([1/10 1 0.1 100 0.01]);
R = 100;
Kn = lqr(An,Bn,Q,R);
K = Kn(1:4)
Ki = Kn(5) 

Qo = diag([0.1 1 100000 10]);
Ro = [10000 0; 0 100];
Ko = lqr(Ao,Bo,Qo,Ro)

% Funciones a utilizar
ref = zeros(1,n);
masa = zeros(1,n);
for j=1:1:n-1
    if (j < (n-1)/2 )
        ref(j) = 2;
        masa(j) = m;
    else (j >= (n-1)/2 )
        ref(j) = 0;
        masa(j) = m*10;
    end
end

% Condidiones iniciales
X_op = [0 ; 0 ; pi ; 0];
X = zeros(4,n);
X(1,1) = 0;   %delta    inicial
X(2,1) = 0;   %delta_p  inicial
X(3,1) = pi;   %phi      inicial
X(4,1) = 0;   %phi_p    inicial
psi(1) = 0;         %psi      inicial
U(1)=0;

Xhat = zeros(4,n);
Xhat(1,1) = 0;   %deltahat    inicial
Xhat(2,1) = 0;   %deltahat_p  inicial
Xhat(3,1) = pi;   %phihat      inicial
Xhat(4,1) = 0;   %phihat_p    inicial

% Iteración
for i=1:1:n-1
    X_a = X(:,i);         %[delta ; delta_p ; phi ; phi_p ]
    Xhat_a = Xhat(:,i);   %[deltahat ; deltahat_p ; phihat ; phihat_p ]
    Y = C*X_a;
    Yhat = C*Xhat_a;
    err = Y-Yhat;
    psi_p = ref(i)-(C(1,:))*(X_a);
    psi(i+1) = psi(i)+dt*psi_p;
    Ua = -K*(Xhat_a-X_op)-Ki*psi(i+1);    % con observador
%     Ua = -K*(X_a-X_op)-Ki*psi(i+1);     % sin observador
    U = [U Ua];

    Xp_a1 = X_a(2);
    Xp_a2 = -F*X_a(2)/M-masa(i)*g*(X_a(3)-pi)/M+Ua/M;
    Xp_a3 = X_a(4);
    Xp_a4 = -F*X_a(2)/(M*l)-g*(M+masa(i))*(X_a(3)-pi)/(M*l)-Ua/(M*l);
    Xp_a = [Xp_a1 ; Xp_a2 ; Xp_a3 ; Xp_a4];
    Xf = X_a+ dt*(Xp_a); % Realizamos la integracion de euler y actualizamos matriz X
    X(:,i+1) = Xf;
    
    % Observador
    Xhatp_a = Ko'*err+A*(X_a-X_op)+B*Ua;
    Xhatf = Xhat_a+dt*Xhatp_a;
    Xhat(:,i+1) = Xhatf;
    
end

% Gráficas
figure(1);hold on;
plot(t,ref,'b');
plot(t,X(1,:),'k'); grid on;
title('Desplazamiento del carro');
xlabel('Tiempo [s]');
ylabel('Posición [m]');
legend('referencia','sin observador','con observador');

figure(2);hold on;
plot(t,X(3,:),'k'); grid on;
title('Ángulo del péndulo');
xlabel('Tiempo [s]');
ylabel('Ángulo [rad]');
legend('sin observador','con observador');
ylim([3.13 3.15]);

disp('Terminado');