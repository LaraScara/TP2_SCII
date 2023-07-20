% Caso de estudio 2 - Péndulo invertido apartado 5
clc; clear ; close all

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

% LQR
Q = diag([1/10 1 0.1 100 0.01]);
R = 100;
Kn = lqr(An,Bn,Q,R);
K = Kn(1:4);
Ki = Kn(5);

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
U(1) = 0; 

% Iteración
for i=1:1:n-1
    X_a = X(:,i);%[delta ; delta_p ; phi ; phi_p ]
    
    psi_p = ref(i)-(C(1,:))*(X_a);
    psi(i+1) = psi(i)+dt*psi_p;
    Ua = -K*(X_a-X_op)-Ki*psi(i+1);
    U = [U Ua];

    Xp_a1 = X_a(2);
    Xp_a2 = -F*X_a(2)/M-masa(i)*g*(X_a(3)-pi)/M+Ua/M;
    Xp_a3 = X_a(4);
    Xp_a4 = -F*X_a(2)/(M*l)-g*(M+masa(i))*(X_a(3)-pi)/(M*l)-Ua/(M*l);
    
    Xp_a = [Xp_a1 ; Xp_a2 ; Xp_a3 ; Xp_a4];
    Xf = X_a+ dt*(Xp_a);
    X(:,i+1) = Xf;
end

% Gráficas
figure
plot(t,ref,'b');
hold on;
plot(t,X(1,:),'r');
title('Desplazamiento del carro');
xlabel('Tiempo [s]');
ylabel('Posición [m]');
legend('referencia','?(t)');
grid on;


figure
hold on;
plot(t,X(3,:),'r');
axis([0 tf 3.13 3.15]);
title('Ángulo del péndulo');
xlabel('Tiempo [s]');
ylabel('Ángulo [rad]');
legend('?(t)');
grid on;

figure
plot(t,U,'r');
title('Acción de control');
xlabel('Tiempo [s]');
legend('u(t)');
grid on;

disp("Terminado")
