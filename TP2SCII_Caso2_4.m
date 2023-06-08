% Caso de estudio 2 - Péndulo invertido apartado 4
clc; clear ; close all

% Parámetros
m = 0.1;
F = 0.1;
l = 1.6;
g = 9.8;
M = 1.5; 

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

Qo = diag([1 10/1 1 10/1]); 
Ro = 5;
Ko = lqr(Ao,Bo,Qo,Ro);

%Ganancia de prealimentacion
G = -inv(C*inv(A-B*K)*B);

% Condidiones iniciales
tf = 30; dt = 1*10^-3; t = 0:dt:(tf-dt); 
n = round(tf/dt);

X = zeros(4,n);
X(1,1) = 0;   %delta    inicial
X(2,1) = 0;   %delta_p  inicial
X(3,1) = 0.1;   %phi      inicial
X(4,1) = 0;   %phi_p    inicial

Xhat = zeros(4,n); 
Xhat(1,1) = 0;   %deltahat       inicial
Xhat(2,1) = 0;   %deltahat_p     inicial
Xhat(3,1) = 0;   %phihat         inicial
Xhat(4,1) = 0;   %phihat_p       inicial

Xc = zeros(4,n);
Xc(1,1) = 0;   %deltac     inicial
Xc(2,1) = 0;   %deltac_p   inicial
Xc(3,1) = 0.1;   %phic       inicial
Xc(4,1) = 0;   %phic_p     inicial

Ua(1) = 0;
Uc(1) = 0;

% Referencia
ref = -10*ones(1,n);

% Iteración
for i=1:1:n-1
    X_a = X(:,i);         %[delta ; delta_p ; phi ; phi_p ]
    Xhat_a = Xhat(:,i);   %[deltahat deltahat_p phihat phihat_p] 
    Y = X_a*C;
    Yhat = Xhat_a*C;
    err = Y-Yhat;
    Ua(i+1)=-K*Xhat_a+ref(i)*G;
    
    Xp_a = A*X_a+B*Ua(i);
    X(:,i+1) = X_a+ dt*Xp_a;
  
    % Observador
    Xhat_p = err*Ko'+A*Xhat_a+Ua(i)*B;
    Xhat(:,i+1) = Xhat_a+dt*Xhat_p;
    
    % Sistema sin observador
    Xc_a = Xc(:,i);       %[deltac ; deltac_p ; phic ; phic_p ]
    Uc_a = -K*Xc_a+ref(i)*G;
    Xc_p = A*Xc_a+B*Uc_a;
    Xcf = Xc_a+ dt*Xc_p;
    Xc(:,i+1) = Xcf;
end

% Gráficas
figure
plot(t,ref,'color',[0 0.4470 0.7410]);
grid on
hold on
plot(t,Xc(1,:),'r');
hold on
plot(t,X(1,:),'color',[0.4660 0.6740 0.1880]);
title('Desplazamiento del carro');
xlabel('Tiempo [s]');
ylabel('Posición [m]');
legend('referencia','sin observador','con observador');

figure
grid on
hold on
plot(t,Xc(3,:),'r');
hold on
plot(t,X(3,:),'color',[0.4660 0.6740 0.1880]);
title('Ángulo del péndulo');
xlabel('Tiempo [s]');
ylabel('Ángulo [rad]');
legend('sin observador','con observador');

figure
plot(t,Ua,'r');
title('Acción de control');
xlabel('Tiempo [s]');
legend('u(t)');
grid on;

disp("Terminado")
