% Caso de estudio 1 - Motor CC apartado 2
clc; clear ; close all

% Parámetros
Laa = 5*10^-3;
J = 0.004;
Ra = 0.2;
Bm = 0.005;
Ki = 6.5*10^-5;
Km = 0.055;
%Laa=366e-6; J=5e-9; Ra=55.6; Bm=0; Ki=6.49e-3; Km=6.53e-3;

% Matrices ; X = [ia ; tita ; w];
A = [-Ra/Laa 0 -Km/Laa  ; 0 0 1 ; Ki/J 0 -Bm/J];
B = [1/Laa; 0; 0];
C = [0 1 0];
D = [0];

An = [A zeros(3,1); -C 0];
Bn = [B ; 0];
Cn = [1 0];

% Controlabilidad y observabilidad
Co = ctrb(A, B);
rank(Co); % = 3 por ende es controlable
Ob = obsv(A,C);
rank(Ob); % = 3 por ende es observable

% Matrices del sistema observador
Ao = A';
Bo = C';
Co = B';

% LQR	
Q = diag([10 1/1 1/1 1000000/1]);
R = 500;
K4 = lqr(An,Bn,Q,R);
K = K4(1:3);
K_i = -K4(4);

Qo = diag([1 100000/1 1/1]);
Ro = 5;
Ko = lqr(Ao,Bo,Qo,Ro);

% Implementación de funciones a usar
tf = 30; dt = 1*10^-5; t = 0:dt:(tf-dt); per = 15; %[seg]
Tl = 1.15*10^-3;
ref = pi/2*square(2*pi*t/per); % Función de referencia que varia entre pi/2 y -pi/2
fTl = Tl/2*square(2*pi*t/per)+Tl/2; % Función de torque que varia entre 0 y 1.15*10^-3

% Condiciones iniciales
n = round(tf/dt);
X = zeros(3,n);
X(1,1) = 0; %ia inicial
X(2,1) = 0; %tita inicial
X(3,1) = 0; %w inicial
psi(1) = 0; %psi inicial
Xhat = zeros(3,n);
Xhat(1,1) = 0; %ia_hat inicial
Xhat(2,1) = 0; %tita_hat inicial
Xhat(3,1) = 0; %wr_hat inicial
%Up(1) = 0;

% Iteración
for i=1:1:n-1
    X_a = [X(1,i); X(2,i) ; X(3,i)]; %[ia ; tita ; w]
    Xhat_a = [Xhat(1,i) ; Xhat(2,i) ; Xhat(3,i)]; %[ia_hat ; tita_hat ; w_hat]
    Y = C*X_a;
    Yhat = Co*Xhat_a;
    psi_p = ref(i)-Y;
    psi(i+1) = psi(i)+psi_p*dt;
    U = -K(2:3)*X_a(2:3)-K(1)*Xhat_a(1)+K_i*psi(i+1);
    
    Xp_1 = -Ra/Laa*X_a(1)-Km/Laa*X_a(3)+1/Laa*U;  %ia_p
    Xp_2 = X_a(3);                               %tita_p
    Xp_3 = Ki/J*X_a(1)-Bm/J*X_a(3)-1/J*fTl(i);%    %W_p
    Xp_a = [Xp_1 ; Xp_2 ; Xp_3];
    % Realizamos la integracion de euler y actualizamos matriz X
    Xf = X_a+ dt*Xp_a; 
    
    X(1,i+1) = Xf(1);
    X(2,i+1) = Xf(2);
    X(3,i+1) = Xf(3);
    
    % Observador
    err = (Y-Yhat);
    Xhat_p = U*B+Ko'*err+A*Xhat_a;
    Xhatf = Xhat_a + dt*Xhat_p;
    
    Xhat(1,i+1) = Xhatf(1);
    Xhat(2,i+1) = Xhatf(2);
    Xhat(3,i+1) = Xhatf(3);
end

% Gráficas
figure
%subplot(2,1,1);
plot(t,ref);
grid on
hold on
plot(t,X(2,:),'r');
title('Posición angular del motor');
xlabel('Tiempo [s]');
ylabel('Posición angular [rad]');
legend('referencia','θ(t)')

figure
%subplot(2,1,2);
grid on
hold on
plot(t,X(1,:),'r');
title('Corriente de armadura');
xlabel('Tiempo [s]');
ylabel('Corriente [A]');
legend('ia(t)')

disp("Terminado")
