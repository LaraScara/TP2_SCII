% Caso de estudio 1 - Motor CC apartado 1
clc; clear ; close all

% Par�metros
Laa = 5*10^-6;
J = 0.004;
Ra = 0.2;
Bm = 0.005;
Ki = 6.5*10^-5;
Km = 0.055;

% Matrices ; X=[ia ; tita ; w];
A = [-Ra/Laa 0 -Km/Laa  ; 0 0 1 ; Ki/J 0 -Bm/J];
B = [1/Laa; 0; 0];
C = [0 1 0];
D = [0];

% Controlabilidad y observabilidad
Co = ctrb(A, B);
rank(Co); % = 3 por ende es controlable
Ob = obsv(A,C);
rank(Ob); % = 3 por ende es observable

% Implementaci�n de funciones a usar
tf = 30; dt = 1*10^-5; t = 0:dt:(tf-dt); per = 15; %[seg]
Tl = 1.15*10^-3;
ref = pi/2*square(2*pi*t/per); % Funci�n de referencia que varia entre pi/2 y -pi/2
fTl = Tl/2*square(2*pi*t/per)+Tl/2; % Funci�n de torque que varia entre 0 y 1.15*10^-3

% LQR
Q = diag([1 15000 1]);
R = 8;
K = lqr(A,B,Q,R);

% Ganancia de prealimentaci�n
G = -inv(C*inv(A-B*K)*B);

% Condiciones iniciales
n = round(tf/dt);
X = zeros(3,n);
X(1,1) = 0; %ia inicial
X(2,1) = 0; %tita inicial

% Iteraci�n
for i=1:1:n-1
    X_a = [X(1,i); X(2,i) ; X(3,i)];%[ia ; w ; tita]
    U = -K*X_a+ref(i)*G;
    
    Xp_1=-Ra/Laa*X_a(1)-Km/Laa*X_a(3)+1/Laa*U;  %ia_p
    Xp_2= X_a(3);                               %tita_p
    Xp_3=Ki/J*X_a(1)-Bm/J*X_a(3)-1/J*fTl(i);%    %W_p
    
    Xp_a = [Xp_1 ; Xp_2 ; Xp_3];
    
    Xf = X_a + dt*Xp_a; % Realizamos la integracion de euler y actualizamos matriz X
    
    X(1,i+1) = Xf(1);
    X(2,i+1) = Xf(2);
    X(3,i+1) = Xf(3);
    
    %X(:,i+1)=Xf;
end

% Gr�ficas
figure
hold on; grid on;
plot(t,ref);
title('Referencia de entrada');
xlabel('Tiempo [s]');
ylabel('�ngulo [rad]');

figure
hold on; grid on;
plot(t,fTl);
title('Torque de perturbaci�n');
xlabel('Tiempo [s]');
ylabel('Torque [Nm]');

figure
%subplot(2,1,1);
plot(t,ref);
hold on; grid on;
plot(t,X(2,:),'r');
title('�ngulo del motor');
xlabel('Tiempo [s]');
ylabel('�ngulo [rad]');
legend('referencia','?(t)')

figure
%subplot(2,1,2);
hold on; grid on;
plot(t,X(1,:),'r');
title('Corriente de armadura');
xlabel('Tiempo [s]');
ylabel('Corriente [A]');
legend('ia(t)')


 disp("Terminado")