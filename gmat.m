%Método de la matriz global

clear; close all; 
clc

escala = 2; escs = 8;
%Iniciamos definiendo las propiedades elásticas del medio

%Estrato 1 ----------E1
vsE1 = 50; %m/s
vpE1 = sqrt(3)*vsE1; %m/s
Esp1 = 40; %m
rho1 = 2000; %kg/m3
%Estrato 2 ----------E2
vsE2 = 200; %m/s
vpE2 = sqrt(3)*vsE1; %m/s
Esp2 = 40; %m
rho2 = 2000; %kg/m3
%Semiespacio ----------SS
vsSS = 300; %m/s
vpSS = sqrt(3)*vsE1; %m/s
rhoSS = 2000; %kg/m3
%interfaces
z0 = 0;
z1 = Esp1;
z2 = Esp1+Esp2;
% Datos para las estaciones
dz = 10;
nz = 10;
z = 0:dz:dz*(nz-1);
x = 0;

%mu
mu1 = (vsE1^2)*rho1;
mu2 = (vsE2^2)*rho2;
muSS = (vsSS^2)*rhoSS;

%Ángulo de incidencia
gm = 10; gama = gm*pi/180;

N = 1024;   % Número de muestras para el tiempo
dt = 0.02;   % Paso en el tiempo (s)
duracion = N*dt;  % duración de la simulación
t = linspace(0,(N-1)*dt,N);  % tiempo (s)
df = 1/N/dt;
f = linspace(0,(N-1)*df,N);  % Frecuencias  (hz)
f(1,1) = 0.001;
%Fuente
tp = 0.5;
ts = 5;
a = (pi/tp).*(t- ts);
src= (a.^2 - 1/2).*exp(-a.^2); % Onda incidente como pulso de Ricker en t
%src = exp(-a.^2); % Onda incidente como gaussiana en t
figure()
plot(t,src,'-b'); % Plot para la onda incidente 
title('Onda incidente')
xlabel('Tiempo (s)')
%hold on;
% Transformada de Fourier de la onda incidente
ondaf = fft(src)*dt;
figure()
plot(f(1:N/2),abs(ondaf(1:N/2)))

nu1 = [];
nu2 = [];
nuSS = [];
kSS = [];
ww = [];
fdt = [];
C = [];
figure()
for c2 = 1:length(z)
    dp = z(c2);
    for c1 = 1:N/2+1
    ww = 2*pi*f(c1);
    kSS = (ww/vsSS)*sin(gama);        % Número de onda horizontal (semiespacio)
    nu1 = sqrt((ww^2/vsE1^2)-kSS^2); %ww*cos(gama)/vsE1; %
    nu2 = sqrt((ww^2/vsE2^2)-kSS^2); %ww*cos(gama)/vsE2;%
    nuSS = sqrt((ww^2/vsSS^2)-kSS^2);
    A = zeros(5,5);
    B = zeros(5,1);
    A(1,1) = mu1*nu1; A(1,2) = -mu1*nu1;
    A(2,1) = exp(1i*nu1*(z1-z0)); A(2,2) = exp(-1i*nu1*(z1-z0)); A(2,3) = -1; A(2,4) = -1; 
    A(3,1) = -mu1*nu1*exp(1i*nu1*(z1-z0)); A(3,2) = mu1*nu1*exp(-1i*nu1*(z1-z0)); A(3,3) = mu2*nu2; A(3,4) = -mu2*nu2;
    A(4,3) = exp(1i*nu2*(z2-z1)); A(4,4) = exp(-1i*nu2*(z2-z1)); A(4,5) = -1;
    A(5,3) = -mu2*nu2*exp(1i*nu2*(z2-z1)); A(5,4) = mu2*nu2*exp(-1i*nu2*(z2-z1)); A(5,5) = -muSS*nuSS;

    B(4,1) = 1; B(5,1) = -muSS*nuSS;

    C(:,c1) = A\B;
    A_1 = C(1,c1);
    A_2 = C(3,c1);
    B_1 = C(2,c1);
    B_2  = C(4,c1);
    R = C(5,c1);
    if dp<=Esp1
            fdt(c1,c2) = (A_1*exp(1i*nu1*(z(c2)-z0)) + B_1*exp(-1i*nu1*(z(c2)-z0)))*exp(-1i*kSS*x);
        elseif Esp1<dp<=z2
            fdt(c1,c2) = (A_2*exp(1i*nu2*(z(c2)-z1)) + B_2*exp(-1i*nu2*(z(c2)-z1)))*exp(-1i*kSS*x);
        else
            fdt(c1,c2) = (exp(1i*nuSS*(z2)) + R*exp(-1i*nuSS*(z(c2)-(z2))))*exp(-1i*kSS*x);
        end   
        afdt(c1)=fdt(c1,c2); 
    end
    plot(f(1:N/2+1), abs(afdt(1:N/2+1))*escala-(c2-1)*dz,'k', 'linewidth', 1);
    hold on;
end


%% Sismogramas
V = zeros(length(f),1);
figure;
for ii = 1:nz
    V(1)=0+1i*0;                                % frec inicial
    for jj=2:N/2
        V(jj)=ondaf(jj)*fdt(jj,ii);                     % Convolucion
        V(N-jj+2)=conj(V(jj));                    % Crepa
    end
    V(N/2+1)=real(ondaf(N/2+1)*fdt(N/2+1,ii));         % Nyquist 
    v=ifft(V)/dt;
    plot(t, v*escs-dz*(ii-1), 'k','linewidth',1); % Sismogramas
    hold on;
end
grid on;
title ('Registros V(t, z)')
xlabel ('t (sec)')
ylabel ('z (m)')
axis ([0 +duracion  -dz*nz dz ]);