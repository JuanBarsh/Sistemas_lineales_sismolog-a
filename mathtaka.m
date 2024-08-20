close all
clear all
clc
escala = 2; escs = 4;
%Resolución de una onda SH incidente en un modelo con un estrato dentro de
%un semiespacio

%%Comenzamos definiendo el dominio, en este caso el estrato
h = 30; %m, espesor del estrato

%Propiedades elásticas del estrato, este estrato es homogéneo

betaR = 100; %m/s, velocidad de la onda S
rhoR  = 1100; %kg/m3, densidad del medio

%Propiedades elásticas del semiespacio

betaE = 300; %m/s, velocidad de la onda S
rhoE = 1100; %kg/m3, densidad del semiespacio


% Datos para las estaciones
dz = 10;
nz = 10;
z = 0:dz:dz*(nz-1);
x = 0;


gm = 10; gamma = gm*pi/180; %grados, el ángulo con el que incide la onda SH

%Del cortante:
muR  = rhoR*betaR^2;
muE = rhoE*betaE^2;

N = 1024;   % Número de muestras para el tiempo
dt = 0.02;   % Paso en el tiempo (s)
duracion = N*dt;  % duración de la simulación
t = linspace(0,(N-1)*dt,N);  % tiempo (s)
df = 1/N/dt;
f = linspace(0,(N-1)*df,N);  % Frecuencias  (hz)
f(1,1) = 0.001;

%Parámetros de la señal
%Fuente
tp = 0.5;
ts = 5;
a = (pi/tp).*(t- ts);
% src= (a.^2 - 1/2).*exp(-a.^2); % Onda incidente como pulso de Ricker en t
src = exp(-a.^2); % Onda incidente como gaussiana en t
figure()
plot(t,src,'-b'); % Plot para la onda incidente 
title('Onda incidente')
xlabel('Tiempo (s)')
%hold on;
% Transformada de Fourier de la onda incidente
ondaf = fft(src)*dt;
figure()
plot(f(1:N/2),abs(ondaf(1:N/2)))

nuE = [];
nuR = [];
k = [];
w = [];
fdt = [];

for c1=1:length(z)
    dp = z(c1);
    for c2 = 1:N/2+1
        w = 2*pi*f(c2);
        k = (w/betaE)*sin(gamma);
        nuE = sqrt((w^2/betaE^2)-k^2);
        nuR = sqrt((w^2/betaR^2)-k^2);
        eta = (muR*nuR)/(muE*nuE);
        A = 2*exp(1i*nuE*h)/(cos(nuR*h)+1i*eta*sin(nuR*h));
        if dp < h
            fdt(c2,c1) = A*cos(nuR*dp)*exp(-1i*k*x);
        else
            fdt(c2,c1) = exp(1i*nuE*dp) + (exp(1i*nuE*h)*((cos(nuR*h)-1i*eta*sin(nuR*h))/(cos(nuR*h)+1i*eta*sin(nuR*h))))*exp(-1i*nuE*(dp-h));
        end  
        afdt(c2)=fdt(c2,c1); 
    end
    plot(f(1:N/2+1), abs(afdt(1:N/2+1))*escala-(c1-1)*dz,'k', 'linewidth', 1);
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

