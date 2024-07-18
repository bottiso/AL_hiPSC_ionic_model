% Software implementation of the Botti2024 model of the action potential 
% of human induced pluripotent stem cell-derived cardiomyocytes, 
% used in 10.1016/j.bpj.2020.03.018
%
% This software is provided for NON-COMMERCIAL USE ONLY 
% (read the license included in the zip file).

clear, close all
options = odeset('MaxStep',1e-3,'InitialStep',2e-5);

%% Original initial state - run the simulations for 800s
% Obtained with stimFlag = 1;
Y(1) = -0.0908810000000000;          % Vm       'V'
Y(2) = 0.0374480000000000;           % Ca_SR    'mM'
Y(3) = 3.65880000000000e-05;         % Ca_i     'mM'
Y(4) = 5.33930000000000e-06;         % d        '-'
Y(5) = 0.925250000000000;            % f1       '-'
Y(6) = 1;                            % f2       '-'
Y(7) = 0.997840000000000;            % fCa      '-'
Y(8) = 4.10050000000000e-07;         % Xr1      '-'
Y(9) = 0.514850000000000;            % Xr2      '-'
Y(10) = 0.0117060000000000;          % Xs       '-'
Y(11) = 0.973370000000000;           % h        '-'
Y(12) = 0.974170000000000;           % j        '-'
Y(13) = 0.00958940000000000;         % m        '-'
Y(14) = 0.512630000000000;           % Xf       '-'
Y(15) = 0.948970000000000;           % q        '-'
Y(16) = 0.00237220000000000;         % r        '-'
Y(17) = 14.9650000000000;            % Nai      'mM'
Y(18) = 0.000108850000000000;        % mL       '-'
Y(19) = 0.625550000000000;           % hL       '-'
Y(20) = 0.0981820000000000;          % RyRa     '-'
Y(21) = 3.62930000000000e-11;        % RyRo     '-'
Y(22) = 0.984730000000000;           % RyRc     '-'
Y(23) = 0.297490000000000;           % u_a       '-'
Y(24) = 0.995440000000000;           % u_i       '-'
Y(25) = 0.00499310000000000;         % o        '-'

%% Current blockers
tDrugApplication = 10000;
INaFRedMed = 1;
INaLRedMed = 1;
ICaLRedMed = 1;
IKrRedMed  = 1;
IKsRedMed  = 1;
INaCaRedMed = 1;
IfRedMed = 1;

%% 0: spontaneous activity. 1: paced.
stimFlag = 1;           
stimFrequency = 60;     %frequency per minute
T_end = 10;

[t,Yc] = ode15s(@Botti2024,[0 T_end],Y, options, stimFlag, stimFrequency, tDrugApplication, INaFRedMed, INaLRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed);

Vm   = Yc(:,1);
dVm  = [0; diff(Vm)./diff(t)];

for i= 1:size(Yc,1)
[~, dati]    = Botti2024(t(i), Yc(i,:), stimFlag, stimFrequency, tDrugApplication, INaFRedMed, INaLRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed);
    INa(i)   = dati(1);
    If(i)    = dati(2);
    ICaL(i)   = dati(3);
    Ito(i)   = dati(4);
    IKs(i)   = dati(5);
    IKr(i)   = dati(6);
    IK1(i)   = dati(7);
    INaCa(i) = dati(8);
    INaK(i)  = dati(9);
    IpCa(i)  = dati(10);
    IbNa(i)  = dati(11);
    IbCa(i)  = dati(12);
    Irel(i)  = dati(13);
    Iup(i)   = dati(14);
    Ileak(i) = dati(15); 
    Istim(i) = dati(16);
    E_K(i)   = dati(17);
    E_Na(i)  = dati(18);
    INaL(i)  = dati(19);
    IKur(i)  = dati(20);
    IKCa(i)  = dati(21);
end


%% post-processing (plot=
plot(t,Vm*1e3);
xlabel('time [s]')
ylabel('Vm [mV]')