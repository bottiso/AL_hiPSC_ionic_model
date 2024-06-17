% Software implementation of the Paci2020 model of the action potential 
% of human induced pluripotent stem cell-derived cardiomyocytes, 
% used in 10.1016/j.bpj.2020.03.018
%
% This software is provided for NON-COMMERCIAL USE ONLY 
% (read the license included in the zip file).

% clear, close all
options = odeset('MaxStep',1e-3,'InitialStep',2e-5);

%% Original initial state - run the simulations for 800s
%Y=[       -0.070    0.32   0.0002  0    0    1     1     1      0      1      0   0.75  0.75  0   0.1    1    0    9.2    0     0.75    0.3     0.9     0.1     3.67e-4     0.9673    0.00496     0.999      0.10637559     0 ];
% YNames = {'Vm', 'Ca_SR', 'Cai', 'g', 'd', 'f1', 'f2', 'fCa', 'Xr1', 'Xr2', 'Xs', 'h', 'j', 'm', 'Xf', 'q', 'r', 'Nai', 'm_L', 'h_L', 'RyRa', 'RyRo', 'RyRc', 'aur', 'iur' 'u_a', 'u_i',  'o', 'ok2p'};
% YUnits = {'V',   'mM',   'mM',  '-', '-', '-',  '-',  '-',   '-',   '-',   '-',  '-', '-', '-', '-',  '-', '-', 'mM',   '-',   '-',    '-',    '-',    '-', '-', '-', '-', '-'};

% %% SS after 1500s, PACED ACTIVITY 60 bpm
% Obtained with stimFlag = 1;

Y_read=readmatrix('IC_SS.csv'); 
Y=Y_read(end,:);

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
i_stim_frequency = 570; %frequency per minute
T_end = 1000;

load("coefficients.mat");
coefficients(13) = 3.5;

[t,Yc] = ode15s(@Paci2020,[0 T_end],Y, options, stimFlag, i_stim_frequency, tDrugApplication, INaFRedMed, INaLRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed, coefficients);

Vm   = Yc(:,1);
dVm  = [0; diff(Vm)./diff(t)];
caSR = Yc(:,2);
Cai  = Yc(:,3);
Nai  = Yc(:,18);
RyRo = Yc(:, 22);
RyRc = Yc(:, 23);
RyRada = Yc(:, 21);
f1=Yc(:,6);
f2=Yc(:,7);

for i= 1:size(Yc,1)
[~, dati]    = Paci2020(t(i), Yc(i,:), stimFlag, i_stim_frequency, tDrugApplication, INaFRedMed, INaLRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed,coefficients);
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
    tau_f1(i)  = dati(22);
    tau_f2(i)  = dati(23);
end
result       = [INa; If; ICaL; Ito; IKs; IKr; IK1; INaCa; INaK; IpCa; IbNa; IbCa; Irel; Iup; Ileak; Istim; E_K; E_Na; INaL; IKur; IKCa];
mat_correnti = [INa; If; ICaL; Ito; IKs; IKr; IK1; INaCa; INaK; IpCa; IbNa; IbCa; Istim; INaL; IKur; IKCa];
I_tot=sum(mat_correnti);

cd biomarkers
% compute_features_rate;
% 
% indices = (t >= 500) & (t <= 500+APD90mean/1000);
% t_limited = t(indices);
% v_limited = Vm(indices)*1e3-MDPmean;
% I = trapz(t_limited, v_limited);
% disp(['IAPD90 = ', num2str(I)])

tstar=T_end-1.3330; %1.0667; %2;
indices=t-tstar>=0;
csvwrite('PV_9_5.csv', [t(indices)-tstar, Vm(indices)*1000, Cai(indices)*1000])
