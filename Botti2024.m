% Software implementation of the Botti2024 model of the action potential 
% of human induced pluripotent stem cell-derived cardiomyocytes, 
% used in 10.1016/j.bpj.2020.03.018
%
% This software is provided for NON-COMMERCIAL USE ONLY 
% (read the license included in the zip file).

function [dY, dati] = Botti2024(time, Y, stimFlag, stimFrequency, tDrugApplication, INaFRedMed, INaLRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, IfRedMed)

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% 1: Vm (volt) (in Membrane)
% 2: Ca_SR (millimolar) (in calcium_dynamics)
% 3: Cai (millimolar) (in calcium_dynamics)
% 4: d (dimensionless) (in i_CaL_d_gate)
% 5: f1 (dimensionless) (in i_CaL_f1_gate)
% 6: f2 (dimensionless) (in i_CaL_f2_gate)
% 7: fCa (dimensionless) (in i_CaL_fCa_gate)
% 8: Xr1 (dimensionless) (in i_Kr_Xr1_gate)
% 9: Xr2 (dimensionless) (in i_Kr_Xr2_gate)
% 10: Xs (dimensionless) (in i_Ks_Xs_gate)
% 11: h (dimensionless) (in i_Na_h_gate)
% 12: j (dimensionless) (in i_Na_j_gate)
% 13: m (dimensionless) (in i_Na_m_gate)
% 14: Xf (dimensionless) (in i_f_Xf_gate)
% 15: q (dimensionless) (in i_to_q_gate)
% 16: r (dimensionless) (in i_to_r_gate)
% 17: Nai (millimolar) (in sodium_dynamics)
% 18: m_L (dimensionless) (in i_NaL_m_gate)
% 19: h_L (dimensionless) (in i_NaL_h_gate)
% 20: RyRa (dimensionless) (in calcium_dynamics)
% 21: RyRo (dimensionless) (in calcium_dynamics)
% 22: RyRc (dimensionless) (in calcium_dynamics)
% 23: u_a (dimensionless) (in Courtemanche I_Kur dynamics)
% 24: u_i (dimensionless) (in CourtemancheI_Kur dynamics)
% 25: O (dimensionless) (in  dynamics)

%% Constants
F = 96485.3415;     % coulomb_per_mole (in model_parameters)
R = 8.314472;       % joule_per_mole_kelvin (in model_parameters)
T = 310.0;          % kelvin (in model_parameters) %37°C

%% Cell geometry
V_SR = 465.199;
Vc = 7012.0;
Cm_ventr   = 9.87109e-11;   % farad (in model_parameters)
Cm_bas = Cm_ventr / 1.113;
Cm_atrial = Cm_bas * 0.887;
Cm = Cm_atrial;

%% Extracellular concentrations
Nao = 154.0; %151 % millimolar (in model_parameters)
if stimFlag == 0
    Ko  = 5.4;   % millimolar (in model_parameters)
elseif stimFlag == 1
    Ko  = 4;   % millimolar (IN PACED DC MODEL)
end
Cao = 2.0;   % millimolar (in model_parameters)

%% Intracellular concentrations
Ki = 150.0;   % millimolar (in model_parameters)

%% Nernst potential
E_Na = R*T/F*log(Nao/Y(17));
E_Ca = 0.5*R*T/F*log(Cao/Y(3));
E_K  = R*T/F*log(Ko/Ki);
PkNa = 0.03;   % dimensionless (in electric_potentials)
E_Ks = R*T/F*log((Ko+PkNa*Nao)/(Ki+PkNa*Y(17)));

%% INa adapted from DOI:10.3389/fphys.2018.00080
g_Na_atrial = 9800.86868852711;
i_Na        =  ((time<tDrugApplication)*1+(time >= tDrugApplication)*INaFRedMed)*g_Na_atrial * Y(13)^3.0*Y(11)*Y(12)*(Y(1) - E_Na);

m_inf       = 1 / (1 + exp((Y(1)*1000 + 39)/-11.2));
tau_m       = (0.00001 + 0.00013*exp(-((Y(1)*1000 + 48)/15)^2) + 0.000045 / (1 + exp((Y(1)*1000 + 42)/-5)));
dY(13, 1)   = (m_inf-Y(13))/tau_m;

h_inf       = 1 / (1 + exp((Y(1)*1000 + 66.5)/6.8));
tau_h       = (0.00007 + 0.034 / (1 + exp((Y(1)*1000 + 41)/5.5) + exp(-(Y(1)*1000 + 41)/14)) + 0.0002 / (1 + exp(-(Y(1)*1000 + 79)/14)));
dY(11, 1)   = (h_inf-Y(11))/tau_h;

j_inf       = h_inf;
tau_j       = 10*(0.0007 + 0.15 / (1 + exp((Y(1)*1000 + 41)/5.5) + exp(-(Y(1)*1000 + 41)/14)) + 0.002 / (1 + exp(-(Y(1)*1000 + 79)/14)));
dY(12, 1)   = (j_inf-Y(12))/tau_j;


%% INaL
myCoefTauM  = 1;
tauINaL     = 200; %ms
GNaLmax     = 13.3509157535054; %(S/F)
Vh_hLate    = 87.61;
i_NaL       = ((time<tDrugApplication)*1+(time >= tDrugApplication)*INaLRedMed)*GNaLmax* Y(18)^(3)*Y(19)*(Y(1)-E_Na);

m_inf_L     = 1/(1+exp(-(Y(1)*1000+42.85)/(5.264)));
alpha_m_L   = 1/(1+exp((-60-Y(1)*1000)/5));
beta_m_L    = 0.1/(1+exp((Y(1)*1000+35)/5))+0.1/(1+exp((Y(1)*1000-50)/200));
tau_m_L     = 1/1000 * myCoefTauM*alpha_m_L*beta_m_L;
dY(18, 1)   = (m_inf_L-Y(18))/tau_m_L;

h_inf_L     = 1/(1+exp((Y(1)*1000+Vh_hLate)/(7.488)));
tau_h_L     = 1/1000 * tauINaL;
dY(19, 1)   = (h_inf_L-Y(19))/tau_h_L;

%% If adapted from DOI:10.3389/fphys.2018.00080
g_f         = ((time<tDrugApplication)*1+(time >= tDrugApplication)*IfRedMed)*25.9541335109086;
fNa         = 0.37;
fK          = 1 - fNa;
i_fK        = fK*g_f*Y(14)*(Y(1) - E_K);
i_fNa       = fNa*g_f*Y(14)*(Y(1) - E_Na);
i_f         = i_fK + i_fNa;

Xf_infinity = 1.0/(1.0 + exp((Y(1)*1000 + 69)/8));
tau_Xf      = 5600 / (1 + exp((Y(1)*1000 + 65)/7) + exp(-(Y(1)*1000 + 65)/19));
dY(14, 1)   = 1000*(Xf_infinity-Y(14))/tau_Xf;

%% ICaL
g_CaL       = 7.44995532741049e-05;   % metre_cube_per_F_per_s (in i_CaL)
i_CaL       = ((time<tDrugApplication)*1+(time >= tDrugApplication)*ICaLRedMed)*g_CaL*4.0*Y(1)*F^2.0/(R*T)*(Y(3)*exp(2.0*Y(1)*F/(R*T))-0.341*Cao)/(exp(2.0*Y(1)*F/(R*T))-1.0)*Y(4)*Y(5)*Y(6)*Y(7);

Vd_ventr = -9.1;
Vd_atrial = Vd_ventr + 3.114; 

d_infinity  = 1.0/(1.0+exp(-(Y(1)*1000.0-Vd_atrial)/7.0));

alpha_d     = 0.25+1.4/(1.0+exp((-Y(1)*1000.0-35.0)/13.0));
beta_d      = 1.4/(1.0+exp((Y(1)*1000.0+5.0)/5.0));
gamma_d     = 1.0/(1.0+exp((-Y(1)*1000.0+50.0)/20.0));
tau_d       = (alpha_d*beta_d+gamma_d)*1.0/1000.0;
dY(4, 1)    = (d_infinity-Y(4))/tau_d;

Vf_ventr = -26.0;
Vf_atrial = Vf_ventr +0.774;

f1_inf      = 1.0/(1.0+exp((Y(1)*1000.0-Vf_atrial)/3.0));

if (f1_inf-Y(5) > 0.0)
    constf1 = 1.0+1433.0*(Y(3)-50.0*1.0e-6);
    constf1 = constf1 * 1.35;
else
    constf1 = 1.0;
end
tau_f1      = (20.0+1102.5*exp(-((Y(1)*1000.0+50.0)/15.0)^2.0)+200.0/(1.0+exp((13.0-Y(1)*1000.0)/10.0))+280.0/(1.0+exp((30.0+Y(1)*1000.0)/10.0)))*constf1/1000.0;
dY(5, 1)    = (f1_inf-Y(5))/tau_f1;

Vf2_ventr = -32.0;
Vf2_atrial = Vf2_ventr + 0.774; 

f2_inf      = 0.33+0.67/(1.0+exp((Y(1)*1000.0-Vf2_atrial)/4.0));

constf2     = 1.0;
tau_f2      = (600.0*exp(-(Y(1)*1000.0+50.0)^2.0/400.0)+31.0/(1.0+exp((25.0-Y(1)*1000.0)/10.0))+1.0/(1.0+exp((30.0+Y(1)*1000.0)/10.0)))*constf2/1000.0;
tau_f2_atrial = 2 * tau_f2;
dY(6, 1)    = (f2_inf-Y(6))/tau_f2_atrial;

alpha_fCa   = 1.0/(1.0+(Y(3)/0.0006)^8.0);
beta_fCa    = 0.1/(1.0+exp((Y(3)-0.0009)/0.0001));
gamma_fCa   = 0.3/(1.0+exp((Y(3)-0.00075)/0.0008));
fCa_inf     = (alpha_fCa+beta_fCa+gamma_fCa)/1.3156;
if ((Y(1) > -0.06) && (fCa_inf > Y(7)))
    constfCa = 0.0;
else
    constfCa = 1.0;
end
tau_fCa     = 0.002;   % second (in i_CaL_fCa_gate)
dY(7, 1)    = constfCa*(fCa_inf-Y(7))/tau_fCa;

%% Ito
g_to_atrial = 59.4564976228642;

i_to        = g_to_atrial*(Y(1)-E_K)*Y(15)*Y(16);

q_inf       = 1.0/(1.0+exp((Y(1)*1000.0+53.0)/13.0));
tau_q       = (6.06+39.102/(0.57*exp(-0.08*(Y(1)*1000.0+44.0))+0.065*exp(0.1*(Y(1)*1000.0+45.93))))/1000.0;
dY(15, 1)   = (q_inf-Y(15))/tau_q;

r_inf       = 1.0/(1.0+exp(-(Y(1)*1000.0-22.3)/18.75));
tau_r       = (2.75352+14.40516/(1.037*exp(0.09*(Y(1)*1000.0+30.61))+0.369*exp(-0.12*(Y(1)*1000.0+23.84))))/1000.0;
dY(16, 1)   = (r_inf-Y(16))/tau_r;

%% IKs
g_Ks        = 2.08556136889656;   % S_per_F (in i_Ks)
i_Ks        = ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKsRedMed)*g_Ks*(Y(1)-E_Ks)*Y(10)^2.0*(1.0+0.6/(1.0+(3.8*0.00001/Y(3))^1.4));

Xs_infinity = 1.0/(1.0+exp((-Y(1)*1000.0-20.0)/16.0));
alpha_Xs    = 1100.0/sqrt(1.0+exp((-10.0-Y(1)*1000.0)/6.0));
beta_Xs     = 1.0/(1.0+exp((-60.0+Y(1)*1000.0)/20.0));
tau_Xs      = 1.0*alpha_Xs*beta_Xs/1000.0;
dY(10, 1)   = (Xs_infinity-Y(10))/tau_Xs;

%% IKr
L0           = 0.025;   % dimensionless (in i_Kr_Xr1_gate)
Q            = 2.3;     % dimensionless (in i_Kr_Xr1_gate)
g_Kr         = 33.7444662924660/2.0;   % S_per_F (in i_Kr)
i_Kr         = ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKrRedMed)*g_Kr*(Y(1)-E_K)*Y(8)*Y(9)*sqrt(Ko/5.4);

V_half       = 1000.0*(-R*T/(F*Q)*log((1.0+Cao/2.6)^4.0/(L0*(1.0+Cao/0.58)^4.0))-0.019);

Xr1_inf      = 1.0/(1.0+exp((V_half-Y(1)*1000.0)/4.9));
alpha_Xr1    = 450.0/(1.0+exp((-45.0-Y(1)*1000.0)/10.0));
beta_Xr1     = 6.0/(1.0+exp((30.0+Y(1)*1000.0)/11.5));
tau_Xr1      = 1.0*alpha_Xr1*beta_Xr1/1000.0;
dY(8, 1)     = (Xr1_inf-Y(8))/tau_Xr1;

Xr2_infinity = 1.0/(1.0+exp((Y(1)*1000.0+88.0)/50.0));
alpha_Xr2    = 3.0/(1.0+exp((-60.0-Y(1)*1000.0)/20.0));
beta_Xr2     = 1.12/(1.0+exp((-60.0+Y(1)*1000.0)/20.0));
tau_Xr2      = 1.0*alpha_Xr2*beta_Xr2/1000.0;
dY(9, 1)    = (Xr2_infinity-Y(9))/tau_Xr2;


%% IK1_Koivumaki
g_K1 = 0.169261259039458;
i_K1 = g_K1*Ko^(0.4457)*1000.0*(Y(1)-E_K)/(1.0+exp(1.5*(Y(1)*1000.0-E_K*1000.0+3.6)*F/(R*1000.0*T)));

%% INaCa
KmCa        = 1.38;   % millimolar (in i_NaCa)
KmNai       = 87.5;   % millimolar (in i_NaCa)
Ksat        = 0.1;    % dimensionless (in i_NaCa)
gamma       = 0.35;   % dimensionless (in i_NaCa)
alpha       = 2.16659;
kNaCa_atrial = 3450.73331808885;
i_NaCa      = kNaCa_atrial * (exp(gamma*Y(1)*F/(R*T))*Y(17)^3.0*Cao-exp((gamma-1.0)*Y(1)*F/(R*T))*Nao^3.0*Y(3)*alpha)/((KmNai^3.0+Nao^3.0)*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*Y(1)*F/(R*T))));

%% INaK
Km_K        = 1.0;    % millimolar (in i_NaK)
Km_Na       = 40.0;   % millimolar (in i_NaK)
PNaK_atrial = 2.27178916268546;
i_NaK       = PNaK_atrial * Ko/(Ko+Km_K)*Y(17)/(Y(17)+Km_Na)/(1.0+0.1245*exp(-0.1*Y(1)*F/(R*T))+0.0353*exp(-Y(1)*F/(R*T)));
i_NaK       = i_NaK*1.3; 

%% IpCa
KPCa        = 0.0005;   % millimolar (in i_PCa)
g_PCa       = 0.456955096510698;   % A_per_F (in i_PCa)
i_PCa       = g_PCa*Y(3)/(Y(3)+KPCa);

%% IKur (Courtemanche formulation)
KQ10 = 3.0;
coeff_Kur = 3.5;
g_Kur = coeff_Kur * (0.005 + 0.05/(1+exp((Y(1)*1000.0-15.0)/-13.0))); % nS/pF
i_Kur =  g_Kur * Y(23)^3.0 * Y(24) * (1000.0*Y(1)-1000.0*E_K);

u_a_inf = 1.0 / (1.0+exp(-(1000.0*Y(1)+30.3)/9.6));
u_i_inf = 1.0 / (1.0+exp((1000.0*Y(1)-99.45)/27.48));
alpha_ua = 0.65 / (exp(-(1000.0*Y(1)+10.0)/8.5)+exp(-(1000.0*Y(1)-30.0)/59.0));
beta_ua = 0.65 / (exp(-(1000.0*Y(1)+82.0)/17.0)+2.5);
alpha_ui = 1.0 / (21.0+exp((1000.0*Y(1)-185.0)/28.0));
beta_ui = exp((Y(1)*1000.0-158.0)/16.0);
tau_ua = KQ10/(alpha_ua+beta_ua);
tau_ui = KQ10/(alpha_ui+beta_ui);
dY(23, 1) = (u_a_inf - Y(23))/tau_ua;
dY(24, 1) = (u_i_inf - Y(24))/tau_ui;

%% IKCa (Skibsbye formulation)
g_KCa = 0.0753851103617082;
KCa_on = 47.0e6;
KCa_off = 13.0;
dY(25, 1)= (1.0-Y(25)) * KCa_on * Y(3)^2.0 - Y(25) * KCa_off;
i_KCa = g_KCa * Y(25) *(1.0/(1.0 + exp((Y(1)*1000.0 - E_K*1000.0 + 120.0)/45.0))) * (Y(1)*1000.0 - E_K*1000.0);

%% Background currents
g_b_Na      = 1.14;         % S_per_F (in i_b_Na)
i_b_Na      = g_b_Na*(Y(1)-E_Na);
i_b_Na       = 0.1 * i_b_Na;

g_b_Ca      = 0.8727264;    % S_per_F (in i_b_Ca)
i_b_Ca      = g_b_Ca*(Y(1)-E_Ca);

%% Sarcoplasmic reticulum
VmaxUp		= 0.82205;
VmaxUp_atrial = 0.3924*VmaxUp;
Kup			=  4.40435e-4;
i_up        = VmaxUp_atrial/(1.0+Kup^2.0/Y(3)^2.0);

V_leak		= 4.48209e-4;
i_leak      = (Y(2)-Y(3))*V_leak;

% RyR
g_irel_max	= 55.808061;
g_irel_max_atrial = g_irel_max * 1.3514;
RyRa1       = 0.05169*2.0;
RyRa2       = 0.050001;
RyRahalf    = 0.02632;
RyRohalf    = 0.00944;
RyRchalf    = 0.00167;

RyRSRCass   = (1 - 1/(1 +  exp((Y(2)-0.3)/0.1)));
i_rel       = g_irel_max_atrial * RyRSRCass*Y(21)*Y(22)*(Y(2)-Y(3));

RyRainfss   = RyRa1-RyRa2/(1 + exp((1000*Y(3)-(RyRahalf))/0.0082));
RyRtauadapt = 1; %s
dY(20, 1)    = (RyRainfss- Y(20))/RyRtauadapt;

RyRoinfss   = (1 - 1/(1 +  exp((1000*Y(3)-(Y(20)+ RyRohalf))/0.003)));
if (RyRoinfss>= Y(21))
    RyRtauact = 18.75e-3;       %s
else
    RyRtauact = 0.1*18.75e-3;   %s
end
dY(21, 1)    = (RyRoinfss- Y(21))/(RyRtauact);

RyRcinfss   = (1/(1 + exp((1000*Y(3)-(Y(20)+RyRchalf))/0.001)));
if (RyRcinfss>= Y(22))
    RyRtauinact = 2*87.5e-3;    %s
else
    RyRtauinact = 87.5e-3;      %s
end
dY(22, 1)    = (RyRcinfss- Y(22))/(RyRtauinact);

%% Ca2+ buffering
Buf_C       = 0.25;   % millimolar (in calcium_dynamics)
Buf_SR      = 10.0;   % millimolar (in calcium_dynamics)
Kbuf_C      = 0.001;   % millimolar (in calcium_dynamics)
Kbuf_SR     = 0.3;   % millimolar (in calcium_dynamics)
Cai_bufc    = 1.0/(1.0+Buf_C*Kbuf_C/(Y(3)+Kbuf_C)^2.0);
Ca_SR_bufSR = 1.0/(1.0+Buf_SR*Kbuf_SR/(Y(2)+Kbuf_SR)^2.0);

%% Ionic concentrations
%Nai
dY(17, 1)   = -Cm*(i_Na+i_NaL+i_b_Na+3.0*i_NaK+3.0*i_NaCa+i_fNa)/(F*Vc*1.0e-18);
%Cai
dY(3, 1)    = Cai_bufc*(i_leak-i_up+i_rel-(i_CaL+i_b_Ca+i_PCa-2.0*i_NaCa)*Cm/(2.0*Vc*F*1.0e-18));
%caSR
dY(2, 1)    = Ca_SR_bufSR*Vc/V_SR*(i_up-(i_rel+i_leak));

%% Stimulation
i_stim_Amplitude 		= 14.1e-10; % ampere (in stim_mode)
i_stim_End 				= 1000000.0;   % second (in stim_mode)
i_stim_PulseDuration	= 0.002; %0.005;   % 0.005 second (in stim_mode)
i_stim_Start 			= 0; %493.8;   % 493.8 second (in stim_mode)
i_stim_frequency        = stimFrequency;   % per_second (in stim_mode)
stim_flag 				= stimFlag;   % dimensionless (in stim_mode)
i_stim_Period 			= 60.0/i_stim_frequency;


if stim_flag~=0 && stim_flag~=1
error('Paci2020: wrong pacing! stimFlag can be only 0 (spontaneous) or 1 (paced)');
end

if ((time >= i_stim_Start) && (time <= i_stim_End) && (time-i_stim_Start-floor((time-i_stim_Start)/i_stim_Period)*i_stim_Period <= i_stim_PulseDuration))
    i_stim = stim_flag*i_stim_Amplitude/Cm;
else
    i_stim = 0.0;
end

%% Membrane potential
dY(1, 1) = -(i_KCa+i_Kur+i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_NaL+i_NaCa+i_PCa+i_f+i_b_Na+i_b_Ca-i_stim);

%% Output variables
IKur    = i_Kur;
IK1     = i_K1;
Ito     = i_to;
IKr     = i_Kr;
IKs     = i_Ks;
ICaL    = i_CaL;
INaK    = i_NaK;
INa     = i_Na;
INaCa   = i_NaCa;
IpCa    = i_PCa;
If      = i_f;
IbNa    = i_b_Na;
IbCa    = i_b_Ca;
Irel    = i_rel;
Iup     = i_up;
Ileak   = i_leak;
Istim   = i_stim;
INaL    = i_NaL;
IKCa    = i_KCa;

dati = [INa, If, ICaL, Ito, IKs, IKr, IK1, INaCa, INaK, IpCa, IbNa, IbCa, Irel, Iup, Ileak, Istim, E_K, E_Na, INaL, IKur, IKCa];
