%%   Box model for DOC 
% ---- The unit of time is **Ma** ------
% ---- The unit of amount is **Gmol** -----
% ---- The unit of distance is **meter** ----
 

% -------  SET PARAMETERS -----------
% clear wF1_DOC

Friverine_low    = 0.5;       % riverine DOC flux 200 Tmol/yr (Gmol/10^6 year)
Friverine_high   = 2;
Fhydrotherm_low  = 0.5;        % hydrothermal input of S to the ocean 720 Gmol/yr
Fhydrotherm_high = 1;
Fphoto_low       = 10.*0.1;       % riverine DOC flux 200 Tmol/yr (Gmol/10^6 year)
Fphoto_high      = 10.*0.2;

U24 = 60*1E6*60*60*24*365*1E6;  % Exchange and circulation fluxes (Sv: 10^6 m3 s-1) m3/10ky
T24 = 0.4*1E6*60*60*24*365*1E6; 
U12 = 20*1E6*60*60*24*365*1E6; 
T12 = 10*1E6*60*60*24*365*1E6; 
U45 = 1*1E6*60*60*24*365*1E6; 
T45 = 0.4*1E6*60*60*24*365*1E6; 
T56 = 1*1E6*60*60*24*365*1E6; 
U62 = 5*1E6*60*60*24*365*1E6; 
T62 = 1*1E6*60*60*24*365*1E6; 
T47 = 19*1E6*60*60*24*365*1E6; 
U43 = 38*1E6*60*60*24*365*1E6; 
T73 = 19*1E6*60*60*24*365*1E6; 
U38 = 48*1E6*60*60*24*365*1E6; 
T38 = 19*1E6*60*60*24*365*1E6; 
T84 = 19*1E6*60*60*24*365*1E6; 


R_DOCdeep_low   = 5;   % DOC efflux in deep sediment 50 umol/m2/day
R_DOCdeep_high  = 50;   

R_DOCcoast_low   = 100;  % DOC efflux in coastal sediment 500 umol/m2/day
R_DOCcoast_high  = 1000;  
fprod            = 1;%0.04;                % fraction of modern primary production for SR

K_DOC_low  = 50;             % Monod for DOC 230 uM  
K_DOC_high = 500;             % Monod for DOC 230 uM  


Q10            = 2.5;
T_seawater     = 25;
T_seawater_ref = 25;
K_temp         = (Q10.^((T_seawater-T_seawater_ref)./10));

tmax     = 1000;      % final time in 10^6a  
maxstep  = 0.1;        % maximum step for integration (10^4a)

% ------- Other parameters ------ 

A1 = 3.8*1e7*1e6;   % m2
A2 = 3.1*1e8*1e6;   
A6 = 1.9*1e7*1e6;
A5 = A6;
A4 = A2;
A7 = 2.9*1e7*1e6;
A8 = 2.9*1e7*1e6;

AREAtot = A1 + A2 + A6 + A7 + A8;

A3 = (AREAtot-A1);

V1 = 3.8*1e6*1e9;  % volume of continental m3
V2 = 3.3*1e7*1e9;  % volume of surface ocean
V3 = 1.3*1e9*1e9;  % volume of deep ocean
V4 = 3.4*1e8*1e9;  % voume of intermediate ocean
V5 = 550000*1e9;   % volume of upwelling surface
V6 = 550000*1e9;   % volume of upwelling slope
V7 = 2.9*1e7*1e9;  % volume of high lat downwelling
V8 = 2.9*1e7*1e9;  % volume of high lat upwelling

Vtotal = V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8;


% ----- SIMULATION  ---------

t0=0;                 % initial time  

DOC10=5E-6*1E-9*1E3  *V1; % 8E-6*1E-9*1E3, 27E-6*1E-9*1E6 % initial amounts in Gmol
DOC20=5E-6*1E-9*1E3  *V2;
DOC30=5E-6*1E-9*1E3  *V3;
DOC40=5E-6*1E-9*1E3  *V4;
DOC50=5E-6*1E-9*1E3  *V5;
DOC60=5E-6*1E-9*1E3  *V6;
DOC70=5E-6*1E-9*1E3  *V7;
DOC80=5E-6*1E-9*1E3  *V8;


W0   = [DOC10 DOC20 DOC30 DOC40 DOC50 DOC60 DOC70 DOC80]';

DOC1conc = @(W) W(1)/V1;   % DOC concentration Gmol/m3
DOC2conc = @(W) W(2)/V2;
DOC3conc = @(W) W(3)/V3;
DOC4conc = @(W) W(4)/V4;
DOC5conc = @(W) W(5)/V5;
DOC6conc = @(W) W(6)/V6;
DOC7conc = @(W) W(7)/V7;
DOC8conc = @(W) W(8)/V8;


%  ----------  EXCHANGE AND CIRCULATION FLUXES ------------


SFE24 = @ (W,t) U24*(DOC2conc(W)-DOC4conc(W));  
SFM24 = @ (W,t) T24*(DOC2conc(W)); 
SFE12 = @ (W,t) U12*(DOC1conc(W)-DOC2conc(W)); 
SFM12 = @ (W,t) T12*(DOC1conc(W)); 
SFE45 = @ (W,t) U45*(DOC4conc(W)-DOC5conc(W)); 
SFM45 = @ (W,t) T45*(DOC4conc(W)); 
SFM56 = @ (W,t) T56*(DOC5conc(W)); 
SFE62 = @ (W,t) U62*(DOC6conc(W)-DOC2conc(W)); 
SFM62 = @ (W,t) T62*(DOC6conc(W)); 
SFM47 = @ (W,t) T47*(DOC4conc(W)); 
SFE43 = @ (W,t) U43*(DOC4conc(W)-DOC3conc(W)); 
SFM73 = @ (W,t) T73*(DOC7conc(W)); 
SFE38 = @ (W,t) U38*(DOC3conc(W)-DOC8conc(W)); 
SFM38 = @ (W,t) T38*(DOC3conc(W)); 
SFM84 = @ (W,t) T84*(DOC8conc(W)); 

K_integr_low  = 0.5;
K_integr_high = 2;

alpha_low  = 4;
alpha_high = 4;

for i=1:300
    
Friverine   = 1E5 *1E6 * unifrnd(Friverine_low ,Friverine_high);
Fhydrotherm = 1E4 *1E6 * unifrnd(Fhydrotherm_low ,Fhydrotherm_high);
Fphoto      = 1E5 *1E6 * unifrnd(Fphoto_low ,Fphoto_high);
R_DOCdeep   = 1E-6*1E-9 *365 * 1E6 * unifrnd(R_DOCdeep_low ,R_DOCdeep_high);
R_DOCcoast  = 1E-6*1E-9 *365 * 1E6 * unifrnd(R_DOCcoast_low ,R_DOCcoast_high);
K_DOC       = 1E-6*1E-9*1E3 * unifrnd(K_DOC_low ,K_DOC_high);
K_integr    = unifrnd(K_integr_low ,K_integr_high);
alpha       = unifrnd(alpha_low ,alpha_high);



% --------- DOC production rate ------------
Integ_R0 = 100 *ones(400,1);

Integ_R0_MC = K_integr.* Integ_R0;
Integ_R     = fprod.*Integ_R0_MC.* 1E-9 * 1E6; %Gmol/m2/Ma Note: Integ_R0 is the depth integrated mineralization rate from coupled particle 1D model

Integ_R1 = Integ_R(20,1);
Integ_R2 = Integ_R(20,1);
Integ_R3 = Integ_R(400,1)-Integ_R(100,1);
Integ_R4 = Integ_R(100,1)-Integ_R(20,1);
Integ_R5 = Integ_R(100,1)-Integ_R(20,1);
Integ_R6 = Integ_R(20,1);
Integ_R7 = 0.1 * Integ_R(100,1);
Integ_R8 = 0.1 * Integ_R(100,1);


F_DOC1 = @(W,t) Integ_R1 * A1;
F_DOC2 = @(W,t) Integ_R2 * A2;
F_DOC3 = @(W,t) Integ_R3 * A3;  
F_DOC4 = @(W,t) Integ_R4 * A4;
F_DOC5 = @(W,t) Integ_R5 * A5;
F_DOC6 = @(W,t) Integ_R6 * A6;
F_DOC7 = @(W,t) Integ_R7 * A7;
F_DOC8 = @(W,t) Integ_R8 * A8;

% --------- DOC flux from sediment ------------

Fsed1 = @(W,t) A1 * fprod * K_temp * R_DOCcoast;
Fsed3 = @(W,t) A3 * fprod * K_temp * R_DOCdeep;
Fsed5 = @(W,t) A6 * fprod * K_temp * R_DOCcoast;

% ------ DOC to DIC ------------

F_DIC1 = @(W,t) alpha * Integ_R1 * A1 * DOC1conc(W)/(DOC1conc(W)+K_DOC);
F_DIC2 = @(W,t) alpha * Integ_R2 * A2 * DOC2conc(W)/(DOC2conc(W)+K_DOC);
F_DIC3 = @(W,t) alpha * Integ_R3 * A3 * DOC3conc(W)/(DOC3conc(W)+K_DOC);  
F_DIC4 = @(W,t) alpha * Integ_R4 * A4 * DOC4conc(W)/(DOC4conc(W)+K_DOC);
F_DIC5 = @(W,t) alpha * Integ_R5 * A5 * DOC5conc(W)/(DOC5conc(W)+K_DOC);
F_DIC6 = @(W,t) alpha * Integ_R6 * A6 * DOC6conc(W)/(DOC6conc(W)+K_DOC);
F_DIC7 = @(W,t) alpha * Integ_R7 * A7 * DOC7conc(W)/(DOC7conc(W)+K_DOC);
F_DIC8 = @(W,t) alpha * Integ_R8 * A8 * DOC8conc(W)/(DOC8conc(W)+K_DOC);


% ------------- SETTING UP THE ODEs ------------------

% DOC

dDOC1dt  = @(t,W) F_DOC1(W,t) - F_DIC1(W,t) + Friverine  - SFE12(W,t) - SFM12(W,t) + Fsed1(W,t) + Fphoto;  % time derivative for concentration in A
dDOC2dt  = @(t,W) F_DOC2(W,t) - F_DIC2(W,t) - SFE24(W,t) - SFM24(W,t) + SFE12(W,t) + SFM12(W,t) + SFE62(W,t) + SFM62(W,t) + Fphoto;
dDOC3dt  = @(t,W) F_DOC3(W,t) - F_DIC3(W,t) + SFM73(W,t) + Fsed3(W,t) - SFM38(W,t) - SFE38(W,t) + SFE43(W,t) - Fhydrotherm;
dDOC4dt  = @(t,W) F_DOC3(W,t) - F_DIC3(W,t) + SFM84(W,t) + SFE24(W,t) + SFM24(W,t) - SFE45(W,t) - SFM45(W,t) - SFM47(W,t) - SFE43(W,t);
dDOC5dt  = @(t,W) F_DOC5(W,t) - F_DIC5(W,t) + SFE45(W,t) + SFM45(W,t) - SFM56(W,t) + Fsed5(W,t);
dDOC6dt  = @(t,W) F_DOC6(W,t) - F_DIC6(W,t) + SFM56(W,t) - SFE62(W,t) - SFM62(W,t);
dDOC7dt  = @(t,W) F_DOC7(W,t) - F_DIC7(W,t) + SFM47(W,t) - SFM73(W,t);
dDOC8dt  = @(t,W) F_DOC8(W,t) - F_DIC8(W,t) + SFM38(W,t) + SFE38(W,t) - SFM84(W,t);


% ----- Call the numerical solver  -----------

dWdt    = @(t,W) [dDOC1dt(t,W) dDOC2dt(t,W) dDOC3dt(t,W) dDOC4dt(t,W) dDOC5dt(t,W) dDOC6dt(t,W) dDOC7dt(t,W) dDOC8dt(t,W)]' ;   % define the derivative 
options = odeset('MaxStep',maxstep);                                           % set solver options
[T,Y]   = ode15s(dWdt,[t0 tmax],W0,options);                                   % call the solver 



DOC1c = Y(:,1)*1E6*1E9/1E3 /V1;  % plot concentrations in uM
DOC2c = Y(:,2)*1E6*1E9/1E3 /V2;
DOC3c = Y(:,3)*1E6*1E9/1E3 /V3;
DOC4c = Y(:,4)*1E6*1E9/1E3 /V4;
DOC5c = Y(:,5)*1E6*1E9/1E3 /V5;
DOC6c = Y(:,6)*1E6*1E9/1E3 /V6;
DOC7c = Y(:,7)*1E6*1E9/1E3 /V7;
DOC8c = Y(:,8)*1E6*1E9/1E3 /V8;

DOC1c_mean(i,1) = mean(DOC1c);
DOC2c_mean(i,1) = mean(DOC2c);
DOC3c_mean(i,1) = mean(DOC3c);
DOC4c_mean(i,1) = mean(DOC4c);
DOC5c_mean(i,1) = mean(DOC5c);
DOC6c_mean(i,1) = mean(DOC6c);
DOC7c_mean(i,1) = mean(DOC7c);
DOC8c_mean(i,1) = mean(DOC8c);


end

M4 = [DOC1c_mean DOC2c_mean DOC3c_mean DOC4c_mean DOC5c_mean DOC6c_mean DOC7c_mean DOC8c_mean];

% 
% for j=1:size(T,1)
% 
%     wF1_DOC(j)=F_DOC2(Y(j,:),T(j))-F_DIC2(Y(j,:),T(j));
% 
% end
% 
% subplot(1,2,1)
% 
% plot(T,DOC2c);
% 
% subplot(1,2,2)
% 
% plot(T,wF1_DOC'.*1E-6.*1E-3);


% -------------- PLOT RESULTS  -------------

subplot(2,4,1);
hist(DOC1c_mean);

subplot(2,4,2);
hist(DOC2c_mean);

subplot(2,4,3);
hist(DOC3c_mean);

subplot(2,4,4);
hist(DOC4c_mean);

subplot(2,4,5);
hist(DOC5c_mean);

subplot(2,4,6);
hist(DOC6c_mean);

subplot(2,4,7);
hist(DOC7c_mean);

subplot(2,4,8);
hist(DOC8c_mean);


