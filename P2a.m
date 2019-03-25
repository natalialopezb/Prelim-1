%% Prelim 1.2a
close all;
clear all;
clc;

%Global variables
global wI1 w11 w12 w13 w23 rx1 rx2 rx3 Kl1 Kl2 Kl3 tl1 tl2 tl3 Sxp Rlt Dx Dl mu

%Basic Parameters
Dt = 40; %Doubling time (min)
DW = 0.3; %Percentage of dry mass per cell
Gc = 200; %Copies per cell
mRNA_h = 2.1; %mRNA half-life (min)
prot_h = 24*60; %Protein half-life (min)
Cv = 9*10^(-17); %Cell volume
Cm = 2.8*10^(-13); %Cell mass
Cc = 5*10^7; %Cell concentration (cell/mL)
Kep = 60*60; %Elongation rate (nts/min)
Klp = 16.5*60; %Translation rate (aa/min)
Ki = 0.024*60; %Initiation rate (1/min)
Rxt = 1150; %Total RNAP active copies per cell
Lx1 = 1200; %gene1 length (nts)
Lx2 = 2400; %gene2 length (nts)
Lx3 = 600; %gene3 length (nts)
Ll1 = 400; %Protein 1 length (AA)
Ll2 = 800; %Protein 2 length (AA)
Ll3 = 200; %Protein 3 length (AA)
wI1 = 100; %Inducer weight
w11 = 0.000001; %Background weight
w12 = 10; %Weight protein 1 in 2
w13 = 5; %Weight protein 1 in 3
w23 = 25; %Weight protein 2 in 3
Av = 6.023*10^23; %Avogadro number
Rib = 45000; %Number of ribosomes per cell

%Compound Parameters
Rxt = (Rxt*Cc*1000/Av)*(Cv*DW/Cm); %RNAP concentration (mol/gDW)
Rlt = (Rib*Cc*1000/Av)*(Cv*DW/Cm); %Ribosomes concentration (mol/gDW)
Dx = log(2)/mRNA_h; %Degradation rate of mRNA
Dl = log(2)/prot_h; %Degradation rate of proteins
mu = log(2)/Dt; %Dilution factor

Sxp = (1.04*Ki*Cv*DW/Cm)*10^(-6); %Saturation constant (mol/gDW)
Ke1 = Kep/Lx1; %Elongation rate for mRNA1
Ke2 = Kep/Lx2; %Elongation rate for mRNA2
Ke3 = Kep/Lx3; %Elongation rate for mRNA3
tx1 = Ke1/Ki; %Tau for mRNA1
tx2 = Ke2/Ki; %Tau for mRNA2
tx3 = Ke3/Ki; %Tau for mRNA3
MW1 = Lx1*607.4+157.9; %Molecular Weight of gene 1
MW2 = Lx2*607.4+157.9; %Molecular Weight of gene 2
MW3 = Lx3*607.4+157.9; %Molecular Weight of gene 3
Gp1 = (Cc*1000*Gc*50*10^(-9)/MW1)*(Cv*DW/Cm); %Gp for gene 1
Gp2 = (Cc*1000*Gc*50*10^(-9)/MW2)*(Cv*DW/Cm); %Gp for gene 2
Gp3 = (Cc*1000*Gc*50*10^(-9)/MW3)*(Cv*DW/Cm); %Gp for gene 3

rx1 = Ke1*Rxt*(Gp1/(Sxp*tx1+Gp1*tx1+Gp1)); %Transcription rate for gene 1
rx2 = Ke2*Rxt*(Gp2/(Sxp*tx2+Gp2*tx2+Gp2)); %Transcription rate for gene 2
rx3 = Ke3*Rxt*(Gp3/(Sxp*tx3+Gp3*tx3+Gp3)); %Transcription rate for gene 3

Kl1 = Klp/Ll1; %Translation rate protein 1
Kl2 = Klp/Ll2; %Translation rate protein 2
Kl3 = Klp/Ll3; %Translation rate protein 3
tl1 = Kl1/Ki; %Tau for protein 1
tl2 = Kl2/Ki; %Tau for protein 1
tl3 = Kl3/Ki; %Tau for protein 1

%Initial conditions
t_i = 0;
t_f = 460; %Final time (min)
step = 1;
t_span = t_i:step:t_f; %Time vector (min)
[m,n] = size(t_span); %Size of time
I = zeros(n+1,1);
%I(1:160,1) = 10*10^(-3); %Inducer initial concentration (mol/gDW)
x0 = [0;0;0;0;0;0]; %Initial conditions for x vector

%% Case 1 - Incoherent loop steady-state without inducer

[t,X] = ode15s(@(t,x) sys(t,x,I),t_span,x0);
X = X.*(10^(12));

figure(1)
q = plot(t_span,X(:,1),t_span,X(:,2),t_span,X(:,3));
q(1).LineWidth = 1.2;
q(1).Color = 'black';
q(1).LineStyle = '-';
q(2).LineWidth = 1.2;
q(2).Color = [0.4 0.4 0.5];
q(2).LineStyle = '--';
q(3).LineWidth = 1.2;
q(3).Color = [0.6 0.6 0.6];
q(3).LineStyle = '-.';
xlabel('Time [min]','fontweight','bold')
ylabel('mRNA concentration [nmol/gDW]','fontweight','bold')
legend('m1','m2','m3')

figure(2)
q = plot(t_span,X(:,4),t_span,X(:,5),t_span,X(:,6));
q(1).LineWidth = 1.2;
q(1).Color = 'black';
q(1).LineStyle = '-';
q(2).LineWidth = 1.2;
q(2).Color = [0.4 0.4 0.5];
q(2).LineStyle = '--';
q(3).LineWidth = 1.2;
q(3).Color = [0.6 0.6 0.6];
q(3).LineStyle = '-.';
xlabel('Time [min]','fontweight','bold')
ylabel('Protein concentration [nmol/gDW]','fontweight','bold')
legend('p1','p2','p3')




%% Case 2 - Adding inducer

%Initial conditions
t_i = 0;
t_f = 300; %Final time (min)
step = 1;
t_span = t_i:step:t_f; %Time vector (min)
[m,n] = size(t_span); %Size of time
I = zeros(n,1);
I(60:n,1) = 10*10^(-3); %Inducer initial concentration (mol/gDW)

X = X.*(10^(-12));
[i,j] = size(X);
x0 = [X(i,1);X(i,2);X(i,3);X(i,4);X(i,5);X(i,6)]; %Initial conditions for x vector

[t,Y] = ode15s(@(t,x) sys(t,x,I),t_span,x0);
Y = Y.*(10^(12));

figure(3)
q = plot(t_span,Y(:,1),t_span,Y(:,2),t_span,Y(:,3));
q(1).LineWidth = 1.2;
q(1).Color = 'black';
q(1).LineStyle = '-';
q(2).LineWidth = 1.2;
q(2).Color = [0.4 0.4 0.5];
q(2).LineStyle = '--';
q(3).LineWidth = 1.2;
q(3).Color = [0.6 0.6 0.6];
q(3).LineStyle = '-.';
xlabel('Time [min]','fontweight','bold')
ylabel('mRNA concentration [nmol/gDW]','fontweight','bold')
legend('m1','m2','m3')

figure(4)
q = plot(t_span,Y(:,4),t_span,Y(:,5),t_span,Y(:,6));
q(1).LineWidth = 1.2;
q(1).Color = 'black';
q(1).LineStyle = '-';
q(2).LineWidth = 1.2;
q(2).Color = [0.4 0.4 0.5];
q(2).LineStyle = '--';
q(3).LineWidth = 1.2;
q(3).Color = [0.6 0.6 0.6];
q(3).LineStyle = '-.';
xlabel('Time [min]','fontweight','bold')
ylabel('Protein concentration [nmol/gDW]','fontweight','bold')
legend('p1','p2','p3')

