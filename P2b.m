%% Prelim 1.2b
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

Par = [Dt DW Gc mRNA_h prot_h Cv Cm Cc Kep Klp Ki Rxt Lx1 Lx2 Lx3 Ll1 ...
    Ll2 Ll3 wI1 w11 w12 w13 w23 Rib]; %Vector with basic parameters

%% Case 1 - Incoherent loop steady-state without inducer

%Let's loop!
[a,stop] = size(Par);
z = 0; %For me

for i=1:stop

    %Compound Parameters
    tempPar = Par(1,i);
    Par(1,i) = Par(1,i)/10; %Decrease parameter 1 order of magnitude
    Rxt = (Par(1,12)*Par(1,8)*1000/Av)*(Par(1,6)*Par(1,2)/Par(1,7)); %RNAP concentration (mol/gDW)
    Rlt = (Par(1,24)*Par(1,8)*1000/Av)*(Par(1,6)*Par(1,2)/Par(1,7)); %Ribosomes concentration (mol/gDW)
    Dx = log(2)/Par(1,4); %Degradation rate of mRNA
    Dl = log(2)/Par(1,5); %Degradation rate of proteins
    mu = log(2)/Par(1,1); %Dilution factor

    Sxp = (1.04*Par(1,11)*Par(1,6)*Par(1,2)/Par(1,7))*10^(-6); %Saturation constant (mol/gDW)
    Ke1 = Par(1,9)/Par(1,13); %Elongation rate for mRNA1
    Ke2 = Par(1,9)/Par(1,14); %Elongation rate for mRNA2
    Ke3 = Par(1,9)/Par(1,15); %Elongation rate for mRNA3
    tx1 = Ke1/Par(1,11); %Tau for mRNA1
    tx2 = Ke2/Par(1,11); %Tau for mRNA2
    tx3 = Ke3/Par(1,11); %Tau for mRNA3
    MW1 = Par(1,13)*607.4+157.9; %Molecular Weight of gene 1
    MW2 = Par(1,14)*607.4+157.9; %Molecular Weight of gene 2
    MW3 = Par(1,15)*607.4+157.9; %Molecular Weight of gene 3
    Gp1 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW1)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 1
    Gp2 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW2)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 2
    Gp3 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW3)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 3

    rx1 = Ke1*Rxt*(Gp1/(Sxp*tx1+Gp1*tx1+Gp1)); %Transcription rate for gene 1
    rx2 = Ke2*Rxt*(Gp2/(Sxp*tx2+Gp2*tx2+Gp2)); %Transcription rate for gene 2
    rx3 = Ke3*Rxt*(Gp3/(Sxp*tx3+Gp3*tx3+Gp3)); %Transcription rate for gene 3

    Kl1 = Par(1,10)/Par(1,16); %Translation rate protein 1
    Kl2 = Par(1,10)/Par(1,17); %Translation rate protein 2
    Kl3 = Par(1,10)/Par(1,18); %Translation rate protein 3
    tl1 = Kl1/Par(1,11); %Tau for protein 1
    tl2 = Kl2/Par(1,11); %Tau for protein 1
    tl3 = Kl3/Par(1,11); %Tau for protein 1
    
    wI1 = Par(1,19);
    w11 = Par(1,20);
    w12 = Par(1,21);
    w13 = Par(1,22);
    w23 = Par(1,23);

    %Initial conditions
    t_i = 0;
    t_f = 20; %Final time (min)
    step = 1;
    t_span = t_i:step:t_f; %Time vector (min)
    [m,n] = size(t_span); %Size of time
    I = zeros(n+1,1);
    %I(1:160,1) = 10*10^(-3); %Inducer initial concentration (mol/gDW)
    x0 = [0;0;0;0;0;0]; %Initial conditions for x vector

    [t,X] = ode15s(@(t,x) sys(t,x,I),t_span,x0);
    X = X.*(10^(12));
    x_inf = X';
    
    Par(1,i) = Par(1,i)*100; %Increase parameter 1 order of magnitude
    Rxt = (Par(1,12)*Par(1,8)*1000/Av)*(Par(1,6)*Par(1,2)/Par(1,7)); %RNAP concentration (mol/gDW)
    Rlt = (Par(1,24)*Par(1,8)*1000/Av)*(Par(1,6)*Par(1,2)/Par(1,7)); %Ribosomes concentration (mol/gDW)
    Dx = log(2)/Par(1,4); %Degradation rate of mRNA
    Dl = log(2)/Par(1,5); %Degradation rate of proteins
    mu = log(2)/Par(1,1); %Dilution factor

    Sxp = (1.04*Par(1,11)*Par(1,6)*Par(1,2)/Par(1,7))*10^(-6); %Saturation constant (mol/gDW)
    Ke1 = Par(1,9)/Par(1,13); %Elongation rate for mRNA1
    Ke2 = Par(1,9)/Par(1,14); %Elongation rate for mRNA2
    Ke3 = Par(1,9)/Par(1,15); %Elongation rate for mRNA3
    tx1 = Ke1/Par(1,11); %Tau for mRNA1
    tx2 = Ke2/Par(1,11); %Tau for mRNA2
    tx3 = Ke3/Par(1,11); %Tau for mRNA3
    MW1 = Par(1,13)*607.4+157.9; %Molecular Weight of gene 1
    MW2 = Par(1,14)*607.4+157.9; %Molecular Weight of gene 2
    MW3 = Par(1,15)*607.4+157.9; %Molecular Weight of gene 3
    Gp1 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW1)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 1
    Gp2 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW2)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 2
    Gp3 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW3)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 3

    rx1 = Ke1*Rxt*(Gp1/(Sxp*tx1+Gp1*tx1+Gp1)); %Transcription rate for gene 1
    rx2 = Ke2*Rxt*(Gp2/(Sxp*tx2+Gp2*tx2+Gp2)); %Transcription rate for gene 2
    rx3 = Ke3*Rxt*(Gp3/(Sxp*tx3+Gp3*tx3+Gp3)); %Transcription rate for gene 3

    Kl1 = Par(1,10)/Par(1,16); %Translation rate protein 1
    Kl2 = Par(1,10)/Par(1,17); %Translation rate protein 2
    Kl3 = Par(1,10)/Par(1,18); %Translation rate protein 3
    tl1 = Kl1/Par(1,11); %Tau for protein 1
    tl2 = Kl2/Par(1,11); %Tau for protein 1
    tl3 = Kl3/Par(1,11); %Tau for protein 1
    
    wI1 = Par(1,19);
    w11 = Par(1,20);
    w12 = Par(1,21);
    w13 = Par(1,22);
    w23 = Par(1,23);

    %Initial conditions
    t_i = 0;
    t_f = 20; %Final time (min)
    step = 1;
    t_span = t_i:step:t_f; %Time vector (min)
    [m,n] = size(t_span); %Size of time
    I = zeros(n+1,1);
    %I(1:160,1) = 10*10^(-3); %Inducer initial concentration (mol/gDW)
    x0 = [0;0;0;0;0;0]; %Initial conditions for x vector

    [t,X] = ode15s(@(t,x) sys(t,x,I),t_span,x0);
    X = X.*(10^(12));
    x_sup = X';
    
    for j=1:6
        s1(j+z,:) = (x_sup(j,:)-x_inf(j,:))/20;
    end
    
    
    Par(1,i) = tempPar;
    z = z+6;

end


%% Case 2 - Early inducer

%Let's loop!
[a,stop] = size(Par);
z = 0; %For me

for i=1:stop

    %Compound Parameters
    tempPar = Par(1,i);
    Par(1,i) = Par(1,i)/10; %Decrease parameter 1 order of magnitude
    Rxt = (Par(1,12)*Par(1,8)*1000/Av)*(Par(1,6)*Par(1,2)/Par(1,7)); %RNAP concentration (mol/gDW)
    Rlt = (Par(1,24)*Par(1,8)*1000/Av)*(Par(1,6)*Par(1,2)/Par(1,7)); %Ribosomes concentration (mol/gDW)
    Dx = log(2)/Par(1,4); %Degradation rate of mRNA
    Dl = log(2)/Par(1,5); %Degradation rate of proteins
    mu = log(2)/Par(1,1); %Dilution factor

    Sxp = (1.04*Par(1,11)*Par(1,6)*Par(1,2)/Par(1,7))*10^(-6); %Saturation constant (mol/gDW)
    Ke1 = Par(1,9)/Par(1,13); %Elongation rate for mRNA1
    Ke2 = Par(1,9)/Par(1,14); %Elongation rate for mRNA2
    Ke3 = Par(1,9)/Par(1,15); %Elongation rate for mRNA3
    tx1 = Ke1/Par(1,11); %Tau for mRNA1
    tx2 = Ke2/Par(1,11); %Tau for mRNA2
    tx3 = Ke3/Par(1,11); %Tau for mRNA3
    MW1 = Par(1,13)*607.4+157.9; %Molecular Weight of gene 1
    MW2 = Par(1,14)*607.4+157.9; %Molecular Weight of gene 2
    MW3 = Par(1,15)*607.4+157.9; %Molecular Weight of gene 3
    Gp1 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW1)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 1
    Gp2 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW2)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 2
    Gp3 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW3)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 3

    rx1 = Ke1*Rxt*(Gp1/(Sxp*tx1+Gp1*tx1+Gp1)); %Transcription rate for gene 1
    rx2 = Ke2*Rxt*(Gp2/(Sxp*tx2+Gp2*tx2+Gp2)); %Transcription rate for gene 2
    rx3 = Ke3*Rxt*(Gp3/(Sxp*tx3+Gp3*tx3+Gp3)); %Transcription rate for gene 3

    Kl1 = Par(1,10)/Par(1,16); %Translation rate protein 1
    Kl2 = Par(1,10)/Par(1,17); %Translation rate protein 2
    Kl3 = Par(1,10)/Par(1,18); %Translation rate protein 3
    tl1 = Kl1/Par(1,11); %Tau for protein 1
    tl2 = Kl2/Par(1,11); %Tau for protein 1
    tl3 = Kl3/Par(1,11); %Tau for protein 1
    
    wI1 = Par(1,19);
    w11 = Par(1,20);
    w12 = Par(1,21);
    w13 = Par(1,22);
    w23 = Par(1,23);

    %Initial conditions
    t_i = 0;
    t_f = 20; %Final time (min)
    step = 1;
    t_span = t_i:step:t_f; %Time vector (min)
    [m,n] = size(t_span); %Size of time
    I = zeros(n+1,1);
    I(:,1) = 10*10^(-3); %Inducer initial concentration (mol/gDW)
    x0 = [0;0;0;0;0;0]; %Initial conditions for x vector

    [t,X] = ode15s(@(t,x) sys(t,x,I),t_span,x0);
    X = X.*(10^(12));
    x_inf = X';
    
    Par(1,i) = Par(1,i)*100; %Increase parameter 1 order of magnitude
    Rxt = (Par(1,12)*Par(1,8)*1000/Av)*(Par(1,6)*Par(1,2)/Par(1,7)); %RNAP concentration (mol/gDW)
    Rlt = (Par(1,24)*Par(1,8)*1000/Av)*(Par(1,6)*Par(1,2)/Par(1,7)); %Ribosomes concentration (mol/gDW)
    Dx = log(2)/Par(1,4); %Degradation rate of mRNA
    Dl = log(2)/Par(1,5); %Degradation rate of proteins
    mu = log(2)/Par(1,1); %Dilution factor

    Sxp = (1.04*Par(1,11)*Par(1,6)*Par(1,2)/Par(1,7))*10^(-6); %Saturation constant (mol/gDW)
    Ke1 = Par(1,9)/Par(1,13); %Elongation rate for mRNA1
    Ke2 = Par(1,9)/Par(1,14); %Elongation rate for mRNA2
    Ke3 = Par(1,9)/Par(1,15); %Elongation rate for mRNA3
    tx1 = Ke1/Par(1,11); %Tau for mRNA1
    tx2 = Ke2/Par(1,11); %Tau for mRNA2
    tx3 = Ke3/Par(1,11); %Tau for mRNA3
    MW1 = Par(1,13)*607.4+157.9; %Molecular Weight of gene 1
    MW2 = Par(1,14)*607.4+157.9; %Molecular Weight of gene 2
    MW3 = Par(1,15)*607.4+157.9; %Molecular Weight of gene 3
    Gp1 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW1)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 1
    Gp2 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW2)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 2
    Gp3 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW3)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 3

    rx1 = Ke1*Rxt*(Gp1/(Sxp*tx1+Gp1*tx1+Gp1)); %Transcription rate for gene 1
    rx2 = Ke2*Rxt*(Gp2/(Sxp*tx2+Gp2*tx2+Gp2)); %Transcription rate for gene 2
    rx3 = Ke3*Rxt*(Gp3/(Sxp*tx3+Gp3*tx3+Gp3)); %Transcription rate for gene 3

    Kl1 = Par(1,10)/Par(1,16); %Translation rate protein 1
    Kl2 = Par(1,10)/Par(1,17); %Translation rate protein 2
    Kl3 = Par(1,10)/Par(1,18); %Translation rate protein 3
    tl1 = Kl1/Par(1,11); %Tau for protein 1
    tl2 = Kl2/Par(1,11); %Tau for protein 1
    tl3 = Kl3/Par(1,11); %Tau for protein 1
    
    wI1 = Par(1,19);
    w11 = Par(1,20);
    w12 = Par(1,21);
    w13 = Par(1,22);
    w23 = Par(1,23);

    %Initial conditions
    t_i = 0;
    t_f = 20; %Final time (min)
    step = 1;
    t_span = t_i:step:t_f; %Time vector (min)
    [m,n] = size(t_span); %Size of time
    I = zeros(n+1,1);
    I(:,1) = 10*10^(-3); %Inducer initial concentration (mol/gDW)
    x0 = [0;0;0;0;0;0]; %Initial conditions for x vector

    [t,X] = ode15s(@(t,x) sys(t,x,I),t_span,x0);
    X = X.*(10^(12));
    x_sup = X';
    
    for j=1:6
        s2(j+z,:) = (x_sup(j,:)-x_inf(j,:))/20;
    end
    
    
    Par(1,i) = tempPar;
    z = z+6;

end

%% Case 3 - Late inducer

%Let's loop!
[a,stop] = size(Par);
z = 0; %For me

for i=1:stop

    %Compound Parameters
    tempPar = Par(1,i);
    Par(1,i) = Par(1,i)/10; %Decrease parameter 1 order of magnitude
    Rxt = (Par(1,12)*Par(1,8)*1000/Av)*(Par(1,6)*Par(1,2)/Par(1,7)); %RNAP concentration (mol/gDW)
    Rlt = (Par(1,24)*Par(1,8)*1000/Av)*(Par(1,6)*Par(1,2)/Par(1,7)); %Ribosomes concentration (mol/gDW)
    Dx = log(2)/Par(1,4); %Degradation rate of mRNA
    Dl = log(2)/Par(1,5); %Degradation rate of proteins
    mu = log(2)/Par(1,1); %Dilution factor

    Sxp = (1.04*Par(1,11)*Par(1,6)*Par(1,2)/Par(1,7))*10^(-6); %Saturation constant (mol/gDW)
    Ke1 = Par(1,9)/Par(1,13); %Elongation rate for mRNA1
    Ke2 = Par(1,9)/Par(1,14); %Elongation rate for mRNA2
    Ke3 = Par(1,9)/Par(1,15); %Elongation rate for mRNA3
    tx1 = Ke1/Par(1,11); %Tau for mRNA1
    tx2 = Ke2/Par(1,11); %Tau for mRNA2
    tx3 = Ke3/Par(1,11); %Tau for mRNA3
    MW1 = Par(1,13)*607.4+157.9; %Molecular Weight of gene 1
    MW2 = Par(1,14)*607.4+157.9; %Molecular Weight of gene 2
    MW3 = Par(1,15)*607.4+157.9; %Molecular Weight of gene 3
    Gp1 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW1)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 1
    Gp2 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW2)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 2
    Gp3 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW3)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 3

    rx1 = Ke1*Rxt*(Gp1/(Sxp*tx1+Gp1*tx1+Gp1)); %Transcription rate for gene 1
    rx2 = Ke2*Rxt*(Gp2/(Sxp*tx2+Gp2*tx2+Gp2)); %Transcription rate for gene 2
    rx3 = Ke3*Rxt*(Gp3/(Sxp*tx3+Gp3*tx3+Gp3)); %Transcription rate for gene 3

    Kl1 = Par(1,10)/Par(1,16); %Translation rate protein 1
    Kl2 = Par(1,10)/Par(1,17); %Translation rate protein 2
    Kl3 = Par(1,10)/Par(1,18); %Translation rate protein 3
    tl1 = Kl1/Par(1,11); %Tau for protein 1
    tl2 = Kl2/Par(1,11); %Tau for protein 1
    tl3 = Kl3/Par(1,11); %Tau for protein 1
    
    wI1 = Par(1,19);
    w11 = Par(1,20);
    w12 = Par(1,21);
    w13 = Par(1,22);
    w23 = Par(1,23);

    %Initial conditions
    t_i = 0;
    t_f = 300; %Final time (min)
    step = 1;
    t_span = t_i:step:t_f; %Time vector (min)
    [m,n] = size(t_span); %Size of time
    I = zeros(n+1,1);
    I(60:n,1) = 10*10^(-3); %Inducer initial concentration (mol/gDW)
    x0 = [0;0;0;0;0;0]; %Initial conditions for x vector

    [t,X] = ode15s(@(t,x) sys(t,x,I),t_span,x0);
    X = X.*(10^(12));
    x_inf = X';
    
    Par(1,i) = Par(1,i)*100; %Increase parameter 1 order of magnitude
    Rxt = (Par(1,12)*Par(1,8)*1000/Av)*(Par(1,6)*Par(1,2)/Par(1,7)); %RNAP concentration (mol/gDW)
    Rlt = (Par(1,24)*Par(1,8)*1000/Av)*(Par(1,6)*Par(1,2)/Par(1,7)); %Ribosomes concentration (mol/gDW)
    Dx = log(2)/Par(1,4); %Degradation rate of mRNA
    Dl = log(2)/Par(1,5); %Degradation rate of proteins
    mu = log(2)/Par(1,1); %Dilution factor

    Sxp = (1.04*Par(1,11)*Par(1,6)*Par(1,2)/Par(1,7))*10^(-6); %Saturation constant (mol/gDW)
    Ke1 = Par(1,9)/Par(1,13); %Elongation rate for mRNA1
    Ke2 = Par(1,9)/Par(1,14); %Elongation rate for mRNA2
    Ke3 = Par(1,9)/Par(1,15); %Elongation rate for mRNA3
    tx1 = Ke1/Par(1,11); %Tau for mRNA1
    tx2 = Ke2/Par(1,11); %Tau for mRNA2
    tx3 = Ke3/Par(1,11); %Tau for mRNA3
    MW1 = Par(1,13)*607.4+157.9; %Molecular Weight of gene 1
    MW2 = Par(1,14)*607.4+157.9; %Molecular Weight of gene 2
    MW3 = Par(1,15)*607.4+157.9; %Molecular Weight of gene 3
    Gp1 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW1)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 1
    Gp2 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW2)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 2
    Gp3 = (Par(1,8)*1000*Par(1,3)*50*10^(-9)/MW3)*(Par(1,6)*Par(1,2)/Par(1,7)); %Gp for gene 3

    rx1 = Ke1*Rxt*(Gp1/(Sxp*tx1+Gp1*tx1+Gp1)); %Transcription rate for gene 1
    rx2 = Ke2*Rxt*(Gp2/(Sxp*tx2+Gp2*tx2+Gp2)); %Transcription rate for gene 2
    rx3 = Ke3*Rxt*(Gp3/(Sxp*tx3+Gp3*tx3+Gp3)); %Transcription rate for gene 3

    Kl1 = Par(1,10)/Par(1,16); %Translation rate protein 1
    Kl2 = Par(1,10)/Par(1,17); %Translation rate protein 2
    Kl3 = Par(1,10)/Par(1,18); %Translation rate protein 3
    tl1 = Kl1/Par(1,11); %Tau for protein 1
    tl2 = Kl2/Par(1,11); %Tau for protein 1
    tl3 = Kl3/Par(1,11); %Tau for protein 1
    
    wI1 = Par(1,19);
    w11 = Par(1,20);
    w12 = Par(1,21);
    w13 = Par(1,22);
    w23 = Par(1,23);

    %Initial conditions
    t_i = 0;
    t_f = 300; %Final time (min)
    step = 1;
    t_span = t_i:step:t_f; %Time vector (min)
    [m,n] = size(t_span); %Size of time
    I = zeros(n+1,1);
    I(60:n,1) = 10*10^(-3); %Inducer initial concentration (mol/gDW)
    x0 = [0;0;0;0;0;0]; %Initial conditions for x vector

    [t,X] = ode15s(@(t,x) sys(t,x,I),t_span,x0);
    X = X.*(10^(12));
    x_sup = X';
    
    for j=1:6
        s3_temp(j+z,:) = (x_sup(j,:)-x_inf(j,:))/20;
    end
    
    
    Par(1,i) = tempPar;
    z = z+6;

end

s3 = zeros(144,21);
[a,b] = size(s3_temp);
s3 = s3_temp(:,b-20:b);
