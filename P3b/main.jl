#Script for Prelim 1 P3b

#Magic function
include("Flux.jl");

#Plot plot plot
using PyPlot

#Parameters
Ll = 308; #Length of protein [aa]
Lx = 924; #Length of gene [nts]
Gp = 5/10^3; #Plasmid concentration [uM]
Cv = 15; #Reaction volume [uL]
Rx = 0.15; #RNAP concentration [uM]
Rl = 1.6; #Ribosomes concentration [uM]
vx = 60; #Elongation rate [nts/s]
vl = 16.5; #Transcription rate [aa/s]
Kx = 0.3; #Saturation constant  translation [uM]
Kl = 57; #Saturation constant Transcription [uM]
tx = 2.7; #Time constant translation
tl = 0.8; #Time constant Transcription
kdx = 8.35/3600; #Degradation rate translation [1/s]
kdl = (9.9*10^(-3))/3600; #Degradation rate Transcription [1/s]

w1 = 0.26;
w2 = 300;
K = 0.3*10^3; #[uM]
n = 1.5;

#Bounds

v1 = 100000/3600; #[uM/s]
v3 = 100000/3600; #[uM/s]
v4 = 100000/3600; #[uM/s]
v6 = 100000/3600; #[uM/s]
b1 = 100000/3600; #[uM/s]
b2 = 100000/3600; #[uM/s]
b3 = 100000/3600; #[uM/s]
b4 = 100000/3600; #[uM/s]
b5 = 100000/3600; #[uM/s]
b6 = 100000/3600; #[uM/s]
b7 = 100000/3600; #[uM/s]
b8 = 100000/3600; #[uM/s]
b9 = 100000/3600; #[uM/s]

#Voigt type model
I_i = 0.0001*10^3; #Inducer initial concentration [uM]
I_f = 10*10^3; #Inducer final concentration [uM]
step = 0.01*10^3; #step size
I = I_i:step:I_f; #Inducer array

f = (I.^n)./(K^n.+I.^n); #Inducer function
v2 = f.*(vx/Lx)*Rx*((Gp)/(Kx*tx+(1+tx)*Gp)); #[uM/s]

mRNA = v2./kdx; #mRNA concentration [uM]

v5 = (vl/Ll)*Rl*((mRNA)./(Kl*tl.+(1+tl).*mRNA)); #[uM/s]

#Stochiometric matrix
#Columns: v1,v2,v3,v4,v5,v6,b1,b2,b3,b4,b5,b6,b7,b8,b9
#rows: G,RNAP,G*,NTP,mRNA,Pi,NMP,rib,rib*,AAtRNA,GTP,tRNA,GDP,Protein,AA,ATP,AMP
stoichiometric_matrix = [[-1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
                         [-1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
                         [1.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
                         [0.0 -1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
                         [0.0 1.0 -1.0 -1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
                         [0.0 2.0 0.0 0.0 2.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0];
                         [0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0];
                         [0.0 0.0 0.0 -1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
                         [0.0 0.0 0.0 1.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
                         [0.0 0.0 0.0 0.0 -1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
                         [0.0 0.0 0.0 0.0 -2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0];
                         [0.0 0.0 0.0 0.0 1.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
                         [0.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0];
                         [0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0];
                         [0.0 0.0 0.0 0.0 0.0 -1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
                         [0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0];
                         [0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0]];


m = size(v5,1);
opv5 = zeros(m,1);

for i=1:m
        #Bounds array
        default_bounds_array = [[0 v1];
                                [v2[i,1] v2[i,1]];
                                [0 v3];
                                [0 v4];
                                [0 v5[i,1]];
                                [0 v6];
                                [-b1 b1];
                                [-b2 b2];
                                [-b3 b3];
                                [-b4 b4];
                                [-b5 b5];
                                [-b6 b6];
                                [-b7 b7];
                                [-b8 b8];
                                [-b9 b9]];
        #Species bounds array
        species_bounds_array = zeros(Float64, 17, 2);
        #Objective array
        objective_coefficient_array = [0.0; 0.0; 0.0; 0.0; -1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0];
        flux, n1, n2, n3, n4, n5 = calculate_optimal_flux_distribution(stoichiometric_matrix,default_bounds_array,species_bounds_array,objective_coefficient_array);
        opv5[i,1] = flux;
end

ProteinC = (abs.(opv5))./kdl; #Protein concentration [uM]
I = I./10^3; #Inducer concentration [mM]

#Time to plot
figure()
semilogx(I[:,1],ProteinC[:,1],color = "black", lw = 2)
xlabel("Inducer [mM]", fontsize = 16)
ylabel("Protein [uM]", fontsize = 16)
savefig("Protein.pdf")
