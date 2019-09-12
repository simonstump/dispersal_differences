%This analyzes the dispersal data by Muller-Landau et al. (2008), to estimate the negative impact that seed dispersal limitation is having on the stabilizing mechanism.

clc, clear

%These are the alpha values listed in Muller-Landau et al. (2008).
alpha=[4.9	% Anacardium excelsum
3.84	% Beilschmiedia pendula
7.06	% Calophyllum longifolium
5.6	% Chrysophyllum cainito
6.51	% Cupania rufescens
5.73	% Dendropanax arboreus
4.74	% Dipteryx oleifera
4.61	% Drypetes standleyi
6.47	% Guapira standleyana
6.33	% Guatteria dumetorum
5.8	% Jacaranda copaia*
5.67	% Luehea seemannii*
4.96	% Platypodium elegans
4.5	% Platymiscium pinnatum
5.15	% Poulsenia armata
5.2	% Pouteria reticulata
5.19	% Pterocarpus rohrii
4.24	% Quararibea asterolepis
6.02	% Simarouba amara
4.61	% Tabebuia guayacan*
6.3	% Tabebuia rosea*
8.24	% Terminalia amazonia
4.8	% Tetragastris panamensis
8.75	% Poulsenia armata
3.93	% Trichilia tuberculata
8.4];	% Zanthoxylum ekmanii


%this is the number of sites a distance (x-1) away.
this=[.25,(1:5)]*4;

%Here we calculate the dispersal kernel for all species.
for i=1:length(alpha)
   disp(i,:)=dispersal_kern(exp(alpha(i)),10,5);
end


%The 7th row is long-distance dispersal; we remove it.
disp(:,7)=[]

%dd is value of \tilde{d_j d_j^P} we report in Table A2.
dd=(disp.^2)*this'

mean(dd)