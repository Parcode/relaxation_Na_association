# relaxation_Na_association
Analyzing Na-ion relaxation and diffusion data to retrieve ion-association to surfactants.

The use is discussed in

manuscript , (please cite if used): 

Pär Håkansson, Pau Mayorga Delgado, Anne Selent, Ritu Ghanghas, Ilari Ainasoja, Sanna Komulainen, Jiří Mareš, Perttu Lantto, Nønne L. Prisle and Ville-Veikko Telkki

Association and dissociation of Na+ between bulk, cluster and micelle sites in aqueous sodium decanoate solutions elucidated by 23Na NMR relaxation experiments and quadrupolar relaxation modelling, submitted for publication, 2026. (will update with link soon)

The main script is FITSCRIPT_Na_v14.m that takes inputparameters from MD simulation, experimental data and runs Markov chain Monte Carlo (MCMC)

to determine the ion association and some additional model parameters. The actual model is assembled in R1_R2model_LM_v14.m. Similar for v9-model. We will keep this code relevant for above manuscript. Currently code is a bit hard to read, if needed we will update with some more comments and remove some variables not used.
