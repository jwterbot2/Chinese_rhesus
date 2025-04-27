## Usage: python pubrhe.msprimeSimulations_v0.03.py ModelName ReplicateNumber

#import libraries to be used.
import sys
import os
import msprime
import tskit

DEBUG=True;

model=sys.argv[1];
rep=int(sys.argv[2]);
nrep = 1;
samp = 79;

#Define model based on command line argument
if DEBUG: print("Setting Demography based on " + model);
demog = msprime.Demography()

recombRate = 1e-8;

if ( model == "fsc2-3" ):
	NCur = 126626/2;
	N_T1 = 374873/2;
	N_T2 = 219439/2;
	N_T3 = 65196/2;
	
	T1 = 562;
	T2 = 31130;
	T3 = 68368;
	
	demog.add_population(name="Rhesus", description="Single Rhesus Population", initial_size=NCur, growth_rate=0);
	demog.add_population_parameters_change(time=T1, initial_size=N_T1, population="Rhesus");
	demog.add_population_parameters_change(time=T2, initial_size=N_T2, population="Rhesus");
	demog.add_population_parameters_change(time=T3, initial_size=N_T3, population="Rhesus");

elif ( model == "fsc2-4" ):
	NCur = 160646/2;
	N_T1 = 640374/2;
	N_T2 = 275692/2;
	N_T3 = 393753/2;
	N_T4 = 60170/2;
	
	T1 = 1707;
	T2 = 9494;
	T3 = 42956;
	T4 = 67421;
	
	demog.add_population(name="Rhesus", description="Single Rhesus Population", initial_size=NCur, growth_rate=0);
	demog.add_population_parameters_change(time=T1, initial_size=N_T1, population="Rhesus");
	demog.add_population_parameters_change(time=T2, initial_size=N_T2, population="Rhesus");
	demog.add_population_parameters_change(time=T3, initial_size=N_T3, population="Rhesus");
	demog.add_population_parameters_change(time=T4, initial_size=N_T4, population="Rhesus");
	
else:
	print("Model " + model + " does not exist. Exiting.");
	exit();

demog.sort_events();

#Make Tree(s)
if DEBUG: print("Making Tree");
treeRep = msprime.sim_ancestry(samp, demography=demog, recombination_rate = recombRate, sequence_length = 223000000, num_replicates=nrep);
if DEBUG: print(treeRep);

#Mutate Tree(s)
ntree=0
for tree in treeRep:
	#Mutate each tree generated individually
	if DEBUG: print(tree);
	if DEBUG: print("Mutating Tree");
	mutated_tree = msprime.sim_mutations(tree, rate=1.08e-8);
	if DEBUG: print(mutated_tree);
	
	#Output vcf to scratch folder.
	out = open("[/path/to/output/folder/file-header]." + model + "." + str(rep) + "-" + str(ntree) + ".SNPs.vcf", "w");
	tskit.TreeSequence.write_vcf(mutated_tree, out);
	ntree = ntree + 1;

if DEBUG: print("done.")
