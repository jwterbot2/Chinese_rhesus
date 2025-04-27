//Number of population samples (demes)
3
//Population effective sizes (number of genomes)
Nbrv
Ntch
Nml
//Sample sizes (haploid size of the sample, sampling time in generations back)
10
10
138
//Growth rates
0
0
0
//Number of migration matrices : If 0 : No migration between demes
0
//Historical event: time, source, sink, migrants, new deme size, new growth rate, new migration matrix
4 historical events
$T1 2 1 1 NA1 0 0 nomig absoluteResize
$T1 2 2 0 0 0 0 // kill ml
$T2 1 0 1 NA2 0 0 absoluteResize
$T2 1 1 0 0 0 0 // kill tch
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block: data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 1e-8 1.08e-8 OUTEXP