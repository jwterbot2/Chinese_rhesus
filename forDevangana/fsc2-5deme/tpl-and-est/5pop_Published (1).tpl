//Number of population samples (demes)
5 
//Population effective sizes (number of genes)
Ntc
Nli
Nbr
Nla
Nmu
//Sample sizes (haploid size of the sample, sampling time in generations back)
10
56
10
62
20
//Growth rates (negative growth =population expansion, because its measured back in time)
0
0
0
0
0
//Number of migration matrices : If 0 : No migration between demes
5
//Migration matrix 0
0 	Mtc2li 	0 	Mtc2la 	0
Mli2tc 	0 	Mli2br	Mli2la 	Mli2mu
0	Mbr2li 	0 	0 	0
Mla2tc 	Mla2li  0 	0 	Mla2mu
0 	Mmu2li  0 	Mmu2la 	0
//Migration matrix 1
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
//Migration matrix 2
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
//Migration matrix 3
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
//Migration matrix 4
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
//Historical event: time, source, sink, migrants, new deme size, new growth rate, new migration matrix
8 historical events
t1 0 1 1 NA1 0 1 absoluteResize
t1 0 0 0 0 0 1 //kill tc
t2 1 2 1 NA2 0 2 absoluteResize
t2 1 1 0 0 0 2 //kill li
t3 2 3 1 NA3 0 3 absoluteResize
t3 2 2 0 0 0 3 //kill br
t4 3 4 1 NA4 0 4 absoluteResize
t4 3 3 0 0 0 4 //kill la
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 1.00e-8 OUTEXP
