# Chinese_rhesus

Scripts and associated files for [CITATION]().

The "fastsimcoal2" folder contains the files used to analyze the 5-deme, 3-deme, and 1-deme (single deme) 
models using fastsimcoal2. Each sub-folder contains a folder named "tpl-and-est" containing the .tpl and 
.est files used for parameter estimation, a folder named "sfs-files" containing the site frequency 
spectra for the corresponding number of demes, and a folder named "best-parameters" containing the .par
files used to run simulations of the best models. For the 5-deme folder, there are files pertaining to
the "Liu-et-al_2018" model and the "re-estimate" model. For the 3-deme folder, there are files for the
model with *M. m. brevicaudus* ("\*brvFirst\*") or *M. m. tcheliensis* ("\*tchFirst\*") branching first.
For the 1-deme folder, there are files for 0 ("\*_0\*") to 4 ("\*_4\*") population size change event models.

The "dadi" folder contains the files used to analyze the 1-deme models using dadi. The file 
"dadi_data_prep.py" is a script used to generate a dadi-formatted site-frequency spectrum from an input
.vcf file. The file "dadi_demography.py" is the script used to perform the model parameterization using
dadi described in the manuscript. The final file, "dadi_demography_aux.py" is an example script for 
running dadi analyses on additional 1-deme models not discussed in the manuscript.

The "msprime" folder contains a single file "pubrhe_msprimeSimulations_v0.03.py" which is the script used
to simulate .vcf files under the 3-event or 4-event model parameterized by fastsimcoal2 for recombination
rate comparisons. The script requires as command-line input the model to run (either "fsc2-3" or "fsc2-4")
and the replicate number being run, both of which are included in the output name.

The "demes" folder contains the parameterized models in the "Demes" specification (Gower et al. [2022](https://doi.org/10.1093/genetics/iyac131)).