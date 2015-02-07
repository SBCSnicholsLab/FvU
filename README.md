# FvU
Fused versus Unfused - Modelling the behaviour of chromosomal clines

Veltsos et al. 2008 modelled two things. (1) the selection conditions for which a neo-X/neo-Y sex determination system can become fixed in a population with X0-system, and (2) which conditions favour the spread of the neo-sex chromosome system into other populations. This repo recapitulates and extends (2).


## Files
###**FvsU.R**
recapitulates what was done in Veltsos et al. 2008 to model the spread of a sex determination system among populations (along a gradient).

###**FvU3D_apo_glm.R**
extends **FvsU.R** from looking at a gradient to include an area. It iterates dispersal, calculation of gamete ratios, and zygote formation by randomly drawing from the gamete pool until there is no change in the distribution of the fused X and the Y-chromosome. It will not model more than 10000 generations, though. Every tenth generation, two glms with logit link-function are being fitted, one for the fused X and one for the Y-chromosome, to infer centre and width of their clines. Then, a list containing the allele matrix (genotypes\$g),  the model coefficients and inferred cline centres (genotypes\$coefs), and the information whether the cline centres transcend (or previously have transcended) the model margins (genotypes\$trans) will be saved using the save() function. The results can be visualised using **get_coefs_apo.R**

The script is meant to be run from the command line (rather than interactively in an R editor). It requires nine arguments. The first seven arguments are: (1) selection against the b allele in males, (2) dominance of the a allele in males, (3) selection against the a allele in females, (4) dominance of the a allele in females, (5) heterozygote disadvantage, (6) dosage compensation, and (7) sex chromosome coadaption, each as used in Veltsos et al. 2008. (8) defines what fraction of each population disperses into each cardinal direction (i.e. total dispersal is (8) times four). (9) must be a one-digit integer. It serves to tell apart replicates with otherwise identical parameters.

USAGE: Rscript FvU3D_apo_glm.R \<sm\> \<dm\> \<sf\> \<df\> \<sh\> \<sd\> \<sc\> \<disp\> \<no of rep\> 

(1-8) each must be larger or equal 0 and smaller or equal 1 with up to three decimal places. (9) must be a one-digit integer.

For a quick test-run "Rscript FvU3D_apo_glm.R 0.9 0.5 0.9 0.5 0.01 0 0 0.03 1" can be used. On a recent machine it should take five minute or less.

###**get_coefs_apo.R**
parses the output files of **FvU3D_apo_glm.R** and extracts genotypes\$coefs and genotypes$trans (coefficients of the glms, cline centres, and whether the cline centres have transcended). These information will be written to a file and a PNG image showing the movement of the clines and their widths will be produced.

**get_coefs_apo.R** will list all files beginning with 'sm' in a directory. It will then group files according to the selection parameters used and extract the data as explained above for each generation. Errors will be raised if in a directory (1) there are files beginning with 'sm' that are not the output of **FvU3D_apo_glm.R**, or (2) there are generations missing.

USAGE: Rscript FvU3D_apo_glm.R

The script requires a subfolder called "dances/" (in which the output will be saved). This folder should be empty.

###**gametes.csv**
TAB-separated table of gamete ratios obtained from ancestral genotypes


###**genotypes.csv**
TAB-separated table of zygotic genotypes resulting from all possible combinations of gametes
