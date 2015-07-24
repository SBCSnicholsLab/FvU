

# FvU
Fused versus Unfused - Modelling the behaviour of chromosomal clines with sexually antagonistic selection acting

Veltsos et al. 2008 modelled two things. (1) the selection conditions for which a neo-X/neo-Y sex determination system can become fixed in a population with X0-system, and (2) which conditions favour the spread of the neo-sex chromosome system into other populations. This repo recapitulates and extends (2).


## Files

###FvU3D_sizelayer.py
Models hybrid zone with directional gene flow due to a population size gradient.

####Usage
python FvU3D.py \<sm\> \<dm\> \<sf\> \<df\> \<sh\> \<sd\> \<sc\> \<disp\> \<grad\> \<run no\>  
EXAMPLE: python FvU3D_sizelayer.py .01 1 .01 1 0 0 0 .03 4 1  
This will use a gradient of 6% effectively stopping the clines' forward-movement.
The parameters are:  
**sm** (selection against the Y in males),  
**dm** (dominance of the Y in males),  
**sf** (selection against the Y in males),  
**df** (dominance of the Y in males),  
**sh** (selection against heterozygote females, heterozygote disadvantage),  
**sd** (dosage compensation),  
**sc** (sex chromosome co-adaption),  
**disp** (fraction of a deme that disperse to each direction i.e. the actual dispersal is dispx4),  
**grad** (population size gradient to be used; possible options are 0-8 for	speepnesses of 0, 1, 3, 5, 6, 7, 9, 11, and 13%), and  
**run** (a number that will be part of the output file names)

####Outputs
1. one tab-delimited TXT file with GLM parameters (intercept, slope, and centre of the Y cline plus intercept, slope, and centre of the fused-X cline),  
2. every ten generations: one GZ-compressed pickle file of the genotype dictionary, and   
3. every ten generations: one GZ-compressed TXT table with hits and misses for F and Y (can be imported into R to visualise it with rgl or to do other funny things)  
The standard array shape is 10x100 populations with cline inflection points between 60 and 61.  

####Requirements
 
Requires the python libraries numpy, pickle, and statsmodels.

###FvU3D_sizelayer_wind.py

Very similar to FvU3D_sizelayer.py, but there is no population size gradient. The simulation will be run on a 10x120 array with 40 individuals in each generation. Directional gene flow is implemented as a fraction (dispmod) that is taken away from rightward dispersal and added to leftward dispersal.

####Usage

python FvU3D.py \<sm\> \<dm\> \<sf\> \<df\> \<sh\> \<sd\> \<sc\> \<disp\> \<dispmod\> \<run no\>  
EXAMPLE: python FvU3D_sizelayer.py .01 1 .01 1 0 0 0 .03 4 1  
This will use a gradient of 6% effectively stopping the clines' forward-movement.
The parameters are:  
**sm** (selection against the Y in males),  
**dm** (dominance of the Y in males),  
**sf** (selection against the Y in males),  
**df** (dominance of the Y in males),  
**sh** (selection against heterozygote females, heterozygote disadvantage),  
**sd** (dosage compensation),  
**sc** (sex chromosome co-adaption),  
**disp** (fraction of a deme that disperse to each direction i.e. the actual dispersal is dispx4),  
**dispmod** (percentage that is removed from rightward dispersal and added to leftward dispersal), and  
**run** (a number that will be part of the output file names)

####Outputs
1. one tab-delimited TXT file with GLM parameters (intercept, slope, and centre of the Y cline plus intercept, slope, and centre of the fused-X cline),  
2. every ten generations: one GZ-compressed pickle file of the genotype dictionary, and   
3. every ten generations: one GZ-compressed TXT table with hits and misses for F and Y (can be imported into R to visualise it with rgl or to do other funny things)  
The standard array shape is 10x100 populations with cline inflection points between 60 and 61.  

####Requirements
 
Requires the python libraries numpy, pickle, and statsmodels.
