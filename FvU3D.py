#!/usr/bin/python


#FvU3D.py currently takes three pramameters all of which have to be supplied when calling the programme.
#These are sm (selection against the)
#EXAMPLE: python FvU3D.py .05 1 .1 1 0 0 0 .03 0 1



import numpy as np
import sys
import pickle
import gzip
import statsmodels.api as stm
#import time



#The 11 genotypes are: 'Ubb', 'Fa', 'Uab', 'Fb', 'Uaa', 'UUbb', 'FF', 'FUb', 'UUab', 'FUa', and 'UUaa' (0-10)
#The 8 gametes (5 male + 3 female) are: 'mF', 'ma', 'mb', 'mUa', 'mUb', 'fF', 'fUa', and 'fUb' (0-7) 

#number of Fs containd in genotypes
Fgenotypes = np.array([0,1,0,1,0, 0,2,1,0,1,0])
#number of potential Fs (depending on sex)
Fpot = np.array([1,1,1,1,1, 2,2,2,2,2,2])
Fabs = np.array([1,0,1,0,1, 2,0,1,2,1,2])
#number of Ys containd in genotypes
Ygenotypes = np.array([0,1,1,0,2, 0,0,0,1,1,2])
#number of potential Ys (depending on wheter and how many Fs are present)
Ypot = np.array([2,1,2,1,2, 2,0,1,2,1,2])
Yabs = np.array([2,0,1,1,0, 2,0,1,1,0,0])

def Fratio(x):
	return np.sum(x*Fgenotypes, axis=2)/np.sum(x*Fpot, axis=2)



def YFratio(x):
	return np.sum(x*Ygenotypes, axis=2)/np.sum(x*Ypot, axis=2)


def Fhits(x):
	return np.sum(x*Fgenotypes, axis=2)

def Fmisses(x):
	return np.sum(x*Fabs, axis=2)

def Yhits(x):
	return np.sum(x*Ygenotypes, axis=2)

def Ymisses(x):
	return np.sum(x*Yabs, axis=2)

genconv=np.array([6,9,7,1,4,2,3,2,0,9,10,8,7,8,5]).reshape(5,3)

gamconv=np.array([0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 2, 0, 0, 0,\
	0, 0, 0, 2, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0,\
	0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2,\
	0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 4, 0, 0, 0, 0, 0, 4, 0, 2, 2, 0, 0]).reshape(8,11)


def params():
	
	if len(sys.argv) <> 11:
		print 'Wrong number of parameters supplied. Exit.'
		exit()

	global 	sm, dm, sf, df, sh, sd, sc, disp, grad, run

	sm = float(sys.argv[1])
	dm = float(sys.argv[2])
	sf = float(sys.argv[3])
	df = float(sys.argv[4])
	sh = float(sys.argv[5])
	sd = float(sys.argv[6])
	sc = float(sys.argv[7])
	
	disp = float(sys.argv[8])
	grad = int(sys.argv[9])
	run = int(sys.argv[10])


	#vector of selection acting on genotypes (SELECTion vectOR)
	global selector
	selector=np.array([[1-sm, 1, (1-sm*dm)*(1-.5*sd)*(1-.5*sc), (1-sm)*(1-.5*sd)*(1-sc),\
		(1-.5*sd)*(1-sc), 1, 1, (1-sh)*(1-.5*sd)*(1-.5*sc), (1-sf*df)*(1-.5*sd)*(1-.5*sc),\
		(1-sf*df)*(1-sh)*(1-.5*sc), (1-sf)*(1-sd)*(1-sc)]])


#Intialises a 3D genotype array ['g']
#(plus a counter['c'], an abundance matrix ['abu'] and other stuff
#that may become important later).
#Returns a dictionary containing all that stuff
def initialise(wi=10, leng=120, popsize=40):

	a=np.zeros([wi, leng, 11])

	abu=np.repeat(popsize , wi*leng).reshape(wi, leng)


	a[:, 0:20, [1,6]]=.5	#fused
	a[:, 20:, [0,5]]=.5	#unfused

	a=a*abu[:,:,np.newaxis]

	c=0

	trans=np.array([0,0])

	print "Matrix initialised with %dx%d populations of %d individuals each." % (wi, leng, popsize)
	return {'g': a, 'abu': abu, 'c': c, 'fm': a[:,[0],:], 'um': a[:,[leng-1],:], 'trans':trans}




#Get Gamete Ratios, takes one 1D np array (genotype vector) as argument and
#returns a 1D 8-item gamete vector
def ggr(x):	

	return np.sum(x*gamconv, axis=1)

#GAmete SAample
#samples gametes according to gamete ratios, does sex and
#returns the new generations genotypes
def gasa(x):	

	#Male gAMETES
	mametes=np.repeat(np.arange(5),np.random.multinomial(40, x[0:5]/sum(x[0:5]), 1).flatten())
	np.random.shuffle(mametes)

	#Female gAMETES
	fametes=np.repeat(np.arange(3),np.random.multinomial(40, x[5:8]/sum(x[5:8]), 1).flatten())
	np.random.shuffle(fametes)

	##FASTAER SEX NOW - MULTIPLE SELECTION FROM ARRAY BY "ADVANCED" INDEXING
	a = np.bincount(genconv[mametes, fametes])
	#topping upt the array in case not all 11 genotypes were produced
	a = np.hstack((a, np.zeros(11-a.shape[0])))
	return a

#the heart of this programme, the Time Goes By function
def tgb(genot, gens=10):

	#dimensions of the genotype matrix	
	dimensions=np.array(genot['g'].shape)
	
	
	#intgen is the tbg-internal version of the genotype array
	#it is two elements wirde in the firs two dimensions to allow
	#distribution to be modelled by addition of an offset vrsion of the array
	
	intgen=np.zeros(dimensions+np.array([2,2,0]))
	
	intgen[1:(dimensions[0]+1), 1:(dimensions[1]+1), :]=genot['g']

	#DISPersal FACTor is an arry with the same dimensions as the genotype array
	#it is being multiplid with the genotypes to account for margin effects
#	dispfact=np.repeat(1-disp*4, np.prod(dimensions)).reshape(dimensions)
#	dispfact[[0,dimensions[0]-1],:,:]=1-disp*3
#	dispfact[:,[0,dimensions[1]-1],:]=1-disp*3
#	dispfact[[0, 0, dimensions[0]-1, dimensions[0]-1],[0, dimensions[1]-1, 0, dimensions[1]-1],:]=1-disp*2
#	print dispfact

	for i in range(gens):
		#dispersal
		intgen[1:(dimensions[0]+1), 1:(dimensions[1]+1), :]=\
			intgen[1:(dimensions[0]+1), 1:(dimensions[1]+1), :]*(1-4*disp)+\
			intgen[0:(dimensions[0]), 1:(dimensions[1]+1), :]*disp+\
			intgen[2:(dimensions[0]+2), 1:(dimensions[1]+1), :]*disp+\
			intgen[1:(dimensions[0]+1), 0:(dimensions[1]), :]*disp+\
			intgen[1:(dimensions[0]+1), 2:(dimensions[1]+2), :]*disp
		#print intgen[:,18:22,:]
		

		#selection acts
		#(selection vector is being multiplied with the genotype array)
		gametes=intgen[1:(dimensions[0]+1), 1:(dimensions[1]+1), :]*selector


		#gamete ratios, from the genotypes arry every vector along axis 2 (= genotype axis)
		#is being multiplied with the genotype conversion matrix
		#subsequently, the rowsums are being calculated resulting in a vector of gamete ratios
		#"apply_along_axis(function, 2, 3Darray)" is similar to "apply(3D array, c(1,2), function)" in R
		gametes=np.apply_along_axis(ggr, 2, gametes)
#		print gametes[:,19:23,:]

		#back to genotypes. this where drift is acting, cf. gasa function
		intgen[1:(dimensions[0]+1), 1:(dimensions[1]+1), :]=np.apply_along_axis(gasa, 2, gametes)
		genot['c']+=1
		#print genot['c']
	genot['g']=intgen[1:(dimensions[0]+1), 1:(dimensions[1]+1), :]

	return genot
		
		


def __main__():

	#get global variable "selector"
	params()
	
	print 'Selection vector', selector


	genotypes=initialise()

	paramtab = open('sm%1.3f_dm%1.3f_sf%1.3f_df%1.3f_sh%1.3f_sd%1.3f_sc%1.3f_disp%1.3f_grad%1.3f_run%d_gens%05d.tab'%\
			(sm, dm, sf, df, sh, sd, sc, disp, grad, run, genotypes['c']), 'w')
	paramtab.write('Yint\tYslope\tYcent\tFint\tFslope\tFcent\n')

	for i in range(10):
		genotypes=tgb(genotypes, 10)


		Fh=Fhits(genotypes['g']).flatten()
		Fm=Fmisses(genotypes['g']).flatten()
		Yh=Yhits(genotypes['g']).flatten()
		Ym=Ymisses(genotypes['g']).flatten()

		yvals=stm.tools.add_constant(np.tile(np.arange(genotypes['g'].shape[1]),genotypes['g'].shape[0]))


		genotypes['modelparams']=[]
		glmy=stm.GLM(np.column_stack((Yh, Ym)), yvals, family=stm.families.Binomial())
		glmyp=glmy.fit().params
		genotypes['modelparams'] += [glmyp[0], glmyp[1], -glmyp[0]/glmyp[1]]

		glmf=stm.GLM(np.column_stack((Fh, Fm)), yvals, family=stm.families.Binomial())
		glmfp=glmf.fit().params
		genotypes['modelparams'] += [glmfp[0], glmfp[1], -glmfp[0]/glmfp[1]]

		pickle.dump(genotypes, gzip.open('sm%1.3f_dm%1.3f_sf%1.3f_df%1.3f_sh%1.3f_sd%1.3f_sc%1.3f_disp%1.3f_grad%1.3f_run%d_gens%05d.pyd.gz'%\
			(sm, dm, sf, df, sh, sd, sc, disp, grad, run, genotypes['c']), 'w'))


		out = gzip.open('sm%1.3f_dm%1.3f_sf%1.3f_df%1.3f_sh%1.3f_sd%1.3f_sc%1.3f_disp%1.3f_grad%1.3f_run%d_gens%05d.csv.gz'%\
			(sm, dm, sf, df, sh, sd, sc, disp, grad, run, genotypes['c']), 'w')		
		out.write(str(genotypes['g'].shape)[1:-1]+'\ny\tFhit\tFmis\tYhit\tYmis\n')
		np.savetxt(out, np.column_stack((yvals[:,1], Fh, Fm, Yh, Ym)), "%d", delimiter='\t')
		out.close()

		if genotypes['trans'][0]==0:
			if genotypes['modelparams'][2] < 0 or genotypes['modelparams'][2] > genotypes['g'].shape[1]:
				genotypes['trans'][0]=1

		if genotypes['trans'][1]==0:
			if genotypes['modelparams'][5] < 0 or genotypes['modelparams'][5] > genotypes['g'].shape[1]:
				genotypes['trans'][1]=1

		print "Generation %d, cline centres: %3.2f and %3.2f" % (genotypes['c'], genotypes['modelparams'][2], genotypes['modelparams'][5])
		paramtab.write('\t'.join(str(genotypes['modelparams'])[1:-1].split(', '))+'\n')		

		if 0 not in genotypes['trans']:
			print 'Both variants fixed. Break.'
			break

	paramtab.close()

## Plot using matplotlib
#	from mpl_toolkits.mplot3d import axes3d
#	import matplotlib.pyplot as plt
#	from matplotlib import cm
#
#	fig = plt.figure()
#	ax = fig.add_subplot(111, projection='3d')
#	Y = range(10)
#	X = range(90)
#	X, Y = np.meshgrid(X, Y)
#	Z = Fratio(genotypes['g'])
#	ax.plot_wireframe(X, Y, Z)
#	plt.show()



__main__()

