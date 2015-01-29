
##Adjust wd to location of the CSV tables
setwd('/home/hannes/Dropbox/scripts/FvsU/')



#Selection parameters (cf. Veltsos et al. 2008 Tab. 1)

sm<-0.01         #selection against b in males
dm<-1           #dominance of a in males
sf<-0.02           #selection against a in females
df<-1           #dominance of a in females
sh<-0           #"heterozygote disadvantage"
sd<-0           #"dosage compensation"
sc<-0           #"sex chromosome coadaption"

disp<-0.08      #dispersal in each direction

#vector of selection acting on genotypes (SELECTion vectOR)
selector<-c(1-sm, 1, (1-sm*dm)*(1-.5*sd)*(1-.5*sc), (1-sm)*(1-.5*sd)*(1-sc), (1-.5*sd)*(1-sc),
            1, 1, (1-sh)*(1-.5*sd)*(1-.5*sc), (1-sf*df)*(1-.5*sd)*(1-.5*sc), (1-sf*df)*(1-sh)*(1-.5*sc), (1-sf)*(1-sd)*(1-sc))
names(selector)<-c('Ubb', 'Fa', 'Uab', 'Fb', 'Uaa',
                   'UUbb', 'FF', 'FUb', 'UUab', 'FUa', 'UUaa')

selector

##What gamete ratios (rows) result from which genotype (columns)? This table will tell us:
##ajust path to file
gamconv<-as.matrix(read.table('gametes.csv', sep='\t', row.names=1, head=T))
gamconv

#What genotype results from which combination of gametes?
##ajust path to file
genconv<-as.matrix(read.table('genotypes.csv', sep='\t', row.names=1, head=T))
genconv

#gradient of carrying capacities
abu<-rep(40,40)

#abu<-rep(c(80, 60, 40, 60, 80), c(2, 2, 32, 2, 2))
length(abu)
abu

#zygotes from gametes - creates n zygotes from a population's gamete vector x
zfg<-function(x,n){             
        #male gemetes are x[1:5]
        mametes<-sapply(runif(n, 0, sum(x[1:5])), function(y) min(which(y<cumsum(x[1:5]))))
        #female gametes are in x[6:8]
        fametes<-sapply(runif(n, 0, sum(x[6:8])), function(z) min(which(z<cumsum(x[6:8]))))
        #find zygote genotype in imported table
        zygotes<-genconv[cbind(mametes, fametes)]
        zygotetab<-vector('numeric',11)
        names(zygotetab)<-c('Ubb', 'Fa', 'Uab', 'Fb', 'Uaa','UUbb', 'FF', 'FUb', 'UUab', 'FUa', 'UUaa')
        for(v in zygotes){zygotetab[which(names(zygotetab)==v)]<-zygotetab[which(names(zygotetab)==v)]+1}
        return(zygotetab)
}


#calculates a population's sex ratio
sexratio<-function(pop){
        #print(pop)
        sum(pop[1:5])/sum(pop)
}

#calculates a population's ratio of present versus theoretically possible fused Xes
Fratio<-function(pop){
        sum(pop*c(0,1,0,1,0, 0,2,1,0,1,0))/sum(pop*c(1,1,1,1,1, 2,2,2,2,2,2))
}

#calculates a population's ratio of present versus theoretically possible fused Xes (in males only)
Fmratio<-function(pop){
        sum(pop*c(0,1,0,1,0, 0,0,0,0,0,0))/sum(pop*c(1,1,1,1,1, 0,0,0,0,0,0))
}

#calculates a population's ratio of present versus theoretically possible fused Xes (in females only)
Ffratio<-function(pop){
        sum(pop*c(0,0,0,0,0, 0,2,1,0,1,0))/sum(pop*c(0,0,0,0,0, 2,2,2,2,2,2))
}


#calculates a population's ratio of present versus theoretically possible Ys (where there is a F no Y can be)
Yratio<-function(pop){
        sum(pop*c(0,1,1,0,2, 0,0,0,1,1,2))/sum(pop*c(2,1,2,1,2, 2,0,1,2,1,2))
}

#calculates a population's ratio of present versus theoretically possible Ys (in males only)
Ymratio<-function(pop){
        sum(pop*c(0,1,1,0,2, 0,0,0,0,0,0))/sum(pop*c(2,1,2,1,2, 0,0,0,0,0,0))
}

#calculates a population's ratio of present versus theoretically possible Ys (in females only) fixation of F causes NAs
Yfratio<-function(pop){
        sum(pop*c(0,0,0,0,0, 0,0,0,1,1,2))/sum(pop*c(0,0,0,0,0, 2,0,1,2,1,2))
}



#plots the ratios defined above for all populations
#does not take arguments, always plots what is in 'genotypes'
#the is a spare panel (bottom right)
genoplot<-function(){
        par(mfrow=c(2,2))

        plot(apply(genotypes$g, 2, Fratio), type='n',ylim=c(0,1), main= 'Chromosomes: F in m - blue, F in f - red,\n Y in m - green, Y in m - orange'); abline(h=.5, lty=2)
        lines(apply(genotypes$g, 2, Fmratio), col='blue', lwd=2)
        lines(apply(genotypes$g, 2, Ffratio), col='red', lwd=2)
        lines(apply(genotypes$g, 2, Ymratio), col='darkgreen', lwd=2)
        lines(apply(genotypes$g, 2, Yfratio), col='orange', lwd=2)
        plot(apply(genotypes$g, 2, sexratio), ylim=c(0,1), main='Sex ratios (male/all)'); abline(h=.5, lty=2)
        
        
        
        plot(abu, ylim=c(0, max(abu)), main='Population sizes')
        
        plot(1:10, 10:1, main=genotypes$c)
}

#matrix of genotypes (40 columns for 11 potential genotypes)
initialise<-function(){
        genotypes<-matrix(0,11,40, dimnames=list(c('Ubb', 'Fa', 'Uab', 'Fb', 'Uaa','UUbb', 'FF', 'FUb', 'UUab', 'FUa', 'UUaa'),NULL))
        
        #populating the matrix of genotypes
        genotypes[c(1,6),1:20]<-.5
        genotypes[c(2,7),21:40]<-.5
        genotypes<-genotypes* matrix(abu, 11,40, byrow=T)
        return(list(g=genotypes,c=0, left=genotypes[,1], right=genotypes[,40]))
}


#this is where the loop/function should start

##tgb (time goes by) takes three arguments:
#the initial genotype matrix, the number of generations, and
#whether there is reinforcemnt from outside the matrix (open=T) or not (open=F, standard)

tgb<-function(genot=genotypes, gens=100, open=F){
        
        for(i in 1:gens){
                
                
                #dispersal 
                if(open==F){
                        genot$g<-genot$g * rep(c((1-disp),(1-disp*2),(1-disp)), c(1,length(abu)-2,1)) + cbind(genot$g,rep(0, 11))[,2:41] * disp + cbind(rep(0, 11), genot$g)[,1:40] * disp
                } else {
                        genot$g<-genot$g * (1-disp*2) + cbind(genot$g,genot$right)[,2:41] * disp + cbind(genot$left , genot$g)[,1:40] * disp 
                }        
                #get matrix of gamete ratios
                gametes<-apply(genot$g,2, function(x) rowSums(matrix(x*selector,8,11,byrow=T)*gamconv))
                
                #get matrix of genotypes (involves drift by samplig from a 'uniform' distribution)
                for(i in 1:length(abu)){
                        genot$g[,i]<-zfg(gametes[,i], abu[i])
                }
                #genot$g<-apply(gametes,2,function(x) zfg(x,80))
                genot$c<-genot$c+1
        }
        #print(paste(counter, 'generations'))
        return(genot)
}

#selector

##on hardy, 100 generations take 3 s, 1000 gens need about 22 s (37 s on curie, and [woohoo] 26 s on weinberg)







#Initialise a genotypes object. That is a list containing $g the genotype array, $c a variable counting the generations,
# $left and $right, the marginal genotypes that will be used in case of open=T.
genotypes<-initialise()


genoplot()


#Run 100 iteration of the time goes by function (each of which simulates 100 generations). Plot after each iteretion.


for(i in 1:100){
        genotypes<-tgb(genotypes, gens=100, open=F)
        genoplot()
}
#genotypes
