#now includes glm and genotypes$trans which tells whether
#the cline centres have transcended the array borders and thus will not be plotted
#to avoid vertical bars

#options(echo=TRUE)
args <- as.numeric(commandArgs(trailingOnly = TRUE))
print (args)
if(length(args)==9){
        
        writeLines('Correct number of arguments supplied.')
        if(NA %in% args){
                writeLines('No all arguments seem to be numeric.')
                stop('Exit.')
        } else {
                sm<-args[1]
                dm<-args[2]
                sf<-args[3]
                df<-args[4]
                sh<-args[5]
                sd<-args[6]
                sc<-args[7]
                disp<-args[8]
                run<-args[9]
                writeLines('The assigned values are:')
                writeLines(sprintf('sm = %1.3f, dm= %1.3f, sf = %1.3f, df = %1.3f', args[1], args[2], args[3], args[4]))
                writeLines(sprintf('sh = %1.3f, sd= %1.3f, sc = %1.3f', args[5], args[6], args[7]))
                writeLines(sprintf('disp = %1.3f , run = %d', args[8], args[9]))
        }
} else {
        writeLines('Too few arguments supplied.')
        stop('Exit.')
}




##Adjust working directory to the location of genotypes.csv and gametes.csv
#setwd('/home/hannes/Dropbox/scripts/FvU3D/')




#defines the dimensions of the genotype array

initialise<-function(wi=10, len=40, popsize=40){
        genotypenr<-11
        genotypenames<-c('Ubb', 'Fa', 'Uab', 'Fb', 'Uaa', 'UUbb', 'FF', 'FUb', 'UUab', 'FUa', 'UUaa')
        a<-array(0, dim=c(wi, len, genotypenr),
                dimnames=list(paste0('w', 1:wi), paste0('l', 1: len), genotypenames))
        
        abu<-matrix(popsize, wi, len)
        
        a[,1:20, c(1,6)]<-0.5
        a[,21:40, c(2,7)]<-0.5
        for(i in 1:dim(a)[3]){
                a[,,i]<-a[,,i]*abu
        }
        
        trans=c(0,0)
        names(trans)<=c('Fout', 'Yout')
        
        return(list(g=a, c=0, left=a[,1,], right=a[,len,], to=a[1,,], bo=a[wi,,], abu=abu, trans=trans))
}

#genotypes<-initialise()
#genotypes$g
#genotypes$c
#genotypes$left
#genotypes$right
#genotypes$to
#genotypes$bo
#genotypes$abu

#Selection parameters (cf. Veltsos et al. 2008 Tab. 1)

#sm<-0.05         #selection against b in males
#dm<-1           #dominance of a in males
#sf<-0.02           #selection against a in females
#df<-1           #dominance of a in females
#sh<-0           #"heterozygote disadvantage"
#sd<-0           #"dosage compensation"
#sc<-0           #"sex chromosome coadaption"

#disp<-0.03      #dispersal in each direction

#vector of selection acting on genotypes (SELECTion vectOR)
selector<-c(1-sm, 1, (1-sm*dm)*(1-.5*sd)*(1-.5*sc), (1-sm)*(1-.5*sd)*(1-sc), (1-.5*sd)*(1-sc),
            1, 1, (1-sh)*(1-.5*sd)*(1-.5*sc), (1-sf*df)*(1-.5*sd)*(1-.5*sc), (1-sf*df)*(1-sh)*(1-.5*sc), (1-sf)*(1-sd)*(1-sc))
names(selector)<-c('Ubb', 'Fa', 'Uab', 'Fb', 'Uaa', 'UUbb', 'FF', 'FUb', 'UUab', 'FUa', 'UUaa')

writeLines("\nGenotypes and selection acting on them:")
print(selector)


##What gamete ratios (rows) result from which genotype (columns)? This table will tell us:
##ajust path to file
gamconv<-as.matrix(read.table('gametes.csv', sep='\t', row.names=1, head=T))

writeLines('\n Gamete ratios:\n')
gamconv

#What genotype results from which combination of gametes?
##ajust path to file
genconv<-as.matrix(read.table('genotypes.csv', sep='\t', row.names=1, head=T))

writeLines('\n Genotypes from gametes:\n')
genconv


#'zygotes from gametes' creates n zygotes from a population's gamete vector x
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

##########
##########Find a way to visualise the results!

Fratio<-function(genot=genotypes$g){
        apply(genot, c(1,2), function(x) sum(x*c(0,1,0,1,0, 0,2,1,0,1,0))/sum((x * c(1,1,1,1,1, 2,2,2,2,2,2))))
}

YFratio<-function(genot=genotypes$g){
        apply(genot, c(1,2), function(x) sum(x*c(0,1,1,0,2, 0,0,0,1,1,2))/sum((x * c(2,1,2,1,2, 2,0,1,2,1,2))))
}


#define some utility functions 

Fpot<-function(genot){
  apply(genot, c(1,2), function(x) sum(x* c(1,1,1,1,1, 2,2,2,2,2,2)))
}

Fhits<-function(genot){
  apply(genot, c(1,2), function(x) sum(x*c(0,1,0,1,0, 0,2,1,0,1,0)))
}



Fmisses<-function(genot){
  apply(genot, c(1,2), function(x) sum(x*c(1,0,1,0,1, 2,0,1,2,1,2)))
}


#get do the same for the Y
Ypot<-function(genot){
  apply(genot, c(1,2), function(x) sum(x* c(2,1,2,1,2, 2,0,1,2,1,2)))
}

Yhits<-function(genot){
  apply(genot, c(1,2), function(x) sum(x*c(0,1,1,0,2, 0,0,0,1,1,2)))
}



Ymisses<-function(genot){
  apply(genot, c(1,2), function(x) sum(x*c(2,0,1,1,0, 2,0,1,1,0,0)))
}



#####Print results to png
print_to_png<-function(c, data){
        
        if(c %% 100 ==0){
                
                png(paste('sm', sm, 'dm', dm, 'sf', sf, 'df', df, 'sh', sh, 'sd', sd, 'sc', sc, 'disp', disp,
                          sprintf('gen_%05d.png', c), sep="_"), width=1400, height=800)
                
                par(mfrow=c(1,2))
                print(dim(data))
                image(Fratio(data), asp=dim(data)[2]/dim(data)[1], frame=F, col=heat.colors(100), breaks=seq(0,1,0.01), axes=F, main=paste('F,', c, 'generations'))
                #image(abc$g[,,6], col=heat.colors(40), breaks=seq(0,40,1), main=paste(abc$c, 'generations'), asp=4, frame=F)
                contour(Fratio(data), add=T, levels=c(0, .05, .5, .95,1))
                abline(h=0.5, lty=2)
                image(YFratio(data), asp=dim(data)[2]/dim(data)[1], frame=F, col=heat.colors(100), breaks=seq(0,1,0.01), axes=F, main=paste('Y|not F,', c, 'generations'))
                contour(YFratio(data), add=T, levels=c(0, .05, .5, .95,1))
                abline(h=0.5, lty=2)
                
                dev.off()
                               
        }
}


#########
#########The new TIME GOES BY function
tgb<-function(genot=genotypes, gens=10, open=F, printout=F){
        
        #defining and internal object representing the genotypes, 'intgen'
        
        dimensions<-dim(genot$g)

        intgen<-array(0, c(dimensions[1]+2,dimensions[2]+2,dimensions[3]))
        
        intgen[2:(dimensions[1]+1), 2:(dimensions[2]+1),]<-
                intgen[2:(dimensions[1]+1), 2:(dimensions[2]+1),]+genot$g
        
        #if it is closed, create a matrix, closedfact
        if(open==F){
                
                f4<-1-4*disp
                f3<-1-3*disp
                f2<-1-2*disp
                
                closedfact<-array(0, dimensions)
                closedfact[2:(dimensions[1]-1), 2:(dimensions[2]-1),]<-f4
                closedfact[1:(dimensions[1]), c(1, dimensions[2]),]<-f3
                closedfact[c(1, dimensions[1]), 1:(dimensions[2]),]<-f3
                closedfact[1,1,1:dimensions[3]]<-f2
                closedfact[dimensions[1],1,1:dimensions[3]]<-f2
                closedfact[1,dimensions[2],1:dimensions[3]]<-f2
                closedfact[dimensions[1],dimensions[2],1:dimensions[3]]<-f2
                #print(closedfact)
                
                #if it is not closed, intgen's margins will be set equal to their neighbouring elements        
        } else {               
                
                intgen[2:(dimensions[1]+1),   1,]<-genot$left
                intgen[2:(dimensions[1]+1), dimensions[2]+2,]<-genot$right
                intgen[1,  2:(dimensions[2]+1),]<-genot$to
                intgen[dimensions[1]+2, 2:(dimensions[2]+1),]<-genot$bo
                
        }
     
        
        #the actual loop starts here
        if(open==F){
                for(i in 1:gens){
                        
                        #dispersal
                        intgen[2:(dimensions[1]+1), 2:(dimensions[2]+1),]<-
                                intgen[2:(dimensions[1]+1), 2:(dimensions[2]+1),]*closedfact+
                                intgen[1:(dimensions[1]),   2:(dimensions[2]+1),]*disp+
                                intgen[3:(dimensions[1]+2), 2:(dimensions[2]+1),]*disp+
                                intgen[2:(dimensions[1]+1), 1:(dimensions[2]),]*disp+
                                intgen[2:(dimensions[1]+1), 3:(dimensions[2]+2),]*disp
                        
                        #gametes
                        #get array of gamete ratios
                        gametes<-apply(intgen[2:(dimensions[1]+1), 2:(dimensions[2]+1),], c(1,2),
                                       function(x) rowSums(matrix(x*selector,8,11,byrow=T)*gamconv))                        
                        
                        #get matrix of genotypes (involves drift by samplig from a 'uniform' distribution)
                        
                        for(k in 1:dimensions[1]){
                                #print(k)
                                for(l in 1:dimensions[2]){
                                        #print(l)
                                        intgen[k+1, l+1,]<-zfg(gametes[,k,l], genot$abu[k,l])        
                                }
                        }
                        
                        
                        
                        
                        #counter
                        genot$c<-genot$c+1
                        print(paste(genot$c, 'generations'))
                        
                        if(printout==T){print_to_png(genot$c, intgen[2:(dimensions[1]+1),2:(dimensions[2]+1),])}
                        
                        
                }
        } else {        #i.e., if not open, do this loop:
                for(i in 1:gens){
                        
                        #dispersal
                        intgen[2:(dimensions[1]+1), 2:(dimensions[2]+1),]<-
                                intgen[2:(dimensions[1]+1), 2:(dimensions[2]+1),]*(1-4*disp)+
                                intgen[1:(dimensions[1]), 2:(dimensions[2]+1),]*disp+
                                intgen[3:(dimensions[1]+2), 2:(dimensions[2]+1),]*disp+
                                intgen[2:(dimensions[1]+1), 1:(dimensions[2]),]*disp+
                                intgen[2:(dimensions[1]+1), 3:(dimensions[2]+2),]*disp
                        
                        
                        
                        #gametes
                        #get array of gamete ratios
                        gametes<-apply(intgen[2:(dimensions[1]+1), 2:(dimensions[2]+1),], c(1,2),
                                       function(x) rowSums(matrix(x*selector,8,11,byrow=T)*gamconv))                        
                        
                        #get matrix of genotypes (involves drift by samplig from a 'uniform' distribution)
                        
                        for(k in 1:dimensions[1]){
                                #print(k)
                                for(l in 1:dimensions[2]){
                                        #print(l)
                                        intgen[k+1, l+1,]<-zfg(gametes[,k,l], genot$abu[k,l])        
                                }
                        }
                        
                        
                        
                        
                        #counter
                        genot$c<-genot$c+1
                        print(paste(genot$c, 'generations'))
                        if(printout==T){print_to_png(genot$c, intgen[2:(dimensions[1]+1),2:(dimensions[2]+1),])}
                        
                }
        }
        
        
        
        return(list(g=intgen[2:(dimensions[1]+1), 2:(dimensions[2]+1),], c=genot$c, left=genot$left, right=genot$right,
                    to=genot$to, bo=genot$bo, abu=genot$abu, trans=genot$trans))
        
}







genotypes<-initialise(40,40,40)





#running tgb in a loop plotting graphs every 10 generatins
for (i in 1:1000){
        genotypes<-tgb(genotypes,gens=10,open=F)
        
        
        #glms to estimate cline positions
        
        genodims<-dim(genotypes$g)
        #print(genodims)
        Fh<-Fhits(genotypes$g)
        Fm<-Fmisses(genotypes$g)
        #Fp<-Fpot()
        
        Yh<-Yhits(genotypes$g)
        Ym<-Ymisses(genotypes$g)
        #Yp<-Ypot()
        
        glma <- glm(cbind(as.vector(Fh), as.vector(Fm))~rep(1:genodims[2], each=genodims[1]), family=binomial)
        glmb <- glm(cbind(as.vector(Yh), as.vector(Ym))~rep(1:genodims[2], each=genodims[1]), family=binomial)
        
        
        IntF<-coefficients(glma)[1]
        SlF<-coefficients(glma)[2]
        #centre of cline seems to be minus intercept over slope
        CentF<- -IntF/SlF
        
        IntY<-coefficients(glmb)[1]
        SlY<-coefficients(glmb)[2]
        CentY<- -IntY/SlY
        
        if(genotypes$trans[1]==0){
                if (CentF<0 | CentF>40){
                        genotypes$trans[1]<-1
                }
        }
        
        if(genotypes$trans[2]==0){
                if (CentY<0 | CentY>40){
                        genotypes$trans[2]<-1
                }
        }
        
        
        outvec<-c(IntF, SlF, CentF, IntY, SlY, CentY, genotypes$c)
        names(outvec)<-c('IntF', 'SlF', 'CentF', 'IntY', 'SlY', 'CentY', 'gen')
        
        #outvec, the vector containing the glm coefficients and cline centres, will be added to genotypes
        #just before the save() step, it will be lost during the time goes by function
        genotypes$coefs <- outvec
          
          
        save(genotypes, file=sprintf('sm%1.3f_dm%1.3f_sf%1.3f_df%1.3f_sh%1.3f_sd%1.3f_sc%1.3f_disp%1.3f_run%d_gens%05d.RData',
                                     sm, dm, sf, df, sh, sd, sc, disp, run, genotypes$c))
        
        #check if there is still change, break if not  (Check if square diffs to mean are zero. If so, break.)
        Fr<-Fratio(genotypes$g)
        Yr<-YFratio(genotypes$g)
        
        
        if((sum((Fr-mean(Fr))**2)==0) & (sum((Yr-mean(Yr))**2)==0)){
                writeLines(sprintf('No change in F and Y ratios in generation %d.', genotypes$c))
                break
        }
}



#tbg can be run in printout mode saving heat map and contour every 100 generations
#files will be saved to the working directory and will be names with the parameters chosen

#tgb(genotypes, gens=20000, printout=T)



######################
######################
######################
##TEST AREA


#100 gens need approx. 22 s on hardy
