#should now avoid plotting vertical bars (which were the result of crazy glms [after transcendence of the cline centres])
#Since the glms will now be calculated by the main script, FvU3D_apo_glm.R,
#this is only to extract $trans and $coefs from genotypes and to do both, draw the dance and save a table.


#what have we got in this directory (that is relevant)?

whwg<- function (dir=getwd()) {
        liste<-list.files(path = dir, pattern="^sm.*")
        tab <- table(sapply(liste, function(x) substring(x, 1, 70), USE.NAMES=F))
        extr<-t(sapply(names(tab), function(b) as.numeric(c(substr(b,3, 7), substr(b,11, 15), substr(b,19, 23), substr(b,27, 31),
                                                            substr(b,35, 39), substr(b, 43, 47), substr(b,51, 55), substr(b,61, 65),
                                                            substr(b,70, 70))), USE.NAMES=F))
        extr<-cbind(extr, tab)
        dimnames(extr)<-list(NULL, c('sm', 'dm', 'sf', 'df', 'sh', 'sd', 'sc', 'disp', 'run', 'gen'))
        
        return(extr)
        
        
}


#uses 'what have we got' function from list_output.R to get a table of parameters
all<-whwg()

#remove single files, get_coefs() needs at least two files

todo<-all[all[,10]>1,]

#function to get the coefficients and trancendence data

getcoeftrans <- function(file){
        load(file)
        return(c(genotypes$coefs, genotypes$trans))
}




#Funtion to create a file name from selection parameters
#the input needs to be a numeric vector of length 10
nutona <- function (b) {
        sprintf('sm%1.3f_dm%1.3f_sf%1.3f_df%1.3f_sh%1.3f_sd%1.3f_sc%1.3f_disp%1.3f_run%d_gens%05d.RData', b[1], b[2],
                b[3], b[4], b[5], b[6], b[7], b[8], as.integer(b[9]), as.integer(b[10])*10)
}

#in case there is only one set of parameters
if(is.null(dim(todo)[1])){
       
        
        t0<-Sys.time()
        
        abc<-sapply(list.files(path = getwd(), pattern=paste0(substring(nutona(todo), 1, 74), '.*')),
                    function(x) getcoeftrans(x), USE.NAMES=F)
        
        abc<- abc[,order(abc[7,])]
        
        if(1 %in% abc[8,]){
                Flim<-min(which(abc[8,]==1))
        } else {Flim <- 1000}
        if(1 %in% abc[9,]){
                Ylim<-min(which(abc[9,]==1))
        } else {Ylim <- 1000}
        if(Flim<Ylim){abslim<-Ylim} else {abslim<-Flim}
        
        png(paste0('dances/dance_', substr(nutona(todo), 1, 79), '.png'), width= 1400, height= 800)
        plot(abc[7,1:Flim], abc[3,1:Flim], type='l', col='blue', ylim=c(0,40), xlim=c(0,10000), main=substr(nutona(todo), 1, 79), frame=F)
        points(abc[7,1:Flim], abc[3,1:Flim]- abs(1/abc[2,0:Flim]/2), type='l', col='lightblue')
        points(abc[7,1:Flim], abc[3,1:Flim]+ abs(1/abc[2,0:Flim]/2), type='l', col='lightblue')
        points(abc[7,1:Ylim], abc[6,1:Ylim], type='l', col='red')
        points(abc[7,1:Ylim], abc[6,1:Ylim]+ abs(1/abc[5,0:Ylim]/2), type='l', col='pink')
        points(abc[7,1:Ylim], abc[6,1:Ylim]- abs(1/abc[5,0:Ylim]/2), type='l', col='pink')
        segments(0,0,10000,0)
        segments(0,40,10000,40)
        dev.off()
        
        write.table(t(abc), paste0('dances/dance_', substr(nutona(todo), 1, 79)),
                    sep='\t', row.names=F)
        
        print(Sys.time()-t0)
} else {
        for(i in 1:dim(todo)[1]){
                
                print(paste0('Parameter set ', i, ' of ', nrow(todo), ':'))
                print(todo[i,])
                
                t0<-Sys.time()
                
                abc<-sapply(list.files(path = getwd(), pattern=paste0(substring(nutona(todo[i,]), 1, 74), '.*')),
                            function(x) getcoeftrans(x), USE.NAMES=F)
                
                abc<- abc[,order(abc[7,])]
                
                if(1 %in% abc[8,]){
                        Flim<-min(which(abc[8,]==1))
                } else {Flim <- 1000}
                if(1 %in% abc[9,]){
                        Ylim<-min(which(abc[9,]==1))
                } else {Ylim <- 1000}
                if(Flim<Ylim){abslim<-Ylim} else {abslim<-Flim}
                
                png(paste0('dances/dance_', substr(nutona(todo[i,]), 1, 79), '.png'), width= 1400, height= 800)
                plot(abc[7,1:Flim], abc[3,1:Flim], type='l', col='blue', ylim=c(0,40), xlim=c(0,10000), main=substr(nutona(todo[i,]), 1, 79), frame=F)
                points(abc[7,1:Flim], abc[3,1:Flim]- abs(1/abc[2,0:Flim]/2), type='l', col='lightblue')
                points(abc[7,1:Flim], abc[3,1:Flim]+ abs(1/abc[2,0:Flim]/2), type='l', col='lightblue')
                points(abc[7,1:Ylim], abc[6,1:Ylim], type='l', col='red')
                points(abc[7,1:Ylim], abc[6,1:Ylim]+ abs(1/abc[5,0:Ylim]/2), type='l', col='pink')
                points(abc[7,1:Ylim], abc[6,1:Ylim]- abs(1/abc[5,0:Ylim]/2), type='l', col='pink')
                segments(0,0,10000,0)
                segments(0,40,10000,40)
                dev.off()
                
                write.table(t(abc), paste0('dances/dance_', substr(nutona(todo[i,]), 1, 79)),
                            sep='\t', row.names=F)
                
                print(Sys.time()-t0)
        }
        
}