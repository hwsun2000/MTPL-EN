
library(glmnet)
library(survival)
library(mvtnorm)
library(glmnetUtils)
library(SIS)
library(survival)
library(glmpath)
library(MASS)
library(pROC)


dyn.load('coxmart.dll')

n=300  #sample size
rho=0.9 #correlation coefficients
sample=100 #simulation times
q <- 12
p <- 1000# covarites
CL<- 1

CU<- 3

epsilon=0.1

block1 <- c(rep(0.3,5),rep(0,45))
block2 <- c(rep(-0.4,5),rep(0,45))
block3 <- c(rep(0.5,5),rep(0,45))
block4 <- c(rep(-0.6,5),rep(0,45))
block5 <- c(rep(0.7,5),rep(0,45))
block6 <- c(rep(-0.8,5),rep(0,45))
block7 <- c(rep(0.8,5),rep(0,45))
block8 <- c(rep(-0.7,5),rep(0,45))
block9 <- c(rep(0.6,5),rep(0,45))
block10 <- c(rep(-0.5,5),rep(0,45))
block11 <- c(rep(0.4,5),rep(0,45))
block12 <- c(rep(-0.3,5),rep(0,45))

block13 <- c(rep(0,p-600))
beta <- matrix(c(block1,block2,block3,block4,block5,block6,block7,block8,block9,block10,block11,block12,block13))
beta1 <- which(beta!=0)##index variable
beta0 <-which(beta==0)

#generate correlation matrix
size=50 #colum number of sigma_rho
num=20  #diagnal blocks in sigma
sigma <- matrix(rep(0, times =size*size*num*num),nrow=size*num)

for (rr in 1:num){
    for (ss in 1:size){
        for (ks in 1:size){
            if (rr%%2) sigma[(rr-1)*size+ss,(rr-1)*size+ks]=rho^abs(ss-ks)
            else sigma[(rr-1)*size+ss,(rr-1)*size+ks]=(-rho)^abs(ss-ks)
        }
    }
}


nn=30   #replicates 

  
    lfdr=rep(0,nn)     #EN
    ltfdr=rep(0,nn)   #rwt
    lt2fdr=rep(0,nn)   #raw
    
     lpsr=rep(0,nn)    
     ltpsr=rep(0,nn)  
     lt2psr=rep(0,nn)  
    
     lnum=rep(0,nn)     
     ltnum=rep(0,nn) #rwt
     lt2num=rep(0,nn)   #raw

     lh=rep(0,nn)     ###likelihood
     ltlh=rep(0,nn) 
     lt2lh=rep(0,nn)

     PI=rep(0,nn)   ##prognostic index
     tPI=rep(0,nn)   ##prognostic index 
     t2PI=rep(0,nn)   ##prognostic index 


     onum=rep(0,nn)     ##EN 
     otnum=rep(0,nn)     ### rwt
     ot2num=rep(0,nn)    ###raq
  

     power0=rep(0,nn)  #EN,proportion of correctly identified outliers
     Ierror0=rep(0,nn)  #EN,proportion of false discoveries

     power=rep(0,nn)  #rwt,proportion of correctly identified outliers
     Ierror=rep(0,nn)  #rwt,proportion of false discoveries

     otnum_2=rep(0,nn)   
      power_2=rep(0,nn)  
     Ierror_2=rep(0,nn)  

      power2=rep(0,nn)  #raw,roportion of correctly identified outliers
     Ierror2=rep(0,nn)  #raw,proportion of false discoveries

     devresmat0<-matrix(0,nn,n)
      AUC0<-rep(0,nn)
     

      devresmat<-matrix(0,nn,n)
      AUCMTL<-rep(0,nn)

      outliers_indi<-matrix(0,nn,n) 
      rwtindicesjl<-matrix(0,nn,n)
      rawindicesjl<-matrix(0,nn,n)

      ht_mat<-matrix(0,nn,n)
    
      betajl<-matrix(0,nn,p)
      tbetajl<-matrix(0,nn,p)
      t2betajl<-matrix(0,nn,p)
 
      CensorRate=rep(0,nn)  #记录截尾比例

   for(i in 1:nn){
   seed　<- 8601+i
   set.seed(seed)

   x=rmvnorm(n,rep(0,p),sigma)
    
       ht=exp(x%*%beta)
       maxht=max(ht)
       minht=min(ht)
       
        kk=sample(n,floor(n*epsilon)) 
     
     seed　<- 8007+i
     set.seed(seed)
  
     flag=rbinom(length(kk),1,0.5)
     ht[kk]=minht*flag+(maxht)*(1-flag)
        
       routliers=kk
       rinliers=(1:n)[-routliers]
  
       outliers_indi[i,routliers]=1

	 #observed survivial time
	
	 logS=log(matrix(runif(n,0,1),n,1)) #log[S(t)]
	 T=-logS/ht  #survival time

	 #censored time
	 myrate <- matrix(runif(n,CL,CU),n,1)/ht
	 C <- rexp(n, rate = 1/myrate)
	
	#survival time and state
	 time <- apply(cbind(C,T),1,min)
 
      t<-time
	stat <- (T<C)
      surv=Surv(time,stat)

      CensorRate[i]<- length(which(as.numeric(stat) == 0))/n*100
  
     ##generate test data

         seed　<- 8103+i
        set.seed(seed) 
      x1=rmvnorm(n,rep(0,p),sigma)  
       ht1=exp(x1%*%beta)
      
       #observed survivial time
	
      logS=log(matrix(runif(n,0,1),n,1)) #log[S(t)]
      T1=-logS/ht1  #survival time

	 #censored time
	 myrate1 <- matrix(runif(n,CL,CU),n,1)/ht1
	 C1<- rexp(n, rate = 1/myrate1)
	
	#survival time and state
	 time1 <- apply(cbind(C1,T1),1,min)

	stat1<-(T1<C1)
      surv1=Surv(time1,stat1)     ###test data
   
       center=apply(x,2,mean)    
      sd=apply(x,2,sd)                
       x.c = sweep(x, 2, center)
       xsd=sweep(x.c,2,sd,"/") 


       center=apply(x1,2,mean)    
      sd=apply(x1,2,sd)            
       x.c = sweep(x1, 2, center)
       xsd1=sweep(x.c,2,sd,"/") 

     #####glmnet
     
      my.alpha <- seq(0.1,1,0.1)     
      cvm=rep(0,length(my.alpha))
      var.selected.EN <- matrix(0,p,length(my.alpha)) 
      link<- matrix(0,n,length(my.alpha))     
      
     fit.EN.cv=cva.glmnet(xsd,y=surv,alpha=my.alpha,family="cox")

      for (j in 1:length(my.alpha)){

           cvm[j]=fit.EN.cv$modlist[[j]]$cvm[which(fit.EN.cv$modlist[[j]]$cvm==min(fit.EN.cv$modlist[[j]]$cvm))]

          var.selected.EN[,j]=as.matrix(coef(fit.EN.cv$modlist[[j]],s="lambda.min"))
          link[,j]=xsd%*%var.selected.EN[,j]
 
       }


        l.cvl=var.selected.EN[,which(cvm==min(cvm))]
   
        ww=as.matrix(l.cvl)

         betajl[i,]<-ww

        lnum[i]=length(which(ww!=0))   

        shi1=which(beta!=0)    
        TP=length(which(ww[shi1]!=0))    
        FN=length(which(ww[shi1]==0))   

        shi2=which(beta==0)  
        TN=length(which(ww[shi2]==0)) 
        FP=length(which(ww[shi2]!=0))    

        lpsr[i]=TP/(TP+FN)
      
         if((TP+FP)==0)
           {lfdr[i]=0} else
          {lfdr[i]=FP/(TP+FP)}
 
         lh[i]=logplik(xsd1,time1,stat1,ww)  

           PI[i]=mean((xsd1%*%ww-xsd1%*%beta)^2)   ####prognostic index
        

       link_EN=link[,which(cvm==min(cvm))]
   
       strata=rep(1,n)    

       time <- time          
       status <- stat   
  
	 sorted <- order(strata, time)      
	 strata <- strata[sorted]        
	 newstrat <- as.integer(c(1*(diff(as.numeric(strata))!=0), 1))  
                                                                     
       stime <- as.double(time[sorted])
       sstat <- as.integer(status[sorted])
     
     ### method=efron
   
       score=exp(link_EN)[sorted]
       weights=rep(1,n)  
      
      coxres <- .C("coxmart", as.integer(n),
				1,
				stime,
				sstat,
				newstrat,
				as.double(score),
				as.double(weights),
				resid=double(n))

            resid <- double(n)
            resid[sorted] <- coxres$resid     ##Matingaile residual 得到

    	     devres=sign(resid) *sqrt(-2* (resid+
			     ifelse(status==0, 0, status*log(status-resid))))

 
         devresmat0[i,]<-devres
  
        

       AUC0[i]<-auc(roc(outliers_indi[i,],abs(devres)))
   

        q1 <- qnorm(0.975)
        q2 <- qnorm(0.975)
 
       ok1<- which(devres<=q1&devres>=0)

      ok2<- which(devres>=-q2&devres<0)


       ok=c(ok1,ok2)


        rwt.indices<-sort(ok)
        outliers_EN=(1:n)[-ok]

      odent_EN=outliers_EN
      onum[i]=length(odent_EN)
        ident_EN=(1:n)[-odent_EN]

      TPresult<- length(intersect(odent_EN,routliers)) #TP
	FPresult<- length(intersect(odent_EN,rinliers)) #FP
      TNresult <- length(intersect(ident_EN,rinliers)) #TN
	FNresult <- length(intersect(ident_EN,routliers)) #FN
    
      power0[i]=TPresult/length(routliers)
      Ierror0[i]=FPresult/length(rinliers)



    ####mtl lasso
    
        mtlcox<-fitter(time,stat,x,gamma=0.05,iter="EN",method="EN",n.multistart=1)    
        l3=mtlcox$rwt.betas
          tbetajl[i,]<-l3
   
       ltnum[i]=length(which(l3!=0)) 
     
       TP=length(which(l3[shi1]!=0)) 
       FN=length(which(l3[shi1]==0)) 
       TN=length(which(l3[shi2]==0))    
       FP=length(which(l3[shi2]!=0))   
       ltpsr[i]=TP/(TP+FN)
      
        if((TP+FP)==0)
        {ltfdr[i]=0} else
      {ltfdr[i]=FP/(TP+FP)}
         

       ltlh[i]=logplik(xsd1,time1,stat1,l3)  
    
   tPI[i]=mean((xsd1%*%l3-xsd1%*%beta)^2)   ####prognostic index


         ###C-index
  
      ## link=xsd1%*%l3   

    ##  AUC.CC0=risksetAUC(Stime=time,status=stat,marker=link,method="Cox", tmax=30)

     ### tcdx[i]=AUC.CC0$Cindex
    
     
       odent=mtlcox$rwt.outliers
      otnum[i]=length(odent)
        ident=(1:n)[-odent]

         rwtindicesjl[i,odent]<-1
     


      TPresult<- length(intersect(odent,routliers)) #TP
	FPresult<- length(intersect(odent,rinliers)) #FP
      TNresult <- length(intersect(ident,rinliers)) #TN
	FNresult <- length(intersect(ident,routliers)) #FN
    
      power[i]=TPresult/length(routliers)
      Ierror[i]=FPresult/length(rinliers)
  

  odent=mtlcox$rwt.outliers2       
      otnum_2[i]=length(odent)    
        ident=(1:n)[-odent]

  
      TPresult<- length(intersect(odent,routliers)) #TP
	FPresult<- length(intersect(odent,rinliers)) #FP
      TNresult <- length(intersect(ident,rinliers)) #TN
	FNresult <- length(intersect(ident,routliers)) #FN
    
      power_2[i]=TPresult/length(routliers)     
      Ierror_2[i]=FPresult/length(rinliers)    


    ####raw

       l2=mtlcox$max.beta
        t2betajl[i,]<-l2

       lt2num[i]=length(which(l2!=0))   
     
       TP=length(which(l2[shi1]!=0))   
       FN=length(which(l2[shi1]==0))   
       TN=length(which(l2[shi2]==0))    
       FP=length(which(l2[shi2]!=0))   
       lt2psr[i]=TP/(TP+FN)
      
        if((TP+FP)==0)
        {lt2fdr[i]=0} else
      {lt2fdr[i]=FP/(TP+FP)}
         

       lt2lh[i]=logplik(xsd1,time1,stat1,l2) 

    t2PI[i]=mean((xsd1%*%l2-xsd1%*%beta)^2)   ####prognostic index

  


     
       ident2=mtlcox$units
      odent2=(1:n)[-ident2]
      ot2num[i]=length(odent2)
       
      rawindicesjl[i,odent2]<-1

      TPresult<- length(intersect(odent2,routliers)) #TP
	FPresult<- length(intersect(odent2,rinliers)) #FP
      TNresult <- length(intersect(ident2,rinliers)) #TN
	FNresult <- length(intersect(ident2,routliers)) #FN
    
      power2[i]=TPresult/length(routliers)
      Ierror2[i]=FPresult/length(rinliers)

      devresmat[i,]<-mtlcox$devres
     
            ht_mat[i,]<-ht

      AUCMTL[i]<-auc(roc(outliers_indi[i,],abs(mtlcox$devres)))

   }

   ##EN

   lnum
   lpsr
   lfdr
   sqrt(lpsr*(1-lfdr))
   lh
   
   onum
    power0
   Ierror0



  ###均数
      Zlnum=mean(lnum) 
      Zlpsr=mean(lpsr)
      Zlfdr=mean(lfdr)
      ZlGM=mean(sqrt(lpsr*(1-lfdr)))    
      Zllh=mean(lh)
      ZPI=mean(PI)

      Zonum=mean(onum)
      Zpower0=mean(power0)
      ZIerror0=mean(Ierror0)
         ZAUC0=mean(AUC0)
  ###中位数
      Mlnum=median(lnum) 
      Mlpsr=median(lpsr)
      Mlfdr=median(lfdr)
      MlGM=median(sqrt(lpsr*(1-lfdr)))    
      Mllh=median(lh)
     MPI=median(PI)

      Monum=median(onum)
      Mpower0=median(power0)
      MIerror0=median(Ierror0)
      MAUC0=median(AUC0)

   ##  mad/0.6749

        Malnum=mad(lnum) /0.6749
      Malpsr=mad(lpsr) /0.6749
      Malfdr=mad(lfdr) /0.6749
      MalGM=mad(sqrt(lpsr*(1-lfdr))) /0.6749    
      Mallh=mad(lh) /0.6749
     MaPI=mad(PI)/0.6749
      Maonum=mad(onum)/0.6749    
      Mapower0=mad(power0)/0.6749    
      MaIerror0=mad(Ierror0)/0.6749    
       MaAUC0=mad(AUC0)/0.6749




 #Rwt MTL EN


    ltnum
    ltpsr
    ltfdr
    sqrt(ltpsr*(1-ltfdr))
    ltlh

    otnum
    power
    Ierror

  

     Zltnum=mean(ltnum) 
     Zltpsr=mean(ltpsr)
     Zltfdr=mean(ltfdr)
     ZltGM=mean(sqrt(ltpsr*(1-ltfdr)))   
     Zltlh=mean(ltlh)
     ZtPI=mean(tPI)
  
    Zotnum=mean(otnum)      
     Zpower=mean(power)
     ZIerror=mean(Ierror)
      ZAUCMTL=mean(AUCMTL)

   
    Zotnum_2=mean(otnum_2)

   Zpower_2=mean(power_2)
      ZIerror_2=mean(Ierror_2)



 ###median
      Mltnum=median(ltnum) 
      Mltpsr=median(ltpsr)
      Mltfdr=median(ltfdr)
      MltGM=median(sqrt(ltpsr*(1-ltfdr)))    
      Mltlh=median(ltlh)
      MtPI=median(tPI)
  
      Motnum=median(otnum)
      Mpower=median(power)
      MIerror=median(Ierror)
          MAUCMTL=median(AUCMTL)

   Motnum_2=median(otnum_2)
  
         Mpower_2=median(power_2)
      MIerror_2=median(Ierror_2)




   ##  mad/0.6749

       Maltnum=mad(ltnum) /0.6749
      Maltpsr=mad(ltpsr) /0.6749
      Maltfdr=mad(ltfdr) /0.6749
      MaltGM=mad(sqrt(ltpsr*(1-ltfdr))) /0.6749    
      Maltlh=mad(ltlh) /0.6749
      MatPI=mad(tPI)/0.6749

      Maotnum=mad(otnum)/0.6749    
      Mapower=mad(power)/0.6749    
      MaIerror=mad(Ierror)/0.6749  
      MaAUCMTL=mad(AUCMTL)/0.6749
  
    Maotnum_2=mad(otnum_2)/0.6749  
  
         Mapower_2=mad(power_2)/0.6749  
      MaIerror_2=mad(Ierror_2)/0.6749  

   
 #Raw MTL EN
   lt2num
   lt2psr
   lt2fdr
   lt2fdr
   sqrt(lt2psr*(1-lt2fdr))
   lt2lh

   ot2num
   power2
   Ierror2


  CensorRate



 ####average

     Zlt2num=mean(lt2num) 
     Zlt2psr=mean(lt2psr)
     Zlt2fdr=mean(lt2fdr)
     Zlt2GM=mean(sqrt(lt2psr*(1-lt2fdr)))  
     Zlt2lh=mean(lt2lh)  
       Zt2PI=mean(t2PI)  
  
     Zot2num=mean(ot2num)   
     Zpower2=mean(power2)
     ZIerror2=mean(Ierror2)

   ###median
      Mlt2num=median(lt2num) 
      Mlt2psr=median(lt2psr)
      Mlt2fdr=median(lt2fdr)
      Mlt2GM=median(sqrt(lt2psr*(1-lt2fdr)))    
      Mlt2lh=median(lt2lh)
      Mt2PI=median(t2PI) 
  
      Mot2num=median(ot2num)
      Mpower2=median(power2)
      MIerror2=median(Ierror2)

   ##  mad/0.6749

       Malt2num=mad(lt2num) /0.6749
      Malt2psr=mad(lt2psr) /0.6749
      Malt2fdr=mad(lt2fdr) /0.6749
      Malt2GM=mad(sqrt(lt2psr*(1-lt2fdr))) /0.6749    
      Malt2lh=mad(lt2lh) /0.6749
         Mat2PI=mad(t2PI)/0.6749 
      Maot2num=mad(ot2num)/0.6749    
      Mapower2=mad(power2)/0.6749    
      MaIerror2=mad(Ierror2)/0.6749  

ZCensorRate=mean(CensorRate)
   MCensorRate=median(CensorRate)
  ### 平均值
    ENoutput=c(Zlnum,Zlpsr,Zlfdr,ZlGM,Zllh,ZPI,Zonum,Zpower0,ZIerror0,ZAUC0)
    ENoutput
    RwtMTLENoutput=c(Zltnum,Zltpsr,Zltfdr,ZltGM,Zltlh,ZtPI,Zotnum,Zpower,ZIerror,ZAUCMTL,Zotnum_2,Zpower_2,ZIerror_2)
    RwtMTLENoutput
    RawMTLENoutput=c(Zlt2num,Zlt2psr,Zlt2fdr,Zlt2GM,Zlt2lh,Zt2PI,Zot2num,Zpower2,ZIerror2)
    RawMTLENoutput


     ZCensorRate

### median

     MENoutput=c(Mlnum,Mlpsr,Mlfdr,MlGM,Mllh,MPI,Monum,Mpower0,MIerror0,MAUC0)
    MENoutput
    MRwtMTLENoutput=c(Mltnum,Mltpsr,Mltfdr,MltGM,Mltlh,MtPI,Motnum,Mpower,MIerror,MAUCMTL,Motnum_2,Mpower_2,MIerror_2)
    MRwtMTLENoutput
    MRawMTLENoutput=c(Mlt2num,Mlt2psr,Mlt2fdr,Mlt2GM,Mlt2lh,Mt2PI,Mot2num,Mpower2,MIerror2)
    MRawMTLENoutput


   MCensorRate

    ### mad/0.6749

 MaENoutput=c(Malnum,Malpsr,Malfdr,MalGM,Mallh,MaPI,Maonum,Mapower0,MaIerror0,MaAUC0)
    MaENoutput
    MaRwtMTLENoutput=c(Maltnum,Maltpsr,Maltfdr,MaltGM,Maltlh,MatPI,Maotnum,Mapower,MaIerror,MaAUCMTL,Maotnum_2,Mapower_2,MaIerror_2)
    MaRwtMTLENoutput
    MaRawMTLENoutput=c(Malt2num,Malt2psr,Malt2fdr,Malt2GM,Malt2lh,Mat2PI,Maot2num,Mapower2,MaIerror2)
    MaRawMTLENoutput

