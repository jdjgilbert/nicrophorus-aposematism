

TRAIT     		HERITABILITY  HERITABILITY LCI  HERITABILITY HCI
stripe1   		0.9311069	0.8070718	0.9691263
stripe2	  		0.9693553	0.8599488	0.9884443
orange_total		0.9693553	0.8599488 	0.9884443
elytrasize		0.9920369	0.8478869	0.9990856  
Lum_2stripe		0.5812081	0.2647815	0.7716477
UV_2stripe		0.4117868	0.2145945	0.6648513
SW_2stripe		0.4117868	0.2145945	0.6648513
MW_2stripe		0.5830845	0.3994038	0.9937514
LW_2stripe		0.4864544	0.2355224	0.6697852
DBL_2stripe		0.4810725	0.249275	0.8066984
Meanref_2stripe		0.4379333	0.2436059	0.7265899
bodysize_mg		0.9944337	0.9258554	0.998655
eclosionfluid		0.4744981	0.2857645	0.7276357
orange_prop		0.9361451	0.820006	0.9651244
firststripe_prop	0.8523544	0.7222129	0.9146386
secondstripe_prop	0.8314307	0.6978735	0.8975755



## Carita colour vs. exudate analysis - heritability

library(MCMCglmm)
library(MasterBayes)

nrow(col <- read.csv('data/Carita2016_colour.csv'))
# 909
nrow(ex  <- read.csv('data/Carita2016_exudate.csv'))
# 224

col[1:5, 1:5]
##      Family Collection treatment individual    sex
##  1    102        F30   0_hours          1 female
##  2    102        F30   0_hours          3 female
##  3    105        F30   0_hours          1 female
##  4    105        F30   0_hours          2   male
##  5    105        F30   0_hours          3   male

ex[1:5, 1:5]
##      Family Collection treatment individual    sex
##  1    107        F30   control          1 female
##  2    107        F30   control          2   male
##  3    107        F30   control          3 female
##  4    107        F30   control          4 female
##  5    107        F30   control          5 female

col$animal <- paste(col$Family, col$individual, sep='-')

head(ped <- cbind.data.frame(id=paste(col$Family, col$individual, sep='-'), dam=col$Family))
##      id dam
##  1 102-1 102
##  2 102-3 102
##  3 105-1 105
##  4 105-2 105
##  5 105-3 105
##  6 105-4 105

#ped$sire <- rep('sire1', nrow(ped))  ## all indivs have the same father
ped$sire <- paste('M', ped$dam, sep='') ## all families have distinct fathers

nrow(ped <- insertPed(ped))
## [1] 1007
nrow(ped <- prunePed(ped, keep=col$animal))
## [1] 929

### select only individuals with firststripe data

nrow(dat.ped <- subset(col, !is.na(firststripe_.mm.)))
## 417
dat.ped$Collection <- factor(dat.ped$Collection)

nrow(ped.dat <- prunePed(ped, keep=dat.ped$animal))
## 518

prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m0<-MCMCglmm(firststripe_.mm. ~ 1,
			random=~animal,
			family="gaussian",
			prior=prior,
			pedigree=ped.dat,
			data=dat.ped,
			nitt=250000,
			burnin=1000,
			thin=200)
### runs fine

autocorr(m0$Sol)  ## fine
##    , , (Intercept)
##    
##    (Intercept)
##    Lag 0      1.000000e+00
##    Lag 200   -1.170837e-02
##    Lag 1000   1.965497e-05
##    Lag 2000   8.177264e-03
##    Lag 10000 -3.671339e-03

plot(m0$Sol)
plot(m0$VCV)

herit <- m0$VCV[, "animal"]/(m0$VCV[, "animal"] + m0$VCV[, "units"])

effectiveSize(herit)  ## 1686
autocorr(herit) ## fine

posterior.mode(herit)
##  var1
##  0.9311069

HPDinterval(herit)
##      lower     upper
##      var1 0.8070718 0.9691263
##      attr(,"Probability")
##      [1] 0.9502008


plot(herit)

## Add generation effect
##    
##    prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
##    m1<-MCMCglmm(firststripe_.mm. ~ Collection,
##    			random=~animal,
##    			family="gaussian",
##    			prior=prior,
##    			pedigree=ped.dat,
##    			data=dat.ped,
##    			nitt=50000,
##    			burnin=1000,
##    			thin=200)
##    
##    autocorr(m1$Sol)  ## fine
##    
##    herit <- m1$VCV[, "animal"]/(m1$VCV[, "animal"] + m1$VCV[, "units"])
##    
##    effectiveSize(herit)  ## 245
##    autocorr(herit) ## fine
##    
##    posterior.mode(herit)
##    ##  var1 
##    ##  0.9454687
##    HPDinterval(herit)
##    ##    lower     upper
##    ##    var1 0.8777721 0.9755327
##    ##    attr(,"Probability")
##    ##    [1] 0.9510204
##    
##    m1$DIC ## 1427.71
##    m0$DIC ## 1428.76
##    
### Take simpler model m0

plot(m0$Sol)
plot(m0$VCV)

### Now do all the other traits - pasted into the table at the top
names(dat.ped)
##   [1] "Family"            "Collection"        "treatment"         "individual"       
##   [5] "sex"               "firststripe_.mm."  "secondstripe_.mm." "orange_total"     
##   [9] "elytrasize"        "Lum_2stripe"       "UV_2stripe"        "SW_2stripe"       
##   [13] "MW_2stripe"        "LW_2stripe"        "DBL_2stripe"       "Meanref_2stripe"  
##   [17] "Lum_black"         "UV_black"          "SW_black"          "MW_black"         
##   [21] "LW_black"          "DBL_black"         "Meanref_black"     "Lum_1stripe"      
##   [25] "UV_1stripe"        "SW_1stripe"        "MW_1stripe"        "LW_1stripe"       
##   [29] "DBL_1stripe"       "Meanref_1stripe"   "bodysize_mg"       "eclosionfluid"    
##   [33] "care_low_high"     "orange_prop"       "firststripe_prop"  "secondstripe_prop"
##   [37] "animal"  


prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m0<-MCMCglmm(secondstripe_.mm. ~ 1,
			random=~animal,
			family="gaussian",
			prior=prior,
			pedigree=ped.dat,
			data=dat.ped,
			nitt=500000,
			burnin=1000,
			thin=1000)

autocorr(m0$Sol)  ## fine
herit <- m0$VCV[, "animal"]/(m0$VCV[, "animal"] + m0$VCV[, "units"])
effectiveSize(herit)  ## 1686
autocorr(herit) ## fine
posterior.mode(herit)
HPDinterval(herit)

##    > posterior.mode(herit)
##    var1 
##    0.986587 
##    > HPDinterval(herit)
##    lower     upper
##    var1 0.8513733 0.9971242
##    attr(,"Probability")
##    [1] 0.9498998


prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m1<-MCMCglmm(orange_total ~ 1,
			random=~animal,
			family="gaussian",
			prior=prior,
			pedigree=ped.dat,
			data=dat.ped,
			nitt=500000,
			burnin=1000,
			thin=1000)

autocorr(m1$Sol)  ## fine
herit <- m1$VCV[, "animal"]/(m1$VCV[, "animal"] + m1$VCV[, "units"])
effectiveSize(herit)  ## 1686
autocorr(herit) ## fine
posterior.mode(herit)
HPDinterval(herit)

##    > posterior.mode(herit)
##    var1 
##    0.9693553 
##    > HPDinterval(herit)
##    lower     upper
##    var1 0.8599488 0.9884443
##    attr(,"Probability")
##    [1] 0.9498998

prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m2<-MCMCglmm(elytrasize ~ 1,
			random=~animal,
			family="gaussian",
			prior=prior,
			pedigree=ped.dat,
			data=dat.ped,
			nitt=500000,
			burnin=1000,
			thin=1000)

autocorr(m2$Sol)  ## fine
herit <- m2$VCV[, "animal"]/(m2$VCV[, "animal"] + m2$VCV[, "units"])
effectiveSize(herit)  ## 1686
autocorr(herit) ## fine
posterior.mode(herit)
HPDinterval(herit)

###     > posterior.mode(herit)
###     var1 
###     0.9920369 
###     > HPDinterval(herit)
###     lower     upper
###     var1 0.8478869 0.9990856
###     attr(,"Probability")
###     [1] 0.9498998
###     > 


prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m3<-MCMCglmm(Lum_2stripe ~ 1,
             random=~animal,
             family="gaussian",
             prior=prior,
             pedigree=ped.dat,
             data=dat.ped,
             nitt=500000,
             burnin=1000,
             thin=1000)

autocorr(m3$Sol)  ## fine
herit <- m3$VCV[, "animal"]/(m3$VCV[, "animal"] + m3$VCV[, "units"])
effectiveSize(herit)  ## 1686
autocorr(herit) ## fine
posterior.mode(herit)
HPDinterval(herit)

##      var1 
##      0.5812081 
##      > HPDinterval(herit)
##      lower     upper
##      var1 0.2647815 0.7716477
##      attr(,"Probability")
##      [1] 0.9498998


prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m4<-MCMCglmm(UV_2stripe ~ 1,
             random=~animal,
             family="gaussian",
             prior=prior,
             pedigree=ped.dat,
             data=dat.ped,
             nitt=500000,
             burnin=1000,
             thin=1000)

autocorr(m4$Sol)  ## fine
herit <- m4$VCV[, "animal"]/(m4$VCV[, "animal"] + m4$VCV[, "units"])
effectiveSize(herit)  ## 499
autocorr(herit) ## fine
posterior.mode(herit)
HPDinterval(herit)

##    > posterior.mode(herit)
##    var1 
##    0.4117868 
##    > HPDinterval(herit)
##    lower     upper
##    var1 0.2145945 0.6648513
##    attr(,"Probability")
##    [1] 0.9498998



prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m5<-MCMCglmm(SW_2stripe ~ 1,
             random=~animal,
             family="gaussian",
             prior=prior,
             pedigree=ped.dat,
             data=dat.ped,
             nitt=500000,
             burnin=1000,
             thin=1000)

autocorr(m5$Sol)  ## fine
herit <- m5$VCV[, "animal"]/(m5$VCV[, "animal"] + m5$VCV[, "units"])
effectiveSize(herit)  ## 499
autocorr(herit) ## fine
posterior.mode(herit)
HPDinterval(herit)

###   > posterior.mode(herit)
###   var1 
###   0.4117868 
###   > HPDinterval(herit)
###   lower     upper
###   var1 0.2145945 0.6648513
###   attr(,"Probability")
###   [1] 0.9498998


prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m6<-MCMCglmm(MW_2stripe ~ 1,
             random=~animal,
             family="gaussian",
             prior=prior,
             pedigree=ped.dat,
             data=dat.ped,
             nitt=500000,
             burnin=1000,
             thin=1000)

autocorr(m6$Sol)  ## fine
herit <- m6$VCV[, "animal"]/(m6$VCV[, "animal"] + m6$VCV[, "units"])
effectiveSize(herit)  ## 388.8679 
autocorr(herit) ## fine
posterior.mode(herit)
HPDinterval(herit)

###     > posterior.mode(herit)
###     var1 
###     0.5830845 
###     > HPDinterval(herit)
###     lower     upper
###     var1 0.3994038 0.9937514
###     attr(,"Probability")
###     [1] 0.9498998


prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m7<-MCMCglmm(LW_2stripe ~ 1,
             random=~animal,
             family="gaussian",
             prior=prior,
             pedigree=ped.dat,
             data=dat.ped,
             nitt=500000,
             burnin=1000,
             thin=1000)

autocorr(m7$Sol)  ## fine
herit <- m7$VCV[, "animal"]/(m7$VCV[, "animal"] + m7$VCV[, "units"])
effectiveSize(herit)  ## 499
autocorr(herit) ## fine
posterior.mode(herit)
HPDinterval(herit)

##      > posterior.mode(herit)
##      var1 
##      0.4864544 
##      > HPDinterval(herit)
##      lower     upper
##      var1 0.2355224 0.6697852
##      attr(,"Probability")
##      [1] 0.9498998


prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m8<-MCMCglmm(DBL_2stripe ~ 1,
             random=~animal,
             family="gaussian",
             prior=prior,
             pedigree=ped.dat,
             data=dat.ped,
             nitt=500000,
             burnin=1000,
             thin=1000)

autocorr(m8$Sol)  ## fine
herit <- m8$VCV[, "animal"]/(m8$VCV[, "animal"] + m8$VCV[, "units"])
effectiveSize(herit)  ## 499
autocorr(herit) ## fine
posterior.mode(herit)
HPDinterval(herit)

#     > posterior.mode(herit)
#     var1 
#     0.4810725 
#     > HPDinterval(herit)
#     lower     upper
#     var1 0.249275 0.8066984
#     attr(,"Probability")
#     [1] 0.9498998


prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m9<-MCMCglmm(Meanref_2stripe ~ 1,
             random=~animal,
             family="gaussian",
             prior=prior,
             pedigree=ped.dat,
             data=dat.ped,
             nitt=500000,
             burnin=1000,
             thin=1000)

autocorr(m9$Sol)  ## fine
herit <- m9$VCV[, "animal"]/(m9$VCV[, "animal"] + m9$VCV[, "units"])
effectiveSize(herit)  ## 499
autocorr(herit) ## fine
posterior.mode(herit)
HPDinterval(herit)
###   > posterior.mode(herit)
###   var1 
###   0.4379333 
###   > HPDinterval(herit)
###   lower     upper
###   var1 0.2436059 0.7265899
###   attr(,"Probability")
###   [1] 0.9498998


prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m10<-MCMCglmm(bodysize_mg ~ 1,
             random=~animal,
             family="gaussian",
             prior=prior,
             pedigree=ped.dat,
             data=dat.ped,
             nitt=500000,
             burnin=1000,
             thin=1000)

autocorr(m10$Sol)  ## fine
herit <- m10$VCV[, "animal"]/(m10$VCV[, "animal"] + m10$VCV[, "units"])
effectiveSize(herit)  ## 499
autocorr(herit) ## fine
posterior.mode(herit)
HPDinterval(herit)

##  > posterior.mode(herit)
##  var1 
##  0.9944337 
##  > HPDinterval(herit)
##  lower    upper
##  var1 0.9258554 0.998655
##  attr(,"Probability")
##  [1] 0.9498998


prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m11<-MCMCglmm(eclosionfluid ~ 1,
              random=~animal,
              family="gaussian",
              prior=prior,
              pedigree=ped.dat,
              data=dat.ped,
              nitt=500000,
              burnin=1000,
              thin=1000)

autocorr(m11$Sol)  ## fine
herit <- m11$VCV[, "animal"]/(m11$VCV[, "animal"] + m11$VCV[, "units"])
effectiveSize(herit)  ## 499
autocorr(herit) ## fine
posterior.mode(herit)
HPDinterval(herit)

##    > posterior.mode(herit)
##    var1 
##    0.4744981 
##    > HPDinterval(herit)
##    lower     upper
##    var1 0.2857645 0.7276357
##    attr(,"Probability")
##    [1] 0.9498998



prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m12<-MCMCglmm(orange_prop ~ 1,
              random=~animal,
              family="gaussian",
              prior=prior,
              pedigree=ped.dat,
              data=dat.ped,
              nitt=500000,
              burnin=1000,
              thin=1000)

autocorr(m12$Sol)  ## fine
herit <- m12$VCV[, "animal"]/(m12$VCV[, "animal"] + m12$VCV[, "units"])
effectiveSize(herit)  ## 499
autocorr(herit) ## fine
posterior.mode(herit)
HPDinterval(herit)

###   > posterior.mode(herit)
###   var1 
###   0.9361451 
###   > HPDinterval(herit)
###   lower     upper
###   var1 0.820006 0.9651244
###   attr(,"Probability")
###   [1] 0.9498998



prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m13<-MCMCglmm(firststripe_prop ~ 1,
              random=~animal,
              family="gaussian",
              prior=prior,
              pedigree=ped.dat,
              data=dat.ped,
              nitt=500000,
              burnin=1000,
              thin=1000)

autocorr(m13$Sol)  ## fine
herit <- m13$VCV[, "animal"]/(m13$VCV[, "animal"] + m13$VCV[, "units"])
effectiveSize(herit)  ## 499
autocorr(herit) ## fine
posterior.mode(herit)
HPDinterval(herit)

##    > posterior.mode(herit)
##    var1 
##    0.8523544 
##    > HPDinterval(herit)
##    lower     upper
##    var1 0.7222129 0.9146386
##    attr(,"Probability")
##    [1] 0.9498998


prior<-list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
m14<-MCMCglmm(secondstripe_prop ~ 1,
              random=~animal,
              family="gaussian",
              prior=prior,
              pedigree=ped.dat,
              data=dat.ped,
              nitt=500000,
              burnin=1000,
              thin=1000)

autocorr(m14$Sol)  ## fine
herit <- m14$VCV[, "animal"]/(m14$VCV[, "animal"] + m14$VCV[, "units"])
effectiveSize(herit)  ## 499
autocorr(herit) ## fine
posterior.mode(herit)
HPDinterval(herit)

##    > posterior.mode(herit)
##    var1 
##    0.8314307 
##    > HPDinterval(herit)
##    lower     upper
##    var1 0.6978735 0.8975755
##    attr(,"Probability")
##    [1] 0.9498998

