#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
#> Supporting Material: EQSIM analysis 
#> Authors: H. Winker, M. Cardinale, R. Sharma, L.T. Kell, I. Mosqueira, C. Griffiths
#> Title: "Assessing the progress of stock rebuilding in the Northeast Atlantic 
#>           against levels that can produce Maximum Sustainable Yield"
#> Published in Fish and Fishery 
#> Code developed by: M. Cardinale (SLU)  
#>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

#Load the library
#devtools::install_github("flr/FLCore")
#devtools::install_github("ices-tools-prod/msy")

library(devtools)
library(FLCore)
library(msy) # EQSIM


stks <- FLStocks(lapply(stksR,function(x){
  stk <- as(x,"FLStock")
  attr(stk,"stockinfo") <- x@stockinfo
  attr(stk,"benchmark") <- x@benchmark
  return(stk)
}))
save(stks, file="../data/ICES_stksinp_n82.Rdata")

# Load stocks input as FLR object
# setwd to working directory to source file ~ICES-Bmsy-FaF/data/ICES_stksinp_n82.Rdata
load("../data/ICES_stksinp_n82.Rdata",verbose = T)


####Create a list from single model for the Baltic Sea example
stk.ls <- c("her.27.3031","cod.27.24-32","cod.27.22-24","her.27.25-2932","her.27.28","her.27.25-2932","ple.27.21-23","spr.27.22-32")
stks.bs <- stks[stk.ls]


# Run Baltic Sea subset of stocks through EQSIM in loop 
for (i in 1:length(stks.bs)) {
  
  stk <- stks.bs[[i]]
  
  ###Set reference points based on official benchmarks
  Fmsy <- an(stk@benchmark["Fmsy"])
  Blim <- an(stk@benchmark["Blim"])
  Bpa <- an(stk@benchmark["Bpa"])
  Btrigger <- an(stk@benchmark["Btrigger"])
  
  ####Create the Hockey stick SR with Blim
  stk <- trim(stk, year=c(range(stk)["minyear"]:range(stk)["maxyear"]))
  segreg3  <- function(ab, ssb) log(ifelse(ssb >= Blim, ab$a * Blim, ab$a * ssb))
  
  #Fit the SR data using only segreg3 model (HS at blim)
  FIT <- eqsr_fit(stk, nsamp = 1000, models = c("segreg3"))
  
  #Fit the SR data with the 3 models combined 
  #FIT <- eqsr_fit(stk, nsamp = 1000, models = c("Bevholt", "Ricker","Segreg"))
  #FIT <- eqsr_fit(Herring, nsamp = 1000, models = c("Bevholt"))
  
  #Define the F range to which run the simulations
  Fscan = seq(0,2,len=60)
  
  #Setting the the biology and selectivity
  bio.years = c(range(stk)["maxyear"]-2,range(stk)["maxyear"])
  sel.years = c(range(stk)["maxyear"]-2,range(stk)["maxyear"])
  
##STEP 1 with Segmented Regression with breakpoint at Blim
SIM <- eqsim_run(FIT, bio.years=bio.years, bio.const=FALSE, sel.years=sel.years, sel.const=FALSE, Fcv=0.212, Fphi=0.423, Blim=Blim, Bpa=Bpa, Fscan = Fscan, verbose=TRUE, process.error=TRUE, Nrun=200, Btrigger=Btrigger, rhologRec=TRUE, SSBcv=0, extreme.trim=c(0.05,0.95))
  
  #Estimate the catch at the equilibrium for a given F
  SIM$rbp$p50[SIM$rbp$variable=="Catch"]
  
  ##Interpolate to get equilbrium catches
  Catch_50perc <- approx(Fscan, SIM$rbp$p50[SIM$rbp$variable=="Catch"], xout=seq(min(Fscan),max(Fscan),length=40))
  
  ##The following matrix now gives you the 50th percentile for a whole long list of F values from which then you can extract the one corresponding to Fmsy
  cbind(Catch_50perc$x, Catch_50perc$y)
  
  Catchequi <- Catch_50perc$y[which.min(abs(Catch_50perc$x - Fmsy))]
  ## where Fmsy is the preliminary Fmsy value identified in Step 1 above
  
  ##The following lines gives you the 50th percentile of the long-term equilibrium distribution of SSB (BMSY) when fishing at each of the values F in the vector Fscan (the vector of F values you used to run EqSim):
  SIM$rbp$p50[SIM$rbp$variable=="Spawning stock biomass"]
  
  ##Interpolate to get this for more F values:
  BF_50perc1 <- approx(Fscan, SIM$rbp$p50[SIM$rbp$variable=="Spawning stock biomass"], xout=seq(min(Fscan),max(Fscan),length=1000))
  
  ##The following matrix now gives you the 50th percentile for a whole long list of F values from which then you can extract the one corresponding to Fmsy and F0
  cbind(BF_50perc1$x, BF_50perc1$y)
  
  #Estimate B0 and BMSY
  F.0 = 0
  BMSY <- BF_50perc1$y[which.min(abs(BF_50perc1$x - Fmsy))]
  B0 <- BF_50perc1$y[which.min(abs(BF_50perc1$x - F.0))]
  BMSY/B0
  Flim <- BF_50perc1$x[which.min(abs(BF_50perc1$y-Blim))]
  
  frp <- data.frame(benchmarks[i,3:8])
  eqsim = data.frame(stock= stk@name)
  
  eqsim$Catchequi <- Catchequi
  
  eqsim$BMSY <-BMSY
  
  eqsim$B0 <-B0
  eqsim$FmsyMedianC <- SIM$refs_interval["FmsyMedianC"]
  eqsim$FmsyMedianL <- SIM$refs_interval["FmsyMedianL"]
  eqsim$F5percRiskBlim <- SIM$refs_interval["F5percRiskBlim"]
  # assign reference points from EQSIM
  attr(stks.bs[[i]],"eqsim") = eqsim
  
}

# Test if EQSIM outcomes are assigned
stks.bs$her.27.3031@eqsim


