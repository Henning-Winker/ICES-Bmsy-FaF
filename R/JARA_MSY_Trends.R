#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
#> Supporting Material: EQSIM analysis 
#> Authors: H. Winker, M. Cardinale, R. Sharma, L.T. Kell, I. Mosqueira, C. Griffiths
#> Title: "Assessing the progress of stock rebuilding in the Northeast Atlantic 
#>           against levels that can produce Maximum Sustainable Yield"
#> Published in Fish and Fishery 
#> Code developed by: H. Winker (FAO), henning.winker@fao.org  
#>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

library(JARA)
library(reshape2)
library(FLCore)
library(ggplot2)
library(FLRef)
library(gridExtra)
library(grid)
library(ggpubr)


# Use equisim 
load("../data/ICES_stksR_n82.Rdata",verbose=T)
th <- theme(axis.text.x  = element_text(angle=90, vjust=0.5), panel.grid.minor = element_blank())
saveplot = TRUE
stks = stksR

# Select stocks with Fmsy
nt = do.call(c,lapply(stks, function(x)x@eqsim[["BMSY"]]))
nam = stks@names[!is.na(nt)]
stks1 =stks[nam]

refyr = do.call(c,lapply(stks, function(x)range(x)["maxyear"]))

# endyr
stks= FLStocks(lapply(stks1,function(x){
  if(range(x)["maxyear"]>2020) x =window(x,end=2020)
  if(range(x)["minyear"]<1970) x =window(x,start=1970)
  return(x)
}))


Tab1 = do.call(rbind,lapply(stks,function(x){
  data.frame(x@stockinfo,productivity=x@productivity@desc,Fmsy=an(x@refpts["Fmsy"]),Bmsy=an(x@refpts["Bmsy"]),Btrigger=an(x@refpts["Btrigger"]),Blim=an(x@refpts["Blim"]),B0=an(x@refpts["B0"]),Bcur=an(tail(ssb(x),1)),Fcur=an(tail(fbar(x),1)))
  
}))

# output
dir.create("../output",showWarnings  =F)
write.csv(Tab1,"../output/Tab1.csv",row.names = F)
yrs= do.call(c,lapply(stks,function(x){range(x)["minyear"]}))


#### Overview plots
# Make stock long table
dat = do.call(rbind,lapply(stks, function(x){
  data.frame(data.frame(x@stockinfo,productivity=x@productivity@desc),Year = dims(x)$minyear:dims(x)$maxyear ,
             FFmsy = c(fbar(x)/x@refpts["Fmsy"]), 
             BBmsy=c(ssb(x)/x@refpts["Bmsy"]),
             BBtri = c(ssb(x)/x@refpts["Btrigger"]),
             BBlim = c(ssb(x)/x@refpts["Blim"]),
             # Compute ABItgt  
             ABImsy = c(ABItgt(as(x,"FLStock"),ftgt=x@refpts[[1]])))
  
}))

write.csv(dat,"../output/ices.status.ts.csv",row.names = F)

# Create ratios of stock trajectories
stky = aggregate(FFmsy~Year,dat,length)
stky$Coverage = stky$FFmsy/max(stky$FFmsy)*100

df = stky
if(saveplot) png("../output/jara_ts_coverage.jpg", 2000, 1200, res=300)
pcover = ggplot(df, aes(x=Year)) + 
  expand_limits(y=c(0,2)) + 
  ylab("Number of stocks") + xlab("Year")+ 
  geom_point(data=df[df$Year==2020,],aes(x=Year,y=FFmsy),bg="darkgrey",pch=21,size=3.5)+
  geom_point(data=df,aes(Year,FFmsy,color=Coverage),pch=16)+
  scale_x_continuous(breaks = seq(1970,2020,5))+
  scale_y_continuous(breaks = seq(0,75,5))+
  theme_bw()+
  th+ theme(legend.position = "right",
            legend.box.spacing = unit(0.1, 'cm'))+
  geom_hline(yintercept = 73,lty=2,col=2)
pcover 
dev.off()

# Make StockInfo Table
info = Tab1 
info = info[info$stock%in%stks@names,]

stkn = aggregate(Fmsy~info$common.name+info$order,info,length)
colnames(stkn) = c("species","order","n")
stkn = stkn[order(factor(stkn$species),decreasing=T),]
stkn = stkn[order(stkn$n),]
stkn$species = factor(stkn$species,levels=unique(stkn$species)) 

stka = aggregate(Fmsy~info$ecoregion,info,length)
colnames(stka) = c("area","n")
stka = stka[order(stka$n),]
stka$area = factor(stka$area,levels=unique(stka$area)) 

stkm = aggregate(Fmsy~model,info,length)
colnames(stkm) = c("model","n")
stkm = stkm[order(stkm$n),]
stkm$model = factor(stkm$model,levels=unique(stkm$model)) 


ps = ggplot(stkn,aes(species,n))+
  geom_bar(stat = 'identity', alpha = 0.6, fill = 'grey', col = 'black')+
  theme_classic()+
  coord_flip()+
  ylab('No. of stocks')+
  xlab('')+
  scale_y_continuous(breaks = seq(0,12,1), position = 'right')+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(hjust = 0))

pm = ggplot(stkm,aes(model,n))+
  geom_bar(stat = 'identity', alpha = 0.6, fill = 'grey', col = 'black')+
  theme_classic()+coord_flip()+
  ylab('No. of stocks')+
  xlab('')+
  scale_y_continuous(breaks = seq(0,100,2), position = 'right')+
  theme(text = element_text(size = 10),
  axis.text.y = element_text(hjust = 0,margin = margin(40, 40, 40, 40)*0.25))

pa = ggplot(stka,aes(area,n))+
  geom_bar(stat = 'identity', alpha = 0.6, fill = 'grey', col = 'black')+
  theme_classic()+coord_flip()+
  ylab('No. of stocks')+
  xlab('')+
  scale_y_continuous(breaks = seq(0,20,1), position = 'right')+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(hjust = 0))

if(saveplot) png("../output/jara_Status_dstat.jpg", 2200, 1800, res=300)
ggarrange(ps,                                                 # First row with scatter plot
          ggarrange(pm, pa, ncol = 1, labels = c("(b)", "(c)"), font.label=list(size=12)), # Second row with box and dot plots
          nrow = 1, 
          labels = "(a)",
          font.label=list(size=12),
          widths = c(1, 1.6)
) 
dev.off()


### Make Ratio plot
x = stks[[1]]
stkr= do.call(rbind,lapply(stks,function(x){
  i = tail(dat[dat$stock%in%x@name,],1)
  j = info[info$stock%in%x@name,]
  rbind(
  data.frame(stock=i$stock[1],region = i$ecoregion,name=i$common.name,order=i$order,
             quant="F[2020]/F[MSY]",y=i$FFmsy),
  data.frame(stock=x@name,region = i$ecoregion,name=i$common.name,order=i$order,
             quant="B[2020]/B[lim]",y=i$BBlim),
  data.frame(stock=x@name,region = i$ecoregion,name=i$common.name,order=i$order,
             quant="B[2020]/B[trigger]",y=i$BBtri),
  data.frame(stock=x@name,region = i$ecoregion,name=i$common.name,order=i$order,
             quant="B[2020]/B[MSY]",y=i$BBmsy),
  data.frame(stock=x@name,region = i$ecoregion,name=i$common.name,order=i$order,
             quant="ABI[MSY]",y=i$ABImsy),
  data.frame(stock=x@name,region = i$ecoregion,name=i$common.name,order=i$order,
             quant="B[MSY]/B[trigger]",y=j$Bmsy/j$Btrigger),
  data.frame(stock=x@name,region = i$ecoregion,name=i$common.name,order=i$order,
             quant="B[MSY]/B[lim]",y=j$Bmsy/j$Blim),
  data.frame(stock=x@name,region = i$ecoregion,name=i$common.name,order=i$order,
             quant="B[trigger]/B[lim]",y=j$Btrigger/j$Blim),
  data.frame(stock=x@name,region = i$ecoregion,name=i$common.name,order=i$order,
             quant="B[trigger]/B[0]",y=j$Btrigger/j$B0),
  data.frame(stock=x@name,region = i$ecoregion,name=i$common.name,order=i$order,
             quant="B[MSY]/B[0]",y=j$Bmsy/j$B0)
  )
  
}))

stkr = stkr[order(stkr$order,stkr$name,decreasing=T),]

stkr$stock = factor(stkr$stock,levels=unique(stkr$stock)) 

if(saveplot) png("../output/jara_Status_S1.jpg", 2200, 2000, res=300)
ggplot(stkr[stkr$quant%in%c("B[2020]/B[MSY]","F[2020]/F[MSY]"),],aes(stock,y))+
  geom_segment(aes(x=stock,xend=stock, y= 1, yend = y))+
  geom_point(aes(col=region))+
  theme_bw()+facet_grid(~quant,labeller = label_parsed)+
  coord_flip()+
  ylab('Ratio')+
  xlab('')+geom_hline(yintercept = 1)+
  scale_y_continuous(breaks = c(0.2,c(seq(0,1,0.5),2,4,6,8)), position = 'right',trans="sqrt")+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(hjust = 0),
        legend.title = element_blank(),
        axis.text.x=element_text(size=6))
dev.off()

df = stkr[stkr$quant%in%c("B[MSY]/B[trigger]","B[MSY]/B[lim]","B[trigger]/B[lim]"),]
nd =aggregate(Bcur~order,info,length)
nd = nd[order(factor(nd$Bcur),decreasing=T),]
df1 = df
df1$order="Combined"
df=df[df$order%in%nd[1:5,1],]
df2 = rbind(df,df1)
df2$order = factor(df2$order,levels=unique(c(nd[,1],"Combined")))
df2$quant = factor(df2$quant ,levels=c("B[MSY]/B[trigger]","B[MSY]/B[lim]","B[trigger]/B[lim]"))

dir.create("../results",showWarnings = F)
if(saveplot) png("../results/Appendix_jara_ViolinBB.jpg", 2500, 1500, res=300)
ggplot(df2,aes(x =order,  y=y,fill=order))+ 
  geom_violin(width=1.1,fill="darkgrey")+theme_bw()+
  facet_grid(~quant, labeller = label_parsed)+
  theme(legend.position = "none")+ylab("Ratio")+xlab("Order")+
  geom_boxplot(size=0.5,position=position_dodge(width=1),width=0.2,fill="lightgrey")+
 geom_point(aes(group=order),col="black",fill="white",pch=21,size=1.5,alpha=1)+
 geom_hline(yintercept = 1,linetype="dashed")+th+
scale_y_continuous(breaks=seq(0,5,0.25),limits=c(0,5))
 
dev.off()


orders = c("Gadiformes","Pleuronectiformes","Clupeiformes")


#-----------------------------------------------------------
# JARA Analysis
#-----------------------------------------------------------

refit = FALSE
saveplot = TRUE

dat = dat[dat$Year>=1970,]
rownames(dat)= 1:nrow(dat)
dat$lBBmsy = dat$BBmsy/0.8
dat$uBBmsy = dat$BBmsy/1.2


# Make input data sets for JARA 
d. = reshape2::dcast(dat,Year~stock,value.var="BBmsy",fun.aggregate=sum)
d.[d.==0] <- NA
bbmsy=d.
ncol(d.)

d. = reshape2::dcast(dat,Year~stock,value.var="lBBmsy",fun.aggregate=sum)
d.[d.==0] <- NA
lbbmsy=d.
ncol(d.)

d. = reshape2::dcast(dat,Year~stock,value.var="uBBmsy",fun.aggregate=sum)
d.[d.==0] <- NA
ubbmsy=d.
ncol(d.)



d. = reshape2::dcast(dat,Year~stock,value.var="FFmsy",fun.aggregate=sum)
d.[d.==0] <- NA
ffmsy=d.
ffmsy[ffmsy$Year==1990,"aru.27.5a14"] =NA

d. = reshape2::dcast(dat,Year~stock,value.var="BBtri",fun.aggregate=sum)
d.[d.==0] <- NA
bbtri=d.

d. = reshape2::dcast(dat,Year~stock,value.var="BBlim",fun.aggregate=sum)
d.[d.==0] <- NA
bblim=d.

d. = reshape2::dcast(dat,Year~stock,value.var="ABImsy",fun.aggregate=sum)
d.[d.==0] <- NA
abimsy=d.

write.csv(dat,"../data/ices.ts.cat1.csv",row.names = F)
# General settings
# minimum obs error
obsmin=0.05
# Time block for trends 2000
tbf = 2015
savegraph = TRUE


dir.create("../rdata",showWarnings = F)
set.seed(1507)
do.call(rbind,lapply(stks,function(x){
  max(an(dimnames(x)$year))
}))

if(refit){
# Fit B/Bmsy trends
bb =build_jara (I = bbmsy,model.type="mixed.trends",assessment="ices",scenario = "BBmsy",fixed.obsE =obsmin,timeblock=2015,backcast.na = 2000)
fitb = fit_jara(bb,quickmcmc = T)
# Fit F/Fmsy trends
bf =build_jara (I = ffmsy,model.type="mixed.trends",assessment="ices",scenario = "FFmsy",fixed.obsE =obsmin,timeblock=2015,backcast.na = 2000)
fitf = fit_jara(bf,quickmcmc = T)
# Fit B/Btri trends
bt =build_jara (I = bbtri,model.type="mixed.trends",assessment="ices",scenario = "BBtri",fixed.obsE =obsmin,timeblock=2015,backcast.na = 2000)
fitt = fit_jara(bt,quickmcmc = T)
# Fit AMBmsy
abi =build_jara (I = abimsy,model.type="mixed.trends",assessment="ices",scenario = "ABImsy", fixed.obsE = 0.2,proc.pen = 0.5,timeblock=2015,backcast.na = 2000)
#abimsy[abimsy<0.05] = NA
fita = fit_jara(abi,quickmcmc = T)

save(fitb,fitf,fitt,fita,file="../rdata/jarafits.ices.rdata")


} else {
  load(file="../rdata/jarafits.ices.rdata")
}
stn = names(stks) 
# Check fits
if(saveplot) png("../results/jara_fits_abi.jpg", 2000, 2000, res=300)
df = dfidx(fita)
ggplot(df$i,aes(x=yr,y=mu,group=name))+ 
  geom_ribbon(aes(ymin=lci, ymax=uci,fill=name),alpha=0.3, linetype=1,color=NA)+theme_bw()+
  geom_line(cex=0.2)+theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y = element_text(size = 7),strip.text.x = element_text(size = 5))+
  geom_point(data=df$o,aes(year,obs,color=name))+ylab(expression(ABI/ABI[MSY]))+
  facet_wrap(~name,ncol=6,scale="free_y")+xlab("Year")
  geom_hline(yintercept = 1,linetype=3)
if(saveplot) dev.off()

if(saveplot) png("../results/jara_fits_bbmsy.jpg", 2000, 2000, res=300)
df = dfidx(fitb)

ggplot(df$i,aes(x=yr,y=mu,group=name))+ 
  geom_ribbon(aes(ymin=lci, ymax=uci,fill=name),alpha=0.3, linetype=1,color=NA)+theme_bw()+
  geom_line(cex=0.2)+theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y = element_text(size = 7),strip.text.x = element_text(size = 5))+
  geom_point(data=df$o,aes(year,obs,color=name))+ylab(expression(B/B[MSY]))+
  facet_wrap(~name,ncol=6,scale="free_y")+xlab("Year")+
  geom_hline(yintercept = 1,linetype=3)
dev.off()

if(saveplot) png("../results/jara_fits_ffmsy.jpg", 2000, 2000, res=300)
df = dfidx(fitf)
df$i = df$i[df$i$name%in%stn,]
df$o = df$o[df$o$name%in%stn,]
ggplot(df$i,aes(x=yr,y=mu,group=name))+ 
  geom_ribbon(aes(ymin=lci, ymax=uci,fill=name),alpha=0.3, linetype=1,color=NA)+theme_bw()+
  geom_line(cex=0.2)+theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y = element_text(size = 7),strip.text.x = element_text(size = 5))+
  geom_point(data=df$o,aes(year,obs,color=name))+ylab(expression(F/F[MSY]))+
  facet_wrap(~name,ncol=6,scale="free_y")+xlab("Year")+
  geom_hline(yintercept = 1,linetype=3)
if(saveplot) dev.off()

if(saveplot) png("../results/jara_fits_bbtri.jpg", 2000, 2000, res=300)
df = dfidx(fitt)
df$i = df$i[df$i$name%in%stn,]
df$o = df$o[df$o$name%in%stn,]
ggplot(df$i,aes(x=yr,y=mu,group=name))+ 
  geom_ribbon(aes(ymin=lci, ymax=uci,fill=name),alpha=0.3, linetype=1,color=NA)+theme_bw()+
  geom_line(cex=0.2)+theme(legend.position = "none", axis.text.x=element_blank(), axis.text.y = element_text(size = 7),strip.text.x = element_text(size = 5))+
  geom_point(data=df$o,aes(year,obs,color=name))+ylab(expression(B/B[trigger]))+
  facet_wrap(~name,ncol=6,scale="free_y")+xlab("Year")+
  geom_hline(yintercept = 1,linetype=3)
if(saveplot) dev.off()

# Main Orders
if(refit){
  ordB = lapply(orders, function(x){
    inp = bbmsy[,c(1,which(do.call(rbind,lapply(stks,function(y)y@stockinfo[["order"]]))%in%x)+1)]
    # Fit B/Bmsy trends
    bb =build_jara (I = inp,model.type="mixed.trends",assessment="ices",scenario = "BBmsy",fixed.obsE =obsmin,timeblock=2015,backcast.na = 2000)
    fit_jara(bb,quickmcmc = T)
  })
  ordF = lapply(orders, function(x){
    inp = ffmsy[,c(1,which(do.call(rbind,lapply(stks,function(y)y@stockinfo[["order"]]))%in%x)+1)]
    # Fit B/Bmsy trends
    bb =build_jara (I = inp,model.type="mixed.trends",assessment="ices",scenario = "FFmsy",fixed.obsE =obsmin,timeblock=2015,backcast.na = 2000)
    fit_jara(bb,quickmcmc = T)
  })
  
  ordBt = lapply(orders, function(x){
    inp = bbtri[,c(1,which(do.call(rbind,lapply(stks,function(y)y@stockinfo[["order"]]))%in%x)+1)]
    # Fit B/Bmsy trends
    bb =build_jara (I = inp,model.type="mixed.trends",assessment="ices",scenario = "FFmsy",fixed.obsE =obsmin,timeblock=2015,backcast.na = 2000)
    fit_jara(bb,quickmcmc = T)
  })
  
  
  ordABI = lapply(orders, function(x){
    inp = abimsy[,c(1,which(do.call(rbind,lapply(stks,function(y)y@stockinfo[["order"]]))%in%x)+1)]
    # Fit B/Bmsy trends
    bb =build_jara (I = inp,model.type="mixed.trends",assessment="ices",scenario = "ABImsy",fixed.obsE =0.2,proc.pen = 0.5,timeblock=2015,backcast.na = 2000)
    fit_jara(bb,quickmcmc = T)
  })
  # Fit B/Btri trends
  save(ordB,ordF,ordABI,ordBt,file="rdata/jarafits.orders.rdata")
} else {
  load("../rdata/jarafits.orders.rdata",verbose=T)
  
}

mb =cbind(mixed.trend(fitb,run="B/Bmsy"),quant="Spawning Biomass") # expected mean
mf =cbind(mixed.trend(fitf,run="F/Fmsy"),quant="Fishing Mortality") # expected mean
mt =cbind(mixed.trend(fitt,run="B/Btrigger"),quant="Spawning Biomass") # expected mean
ma =cbind(mixed.trend(fita,run="ABImsy"),quant="Spawning Biomass") # expected mean
md = rbind(mb,mf,mt,ma)
md$run = factor(md$run,levels=c("F/Fmsy","B/Bmsy","B/Btrigger","ABImsy"))
# Probability
pb =cbind(mixed.trend(fitb,run="B>Bmsy",type="pr"),quant="Spawning Biomass") # expected mean
pf =cbind(mixed.trend(fitf,run="F<Fmsy",type="pr"),quant="Fishing Mortality") # expected mean
pf[,4:8] =1- pf[,4:8]  
pt =cbind(mixed.trend(fitt,run="B>Btrigger",type="pr"),quant="Spawning Biomass") # expected mean
pa =cbind(mixed.trend(fita,run="ABImsy>1",type="pr"),quant="Spawning Biomass") # expected mean

# Fits
idxb = dfidx(fitb,run="idxb") # B/Bmsy
idxf = dfidx(fitf,run="idxf") # F/Fmsy
idxt = dfidx(fitt,run="idxt") # Trigger
idxa = dfidx(fita,run="idxa") # ABImsy

idxb$o$pr = ifelse(idxb$o$obs>=1,1,0)
idxf$o$pr = ifelse(idxf$o$obs<=1,1,0)
idxt$o$pr = ifelse(idxt$o$obs>=1,1,0)
idxa$o$pr = ifelse(idxa$o$obs>=1,1,0)

pb$Obs = aggregate(pr~year,idxb$o,mean)[,2]
pf$Obs = aggregate(pr~year,idxf$o,mean)[,2]
pt$Obs = aggregate(pr~year,idxt$o,mean)[,2]
pa$Obs = aggregate(pr~year,idxa$o,mean)[,2]

pd = rbind(pb,pf,pt,pa)
pd$run = factor(pd$run,levels=c("F<Fmsy","B>Bmsy","B>Btrigger","ABImsy>1"))


dir.create("../results",showWarnings = F)



# Color shading
df = pd
pstpc1 = ggplot(df, aes(x=Year,group=run)) + 
  geom_ribbon(data=df,aes(ymin = `5%`, ymax = `95%`,fill=run), alpha=0.15) + 
  geom_ribbon(data=df,aes(ymin = `25%`, ymax = `75%`,fill=run), alpha=0.2)+
  geom_line(aes(y=`50%`,color=run),linewidth=0.6) + 
  expand_limits(y=c(0,1)) + 
  geom_hline(yintercept = c(0.25,0.5,0.75), linetype=2,col="darkgrey")+ #+ facet_wrap(~quant)+
  ylab("Proportion") + xlab("Year")+ 
  geom_point(aes(Year,Obs,color=run),pch=21,bg="white",cex=1.5)+
  geom_line(aes(y=`50%`,color=run))+
  theme(legend.position = "right") + theme_bw()+
  scale_x_continuous(breaks=seq(min(df$Year), max(df$Year), 5))+ #, labels=as.character(vp)
  th + 
  scale_color_manual(values=c("black", "darkgreen","darkorange","blue"))+
  scale_fill_manual(values=c("black", "darkgreen","darkorange","blue"))+
  theme(legend.title = element_blank(),legend.position = "right",
        legend.box.spacing = unit(0.1, 'cm'))


df = md
pstmc1 = ggplot(df, aes(x=Year,group=run)) + 
  geom_ribbon(data=df,aes(ymin = `5%`, ymax = `95%`,fill=run), alpha=0.15) + 
  geom_ribbon(data=df,aes(ymin = `25%`, ymax = `75%`,fill=run), alpha=0.2)+
  geom_line(aes(y=`50%`,color=run),linewidth=0.6) + 
  expand_limits(y=c(0,1)) + 
  geom_hline(yintercept = 1, linetype=2,col="black")+ #+ facet_wrap(~quant)+
  ylab("Ratio") + xlab("Year")+ 
  geom_point(aes(Year,Obs,color=run),pch=21,bg="white",cex=1.5)+
  geom_line(aes(y=`50%`,color=run))+
  theme(legend.position = "right") + theme_bw()+
  scale_x_continuous(breaks=seq(min(df$Year), max(df$Year), 5))+ #, labels=as.character(vp)
  scale_y_continuous(breaks=seq(0, 3, .5))+ #, labels=as.character(vp)
  th + 
  scale_color_manual(values=c("black", "darkgreen","darkorange","blue"))+
  scale_fill_manual(values=c("black", "darkgreen","darkorange","blue"))+
  theme(legend.title = element_blank(),legend.position = "right",
        legend.box.spacing = unit(0.1, 'cm'))



if(saveplot) png("../results/jara_Trends_2x1_col.jpg", 2000, 2100, res=300)
ggarrange(pstmc1, pstpc1,ncol=1,align="v")
dev.off()



pf =cbind(mixed.trend(fitf,run="F<Fmsy",type="pr"),quant="Fishing Mortality") # expected mean
pf[,4:8] =1- pf[,4:8]  
pt =cbind(mixed.trend(fitt,run="B>Btrigger",type="pr"),quant="Spawning Biomass") # expected mean
pa =cbind(mixed.trend(fita,run="ABImsy>1",type="pr"),quant="Spawning Biomass") # expected mean

# Fits
idxb = dfidx(fitb,run="idxb",backcast.na = 2000) # B/Bmsy
idxf = dfidx(fitf,run="idxf",backcast.na = 2000) # F/Fmsy
idxt = dfidx(fitt,run="idxt",backcast.na = 2000) # Trigger
idxa = dfidx(fita,run="idxa",backcast.na = 2000) # ABImsy

idxb$o$pr = ifelse(idxb$o$obs>=1,1,0)
idxf$o$pr = ifelse(idxf$o$obs<=1,1,0)
idxt$o$pr = ifelse(idxt$o$obs>=1,1,0)
idxa$o$pr = ifelse(idxa$o$obs>=1,1,0)

pb$Obs = aggregate(pr~year,idxb$o,mean)[,2]
pf$Obs = aggregate(pr~year,idxf$o,mean)[,2]
pt$Obs = aggregate(pr~year,idxt$o,mean)[,2]
pa$Obs = aggregate(pr~year,idxa$o,mean)[,2]


# Orders probs
mob = do.call(rbind,Map(function(x,y){
  pr = cbind(mixed.trend(x,run="B>B[MSY]",type="pr"),quant=paste0(y, " (n ~`=`~",length(x$indices),")"))
  idx = dfidx(x) # B/Bmsy
  idx$o$pr = ifelse(idx$o$obs>=1,1,0) 
  pr$Obs = aggregate(pr~year,idx$o,mean)[,2]
  pr
},x=ordB,y=orders))
mof = do.call(rbind,Map(function(x,y){
  pr = cbind(mixed.trend(x,run="F<F[MSY]",type="pr"),quant=paste0(y, " (n ~`=`~",length(x$indices),")"))
  pr[,4:8] =1- pr[,4:8]  
  idx = dfidx(x) # B/Bmsy
  idx$o$pr = ifelse(idx$o$obs<=1,1,0) 
  pr$Obs = aggregate(pr~year,idx$o,mean)[,2]
  pr
},x=ordF,y=orders))

mot = do.call(rbind,Map(function(x,y){
  pr = cbind(mixed.trend(x,run="B>B[trigger]",type="pr"),quant=paste0(y, " (n ~`=`~",length(x$indices),")"))
  idx = dfidx(x) # B/Bmsy
  idx$o$pr = ifelse(idx$o$obs>=1,1,0) 
  pr$Obs = aggregate(pr~year,idx$o,mean)[,2]
  pr
  
},x=ordBt,y=orders))
moa = do.call(rbind,Map(function(x,y){
  pr = cbind(mixed.trend(x,run="ABI[MSY]>1",type="pr"),quant=paste0(y, " (n ~`=`~",length(x$indices),")"))
  idx = dfidx(x) # B/Bmsy
  idx$o$pr = ifelse(idx$o$obs>=1,1,0) 
  pr$Obs = aggregate(pr~year,idx$o,mean)[,2]
  pr
},x=ordABI,y=orders))

modp = rbind(mof,mot,mob,moa)
modp$quant = factor(modp$quant,levels=unique(modp$quant))
modp$run = factor(modp$run,levels=c("F<F[MSY]","B>B[trigger]","B>B[MSY]","ABI[MSY]>1"))


# Orders geomean
mob = do.call(rbind,Map(function(x,y){
  cbind(mixed.trend(x,run="B/B[MSY]"),quant=paste0(y, " (n ~`=`~",length(x$indices),")"))
  
},x=ordB,y=orders))
mof = do.call(rbind,Map(function(x,y){
   cbind(mixed.trend(x,run="F/F[MSY]"),quant=paste0(y, " (n ~`=`~",length(x$indices),")"))
  
},x=ordF,y=orders))

mot = do.call(rbind,Map(function(x,y){
  cbind(mixed.trend(x,run="B/B[trigger]"),quant=paste0(y, " (n ~`=`~",length(x$indices),")"))
  },x=ordBt,y=orders))
moa = do.call(rbind,Map(function(x,y){
  cbind(mixed.trend(x,run="ABI[MSY]"),quant=paste0(y, " (n ~`=`~",length(x$indices),")"))
  
},x=ordABI,y=orders))

mod = rbind(mof,mot,mob,moa)
mod$quant = factor(mod$quant,levels=unique(mod$quant))
mod$run = factor(mod$run,levels=c("F/F[MSY]","B/B[trigger]","B/B[MSY]","ABI[MSY]"))


df = mod
if(saveplot) png("../results/jara_Status_Order.png", 2500, 2000, res=300)
pord = ggplot(df, aes(x=Year)) + facet_grid(run~quant,scales="free_y",labeller = label_parsed)+
  geom_ribbon(data=df,aes(ymin = `5%`, ymax = `95%`), fill="gray", alpha=0.5) + 
  geom_ribbon(data=df,aes(ymin = `25%`, ymax = `75%`), fill="gray", alpha=0.8)+
  expand_limits(y=c(0,1)) + 
  geom_hline(yintercept = c(1), linetype=2) +
  ylab("Ratio") + xlab("Year")+ 
  geom_point(aes(Year,Obs,color=n),pch=16)+
  geom_line(aes(y=`50%`)) +
  geom_point(data=df[df$Year==2020,],aes(y=`50%`),bg="white",pch=21,size=1)+
  theme(legend.position = "right") + theme_bw()+
  th+ theme(legend.title = element_blank(),legend.position = "right",
            legend.box.spacing = unit(0.1, 'cm'))
pord  
dev.off()

## Create plot with fits showing marginal distributions

# Fits
idxb = dfidx(fitb,run="idxb",backcast.na = 2000) # B/Bmsy
idxf = dfidx(fitf,run="idxf",backcast.na = 2000) # F/Fmsy
idxt = dfidx(fitt,run="idxt",backcast.na = 2000) # Trigger
idxa = dfidx(fita,run="idxt",backcast.na = 2000) # Trigger

df = idxf
# Time-Series
pF= ggplot(df$i,aes(x=yr,y=mu,group=name))+ 
  geom_ribbon(aes(ymin=lci, ymax=pmin(10,uci),fill=name),alpha=0.2, linetype=1,color=NA)+theme_bw()+
  geom_line(aes(color=name),cex=0.4)+theme(legend.position = "none")+
  geom_point(data=df$o,aes(year,obs,color=name),cex=0.6)+ylab(expression(F/F[MSY]))+
  xlab("")+scale_y_continuous(expand = c(0, 0),breaks=seq(0,10,1))+scale_x_continuous(expand = c(0, 0))+
  geom_hline(yintercept = 1,linetype=1,col="black")+
  coord_cartesian(ylim=c(0,10))+theme( plot.margin = unit(c(3.5, 5.5, 0, 0), "points"))

dd = df$i[df$i$yr==max(df$i$yr),]
dF<-ggplot(dd) + 
  geom_histogram(aes(x = mu, y =  ..count..),alpha=0.5,fill="darkgrey",col="gray",position = "stack")+ 
  geom_vline(xintercept=1,col="black")+
  scale_y_continuous(expand = c(0, 0))+scale_x_continuous(expand = c(0, 0))+ coord_flip(xlim=c(0,10))+
  theme(legend.position = "none", 
        axis.title.x = element_text(colour ='NA'), 
        axis.text.x  = element_text(colour ="NA"), 
        axis.ticks.x = element_line(colour ="NA"),
        axis.ticks =   element_line(colour ="NA"),
        axis.title.y = element_blank(), 
        axis.text.y  = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(3.5, 5.5, 0, 0), "points"),
        panel.background = element_rect(fill   ="NA", colour ="NA"), 
        panel.border     = element_rect(fill   ="NA", colour ="NA"), 
        panel.grid.major = element_line(colour ="NA"), 
        panel.grid.minor = element_line(colour ="NA")                    
  )

df = idxb
# Time-Series B/Bmsy
pB= ggplot(df$i,aes(x=yr,y=mu,group=name))+ 
  geom_ribbon(aes(ymin=lci, ymax=pmin(5,uci),fill=name),alpha=0.2, linetype=1,color=NA)+theme_bw()+
  geom_line(aes(color=name),cex=0.4)+theme(legend.position = "none")+
  geom_point(data=df$o,aes(year,obs,color=name),cex=0.6)+ylab(expression(B/B[MSY]))+
  xlab("Year")+scale_y_continuous(expand = c(0, 0))+scale_x_continuous(expand = c(0, 0))+
  geom_hline(yintercept = 1,linetype=1,col="black")+
  coord_cartesian(ylim=c(0,5))+theme( plot.margin = unit(c(3.5, 5.5, 0, 0), "points"))

dd = df$i[df$i$yr==max(df$i$yr),]
dB<-ggplot(dd) + 
  geom_histogram(aes(x = mu, y =  ..count..),fill="darkgreen",alpha=0.5,col="gray",position = "stack")+ 
  geom_vline(xintercept=1,col=1)+
  scale_y_continuous(expand = c(0, 0))+scale_x_continuous(expand = c(0, 0))+ coord_flip(xlim=c(0,5))+
  theme(legend.position = "none", 
        axis.title.x = element_text(colour ='NA'), 
        axis.text.x  = element_text(colour ="NA"), 
        axis.ticks.x = element_line(colour ="NA"),
        axis.ticks =   element_line(colour ="NA"),
        axis.title.y = element_blank(), 
        axis.text.y  = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(3.5, 5.5, 0, 0), "points"),
        panel.background = element_rect(fill   ="NA", colour ="NA"), 
        panel.border     = element_rect(fill   ="NA", colour ="NA"), 
        panel.grid.major = element_line(colour ="NA"), 
        panel.grid.minor = element_line(colour ="NA")                    
  )


df = idxt
# Time-Series B/Btrig
# Time-Series
pT= ggplot(df$i,aes(x=yr,y=mu,group=name))+ 
  geom_ribbon(aes(ymin=lci, ymax=pmin(10,uci),fill=name),alpha=0.2, linetype=1,color=NA)+theme_bw()+
  geom_line(aes(color=name),cex=0.4)+theme(legend.position = "none")+
  geom_point(data=df$o,aes(year,obs,color=name),cex=0.6)+ylab(expression(B/B[trigger]))+
  xlab("")+scale_y_continuous(expand = c(0, 0),breaks=seq(0,10,1))+scale_x_continuous(expand = c(0, 0))+
  geom_hline(yintercept = 1,linetype=1,col="black")+
  coord_cartesian(ylim=c(0,10))+theme( plot.margin = unit(c(3.5, 5.5, 0, 0), "points"))

dd = df$i[df$i$yr==max(df$i$yr),]
dT<-ggplot(dd) + 
  geom_histogram(aes(x = mu, y =  ..count..),alpha=0.5,fill="orange",col="gray",position = "stack")+ 
  geom_vline(xintercept=1,col="black")+
  scale_y_continuous(expand = c(0, 0))+scale_x_continuous(expand = c(0, 0))+ coord_flip(xlim=c(0,10))+
  theme(legend.position = "none", 
        axis.title.x = element_text(colour ='NA'), 
        axis.text.x  = element_text(colour ="NA"), 
        axis.ticks.x = element_line(colour ="NA"),
        axis.ticks =   element_line(colour ="NA"),
        axis.title.y = element_blank(), 
        axis.text.y  = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(3.5, 5.5, 0, 0), "points"),
        panel.background = element_rect(fill   ="NA", colour ="NA"), 
        panel.border     = element_rect(fill   ="NA", colour ="NA"), 
        panel.grid.major = element_line(colour ="NA"), 
        panel.grid.minor = element_line(colour ="NA")                    
  )


df = idxa
# Time-Series B/Btrig
# Time-Series
pA= ggplot(df$i,aes(x=yr,y=mu,group=name))+ 
  geom_ribbon(aes(ymin=lci, ymax=pmin(5,uci),fill=name),alpha=0.2, linetype=1,color=NA)+theme_bw()+
  geom_line(aes(color=name),cex=0.4)+theme(legend.position = "none")+
  geom_point(data=df$o,aes(year,obs,color=name),cex=0.6)+ylab(expression(ABI[MSY]))+
  xlab("Year")+scale_y_continuous(expand = c(0, 0))+scale_x_continuous(expand = c(0, 0))+
  geom_hline(yintercept = 1,linetype=1,col="black")+
  coord_cartesian(ylim=c(0,5))+theme( plot.margin = unit(c(3.5, 5.5, 0, 0), "points"))

dd = df$i[df$i$yr==max(df$i$yr),]
dA<-ggplot(dd) + 
  geom_histogram(aes(x = mu, y =  ..count..),fill="blue",alpha=0.5,col="gray",position = "stack")+ 
  geom_vline(xintercept=1,col=1)+
  scale_y_continuous(expand = c(0, 0))+scale_x_continuous(expand = c(0, 0))+ coord_flip(xlim=c(0,5))+
  theme(legend.position = "none", 
        axis.title.x = element_text(colour ='NA'), 
        axis.text.x  = element_text(colour ="NA"), 
        axis.ticks.x = element_line(colour ="NA"),
        axis.ticks =   element_line(colour ="NA"),
        axis.title.y = element_blank(), 
        axis.text.y  = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(3.5, 5.5, 0, 0), "points"),
        panel.background = element_rect(fill   ="NA", colour ="NA"), 
        panel.border     = element_rect(fill   ="NA", colour ="NA"), 
        panel.grid.major = element_line(colour ="NA"), 
        panel.grid.minor = element_line(colour ="NA")                     
  )




fnVP=function(pF,dF,pT,dT,pB,dB,pA,dA){
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(10, 8)))  # 5 by 5 grid
  print(dF , vp=vplayout(1:5,4))         # 2nd to the left +opts(legend.position = c(0,1.05)) + opts(legend.text = theme_text(colour = "black")) 
  print(pF, vp=vplayout(1:5,1:3))   
  print(dT , vp=vplayout(1:5,8))         # 2nd to the left +opts(legend.position = c(0,1.05)) + opts(legend.text = theme_text(colour = "black")) 
  print(pT, vp=vplayout(1:5,5:7)) 
  print(dB , vp=vplayout(6:10,4))         # 2nd to the left +opts(legend.position = c(0,1.05)) + opts(legend.text = theme_text(colour = "black")) 
  print(pB, vp=vplayout(6:10,1:3))   
  print(dA , vp=vplayout(6:10,8))         # 2nd to the left +opts(legend.position = c(0,1.05)) + opts(legend.text = theme_text(colour = "black")) 
  print(pA, vp=vplayout(6:10,5:7)) 
}
# Plot trends with Marginal
if(saveplot) png("../results/jara_fitsMarg.png", 2500,2100, res=300)
fnVP(pF,dF,pT,dT,pB,dB,pA,dA)
dev.off()

#><>><>><>><>><>><>><>><>><>><>><>
# Stock Status Classification
#><>><>><>><>><>><>><>><>><>><>><>

# Kobe plot
dkb = rbind(mb,mt)
dkb$x = dkb$`50%`
dkb$y = (rep(mf$`50%`,2))


dkb$Obs
dkb$run = ifelse(dkb$run=="B/Bmsy","B[MSY]","B[MSY]~`=`~B[trigger]")
dkb$run = factor(dkb$run,levels=c("B[MSY]~`=`~B[trigger]","B[MSY]"))
if(saveplot) png("../results/jara_GoldiLocks.jpg", 3000,2000, res=300)

kb <- ggplot(dkb, aes(x = x, y = y)) + 
    geom_rect(aes(xmin = 1, xmax = Inf, ymin = 0, ymax = 1), 
              colour = "green", fill = "green") + 
     geom_rect(aes(xmin = 0,xmax = 1, ymin = 0, ymax = 1), 
               colour = "yellow",  fill = "yellow") + 
     geom_rect(aes(xmin = 1, xmax = Inf,                                                    ymin = 1, ymax = Inf), colour = "orange", fill = "orange") + geom_rect(aes(xmin = 0, xmax = 1, ymin = 1, ymax = Inf), 
              colour = "red", fill = "red") + 
    ylab(expression(F/F[MSY]))+xlab(expression(B/B[MSY]))+
    geom_hline(aes(yintercept = 1)) + 
    geom_vline(aes(xintercept = 1)) + 
     scale_x_continuous(expand = c(0,0),lim=c(0,2.05)) + 
     scale_y_continuous(expand = c(0, 0),lim=c(0,2.))+
     facet_grid(~run,labeller = label_parsed)+geom_path()+theme(legend.title = element_blank())+
       geom_point(aes(colour = Year)) + scale_colour_gradient(low="lightblue",high="darkblue")+
  theme(legend.position = "right",
        legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_blank(), #change legend title font size
        legend.text = element_text(size=7),
        plot.title =  element_text(size=10))
kb
dev.off()    

d = idxf$i
d$ffmsy = d$mu
d$Year = d$yr
ices.status=rbind(
data.frame(d,bbmsy=idxb$i$mu,ref="B[MSY]"),
data.frame(d,bbmsy=idxt$i$mu,ref="B[MSY]~`=`~B[trigger]"))


# FAO status (Bmsy)
ices.status$fao.blue = ifelse(ices.status$bbmsy>1.2,1,0)
ices.status$fao.green = ifelse(ices.status$bbmsy<=1.2 & ices.status$bbmsy>0.8,1,0)
ices.status$fao.red = ifelse(ices.status$bbmsy<0.8,1,0)
# Kobe status (Bmsy)
ices.status$kobe.green = ifelse(ices.status$bbmsy>=1 & ices.status$ffmsy<1,1,0)
ices.status$kobe.yellow = ifelse(ices.status$bbmsy<1 & ices.status$ffmsy<1,1,0)
ices.status$kobe.orange = ifelse(ices.status$bbmsy>=1 & ices.status$ffmsy>1,1,0)
ices.status$kobe.red = ifelse(ices.status$bbmsy<1 & ices.status$ffmsy>1,1,0)

# Prepare Kobe time-series
df=  reshape2::melt(ices.status, id.vars = c("Year","ref"), measure.vars =  c("kobe.green", "kobe.yellow","kobe.orange","kobe.red"),na.rm = T)
df=df[df$ref=="B[MSY]",]
tdf=data.frame(year=reshape2::dcast(df, formula = Year ~ variable,sum)[,1],reshape2::dcast(df, formula = Year ~ variable,sum)[,-1]/apply(reshape2::dcast(df, formula = Year ~ variable,sum)[,-1],1,sum))
ts = reshape2::melt(tdf, id.vars = c("year")) 
ts$Status =  c("Sustainable","Rebuilding","Overfishing","Overfished")[match(ts$variable,unique(ts$variable))]
ts$Status = factor(ts$Status,levels=unique(ts$Status))
ts.msy = ts
df=  reshape2::melt(ices.status, id.vars = c("Year","ref"), measure.vars =  c("kobe.green", "kobe.yellow","kobe.orange","kobe.red"),na.rm=T)
df=df[df$ref=="B[MSY]~`=`~B[trigger]",]
tdf=data.frame(year=reshape2::dcast(df, formula = Year ~ variable,sum)[,1],reshape2::dcast(df, formula = Year ~ variable,sum)[,-1]/apply(reshape2::dcast(df, formula = Year ~ variable,sum)[,-1],1,sum))
ts = reshape2::melt(tdf, id.vars = c("year")) 
ts$Status =  c("Sustainable","Rebuilding","Overfishing","Overfished")[match(ts$variable,unique(ts$variable))]
ts$Status = factor(ts$Status,levels=unique(ts$Status))
ts.tri = ts

kbts = rbind(data.frame(ts.msy,ref="B[MSY]"),
           data.frame(ts.tri,ref="B[MSY]~`=`~B[trigger]"))
kbts$ref = factor(kbts$ref,levels=c("B[MSY]~`=`~B[trigger]","B[MSY]"))

# Plot KB time series
pkbts = ggplot(kbts, aes(x=year, y=value, fill=Status)) + 
  facet_grid(~ref,labeller = label_parsed)+
  geom_area()+ylab("Proportions of Stocks")+xlab("Year")+ggtitle("(b)")+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1))+ 
  scale_x_continuous(expand = c(0,0),breaks=seq(1975,2020,5))+
  scale_fill_manual(values =c("green","yellow","orange","red"))+#geom_line(aes(year,catch),col=1,size=0.8)+
  theme(legend.position = "right",
        legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_blank(), #change legend title font size
        legend.text = element_text(size=7),
        plot.title =  element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5,size=10))+ #change legend text font size
        geom_hline(yintercept = c(0.25,0.5,0.75),size=0.3,col="grey",linetype=2)


# Prepare FAO time series
df=  reshape2::melt(ices.status, id.vars = c("Year","ref"), measure.vars = c("fao.blue", "fao.green","fao.red"),na.rm=T)
df=df[df$ref=="B[MSY]",]
tdf=data.frame(year=reshape2::dcast(df, formula = Year ~ variable,sum)[,1],reshape2::dcast(df, formula = Year ~ variable,sum)[,-1]/apply(reshape2::dcast(df, formula = Year ~ variable,sum)[,-1],1,sum))
ts = reshape2::melt(tdf, id.vars = c("year")) 
ts$Status =  c("Underfished","Sustainable","Overfished")[match(ts$variable,unique(ts$variable))]
ts$Status = factor(ts$Status,levels=unique(ts$Status))
ts.msy = ts
df=  reshape2::melt(ices.status, id.vars = c("Year","ref"), measure.vars =  c("fao.blue", "fao.green","fao.red"),na.rm=T)
df=df[df$ref=="B[MSY]~`=`~B[trigger]",]
tdf=data.frame(year=reshape2::dcast(df, formula = Year ~ variable,sum)[,1],reshape2::dcast(df, formula = Year ~ variable,sum)[,-1]/apply(reshape2::dcast(df, formula = Year ~ variable,sum)[,-1],1,sum))
ts = reshape2::melt(tdf, id.vars = c("year")) 
ts$Status =  c("Underfished","Sustainable","Overfished")[match(ts$variable,unique(ts$variable))]
ts$Status = factor(ts$Status,levels=unique(ts$Status))
ts.tri = ts

faots = rbind(data.frame(ts.msy,ref="B[MSY]"),
             data.frame(ts.tri,ref="B[MSY]~`=`~B[trigger]"))
faots$ref = factor(faots$ref,levels=c("B[MSY]~`=`~B[trigger]","B[MSY]"))

pfao = ggplot(faots, aes(x=year, y=value, fill=Status)) + 
  facet_grid(~ref,labeller = label_parsed)+
  geom_area()+ylab("Proportions of Stocks")+xlab("Year")+ggtitle("(a)")+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1))+ 
  scale_x_continuous(expand = c(0,0),breaks=seq(1975,2020,5))+
  scale_fill_manual(values =c("blue","green","red"))+
  theme(legend.position = "right",
        legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_blank(), #change legend title font size
        legend.text = element_text(size=7),
        plot.title =  element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5,size=10))+ #change legend text font size
  geom_hline(yintercept = c(0.25,0.5,0.75),size=0.3,col="grey",linetype=2)


if(saveplot) png("../results/jara_ices_fao.jpg", 2500,2500, res=300)
grid.arrange(pfao,pkbts)
dev.off()


kbc= kb+theme(legend.position = "right",
      legend.key.size = unit(1, 'cm'), #change legend key size
      legend.key.height = unit(0.5, 'cm'), #change legend key height
      legend.key.width = unit(0.5, 'cm'), #change legend key width
      legend.title = element_blank(), #change legend title font size
      legend.text = element_text(size=7),
      plot.title =  element_text(size=10),
      legend.box.spacing = unit(1, 'cm'))+ggtitle("(c) Goldilocks")


if(saveplot) png("../results/jara_ices_fao_goldi.jpg", 2500,3400, res=300)
ggarrange(pfao,pkbts,kb+ggtitle("(c)"),  align="v",ncol=1)
dev.off()
