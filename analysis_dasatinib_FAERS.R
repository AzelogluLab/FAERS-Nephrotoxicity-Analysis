######################################################################
# Script for analysis of FAERS dataset for Dasatinib renal toxicity  #
######################################################################

# Import FAERS frequency counts from the AERSMINE resource
dat<-read.table("data/aersMineExploreDataSet_4372 additional renal ADRS.tsv",skip=5,header=TRUE,sep="\t")

# Remove irrelevant columns
dat<-dat[,-grep("FDR|Fishers",names(dat)),]
dat<-dat[,13:24]

# Define renal ADRs of interest and assign to dataset
ADRs<-c("glomnephritisNS",
        "nephropathies",
        "hypertension",
        "RenalDisorderExclNephr",
        "Tdisorder",
        "Tnecrose",
        "Tacidose",
        "Tatrophy",
        "Tdysfun",
        "Tinjury"                       )
names(dat)[3:12]<-ADRs

# Clean up dataset
remComma<-function(x){return(as.numeric(gsub(pattern = "\\,","",x)))}

dat[,3]<-remComma(dat[,3])
dat[,4]<-remComma(dat[,4])
dat[,5]<-remComma(dat[,5])
dat[,6]<-remComma(dat[,6])
dat[,7]<-remComma(dat[,7])
dat[,8]<-remComma(dat[,8])
dat[,9]<-remComma(dat[,9])
dat[,10]<-remComma(dat[,10])
dat[,11]<-remComma(dat[,11])
dat[,12]<-remComma(dat[,12])
dat$Total_Drugs_Reports<-remComma(dat$Total_Drugs_Reports)

# Select kinase inhibitors
drugsOI<-as.character(dat$Drugs[grep("nib",dat$Drugs)])
dat<-dat[dat$Drugs %in% c(drugsOI,"Unique Patients"),]

# Compute RORs.

library(dplyr)

res<-data.frame()
for(tox in ADRs){
  
  for(d in drugsOI){
    
    dycy = dat[dat$Drugs==d,tox]  
    dycn = dat$Total_Drugs_Reports[dat$Drugs==d]   - dycy
    
    if(tox == "glomnephritisNS")        dncy <- sum(dat$glomnephritisNS[dat$Drugs != d & dat$Drugs != "Unique Patients"],na.rm=T) 
    if(tox == "nephropathies")          dncy <- sum(dat$nephropathies[dat$Drugs != d & dat$Drugs != "Unique Patients"],na.rm=T)  
    if(tox == "hypertension")           dncy <- sum(dat$hypertension[dat$Drugs != d & dat$Drugs != "Unique Patients"],na.rm=T) 
    if(tox == "RenalDisorderExclNephr") dncy <- sum(dat$RenalDisorderExclNephr[dat$Drugs != d & dat$Drugs != "Unique Patients"],na.rm=T) 
    if(tox == "Tdisorder")              dncy <- sum(dat$Tdisorder[dat$Drugs != d & dat$Drugs != "Unique Patients"],na.rm=T) 
    if(tox == "Tnecrose")               dncy <- sum(dat$Tnecrose[dat$Drugs != d & dat$Drugs != "Unique Patients"],na.rm=T) 
    if(tox == "Tacidose")               dncy <- sum(dat$Tacidose[dat$Drugs != d & dat$Drugs != "Unique Patients"],na.rm=T) 
    if(tox == "Tatrophy")               dncy <- sum(dat$Tatrophy[dat$Drugs != d & dat$Drugs != "Unique Patients"],na.rm=T) 
    if(tox == "Tdysfun")                dncy <- sum(dat$Tdysfun[dat$Drugs != d & dat$Drugs != "Unique Patients"],na.rm=T) 
    if(tox == "Tinjury")                dncy <- sum(dat$Tinjury[dat$Drugs != d & dat$Drugs != "Unique Patients"],na.rm=T) 
    
    dncn = sum(dat$Total_Drugs_Reports[dat$Drugs != d],na.rm=T) - dncy
    
    orr=(dycy / dycn)/ (dncy / dncn)
    ln_orr<-log(orr)
    SE_ln_orr <- sqrt((1/dycy) + (1/dycn) + (1/dncy) + (1/dncn))
    lnCIup<-ln_orr+1.96 * SE_ln_orr
    lnCIdown<-ln_orr-1.96 * SE_ln_orr
    CIup<-exp(lnCIup)
    CIdown<-exp(lnCIdown)
    
    res<-rbind(res,data.frame(tox=tox,drug=d,ror=orr,CIup=CIup,CIdown=CIdown,dycy=dycy,dycn=dycn))
  }
}

# Define specific KIs of interest
drugsOI<-c("axitinib" ,    "crizotinib" , "dasatinib",   "erlotinib"  ,
           "gefitinib"  ,"imatinib"  ,  "ruxolitinib" ,"lapatinib"  , "nilotinib" ,
           "pazopanib" ,  "sorafenib"  , "sunitinib" ,  "vandetanib"  ,"tofacitinib",
           "bosutinib" )
res$drug<-as.character(res$drug)

# Select from ROR analysis only KIs of interest
res<-res[res$drug %in% drugsOI,]
res$ror[is.infinite(res$ror)]<-NA


# Generate figures
library(ggplot2)

# Figure all tubular
pdf("figures_paper_submission2/figure_Supplemental_all_Tubular_ADRs.pdf", width=10,height=8)
ggplot(data=res[res$tox %in% c( "Tdisorder",
                                "Tnecrose",
                                "Tacidose",
                                "Tatrophy",
                                "Tdysfun",
                                "Tinjury"     ),])+
  geom_point(aes(x=drug,y=ror))+
  geom_errorbar(aes(x=drug,ymin=CIdown,ymax=CIup))+
  ggtitle("All tubular ADRs")+
  ylab("Reporting odds ratio")+xlab("")+
  coord_flip()+facet_wrap(~tox)+
  theme(legend.position="none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))
dev.off()

postscript("figures_paper_submission2/figure_Supplemental_all_Tubular_ADRs.eps", width=10,height=8)
ggplot(data=res[res$tox %in% c( "Tdisorder",
                                "Tnecrose",
                                "Tacidose",
                                "Tatrophy",
                                "Tdysfun",
                                "Tinjury"     ),])+
  geom_point(aes(x=drug,y=ror))+
  geom_errorbar(aes(x=drug,ymin=CIdown,ymax=CIup))+
  ggtitle("All tubular ADRs")+
  ylab("Reporting odds ratio")+xlab("")+
  coord_flip()+facet_wrap(~tox)+
  theme(legend.position="none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))
dev.off()


# Nephropathies figure
resNP<-res[res$tox=="nephropathies" ,]
resNP<-resNP[order(resNP$ror),]
resNP$drug<-factor(resNP$drug,levels = unique(resNP$drug))


colors<-data.frame(drug=c("dasatinib","imatinib","nilotinib","vandetanib","erlotinib","bosutinib"),
                   col=c(rgb(maxColorValue = 1, red=.6,green=.1,blue=.0),
                         rgb(maxColorValue = 1, red=.8,green=.4,blue=.7),
                         rgb(maxColorValue = 1, red=.2,green=.5,blue=.4),
                         rgb(maxColorValue = 1, red=.4,green=.5,blue=.2),
                         rgb(maxColorValue = 1, red=.9,green=.6,blue=.3),
                         rgb(maxColorValue = 1, red=.2,green=.4,blue=.6)))
resNP<-merge(resNP,colors,by="drug",all=T)
resNP$col<-as.character(resNP$col)
resNP$col[is.na(resNP$col)]<-rgb(maxColorValue = 1, red=0,green=0,blue=0)




pdf("figures_paper_submission2/figure_1_Nephropathies.pdf")
ggplot(data=resNP)+
  geom_point(aes(x=drug,y=ror),col="black")+
  theme_classic()+
  geom_errorbar(aes(x=drug,ymin=CIdown,ymax=CIup),col="black")+
  ggtitle("Nephropathies")+
  ylab("Reporting odds ratio")+xlab("")+
  coord_flip()+
  scale_color_gradient(low="blue", high="red")+
  theme(legend.position="none")+
  theme(legend.position="none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))
dev.off()

postscript("figures_paper_submission2/figure_1_Nephropathies.eps")
ggplot(data=resNP)+
  geom_point(aes(x=drug,y=ror,col=ror))+
  theme_classic()+
  geom_errorbar(aes(x=drug,ymin=CIdown,ymax=CIup,col=ror))+
  ggtitle("Nephropathies")+
  ylab("Reporting odds ratio")+xlab("")+
  coord_flip()+
  scale_color_gradient(low="blue", high="red")+
  theme(legend.position="none")+
  theme(legend.position="none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))
dev.off()






resTUB<-res[res$tox %in% c("Tdisorder"),]
resTUB<-resTUB[order(resTUB$ror),]
resTUB<-resTUB[c(15,1:14),]
resTUB$drug<-factor(resTUB$drug,levels = unique(resTUB$drug))


pdf("figures_paper_submission2/figure_Supplemental_TubularDisorder.pdf")
ggplot(data=resTUB)+
  geom_point(aes(x=drug,y=ror))+
  theme_classic()+
  geom_errorbar(aes(x=drug,ymin=CIdown,ymax=CIup))+
  ggtitle("Tubular ADRs")+
  ylab("Reporting odds ratio")+xlab("")+
  coord_flip()+
  theme(legend.position="none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))
dev.off()

postscript("figures_paper_submission2/figure_Supplemental_TubularDisorder.eps")
ggplot(data=resTUB)+
  geom_point(aes(x=drug,y=ror))+
  theme_classic()+
  geom_errorbar(aes(x=drug,ymin=CIdown,ymax=CIup))+
  ggtitle("Tubular ADRs")+
  ylab("Reporting odds ratio")+xlab("")+
  coord_flip()+
  theme(legend.position="none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))
dev.off()


# Figure GN
resGN<-res[res$tox=="glomnephritisNS" ,]
resGN<-resGN[order(resGN$ror),]
resGN$drug<-factor(resGN$drug,levels = unique(resGN$drug))

pdf("figures_paper_submission2/figure_Supplemental_Glomerulonephritis_Nephrotic syndrome.pdf")
ggplot(data=resGN)+
  geom_point(aes(x=drug,y=ror,col=ror))+
  theme_classic()+
  geom_errorbar(aes(x=drug,ymin=CIdown,ymax=CIup,col=ror))+
  ggtitle("Glomerulonephritis & Nephrotic syndrome ")+
  ylab("Reporting odds ratio")+xlab("")+
  scale_color_gradient(low="blue", high="red")+
  coord_flip()+
  theme(legend.position="none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))
dev.off()


postscript("figures_paper_submission2/figure_Supplemental_Glomerulonephritis_Nephrotic syndrome.eps")
ggplot(data=resGN)+
  geom_point(aes(x=drug,y=ror,col=ror))+
  theme_classic()+
  geom_errorbar(aes(x=drug,ymin=CIdown,ymax=CIup,col=ror))+
  ggtitle("Glomerulonephritis & Nephrotic syndrome ")+
  ylab("Reporting odds ratio")+xlab("")+
  scale_color_gradient(low="blue", high="red")+
  coord_flip()+
  theme(legend.position="none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))
dev.off()


### Figure GN vs HT
resGN<-res[res$tox%in%c("glomnephritisNS") ,-c(1,6,7)]
names(resGN)[2:4]<-c("ROR_GN","CIup_GN","CIdown_GN")

resHT<-res[res$tox%in%c("hypertension") ,-c(1,6,7)]
names(resHT)[2:4]<-c("ROR_HT","CIup_HT","CIdown_HT")

resGNHT<-cbind(resHT,resGN)
resGNHT[,5]<-NULL

library(ggrepel)



pdf("figures_paper_submission2/figure_1_GlomNeph_vs_Hypertension_color.pdf")
ggplot(resGNHT,aes(x=ROR_HT,y=ROR_GN))+
  geom_errorbar(aes(x=ROR_HT,ymin=CIdown_GN,ymax=CIup_GN),col="grey",alpha=.5)+
  geom_errorbarh(aes(y=ROR_GN,xmin=CIdown_HT,xmax=CIup_HT),col="grey",alpha=.5)+
  geom_point(aes(col=drug))+
  geom_text_repel(aes(label=drug,col=drug),nudge_y=.05,nudge_x=.05,size=6)+
  theme_classic()+
  #coord_cartesian(ylim=c(0,1))+
  xlab("Hypertension ROR")+ylab("Glomerulonephritis & Nephrotic syndrome ROR")+
  theme(legend.position="none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))
dev.off()

colors<-data.frame(drug=c("dasatinib","imatinib","nilotinib","vandetanib","erlotinib","bosutinib"),
                   col=c(rgb(maxColorValue = 1, red=.6,green=.1,blue=.0),
                         rgb(maxColorValue = 1, red=.8,green=.4,blue=.7),
                         rgb(maxColorValue = 1, red=.2,green=.5,blue=.4),
                         rgb(maxColorValue = 1, red=.4,green=.5,blue=.2),
                         rgb(maxColorValue = 1, red=.9,green=.6,blue=.3),
                         rgb(maxColorValue = 1, red=.2,green=.4,blue=.6)))
resGNHT2<-merge(resGNHT,colors,by="drug",all=T)
resGNHT2$col<-as.character(resGNHT2$col)
resGNHT2$col[is.na(resGNHT2$col)]<-rgb(maxColorValue = 1, red=0,green=0,blue=0)

pdf("figures_paper_submission2/figure_1_GlomNeph_vs_Hypertension_col.pdf")
ggplot(resGNHT2,aes(x=ROR_HT,y=ROR_GN))+
  geom_errorbar(aes(x=ROR_HT,ymin=CIdown_GN,ymax=CIup_GN),col="grey",alpha=.5)+
  geom_errorbarh(aes(y=ROR_GN,xmin=CIdown_HT,xmax=CIup_HT),col="grey",alpha=.5)+
  geom_point(aes(col=drug))+
  geom_text_repel(aes(label=drug),nudge_y=.05,nudge_x=.05,size=6,col=resGNHT2$col)+
  theme_classic()+
  scale_color_manual(values=resGNHT2$col) +
  #coord_cartesian(ylim=c(0,1))+
  xlab("Hypertension ROR")+ylab("Glomerulonephritis & Nephrotic syndrome ROR")+
  theme(legend.position="none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))
dev.off()

postscript("figures_paper_submission2/figure_1_GlomNeph_vs_Hypertension_col.eps")
ggplot(resGNHT2,aes(x=ROR_HT,y=ROR_GN))+
  geom_errorbar(aes(x=ROR_HT,ymin=CIdown_GN,ymax=CIup_GN),col="grey",alpha=.5)+
  geom_errorbarh(aes(y=ROR_GN,xmin=CIdown_HT,xmax=CIup_HT),col="grey",alpha=.5)+
  geom_point(aes(col=drug))+
  geom_text_repel(aes(label=drug),nudge_y=.05,nudge_x=.05,size=6,col=resGNHT2$col)+
  theme_classic()+
  scale_color_manual(values=resGNHT2$col) +
  #coord_cartesian(ylim=c(0,1))+
  xlab("Hypertension ROR")+ylab("Glomerulonephritis & Nephrotic syndrome ROR")+
  theme(legend.position="none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))
dev.off()





# Tubular vs NP
resTD<-res[res$tox=="Tdisorder",]

resNP<-res[res$tox=="nephropathies" ,]
resNP<-resNP[order(resNP$ror),]
resNP$drug<-factor(resNP$drug,levels = unique(resNP$drug))

resTDNP<-merge(resTD,resNP,by="drug")
names(resTDNP)[3]<-"rorTubDis"
names(resTDNP)[4]<-"ciupTubDis"
names(resTDNP)[5]<-"cidownTubDis"

names(resTDNP)[9]<-"rorNP"
names(resTDNP)[10]<-"ciupNP"
names(resTDNP)[11]<-"cidownNP"

ggplot(resTDNP,aes(x=rorTubDis,y=rorNP))+
  geom_errorbar(aes(x=rorTubDis,ymin=cidownNP,ymax=ciupNP),col="grey",alpha=.5)+
  geom_errorbarh(aes(y=rorNP,xmin=cidownTubDis,xmax=ciupTubDis),col="grey",alpha=.5)+
  geom_point(aes(col=drug))+
  geom_text_repel(aes(label=drug,col=drug),nudge_y=.05,nudge_x=.05,size=4)+
  theme_classic()+
  coord_cartesian(xlim=c(0,4))+
  xlab("Nephropathies ROR")+ylab("TubDis & Nephropathies ROR")+
  theme(legend.position="none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))

# Write numbers to tables
write.csv(x=resGN,file="figures_paper_submission2/table_GlomNep.csv",quote=F,row.names = F)
write.csv(x=resNP,file="figures_paper_submission2/table_Nephropath.csv",quote=F,row.names = F)
write.csv(x=resTD,file="figures_paper_submission2/table_TubularDysfun.csv",quote=F,row.names = F)




