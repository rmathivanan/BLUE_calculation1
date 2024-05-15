# BLUE_calculation1

phe <- read.table("PHENOTYPE.csv", head = TRUE, sep = ",")

blues=data.frame(unique(as.factor(phe$taxa)))
colnames(blues)[1]='taxa'

traits=colnames(phe)[4:29]

##BLUEs calculation

count=0
for (i in 1:length(traits)){
  tempdf=phe[,c('taxa','Batch','Run.Order',traits[i])]
  colnames(tempdf)[4]="ResponseVariable"
  NAdf=tempdf[is.na(tempdf$ResponseVariable),]
  NAtaxa=data.frame(unique(NAdf$taxa))
  print(nrow(NAtaxa))
  tempdf=na.omit(tempdf)
  if (nrow(NAtaxa!=0)){
    NAtaxa$effect=NA}
  colnames(NAtaxa)[1]='taxa'
  bluemodel=lmer(ResponseVariable~taxa+(1|Batch)+(1|Run.Order), data=tempdf)
  bluestemp=data.frame(fixef(bluemodel))
  bluestemp=data.frame(rownames(bluestemp),bluestemp)
  colnames(bluestemp)=c('taxa', 'effect')
  bluestemp$taxa=gsub('taxa','',bluestemp$taxa)
  rownames(bluestemp)=NULL
  intercept=bluestemp[1,2]
  print(nrow(bluestemp))
  NAtaxa=NAtaxa %>%
    filter(!taxa %in% unique(bluestemp$taxa))
  if (nrow(NAtaxa!=0)){
    bluestemp=rbind(bluestemp,NAtaxa)
  }
  print(nrow(bluestemp))
  name=setdiff(blues$taxa, bluestemp$taxa)
  bluestemp[1,1]=name
  bluestemp[1,2]=0
  bluestemp$effect=bluestemp$effect+intercept
  colnames(bluestemp)[2]=traits[i]
  blues <- merge(blues, bluestemp, by = 'taxa', all.x = TRUE)
  count=count+1
  print(count)
}

maizegenotype=blues$genotype
write.csv(blues, "Maize_blues.csv", row.names = F)
write.table(maizegenotype, "maizegenotype.csv", row.names = F, col.names = FALSE)






library(ggplot2)

data<-read.csv("PHENOTYPE.csv")

ggplot(data, aes(x=Raffinose))+geom_histogram(bins = 30)
