#Code for calculation of ribosomal distances relative to either the translation start site or the stop codon. 

library(dplyr)
library(stringr)
library(ggplot2)

#Part A : preparation of startcodon file from a gtf file containing coordinates of start codons.

tss <- read.delim("startsites.gtf", header=FALSE)
tss <- select(tss, V4, V5, V9)
a <- str_split_fixed(tss$V9, ";", 3)
a <- a[,1]
a <- str_split_fixed(a, " ", 3)
a <- a[,2]
tss <- select(tss, V4, V5)
tss <- cbind(tss, a)
tss <- mutate(tss, startplusone=round((V4+V5)/2))
tss <- select(tss, a, startplusone)

#Part B: Preparation of stopcodon file from a gtf file containing coordinates of stop codons.

sc <- read.delim("Galaxy378-[Filter_on_data_353].gtf", header=FALSE)
sc <- select(sc, V4, V5, V9)
a <- str_split_fixed(sc$V9, ";", 3)
a <- a[,2]
a <- str_split_fixed(a, " ", 3)
a <- a[,3]
sc <- select(sc, V4, V5)
sc <- cbind(sc, a)
sc <- mutate(sc, stopminusone=round((V4+V5)/2))
sc <- select(sc, a, stopminusone)

#Part C: Read in a bed file of the aligned reads with an interval corresponding to mDNA coding sequences. 
#Create separate dataframes with reads on positive or negative strand. Extract the name of the 
#corresponding mRNA and calculate the middle position of the interval. 
inter <- read.delim("D2.bed", header=FALSE)
pos <- inter[inter$V7=="+",]
neg <- inter[inter$V7=="-",]

pos <- select(pos, V1,V4,V5,V7,V9)
pos <- mutate(pos, base=round((V4+V5)/2))

a <- str_split_fixed(pos$V9, ";", 5)
a <- a[,2]
a <- str_split_fixed(a, " ", 2)
a <- a[,2]
a <- str_split_fixed(a, " ", 2)
a <- a[,2]

pos <- cbind(pos, a)

neg <- select(neg, V1,V4,V5,V7,V9)
neg <- mutate(neg, base=round((V4+V5)/2))

a <- str_split_fixed(neg$V9, ";", 5)
a <- a[,2]
a <- str_split_fixed(a, " ", 2)
a <- a[,2]
a <- str_split_fixed(a, " ", 2)
a <- a[,2]

neg <- cbind(neg, a)

#Part D: Add information about position of translation start site and calculate distances of ribosome from translation
#start site. 
#Combine distances of negative and positive strand-based calculations.

pos<- inner_join(pos, tss)
pos <- mutate(pos, dist=base-startplusone+1)

neg <- inner_join(neg, tss)
neg <- mutate(neg, dist=startplusone+1-base)
D2<- rbind(pos, neg)
write.csv(D2, "D2DistanceFromTSS.csv")

condition1<- c(pos$dist, neg$dist)

#Part E: Part E is an alternative to part D. Add information about position of stop codon and calculate distances of ribosome
#from stop codon. 
#Combine distances of negative and positive strand-based calculations.

pos <- inner_join(pos, sc)
pos <- mutate(pos, dist=stopminusone-base+1)

neg <- inner_join(neg, sc)
neg <- mutate(neg, dist=base-stopminusone+1)

D2<- rbind(pos, neg)
write.csv(D2, "D2DistanceFromStopCodon.csv")

condition1<- c(pos$dist, neg$dist)

#Part C combined with either part D or E should be repeated for other bed files corresponding to separate experiments of 
#different conditions. 

#Part D: Plotting a coloured single density plot of ribosomal positions for all conditions. 

a <- rep("cond1", length(condition1))
b <- rep("cond2", length(condition2))
c <- rep("cond3", length(condition3))
d <- rep("cond4", length(condition4))
fac <- c(a,b,c, d)
all <- c(condition1, condition2, condition3, condition4)
data <-data.frame(dens=all, lines=fac)

ggplot(data, aes(x=dens, colour=lines)) +geom_density() +coord_cartesian(xlim = c(-300, 3700)) +xlab("Distance of ribosome from start/stop codon")

#Alternatively, a black and white density plot can be plot.  

ggplot(data, aes(x=dens, linetype=lines)) +geom_density() +coord_cartesian(xlim = c(-300, 3700)) +xlab("Distance of ribosome from start/stop codon")

#Statistical analysis: I excluded mRNAs of intron-containing genes from the analysis. 329 mRNAs out of 6692 were excluded. 
#(6363 mRNAs were included in cds2.) 
#Minimal length of mRNA=47. In cds3, all transcripts shorter than 300 from cds2 were excluded.
cds <- read.csv("cds.csv")
cds2 <- cds[duplicated(cds$a),]
list <- cds2$a
cds2 <- cds[!cds$a %in% list, ]
cds2 <- mutate(cds2, dist=abs(V4-V5))
cds3 <- cds2[cds2$dist < 300, ]

inter <- read.delim("libD2high.bed", header=FALSE)
pos <- inter[inter$V7=="+",]
neg <- inter[inter$V7=="-",]
 
pos <- select(pos, V4,V5,V9)
pos <- mutate(pos, base=round((V4+V5)/2))
a <- str_split_fixed(pos$V9, ";", 5)
a <- a[,2]
a <- str_split_fixed(a, " ", 2)
a <- a[,2]
a <- str_split_fixed(a, " ", 2)
a <- a[,2]
pos <- select(pos, base)
pos <- cbind(pos, a)
neg <- select(neg, V4,V5,V9)
neg <- mutate(neg, base=round((V4+V5)/2))
a <- str_split_fixed(neg$V9, ";", 5)
a <- a[,2]
a <- str_split_fixed(a, " ", 2)
a <- a[,2]
a <- str_split_fixed(a, " ", 2)
a <- a[,2]
neg <- select(neg, base)
neg <- cbind(neg, a)
neg2<- neg[neg$a %in% cds2$a,]
pos2<- pos[pos$a %in% cds2$a,]
neg3<- neg[neg$a %in% cds3$a,]
pos3<- pos[pos$a %in% cds3$a,]
pos<- inner_join(pos, tss)

pos <- mutate(pos, dist=base-startplusone+1)
neg <- inner_join(neg, tss)
neg <- mutate(neg, dist=startplusone+1-base)
pos2<- inner_join(pos2, tss)
pos2 <- mutate(pos2, dist=base-startplusone+1)
neg2 <- inner_join(neg2, tss)
neg2 <- mutate(neg2, dist=startplusone+1-base)
pos3<- inner_join(pos3, tss)
pos3 <- mutate(pos3, dist=base-startplusone+1)
neg3 <- inner_join(neg3, tss)
neg3 <- mutate(neg3, dist=startplusone+1-base)
libD2all <-c(pos$dist, neg$dist)
libD2nointron <- c(pos2$dist, neg2$dist)
libD2long <- c(pos3$dist, neg3$dist)
save(libD2all, libD2nointron, libD2long, file="libD2test")
