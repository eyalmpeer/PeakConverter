library(ggplot2)
args<-commandArgs()
if (length(args)!=3) { 
	stop("One argument must be supplied (input file)", call.=FALSE)
}

peaks<-read.table(args[3])
colnames(peaks)<-c('Tx_ID', 'Peak_Start', 'Peak_End', 'Peak_Names', 'Peak_Middle', 'ATG', 'Stop', 'Length', 'FirstSS', 'LastSS')

#For Metagene and ATG window analysis use the subset of coding transcripts:
x<- subset(peaks,ATG>-1,select=c(Tx_ID,Peak_Middle,ATG,Stop,Length))

#Add columns of segment lengths:
x['5UTR']<-x['ATG']
x['CDS']<-x['Stop']-x['ATG']
x['3UTR']<-x['Length']-x['Stop']

#Get the mean length fractions:
m5len<-(mean(x$'5UTR')/mean(x$Length))
mcdslen<-(mean(x$CDS)/mean(x$Length))
m3len<- (mean(x$'3UTR')/mean(x$Length))

#Add a column of normalized peak middle location in transcript:
x$NormLoc <- ifelse(x$Peak_Middle<x$ATG, (x$Peak_Middle/x$'5UTR')*m5len, 
ifelse(x$Peak_Middle<x$Stop, ((x$Peak_Middle-x$ATG)/x$CDS)*mcdslen+m5len, 
((x$Peak_Middle-x$Stop)/x$'3UTR')*m3len+mcdslen+m5len))

#Create PDF file of metagene plot:
pdf('m1A_peaks_metagene.pdf',width=7,height=5)
ggplot(x, aes(x=NormLoc, fill="m1A Peaks")) + 
geom_density(alpha = 0.1)+ theme(panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
labs(title = "Distribution of m1A peaks in transcripts", x = NULL, y = "Peak Frequency") + 
geom_vline(xintercept = c(m5len,mcdslen+m5len))+scale_x_continuous(breaks=c(m5len,mcdslen+m5len), labels=c("ATG","Stop")) + 
scale_fill_manual(values = c("blue"))+guides(fill=FALSE)
dev.off()

#Calculate peak middle distance from ATG:
x$ATGDist<-x$Peak_Middle-x$ATG
#Take peaks that lie inside 600 nt window around the ATG:
atgwindow<-subset(x,x$ATGDist<=300 & x$ATGDist>=-300)

#Create PDF file of AUG window plot
pdf('m1A_peaks_AUG.pdf',width=4,height=4)
ggplot(atgwindow, aes(x=ATGDist, fill="m1A Peaks")) + 
geom_density(alpha = 0.1)+ theme(panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
labs(title = "m1A peaks around canonical AUG", x = NULL, y = "Peak Frequency") + scale_fill_manual(values = c("blue"))+guides(fill=FALSE)
dev.off()

#Create subset of spliced transcripts:
splice<-subset(peaks,peaks$FirstSS>-1)[c('Tx_ID', 'Peak_Middle', 'FirstSS')]
#Calculate peak middle distance from the first splice site:
splice$FSSDist<-splice$Peak_Middle-splice$FirstSS
#Take peaks that lie inside 600 nt window around the first splice site:
splicewindow<-subset(splice,splice$FSSDist<=300 & splice$FSSDist>=-300)
#Create PDF file of first splice site window plot:
pdf('m1A_peaks_splicesite.pdf',width=4,height=4)
ggplot(splicewindow, aes(x=FSSDist, fill="m1A Peaks")) + 
geom_density(alpha = 0.1)+ theme(panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
labs(title = "m1A peaks around 1st splice site", x = NULL, y = "Peak Frequency") + scale_fill_manual(values = c("blue"))+guides(fill=FALSE)
dev.off()
