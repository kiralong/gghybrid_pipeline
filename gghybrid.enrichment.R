##############################################################################
###########gghybrid example code##############################################
##############################################################################

#Original code written by Richard Ian Bailey 17 April 2018#
#Adapted code viewed here written by Kira Long and Angel Rivera-Colon

#Install the package from GitHub#
#install.packages("devtools")
#library(devtools)
#devtools::install_github("ribailey/gghybrid")

#Attach it#
library(gghybrid)
library(data.table)

#Take a look at the function help files#
#?read.data	#Read in a data file in structure format or similar#
#?data.prep	#Prepare the data for analysis#
#?esth		#Run hybrid index estimation#
#?plot_h	#Plot estimated hybrid indices and credible intervals#
#?ggcline	#Run genomic cline estimation#
#?plot_clinecurve	#Plot one or more fitted cline curves, with or without individual data points#
#?compare.models	#compare two models run on the same data set using the widely applicable information criterion#

#(Note: gghybrid relies on the data.table package for data manipulation)#


#Example: Using a data file in structure format but with a complete header row.
#The file contains 1 column per marker (equivalent to ONEROW=0 in structure) for a set
#of diploid markers, plus haploid (mitochondrial) ND2. For ND2, the second (non-existent)
#allele copy is coded as missing data. Other file formats would require more information
#when running 'read.data' (see ?read.data).

#The data file contains two mandatory columns to the left of the marker columns, named
#INDLABEL (the individual references) and POPID (the sample population). These columns
#don't need to be named in the header row (or could have different names), but it makes it easier if they are.

#I downloaded the data (and then created a single input file in the right format)
#from Dryad here: https://www.datadryad.org/resource/doi:10.5061/dryad.v6f4d

#Set working directory
#Kira's labtop
#setwd("/Users/kira/Google Drive/Bioinformatics/Candidate_Genes/GGhybrid/Enrichment_analysis/")
#UIUC Campus Cluster
setwd("/projects/aces/kira/stacks_ref_map/ref_map_Mar30_2020/cline_analysis/gghybrid/2020_October_15_enrichment")

#Read in the data file (This is the simplest file format for reading in)#

dat=read.data("populations.p6.r50.sSNP.enrichment.structure.csv",
              nprecol=2,MISSINGVAL=NA)

#Take a look at the help file to get a handle on what the object 'dat' now contains.

###

#Data preparation and filtering. Here I'm filtering out loci that have a minor allele
#frequency greater than 0.1 in both parental reference sets. There are also options for
#filtering by difference in parental allele frequencies, and for number of allele copies
#in each parental reference set (variable among loci due to missing data).

#The function uses objects produced by 'read.data'#

prepdata=data.prep(data=dat$data,
                   loci=dat$loci,
                   alleles=dat$alleles,
                   S0=c("02SS"), #POPID names for the first parental reference set#
                   S1="10CG", #POPID names for the second parental reference set#
                   precols=dat$precols,
                   max.S.MAF = 0.05,	#Filtering by parental minor allele frequency# @KML: Original cutoff was at 10%, changed it to 5%
                   return.genotype.table=T,
                   return.locus.table=T)

#'return.genotype.table=T' makes an optional table of genotypes, where for each locus
#an individual's genotype (assuming diploidy) will be 0 (two copies of the allele with
#relatively higher frequency in the 'S0' parental set), 1 (heterozygote), 2 (two copies
#of the designated 'S1' allele). This table isn't needed in downstream functions, but
#could be useful e.g. for estimating parental linkage disequilibria (associations of
#alleles from the same parent species).

#'return.locus.table=T' is also optional and not needed downstream. It's just a table
#with one row per marker, giving some information on parental allele frequencies, sample
#sizes etc.

###
#Create a pdf for the outputs
pdf('./gghybrid_graphs_r50.pdf',8,8) #To make a single pdf with all graphs together

#Next, run hybrid index estimation#

#This function uses objects produced by both the previous functions#

hindlabel= esth(data.prep.object = prepdata$data.prep,
                read.data.precols = dat$precols,
                include.Source = TRUE,	#Set to TRUE if you want hybrid indices for the parental reference individuals#
                plot.ind = c("pr_2017_095","ru_2017_001","qp_2017_023","ro_2017_036",
                             "pr_2017_122","ru_2017_014"),
                plot.col = c("blue","green","cyan","purple","magenta","red"),
                nitt=50000,burnin=10000) #high run
                #nitt=10000,burnin=5000) #original
                #nitt=100,burnin=10) #For testing that everything is running

#The plots ('plot.ind' and 'plot.col' options) are optional. They plot accepted posterior hybrid
#index values in real time, in this case for 5 randomly chosen individuals.

#'esth' has more functionality - this above just shows the basics#

#Take a look at the results#

#hindlabel

#data.tables sometimes have a strange habit of not showing up the first
#time - if that happens just run the above line again.

###

#Plot a subset of the estimated hybrid indices (the resulting object 'abc' is useful for making a legend)#
#pdf('./hybrid_index.pdf',8,8) #To make a pdf of the graph
#par(mar=c(5,5,4,2)) #To resize the margins of the graph on the pdf

setkey(hindlabel$hi,POPID)	#function from data.table, for rapid sorting and subsetting#

par(mar=c(5,5,4,5)) #To resize the margins of the graph for hybrid index
#
abc = plot_h(data=hindlabel$hi[c("02SS","05RO","06QP","08RU","09PR","10CG")],#Subset of POPIDs#
             test.subject=hindlabel$test.subject,
             mean.h.by="POPID",			#Calculate the mean hybrid index for each value of the "POPID" column#
             sort.by=c("mean_h","POPID","h_posterior_mode"),	#Order test subjects along the x axis by the mean hybrid
             #index calculated above and also by individual hybrid index#
             col.group="POPID",
             group.sep="POPID",
             fill.source=TRUE,
             basic.lines=FALSE,
             source.col=c("dodgerblue","firebrick"),
             source.limits=c("dodgerblue","firebrick"),
             cex=1,pch=16,
             cex.lab=1.5,cex.main=1.5,ylim=c(0,1),
             las=1)
#

#Reshape the plot window as you want#

#Add a legend using the 'plot_h' object 'abc'#

setkey(abc,rn)		#Order data by row number#
legend(72,0.70,	#Place the legend in the top left of the figure# @KML: Overwrote legend position to outside the plot
       c("2","5","6","8","9","10"),   # @KML: Override the legend names to have the right IDs for CG and SS. Careful as it overrides the data itself.
       # abc[,POPID], 		#Name of the field by which data point colours are grouped# (@KML: These are the original legend labels)
       bg="white",			#Background colour#
       text.col=c("black"), #Text colour#
       pch=22, 				#Text size#
       col=abc[,col.Dark2], #Name of the field containing colour information#
       pt.bg=abc[,col.Dark2],	#Name of the field containing colour information#
       ncol=1,				#Number of columns for the group names#
       cex=1, pt.cex=1,
       xpd = TRUE,
       title = "Population")
#dev.off() #to turn off the pdf, if we were separating the graphs to their own pdfs
###

#Run genomic cline analysis#

gc1=ggcline(
  data.prep.object=prepdata$data.prep,
  esth.object=hindlabel,
  read.data.precols=dat$precols,
  #plot.test.subject=c("68177_27"),	#Optional (This is the loci ID for SNP in BCO2)
  #plot.col=c("orange"),		#Optional
  #plot.ylim=c(-2,3),			#Optional
  #nitt = 200, burnin = 20, print.k = 50) #For testing that everything is running
  #nitt = 10000, burnin = 5000, print.k = 50) #original
  nitt = 50000, burnin = 10000, print.k = 50)

#Again the real-time plots are optional. Open circles are parameter v (width),
#'plus' symbols are logit(centre). centre ranges from 0 to 1 (the value of the hybrid index
#at which allele frequencies for the locus are halfway between the parental values).
#The null value for centre is 0.5, but logit(centre) is plotted here to
#try and improve clarity. logit(centre) ranges from -Inf to Inf and the null value is 0.

#Take a look#

#gc1	#Run again if it doesn't show up first time#

###
# @KML created a loop to plot all the loci based on a list of selected IDs

#Plot the cline curve for one locus, and add a title and axis labels#
SNP_labels <- read.delim("./snp_chromosome_coordinates_r50.tsv",header = FALSE,stringsAsFactors = FALSE)

#Loop over the desired SNPs you want to plot, this will plot them all in a single pdf
for (row_index in 1:nrow(SNP_labels)){
  locus_of_interest <- SNP_labels[row_index,]$V1 #pulling out all the locus IDs
  chromosome <- SNP_labels[row_index,]$V2 #Pulling out all the region names for plotting
  basepair <- SNP_labels[row_index,]$V3 #Pulling out the basepair for labeling
  par(mar=c(5,5,4,1)) #set margins of graphing space of cline graph
  # Verify if the data is empty for a locus, if yes skip
  n = nrow(subset(prepdata$data.prep, locus == locus_of_interest))
  if (n==0){
    next
  }
  # Make the Plot for the locus
  plot_clinecurve(
    ggcline.object=gc1$gc,
    data.prep.object=prepdata$data.prep,
    esth.object=hindlabel$hi,
    cline.locus=locus_of_interest,
    cline.col="firebrick",
    null.line.locus=locus_of_interest,	#A diagonal dashed line representing the null expectation#
    null.line.col="black",
    plot.data=locus_of_interest,		#Plot the data too#
    data.col="black",
    plot.genotype=TRUE,		#Plot genotypes, scaled from 0 to 1, rather than alleles (which would all be 0 or 1)#
    PLOIDY=2,			#Needed when plotting genotypes#
    cline.centre.line=locus_of_interest,#Plot the cline centre#
    cline.centre.col="firebrick",
    las=1,
    ylim=c(0,1))
  title(main = paste("SNP ID:",locus_of_interest,"at",chromosome, basepair),xlab="Hybrid index",
        ylab="Genotype frequency",cex.main=1.5,cex.lab=1.5)
}
dev.off()


#@KML This next section is to extract the data tables that have 95 CI, the cline parameter estimates etc.
write.table(prepdata$geno.data,"genotype_clines_r50.txt",quote = FALSE, row.names = FALSE)
write.table(prepdata$locus.data,"locus_clines_r50.txt",quote = FALSE, row.names = FALSE)
write.table(hindlabel,file = "HI_table_r50.txt",quote = FALSE, row.names = FALSE)
write.table(gc1,file = "cline_table_r50.txt",quote = FALSE, row.names = FALSE)

#These functions have more options than shown here, and there are also possibilities for
#pooling or fixing parameters in a variety of different ways,
#followed by model comparison of pairs of models using the 'compare.models' function.
