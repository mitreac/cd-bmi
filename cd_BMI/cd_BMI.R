library(org.Hs.eg.db)
library(devtools)
# library(arcdiagram) # devtools::install_github('gastonstat/arcdiagram')
library(ROntoTools)
library(plyr)
library(mixtools)
library(ggplot2)
library(RColorBrewer)
#install.packages("gplots")
library(gplots)

# script that contains the implementation of the Change Detection method
# not released yet by the authors
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7235093/
source(ChangeDetector.R)


# FUNCTIONS

# __________________________________________________________________________

# function used to get the graph data for all KEGG pathways
downloadKeggPathwayNames = function() {
  require(ROntoTools)
  kpg <- keggPathwayGraphs("hsa",updateCache=T)
  kpg <- setEdgeWeights(kpg)
  kpg <- setNodeWeights(kpg, defaultWeight = 1)
  
  kpgnames <- keggPathwayNames(organism = "hsa", updateCache = T, verbose = TRUE)
  cache <- list(kpgnames=kpgnames, kpg=kpg)
  saveRDS(cache,file=paste("keggData", format(Sys.time(), "%Y-%b-%d"), ".rds",sep=""))
  return (cache)
}

# __________________________________________________________________________


# function used to get the gene expression data for a pathway
get_mydata = function(all_data, mysystem){
  # map the gene symbols to entrez ids to map them to KEGG pathways later
  mymapping = toTable(org.Hs.egSYMBOL)
  mymappingEntrezID = mymapping$gene_id
  names(mymappingEntrezID) = mymapping$symbol
  
  myEntrezIds=mymappingEntrezID[rownames(all_data)]
  
  all_data1 = myQBCdata[!is.na(myEntrezIds),]
  rownames(all_data1) = paste("hsa:",as.character(myEntrezIds[!is.na(myEntrezIds)]),sep="")
  
  mydata <- all_data1[rownames(all_data1) %in% mysystem[[1]]@nodes,]
  return(mydata)
}

# __________________________________________________________________________

# function used to get the pathway graph for a pathway
get_mysystem = function(keggData, pathwayID){
  mysystem <- keggData$kpg[pathwayID] 
  return(mysystem)
}

# __________________________________________________________________________

# function used to plot the pathway PF distribution and the mixture of the gamma distributions
plot_distributions = function(myPFs, mymix, myxCoord, maxy, file_name){
  tiff(file_name, units="in", width=9, height=5, res=120)
  plot(density(myPFs, from = 0), ylim = c(0,maxy), lwd = 4, col = "black", main = "", xlab = "system perturbation")
  curve(dgamma(x, shape = mymix$gamma.pars[1,1], scale=mymix$gamma.pars[2,1]),  lwd=2, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.75), add=TRUE, n=500, lty = 3)    
  curve(dgamma(x, shape = mymix$gamma.pars[1,2], scale=mymix$gamma.pars[2,2]),  lwd=2, col=rgb(red = 1, green = 0, blue = 0, alpha = 0.75), add=TRUE, n=500, lty = 3)    
  curve(mymix$lambda[1]*dgamma(x, shape = mymix$gamma.pars[1,1], scale=mymix$gamma.pars[2,1]) + mymix$lambda[2]*dgamma(x, shape = mymix$gamma.pars[1,2], scale=mymix$gamma.pars[2,2]),  lwd=4, col=rgb(red = 1, green = 0, blue = 1, alpha = 0.75), add=TRUE, n=1000)    
  segments(x0 = myxCoord, y0 = -0.01, x1 = myxCoord, y1 = maxy, col = "yellow", lwd = 5)
  dev.off()
}

# __________________________________________________________________________

# function used to plot arcplots
plot_arcplot = function(file_name, plot_title, myCDs, BMI_groups_labels, pathways){
  pcolors = c("#009E73","turquoise","#56B4E9","#0072B2","#666699","#CC79A7","hotpink","tomato","red")
  names(pcolors) = sort(pathways)
  tiff(file_name, units="in", width=9, height=5, res=120)
  ggplot(data = data.frame(x = 1:7, y = rep(0,7))) +  
    ylim(-0.2,0.5) +
    geom_line(mapping = aes(x = x, y = y), data = data.frame(x = 1:7, y = rep(0,7)), color = "gray") +
    geom_curve(data = data.frame(myCDs[1,, drop = FALSE]), aes(x = X1, y = rep(0,1), xend = X2, yend = rep(0,1), color = pathways[1]), lwd = 1.5, curvature = -0.55) +
    geom_curve(data = data.frame(myCDs[2:4,]), aes(x = X1, y = rep(0,3), xend = X2, yend = rep(0,3), color = pathways[2:4]), lwd = 1.5, curvature = 0.55) + 
    geom_curve(data = data.frame(myCDs[5,, drop = FALSE]), aes(x = X1, y = rep(0,1), xend = X2, yend = rep(0,1), color = pathways[5]), lwd = 1.5, curvature = 0.3) +
    geom_curve(data = data.frame(myCDs[6:9,]), aes(x = X1, y = rep(0,4), xend = X2, yend = rep(0,4), color = pathways[6:9]), lwd = 1.5, curvature = -0.45) + 
    scale_colour_manual(values=pcolors) +
    ggtitle(plot_title) + 
    theme_minimal() + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
          axis.title.x = element_text(margin = margin(b = 0), vjust = 2, size = 14),
          axis.text.x = element_text(angle = 0, vjust = 4, hjust=0.25, margin = margin(b = 0), size = 14),
          panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.margin = unit(c(0,0,0,0), unit = "pt"), legend.position = c(0.5, 0.75), legend.text=element_text(size = 12),
          plot.title = element_text(hjust = 0.5, vjust = -5, size = 18), legend.title = element_text(size = 14)) +
    xlab("BMI range") + 
    scale_x_continuous(breaks=c(1:7), labels=BMI_groups_labels) +
    labs(colour="Pathway") +
    guides(color=guide_legend(ncol=3)) +
    geom_point(shape = 15, size = 3, mapping = aes(x = x, y = y), data = data.frame(x = 1:7, y = rep(0,7))) 
  dev.off()
  
}

# __________________________________________________________________________

# function used to plot heatmap
plot_heatmap = function(mydata, file_name){
  # from symbols to entrezids
  mymapping = toTable(org.Hs.egSYMBOL)
  mymappingEntrezID = mymapping$gene_id
  names(mymappingEntrezID) = mymapping$symbol
  # from entrezids to symbols
  mymappingEntrezID1 = mymappingEntrezID
  mymappingEntrezID1 = paste("hsa", mymappingEntrezID, sep = ":")
  names(mymappingEntrezID1) = names(mymappingEntrezID)
  g_symb = names(mymappingEntrezID1[mymappingEntrezID1 %in% rownames(mydata)]) 
  names(g_symb) = mymappingEntrezID1[mymappingEntrezID1 %in% rownames(mydata)]
  rownames(mydata1) = g_symb[rownames(mydata)]
  
  sd_data = apply(X = mydata, MARGIN = 1, FUN = sd)
  mydata2 = apply(X = mydata[sd_data>quantile(sd_data)[2],], MARGIN = 1, FUN = smooth)
  mydata2 = t(mydata2)
  colnames(mydata2) = colnames(mydata)
  tiff(file_name, units="in", width=9, height=15, res=120)
  heatmap.2(gene_data,Colv = NA, col=rev(brewer.pal(11, "RdBu")), dendrogram = "none", scale = "row", 
            density.info = "none",trace = "none", margins = c(7,7), lhei = c(1,12))
  dev.off()
}

# __________________________________________________________________________

#function to compute the change intervals for a pathway
process_pathway = function(mydata, mysystem){

  myRESCD_QBC=changeDetectorAux(data = mydata, system = mysystem, mm = 1, tsd = 0, tcutoff = 1)
  myPFs=myRESCD_QBC$pPFs
  
  mymean1 = min(myPFs)
  mymean2 = max(myPFs)
  mysd = sd(myPFs)/2
  
  mya1 = mymean1^2/mysd^2
  myb1 = mysd^2/mymean1
  mya2 = mymean2^2/mysd^2
  myb2 = mysd^2/mymean2
  
  set.seed(42)
  mymix = gammamixEM(x = myPFs, maxit=100, k=2, alpha = c(mya1,mya2), beta = c(myb1,myb2))
  
  myxval = c(0:(max(myPFs)*1000))/1000
  myd1 = dgamma(x=myxval, shape = mymix$gamma.pars[1,1], scale=mymix$gamma.pars[2,1])
  myd2 = dgamma(x=myxval, shape = mymix$gamma.pars[1,2], scale=mymix$gamma.pars[2,2])
  
  myxCoord = which(myd2>myd1 & myd2 > 1)[1]/1000
  myxCoord = which(myd2>myd1)[1]/1000
  
  myCDs=getChangeIntervals(stateLabels = colnames(mydata), pPFs = myPFs, pfcutoff = myxCoord)
  return(list(myCDs = myCDs, myPFs = myPFs, myxCoord = myxCoord, mymix = mymix, 
              maxy = max(c(max(myd1), max(myd2), max(myPFs))) ))
}

#================================================================================
#================================================================================

# CODE TO RUN FUNCTIONS - help generate figures from the manuscript
#================================================================================

# gene expression data summarized by BMI group form GSE142102
TNBC_genedata = read.table(file = "TNBC_data_BMI.txt", header = T, row.names = 1, sep = "\t", check.names =  F)


# processed file from GSE189757
ERpBCdata=read.table(file = "Processed_data.txt", header = T, sep = "\t", row.names = NULL)
# create a text file with the BMI group for each sample from GSE142102 
groups_data=read.table(file = "groups_data.txt", header = T, sep = "\t", row.names = 1)
bmi_groupNames = c("18.5-24.9",  "25-27.9", "28-29.9",  "30-31.9",  "32-34",  "34.01-39.9",  ">=40")


myERpBCdata_matrix = ERpBCdata[,2:dim(ERpBCdata)[2]]
myERpBCdata_gene_level = aggregate(x = myERpBCdata_matrix, by = list(ERpBCdata[,1]), FUN = "mean", na.rm=TRUE, na.action=NULL)
myERpBCdata_matrix_temp = myERpBCdata_gene_level[,2:45]
myERpBCdata_matrix_BMI_group = aggregate(t(myERpBCdata_matrix_temp), by = list(groups_data), FUN = "mean", na.rm=TRUE, na.action=NULL)
ERpBC_genedata = t(myERpBCdata_matrix_BMI_group[,2:dim(myERpBCdata_matrix_BMI_group)[2]])
rownames(ERpBC_genedata) = myERpBCdata_matrix_BMI_group[,1]
colnames(ERpBC_genedata) = bmi_groupNames

keggData = downloadKeggPathwayNames()

KEGG_pathways = c(
  "path:hsa04920", #Adipocytokine signaling pathway
  "path:hsa04066", #Hif1a signaling
  "path:hsa04930", #type 2 diabetes
  "path:hsa03320", #PPAR signaling
  "path:hsa04910", #insulin signaling 
  "path:hsa04060", #cytokine-cytokine interaction
  "path:hsa04150", #mTOR
  "path:hsa04151", #pik3ca-akt signaling
  "path:hsa04915"  #estrogen signaling
)

CDintervals_ERpBC = NULL
CDintervals_TNC = NULL
for (pathway_ID in KEGG_pathways){
  mysystem = get_mysystem(keggData = keggData, pathwayID = pathway_ID)
  mydaya_ERpBC = get_mydata(all_data = ERpBC_genedata, mysystem = mysystem)
  mydata_TNBC = get_mydata(all_data = TNBC_genedata, mysystem = mysystem)
  ERp_res = process_pathway(mydata = mydaya_ERpBC, mysystem = mysystem)
  TNBC_res = process_pathway(mydata = mydaya_TNBC, mysystem = mysystem)
  CDintervals_ERpBC = rbind(CDintervals_ERpBC, ERp_res$myCDs)
  CDintervals_TNBC = rbind(CDintervals_TNBC, TNBC_res$myCDs)
}

colnames(CDintervals_ERpBC) = c("X1","X2")
colnames(CDintervals_TNBC) = c("X1","X2")

pathways = keggData$kpgnames[KEGG_pathways]

plot_arcplot("ERpBC_CD_analysis.tiff", 
             "ER+ breast cancer pathway activation analysis", 
             CDintervals_ERpBC, bmi_groupNames, pathways)

plot_arcplot("TNBC_CD_analysis.tiff", 
             "Triple negative breast cancer pathway activation analysis", 
             CDintervals_TNC, bmi_groupNames, pathways)

pathway_ID = "path:hsa04920" #Adipocytokine signaling pathway
mysystem = get_mysystem(keggData = keggData, pathwayID = pathway_ID)
mydata_AdipPath = get_mydata(all_data = ERpBC_genedata, mysystem = mysystem)
AdipPath_res = process_pathway(mydata = mydaya_AdipPath, mysystem = mysystem)

plot_distributions(myPFs = AdipPath_res$myPFs, 
                   mymix = AdipPath_res$mymix, 
                   myxCoord = AdipPath_res$myxCoord, 
                   maxy = AdipPath_res$maxy,
                   file_name = "ERp_Adipocytokine_pathway_distr.tiff")

plot_heatmap(mydata = mydata_AdipPath, file_name = "ERp_Adipocytokine_pathway_heatmap.tiff")

pathway_ID = "path:hsa04066" #Hif1a signaling
mysystem = get_mysystem(keggData = keggData, pathwayID = pathway_ID)
mydata_Hif1aPath = get_mydata(all_data = ERpBC_genedata, mysystem = mysystem)

plot_heatmap(mydata = mydata_Hif1aPath, file_name = "ERp_Hif1a_pathway_heatmap.tiff")



