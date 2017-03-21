library(microbiome)
biompath="C:/Users/cfmasarweh/Google Drive/ChadPhD/SimmonsRotation/qingqing_otu_mfc.txt"
biom=read.biom(biompath,new=T) 
Label = rownames(biom$RA.Otus) 
Label2 = rownames(biom$RA.Otus)
meta = data.frame(Label,Label2) #creates a mock metadata file in the absence of one
xorder=c("S1.1","S2.1","X1cc","X1ct","X1Ic","X1It","X2b","X3b","X4b","X5b","X6b","X7b","X2t","X3t","X4t","X5t","X6t","X7t","X3c","X4c","X6c","X7c","X8.1","X9.1")
meta$Label = factor(Label, levels = xorder) #order xaxis the way you want, step 1
meta=meta[order(meta$Label),] #order xaxis the way you want, step 2

setwd("C:/Users/cfmasarweh/Google Drive/ChadPhD/SimmonsRotation/microbiome-master/R")
for (f in list.files(pattern="*.R")) {
  source(f)
} #loads Josh's functions


###Generates network graphs and corresponding stacked bar charts over a range of correlation strengths and OTU incidence values###
#must change the title name and cluster type depending on the desired algorithm. 
#I modified josh's barplot_module to allow for variable plot title
#I have decided to eliminate 
for(i in c(.3,.5,.8)){ #these values represent the % presence in all samples to filter OTUs
  step=i
  setwd("C:/Users/cfmasarweh/Google Drive/ChadPhD/SimmonsRotation/igraphNetworkClusters/")
  for(j in c(.6)){ #these values represent the correlation strength
    stepp=j
    
    ##Generation of network and cluster objects##
    biom_fil=cooccur_filter(RA=biom$RA.Otus,co_per=step) #% incidence filtering
    biom_netw=cooccurrence(biom_fil,taxon = biom$taxon,cor=stepp) #cooccurrence network generation 
    biom_info=cluster_infomap(biom_netw$netw) #clustering algorithm
    
    ##ZiPi plots##
    zipiname= paste("ZiPi incidence=", step, "strength=", stepp, ".pdf")
    zplottitle=paste("Interactions within and between Ecological Modules\nincidence=", step, "strength=", stepp)
    biom_zipi=ZiPi(biom_netw$netw,modules=biom_info$membership)
    pdf(file=zipiname, width=7, height=7)
    par(xpd=F) #required to allow editing of the plot after creation and before export
    plot(biom_zipi$P,biom_zipi$Z,
         ylim = c(min(biom_zipi$Z, na.rm = T)-0.5, max(biom_zipi$Z, na.rm = T)+0.5), #value of ylim depends on min/max of zi 
         ylab="Within-module degree (hubs)",
         xlab="Between-module connectivity(connectors)",
         main=zplottitle
    )
    abline(v=0.62)
    abline(h=2.5)
    points(biom_zipi$P[biom_zipi$P>=0.62],biom_zipi$Z[biom_zipi$P>=0.62],col="red",pch=1)
    points(biom_zipi$P[biom_zipi$Z>=2.5],biom_zipi$Z[biom_zipi$Z>=2.5],col="blue",pch=1)
    text(0,max(biom_zipi$Z, na.rm = T)+0.4, labels="Module hubs", adj=c(0,0)) #upper  left, left/bottom justified
    text(max(biom_zipi$P, na.rm = T), max(biom_zipi$Z, na.rm = T)+0.4, labels="Network hubs", adj=c(1,0)) #upper right, right/top justified
    text(0, min(biom_zipi$Z, na.rm = T)-0.5, labels="Peripherals", adj=c(0,0)) #lower left, left/bottom justified
    text(max(biom_zipi$P, na.rm = T), min(biom_zipi$Z, na.rm = T)-0.5, labels="Connectors", adj=c(1,0)) #lower right
    dev.off
    graphics.off() #for some reason, dev.off doesn't really turn it off, I need graphics.off
    
    ##ZiPi tables (tables of keystone taxa)##
    biom_zipi_edited=biom_zipi[biom_zipi$Z != "NaN",]
    nethub_a=as.vector(biom_zipi_edited$names[biom_zipi_edited$P>.62 & biom_zipi_edited$Z>2.5]) #the result is a factor, turn it iinto a vector
    modhub_a=as.vector(biom_zipi_edited$names[biom_zipi_edited$P<.62 & biom_zipi_edited$Z>2.5]) #the result is a factor, turn it iinto a vector
    connect_a=as.vector(biom_zipi_edited$names[biom_zipi_edited$P>.62 & biom_zipi_edited$Z<2.5]) 
    connect_a=connect_a[!is.na(connect_a)]
    modhub_a=modhub_a[!is.na(modhub_a)]
    nethub_a=nethub_a[!is.na(nethub_a)]
    write.table(biom$taxon[biom$taxon$otus %in% modhub_a,], 
                file = paste("C:/Users/cfmasarweh/Google Drive/ChadPhD/SimmonsRotation/igraphNetworkClusters/ZiPi_tables/ModuleHubs",
                             "incidence=",step,"strength=",stepp,".txt"), sep="\t", col.names = NA)
    write.table(biom$taxon[biom$taxon$otus %in% nethub_a,], 
                file = paste("C:/Users/cfmasarweh/Google Drive/ChadPhD/SimmonsRotation/igraphNetworkClusters/ZiPi_tables/NetworkHubs",
                             "incidence=",step,"strength=",stepp,".txt"), sep="\t", col.names = NA)
    write.table(biom$taxon[biom$taxon$otus %in% connect_a,], 
                file = paste("C:/Users/cfmasarweh/Google Drive/ChadPhD/SimmonsRotation/igraphNetworkClusters/ZiPi_tables/Connectors",
                             "incidence=",step,"strength=",stepp,".txt"), sep="\t", col.names = NA)
    
    ## Objects for custom Network graphs##
  #file name and graph title names
    graphtitleP = paste("Infomap graph, incidence=", step, "strength=", stepp, "\nColored by Phylum, sized by degree") #automatic title generation
    graphtitleF = paste("Infomap graph, incidence=", step, "strength=", stepp, "\nColored by Family, sized by degree")
    graphtitleRaw = paste("Infomap graph, incidence=", step, "strength=", stepp, "\nRaw")
    graphtitleZiPi = paste("Infomap graph, incidence=", step, "strength=", stepp, 
                           "\nBlue = Connector, Red = Module Hub, Yellow = Network Hub, sized by degree")
    filenameP = paste("Infomap graph, incidence=", step, "strength=", stepp,"Phylum",".pdf") #automatic filename generation
    filenameF = paste("Infomap graph, incidence=", step, "strength=", stepp,"Family",".pdf")
    filenameRaw = paste("Infomap graph, incidence=", step, "strength=", stepp,"Raw",".pdf")
    filenameZiPi = paste("Infomap graph, incidence=", step, "strength=", stepp,"ZiPi",".pdf")
    
  #sizing nodes by degree of connectivity 
    V(biom_netw$netw)$vertex_degree <-  degree(biom_netw$netw) #defining the object to be used for node size
  #Object options for drawing polygons around clusters
    cluster_list = vector(mode="list",length = max(biom_info$membership)) #defining the object to be used for cluster highlighting
    for( i in 1:max(biom_info$membership)){ 
      cluster_list[[i]]=which(biom_info$membership == i)
    } #completes creation of list object for highlighting clusters with polygon
    cluster_list_short_NULLs = vector(mode="list",length = max(biom_info$membership)) #defining the object to be used for modified cluster highlighting
    for( i in c(1:max(biom_info$membership))){
      if(length(which(biom_info$membership==i)) > 3){ #set the number of desired elements here
        cluster_list_short_NULLs[[i]] = which(biom_info$membership == i)
      }
    } #completes creation of list object for highlighting clusters with set number of elements
    is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null)) #creates fxn for identifying Nulls in the cluster list
    rmNullObs <- function(x) {
      x <- Filter(Negate(is.NullOb), x)
      lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
    } #creates fxn for removing NULLs for the cluster list
    cluster_list_short=rmNullObs(cluster_list_short_NULLs) #creates new object with NULLs removed from cluster list
  #Creates objects for coloring the network by keystone taxa
    netw=V(biom_netw$netw)$name
    connectors_in_netw=vector(length = length(netw))
    connectors_in_netw[which(netw %in% connect_a)]=connect_a
    connectors_in_netw[which(connectors_in_netw==FALSE)]=NA
    modhubs_in_netw=vector(length = length(netw))
    modhubs_in_netw[which(netw %in% modhub_a)]=modhub_a
    modhubs_in_netw[which(modhubs_in_netw==FALSE)]=NA
    nethubs_in_netw=vector(length = length(netw))
    nethubs_in_netw[which(netw %in% nethub_a)]=nethub_a
    nethubs_in_netw[which(nethubs_in_netw==FALSE)]=NA
    #alternate approach
    #modhub_b = rep(NA, length(netw))
    #modhub_b[which(netw %in% modhub_a)] <- which(netw%in%modhub_a)
    #modhub_in_new=netw[modhub_b]
    
    intermed1 = pmin(connectors_in_netw, modhubs_in_netw, na.rm = TRUE)
    keystone_in_netw_names = pmin(intermed1, nethubs_in_netw, na.rm = TRUE)
    
    connectors_in_netw[!is.na(connectors_in_netw)] = "red"
    modhubs_in_netw[!is.na(modhubs_in_netw)] = "blue"
    nethubs_in_netw[!is.na(nethubs_in_netw)] = "yellow"
    intermed2 = pmin(connectors_in_netw, modhubs_in_netw, na.rm = TRUE)
    keystone_in_netw_colors = pmin(intermed2, nethubs_in_netw, na.rm = TRUE)
    
    keystone_in_netw_colors[is.na(keystone_in_netw_colors)] = "grey"
    
    keystone_in_netw_Family = as.character(biom$taxon$Family[match(keystone_in_netw_names,biom$taxon$otu)])
    
    #Creating network graph and exporting
    pdf(file=filenameP,width=7, height=7) #export files
    set.seed(4) #keeps the graphs from randomly generating, value is arbitrary
    plot(biom_netw$netw,
         vertex.color=as.factor(biom$taxon$Phylum[match(V(biom_netw$netw)$name,biom$taxon$otu)]),
         mark.groups = cluster_list_short, #circles certain clusters
         #try setting this to circle clusters involved in nodes with high degrees of connectivity
         edge.color="grey22", 
         vertex.size=V(biom_netw$netw)$vertex_degree/3, #sizes nodes by degree of connectivity 
         vertex.label=NA, #turns off node label
         main= graphtitleP #automatic title generation
    )
    dev.off()
    pdf(file=filenameF,width=7, height=7) #export files
    set.seed(4) #keeps the graphs from randomly generating, value is arbitrary
    plot(biom_netw$netw,
         vertex.color=as.factor(biom$taxon$Family[match(V(biom_netw$netw)$name,biom$taxon$otu)]),
         mark.groups = cluster_list_short, #circles certain clusters
         #try setting this to circle clusters involved in nodes with high degrees of connectivity
         edge.color="grey22", 
         vertex.size=V(biom_netw$netw)$vertex_degree/3, #sizes nodes by degree of connectivity 
         vertex.label=NA, #turns off node label
         main= graphtitleF #automatic title generation
    )
    dev.off()
    pdf(file=filenameRaw,width=7, height=7) #export files
    set.seed(4) #keeps the graphs from randomly generating, value is arbitrary
    plot(biom_netw$netw,
         vertex.color="orange",
         #try setting this to circle clusters involved in nodes with high degrees of connectivity
         edge.color="grey22", 
         vertex.size=3, #sizes nodes by degree of connectivity 
         vertex.label=NA, #turns off node label
         main= graphtitleRaw #automatic title generation
    )
    dev.off()
    pdf(file=filenameRaw,width=7, height=7) #export files
    set.seed(4) #keeps the graphs from randomly generating, value is arbitrary
    plot(biom_netw$netw,
         vertex.color="orange",
         #try setting this to circle clusters involved in nodes with high degrees of connectivity
         edge.color="grey22", 
         vertex.size=3,
         vertex.label=NA, #turns off node label
         main= graphtitleRaw #automatic title generation
    )
    dev.off()
    pdf(file=filenameZiPi,width=7, height=7) #export files
    set.seed(4) #keeps the graphs from randomly generating, value is arbitrary
    plot(biom_netw$netw,
         vertex.color=keystone_in_netw_colors,
         edge.color="grey22", 
         vertex.size=V(biom_netw$netw)$vertex_degree/3, #sizes nodes by degree of connectivity 
         vertex.label=keystone_in_netw_Family, #Gives Family name to keystone taxa
         vertex.label.dist=0.25, #offsets the Family names by a little bit
         vertex.label.font=2,
         vertex.label.color=keystone_in_netw_colors,
         main= graphtitleZiPi #automatic title generation
    )
    dev.off()
    
    ###Stacked Bar Plots of cluster presence per sample##
    plottitle = paste("Infomap clustering by sample\n, incidence=", step, "strength=", stepp) #automatic title generation
    filename2 = paste("Infomap clustering by sample, incidence=", step, "strength=", stepp,".pdf") #automatic filename generation
    plot_module=barplot_module(data=biom$RA.Otus,niche = biom_info,meta = meta,categories = "Label",plottitle = plottitle)
    ggsave(plot_module$plot, filename=filename2, width=7, height=7) #export files
  }
}






################################################################################################################################

###For generating cluster RA stacked bar plot###
#to generate the barplot_module plots, one must manually generate each one,
#depending on the desired correlation strengths and OTU incidence. 
#select parameters based on the results of the above forloop
biom_fil=cooccur_filter(RA=biom$RA.Otus,co_per=0.5)
biom_netw=cooccurrence(biom_fil,taxon = biom$taxon,cor=0.5,conn=1)
biom_info=cluster_infomap(biom_netw$netw)
#plot_module=barplot_module(data=biom$RA.Otus,niche = biom_info,meta = meta,categories = "Label")
#plot_module$plot
set.seed(4)
plot(biom_netw$netw,vertex.color=as.factor(biom_info$membership), edge.color="black", vertex.size=4, vertex.label="", main= "Infomap, cor=0.8, pvalue=0.01")


###############################################################################################################################

###For troubleshooting cooccur_filter, cooccurrence, and barplot_module###
#this set will use the tutorial files to generate plots, so if it doesn't work, something's wrong
setwd("C:/Users/cfmasarweh/Google Drive/ChadPhD/SimmonsRotation/microbiome-master/data")
load("metadata.rda")
biom = read.biom(fir_data, new = F)
setwd("C:/Users/cfmasarweh/Google Drive/ChadPhD/SimmonsRotation/Troubleshooting")
biom_fil=cooccur_filter(RA=biom$RA.Otus,co_per=.5)
biom_netw=cooccurrence(biom_fil,taxon = biom$taxon,cor=.6)
biom_info=infomap.community(biom_netw$netw)
V(biom_netw$netw)$vertex_degree <-  degree(biom_netw$netw) #defining the object to be used for node size
a_list = vector(mode="list",length = max(biom_info$membership))
for( i in 1:max(biom_info$membership)){
  a_list[[i]]=which(biom_info$membership==i)}


#Labeling only those OTUs that are keystone taxa according to zipi
b=as.integer(V(biom_netw$netw)$name) #vector of OTU names, the order of which reflects the order in which nodes are plotted
a=as.vector(biom_zipi$names[biom_zipi$P>.62 & biom_zipi$Z>1]) #the result is a factor, turn it iinto a vector
a#character class and character type
a=a[!is.na(a)] #same class and type
a=as.integer(a)#OTUs where the Zi and Pi parameters above are true
c=biom$taxon[biom$taxon$otus %in% a,] #returns biom$taxon info where otus meet certain zipi values
c=as.data.frame(c)#creates list type, dataframe class of c
d=b[match(b,a,nomatch=NA)] #returns vector of length equal to b, where all non-matches to a are NA
set.seed(4)
plot.igraph(biom_netw$netw, 
     vertex.color=as.factor(biom$taxon$Family[match(V(biom_netw$netw)$name,biom$taxon$otu)]), 
     edge.color="grey22", 
     mark.groups =a_list, 
     vertex.label=NA, #labels only the keystone taxa, as defined by b
     vertex.label.dist=0,
     vertex.size=3, 
     main= "bla"
     )
c=biom$taxon$Family
set.seed(4)
plot.igraph(biom_netw$netw, 
            vertex.color="grey", 
            edge.color="grey22", 
            #mark.groups =b_list, 
            vertex.label=NA, #labels only the keystone taxa, as defined by b
            vertex.size=3, 
            main= "bla"
)
plot_module=barplot_module(data=biom$RA.Otus,niche = biom_info,meta = meta,categories = "Label")
plot_module$plot

biom_zipi=ZiPi(biom_netw$netw,modules=biom_info$membership)
par(xpd=F)
plot(biom_zipi$P,biom_zipi$Z,
     ylim = c(min(biom_zipi$Z, na.rm = T)-0.5, max(biom_zipi$Z, na.rm = T)+0.5), #value of ylim depends on min/max of zi 
     ylab="Zi",
     xlab="Pi"
)
abline(v=0.62)
abline(h=2.5)
points(biom_zipi$P[biom_zipi$P>=0.62],biom_zipi$Z[biom_zipi$P>=0.62],col="red",pch=1)
points(biom_zipi$P[biom_zipi$Z>=2.5],biom_zipi$Z[biom_zipi$Z>=2.5],col="blue",pch=1)



######################################################################################################################################

##Below creates a bunch of files into the WD, which display the info in key objects created in this pipeline##
#-->I want to plot the relative abundance of each OTU within each cluster. 
#List of OTUs is biom_info$names; list of cluster membership of OTUs is biom_info$membership
write.table(biom$taxon, file = "biom$taxon.txt",  sep="\t", row.names=F) #I think this is each taxon and the number of OTUs it represents
write.table(biom$RA.Otus, file = "biom$RA.Otus.txt",  sep="\t") #RA of each oTU (colname) by sample (rowname)
write.table(biom$biom_tab, file = "biom$biom_tab.txt",  sep="\t") #a table with samples as header, OTUs as rowname, and taxon as terminal row cell, with OTU counts as values
#write.table(biom_info, file = "biom_info.txt",  sep="\t") #igraph object type, not sure what it is
#write.table(biom_netw$corr, file = "biom_netw$corr.txt",  sep="\t") #igraph object type, not sure what it is
write.table(biom_netw$corrMin, file = "biom_netw$corrMin.txt",  sep="\t") #matrix of OTU x OTU, with 1 or 0 as value. dunno what it is
write.table(biom_netw$taxon.netw, file = "biom_netw$taxon.netw.txt",  sep="\t") #I think this is each taxon and the number of OTUs it represents, filtered
write.table(biom_netw$pAdjusted, file = "biom_netw$pAdjusted.txt",  sep="\t") #dunno what this is
write.table(meta, file = "meta.txt", sep = "\t")
write.table(biom_info$membership, file = "biom_info$membership.txt",  sep="\t") #right column is cluster number, left is vcount
write.table(biom_info$names, file = "biom_info$names.txt",  sep="\t") #left column is vcount, right is OTU number






###The following is explanations of each aspect of the workflow for generating a network graph###
##Large, illegible network
#filter OTU's to 50% presence in all samples. RA_otus is the list object with the numbers in the final plot, and is known in this fxn as RA
#final product of cooccur_filter is biom_fil, which is the filtered RA_otus
biom_fil=cooccur_filter(RA=biom$RA.Otus,co_per=0.2)

#run co-occurence. Taxon can be excluded and identified later if desired. 
#this takes biom_fil (i.e. filtered TA_otus; in this case, it's the fxn argument "data") and does stuff with it
biom_netw=cooccurrence(biom_fil,taxon = biom$taxon,cor=0.8)
#this plot will change every time you run it
plot(biom_netw$netw,  edge.color="black", vertex.size=4, vertex.label="", main= "Base graph, cor=0.8, pvalue=0.01")
##

#we've generated our filtered, modified, plottable matrix, now we can do other things with it using, again, igraph
#we have the option of using different matrices 
#at some point, try to forloop these commands, maybe even with a series of correlation strengths
biom_info=cluster_louvain(biom_netw$netw)
#biom_info$names are the OTU numbers

#these plots will change every time you run them
#biom_info$membership, colors by cluster membership
#biom_netw$netw is an igraph object, which seems to be a list of connections
plot(biom_netw$netw,vertex.color=as.factor(biom_info$membership), edge.color="black", vertex.size=4, vertex.label="", main= "Infomap, cor=0.8, pvalue=0.01")