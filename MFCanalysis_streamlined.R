library(microbiome) #loads Josh's microbiome package
for (f in list.files(path = "C:/Users/cfmasarweh/Google Drive/ChadPhD/SimmonsRotation/microbiome-master/R", pattern="*.R")) {
  setwd("C:/Users/cfmasarweh/Google Drive/ChadPhD/SimmonsRotation/microbiome-master/R")
  source(f)
} #loads Josh's functions. Be sure to specify the path in which Josh's functionss reside.

#Loading the input files 
biom=read.biom("C:/Users/cfmasarweh/Google Drive/ChadPhD/SimmonsRotation/qingqing_otu_mfc.txt",new=T) 
meta=read.table("C:/Users/cfmasarweh/Google Drive/ChadPhD/SimmonsRotation/meta_MFC_fake.txt")

#Parameters for you to set
xorder=c("S1.1","S2.1","X1cc","X1ct","X1Ic","X1It","X2b","X3b","X4b","X5b","X6b","X7b","X2t","X3t","X4t","X5t","X6t","X7t","X3c","X4c","X6c","X7c","X8.1","X9.1") #order of the x axis
meta=meta[order(meta$Label),] #order xaxis the way you want, step 2
clusteringAlgorithm = "cluster_optimal" #write the igraph clustering fxn in quotes. See documentation for full list of choices
incidence_parameters = c(0.8) #example: c(0.3,0.5,0.6,0.8,0.9)
correlation_strength_parameters = c(0.6)
zipitablefilepath = "C:/Users/cfmasarweh/Google Drive/ChadPhD/SimmonsRotation/igraphNetworkClusters/ZiPi_tables/" #gotta have that last forward slash
plotsfilepath = "C:/Users/cfmasarweh/Google Drive/ChadPhD/SimmonsRotation/igraphNetworkClusters/" #gotta have that last forward slash

#must change the title name and cluster type depending on the desired algorithm. 
#I modified josh's barplot_module to allow for variable plot title

#The follow generates the tables and plots#
for(i in incidence_parameters){ #these values represent the % presence in all samples to filter OTUs
  incparam=i
  for(j in correlation_strength_parameters){ #these values represent the correlation strength
    strengthparam=j
    
    ##Generation of network and cluster objects##
    biom_fil=cooccur_filter(RA=biom$RA.Otus,co_per=incparam) # % incidence filtering
    biom_netw=cooccurrence(biom_fil,taxon = biom$taxon,cor=strengthparam) #cooccurrence network generation 
    if(clusteringAlgorithm == "cluster_infomap"){
      biom_netw_clust=cluster_infomap(biom_netw$netw)
    } else if(clusteringAlgorithm == "cluster_louvain"){ #louvain, aka multilevel optimization
      biom_netw_clust=cluster_louvain(biom_netw$netw)
    } else if(clusteringAlgorithm == "cluster_fast_greedy"){ #fast greedy
      biom_netw_clust=cluster_fast_greedy(biom_netw$netw)
    } else if(clusteringAlgorithm == "cluster_optimal"){ #optimal
      if(winDialog(type = "yesno", message = "The cluster_optimal algorithm uses a ton of memory and may crash your machine.\nDo you want to proceed?")=="YES"){
        biom_netw_clust=cluster_optimal(biom_netw$netw)
      } else {winDialog(type = "ok", message = "Very good, Sir.\nYou may choose a new algorithm at your leisure, Sir")}
    } else if(clusteringAlgorithm == "cluster_edge_betweenness"){ #edge betweenness
      biom_netw_clust=cluster_edge_betweenness(biom_netw$netw)
    } else if(clusteringAlgorithm == "cluster_label_prop"){ #propagating labels
      biom_netw_clust=cluster_label_prop(biom_netw$netw)
    } else if(clusteringAlgorithm == "cluster_leading_eigen"){ #leading eigenvector
      biom_netw_clust=cluster_leading_eigen(biom_netw$netw)
    } else if(clusteringAlgorithm == "cluster_spinglass"){ #spinglass model and simulated annealing
      biom_netw_clust=cluster_spinglass(biom_netw$netw)
    } else if(clusteringAlgorithm == "cluster_walktrap"){ #short, random walks
      biom_netw_clust=cluster_walktrap(biom_netw$netw)
    } else {
      print("Please check your cluster algorithm spelling")
    }
    
    ##ZiPi plots##
    zipifilename= paste("ZiPi incidence=", incparam, "strength=", strengthparam, ".pdf")
    zipiplottitle=paste("Interactions within and between Ecological Modules\nincidence=", incparam, "strength=", strengthparam)
    biom_zipi=ZiPi(biom_netw$netw,modules=biom_netw_clust$membership)
    pdf(file=paste(plotsfilepath,zipifilename), width=7, height=7)
    par(xpd=F) #required to allow editing of the plot after creation and before export
    plot(biom_zipi$P,biom_zipi$Z,
         ylim = c(min(biom_zipi$Z, na.rm = T)-0.5, max(biom_zipi$Z, na.rm = T)+0.5), #value of ylim depends on min/max of zi 
         ylab="Within-module degree (hubs)",
         xlab="Between-module connectivity(connectors)",
         main=zipiplottitle
    )
    abline(v=0.62) #visually demarcate the boundary between connectors and peripherals
    abline(h=2.5) #visually demarcate the boundary between hubs and peripherals
    points(biom_zipi$P[biom_zipi$P>=0.62],biom_zipi$Z[biom_zipi$P>=0.62],col="blue",pch=1) #color connectors blue 
    points(biom_zipi$P[biom_zipi$Z>=2.5],biom_zipi$Z[biom_zipi$Z>=2.5],col="red",pch=1) #color hubs red
    text(0,max(biom_zipi$Z, na.rm = T)+0.4, labels="Module hubs", adj=c(0,0)) #upper  left, left/bottom justified
    text(max(biom_zipi$P, na.rm = T), max(biom_zipi$Z, na.rm = T)+0.4, labels="Network hubs", adj=c(1,0)) #upper right, right/top justified
    text(0, min(biom_zipi$Z, na.rm = T)-0.5, labels="Peripherals", adj=c(0,0)) #lower left, left/bottom justified
    text(max(biom_zipi$P, na.rm = T), min(biom_zipi$Z, na.rm = T)-0.5, labels="Connectors", adj=c(1,0)) #lower right
    dev.off
    graphics.off() #for some reason, dev.off doesn't really turn it off, I need graphics.off, too
    
    ##ZiPi tables (tables of keystone taxa)##
    #zipitablefilepath = "C:/Users/cfmasarweh/Google Drive/ChadPhD/SimmonsRotation/igraphNetworkClusters/ZiPi_tables"
    biom_zipi_edited=biom_zipi[biom_zipi$Z != "NaN",] #gets rid of lines with non-number Z values
    nethub_a=as.vector(biom_zipi_edited$names[biom_zipi_edited$P>=.62 & biom_zipi_edited$Z>=2.5]) #vector of names of OTUs that are network hubs
    modhub_a=as.vector(biom_zipi_edited$names[biom_zipi_edited$P<=.62 & biom_zipi_edited$Z>=2.5]) #vector of names of OTUs that are module hubs
    connect_a=as.vector(biom_zipi_edited$names[biom_zipi_edited$P>=.62 & biom_zipi_edited$Z<=2.5]) ##vector of names of OTUs that are connectors
    connect_a=connect_a[!is.na(connect_a)] #gets rid of NAs this might be uncecessary given biom_zipi_edited
    modhub_a=modhub_a[!is.na(modhub_a)] #gets rid of NAs
    nethub_a=nethub_a[!is.na(nethub_a)] #gets rid of NAs
    write.table(biom$taxon[biom$taxon$otus %in% modhub_a,], 
                file = paste(zipitablefilepath,
                             "modulehub_incidence=",incparam,"_strength=",strengthparam,".txt", sep = ""), 
                sep="\t", col.names = NA)
    write.table(biom$taxon[biom$taxon$otus %in% nethub_a,], 
                file = paste(zipitablefilepath,
                             "netwokhub_incidence=",incparam,"_strength=",strengthparam,".txt", sep = ""), 
                sep="\t", col.names = NA)
    write.table(biom$taxon[biom$taxon$otus %in% connect_a,], 
                file = paste(zipitablefilepath,
                             "connectors_incidence=",incparam,"_strength=",strengthparam,".txt", sep = ""), 
                sep="\t", col.names = NA)
    
    ## Objects for custom Network graphs##
    #file name and graph title names
    graphtitleP = paste("Infomap graph, incidence=", incparam, "strength=", strengthparam, "\nColored by Phylum\nsized by degree") #automatic title generation
    graphtitleF = paste("Infomap graph, incidence=", incparam, "strength=", strengthparam, "\nColored by Family\nsized by degree")
    graphtitleRaw = paste("Infomap graph, incidence=", incparam, "strength=", strengthparam, "\nRaw")
    graphtitleZiPi = paste("Infomap graph, incidence=", incparam, "strength=", strengthparam, 
                           "\nBlue = Connector, Red = Module Hub, Yellow = Network Hub, sized by degree")
    filenameP = paste("Infomap graph, incidence=", incparam, "strength=", strengthparam,"Phylum",".pdf") #automatic filename generation
    filenameF = paste("Infomap graph, incidence=", incparam, "strength=", strengthparam,"Family",".pdf")
    filenameRaw = paste("Infomap graph, incidence=", incparam, "strength=", strengthparam,"Raw",".pdf")
    filenameZiPi = paste("Infomap graph, incidence=", incparam, "strength=", strengthparam,"ZiPi",".pdf")
    
    #sizing nodes by degree of connectivity 
    V(biom_netw$netw)$vertex_degree <-  degree(biom_netw$netw) #defining the object to be used for node size
    #Object options for drawing polygons around clusters
    cluster_list = vector(mode="list",length = max(biom_netw_clust$membership)) #defining the object to be used for cluster highlighting
    for( i in 1:max(biom_netw_clust$membership)){ 
      cluster_list[[i]]=which(biom_netw_clust$membership == i)
    } #completes creation of list object for highlighting clusters with polygon
    cluster_list_short_NULLs = vector(mode="list",length = max(biom_netw_clust$membership)) #defining the object to be used for modified cluster highlighting
    for( i in c(1:max(biom_netw_clust$membership))){
      if(length(which(biom_netw_clust$membership==i)) > 3){ #set the number of desired elements here
        cluster_list_short_NULLs[[i]] = which(biom_netw_clust$membership == i)
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
    
    connectors_in_netw[!is.na(connectors_in_netw)] = "blue"
    modhubs_in_netw[!is.na(modhubs_in_netw)] = "red"
    nethubs_in_netw[!is.na(nethubs_in_netw)] = "yellow"
    intermed2 = pmin(connectors_in_netw, modhubs_in_netw, na.rm = TRUE)
    keystone_in_netw_colors = pmin(intermed2, nethubs_in_netw, na.rm = TRUE)
    
    keystone_in_netw_colors[is.na(keystone_in_netw_colors)] = "grey"
    
    keystone_in_netw_Family = as.character(biom$taxon$Family[match(keystone_in_netw_names,biom$taxon$otu)])
    
    #Creating network graph and exporting
    pdf(file=paste(plotsfilepath,filenameP),width=7, height=7) #export files
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
    pdf(file=paste(plotsfilepath,filenameF),width=7, height=7) #export files
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
    pdf(file=paste(plotsfilepath,filenameRaw),width=7, height=7) #export files
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
    pdf(file=paste(plotsfilepath,filenameRaw),width=7, height=7) #export files
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
    pdf(file=paste(plotsfilepath,filenameZiPi),width=7, height=7) #export files
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
    plottitle = paste("Infomap clustering by sample\n, incidence=", incparam, "strength=", strengthparam) #automatic title generation
    filename2 = paste("Infomap clustering by sample, incidence=", incparam, "strength=", strengthparam,".pdf") #automatic filename generation
    plot_module=barplot_module(data=biom$RA.Otus,niche = biom_netw_clust,meta = meta,categories = "Label",plottitle = plottitle)
    ggsave(plot_module$plot, filename=paste(plotsfilepath,filename2), width=7, height=7) #export files
  }
}