#dropbox -> Kinabalu community -> data

source("./R/bip_ggnet.R")
source("./R/bip_igplot.R")
source("./R/bip_qtplot.R")
source("./R/bip_railway.R")
source("./R/vectorize.R")
source("./R/bip_edgewt.R")
require(ggplot2)
require(network)
require(igraph)
require(sna)
require(ggnetwork)
require(ergm)
require(intergraph)
require(RColorBrewer)
library('igraph')
library(reshape)



############################
sp1929 <-read.csv("Smith_1929_species_wide.csv",head=T)
sp2007 <-read.csv("Bickford_2007_species_wide_NA_removed.csv",head=T) 
sp2010 <-read.csv("Bickford_2010_species_wide_NA_removed.csv",head=T)
sp2016 <-read.csv("Karlsson_2016_species_wide_KED_removed.csv",head=T)
allyearzone <-read.csv("network_yearzone.csv", head = T)
allyearzone <- allyearzone[,-1]
allyearzone1 <- melt(allyearzone, id=c("Yearzone"))
allyearzone2 <- cast(allyearzone1, variable~Yearzone, sum)
colnames(allyearzone2)[1] <- "Species"
allyearzone2 <- as.data.frame(allyearzone2)


#set up datasets as matrix
g29 <- as.matrix(sp1929[,-c(1,2,31,32)]); row.names(g29) = sp1929$Site
g29[g29 > 0.001] <- 1
i29 <- graph.incidence(g29, mode=c('all') )
#to get node values
V(i29)
#set colour of nodes, first the group nodes
V(i29)$color[2:8] <- rgb(1,0,0,.5)
V(i29)$color[9:36] <- rgb(0,1,0,.5)
V(i29)$label <- V(i29)$name
V(i29)$label.color <- rgb(0,0,.2,.5)
V(i29)$label.cex <- .4
V(i29)$size <- 6
V(i29)$frame.color <- NA
E(i29)$color <- rgb(.5,.5,0,.2)
pdf('i29.pdf')
plot(i29, layout=layout.fruchterman.reingold)
dev.off()

g7 <- as.matrix(sp2007[,-c(1,2,29,30,31,32)]); row.names(g7) = sp2007$Site
g7[g7 > 0.001] <- 1
i7 <- graph.incidence(g7, mode=c('all') )
V(i7)
#set colour of nodes, first the group nodes
V(i7)$color[2:14] <- rgb(1,0,0,.5)
V(i7)$color[15:40] <- rgb(0,1,0,.5)
V(i7)$label <- V(i7)$name
V(i7)$label.color <- rgb(0,0,.2,.5)
V(i7)$label.cex <- .4
V(i7)$size <- 6
V(i7)$frame.color <- NA
E(i7)$color <- rgb(.5,.5,0,.2)
pdf('i7.pdf')
plot(i7, layout=layout.fruchterman.reingold)
dev.off()

g10 <- as.matrix(sp2010[,-c(1,2,35,36,37,38,39)]); row.names(g10) = sp2010$Site
g10[g10 > 0.001] <- 1
i10 <- graph.incidence(g10, mode=c('all') )
V(i10)
#set colour of nodes, first the group nodes
V(i10)$color[1:18] <- rgb(1,0,0,.5)
V(i10)$color[19:58] <- rgb(0,1,0,.5)
V(i10)$label <- V(i10)$name
V(i10)$label.color <- rgb(0,0,.2,.5)
V(i10)$label.cex <- 0.6
V(i10)$size <- 6
V(i10)$frame.color <- NA
E(i10)$color <- rgb(.5,.5,0,.2)
pdf('i10.pdf')
plot(i10, layout=layout.fruchterman.reingold)
dev.off()

g16 <- as.matrix(sp2016[,-c(1,2,31,32,33,34)]); row.names(g16) = sp2016$Site
g16[g16 > 0.001] <- 1
i16 <- graph.incidence(g16, mode=c('all') )
V(i16)
#set colour of nodes, first the group nodes
V(i16)$color[1:27] <- rgb(1,0,0,.5)
V(i16)$color[28:57] <- rgb(0,1,0,.5)
V(i16)$label <- V(i16)$name
V(i16)$label.color <- rgb(0,0,.2,.5)
V(i16)$label.cex <- .4
V(i16)$size <- 6
V(i16)$frame.color <- NA
E(i16)$color <- rgb(.5,.5,0,.2)
pdf('i16.pdf')
plot(i16, layout=layout.fruchterman.reingold)
#plot(i16, layout=layout.kamada.kawai)
#plot(i16, layout=layout.fruchterman.reingold.grid)
dev.off()

#Yearzone
gzone <- as.matrix(allyearzone2[,-c(1)]); row.names(gzone) = allyearzone2$Species
gzone[gzone > 0.001] <- 1
izone <- graph.incidence(gzone, mode=c('all') )
V(izone)
#set colour of nodes, first the group nodes
V(izone)$color[2:74] <- rgb(1,0,0,.5)
V(izone)$color[75:86] <- rgb(0,1,0,.5)
V(izone)$label <- V(izone)$name
V(izone)$label.color <- rgb(0,0,.2,.5)
V(izone)$label.cex <- .4
V(izone)$size <- 6
V(izone)$frame.color <- NA
E(izone)$color <- rgb(.5,.5,0,.2)
pdf('izone.pdf')
plot(izone, layout=layout.fruchterman.reingold)
#plot(i16, layout=layout.kamada.kawai)
#plot(i16, layout=layout.fruchterman.reingold.grid)
dev.off()
izone$layout <- layout.fruchterman.reingold(izone)
V(izone)$label <- V(izone)$name
tkplot(izone)
izone$layout <- tkplot.getcoords(17)
pdf("izonecorrected_3.pdf")
plot(izone)
dev.off()

adjm <- matrix(sample(0:1, 100, replace=TRUE, prob=c(0.9,0.1)), nc=10)

g8 <- graph_from_adjacency_matrix(gzone, weighted=TRUE, mode="plus")
d2 <- function(x) { diag(x) <- diag(x)/2; x }
all(nzs((d2(adjm+t(adjm)))[lower.tri(adjm)]) == sort(E(g8)$weight))

g9 <- graph_from_adjacency_matrix(gzone, weighted=TRUE, mode="plus", diag=FALSE)
d0 <- function(x) { diag(x) <- 0 }
all(nzs((d0(adjm+t(adjm)))[lower.tri(adjm)]) == sort(E(g9)$weight))

#We???ve emphasized groups in this visualization so much, that we might want to just create a network consisting of group co-membership. First we need to create a new network object. We???ll do that the same way for this network as for our example at the top of this page:
g29e <- t(g29) %*% g29
g7e <- t(g7) %*% g7
g10e <- t(g10) %*% g10
g16e <- t(g16) %*% g16
gzonee <- t(gzone) %*% gzone

i29e <- graph.adjacency(g29e, mode = "undirected")
i7e <- graph.adjacency(g7e, mode = "undirected")
i10e <- graph.adjacency(g10e, mode = "undirected")
i16e <- graph.adjacency(g16e, mode = "undirected")
izonee <- graph.adjacency(gzonee, mode = "undirected")

#Now we need to tansform the graph so that multiple edges become an attribute ( E(g)$weight ) of each unique edge:
E(i29e)$weight <- count.multiple(i29e)
i29e <- simplify(i29e)

E(i7e)$weight <- count.multiple(i7e)
i7e <- simplify(i7e)

E(i10e)$weight <- count.multiple(i10e)
i10e <- simplify(i10e)

E(i16e)$weight <- count.multiple(i16e)
i16e <- simplify(i16e)

E(izonee)$weight <- count.multiple(izonee)
izonee <- simplify(izonee)

#Now we???ll set the other plotting parameters as we did above:
# Set vertex attributes
V(i29e)$label <- V(i29e)$name
V(i29e)$label.color <- rgb(0,0,.2,.8)
V(i29e)$label.cex <- .6
V(i29e)$size <- 6
V(i29e)$frame.color <- NA
V(i29e)$color <- rgb(0,0,1,.5)

V(i7e)$label <- V(i7e)$name
V(i7e)$label.color <- rgb(0,0,.2,.8)
V(i7e)$label.cex <- .6
V(i7e)$size <- 6
V(i7e)$frame.color <- NA
V(i7e)$color <- rgb(0,0,1,.5)

V(i10e)$label <- V(i10e)$name
V(i10e)$label.color <- rgb(0,0,.2,.8)
V(i10e)$label.cex <- .6
V(i10e)$size <- 6
V(i10e)$frame.color <- NA
V(i10e)$color <- rgb(0,0,1,.5)

V(i16e)$label <- V(i16e)$name
V(i16e)$label.color <- rgb(0,0,.2,.8)
V(i16e)$label.cex <- .6
V(i16e)$size <- 6
V(i16e)$frame.color <- NA
V(i16e)$color <- rgb(0,0,1,.5)

V(izonee)$label <- V(izonee)$name
V(izonee)$label.color <- rgb(0,0,.2,.8)
V(izonee)$label.cex <- .6
V(izonee)$size <- 6
V(izonee)$frame.color <- NA
V(izonee)$color <- rgb(0,0,1,.5)


# Set edge gamma according to edge weight
egam29 <- (log(E(i29e)$weight)+.3)/max(log(E(i29e)$weight)+.3)
E(i29e)$color <- rgb(.5,.5,0,egam29)

egam7 <- (log(E(i7e)$weight)+.3)/max(log(E(i7e)$weight)+.3)
E(i7e)$color <- rgb(.5,.5,0,egam7)

egam10 <- (log(E(i10e)$weight)+.3)/max(log(E(i10e)$weight)+.3)
E(i10e)$color <- rgb(.5,.5,0,egam10)

egam16 <- (log(E(i16e)$weight)+.3)/max(log(E(i16e)$weight)+.3)
E(i16e)$color <- rgb(.5,.5,0,egam16)

egamzone <- (log(E(izonee)$weight)+.3)/max(log(E(izonee)$weight)+.3)
E(izonee)$color <- rgb(.5,.5,0,egamzone)

izonee$layout <- layout.kamada.kawai(izonee)
V(izonee)$label <- V(izonee)$name
tkplot(izonee)
izonee$layout <- tkplot.getcoords(18)
pdf("izonecorrected2.pdf")
plot(izonee)
dev.off()

#We set edge gamma as a function of how many edges exist between two nodes, or in this case, how many species each site has in common. For illustrative purposes, let???s compare how the Kamada-Kawai and Fruchterman-Reingold algorithms render this graph:
pdf("i29e.pdf")
plot(i29e, main = "layout.kamada.kawai", layout=layout.kamada.kawai)
plot(i29e, main = "layout.fruchterman.reingold", layout=layout.fruchterman.reingold)
dev.off()

pdf("i7e.pdf")
plot(i7e, main = "layout.kamada.kawai", layout=layout.kamada.kawai)
plot(i7e, main = "layout.fruchterman.reingold", layout=layout.fruchterman.reingold)
dev.off()

pdf("i10e.pdf")
plot(i10e, main = "layout.kamada.kawai", layout=layout.kamada.kawai)
plot(i10e, main = "layout.fruchterman.reingold", layout=layout.fruchterman.reingold)
dev.off()


pdf("i16e.pdf")
plot(i16e, main = "layout.kamada.kawai", layout=layout.kamada.kawai)
plot(i16e, main = "layout.fruchterman.reingold", layout=layout.fruchterman.reingold)
dev.off()

pdf("izonee.pdf")
plot(izonee, main = "layout.kamada.kawai", layout=layout.kamada.kawai)
plot(izonee, main = "layout.fruchterman.reingold", layout=layout.fruchterman.reingold)
dev.off()


#########asess percent overlap of groups
#irst we???ll need to create a percent overlap graph. We start by dividing each row by the diagonal (this is really easy in R):
ol29 <- g29e/diag(g29e)
ol7 <- g7e/diag(g7e)
ol10 <- g10e/diag(g10e)
ol16 <- g16e/diag(g16e)
olzone <- gzonee/diag(gzonee)

#Next, sum the matricies and set any NA cells (caused by dividing by zero in the step above) to zero:
#magall <- ol29 + ol7 + ol10 + ol16, only works if array looks the same
magall29 <-ol29 
magall7 <-ol7 
magall10 <-ol10
magall16 <- ol16 
magallzone <- olzone 
magall29[is.na(magall29)] <- 0
magall7[is.na(magall7)] <- 0
magall10[is.na(magall10)] <- 0
magall16[is.na(magall16)] <- 0
magallzone[is.na(magallzone)] <- 0

#Note that magall now consists of a percent overlap matrix, if it is summed over years, the maximun is now the number of years.
#Let???s compute average number of a species, by taking the mean across each value in each diagonal:
magdiag29 <- apply(cbind(diag(g29e)), 1, mean )
magdiag7 <- apply(cbind(diag(g7e)), 1, mean )
magdiag10 <- apply(cbind(diag(g10e)), 1, mean )
magdiag16 <- apply(cbind(diag(g16e)), 1, mean )
magdiagzone <- apply(cbind(diag(gzonee)), 1, mean )

#Finally, we???ll generate centrality measures for magall. When we create the igraph object from our matrix, we need to set weighted=T because otherwise igraph dichotomizes edges at 1. This can distort our centrality measures because now edges represent  more than binary connections???they represent the percent of membership overlap.
magallg29 <- graph.adjacency(magall29, weighted=T)
magallg7 <- graph.adjacency(magall7, weighted=T)
magallg10 <- graph.adjacency(magall10, weighted=T)
magallg16 <- graph.adjacency(magall16, weighted=T)
magallgzone <- graph.adjacency(magallzone, weighted=T)

# Degree
V(magallg29)$degree <- degree(magallg29)
V(magallg7)$degree <- degree(magallg7)
V(magallg10)$degree <- degree(magallg10)
V(magallg16)$degree <- degree(magallg16)
V(magallgzone)$degree <- degree(magallgzone)

# Betweenness centrality
V(magallg29)$btwcnt <- betweenness(magallg29)
V(magallg7)$btwcnt <- betweenness(magallg7)
V(magallg10)$btwcnt <- betweenness(magallg10)
V(magallg16)$btwcnt <- betweenness(magallg16)
V(magallgzone)$btwcnt <- betweenness(magallgzone)

#Before we plot this, we should probably filter some of the edges, otherwise our graph will probably be too busy to make sense of visually.  Take a look at the distribution of connection strength by plotting the density of the magall matrix:
plot(density(magall29))
plot(density(magall7))
plot(density(magall10))
plot(density(magall16))
plot(density(magallzone))

#Assess the edge weights and the percent overlap for most species (most below one is an overlap of less than 1/3). If you then filter at 1, an edge will consists of group overlap of more than 1/3 of the group???s members in question.
magallgt129 <- magall29
magallgt129[magallgt129<10] <- 0
magallggt129 <- graph.adjacency(magallgt129, weighted=T)

magallgt17 <- magall7
magallgt17[magallgt17<0.2] <- 0
magallggt17 <- graph.adjacency(magallgt17, weighted=T)

magallgt110 <- magall10
magallgt110[magallgt110<0.2] <- 0
magallggt110 <- graph.adjacency(magallgt110, weighted=T)

magallgt116 <- magall16
magallgt116[magallgt116<0.5] <- 0
magallggt116 <- graph.adjacency(magallgt116, weighted=T)

magallgt1zone <- magallzone
magallgt1zone[magallgt1zone<0.8] <- 0
magallggt1zone <- graph.adjacency(magallgt1zone, weighted=T)

# Removes loops:
magallggt129 <- simplify(magallggt129, remove.multiple=FALSE, remove.loops=TRUE)
magallggt17 <- simplify(magallggt17, remove.multiple=FALSE, remove.loops=TRUE)
magallggt110 <- simplify(magallggt110, remove.multiple=FALSE, remove.loops=TRUE)
magallggt116 <- simplify(magallggt116, remove.multiple=FALSE, remove.loops=TRUE)
magallggt1zone <- simplify(magallggt1zone, remove.multiple=FALSE, remove.loops=TRUE)

#Before we do anything else, we???ll create a custom layout based on Fruchterman.-Ringold wherein we adjust the coordates by hand using the tkplot gui tool to make sure all of the labels are visible. This is very useful if you want to create a really sharp-looking network visualization for publication.
#Let the plot load, then maximize the window, and select to View -> Fit to Screen so that you get maximum resolution for this large graph. Now hand-place the nodes, making sure no labels overlap: Use the coordinate number for the number tk object that you have called (which can be found in the top right hand corner of the graph, i.e. graph object 7 should be saved as "7" not 1. )
magallggt129$layout <- layout.fruchterman.reingold(magallggt129)
V(magallggt129)$label <- V(magallggt129)$name
tkplot(magallggt129)
magallggt129$layout <- tkplot.getcoords(1)

magallggt17$layout <- layout.fruchterman.reingold(magallggt17)
V(magallggt17)$label <- V(magallggt17)$name
tkplot(magallggt17)
magallggt17$layout <- tkplot.getcoords(9)

magallggt110$layout <- layout.fruchterman.reingold(magallggt110)
V(magallggt110)$label <- V(magallggt110)$name
tkplot(magallggt110)
magallggt110$layout <- tkplot.getcoords(10)

magallggt116$layout <- layout.fruchterman.reingold(magallggt116)
V(magallggt116)$label <- V(magallggt116)$name
tkplot(magallggt116)
magallggt110$layout <- tkplot.getcoords(11)


magallggt1zone$layout <- layout.fruchterman.reingold(magallggt1zone)
V(magallggt1zone)$label <- V(magallggt1zone)$name
tkplot(magallggt1zone)
magallggt1zone$layout <- tkplot.getcoords(16)

#Let the plot load, then maximize the window, and select to View -> Fit to Screen so that you get maximum resolution for this large graph. Now hand-place the nodes, making sure no labels overlap: Use the coordinate number for the number tk object that you have called (which can be found in the top right hand corner of the graph, i.e. graph object 7 should be saved as "7" not 1. )

# Set vertex attributes
V(magallggt129)$label <- V(magallggt129)$name
V(magallggt129)$label.color <- rgb(0,0,.2,.6)
V(magallggt129)$size <- 6
V(magallggt129)$frame.color <- NA
V(magallggt129)$color <- rgb(0,0,1,.5)

V(magallggt17)$label <- V(magallggt17)$name
V(magallggt17)$label.color <- rgb(0,0,.2,.6)
V(magallggt17)$size <- 6
V(magallggt17)$frame.color <- NA
V(magallggt17)$color <- rgb(0,0,1,.5)

V(magallggt110)$label <- V(magallggt110)$name
V(magallggt110)$label.color <- rgb(0,0,.2,.6)
V(magallggt110)$size <- 6
V(magallggt110)$frame.color <- NA
V(magallggt110)$color <- rgb(0,0,1,.5)

V(magallggt116)$label <- V(magallggt116)$name
V(magallggt116)$label.color <- rgb(0,0,.2,.6)
V(magallggt116)$size <- 6
V(magallggt116)$frame.color <- NA
V(magallggt116)$color <- rgb(0,0,1,.5)

V(magallggt1zone)$label <- V(magallggt1zone)$name
V(magallggt1zone)$label.color <- rgb(0,0,.2,.6)
V(magallggt1zone)$size <- 6
V(magallggt1zone)$frame.color <- NA
V(magallggt1zone)$color <- rgb(0,0,1,.5)

#   One thing that we can do with this graph is to set label size as a function of degree, which adds a ???tag-cloud???-like element to the visualization:

V(magallggt129)$label.cex <- V(magallggt129)$degree/(max(V(magallggt129)$degree)/2)+ .3
V(magallggt17)$label.cex <- V(magallggt17)$degree/(max(V(magallggt17)$degree)/2)+ .3
V(magallggt110)$label.cex <- V(magallggt110)$degree/(max(V(magallggt110)$degree)/2)+ .3
V(magallggt116)$label.cex <- V(magallggt116)$degree/(max(V(magallggt116)$degree)/2)+ .3
#note, unfortunately one must play with the formula above to get the
#ratio just right

pdf("magallggt1customlayout.pdf")
plot(magallggt129)
dev.off()

pdf("magallggt1customlayout.pdf")
plot(magallggt17)
dev.off()

pdf("magallggt1customlayout.pdf")
plot(magallggt110)
dev.off()

pdf("magallggt1customlayout.pdf")
plot(magallggt116)
dev.off()

pdf("magallggt1customlayout.pdf")
plot(magallggt1zone)
dev.off()






