############################################################################.
############  Stepping-Stone approach to larval connectivity  ##############.
############################################################################.

############  Script by Emilie Boulanger and Marco Andrello  ###############.
############  for (Ecography ref)                            ###############.

# do ctrl + shift + F10 to restart R session

# load libraries
library(igraph)
library(maps)
library(mapdata)
library(Rfast)   # floyd algorithm
library(fields)  # for colour ramps

# load data
# coordinates of 100 grid cells
coord100 <- read.table("data/coord_mullus_100gridCells.txt")

# larval dispersal probability matrix
dispersal100 <- read.table("data/larvalDispersal_mullus_100gridCells.txt") # -> the C100 matrix
colnames(dispersal100) <- rownames(dispersal100) # 

#  47 sampled population sites
pool_frequencies <- read.table("data/popFreq_mullus_47sites_subset.txt") # this is a subset of the pool frequencies datafile available on dryad
# get the cell names corresponding to the 100 grid cells
pop_cell <- sub(".*?p(.*?)f.*", "\\1", rownames(pool_frequencies)) 
pop_cell <- sub(".*?_.*", "\\1", pop_cell)

##########################  plot spatial graphs  ##########################
############################## .. 100 SITES ###############################

# remove self-recruitment
diag(dispersal100) <- 0.0
# prepare data
dispersal100 <- as.matrix(dispersal100)
# simple spatial graph
g100 <-graph.adjacency(dispersal100, mode="directed", weighted=TRUE)
length(E(g100))
coord100 <- coord100[match(V(g100)$name, rownames(coord100)), ]

V(g100)$Latitude  <- coord100[,"Latitude"] ## V(g): nodes (=sampling sites) of the graph
V(g100)$Longitude <- coord100[, "Longitude"]

# plot the map with spatial graph
map("worldHires", xlim=c(-8,37), ylim=c(29.5,47), col="gray90", fill=TRUE)
plot(g100, add=TRUE, rescale=FALSE,layout=coord100,
     vertex.size=40,        
     vertex.color="black",  
     edge.color="red",      
     vertex.label=NA)       

# the mediterranean sea is fully connected when we consider the 100 grid cell centroids

########################### .. 47 POPULATIONS #############################

# subset larval dispersal data
dispersal47 <- dispersal100[pop_cell, pop_cell] # -> the C47 matrix
# get the coordinates for the 47 sites
coord47 <- coord100[pop_cell,]
# simple spatial graph
g47 <-graph.adjacency(dispersal47, mode="undirected", weighted=TRUE)
length(E(g47))
coord47 <- coord47[match(V(g47)$name, rownames(coord47)), ]

V(g47)$Latitude  <- coord47[,"Latitude"] ## V(g): nodes (=sampling sites) of the graph
V(g47)$Longitude <- coord47[, "Longitude"]

# plot the map with spatial graph
map("worldHires", xlim=c(-8,37), ylim=c(29.5,47), col="gray90", fill=TRUE)
plot(g47, add=TRUE, rescale=FALSE,layout=coord47,
     vertex.size=40,        
     vertex.color="black",  
     edge.color="red",      
     vertex.label=NA) 

# the 47 sampled popualtion sites are no fully connected but form three clusters

#################  calculate the shortest-path distance  ##################

####################  .. Floyd-Warshal shortest-paths  ####################
# to correlate with genetic distances

# first transform larval dispersal probabilities to larval dispersal distances
distance100 <- log(1/dispersal100)
diag(distance100) <- 0
# apply the Floyd-Warshall algorithm for shortest paths in a directed graph
sp.fw.100 <- floyd(distance100)
rownames(sp.fw.100) <- rownames(distance100)
colnames(sp.fw.100) <- colnames(distance100)
# subset 47 sampled population sites
sp.fw.47 <- sp.fw.100[pop_cell, pop_cell]
# plot the spatial graph
gsp47 <-graph.adjacency(sp.fw.47, mode="undirected", weighted=TRUE)
length(E(gsp47))
coord47 <- coord47[match(V(gsp47)$name, rownames(coord47)), ]

V(gsp47)$Latitude  <- coord47[,"Latitude"] ## V(g): nodes (=sampling sites) of the graph
V(gsp47)$Longitude <- coord47[, "Longitude"]

# plot the map with spatial graph
map("worldHires", xlim=c(-8,37), ylim=c(29.5,47), col="gray90", fill=TRUE)
plot(gsp47, add=TRUE, rescale=FALSE,layout=coord47,
     vertex.size=40,        
     vertex.color="black",  
     edge.color="red",      
     vertex.label=NA) 
# the larval distance matrix is now filled up 
# and all sites are connected through shortest paths representing multiple generations
# the sp.fw.47 is then used to compare multi-generational propagule dispersal to multi-generational genetic connectivity

###################  .. Explicit shortest path lengths  ###################

# for visualization & to infer temporal scale of connectivity

# create a spatial graph based on the larval distance values between the 100 grid cells
gdist100 <- graph.adjacency(distance100, mode = "directed", weighted = TRUE)

# get specific shortest path from one source
get.shortest.paths(gdist100, from = 5, to=100, output = "both")
   # $vpath
   # + 4/100 vertices:
   #   [1] 5   65  99  100
   # $epath
   # + 3/9900 edges:
   #   [1] 5 ->65  65->99  99->100
get.shortest.paths(gdist100, from = 100, to=5, output = "both")
   # $vpath
   # + 3/100 vertices:
   #   [1] 100 3   5  
   # $epath
   # + 2/9900 edges:
   #   [1] 100->3 3  ->5

# loop this: get all specific shortest paths from each source node to the whole graph
x <- 1
mylist <- list()
for (i in 1:length(V(gdist100))) {
  sp <- shortest_paths(gdist100, from = i, to = V(gdist100), output = "both")
  # extract the number of shortest paths and put it in vector x
  for (j in 1:length(sp$epath)){
    x[j] <- length(sp$epath[[j]])
    }
  # fill list with generated vectors
  mylist[[i]] <- x
  }
# matrix of the explicit shortest path lengths between all pairs of nodes
sp_matrix <- do.call("rbind",mylist)
rownames(sp_matrix) <- V(gdist100)
colnames(sp_matrix) <- V(gdist100)
# verify if correct
shortest_paths(gdist100, from = 20, to=74)
   # 24/100 vertices: -> this equals to 23 edges connecting teh 24 vertices or nodes
   # [1] 20 4  5  3  45 46 49 21 53 52 51 25 14 75 21 32 76 36 35 40 54 56 73 74
sp_matrix[20,74]
   # 23
summary(as.vector(sp_matrix))

## subset 47 sampled population sites
sp_matrix_47 <- sp_matrix[pop_cell, pop_cell]
summary(as.vector(sp_matrix_47))

################### .... represent in spatial graph map: colour levels ----
# to recreate MS Figure 3a

# Create matrices lon and lat
lon <- coord47$Longitude
lat <- coord47$Latitude

idmat		<- matrix(data=1,nrow=1,ncol=ncol(sp_matrix_47))

lonmat	<- matrix(data=lon,nrow=nrow(sp_matrix_47),ncol=1)
lon_dest	<- lonmat%*%idmat
lon_or	<- t(lon_dest)

latmat	<- matrix(data=lat,nrow=nrow(sp_matrix_47),ncol=1)
lat_dest	<- latmat%*%idmat
lat_or	<- t(lat_dest)

# Read data connectivity
c <- sp_matrix_47                       

# Finds nonzero edges and sort them
a 		<- sort(c, decreasing = T, index.return=T) #decreasing because work with distances -> want smallest distances on top when plotting
l0 		<- length(which(a$x==0))
c_strength 	<- a$x[1:c(length(a$x)-(l0+1))] # to include all but 47 0's that are as last
c_id 		<- a$ix[(l0+1):length(a$ix)]

# Find lon and lat of segment origin
lon_1 <- lon_or[c_id]
lat_1 <- lat_or[c_id]

# Find lon and lat of segment destination
lon_2 <- lon_dest[c_id]
lat_2 <- lat_dest[c_id]

# Colourblind friendly spatial graph 
# define colour bins
# bins = seq(0,max(c),1)       # if want one colour for each value
bins = c(seq(0,7),10,15,25)  # for a map that is more easily readable
colours  <- cut(c_strength, breaks=bins,include.lowest=T)
table(colours)
# define colour palette
    # green to brown
colfunc <- colorRampPalette(c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))
col_palette <- colfunc(length(bins))
col_palette <- col_palette[length(col_palette):1] # invert colour palette
palette(col_palette)

# represent the spatial graph
map("worldHires", xlim=c(-8,37), ylim=c(29.5,47),col = "gray85", boundary = TRUE, interior = FALSE, fill = TRUE, border = NA)
points(lon, lat, pch=19, col="black", cex=1)
segments(lon_1,lat_1,lon_2,lat_2,lty=1,lwd=2.5,col=colours)

legend("bottomleft",
       legend = c(1:7,"8-10","11-15","16-25"),
       pch = 19,
       col = col_palette,
       ncol = 3,
       cex = 1.4,
       pt.cex = 2.2,
       bty = "n")
map.axes(cex.axis=1)
map.scale(3, 31, ratio=FALSE, relwidth=0.15, cex=1)

################### highlight specific path in network ####################
# to recreate Supplementary Material Figure A3
# adapted from http://kateto.net/network-visualization

# graph of direct larval distances
distance100noInf <- distance100
distance100noInf[is.infinite(distance100noInf)] <- 0

gt100 <-graph.adjacency(distance100noInf, mode="directed", weighted=TRUE)

# visualize specific path 
path_91_74 <- shortest_paths(gt100, 
                             from = V(gt100)[91], 
                             to  = V(gt100)[74],
                             output = "both") # both path nodes and edges
shortest.paths(gt100, V(gt100)[91], to=V(gt100)[74])
get.shortest.paths(gt100, V(gt100)[91], to=V(gt100)[74], output = "both")

# Generate edge color variable to plot the path:
ecol <- rep("gray80", ecount(g100))
ecol[unlist(path_91_74$epath)] <- "orange"
# Generate edge width variable to plot the path:
ew <- rep(2, ecount(g100))
ew[unlist(path_91_74$epath)] <- 10
# Generate node color variable to plot the path:
vcol <- rep("gray40", vcount(g100))
vcol[unlist(path_91_74$vpath)] <- "gold"
# generate node weight variable to plot the path:
vw <- rep(20, vcount(g100))
vw[unlist(path_91_74$vpath)] <- 40

# plot the path on the spatial graph
#jpeg(file="SpatialGraph_LarvalConn_100_path_91-74_47coords.jpeg", width = 1700, height = 900, units = "px")
map("worldHires", xlim=c(-8,37), ylim=c(29.5,47), col="gray90", fill=TRUE)
plot(gt100, add = TRUE, rescale= FALSE, layout = coord100,
     vertex.color=vcol, edge.color=ecol, 
     edge.width=ew, edge.arrow.mode=0,
     vertex.size = vw,
     vertex.label = NA)
# add 47 sampling points to show interest of going outside sampled populations to calculate path
points(coord47$Longitude, coord47$Latitude, pch=23, col = "black", bg = "black", cex=1.5)
#map.scale(-3.5, 30.5, metric = TRUE, ratio = FALSE, cex = 2)
map.axes(cex.axis=2)
legend("bottomleft",
       legend = c("100 grid cell centroids",
                  "47 sampling sites", 
                  "shortest path stepping-stones"),
       pch = c(19, 18, 21),
       col = c("black", "black", "black"),
       pt.bg = c("black", "black", "gold"),
       cex = 2)
#dev.off()

 