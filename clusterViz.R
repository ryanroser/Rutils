# Ryan Roser
# ryanroser@gmail.com
# A collection of functions for visualizing clusters

library(igraph)
library(lsa)
library(scales)

clusterMatrix = function(m, minK=5, maxK=20, nstart=5) {
  # This helps pick K for k-means clustering

  out = c()
  for( k in minK:maxK ) {
    model = kmeans(m,centers=k,nstart=nstart)
    out = rbind(out,c(k, model$betweenss, model$totss, mean(model$withinss)))
  }
  
  # see where things level off and pick the number of clusters to use
  cl_thresh = 0.1*max(out[,4]) + min(out[,4])*0.9
  n_clusters = min((minK:maxK)[out[,4] <=  cl_thresh])
  
  plot(data.frame(k=out[,1],
                  `Mean WithinSS`=out[,4]),
       type='l',
       main=sprintf("K-Means Clustering - 90%% Reduction in Average WithinSS at k=%.f",n_clusters))
  points(out[,1],rep(cl_thresh,nrow(out)),type='l',col='red')
  points(n_clusters,out[out[,1]==n_clusters,4],col="red",pch=4)
  
  model = kmeans(m,centers=n_clusters,nstart=nstart)
  return(model)
}

pickClusterColors = function(n_clusters) {
  # Helper function to get nice colors for the clusters
  rainbow(n_clusters)
}

loadGraph = function(edges, clusters, min_edge_value=NA, max_edge_value=NA) {
  # creates a graph object from a matrix of edges
  # provides functionallity to prune edges if necessary
  
  if(!is.na(min_edge_value) | !is.na(max_edge_value)) {
    if( is.na(min_edge_value) ) {
      min_edge_value = 1
    }
    if( is.na(max_edge_value) ) {
      max_edge_value = exp(max(edges$z))+1
    }
    # View the initial distribution
    hist(edges$z, main="Distribution of edge values prior to pruning")
    edge.net = graph.data.frame(edges[edges$z >= log(min_edge_value) & edges$z <= log(max_edge_value),], directed=F)
    hist(E(edge.net)$z, main="Distribution of edge values after pruning")
  } else {
    edge.net = graph.data.frame(edges, directed=F)
  }
  
  # assign some attributes to the verticies
  cluster_colors = pickClusterColors(length(unique(clusters)))
  V(edge.net)$size = log(degree(edge.net))
  V(edge.net)$cluster = clusters[V(edge.net)$name]
  V(edge.net)$color = cluster_colors[V(edge.net)$cluster]
  
  return(edge.net)
}

plotCluster = function(edge.net, cluster, cluster_id, max_cluster_size=50) {
  # plots the nodes in a cluster within a graph object
  
  vnames = V(edge.net)$name
  cluster_names = names(cluster[cluster==cluster_id])
  cluster.vs = V(edge.net)[vnames %in% cluster_names] 
  cluster.net = subgraph(edge.net, cluster.vs)

  # scale the cluster sizes based on the weights of connected edges
  total_wts = sapply(1:length(V(cluster.net)),function(i){
    sum(E(cluster.net)[ from(i) ]$z)
  })
  
  layout = layout.fruchterman.reingold(cluster.net,weight=V(cluster.net)$z)
  #layout = layout.spring(cluster.net,repulse=T)
  #layout = layout.drl(cluster.net,weight=E(cluster.net)$z)

  # nicely calculate the label locations based on the layout
  y = rescale(layout[,2],c(-1,1),range(layout[,2]))
  x = rescale(layout[,1],c(-1,1),range(layout[,1]))
  lab.locs = atan( y / abs(x) )
  lab.locs = -lab.locs
  lab.locs = ifelse(x>0, lab.locs, pi - lab.locs)
  
  plot(cluster.net, layout=layout,
       vertex.size=pmin(max_cluster_size,pmax(total_wts,1)),
       edge.width=E(cluster.net)$z*2,
       vertex.label.dist=0.5,
       vertex.label.color='black',
       vertex.label.degree=lab.locs,
       main=sprintf("Graph for Cluster %s (%i members)",cluster_id, length(V(cluster.net))),
       vertex.label.cex=1)
  return(cluster.net)
}

plotEdges = function(edges, clusters, percentile=0, fontSize=0.8, min_n=1) {
  # creates a graph based on a matrix of edges
  # clusters provides the assignments for the verticies
  
  # make graph
  prunedEdges = pruneObservations(as.matrix(edges), percentile, min_n)
  edge.net = graph.data.frame(prunedEdges, directed=F)
  
  # assign some attributes to the verticies
  n_clusters = length(unique(clusters))
  cluster_colors = pickClusterColors(n_clusters)
  V(edge.net)$size = log(degree(edge.net))
  V(edge.net)$cluster = clusters[V(edge.net)$name]
  V(edge.net)$color = cluster_colors[V(edge.net)$cluster]  
  
  layout = layout.fruchterman.reingold(edge.net,weight=E(edge.net)$z)
  #layout = layout.drl(edge.net,weight=E(edge.net)$z)
  #layout = layout.spring(edge.net,repulse=T)
  
  total_wts = sapply(1:length(V(edge.net)),function(i){
    sum(E(edge.net)[ from(i) ]$z)
  })
  total_wts = 10 * rank(total_wts)/length(total_wts) + 1
  
  # normalize labels
  y = rescale(layout[,2],c(-1,1),range(layout[,2]))
  x = rescale(layout[,1],c(-1,1),range(layout[,1]))
  lab.locs = atan( y / abs(x) )
  lab.locs = -lab.locs
  lab.locs = ifelse(x>0, lab.locs, pi - lab.locs)
  
  plot(edge.net, layout=layout, 
       vertex.size=total_wts,
       edge.width=E(edge.net)$z*2,
       vertex.label.dist=0.1,
       vertex.label.cex=fontSize,
       vertex.label.degree=lab.locs,
       asp=1,
       main="Full Graph")
  legend("topleft",legend=paste("C",1:length(cluster_colors),sep=""),fill=cluster_colors)
  return(edge.net)
}

plotSubsetOfClusters = function(edge.net, mask) {
  # plots a subset of the graph verticies based on mask
  
  bad.vs = V(edge.net)[mask] #the ones to remove
  clean.net = delete.vertices(edge.net, bad.vs) 
  
  #layout = layout.spring(clean.net,repulse=T)
  layout = layout.drl(clean.net,weight=E(clean.net)$z)
  #layout = layout.auto
  
  total_wts = sapply(1:length(V(clean.net)),function(i){
    sum(E(clean.net)[ from(i) ]$z)
  })
  total_wts = 10 * rank(total_wts)/length(total_wts) + 1
  
  plot(clean.net, layout=layout, 
       vertex.size=total_wts,
       edge.width=E(clean.net)$z*2,
       vertex.label.dist=0.5,
       vertex.label.cex=0.5,
       asp=0)
  legend("topleft",legend=paste("Cluster",1:n_clusters),fill=cluster_colors)
  return(clean.net)
}

pruneObservations = function(allObs, percentile, min_n=1) {
  # removes edges from a matrix based on the percentile cutoff
  
  x = allObs[,1]
  y = allObs[,2]
  z = as.numeric(allObs[,3])
  min_z = quantile(as.numeric(z), percentile)
  
  obs = c()
  for( i in 1:nrow(allObs) ) {
    if( !is.na(min_z) & z[i] >= min_z )
      obs = rbind(obs, c(x[i],y[i],z[i]))
  }
  
  clusters_with_obs = sort(unique(c(obs[,1],obs[,2])))
  
  # get all cluster names, ordered by max(z) for name, decreasing
  #maxZ = tapply(c(productEdges$z,productEdges$z),c(productEdges$x,productEdges$y),max,na.rm=T)
  #maxZ = tapply(c(z,z),c(x,y),max,na.rm=T)
  
  all_clusters = sort(unique(c(x,y)))
  for( i in all_clusters ) {
    # there arent any observations, so pick the observation with the greatest z
    if(!(i %in% clusters_with_obs)) {
      i_mask = (x==i | y==i)
      z_order = order(z[i_mask], decreasing=T)
      rng = 1:min(length(z_order),min_n)
      #best_mask = which.max(as.numeric(z[i_mask]))
      best_x = x[i_mask][z_order][rng]
      best_y = y[i_mask][z_order][rng]
      best_z = z[i_mask][z_order][rng]
      obs = rbind(obs,cbind(best_x, best_y, best_z))
      clusters_with_obs = sort(unique(c(clusters_with_obs, best_x, best_y)))
    }
  }
  return(data.frame(x=obs[,1], y=obs[,2], z=as.numeric(obs[,3])))
}

plotCentroids = function(model,percentile=0.9) {
  # make cluster graph for the k-means model centroids

  centroids = t(model$centers)
  n_clusters = ncol(centroids)
  cluster_colors = pickClusterColors(n_clusters)
  
  # assemble edges
  cldata = c()
  cldist = cosine(centroids)
  for(i in 1:(nrow(cldist)-1)){
    for(j in (i+1):nrow(cldist)){
      cldata = rbind(cldata, c(i,j,cldist[i,j]))
    }
  }
  colnames(cldata) = c("x","y","z")
  
  # prune graph
  obs = pruneObservations(cldata,percentile)
  
  centroid.net = graph.data.frame(obs, directed=F)
  V(centroid.net)$size = log(degree(centroid.net))
  V(centroid.net)$cluster = as.numeric(V(centroid.net)$name)
  V(centroid.net)$color = cluster_colors[V(centroid.net)$cluster]
  
  # layout graph
  layout = layout.auto
  
  total_wts = sapply(1:length(V(centroid.net)),function(i){
    sum(E(centroid.net)[ from(i) ]$z)
  })
  
  plot(centroid.net, layout=layout, 
       vertex.size=pmax(5,pmax(2,total_wts)*1.5),
       edge.width= qnorm(E(centroid.net)$z)*3,
       vertex.label.dist=0,
       vertex.label.cex=1,
       asp=0,
       main="Relationships between clusters")
  return(centroid.net)
}

createKmeansReport = function(mtx, edges, maxK=30, filename=NA) {
  # creates a report showing the clustering for a matrix and edges
  # using k-means to generate the clusters
  
  # if filename is not NA, a pdf is created
  if(!is.na(filename)) pdf(filename, width=11, height=8.5, pointsize=10,compress=F)
  
  # create kmeans model
  model = clusterMatrix(mtx,maxK=maxK)
  cluster_memberships = model$cluster
  n_clusters = length(unique(cluster_memberships))
  
  # create a graph using the cluster memberships and edges
  # TODO: edges can be determined from the matrix
  edge.net = loadGraph(edges, cluster_memberships)
  
  plotCentroids(model)
  for(cluster_id in 1:n_clusters) {
    plotCluster(edge.net, cluster_memberships, cluster_id, max_cluster_size=20)
  }
  
  if(!is.na(filename)) dev.off()
}

createGraphReport = function(edges, filename=NA, cutoff=0.9, min_n=1, community_function, ...) {
  # create a report that shows the clustering for a matrix of edges
  # using igraph to create the clusters
  
  # prepare graph
  edge.net = plotEdges(edges, 1, percentile=cutoff, min_n=min_n)
  
  if(!is.na(filename)) pdf(filename, width=11, height=8.5, pointsize=10)
  
  comm = community_function(edge.net, ...)
  cluster_memberships = membership(comm)
  edge.net = plotEdges(edges, cluster_memberships, percentile=cutoff, min_n=min_n)
  table(cluster_memberships)
  
  n_clusters = length(unique(cluster_memberships))
  for(cluster_id in 1:n_clusters) {
    plotCluster(edge.net, cluster_memberships, cluster_id, max_cluster_size=20)
  }
  
  if(!is.na(filename)) dev.off()
}

createGraphReportSVG = function(edges, onefile=F, cutoff=0.9, min_n=1, community_function, ...) {
  # same as createGraphReport, but writes to SVG files
  
  # prepare graph
  edge.net = plotEdges(edges, 1, percentile=cutoff, min_n=min_n)
  
  svg(onefile=onefile, width=11, height=8.5, pointsize=10)
  
  comm = community_function(edge.net, ...)
  cluster_memberships = membership(comm)
  edge.net = plotEdges(edges, cluster_memberships, percentile=cutoff, min_n=min_n)
  table(cluster_memberships)
  
  n_clusters = length(unique(cluster_memberships))
  for(cluster_id in 1:n_clusters) {
    plotCluster(edge.net, cluster_memberships, cluster_id, max_cluster_size=20)
  }
  
  dev.off()
}

