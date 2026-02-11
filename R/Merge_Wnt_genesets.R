#' Removal of gene sets containing few genes after purification, and integration of biologically similar Wnt-related gene sets through construction of Jaccard index networks.
#' 
#' @title Merge_Wnt_genesets
#' @param file_paths Working directory path.
#' @param activation_geneset Wnt pathway activation gene sets to be merged. (Place in the working directory)
#' @param inhibition_geneset Wnt pathway inhibition gene sets to be merged. (Place in the working directory)
#' @param delete_GN TRUE: Remove gene sets containing insufficient numbers of Wnt-related genes. FALSE: Retain all gene sets. Default parameter is "TRUE".
#' @param min_GN Gene sets with fewer genes than min_GN will be removed from the Wnt pathway collection. Default parameter is 5.
#' @param similarity When the Jaccard index similarity between two Wnt-related gene sets exceeds this threshold, the gene sets will be integrated. Default parameter is 0.3.
#' @author Dingkang Zhao
#' @examples
#' @return
#' @export
Merge_Wnt_genesets=function(file_paths,
                        activation_geneset,
                        inhibition_geneset,
                        delete_GN=T,
                        delete_time=c("before","after"),
                        min_GN=5,
                        integration_method=c("Jaccard","Sorensen-Dice","Hub-Promoted","Hub-Depressed"),
                        de_redundant_basis=c("igraph_component","agglomeration_ward.D",
                                             "agglomeration_ward.D2","agglomeration_single",
                                             "agglomeration_complete","agglomeration_average",
                                             "agglomeration_mcquitty","agglomeration_median",
                                             "agglomeration_centroid"),
                        similarity=0.3,
                        export_file=F
                          ){
  setwd(file_paths)
  library(reshape2)
  library(igraph)
  GENE_SETS_UP=activation_geneset[,colSums(activation_geneset==""|is.na(activation_geneset)|activation_geneset==" ")<nrow(activation_geneset)]
  GENE_SETS_DN=inhibition_geneset[,colSums(inhibition_geneset==""|is.na(inhibition_geneset)|inhibition_geneset==" ")<nrow(inhibition_geneset)]
  update_GENE_SETS_UP=list()
  update_GENE_SETS_DN=list()
  for (i in 1:length(colnames(GENE_SETS_UP))){
    colup=GENE_SETS_UP[,i]
    colup=colup[colup != ""]
    update_GENE_SETS_UP[[i]]=colup
    names(update_GENE_SETS_UP)[i]=paste0(colnames(GENE_SETS_UP)[i])
  }
  for (i in 1:length(colnames(GENE_SETS_DN))){
    coldn=GENE_SETS_DN[,i]
    coldn=coldn[coldn != ""]
    update_GENE_SETS_DN[[i]]=coldn
    names(update_GENE_SETS_DN)[i]=paste0(colnames(GENE_SETS_DN)[i])
  }
  if(delete_time=="before"){
    if(delete_GN){
      update_GENE_SETS_UP=update_GENE_SETS_UP[sapply(update_GENE_SETS_UP,length)>=min_GN]
      update_GENE_SETS_DN=update_GENE_SETS_DN[sapply(update_GENE_SETS_DN,length)>=min_GN]
    }else{
      print("No genesets are deleted.")
    }
  }
  print("Start calculating correlations between genesets using integration coefficients.")
  if(de_redundant_basis=="igraph_component"){
    if(integration_method=="Jaccard"){
      integration_ratio=function(set1, set2) {
        intersection=length(intersect(set1, set2))
        union=length(union(set1, set2))
        ratio=intersection/union
        return(ratio)
      }
    }else if(integration_method=="Sorensen-Dice"){
      integration_ratio=function(set1, set2) {
        intersection=length(intersect(set1, set2))
        Nset1=length(set1)
        Nset2=length(set2)
        ratio=(intersection*2)/(Nset1+Nset2)
        return(ratio)
      }
    }else if(integration_method=="Hub-Promoted"){
      integration_ratio=function(set1, set2) {
        intersection=length(intersect(set1, set2))
        Nset1=length(set1)
        Nset2=length(set2)
        ratio=intersection/min(Nset1,Nset2)
        return(ratio)
      }
    }else if(integration_method=="Hub-Depressed"){
      integration_ratio=function(set1, set2) {
        intersection=length(intersect(set1, set2))
        Nset1=length(set1)
        Nset2=length(set2)
        ratio=intersection/max(Nset1,Nset2)
        return(ratio)
      }
    }
    matrix_UP=outer(
      names(update_GENE_SETS_UP), names(update_GENE_SETS_UP), 
      Vectorize(function(x, y) integration_ratio(update_GENE_SETS_UP[[x]], update_GENE_SETS_UP[[y]]))
    )
    dimnames(matrix_UP)=list(names(update_GENE_SETS_UP), names(update_GENE_SETS_UP))
    matrix_DN=outer(
      names(update_GENE_SETS_DN), names(update_GENE_SETS_DN), 
      Vectorize(function(x, y) integration_ratio(update_GENE_SETS_DN[[x]], update_GENE_SETS_DN[[y]]))
    )
    dimnames(matrix_DN)=list(names(update_GENE_SETS_DN), names(update_GENE_SETS_DN))
    matrix_UP[upper.tri(matrix_UP)]=NA
    matrix_DN[upper.tri(matrix_DN)]=NA
    matrix_UP=reshape2::melt(as.matrix(matrix_UP))
    colnames(matrix_UP)=c("rownames", "colnames", "value")
    matrix_DN=reshape2::melt(as.matrix(matrix_DN))
    colnames(matrix_DN)=c("rownames", "colnames", "value")
    matrix_UP=matrix_UP[matrix_UP$rownames != matrix_UP$colnames, ]
    matrix_DN=matrix_DN[matrix_DN$rownames != matrix_DN$colnames, ]
    matrix_UP=na.omit(matrix_UP)
    matrix_DN=na.omit(matrix_DN)
    matrix_UP0=matrix_UP[matrix_UP$value>=similarity,]
    matrix_DN0=matrix_DN[matrix_DN$value>=similarity,]
    if (is.null(matrix_UP0)) {
      print("Low correlation between genesets, no need for joint genesets(UP).")
    } else {
      g=graph_from_data_frame(matrix_UP0, directed = FALSE)
      clusters_UP=components(g)
      membership_UP=as.data.frame(clusters_UP$membership)
      colnames(membership_UP)="cluster"
      num_cluster=unique(membership_UP$cluster)
      for(z in 1:length(num_cluster)){
        membershipup=membership_UP[membership_UP$cluster==num_cluster[z],,drop=F]
        union_vector=c()
        for (name in rownames(membershipup)) {
          union_vector=union(union_vector, update_GENE_SETS_UP[[name]])
          update_GENE_SETS_UP[[name]]=NULL
        }
        update_GENE_SETS_UP[[length(update_GENE_SETS_UP)+1]]=union_vector
        #names(update_GENE_SETS_UP)[length(update_GENE_SETS_UP)]=paste0("Integration_",num_cluster[z],"_UP.",paste(rownames(membershipup),collapse = "."),".")
        names(update_GENE_SETS_UP)[length(update_GENE_SETS_UP)]=paste0("Integration_",num_cluster[z],"_UP")
      }
    }
    if (is.null(matrix_DN0)) {
      print("Low correlation between genesets, no need for joint genesets(DN).")
    } else {
      g=graph_from_data_frame(matrix_DN0, directed = FALSE)
      clusters_DN=components(g)
      membership_DN=as.data.frame(clusters_DN$membership)
      colnames(membership_DN)="cluster"
      num_cluster=unique(membership_DN$cluster)
      for(z in 1:length(num_cluster)){
        membershipdn=membership_DN[membership_DN$cluster==num_cluster[z],,drop=F]
        union_vector=c()
        for (name in rownames(membershipdn)) {
          union_vector=union(union_vector, update_GENE_SETS_DN[[name]])
          update_GENE_SETS_DN[[name]]=NULL
        }
        update_GENE_SETS_DN[[length(update_GENE_SETS_DN)+1]]=union_vector
        #names(update_GENE_SETS_DN)[length(update_GENE_SETS_DN)]=paste0("Integration_",num_cluster[z],"_DN.",paste(rownames(membershipdn),collapse = "."),".")
        names(update_GENE_SETS_DN)[length(update_GENE_SETS_DN)]=paste0("Integration_",num_cluster[z],"_DN")
      }
    }
  }else if(de_redundant_basis%in%c("agglomeration_ward.D",
                                   "agglomeration_ward.D2","agglomeration_single",
                                   "agglomeration_complete","agglomeration_average",
                                   "agglomeration_mcquitty","agglomeration_median",
                                   "agglomeration_centroid")){
    aux.merge = function(gsets, minsize = 5, minjac=0.2){
      if(integration_method=="Jaccard"){
        integration_ratio=function(set1, set2) {
          intersection=length(intersect(set1, set2))
          union=length(union(set1, set2))
          ratio=intersection/union
          return(ratio)
        }
      }else if(integration_method=="Sorensen-Dice"){
        integration_ratio=function(set1, set2) {
          intersection=length(intersect(set1, set2))
          Nset1=length(set1)
          Nset2=length(set2)
          ratio=(intersection*2)/(Nset1+Nset2)
          return(ratio)
        }
      }else if(integration_method=="Hub-Promoted"){
        integration_ratio=function(set1, set2) {
          intersection=length(intersect(set1, set2))
          Nset1=length(set1)
          Nset2=length(set2)
          ratio=intersection/min(Nset1,Nset2)
          return(ratio)
        }
      }else if(integration_method=="Hub-Depressed"){
        integration_ratio=function(set1, set2) {
          intersection=length(intersect(set1, set2))
          Nset1=length(set1)
          Nset2=length(set2)
          ratio=intersection/max(Nset1,Nset2)
          return(ratio)
        }
      }
      integration_matrix = outer(1:length(gsets), 1:length(gsets), 
                             Vectorize(function(x, y) integration_ratio(gsets[[x]], 
                                                                    gsets[[y]])))
      dimnames(integration_matrix) = list(names(gsets),names(gsets))
      
      clust = cutree(hclust(as.dist(1-integration_matrix),method=sub(".*_", "", de_redundant_basis)),h=1-minjac)
      tapply(1:length(gsets),clust,function(i){ sort(unique(unlist(gsets[i]))) })	
    }
    update_GENE_SETS_UP = aux.merge(update_GENE_SETS_UP, minsize = min_GN, minjac=similarity)
    update_GENE_SETS_DN = aux.merge(update_GENE_SETS_DN, minsize = min_GN, minjac=similarity)
  }
  if(delete_time=="after"){
    if(delete_GN){
      update_GENE_SETS_UP=update_GENE_SETS_UP[sapply(update_GENE_SETS_UP,length)>=min_GN]
      update_GENE_SETS_DN=update_GENE_SETS_DN[sapply(update_GENE_SETS_DN,length)>=min_GN]
    }else{
      print("No genesets are deleted.")
    }
  }
  max_length_UP= max(sapply(update_GENE_SETS_UP,length))
  update_GENE_SETS_UP=lapply(update_GENE_SETS_UP, function(x) {
    if (length(x) < max_length_UP) {
      c(x, rep("", max_length_UP - length(x)))
    } else {
      x
    }
  })
  max_length_DN=max(sapply(update_GENE_SETS_DN,length))
  update_GENE_SETS_DN=lapply(update_GENE_SETS_DN, function(x) {
    if (length(x) < max_length_DN) {
      c(x, rep("", max_length_DN - length(x)))
    } else {
      x
    }
  })
  if(de_redundant_basis=="igraph_component"){
    joint_GENE_SETS_UP=as.data.frame(update_GENE_SETS_UP)
    joint_GENE_SETS_DN=as.data.frame(update_GENE_SETS_DN)
    joint_geneset=list(joint_GENE_SETS_UP,joint_GENE_SETS_DN,matrix_UP,matrix_DN)
    names(joint_geneset)=c("joint_activation_geneset","joint_inhibition_geneset","jaccard_network_UP","jaccard_network_DN")
  }else if(de_redundant_basis%in%c("agglomeration_ward.D",
                                   "agglomeration_ward.D2","agglomeration_single",
                                   "agglomeration_complete","agglomeration_average",
                                   "agglomeration_mcquitty","agglomeration_median",
                                   "agglomeration_centroid")){
    joint_GENE_SETS_UP=as.data.frame(update_GENE_SETS_UP)
    joint_GENE_SETS_DN=as.data.frame(update_GENE_SETS_DN)
    joint_geneset=list(joint_GENE_SETS_UP,joint_GENE_SETS_DN)
    names(joint_geneset)=c("joint_activation_geneset","joint_inhibition_geneset")
  }
  if(export_file==T){
  write.table(joint_GENE_SETS_UP,file="joint_activation_geneset.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
  write.table(joint_GENE_SETS_DN,file="joint_inhibition_geneset.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
  }
  print("This calculation is over.")
  return(joint_geneset) 
}










