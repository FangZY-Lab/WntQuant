#' The identified Wnt-associated gene sets undergo a precision cleaning procedure based on the Score2 metric.
#' Score2 can also be computed for specific Wnt-related biological contexts to validate the accuracy of external Wnt pathway gene sets.
#' 
#' @title Wnt_purification_system
#' @param file_paths Working directory path.
#' @param expression_accession_vector A vector of all sub-dataset names to be used for Wnt pathway analysis.
#' @param group_HL A file required for each sub-dataset to distinguish between high (H) and low (L) Wnt signaling activity groups. File column 1: Accession, column 2: Group1, column 3: Group2, column 4: Group1_Status, column 5: Group2_Status.
#' @param gene_difference_method Method for identifying differentially expressed genes associated with Wnt pathway activity: "limma", "t_test", "wilcox_test".
#' @param alternative When using t.test/wilcox.test/limma, select "two.sided", "less", or "greater" to assess Wnt activity-associated expression changes.
#' @param p_combine_method Method for integrating p-values across multiple sub-datasets within the same Wnt-related study, options include "fisher", "z.transform", "logit", "ctt".
#' @param using_FC Under "two.sided" condition. TRUE: weight p-values with directionality based on log2FC. FALSE: weight p-values with directionality based on sign(log2FC). Default is "FALSE".
#' @param statistics Method for computing Score2 through specific p-value integration approaches: "arithmetic_mean", "median", "geometric_mean".
#' @param na_ratio Maximum percentage of missing expression values allowed per gene across samples (default 50%). Genes with missing percentage exceeding na_ratio will be removed from Wnt pathway analysis.
#' @param purpose Select "cleaned" for precision refinement of Wnt-related gene sets, or select "validated" for evaluation of external Wnt pathway gene sets.
#' @param quantile_threshold Based on Score2 ranking. Under "two.sided" condition, genes above the (1-quantile_threshold)×100% percentile are retained as Wnt-activated genes, while genes below the (1-quantile_threshold)×100% percentile are retained as Wnt-inhibited genes. Under "greater" and "less" conditions, genes in the upper (1-quantile_threshold)×100% percentile are retained.
#' @param activation_geneset Wnt pathway activation gene sets to be cleaned. (Place in the working directory)
#' @param inhibition_geneset Wnt pathway inhibition gene sets to be cleaned. (Place in the working directory)
#' @param geneSets_gmt Wnt-related gene sets to be validated. (GMT file format, read using the 'read.gmt' function from the 'clusterProfiler' package)
#' @param min.sz Minimum gene set size for fGSEA analysis of Wnt pathway gene sets (default = 1).
#' @param max.sz Maximum gene set size for fGSEA analysis of Wnt pathway gene sets (default = 10000).
#' @author Dingkang Zhao
#' @examples
#' @return
#' @export
Wnt_purification_system=function(file_paths,
                        expression_accession_vector,
                        group_HL,
                        gene_difference_method=c("limma","t_test","wilcox_test"),
                        alternative=c("two.sided","less","greater"),
                        p_combine_method=c("fisher", "z.transform", "logit", "cct", "sumz","geometric_mean"),
                        using_FC=F,
                        na_ratio=0.5,
                        using_KNN=T,
                        statistics=c("arithmetic_mean","median","geometric_mean","sumz"),
                        purpose=c("cleaned","validated"),
                        threshold_type=c("quantile","rank"),
                        rank_threshold=100,
                        quantile_threshold=0.99,
                        activation_geneset=NA,
                        inhibition_geneset=NA,
                        geneSets_gmt=NA,
                        min.sz=1,
                        max.sz=10000,
                        export_file=F
                        
){
  setwd(file_paths)
  library(doBy)
  library(limma)
  library(metap)
  ##########################################Limma one-sided test
  #Artmann S, Jung K, Bleckmann A, Beissbarth T. Detection of simultaneous group effects in microRNA expression and related target gene sets. PLoS One. 2012;7(6):e38365. doi: 10.1371/journal.pone.0038365. Epub 2012 Jun 19. PMID: 22723856; PMCID: PMC3378551.
  limma.one.sided = function(fit, lower = FALSE){
    se.coef=sqrt(fit$s2.post) * fit$stdev.unscaled
    df.total=fit$df.prior + fit$df.residual
    pt(fit$t, df = df.total, lower.tail = lower)
  }
  ##########################################
  library(impute)
  library(survcomp)
  library(clusterProfiler)
  library(fgsea)
  ##########################################CTT: a p-value integration method
  #Liu Y, Xie J. Cauchy combination test: a powerful test with analytic p-value calculation under arbitrary dependency structures. J Am Stat Assoc. 2020;115(529):393-402. doi: 10.1080/01621459.2018.1554485. Epub 2019 Apr 25. PMID: 33012899; PMCID: PMC7531765.
  #Liu Y, Chen S, Li Z, Morrison AC, Boerwinkle E, Lin X. ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant Analysis in Sequencing Studies. Am J Hum Genet. 2019 Mar 7;104(3):410-421. doi: 10.1016/j.ajhg.2019.01.002. PMID: 30849328; PMCID: PMC6407498.
  CCT = function(pvals, weights=NULL){
    if(sum(is.na(pvals)) > 0){
      stop("Cannot have NAs in the p-values!")
    }
    if((sum(pvals<0) + sum(pvals>1)) > 0){
      stop("All p-values must be between 0 and 1!")
    }
    is.zero = (sum(pvals==0)>=1)
    is.one = (sum(pvals==1)>=1)
    if(is.zero && is.one){
      stop("Cannot have both 0 and 1 p-values!")
    }
    if(is.zero){
      return(0)
    }
    if(is.one){
      warning("There are p-values that are exactly 1!")
      return(1)
    }
    if(is.null(weights)){
      weights = rep(1/length(pvals),length(pvals))
    }else if(length(weights)!=length(pvals)){
      stop("The length of weights should be the same as that of the p-values!")
    }else if(sum(weights < 0) > 0){
      stop("All the weights must be positive!")
    }else{
      weights = weights/sum(weights)
    }
    is.small = (pvals < 1e-16)
    if (sum(is.small) == 0){
      cct.stat = sum(weights*tan((0.5-pvals)*pi))
    }else{
      cct.stat = sum((weights[is.small]/pvals[is.small])/pi)
      cct.stat = cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
    }
    if(cct.stat > 1e+15){
      pval = (1/cct.stat)/pi
    }else{
      pval = 1-pcauchy(cct.stat)
    }
    return(pval)
  }
  ##########################################
  vector=expression_accession_vector
  list_data=list()
  list_G_data=list()
  print("Start loading expression data and grouping information.")
  for (q in 1:length(vector)) {
    list_data[[q]]=get(vector[q])
    list_G_data[[q]]=get(paste0(vector[q],"_G"))
    names(list_data)[q]=paste0(vector[q])
    names(list_G_data)[q]=paste0(vector[q],"_G")
  }
  modified_vector=lapply(vector, function(x) {
    ifelse(grepl("_",x),substr(x,1,regexpr("_", x)-1),x)
  })
  modified_vector=unique(as.character(modified_vector))
  large_list=c()
  for (item in 1:length(modified_vector)) {
    assign(paste0("list_", modified_vector[item]),list_data[grep(paste0("^",modified_vector[item]), names(list_data))])
    large_list=c(large_list,paste0("list_", modified_vector[item]))
  }
  print("Start the differential analysis of genes.")
  for (j in 1:length(large_list)) {
    expList=get(large_list[j])
    for (i in 1:length(expList)){
      exp=expList[[i]]
      dimnames=list(rownames(exp),colnames(exp))
      exp=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
      exp=as.data.frame(exp)
      group=list_G_data[[paste0(names(expList)[i],"_G")]]
      group=as.data.frame(group)
      Group_HL=group_HL
      rownames(Group_HL)=Group_HL$Accession
      Group_HL=Group_HL[,-1]
      group_HL_selected=Group_HL[paste0(names(expList)[i]),,drop=F]
      group$group=ifelse(group$group==group_HL_selected$Group1,group_HL_selected$Group1_Status,group_HL_selected$Group2_Status)
      if(gene_difference_method=="t_test"){
        Get_T_test=function(gene=gene,group=group){
          result=NULL
          for (n in 1:nrow(gene)) {
            gene_n=data.frame(t(gene[n,subset(group, group %in% c("H","L"))$Tag]))
            gene_id=names(gene_n)[1]
            names(gene_n)[1]='gene'
            gene_n$Tag=rownames(gene_n)
            gene_n=merge(gene_n, group, by = 'Tag', all.x = TRUE)
            gene_n$group=factor(gene_n$group,levels = c("H","L"))
            P.Value=t.test(gene~group, gene_n,alternative=alternative)$p.value
            if (!is.na(P.Value)) {
              stat=summaryBy(gene~group,gene_n, FUN = c(mean, median))
              result=rbind(result, c(gene_id, as.character(stat[1,1]), stat[1,2], stat[1,3], as.character(stat[2,1]), stat[2,2], stat[2,3], P.Value))
            }
          }
          result=data.frame(result)
          names(result)=c('gene_id', 'group1', 'mean1', 'median1', 'group2', 'mean2', 'median2', 'P.Value')
          result$mean1=as.numeric(result$mean1)
          result$mean2=as.numeric(result$mean2)
          result=result %>% mutate(logFC=mean1-mean2)
          result=as.data.frame(result)
          result$P.Value=as.numeric(result$P.Value)
          result$logFC=as.numeric(result$logFC)
          return(result) 
        } 
        result=Get_T_test(gene=exp,group=group)
      }else if(gene_difference_method=="wilcox_test"){
        Get_wilcox_test=function(gene=gene,group=group){
          result=NULL
          for (n in 1:nrow(gene)) {
            gene_n=data.frame(t(gene[n,subset(group, group %in% c("H","L"))$Tag]))
            gene_id=names(gene_n)[1]
            names(gene_n)[1]='gene'
            gene_n$Tag=rownames(gene_n)
            gene_n=merge(gene_n, group, by = 'Tag', all.x = TRUE)
            gene_n$group=factor(gene_n$group,levels = c("H","L"))
            P.Value=wilcox.test(gene~group, gene_n,alternative=alternative)$p.value
            if (!is.na(P.Value)) {
              stat=summaryBy(gene~group, gene_n, FUN = c(mean, median))
              result=rbind(result, c(gene_id, as.character(stat[1,1]), stat[1,2], stat[1,3], as.character(stat[2,1]), stat[2,2], stat[2,3], P.Value))
            }
          }
          result=data.frame(result)
          names(result)=c('gene_id', 'group1', 'mean1', 'median1', 'group2', 'mean2', 'median2', 'P.Value')
          result$mean1=as.numeric(result$mean1)
          result$mean2=as.numeric(result$mean2)
          result=result %>% mutate(logFC=mean1-mean2)
          result=as.data.frame(result)
          result$P.Value=as.numeric(result$P.Value)
          result$logFC=as.numeric(result$logFC)
          return(result)
        }
        result=Get_wilcox_test(gene=exp,group=group)
      }else if(gene_difference_method=="limma"){
        Get_limma=function(gene=gene,group=group,alternative_limma=alternative){
          design=model.matrix(~0+factor(group$group))
          colnames(design)=levels(factor(group$group))
          rownames(design)=colnames(gene)
          compare=makeContrasts(H-L, levels=design)
          fit=lmFit(gene, design)
          fit=contrasts.fit(fit, compare)
          fit=eBayes(fit)
          if(alternative_limma=="greater"){
            result=limma.one.sided(fit, lower=FALSE)
            result=as.data.frame(result)
            colnames(result)="P.Value"
            result$gene_id=rownames(result)
            topTable=topTable(fit, coef=1, number=Inf)
            topTable$gene_id=rownames(topTable)
            topTable=topTable[,c("gene_id","logFC")]
            result=merge(result,topTable,by="gene_id")
          }else if(alternative_limma=="less"){
            result=limma.one.sided(fit, lower=TRUE)
            result=as.data.frame(result)
            colnames(result)="P.Value"
            result$gene_id=rownames(result)
            topTable=topTable(fit, coef=1, number=Inf)
            topTable$gene_id=rownames(topTable)
            topTable=topTable[,c("gene_id","logFC")]
            result=merge(result,topTable,by="gene_id")
          }else if(alternative_limma=="two.sided"){
            result=topTable(fit, coef=1, number=Inf)
            result=na.omit(result)
            result$gene_id=rownames(result)
          }else{
            print("Please enter a correct method.(greater/less/two.sided)")
          }
          return(result)
        }
        result=Get_limma(gene=exp,group=group,alternative_limma=alternative)
      }else{
        print("Please enter a correct method.(limma/t_test/wilcox_test)")
      }
      print("Start screening genes.")
      result=result[,c("gene_id","logFC","P.Value")]
      assign(paste0(names(expList)[i],"_result"),result)
    }
    if (length(expList)>1){
      for (mp in 1:length(expList)){
        if (mp == 1) { 
          result_single=get(paste0(names(expList)[mp],"_result"))
          result_p=result_single[,c("gene_id","P.Value")]
          colnames(result_p)[2]=c("P.Value1")
          result_FC=result_single[,c("gene_id","logFC")]
          colnames(result_FC)[2]=c("logFC1")
          result_p_all=result_p
          result_FC_all=result_FC
        } else {
          result_single=get(paste0(names(expList)[mp],"_result"))
          result_p=result_single[,c("gene_id","P.Value")]
          colnames(result_p)[2]=paste0("P.Value",mp)
          result_FC=result_single[,c("gene_id","logFC")]
          colnames(result_FC)[2]=paste0("logFC",mp)
          result_p_all=merge(result_p_all,result_p,by="gene_id")
          result_FC_all=merge(result_FC_all,result_FC,by="gene_id")
        }
      }
      rownames(result_p_all)=result_p_all[,1]
      result_p_all=result_p_all[,-1]
      rownames(result_FC_all)=result_FC_all[,1]
      result_FC_all=result_FC_all[,-1]
      if(p_combine_method %in% c("fisher", "z.transform", "logit")){
        result_p_all_combine=apply(result_p_all,1,function(row) combine.test(p=na.omit(row), method=p_combine_method))
      }else if(p_combine_method=="cct"){
        result_p_all_combine=apply(result_p_all,1,function(row) CCT(p=na.omit(row)))
      }else if(p_combine_method=="sumz"){
        result_p_all_combine=apply(result_p_all,1,function(row) as.numeric(sumz(na.omit(row))$p))
      }else if(p_combine_method=="geometric_mean"){
        result_p_all_combine=apply(result_p_all,1,function(row) exp(mean(log(row),na.rm=T)))
      }
      result_p_all_combine=as.data.frame(result_p_all_combine)
      result_p_all_combine$result_p_all_combine[result_p_all_combine$result_p_all_combine == 0]=min(result_p_all_combine$result_p_all_combine[result_p_all_combine$result_p_all_combine!=0])
      result_p_all$p_combine=result_p_all_combine$result_p_all_combine
      result_p_all$rowname=rownames(result_p_all)
      result_p_all=result_p_all[,c("rowname","p_combine")]
      colnames(result_p_all)=c("gene_id","P.Value")
      result_FC_all$FC_mean=rowMeans(result_FC_all)
      result_FC_all$gene_id=rownames(result_FC_all)
      result_FC_all=result_FC_all[,c("gene_id","FC_mean")]
      colnames(result_FC_all)=c("gene_id","logFC")
      result=result_FC_all
      result$P.Value=result_p_all$P.Value
      assign(paste0(gsub("_.*$","",names(expList)[1]),"_result"),result)
    }else{
      print("Genesets are not processed for merging.")
    }
  }
  print("Start building the p-value matrix.")
  for (mn in 1:length(modified_vector)) {
    if (mn == 1){
      pmatrix=get(paste0(modified_vector[mn],"_result"))[,c("gene_id","P.Value")]
      colnames(pmatrix)[mn+1]=paste0(modified_vector[mn])
      FCmatrix=get(paste0(modified_vector[mn],"_result"))[,c("gene_id","logFC")]
      colnames(FCmatrix)[mn+1]=paste0(modified_vector[mn])
    }else{
      pmatrix=merge(pmatrix,get(paste0(modified_vector[mn],"_result"))[,c("gene_id","P.Value")],by="gene_id",all=T)
      colnames(pmatrix)[mn+1]=paste0(modified_vector[mn])
      FCmatrix=merge(FCmatrix,get(paste0(modified_vector[mn],"_result"))[,c("gene_id","logFC")],by="gene_id",all=T)
      colnames(FCmatrix)[mn+1]=paste0(modified_vector[mn])
    }
  }
  rownames(pmatrix)=pmatrix[,1]
  pmatrix=pmatrix[,-1]
  pmatrix$na_count=rowSums(is.na(pmatrix))/ncol(pmatrix)
  pmatrix=pmatrix[pmatrix$na_count <= na_ratio,]
  pmatrix=pmatrix[,!names(pmatrix) %in% "na_count"]
  rownames(FCmatrix)=FCmatrix[,1]
  FCmatrix=FCmatrix[,-1]
  FCmatrix=FCmatrix[rownames(pmatrix),]
  if(using_KNN){
    pmatrix=impute.knn(as.matrix(pmatrix),k=10,rowmax=1,colmax=1)
    pmatrix=pmatrix$data
    pmatrix=as.data.frame(pmatrix)
    FCmatrix=impute.knn(as.matrix(FCmatrix),k=10,rowmax=1,colmax=1)
    FCmatrix=FCmatrix$data
    FCmatrix=as.data.frame(FCmatrix)
  }else{
    print("Completing missing data without using the KNN algorithm.")
  }
  if(statistics=="arithmetic_mean"){
    if(alternative=="two.sided"){
      if(using_FC){
        FCmatrix=FCmatrix
      }else{
        FCmatrix=lapply(FCmatrix, function(x) ifelse(x != 0, sign(x), 0))
      }
    }else{
      FCmatrix[!is.na(FCmatrix)]=1
    }
    pmatrix=(-log10(pmatrix))*FCmatrix
    pmatrix$p_score=apply(pmatrix,1,function(x) mean(x,na.rm=T))
  }else if(statistics=="median"){
    if(alternative=="two.sided"){
      if(using_FC){
        FCmatrix=FCmatrix
      }else{
        FCmatrix=lapply(FCmatrix, function(x) ifelse(x != 0, sign(x), 0))
      }
    }else{
      FCmatrix[!is.na(FCmatrix)]=1
    }
    pmatrix=(-log10(pmatrix))*FCmatrix
    pmatrix$p_score=apply(pmatrix,1,function(x) median(x,na.rm=T))
  }else if(statistics=="geometric_mean"){
    if(alternative=="two.sided"){
      if(using_FC){
        FCmatrix$FC_score=apply(FCmatrix,1,function(x) mean(x,na.rm=T))
      }else{
        FCmatrix$FC_score=apply(FCmatrix,1,function(x) mean(x,na.rm=T))
        FCmatrix$FC_score=sign(FCmatrix$FC_score)
      }
    }else{
      FCmatrix$FC_score=1
    }
    pmatrix$p_score=apply(pmatrix,1,function(x) exp(mean(log(x),na.rm=T)))
    pmatrix$p_score=-log10(pmatrix$p_score)*FCmatrix$FC_score
  }else if(statistics=="sumz"){
    if(alternative=="two.sided"){
      if(using_FC){
        FCmatrix$FC_score=apply(FCmatrix,1,function(x) mean(x,na.rm=T))
      }else{
        FCmatrix$FC_score=apply(FCmatrix,1,function(x) mean(x,na.rm=T))
        FCmatrix$FC_score=sign(FCmatrix$FC_score)
      }
    }else{
      FCmatrix$FC_score=1
    }
    pmatrix$p_score=apply(pmatrix,1,function(x) as.numeric(sumz(na.omit(x))$p))
    pmatrix$p_score=-log10(pmatrix$p_score)*FCmatrix$FC_score
  }else{
    print("Please enter a correct method.(arithmetic_mean/median/geometric_mean/sumz)")
  }
  pmatrix=pmatrix[order(pmatrix$p_score,decreasing = T),]
  if (purpose=="cleaned"){
    if(alternative=="less"){
      if(threshold_type=="quantile"){
        pmatrix_DN=pmatrix[pmatrix$p_score >= quantile(pmatrix$p_score, quantile_threshold),]
      }else if(threshold_type=="rank"){
        pmatrix_DN=pmatrix[1:rank_threshold,]
      }
      Element_DN=as.character(rownames(pmatrix_DN))
      update_inhibition_geneset=list()
      for (mi in 1:length(colnames(inhibition_geneset))){
        coldn=inhibition_geneset[,mi]
        coldn=coldn[coldn != ""]
        col_intersect_dn=intersect(coldn,Element_DN)
        update_inhibition_geneset[[mi]]=col_intersect_dn
        names(update_inhibition_geneset)[mi]=paste0(colnames(inhibition_geneset)[mi],"_PC")
      }
      max_length_DN=max(sapply(update_inhibition_geneset,length))
      update_inhibition_geneset=lapply(update_inhibition_geneset, function(x) {
        if (length(x) < max_length_DN) {
          c(x, rep("", max_length_DN - length(x)))
        } else {
          x
        }
      })
      update_inhibition_geneset=as.data.frame(update_inhibition_geneset)
      update_geneset_I=list(update_inhibition_geneset,pmatrix)
      names(update_geneset_I)=c("update_inhibition_geneset","cleaning_system")
      if(export_file==T){
      write.table(update_inhibition_geneset,file="update_inhibition_geneset_precise.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
      }
      aim=update_geneset_I
    }else if(alternative=="greater"){
      if(threshold_type=="quantile"){
        pmatrix_UP=pmatrix[pmatrix$p_score >= quantile(pmatrix$p_score, quantile_threshold),]
      }else if(threshold_type=="rank"){
        pmatrix_UP=pmatrix[1:rank_threshold,]
      }
      Element_UP=as.character(rownames(pmatrix_UP))
      update_activation_geneset=list()
      for (ma in 1:length(colnames(activation_geneset))){
        colup=activation_geneset[,ma]
        colup=colup[colup != ""]
        col_intersect_up=intersect(colup,Element_UP)
        update_activation_geneset[[ma]]=col_intersect_up
        names(update_activation_geneset)[ma]=paste0(colnames(activation_geneset)[ma],"_PC")
      }
      max_length_UP= max(sapply(update_activation_geneset,length))
      update_activation_geneset=lapply(update_activation_geneset, function(x) {
        if (length(x) < max_length_UP) {
          c(x, rep("", max_length_UP - length(x)))
        } else {
          x
        }
      })
      update_activation_geneset=as.data.frame(update_activation_geneset)
      update_geneset_A=list(update_activation_geneset,pmatrix)
      names(update_geneset_A)=c("update_activation_geneset","cleaning_system")
      if(export_file==T){
      write.table(update_activation_geneset,file="update_activation_geneset_precise.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
      }
      aim=update_geneset_A
    }else if(alternative=="two.sided"){
      if(threshold_type=="quantile"){
        pmatrix_UP=pmatrix[pmatrix$p_score >= quantile(pmatrix$p_score, quantile_threshold),]
        pmatrix_DN=pmatrix[pmatrix$p_score < quantile(pmatrix$p_score, 1-quantile_threshold),]
      }else if(threshold_type=="rank"){
        pmatrix_UP=pmatrix[1:rank_threshold,]
        pmatrix_DN=pmatrix[(length(pmatrix$p_score)-(rank_threshold-1)):length(pmatrix$p_score),]
      }
      Element_UP=as.character(rownames(pmatrix_UP))
      Element_DN=as.character(rownames(pmatrix_DN))
      update_activation_geneset=list()
      update_inhibition_geneset=list()
      for (ma in 1:length(colnames(activation_geneset))){
        colup=activation_geneset[,ma]
        colup=colup[colup != ""]
        col_intersect_up=intersect(colup,Element_UP)
        update_activation_geneset[[ma]]=col_intersect_up
        names(update_activation_geneset)[ma]=paste0(colnames(activation_geneset)[ma],"_PC")
      }
      for (mi in 1:length(colnames(inhibition_geneset))){
        coldn=inhibition_geneset[,mi]
        coldn=coldn[coldn != ""]
        col_intersect_dn=intersect(coldn,Element_DN)
        update_inhibition_geneset[[mi]]=col_intersect_dn
        names(update_inhibition_geneset)[mi]=paste0(colnames(inhibition_geneset)[mi],"_PC")
      }
      max_length_UP= max(sapply(update_activation_geneset,length))
      update_activation_geneset=lapply(update_activation_geneset, function(x) {
        if (length(x) < max_length_UP) {
          c(x, rep("", max_length_UP - length(x)))
        } else {
          x
        }
      })
      max_length_DN=max(sapply(update_inhibition_geneset,length))
      update_inhibition_geneset=lapply(update_inhibition_geneset, function(x) {
        if (length(x) < max_length_DN) {
          c(x, rep("", max_length_DN - length(x)))
        } else {
          x
        }
      })
      update_activation_geneset=as.data.frame(update_activation_geneset)
      update_inhibition_geneset=as.data.frame(update_inhibition_geneset)
      update_geneset=list(update_activation_geneset,update_inhibition_geneset,pmatrix)
      names(update_geneset)=c("update_activation_geneset","update_inhibition_geneset","cleaning_system")
      if(export_file==T){
      write.table(update_activation_geneset,file="update_activation_geneset_precise.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
      write.table(update_inhibition_geneset,file="update_inhibition_geneset_precise.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
      }
      aim=update_geneset
    }else{
      print("Please enter a correct method.(two.sided/greater/less)")
    }
  }else if(purpose=="validated"){
    FCmatrix$FC_score=rowMeans(FCmatrix[,1:(length(colnames(FCmatrix))-1)],na.rm=T)
    FCmatrix=FCmatrix[order(FCmatrix$FC_score,decreasing = T),]
    canonical_pathways=geneSets_gmt
    canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
    canonical_pathways$term=as.character(canonical_pathways$term)
    canonical_pathways=canonical_pathways %>% split(x = .$gene, f = .$term)
    id=FCmatrix$FC_score
    names(id)=rownames(FCmatrix)
    fgseaRes=fgsea(pathways=canonical_pathways,
                   stats=id,
                   minSize=min.sz,
                   maxSize=max.sz)
    fgseaRes=data.frame(fgseaRes)
    aim=list(fgseaRes,FCmatrix)
    names(aim)=c("fGSEA_result","validated_system")
  }else{
    print("Please enter a correct method.(cleaned/validated)")
  }
  print("This calculation is over.")
  return(aim)
}

 

