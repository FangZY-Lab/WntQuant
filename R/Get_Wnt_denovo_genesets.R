#' To obtain de novo Wnt signaling pathway-related gene sets from multiple bulk RNA-Seq datasets,
#' differentially expressed genes associated with Wnt pathway activity are first identified via limma/t.test/wilcox.test,
#' then up-regulated and down-regulated differentially expressed genes are designated as two de novo gene sets representing Wnt pathway activation and inhibition, respectively.
#' 
#' @title Get_Wnt_denovo_genesets
#' @param file_paths Working directory path.
#' @param expression_accession_vector A vector of all sub-dataset names to be used for de novo Wnt pathway-related gene set identification.
#' @param group_HL A file required for each sub-dataset to distinguish between high (H) and low (L) Wnt signaling activity groups. File column 1: Accession, column 2: Group1, column 3: Group2, column 4: Group1_Status, column 5: Group2_Status.
#' @param gene_difference_method Method for identifying differentially expressed genes associated with Wnt pathway activity: "limma", "t_test", "wilcox_test".
#' @param alternative When using t.test/wilcox.test, select "two.sided", "less", or "greater" to assess Wnt activity-associated expression changes.
#' @param p_combine_method Method for integrating p-values across multiple sub-datasets within the same Wnt-related study, options include "fisher", "z.transform", "logit", "ctt".
#' @param threshold_or_rank Strategy for selecting Wnt pathway-associated differentially expressed genes: "rank" for ranking by p-value, or "threshold" for applying significance cutoffs.
#' @param top_genes If "rank" is selected for threshold_or_rank, the top-ranked genes associated with differential Wnt pathway activity will be retained (default = 500).
#' @param gene_Pfilter If "threshold" is selected for threshold_or_rank, this defines the p-value cutoff for Wnt-related differential expression (default = 0.05).
#' @param gene_FCfilter If "threshold" is selected for threshold_or_rank, this defines the log2 fold change threshold for Wnt pathway activity-associated differential expression (default = 1).
#' @author Dingkang Zhao
#' @examples
#' @return A list of de novo gene sets representing Wnt pathway activation (up-regulated genes) and inhibition (down-regulated genes).
#' @export
Get_Wnt_denovo_genesets=function(file_paths,
                            expression_accession_vector,
                            group_HL,
                            gene_difference_method=c("limma","t_test","wilcox_test"),
                            alternative=c("two.sided","less","greater"),
                            p_combine_method=c("fisher", "z.transform", "logit", "cct","sumz","geometric_mean"),
                            threshold_or_rank=c("rank","threshold"),
                            top_genes=500,
                            gene_Pfilter=0.05,
                            gene_FCfilter=1,
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
  library(pheatmap)
  library(ggplot2)
  library(svglite)
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
        result_p_all_combine=apply(result_p_all,1,function(row)  exp(mean(log(row),na.rm=T)))
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
  print("Start screening genes.")
  for (mq in 1:length(modified_vector)) {
    result=get(paste0(modified_vector[mq],"_result"))
    if (threshold_or_rank=="rank"){
      if(alternative=="two.sided"){
        upregulated_gene=subset(result,logFC>0)
        downregulated_gene=subset(result,logFC<0)
        upregulated_gene=upregulated_gene[order(upregulated_gene$P.Value), ]
        upregulated_gene=upregulated_gene[1:top_genes, ]
        upregulated_gene=upregulated_gene$gene_id
        downregulated_gene=downregulated_gene[order(downregulated_gene$P.Value), ]
        downregulated_gene=downregulated_gene[1:top_genes, ]
        downregulated_gene=downregulated_gene$gene_id
        assign(paste0(modified_vector[mq],"_UP"),upregulated_gene)
        assign(paste0(modified_vector[mq],"_DN"),downregulated_gene)
      }else if(alternative=="less"){
        downregulated_gene=subset(result,logFC<0)
        downregulated_gene=downregulated_gene[order(downregulated_gene$P.Value), ]
        downregulated_gene=downregulated_gene[1:top_genes, ]
        downregulated_gene=downregulated_gene$gene_id
        assign(paste0(modified_vector[mq],"_DN"),downregulated_gene)
      }else if(alternative=="greater"){
        upregulated_gene=subset(result,logFC>0)
        upregulated_gene=upregulated_gene[order(upregulated_gene$P.Value), ]
        upregulated_gene=upregulated_gene[1:top_genes, ]
        upregulated_gene=upregulated_gene$gene_id
        assign(paste0(modified_vector[mq],"_UP"),upregulated_gene)
      }else{
        print("Please enter a correct method.(greater/less/two.sided)")
      }
    }else if(threshold_or_rank=="threshold"){
      if(alternative=="two.sided"){
        Pfilter=gene_Pfilter
        FCfilter=log2(gene_FCfilter)
        result$Significant=ifelse((result$P.Value<Pfilter & abs(result$logFC)>FCfilter),ifelse(result$logFC>FCfilter,"Up","Down"),"No")
        upregulated_gene=subset(result, P.Value < Pfilter & logFC > FCfilter)
        upregulated_gene=upregulated_gene$gene_id
        downregulated_gene=subset(result, P.Value < Pfilter & logFC < -FCfilter)
        downregulated_gene=downregulated_gene$gene_id
        assign(paste0(modified_vector[mq],"_UP"),upregulated_gene)
        assign(paste0(modified_vector[mq],"_DN"),downregulated_gene)
      }else if(alternative=="less"){
        Pfilter=gene_Pfilter
        FCfilter=log2(gene_FCfilter)
        result$Significant=ifelse((result$P.Value<Pfilter & abs(result$logFC)>FCfilter),ifelse(result$logFC>FCfilter,"Up","Down"),"No")
        downregulated_gene=subset(result, P.Value < Pfilter & logFC < -FCfilter)
        downregulated_gene=downregulated_gene$gene_id
        assign(paste0(modified_vector[mq],"_DN"),downregulated_gene)
      }else if(alternative=="greater"){
        Pfilter=gene_Pfilter
        FCfilter=log2(gene_FCfilter)
        result$Significant=ifelse((result$P.Value<Pfilter & abs(result$logFC)>FCfilter),ifelse(result$logFC>FCfilter,"Up","Down"),"No")
        upregulated_gene=subset(result, P.Value < Pfilter & logFC > FCfilter)
        upregulated_gene=upregulated_gene$gene_id
        assign(paste0(modified_vector[mq],"_UP"),upregulated_gene)
      }else{
        print("Please enter a correct method.(greater/less/two.sided)")
      }
    }else{
      print("Please enter a correct method.(rank/threshold)")
    }
  }
  if(alternative=="two.sided"){
    combined_UP=list()
    combined_DN=list()
    for (mn in 1:length(modified_vector)) {
      combined_UP[[mn]]=get(paste0(modified_vector[mn],"_UP"))
      names(combined_UP)[mn]=paste0(modified_vector[mn], "_UP")
      combined_DN[[mn]]=get(paste0(modified_vector[mn],"_DN"))
      names(combined_DN)[mn]=paste0(modified_vector[mn], "_DN")
    }
    max_length_UP=max(sapply(combined_UP, length))
    denovo_genesets_UP=lapply(combined_UP, function(x) {
      if (length(x) < max_length_UP) {
        c(x, rep("", max_length_UP - length(x)))
      } else {
        x
      }
    })
    max_length_DN=max(sapply(combined_DN, length))
    denovo_genesets_DN=lapply(combined_DN, function(x) {
      if (length(x) < max_length_DN) {
        c(x, rep("", max_length_DN - length(x)))
      } else {
        x
      }
    })
    denovo_genesets_UP=as.data.frame(denovo_genesets_UP)
    denovo_genesets_DN=as.data.frame(denovo_genesets_DN)
    if(export_file==T){
    write.table(denovo_genesets_UP, file = "denovo_genesets_UP.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(denovo_genesets_DN, file = "denovo_genesets_DN.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    }
    padded_list=list(denovo_genesets_UP,denovo_genesets_DN)
    names(padded_list)=c("activation","inhibition")
  }else if(alternative=="less"){
    combined_DN=list()
    for (mn in 1:length(modified_vector)) {
      combined_DN[[mn]]=get(paste0(modified_vector[mn],"_DN"))
      names(combined_DN)[mn]=paste0(modified_vector[mn], "_DN")
    }
    max_length_DN=max(sapply(combined_DN, length))
    denovo_genesets_DN=lapply(combined_DN, function(x) {
      if (length(x) < max_length_DN) {
        c(x, rep("", max_length_DN - length(x)))
      } else {
        x
      }
    })
    denovo_genesets_DN=as.data.frame(denovo_genesets_DN)
    if(export_file==T){
    write.table(denovo_genesets_DN, file = "denovo_genesets_DN.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    }
    padded_list=denovo_genesets_DN
  }else if(alternative=="greater"){
    combined_UP=list()
    for (mn in 1:length(modified_vector)) {
      combined_UP[[mn]]=get(paste0(modified_vector[mn],"_UP"))
      names(combined_UP)[mn]=paste0(modified_vector[mn], "_UP")
    }
    max_length_UP=max(sapply(combined_UP, length))
    denovo_genesets_UP=lapply(combined_UP, function(x) {
      if (length(x) < max_length_UP) {
        c(x, rep("", max_length_UP - length(x)))
      } else {
        x
      }
    })
    denovo_genesets_UP=as.data.frame(denovo_genesets_UP)
    if(export_file==T){
    write.table(denovo_genesets_UP, file = "denovo_genesets_UP.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    }
    padded_list=denovo_genesets_UP
  }else{
    print("Please enter a correct method.(greater/less/two.sided)")
  }
  print("This calculation is over.")
  return(padded_list)
}
