# 固定PAM评分及五个cfRNA与骨参数的关联分析。
options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(DESeq2)
  library(readxl)
  library(ggplot2)
  library(ragg)
  library(jsonlite)
})

project_dir <- 'F:/DE-analysis/final/cfRNA_project'
supp_dir <- 'F:/DE-analysis/biomedicines/Supplementary Data and Fig 2'
output_dir <- file.path('F:/DE-analysis/biomedicines', '\u8fd4\u4fee', 'Reviewer2_Comment1')
seed <- 20260717L
n_permutations <- 99999L
n_bootstrap <- 10000L
confidence_level <- 0.95
permutation_tolerance <- 1e-12
six_classes <- c('mRNA','lncRNA','snoRNA','snRNA','miRNA','piRNA')
sample_ids <- c(paste0('N',1:10),paste0('OP',1:10))
expected_features <- c('ENSRNOG00000000047','ENSRNOG00000036688',
                       'piR-rno-2281228','ENSRNOG00000054973','piR-rno-2096851')
feature_rds <- file.path(project_dir,'03_results/03_ML/02_strict_model/01_feature_selection/rds/strict_feature_selection.rds')
de_rds <- file.path(project_dir,'03_results/02_DE/rds/DE_objects.rds')
prediction_file <- file.path(project_dir,'03_results/03_ML/02_strict_model/03_model_evaluation/tables/Model_predictions.csv')
feature_file <- file.path(project_dir,'03_results/03_ML/02_strict_model/01_feature_selection/tables/final_features.csv')
split_file <- file.path(project_dir,'03_results/03_ML/02_strict_model/01_feature_selection/tables/sample_split.csv')
ct_file <- file.path(supp_dir,'Supplementary Data 1.xlsx')
score_file <- file.path(supp_dir,'Supplementary Data 6.xlsx')
stopifnot(all(file.exists(c(feature_rds,de_rds,prediction_file,feature_file,split_file,ct_file,score_file))))

main <- function() {
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  original_workdir <- getwd()
  setwd(output_dir)
  on.exit(setwd(original_workdir), add=TRUE)
  output_dir <- '.'
  if (file.exists(file.path(output_dir,'correlation_results.csv'))) stop('Analysis results already exist; do not overwrite.')
  sink(file.path(output_dir,'analysis_log.txt'), split=TRUE)
  on.exit(sink(), add=TRUE)
  cat('R:',R.version.string,'\n')
  cat('Seed:',seed,' Permutations:',n_permutations,' Bootstrap:',n_bootstrap,'\n')
  feature_object <- readRDS(feature_rds)
  cat('\nFeature object structure:\n'); print(names(feature_object)); str(feature_object,max.level=2)
  de_object <- readRDS(de_rds)
  cat('\nDE object structure:\n'); print(names(de_object)); str(de_object,max.level=2)
  stopifnot(all(c('final_features','train_ids','test_ids') %in% names(feature_object)),
            identical(names(de_object),c('dds','result','metadata')))

  features <- read.csv(feature_file,check.names=FALSE)
  split <- read.csv(split_file,check.names=FALSE)
  stopifnot(identical(features$feature_id,expected_features),
            identical(feature_object$final_features$feature_id,expected_features),
            !anyDuplicated(split$sample_id),setequal(split$sample_id,sample_ids))
  split <- split[match(sample_ids,split$sample_id),]
  stopifnot(identical(split$group,c(rep('N',10),rep('OP',10))),
            sum(split$dataset=='Training')==12L,sum(split$dataset=='Validation')==8L,
            setequal(split$sample_id[split$dataset=='Training'],feature_object$train_ids),
            setequal(split$sample_id[split$dataset=='Validation'],feature_object$test_ids))

  predictions <- read.csv(prediction_file,check.names=FALSE)
  predictions <- predictions[predictions$Model=='PAM',]
  stopifnot(nrow(predictions)==20L,!anyDuplicated(predictions$Sample_ID),
            setequal(predictions$Sample_ID,sample_ids))
  predictions <- predictions[match(sample_ids,predictions$Sample_ID),]
  stopifnot(all(is.finite(predictions$Probability)),all(predictions$Probability>=0 & predictions$Probability<=1),
            identical(predictions$Observed,ifelse(split$group=='OP','OP_like','N_like')),
            identical(predictions$Set,paste(split$dataset,'set')))
  published_scores <- as.data.frame(read_excel(score_file,sheet='sample_predictions',skip=1,.name_repair='minimal'))
  published_scores <- published_scores[published_scores$Model=='PAM',]
  stopifnot(nrow(published_scores)==20L,!anyDuplicated(published_scores$Sample_ID))
  published_scores <- published_scores[match(sample_ids,published_scores$Sample_ID),]
  stopifnot(max(abs(as.numeric(published_scores$Probability)-predictions$Probability))<1e-14)

  ct <- as.data.frame(read_excel(ct_file,sheet='micro_CT',range='A2:H22',
                               col_types=c('text','text',rep('numeric',6)),.name_repair='minimal'))
  stopifnot(nrow(ct)==20L,!anyDuplicated(ct$sample_id),setequal(ct$sample_id,sample_ids),
            all(c('sample_id','group','BV/TV (%)','BMD (mg HA/cm\u00b3)') %in% names(ct)))
  ct <- ct[match(sample_ids,ct$sample_id),]
  stopifnot(identical(ct$group,split$group))

  dds <- de_object$dds
  stopifnot(setequal(colnames(dds),sample_ids),all(expected_features %in% rownames(dds)))
  sf <- sizeFactors(dds)
  stopifnot(length(sf)==20L,all(is.finite(sf)),all(sf>0),is.null(normalizationFactors(dds)))
  expression <- counts(dds,normalized=TRUE)[expected_features,sample_ids,drop=FALSE]
  raw_expression <- counts(dds,normalized=FALSE)[expected_features,sample_ids,drop=FALSE]
  stopifnot(max(abs(expression-sweep(raw_expression,2L,sf[sample_ids],'/')))<1e-10)
  expression <- as.data.frame(t(expression),check.names=FALSE)
  master <- data.frame(sample_id=sample_ids,group=split$group,dataset=split$dataset,
                       PAM_probability=predictions$Probability,
                       BMD_mgHA_cm3=ct[['BMD (mg HA/cm\u00b3)']],BVTV_percent=ct[['BV/TV (%)']],
                       expression,check.names=FALSE)
  stopifnot(identical(master$sample_id,sample_ids),all(is.finite(as.matrix(master[,4:11]))))
  variables <- c('PAM_probability',expected_features)
  outcomes <- c('BMD_mgHA_cm3','BVTV_percent')
  scopes <- list(All=seq_len(20L),Sham=which(master$group=='N'),OVX=which(master$group=='OP'))
  for (scope in names(scopes)) {
    for (variable in c(variables,outcomes)) {
      if (length(unique(master[scopes[[scope]],variable]))<2L) stop('Constant variable: ',scope,'/',variable)
    }
  }
  write.csv(master,file.path(output_dir,'sample_analysis_table.csv'),row.names=FALSE)
  write.csv(data.frame(sample_id=sample_ids,DESeq2_size_factor=sf[sample_ids]),
            file.path(output_dir,'expression_size_factors.csv'),row.names=FALSE)
  write.csv(features[,c('feature_id','RNA_type')],file.path(output_dir,'feature_identity.csv'),row.names=FALSE)
  cat('\nMatched 20 animals, 5 cfRNAs, 6 score/expression variables, and 2 CT outcomes.\n')
  print(table(master$group,master$dataset))

  set.seed(seed)
  results <- list()
  bootstrap_distributions <- list()
  result_index <- 0L
  for (scope in names(scopes)) {
    dat <- master[scopes[[scope]],,drop=FALSE]
    n <- nrow(dat)
    x <- as.matrix(dat[,variables,drop=FALSE])
    y <- as.matrix(dat[,outcomes,drop=FALSE])
    xr <- apply(x,2L,rank,ties.method='average')
    yr <- apply(y,2L,rank,ties.method='average')
    xc <- sweep(xr,2L,colMeans(xr),'-')
    yc <- sweep(yr,2L,colMeans(yr),'-')
    denominator <- outer(sqrt(colSums(xc^2)),sqrt(colSums(yc^2)))
    observed <- crossprod(xc,yc)/denominator
    stopifnot(max(abs(observed-cor(x,y,method='spearman')))<1e-12)
    permutation_index <- replicate(n_permutations,sample.int(n))
    exceedances <- matrix(0L,nrow=length(variables),ncol=length(outcomes))
    for (j in seq_along(outcomes)) {
      yp <- matrix(yc[permutation_index,j],nrow=n,ncol=n_permutations)
      rho_perm <- sweep(crossprod(xc,yp),1L,denominator[,j],'/')
      exceedances[,j] <- rowSums(abs(rho_perm)>=abs(observed[,j])-permutation_tolerance)
    }
    p_values <- (exceedances+1)/(n_permutations+1)

    strata <- split(seq_len(n),dat$group)
    bootstrap_index <- replicate(n_bootstrap,unlist(lapply(strata,function(indices) {
      indices[sample.int(length(indices),length(indices),replace=TRUE)]
    }),use.names=FALSE))
    bootstrap_rho <- vapply(seq_len(n_bootstrap),function(b) {
      ii <- bootstrap_index[,b]
      as.vector(cor(x[ii,,drop=FALSE],y[ii,,drop=FALSE],method='spearman'))
    },numeric(length(variables)*length(outcomes)))
    if (any(!is.finite(bootstrap_rho))) stop('Undefined bootstrap correlation in ',scope,'; no replicates were discarded.')
    alpha <- (1-confidence_level)/2
    lower <- apply(bootstrap_rho,1L,quantile,probs=alpha,names=FALSE,type=7)
    upper <- apply(bootstrap_rho,1L,quantile,probs=1-alpha,names=FALSE,type=7)
    bootstrap_distributions[[scope]] <- bootstrap_rho
    grid <- expand.grid(feature=variables,bone_parameter=outcomes,stringsAsFactors=FALSE)
    result_index <- result_index+1L
    results[[result_index]] <- data.frame(analysis_group=scope,grid,n=n,
      Spearman_rho=as.vector(observed),CI95_low=lower,CI95_high=upper,
      P_permutation=as.vector(p_values),permutation_exceedances=as.vector(exceedances),
      bootstrap_valid=n_bootstrap,check.names=FALSE)
    cat('Completed:',scope,'n=',n,'12 associations; all bootstrap replicates finite.\n')
  }
  result <- do.call(rbind,results)
  result$BH_FDR <- p.adjust(result$P_permutation,method='BH')
  stopifnot(nrow(result)==36L,all(is.finite(as.matrix(result[,c('Spearman_rho','CI95_low','CI95_high','P_permutation','BH_FDR')]))))
  result <- result[,c('analysis_group','feature','bone_parameter','n','Spearman_rho',
                      'CI95_low','CI95_high','P_permutation','BH_FDR','permutation_exceedances','bootstrap_valid')]
  write.csv(result,file.path(output_dir,'correlation_results.csv'),row.names=FALSE)
  write.csv(result[result$feature=='PAM_probability',],file.path(output_dir,'PAM_score_correlations.csv'),row.names=FALSE)
  write.csv(result[result$feature!='PAM_probability',],file.path(output_dir,'cfRNA_correlations.csv'),row.names=FALSE)

  parameters <- list(seed=seed,n_permutations=n_permutations,n_bootstrap=n_bootstrap,
    correlation='Spearman; average ranks; no probability rounding',
    p_value='(1 + number of permutations with absolute rho >= observed absolute rho - 1e-12)/(99999 + 1)',
    CI='95% paired percentile bootstrap; type-7 quantiles; pooled bootstrap stratified by N/OP',
    multiplicity='One BH family comprising all 36 pre-specified tests',
    expression='Existing DESeq2 normalized counts from the full-cohort DE_objects.rds; no parameter refitting',
    score='Previously saved PAM OP_like probability; no model retraining or feature selection',
    inference='Exploratory associations conditional on the existing selected features/model; pooled data include training animals',
    R_version=R.version.string,packages=lapply(c('DESeq2','ggplot2','readxl','ragg'),function(p)list(package=p,version=as.character(packageVersion(p)))))
  write_json(parameters,file.path(output_dir,'analysis_parameters.json'),pretty=TRUE,auto_unbox=TRUE)
  saveRDS(list(parameters=parameters,master=master,features=features[,c('feature_id','RNA_type')],
               results=result,bootstrap_distributions=bootstrap_distributions),file.path(output_dir,'association_objects.rds'))
  writeLines(c(
    'Methods',
    'The locked PAM probabilities and five pre-selected cfRNAs were matched to BMD and BV/TV by confirmed animal IDs.',
    'Individual cfRNA expression used normalized counts from the existing full-cohort DESeq2 object and its saved size factors. No normalization parameters, feature selection, or predictive models were re-fitted.',
    'Spearman correlations were evaluated in all 20 rats, 10 Sham rats, and 10 OVX rats.',
    'Two-sided Monte Carlo pairing-permutation P values used 99,999 permutations and absolute rho, with a plus-one correction and 1e-12 comparison tolerance.',
    '95% paired percentile bootstrap confidence intervals used 10,000 resamples; the pooled analysis retained 10 animals per group in each resample.',
    'BH correction was applied once across all 36 comparisons. The random seed was 20260717.',
    'Scatterplots show the saved PAM probability or log2(normalized cfRNA count + 1); correlations use the original continuous values.',
    'Intervals condition on the already selected features and fitted model and do not account for prior model-selection uncertainty. All-sample associations include training observations.',
    'Training and validation VST values were not concatenated because their original transformation procedures differ.',
    'References: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html',
    'https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html',
    'https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html',
    'https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html'
  ),file.path(output_dir,'analysis_methods.txt'))

  figure_dir <- file.path(output_dir,'figures')
  dir.create(figure_dir,showWarnings=FALSE)
  theme_set(theme_bw(base_size=10,base_family='Arial')+
      theme(panel.grid.minor=element_blank(),strip.background=element_rect(fill='grey95'),
            legend.position='bottom',plot.title=element_text(face='bold')))
  save_plot <- function(plot,name,width,height) {
    ggsave(file.path(figure_dir,paste0(name,'.png')),plot=plot,device=ragg::agg_png,width=width,height=height,dpi=300,bg='white')
    ggsave(file.path(figure_dir,paste0(name,'.pdf')),plot=plot,device=grDevices::cairo_pdf,width=width,height=height,bg='white')
  }
  outcome_labels <- c(BMD_mgHA_cm3='BMD (mg HA/cm3)',BVTV_percent='BV/TV (%)')
  scope_labels <- c(All='All rats (n=20)',Sham='Sham (n=10)',OVX='OVX (n=10)')
  plot_result <- result
  plot_result$analysis_group <- factor(plot_result$analysis_group,levels=names(scopes),labels=scope_labels)
  plot_result$bone_parameter <- factor(plot_result$bone_parameter,levels=outcomes,labels=outcome_labels)
  plot_result$feature <- factor(plot_result$feature,levels=rev(variables),labels=rev(c('PAM score',expected_features)))
  plot_result$FDR <- factor(ifelse(plot_result$BH_FDR<0.05,'< 0.05','>= 0.05'),levels=c('< 0.05','>= 0.05'))
  forest <- ggplot(plot_result,aes(x=Spearman_rho,y=feature,color=FDR))+
    geom_vline(xintercept=0,color='grey65',linetype='dashed')+
    geom_segment(aes(x=CI95_low,xend=CI95_high,yend=feature),linewidth=0.55)+
    geom_point(size=2)+facet_grid(bone_parameter~analysis_group)+
    scale_color_manual(values=c('< 0.05'='#1F4E79','>= 0.05'='#777777'),name='BH FDR')+
    scale_x_continuous(limits=c(-1.02,1.02),breaks=c(-1,-0.5,0,0.5,1))+
    labs(x='Spearman rho (95% bootstrap CI)',y=NULL,title='Associations with bone measurements')
  save_plot(forest,'association_summary',13,7.5)

  format_p <- function(p) formatC(p,format='g',digits=3)
  annotation <- function(row) sprintf('rho = %.2f; P = %s\nBH FDR = %s',row$Spearman_rho,format_p(row$P_permutation),format_p(row$BH_FDR))
  score_points <- list(); score_annotations <- list(); labels <- character(); k <- 0L
  for (outcome in outcomes) for (scope in names(scopes)) {
    k <- k+1L; dat <- master[scopes[[scope]],]
    label <- paste(outcome_labels[[outcome]],scope_labels[[scope]],sep=' | ')
    labels <- c(labels,label)
    score_points[[k]] <- data.frame(panel=label,x=dat$PAM_probability,y=dat[[outcome]],group=dat$group,dataset=dat$dataset)
    r <- result[result$feature=='PAM_probability' & result$bone_parameter==outcome & result$analysis_group==scope,]
    score_annotations[[k]] <- data.frame(panel=label,label=annotation(r))
  }
  points <- do.call(rbind,score_points); notes <- do.call(rbind,score_annotations)
  points$panel <- factor(points$panel,levels=labels); notes$panel <- factor(notes$panel,levels=labels)
  score_plot <- ggplot(points,aes(x=x,y=y,color=group,shape=dataset))+geom_point(size=2.4)+
    geom_text(data=notes,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-0.03,vjust=1.15,size=3)+
    facet_wrap(~panel,ncol=3,scales='free')+
    scale_color_manual(values=c(N='#D73027',OP='#4575B4'),labels=c(N='Sham',OP='OVX'),name='Group')+
    scale_shape_manual(values=c(Training=16,Validation=17),name='Original split')+
    scale_y_continuous(expand=expansion(mult=c(0.06,0.3)))+
    labs(x='PAM probability of OVX',y=NULL,title='Locked PAM score and bone measurements')
  save_plot(score_plot,'PAM_score_scatter',12,7)

  for (scope in names(scopes)) {
    dat <- master[scopes[[scope]],]
    pieces <- list(); anns <- list(); labels <- character(); k <- 0L
    for (feature in expected_features) for (outcome in outcomes) {
      k <- k+1L
      label <- paste(feature,outcome_labels[[outcome]],sep=' | ')
      labels <- c(labels,label)
      pieces[[k]] <- data.frame(panel=label,x=log2(dat[[feature]]+1),y=dat[[outcome]],group=dat$group,dataset=dat$dataset)
      r <- result[result$feature==feature & result$bone_parameter==outcome & result$analysis_group==scope,]
      anns[[k]] <- data.frame(panel=label,label=annotation(r))
    }
    points <- do.call(rbind,pieces); notes <- do.call(rbind,anns)
    points$panel <- factor(points$panel,levels=labels); notes$panel <- factor(notes$panel,levels=labels)
    plot <- ggplot(points,aes(x=x,y=y,color=group,shape=dataset))+geom_point(size=2)+
      geom_text(data=notes,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-0.03,vjust=1.15,size=2.7)+
      facet_wrap(~panel,ncol=2,scales='free')+
      scale_color_manual(values=c(N='#D73027',OP='#4575B4'),labels=c(N='Sham',OP='OVX'),name='Group')+
      scale_shape_manual(values=c(Training=16,Validation=17),name='Original split')+
      scale_y_continuous(expand=expansion(mult=c(0.06,0.35)))+
      labs(x='Log2(normalized count + 1)',y=NULL,title=paste('Five cfRNAs and bone measurements:',scope_labels[[scope]]))
    save_plot(plot,paste0('cfRNA_scatter_',scope),10,13)
  }
  cat('\nPAM score results:\n')
  print(result[result$feature=='PAM_probability',c('analysis_group','bone_parameter','n','Spearman_rho','CI95_low','CI95_high','P_permutation','BH_FDR')],row.names=FALSE)
  cat('\nWithin-group results with BH FDR < 0.05:\n')
  print(result[result$analysis_group!='All' & result$BH_FDR<0.05,],row.names=FALSE)
  cat('\nFinished. Output:',output_dir,'\n')
  writeLines(capture.output(sessionInfo()),file.path(output_dir,'sessionInfo.txt'))
}
main()
