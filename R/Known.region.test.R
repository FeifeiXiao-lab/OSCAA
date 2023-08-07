# library("emdbook") #dmvnorm function
# library("cluster")
# library("factoextra")
# library("Rcpp")
# library("ClusterR")
# library("umap")

# create dummy variables for copy number states( apply Two-Dimensional GMM with first two PCs)
dummy.var<-function(s.vec,nCNV){
  # intercept for cnv1
  intercept<-rep(1,length(s.vec))
  # comparsion of other cnvs with cnv1
  dumm.var<-matrix(rep(0,(nCNV-1)*length(s.vec)),ncol =nCNV-1,byrow = T )
  
  colnames(dumm.var)<-paste0("CN",sort(unique(s.vec))[-1],"-",min(unique(s.vec)))
  
  for (i in 1:length(sort(unique(s.vec))[-1]) ) {
    dumm.var[,i]<-ifelse(s.vec==sort(unique(s.vec))[i+1],1,0)
  }
  
  dummy.vars<-cbind(intercept,dumm.var)
  
  return(dummy.vars=dummy.vars)
}

# ----------------------------------------------------------- #
# data augmentation - prepare the input for gaussian mixture model

ExpandData.2D<-function(trait, # Each element of the list can be either a vector of quantitative traits or a vector of 0 and 1 in a case/control framework
                        subID=NULL, # individual IDs
                        signal, #vector or signal matrix
                        covariates=NULL, # covariates to be adjusted
                        nCNV,  # number of CNV states
                        CNVs  # unique CN vector
){
  # signal is a vector
  if ( !is.matrix(signal) ){
    signal<-as.matrix(signal)
  }
  
  if (is.null(covariates)){
    # signal is a matrix
    if (is.null(subID)){
      subID<-seq(1:nrow(signal))
    }
    sortID<-seq(1:nrow(signal))
    data.for.GM<-c()
    for (i in 1:nCNV){
      tmp.data<-data.frame(subID=subID,
                           sortID=sortID,
                           trait=trait,
                           signal=signal,
                           CNV.state=CNVs[i])
      data.for.GM<-rbind(data.for.GM,tmp.data)
    }
    
  }
  
  if (!is.null(covariates)){
    if ( !is.matrix(covariates) ){
      covariates<-as.matrix(covariates)
    }
    # signal is a matrix
    if (is.null(subID)){
      subID<-seq(1:nrow(signal))
    }
    sortID<-seq(1:nrow(signal))
    data.for.GM<-c()
    for (i in 1:nCNV){
      tmp.data<-data.frame(subID=subID,
                           sortID=sortID,
                           trait=trait,
                           signal=signal,
                           CNV.state=CNVs[i],
                           covariates)
      data.for.GM<-rbind(data.for.GM,tmp.data)
    }
    
  }
  
  data.for.GM<-data.for.GM[order(data.for.GM$sortID,data.for.GM$CNV.state),]
  rownames(data.for.GM)<-seq(1,nrow(data.for.GM))
  return(data.for.GM=data.for.GM)
}
# data<-ExpandData(batch=batch,
#                  trait=trait, 
#                  subID=sample,
#                  signal=pca.signal,
#                  nCNV=3,
#                  association.strata = NULL)



# ---------------------------------------------------------------------- #
# calculate covariance for each covariates level
cov.cal<-function(design.mat,
                  unique.design.mat,
                  signal1,
                  mean1,
                  signal2,
                  mean2,
                  weights){
  cov.prod<-weights*(signal1-mean1)*(signal2-mean2)
  cov.list<-rep(1,length(cov.prod))
  for (i in 1:nrow(unique.design.mat)){
    tmp.seq<-which(apply(design.mat, 1, function(x) return(all(x == unique.design.mat[i,]))))
    cov.list[tmp.seq]<-sum(cov.prod[tmp.seq])/sum(weights[tmp.seq])
  }
  return(cov.list)
}


# ----------------------------------------------------------------------------- #
# calculate likelihood for agumented data based on current parameters in the Gaussian mixture model
ll.cal<-function(design.mat,
                 unique.design.mat,
                 signal1,mean1,var1,
                 signal2,mean2,var2,
                 covar,
                 posterior.prob,
                 prob.phenotype,
                 alpha,
                 nCNV){
  signal.mat.ll<-cbind(signal1,signal2)
  ll.list<-posterior.prob
  alpha.list<-rep(alpha,length(posterior.prob)/nCNV)
  for (i in 1:nrow(unique.design.mat)){
    tmp.seq<-which(apply(design.mat, 1, function(x) return(all(x == unique.design.mat[i,]))))
    tmp.cov.mat<-matrix(c(unique(var1[tmp.seq]),unique(covar[tmp.seq]),
                          unique(covar[tmp.seq]),unique(var2[tmp.seq])),nrow = 2)
    ll.list[tmp.seq]<-dmvnorm(signal.mat.ll[tmp.seq,],  mu = c(unique(mean1[tmp.seq]),unique(mean2[tmp.seq])), Sigma = tmp.cov.mat)*
      prob.phenotype[tmp.seq]* #phenotype model
      alpha.list[tmp.seq]
    
  }
  
  return(ll.list)
  
}

# ll.cal(design.mat=design.signal,unique.design.mat=unique(design.signal),
#   signal1=signal.x1,mean1 = mean.signal1,var1=var.signal1,
#   signal2=signal.x2,mean2 = mean.signal2,var2=var.signal2,
#   posterior.prob=pis,nCNV=3)


# first part, EM algorithm for one location
# ============================================================================================ #
# fit gaussian mixture model for one CNVR (copy number variation region)

OSA.2D.ECM<-function(
  signal.mat,           # n by k matrix with sample names
  phenotype,
  nCNV=NULL,
  Dim.reduction="PCA", # dimension reduction technique, PCA or UMAP
  dt.nCNV=1,           # the number of PCs used to determine the best number of CNVs.1: 1 PC, 2: 2 PCs
  assumption="binary", #
  label.initial="hclust", # initial labels by hierarchical clustering "hclust" or "GMM"
  max.clusters=5,
  covariates=NULL,
  covariates.signal.names=NULL,# covariates  to be adjusted in signal model
  covariates.pheno.names=NULL,  # covariates to be adjusted in phenotype model
  max.iter=200,    # maximum number of iterations
  tol=1,          # convergence threshold
  # ldf to get signal with larger difference among clusters
  applyldf=FALSE,
  post=NULL,
  # warm start: set initial values
  a1=NULL,
  d1=NULL,
  a2=NULL,
  d2=NULL,
  initial.pis=NULL
){
  
  if (is.null(colnames(signal.mat)) ){
    warning("probe names in signal.mat is missing, assigning probe names")
    colnames(signal.mat)<-paste0("probe",seq(1,ncol(signal.mat)))
  }
  if (is.null(rownames(signal.mat)) ){
    warning("sample names in signal.mat is missing, assigning sample names")
    rownames(signal.mat)<-paste0("sampele",seq(1,nrow(signal.mat)))
  }
  
  # sample size
  n      <- nrow(signal.mat)
  
  # two dimension signal
  if (applyldf==FALSE){
    
    # PCA
    if(Dim.reduction=="PCA"){
      tmp.signal.mat<-scale(signal.mat,center = TRUE,scale = TRUE)
      log2R.pca<-prcomp(tmp.signal.mat)
      signal.dred<-log2R.pca[["x"]][,c(1,2)]
    }
    # UMAP
    if(Dim.reduction=="UMAP"){
      signal.dred<-umap(signal.mat)[["layout"]]
    }
    
  }else if (isTRUE(applyldf)){
    signal.dred<-apply.ldf(signal.mat,post)
  }
  
  rownames(signal.dred)<-rownames(signal.mat)
  colnames(signal.dred)<-c("D1","D2")
  
  if(is.null(nCNV)){
      if (dt.nCNV==1){
        tmp.op<-Optimal_Clusters_GMM(data.frame(signal.dred[,1]),max_clusters = max.clusters,criterion = "BIC",
                                     dist_mode="eucl_dist",verbose = FALSE,plot_data = FALSE)
      }
      if (dt.nCNV==2){
        tmp.op<-Optimal_Clusters_GMM(data.frame(signal.dred),max_clusters = max.clusters,criterion = "BIC",
                                     dist_mode="eucl_dist",verbose = FALSE,plot_data = FALSE)
      }
      nCNV<-which.min(tmp.op)
  }

  # initial labels
  initial.pis<-c()
  if(Dim.reduction=="UMAP"){
    hcluster<-hclust(dist(signal.dred),"ave")
    initial.pis<-cutree(hcluster,k=nCNV)
  }
  if (is.null(initial.pis) |Dim.reduction=="PCA"){
    # initial labels using GMM method
    if (label.initial=="GMM"){
      tmpll.list<-c()
      for (tmp.seed1 in 1:100) {
        tmp.label<-GMM(data.frame(signal.dred[,1]),gaussian_comps =nCNV,dist_mode = "eucl_dist",verbose = FALSE,seed=tmp.seed1)$Log_likelihood
        tmpll.list<-c(tmpll.list,sum(apply(tmp.label,1,max)))
      }
      tmp.seed2<-which.max(tmpll.list)
      tmp.label<-GMM(data.frame(signal.dred[,1]),gaussian_comps =nCNV,dist_mode = "eucl_dist",verbose = FALSE,seed=tmp.seed2)$Log_likelihood
      initial.pis<-apply(tmp.label,1,which.max)
      names(initial.pis)<-names(signal.dred[,1])
    }
    # initial labels using hierarchical clustering
    if (label.initial=="hclust"){
      hcluster<-hclust(dist(signal.dred),"ave")
      initial.pis<-cutree(hcluster,k=nCNV)
      names(initial.pis)<-names(signal.dred[,1])
    }
    
  }
  if ( min(table(initial.pis))==1){
    nCNV<-nCNV-1
    
    # initial labels using GMM method
    if (label.initial=="GMM"){
      tmpll.list<-c()
      for (tmp.seed1 in 1:100) {
        tmp.label<-GMM(data.frame(signal.dred[,1]),gaussian_comps =nCNV,dist_mode = "eucl_dist",verbose = FALSE,seed=tmp.seed1)$Log_likelihood
        tmpll.list<-c(tmpll.list,sum(apply(tmp.label,1,max)))
      }
      tmp.seed2<-which.max(tmpll.list)
      tmp.label<-GMM(data.frame(signal.dred[,1]),gaussian_comps =nCNV,dist_mode = "eucl_dist",verbose = FALSE,seed=tmp.seed2)$Log_likelihood
      initial.pis<-apply(tmp.label,1,which.max)
      names(initial.pis)<-names(signal.dred[,1])
    }
    
    # initial labels using hierarchical clustering
    if (label.initial=="hclust"){
      hcluster<-hclust(dist(signal.dred),"ave")
      initial.pis<-cutree(hcluster,k=nCNV)
      names(initial.pis)<-names(signal.dred[,1])
    }
  }
  
  # define CNV states according to the raw intensity, 5 states,0,1,2,3,4
  tmp.center<-rep(0,nCNV)
  for (c in 1:nCNV) {
    tmp.center[c]<-mean(rowMeans(signal.mat)[which(initial.pis==c)])
  }
  center.dat<-data.frame(raw.center=tmp.center,cluster=seq(1,nCNV))
  center.dat<-center.dat[order(center.dat$raw.center),]
  tmp.center<-sort(tmp.center)
  
  center.dat$sort.cluster<-seq(1,nCNV)
  CNV<-center.dat$sort.cluster+(2-which.min(abs(tmp.center)))
  tmp.seq<-data.frame(raw.val=center.dat$raw.center,cluster=center.dat$cluster,CNV=CNV)
  cur.val<-data.frame(raw.val=center.dat$raw.center[match(initial.pis,center.dat$cluster)],sample=names(initial.pis))
  tmp.seq<-merge(cur.val,tmp.seq,by="raw.val")
  tmp.seq<-tmp.seq[match(rownames(signal.mat),tmp.seq$sample),]
  
  # data augmentation
  data<-ExpandData.2D(trait=phenotype,
                      subID=rownames(signal.dred),
                      signal=signal.dred,
                      CNVs=CNV,
                      nCNV=nCNV,
                      covariates=covariates)
  
  
  # initiate poesterior probabilities
  pis<-c()
  for (i in 1:n){
    tmp.p1                         <- rep(0,nCNV)
    tmp.p1[match(tmp.seq$CNV[i],CNV)]  <- 1
    pis                            <- c(pis,tmp.p1)
  }
  data$pis=pis
  
  
  # create design matrix for signal models
  s.vec  <- data$CNV.state
  # dummy variables created for CNV states
  design.signal<-dummy.var(s.vec,nCNV)
  if ( !is.null(covariates.signal.names) ){
    design.signal<-cbind(design.signal,data[,covariates.signal.names])
    colnames(design.signal)[(ncol(design.signal)+1-length(covariates.signal.names)):ncol(design.signal)]<-covariates.signal.names
    design.signal<-as.matrix(design.signal)
  }
  
  
  # create design matrix for phenotype model
  phe.CNV.state <- data$CNV.state
  phe.CNV.state[which(data$CNV.state==2)] <- 0
  phe.CNV.state[which(data$CNV.state!=2)] <- 1
  design.phe<-cbind(intercept=rep(1,nrow(data)),phe.CNV.state=phe.CNV.state)
  if ( !is.null(covariates.pheno.names) ){
    design.phe<-cbind(design.phe,data[,covariates.pheno.names])
    colnames(design.phe)[(ncol(design.phe)+1-length(covariates.pheno.names)):ncol(design.phe)]<-covariates.pheno.names
    design.phe<-as.matrix(design.phe)
  }
  
  
  
  if(is.null(a1)){
    a1<-rep(1,ncol(design.signal))
    a2<-rep(1,ncol(design.signal))
    d1<-rep(1,ncol(design.signal))
    d2<-rep(1,ncol(design.signal))
  }
  
  # signal model
  # design.signal<-cbind(design.signal,tmp.phenotype.y)
  signal.x1 <- data$signal.D1
  signal.x2 <- data$signal.D2
  trait    <- data$trait
  # phenotype variable
  tmp.phenotype.y<-data$trait
  
  distance.list <- c()
  alpha         <- rep(1/nCNV,nCNV)
  num_iter      <- 0
  log.ll.int    <- log(0)
  alpha.int      <- rep(1/nCNV,nCNV)
  
  repeat{
    num_iter<-num_iter+1
    # print(paste0("iteration=",num_iter))
    # -------------------------------------------------------------------------------------- #
    # -------------------------------------------------------------------------------------- #
    # M STEP
    # fit signal mean model
    a1<-glm(signal.x1~0+design.signal,family = "gaussian"(link = "identity"),weights = as.vector(pis/exp(design.signal%*%d1)))$coefficients
    a2<-glm(signal.x2~0+design.signal,family = "gaussian"(link = "identity"),weights = as.vector(pis/exp(design.signal%*%d2)))$coefficients
    # updated mean for nCNV components
    mean.signal1   <- as.vector(design.signal%*%a1) #
    mean.signal2   <- as.vector(design.signal%*%a2) #
    # fit signal variance model, update signal for each component
    tmp.var.mod.y1 <- (signal.x1-(design.signal%*%a1))**2
    d1             <- glm(tmp.var.mod.y1 ~ 0+design.signal, family = Gamma(link = "log"),weights = pis)$coefficients
    var.signal1     <- exp(design.signal%*%d1) #
    tmp.var.mod.y2 <- (signal.x2-(design.signal%*%a2))**2
    d2             <- glm(tmp.var.mod.y2 ~ 0+design.signal, family = Gamma(link = "log"),weights = pis)$coefficients
    var.signal2     <- exp(design.signal%*%d2) #
    # updated correlation within each level of covariates
    cov.list<-cov.cal(design.mat=design.signal,unique.design.mat=unique(design.signal),
                      signal1=signal.x1,mean1 = mean.signal1,
                      signal2=signal.x2,mean2 = mean.signal2,
                      weights=pis)
    # -------------------------------------------------------------------------------------- #
    # fit phenotype model:assume constant effect from del-normal-dup
    phe.model<-glm(tmp.phenotype.y ~ design.phe, family = binomial(link = "logit"),weights = pis)
    prob.phe <-predict(phe.model,type = "response")

    # -------------------------------------------------------------------------------------- #
    # fit copy number model
    for (j in 1:nCNV){
      alpha[j] <- sum(pis[seq(0,n-1)*nCNV+j])/sum(pis)
    }
    
    # -------------------------------------------------------------------------------------- #
    # -------------------------------------------------------------------------------------- #
    # E STEP
    # update posterior probability pis
    ll.list<-ll.cal(design.mat=design.signal,unique.design.mat=unique(design.signal),
                    signal1=signal.x1,mean1 = mean.signal1,var1=var.signal1,
                    signal2=signal.x2,mean2 = mean.signal2,var2=var.signal2,
                    covar=cov.list,
                    posterior.prob=pis,prob.phenotype=prob.phe,alpha=alpha,nCNV=nCNV)
    # update pis
    new.pis<-pis
    log.ll <-0
    for (k in 1:n) {
      tmp.likelihood<-ll.list[seq(1,nCNV)+(k-1)*nCNV]                                                                                      #copy number model
      new.pis[seq(1,nCNV)+(k-1)*nCNV]<-tmp.likelihood/sum(tmp.likelihood)
      # log.ll<-log.ll+sum(log(tmp.likelihood))
      log.ll<-log.ll+sum(log(tmp.likelihood+10**(-10)))
    }
    
    distance.list  <- c(distance.list,log.ll)
    
    # use latest estimated parameters
    pis            <- new.pis
    alpha.int      <- alpha
    mean.signal.int1<-mean.signal1
    mean.signal.int2<-mean.signal2
    
    if (abs(log.ll-log.ll.int)<tol | num_iter>max.iter ){break}
    # aim: keep the second last pis if last.dist.diff>=0
    log.ll.int     <- log.ll
  }
  # update when only one iteration run
  if(num_iter==1){
    pis            <- new.pis
    log.ll.int     <- log.ll
    alpha.int      <- alpha
    mean.signal.int1<-mean.signal1
    mean.signal.int2<-mean.signal2
  }
  
  
  
  # ========================================================================= #
  # relocate copy number states
  new.cnv<-c()
  for (i in 1:(length(pis)/nCNV)){
    new.cnv<-c(new.cnv,which.max(pis[seq(1,nCNV)+(i-1)*nCNV]))
  }
  tmp.center2<-rep(0,nCNV)
  for (c in 1:nCNV) {
    tmp.center2[c]<-mean(rowMeans(signal.mat)[new.cnv==c])
  }
  center.dat1<-data.frame(raw.center=tmp.center2,cluster=seq(1,nCNV))
  center.dat1<-center.dat1[order(center.dat1$raw.center),]
  center.dat1$new.cluster<-seq(1,nCNV)
  tmp.center2<-sort(tmp.center2)
  center.dat1$new.cluster<-center.dat1$new.cluster+(2-which.min(abs(tmp.center2)))
  center.dat1<-center.dat1[order(center.dat1$cluster),]
  
  phe.CNV.state<-rep(center.dat1$new.cluster,n)
  # association test after convergence 
  if (assumption=="constant"){
    if ( !is.null(covariates.pheno.names) ){
      tmp.design<-cbind(phe.CNV.state,data[,covariates.pheno.names])
      colnames(tmp.design)[(ncol(tmp.design)+1-length(covariates.pheno.names)):ncol(tmp.design)]<-covariates.pheno.names
      tmp.design<-as.matrix(tmp.design)
    } else 
    {tmp.design<-phe.CNV.state}
    p1<-glm(tmp.phenotype.y ~ tmp.design, family = binomial(link = "logit"),weights = pis)
    p.val<-summary(p1)$coefficients[2,"Pr(>|z|)"]
    beta<-summary(p1)$coefficients[2,"Estimate"]
    std <-summary(p1)$coefficients[2,"Std. Error"]
  }
  
  if (assumption=="binary"){
    phe.CNV.state.bin<-rep(1,length(phe.CNV.state))
    phe.CNV.state.bin[which(phe.CNV.state==2)]<-0
    
    if ( !is.null(covariates.pheno.names) ){
      tmp.design<-cbind(phe.CNV.state.bin,data[,covariates.pheno.names])
      colnames(tmp.design)[(ncol(tmp.design)+1-length(covariates.pheno.names)):ncol(tmp.design)]<-covariates.pheno.names
      tmp.design<-as.matrix(tmp.design)
    }else 
    {tmp.design<-phe.CNV.state.bin}
    p1<-glm(tmp.phenotype.y ~ tmp.design, family = binomial(link = "logit"),weights = pis)
    p.val<-summary(p1)$coefficients[2,"Pr(>|z|)"]
    beta<-summary(p1)$coefficients[2,"Estimate"]
    std <-summary(p1)$coefficients[2,"Std. Error"]
    # beta<-summary(p1)$coefficients["tmp.designphe.CNV.state.bin","Estimate"]
  }
  
  if (assumption=="deletion"){
    # table(phe.CNV.state)
    phe.CNV.state.del<-rep(0,length(phe.CNV.state))
    phe.CNV.state.del[which(phe.CNV.state<=1)]<-1
    if (any(phe.CNV.state.del==1)){
      if ( !is.null(covariates.pheno.names) ){
        tmp.design<-cbind(phe.CNV.state.del,data[,covariates.pheno.names])
        colnames(tmp.design)[(ncol(tmp.design)+1-length(covariates.pheno.names)):ncol(tmp.design)]<-covariates.pheno.names
        tmp.design<-as.matrix(tmp.design)
      }else 
      {tmp.design<-phe.CNV.state.del}
      p1<-glm(tmp.phenotype.y ~ tmp.design, family = binomial(link = "logit"),weights = pis)
      p.val<-summary(p1)$coefficients[2,"Pr(>|z|)"]
      beta<-summary(p1)$coefficients[2,"Estimate"]
      std <-summary(p1)$coefficients[2,"Std. Error"]
      # beta<-summary(p1)$coefficients["tmp.designphe.CNV.state.bin","Estimate"]
    }else {
      p.val<-1
      beta<-0
      std <-100
    }
    
  }
  
  # ========================================================================= #
  # to see the dominant CN in this probe
  center.dat1<-center.dat1[order(center.dat1$cluster),]
  CN.proportion<-c(0,1,0)
  CN.proportion[1]<-sum(alpha.int[which(center.dat1$new.cluster<2)])
  CN.proportion[2]<-sum(alpha.int[which(center.dat1$new.cluster==2)])
  CN.proportion[3]<-sum(alpha.int[which(center.dat1$new.cluster>2)])
  names(CN.proportion)<-c("deletion","normal","duplication")
  
  if (!is.null(covariates)){
    tmp.centroids<-data.frame()
    tmp.center3<-rep(1,nCNV)
    raw.rowmean<-rowMeans(signal.mat)
    unique.centroids.design.mat<-as.matrix(unique(covariates[,covariates.signal.names]))
    for (i in 1:nrow(unique.centroids.design.mat) ) {
      tmp.seq<-which(apply(as.matrix(covariates[,covariates.signal.names]), 1, 
                           function(x) return(all(x == unique.centroids.design.mat[i,]))))
      for (c in 1:nCNV) {
        tmp.center3[c]<-mean(raw.rowmean[which(new.cnv[tmp.seq]==c)])
      }
      tmp.centroids<-rbind(tmp.centroids,tmp.center3)
    }
    tmp.centroids<-cbind(unique.centroids.design.mat,tmp.centroids)
    colnames(tmp.centroids)<-c(covariates.signal.names,paste0("CNV",center.dat1$new.cluster))
  }else{
    tmp.centroids<-center.dat1
  }
  
  
  # final cluster staus
  pp.mat        <- matrix(new.pis,ncol = nCNV,byrow = T)
  clusters      <- apply(pp.mat, 1, which.max)
  # <2 indicate deletion; =2 indicate normal; >2 indicate duplication
  final.cluster <- center.dat1$new.cluster[clusters]
  
  
  return(list(weights         = alpha.int,
              CN.cluster      = center.dat1$new.cluster,
              raw.overall.mean= center.dat1$raw.center,
              sample.cluster  = final.cluster,
              centroids       = tmp.centroids,
              # covariance.mat  = cov.signal,
              posterior.prob  = new.pis,
              beta            = beta,
              std             = std,
              p.val           = p.val,
              nCNV            = nCNV,
              CN.proportion   = CN.proportion,
              signal.dred     = signal.dred,
              joint.log.ll    = log.ll.int,
              distance.list   = distance.list))
  
}





# fit gaussian mixture models for all CNVRs
Known.region.test<-function(position,
                            signal.mat,
                            phenotype,
                            Dim.reduction="PCA", # dimension reduction technique, PCA or UMAP
                            dt.nCNV=1,           # the number of PCs used to determine the best number of CNVs.1: 1 PC, 2: 2 PCs
                            assumption="binary",
                            label.initial="hclust", # initial labels by hierarchical clustering "hclust" (recommended for correlated PCs) or  "GMM" 
                            nCNV=NULL,
                            max.clusters=5,
                            covariates=NULL,
                            covariates.signal.names=NULL,# covariates  to be adjusted in phenotype model
                            covariates.pheno.names=NULL,  # covariates to be adjusted in phenotype model
                            overall.alpha=0.1,
                            applyldf=FALSE,
                            smooth=FALSE,
                            max.iter=200,
                            tol=0.01
){
  # ---------------------------------------------------------------------------#
  
  # ---------------------------------------------------------------------------#
  # apply 2D GMM based on known CNV regions
  step2.res.list<-list()
  step3.res.list<-list()
  f<-1
  k<-1
  # h<-1
  for (h in 1:nrow(position)) {
    print(paste0("Analyzing # ",h, " CNV"))
    if (smooth){
      tmp.signal.mat.pre<-signal.mat[,position[h,1]:position[h,2]]
      tmp.signal.mat<-modSaRa::smooth(tmp.signal.mat.pre)
      rownames(tmp.signal.mat)<-rownames(tmp.signal.mat.pre)
      colnames(tmp.signal.mat)<-colnames(tmp.signal.mat.pre)
    }
    if (!smooth){
      tmp.signal.mat<-signal.mat[,position[h,1]:position[h,2]]
    }
    
    # go to next iteration if nCNV is 1
    # PCA
    if(Dim.reduction=="PCA"){
      log2R.pca<-prcomp(tmp.signal.mat,center=TRUE,scale.=TRUE)
      tmp.signal.dred2d<-log2R.pca[["x"]][,c(1,2)]
      colnames(tmp.signal.dred2d)<-c("D1","D2")
    }
    # UMAP
    if(Dim.reduction=="UMAP"){
      tmp.signal.dred2d<-umap(tmp.signal.mat)[["layout"]]
      colnames(tmp.signal.dred2d)<-c("D1","D2")
    }
    
    
    # determine number of clusters
    # # GMM method
    if (is.null(nCNV) ){
      if (dt.nCNV==1){
        tmp.op<-Optimal_Clusters_GMM(data.frame(tmp.signal.dred2d[,1]),max_clusters = max.clusters,criterion = "BIC",
                                     dist_mode="eucl_dist",verbose = FALSE,plot_data = FALSE)
      }
      if (dt.nCNV==2){
        tmp.op<-Optimal_Clusters_GMM(data.frame(tmp.signal.dred2d),max_clusters = max.clusters,criterion = "BIC",
                                     dist_mode="eucl_dist",verbose = FALSE,plot_data = FALSE)
      }
      tmp.nCNV<-which.min(tmp.op)
    }
    
    if (!is.null(nCNV) ){
      tmp.nCNV<-nCNV
    }

    
    # ==================================================================================== #
    # directly give the true number of clusters
    tmp.2D.OSA<-list(p.val=1,beta=NA,posterior.prob=NA,CN.proportion=NA,nCNV=tmp.nCNV,centroids=NA,
                     signal.dred=tmp.signal.dred2d,sample.cluster=NA)
    if(tmp.nCNV!=1){
      # if error occurs, exit current loop
      tryCatch(tmp.2D.OSA<-OSA.2D.ECM(signal.mat = tmp.signal.mat,
                                      phenotype               = phenotype,
                                      label.initial           = label.initial,
                                      dt.nCNV                 = dt.nCNV,
                                      covariates              = covariates,
                                      covariates.pheno.names        = covariates.pheno.names,
                                      covariates.signal.names = covariates.signal.names,
                                      nCNV                    = tmp.nCNV, 
                                      Dim.reduction           = Dim.reduction, 
                                      max.iter=max.iter, tol=tol,assumption=assumption, max.clusters=max.clusters),
               error=function(e){ 
                 # if error occurs, exit current loop 
                 
               } 
      )
      
      if (isTRUE(applyldf)){
        post<-tmp.2D.OSA$posterior.prob[tmp.nCNV*c(0:(nrow(tmp.signal.dred2d)-1))+1]
        for (j in 2:tmp.nCNV){
          post<-cbind(post,tmp.2D.OSA$posterior.prob[tmp.nCNV*c(0:(nrow(tmp.signal.dred2d)-1))+j])
        }
        # post.prob<-apply(post, 1, max)
        singal.ldf<-apply.ldf(tmp.signal.mat,post)
        if (dt.nCNV==1){
          tmp.op<-Optimal_Clusters_GMM(data.frame(singal.ldf[,1]),max_clusters = max.clusters,criterion = "BIC",
                                       dist_mode="eucl_dist",verbose = FALSE,plot_data = FALSE)
        }
        if (dt.nCNV==2){
          tmp.op<-Optimal_Clusters_GMM(data.frame(singal.ldf),max_clusters = max.clusters,criterion = "BIC",
                                       dist_mode="eucl_dist",verbose = FALSE,plot_data = FALSE)
        }
        
        tmp.nCNV<-which.min(tmp.op)
        
        tryCatch(tmp.2D.OSA<-OSA.2D.ECM(signal.mat = tmp.signal.mat,
                                        phenotype               = phenotype,
                                        label.initial           = label.initial,
                                        dt.nCNV                 = dt.nCNV,
                                        covariates              = covariates,
                                        covariates.pheno.names        = covariates.pheno.names,
                                        covariates.signal.names = covariates.signal.names,
                                        nCNV = tmp.nCNV,  
                                        post=post,
                                        applyldf=TRUE,
                                        Dim.reduction=Dim.reduction, 
                                        max.iter=max.iter, tol=tol,assumption=assumption, max.clusters=max.clusters),
                 error=function(e){ 
                   # if error occurs, exit current loop 
                 } 
        )
        # end apply ldf
      }
      
      step2.res.list[[f]]<-list(CNVRind=h,
                                region.length=position[h,2]-position[h,1],
                                range=c(start=position[h,1],
                                        end=position[h,2]),
                                nCNV=tmp.nCNV,
                                signal.dred     = tmp.2D.OSA$signal.dred,
                                sample.cluster  = tmp.2D.OSA$sample.cluster,
                                posterior.prob  = tmp.2D.OSA$posterior.prob,
                                CNV.prop        = tmp.2D.OSA$CN.proportion,
                                centroids       = tmp.2D.OSA$centroids,
                                overall.p.val   = tmp.2D.OSA$p.val,
                                overall.beta    = tmp.2D.OSA$beta,
                                overall.std     = tmp.2D.OSA$std)
      f<-f+1
      
      # store significant results
      if (tmp.2D.OSA$p.val<=overall.alpha){
        step3.res.list[[k]]<-list(CNVRind=h,
                                  region.length=position[h,2] - position[h,1],
                                  range=c(start=position[h,1],  end=position[h,2]),
                                  nCNV=tmp.nCNV,
                                  signal.dred     = tmp.2D.OSA$signal.dred,
                                  sample.cluster  = tmp.2D.OSA$sample.cluster,
                                  posterior.prob  = tmp.2D.OSA$posterior.prob,
                                  CNV.prop        = tmp.2D.OSA$CN.proportion,
                                  centroids       = tmp.2D.OSA$centroids,
                                  overall.p.val   = tmp.2D.OSA$p.val,
                                  overall.beta    = tmp.2D.OSA$beta,
                                  overall.std     = tmp.2D.OSA$std)
        k<-k+1
      }
      
      
      
    }
    
    # end step 2 loop
  }
  # end step 2 loop
  
  return(list(candidate.CNR.test=step2.res.list,
              sig.CNR.test=step3.res.list))
}

# ------------------------------------------- #
# apply ldf
# posterior is called using posterior from GMM with PCA
apply.ldf<-function (full.signal, posterior) 
{
  # posterior <- posterior[dimnames(full.signal)[[1]], ]
  # if (sum(dimnames(full.signal[[1]]) != dimnames(posterior)[[1]]) > 
  #     0) 
  #   stop("Names are not matching")
  can <- cancor(x = full.signal, y = posterior)
  ldf <- full.signal %*% can$xcoef[, 1:2]
  # ldf <- as.numeric(full.signal %*% can$xcoef[, 1:2])
  ldf.signal<-apply(ldf,2,function(x){x/sd(x)})
  # sd(ldf.signal[,1])
  return(ldf.signal)
}

