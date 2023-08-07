#library("ggplot2")
#library("gridExtra")
# plot clusters based on 1st and 2nd PCs
OSCAA.plot<-function(OSA.PCA.res){
  
  clusters.num<-OSA.PCA.res$sample.cluster
  CNV<-c("double deletion","single deletion","normal","single duplication","double duplication")[clusters.num+1]
  clusters.color<-c("springgreen4","darkgreen","red","gold4","brown")[clusters.num+1]
  
  
  # p1<-ggplot(as.data.frame(OSA.PCA.res[["signal.dred"]]),aes(x=D1, y=D2)) + geom_point() +
  #   ggtitle("raw data without clustering")
  
  p2<-ggplot(as.data.frame(OSA.PCA.res[["signal.dred"]]),aes(x=D1, y=D2,color=CNV)) + geom_point() +
    scale_color_manual(values = c("double deletion" = "pink3", "single deletion" = "darkgreen",
                                  "normal" = "firebrick1", 
                                  "single duplication" = "gold4", "double duplication" = "brown")[1+sort(unique(clusters.num))]) +
    ggtitle(paste0(" length: ",OSA.PCA.res[["region.length"]],
                   "; postion: ",OSA.PCA.res[["range"]][1], "-",OSA.PCA.res[["range"]][2],
                   "; p val: ", round(OSA.PCA.res[["overall.p.val"]],4)))
  print(p2)
  # grid.arrange(p1, p2, nrow = 1)
  
}


OSCAA.plot(OSA.PCA.res=CNVR.res[["sig.CNR.test"]][[1]])
