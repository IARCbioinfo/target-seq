nb_mut_thr <- function(nbmut){
  fit_NB_rob=glmrob.nb(nbmut,rep(100,length(nbmut)))$coef
  x2=0:max(nbmut)
  #Create a pdf with histogram
  pdf("Hist_nb_mut.pdf",6.5,5)
  hist(nbmut,col="grey",br=100,freq=F,main="")
  lines(x2,dnbinom(x2,size = 1/fit_NB_rob["sigma"],mu=fit_NB_rob["slope"]*100),lwd=2,col="red")
  threshold_rob=qnbinom(1-0.05/length(nbmut),size = 1/fit_NB_rob["sigma"],mu=fit_NB_rob["slope"]*100)
  abline(v=threshold_rob,col="red",lwd=2)
  dev.off()
  return(threshold_rob)
}
