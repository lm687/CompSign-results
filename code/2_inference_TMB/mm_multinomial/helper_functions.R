plot_pairs_with_identity = function(x, ...){
  #' plots scatterplots for any two columns of x and adds a identity line to 
  #' compare the values and not simply the correlation
  par(mfrow=c(ncol(x),ncol(x)), mar=c(0.1,0.1,0.1,0.1))
  for(i in 1:ncol(x)){
    for(j in 1:ncol(x)){
      if(j > i){
        plot(x[,i], x[,j], xaxt='n',yaxt='n', ann=FALSE, ...)
        abline(coef=c(0,1))
      }else if(i==j){
        plot(0, xlim=c(0, 2), ylim=c(0, 2), col='white'); text(1,1, labels = colnames(x)[i])
      }else{
        plot.new()
        # plot(0, xaxt='n',yaxt='n', ann=FALSE, ...)
      }
    }
  }
}