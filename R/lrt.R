#' Likelihood ratio test
#'
#' @param obj1 fit results model 1
#' @param obj2 fit results model 2
#' @return test statistic, df, and p-value
#' @export
lrt <- function(obj1,obj2) {
  stat <- 2*abs(obj1$loglik-obj2$loglik)
  df <- abs(obj1$npars-obj2$npars)
  pv <- 1 - stats::pchisq(stat,df)
  cat(paste('Likelihood ratio test:'),'\n\n')
  cat(paste('Test statistic =',round(stat,3),', df =',df,', p-value =',round(pv,3)))
}
