#' Accurately sum a protected binary occupation-time total
#' @noRd
# Inputs are ordered tau, -time, accrued so valid occupation inputs do not
# overflow an intermediate tau + accrued before their cancellation with time.
markov_sum3 <- function(tau,time,accrued) {
  partials<-numeric(3L);n<-0L
  for(x in c(tau,-time,accrued)) {
    i<-0L
    if(n>0L)for(j in seq_len(n)) {
      y<-partials[j]
      if(abs(x)<abs(y)) { temporary<-x;x<-y;y<-temporary }
      hi<-x+y
      if(!is.finite(hi))stop('Protected total overflowed; check occupation inputs.')
      lo<-y-(hi-x)
      if(lo!=0) {i<-i+1L;partials[i]<-lo}
      x<-hi
    }
    n<-i
    if(x!=0) {n<-n+1L;partials[n]<-x}
  }
  if(n==0L)return(0)
  hi<-partials[n];n<-n-1L;lo<-0
  while(n>0L) {
    x<-hi;y<-partials[n];n<-n-1L
    hi<-x+y
    lo<-y-(hi-x)
    if(lo!=0)break
  }
  # A remaining low-order component can decide an exact halfway case. Round
  # away from the halfway endpoint only when the residual components share its
  # sign and twice the leading residual is represented by the next addition.
  if(n>0L && ((lo<0 && partials[n]<0)||(lo>0 && partials[n]>0))) {
    y<-lo*2;x<-hi+y
    if(y==x-hi)hi<-x
  }
  hi
}

#' Evaluate complete protected rewards before frozen-score centering
#' @noRd
markov_protected_total <- function(tau,time,reward,accrued) {
  stopifnot(length(tau)==1L,length(time)==1L,length(reward)==length(accrued))
  # tau/time are scalar; reward/accrued are matching component vectors. The
  # public reward contract is 0/1. Other internal numeric rewards deliberately
  # retain the existing expression, without claiming exact product rounding.
  result<-accrued
  unit<-which(reward==1)
  if(length(unit))result[unit]<-vapply(accrued[unit],function(x)
    markov_sum3(tau,time,x),numeric(1))
  other<-which(reward!=0 & reward!=1)
  if(length(other))result[other]<-accrued[other]+(tau-time)*reward[other]
  result
}
