alive_conditional_fixture <- function() {
  list(jump_times=c(.25,.5,.75,1),kernels=NULL,rewards=cbind(c(1,0,0),c(1,1,0)),tau=1,
    query_times=c(0,.25,.5),query_state=c(1,2,3),accrued=cbind(c(0,.125,.25),c(0,.25,.4)),
    center=cbind(rep(.48,3),rep(.64,3)),scale=matrix(1,3,2),q=c(.07,.19,.33),
    side="right",grid_size=128L,max_nodes=10000L,tolerance=1,strict=FALSE,
    max_refinements=0L,exact_limit=256L,point_mass_threshold=.001,
    kernel_factors=list(base_increments=rbind(c(.12,.04,.22,.11),c(.08,.21,.09,.13),
      c(.2,.05,.12,.07),c(.15,.1,.22,.18)),relative_risk=rep(1,4),
      from=c(1,2,1,2),to=c(2,1,3,3),absorbing=3))
}

test_that("conditional scalar envelopes contain the exact joint law for every state order", {
  a <- alive_conditional_fixture()
  ref <- do.call(markov_conditional_coverage_grid,a)
  expect_equal(ref$lower,ref$upper,tolerance=1e-13)
  for (s in list(1:3,c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),3:1)) {
    z <- a
    z$rewards <- a$rewards[s,]
    z$query_state <- match(a$query_state,s)
    for (nm in c("from","to","absorbing"))
      z$kernel_factors[[nm]] <- match(a$kernel_factors[[nm]],s)
    for (screen in c(FALSE,TRUE)) {
      val <- alive_conditional_coverage(z,markov_conditional_coverage_grid,screen)
      expect_true(all(val$lower<=ref$probability+1e-14 & val$upper>=ref$probability-1e-14))
      expect_null(val$diagnostics$alive_fallback)
      if(!screen) expect_equal(val$probability,ref$probability,tolerance=1e-12)
    }
  }
})

test_that("atomic ties and inconsistent living clocks retain original query arithmetic", {
  a <- alive_conditional_fixture()
  q0 <- abs(.5-a$center[1,2])
  a$q <- c(q0-.Machine$double.eps,q0,q0+.Machine$double.eps)
  a$query_weights <- c(1.3,-.8,.5)
  for(screen in c(FALSE,TRUE)) {
    val <- alive_conditional_coverage(a,markov_conditional_coverage_grid,screen)
    expect_identical(val$diagnostics$alive_fallback,"terminal atom boundary")
    expect_identical(val$reduced,do.call(markov_conditional_coverage_grid,a)$reduced)
  }
  a$q <- c(.07,.19,.33)
  a$accrued[1,2] <- .Machine$double.eps
  val <- alive_conditional_coverage(a,markov_conditional_coverage_grid,FALSE)
  expect_identical(val$diagnostics$alive_fallback,"unsupported profile")
  expect_identical(val$reduced,do.call(markov_conditional_coverage_grid,a)$reduced)
})

test_that("signed known memberships are exact and unknown contributions remain enclosed", {
  a <- alive_conditional_fixture()
  a$query_weights <- c(-.5,.25,1.25)
  a$query_known <- matrix(NA_real_,3,3)
  a$query_known[1,2] <- 1
  ref <- do.call(markov_conditional_coverage_grid,a)$reduced
  for(screen in c(FALSE,TRUE)) {
    z <- alive_conditional_coverage(a,markov_conditional_coverage_grid,screen)$reduced
    expect_true(all(z$lower <= ref$probability & ref$probability <= z$upper))
  }
  a$query_known[,] <- 1
  z <- alive_conditional_coverage(a,markov_conditional_coverage_grid,FALSE)$reduced
  expect_identical(z$lower,z$upper)
  expect_identical(z$probability,rep(1,3))
})

test_that("public calibration selects the same first crossing with an unchanged fitted score", {
  records <- do.call(rbind,lapply(1:12,function(i) data.frame(id=paste0("t",i),
    time=c(i/24,.5+i/24,1+i/24),status=c(1,7,2))))
  fit <- life_train(LifeHist(time,status)~1,records,absorbing=2,
    status_map=c("1"=1,"7"=0,"2"=2),targets=list(home=0,alive=c(0,1)),tau=2)
  data <- data.frame(id=paste0("c",1:12),time=1/32,status=0)
  q <- seq(.015,2.015,length.out=21)
  ref <- life_calibrate(fit,data,alpha=.3,q=q)
  z <- life_calibrate(fit,data,alpha=.3,q=q,control=list(markov_occupation="condition_on_alive"))
  expect_identical(z$threshold,ref$threshold)
  expect_true(z$diagnostics$numerical$threshold_stable)
  expect_identical(z$diagnostics$numerical$occupation_engine,"condition_on_alive")
  legacy <- fit
  legacy$control$markov_occupation <- NULL
  expect_identical(life_calibrate(legacy,data,alpha=.3,q=q)$threshold,ref$threshold)
  expect_error(life_calibrate(fit,data,control=list(markov_occupation="invalid")),
               "markov_occupation")
})
