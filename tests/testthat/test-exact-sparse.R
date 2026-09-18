test_that("sparse exact is opt-in and preserves the unregularized objective", {
  expect_identical(eihrms_default_params()$ot_method, "exact")
  p <- eihrms_default_params(); p$ot_method <- "exact_sparse"
  expect_silent(validate_params(p))
  expect_true(distance_method_is_symmetric("ppm_wasserstein", "exact_sparse"))
  set.seed(908)
  for (width in c(.001, 45, 5000, 3e6)) for (i in 1:15) {
    a <- cbind(sample(c(50,100,200),20,TRUE)*(1+runif(20)*min(width,1e5)/1e6),rexp(20))
    b <- cbind(sample(c(50,100,200),17,TRUE)*(1+runif(17)*min(width,1e5)/1e6),rexp(17))
    a[1,2] <- 0
    d <- ppm_wasserstein_distance(a,b,ppm=width/3,ot_method="exact_sparse")
    expect_equal(d,ppm_wasserstein_distance(a,b,ppm=width/3),tolerance=1e-10)
    expect_equal(d,ppm_wasserstein_distance(b,a,ppm=width/3,ot_method="exact_sparse"),tolerance=1e-10)
  }
  # Capped transport need not follow the monotone 1D coupling.
  a <- cbind(c(100,100.01),c(.5,.5)); b <- cbind(c(100.01,100.02),c(.5,.5))
  expect_equal(ppm_wasserstein_distance(a,b,ot_method="exact_sparse"),.5)
})

test_that("sparse exact preserves boundary, alignment and invalid-input handling", {
  for(w in c(1,45,5000)) for(sign in c(-1,1)) for(eps in c(-1e-12,0,1e-12)) {
    y <- 100*(1+sign*w/2e6)/(1-sign*w/2e6)*(1+eps)
    a <- cbind(c(100,100),c(.3,.7)); b <- cbind(c(y,300),c(.8,.2))
    for(align in c(FALSE,TRUE)) expect_equal(
      ppm_wasserstein_distance(a,b,ppm=w/3,align=align,ot_method="exact_sparse"),
      ppm_wasserstein_distance(a,b,ppm=w/3,align=align),tolerance=1e-10)
  }
  a <- cbind(c(100,100.001),c(1e-10,1))
  for(b in list(matrix(numeric(),ncol=2),cbind(100,0),cbind(100,NA_real_),cbind(-1,1),cbind(Inf,1))) {
    expect_identical(ppm_wasserstein_distance(a,b,ot_method="exact_sparse"),
                     ppm_wasserstein_distance(a,b))
  }
})

test_that("sparse diagnostics distinguish analytical and validated solver components", {
  a <- cbind(c(100,100.001),c(.2,.8)); b <- cbind(c(100.0005,100.002),c(.7,.3))
  collector <- new.env()
  d <- ppm_wasserstein_distance(a,b,ot_method="exact_sparse",.diagnostics=collector)
  r <- collector$records[[1]]
  expect_identical(r$requested_method,"exact_sparse")
  expect_identical(r$selected_path,"exact_sparse_primary")
  expect_false(r$fallback_used)
  expect_true(all(r$attempts$accepted))
  expect_equal(r$sparse_details$general_components,1)
  expect_equal(r$selected_total,d)
  bad <- function(...) data.frame(from=1,to=1,mass=0)
  expect_error(ppm_wasserstein_distance(a,b,ot_method="exact_sparse",
               .exact_solver=bad,.diagnostics=collector),"failed plan validation")
  expect_identical(tail(collector$records,1)[[1]]$status,"error")
  broken <- function(...) stop("injected solver failure")
  expect_error(ppm_wasserstein_distance(a,b,ot_method="exact_sparse",
               .exact_solver=broken,.diagnostics=collector),"failed plan validation")
  expect_identical(tail(collector$records,1)[[1]]$original_solver_error,"injected solver failure")
  ppm_wasserstein_distance(cbind(100,1),cbind(100.001,1),ot_method="exact_sparse",.diagnostics=collector)
  expect_equal(tail(collector$records,1)[[1]]$sparse_details$analytic_components,1)
})

test_that("matrix, detailed and parallel routes accept the sparse option", {
  f <- list(a=cbind(c(100,150),c(.4,.6)),b=cbind(c(100.001,150.004),c(.7,.3)),c=cbind(c(100.002,200),c(.3,.7)))
  l <- lapply(f,function(x) {x[,1]<-x[,1]-50;x})
  p <- eihrms_default_params(); p$distance_method <- "ppm_wasserstein"; p$tol_ppm <- 15
  expected <- compute_distance_matrix(f,l,p,progress=FALSE)
  expected_search <- compute_distance_matrix_search(f,l,f,l,p,progress=FALSE)
  p$ot_method <- "exact_sparse"
  expect_equal(compute_distance_matrix(f,l,p,progress=FALSE),expected,tolerance=1e-10)
  expect_equal(compute_distance_matrix_search(f,l,f,l,p,progress=FALSE),expected_search,tolerance=1e-10)
  p$use_split_loss <- TRUE; p$use_typical_loss <- TRUE
  q <- p; q$ot_method <- "exact"
  extra <- list(loss_anchor_list=l,loss_pair_list=l,loss_anchor_typ_list=l,loss_pair_typ_list=l)
  expect_equal(do.call(compute_distance_matrix,c(list(f,l,p,progress=FALSE),extra)),
               do.call(compute_distance_matrix,c(list(f,l,q,progress=FALSE),extra)),tolerance=1e-10)
  p$return_distance_components <- TRUE; q$return_distance_components <- TRUE
  expect_equal(do.call(compute_distance_matrix,c(list(f,l,p,progress=FALSE),extra)),
               do.call(compute_distance_matrix,c(list(f,l,q,progress=FALSE),extra)),tolerance=1e-10)
  if (.Platform$OS.type != "windows") {
    p$use_parallel <- TRUE; p$n_cores <- 2L
    expect_equal(do.call(compute_distance_matrix,c(list(f,l,p,progress=FALSE),extra)),
                 do.call(compute_distance_matrix,c(list(f,l,q,progress=FALSE),extra)),tolerance=1e-10)
  }
})
