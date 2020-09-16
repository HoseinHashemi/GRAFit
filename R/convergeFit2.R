convergeFit2 <- function (Data, init = Data$init, initalg = "CHARM", initspecs = list(alpha.star = 0.44), 
          finalalg = initalg, finalspecs = initspecs, inititer = 500, 
          finaliter = inititer * 4, maxruns = Inf, inititerharm = inititer, 
          inititermcmc = inititer, maxeval = inititer, ftol_rel = 0, 
          domlfit = TRUE, domcmc = TRUE, docma = FALSE, cmasigma = NULL, 
          cmasigmamult = c(2, 1, 0.5), cmatolfrac = 0.05, cmaiter = inititer, 
          cmaresetmaxruns = TRUE, maxwalltime = Inf) 
{
  tfinal = timeinmins() + maxwalltime
  overtime = FALSE
  bestLP = profitLikeModel(init, Data, makeplots = F)$LP
  converged = FALSE
  run = 0
  if (docma) {
    stopifnot(!domcmc)
    stopifnot(is.numeric(cmasigma) && (length(cmasigma) == 
                                         length(init)))
    lims = unlist(Data$intervals)
    whichfit = which(unlist(Data$tofit))
    whichlog = which(unlist(Data$tolog))
    lims = list(lower = lims[endsWith(names(lims), "lim1")], 
                upper = lims[endsWith(names(lims), "lim2")])
    for (lim in names(lims)) {
      lims[[lim]][whichlog] = log10(lims[[lim]][whichlog])
      lims[[lim]] = lims[[lim]][whichfit]
    }
  }
  if (domlfit) {
    algofuncs = c("NEWUOA", "BOBYQA", "NM")
    fits = list()
    for (func in algofuncs) {
      fits[[func]] = mlFit(Data, init = init, maxeval = maxeval, 
                           algo.func = func, maxwalltime = tfinal - timeinmins())
      if (fits[[func]]$overtime) 
        break
    }
    best = which.max(unlist(lapply(fits, function(x) {
      x$value
    })))
    fit = fits[[best]]
    init = fit$par
    bestLP = fit$value
    mlfunc = algofuncs[best]
    print(sprintf(paste0("Got bestfunc: ", mlfunc, "; value: %.6e; par:"), 
                  bestLP))
    print(init)
  }
  overtime = exceededmaxtime(timeinmins(), tfinal)
  while (!converged && !overtime) {
    if (domlfit) {
      fit = mlFit(Data, init = init, algo.func = mlfunc, 
                  maxeval = maxeval, ftol_rel = ftol_rel, maxwalltime = tfinal - 
                    timeinmins())
      newLP = fit$value
      if (newLP > bestLP) {
        init = profitRemakeModellist(fit$par, Data = Data)$parm
      }
      print(paste0("mlfit from LP=", sprintf("%.5e", bestLP), 
                   " to LP=", sprintf("%.5e", newLP), " deltaLP=", 
                   sprintf("%.3e", newLP - bestLP), "; parm:"))
      print(fit$par)
    }
    else {
      fit = mlFit(Data = Data, init = init, maxiter = inititer, 
                  algo.func = "HAR", maxwalltime = tfinal - timeinmins())
      newLP = fit$value
      if (newLP > bestLP) 
        init = fit$par
    }
    run = run + 1
    converged = (newLP - bestLP) < exp(1) || (run >= maxruns)
    if (fit$overtime) 
      break
    if (converged && !(domlfit && mlfunc == "NM")) {
      converged = FALSE
      domlfit = TRUE
      mlfunc = "NM"
    }
    if (converged) {
      if (newLP > bestLP) 
        bestLP = newLP
      print(sprintf(paste0("MLFit converged at value: %.6e; par:"), 
                    bestLP))
      print(init)
      if (docma) {
        algofunc = Data$algo.func
        Data$algo.func = "CMA"
        for (mult in cmasigmamult) {
          fit = cmaeshpc(init, profitLikeModel, Data = Data, 
                         control = list(maxit = cmaiter, fnscale = -1, 
                                        sigma = mult * cmasigma, diag.sigma = TRUE, 
                                        diag.eigen = TRUE, diag.pop = TRUE, diag.value = TRUE, 
                                        maxwalltime = tfinal - timeinmins(), trace = TRUE, 
                                        stopfitness = Inf, stop.tolx = cmatolfrac * 
                                          mult * cmasigma, lower = lims$lower, 
                                        upper = lims$upper))
          if ((fit$value - bestLP) > exp(1)) {
            converged = FALSE
            fit$par = profitRemakeModellist(fit$par, 
                                            Data = Data)$parm
            init = fit$par
            bestLP = fit$value
            print(sprintf(paste0("CMA converged at value: %.6e; par:"), 
                          bestLP))
            print(init)
            if (cmaresetmaxruns) 
              run = 0
          }
          else docma = FALSE
          overtime = exceededmaxtime(timeinmins(), tfinal)
          if (overtime) 
            break
        }
        Data$algo.func = algofunc
      }
      if (domcmc && !overtime) {
        print(paste0("LD convergence starting with params: c(", 
                     paste0(sprintf("%.6f", init), collapse = ","), 
                     ")"))
        fit = LaplacesDemon(profitLikeModel, Initial.Values = init, 
                            Data = Data, Iterations = inititerharm, Algorithm = "HARM", 
                            Thinning = 1, Specs = list(alpha.star = 0.234, 
                                                       B = NULL), CheckDataMatrixRanks = FALSE, 
                            MaxWalltime = tfinal - timeinmins())
        best = which.max(fit$Monitor[, "LP"])
        newLP = fit$Monitor[best, "LP"]
        converged = (newLP - bestLP) < exp(1)
        if (newLP > bestLP) {
          init = fit$Posterior1[best, ]
          bestLP = newLP
        }
        overtime = exceededmaxtime(timeinmins(), tfinal)
        if (overtime) 
          break
        if (converged) {
          print(paste0("LD MCMC starting with params: c(", 
                       paste0(sprintf("%.6f", init), collapse = ","), 
                       ")"))
          if (finalalg != "HARM") {
            fit = LaplacesDemon(profitLikeModel, Initial.Values = init, 
                                Data = Data, Iterations = inititermcmc, 
                                Algorithm = finalalg, Thinning = 1, Specs = finalspecs, 
                                CheckDataMatrixRanks = FALSE, MaxWalltime = tfinal - 
                                  timeinmins())
            best = which.max(fit$Monitor[, "LP"])
            newLP = fit$Monitor[best, "LP"]
          }
          converged = (newLP - bestLP) < exp(1)
          if (newLP > bestLP) {
            init = fit$Posterior1[best, ]
            bestLP = newLP
          }
          overtime = exceededmaxtime(timeinmins(), tfinal)
          if (overtime) 
            break
          if (converged) {
            fit = LaplacesDemon(profitLikeModel, Initial.Values = init, 
                                Data = Data, Iterations = finaliter, Algorithm = finalalg, 
                                Thinning = 1, Specs = finalspecs, CheckDataMatrixRanks = FALSE, 
                                MaxWalltime = tfinal - timeinmins())
            best = which.max(fit$Monitor[, "LP"])
            newLP = fit$Monitor[best, "LP"]
            converged = (newLP - bestLP) < exp(1)
            if (newLP > bestLP) {
              init = fit$Posterior1[best, ]
              bestLP = newLP
            }
            overtime = exceededmaxtime(timeinmins(), 
                                       tfinal)
            if (overtime) 
              break
          }
          print(paste0("LD MCMC finished (converged=", 
                       converged, ") with params: c(", paste0(sprintf("%.6f", 
                                                                      init), collapse = ","), ")"))
        }
      }
    }
    if (newLP > bestLP) 
      bestLP = newLP
    overtime = exceededmaxtime(timeinmins(), tfinal)
  }
  if (!exists("fit")) {
    fit = list(par = init)
    Data$algo.func = "LD"
    fit$value = profitLikeModel(fit$par, Data = Data)$LP
  }
  fit$converged = converged
  fit$overtime = overtime
  if (is.demonoid(fit) || is.laplace(fit)) {
    fit$par = init
    fit$value = bestLP
  }
  fit$par = profitRemakeModellist(fit$par, Data = Data)$parm
  return(fit)
}