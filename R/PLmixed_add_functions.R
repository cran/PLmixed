#' summary.PLmod
#'
#' Obtain key output for a fitted PLmixed model object of class PLmod.
#' @rdname summary.PLmod
#' @param object an object of class PLmod
#' @return An object containing all parameter estimates and model characteristics.
#' @param ... Additional arguments.
#' @method summary PLmod
#' @export
#' @keywords GLMM GLLAMM IRT
#'

summary.PLmod <- function(object, ...){
  details <- list("nobs" = nobs(object$'lme4 Model'),
                  "ngrps" = ngrps(object$'lme4 Model'))
  if (object$'REML' == F){
    fit.stat <- list("AIC" = AIC(object$'lme4 Model') + 2*length(object$'Param'),
                     "BIC" = (BIC(object$'lme4 Model')
                              + length(object$'Param')*(log(nobs(object$'lme4 Model')))),
                     "logLik" = object$'Log-Likelihood',
                     "deviance" = deviance(object$'lme4 Model'),
                     "df.resid" = df.residual(object$'lme4 Model') - length(object$'Param'))
  }
  else{
    fit.stat <- list("REML" = lme4::REMLcrit(object$'lme4 Model')[1])
  }
  return.object <- list("Formula" = Reduce(paste, deparse(object$'Model')),
                        "Family" = family(object$'lme4 Model'),
                        "Data" = object$'Data',
                        "Fit" = fit.stat,
                        "Optim Iterations" = object$'Total Iterations',
                        "Estimation Time" = object$'Estimation Time',
                        "Lambda" = object$'Lambda',
                        "load.var" = object$'Load.Var',
                        "Random Effects" = object$'Random Effects',
                        "Fixed Effects" = object$'Fixed Effects',
                        "Details" = details,
                        "Residuals" = residuals(object),
                        "Param" = object$'Param',
                        "Optimizer" = object$'Optimizer')
  class(return.object) <- append("summary.PLmod", class(return.object))
  return.object
}

#####


#' print.summary.PLmod
#'
#' Print the output for a PLmixed model object of class PLmod.
#' @rdname print.summary.PLmod
#' @param x an object of class PLmod
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param ... Additional arguments.
#' @method print summary.PLmod
#' @export
#' @keywords GLMM GLLAMM IRT
#'

print.summary.PLmod <- function(x, digits = 4, ...){
  object <- x
  cat("Profile-based Mixed Effect Model Fit With PLmixed Using lme4 \n")
  cat("Formula: ", object$'Formula', "\n")
  cat("Data: ", object$'Data', "\n")
  .prt.family(object$'Family')
  cat("\n")
  aictab <- unlist(object$'Fit')
  .prt.aictab(aictab, digits = 2)
  cat("\n")
  .prt.resids(object$'Residuals', digits = digits)
  for (i in 1:length(object$'Lambda')){
    lam <- as.data.frame(object$'Lambda'[[i]])
    cat("Lambda: ", object$'load.var'[i], "\n")
    print(lam, digits = digits)
    cat("\n")
  }
  .prt.VC(object$'Random Effects', comp = c("Var", "Std.Dev."), digits = digits)
  .prt.grps(ngrps = object$'Details'$'ngrps', nobs = object$'Details'$'nobs')
  cat(" \nFixed effects: \n")
  print(object$'Fixed Effects', digits = digits)
  cat("\n")
  cat("lme4 Optimizer: ", unlist(object$'Optimizer'['lme4 Optimizer']), "\n")
  cat("Optim Optimizer: ", unlist(object$'Optimizer'['Optim Optimizer']), "\n")
  cat("Optim Iterations: ", object$'Optim Iterations', "\n")
  cat("Estimation Time: ", round(object$'Estimation Time', digits = 2), "minutes \n")
  return(NULL)
}


#####


#' print.PLmod
#'
#' Print the fitted PLmixed model object of class PLmod.
#' @rdname print.PLmod
#' @param x an object of class PLmod
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param ... Additional arguments.
#' @method print PLmod
#' @export
#' @keywords GLMM GLLAMM IRT
#'

print.PLmod <- function(x, digits = 4, ...){
  object <- x
  cat("Profile-based Mixed Effect Model Fit With PLmixed Using lme4 \n")
  cat("Formula: ", Reduce(paste, deparse(object$'Model')), "\n")
  cat("Data: ", object$'Data', "\n")
  .prt.family(family(object$'lme4 Model'))
  if (object$'REML' == F){
    fit.stat <- list("AIC" = AIC(object$'lme4 Model') + 2*length(object$'Param'),
                     "BIC" = (BIC(object$'lme4 Model')
                              + length(object$'Param')*(log(nobs(object$'lme4 Model')) - log(2*pi))),
                     "logLik" = object$'Log-Likelihood',
                     "deviance" = deviance(object$'lme4 Model'),
                     "df.resid" = df.residual(object$'lme4 Model') - 2*length(object$'Param'))
  }
  else{
    fit.stat <- list("REML" = REMLcrit(object$'lme4 Model')[1])
  }
  aictab <- unlist(fit.stat)
  .prt.aictab(aictab, digits = 2)
  cat("Lambda Estimates: \n")
  lam.est <- round(object$'Param', digits = digits)
  cat(lam.est, "\n")
  .prt.VC(object$'Random Effects', comp = "Std.Dev.", digits = digits)
  .prt.grps(ngrps = object$'Details'$'ngrps', nobs = object$'Details'$'nobs')
  cat(" \rFixed effects: \n")
  print(object$'Fixed Effects'[,1], digits = digits)
  return(cat(""))
}

#####

#' coef.PLmod
#'
#' Obtain coefficients for a model of class PLmod.
#' @param object an object of class PLmod
#' @param ... Additional arguments from \code{\link{coef.merMod}}.
#' @return sum of the random and fixed effects coefficients for each explanatory variable for
#' each level of the grouping factor.
#' @keywords GLMM GLLAMM IRT
#' @export

coef.PLmod <- function(object,...){
  coef(object$'lme4 Model',...)
}

#####


#' residuals.PLmod
#'
#' Obtain residuals for a model of class PLmod.
#' @param object an object of class PLmod
#' @param ... Additional arguments from \code{\link{residuals.merMod}}.
#' @keywords GLMM GLLAMM IRT
#' @export

residuals.PLmod <- function(object,...){
  residuals(object$'lme4 Model',...)
}

#####


#' ranef.PLmod
#'
#' Obtain conditional modes of the random effects for a model of class PLmod.
#' @param object an object of class PLmod
#' @param ... Additional arguments from \code{\link{ranef.merMod}}.
#' @keywords GLMM GLLAMM IRT
#' @export

ranef.PLmod <- function(object,...){
  lme4::ranef(object$'lme4 Model',...)
}

#####


#' fixef.PLmod
#'
#' Obtain fixed effect estimates for a model of class PLmod.
#' @param object an object of class PLmod
#' @param ... Additional arguments from \code{\link{fixef.merMod}}.
#' @keywords GLMM GLLAMM IRT
#' @export

fixef.PLmod <- function(object,...){
  lme4::fixef(object$'lme4 Model',...)
}

#####


#' fitted.PLmod
#'
#' Obtain fitted values for a model of class PLmod.
#' @param object an object of class PLmod
#' @param ... Additional arguments from \code{\link{fitted.merMod}}.
#' @keywords GLMM GLLAMM IRT
#' @export
#'

fitted.PLmod <- function(object,...){
  fitted(object$'lme4 Model',...)
}

#####

#' simulate.PLmod
#'
#' Simulate responses from a model of class PLmod.
#' @param object an object of class PLmod
#' @param ... Additional arguments from \code{\link{simulate.merMod}}.
#' @keywords GLMM GLLAMM IRT
#' @export
#'

simulate.PLmod <- function(object,...){
  simulate(object$'lme4 Model',...)
}

#####


#' predict.PLmod
#'
#' Predict response values from a model of class PLmod.
#' @param object an object of class PLmod
#' @param newdata data frame to obtain predictions for
#' @param ... Additional arguments from \code{\link{predict.merMod}}.
#' @keywords GLMM GLLAMM IRT
#' @export
#'

predict.PLmod <- function(object, newdata = NULL,...){

  if (is.null(newdata) == T){
    predict(object$'lme4 Model',...)
  }
  else{
    obs <- newdata
    estimates <- object$'Param'
    num.est <- 1
    for (h in 1:length(object$'Load.Var')){
      col <- which(colnames(obs) == object$'Load.Var'[h])
      uniq <- unique(obs[,col])
      num <- length(unique(obs[,col]))
      num.fac <- length(object$'Factor'[[h]])
      factor <- object$'Factor'[[h]]
      consts <- object$'Lambda.raw'[[h]]
      for(q in 1:num.fac){
        for(i in 1:num){
          if(is.na(consts[i,q]) == 0){
            obs$weighted.var[obs[,col]==uniq[i]] <- consts[i,q]
          }
          else {
            obs$weighted.var[obs[,col]==uniq[i]] <- estimates[num.est]
            num.est <- num.est + 1
          }
        }
        names(obs)[names(obs) == "weighted.var"] <- factor[q]
      }
    }
    predict(object$'lme4 Model', newdata = obs, ...)
  }
}

#####

#' plot.PLmod
#'
#' Diagnostic plots for a model of class PLmod.
#' @param x an object of class PLmod
#' @param ... Additional arguments from \code{\link{plot.merMod}}.
#' @keywords GLMM GLLAMM IRT
#' @export
#'

plot.PLmod <- function(x,...){
  object <- x
  plot(object$'lme4 Model',...)
}

#####

#' iterPlot
#'
#' Plot parameter estimates at each \code{\link{optim}} iteration.
#' @param object an object of class PLmod
#' @keywords GLMM GLLAMM IRT
#' @export
#'

iterPlot <- function(object){
  par(mfrow = c(3, 2))
  plot(1:nrow(object$'Iteration Summary'$'Lambda'),
       object$'Iteration Summary'$'Lambda'[,1],
       type = 'l', ylim = c(min(object$'Iteration Summary'$'Lambda'),
                            max(object$'Iteration Summary'$'Lambda')),
       xlab = "Iteration", ylab = "Lambda Estimate")
  for (i in 2:ncol(object$'Iteration Summary'$'Lambda')){
    points(object$'Iteration Summary'$'Lambda'[,i], type = 'l', col = i)
  }

  plot(1:nrow(object$'Iteration Summary'$'Random Effects'),
       object$'Iteration Summary'$'Random Effects'[,1],
       type = 'l', ylim = c(min(object$'Iteration Summary'$'Random Effects'),
                            max(object$'Iteration Summary'$'Random Effects')),
       xlab = "Iteration", ylab = "Random Effect Estimate")
  for (i in 2:ncol(object$'Iteration Summary'$'Random Effects')){
    points(object$'Iteration Summary'$'Random Effects'[,i], type = 'l', col = i)
  }

  plot(1:nrow(object$'Iteration Summary'$'Fixed Effects'),
       object$'Iteration Summary'$'Fixed Effects'[,1],
       type = 'l', ylim = c(min(object$'Iteration Summary'$'Fixed Effects'),
                            max(object$'Iteration Summary'$'Fixed Effects')),
       xlab = "Iteration", ylab = "Fixed Effect Estimate")
  for (i in 2:ncol(object$'Iteration Summary'$'Fixed Effects')){
    points(object$'Iteration Summary'$'Fixed Effects'[,i], type = 'l', col = i)
  }

  plot(1:nrow(object$'Iteration Summary'$'Log-Likelihood'),
       object$'Iteration Summary'$'Log-Likelihood', type = 'l',
       xlab = "Iteration", ylab = "Log-Likelihood")

  plot(1:nrow(object$'Iteration Summary'$'Time'),
       object$'Iteration Summary'$'Time', type = 'l', ylim = c(0, max(object$'Iteration Summary'$'Time')),
       xlab = "Iteration", ylab = "Iteration Time (s)")

  predictions <- predict(object$'lme4 Model')
  residuals <- residuals(object$'lme4 Model')
  plot(predictions, residuals, cex = .5, xlab = "Fitted", ylab = "Residual")
  par(mfrow = c(1, 1))
}







#' Simulated KYPS dataset.
#'
#' A simulated dataset that replicates the dataset from KYPS.
#'
#' @format A data frame with 11494 rows and 5 variables:
#' \describe{
#'   \item{mid}{Middle School ID}
#'   \item{hid}{High School ID}
#'   \item{sid}{Student ID}
#'   \item{time}{Time Identifier}
#'   \item{esteem}{Self Esteem}
#' }
"KYPSsim"

#' Simulated KYPS item-level dataset.
#'
#' A simulated dataset that replicates the dataset item-level
#' data from KYPS.
#'
#' @format A data frame with 66947 rows and 6 variables:
#' \describe{
#'   \item{id}{Student ID}
#'   \item{time}{Time Identifier}
#'   \item{item}{Item ID}
#'   \item{mid}{Middle School ID}
#'   \item{hid}{High School ID}
#'   \item{response}{Item Response}
#' }
"KYPSitemsim"

#' Simulated multilevel IRT dataset.
#'
#' A simulated dataset that replicates the dataset from CITO.
#'
#' @format A data frame with 2500 rows and 4 variables:
#' \describe{
#'   \item{sid}{Student ID}
#'   \item{school}{School ID}
#'   \item{item}{Item ID}
#'   \item{y}{Response}
#' }
"IRTsim"


#' Simulated Multi-rater Multi-response dataset.
#'
#' A simulated dataset that replicates the dataset from a multi-rater
#' mult-reponse study where teachers and students provided responses
#' about two student traits.
#'
#' @format A data frame with 54462 rows and 7 variables:
#' \describe{
#'   \item{item}{Item ID}
#'   \item{method}{1 = teacher response, 2 = student response}
#'   \item{trait}{1 = trait 1, 2 = trait 2}
#'   \item{stu}{Student ID}
#'   \item{class}{Classroom ID}
#'   \item{tch}{Teacher ID}
#'   \item{response}{Item response}
#' }
"JUDGEsim"