#' PLmixed: A package for estimating GLMMs with factor structures.
#'
#' The \code{PLmixed} package's main function is \code{\link{PLmixed}}, which estimates
#' the model through nested maximizations using the \pkg{\link{lme4}} and \pkg{\link{optimx}} packages
#' (previously the \code{\link{optim}} function). This extends the capabilities of \pkg{\link{lme4}}
#' to allow for estimated factor structures, making it useful for estimating multilevel
#' factor analysis and item response theory models with an arbitrary number of hierarchical
#' levels or crossed random effects.
#' @docType package
#' @name PLmixed-package
#' @references
#' Rockwood, N. J., & Jeon, M. (2019). Estimating complex measurement and growth models
#' using the R package PLmixed.\emph{Multivariate Behavioral Research, 54}(2), 288-306.
#'
#' Jeon, M., & Rabe-Hesketh, S. (2012). Profile-likelihood approach for estimating
#' generalized linear mixed models with factor structures. \emph{Journal of Educational
#' and Behavioral Statistics, 37}(4), 518-542.
NULL


#' Fit GLMM with Factor Structure
#'
#' Fit a (generalized) linear mixed effects model (GLMM) with factor structures. Utilizes both the
#' \pkg{\link{lme4}} package and \code{\link{optim}} function for estimation using a profile-likelihood based
#' approach.
#' @param formula A formula following that of \pkg{\link{lme4}}, with the addition that factors can be specified
#' as random effects. Factor names should not be names of variables in the data set, and are instead
#' defined with the \code{factor} argument.
#' @seealso \pkg{\link{lme4}}
#' @param data A data frame containing the variables used in the model (but not factor names).
#' @param family A GLM family, see \code{\link{glm}} and \code{\link{family}}.
#' @param load.var A variable in the dataframe identifying what the factors load onto. Each unique element in \code{load.var} will have
#' a unique factor loading. All rows in the dataset with the same value for \code{load.var} will have the same factor loading.
#' @param lambda A matrix or list of matrices corresponding to the loading matrices. A value of NA
#' indicates the loading is freely estimated, while a numeric entry indicates a constraint.
#' @param factor A list of factors corresponding to the loading matrices and factors specified in model.
#' @param init A scalar (default = \code{1}) or vector of initial lambda values. If a scalar, the value is applied to all lambda parameters.
#' If a vector, the values apply in row by column by matrix order.
#' @param nlp A character vector containing the names of additional nonlinear parameters that are in the model formula.
#' @param init.nlp A scalar (default = \code{1}) or vector of initial nlp values. If a scalar, the value is applied to all nlp parameters.
#' If a vector, the values apply in the order listed.
#' @param nAGQ If family is non-gaussian, the number of points per axis for evaluating the adaptive
#' Gauss-Hermite approximation to the log-likelihood. Defaults to \code{1}, corresponding to the Laplace approximation.
#' See \code{\link{glmer}}.
#' @seealso \code{\link{glmer}}
#' @seealso \code{\link{lmer}}
#' @param method The \code{\link{optimx}} optimization method. Defaults to \code{L-BFGS-B}.
#' @param lower Lower bound on lambda parameters if applicable.
#' @param upper Upper bound on lambda parameters if applicable.
#' @param lme4.optimizer The \pkg{\link{lme4}} optimization method.
#' @param lme4.start Start values used for \pkg{\link{lme4}}.
#' @param lme4.optCtrl A list controlling the lme4 optimization. See \code{\link{lmerControl}}
#' or \code{\link{glmerControl}}
#' @param opt.control Controls for the \code{\link{optimx}} optimization.
#' @param REML Use REML if model is linear? Defaults to \code{FALSE}.
#' @param SE Method of calculating standard errors for fixed effects.
#' @param ND.method Method of calculating numerical derivatives.
#' @param check Check number of observations vs. levels and number of observations vd. random effects.
#' @param est Return parameter estimates.
#' @param iter.count Print the iteration counter during optimization.
#' @return An object of class \code{PLmod}, which contains an object of class \code{merMod} as one of its elements.
#' Some functions for class \code{merMod} have been adapted to work with class \code{PLmod}. Others can be utilized
#' using \code{object$'lme4 Model'}, where \code{object} is an object of class \code{PLmod}.
#' @details Factors are listed within the \code{formula} in the same way that random effects are specified
#' in \pkg{\link{lme4}}. The grouping variable listed after \code{|} defines what the factor values randomly
#' vary over, just as \code{|} does for other random effects. The names of factors and other random
#' effect terms can be listed within the same set of parentheses, allowing the covariance between the
#' factor(s) and random effect(s) to be estimated. The same factor may be specified for multiple grouping
#' variables, allowing for multilevel or crossed effects.
#'
#' The \code{factor} argument must list any factor that appears in the \code{formula}. The ordering will
#' depend on the ordering of the matrices listed within \code{lambda}. The matrices in \code{lambda}
#' specify the factor loading matrices. The number of matrices in \code{lambda} should equal the number
#' of character vectors in \code{factor} and the number of elements in \code{load.var}. The number of
#' rows in the \emph{k}th matrix listed in \code{lambda} should correspond to the number of unique elements
#' in the dataset for the \emph{k}th variable listed in \code{load.var}, and the number of columns in the \emph{k}th
#' matrix should correspond to the number of factors listed in the \emph{k}th character vector of \code{factor}.
#'
#' Within the \emph{k}th matrix, the \emph{(i, j)} cell corresponds to the factor loading for the \emph{i}th unique element
#' of the \emph{k}th variable listed in \code{load.var} on the \emph{j}th factor listed in the \emph{k}th character vector
#' of \code{factor}. Each element of the matrix should be either a number or \code{NA}. If the element is a
#' number, the loading will be constrained to that value. If the element is an \code{NA}, the loading will
#' be freely estimated. For identification, it is necessary (but not sufficient) for at least one element in
#' each column to be constrained.
#'
#' The \code{nlp} argument can be viewed as a special case of the \code{factor} argument, where the character vector
#' listed in \code{nlp} is automatically linked to 1 x p lambda matrix, where p is the number of elements in \code{nlp}.
#' The \code{load.var} for these parameters is viewed as a constant, so that the \code{nlp} parameters are equivalent for
#' all rows in the dataset. Thus, \code{nlp} simplifies the process of adding additional nonlinear parameters to the model
#' without having to specify corresponding \code{lambda} and \code{load.var} values.
#'
#' @references
#'
#' Rockwood, N. J., & Jeon, M. (2019). Estimating complex measurement and growth models
#' using the R package PLmixed.\emph{Multivariate Behavioral Research, 54}(2), 288-306.
#'
#' Jeon, M., & Rabe-Hesketh, S. (2012). Profile-likelihood approach for estimating
#' generalized linear mixed models with factor structures. \emph{Journal of Educational
#' and Behavioral Statistics, 37}(4), 518-542.
#'
#' @keywords GLMM GLLAMM IRT Factor
#' @export
#' @import lme4 numDeriv Matrix optimx
#' @importFrom stats df.residual gaussian logLik optim vcov AIC BIC coef deviance family fitted
#' nobs pnorm predict residuals simulate
#' @importFrom graphics par plot points
#' @examples
#' data("IRTsim") # Load the IRTsim data
#'
#' IRTsub <- IRTsim[IRTsim$item < 4, ] # Select items 1-3
#' set.seed(12345)
#' IRTsub <- IRTsub[sample(nrow(IRTsub), 300), ] # Randomly sample 300 responses
#'
#' IRTsub <- IRTsub[order(IRTsub$item), ] # Order by item
#' irt.lam = c(1, NA, NA) # Specify the lambda matrix
#'
#' # Below, the # in front of family = binomial can be removed to change the response distribution
#' # to binomial, where the default link function is logit.
#'
#' irt.model <- PLmixed(y ~ 0 + as.factor(item) + (0 + abil.sid |sid) +(0 + abil.sid |school),
#'                      data = IRTsub, load.var = c("item"), # family = binomial,
#'                      factor = list(c("abil.sid")), lambda = list(irt.lam))
#' summary(irt.model)
#'
#' \dontrun{
#' # A more time-consuming example.
#' # ~ 5-10 minutes
#'
#' data("KYPSsim") # Load the KYPSsim data
#'
#' kyps.lam <- rbind(c( 1,  0),  # Specify the lambda matrix
#'                   c(NA,  0),
#'                   c(NA,  1),
#'                   c(NA, NA))
#'
#' kyps.model <- PLmixed(esteem ~ as.factor(time) +  (0 + hs | hid)
#'                       + (0 + ms | mid) + (1 | sid), data = KYPSsim,
#'                       factor = list(c("ms", "hs")), load.var = c("time"),
#'                       lambda = list(kyps.lam))
#' summary(kyps.model)
#'
#' data("JUDGEsim")
#' JUDGEsim <- JUDGEsim[order(JUDGEsim$item), ] # Order by item
#' unique(JUDGEsim$item)
#'
#' # Specify Lambda matrix
#' judge.lam <- rbind(c( 1,  0,  1,  0,  0,  0),
#'                    c(NA,  0, NA,  0,  0,  0),
#'                    c(NA,  0, NA,  0,  0,  0),
#'                    c( 0,  1,  0,  1,  0,  0),
#'                    c( 0, NA,  0, NA,  0,  0),
#'                    c( 0, NA,  0, NA,  0,  0),
#'                    c( 0,  0,  0,  0,  1,  0),
#'                    c( 0,  0,  0,  0, NA,  0),
#'                    c( 0,  0,  0,  0, NA,  0),
#'                    c( 0,  0,  0,  0,  0,  1),
#'                    c( 0,  0,  0,  0,  0, NA),
#'                    c( 0,  0,  0,  0,  0, NA))
#'
#' # Conduct analysis
#' judge.example <- PLmixed(response ~ 0 + as.factor(item) + (1 | class)
#'                          + (0 + trait1.t + trait2.t + trait1.s + trait2.s | stu)
#'                          + (0 + teacher1 + teacher2 | tch), data = JUDGEsim,
#'                          lambda = list(judge.lam), load.var = "item",
#'                          factor = list(c("teacher1", "teacher2", "trait1.t",
#'                                          "trait2.t", "trait1.s", "trait2.s")))
#'
#' summary(judge.example)
#'
#' data("KYPSitemsim")
#'
#' time.lam <- rbind(c( 1,  0),  # Specify time lambda matrix
#'                   c(NA,  0),
#'                   c(NA,  1),
#'                   c(NA, NA))
#'
#' item.lam <- c(1, NA, NA, NA, NA, NA) # Specify item lambda matrix
#'
#' KYPSitemsim$time2 <- (KYPSitemsim$time == 2) * 1
#' KYPSitemsim$time3 <- (KYPSitemsim$time == 3) * 1
#' KYPSitemsim$time4 <- (KYPSitemsim$time == 4) * 1
#'
#' kyps.item.model <- PLmixed(response ~ 0 + as.factor(item) + lat.var:time2
#'                            + lat.var:time3 + lat.var:time4 + (0 + hs:lat.var | hid)
#'                            + (0 + ms:lat.var | mid) + (0 + lat.var:as.factor(time) | id),
#'                            data = KYPSitemsim, lambda = list(time.lam, item.lam),
#'                            factor = list(c("ms", "hs"), "lat.var"),
#'                            load.var = c("time", "item"))
#'
#' summary(kyps.item.model)
#'
#' }
#'

################## PLmixed Function ######################

PLmixed <- function(formula, data, family = gaussian, load.var = NULL, lambda = NULL, factor = NULL, init = 1,
                    nlp = NULL, init.nlp = 1, nAGQ = 1, method = "L-BFGS-B", lower = -Inf, upper = Inf,
                    lme4.optimizer = "bobyqa", lme4.start = NULL, lme4.optCtrl = list(), opt.control = NULL,
                    REML = FALSE, SE = 1, ND.method = "simple", check = "stop", est = TRUE, iter.count = TRUE) {

  start.time <- proc.time()
  iter.counter.global <- 0
  lme4.initial.theta.values <- lme4.start
  stop.iter.counter <- NULL
  model <- formula

  df.name <- deparse(substitute(data))
  if(inherits(data, "tbl_df")) data <- as.data.frame(data)

  # For Iteration Summary
  ll.list.global <- NULL
  lambda.est.global <- NULL
  fixed.effect.est.global <- NULL
  random.effect.est.global <- NULL
  time.global <- NULL

  new.lambda <- vector("list", length(load.var))
  new.factor <- vector("list", length(load.var))

  if(is.null(lambda)){
    lambda <- list(NA)
  }

  if(length(load.var) > length(lambda)){
    warning("More load.var listed than lambda matrices.",
            immediate. = TRUE)
    for (i in ((length(lambda) + 1):length(load.var))){
      lambda[i] = NA
    }
  }

  if(is.null(factor)){
    factor <- list(NA)
  }

  factors <- rep(NA, length(load.var))
  col.vec <- rep(NA, length(load.var))
  uniq.list <- vector("list", length(load.var))
  num.vec <- rep(NA, length(load.var))

  tot.est <- 0

  if (!is.null(load.var)){
    for (l in 1:length(load.var)){

      # Determine column of load.var, which are unique, and number of unique

      col <- which(colnames(data) == load.var[l])
      col.vec[l] <- col
      uniq <- unique(data[, col])
      uniq.list[[l]] <- uniq
      num <- length(uniq)
      num.vec[l] <- num


      # If no lambda is given, and no factor names are given,
      # Lambda will be a 1 by q matrix where q is the number
      # of unique values in load.var.
      # Lambda[1,1] is constrained to 1.

      # If no lambda is given, but factor names are given,
      # Lambda will be a f by q matrix were f is the number
      # of factor names given, and q is the number of unique
      # values in load.var.
      # For each factor, one value in Lambda is constrained
      # to 1.

      # The exception is if the load.var contains only 1 unique element.
      # This allows for constraints and nonlinear fixed effect parameters
      # to be included in the model using a single element lambda matrix
      # corresponding to a load.var with 1 unique element and a nonlinear parameter
      # name given in factor, which is included in the model formula.

      if(is.na(lambda[l]) & num > 1){
        tmp.lambda <- NULL
      }
      else{
        tmp.lambda <- lambda[[l]]
      }

      if(is.na(factor[l])){
        tmp.factor <- NULL
      }
      else{
        tmp.factor <- factor[[l]]
      }

      if(is.null(tmp.lambda)){
        if(is.null(tmp.factor)){
          warning("For load.var = ", paste(load.var[l]), "\nNo lambda and factor arguments specified.\n",
                  "Lambda will be a 1 x q matrix where q = length(unique(load.var)).\n",
                  "lambda[1,1] is constrained to 1.",
                  immediate. = TRUE)
          tmp.lambda <- matrix(NA, nrow = num, ncol = 1)
          tmp.lambda[1, 1] <- 1
        }
        else{
          warning("For load.var = ", paste(load.var[l]), "\nNo lambda matrix provided.\nLambda will be ",
                  "an f x q matrix where f is the number of factor names provided and q = length(unique(load.var)).\n",
                  "One element of each column in lambda is constrained to 1.",
                  immediate. = TRUE)
          tmp.lambda <- matrix(NA, nrow = num, ncol = length(tmp.factor))
          for(i in 1:length(tmp.factor)){
            tmp.lambda[i, i] <- 1
          }
        }
      }



      num.est.par <- sum(is.na(tmp.lambda))
      tot.est <- c(tot.est, num.est.par)


      # Determines the number of factors specified. If factor
      # names are provided, the number of factors is set equal to
      # the number of names provided.
      # If no names are provided, the number of factors is set to the
      # number of columns in Lambda. Factor names are assigned column
      # by column from lambda to lambda as f1.1, f2.1, f1.2, f2.2, etc.

      if (!is.null(tmp.factor)){
        factors[l] <- length(tmp.factor)
      }
      else{
        warning("For load.var = ", paste(load.var[l]), "\nNo factor names provided.\n",
                "Factors named as fc.v, where c is the column number of the vth load.var.",
                immediate. = TRUE)
        factors[l] <- ncol(as.matrix(tmp.lambda))
        tmp.factor <- paste("f", 1:factors[l],".",l, sep="")
      }

      tmp.lambda <- as.matrix(tmp.lambda)
      new.lambda[[l]] <- tmp.lambda
      new.factor[[l]] <- tmp.factor

    }
  }


  # Determine how many parameters in Lambda must be estimated.
  # For those that must be estimated, start values are assigned
  # based on those provided by the user. If none are provided,
  # the start values used are 1.

  tot.est <- sum(tot.est)

  if (length(init) == 1){
    init.tot <- rep(init, tot.est)
  }
  else{
    init.tot <- c(init, rep(1, (tot.est - length(init))))
  }

  # Determine the number of parameters in nlp. Start values are assigned
  # based on those provided by the user. If none are provided,
  # the start values used are 1.

  tot.est.nlp <- length(nlp)

  if(tot.est.nlp > 0){
    if (length(init.nlp) == 1){
      init.tot.nlp <- rep(init.nlp, tot.est.nlp)
    }
    else{
      init.tot.nlp <- c(init.nlp, rep(1, (tot.est.nlp - length(init.nlp))))
    }
  } else {
    init.tot.nlp <- NULL
  }

  init.tot <- c(init.tot.nlp, init.tot)

  lambda <- new.lambda
  factor <- new.factor

  # Prepare family argument for fitting with lme4
  # Taken from glmer code to format the family
  # argument correctly. Will be used to determine if
  # lmer or glmer is used.
  if (is.character(family)){
    family <- get(family, mode = "function", envir = parent.frame(2))
  }
  if( is.function(family)){
    family <- family()
  }

  ################ La.load Function ########################

  La.load <- function(start, model, data, load.var, factor = NULL, consts = NULL,
                      est = FALSE, lik = TRUE, delta = FALSE, k = 1, derivs = TRUE, nAGQ = 1) {

    obs <- data
    time.local <- proc.time()
    ### num.est is used to fill in the NAs from the start values in the next loop
    num.est <- 1

    if (is.null(stop.iter.counter)){
      assign('iter.counter.global', iter.counter.global + 1, inherits = TRUE)
      if (iter.count){
        cat('\r',paste0('Iteration Number: ', iter.counter.global))
      }
    }

    if (tot.est.nlp > 0){
      nlp.vals <- start[1:tot.est.nlp]
      for (p in 1:tot.est.nlp){
        obs$place.holder.variable.for.nlp <- nlp.vals[p]
        names(obs)[names(obs) == "place.holder.variable.for.nlp"] <- nlp[p]
      }

      if (tot.est > 0){
        start <- start[(tot.est.nlp + 1): length(start)]
      }
    }

    if (!is.null(load.var)){
      for (h in 1:length(load.var)){
        num.fac <- factors[h]
        num <- num.vec[h]
        col <- col.vec[h]
        uniq <- uniq.list[[h]]
        factor.2 <- factor[[h]]
        consts.2 <- consts[[h]]



        ### Cycle through each column of the loading matrix
        for(q in 1:num.fac){

          ### Apply contstraints and lambda values to weighted.var
          for(i in 1:num){
            if(is.na(consts.2[i, q]) == 0){
              obs$weighted.var[obs[, col] == uniq[i]] <- consts.2[i, q]
            }
            else {
              obs$weighted.var[obs[, col] == uniq[i]] <- start[num.est]
              num.est <- num.est + 1
            }
          }
          names(obs)[names(obs) == "weighted.var"] <- factor.2[q]
        }

      }
    }
    #print(start)
    #print(num.est)
    #print(lme4.initial.theta.values)

    if (is.null(lme4.initial.theta.values) == 0){
      lme4.start.val <- lme4.initial.theta.values
    }
    else {
      lme4.start.val <- NULL
    }


    ### Estimate GLMM - Use lmer if family = gaussian(link = "identity") and glmer for others
    if (isTRUE(all.equal(family, gaussian()))) {
      lmer.result <- lme4::lmer(model, data = obs, REML = REML, start = lme4.start.val,
                                control = lmerControl(calc.derivs = derivs,
                                                      check.nobs.vs.nlev = check,
                                                      check.nobs.vs.nRE = check,
                                                      check.nlev.gtr.1 = check,
                                                      optimizer = lme4.optimizer,
                                                      optCtrl = lme4.optCtrl))

      # Obtain estimates to use as lme4 starting values for next optim iteration
      assign('lme4.initial.theta.values', lme4::getME(lmer.result, name = "theta"),
             inherits = TRUE)
    }
    else {
      lmer.result <- lme4::glmer(model, data = obs, family = family, start = lme4.start.val,
                                 control = glmerControl(calc.derivs = derivs,
                                                        check.nobs.vs.nlev = check,
                                                        check.nobs.vs.nRE = check,
                                                        check.nlev.gtr.1 = check,
                                                        optimizer = lme4.optimizer,
                                                        optCtrl = lme4.optCtrl),
                                                        nAGQ = nAGQ)

      # Obtain estimates to use as lme4 starting values for next optim iteration
      if (nAGQ == 0){
        assign('lme4.initial.theta.values',lme4::getME(lmer.result, name = "theta"),
               inherits = TRUE)
      }
      else{
        assign('lme4.initial.theta.values', list(fixef = lme4::fixef(lmer.result),
                                                 theta=lme4::getME(lmer.result, name = "theta")),
               inherits = TRUE)
      }
    }


    #likelihood
    ll <- as.numeric(logLik(lmer.result))
    #print(ll)

    #print(summary(lmer.result))

    time.local.2 <- proc.time()
    time.dif <- time.local.2 - time.local

    assign('ll.list.global', rbind(ll.list.global,ll), inherits = TRUE)
    assign('lambda.est.global', rbind(lambda.est.global,start), inherits = TRUE)
    assign('fixed.effect.est.global', rbind(fixed.effect.est.global,lme4::fixef(lmer.result)),
           inherits = TRUE)
    assign('random.effect.est.global', rbind(random.effect.est.global,
                                             lme4::getME(lmer.result, name = "theta")), inherits = TRUE)
    assign('time.global', rbind(time.global,(time.dif[3])), inherits = TRUE)

    if(lik == TRUE){
      return(-ll)
    }

    if(est == TRUE){
      return(list(ll=-ll, total=lmer.result, degfree = df.residual(lmer.result)))
    }

    if(delta == TRUE){
      fix.b <- fixef(lmer.result)[k]
      return(fix.b)
    }

  }

  ########## End La.Load Function ##########

  if (length(init.tot) > 0) {
    opt.result <- optimx::optimx(par = init.tot, fn = La.load, hessian = TRUE, method = method, control = opt.control,
                          consts = lambda, factor = factor, data = data, load.var = load.var, lower = lower, upper = upper,
                          model = model, est = FALSE, lik = TRUE, delta = FALSE, derivs = FALSE, nAGQ = nAGQ)
    # value 1
    Est <- coef(opt.result)[1, ]
    final.lik <- (opt.result[1, "value"]) * (-1)
    #cov.mat <- solve(opt.result$hessian)
    hess <- attr(opt.result, "details")[[3]]
    cov.mat <- solve(hess)
    st.er <- sqrt(diag(cov.mat))
    table <- cbind(Est,st.er)
    code <- opt.result[1, "convcode"]
  }
  else {
    opt.result <- optim(par = init.tot, fn = La.load, hessian = TRUE, method = "L-BFGS-B", control = opt.control,
                         consts = lambda, factor = factor, data = data, load.var = load.var, lower = lower, upper = upper,
                         model = model, est = FALSE, lik = TRUE, delta = FALSE, derivs = FALSE, nAGQ = nAGQ)
    final.lik <- (opt.result$value) * (-1)
    Est <- opt.result$par
    table <- Est
    cov.mat <- NULL
    code <- opt.result$convergence
    method <- NA
  }


  tot.est.par <- length(Est)
  total.iters <- iter.counter.global


  assign('stop.iter.counter', 1, inherits = TRUE)

  # value 4
  if (est == TRUE) {
    final.model <- La.load(start = Est, consts = lambda, data = data,
                           factor = factor, load.var = load.var, model = model,
                           est = TRUE, lik = FALSE, nAGQ = nAGQ)

    final <- final.model$total

    # Needed for SE calculations
    naive <- diag(lme4::vcov.merMod(final))
    beta <- lme4::fixef(final)
    num.fix.ef <- length(beta)

    adj.se.b <- rep(NA, num.fix.ef)

    if(SE == 1 & tot.est.par > 0){
      ### Option 1 ###
      #Direct numerical derivative
      eps <- 10^(-4)
      iden <- diag(1, nrow = tot.est.par, ncol = tot.est.par)
      beta.ep <- matrix(0, num.fix.ef, tot.est.par)

      rep.lam <- vector("list", length(lambda))
      for(i in 1:tot.est.par){
        Est.new <- Est + iden[, i]*eps

        if (!is.null(load.var)){
          if (tot.est.nlp > 0){
            lam.new <- Est.new[(tot.est.nlp + 1):length(Est.new)]
          } else {
            lam.new <- Est.new
          }

          cycle <- 1
          for(c in 1:length(lambda)){
            n.lam <- lambda[[c]]
            full.lam <- matrix(0,nrow = nrow(n.lam), ncol = ncol(n.lam))

            for(p in 1:ncol(n.lam)){
              for(q in 1:nrow(n.lam)){
                if(is.na(n.lam[q, p]) == 0){
                  full.lam[q, p] <- n.lam[q, p]
                }
                else{
                  full.lam[q, p] <- lam.new[cycle]
                  cycle <- cycle + 1
                }
              }
            }
            rep.lam[[c]] <- full.lam
          }
        }
        final.new <- La.load(start = Est.new, consts = rep.lam, data = data,
                             factor = factor, load.var = load.var, model = model,
                             est = TRUE, lik = FALSE, delta = FALSE, nAGQ = nAGQ)
        beta.ep[, i] <- lme4::fixef(final.new$total)
      }

      grad.b <- (matrix(rep(beta, tot.est.par), num.fix.ef, tot.est.par, byrow = FALSE) - beta.ep)/ eps

      #Adjusted standard errors for beta
      adj.se.b <- rep(0, num.fix.ef)
      for(i in 1:num.fix.ef){
        adj.se.b[i] <- sqrt(naive[i] + t(grad.b[i, ]) %*% cov.mat %*% grad.b[i, ])
      }

    }

    else if(SE == 2 & tot.est.par > 0){

      ### Option 2 ###
      # Use R package numDeriv for numerical derivatives

      grad.b2 <- matrix(0, num.fix.ef, tot.est.par)
      for(i in 1:num.fix.ef){
        grad.b2[i, ] <- numDeriv::grad(La.load, method = ND.method, x = Est, consts = lambda, data = data,
                                      factor = factor, load.var = load.var, model = model, lik = FALSE,
                                      est = FALSE, delta = TRUE, k = i, nAGQ = nAGQ)
      }

      # Adjust standard errors
      adj.se.b <- rep(0,num.fix.ef)
      for (i in 1:num.fix.ef){
        adj.se.b[i] <- sqrt(naive[i] + t(grad.b2[i, ]) %*% cov.mat %*% grad.b2[i, ])
      }

    }
    else if(tot.est.par == 0){
      adj.se.b <- sqrt(naive)
    }



    ### Organize Output

    cat('\r', "                                 ","\n")
    if (inherits(final, "lmerMod")){
      fix <- cbind(beta, adj.se.b, beta/adj.se.b)
      colnames(fix) <- c("Beta", "SE", "t value")
    }
    else{
      fix <- cbind(beta, adj.se.b, beta/adj.se.b,
                   2*pnorm(abs(beta/adj.se.b), lower.tail = F))
      colnames(fix) <- c("Beta", "SE", "z value", "Pr(>|z|)")
    }
    rand.ef <- lme4::VarCorr(final)

    fin.lam <- vector("list", length(load.var))

    if (tot.est.nlp > 0){
      nlp.Est <- Est[1:tot.est.nlp]
      nlp.st.er <- st.er[1:tot.est.nlp]
      nlp.tab <- cbind(nlp.Est, nlp.st.er)
      row.names(nlp.tab) <- nlp
      colnames(nlp.tab) <- c("Estimate", "SE")
      if (tot.est > 0){
        Est <- Est[(tot.est.nlp + 1): length(Est)]
        st.er <- st.er[(tot.est.nlp + 1): length(st.er)]
      } else {
        Est <- NULL
        st.er <- NULL
      }
    } else{
      nlp.Est <- NULL
      nlp.tab <- NULL
    }


    fill <- 1

    if(!is.null(load.var)){
      for (h in 1:length(load.var)){

        num.fac <- factors[h]
        num <- num.vec[h]
        col <- col.vec[h]
        uniq <- uniq.list[[h]]
        factor.2 <- factor[[h]]
        lambda.1 <- lambda[[h]]

        final.lambda <- matrix(NA, nrow = num.vec[h], ncol = 2 * (length(factor.2)))
        se.names <- rep("SE", length(factor.2))
        comb <- rbind(factor.2, se.names)
        comb <- matrix(comb, nrow = 1, ncol = 2 * length(factor.2), byrow = F)
        rownames(final.lambda) <- uniq
        colnames(final.lambda) <- comb


        for(i in 1:nrow(as.matrix(lambda.1))){
          for (j in 1:ncol(as.matrix(lambda.1))){
            if(is.na(lambda.1[i, j]) == 1){
              final.lambda[i, (j * 2 - 1)] <- Est[fill]
              final.lambda[i, j * 2] <- st.er[fill]
              fill <- fill + 1
            }
            else{
              final.lambda[i, (j * 2 - 1)] <- lambda.1[i, j]
            }
          }
        }
        fin.lam[[h]] <- final.lambda
      }


      lam.names <- paste0("lambda.", load.var)
      names(fin.lam) <- lam.names

    } else{
      fin.lam <- NULL
    }


    if (inherits(final, "lmerMod")){
      reml <- REML
    }
    else{
      reml <- FALSE
    }
    optimizers <- list("lme4 Optimizer" = lme4.optimizer, "Optim Optimizer" = method)

    finish.time <- proc.time()
    total.estimation.time <- (finish.time-start.time)/60
    total.estimation.time <- total.estimation.time[3]

    iter.summary <- list("Log-Likelihood" = ll.list.global, "Lambda" = lambda.est.global,
                         "Fixed Effects" = fixed.effect.est.global, "Random Effects" = random.effect.est.global,
                         "Time" = time.global)

    final.object <- list("Log-Likelihood" = final.lik, "Fixed Effects" = fix, "Random Effects" = rand.ef, "Lambda" = fin.lam,
                         "nlp" = nlp.tab, "Cov Matrix" = cov.mat, "Error code" = code, "Total Iterations" = total.iters,
                         "lme4"= final, "Iteration Summary" = iter.summary, "Model" = model, "Family" = family(final),
                         "Data" = df.name, "Load.Var" = load.var, "Factor" = factor, "nAGQ" = nAGQ,
                         "Lambda.raw" = new.lambda, "Param" = c(nlp.Est, Est), "Estimation Time" = total.estimation.time, "REML" = reml,
                         "Optimizer" = optimizers, "nlp.vars" = nlp)

    class(final.object) <- append("PLmod", class(final.object))
    return(final.object)

  }
  else{
    return(list("Log-Likelihood" = final.lik, "Estimates"=table, "Cov matrix"=cov.mat,
                "Error code"=code,"Optim Iterations"=total.iters))
  }

}
########## End of PLmixed Function ##########

