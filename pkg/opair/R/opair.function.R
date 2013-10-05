
opair <-
  function(descriptors, products, d.equiv=.5, conf.level=0.95,
           abbreviate.names=FALSE, ...) 
{
  ## Should descriptors and products be data.frames/matrices/variables
  ## or names of variables to be interpreted in a data(.frame)
  ## argument?

  ## Assuming descriptors is a data.frame and products is a
  ## variable/factor with default level equal to the reference; one
  ## reference product and one-or-more test-products; at least two
  ## levels

  ## Various checks:
  stopifnot(is.factor(products) && nlevels(products) >= 2)
  stopifnot(is.data.frame(descriptors))
  ## conf.level checked in confint.clm

  ## Compute various variables:
  nprod <- nlevels(products)
  ndscr <- ncol(descriptors)
  testprodnames <- levels(products)[-1]
  attr.names <- names(descriptors)
  d.equiv <- as.numeric(d.equiv)[1]
  abbreviate.names <- as.logical(abbreviate.names)[1]

  ## Ensure all descriptor/attribute variables are ordered factors: 
  attr.df <- data.frame(lapply(descriptors, num2ord))
  ## Compute CLMs:
  fit.list <- lapply(attr.df, function(y) {
    eval.parent(clm(y ~ products, link="probit")) })
  ## Extract list of d-primes:
  dp.list <- lapply(fit.list, function(fl) {
    b <- fl$beta
    names(b) <- testprodnames
    b * sqrt(2)
  })
  ## Extract summmary of CLM objects:
  summ.list <- lapply(fit.list, summary)
  ## Extract list of standard errors of d-primes:
  se.list <- lapply(summ.list, function(s) {
    b <- coef(s)
    nalpha <- length(s$alpha)
    b[-(1:nalpha), 2] * sqrt(2)
  })
  ## Extract list of CIs:
  ci.list <- lapply(summ.list, function(s) {
    b <- coef(s)
    nalpha <- length(s$alpha)
    b <- b[-(1:nalpha), 1:2, drop=FALSE] * sqrt(2)
    a <- 1 - conf.level
    z <- qnorm(1 - a/2)
    ## b[, 1] + b[, 2] %*% t(c(-z, z))
    cbind(b[, 1] - z * b[, 2],
          b[, 1] + z * b[, 2])
  })
  ## ci.list <- lapply(fit.list, function(fl) {
  ##   ci <- confint(fl, type="Wald", level=conf.level)
  ##   from <- nrow(ci) - nprod + 2
  ##   ci[from:nrow(ci), ] * sqrt(2)
  ## })
  ci.mat <- do.call(rbind, ci.list)
  ## Extract p-values for difference tests:
  pdiff.list <- lapply(summ.list, function(s) {
    b <- coef(s)
    nalpha <- length(s$alpha)
    b[-(1:nalpha), 4]
  })
  ## Compute p-values for equivalence tests:
  pmat.list <- lapply(summ.list, function(s) {
    b <- coef(s)
    nalpha <- length(s$alpha)
    z1 <- (b[-(1:nalpha), 1] + d.equiv/sqrt(2)) / b[-(1:nalpha), 2]
    z2 <- (b[-(1:nalpha), 1] - d.equiv/sqrt(2)) / b[-(1:nalpha), 2]
    cbind("low"=pnorm(z1, lower.tail=FALSE),
          "upr"=pnorm(z2, lower.tail=TRUE))
  })
  pequiv <- do.call(c, lapply(pmat.list, function(pm)  {
    apply(pm, 1, max) }))

  ## Optionally abbreviate Descriptor and Product names:
  abbr.attr.names <-
    if(abbreviate.names && any(nchar(attr.names) > 10))
      abbreviate(attr.names, minlength=10) else attr.names
  abbr.testprodnames <-
    if(abbreviate.names && any(nchar(testprodnames) > 7))
      abbreviate(testprodnames, minlength=7) else testprodnames
  ## Compute coefficient matrix (data.frame):
  coefmat <- data.frame("Descriptor"=rep(abbr.attr.names, each=nprod-1),
                        "Product"=rep(abbr.testprodnames,
                          length(attr.names)),
                        "d-prime"=unlist(dp.list),
                        "std.err"=unlist(se.list),
                        "lower"=ci.mat[, 1],
                        "upper"=ci.mat[, 2],
                        "p-diff"=do.call(c, pdiff.list),
                        "p-equiv"=pequiv)
  rownames(coefmat) <- NULL
  res <- list(coefficients=coefmat)

  ## check convergence
  res$convergence <- sapply(fit.list, function(f) f$convergence)
  ## More thorough convergence tests:
  ## maxabsgrad <- sapply(fit.list, function(f) max(abs(f$gradient)))
  ## step <- sapply(fit.list, function(f)
  ##                max(abs(solve(f$Hessian, f$gradient))))
  ## condHess <- sapply(summ.list, function(f) f$condHess)

  ## Assign class and return object:
  class(res) <- "opair"
  res
}

print.opair <-
  function(x, digits=max(3, getOption("digits") - 3), ...)
{
  ## Print title line:
  cat("\nThurstonian model for ordinal paired comparisons\n")
  ## Format and print table of coefficients:
  cat("\nd-prime estimates:\n")
  tab <- x$coefficients
  ## Format d.prime correctly:
  tab$d.prime <- formatC(tab$d.prime, format="f", digits=digits)
  ## Remove std.err from printed table:
  tab <- tab[, -4]
  ## Print table:
  print.data.frame(tab, digits=digits, quote=FALSE,
                   row.names=FALSE)
  ## Return object invisibly:
  invisible(x)
}

num2ord <- function(x) {
  ## Result is an ordered factor with numerical levels. If x is an
  ## ordered factor, the resulting factor might have a different
  ## ordering. 

  ## Convert x to numeric if it isn't already:
  x <- if(is.factor(x)) as.numeric(as.character(x)) else as.numeric(x)
  stopifnot(is.numeric(x)) ## Cannot convert x to numeric
  ux <- sort(unique(x), decreasing=FALSE)
  factor(x, levels=ux, ordered=TRUE)
}

plot.opair <- function(x, type=1, ...) {
  type <- round(as.numeric(type))
  if(!type %in% 1:2) stop("'type' has to be one of (1, 2)")
  b <- coef(x)
  d <- b$d.prime
  y <- range(d)
  ## Make better xlimits:
  drange <- range(c(floor(y * 2) / 2, ceiling(y * 2) / 2 ))
  if(type == 1) { ## dot-plot:
    dotchart(b$d.prime,
             labels=paste(b$Descriptor, ", ", b$Product, sep=""),
             xlim=drange, color=c("blue", "red")[(d > 0) + 1],
             xlab="d-prime", ...) 
    abline(v=0)
  }
  if(type == 2) {
    plot(d, type="h", col=c("blue", "red")[(d > 0) + 1],
         frame.plot=FALSE, ylim=drange, axes=FALSE, xlab="",
         ylab="d-prime", ...)
    axis(1, at=1:length(d),
         labels=as.character(1:length(d)))
    ## axis(1, at=1:length(dest), labels=as.character(1:length(dest)))
    axis(2, las=1)
    abline(h=0)
  }
  invisible(x)
}

save.opair <-
  function(fit, file=stop("'file' must be specified"), sep=",",
           dec=".", ...)
{
  name <- paste0(file, ".csv")
  write.table(coef(fit), file=name, quote=FALSE,
              sep=sep, row.names=FALSE, dec=dec, ...)
}
