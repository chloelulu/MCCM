
perm.fdr.adj <- function (F0, Fp) {
  ord <- order(F0, decreasing = T)
  F0 <- F0[ord]
  perm.no <- ncol(Fp)
  Fp <- as.vector(Fp)
  Fp <- Fp[!is.na(Fp)]
  Fp <- sort(c(Fp, F0), decreasing = F)
  n <- length(Fp)
  m <- length(F0)
  FPN <- (n + 1) - match(F0, Fp) - 1:m
  p.adj.fdr <- FPN / perm.no / (1:m)
  
  # Impute 0s - pseudo-ct 0.5
  p.adj.fdr[p.adj.fdr == 0] <- 0.5 / perm.no / which(p.adj.fdr == 0)
  
  p.adj.fdr <- pmin(1, rev(cummin(rev(p.adj.fdr))))[order(ord)]
  
  return(p.adj.fdr)
}

perm.fwer.adj <- function (F0, Fp) {
  ord <- order(F0, decreasing = T)
  m <- length(F0)
  F0 <- F0[ord]
  Fp <- Fp[ord, , drop=FALSE]
  col.max <- Fp[m, ]
  p.adj.fwer <- sapply(m:1, function(i) {
    x <- F0[i]
    y <- Fp[i, ]
    col.max <<- ifelse(y > col.max, y, col.max)
    # Impute 0s
    mean(c(1, col.max >= x))
  })
  p.adj.fwer <- rev(p.adj.fwer)
  
  p.adj.fwer <- pmin(1, rev(cummin(rev(p.adj.fwer))))[order(ord)]
  return(p.adj.fwer)
}

na.pad <- function (vec, ind) {
  vec0 <- numeric(length(ind))
  vec0[!ind] <- vec
  vec0[ind] <- NA
  return(vec0)
}

find.ref.z1 <- function(otu.tab, M,  ref.pct = 0.50, feature.dat.type, p.max = 500) {
  
  # p.max to reduce the computation time
  p <- nrow(otu.tab)
  n <- ncol(otu.tab)
  
  if (p > p.max) {
    ord <- order(rowSums(otu.tab), decreasing = TRUE)
    otu.tab <- otu.tab[ord[1 : p.max], ]
    p <- p.max
  } else {
    ord <- 1 : p
  }
  
  if (feature.dat.type == 'count') {
    otu.tab <- otu.tab + 1
  }
  
  if (feature.dat.type == 'proportion') {
    otu.tab <- t(apply(otu.tab, 1, function (x) {
      x[x == 0] <- min(x[x != 0]) / 2
      return(x)
    }))
  }
  
  otu.tab <- t(otu.tab)
  res <- matrix(NA, p, p)
  
  XXI <- solve(t(M) %*% M)
  dof <- n - ncol(M)
  XXIM <- XXI %*% t(M)
  
  for (i in 1 : (p - 1)) {
    
    Y <- (as.matrix(log((otu.tab[, i]) / (otu.tab[, (i + 1) : p]))))
    
    est <- XXIM %*% Y
    resid <- Y - M %*% est
    sigma <- sqrt(colSums(resid^2) / dof)
    
    res[i, (i + 1) : p] <- res[(i + 1) : p, i] <- sigma
  }
  ref.ind <- order(rowMedians(res, na.rm = TRUE))[1 : round(ref.pct * p)]
  return(ord[ref.ind])
}

bbmix.fit.MM <- function (ct, dep,  nIter = 10, winsor.qt = 1.0) {
  
  if (mean(ct == 0) >= 0.95 | sum(ct != 0) < 5) {
    stop('The number of nonzero is too small to fit the model! Consider removing it from testing or please do not enable posterior sampling!\n')
  }
  
  # Initialization
  prop0 <- ct / dep
  qt <- quantile(prop0, winsor.qt)
  prop0[prop0 >= qt] <- qt
  
  var1 <- var(prop0) / 2
  mean1 <- mean(prop0) / 2
  
  var2 <- var(prop0) / 2
  mean2 <- 3 * mean(prop0) / 2
  
  pi <- 0.5
  
  for (i in 1 : nIter) {
    
    
    shape1.1 <- ((1 - mean1) / var1 - 1 / mean1) * mean1 ^ 2
    shape1.2 <- shape1.1 * (1 / mean1 - 1)
    
    shape2.1 <- ((1 - mean2) / var2 - 1 / mean2) * mean2 ^ 2
    shape2.2 <- shape2.1 * (1 / mean2 - 1)
    
    m1 <- shape1.1 / (shape1.1 + shape1.2)
    s1 <- shape1.1 + shape1.2
    
    m2 <- shape2.1 / (shape2.1 + shape2.2)
    s2 <- shape2.1 + shape2.2
    
    f1 <- dbetabinom(ct, dep, m1, s1)
    f2 <- dbetabinom(ct, dep, m2, s2)
    
    q1 <-  pi * f1 /  (pi * f1 + (1 - pi) * f2)
    q2 <-  1 - q1
    
    pi <- mean(q1)
    
    # Rough estimation
    mean1 <- sum(prop0 * q1) / sum(q1)
    var1 <- sum((prop0 - mean1)^2 * q1) / sum(q1)
    
    mean2 <- sum(prop0 * q2) / sum(q2)
    var2 <- sum((prop0 - mean2)^2 * q2) / sum(q2)
  }
  
  return(list(shape1.1 = shape1.1, shape1.2 = shape1.2, shape2.1 = shape2.1, shape2.2 = shape2.2, pi = pi, q1 = q1))
  
}

#' @title Permutation-based differential abundance analysis
#' @param meta.dat a data frame containing the sample information
#' @param feature.dat a matrix of counts, row - features (OTUs, genes, etc) , column - samples.
#' @param feature.dat.type  the type of the feature data. It could be "count", "proportion" or "other". For "proportion" data type, posterior sampling will not be performed,
#' but the reference-based ratio approach will still be used to address compositional effects. For "other" data type, neither posterior sampling or reference-base ratio approach
#' will be used.
#' @param grp.name the name for variable of interest. It could be numeric or categorical; should be in "meta.dat".
#' @param adj.name  the name for the variable(s) to be adjusted. They could be numeric or categorical; should be in "meta.dat."
#' @param prev.filter the prevalence cutoff, under which the features will  be filtered. 
#' @param mean.abund.filter the mean relative abundance cutoff, under which the features will  be filtered. 
#' @param max.abund.filter the max relative abundance cutoff, under which the features will  be filtered. 
#' @param min.prop  proportions less than this value will be replaced with this value. Only relevant when log transformation is used. The default is 0. 
#' @param is.winsor a logical value indicating whether winsorization should be performed to replace outliers. The default is TRUE.
#' @param winsor.end a character indicating whether the outliers at the top, bottom or both will be winsorized.
#' @param outlier.pct the percentage of expected outliers. These outliers will be winsorized.
#' @param is.post.sample a logical value indicating whether to perform posterior sampling of the underlying proportions. Only relevant when the feature data are counts.
#' @param post.sample.no the number of posterior samples if posterior sampling is used.
#' @param link.func  a list of transformation functions for the feautre data.
#' @param perm.no the number of permutations; If the raw p values are of the major interest, set "perm.no" to at least 999.
#' @param strata  a factor indicating the permutation strata; permutation will be confined to each stratum.
#' @param stats.combine.func function to combine the F-statistic for the omnibus test.
#' @param ref.pct percentage of reference taxa.
#' @param stage.no the number of stages if multiple-stage normalization is used.
#' @param excl.pct the maximum percentage of significant features in the reference set that can be removed. Only relevant when multiple-stage normalization is used.
#' @param p.max the maximum number of most abundant features to be considered in reference selection; only relevant for large data sets 
#' @param is.fwer a logical value indicating whether the family-wise error rate control (West-Young) will be performed.
#' @param verbose a logical value indicating whether the trace information should be printed out.
#' @param return.feature.dat a logical value indicating whether the wisorized, filtered "feature.dat" matrix should be returned.
#' @return A list with the elements
#' \item{call}{the call}
#' \item{feature.dat}{the wisorized, filtered "feature.dat" matrix}
#' \item{meta.dat}{"meta.dat" used}
#' \item{grp.name}{"grp.name" used}
#' \item{filter.ind}{a vector of logical values indicating which features are filtered}
#' \item{R2}{a matrix of percent explained variance (number of features by number of transformation functions)}
#' \item{F0}{a matrix of observed F-statistics (number of features by number of functions)}
#' \item{RSS}{a matrix of residual sum squares (number of features by number of functions)}
#' \item{df.model, df.residual}{degrees of freedom for the model and residual space}
#' \item{p.raw}{the raw p-values based on permutations (not accurate if "perm.no" is small)}
#' \item{p.adj.fdr}{permutation-based FDR-adjusted p-values}
#' \item{p.adj.fwer}{permutation-based FWER-adjusted (West-Young) p-values}
#' @import stats
#' @import nlme
#' @import matrixStats
#' @import vegan
#' @import ape
#' @import statmod
#' @importFrom rmutil dbetabinom
#' @rdname ZicoSeq
#' @export

ZicoSeq2 <- function (meta.dat, feature.dat, grp.name, adj.name = NULL, feature.dat.type = c('count', 'proportion', 'other'),
    prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0, min.prop = 0,
    is.winsor = TRUE, outlier.pct = 0.03, winsor.end = c('top', 'bottom', 'both'),
    is.post.sample = TRUE,  post.sample.no = 25,
    link.func = list(function (x) sign(x) * (abs(x))^0.5), stats.combine.func = max,
    perm.no = 99,  strata = NULL, 
    ref.pct = 0.50, stage.no = 6, excl.pct = 0.20, p.max = 500,
    is.fwer = FALSE, verbose = TRUE, return.feature.dat = TRUE) {
  
  this.call <- match.call()
  feature.dat.type <- match.arg(feature.dat.type)
  winsor.end <- match.arg(winsor.end)
  
  if (!is.null(strata)) {
    if (length(strata) == 1 & is.character(strata)) {
      strata <- factor(meta.dat[, strata])
    } else {
      strata <- factor(strata)
    }
    if (length(strata) != nrow(meta.dat)) {
      stop("'strata' parameter is not correct. Dimension does not match!\n")
    }
  }
  
  if(is.null(rownames(feature.dat)[1])) {
    rownames(feature.dat) <- paste0("S", 1 : nrow(feature.dat))
  }
  
  sds <- matrixStats::rowSds(feature.dat)
  if (sum(sds == 0) != 0) {
    stop(paste('Feature ',  paste(which(sds == 0), collapse = ','), 
               'have identical values (e.g. all 0s)! Please remove them!\n'))
  }
  
  
  if (perm.no < 99) {
    warning('To obtain stable results, number of permutations should be at least 99!\n')
  }
  
  if (feature.dat.type == 'count') {
    if (sum(feature.dat < 0 | (feature.dat < 1 & feature.dat > 0)) != 0) {
      stop('It seems the feature data are not counts. Please check!\n')
    }
  }
  
  if (feature.dat.type == 'proportion') {
    if (sum(feature.dat < 0 | feature.dat > 1) != 0) {
      stop('It seems the feature data are not proportions. Please check!\n')
    }
  }
  
  if (feature.dat.type %in% c('proportion', 'other') & is.post.sample == TRUE) {
    cat('For proportion and other data types,  posterior sampling will not be performed!\n')
    is.post.sample <- FALSE
  }
  
  if (feature.dat.type == 'other') {
    stage.no <- 1
  }
  
  if (feature.dat.type == 'count' & nrow(meta.dat) <= 40) {
    cat('For sample size less than 40, posterior sampling will not be used!\n')
    is.post.sample <- FALSE
  }
  
  if (feature.dat.type == 'count' & is.post.sample == FALSE){
    if (min(coef(summary(lm(sqrt(colSums(feature.dat)) ~  meta.dat[, grp.name])))[-1, 'Pr(>|t|)']) < 0.05) {
      warning('The sequencing depth is correlated with the variable of interest!
							Rarefaction is needed to control for false positves!\n')
      
    }
  }
  
  ###############################################################################
  # Filter to remove very rare taxa
  if (feature.dat.type %in% c('count', 'proportion')) {
    temp <- t(t(feature.dat) / colSums(feature.dat))
    filter.ind <- rowMeans(temp != 0) >= prev.filter & rowMeans(temp) >= mean.abund.filter & matrixStats::rowMaxs(temp) >= max.abund.filter
    names(filter.ind) <- rownames(feature.dat)
    if (verbose)  cat(sum(!filter.ind), ' features are filtered!\n')
    filter.features <- rownames(feature.dat)[!filter.ind]
    feature.dat <- feature.dat[filter.ind, ]
  } else {
    filter.ind <- rep(TRUE, ncol(feature.dat))
    filter.features <- NULL
  }
  
  cat('Filter finished \n')
  
  # After filtering, some samples could have all zeros, will need to remove them.
  if (feature.dat.type %in% c('count', 'proportion')) {
    ind <- colSums(feature.dat) != 0
    if (sum(!ind) != 0) {
      warning(sum(!ind), ' samples have all zeros after filtering! They will be removed! \n')
      feature.dat <- feature.dat[, ind, drop = FALSE]
      meta.dat <- meta.dat[ind, ]
      if (!is.null(strata)) {
        strata <- strata[ind]
      }
    }
  }
  
  cat('Remove 0s finished \n')
  
  if (feature.dat.type == 'proportion') {
    # Renormalization
    feature.dat <- t(t(feature.dat) / colSums(feature.dat))
  }
  
  cat('Renormalization finished \n')
  sample.no <- ncol(feature.dat)
  otu.no <- nrow(feature.dat)
  row.names <- rownames(feature.dat)
  depth <- colSums(feature.dat)
  
  if (verbose) cat('The data has ', sample.no, ' samples and ', otu.no, ' features will be tested!\n' )
  
  if (feature.dat.type %in% c('count', 'proportion') & sum(rowSums(feature.dat != 0) <= 2) != 0) {
    warning('Some features have less than 3 nonzero values! 
               They have virtually no statistical power. Please consider removing them in the analysis!\n')
  }
  
  ###############################################################################
  # Winsorization to reduce the influence of outliers
  
  if (is.winsor == TRUE) {
    depth <- colSums(feature.dat)
    
    cat("On average, ", round(length(depth) * outlier.pct), 
        " outlier counts will be replaced for each feature!\n") 
    
    feature.dat <- apply(feature.dat, 1, function (x) {
      
      if (feature.dat.type == 'count') {
        x <- x / depth
      }
      
      if (winsor.end == 'top') {
        qt <- quantile(x, 1 - outlier.pct)
        x[x >= qt] <- qt
      }
      
      if (winsor.end == 'bottom') {
        qt <- quantile(x, outlier.pct)
        x[x <= qt] <- qt
      }
      
      if (winsor.end == 'both') {
        qt1 <- quantile(x, 1 - outlier.pct / 2)
        qt2 <- quantile(x, outlier.pct / 2)
        x[x >= qt1] <- qt1
        x[x <= qt2] <- qt2
      }
      
      if (feature.dat.type == 'count') {
        return(round(x * depth))
      } else {
        return(x)
      }				
    })
    
    feature.dat <- t(feature.dat)
    sds <- matrixStats::rowSds(feature.dat)
    if (sum(sds == 0) != 0) {
      stop(paste('After winsorization, feature ',  paste(which(sds == 0), collapse = ','),  'have identical values! Did you set the "outlier.pct" too high or some features are extremely sparse?\n'))
    }
  }
  
  
  if (feature.dat.type == 'proportion') {
    # Renormalization
    feature.dat <- t(t(feature.dat) / colSums(feature.dat))
  }
  
  ###############################################################################
  # Generate samples from posterior distribution (stacking it)
  if (is.post.sample == TRUE) {
    
    if (verbose) cat('Fitting beta mixture ...\n')
    feature.dat.p <- apply(feature.dat, 1, function (x) {
      err1 <- try(res <- bbmix.fit.MM(x, depth), silent = TRUE)
      
      # Handle error
      if (!inherits(err1, 'try-error')) {
        prop1.1 <- rbeta(sample.no * post.sample.no, shape1 = x + res$shape1.1, shape2 = res$shape1.2 + depth - x)
        prop1.2 <- rbeta(sample.no * post.sample.no, shape1 = x + res$shape2.1, shape2 = res$shape2.2 + depth - x)
        prop <- ifelse(runif(sample.no * post.sample.no) <= res$q1, prop1.1, prop1.2)
      } else {
        prop <- x / depth
        v <- var(prop)
        m <- mean(prop)
        
        a1 <- ((1 - m) / v - 1 / m) * m ^ 2
        a2 <- a1 * (1 / m - 1)
        
        if (is.na(a1) | a1 < 0) {
          # uniform prior
          prop <- rbeta(sample.no * post.sample.no, shape1 = x + 1, shape2 = otu.no + depth - x)
        } else {
          # beta prior
          prop <- rbeta(sample.no * post.sample.no, shape1 = x + a1, shape2 = a2 + depth - x)
        }
      }
      return(prop)
    })
    
    feature.dat.p <- t(feature.dat.p)
    feature.dat.p.list <- list()
    st <- 1
    end <- sample.no
    for (i in 1 : post.sample.no) {
      feature.dat.p.list[[i]] <- feature.dat.p[, st : end]
      # Normalization
      #    feature.dat.p.list[[i]] <- t(t(feature.dat.p[, st : end]) / colSums(feature.dat.p[, st : end]))
      st <- st + sample.no
      end <- end + sample.no
    }
  } else {
    feature.dat.p.list <- list()
    if (feature.dat.type == 'count') {
      feature.dat.p.list[[1]] <- t(t(feature.dat) / depth)
    } else{
      feature.dat.p.list[[1]] <- feature.dat
    }
    
    post.sample.no <- 1
  }
  
  # Replace zeros or extremely small values for log calculation
  if (feature.dat.type %in% c('count', 'proportion')) {
    for (i in 1 : post.sample.no) {
      temp <- feature.dat.p.list[[i]]
      temp[temp <= min.prop] <- min.prop
      feature.dat.p.list[[i]] <- temp
    }
  }
  #return(feature.dat.p.list)
  ###############################################################################
  # Covariate space (including intercept)
  
  
  if (is.null(adj.name)) {
    M0 <- model.matrix(~ 1, meta.dat)
  } else {
    data0 <- meta.dat[, c(adj.name), drop = FALSE]
    if (sum(is.na(data0)) != 0) {
      stop("Please remove or impute NAs in the variables to be adjusted!\n")
    }
    M0 <- model.matrix( ~ ., data0)
  }
  
  data1 <- meta.dat[, c(grp.name), drop = FALSE]
  
  if(sum(is.na(data1)) != 0) {
    stop("Please remove or impute NAs in the variable of interest!\n")
  }
  
  M1 <-  model.matrix( ~ ., data1)[, -1, drop = FALSE]  # No intercept
  
  M01 <- cbind(M0, M1)
  
  # QR decompostion
  qrX0 <- qr(M0, tol = 1e-07)
  Q0 <- qr.Q(qrX0)
  Q0 <- Q0[, 1:qrX0$rank, drop = FALSE]
  H0 <- (Q0 %*% t(Q0))
  
  qrX1 <- qr(M1, tol = 1e-07)
  Q1 <- qr.Q(qrX1)
  Q1 <- Q1[, 1:qrX1$rank, drop = FALSE]
  
  qrX01 <- qr(M01, tol = 1e-07)
  Q01 <- qr.Q(qrX01)
  Q01 <- Q01[, 1:qrX01$rank, drop = FALSE]
  
  R0 <- as.matrix(resid(lm(Q1 ~ Q0 - 1)))
  
  pX0 <- ncol(Q0)
  pX1 <- ncol(Q1)
  pX01 <- ncol(Q01)
  
  df.model <- pX01 - pX0
  df.residual <- sample.no - pX01
  
  func.no <- length(link.func)
  
  ###############################################################################
  if (feature.dat.type %in% c('count', 'proportion')) {
    if (verbose)  cat('Finding the references ...\n')
    ref.ind <- tryCatch({
       find.ref.z1(feature.dat, M0, ref.pct = ref.pct, feature.dat.type = feature.dat.type, p.max = p.max)
    }, error = function(e) {
      linda.obj <- MicrobiomeStat::linda(feature.dat, meta.dat, formula = paste0('~', grp.name, '+', paste0(adj.name, collapse = '+')))
      order(abs(linda.obj$output[[1]]$log2FoldChange), decreasing = FALSE)[1:round(ref.pct * nrow(feature.dat))]
    })
    size.factor <- colSums(feature.dat[ref.ind, ])
  } else {
    ref.ind <- NULL
  }
  ###############################################################################
  # Perform multiple stage normalization
  if (verbose) cat('Permutation testing ...\n')
  # norm.ind <- NULL
  
  for (i in 1 : stage.no) {
    # Reference proportion
    
    if (feature.dat.type == 'other') {
      divisor <- 1
    } else {
      divisor <- size.factor / depth
      divisor[divisor == 0] <- 0.5 * min(divisor[divisor != 0])
    }
    
    # Create the giant Y matrix
    Y <- matrix(NA, sample.no, func.no * otu.no * post.sample.no)
    
    # Change order - func.no * otu.no
    for (k in 1 : post.sample.no) {
      for (j in 1 : func.no) {
        func <- link.func[[j]]
        feature.dat.p <- feature.dat.p.list[[k]]
        
        Y[, (k - 1) * func.no * otu.no + func.no * (0 : (otu.no - 1)) + j] <-
          func(t(feature.dat.p) / divisor)  # No scaling first
        
      }
    }
    
    Y <- t(Y)
    TSS <- rowSums(Y^2)
    MSS01 <- rowSums((Y %*% Q01)^2)
    MSS0 <- rowSums((Y %*% Q0)^2)
    
    MSS <- (MSS01 - MSS0)
    RSS <- (TSS - MSS01)
    
    
    getPermuteMatrix <- getFromNamespace("getPermuteMatrix", "vegan")
    perm.ind <- getPermuteMatrix(perm.no, sample.no, strata = strata)
    perm.no <- nrow(perm.ind)
    
    MRSSp <- sapply(1 : perm.no, function(ii) {
      if (verbose) {
        if (ii %% 10 == 0) cat('.')
      }
      
      Rp <- R0[perm.ind[ii, ], , drop = FALSE]
      # Project to the reisdual space
      Rp <- Rp - H0 %*% Rp
      
      qrRp <- qr(Rp, tol = 1e-07)
      Q1p <- qr.Q(qrRp)
      Q1p <- Q1p[, 1:qrRp$rank, drop = FALSE]
      
      MSS01p <- MSS0 + rowSums((Y %*% Q1p)^2)
      
      MSSp <- (MSS01p - MSS0)
      RSSp <- (TSS - MSS01p)
      
      c(MSSp, RSSp)
      
    })
    
    unit <- func.no * otu.no * post.sample.no
    MSSp <- MRSSp[1 : unit, ]
    RSSp <- MRSSp[(unit + 1) : (2 * unit), ]
    
    # EB is based on the aggregated RSS
    RSS.m <- array(RSS, c(func.no, otu.no,  post.sample.no))
    RSS.m <- t(apply(RSS.m, c(1, 2), mean))  # otu.no by func.no
    
    F0.m <- array((MSS / df.model) / (RSS / df.residual), c(func.no, otu.no,  post.sample.no))
    F0.m <-  t(apply(F0.m, c(1, 2), mean))  # otu.no by func.no
    
    R2.m <- array(MSS / TSS, c(func.no, otu.no, post.sample.no))
    R2.m <- t(apply(R2.m, c(1, 2), mean))  # otu.no by func.no
    
    F0 <- (MSS / df.model)  /  (RSS  / df.residual)
    Fp <- (MSSp / df.model)  /  (RSSp  / df.residual)
    
    # Expectation of F0 and Fp
    F0 <- array(F0, c(func.no, otu.no, post.sample.no))
    Fp <- array(Fp, c(func.no, otu.no, post.sample.no, perm.no))
    
    F0 <- apply(F0, c(1, 2), mean)    # func.no * otu.no
    Fp <- apply(Fp, c(1, 2, 4), mean) # func.no * otu.no * perm.no
    
    ###############################################################################
    # Omnibus test by taking maximum
    F0 <- apply(F0, 2, stats.combine.func)
    Fp <- apply(Fp, c(2, 3), stats.combine.func)  # otu.no by perm.no
    
    if (verbose) cat('\n')
    
    if (mean(is.na(F0)) >= 0.1) {
      warning('More than 10% observed F stats have NA! Please check! \n')
    }
    
    if (mean(is.na(Fp)) >= 0.1) {
      warning('More than 10% permuted F stats have NA! Please check! \n')
    }
    
    na.ind <- is.na(F0)
    F0 <- F0[!na.ind]
    Fp <- Fp[!na.ind, ]
    
    which.nan.ind <- which(!na.ind)
    ###############################################################################
    p.raw <- rowMeans(cbind(Fp, F0) >= F0)
    
    if (i == stage.no) {
      break
    } else {
      if (mean(p.raw <= 0.05) >= excl.pct) {
        ind <- order(p.raw)[1 : round(length(p.raw) * excl.pct)]
      } else {
        ind <- which(p.raw <= 0.05)
      }
      size.factor <- colSums(feature.dat[setdiff(ref.ind, which.nan.ind[ind]), ])
      # norm.ind <- cbind(norm.ind, !(1:nrow(feature.dat) %in% setdiff(ref.ind, which.nan.ind[ind])))
    }
  }
  
  
  
  p.adj.fdr <- perm.fdr.adj(F0, Fp)
  p.raw <- na.pad(p.raw, na.ind)
  p.adj.fdr <- na.pad(p.adj.fdr, na.ind)
  names(p.raw) <- names(p.adj.fdr) <- rownames(R2.m) <- rownames(RSS.m) <- rownames(F0.m) <- row.names
  colnames(R2.m) <- colnames(F0.m) <- colnames(RSS.m) <- paste0('Func', 1 : func.no)
  
  if (is.fwer) {
    p.adj.fwer <- perm.fwer.adj(F0, Fp)
    p.adj.fwer <- na.pad(p.adj.fwer, na.ind)
    names(p.adj.fwer)  <- row.names
  } else {
    p.adj.fwer <- NULL
  }
  
  # Finally, get the direction of change under specific transformations
  if (length(ref.ind) != 0) {
    divisor <- colSums(feature.dat[ref.ind, ])
    divisor[divisor == 0] <- 1
    temp.dat <- t(t(feature.dat) / divisor)
  } else {
    temp.dat <- feature.dat
  }
  
  coef.list <- NULL
  for (j in 1 : func.no) {
    func <- link.func[[j]]
    coef.tmp <- tryCatch({
     solve(t(M01) %*% M01) %*% t(M01) %*% t(func(temp.dat))
    }, error = function(e) {
      ginv(t(M01) %*% M01) %*% t(M01) %*% t(func(temp.dat))
      
    })
    
    if(is.null(rownames(coef.tmp))){
      rownames(coef.tmp) <- colnames(M01)
    }
    coef.list[[j]] <- coef.tmp
    
  }
  
  
  if (!return.feature.dat) {
    feature.dat <- NULL
  }
  if (verbose) cat('Completed!\n')
  
  return(list(call = this.call, feature.dat = feature.dat, meta.dat = meta.dat, grp.name = grp.name,
              filter.features = filter.features, ref.features = row.names[ref.ind],
              R2 = R2.m, F0 = F0.m, RSS = RSS.m, df.model = df.model, df.residual = df.residual, coef.list = coef.list,
              p.raw = p.raw, p.adj.fdr = p.adj.fdr,  p.adj.fwer = p.adj.fwer))
  
}




#' Plot ZicoSeq results
#' The function plots the association strength and direction based on the output from \code{ZicoSeq}.
#' @param ZicoSeq.obj return from function \code{ZicoSeq}.
#' @param pvalue.type character; It could be 'p.adj.fdr','p.raw' or 'p.adj.fwer'.
#' @param cutoff a real value between 0 and 1; cutoff for pvalue.type.
#' @param text.size text size for the plots.
#' @param out.dir character; the directory to save the figures, e.g., \code{getwd()}. Default is NULL. If NULL, figures will not be saved.
#' @param file.name character; the name of the file.
#' @param width the width of the graphics region in inches. See R function \code{ggsave}.
#' @param height the height of the graphics region in inches. See R function \code{ggsave}.
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @import dplyr
#' @import tibble
#' @importFrom ggpubr ggarrange
#' @rdname ZicoSeq.plot
#' @export
#' 


ZicoSeq.plot <- function(ZicoSeq.obj, pvalue.type = c('p.adj.fdr','p.raw','p.adj.fwer'), 
                         cutoff = 0.1, text.size = 10, out.dir = NULL, file.name = 'ZicoSeq.plot.pdf',  width = 10, height = 6){
  
  #  if(sum(ZicoSeq.obj[[pvalue.type]]<= cutoff) == 0){
  #    cat(paste0('No taxa within ', pvalue.type, ' <= ',cutoff, '!'))
  #  }
  
  grp.name <-  ZicoSeq.obj$grp.name
  meta.dat <- ZicoSeq.obj$meta.dat
  
  # cat('Significant(',pvalue.type,'<=',cutoff,') association strength between taxa and ',grp.name, ' is visualized!\n' )
  
  if(ZicoSeq.obj$call$feature.dat.type == 'proportion'){
    abundance <- ZicoSeq.obj$feature.dat
    prevalence <- apply(ZicoSeq.obj$feature.dat, 1, function(x) mean(x > 0))
  }
  
  if(ZicoSeq.obj$call$feature.dat.type == 'count'){
    abundance <- t(t(ZicoSeq.obj$feature.dat)/colSums(ZicoSeq.obj$feature.dat))
    prevalence <- apply(ZicoSeq.obj$feature.dat, 1, function(x) mean(x > 0))
  }
  
  if(ZicoSeq.obj$call$feature.dat.type == 'other'){
    abundance <- ZicoSeq.obj$feature.dat
    prevalence <- apply(abundance, 1, function(x) sd(x))
  }
  
  # relative abundance
  abundance <- apply(abundance, 1, function(x) mean(x))
  
  if(!(length(unique(meta.dat[,grp.name])) > 2 & !is.numeric(meta.dat[,grp.name]))) {
    coefs <- sapply(ZicoSeq.obj$coef.list, function(x) x[grep(grp.name,rownames(x)),])
    colnames(coefs) <- paste0('coef_Func',1:length(ZicoSeq.obj$coef.list))
  }  
  
  if (!is.numeric(meta.dat[,grp.name]) & length(unique(meta.dat[,grp.name])) == 2) {
    level2 <- rownames(ZicoSeq.obj$coef.list[[1]])[grep(paste0('^',grp.name),rownames(ZicoSeq.obj$coef.list[[1]]))]
    level2 <- gsub(grp.name, '', level2)
    base <- setdiff(unique(meta.dat[,grp.name]), level2)
  }
  
  # Break the ties
  ZicoSeq.obj$R2 <- ZicoSeq.obj$R2 + runif(length(ZicoSeq.obj$R2), 0, 1e-10)
  
  R2 <- matrixStats::rowMaxs(ZicoSeq.obj$R2)
  
  ## pvalue = 0??? -log10(p)
  if(length(unique(meta.dat[,grp.name])) > 2 & !is.numeric(meta.dat[,grp.name])){
    plot.data <- data.frame(pvals = ZicoSeq.obj[[pvalue.type]], 
                            prevalence = prevalence, abundance = abundance, 
                            R2 = R2, 
                            taxa = rownames(ZicoSeq.obj$R2))
    plot.data[plot.data$pvals > cutoff, 'taxa'] <- ''
    
    
  }else{
    
    signs <- t(sign(coefs))[t(ZicoSeq.obj$R2 == R2)]
    signs[signs == 0] <- 1
    
    plot.data <- data.frame(pvals = ZicoSeq.obj[[pvalue.type]], 
                            prevalence = prevalence, abundance = abundance, 
                            R2 = R2 * signs, 
                            taxa = rownames(ZicoSeq.obj$R2))
    
    plot.data[plot.data$pvals > cutoff, 'taxa'] <- ''
    
    
  }
  
  if(is.numeric(meta.dat[,grp.name])) {
    title.name = paste0('Association between ', grp.name, ' and  taxa abundance')
  }
  
  if(!is.numeric(meta.dat[,grp.name]) & length(unique(meta.dat[,grp.name])) > 2){
    title.name = paste0('Association between ', grp.name, ' and  taxa abundance')
  }
  
  if(!is.numeric(meta.dat[,grp.name]) & length(unique(meta.dat[,grp.name])) == 2){
    title.name = paste0('Differential abundance between ', paste0(base, ' (reference) and ', level2))
  }
  
  pvals <- R2 <- taxa <- NULL
  plot.obj <-
    ggplot(plot.data, aes(x = R2, y = -log10(pvals))) +
    geom_point(aes(size = abundance, color = prevalence)) +
    geom_vline(aes(xintercept = 0), color = 'gray', linetype = 'dashed') +
    geom_hline(aes(yintercept = -log10(cutoff)), color = 'gray', linetype = 'dashed') +
    scale_colour_gradient2(low = "white", high = "#006D2C") +
    scale_y_continuous(limits = c(0, max(-log10(plot.data$pvals)) * 1.3)) +
    ggrepel::geom_text_repel(aes(label = taxa), max.overlaps = Inf, color = 'black') +
    labs(x = bquote(R^2), y = paste0('-log10(',pvalue.type,')'),
         color = ifelse(ZicoSeq.obj$call$feature.dat.type == 'other','Standard deviation','Prevalence'),
         size = ifelse(ZicoSeq.obj$call$feature.dat.type == 'other','Mean value','Mean abundance')) +
    theme_bw() +
    theme(axis.text = element_text(color = 'black', size = text.size),
          axis.title = element_text(color = 'black', size = text.size),
          legend.text = element_text(color = 'black', size = text.size),
          legend.title = element_text(color = 'black', size = text.size)) +
    ggtitle(title.name)
  
  
  
  #	plots <- ggarrange(plotlist = plot.list, common.legend = T)
  
  if(!is.null(out.dir)) {
    print(plot.obj)
    ggsave(paste0(out.dir, file.name), width = width, height = height)
  }
  return(plot.obj)
}

















vkey2 <- function (map, title = NA, side = 2, stretch = 1.4, x, y, 
		wh)
{
	opar <- par(xpd = NA)
	on.exit(par(opar))
	n <- length(map$breaks) + 1
	dy <- strheight("A")
	aspect <- diff(grconvertX(1:2, from = "inches"))/diff(grconvertY(1:2,
					from = "inches"))
	dx <- dy * aspect
	if (missing(wh)) {
		
		wh <- 1:(n-1)
		
	}
	labs <- format(map$breaks[wh])
	maxlabwidth <- max(strwidth(labs))
	if (missing(x)) {
		x <- grconvertX(1, from = "nfc") - (2 * dx)
		if (side == 4)
			x <- x - maxlabwidth - dx
	}
	else {
		if (is.list(x)) {
			y <- x$y
			x <- x$x
		}
	}
	if (missing(y))
		y <- par("usr")[3] + dy
	ybord <- y + ((0:(n - 1)) * dy * stretch)
	rect(x, ybord[-n], x + dx, ybord[-1], col = map$colors, border = NA)
	if (side == 4) {
		xtext <- x + dx
		text(x = x, y = ybord[n] + (1.5 * dy), title, adj = c(0,
						0))
	}
	if (side == 2) {
		xtext <- x
		text(x = x + dx, y = ybord[n] + (1.5 * dy), title, adj = c(1,
						0))
	}
	text(x = xtext, y = ybord[wh] + 0.5 * dy, labels = labs, pos = side)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	require(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
				ncol = cols, nrow = ceiling(numPlots/cols), byrow=TRUE)
	}
	
	if (numPlots==1) {
		print(plots[[1]])
		
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
							layout.pos.col = matchidx$col))
		}
	}
}


heatmap.3 <- function(x,
		Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
		distfun = dist,
		hclustfun = hclust,
		dendrogram = c("both","row", "column", "none"),
		symm = FALSE,
		scale = c("none","row", "column"),
		na.rm = TRUE,
		revC = identical(Colv,"Rowv"),
		add.expr,
		breaks,
		symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
		col = "heat.colors",
		colsep,
		rowsep,
		sepcolor = "white",
		sepwidth = c(0.05, 0.05),
		cellnote,
		notecex = 1,
		notecol = "cyan",
		na.color = par("bg"),
		trace = c("none", "column","row", "both"),
		tracecol = "cyan",
		hline = median(breaks),
		vline = median(breaks),
		linecol = tracecol,
		margins = c(5,5),
		ColSideColors = NULL,
		RowSideColors = NULL,
		side.height.fraction=0.3,
		cexRow = 0.2 + 1/log10(nr),
		cexCol = 0.2 + 1/log10(nc),
		labRow = NULL,
		labCol = NULL,
		key = TRUE,
		keysize = 1.5,
		density.info = c("none", "histogram", "density"),
		denscol = tracecol,
		symkey = max(x < 0, na.rm = TRUE) || symbreaks,
		densadj = 0.25,
		main = NULL,
		xlab = NULL,
		ylab = NULL,
		lmat = NULL,
		lhei = NULL,
		lwid = NULL,
		NumColSideColors = 1,
		NumRowSideColors = 1,
		KeyValueName="Value",...){
	
	invalid <- function (x) {
		if (missing(x) || is.null(x) || length(x) == 0)
			return(TRUE)
		if (is.list(x))
			return(all(sapply(x, invalid)))
		else if (is.vector(x))
			return(all(is.na(x)))
		else return(FALSE)
	}
	
	x <- as.matrix(x)
	scale01 <- function(x, low = min(x), high = max(x)) {
		x <- (x - low)/(high - low)
		x
	}
	retval <- list()
	scale <- if (symm && missing(scale))
				"none"
			else match.arg(scale)
	dendrogram <- match.arg(dendrogram)
	trace <- match.arg(trace)
	density.info <- match.arg(density.info)
	if (length(col) == 1 && is.character(col))
		col <- get(col, mode = "function")
	if (!missing(breaks) && (scale != "none"))
		warning("Using scale=\"row\" or scale=\"column\" when breaks are",
				"specified can produce unpredictable results.", "Please consider using only one or the other.")
	if (is.null(Rowv) || is.na(Rowv))
		Rowv <- FALSE
	if (is.null(Colv) || is.na(Colv))
		Colv <- FALSE
	else if (Colv == "Rowv" && !isTRUE(Rowv))
		Colv <- FALSE
	if (length(di <- dim(x)) != 2 || !is.numeric(x))
		stop("`x' must be a numeric matrix")
	nr <- di[1]
	nc <- di[2]
	if (nr <= 1 || nc <= 1)
		stop("`x' must have at least 2 rows and 2 columns")
	if (!is.numeric(margins) || length(margins) != 2)
		stop("`margins' must be a numeric vector of length 2")
	if (missing(cellnote))
		cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
	if (!inherits(Rowv, "dendrogram")) {
		if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
					c("both", "row"))) {
			if (is.logical(Colv) && (Colv))
				dendrogram <- "column"
			else dedrogram <- "none"
			warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
					dendrogram, "'. Omitting row dendogram.")
		}
	}
	if (!inherits(Colv, "dendrogram")) {
		if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
					c("both", "column"))) {
			if (is.logical(Rowv) && (Rowv))
				dendrogram <- "row"
			else dendrogram <- "none"
			warning("Discrepancy: Colv is FALSE, while dendrogram is `",
					dendrogram, "'. Omitting column dendogram.")
		}
	}
	if (inherits(Rowv, "dendrogram")) {
		ddr <- Rowv
		rowInd <- order.dendrogram(ddr)
	}
	else if (is.integer(Rowv)) {
		hcr <- hclustfun(distfun(x))
		ddr <- as.dendrogram(hcr)
		ddr <- reorder(ddr, Rowv)
		rowInd <- order.dendrogram(ddr)
		if (nr != length(rowInd))
			stop("row dendrogram ordering gave index of wrong length")
	}
	else if (isTRUE(Rowv)) {
		Rowv <- rowMeans(x, na.rm = na.rm)
		hcr <- hclustfun(distfun(x))
		ddr <- as.dendrogram(hcr)
		ddr <- reorder(ddr, Rowv)
		rowInd <- order.dendrogram(ddr)
		if (nr != length(rowInd))
			stop("row dendrogram ordering gave index of wrong length")
	}
	else {
		rowInd <- nr:1
	}
	if (inherits(Colv, "dendrogram")) {
		ddc <- Colv
		colInd <- order.dendrogram(ddc)
	}
	else if (identical(Colv, "Rowv")) {
		if (nr != nc)
			stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
		if (exists("ddr")) {
			ddc <- ddr
			colInd <- order.dendrogram(ddc)
		}
		else colInd <- rowInd
	}
	else if (is.integer(Colv)) {
		hcc <- hclustfun(distfun(if (symm)
									x
								else t(x)))
		ddc <- as.dendrogram(hcc)
		ddc <- reorder(ddc, Colv)
		colInd <- order.dendrogram(ddc)
		if (nc != length(colInd))
			stop("column dendrogram ordering gave index of wrong length")
	}
	else if (isTRUE(Colv)) {
		Colv <- colMeans(x, na.rm = na.rm)
		hcc <- hclustfun(distfun(if (symm)
									x
								else t(x)))
		ddc <- as.dendrogram(hcc)
		ddc <- reorder(ddc, Colv)
		colInd <- order.dendrogram(ddc)
		if (nc != length(colInd))
			stop("column dendrogram ordering gave index of wrong length")
	}
	else {
		colInd <- 1:nc
	}
	retval$rowInd <- rowInd
	retval$colInd <- colInd
	retval$call <- match.call()
	x <- x[rowInd, colInd]
	x.unscaled <- x
	cellnote <- cellnote[rowInd, colInd]
	if (is.null(labRow))
		labRow <- if (is.null(rownames(x)))
					(1:nr)[rowInd]
				else rownames(x)
	else labRow <- labRow[rowInd]
	if (is.null(labCol))
		labCol <- if (is.null(colnames(x)))
					(1:nc)[colInd]
				else colnames(x)
	else labCol <- labCol[colInd]
	if (scale == "row") {
		retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
		x <- sweep(x, 1, rm)
		retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
		x <- sweep(x, 1, sx, "/")
	}
	else if (scale == "column") {
		retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
		x <- sweep(x, 2, rm)
		retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
		x <- sweep(x, 2, sx, "/")
	}
	if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
		if (missing(col) || is.function(col))
			breaks <- 16
		else breaks <- length(col) + 1
	}
	if (length(breaks) == 1) {
		if (!symbreaks)
			breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
					length = breaks)
		else {
			extreme <- max(abs(x), na.rm = TRUE)
			breaks <- seq(-extreme, extreme, length = breaks)
		}
	}
	nbr <- length(breaks)
	ncol <- length(breaks) - 1
	if (class(col) == "function")
		col <- col(ncol)
	min.breaks <- min(breaks)
	max.breaks <- max(breaks)
	x[x < min.breaks] <- min.breaks
	x[x > max.breaks] <- max.breaks
	if (missing(lhei) || is.null(lhei))
		lhei <- c(keysize, 4)
	if (missing(lwid) || is.null(lwid))
		lwid <- c(keysize, 4)
	if (missing(lmat) || is.null(lmat)) {
		lmat <- rbind(4:3, 2:1)
		
		if (!is.null(ColSideColors)) {
			#if (!is.matrix(ColSideColors))
			#stop("'ColSideColors' must be a matrix")
			if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
				stop("'ColSideColors' must be a matrix of nrow(x) rows")
			lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
			#lhei <- c(lhei[1], 0.2, lhei[2])
			lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
		}
		
		if (!is.null(RowSideColors)) {
			#if (!is.matrix(RowSideColors))
			#stop("'RowSideColors' must be a matrix")
			if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
				stop("'RowSideColors' must be a matrix of ncol(x) columns")
			lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
			#lwid <- c(lwid[1], 0.2, lwid[2])
			lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
		}
		lmat[is.na(lmat)] <- 0
	}
	
	if (length(lhei) != nrow(lmat))
		stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
	if (length(lwid) != ncol(lmat))
		stop("lwid must have length = ncol(lmat) =", ncol(lmat))
	op <- par(no.readonly = TRUE)
	on.exit(par(op))
	
	layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
	
	if (!is.null(RowSideColors)) {
		if (!is.matrix(RowSideColors)){
			par(mar = c(margins[1], 0, 0, 0.5))
			image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
			box()
		} else {
			par(mar = c(margins[1], 0, 0, 0.5))
			rsc = t(RowSideColors[,rowInd, drop=F])
			rsc.colors = matrix()
			rsc.names = names(table(rsc))
			rsc.i = 1
			for (rsc.name in rsc.names) {
				rsc.colors[rsc.i] = rsc.name
				rsc[rsc == rsc.name] = rsc.i
				rsc.i = rsc.i + 1
			}
			rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
			image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
			
			if (length(rownames(RowSideColors)) > 0) {
				axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), rownames(RowSideColors), las = 2, tick = FALSE)
			}
			box()
		}
	}
	
	if (!is.null(ColSideColors)) {
		
		if (!is.matrix(ColSideColors)){
			par(mar = c(0.5, 0, 0, margins[2]))
			image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
			box()
		} else {
			par(mar = c(0.5, 0, 0, margins[2]))
			csc = ColSideColors[colInd, , drop=F]
			csc.colors = matrix()
			csc.names = names(table(csc))
			csc.i = 1
			for (csc.name in csc.names) {
				csc.colors[csc.i] = csc.name
				csc[csc == csc.name] = csc.i
				csc.i = csc.i + 1
			}
			csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
			image(csc, col = as.vector(csc.colors), axes = FALSE)
			if (length(colnames(ColSideColors)) > 0) {
				axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
			}
			box()
		}
	}
	
	par(mar = c(margins[1], 0, 0, margins[2]))
	x <- t(x)
	cellnote <- t(cellnote)
	if (revC) {
		iy <- nr:1
		if (exists("ddr"))
			ddr <- rev(ddr)
		x <- x[, iy]
		cellnote <- cellnote[, iy]
	}
	else iy <- 1:nr
	image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
	retval$carpet <- x
	if (exists("ddr"))
		retval$rowDendrogram <- ddr
	if (exists("ddc"))
		retval$colDendrogram <- ddc
	retval$breaks <- breaks
	retval$col <- col
	if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
		mmat <- ifelse(is.na(x), 1, NA)
		image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
				col = na.color, add = TRUE)
	}
	axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
			cex.axis = cexCol)
	if (!is.null(xlab))
		mtext(xlab, side = 1, line = margins[1] - 1.25)
	axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
			cex.axis = cexRow)
	if (!is.null(ylab))
		mtext(ylab, side = 4, line = margins[2] - 1.25)
	if (!missing(add.expr))
		eval(substitute(add.expr))
	if (!missing(colsep))
		for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
	if (!missing(rowsep))
		for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
	min.scale <- min(breaks)
	max.scale <- max(breaks)
	x.scaled <- scale01(t(x), min.scale, max.scale)
	if (trace %in% c("both", "column")) {
		retval$vline <- vline
		vline.vals <- scale01(vline, min.scale, max.scale)
		for (i in colInd) {
			if (!is.null(vline)) {
				abline(v = i - 0.5 + vline.vals, col = linecol,
						lty = 2)
			}
			xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
			xv <- c(xv[1], xv)
			yv <- 1:length(xv) - 0.5
			lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
		}
	}
	if (trace %in% c("both", "row")) {
		retval$hline <- hline
		hline.vals <- scale01(hline, min.scale, max.scale)
		for (i in rowInd) {
			if (!is.null(hline)) {
				abline(h = i + hline, col = linecol, lty = 2)
			}
			yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
			yv <- rev(c(yv[1], yv))
			xv <- length(yv):1 - 0.5
			lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
		}
	}
	if (!missing(cellnote))
		text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
				col = notecol, cex = notecex)
	par(mar = c(margins[1], 0, 0, 0))
	if (dendrogram %in% c("both", "row")) {
		plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
	}
	else plot.new()
	par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
	if (dendrogram %in% c("both", "column")) {
		plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
	}
	else plot.new()
	if (!is.null(main))
		title(main, cex.main = 1.5 * op[["cex.main"]])
	if (key) {
		par(mar = c(5, 4, 2, 1), cex = 0.75)
		tmpbreaks <- breaks
		if (symkey) {
			max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
			min.raw <- -max.raw
			tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
			tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
		}
		else {
			min.raw <- min(x, na.rm = TRUE)
			max.raw <- max(x, na.rm = TRUE)
		}
		
		z <- seq(min.raw, max.raw, length = length(col))
		image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
				xaxt = "n", yaxt = "n")
		par(usr = c(0, 1, 0, 1))
		lv <- pretty(breaks)
		xv <- scale01(as.numeric(lv), min.raw, max.raw)
		axis(1, at = xv, labels = lv)
		if (scale == "row")
			mtext(side = 1, "Row Z-Score", line = 2)
		else if (scale == "column")
			mtext(side = 1, "Column Z-Score", line = 2)
		else mtext(side = 1, KeyValueName, line = 2)
		if (density.info == "density") {
			dens <- density(x, adjust = densadj, na.rm = TRUE)
			omit <- dens$x < min(breaks) | dens$x > max(breaks)
			dens$x <- dens$x[-omit]
			dens$y <- dens$y[-omit]
			dens$x <- scale01(dens$x, min.raw, max.raw)
			lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
					lwd = 1)
			axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
			title("Color Key\nand Density Plot")
			par(cex = 0.5)
			mtext(side = 2, "Density", line = 2)
		}
		else if (density.info == "histogram") {
			h <- hist(x, plot = FALSE, breaks = breaks)
			hx <- scale01(breaks, min.raw, max.raw)
			hy <- c(h$counts, h$counts[length(h$counts)])
			lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
					col = denscol)
			axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
			title("Color Key\nand Histogram")
			par(cex = 0.5)
			mtext(side = 2, "Count", line = 2)
		}
		else title("Color Key")
	}
	else plot.new()
	retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
			high = retval$breaks[-1], color = retval$col)
	invisible(retval)
}


load_package <- function (package.names=c('rmutil','MASS',  'ape', 'vegan', 'pscl', 'glmmADMB', 'aod', 'nlme', 'MiRKAT', 'matrixStats', 
				'gplots', 'scales', 'ggplot2', 'rpart', 'qvalue', 'DESeq2', 'lme4', 'dplyr', 'tidyr','tibble',
				 'phyloseq','GUniFrac', 'RColorBrewer', 'squash', 'rhdf5', 'biom', 'reshape', 'randomForest',  'ade4',  'arsenal'), lib.loc=NULL) { #'phangorn', 
	for (package.name in package.names) {
		require(package.name, character.only=T, lib.loc = lib.loc)
	}
}

getPermuteMatrixBlock <- function(permutations, strata) {
	strata <- factor(strata)
	strata <- factor(as.numeric(strata))
	res <- sapply(1:permutations, function(i) {
				strata1 <- strata
				levels(strata1) <- sample(levels(strata1))
				order(as.numeric(as.character(strata1)))
			})
	t(res)
}

# adonis2: permutate the covariate instead of data matrix - validated
# The variable of interest should be put at the last; Only p values of variable of interest will be reported
# Currently it can only be applied to subject-specific variable (nonvarying within subject) but covariates are allowed to vary
# Need to improve the speed
adonis2 <- function (formula, data = NULL, permutations = 999, method = "bray", 
		strata = NULL, block.perm = TRUE, contr.unordered = "contr.sum", contr.ordered = "contr.poly", 
		...) 
{
	TOL <- 1e-07
	Terms <- terms(formula, data = data, keep.order=TRUE)
	lhs <- formula[[2]]
	lhs <- eval(lhs, data, parent.frame())
	formula[[2]] <- NULL
	rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
	op.c <- options()$contrasts
	options(contrasts = c(contr.unordered, contr.ordered))
	rhs <- model.matrix(formula, rhs.frame)
	options(contrasts = op.c)
	grps <- attr(rhs, "assign")
	qrhs <- qr(rhs)
	rhs <- rhs[, qrhs$pivot, drop = FALSE]
	rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
	grps <- grps[qrhs$pivot][1:qrhs$rank]
	u.grps <- unique(grps)
	nterms <- length(u.grps) - 1
	st <- ifelse(nterms == 1, 2, nterms)
	H.s <- lapply(st:(nterms+1), function(j) {
				Xj <- rhs[, grps %in% u.grps[1:j]]
				qrX <- qr(Xj, tol = TOL)
				Q <- qr.Q(qrX)
				tcrossprod(Q[, 1:qrX$rank])
			})
	if (inherits(lhs, "dist")) {
		if (any(lhs < -TOL)) 
			stop("dissimilarities must be non-negative")
		dmat <- as.matrix(lhs^2)
	}
	else {
		dist.lhs <- as.matrix(vegdist(lhs, method = method, ...))
		dmat <- dist.lhs^2
	}
	n <- nrow(dmat)
	G <- -sweep(dmat, 1, rowMeans(dmat))/2
	SS.Exp.comb <- sapply(H.s, function(hat) sum(G * t(hat)))
	if (nterms == 1) {
		SS.Exp.each <- SS.Exp.comb
	} else {
		SS.Exp.each <- SS.Exp.comb[2] - SS.Exp.comb[1]
	}
	
	H.snterm <- H.s[[length(H.s)]]
	tIH.snterm <- t(diag(n) - H.snterm)
	
#	if (length(H.s) > 1) {
#		H.s[[2]] <- H.s[[2]] - H.s[[1]]
#	}
	
	SS.Res <- sum(G * tIH.snterm)
	df.Exp <- sum(grps == u.grps[length(u.grps)])
	df.Res <- n - qrhs$rank
	
	if (inherits(lhs, "dist")) {
		beta.sites <- qr.coef(qrhs, as.matrix(lhs))
		beta.spp <- NULL
	}
	else {
		beta.sites <- qr.coef(qrhs, dist.lhs)
		beta.spp <- qr.coef(qrhs, as.matrix(lhs))
	}
	colnames(beta.spp) <- colnames(lhs)
	colnames(beta.sites) <- rownames(lhs)
	
	F.Mod <- (SS.Exp.each/df.Exp)/(SS.Res/df.Res)
	
	f.test <- function(tH, G, df.Exp, df.Res, tIH.snterm) {
		(sum(G * tH)/df.Exp)/(sum(G * tIH.snterm)/df.Res)
	}
	
	rhs.1 <- rhs[, grps %in% u.grps[1:nterms], drop=F]
	rhs.2 <- rhs[, grps %in% u.grps[(nterms+1)], drop=F]	
	
	if (missing(strata)) strata <- NULL
	
	if (block.perm == FALSE) {
		permute.ind <- vegan:::getPermuteMatrix(permutations, n, strata = strata)
	} else {
		if (is.null(strata)) stop('Block permutation requires strata!\n')
		strata.u <- unique(strata)
		reorder.ind <- unlist(lapply(strata.u, function(x) which(strata == x)))
		expand.ind <- rep(1:length(strata.u), sapply(strata.u, function(x) sum(strata == x)))
		rhs.2.u <- rhs.2[sapply(strata.u, function(x) which(strata == x)[1]), , drop=F]
	}
	
	# To be checked and validated. 
	f.perms <- 
			sapply(1:permutations, function(i) {
						if (block.perm == FALSE) {
							rhs.2 <- rhs[permute.ind[i, ], grps %in% u.grps[(nterms+1)], drop=F]	
						} else {
							rhs.2[reorder.ind, ] <- rhs.2.u[sample(nrow(rhs.2.u)), , drop=F][expand.ind, , drop=F]
						}		
						Xj <- cbind(rhs.1, rhs.2)
						qrX <- qr(Xj, tol = TOL)
						Q <- qr.Q(qrX)
						tH.snterm <- t(tcrossprod(Q[, 1:qrX$rank]))
						tIH.snterm <- diag(n) - tH.snterm
						
						if (nterms > 1) {
							tH.snterm <- tH.snterm - t(H.s[[1]])
						}							
						f.test(tH.snterm, G, df.Exp, df.Res, 
								tIH.snterm)
					}
			)
	
	f.perms <- round(f.perms, 12)
	F.Mod <- round(F.Mod, 12)
	SumsOfSqs = c(SS.Exp.each, SS.Res, SS.Exp.comb[length(SS.Exp.comb)] + SS.Res)
	tab <- data.frame(Df = c(df.Exp, df.Res, n - 1), SumsOfSqs = SumsOfSqs, 
			MeanSqs = c(SS.Exp.each/df.Exp, SS.Res/df.Res, NA), F.Model = c(F.Mod, 
					NA, NA), R2 = SumsOfSqs/SumsOfSqs[length(SumsOfSqs)], 
			P = c((sum(f.perms >= F.Mod) + 1)/(permutations + 
								1), NA, NA))
	rownames(tab) <- c(attr(attr(rhs.frame, "terms"), "term.labels")[u.grps[length(u.grps)]], 
			"Residuals", "Total")
	colnames(tab)[ncol(tab)] <- "Pr(>F)"
	attr(tab, "heading") <- "Terms added sequentially (first to last)\n"
	class(tab) <- c("anova", class(tab))
	out <- list(aov.tab = tab, call = match.call(), coefficients = beta.spp, 
			coef.sites = beta.sites, f.perms = as.matrix(f.perms), model.matrix = rhs, 
			terms = Terms)
	class(out) <- "adonis"
	out
}

######################################
# New: 2018_03_09
set_colour_theme <- function (type = 'stata', pal = palette(RColorBrewer::brewer.pal(5, "Set1"))) {
	library(ggthemes) # Load
	if (type == 'stata') {
		scale_colour_discrete <<- function(...) scale_color_stata()
		scale_fill_discrete <<- function(...) scale_fill_stata()
	} 
	if (type == 'economist') {
		scale_colour_discrete <<- function(...) scale_color_economist()
		scale_fill_discrete <<- function(...) scale_fill_economist()
	}
	if (type == 'wsj') {
		scale_colour_discrete <<- function(...) scale_color_wsj()
		scale_fill_discrete <- function(...) scale_fill_wsj()
	}
	if (type == 'hc') {
		scale_colour_discrete <<- function(...) scale_color_hc()
		scale_fill_discrete <<- function(...) scale_fill_hc()
	}
	if (type == 'default') {
		scale_fill_discrete <<-  ggplot2::scale_fill_discrete
		scale_colour_discrete <<- ggplot2::scale_colour_discrete
	}
	if (type ==  'cb') {
		cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	
		scale_fill_discrete <<- function(...) scale_fill_manual(values=cbPalette)
		scale_colour_discrete <<- function(...) scale_colour_manual(values=cbPalette)
	}
	if (type ==  'cbb') {

		cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
		scale_fill_discrete <<- function(...) scale_fill_manual(values=cbbPalette)
		scale_colour_discrete <<- function(...) scale_colour_manual(values=cbbPalette)
	}
	
	if (type == 'customize') {
		scale_fill_discrete <<- function(...) scale_fill_manual(values=pal)
		scale_colour_discrete <<- function(...) scale_colour_manual(values=pal)
	}

}

########################################
# This generates the matrix columns-wise
# From JnPaulson
generate_matrix <- function(x){
	indptr  = x$sample$matrix$indptr+1
	indices = x$sample$matrix$indices+1
	data    = x$sample$matrix$data
	nr = length(x$observation$ids)
	
	counts = sapply(2:length(indptr),function(i){
				x = rep(0,nr)
				seq = indptr[i-1]:(indptr[i]-1)
				x[indices[seq]] = data[seq]
				x
			})
	rownames(counts) = x$observation$ids
	colnames(counts) = x$sample$ids
	# I wish this next line wasn't necessary
	lapply(1:nrow(counts),function(i){
				counts[i,]
			})
}
generate_metadata <- function(x){
	metadata = x$metadata
	metadata = lapply(1:length(x$ids),function(i){
				id_metadata = lapply(metadata,function(j){
							if(length(dim(j))>1){ as.vector(j[,i,drop=FALSE]) }
							else{ j[i] }
						})
				list(id = x$ids[i],metadata=id_metadata)
			})
	return(metadata)
}
namedList <- function(...) {
	L <- list(...)
	snm <- sapply(substitute(list(...)),deparse)[-1]
	if (is.null(nm <- names(L))) nm <- snm
	if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
	setNames(L,nm)
}

read_hdf5_biom <- function(file_input){
	x = h5read(file_input,"/",read.attributes = TRUE)
	data = generate_matrix(x)
	rows = generate_metadata(x$observation)
	columns = generate_metadata(x$sample)
	shape = c(length(data),length(data[[1]])) # dim(data)
	
	# Experimental -- need to actually load these from file
	id = attr(x,"id")
	vs = attr(x,"format-version")
	format = sprintf("Biological Observation Matrix %s.%s",vs[1],vs[2])
	format_url = attr(x,"format-url")
	type = "OTU table"
	#type=attr(x,"type")
	generated_by = attr(x,"generated-by")
	date = attr(x,"creation-date")
	matrix_type = "dense"
	matrix_element_type = "int"
	
	namedList(id,format,format_url,type,generated_by,date,matrix_type,matrix_element_type,
			rows,columns,shape,data)
}
############################################
# comm are otu counts, row: otus, column samples
# intersect.no: Pairwise ratio calculated on pairs with at least 'intersect.no' common taxa; 
# Rev: 2016_12_12 add prevalence filter, which is not necessary
# Rev: 2017_02_07 
GMPR.old <- function (comm, intersect.no=4, prev.filter=0.0) {
	comm <- comm[rowMeans(comm != 0) >= prev.filter, ]
	
	ind.vec <- numeric(ncol(comm))
	res <- sapply(1:ncol(comm),  function(i) {
				x <- comm[, i]
				pr <- sapply(1:ncol(comm),  function(j) {
							y <- comm[, j]
							ind <- x != 0 & y != 0
							if (sum(ind) >= intersect.no) {
								res <- median(x[ind] / y[ind])
							} else {
								res <- NA
							}
						})
				if (sum(is.na(pr)) != 0) ind.vec[i] <<- 1
				exp(mean(log(pr[!is.na(pr)])))
			}
	
	)
	if (sum(ind.vec)) warnings(paste0('Sample ', paste(which(ind.vec!=0), collapse=' '), ' do not have at least ', 
						intersect.no, ' common taxa with some of the other samples. '))
	if (sum(is.nan(res))) {
		ind <- is.nan(res)
		res[ind] <- 1
		warning(paste0('Sample ', paste(which(ind), collapse=' '), ' do not have at least ', intersect.no, ' common taxa for any of the other samples. ',
						'For these samples, we force the size factors to be 1! ', 
						'You may consider removing these samples since they are very different from other samples!\n'))
	}
	cat('Finished!')
	res
}

# New: 2017_02_07
GMPR <- function (comm, intersect.no=4, ct.min=2, trace=FALSE) {
	# Computes the GMPR size factor
	#
	# Args:
	#   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
	#   intersect.no: the minimum number of shared features between sample pair, where the ratio is calculated
	#   ct.min: the minimum number of counts required to calculate ratios ct.min = 5 has better results
	
	#
	# Returns:
	#   a list that contains:
	#      gmpr： the GMPR size factors for all samples; Samples with distinct sets of features will be output as NA.
	#      nss:   number of samples with significant sharing (> intersect.no) including itself
	
	# mask counts < ct.min
	comm[comm < ct.min] <- 0
	
	if (is.null(colnames(comm))) {
		colnames(comm) <- paste0('S', 1:ncol(comm))
	}
	
	if(trace) cat('Begin GMPR size factor calculation ...\n')
	
	comm.no <- numeric(ncol(comm))
	gmpr <- sapply(1:ncol(comm),  function(i) {		
				if (i %% 50 == 0) {
					if(trace) cat(i, '\n')
				}
				x <- comm[, i]
				# Compute the pairwise ratio
				pr <- x / comm
				# Handling of the NA, NaN, Inf
				pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
				# Counting the number of non-NA, NaN, Inf
				incl.no <- colSums(!is.na(pr))		
				# Calculate the median of PR
				pr.median <- colMedians(pr, na.rm=TRUE)
				# Record the number of samples used for calculating the GMPR
				comm.no[i] <<- sum(incl.no >= intersect.no)
				# Geometric mean of PR median
				if (comm.no[i] > 1) {
					return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
				} else {
					return(NA)
				}
			}
	)
	
	if (sum(is.na(gmpr))) {
		warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'), 
						'\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
						'For these samples, their size factors are set to be NA! \n', 
						'You may consider removing these samples since they are potentially outliers or negative controls!\n',
						'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
	}
	
	if(trace) cat('Completed!\n')
	if(trace) cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
	attr(gmpr, 'NSS') <- comm.no
	# Rev: 2017_09_07
    gmpr <- gmpr * median(colSums(comm))
	names(gmpr) <- colnames(comm)
	return(gmpr)
}


# Rev: 2016_10_25
uniquefy_taxa_names <- function (data.obj) {
	for (level in names(data.obj$abund.list)) {
		obj <- data.obj$abund.list[[level]]
#		rownames(obj) <- gsub('unclassified', paste0('Unclassified', substr(level, 1, 1)), rownames(obj))
#		rownames(obj) <- gsub('unassigned', paste0('Unclassified', substr(level, 1, 1)), rownames(obj))
		rownames(obj) <- paste0(rownames(obj), substr(level, 1, 1))
		data.obj$abund.list[[level]] <- obj
	}
#	suffix <- c('K', 'P', 'C', 'O', 'F', 'G', 'S')
#	for (i in 1:ncol(data.obj$otu.name)) {
#		data.obj$otu.name[, i] <- paste0(data.obj$otu.name[, i], suffix[i])
#	}
	return(data.obj)
}

# New: 2018_02_24
dereplicate <- function (names) {
	names <- factor(names)
	ord <- order(order(names))
	tab <- table(names)
	suffix <- unlist(sapply(tab, function (x) {
				if (x == 1) {
					return('')
				} else {
					return(paste0('_', 1:x))
				}
			}))

	as.character(paste0(names, suffix[ord]))

}

# Rev: 2016_09_22 Add load.map
# Rev: 2016_12_12 Reogranize rarefy, normalize, winsorize
load_data <- function (otu.file, map.file, tree.file=NULL,  load.map=TRUE, parseFunction=parse_taxonomy_greengenes, version='Old', 
		 species=TRUE, filter.no=1, rep.seq=NULL,
		 rff=FALSE, dep=NULL,
	     norm='TSS', level='OTU', intersect.no=4,
		 winsor=FALSE, winsor.qt=0.97,
		 ko.file=NULL, cog.file=NULL, ko.ann.file=NULL,
		 meta.sep='\t', quote="\"", comment="",
		 read.gg=FALSE,  seed=1234, ...) {
	# ko and cog file are not rarefied	
	# filter.no: filter the OTUs with read support less than filter.no (default is filtering singleton); singleton will not be filtered after rarefaction
	# winsorization and GMPR should be further studied. Current default is false and GMPR is on the genus level
	act.seq <- NULL
	max.nchar <- 6
	set.seed(seed)
	if (load.map == TRUE) {
		cat("Load meta file...\n")
		if (grepl("csv$", map.file)) {
			meta.dat <- read.csv(map.file, header=T, check.names=F, row.names=1, comment=comment, quote=quote, ...)
		} else {
			meta.dat <- read.table(map.file, header=T, check.names=F, row.names=1, comment=comment, sep=meta.sep, quote=quote, ...)
		}
	} else {
		meta.dat <- NULL
	}
	
	# Load Tree
	if (!is.null(tree.file)) {
		cat("Load tree file ...\n")
		if (read.gg == F) {
			tree.12 <- read.tree(tree.file)
		} else {
			tree.12 <- read_tree_greengenes(tree.file)
		}
		
		if (is.rooted(tree.12) == F) {
			tree.12 <- midpoint(tree.12)
		}
		
	} else {
		tree.12 <- NULL
	}
	
	cat("Load OTU file...\n")  # Rewrite load new biom file, rev:2016-06-20
	if (version != 'New') {
		biom.obj <-  import_biom(otu.file, parseFunction = parseFunction)  
		
		otu.tab.12 <- otu_table(biom.obj)@.Data
		otu.ind <- rowSums(otu.tab.12) > filter.no  # change otu.tab.12 != 0, rev:2016-06-20
		otu.tab.12 <- otu.tab.12[otu.ind, ]
		# OTU names
		otu.name.full <- as.matrix(biom.obj@tax_table[, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')])
		otu.name.full <- otu.name.full[otu.ind, ]
		
		otu.name.12 <- otu.name.full[, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')]
		otu.name.12[is.na(otu.name.12) | otu.name.12 %in% c('', 'NA')] <- 'unclassified'
		otu.name.12 <- otu.name.12@.Data
		
		otu.name.12[, ] <- gsub('\\[', '', otu.name.12)
		otu.name.12[, ] <- gsub('\\]', '', otu.name.12)
		
		otu.name.12[, ] <- gsub('_unclassified', '', otu.name.12)
		
		otu.name.full <- otu.name.full@.Data
	} else {
		temp <-  read_hdf5_biom(otu.file)
#		otu.file <- paste0(otu.file, '.old')
#		write_biom(temp, otu.file) 
#		biom.obj <- import_biom(otu.file, parseFunction = parseFunction)  
		otu.tab.12 <- matrix(unlist(temp$data), byrow=T, nrow=temp$shape[1], ncol=temp$shape[2])
	
		otu.ids <- sapply(temp$rows, function(x) x[['id']])
		sam.ids <- sapply(temp$columns, function(x) x[['id']])
		
#		res <- t(sapply(temp$rows, function(i) {
#					i$metadata$taxonomy
#				}))
#		rownames(res) <- otu.ids 
	
		# From phyloseq: to be double checked!! Checked!
		if (all(sapply(sapply(temp$rows, function(i) {
									i$metadata
								}), is.null))) {
			otu.name.full <- NULL
		} else {
			taxlist = lapply(temp$rows, function(i) {
						parseFunction(i$metadata$taxonomy)
					})
			names(taxlist) = sapply(temp$rows, function(i) {
						i$id
					})
			otu.name.full = build_tax_table(taxlist)
		}

		# Rev: 2017_04_18
		otu.name.full <- otu.name.full@.Data
		
		rownames(otu.tab.12) <- otu.ids
		colnames(otu.tab.12) <- sam.ids
		
		otu.ind <- rowSums(otu.tab.12) > filter.no  # change otu.tab.12 != 0, rev:2016-06-20
		otu.tab.12 <- otu.tab.12[otu.ind, ]	
		otu.name.full <- otu.name.full[otu.ind, ]
		otu.name.12 <- otu.name.full[, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')]
		otu.name.12[is.na(otu.name.12) | otu.name.12 %in% c('', 'NA')] <- 'unclassified'
		otu.name.12[, ] <- gsub('\\[', '', otu.name.12)
		otu.name.12[, ] <- gsub('\\]', '', otu.name.12)
	}
	
	if (rff == TRUE) {
		cat('Rarefaction ...\n')
		if (is.null(dep)) {
			otu.tab.12 <- t(Rarefy(t(otu.tab.12))$otu.tab.rff)
		} else {
			otu.tab.12 <- t(Rarefy(t(otu.tab.12), dep)$otu.tab.rff)
		}	
		dep <- colSums(otu.tab.12)[1]
		cat("Depth ", dep, '\n')
		# Remove empty OTUs
		otu.ind <- rowSums(otu.tab.12) > 0  # rev:2016-06-28
		otu.tab.12 <- otu.tab.12[otu.ind, ]	
		otu.name.12 <- otu.name.12[otu.ind, ]
		otu.name.full <- otu.name.full[otu.ind, ]
		act.seq <- paste0(act.seq, 'R')
	} else {
		rff <- FALSE
		dep <- NULL
	}

	# Load mapping file
	if (load.map == TRUE) {
		samIDs <- intersect(rownames(meta.dat), colnames(otu.tab.12))
		if (length(samIDs) == 0) {
			stop('Sample names in the meta file and biom are completely different?\n')
		} 
		if (length(samIDs) < length(colnames(otu.tab.12)) | length(samIDs) < length(rownames(meta.dat))) {
			warning('Sample names in the meta file and biom differ! May be due to rarefaction?\n')
		}
		meta.dat <- meta.dat[samIDs, ]
		otu.tab.12 <- otu.tab.12[, samIDs]
	} else {
		samIDs <- colnames(otu.tab.12)
	}

   # Create abundance list
    cat("Create taxa abundance list ...\n")
	abund.list.12 <- list()
	hierachs <- c('Phylum', 'Class', 'Order', 'Family', 'Genus')
	for (hierach in hierachs) {	
		if (hierach != 'Phylum') {
			single.names <- otu.name.12[, hierach]
		#	single.names[grepl('unclassified', single.names, ignore.case=T)] <- paste0('Unclassified',substr(hierach, 1, 1))
			tax.family <- paste(otu.name.12[, 'Phylum'], single.names, sep=";")
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- paste0('Unclassified_', hierach)
		} else {
			tax.family <- otu.name.12[, 'Phylum']
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- 'Unclassified_Phylum'
		}
		family <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
		rownames(family) <- family[, 1]
		family <- as.matrix(family[, -1])
		abund.list.12[[hierach]] <- family
	}
	
	if (species) {
		abund.list.12[['Species']] <- otu.tab.12
		rownames(abund.list.12[['Species']]) <- paste0("A", substr(rownames(otu.tab.12), 1, max.nchar), ":", otu.name.12[, 'Phylum'], ";", otu.name.12[, 'Genus'])
	}
	
    cat('Normalize (size factor) ...\n')
	if (rff == TRUE) {
		cat('For rarefied data, the size factor for samples can still be different!\n')
	}
		
	if (level == 'OTU') {
		data <- otu.tab.12
	} else {
		if (level %in% names(abund.list.12)) {
			data <- abund.list.12[[level]]
		} else {
			data <- otu.tab.12
			level <-'OTU'
			warning('No or wrong level specified! OTU level will be used!\n')
		}
	}
	
	if (norm == 'GMPR') {
		sf <- GMPR(data, intersect.no=intersect.no)
		warning('GMPR is only suitable for samples from the same body location!\n')
		norm <- 'GMPR'
		names(sf) <- colnames(data)
		act.seq <- paste0(act.seq, 'N')
	} else {
		if (norm == 'TSS') {
			sf <- colSums(data)
			norm <- 'TSS'
			act.seq <- paste0(act.seq, 'N')
		} else {
			warning('Normalization method not specified or unknown! TSS is used!\n')
			sf <- colSums(data)
			norm <- 'TSS'
			act.seq <- paste0(act.seq, 'N')
		}
		
	}

	if (winsor == TRUE) {
		act.seq <- paste0(act.seq, 'W')
		if (rff == TRUE) {
			warning('Winsorization after rarefaction will make the data have different total numbers!\n')
		}
		cat('Winsorize ...\n')
		if (is.null(winsor.qt)) {
			winsor.qt <- 0.97
		}
		# Addressing the outlier (97% percent) or at least one outlier
		abund.list.12 <- sapply(abund.list.12, function(genus) {
					genus.p <- t(t(genus) / sf)
					genus.p <- apply(genus.p, 1, function(x) {
								cutoff <- quantile(x, winsor.qt)
								x[x >= cutoff] <- cutoff
								x
							}
					)
					# column/row switch
					genus.w <- t(round(genus.p * sf))
				})
		# OTU table
		otu.tab.12.p <- t(t(otu.tab.12) / sf)
		otu.tab.12.p <- apply(otu.tab.12.p, 1, function(x) {
					cutoff <- quantile(x, winsor.qt)
					x[x >= cutoff] <- cutoff
					x
				}
		)
		# column/row switch
		otu.tab.12 <- t(round(otu.tab.12.p * sf))
	} else {
		winsor <- FALSE
		winsor.qt <- NULL
	}

	# Rarefaction/Normalizing factors are not calculated for functional data
    # Rev: 2017_01_19 the sample IDs for functional data are not ordered! Potential Danger! augment with NA
	if (!is.null(ko.file)) {
		
		cat("Load kegg file...\n")
		ko <- read_biom(ko.file)
		ko.dat <- as.matrix(biom_data(ko))
		
		if (sum(!(samIDs %in% colnames(ko.dat))) != 0) {
			missingIDs <- setdiff(samIDs, colnames(ko.dat))
			aug.mat <- matrix(NA, nrow(ko.dat), length(missingIDs))
			colnames(aug.mat) <- missingIDs
			rownames(aug.mat) <- rownames(ko.dat)
			ko.dat <- cbind(ko.dat, aug.mat) 
		}	
		
		ko.dat <- ko.dat[, samIDs]
		# Rarefaction?
	  
	    if (is.null(ko.ann.file)) {
			# Old - back compatability
			ko.ann <- observation_metadata(ko)
			ko.ann <- cbind(KEGG_Pathways1=sapply(ko.ann, function(x) x['KEGG_Pathways1']), 
					KEGG_Pathways2=sapply(ko.ann, function(x) x['KEGG_Pathways2']), 
					KEGG_Pathways3=sapply(ko.ann, function(x) x['KEGG_Pathways3']))
			rownames(ko.ann) <- rownames(ko.dat)		
			ko.ann[is.na(ko.ann)] <- 'Unclassified'
			
			hierachs <- c("KEGG_Pathways1", "KEGG_Pathways2", "KEGG_Pathways3")
			for (hierach in hierachs) {	
				tax.family <- ko.ann[, hierach]
				family <- aggregate(ko.dat, by=list(tax.family), FUN=sum)
				rownames(family) <- family[, 1]
				family <- as.matrix(family[, -1])
				abund.list.12[[hierach]] <- family
			}
		} else {
			# New
			load(ko.ann.file)
			#
			kos <- rownames(ko.dat)
			abund.list.12[["KEGG_Pathways3"]] <- NULL
			kos.id <- NULL
			for (ko.item in names(kegg.map)) {
				kos.common <- intersect(kos, kegg.map[[ko.item]])
				if (length(kos.common) != 0) {
					abund.list.12[["KEGG_Pathways3"]] <- rbind(abund.list.12[["KEGG_Pathways3"]], colSums(ko.dat[kos.common, , drop=F]))
					kos.id <- c(kos.id, ko.item)
				}
			}
			rownames(abund.list.12[["KEGG_Pathways3"]]) <- kos.id
			
			abund.list.12[["KEGG_Metabolism"]] <- abund.list.12[["KEGG_Pathways3"]][intersect(kos.id, unlist(kegg.ann[['Metabolism']])), ]
			rownames(abund.list.12[["KEGG_Metabolism"]]) <- paste0('M', rownames(abund.list.12[["KEGG_Metabolism"]])) 
			
			abund.list.12[["KEGG_Defense"]] <- NULL
			kos.id <- NULL
			for (ko.item in names(defense.map)) {
				kos.common <- intersect(kos, defense.map[[ko.item]])
				if (length(kos.common) != 0) {
					abund.list.12[["KEGG_Defense"]] <- rbind(abund.list.12[["KEGG_Defense"]], colSums(ko.dat[kos.common, , drop=F]))
					kos.id <- c(kos.id, ko.item)
				}
			}
			rownames(abund.list.12[["KEGG_Defense"]]) <- kos.id
			
			abund.list.12[["KEGG_Toxin"]] <- NULL
			kos.id <- NULL
			for (ko.item in names(toxin.map)) {
				kos.common <- intersect(kos, toxin.map[[ko.item]])
				if (length(kos.common) != 0) {
					abund.list.12[["KEGG_Toxin"]] <- rbind(abund.list.12[["KEGG_Toxin"]], colSums(ko.dat[kos.common, , drop=F]))
					kos.id <- c(kos.id, ko.item)
				}
			}
			rownames(abund.list.12[["KEGG_Toxin"]]) <- kos.id
		}
	
	}
	
	if (!is.null(cog.file)) {
		cat("Load cog file...\n")
		cog <- read_biom(cog.file)
		cog.dat <- as.matrix(biom_data(cog))
		
		if (sum(!(samIDs %in% colnames(cog.dat ))) != 0) {
			missingIDs <- setdiff(samIDs, colnames(cog.dat ))
			aug.mat <- matrix(NA, nrow(cog.dat), length(missingIDs))
			colnames(aug.mat) <- missingIDs
			rownames(aug.mat) <- rownames(cog.dat)
			cog.dat  <- cbind(cog.dat, aug.mat) 
		}	
		
		cog.dat <- cog.dat[, samIDs]
		
		# rarefaction?
		cog.ann <- observation_metadata(cog)
		hierachs <- c("COG_Category1", "COG_Category2")
		for (hierach in hierachs) {	
			tax.family <- sapply(cog.ann, function(x) x[hierach])
			family <- aggregate(cog.dat, by=list(tax.family), FUN=sum)
			rownames(family) <- family[, 1]
			family <- as.matrix(family[, -1])
			abund.list.12[[hierach]] <- family
		}
	}
		
	# Drop tree tips
    if (!is.null(tree.12)) {
		absent <- tree.12$tip.label[!(tree.12$tip.label %in% rownames(otu.tab.12))]
		if (length(absent) != 0) {
			tree.12 <- drop.tip(tree.12, absent)
			warning("The tree has OTUs not in the OTU table!")
		}
	}

	data.obj <- list(otu.tab=otu.tab.12, abund.list=abund.list.12, meta.dat=meta.dat, tree=tree.12,
			otu.name=otu.name.12, otu.name.full=otu.name.full, 
			size.factor=sf, norm.method=norm, norm.level=level, 
			winsor=winsor, winsor.qt=winsor.qt,
			rff=rff, rff.dep=dep, act.seq=act.seq,
			call=match.call())
}


# 2020_01_16: convert phyloseq object to our object
convert_phyloseq_obj <- function (phylo.obj, species = TRUE) {
	
	otu.tab <- otu_table(phylo.obj)@.Data
	meta.dat <- data.frame(sample_data(phylo.obj))
	otu.name <- tax_table(phylo.obj)@.Data
	tree <- phy_tree(phylo.obj)
	if(!is.rooted(tree)) {
		cat('Root the tree by midpointing ... \n')
		tree <- midpoint(tree)
	} 
	
	ind <- rowSums(otu.tab) > 0
	if (sum(!ind) != 0) {
		cat('Remove ', sum(!ind), ' OTUs with all 0s! \n')
		otu.tab <- otu.tab[ind, ] 
		otu.name <- otu.name[ind, ]
	}
	
	absent <- tree$tip.label[!(tree$tip.label %in% rownames(otu.tab))]
	if (length(absent) != 0) {
		cat('Drop OTUs not in the OTU table ...\n')
		tree <- drop.tip(tree, absent)
	}
	
	
	if (sum(rownames(otu.name) != rownames(otu.tab)) != 0) stop('Inconsistent OTU names!\n')
	if (sum(rownames(meta.dat) != colnames(otu.tab)) != 0) stop('Inconsistent sample names!\n')
	
	
# Create abundance list
	cat("Create taxa abundance list ...\n")
	abund.list <- list()
	hierachs <- c('Phylum', 'Class', 'Order', 'Family', 'Genus')
	for (hierach in hierachs) {	
		if (hierach != 'Phylum') {
			single.names <- otu.name[, hierach]
			#	single.names[grepl('unclassified', single.names, ignore.case=T)] <- paste0('Unclassified',substr(hierach, 1, 1))
			tax.family <- paste(otu.name[, 'Phylum'], single.names, sep=";")
			tax.family[grepl('unclassified', tax.family, ignore.case = T)] <- paste0('Unclassified_', hierach)
		} else {
			tax.family <- otu.name[, 'Phylum']
			tax.family[grepl('unclassified', tax.family, ignore.case = T)] <- 'Unclassified_Phylum'
		}
		family <- aggregate(otu.tab, by = list(tax.family), FUN = sum)
		rownames(family) <- family[, 1]
		family <- as.matrix(family[, -1])
		abund.list[[hierach]] <- family
	}
	
	if (species) {
		abund.list[['Species']] <- otu.tab
		rownames(abund.list[['Species']]) <- paste0("OTU", rownames(otu.tab), ":", otu.name[, 'Phylum'], ";", otu.name[, 'Genus'])
	}
	cat('Completed!\n')
	data.obj <- list(otu.tab = otu.tab, abund.list = abund.list, meta.dat = meta.dat, tree = tree,
			otu.name = otu.name, otu.name.full = otu.name, call=match.call())
	return(data.obj)
}


# New: 2016_12_12, separate rarafaction and load data
# Functional data will not be rarefied
rarefy_data <- function (data.obj, dep=NULL) {
	otu.tab.12 <- data.obj$otu.tab
	otu.name.12 <- data.obj$otu.name
	otu.name.full <- data.obj$otu.name.full
#	abund.list.12 <- data.obj$abund.list
	siz.factor <- data.obj$size.factor
	meta.dat <- data.obj$meta.dat
	
	cat('Rarefaction ...\n')
	if (is.null(dep)) {
		otu.tab.12 <- t(Rarefy(t(otu.tab.12))$otu.tab.rff)
	} else {
		otu.tab.12 <- t(Rarefy(t(otu.tab.12), dep)$otu.tab.rff)
	}	
	dep <- colSums(otu.tab.12)[1]
	cat("Depth ", dep, '\n')
	# Remove empty OTUs
	otu.ind <- rowSums(otu.tab.12) > 0  # rev:2016-06-28
	otu.tab.12 <- otu.tab.12[otu.ind, ]	
	otu.name.12 <- otu.name.12[otu.ind, ]
	otu.name.full <- otu.name.full[otu.ind, ]
	
	samIDs <- intersect(rownames(meta.dat), colnames(otu.tab.12))

	if (length(samIDs) < nrow(meta.dat)) {
		warning('Some samples were lost during rarefaction!\n')
	}
	
	meta.dat <- meta.dat[samIDs, ]
	otu.tab.12 <- otu.tab.12[, samIDs]
	
	data.obj$meta.dat <- meta.dat  # Rev: 2017_01_19
	data.obj$otu.tab <- otu.tab.12
	data.obj$otu.name <- otu.name.12
	data.obj$otu.name.full <- otu.name.full
#	data.obj$abund.list <- abund.list.12
	data.obj$rff <- TRUE
	data.obj$rff.depth <- dep
	cat('After rarefaction, the size factor is automatically set to be the rarefaction depth!\n')
	data.obj$norm.method <- 'TSS'
	data.obj$norm.level <- 'OTU'
	data.obj$size.factor <- colSums(otu.tab.12)
	
	# Futher subset if it contains functional data
    data.obj <- subset_data(data.obj, samIDs)
	data.obj$act.seq <- paste0(data.obj$act.seq, 'R')
	
	cat("Recreate taxa abundance list ...\n")
	data.obj <- aggregate_by_taxonomy(data.obj)

	return(data.obj)
}

rarefy_data2 <- function (data.obj, dep=NULL) {
  otu.tab.12 <- data.obj$otu.tab
  otu.name.12 <- data.obj$otu.name
  #otu.name.full <- data.obj$otu.name.full
  #	abund.list.12 <- data.obj$abund.list

  meta.dat <- data.obj$meta.dat
  
  cat('Rarefaction ...\n')
  if (is.null(dep)) {
    otu.tab.12 <- t(Rarefy(t(otu.tab.12))$otu.tab.rff)
  } else {
    otu.tab.12 <- t(Rarefy(t(otu.tab.12), dep)$otu.tab.rff)
  }	
  dep <- colSums(otu.tab.12)[1]
  cat("Depth ", dep, '\n')
  # Remove empty OTUs
  otu.ind <- rowSums(otu.tab.12) > 0  # rev:2016-06-28
  otu.tab.12 <- otu.tab.12[otu.ind, ]	
  otu.name.12 <- otu.name.12[otu.ind, ]
 # otu.name.full <- otu.name.full[otu.ind, ]
  
  samIDs <- intersect(rownames(meta.dat), colnames(otu.tab.12))
  
  if (length(samIDs) < nrow(meta.dat)) {
    cat('Some samples were lost during rarefaction!\n')
  }
  
  meta.dat <- meta.dat[samIDs, ]
  otu.tab.12 <- otu.tab.12[, samIDs]
  
  data.obj$meta.dat <- meta.dat  # Rev: 2017_01_19
  data.obj$otu.tab <- otu.tab.12
  data.obj$otu.name <- otu.name.12
  #data.obj$otu.name.full <- otu.name.full

  # Futher subset if it contains functional data
  data.obj <- subset_data(data.obj, samIDs)

  
  cat("Recreate taxa abundance list ...\n")
  data.obj <- aggregate_by_taxonomy2(data.obj)
  
  return(data.obj)
}

# New: 2018_08_22, rarefaction by strata
# Functional data will not be rarefied
rarefy_data_strata <- function (data.obj, strata1, strata2) {
	otu.tab.12 <- data.obj$otu.tab
	otu.name.12 <- data.obj$otu.name
	otu.name.full <- data.obj$otu.name.full
#	abund.list.12 <- data.obj$abund.list
	siz.factor <- data.obj$size.factor
	meta.dat <- data.obj$meta.dat
	
	cat('Rarefaction by strata ...\n')
	strata1 <- data.obj$meta.dat[, strata1]
	strata2 <- data.obj$meta.dat[, strata2]
	
	strata <- factor(paste(strata1, strata2))
	slevels <- levels(strata)
	
	for (slevel in slevels) {
		ind <- which(strata == slevel)
		if (length(ind) >= 2) {
			otu.tab.12[, ind] <- t(Rarefy(t(otu.tab.12[, ind]))$otu.tab.rff)
		}
	}

	# Remove empty OTUs
	otu.ind <- rowSums(otu.tab.12) > 0  # rev:2016-06-28
	otu.tab.12 <- otu.tab.12[otu.ind, ]	
	otu.name.12 <- otu.name.12[otu.ind, ]
	otu.name.full <- otu.name.full[otu.ind, ]
	
	samIDs <- intersect(rownames(meta.dat), colnames(otu.tab.12))
	
	if (length(samIDs) < nrow(meta.dat)) {
		warning('Some samples were lost during rarefaction!\n')
	}
	
	meta.dat <- meta.dat[samIDs, ]
	otu.tab.12 <- otu.tab.12[, samIDs]
	
	data.obj$meta.dat <- meta.dat  # Rev: 2017_01_19
	data.obj$otu.tab <- otu.tab.12
	data.obj$otu.name <- otu.name.12
	data.obj$otu.name.full <- otu.name.full
#	data.obj$abund.list <- abund.list.12
	data.obj$rff <- TRUE
	data.obj$rff.depth <- dep
	cat('After rarefaction, the size factor is automatically set to be the rarefaction depth!\n')
	data.obj$norm.method <- 'TSS'
	data.obj$norm.level <- 'OTU'
	data.obj$size.factor <- colSums(otu.tab.12)
	
	# Futher subset if it contains functional data
	data.obj <- subset_data(data.obj, samIDs)
	data.obj$act.seq <- paste0(data.obj$act.seq, 'R')

	cat("Recreate taxa abundance list ...\n")
    data.obj <- aggregate_by_taxonomy(data.obj)
	return(data.obj)
}

# New: 2017_04_18
# Rev: 2020_02_20 'species'
remove_otu <- function (data.obj, ind) {

	data.obj$otu.tab <- data.obj$otu.tab[ind, ]
	data.obj$otu.name <- data.obj$otu.name[ind, ]
	data.obj$otu.name.full <- data.obj$otu.name.full[ind, ]
	
# Re-create Abundance 
#	otu.name <- data.obj$otu.name
#	otu.tab <- data.obj$otu.tab
#	abund.list <- list()
#
#	for (hierach in setdiff(hierachs, 'Species')) {	
#		if (hierach != 'Phylum') {
#			single.names <- otu.name[, hierach]
#			#	single.names[grepl('unclassified', single.names, ignore.case=T)] <- paste0('Unclassified',substr(hierach, 1, 1))
#			tax.family <- paste(otu.name[, 'Phylum'], single.names, sep=";")
#			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- paste0('Unclassified_', hierach)
#		} else {
#			tax.family <- otu.name[, 'Phylum']
#			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- 'Unclassified_Phylum'
#		}
#		family <- aggregate(otu.tab, by=list(tax.family), FUN=sum)
#		rownames(family) <- family[, 1]
#		family <- as.matrix(family[, -1])
#		abund.list[[hierach]] <- family
#	}
#	
#	if ('Species' %in% hierachs) {
#		abund.list[['Species']] <- otu.tab
#		rownames(abund.list[['Species']]) <- paste0("OTU", substr(rownames(otu.tab), 1, max.nchar), ":", otu.name[, 'Phylum'], ";", otu.name[, 'Genus'])
#	}
#	data.obj$abund.list <- abund.list
	
	data.obj <- aggregate_by_taxonomy(data.obj, species = 'Species' %in% names(data.obj$abund.list))
	
	return(data.obj)
}

# New: 2017_10_17
update_name <- update_sample_name <- function (data.obj, new.name) {
	
	if (length(new.name) != nrow(data.obj$meta.dat)) stop('The number of sample do not agree!\n')
	if (length(new.name) != length(unique(new.name))) stop ('The new names have duplicates!\n')
	
	rownames(data.obj$meta.dat) <- new.name
	data.obj$abund.list <- lapply(data.obj$abund.list, function (x) {colnames(x) <- new.name; x})
	colnames(data.obj$otu.tab) <- new.name
    if(!is.null(data.obj$size.factor))  names(data.obj$size.factor) <- new.name
	
	return(data.obj)
}

# Rev: 2016_09_22
update_data <- update_meta_data <- function (data.obj, map.file, meta.sep='\t', quote="\"", comment="", ...) {
	cat("Load meta file...\n")
	if (is.character(map.file)) {
		if (grepl("csv$", map.file)) {
			meta.dat <- read.csv(map.file, header=T, check.names=F, row.names=1, comment=comment, quote=quote, ...)
		} else {
			meta.dat <- read.table(map.file, header=T, check.names=F, row.names=1, comment=comment, sep=meta.sep, quote=quote, ...)
		}
	} else {
		meta.dat <- map.file
	}
    
	data.obj$meta.dat <- meta.dat
	
	samIDs <- intersect(rownames(meta.dat), colnames(data.obj$otu.tab))
	if (length(samIDs) == 0)  stop('Sample names in the meta file and biom file differ?\n')
    data.obj <- subset_data(data.obj, samIDs)
	return(data.obj)
}

# When the otu.name has detailed species level inofrmation
aggregate_by_taxonomy <- function (data.obj, species = TRUE) {
	
	max.nchar <- 4

	hierachs <- c('Phylum', 'Class', 'Order', 'Family', 'Genus')
	otu.name.12 <- data.obj$otu.name
	otu.name.12[is.na(otu.name.12) | otu.name.12 %in% c('', 'NA')] <- 'unclassified'
	otu.tab.12 <- data.obj$otu.tab
	
	abund.list.12 <- list()
	for (hierach in hierachs) {	
		if (!(hierach %in% c('Phylum'))) {
			single.names <- otu.name.12[, hierach]
			#	single.names[grepl('unclassified', single.names, ignore.case=T)] <- paste0('Unclassified',substr(hierach, 1, 1))
			tax.family <- paste(otu.name.12[, 'Phylum'], single.names, sep=";")
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- paste0('Unclassified_', hierach)
		} else {
			tax.family <- otu.name.12[, hierach]
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- paste0('Unclassified_', hierach)
		}
		family <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
		rownames(family) <- family[, 1]
		family <- as.matrix(family[, -1])
		abund.list.12[[hierach]] <- family
	}
	
	if (species) {
		abund.list.12[['Species']] <- otu.tab.12
		hierachs <- c(hierachs, 'Species')
		rownames(abund.list.12[['Species']]) <- paste0('A', substr(rownames(otu.tab.12), 1, max.nchar), ":", otu.name.12[, 'Phylum'], ";", otu.name.12[, 'Genus'])
	}
	
	# still keep other levels
	if (is.null(data.obj$abund.list)) {
		data.obj$abund.list <- list()
	}
	for (hierach in hierachs) {
		data.obj$abund.list[[hierach]] <- abund.list.12[[hierach]]
	}
	
	return(data.obj)
}

# When the otu.name has detailed species level inofrmation
aggregate_by_taxonomy2 <- function (data.obj, species = TRUE) {
	if (!species) {
		hierachs <- c('Phylum', 'Class', 'Order', 'Family', 'Genus')
	} else {
		hierachs <- c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
	}
	
	otu.name.12 <- data.obj$otu.name
	otu.name.12[is.na(otu.name.12) | otu.name.12 %in% c('', 'NA')] <- 'unclassified'
	otu.tab.12 <- data.obj$otu.tab
	
	abund.list.12 <- list()
	for (hierach in hierachs) {	
		if (!(hierach %in% c('Phylum', 'Species'))) {
			single.names <- otu.name.12[, hierach]
			#	single.names[grepl('unclassified', single.names, ignore.case=T)] <- paste0('Unclassified',substr(hierach, 1, 1))
			tax.family <- paste(otu.name.12[, 'Phylum'], single.names, sep=";")
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- paste0('Unclassified_', hierach)
		} else {
			tax.family <- otu.name.12[, hierach]
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- paste0('Unclassified_', hierach)
		}
		family <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
		rownames(family) <- family[, 1]
		family <- as.matrix(family[, -1])
		abund.list.12[[hierach]] <- family
	}
	# still keep other levels
	if (is.null(data.obj$abund.list)) {
		data.obj$abund.list <- list()
	}
	for (hierach in hierachs) {
		data.obj$abund.list[[hierach]] <- abund.list.12[[hierach]]
	}
	
	return(data.obj)
}

# New: 2017_02_16  aggregate data
# New: 2017_11_02  add more variable
aggregate_data <- function (data.obj, subject, strata = NULL) {
	cat('Aggregate data by a factor ...\n')
	subject <- (data.obj$meta.dat[,  subject])
	
	if (!is.null(strata)) {
		grp <- (data.obj$meta.dat[,  strata])
	}

	data.obj.new <- list()
	abund.list<- list()
	for (name in names(data.obj$abund.list)) {
		abund <- data.obj$abund.list[[name]]
		if (is.null(strata)) {
			obj <- aggregate(t(abund), by=list(subject), sum)
			abund.list[[name]] <- t(as.matrix(obj[, -1]))
			colnames(abund.list[[name]]) <- obj[, 1]
		} else {
			obj <- aggregate(t(abund), by=list(subject, grp), sum)
			abund.list[[name]] <- t(as.matrix(obj[, -(1:2)]))
			colnames(abund.list[[name]]) <- paste(obj[, 1], obj[, 2], sep="_")
		}
	}
	data.obj.new$abund.list <- abund.list
	
	abund <- data.obj$otu.tab
	if (is.null(strata)) {
		obj <- aggregate(t(abund), by=list(subject), sum)
		otu.tab <- t(as.matrix(obj[, -1]))
		colnames(otu.tab) <- obj[, 1]
	} else {
		obj <- aggregate(t(abund), by=list(subject, grp), sum)
		otu.tab <- t(as.matrix(obj[, -(1:2)]))
		colnames(otu.tab) <- paste(obj[, 1], obj[, 2], sep="_")
	}

	data.obj.new$otu.tab <- otu.tab
	
	data.obj.new$tree <- data.obj$tree
	data.obj.new$otu.name <- data.obj$otu.name
	data.obj.new$otu.name.full <- data.obj$otu.name.full
	
	if (is.null(strata)) {
		unique.ct <- aggregate(data.obj$meta.dat, by=list(subject), function(x) {
					length(unique(x))
				})
		
		unique.ct <- as.matrix(unique.ct[, -1])
		ind <- colSums(unique.ct) == nlevels(factor(subject))
		
		meta.dat <- aggregate(data.obj$meta.dat, by=list(subject), function(x) {
					x[1]
				})
		rownames(meta.dat) <- meta.dat[, 1]
		meta.dat <- meta.dat[, -1]
		meta.dat <- meta.dat[, ind]
	} else {
		unique.ct <- aggregate(data.obj$meta.dat, by=list(subject, grp), function(x) {
					length(unique(x))
				})
		
		unique.ct <- as.matrix(unique.ct[, -(1:2)])
		ind <- colSums(unique.ct) == nrow(unique.ct)
		
		meta.dat <- aggregate(data.obj$meta.dat, by=list(subject, grp), function(x) {
					x[1]
				})
		rownames(meta.dat) <- paste(meta.dat[, 1], meta.dat[, 2], sep="_")
		meta.dat <- meta.dat[, -(1:2)]
		meta.dat <- meta.dat[, ind]
	}

	data.obj.new$meta.dat <- meta.dat
	data.obj.new$winsor <- data.obj$winsor
	data.obj.new$winsor.qt <- data.obj$winsor.qt
	data.obj.new$rff <- data.obj$rff
	data.obj.new$act.seq <- paste0(data.obj$act.seq, 'A')

	cat('Finished!\n')
	return(data.obj.new)
}

# New: 2017_06_01 Remove reads belonging to specific taxa
remove_taxa <- function (data.obj, taxa.level='Genus', taxa.names) {
	# Remove Methanothermococcus
	if (taxa.level == 'OTU') {
		if (is.numeric(taxa.names)) warning('taxa.names should be characters!\n')
		ind <- !(rownames(data.obj$otu.tab) %in% paste(taxa.names))
		if (sum(!ind) == 0) {
			warning('Do not detect the given OTUs! Check the OTU names!\n')
			return(data.obj)
		}
	} else {
		if (!taxa.level %in% colnames(data.obj$otu.name)) {
			stop("Can't find the specified taxa levels!\n")
		}
		ind <- !(data.obj$otu.name[, taxa.level] %in% taxa.names)
		
		if (sum(!ind) == 0) {
			warning('Do not detect the given taxa! Check the taxa level!\n')
			return(data.obj)
		}
	}

	cat(sum(!ind), ' OTUs will be removed!\n')
	
	data.obj$otu.tab <- data.obj$otu.tab[ind, ]
	data.obj$otu.name <- data.obj$otu.name[ind, ]
	data.obj$otu.name.full <- data.obj$otu.name.full[ind, ]
	
# Re-create Abundance 
	otu.name <- data.obj$otu.name
	otu.tab <- data.obj$otu.tab
	abund.list <- list()
	hierachs <- c('Phylum', 'Class', 'Order', 'Family', 'Genus')
	for (hierach in hierachs) {	
		if (hierach != 'Phylum') {
			single.names <- otu.name[, hierach]
			#	single.names[grepl('unclassified', single.names, ignore.case=T)] <- paste0('Unclassified',substr(hierach, 1, 1))
			tax.family <- paste(otu.name[, 'Phylum'], single.names, sep=";")
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- paste0('Unclassified_', hierach)
		} else {
			tax.family <- otu.name[, 'Phylum']
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- 'Unclassified_Phylum'
		}
		family <- aggregate(otu.tab, by=list(tax.family), FUN=sum)
		rownames(family) <- family[, 1]
		family <- as.matrix(family[, -1])
		abund.list[[hierach]] <- family
	}
	if ('Species' %in% names(data.obj$abund.list)) {
		abund.list[['Species']] <- data.obj$abund.list[['Species']][ind, ]
	}
	data.obj$abund.list <- abund.list
	# Erase size factor
	data.obj$size.factor <- NULL
	data.obj$norm.method <- NULL
	data.obj$norm.level <- NULL
	if (!is.na.null(data.obj$rff)) {
		if (data.obj$rff == TRUE) {
			warning('The function is not intended to be applied to rarefied object! Rarefaction should be performed again!\n')
			data.obj$rff <- FALSE
			data.obj$rff.dep <- NULL
		}
	}
    data.obj$act.seq <- NULL
	cat('Remove taxa Finished!\n')
	
	ind <- colSums(data.obj$otu.tab) != 0
	
	if (sum(!ind) != 0) {
		cat(sum(!ind), ' samples contain no reads and they are removed!')
		data.obj <- subset_data(data.obj, ind)
	}
	return(data.obj)
}
# New: 2016_12_12
# intersect.no=4, winsor.qt=0.97
calculate_size <- function (data.obj, level='OTU', norm='GMPR', ...) {
	if (!is.null(data.obj$rff)) {
		if (data.obj$rff)
		warning("The data has been rarefied! Better calculate the size factor for unrarefied data!\n")
	}
	if (!is.null(data.obj$winsor)) {
		if (data.obj$winsor) {
			warning("The data has been winsorized! Better calculate the size factor before winsorization!\n")
		}
	}

	if (level == 'OTU') {
		data <- data.obj$otu.tab
	} else {
		if (level %in% names(data.obj$abund.list)) {
			data <- data.obj$abund.list[[level]]
		} else {
			data <- data.obj$otu.tab
			level <-'OTU'
			warning('No or wrong level specified! OTU level will be used!\n')
		}
	}
	if (norm == 'TSS') {
		sf <- colSums(data)
		norm <- 'TSS'
	}
	if (norm == 'GMPR') {
		sf <- GMPR(data, ...)
#		warning('GMPR is only suitable for samples from a same body location!\n')
		norm <- 'GMPR'
	}
	names(sf) <- colnames(data)
	data.obj$size.factor <- sf
	data.obj$norm.method <- norm
	data.obj$norm.level <- level
	data.obj$act.seq <- paste0(data.obj$act.seq, 'N')
	return(data.obj)
}

# Rev: 2016_12_08
winsor_data <- function (data.obj, winsor.qt=0.97) {
	if (data.obj$rff) {
		warning("The data has been rarefied! Better winsorize without rarefaction!\n")
	}
    cat('Winsorize ...\n')
	if (is.null(winsor.qt)) {
		winsor.qt <- 0.97
	}
	otu.tab.12 <- data.obj$otu.tab
	abund.list.12 <- data.obj$abund.list
	
    # Rev: 2016_12_08
	if (is.null(data.obj$size.factor)) {
		sf <- colSums(otu.tab.12)
		data.obj$norm.method <- 'TSS'
		data.obj$norm.level <- 'OTU'
		data.obj$size.factor <- sf
		cat('Size factor not available! Calculating using total sum! You may consider using better method such as GMPR!\n')
		data.obj$act.seq <- paste0(data.obj$act.seq, 'N')
	} else {
		sf <- data.obj$size.factor
	}
	
	# Addressing the outlier (97% percent) or at least one outlier
	abund.list.12 <- lapply(abund.list.12, function(genus) {
				genus.p <- t(t(genus) / sf)
				genus.p <- apply(genus.p, 1, function(x) {
							cutoff <- quantile(x, winsor.qt)
							x[x >= cutoff] <- cutoff
							x
						}
				)
				# column/row switch
				genus.w <- t(round(genus.p * sf))
			})
	# OTU table
	otu.tab.12.p <- t(t(otu.tab.12) / sf)
	otu.tab.12.p <- apply(otu.tab.12.p, 1, function(x) {
				cutoff <- quantile(x, winsor.qt)
				x[x >= cutoff] <- cutoff
				x
			}
	)
	# column/row switch
	otu.tab.12 <- t(round(otu.tab.12.p * sf))
	data.obj$winsor <- TRUE
	data.obj$winsor.qt <- winsor.qt
    data.obj$otu.tab <- otu.tab.12
	data.obj$abund.list <- abund.list.12
	data.obj$act.seq <- paste0(data.obj$act.seq, 'W')
	#data.obj$size.factor <- sf
	return(data.obj)
}

# Rev: 2016_11_28, Bray-curtis use normalized data.
construct_distance <- function (data.obj, unifrac.file=NULL,  Phylum='All', dist.RData=NULL, save.RData=NULL, 
		filter.no=0, rff=FALSE, dep=NULL, seed=1234) {
	set.seed(seed)
	
	if (!is.null(dist.RData)) {
		load(dist.RData, envir=.GlobalEnv )
	} else {
		dist.list.12 <- list()
		cat("Generalized UniFrac ...\n")
		
		otu.tab <- t(data.obj$otu.tab)
		
		if (rff == TRUE) {
			if (is.null(dep)) {
				otu.tab <- Rarefy(otu.tab)$otu.tab.rff
			} else {
				otu.tab <- Rarefy(otu.tab, dep)$otu.tab.rff
			}	
		}
				
		if (Phylum != 'All') {
			ind <- data.obj$otu.name[, 'Phylum'] == Phylum
			otu.tab <- otu.tab[, ind]
		}
		
		# Filter otus with reads <= filter.no
		otu.tab <- otu.tab[, colSums(otu.tab) > filter.no]
        
		# Remove samples with no reads
		if (sum(rowSums(otu.tab) == 0) >= 1) {
			otu.tab <- otu.tab[rowSums(otu.tab) != 0, ]
			warning('Some samples do not have reads after rarefaction! Please be careful!\n')
		}
		
		# To make sure the OTUs in otu.tab are in the tree (rev:2016-06-19)
		
		if (sum(!(colnames(otu.tab) %in% data.obj$tree$tip.label))) {
			warning('Some OTU names are not in the tree! An intersection set will be used!\n')	
		}
		common.otus <- intersect(colnames(otu.tab), data.obj$tree$tip.label)
		unifrac12 <- GUniFrac(otu.tab[, common.otus], data.obj$tree)$unifracs
	
		dist.list.12[['WUniFrac']] <- unifrac12[, , 'd_1']
		dist.list.12[['GUniFrac']] <- unifrac12[, , 'd_0.5']
		if (is.null(unifrac.file)) {
			dist.list.12[['UniFrac']] <- unifrac12[, , 'd_UW']
		} else {
			# The orders may be different
			dist.list.12[['UniFrac']] <- as.matrix(read.table(unifrac.file, row.names=1, header=T)) # Rarefaction
		}

		# Need speed up
		# Suggest using rarefied counts 
	    # If case/control has different sequencing depth, it will result in false clustering! 
        # Rev: 2016_11_28 Use normalized data for BC distance calculation to reduce noise in case of variable library size
		dist.list.12[['BC']] <-as.matrix(vegdist(otu.tab / rowSums(otu.tab)))
		
		genus <- t(data.obj$abund.list[['Genus']])
		genus <- genus / rowSums(genus)
		dist.list.12[['Euc']] <- as.matrix(dist(genus))
		
		genus <- sqrt(genus)
		dist.list.12[['Hel']] <-as.matrix(dist(genus))
#		dist.list.12[['JS']] <- as.matrix(distance(otu_table(data.obj$abund.list[['Genus']], taxa_are_rows=T), method='jsd'))
		
		if (!is.null(save.RData)) {
			save(dist.list.12, file=save.RData)
		}
	}

    return(dist.list.12)
}


## add 09172023 for if there is no tree
construct_distance2 <- function (data.obj, unifrac.file=NULL,  Phylum='All', dist.RData=NULL, save.RData=NULL, is.tree =T,
                                filter.no=0, rff=FALSE, dep=NULL, seed=1234) {
  set.seed(seed)
  
  if (!is.null(dist.RData)) {
    load(dist.RData, envir=.GlobalEnv )
  } else {
    dist.list.12 <- list()
    cat("Generalized UniFrac ...\n")
    
    otu.tab <- t(data.obj$otu.tab)
    
    if (rff == TRUE) {
      if (is.null(dep)) {
        otu.tab <- Rarefy(otu.tab)$otu.tab.rff
      } else {
        otu.tab <- Rarefy(otu.tab, dep)$otu.tab.rff
      }	
    }
    
    if (Phylum != 'All') {
      ind <- data.obj$otu.name[, 'Phylum'] == Phylum
      otu.tab <- otu.tab[, ind]
    }
    
    # Filter otus with reads <= filter.no
    otu.tab <- otu.tab[, colSums(otu.tab) > filter.no]
    
    # Remove samples with no reads
    if (sum(rowSums(otu.tab) == 0) >= 1) {
      otu.tab <- otu.tab[rowSums(otu.tab) != 0, ]
      warning('Some samples do not have reads after rarefaction! Please be careful!\n')
    }
    
    # To make sure the OTUs in otu.tab are in the tree (rev:2016-06-19)
    
    if (sum(!(colnames(otu.tab) %in% data.obj$tree$tip.label))) {
      warning('Some OTU names are not in the tree! An intersection set will be used!\n')	
    }
    
    if(is.tree==T){
      common.otus <- intersect(colnames(otu.tab), data.obj$tree$tip.label)
      
      unifrac12 <- GUniFrac(otu.tab[, common.otus], data.obj$tree)$unifracs
      
      dist.list.12[['WUniFrac']] <- unifrac12[, , 'd_1']
      dist.list.12[['GUniFrac']] <- unifrac12[, , 'd_0.5']
      if (is.null(unifrac.file)) {
        dist.list.12[['UniFrac']] <- unifrac12[, , 'd_UW']
      } else {
        # The orders may be different
        dist.list.12[['UniFrac']] <- as.matrix(read.table(unifrac.file, row.names=1, header=T)) # Rarefaction
      }
    }

    
    # Need speed up
    # Suggest using rarefied counts 
    # If case/control has different sequencing depth, it will result in false clustering! 
    # Rev: 2016_11_28 Use normalized data for BC distance calculation to reduce noise in case of variable library size
    dist.list.12[['BC']] <-as.matrix(vegdist(otu.tab / rowSums(otu.tab)))
    
    genus <- t(data.obj$abund.list[['Genus']])
    genus <- genus / rowSums(genus)
    dist.list.12[['Euc']] <- as.matrix(dist(genus))
    
    genus <- sqrt(genus)
    dist.list.12[['Hel']] <-as.matrix(dist(genus))
    #		dist.list.12[['JS']] <- as.matrix(distance(otu_table(data.obj$abund.list[['Genus']], taxa_are_rows=T), method='jsd'))
    
    if (!is.null(save.RData)) {
      save(dist.list.12, file=save.RData)
    }
  }
  
  return(dist.list.12)
}



outlier_detect <- function (data.obj, dist.obj, min.dep=2000) {
	# Future development
	samIDs <- colnames(data.obj$otu.tab)[colSums(data.obj$otu.tab) >= 2000]
	return(samIDs)
}

is.na.null <- function (x) {
	if (is.null(x)) {
		return(TRUE)
	} else {
		if (is.na(x)[1]) {
			return(TRUE)
		}  else {
			return(FALSE)
		}
	}
	
}

# Rev: 2016_09_26 remove empty OTUs/taxa
# Rev: 2016_12_01 add more logical controls
subset_data <- function (data.obj, samIDs) {
	
	# Rev: 2016_1_19 to add error protection
	# Transform logical samIDs into characer samIDs
	if (is.logical(samIDs) | is.numeric(samIDs)) {
		samIDs <- rownames(data.obj$meta.dat)[samIDs]
	}
	
	data.obj$meta.dat <- data.obj$meta.dat[samIDs, , drop=FALSE]

	if (!is.na.null(data.obj$otu.tab)) {

		data.obj$otu.tab <- data.obj$otu.tab[, samIDs, drop=FALSE]
		data.obj$otu.tab <- data.obj$otu.tab[rowSums(data.obj$otu.tab) != 0, , drop=FALSE]
		data.obj$otu.name <- data.obj$otu.name[rownames(data.obj$otu.tab), , drop=FALSE]

		if (!is.na.null(data.obj$otu.name.full)) {
			data.obj$otu.name.full <- data.obj$otu.name.full[rownames(data.obj$otu.tab), , drop=FALSE]
		}
	}

	if (!is.na.null(data.obj$abund.list)) {
		data.obj$abund.list <- lapply(data.obj$abund.list, function(x) {
					xx <- x[, samIDs, drop=FALSE]
					xx <- xx[rowSums(xx) != 0, , drop=FALSE]
				})
	}

	if (!is.na.null(data.obj$size.factor)) {
		data.obj$size.factor <- data.obj$size.factor[samIDs]
	}

	if (!is.na.null(data.obj$ko.list)) {
		data.obj$ko.list <- lapply(data.obj$ko.list, function(x) {
					xx <- x[, samIDs, drop=FALSE]
					xx <- xx[rowSums(xx) != 0, , drop=FALSE]
				})
	}

	if (!is.na.null(data.obj$cog.list)) {
		data.obj$cog.list <- lapply(data.obj$cog.list, function(x) {
					xx <- x[, samIDs, drop=FALSE]
					xx <- xx[rowSums(xx) != 0, , drop=FALSE]
				})
	}
	return(data.obj)
}

# Rev: 2016_12_01 add more logical controls
# Rev: 2016_02_01 fix one error
subset_dist <- function (dist.obj, samIDs) {
	
	# Rev: 2016_1_19 to add error protection
	# Transform logical samIDs into character samIDs
	if (is.logical(samIDs) | is.numeric(samIDs)) {
		samIDs <- rownames(dist.obj[[1]])[samIDs]
	}
	
	lapply(dist.obj, function(x) {
				if(!is.na.null(x)){
					x <- x[samIDs, samIDs]
				} else {
					x
				}
				x
			})
}

# New: 2018_06_10
generate_detailed_taxonomy <- function (data.obj, finest = 'Species', unclassified.symbol='unclassified') {
	# Create abundance list

	otu.name.12 <- data.obj$otu.name
	otu.name.12[otu.name.12 == '' | is.na(otu.name.12)] <- unclassified.symbol 
	otu.tab.12 <- data.obj$otu.tab
	
	if (finest == 'Species') {
		hierachs <- rev(c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
	} else {
		hierachs <- rev(c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus'))
	}
	
	abund.list.12 <- NULL
	otu.tab <- otu.tab.12
	otu.name <- otu.name.12
	for (hierach in hierachs) {	
		single.names <- otu.name[, hierach]
		if (hierach == 'Species' ) {
			tax.family <- paste0(otu.name[, 'Phylum'], ';', single.names)
		} 
		if (hierach %in% c('Phylum', 'Kingdom')) {
			tax.family <- single.names
		} 
		
		if (hierach %in% c('Class', 'Order', 'Family', 'Genus')) {
			tax.family <- single.names
			tax.family <- paste0(otu.name[, 'Phylum'], ';', single.names)
		} 
		
		if (hierach == 'Kingdom') {
			unclassified.ind <- rep(FALSE, length(tax.family))
			classified.ind <- !unclassified.ind
		} else {
			unclassified.ind <- grepl(unclassified.symbol, tax.family, ignore.case=T)
			classified.ind <- !unclassified.ind
		}

		otu.tab0 <- otu.tab[classified.ind, , drop=FALSE]
		tax.family0 <- tax.family[classified.ind, drop=FALSE]
		
		family <- aggregate(otu.tab0, by=list(tax.family0), FUN=sum)
		rownames(family) <- family[, 1]
		family <- as.matrix(family[, -1])
		abund.list.12 <- rbind(abund.list.12, family)

		if (sum(unclassified.ind) != 0) {
			otu.tab <- otu.tab[unclassified.ind, , drop=FALSE]
			otu.name <- otu.name[unclassified.ind, , drop=FALSE]
		} else{
			break
		}

	}
	data.obj$abund.list[['Detailed']] <- abund.list.12
	return(data.obj)
}


# 2020_01_25, allows zero-read samples and allows stratification by grp.name
perform_sequence_stat_analysis2 <- function (data.obj,  ann='', grp.name = NULL, prev.cutoff = 0.9, abund.cutoff = 0.5) {
	
	
	grp <- factor(data.obj$meta.dat[, 'grp.name'])
	levs <- levels(grp)
	
	sink(paste0('Sequence_Analysis_Statistics_', ann, '.txt'))
	otu.tab <- data.obj$otu.tab
	
	# Sequencing depth
	otu.abund <- rowSums(otu.tab)
	sam.abund <- colSums(otu.tab)
	otu.prev <- rowSums(otu.tab != 0) 
	sam.prev <- colSums(otu.tab != 0)
	
	# May keep 0 read sample - e.g. control samples withour reads are also informative
	otu.abund <- otu.abund[otu.abund >= 1]
#	sam.abund <- sam.abund[sam.abund >= 1]
	
	cat('This data set contains ', length(sam.abund), ' samples.\n')
	cat('16S rDNA targeted sequencing yields  a median of ', median(sam.abund), 'reads/sample  (range:', min(sam.abund), '-', max(sam.abund), ').\n')
	cat('Clustering of these 16S sequence tags produces ', sum(otu.abund > 0), ' OTUs at 97% similarity level.\n')
	
	
	phy.abund <- data.obj$abund.list[['Phylum']]
	fam.abund <- data.obj$abund.list[['Family']]
	gen.abund <- data.obj$abund.list[['Genus']]
	
	phy.prev <- rowSums(phy.abund != 0) / ncol(phy.abund)
	fam.prev <- rowSums(fam.abund != 0) / ncol(phy.abund)
	gen.prev <- rowSums(gen.abund != 0) / ncol(phy.abund)
	
	phy.abund <- rowMeans(t(t(phy.abund) / sam.abund))
	fam.abund <- rowMeans(t(t(fam.abund) / sam.abund))
	gen.abund <- rowMeans(t(t(gen.abund) / sam.abund))
	
	cat('These OTUs belong to ', sum(phy.abund > 0), ' phyla,', sum(fam.abund > 0), ' families and ', sum(gen.abund > 0), 'genera.\n\n')
	
	phy.prev <- sort(phy.prev, decr=T)
	phy.prev <- round(phy.prev[phy.prev >= prev.cutoff] * 100, 2)
	
	fam.prev <- sort(fam.prev, decr=T)
	fam.prev <- round(fam.prev[fam.prev >= prev.cutoff] * 100, 2)
	
	gen.prev <- sort(gen.prev, decr=T)
	gen.prev <- round(gen.prev[gen.prev >= prev.cutoff] * 100, 2)
	
	# Rev: 2017_02_19 ' ' -> '\n'
	cat('\nThe most prevalent phyla are:\n', paste(paste0(names(phy.prev), '(', phy.prev, '%)'), collapse='\n'), '\n')
	cat('\nThe most prevalent families are:\n', paste(paste0(names(fam.prev), '(', fam.prev, '%)'), collapse='\n'), '\n')
	cat('\nand the most prevalent genera are:\n', paste(paste0(names(gen.prev), '(', gen.prev, '%)'), collapse='\n'), '\n\n')
	
	phy.abund <- sort(phy.abund, decr=T)
	phy.abund <- round(phy.abund[phy.abund >= abund.cutoff] * 100, 2)
	
	fam.abund <- sort(fam.abund, decr=T)
	fam.abund <- round(fam.abund[fam.abund >= abund.cutoff] * 100, 2)
	
	gen.abund <- sort(gen.abund, decr=T)
	gen.abund <- round(gen.abund[gen.abund >= abund.cutoff] * 100, 2)
	
	cat('\nThe most abundant phyla are:\n', paste(paste0(names(phy.abund), '(', phy.abund, '%)'), collapse='\n'), '\n')
	cat('\nThe most abundant families are:\n', paste(paste0(names(fam.abund), '(', fam.abund, '%)'), collapse='\n'), '\n')
	cat('\nand the most abundant genera are:\n', paste(paste0(names(gen.abund), '(', gen.abund, '%)'), collapse='\n'), '\n\n')
	sink()
	
	pdf(paste0('Sequence_Analysis_Statistics_', ann, '.pdf'), height=5, width=5)
	obj <- ggplot2::ggplot(data=data.frame(x=otu.abund), aes(x=x)) + geom_histogram(col='black', fill='gray') + ylab('Frequency') + xlab('Abundance(Total counts)') +
			scale_x_log10(breaks=c(1, 10, 100, 1000, 10000, 100000, 100000))
	print(obj)
	obj <- ggplot2::ggplot(data=data.frame(x=sam.abund), aes(x=x)) + geom_histogram(col='black', fill='gray')  + ylab('Frequency') + xlab('Sequencing depth')
	print(obj)
	obj <- ggplot2::ggplot(data=data.frame(x=otu.prev), aes(x=x))  + ylab('Frequency') + xlab('OTU prevalence(Occurence frequency)') + geom_histogram(col='black', fill='gray')
	print(obj)
	dev.off()
}

perform_sequence_stat_analysis <- function (data.obj,  ann='') {
	sink(paste0('Sequence_Analysis_Statistics_', ann, '.txt'))
	otu.tab <- data.obj$otu.tab
	
	# Sequencing depth
	otu.abund <- rowSums(otu.tab)
	sam.abund <- colSums(otu.tab)
	otu.prev <- rowSums(otu.tab!=0)/ncol(otu.tab)
	
	# May keep 0 read sample - e.g. control samples withour reads are also informative
	otu.abund <- otu.abund[otu.abund >= 1]
	sam.abund <- sam.abund[sam.abund >= 1]
	cat('This data set contains ', length(sam.abund), ' samples after quality controls.\n')
	cat('16S rDNA targeted sequencing yields ', median(sam.abund), 'reads/sample on average (range:', min(sam.abund), '-', max(sam.abund), ').\n')
	cat('Clustering of these 16S sequence tags produces ', sum(otu.abund > 0), ' OTUs at 97% similarity level.\n')
	
	pdf(paste0('Sequence_Analysis_Statistics_', ann, '.pdf'), height=5, width=5)
	obj <- ggplot2::ggplot(data=data.frame(x=otu.abund), aes(x=x)) + geom_histogram(col='black', fill='gray') + ylab('Frequency') + xlab('Abundance(Total counts)') +
			scale_x_log10(breaks=10^(0 : floor(log10(max(otu.abund))) + 1))
	print(obj)
	obj <- ggplot2::ggplot(data=data.frame(x=sam.abund), aes(x=x)) + geom_histogram(col='black', fill='gray')  + ylab('Frequency') + xlab('Sequencing depth') +
	  scale_x_log10(breaks=10^(0 : floor(log10(max(sam.abund))) + 1))
	print(obj)
	obj <- ggplot2::ggplot(data=data.frame(x=otu.prev), aes(x=x))  + ylab('Frequency') + xlab('Prevalence(Occurence frequency)') + geom_histogram(col='black', fill='gray')
	print(obj)
	dev.off()
	
	phy.abund <- data.obj$abund.list[['Phylum']]
	fam.abund <- data.obj$abund.list[['Family']]
	gen.abund <- data.obj$abund.list[['Genus']]
	spe.abund <- data.obj$otu.tab
	
	phy.prev <- rowSums(phy.abund != 0) / ncol(phy.abund)
	fam.prev <- rowSums(fam.abund != 0) / ncol(phy.abund)
	gen.prev <- rowSums(gen.abund != 0) / ncol(phy.abund)
	spe.prev <- rowSums(spe.abund != 0) / ncol(phy.abund)
	
	phy.abund <- rowMeans(t(t(phy.abund) / sam.abund))
	fam.abund <- rowMeans(t(t(fam.abund) / sam.abund))
	gen.abund <- rowMeans(t(t(gen.abund) / sam.abund))
	spe.abund <- rowMeans(t(t(spe.abund) / sam.abund))
	
	cat('These OTUs/ASVs/Species belong to ', sum(phy.abund > 0), ' phyla,', sum(fam.abund > 0), ' families and ', sum(gen.abund > 0), 'genera.\n\n')
	
	phy.prev <- sort(phy.prev, decr=T)
	phy.prev <- round(phy.prev[phy.prev >= 0.90] * 100, 2)
	
	fam.prev <- sort(fam.prev, decr=T)
	fam.prev <- round(fam.prev[fam.prev >= 0.90] * 100, 2)
	
	gen.prev <- sort(gen.prev, decr=T)
	gen.prev <- round(gen.prev[gen.prev >= 0.90] * 100, 2)
	
  spe.prev <- sort(spe.prev, decr=T)
	spe.prev <- round(spe.prev[spe.prev >= 0.90] * 100, 2)
	
	# Rev: 2017_02_19 ' ' -> '\n'
	cat('\nThe most prevalent phyla are:\n', paste(paste0(names(phy.prev), '(', phy.prev, '%)'), collapse='\n'), '\n')
	cat('\nThe most prevalent families are:\n', paste(paste0(names(fam.prev), '(', fam.prev, '%)'), collapse='\n'), '\n')
	cat('\nand the most prevalent genera are:\n', paste(paste0(names(gen.prev), '(', gen.prev, '%)'), collapse='\n'), '\n')
	cat('\nand the most prevalent species are:\n', paste(paste0(names(spe.prev), '(', spe.prev, '%)'), collapse='\n'), '\n\n')
	
	phy.abund <- sort(phy.abund, decr=T)
	phy.abund <- round(phy.abund[phy.abund >= 0.05] * 100, 2)
	
	fam.abund <- sort(fam.abund, decr=T)
	fam.abund <- round(fam.abund[fam.abund >= 0.05] * 100, 2)
	
	gen.abund <- sort(gen.abund, decr=T)
	gen.abund <- round(gen.abund[gen.abund >= 0.05] * 100, 2)
	
	spe.abund <- sort(spe.abund, decr=T)
	spe.abund <- round(spe.abund[spe.abund >= 0.01] * 100, 2)
	
	cat('\nThe most abundant phyla are:\n', paste(paste0(names(phy.abund), '(', phy.abund, '%)'), collapse='\n'), '\n')
	cat('\nThe most abundant families are:\n', paste(paste0(names(fam.abund), '(', fam.abund, '%)'), collapse='\n'), '\n')
	cat('\nand the most abundant genera are:\n', paste(paste0(names(gen.abund), '(', gen.abund, '%)'), collapse='\n'), '\n')
	cat('\nand the most abundant species are:\n', paste(paste0(names(spe.abund), '(', spe.abund, '%)'), collapse='\n'), '\n\n')
	sink()
}

# Retire, using 'arsenal' tableby
perform_demograph_analysis <- function (data.obj, grp.name) {
	obj <- summary(data.obj$meta.dat)
	write.csv(obj, "meta.data.summary.csv", quote=F)

	grp <- data.obj$meta.dat[, grp.name]
	res <- NULL
	pv.vec <- NULL
	if (is.factor(grp)) {
		for (var1 in setdiff(colnames(data.obj$meta.dat), grp.name)) {
			temp <- data.obj$meta.dat[, var1]
			res <- rbind(res, c("", var1, rep("", nlevels(grp) - 1)))
			res <- rbind(res, c("", levels(grp)))
			if (is.factor(temp)) {
				res <- rbind(res, cbind(levels(temp), table(temp, grp)))
				if (nlevels(temp) == 1) {
					pv <- NA
				} else {
					err <- try(
							pv <- formatC(fisher.test(table(temp, grp))$p.value, digit=3)
					)
					if (inherits(err, "try-error")) {
						pv <- NA
					}
				}
				res <- rbind(res, c('Fisher p', pv, rep("", nlevels(grp) - 1)))

			} else {
				res <- rbind(res, c('mean', aggregate(temp, by=list(grp), FUN='mean')[, 2]))
				res <- rbind(res, c('sd', aggregate(temp, by=list(grp), FUN='sd')[, 2]))
				
				err <- try(
						pv <- formatC(summary(aov(temp ~ grp))[[1]][1, 'Pr(>F)'], digit=3)
				)
				if (!inherits(err, "try-error")) {
					res <- rbind(res, c('ANOVA p', pv, rep("", nlevels(grp) - 1)))
				} else {
					res <- rbind(res, c('ANOVA p', 'NA', rep("", nlevels(grp) - 1)))
					pv <- NA
				}

			}
			res <- rbind(res, rep("", nlevels(grp)+1))
			pv.vec <- c(pv.vec, pv)
		}
		
	} else {

		
	}

	write.csv(res, "meta.data.by.grp.csv", row.names=F, quote=F)
	names(pv.vec) <- setdiff(colnames(data.obj$meta.dat), grp.name)
	return(pv.vec)
}


# Rev: 2017_08_23 add automatically create phylo.obj
generate_rarefy_curve <- function (data.obj, phylo.obj=NULL, grp.name, depth=NULL, npoint=10, iter.no=5,
		measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'), ann='', gg.cmd="theme(legend.justification=c(1,0), legend.position=c(1,0))", wid=5, hei=5) {
	cat("Create rarefaction curves!\n")
	# Rev: 2017_08_23
	if (is.null(phylo.obj)) {
		phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows=T), phy_tree(data.obj$tree), 
				tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))
	}
	if (is.null(depth)) {
		depth <- min(sample_sums(phylo.obj))
		phylo.even <- rarefy_even_depth(phylo.obj, rngseed=12345)
	} else {
		if (depth > min(sample_sums(phylo.obj))) {
			ind <- sample_sums(phylo.obj) >= depth
			cat(sum(!ind), " samples do not have sufficient number of reads!\n")
			sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ] 
			data.obj <- subset_data(data.obj, ind)
		}
		phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed=12345)
	}
	
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	# Rev: 2016_12_12
	if (is.character(grp)) {
		grp <- factor(grp)
	}

	if (!is.factor(grp)) {
		stop('Rarefaction curve needs a factor!\n')
	}
	res <- NULL
	incr <- depth %/% npoint
	sink('temp.txt')
	for (dep in c(10, incr*(1:npoint))) {
		x <- 0
		for (i in 1:iter.no) {
			phylo.even <- rarefy_even_depth(phylo.obj, dep, rngseed=12345+i)
			x <- x + estimate_richness(phylo.even, measures=measures)
		}
		
		res <- rbind(res, t(x[, measures, drop=F]/iter.no))
	}
	colnames(res) <- rownames(df)
	sink()
	
	pdf(paste0("Alpha_diversity_Rarefaction_Curve_", ann, ".pdf"), width=wid, height=hei)
	for (i in 1:length(measures)) {
		measure <- measures[i]
		cat("Measure: ", measure, "\n")
		res2 <- res[(0:(npoint))*length(measures)+i, , drop=F]
		m <- t(apply(res2, 1, function(x) tapply(x, grp, mean)))
		se <- t(apply(res2, 1, function(x) tapply(x, grp, function(y) sd(y)/sqrt(length(y)))))
		uci <- m+se
		lci <- m-se
		
		m <- melt(m)
		uci <- melt(uci)
		lci <- melt(lci)
		
		res2 <- cbind(c(10, incr*(1:npoint)), m[, 2:3], uci[, 3], lci[, 3])
		colnames(res2) <- c('Depth', 'Group', 'mean', 'max', 'min')
		
		res2 <- as.data.frame(res2)
		res2$Group <- factor(res2$Group, levels=levels(grp))
		
		#write.table(res2, paste0("Alpha_diversity_Rarefaction_", ann, "_", measure, ".txt"))
		
		obj <- ggplot(res2, aes(x=Depth, y=mean, color=Group, group=Group)) +
				geom_errorbar(aes(ymin=min, ymax=max), alpha=0.5, width=.25, position=position_dodge(.2)) + 
				geom_line() + 
				geom_point(size=3, shape=21, fill="white") +
				labs(y=measure)

		if (!is.null(gg.cmd)) {
			obj <- obj + eval(parse(text=gg.cmd))
		}
		
		print(obj)		
	}
	dev.off()
}

# Rev: 2016_12_25  Add anova, and record the results
# Rev: 2018_02_24  Remove dependence on Phyloseq
perform_alpha_test <- function (data.obj, phylo.obj=NULL, rarefy=TRUE, depth=NULL, iter.no=5, 
		measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'),  model='lm', 
		formula=NULL, grp.name=NULL, adj.name=NULL, ann='', seed=123, ...) {
	# Implement future random effects model
	# Rev: 2017_08_23
	if (is.null(phylo.obj)) {
		phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows=T), phy_tree(data.obj$tree), 
				tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))
	}
	
	result <- list()
	if (rarefy == TRUE) {
		if (is.null(depth)) {
			depth <- min(sample_sums(phylo.obj))
		} else {
			if (depth > min(sample_sums(phylo.obj))) {
				ind <- sample_sums(phylo.obj) >= depth
				cat(sum(!ind), " samples do not have sufficient number of reads!\n")
				sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ] 
				data.obj <- subset_data(data.obj, ind)
			}
		}
			
		x <- 0
		sink('temp.txt')
		for (i in 1:iter.no) {
			phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed=12345+i)
			x <- x + estimate_richness(phylo.even, measures=measures)
		}
		sink()
		x <- x / iter.no
	} else {
		x <- estimate_richness(phylo.obj, measures=measures)
	}

	df <- data.obj$meta.dat
	
	fitted.obj <- list()
	
	if (rarefy == T) {
		sink(paste0('Alpha_diversity_test_results_rarefied_', ann, '.txt'))
	} else {
		sink(paste0('Alpha_diversity_test_results_unrarefied_', ann, '.txt'))
	}
	
	date()
	
	# variables to adjust always come first in anova analyses
	if (is.null(formula)) {
		if (is.null(adj.name)) {
			formula <- paste('~', grp.name)
		} else {
			formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
		}
	}
	
	for (measure in measures) {
		cat("Alpha diversity:", measure, "\n")
		xx <- x[, measure]
		if (model == 'lm') {
			cat('Linear model:\n')

			lm.obj <- lm(as.formula(paste('xx ', formula)), df, ...)
			prmatrix(summary(lm.obj)$coefficients)
			cat('\nANOVA:\n')
			print(anova(lm.obj))
			cat('\n')
			fitted.obj[[measure]] <- lm.obj
		}
		if (model == 'lme') {
			df$xx <- xx
			cat('Linear mixed effects model:\n')
			lm.obj <- lme(as.formula(paste('xx ', formula)), df, method='ML', ...)
			prmatrix(summary(lm.obj)$tTable)
			cat('\nANOVA:\n')
			print(anova(lm.obj))
			cat('\n')
			fitted.obj[[measure]] <- lm.obj
		}
		cat("\n")
	}
	sink()
	
	result$fitted.obj <- fitted.obj
	result$alpha.diversity <- x
	# Rev: 2016_12_25
	return(invisible(result))

}

# New: 2018_02_25
estimate_richness_ <- function (OTU,  measures = NULL) {
	# Modified over Phyloseq
	
	renamevec = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", 
			"InvSimpson", "Fisher", 'Pielou')
	names(renamevec) <- c("S.obs", "S.chao1", "S.ACE", "shannon", 
			"simpson", "invsimpson", "fisher", 'Pielou')
	if (is.null(measures)) {
		measures = as.character(renamevec)
	}
	if (any(measures %in% names(renamevec))) {
		measures[measures %in% names(renamevec)] <- renamevec[names(renamevec) %in% 
						measures]
	}
	if (!any(measures %in% renamevec)) {
		stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
	}
	outlist = vector("list")
	estimRmeas = c("Chao1", "Observed", "ACE")
	if (any(estimRmeas %in% measures)) {
		outlist <- c(outlist, list(t(data.frame(vegan::estimateR(OTU)))))
	}
	if ("Shannon" %in% measures) {
		outlist <- c(outlist, list(shannon = vegan::diversity(OTU, index = "shannon")))
	}
	if ("Simpson" %in% measures) {
		outlist <- c(outlist, list(simpson = vegan::diversity(OTU, index = "simpson")))
	}
	if ("InvSimpson" %in% measures) {
		outlist <- c(outlist, list(invsimpson = vegan::diversity(OTU, 
								index = "invsimpson")))
	}
	
	if ("Fisher" %in% measures) {
		fisher = tryCatch(vegan::fisher.alpha(OTU, se = TRUE), warning = function(w) {
					warning("estimate_richness: Warning in fisher.alpha(). See `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution")
					suppressWarnings(vegan::fisher.alpha(OTU, se = TRUE)[, c("alpha", "se")])
				})
		if (!is.null(dim(fisher))) {
			colnames(fisher)[1:2] <- c("Fisher", "se.fisher")
			outlist <- c(outlist, list(fisher))
		}
		else {
			outlist <- c(outlist, Fisher = list(fisher))
		}
	}
	
	f <- function(ns) {
		sapply(ns, function (n) {
					
					i <- 1:n
					-sum(1 / i * log (1 / i))
				})
	}
	
	cat(names(outlist))
	if ("Pielou" %in% measures) {
		outlist <- c(outlist, list(Pielou = vegan::diversity(OTU, index = "shannon") / f(rowSums(OTU!=0))))
	}
	
	
	out = do.call("cbind", outlist)
	namechange = intersect(colnames(out), names(renamevec))
	colnames(out)[colnames(out) %in% namechange] <- renamevec[namechange]
	colkeep = sapply(paste0("(se\\.){0,}", measures), grep, colnames(out), 
			ignore.case = TRUE)
	out = out[, sort(unique(unlist(colkeep))), drop = FALSE]
	out <- as.data.frame(out)
	return(out)
}

# New: 2018_02_25
generate_alpha_diversity <- function (data.obj,  rarefy=TRUE, depth=NULL, iter.no=5,
		measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'), seed=123) {	
	# Rev: 2017_08_23

	OTU <- t(data.obj$otu.tab)
	depths <- rowSums(OTU)
	result <- list()
	if (rarefy == TRUE) {
		if (is.null(depth)) {
			depth <- min(depths)
		} else {
			if (depth > min(depths)) {
				ind <- depths >= depth
				cat(sum(!ind), " samples do not have sufficient number of reads!\n")
				# data.obj <- subset_data(data.obj, ind)
			}
		}

		set.seed(seed)
		x <- 0 
		for (i in 1:iter.no) {
			OTU2 <- Rarefy(OTU, depth = depth)$otu.tab.rff
			x <- x + estimate_richness_(OTU2, measures=measures)
		}
		x <- x / iter.no
		rownames(x) <- rownames(OTU2)
	} else {
		x <- estimate_richness_(OTU, measures=measures)
		rownames(x) <- rownames(OTU)
	}
	
	return(x)
}
	

# Rev: 2016_09_10
# Rev: 2016_11_28
# Rev: 2017_04_18
# Rev: 2018_03_06 position_jitter(h=0)
generate_alpha_boxplot <- function (data.obj, phylo.obj=NULL, rarefy=TRUE, depth=NULL, grp.name, strata=NULL, 
		measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'), gg.cmd=NULL, ann='', subject=NULL, p.size=2.5, l.size=0.5,
		hei = NULL, wid = NULL) {	
	# Rev: 2017_08_23
	if (is.null(phylo.obj)) {
		phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows=T), phy_tree(data.obj$tree), 
				tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))
	}
	# To be completed - jetter when strata is not null
	if (rarefy == TRUE) {
		if (is.null(depth)) {
			depth <- min(sample_sums(phylo.obj))
		} else {
			if (depth > min(sample_sums(phylo.obj))) {
				ind <- (sample_sums(phylo.obj) >= depth)
				cat(sum(!ind), " samples do not have sufficient number of reads!\n")
				
				sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ] 
				data.obj <- subset_data(data.obj, ind)
			}
		}

		phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed=12345)
		x <- estimate_richness(phylo.even, measures=measures)
	} else {
		x <- estimate_richness(phylo.obj, measures=measures)
	}
	
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	if (!is.null(subject)) {
		ID <- df[, subject]
	}
	
	if (is.null(hei) | is.null(wid))  {
		hei <- 5
		if (is.null(strata)) {
			wid <- 5
		} else {
			wid <- 5.5
		}
	}

	if (rarefy == T) {
		pdf(paste0('Alpha_diversity_boxplot_rarefied_', ann, '.pdf'), height=hei, width=wid)
	} else {
		pdf(paste0('Alpha_diversity_boxplot_unrarefied_', ann, '.pdf'), height=hei, width=wid)
	}
	if (is.null(subject)) {
		if (is.null(strata)) {
			for (measure in measures) {
				cat(measure, '\n')
				xx <- x[, measure]		
				df2 <- data.frame(Value=xx, Group=grp)
				dodge <- position_dodge(width=0.75)
				obj <- ggplot(df2, aes(x=Group, y=Value, col=Group)) +
						geom_boxplot(position=dodge,  outlier.colour = NA) + 
						geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1, h = 0)) +
						labs(y=measure) +
						theme(legend.position="none")
				if (!is.null(gg.cmd)) {
					obj <- obj + eval(parse(text=gg.cmd))
				}
				
				print(obj)
			}	
		} else {
			for (measure in measures) {
				cat(measure, '\n')
				xx <- x[, measure]		
				grp2 <- df[, strata]
				df2 <- data.frame(Value=xx, Group=grp, Strata=grp2)
				
				dodge <- position_dodge(width=0.9)
				obj <- ggplot(df2, aes(x=Strata, y=Value, col=Group, fill=Group)) +
						geom_jitter(alpha=0.6, size=3.0,  position = position_jitterdodge(dodge.width=0.9)) +
						geom_boxplot(position=dodge,  outlier.colour = NA, fill='white', alpha=0.6) + 
						labs(y=measure, x=strata) 
				if (!is.null(gg.cmd)) {
					obj <- obj + eval(parse(text=gg.cmd))
				}
				print(obj)
			}
			
		}
		dev.off()
	} else {
		if (is.null(strata)) {
			for (measure in measures) {
				cat(measure, '\n')
				xx <- x[, measure]		
				df2 <- data.frame(Value=xx, Group=grp, subject=ID)
				dodge <- position_dodge(width=0.75)
			
				obj <- ggplot(df2, aes(x=Group, y=Value, shape=Group,  group=subject)) +
						geom_point(size=p.size) +
						geom_line(size=l.size) +
						labs(y=measure) +
						theme(legend.position="none")
				if (!is.null(gg.cmd)) {
					obj <- obj + eval(parse(text=gg.cmd))
				}
				
				print(obj)
			}	
		} else {
			warning('Stratum is disallowed when subject is specified!\n')
		}
		dev.off()
	}

}

# Rev: 2018_03_09 Removing phyloseq dependency and add 'subject' parameter, and remove 'model' parameter 

generate_rarefy_curve2 <- function (data.obj,  grp.name, 
		measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'), depth=NULL,  iter.no=5, npoint=10, seed=123,
		ann='', gg.cmd="theme(legend.justification=c(1,0), legend.position=c(1,0))", wid=5, hei=5, error = 'se') {
	
	cat("Create rarefaction curves!\n")

	if (is.null(depth)) {
		depth <- min(colSums(data.obj$otu.tab))
	}
	
	if (error == 'se') k <- 1
	if (error == 'ci') k <- 1.96
	
	ind <- colSums(data.obj$otu.tab) >= depth
	if (sum(!ind) != 0) {
		cat(sum(!ind), " samples do not have sufficient number of reads!\n")
		data.obj <- subset_data(data.obj, ind)
	}
	
	df <- data.obj$meta.dat
	grp <- df[, grp.name]

	if (is.character(grp)) {
		grp <- factor(grp)
	}
	
	if (!is.factor(grp)) {
		stop('Rarefaction curve needs a factor!\n')
	}
	res <- NULL
	incr <- depth %/% npoint

	for (dep in c(10, incr*(1:npoint))) {
		x <- generate_alpha_diversity(data.obj,  rarefy=TRUE, depth=dep, iter.no=iter.no, measures=measures, seed=seed)
		res <- rbind(res, t(x[, measures, drop=F]))
	}
	colnames(res) <- rownames(df)

	cat('Wait ...\n')
	pdf(paste0("Alpha_diversity_Rarefaction_Curve_", ann, ".pdf"), width=wid, height=hei)
	for (i in 1:length(measures)) {
		measure <- measures[i]
		cat("Measure: ", measure, "\n")
		res2 <- res[(0:(npoint))*length(measures)+i, , drop=F]
		m <- t(apply(res2, 1, function(x) tapply(x, grp, mean)))
		se <- t(apply(res2, 1, function(x) tapply(x, grp, function(y) k * sd(y) / sqrt(length(y)))))
		uci <- m+se
		lci <- m-se
		
		m <- melt(m)
		uci <- melt(uci)
		lci <- melt(lci)
		
		res2 <- cbind(c(10, incr*(1:npoint)), m[, 2:3], uci[, 3], lci[, 3])
		colnames(res2) <- c('Depth', 'Group', 'mean', 'max', 'min')
		
		res2 <- as.data.frame(res2)
		res2$Group <- factor(res2$Group, levels=levels(grp))
		
		#write.table(res2, paste0("Alpha_diversity_Rarefaction_", ann, "_", measure, ".txt"))
		
		obj <- ggplot(res2, aes(x=Depth, y=mean, color=Group, group=Group)) +
				geom_errorbar(aes(ymin=min, ymax=max), alpha=0.5, width=.25, position=position_dodge(.2)) + 
				geom_line() + 
				geom_point(size=3, shape=21, fill="white") +
				labs(y=measure)
		
		if (!is.null(gg.cmd)) {
			obj <- obj + eval(parse(text=gg.cmd))
		}
		
		print(obj)		
	}
	dev.off()
	cat('Finished!\n')
}

# New: 2018_03_09 Removing phyloseq dependency and add 'subject' parameter, and remove 'model' parameter 
perform_alpha_test2 <- function (data.obj, alpha.obj=NULL, rarefy=TRUE, depth=NULL, iter.no=5, 
		measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'),  model='lm', 
		formula=NULL, grp.name=NULL, adj.name=NULL, subject=NULL, ann='', seed=123, ...) {
	if (is.null(alpha.obj)) {
		alpha.obj <- generate_alpha_diversity(data.obj,  rarefy=rarefy, depth=depth, iter.no=iter.no, measures=measures, seed=seed)
		
	} else {
		if (sum(!(rownames(alpha.obj) %in% rownames(data.obj$meta.dat))) != 0){
			stop("alpha.obj contains samples not in data.obj!\n")
		} 
	}
	
	
	x <- alpha.obj
	df <- data.obj$meta.dat[rownames(alpha.obj), ]
	
	if (is.null(subject)) {
		model <- 'lm'
	} else {
		model <- 'lme'
	}
	
	result <- list()
	fitted.obj <- list()
	
	if (rarefy == T) {
		sink(paste0('Alpha_diversity_test_results_rarefied_', ann, '.txt'))
	} else {
		sink(paste0('Alpha_diversity_test_results_unrarefied_', ann, '.txt'))
	}
	
	date()
	
	# variables to adjust always come first in anova analyses
	if (is.null(formula)) {
		if (is.null(adj.name)) {
			formula <- paste('~', grp.name)
		} else {
			formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
		}
	}
	
	for (measure in measures) {
		cat("Alpha diversity:", measure, "\n")
		xx <- x[, measure]
		if (model == 'lm') {
			cat('Linear model:\n')
			
			lm.obj <- lm(as.formula(paste('xx ', formula)), df, ...)
			prmatrix(summary(lm.obj)$coefficients)
			cat('\nANOVA:\n')
			print(anova(lm.obj))
			cat('\n')
			fitted.obj[[measure]] <- lm.obj
		}
		if (model == 'lme') {
			df$xx <- xx
			cat('Linear mixed effects model:\n')
			lm.obj <- lme(as.formula(paste('xx ', formula)), df, method='ML', random=as.formula(paste0(' ~ 1 | ', subject)), ...)
			prmatrix(summary(lm.obj)$tTable)
			cat('\nANOVA:\n')
			print(anova(lm.obj))
			cat('\n')
			fitted.obj[[measure]] <- lm.obj
		}
		cat("\n")
	}
	sink()
	
	result$fitted.obj <- fitted.obj
	result$alpha.diversity <- x
	# Rev: 2016_12_25
	return(invisible(result))
	
}


perform_alpha_test3 <- function (data.obj, alpha.obj=NULL, rarefy=TRUE, depth=NULL, iter.no=5, 
                                 measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'),  model='lm', 
                                 formula=NULL, grp.name=NULL, adj.name=NULL, subject=NULL, ann='', seed=123, ...) {
  if (is.null(alpha.obj)) {
    alpha.obj <- generate_alpha_diversity(data.obj,  rarefy=rarefy, depth=depth, iter.no=iter.no, measures=measures, seed=seed)
    
  } else {
    if (sum(!(rownames(alpha.obj) %in% rownames(data.obj$meta.dat))) != 0){
      stop("alpha.obj contains samples not in data.obj!\n")
    } 
  }
  
  
  x <- alpha.obj
  df <- data.obj$meta.dat[rownames(alpha.obj), ]
  
  if (is.null(subject)) {
    model <- 'lm'
  } else {
    model <- 'lme'
  }
  
  
  if (rarefy == T) {
    sink(paste0('Alpha_diversity_test_results_rarefied_', ann, '.txt'))
  } else {
    sink(paste0('Alpha_diversity_test_results_unrarefied_', ann, '.txt'))
  }
  
  date()
  
  # variables to adjust always come first in anova analyses
  if (is.null(formula)) {
    if (is.null(adj.name)) {
      formula <- paste('~', grp.name)
    } else {
      formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
    }
  }
  result <- list()
  fitted.obj <- list()
  res <- matrix(NA, nrow = length(measures), ncol = 3, dimnames = list(measures, c('P','R2','coef')))
  for (measure in measures) {
    cat("Alpha diversity:", measure, "\n")
    xx <- x[, measure]
    if (model == 'lm') {
        cat('Linear model:\n')
        lm.obj <- lm(as.formula(paste('xx ', formula)), df, ...)
        prmatrix(summary(lm.obj)$coefficients)
        cat('\nANOVA:\n')
        print(anova(lm.obj))
        cat('\n')
        fitted.obj[[measure]] <- lm.obj
        p <- anova(lm.obj)[grp.name,'Pr(>F)']
        if(!is.null(adj.name)){
          r2 <- partial.r.adj(x=df[,grp.name],z=df[,adj.name],y=xx)
        }else{
          r2 <- partial.r.adj(x=df[,grp.name],z=NULL,y=xx)
        }
        
        tmp <- summary(lm.obj)
        if(length(grep(grp.name,rownames(coef(tmp))))>1){
          coef <- NA
        }else{
          coef <- tmp$coefficients[,'Estimate']
          coef <- coef[grep(grp.name,names(coef))]
        }
        
    }
    if (model == 'lme') {
      df$xx <- xx
      cat('Linear mixed effects model:\n')
      lm.obj <- lme(as.formula(paste('xx ', formula)), df, method='ML', random=as.formula(paste0(' ~ 1 | ', subject)), ...)
      prmatrix(summary(lm.obj)$tTable)
      cat('\nANOVA:\n')
      print(anova(lm.obj))
      cat('\n')
      fitted.obj[[measure]] <- lm.obj
    }
    
    res[measure, 'P'] <- p
    res[measure, 'R2'] <- r2
    res[measure, 'coef'] <- coef
    result$fitted.obj <- fitted.obj
    result$alpha.diversity <- x
    cat("\n")
  }
  sink()
  
  # Rev: 2016_12_25
  return(list(result= result, res = res, fitted.obj = fitted.obj))
  
}

# New: 2018_03_09 Removing phyloseq dependency
generate_alpha_boxplot2 <- function (data.obj, alpha.obj=NULL,
		rarefy=TRUE, depth=NULL, iter.no=5, measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'),  seed  = 123, 
		grp.name, strata=NULL, gg.cmd=NULL, ann='', subject=NULL, p.size=2.5, l.size=0.5,
		hei = NULL, wid = NULL) {	
	# Rev: 2017_08_23
	if (is.null(alpha.obj)) {
		alpha.obj <- generate_alpha_diversity(data.obj,  rarefy=rarefy, depth=depth, iter.no=iter.no, measures=measures, seed=seed)

	} else {
		if (sum(!(rownames(alpha.obj) %in% rownames(data.obj$meta.dat))) != 0){
			stop("alpha.obj contains samples not in data.obj!\n")
		} 
	}

	x <- alpha.obj
    df <- data.obj$meta.dat[rownames(alpha.obj), ]
	grp <- df[, grp.name]
	if (!is.null(subject)) {
		ID <- df[, subject]
	}
	
	if (is.null(hei) | is.null(wid))  {
		hei <- 5
		if (is.null(strata)) {
			wid <- 5
		} else {
			wid <- 5.5
		}
	}
	
	if (rarefy == T) {
		pdf(paste0('Alpha_diversity_boxplot_rarefied_', ann, '.pdf'), height=hei, width=wid)
	} else {
		pdf(paste0('Alpha_diversity_boxplot_unrarefied_', ann, '.pdf'), height=hei, width=wid)
	}
	if (is.null(subject)) {
		if (is.null(strata)) {
			for (measure in measures) {
				cat(measure, '\n')
				xx <- x[, measure]		
				df2 <- data.frame(Value=xx, Group=grp)
				dodge <- position_dodge(width=0.75)
				obj <- ggplot(df2, aes(x=Group, y=Value, col=Group)) +
						geom_boxplot(position=dodge,  outlier.colour = NA, lwd=0.35) + 
						geom_jitter(alpha=0.6, size=p.size,  position = position_jitter(w = 0.1, h = 0)) +
						labs(y=measure) +
						theme(legend.position="none")
				if (!is.null(gg.cmd)) {
					obj <- obj + eval(parse(text=gg.cmd))
				}
				
				print(obj)
			}	
		} else {
			for (measure in measures) {
				cat(measure, '\n')
				xx <- x[, measure]		
				grp2 <- df[, strata]
				df2 <- data.frame(Value=xx, Group=grp, Strata=grp2)
				
				dodge <- position_dodge(width=0.75)
				obj <- ggplot(df2, aes(x=Strata, y=Value, col=Group, fill=Group)) +
						geom_boxplot(position=dodge,  outlier.colour = NA, fill='white', alpha=0.65, lwd=0.35, width=0.65) + 
						geom_point(position=position_jitterdodge(dodge.width=0.75), size=p.size, alpha=0.6) +
						labs(y=measure, x=strata) 
				if (!is.null(gg.cmd)) {
					obj <- obj + eval(parse(text=gg.cmd))
				}
				print(obj)
				
				
			}
			
		}
		dev.off()
	} else {
		if (is.null(strata)) {
			for (measure in measures) {
				cat(measure, '\n')
				xx <- x[, measure]		
				df2 <- data.frame(Value=xx, Group=grp, subject=ID)
				dodge <- position_dodge(width=0.75)
				
				obj <- ggplot(df2, aes(x=Group, y=Value, shape=Group,  group=subject)) +
						geom_point(size=p.size) +
						geom_line(size=l.size) +
						labs(y=measure) +
						theme(legend.position="none")
				if (!is.null(gg.cmd)) {
					obj <- obj + eval(parse(text=gg.cmd))
				}
				
				print(obj)
			}	
		} else {
			warning('Stratum is disallowed when subject is specified!\n')
		}
		dev.off()
	}
	
}

# New: 2018_03_08  Removing phyloseq dependency
# Rev: 2018_07_16  More advanced strata support
generate_alpha_scatterplot2 <- function (data.obj, alpha.obj=NULL, grp.name, subject=NULL,
		rarefy=TRUE, depth=NULL, iter.no=5, seed  = 123,  measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'), 
		strata=NULL, strata2=NULL, combine.strata=FALSE, combine.strata2=FALSE, smooth.method='loess', 
		pt.shape=16, pt.alpha=0.5, pt.size=2, subject.pt=FALSE,
		gg.cmd=NULL, ann='', hei=NULL, wid=NULL, pdf=TRUE) {	

	if (is.null(alpha.obj)) {
		alpha.obj <- generate_alpha_diversity(data.obj,  rarefy=rarefy, depth=depth, iter.no=iter.no, measures=measures, seed=seed)
		
	} else {
		if (sum(!(rownames(alpha.obj) %in% rownames(data.obj$meta.dat))) != 0){
			stop("alpha.obj contains samples not in data.obj!\n")
		} 
	}
	
	x <- alpha.obj
	df <- data.obj$meta.dat[rownames(alpha.obj), ]
	grp <- df[, grp.name]
	
	if(is.null(subject)) {
		ID <- NA
	} else {
		ID <- factor(df[, subject])
	}
	
	
	if (is.null(hei) | is.null(wid)) {
		hei <- 5
		wid <- 6
		
		if (!is.null(strata) & combine.strata == FALSE) {
			wid <- 4 * nlevels(factor(df[, strata]))
		} 
		if (!is.null(strata2) & combine.strata2 == FALSE) {
			hei <- 2.5 * nlevels(factor(df[, strata2]))
		} 
	}
	
	if (pdf == TRUE) {
		if (rarefy == T) {
			pdf(paste0('Alpha_diversity_scatterplot_rarefied_', ann, '.pdf'), height=hei, width=wid)
		} else {
			pdf(paste0('Alpha_diversity_scatterplot_unrarefied_', ann, '.pdf'), height=hei, width=wid)
		}
	}

	
	if (is.null(strata)) {
		for (measure in measures) {
			cat(measure, '\n')
			xx <- x[, measure]		
			df2 <- data.frame(Value=xx, Group=grp, ID=ID)
			obj <- ggplot(df2, aes(x=Group, y=Value, group=ID)) +
#					geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
					geom_smooth(method=smooth.method, aes(group=NA), size=0.75) +
					labs(y='Value', x=grp.name, title=measure) +
					theme(legend.position="none")

			if (!is.null(subject)) {
				if (subject.pt == FALSE) {
					obj <- obj + geom_line(size=0.25, alpha=0.75) 
				}  else {
					obj <- obj + geom_line(size=0.25, alpha=0.75) + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
				}
				
			} else {
				obj <- obj + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
			}
			if (!is.null(gg.cmd)) {
				obj <- obj + eval(parse(text=gg.cmd))
			}
			print(obj)
		}	
	} else {
		for (measure in measures) {
			cat(measure, '\n')
	
			xx <- x[, measure]
			grp2 <- factor(df[, strata])
			if (!is.null(strata2)) {
				grp3 <- factor(df[, strata2])
			} else {
				grp3 <- grp2
			}
			grp4 <- factor(paste(grp2, grp3))
			
			# VERY bad names with only difference in capital letters, variable names are also confusing
			# Rev: 2018_07_16 
			
			df2 <- data.frame(Value=xx, Group=grp, Strata=grp2, Strata2=grp3, Strata3=grp4, ID=ID)
				
			if (is.null(strata2)) {
				if (combine.strata == TRUE) {
					obj <- ggplot(df2, aes(x=Group, y=Value, group=ID, col=Strata)) +
#							geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
							geom_smooth(method=smooth.method, aes(group=Strata), size=0.75) +
							labs(y='Value', x=grp.name, title=measure) 
				} else {
					obj <- ggplot(df2, aes(x=Group, y=Value, group=ID)) +
#							geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
							geom_smooth(method=smooth.method, aes(group=NA), size=0.75) +
							labs(y='Value', x=grp.name, title=measure) + 
							facet_grid(. ~ Strata)
				}
				
			} else {
				if (combine.strata == TRUE) {
					if (combine.strata2 == TRUE) {
						obj <- ggplot(df2, aes(x=Group, y=Value, group=ID, col=Strata3)) +
#								geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
								geom_smooth(method=smooth.method, aes(group=Strata3), size=0.75) +
								labs(y='Value', x=grp.name, title=measure) 	
					} else {
						obj <- ggplot(df2, aes(x=Group, y=Value, group=ID, col=Strata)) +
#								geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
								geom_smooth(method=smooth.method, aes(group=Strata), size=0.75) +
								labs(y='Value', x=grp.name, title=measure) + 
								facet_grid(Strata2 ~ .)	
					}
				} else {
					if (combine.strata2 == TRUE) {
						warning('Currently, combine.strata2 will be forced to be FALSE if combine.strata=FALSE!\n')
					}
					obj <- ggplot(df2, aes(x=Group, y=Value, group=ID)) +
#							geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
							geom_smooth(method=smooth.method, aes(group=NA), size=0.75) +
							labs(y='Value', x=grp.name, title=measure) + 
							facet_grid(Strata2 ~ Strata)
				}
				
			}		
			
			obj <- obj + theme(legend.title=element_blank())	
		
			if (!is.null(subject)) {
				if (subject.pt == FALSE) {
					obj <- obj + geom_line(size=0.25, alpha=0.75) 
				}  else {
					obj <- obj + geom_line(size=0.25, alpha=0.75) + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
				}
				
			} else {
				obj <- obj + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
			}
			
			if (!is.null(gg.cmd)) {
				obj <- obj + eval(parse(text=gg.cmd))
			}
			print(obj)
		}
		
	}
	if(pdf == TRUE) {
		dev.off()
	}
	
	
}


# New: 2017_02_22 
# Depend on phyloseq - outdated. Please see the newer version generate_alpha_scatterplot2
generate_alpha_scatterplot <- function (data.obj, phylo.obj=NULL, rarefy=TRUE, depth=NULL, grp.name, 
		strata=NULL, strata2=NULL, wrap.nrow=2, smooth.method='auto', hei0=NULL, wid0=NULL,
		measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'), gg.cmd=NULL, ann='') {	
	
	# Rev: 2017_08_23
	if (is.null(phylo.obj)) {
		phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows=T), phy_tree(data.obj$tree), 
				tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))
	}
	
	if (rarefy == TRUE) {
		if (is.null(depth)) {
			depth <- min(sample_sums(phylo.obj))
		} else {
			if (depth > min(sample_sums(phylo.obj))) {
				ind <- (sample_sums(phylo.obj) >= depth)
				cat(sum(!ind), " samples do not have sufficient number of reads!\n")
				
				sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ] 
				data.obj <- subset_data(data.obj, ind)
			}
		}
		
		phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed=12345)
		x <- estimate_richness(phylo.even, measures=measures)
	} else {
		x <- estimate_richness(phylo.obj, measures=measures)
	}
	
	df <- data.obj$meta.dat
	grp <- df[, grp.name]

	if (is.null(hei0)) {
		hei <- 5
	} else {
		hei <- hei0
	}
	
	if (is.null(strata)) {
		if (is.null(wid0)) {
			wid <- 5
		} else {
			wid <- wid0
		}
		
	} else {
		if (is.null(wid0)) {
			wid <- 5.5
		} else {
			wid <- wid0
		}
	}
	
	if (rarefy == T) {
		pdf(paste0('Alpha_diversity_scatterplot_rarefied_', ann, '.pdf'), height=hei, width=wid)
	} else {
		pdf(paste0('Alpha_diversity_scatterplot_unrarefied_', ann, '.pdf'), height=hei, width=wid)
	}
	
	if (is.null(strata)) {
		for (measure in measures) {
			cat(measure, '\n')
			xx <- x[, measure]		
			df2 <- data.frame(Value=xx, Group=grp)
			obj <- ggplot(df2, aes(x=Group, y=Value)) +
					geom_point() +
					geom_smooth(method=smooth.method) +
					labs(y='Value', x=grp.name, title=measure) +
					theme(legend.position="none")
			if (!is.null(gg.cmd)) {
				obj <- obj + eval(parse(text=gg.cmd))
			}
			print(obj)
		}	
	} else {
		for (measure in measures) {
			cat(measure, '\n')
			
			xx <- x[, measure]		
			grp2 <- df[, strata]
			if (!is.null(strata2)) {
				grp3 <- df[, strata2]
			} else {
				grp3 <- grp2
			}
			
			df2 <- data.frame(Value=xx, Group=grp, Strata=grp2, Strata2=grp3)
			if (is.null(strata2)) {
				obj <- ggplot(df2, aes(x=Group, y=Value)) +
						geom_point() +
						geom_smooth(method=smooth.method) +
						labs(y='Value', x=grp.name, title=measure) + 
						facet_wrap(~ Strata, nrow=wrap.nrow)
			} else {
				obj <- ggplot(df2, aes(x=Group, y=Value)) +
						geom_point() +
						geom_smooth(method=smooth.method) +
						labs(y='Value', x=grp.name, title=measure) + 
						facet_grid(Strata2 ~ Strata)
			}		
			obj <- obj + theme(legend.title=element_blank())			
			if (!is.null(gg.cmd)) {
				obj <- obj + eval(parse(text=gg.cmd))
			}
			print(obj)
		}
		
	}
	dev.off()
	
}

#perform_alpha_test_otu <- function (data.obj, formula, grp.name, ...) {
#	# Assume two groups
#	df <- data.obj$meta.dat
#	grp <- df[, grp.name]
#	dep <- colSums(data.obj$otu.tab)
#	obs <- colSums(data.obj$otu.tab !=0)
#	obj <- lm(as.formula(paste('log(obs)', formula, '+ log(dep)')), df, ...)
#	sink('Alpha_diversity_test_results_observed_OTU.txt')
#	date()
#	cat('Linear model correcting for differential sequencing depth:\n')
#	prmatrix(summary(obj)$coefficients)
#	sink()
#	
#	pdf("Alpha_diversity_test_results_observed_OTU.pdf")
#	plot(dep, obs, log='xy', xlab="Sequencing Depth", ylab="Observed OTU Number", 
#			bg=c('blue', 'red')[grp], pch=21)
#	grp.levels <- levels(grp)
#	abline(lm(obs[grp==grp.levels[1]] ~ dep[grp==grp.levels[1]]), col='red')
#	abline(lm(obs[grp==grp.levels[2]] ~ dep[grp==grp.levels[2]]), col='blue')
#	legend('left', levels(grp), fill=c('blue', 'red'))
#	dev.off()
#}


# Rev: 2017_02_13 Individual label
# Rev: 2017_05_06 Add p value 
generate_ordination <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		grp.name, adj.name=NULL, emp.lev=NULL, indiv.lab=FALSE, indiv.lab.cex=0.5, is.lightcolor=TRUE,
		strata=NULL,  pc.separate, pca.method='cmd', ann=NULL, sub=NULL,
		clab=1.0, cex.pt=1.25, ellipse=T, cstar= 1, wid=5, hei=5, pdf=TRUE, is.pvalue=FALSE, fast.pvalue = FALSE,...) {
	# Implment strata
	# To be completed, add continuous case
	strata0 <- strata
	
	if (pdf) {
		if (is.null(ann)) {
			pdf(paste0('Beta_diversity_ordination_', pca.method, '_', grp.name, '.pdf'), width=wid, height=hei)
		} else {
			pdf(paste0('Beta_diversity_ordination_', pca.method, '_', ann, '.pdf'), width=wid, height=hei)
		}
	}

	df <- data.obj$meta.dat
	grp <- factor(df[, grp.name])
	
	if (!is.null(emp.lev)) {
		grp <- factor(grp, levels=c(setdiff(levels(grp), emp.lev), emp.lev))
	}
	
	if (!is.null(strata)) {
		strata <- factor(df[, strata])
	} else {
		strata <- factor(grp)
	}
	
	# To be revised
	darkcols <- hue_pal(l=40)(nlevels(grp))
	if (is.lightcolor) {
		lightcols <- hue_pal(c=45, l=80)(nlevels(grp))
	} else {
		lightcols <- darkcols
	}
	
	pchs <- rep(c(21, 22, 23, 24, 25), ceiling(nlevels(strata) / 5))[1:nlevels(strata)]
	
	for (dist.name in dist.names) {
		dist.temp <- dist.obj[[dist.name]]
		if (!is.null(adj.name)) {
			adj <- as.data.frame(df[, adj.name])
			obj <- cmdscale(as.dist(dist.temp), k=ncol(dist.temp)-1)
			dat2 <- apply(obj, 2, function(x) resid(lm(x ~ ., data=adj)))
			dist.temp <- dist(dat2)
		} 
		if (pca.method == 'cmd') {
			obj <- cmdscale(as.dist(dist.temp), k=2, eig=T)
			pve <- round(obj$eig[1:2]/sum(abs(obj$eig))*100, 1)
			y <- cbind(obj$points[, 1], obj$points[, 2])
			
			xlab <- paste0('PC1(', pve[1], '%)')
			ylab <- paste0('PC2(', pve[2], '%)')
		} 

		if (pca.method == 'nmds') {
				obj <- metaMDS(as.dist(dist.temp), k=2)
				y <- cbind(obj$points[, 1], obj$points[, 2])
				xlab <- 'NMDS1'
                ylab <- 'NMDS2'
		} 
		
		if (pca.method == 'pls') {
			require(mixOmics)
			# Test
			obj <- cmdscale(as.dist(dist.temp), k=ncol(dist.temp)-1)
			obj <- plsda(obj, grp, ncomp=2)
			y <- cbind(obj$variates$X[, 1], obj$variates$X[, 2])
			xlab <- 'PLS1'
			ylab <- 'PLS2'
		}
		
		if (is.pvalue == TRUE & nlevels(factor(grp)) > 1) {
		  if (fast.pvalue == TRUE) {
		    pvalue <- paste0('(P=', dmanova(as.dist(dist.temp) ~ grp)$aov.tab[1, 6], ')')
		  } else {
		    pvalue <- paste0('(P=', adonis(as.dist(dist.temp) ~ grp)$aov.tab[1, 6], ')')
		  }
			
		} else {
			pvalue <- ''
		}
		
		colnames(y) <- c("PC1", "PC2")		
		xlim <- c(-max(range(abs(y[, 1]))) * 1.25, max(range(abs(y[, 1]))) * 1.25)
		ylim <- c(-max(range(abs(y[, 2]))) * 1.25, max(range(abs(y[, 2]))) * 1.25)
		plot(y[, 1], y[, 2], type='n', xlim=xlim, ylim=ylim, 
				xlab=xlab, ylab=ylab)
#		points(y[, 1], y[, 2], bg='yellow', col=darkcols[grp], 
#				type='p', pch = pchs[grp], cex=cex.pt)
		col1 <- lightcols[1:nlevels(grp)]
		col2 <- darkcols[1:nlevels(grp)]
		col3 <- darkcols[grp]
		if (!is.null(emp.lev)) {
			
			col1 <- col2
			col1[which(levels(grp) != emp.lev)] <- rgb(0.5, 0.5, 0.5, 0.5)
			col1[which(levels(grp) == emp.lev)] <- rgb(1.0, 0, 0, 1.0)
			col2[which(levels(grp) != emp.lev)] <- rgb(0.5, 0.5, 0.5, 0.5)
			col2[which(levels(grp) == emp.lev)] <- rgb(1.0, 0, 0, 1.0)
            col3[grp !=  emp.lev] <- rgb(0.5, 0.5, 0.5, 0.5)
			col3[grp ==  emp.lev] <- rgb(1.0, 0, 0, 1.0)

		}
		if (ellipse == T) {
			s.class(y, 
					fac = grp,
					cstar = cstar,
					clab = 0,
					cpoint = 0,
					axesell = F,
					col = col1,
					grid = TRUE,
					add.plot=T,
					...
			)
		}
		
		points(y[, 1], y[, 2], bg=col3, col='black',
				type='p', pch = pchs[strata], cex=cex.pt)
		if (indiv.lab) {
			if (is.null(emp.lev)) {
				lab.temp <- rownames(dist.temp)
				text(y[, 1], y[, 2], labels=lab.temp, cex=indiv.lab.cex, col=rgb(1, 0, 0, 0.75))
			} else {
				lab.temp <- rownames(dist.temp)
				lab.temp[grp !=  emp.lev] <- ''
				text(y[, 1], y[, 2], labels=lab.temp, cex=indiv.lab.cex, col=rgb(1, 0, 0, 0.75))
			}
	
		}
		s.class(y, 
				fac = grp,
				cstar =0,
				cellipse = 0,
				clab = clab,
				cpoint = 0,
				axesell = F,
				col = col2,
				grid = TRUE,
				add.plot=T
		)
		if (!is.null(strata0)) {
			legend('topright', legend=(levels(strata)), pch=pchs[1:nlevels(strata)])	
		}

		title(main=paste(dist.name, "distance"), sub=paste0(sub, pvalue))
#		text(-0.25, 0.3, "PERMANOVA p=0.016")
	}
	
	if (pdf) {
		dev.off()
	}
}


# New: 2018_07_16  PCs trend
generate_PC_scatterplot <- function (data.obj, dist.obj=NULL, grp.name, adj.name=NULL, subject=NULL,
		dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), nPC=2, pca.method='cmd',
		strata=NULL, strata2=NULL, combine.strata=FALSE, combine.strata2=FALSE, smooth.method='loess', 
		pt.shape=16, pt.alpha=0.5, pt.size=2, subject.pt = FALSE,
		ann='', hei=NULL, wid=NULL, pdf=TRUE, gg.cmd=NULL) {	
	
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	if(is.null(subject)) {
		ID <- NA
	} else {
		ID <- factor(df[, subject])
	}
	
	if (is.null(hei) | is.null(wid)) {
		hei <- 5
		wid <- 6
		
		if (!is.null(strata) & combine.strata == FALSE) {
			wid <- 4 * nlevels(factor(df[, strata]))
		} 
		if (!is.null(strata2) & combine.strata2 == FALSE) {
			hei <- 2.5 * nlevels(factor(df[, strata2]))
		} 
	}
	
	if (pdf == TRUE) {
		pdf(paste0('Beta_diversity_PC_scatterplot_', ann, '.pdf'), height=hei, width=wid)
	} 
	
	measures <- NULL
	x <- NULL
	for (dist.name in dist.names) {
		dist.temp <- dist.obj[[dist.name]]
		if (!is.null(adj.name)) {
			adj <- as.data.frame(df[, adj.name])
			obj <- cmdscale(as.dist(dist.temp), k=ncol(dist.temp)-1)
			dat2 <- apply(obj, 2, function(x) resid(lm(x ~ ., data=adj)))
			dist.temp <- dist(dat2)
		} 
		if (pca.method == 'cmd') {
			obj <- cmdscale(as.dist(dist.temp), k=nPC, eig=T)
			pve <- round(obj$eig[1:nPC]/sum(abs(obj$eig))*100, 1)
			x <- cbind(x, obj$points)
			measures <- c(measures, paste0(dist.name, ' - MDS PC', 1:nPC, '(', pve, '%)'))
		} 
		
		if (pca.method == 'nmds') {
			obj <- metaMDS(as.dist(dist.temp), k=nPC)
			x <- cbind(x, obj$points)
			measures <- c(measures, paste0(dist.name, ' - NMDS PC', 1:nPC))
		} 
	}
	colnames(x) <- measures
		

	if (is.null(strata)) {
		for (measure in measures) {
			cat(measure, '\n')
			xx <- x[, measure]		
			df2 <- data.frame(Value=xx, Group=grp, ID=ID)
			obj <- ggplot(df2, aes(x=Group, y=Value, group=ID)) +
					geom_smooth(method=smooth.method, aes(group=NA), size=0.75) +
					labs(y='Value', x=grp.name, title=measure) +
					theme(legend.position="none")
			if (!is.null(gg.cmd)) {
				obj <- obj + eval(parse(text=gg.cmd))
			}
			
			if (!is.null(subject)) {
				if (subject.pt == FALSE) {
					obj <- obj + geom_line(size=0.25, alpha=0.75) 
				}  else {
					obj <- obj + geom_line(size=0.25, alpha=0.75) + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
				}
				
			} else {
				obj <- obj + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
			}
			print(obj)
		}	
	} else {
		for (measure in measures) {
			cat(measure, '\n')
			
			xx <- x[, measure]
			grp2 <- factor(df[, strata])
			if (!is.null(strata2)) {
				grp3 <- factor(df[, strata2])
			} else {
				grp3 <- grp2
			}
			grp4 <- factor(paste(grp2, grp3))
			
			# VERY bad names with only difference in capital letters, variable names are also confusing
			# Rev: 2018_07_16 
			
			df2 <- data.frame(Value=xx, Group=grp, Strata=grp2, Strata2=grp3, Strata3=grp4, ID=ID)
			
			if (is.null(strata2)) {
				if (combine.strata == TRUE) {
					obj <- ggplot(df2, aes(x=Group, y=Value, group=ID, col=Strata)) +
					#		geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
							geom_smooth(method=smooth.method, aes(group=Strata), size=0.75) +
							labs(y='Value', x=grp.name, title=measure) 
				} else {
					obj <- ggplot(df2, aes(x=Group, y=Value, group=ID)) +
					#		geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
							geom_smooth(method=smooth.method,  aes(group=NA),  size=0.75) +
							labs(y='Value', x=grp.name, title=measure) + 
							facet_grid(. ~ Strata)
				}
				
			} else {
				if (combine.strata == TRUE) {
					if (combine.strata2 == TRUE) {
						obj <- ggplot(df2, aes(x=Group, y=Value, group=ID, col=Strata3)) +
					#			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
								geom_smooth(method=smooth.method, aes(group=Strata3), size=0.75) +
								labs(y='Value', x=grp.name, title=measure) 	
					} else {
						obj <- ggplot(df2, aes(x=Group, y=Value, group=ID, col=Strata)) +
					#			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
								geom_smooth(method=smooth.method, aes(group=Strata), size=0.75) +
								labs(y='Value', x=grp.name, title=measure) + 
								facet_grid(Strata2 ~ .)	
					}
				} else {
					if (combine.strata2 == TRUE) {
						warning('Currently, combine.strata2 will be forced to be FALSE if combine.strata=FALSE!\n')
					}
					obj <- ggplot(df2, aes(x=Group, y=Value, group=ID)) +
						#	geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
							geom_smooth(method=smooth.method,  aes(group=NA), size=0.75) +
							labs(y='Value', x=grp.name, title=measure) + 
							facet_grid(Strata2 ~ Strata)
				}
				
			}		
			
			obj <- obj + theme(legend.title=element_blank())	
			
			if (!is.null(gg.cmd)) {
				obj <- obj + eval(parse(text=gg.cmd))
			}
			
			if (!is.null(subject)) {
				if (subject.pt == FALSE) {
					obj <- obj + geom_line(size=0.25, alpha=0.75) 
				}  else {
					obj <- obj + geom_line(size=0.25, alpha=0.75) + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
				}
			} else {
				obj <- obj + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
			}
			
			print(obj)
		}
		
	}
	if (pdf == TRUE) {
		dev.off()
	}
	
	
}


# New: 2018_07_16  grp.name codes the day
generate_ordination_trajectory <- function (data.obj, dist.obj=NULL, day.name=NULL, adj.name=NULL, subject=NULL,
		dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),  pca.method='cmd',
		strata=NULL, strata2=NULL, combine.strata=FALSE, combine.strata2=FALSE, smooth.method='auto', 
		pt.shape=16, pt.alpha=0.75, pt.size=2, subject.pt = TRUE,
		ann='', hei=NULL, wid=NULL, pdf=TRUE, gg.cmd=NULL) {	
	
	nPC <- 2
	df <- data.obj$meta.dat
	ind <- order(df[, day.name])
	
	df <- df[ind, ]
	dist.obj <- lapply(dist.obj, function (x) x[ind, ind])
	
	day <- df[, day.name]
	
	
	if(is.null(subject)) {
		ID <- NA
	} else {
		ID <- factor(df[, subject])
	}
	
	if (is.null(hei) | is.null(wid)) {
		hei <- 5
		wid <- 6
		
		if (!is.null(strata) & combine.strata == FALSE) {
			wid <- 4 * nlevels(factor(df[, strata]))
		} 
		if (!is.null(strata2) & combine.strata2 == FALSE) {
			hei <- 2.5 * nlevels(factor(df[, strata2]))
		} 
	}
	
	if (pdf == TRUE) {
		pdf(paste0('Beta_diversity_Subject_Trajectory_', ann, '.pdf'), height=hei, width=wid)
	} 
	
	x1.measures <- NULL
	x2.measures <- NULL
	x1 <- NULL
	x2 <- NULL
	for (dist.name in dist.names) {
		dist.temp <- dist.obj[[dist.name]]
		if (!is.null(adj.name)) {
			adj <- as.data.frame(df[, adj.name])
			obj <- cmdscale(as.dist(dist.temp), k=ncol(dist.temp)-1)
			dat2 <- apply(obj, 2, function(x) resid(lm(x ~ ., data=adj)))
			dist.temp <- dist(dat2)
		} 
		if (pca.method == 'cmd') {
			obj <- cmdscale(as.dist(dist.temp), k=nPC, eig=T)
			pve <- round(obj$eig[1:nPC]/sum(abs(obj$eig))*100, 1)
			x1 <- cbind(x1, obj$points[, 1])
			x2 <- cbind(x2, obj$points[, 2])
			x1.measures <- c(x1.measures, paste0('MDS PC', 1, '(', pve[1], '%)'))
			x2.measures <- c(x2.measures, paste0('MDS PC', 2, '(', pve[2], '%)'))
		} 
		
		if (pca.method == 'nmds') {
			obj <- metaMDS(as.dist(dist.temp), k=nPC)
			x1 <- cbind(x1, obj$points[, 1])
			x2 <- cbind(x2, obj$points[, 2])
			x1.measures <- c(x1.measures, paste0('MDS PC', 1))
			x2.measures <- c(x2.measures, paste0('MDS PC', 2))
		} 
	}
	colnames(x1) <- dist.names
	colnames(x2) <- dist.names
	names(x1.measures) <- dist.names
	names(x2.measures) <- dist.names
	measures <- dist.names
	
	
	if (is.null(strata)) {
		for (measure in measures) {
			cat(measure, '\n')
			xx1 <- x1[, measure]		
			xx2 <- x2[, measure]
			df2 <- data.frame(Day=factor(day), ID=ID, xx1=xx1, xx2=xx2)
			obj <- ggplot(df2, aes(x=xx1, y=xx2, group=ID, col=Day)) +
					#		geom_smooth(method=smooth.method, aes(group=NA), size=0.75) +
					labs(y=x2.measures[measure], x=x1.measures[measure], title=measure) 

			
			if (subject.pt == FALSE) {
				obj <- obj + geom_path(size=0.25, alpha=0.5, colour='black') 
			}  else {
				obj <- obj + geom_path(size=0.25, alpha=0.5, colour='black') + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
			}
			
			if (!is.null(gg.cmd)) {
				obj <- obj + eval(parse(text=gg.cmd))
			}
			
			obj <- obj + theme(legend.title=element_blank())	
			
			print(obj)
		}	
	} else {
		for (measure in measures) {
			cat(measure, '\n')
			
			xx1 <- x1[, measure]		
			xx2 <- x2[, measure]
			grp2 <- factor(df[, strata])
			if (!is.null(strata2)) {
				grp3 <- factor(df[, strata2])
			} else {
				grp3 <- grp2
			}
			grp4 <- factor(paste(grp2, grp3))
			
			# VERY bad names with only difference in capital letters, variable names are also confusing
			# Rev: 2018_07_16 
			
			df2 <- data.frame(xx1=xx1, xx2=xx2, Strata=grp2, Strata2=grp3, Strata3=grp4, ID=ID, Day=factor(day))
			
			if (is.null(strata2)) {
				if (combine.strata == TRUE) {
					obj <- ggplot(df2, aes(x=xx1, y=xx2, group=ID, col=Day, linetype=Strata)) +
							#		geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
							#	geom_smooth(method=smooth.method, aes(group=Strata), size=0.75) +
							labs(y=x2.measures[measure], x=x1.measures[measure], title=measure)
					
					if (subject.pt == FALSE) {
						obj <- obj + geom_path(size=0.25, alpha=0.75, aes(col=Strata)) 
					}  else {
						obj <- obj + geom_path(size=0.25, alpha=0.75, aes(col=Strata)) + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
					}
					
				} else {
					obj <- ggplot(df2, aes(x=xx1, y=xx2, group=ID, col=Day)) +
							#		geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
							# geom_smooth(method=smooth.method,  aes(group=NA),  size=0.75) +
							labs(y=x2.measures[measure], x=x1.measures[measure], title=measure) +
							facet_grid(. ~ Strata)
					if (subject.pt == FALSE) {
						obj <- obj + geom_path(size=0.25, alpha=0.5, col='black') 
					}  else {
						obj <- obj + geom_path(size=0.25, alpha=0.5, col='black') + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
					}
					
				}
				
			} else {
				if (combine.strata == TRUE) {
					if (combine.strata2 == TRUE) {
						obj <- ggplot(df2, aes(x=xx1, y=xx2, group=ID, col=Day, linetype=Strata3)) +
								#			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
								# geom_smooth(method=smooth.method, aes(group=Strata3), size=0.75) +
								labs(y=x2.measures[measure], x=x1.measures[measure], title=measure) 
						if (subject.pt == FALSE) {
							obj <- obj + geom_path(size=0.25, alpha=0.75, aes(col=Strata3)) 
						}  else {
							obj <- obj + geom_path(size=0.25, alpha=0.75, aes(col=Strata3)) + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
						}
					} else {
						obj <- ggplot(df2, aes(x=xx1, y=xx2, group=ID, col=Day, linetype=Strata)) +
								#			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
								# geom_smooth(method=smooth.method, aes(group=Strata), size=0.75) +
								labs(y=x2.measures[measure], x=x1.measures[measure], title=measure) +
								facet_grid(Strata2 ~ .)	
						if (subject.pt == FALSE) {
							obj <- obj + geom_path(size=0.25, alpha=0.75, aes(col=Strata)) 
						}  else {
							obj <- obj + geom_path(size=0.25, alpha=0.75, aes(col=Strata)) + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
						}
			
					}
				} else {
					if (combine.strata2 == TRUE) {
						warning('Currently, combine.strata2 will be forced to be FALSE if combine.strata=FALSE!\n')
					}
					obj <- ggplot(df2, aes(x=xx1, y=xx2, group=ID, col=Day)) +
							#	geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
							#	geom_smooth(method=smooth.method,  aes(group=NA), size=0.75) +
							labs(y=x2.measures[measure], x=x1.measures[measure], title=measure) +
							facet_grid(Strata2 ~ Strata)
					if (subject.pt == FALSE) {
						obj <- obj + geom_path(size=0.25, alpha=0.5, colour='black') 
					}  else {
						obj <- obj + geom_path(size=0.25, alpha=0.5, colour='black') + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
					}
				}
				
			}		
			
			obj <- obj + theme(legend.title=element_blank())	
			
			if (!is.null(gg.cmd)) {
				obj <- obj + eval(parse(text=gg.cmd))
			}
			
			print(obj)
		}
		
	}
	if (pdf == TRUE) {
		dev.off()
	}
	
	
}

# New: 2017_02_21  Separate PC plots for strata2
# Rev: 2017_05_04 Add p value
# Rev: 2017_06_01 Different pages for different distances
generate_ordination_separate <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		grp.name, adj.name=NULL, emp.lev=NULL, indiv.lab=FALSE, indiv.lab.cex=0.5, is.lightcolor=TRUE,
		strata=NULL,  separate=NULL,  layout=NULL, pc.separate, pca.method='cmd', ann=NULL, sub=NULL,
		clab=1.0, cex.pt=1.25, ellipse=T, cstar= 1, wid=NULL, hei=NULL, pdf=TRUE, is.pvalue = FALSE, fast.pvalue = FALSE,...) {
	# Implment strata
	# To be completed, add continuous case
	strata0 <- strata
	if (is.null(separate)){
		stop('The augment separate is required!\n')
	} else {
		data.obj$meta.dat[, separate] <- factor(data.obj$meta.dat[, separate])
		sep.levels <- levels(data.obj$meta.dat[, separate])
	}
	
	if (!is.null(layout)) {
		if (is.null(hei) | is.null(wid)) {
			hei <- 3 * layout[1]
			wid <- 3 * layout[2]
		}

	} else {
		layout <- c(1, length(sep.levels))
		if (is.null(hei) | is.null(wid)) {
			hei <- 3
			wid <- 3 * length(sep.levels)
		}
	}
	
	if (pdf) {

		if (is.null(ann)) {
			pdf(paste0('Beta_diversity_ordination_', pca.method, '_', grp.name, '_separate.pdf'), width=wid, height=hei)
		} else {
			pdf(paste0('Beta_diversity_ordination_', pca.method, '_', ann, '_separate.pdf'), width=wid, height=hei)
		}
	}
	
	par(mfrow=layout)

	df <- data.obj$meta.dat
	grp <- factor(df[, grp.name])
	
	if (!is.null(emp.lev)) {
		grp <- factor(grp, levels=c(setdiff(levels(grp), emp.lev), emp.lev))
	}
	
	if (!is.null(strata)) {
		strata <- factor(df[, strata])
	} else {
		strata <- factor(grp)
	}
	
	# To be revised
	darkcols <- hue_pal(l=40)(nlevels(grp))
	if (is.lightcolor) {
		lightcols <- hue_pal(c=45, l=80)(nlevels(grp))
	} else {
		lightcols <- darkcols
	}
	
	pchs <- rep(c(21, 22, 23, 24, 25), ceiling(nlevels(strata) / 5))[1:nlevels(strata)]
	
	for (dist.name in dist.names) {
		dist.temp <- dist.obj[[dist.name]]
		
		
		for (sep.level in sep.levels) {
			ind <- df[, separate] == sep.level
			df2 <- df[ind, ]
			grp2 <- grp[ind]
			strata2 <- strata[ind]
			dist.temp2 <- dist.temp[ind, ind]
			
			if (!is.null(adj.name)) {
				adj <- as.data.frame(df2[, adj.name])
				obj <- cmdscale(as.dist(dist.temp2), k=ncol(dist.temp2)-1)
				dat2 <- apply(obj, 2, function(x) resid(lm(x ~ ., data=adj)))
				dist.temp2 <- dist(dat2)
			} 
			if (pca.method == 'cmd') {
				obj <- cmdscale(as.dist(dist.temp2), k=2, eig=T)
				pve <- round(obj$eig[1:2]/sum(abs(obj$eig))*100, 1)
				y <- cbind(obj$points[, 1], obj$points[, 2])
				
				xlab <- paste0('PC1(', pve[1], '%)')
				ylab <- paste0('PC2(', pve[2], '%)')
			} 
			
			if (pca.method == 'nmds') {
				obj <- metaMDS(as.dist(dist.temp2), k=2)
				y <- cbind(obj$points[, 1], obj$points[, 2])
				xlab <- 'NMDS1'
				ylab <- 'NMDS2'
			} 
			
			if (pca.method == 'pls') {
				require(mixOmics)
				# Test
				obj <- cmdscale(as.dist(dist.temp2), k=ncol(dist.temp2)-1)
				obj <- plsda(obj, grp2, ncomp=2)
				y <- cbind(obj$variates$X[, 1], obj$variates$X[, 2])
				xlab <- 'PLS1'
				ylab <- 'PLS2'
			}
			
			if (is.pvalue == TRUE & nlevels(factor(grp2)) > 1) {
			  if (fast.pvalue == TRUE) {
			    pvalue <- paste0('(P=', dmanova(as.dist(dist.temp2) ~ grp2)$aov.tab[1, 6], ')')
			  } else {
			    pvalue <- paste0('(P=', adonis(as.dist(dist.temp2) ~ grp2)$aov.tab[1, 6], ')')
			  }

			} else {
				pvalue <- ''
			}
			colnames(y) <- c("PC1", "PC2")	
			xlim <- c(-max(range(abs(y[, 1]))) * 1.25, max(range(abs(y[, 1]))) * 1.25)
			ylim <- c(-max(range(abs(y[, 2]))) * 1.25, max(range(abs(y[, 2]))) * 1.25)
			plot(y[, 1], y[, 2], type='n', xlim=xlim, ylim=ylim, 
					xlab=xlab, ylab=ylab)
#		points(y[, 1], y[, 2], bg='yellow', col=darkcols[grp2], 
#				type='p', pch = pchs[grp2], cex=cex.pt)
			col1 <- lightcols[1:nlevels(grp2)]
			col2 <- darkcols[1:nlevels(grp2)]
			col3 <- darkcols[grp2]
			if (!is.null(emp.lev)) {
				
				col1 <- col2
				col1[which(levels(grp2) != emp.lev)] <- rgb(0.5, 0.5, 0.5, 0.5)
				col1[which(levels(grp2) == emp.lev)] <- rgb(1.0, 0, 0, 1.0)
				col2[which(levels(grp2) != emp.lev)] <- rgb(0.5, 0.5, 0.5, 0.5)
				col2[which(levels(grp2) == emp.lev)] <- rgb(1.0, 0, 0, 1.0)
				col3[grp2 !=  emp.lev] <- rgb(0.5, 0.5, 0.5, 0.5)
				col3[grp2 ==  emp.lev] <- rgb(1.0, 0, 0, 1.0)
				
			}
			if (ellipse == T) {
				s.class(y, 
						fac = grp2,
						cstar = cstar,
						clab = 0,
						cpoint = 0,
						axesell = F,
						col = col1,
						grid = TRUE,
						add.plot=T,
						...
				)
			}
			
			points(y[, 1], y[, 2], bg=col3, col='black',
					type='p', pch = pchs[strata2], cex=cex.pt)
			if (indiv.lab) {
				if (is.null(emp.lev)) {
					lab.temp <- rownames(dist.temp2)
					text(y[, 1], y[, 2], labels=lab.temp, cex=indiv.lab.cex, col=rgb(1, 0, 0, 0.75))
				} else {
					lab.temp <- rownames(dist.temp2)
					lab.temp[grp2 !=  emp.lev] <- ''
					text(y[, 1], y[, 2], labels=lab.temp, cex=indiv.lab.cex, col=rgb(1, 0, 0, 0.75))
				}
				
			}
			s.class(y, 
					fac = grp2,
					cstar =0,
					cellipse = 0,
					clab = clab,
					cpoint = 0,
					axesell = F,
					col = col2,
					grid = TRUE,
					add.plot=T
			)
			if (!is.null(strata0)) {
				legend('topright', legend=(levels(strata2)), pch=pchs[1:nlevels(strata2)])	
			}
			
			title(main=sep.level, sub=paste(dist.name, "distance", pvalue))

		}


#		text(-0.25, 0.3, "PERMANOVA p=0.016")
	}
	
	if (pdf) {
		dev.off()
	}
}


# Rev:2017_04_19 Add bootstrap standard error

generate_distance_barplot <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
		grp.name, strata=NULL, within=T, between=T, bt.no = 100, ann='') {
	strata.name <- strata
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	grp.levels <- levels(grp)
	grp.nlevels <- nlevels(grp)
	grp.btws <- outer(grp.levels, grp.levels, paste, sep="#")
	grp.btws <- grp.btws[lower.tri(grp.btws)]
	
	if (!is.null(strata.name)) {
		strata <- df[, strata.name]
	} else {
		strata <- factor(rep(1, nrow(df))) # pseudo strata
	}
	res.df <- NULL
	for (dist.name in dist.names) {
		for (stratum in levels(strata)) {
			ind <- strata %in% stratum
			dist.sub <- dist.obj[[dist.name]][ind, ind]
			df2 <- df[ind, , drop=FALSE]
			if (between) {
				for (grp.btw in grp.btws) {
					ind1 <- which(df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[1])
					ind2 <- which(df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[2])	
					temp <- as.vector(dist.sub[ind1, ind2])
					sem <- sd(sapply(1:bt.no, function (i) {
								mean(dist.sub[as.numeric(sample(paste(ind1), repl=TRUE)), 
												as.numeric(sample(paste(ind2), repl=TRUE))])
							}))
					res.df <- rbind(res.df, data.frame(DistanceMetric=dist.name, Strata=stratum, Compare='Between', 
									DistanceType=grp.btw, Distance=mean(temp), sd=sem))
				}	
			}

			if (within) {
				for (grp.wth in grp.levels) {
					ind1 <- which(df2[, grp.name] == grp.wth)
					temp <- dist.sub[ind1, ind1]			
					temp <- temp[lower.tri(temp)]
					sem <- sd(sapply(1:bt.no, function (i) {
										ind2 <- as.numeric(sample(paste(ind1), repl = TRUE))
										temp <- dist.sub[ind2, ind2]			
										temp <- temp[lower.tri(temp)]
										mean(temp)
									}))
					res.df <- rbind(res.df, data.frame(DistanceMetric=dist.name, Strata=stratum, Compare='Within',
									DistanceType=grp.wth, Distance=mean(temp), sd=sem))
				}	
			}
		}
		
	}
	if (between & within) {
		res.df$DistanceType <- factor(res.df$DistanceType, levels=c(grp.btws, grp.levels))
		levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'), paste('Within\n', grp.levels))
	} else {
		if (between) {
			res.df$DistanceType <- factor(res.df$DistanceType)
			levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'))
		} else {
			res.df$DistanceType <- factor(res.df$DistanceType)
			levels(res.df$DistanceType) <- c(paste('Within\n', grp.levels))
		}
	}
	
	if (is.null(strata.name)) {
		pdf(paste0("Beta_diversity_btw", between, "_wth", within, "no_strata_barplot_", ann, ".pdf"),
				width=5, height=5)
		
		limits <- aes(ymax = Distance + sd, ymin=Distance - sd)
		dodge <- position_dodge(width=0.9)
		for (dist.name in dist.names) {
			cat(dist.name, "...\n")
			temp <- res.df[res.df$DistanceMetric == dist.name, ]
			obj1 <- ggplot(temp, aes(x=DistanceType, y=Distance, fill=DistanceType)) + 
					geom_bar(position=dodge, stat="identity", width=0.75) + 
					geom_bar(position=dodge, stat="identity", width=0.75, colour="black", size=0.25) +
					geom_errorbar(limits, position=dodge, size=0.25, width=0.25) +
					labs(y=paste(dist.name, "Distance"), x='') +
					theme(legend.position="none") +
					theme(axis.text.x=element_text(angle=90, hjust=1))
					
			print(obj1)
			
		}
		dev.off()
	} else {
		pdf(paste0("Beta_diversity_btw", between, "_wth", within, "_strata", strata.name, "barplot_", ann, ".pdf"),
				width=2.5*(nlevels(strata)-1) + 5, height=5)
		
		limits <- aes(ymax = Distance + sd, ymin=Distance - sd)
		dodge <- position_dodge(width=0.9)
		for (dist.name in dist.names) {
			cat(dist.name, "...\n")
			temp <- res.df[res.df$DistanceMetric == dist.name, ]
			obj1 <- ggplot(temp, aes(x=Strata, y=Distance, fill=DistanceType)) + 
					geom_bar(position=dodge, stat="identity") + 
					geom_bar(position=dodge, stat="identity", colour="black", size=0.25) +
					geom_errorbar(limits, position=dodge, size=0.25, width=0.5) +
					labs(y=paste(dist.name, "Distance"), x=strata.name) +
					theme(axis.text.x=element_text(angle=90, hjust=1))
			print(obj1)
			
		}
		dev.off()
		
	}
	
}

generate_distance_boxplot <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
		grp.name, strata=NULL, within=F, between=T, ann='') {
	strata.name <- strata
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	grp.levels <- levels(grp)
	grp.nlevels <- nlevels(grp)
	grp.btws <- outer(grp.levels, grp.levels, paste, sep="#")
	grp.btws <- grp.btws[lower.tri(grp.btws)]
	
	if (!is.null(strata.name)) {
		strata <- df[, strata.name]
	} else {
		strata <- factor(rep(1, nrow(df))) # pseudo strata
	}
	res.df <- NULL
	for (dist.name in dist.names) {
		for (stratum in levels(strata)) {
			ind <- strata %in% stratum
			dist.sub <- dist.obj[[dist.name]][ind, ind]
			df2 <- df[ind, , drop=FALSE]
			if (between) {
				for (grp.btw in grp.btws) {
					ind1 <- df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[1]
					ind2 <- df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[2]	
					temp <- as.vector(dist.sub[ind1, ind2])
					n <- length(temp)
					res.df <- rbind(res.df, data.frame(DistanceMetric=rep(dist.name, n), Strata=rep(stratum, n), Compare=rep('Between', n),
									DistanceType=rep(grp.btw, n), Distance=temp))
				}	
			}
			
			if (within) {
				for (grp.wth in grp.levels) {
					ind1 <- df2[, grp.name] == grp.wth
					temp <- dist.sub[ind1, ind1]			
					temp <- temp[lower.tri(temp)]
					n <- length(temp)
					res.df <- rbind(res.df, data.frame(DistanceMetric=rep(dist.name, n), Strata=rep(stratum, n), Compare=rep('Within', n),
									DistanceType=rep(grp.wth, n), Distance=temp))
				}	
			}
		}
		
	}
	if (between & within) {
		res.df$DistanceType <- factor(res.df$DistanceType, levels=c(grp.btws, grp.levels))
		levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'), paste('Within\n', grp.levels))
	} else {
		if (between) {
			res.df$DistanceType <- factor(res.df$DistanceType)
			levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'))
		} else {
			res.df$DistanceType <- factor(res.df$DistanceType)
			levels(res.df$DistanceType) <- c(paste('Within\n', grp.levels))
		}
	}
	
	if (is.null(strata.name)) {
		pdf(paste0("Beta_diversity_btw", between, "_wth", within, "no_strata_boxplot_", ann, ".pdf"),
				width=5, height=5)
		
		for (dist.name in dist.names) {
			cat(dist.name, "...\n")
			temp <- res.df[res.df$DistanceMetric == dist.name, ]
			
			dodge <- position_dodge(width=0.95)		
			obj1 <- ggplot(temp, aes(x=DistanceType, y=Distance, col=DistanceType)) + 
					geom_boxplot(position=dodge, outlier.colour = NA) + 
				# geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) +
					labs(y=paste(dist.name, "Distance"), x='') +
			    theme_bw() + 
					theme(legend.position="none", panel.grid= element_blank())
			
			print(obj1)
			
		}
		dev.off()
	} else {
		pdf(paste0("Beta_diversity_btw", between, "_wth", within, "_strata", strata.name, "boxplot_", ann, ".pdf"),
				width=2.5*(nlevels(strata)-1) + 5, height=5)
		
		for (dist.name in dist.names) {
			cat(dist.name, "...\n")
			temp <- res.df[res.df$DistanceMetric == dist.name, ]
			dodge <- position_dodge(width=0.95)		
			obj1 <- ggplot(temp, aes(x=Strata, y=Distance, col=DistanceType)) + 
					geom_boxplot(position=dodge, outlier.colour = NA) + 
		#			geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) +
					labs(y=paste(dist.name, "Distance"), x=strata.name)
			print(obj1)
			
		}
		dev.off()
		
	}
	
}


# Rev:2016_11_25
generate_clustering <- function(data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), meta.info, cluster.method='ward.D', 
		is.labRow=F, cex.lab=NULL, ann="", wid=10, hei=6, colFnsF=NULL, colFnsC=NULL) {
	if (is.null(colFnsC)) {
		colFnsC <- c(colorRampPalette(c('black', 'green')), colorRampPalette(c('black', 'blue')), colorRampPalette(c('black', 'red')))
	} 
	if (is.null(colFnsF)) {
		rainbow3 <- function(x) {
			rainbow(x + 2)[1:x]
		}
		jet3 <- function(x) {
			jet(x + 2)[1:x]
		}
		colFnsF <- c(rainbow3, jet3)
	}
	
	pdf(paste0("Beta_diversity_Hierachical_clustering_", ann, ".pdf"), width=wid, height=hei)
	for (dist.name in dist.names) {

		df <- data.obj$meta.dat
		dist.sub <- dist.obj[[dist.name]]
		dend <- hclust(as.dist(dist.sub), cluster.method)
		

		key.list <- list()
		mat <- NULL
		for (keyID in meta.info) {
			x <- df[, keyID]
			i <- 0
			j <- 0
			if (is.factor(x)) {
				key.list[[keyID]] <- list(breaks=levels(x), colors=(colFnsF[i+1][[1]])(nlevels(x)), base=NA, col.na=NA, right=F, include.lowest=F)	
				i <- (i + 1) %% length(colFnsF)
				mat <- cbind(mat, key.list[[keyID]]$colors[x])
			} else {
				key.list[[keyID]] <- makecmap(x, n=5, colFn=colFnsC[j+1][[1]])
				j <- (j + 1) %% length(colFnsC)
				mat <- cbind(mat, cmap(x, key.list[[keyID]]))
			}
		}
		colnames(mat) <- meta.info

		par(oma = c(1, 2, 1, 2))
		par(mar = c(5,4,4,3)+0.1)  # make space for color keys
		if (is.labRow == F) {
			labRow <- ""
		} else {
			labRow <- rownames(df)
		}
		if (is.null(cex.lab)) {
			cex.lab <- 25 / nrow(df) 
		}
		
		dendromat(dend, mat, labRow=labRow,
				ylab = 'Distance', main = paste(dist.name, "distance"), cex.lab=cex.lab)
		
		par(oma=c(0, 0, 0, 0))
		y.cord <- (1/length(meta.info)) * (0:(length(meta.info) - 1))
		k <- 1
		for (keyID in meta.info) {
			x <- df[, keyID]
			if (is.factor(x)) {
				vkey2(key.list[[keyID]], keyID, y=y.cord[k], stretch=1.2)
			} else {
				vkey(key.list[[keyID]], keyID, y=y.cord[k], stretch=1.2)
			}
			k <- k + 1
		}		
	}
	dev.off()
}


# Rev: 2016_09_10, Implement p value based omnibus test
PermanovaG2 <- function (formula, dat = NULL, ...) 
{
	save.seed <- get(".Random.seed", .GlobalEnv)
	lhs <- formula[[2]]
	lhs <- eval(lhs, dat, parent.frame())
	rhs <- as.character(formula)[3]
	p.perms <- list()
	p.obs <- list()
	for (i in 1:(dim(lhs)[3])) {
		assign(".Random.seed", save.seed, .GlobalEnv)
		Y <- as.dist(lhs[, , i])
		formula2 <- as.formula(paste("Y", "~", rhs))
		obj <- adonis(formula2, dat, ...)
		perm.mat <- obj$f.perms
		p.perms[[i]] <- 1 - (apply(perm.mat, 2, rank) - 1) / nrow(perm.mat)
		p.obs[[i]] <- obj$aov.tab[1:ncol(perm.mat), "Pr(>F)"]

	}
	
	omni.pv <- NULL
	indiv.pv <- NULL
	for (j in 1:ncol(perm.mat)) {
		p.perms.j <- sapply(p.perms, function (x) x[, j])
		p.obj.j <- sapply(p.obs, function (x) x[j])
		omni.pv <- c(omni.pv, mean(c(rowMins(p.perms.j ) <= min(p.obj.j), 1)))
		indiv.pv <- rbind(indiv.pv, p.obj.j)
	}
	colnames(indiv.pv) <- paste0('D', 1:ncol(indiv.pv), '.p.value')
	rownames(indiv.pv) <- 1:nrow(indiv.pv)
	
	aov.tab <- data.frame(indiv.pv, omni.p.value = omni.pv)
	rownames(aov.tab) <- rownames(obj$aov.tab)[1:ncol(perm.mat)]
	list(aov.tab = aov.tab)
}


perform_permanova_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		PermanovaG.dist=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
		formula=NULL,  grp.name=NULL, adj.name=NULL, pairwise=F, block.perm=F, strata=NULL, ann='', ...) {
	# PermanovaG not implemented for block permutation
	result <- list()
	
	df <- data.obj$meta.dat
	if (!is.null(strata)) {
		if (is.character(strata)) {
			strata <- factor(df[, strata])
		}
	}
	if (!is.null(formula)) {
		ind <- apply(df[, gsub("^\\s+|\\s+$", "", strsplit(strsplit(formula, "~")[[1]][2], "\\+")[[1]]), drop=F], 1, function(x) sum(is.na(x))) == 0
		df <- df[ind, ]
		if (!is.null(strata)) {
			strata <- factor(strata[ind])
		}
		sink(paste0('Beta_diversity_PERMANOVA_test_', ann, '.txt'))
		date()
		cat('\nPERMANOVA test: \n')
		permanova.obj <- list()
		for (dist.name in dist.names) {
			cat(dist.name, " distance: \n")
			dist.mat <- as.dist(dist.obj[[dist.name]][ind, ind])
			if (block.perm == F) {
				obj <- adonis(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
			} else {
				obj <- adonis2(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
			}
			
			prmatrix(obj$aov.tab)
			permanova.obj[[dist.name]] <- obj$aov.tab
			cat("\n")
		}
		result$permanova.obj <- permanova.obj
		permanovaG.obj <- NULL
		if (block.perm == F & !is.null(PermanovaG.dist)) {
			cat('\nPERMANOVA G test combining ', paste(PermanovaG.dist, collapse=','), '\n')
			response <- array(NA, c(sum(ind), sum(ind), length(PermanovaG.dist)), dimnames=list(NULL, NULL, PermanovaG.dist))
			for (dist.name in PermanovaG.dist) {
				response[, , dist.name] <- dist.obj[[dist.name]][ind, ind]
			}
			obj <- PermanovaG2(as.formula(paste("response", formula)), df,  strata=strata, ...)
			prmatrix(obj$aov.tab)
			permanovaG.obj <- obj$aov.tab
			cat("\n")
			result$permanovaG.obj <- permanovaG.obj
		}
		cat("\n")
		sink()
	} else {
		if (pairwise == F) {
			if (is.null(adj.name)) {
				formula <- paste('~', grp.name)
			} else {
				formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
			}
			
			ind <- apply(df[, gsub("^\\s+|\\s+$", "", strsplit(strsplit(formula, "~")[[1]][2], "\\+")[[1]]), drop=F], 1, function(x) sum(is.na(x))) == 0
			df <- df[ind, ]
			if (!is.null(strata)) {
				strata <- factor(strata[ind])
			}
			sink(paste0('Beta_diversity_PERMANOVA_test_', ann, '.txt'))
			date()
			cat('\nPERMANOVA test: \n')
			permanova.obj <- list()
			for (dist.name in dist.names) {
				cat(dist.name, " distance: \n")
				dist.mat <- as.dist(dist.obj[[dist.name]][ind, ind])
				if (block.perm == F) {
					obj <- adonis(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
				} else {
					obj <- adonis2(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
				}
				prmatrix(obj$aov.tab)
				permanova.obj[[dist.name]] <- obj$aov.tab
				cat("\n")
			}
			result$permanova.obj <- permanova.obj
			permanovaG.obj <- NULL
			if (block.perm == F & !is.null(PermanovaG.dist)) {
				cat('\nPERMANOVA G test combining ', paste(PermanovaG.dist, collapse=','), '\n')
				response <- array(NA, c(sum(ind), sum(ind), length(PermanovaG.dist)), dimnames=list(NULL, NULL, PermanovaG.dist))
				for (dist.name in PermanovaG.dist) {
					response[, , dist.name] <- dist.obj[[dist.name]][ind, ind]
				}
				obj <- PermanovaG2(as.formula(paste("response", formula)), df,  strata=strata, ...)
				prmatrix(obj$aov.tab)
				permanovaG.obj <- obj$aov.tab
				cat("\n")
				result$permanovaG.obj <- permanovaG.obj
			}
			
			cat("\n")
			sink()
			
		} else {
			sink(paste0('Beta_diversity_PERMANOVA_test_', ann, '_pairwise.txt'))
			date()
			cat('\nPairwise PERMANOVA test: \n')
			grp <- factor(df[, grp.name])
			grp.levels <- levels(grp)
			grp.nlevels <- nlevels(grp)
			pmat.all <- NULL
			rmat.all <- NULL
			pmat.G <- matrix(NA, grp.nlevels, grp.nlevels)
			colnames(pmat.G) <- rownames(pmat.G) <- grp.levels
			for (dist.name in dist.names) {
				cat(dist.name, " distance: \n")
				pmat <- matrix(NA, grp.nlevels, grp.nlevels)
				colnames(pmat) <- rownames(pmat) <- grp.levels
				rmat <- matrix(NA, grp.nlevels, grp.nlevels)
				colnames(rmat) <- rownames(rmat) <- grp.levels
				for (i in 1:(grp.nlevels-1)) {
					grp.level1 <- grp.levels[i]
					for (j in (i+1):grp.nlevels) {
						
						grp.level2 <- grp.levels[j]
						cat(grp.level1, ' vs ', grp.level2, '\n')
						ind <- grp %in% c(grp.level1, grp.level2)
						df2 <- subset(df, ind)
						df2[, grp.name] <- factor(df2[, grp.name])
						dist.mat <- dist.obj[[dist.name]][ind, ind]
						
						if (!is.null(strata)) {
						  strata2 <- factor(strata[ind])
					    } else {
						  strata2 <- NULL
						}
						
						if (is.null(adj.name)) {
							formula <- paste('~', grp.name)
						} else {
							formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
						}
						
						ind2 <- apply(df2[, gsub("^\\s+|\\s+$", "", strsplit(strsplit(formula, "~")[[1]][2], "\\+")[[1]]), drop=F], 1, function(x) sum(is.na(x))) == 0
						df2 <- df2[ind2, ]
						dist.mat2 <- as.dist(dist.mat[ind2, ind2])
						
						if (!is.null(strata2)) {
						    strata2 <- factor(strata2[ind2])
						}
						
						if (block.perm == F) {
							obj <- adonis(as.formula(paste("dist.mat2", formula)), df2,  strata=strata2, ...)
						} else {
							obj <- adonis2(as.formula(paste("dist.mat2", formula)), df2,  strata=strata2,  ...)
						}
						prmatrix(obj$aov.tab)
						cat("\n")
	
						if (block.perm == F) {
							pmat[i, j] <- pmat[j, i] <- obj$aov.tab[length(adj.name)+1, 6]
							rmat[i, j] <- rmat[j, i] <- obj$aov.tab[length(adj.name)+1, 5]
						} else {
							pmat[i, j] <- pmat[j, i] <- obj$aov.tab[1, 6]
							rmat[i, j] <- rmat[j, i] <- obj$aov.tab[1, 5]
						}

						# PERMANOVA G after last distance
						if (block.perm == F & !is.null(PermanovaG.dist)) {
							if (dist.name == dist.names[length(dist.names)]) {
								cat('\nPERMANOVA G test combining ', paste(PermanovaG.dist, collapse=','), '\n')
								response <- array(NA, c(sum(ind2), sum(ind2), length(PermanovaG.dist)), dimnames=list(NULL, NULL, PermanovaG.dist))
								for (dist.name in PermanovaG.dist) {
									# Rev: 2019_05_20 Fixed an error
									response[, , dist.name] <- dist.obj[[dist.name]][ind, ind][ind2, ind2]
								}
								obj <- PermanovaG2(as.formula(paste("response", formula)), df2,  strata=strata2, ...)
								prmatrix(obj$aov.tab)
								cat("\n")
								pmat.G[i, j] <- pmat.G[j, i] <- obj$aov.tab[length(adj.name)+1, 'omni.p.value']
							}
						}
					}
				}
				cat("\n")
				pmat.all <- rbind(pmat.all, c(dist.name, rep('', grp.nlevels-1)), formatC(pmat), rep("", grp.nlevels))
				rmat.all <- rbind(rmat.all, c(dist.name, rep('', grp.nlevels-1)), formatC(rmat), rep("", grp.nlevels))
			}
			cat("\n")
			sink()
			write.csv(pmat.all, paste0('Beta_diversity_PERMANOVA_test_', ann, '_pairwiseP.csv'))
			write.csv(rmat.all, paste0('Beta_diversity_PERMANOVA_test_', ann, '_pairwiseR.csv'))
			write.csv(pmat.G, paste0('Beta_diversity_PERMANOVA_G_test_', ann, '_pairwiseP.csv'))
		}
	}
	
	return(invisible(result))
}

randomForestTest <- function(x, y, perm.no=999, ...) {
	iris.rf <- randomForest(x=x, y=y, importance=FALSE)
	to1 <- mean(abs((as.numeric(y) - 1) - predict(iris.rf,  type='prob')[, 2])^2) #mean squarred error
	to2 <- mean(y !=  predict(iris.rf,   type='response'))  #mean prediction error
	tp <- sapply(1:perm.no, function(i) {
				if (i %% 10 == 0) cat('.')
				y.p <- sample(y)
				iris.rf <- randomForest(x=x, y=y.p, importance=FALSE)
				
				t1 <- mean(abs((as.numeric(y.p) - 1) - predict(iris.rf,  type='prob')[, 2])^2) 
				t2 <- mean(y.p !=  predict(iris.rf,   type='response'))
				c(t1, t2)
			}
	
	)
	pv1 <- (sum(tp[1, ] <= to1) + 1) / (perm.no + 1)
	pv2 <- (sum(tp[2, ] <= to2) + 1) / (perm.no + 1)
	c(pv1=pv1, pv2=pv2)
}

perform_rf_test <- function (data.obj, grp.name, taxa.level='Genus', perm.no=999, prev=0.1, minp=0.000, ann='',...) {
	if (taxa.level == 'Species') {
		if (taxa.level %in% names(data.obj$abund.list)) {
			ct <- data.obj$abund.list[[taxa.level]]
		} else {
			# Accomodate different version
			ct <- data.obj$otu.tab
			rownames(ct) <- paste0("OTU", rownames(ct), ":", data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
			data.obj$abund.list[['Species']] <- ct
		}
	} else {
		ct <- data.obj$abund.list[[taxa.level]]
	}
	
	prop <- t(t(ct) / colSums(ct))
	prop <- prop[matrixStats::rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
	prop <- t(prop)
	
	grp <- data.obj$meta.dat[, grp.name]
	if (!is.factor(grp)) stop('Current random Forest test is designed to deal with binary factor data!\n')

	cat('RF test P values  ...\n')
	obj <- randomForestTest(prop, grp, perm.no, ...)
	sink(paste0('OverallAssoc_RF_test_', taxa.level, '_', ann, '.txt'))
	cat ('Prediction Probability:', obj[1], '\n')
	cat ('Binary Response:', obj[2], '\n')
	sink()
}

# Rev: 2016_12_02, MiKRAT for binary result, add out_type='D'
perform_mirkat_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		grp.name=NULL, adj.name=NULL, pairwise=F,  ann='', ...) {
	
	# MiRKAT not implemented for correlated data
	df <- data.obj$meta.dat
	
	if (pairwise == F) {
		
		ind <- apply(df[, c(grp.name, adj.name), drop=F], 1, function(x) sum(is.na(x))) == 0
		df <- df[ind, ]
		
		grp <- df[, grp.name]
		
		if (is.character(grp)) {
			grp <- factor(grp)
		}
		# Rev: 2016_12_02
		if (is.factor(grp)) {
			if (nlevels(grp) > 2) {
				stop('Currently MiRKAT only supports binary outcome!')
			} else {
				grp <- as.numeric(grp) - 1
				out_type <- 'D'
			}
		} else {
			out_type <- 'C'
		}
		if (!is.null(adj.name)) {
			# No intercept
			adj <- model.matrix(~ ., data.frame(df[, adj.name]))
			# Remove collinear terms
			qadj <- qr(adj, tol = 1e-07)
			adj <- adj[, qadj$pivot, drop = FALSE]
			adj <- adj[, 1:qadj$rank, drop = FALSE]
			# Remove intercept
			adj <- adj[, colSums(adj==1) != nrow(adj)]
		} else {
			adj <- NULL
		}
		
		sink(paste0('Beta_diversity_MiRKAT_test_', ann, '.txt'))
		date()
		cat('\nMiRKAT test combining ', paste(dist.names, collapse=','), '\n')
		Ks <- list()
		for (dist.name in dist.names) {
			Ks[[dist.name]] <- D2K(dist.obj[[dist.name]][ind, ind])
		}
		# Rev: 2016_12_02
		obj <- MiRKAT(grp, X=adj, Ks, out_type=out_type)
		cat('Individual P value: ')
		prmatrix(t(obj$indivP))
		cat('\nOmnibus P value: ')
		cat(obj$omnibus_p)
		cat("\n")
		sink()
		
	} else {
		sink(paste0('Beta_diversity_MiRKAT_test_', ann, '_pairwise.txt'))
		date()
		cat('\nPairwise MiRKAT test: \n')
		grp <- factor(df[, grp.name])
		grp.levels <- levels(grp)
		grp.nlevels <- nlevels(grp)
		pmat.G <- matrix(NA, grp.nlevels, grp.nlevels)
		colnames(pmat.G) <- rownames(pmat.G) <- grp.levels
		parr <- array(NA, c(grp.nlevels, grp.nlevels, length(dist.names)), dimnames=list(grp.levels, grp.levels, dist.names))
		
		for (i in 1:(grp.nlevels-1)) {
			grp.level1 <- grp.levels[i]
			for (j in (i+1):grp.nlevels) {
				
				grp.level2 <- grp.levels[j]
				cat(grp.level1, ' vs ', grp.level2, '\n')
				ind <- grp %in% c(grp.level1, grp.level2)
				df2 <- subset(df, ind)
				df2[, grp.name] <- factor(df2[, grp.name])
				ind2 <- apply(df2[, c(grp.name, adj.name), drop=F], 1, function(x) sum(is.na(x))) == 0
				df2 <- df[ind2, ]
				
				grp2 <- as.numeric(df2[, grp.name]) - 1
				if (!is.null(adj.name)) {
					# No intercept
					adj <- model.matrix(~ ., data.frame(df2[, adj.name]))
					# Remove collinear terms
					qadj <- qr(adj, tol = 1e-07)
					adj <- adj[, qadj$pivot, drop = FALSE]
					adj <- adj[, 1:qadj$rank, drop = FALSE]
					# Remove intercept
					adj <- adj[, colSums(adj==1) != nrow(adj)]
				} else {
					adj <- NULL
				}
				
				Ks <- list()
				for (dist.name in dist.names) {
					Ks[[dist.name]] <- D2K(dist.obj[[dist.name]][rownames(df2), rownames(df2)])
				}
				# Rev: 2016_12_02
				obj <- MiRKAT(grp2, X=adj, Ks, out_type='D')
				pmat.G[i, j] <- pmat.G[j, i] <- obj$omnibus_p
				parr[i, j, ] <- parr[j, i, ] <- obj$indivP
				cat('Individual P value: ')
				prmatrix(t(obj$indivP))
				cat('\nOmnibus P value: ')
				cat(obj$omnibus_p)
				cat("\n")
				
			}
		}
		cat("\n")
		sink()
		
		pmat.all <- NULL
		for (dist.name in dist.names) {
			pmat <- parr[, , dist.name]
			pmat.all <- rbind(pmat.all, c(dist.name, rep('', grp.nlevels-1)), formatC(pmat), rep("", grp.nlevels))
		}
		
		write.csv(pmat.all, paste0('Beta_diversity_MiRKAT_test_', ann, '_pairwiseP.csv'))
		write.csv(pmat.G, paste0('Beta_diversity_MiRKAT_O_test_', ann, '_pairwiseP.csv'))
	}
}


perform_betadisper_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), grp.name) {
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	sink('Beta_diversity_BETADISPER_test.txt')
	date()
	cat('\nBetadisper test: \n')
	for (dist.name in dist.names) {
		cat(dist.name, " distance: \n")
		dist.mat <- as.dist(dist.obj[[dist.name]])
		obj <- betadisper(dist.mat, grp)
		prmatrix(anova(obj))
		cat("\n")
	}
	cat("\n")
	sink()
}

distance_compare_test <- function (dist.mat, ind1, ind2, ind3, ID2=NULL, ID3=NULL, alternative='greater', nperm=999) {
	
	ind23 <- c(ind2, ind3)
	if (!is.null(ID2) & !is.null(ID3)) {
		ID23 <- c(ID2, ID3)
		ID2u <- unique(ID2)
		ID3u <- unique(ID3)
		ID23u <- c(ID2u, ID3u)
		n2u <- length(ID2u)
		n3u <- length(ID3u)
	}
	n2 <- length(ind2)
	n3 <- length(ind3)
	
	dist12 <- dist.mat[ind1, ind2]
	dist13 <- dist.mat[ind1, ind3]
	
	stat.obs <- mean(dist12) - mean(dist13)
	stat.perm <- sapply(1:nperm, function(i) {
				if (is.null(ID2) & is.null(ID3)) {
					ind23.p <- sample(ind23)
					ind2.p <- ind23.p[1:n2]
					ind3.p <- ind23.p[(n2+1):(n2+n3)]
				} else {
					ID23u.p <- sample(ID23u)
					ID2.p <- ID23u.p[1:n2u]
					ID3.p <- ID23u.p[(n2u+1):(n2u+n3u)]
					ind2.p <- ind23[ID23 %in% ID2.p]
					ind3.p <- ind23[ID23 %in% ID3.p]
				}

				dist12.p <- dist.mat[ind1, ind2.p]
				dist13.p <- dist.mat[ind1, ind3.p]
				mean(dist12.p) - mean(dist13.p)
			})
	if (alternative == 'greater') {
		pv <- mean(c(stat.perm >= stat.obs, TRUE))
	} else {
		pv <- mean(c(stat.perm <= stat.obs, TRUE))
	}  
	pv
}

# Note that this is the block-permutation scheme, not within-sbuject permutation
perform_distance_compare_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
		grp.name, level1, level2, level3, subject=NULL, alternative='greater', nperm=999, seed=123, ann='') {
	sink(paste0('Beta_diversity_test_', level1, '-', level2, '_', alternative, '_', level1, '-', level3, '_', ann, '.txt'))
	set.seed(seed)
	cat('Testing the distance ', level1, '-', level2, ' is ', alternative, ' than ', level1, '-', level3, '\n')
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	IDs <- df[, subject]
	grp.levels <- levels(grp)
	grp.nlevels <- nlevels(grp)
	ind1 <- which(grp == level1)
	ind2 <- which(grp == level2)
	ind3 <- which(grp == level3)
	if (is.null(subject)) {
		ID2 <- NULL
		ID3 <- NULL
	} else {
		ID2 <- as.character(IDs[grp == level2])
		ID3 <- as.character(IDs[grp == level3])
	}
	for (dist.name in dist.names) {
		cat(dist.name, ' Distance:p=') 
		pv <- distance_compare_test(dist.obj[[dist.name]], ind1, ind2, ind3, ID2, ID3, alternative, nperm)
		cat(pv, '\n')
	}
	sink()
	
}


perform_taxa_compare_test <- function (data.obj, grp.name, level1, level2, level3, alternative='greater', nperm=1000, seed=123,
		taxa.levels=c('Phylum', 'Family', 'Genus'), taxa.name='All', prev=0.1, minp=0.002, ann='') {
	
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	grp.levels <- levels(grp)
	grp.nlevels <- nlevels(grp)
	ind1 <- which(grp == level1)
	ind2 <- which(grp == level2)
	ind3 <- which(grp == level3)
	
	pv.list <- list()
	
	for (LOI in taxa.levels) {
		cat(LOI, "\n")
		prop <- data.obj$abund.list[[LOI]]
		prop <- t(t(prop) / colSums(prop))
		pv.vec <- NULL
		if ('All' %in% taxa.name) {
			prop <- prop[matrixStats::rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
		} else {
			prop <- prop[taxa.name, , drop=FALSE]
		}
		
		for (taxon in rownames(prop)) {
			cat(".")
			dist.mat <- as.matrix(dist(rank(prop[taxon, ])))
			pv <- distance_compare_test(dist.mat, ind1, ind2, ind3, alternative, nperm)
			pv.vec <- c(pv.vec, pv)
		}
		names(pv.vec) <- rownames(prop)
		pv.list[[LOI]] <- pv.vec
		cat("\n")
	}
	
	for (LOI in taxa.levels) {
		pv.vec <- pv.list[[LOI]]

		qv.vec <- p.adjust(pv.vec, 'fdr')

		
		write.csv(data.frame("p value"=pv.vec, "q value"=qv.vec), paste0("Taxa_TrendAnalysis_", LOI, "_", ann, ".csv"))
	}
	
}

# New: 2017_11_03
# Now default is averaged distance for within or between-subject comparison
perform_site_correlation_test <- function (data.obj, dist.obj, site.name, sites, subject, 
		dist.names = c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), dist.ave = TRUE, nPerm = 999, ann = '') {

	pmat <- array(NA, c(length(dist.names), length(sites), length(sites)), dimnames=list(dist.names, sites, sites))
	for (site1 in sites) {
		for(site2 in sites) {
			if (site1 != site2 & is.na(pmat[1, site1, site2])) {
				cat('.')
				sampling.site <- c(site1, site2)
				samIDs <- data.obj$meta.dat[, site.name] %in% sampling.site 
				data.obj1 <- subset_data(data.obj, samIDs)
				dist.obj1 <- subset_dist(dist.obj, samIDs)
				
# Compute the observed distance within subject
				lIDs <- as.character(data.obj1$meta.dat[, site.name])
				pIDs <- as.character(data.obj1$meta.dat[, subject])
				uIDs <- as.character(unique(data.obj1$meta.dat[, subject]))
				
				
				pdf(paste('Permutation_test', site1, site2, 'Correlation_', ann, '.pdf', sep='_'), width=6, height=6)
				for (dist.name in dist.names) {
					cat(dist.name, '...\n')
					
					dist.mat <- dist.obj1[[dist.name]]					
# To self
					dist.o1 <- NULL
					for (uID in uIDs) {
						ind1 <- pIDs == uID & lIDs == site1
						ind2 <- pIDs == uID & lIDs == site2
						
						if (sum(ind1) != 0 & sum(ind2) != 0)
							if (dist.ave) {
								dist.o1 <- c(dist.o1, mean(dist.mat[ind1, ind2]))
							} else {
								dist.o1 <- c(dist.o1, dist.mat[ind1, ind2])
							}
							
					}
					
# To other people
					dist.o2 <- NULL
					for (uID in uIDs) {
						ind1 <- pIDs == uID & lIDs == site1
						ind2 <- pIDs != uID & lIDs == site2
						
						if (sum(ind1) != 0 & sum(ind2) != 0)
							if (dist.ave) {
								dist.o2 <- c(dist.o2, mean(dist.mat[ind1, ind2]))
							} else {
								dist.o2 <- c(dist.o2, dist.mat[ind1, ind2])
							}
							
					}
					
					if (!is.null(dist.o1) & !is.null(dist.o2)) {
						boxplot(list(ToSelf=dist.o1, ToOther=dist.o2), col='steelblue', ylab=paste(dist.name, 'distance'))
						stat.o <- mean(dist.o2) - mean(dist.o1)
						
						stat.p <- NULL
						for (i in 1:nPerm) {
							# Exchange the uterus between people :D
							pIDs.p <- pIDs
							pIDs.uterus <- pIDs[lIDs == site1]
							temp <- factor(pIDs.uterus)
							levels(temp) <- sample(levels(temp))
							pIDs.uterus.p <- as.character(temp)
							pIDs.p[lIDs == site1] <- pIDs.uterus.p
							
							# To increase the number of permutation
							pIDs.uterus <- pIDs[lIDs == site2]
							temp <- factor(pIDs.uterus)
							levels(temp) <- sample(levels(temp))
							pIDs.uterus.p <- as.character(temp)
							pIDs.p[lIDs == site2] <- pIDs.uterus.p
							
							# To self
							dist.o1 <- NULL
							for (uID in uIDs) {
								ind1 <- pIDs.p == uID & lIDs == site1
								ind2 <- pIDs.p == uID & lIDs == site2
								
								if (sum(ind1) != 0 & sum(ind2) != 0)
									if (dist.ave) {
										dist.o1 <- c(dist.o1, mean(dist.mat[ind1, ind2]))
									} else {
										dist.o1 <- c(dist.o1, dist.mat[ind1, ind2])
									}
							}
							
# To other people
							dist.o2 <- NULL
							for (uID in uIDs) {
								ind1 <- pIDs.p == uID & lIDs == site1
								ind2 <- pIDs.p != uID & lIDs == site2
								
								if (sum(ind1) != 0 & sum(ind2) != 0)
									if (dist.ave) {
										dist.o2 <- c(dist.o2, mean(dist.mat[ind1, ind2]))
									} else {
										dist.o2 <- c(dist.o2, dist.mat[ind1, ind2])
									}
							}
							
							# boxplot(list(ToSelf=dist.o1, ToOther=dist.o2))
							stat.p <- c(stat.p, mean(dist.o2) - mean(dist.o1))
						}
						
						pv <- mean(c(stat.p >= stat.o, TRUE))
						
						hist(stat.p, xlab=paste(dist.name, 'distance'), col='steelblue', xlim=c(min(stat.p) * 2, max(stat.p) * 2),
								main='Average distances observed under permutation', sub=paste('P < ', pv))
						
						abline(v=stat.o, col='red')
						pmat[dist.name, site1, site2] <- pmat[dist.name, site2, site1] <- pv
					}

				}
				dev.off()
			}
			
		}
	}
	cat('Finished!\n')
	
	pmat.csv <- NULL
	for (dist.name in dist.names) {
		temp <- as.matrix(pmat[dist.name, ,  ])
		mode(temp) <- 'character'
		pmat.csv <- rbind(pmat.csv, rbind(Distance=c(dist.name, rep('', ncol(temp) - 1)), temp, ''))
	}
	
	write.csv(pmat.csv, paste0('SiteCorrelationP_', ann, '.csv'))
	return(invisible(pmat))
}

# New: 2017_11_03
perform_site_taxa_correlation_test <- function (data.obj, site.name, sites, subject, 
		taxa.levels = c('Genus'), taxa.names = NULL, dist.func = function (x) {dist(sqrt(x))}, 
		prev=0.1, minp=0.002, medianp=NULL, ann = '', ...) {
	pmat.list <- list()
	for (taxa.level in taxa.levels) {
		cat(taxa.level, '...\n')
		genus <- data.obj$abund.list[[taxa.level]]
		genus <- t(t(genus) / colSums(genus))
		if (is.null(taxa.names)) {
			ind <- rep(TRUE, nrow(genus))
			if (!is.null(prev)) ind <- ind & rowMeans(genus != 0) > prev
			if (!is.null(minp)) ind <- ind & matrixStats::rowMaxs(genus) > minp
			if (!is.null(medianp)) {
				nz.mean <- apply(genus, 1, function(x) median(x[x!=0]))
				ind <- ind & nz.mean > medianp
			}
			genus <- genus[ind, , drop = FALSE]
			gen.names <- rownames(genus)
		} else {
			if (sum(!(taxa.names %in% rownames(genus))) != 0) {
				stop('Could not find the taxa! Please check!\n')
			}
			genus <- genus[taxa.names, , drop = FALSE]
			gen.names <- taxa.names
		}

		# Convert into distance
	    dist.obj <- list()
		for (gen.name in gen.names) {
			dist.obj[[gen.name]] <- as.matrix(dist.func(genus[gen.name, ]))
		}
		
		pmat.list[[taxa.level]] <- perform_site_correlation_test(data.obj, dist.obj, site.name, sites, subject, 
				dist.names = gen.names,  ann = paste(taxa.level, ann, sep = '_'), ...)
		
	}
	return(invisible(pmat.list))
}

# Rev: 2017_02_19 Add hei and wid
# Rev: 2018_01_15 Add qt.outlier=0.097
# Rev: 2018_03_06 position_jitter(h = 0), coord_cartesian
# Rev: 2018_07_16 
generate_taxa_boxplot <- function (data.obj,  grp.name, strata=NULL, scale='P',  taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'), 
		taxa.name='All', pseudo.ct=0.5, rm.outlier=TRUE, qt.outlier=0.97, prev=0.1, minp=0.002, ann='All', 
		subject=NULL, l.size=0.5, p.size=2.5, hei=NULL, wid=NULL) {
	# To be completed
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	if (!is.null(subject)) {
		ID <- df[, subject]
	}
	
	if (is.null(hei)) {
		hei <- 5
	} else {
		hei <- hei
	}
	
	if (is.null(strata)) {
		if (is.null(wid)) {
			wid <- 5
		} else {
			wid <- wid
		}	
	} else {
		if (is.null(wid)) {
			wid <- 5.5
		} else {
			wid <- wid
		}
	}
	
	if (scale == 'P') {
		ylab <- 'Proportion'
	} 
	if (scale == 'logP') {
		ylab <- 'log10(Proportion)'
	}
	if (scale == 'sqrtP') {
		ylab <- 'sqrt(Proportion)'
	}
	if (scale == 'binary') {
		ylab <- 'Count'
	} 
	if (scale == 'none') {
		ylab <- 'Value'
	} 
	
	for (LOI in taxa.levels) {
		if (LOI == 'All') {
			if ('All' %in% taxa.name) stop('Taxa names are not required to be all when taxa level is also all!\n')
			headnames <- NULL
			prop <- NULL
			for (LOI2 in names(data.obj$abund.list)) {
				prop2 <- data.obj$abund.list[[LOI2]]
				taxa.name2 <- taxa.name[taxa.name %in% rownames(prop2)]
				if (length(taxa.name2) != 0) {
					if (scale == 'logP') {
						prop2 <- prop2 + pseudo.ct
					} 
					prop2 <- t(t(prop2) / colSums(prop2))			
					headnames2 <- sapply(strsplit(rownames(prop2), ";"), paste, collapse="\n")
					names(headnames2) <- rownames(prop2)
					prop <- rbind(prop, prop2[taxa.name2, , drop=FALSE])
					headnames <- c(headnames, headnames2)
				}
			}
			
		} else {

			cat(LOI, "\n")
			ct <- data.obj$abund.list[[LOI]]
			prop <- t(t(ct) / colSums(ct))
			ind <- matrixStats::rowMaxs(prop) > minp & rowSums(prop != 0) > prev * ncol(prop)
			
			if (scale == 'logP') {
				ct <- ct + pseudo.ct
				prop <- t(t(ct) / colSums(ct))
			} 
			
			headnames <- sapply(strsplit(rownames(prop), ";"), paste, collapse="\n")
			names(headnames) <- rownames(prop)
			
			if ('All' %in% taxa.name) {
				prop <- prop[ind, , drop=FALSE]
			} else {
				# Rev: 2018_07_16
				prop <- prop[rownames(prop) %in% taxa.name, , drop=FALSE]
			}
			
			
		}
		
		if (scale == 'logP') {
			prop <- log10(prop)
		}
		
		if (scale == 'sqrtP') {
			prop <- sqrt(prop)
		}
		
		if (scale == 'binary') {
			temp <- prop != 0
			prop[temp] <- 'Presence'
			prop[!temp] <- 'Absence'
		}

		
		pdf(paste("Taxa_Boxplot", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)
		
		if (is.null(strata)) {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				if (scale != 'binary') {
					if (scale == 'P') {
						if (rm.outlier == T) {
							ylims <- c(min(taxon.abund),  min(quantile(taxon.abund, qt.outlier) * 1.25, 1))
						} else {
							ylims <- range(taxon.abund)
						}
						
					} else {
						ylims <- range(taxon.abund)
					}
					
					if (is.null(subject)) {
						df2 <- data.frame(Value=taxon.abund, Group=grp)	
						dodge <- position_dodge(width=0.9)
						obj <- ggplot(df2, aes(x=Group, y=Value, col=Group)) +
								geom_boxplot(position=dodge, outlier.colour = NA, alpha=0.75, lwd=0.35) + 
								geom_jitter(alpha=0.6, size=p.size,  position = position_jitter(w = 0.1, h = 0)) +
								labs(y=ylab, title=headnames[taxon]) 
						if (rm.outlier) {
							obj <- obj + coord_cartesian(ylim = ylims)
						}
						#		ylim(ylims[1], ylims[2]) +
						obj <- obj + theme(legend.position="none")
						print(obj)
					} else {						
						df2 <- data.frame(Value=taxon.abund, Group=grp, subject=ID)	
						dodge <- position_dodge(width=0.9)
						obj <- ggplot(df2, aes(x=Group, y=Value, shape=Group, group=subject)) +
								geom_point(size=p.size) +
								geom_line(size=l.size) +
								labs(y=ylab, title=headnames[taxon]) +
							#	ylim(ylims) +
								theme(legend.position="none")
						print(obj)
						
					}

				} else {
					df2 <- data.frame(Value=taxon.abund, Group=grp)	
					
					obj <- ggplot(df2, aes(x=Group, fill=Value)) +
							geom_bar(width=.5) +
							labs(y=ylab, title=headnames[taxon]) +
							theme(legend.title=element_blank())
					
					print(obj)
				}
				
			}	
		} else {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				if (scale != 'binary') {
					if (scale == 'P') {
						if (rm.outlier == T) {
							ylims <- c(min(taxon.abund),  min(quantile(taxon.abund, qt.outlier) * 1.25, 1))
						} else {
							ylims <- range(taxon.abund)
						}
						
					} else {
						ylims <- range(taxon.abund)
					}
					grp2 <- df[, strata]
					df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2)
					dodge <- position_dodge(width=0.75)
					obj <- ggplot(df2, aes(x=Strata, y=Value, col=Group, fill=Group)) +
							geom_boxplot(fill='white',  position=dodge, outlier.colour = NA, alpha=0.75, width=0.65, lwd=0.35) + 
							geom_point(position=position_jitterdodge(dodge.width=0.75), size=p.size, alpha=0.6) +
							labs(y=ylab, x=strata, title=headnames[taxon])
					if (rm.outlier) {
						obj <- obj + coord_cartesian(ylim = ylims)
					}
					#		ylim(ylims[1], ylims[2]) +
					obj <- obj + theme(legend.title=element_blank())
					print(obj)
				} else {
					grp2 <- df[, strata]
					df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2)
					obj <- ggplot(df2, aes(x=Group, fill=Value)) +
							geom_bar(width=.5) +
							labs(y=ylab, title=headnames[taxon]) +
							facet_wrap(~ Strata) + 
							theme(legend.title=element_blank())
					
					print(obj)
				}
			}
			
		}
		dev.off()
	}
	
}

# Rev: 2017_02_19 Add Pseudo.ct, hei0, wid0; strata2; smooth.method
# Rev: 2018_01_15 Add qt.outlier 
# Rev: 2018_07_16 Add an option to combine strata in the same plot, coord_cartesian
# Rev: 2019_01_03 Selecting taxa before adding pseudo-ct
# Note: if rm.outlier=TRUE, the smooth function only acts on the points in the plot region. May not be desirable
generate_taxa_scatterplot <- function (data.obj,  grp.name, subject=NULL,
		taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'), taxa.name='All', prev=0.1, minp=0.002, 
		strata=NULL, strata2=NULL, combine.strata=FALSE, combine.strata2=FALSE, smooth.method='loess', 
		pt.shape=16, pt.alpha=0.5, pt.size=2, subject.pt=FALSE,
		scale='P', pseudo.ct=0.5, is.pvalue=FALSE, rm.outlier=FALSE, qt.outlier=0.97, 
		ann='All', hei=NULL, wid=NULL) {
	# To be completed
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	if(is.null(subject)) {
		ID <- NA
	} else {
		ID <- factor(df[, subject])
	}
	
	if (is.null(hei) | is.null(wid)) {
		hei <- 5
		wid <- 6
		if (scale != 'binary') {
			if (!is.null(strata) & combine.strata == FALSE) {
				wid <- 4 * nlevels(factor(df[, strata]))
			} 
			if (!is.null(strata2) & combine.strata2 == FALSE) {
				hei <- 2.5 * nlevels(factor(df[, strata2]))
			} 
		} else {
			if (!is.null(strata)) {
				wid <- 2.5 * nlevels(factor(df[, strata]))
			} 
			if (!is.null(strata2)) {
				hei <- 2.5 * nlevels(factor(df[, strata2]))
			} 
		}

	}

	if (scale == 'P') {
		ylab <- 'Proportion'
	} 
	if (scale == 'none') {
		ylab <- 'Value'
	} 
	if (scale == 'logP') {
		ylab <- 'log10(Proportion)'
	}
	if (scale == 'sqrtP') {
		ylab <- 'sqrt(Proportion)'
	}
	if (scale == 'binary') {
		ylab <- 'Count'
	} 
	
	for (LOI in taxa.levels) {
		if (LOI == 'All') {
			if ('All' %in% taxa.name) stop('Taxa names need to be specified when taxa level is  all!\n')
			headnames <- NULL
			prop <- NULL
			for (LOI2 in names(data.obj$abund.list)) {
				prop2 <- data.obj$abund.list[[LOI2]]
				taxa.name2 <- taxa.name[taxa.name %in% rownames(prop2)]
				if (length(taxa.name2) != 0) {
					if (scale == 'logP') {
						prop2 <- prop2 + pseudo.ct
					} 
					prop2 <- t(t(prop2) / colSums(prop2))			
					headnames2 <- sapply(strsplit(rownames(prop2), ";"), paste, collapse="\n")
					names(headnames2) <- rownames(prop2)
					prop <- rbind(prop, prop2[taxa.name2, , drop=FALSE])
					headnames <- c(headnames, headnames2)
				}
			}
			
		} else {
			
			cat(LOI, "\n")
			ct <- data.obj$abund.list[[LOI]]
			prop <- t(t(ct) / colSums(ct))
			ind <- matrixStats::rowMaxs(prop) > minp & rowSums(prop != 0) > prev * ncol(prop)
			
			if (scale == 'logP') {
				ct <- ct + pseudo.ct
				prop <- t(t(ct) / colSums(ct))
			} 
			
			headnames <- sapply(strsplit(rownames(prop), ";"), paste, collapse="\n")
			names(headnames) <- rownames(prop)
			
			if ('All' %in% taxa.name) {
				prop <- prop[ind, , drop=FALSE]
			} else {
				prop <- prop[rownames(prop) %in% taxa.name, , drop=FALSE]
			}
			
		}
		
		if (scale == 'logP') {
			prop <- log10(prop)
		}
		
		if (scale == 'sqrtP') {
			prop <- sqrt(prop)
		}
		if (scale == 'asinsqrtP') {
			prop <- asin(sqrt(prop))
		}
		
		if (scale == 'binary') {
			temp <- prop != 0
			prop[temp] <- 'Presence'
			prop[!temp] <- 'Absence'
		}
		
		
		pdf(paste("Taxa_Scatterplot", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)
		if (is.null(strata)) {
			if (is.pvalue == TRUE) {
				if (scale != 'binary') {
					pvalues <- apply(prop, 1, function (x) {
								formatC(cor.test(as.numeric(x), grp, method='spearman')$p.value, digits=3)
							})
				} else {
					pvalues <- apply(prop, 1, function (x) {
								if (nlevels(factor(x)) > 1) {
									return(formatC(wilcox.test(grp ~ x)$p.value, digits=3))
								} else {
									return('1.0')
								}
								
							})
				}

			} 
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				if (scale != 'binary') {
					if (scale == 'P') {
						if (rm.outlier == T) {
							ylims <- c(min(taxon.abund),  min(quantile(taxon.abund, qt.outlier) * 1.25, 1))
						} else {
							ylims <- range(taxon.abund)
						}
					} else {
						ylims <- range(taxon.abund)
					}
					df2 <- data.frame(Value=taxon.abund, Group=grp, ID=ID)	
					dodge <- position_dodge(width=0.9)
					obj <- ggplot(df2, aes(x=Group, y=Value, group=ID)) +
					#		geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
							geom_smooth(method=smooth.method, aes(group=NA), size=0.75) 
					if (is.pvalue == TRUE) {
						obj <- obj + labs(y=ylab, x=grp.name, title=paste0(headnames[taxon], '(P=', pvalues[taxon], ')'))
					} else {
						obj <- obj + labs(y=ylab, x=grp.name, title=paste0(headnames[taxon]))
					}
					if (rm.outlier) {
						obj <- obj + coord_cartesian(ylim = ylims)
					}
                    obj <- obj + theme(legend.position="none")
					if (!is.null(subject)) {
						if (subject.pt == FALSE) {
							obj <- obj + geom_line(size=0.25, alpha=0.75) 
						}  else {
							obj <- obj + geom_line(size=0.25, alpha=0.75) + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
						}
						
					} else {
						obj <- obj + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
					}
					print(obj)
				} else {
					df2 <- data.frame(Value=taxon.abund, Group=grp, ID=ID)	
					dodge <- 
							obj <- ggplot(df2, aes(y=Group, x=Value)) +
							geom_boxplot(position=position_dodge(width=0.9), outlier.colour = NA, alpha=0.75) + 
							geom_jitter(alpha=pt.alpha, size=pt.size,  position = position_jitter(w = 0.1, h = 0))
					
					if (is.pvalue == TRUE) {
						obj <- obj + labs(y=grp.name, x='Presence/Absence', title=paste0(headnames[taxon], '(P=', pvalues[taxon], ')'))
					} else {
						obj <- obj + labs(y=grp.name, x='Presence/Absence', title=paste0(headnames[taxon]))
					}
					obj <- obj + theme(legend.position="none")
					
					
					print(obj)
				}
				
			}	
		} else {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				if (scale != 'binary') {
					if (scale == 'P') {
						if (rm.outlier == T) {
							ylims <- c(min(taxon.abund), min(quantile(taxon.abund, qt.outlier) * 1.25, 1))
						} else {
							ylims <- range(taxon.abund)
						}
						
					} else {
						ylims <- range(taxon.abund)
					}
					grp2 <- factor(df[, strata])
					if (!is.null(strata2)) {
						grp3 <- factor(df[, strata2])
					} else {
						grp3 <- grp2
					}
					grp4 <- factor(paste(grp2, grp3))
					
					# VERY bad names with only difference in capital letters, variable names are also confusing
			        # Rev: 2018_07_16 
			
					df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2, Strata2=grp3, Strata3=grp4, ID=ID)
					
					
					if (is.null(strata2)) {
						if (combine.strata == TRUE) {
							obj <- ggplot(df2, aes(x=Group, y=Value, group=ID, col=Strata)) +
						#			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
									geom_smooth(method=smooth.method, aes(group=Strata), size=0.75) +
									labs(y=ylab, x=grp.name, title=headnames[taxon])
						} else {
							obj <- ggplot(df2, aes(x=Group, y=Value, group=ID)) +
						#			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
									geom_smooth(method=smooth.method, aes(group=NA), size=0.75) +
									labs(y=ylab, x=grp.name, title=headnames[taxon]) + 
									facet_grid(. ~ Strata)
						}

					} else {
						if (combine.strata == TRUE) {
							if (combine.strata2 == TRUE) {
								obj <- ggplot(df2, aes(x=Group, y=Value, group=ID, col=Strata3)) +
							#			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
										geom_smooth(method=smooth.method, aes(group=Strata3), size=0.75) +
										labs(y=ylab, x=grp.name, title=headnames[taxon])	
							} else {
								obj <- ggplot(df2, aes(x=Group, y=Value, group=ID, col=Strata)) +
							#			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
										geom_smooth(method=smooth.method, aes(group=Strata), size=0.75) +
										labs(y=ylab, x=grp.name, title=headnames[taxon]) +
										facet_grid(Strata2 ~ .)	
							}
						} else {
							if (combine.strata2 == TRUE) {
								warning('Currently, combine.strata2 will be forced to be FALSE if combine.strata=FALSE!\n')
							}
							obj <- ggplot(df2, aes(x=Group, y=Value, group=ID)) +
							#		geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
									geom_smooth(method=smooth.method, aes(group=NA), size=0.75) +
									labs(y=ylab, x=grp.name, title=headnames[taxon]) + 
									facet_grid(Strata2 ~ Strata)
						}

					}		
					
					if (rm.outlier) {
						obj <- obj + coord_cartesian(ylim = ylims)
					}
					
					obj <- obj + theme(legend.title=element_blank())	
					
					if (!is.null(subject)) {
						if (subject.pt == FALSE) {
							obj <- obj + geom_line(size=0.25, alpha=0.75) 
						}  else {
							obj <- obj + geom_line(size=0.25, alpha=0.75) + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
						}
						
					} else {
						obj <- obj + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) 
					}
					print(obj)
					
				} else {				
					grp2 <- df[, strata]
					if (!is.null(strata2)) {
						grp3 <- df[, strata2]
					} else {
						grp3 <- grp2
					}
					df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2, Strata2=grp3)
					dodge <- position_dodge(width=0.9)
							
					if (is.null(strata2)) {
						obj <- ggplot(df2, aes(x=Strata, y=Group, col=Value, fill=Value)) +
								geom_boxplot(fill='white', position=dodge, outlier.colour = NA, alpha=0.75) +
								geom_point(position=position_jitterdodge(dodge.width=0.9), size=pt.size, alpha=pt.alpha) +
								labs(y=grp.name, x=strata, title=headnames[taxon]) 
					}  else {
						obj <- ggplot(df2, aes(x=Strata, y=Group, col=Value, fill=Value)) +
								geom_boxplot(fill='white', position=dodge, outlier.colour = NA, alpha=0.75) +
								geom_point(position=position_jitterdodge(dodge.width=0.9), size=pt.size, alpha=pt.alpha) +
								facet_grid(Strata2 ~ .) +
						        labs(y=grp.name, x=strata, title=headnames[taxon]) 
					}

					print(obj)
				}
			}
			
		}
		dev.off()
	}		
}

# New: 2019_01_03 - e.g., study the correlation between different sampling methods
# Restrict to two-level factor and pairs - labelling outlier has not been completed
generate_taxa_correlationplot <- function(data.obj,  grp.name, subject.name = NULL, strata = NULL, 
		scale='sqrt', pseudo.ct = 0.5, error = 'se', n.bt = 50,  type = c('logP', 'Cor'), p.cutoff = 0.01,
		taxa.levels=c('Species'), taxa.name = 'All', prev = 0.1, minp = 0.002,
		k.outlier = 5, err.wid = 0, err.hei = 0, xlab = 'Grp1.mean', ylab = 'Grp2.mean', 
		wid = 7.5, hei = 6, pt.size = 2,   ann = '') {
	
	type <- match.arg(type)
	# Strata = to color strata
	# We will only focus on pairs
	if (type == 'Cor' & is.null(subject.name)) stop('Correlation measure requires paired samples!\n')
	df <- data.obj$meta.dat
	grp <- factor(df[, grp.name])
	grp2 <- as.numeric(grp)
	ord <- 1:length(grp)
	
	if (!is.null(subject.name)) {
		subject <- factor(df[, subject.name])
		
		# Get one sample for each subject in each group
		ind <- as.vector(tapply(rownames(df), list(grp, subject), function (x) {
							if (length(x) < 1) {
								return(NULL)
							} else {
								sample(x, 1)
							}
						}))
		ind <- ind[!is.na(ind)]
		data.obj <- subset_data(data.obj, ind)
		
		df <- data.obj$meta.dat
		grp <- factor(df[, grp.name])
		subject <- factor(df[, subject.name])
		
		tab <- table(subject)
		ind <- subject %in% names(tab)[tab == 2]
		
		data.obj <- subset_data(data.obj, ind)
		df <- data.obj$meta.dat
		grp <- factor(df[, grp.name])
		subject <- factor(df[, subject.name])
		
		ord <- order(subject)
		df <- df[ord, ]
		grp <- grp[ord]
		grp2 <- as.numeric(grp)
		subject <- subject[ord]
	}
	
	plot.data.list <- list()
	for (i in 1:length(taxa.levels)) {
		LOI <- taxa.levels[i]
		cat(LOI, "\n")
		
		ct <- data.obj$abund.list[[LOI]][, ord]
		
		prop <- t(t(ct) / colSums(ct))
		ind <- matrixStats::rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop)
		
		if (scale == 'log') {
			ct <- ct + pseudo.ct
			prop <- t(t(ct / colSums(ct)))
		}
		
		if ('All' %in% taxa.name) {
			prop <- prop[ind, , drop=FALSE]
		} else {
			# Rev: 2018_07_16
			prop <- prop[rownames(prop) %in% taxa.name, , drop=FALSE]
		}
		
		if (is.null(strata)) {
			
			# Calculate spearman correlation
			df2 <- t(apply(prop, 1, function (x) {
								x1 <- x[grp2 == 1]
								x2 <- x[grp2 == 2]
								if (error == 'se') scale.factor <- 1
								if (error == 'ci') scale.factor <- 1.96	
								se.x1 <- sd(sapply(1:n.bt, function (i) mean(sample(x1, repl=TRUE)))) * scale.factor
								se.x2 <- sd(sapply(1:n.bt, function (i) mean(sample(x2, repl=TRUE)))) * scale.factor
								m.x1 <- mean(x1)
								m.x2 <- mean(x2)
								if (!is.null(subject.name)) {
									cor.val <- cor(x1, x2, method = 'spearman')
									p.val <- wilcox.test(x1, x2, pair = TRUE)$p.value
								} else {
									cor.val <- NA
									p.val <- wilcox.test(x1, x2)$p.value
								}
								xmin <- m.x1 - se.x1
								xmax <- m.x1 + se.x1
								ymin <- m.x2 - se.x2
								ymax <- m.x2 + se.x2
								if (xmin <= 0) xmin <- 0.00001
								if (xmax >= 1) xmax <- 1 - 0.00001
								if (ymin <= 0) ymin <- 0.00001
								if (ymax >= 1) ymax <- 1 - 0.00001
								p.val.b <- as.numeric(p.val <= p.cutoff)
								
								c(m.x1, m.x2, xmin, xmax, ymin, ymax, cor.val, log10(p.val), p.val.b)
							}))
			
			df2 <- data.frame(rownames(prop), df2)
			colnames(df2) <- c('Taxa', 'x', 'y', 'xmin', 'xmax', 'ymin', 'ymax', 'Cor', 'Pval', 'Pval.b')
			df2[, 'Pval.b'] <- factor(ifelse(df2[, 'Pval.b'] == 1, paste0('p <=', p.cutoff), paste0('p >', p.cutoff)),
					levels = c(paste0('p >', p.cutoff), paste0('p <=', p.cutoff)))
			#df2[, 'Pval.b'] <- factor(df2[, 'Pval.b'] )
			#df2 <- df2[!is.na(df2$Cor), ]
			
			obj1 <-  ggplot(data = df2, aes(x = x,y = y)) + 
					geom_errorbar(aes(ymin = ymin, ymax = ymax), width = err.wid, color = 'gray') + 
					geom_errorbarh(aes(xmin = xmin, xmax = xmax), height = err.hei, color = 'gray') 
			if (type == "Cor" & !is.null(subject.name)) {
				obj1 <- obj1 + geom_point(aes(fill = Cor, shape = Pval.b), size = pt.size) 
			} 
			if (type == "logP") {
				obj1 <- obj1 + geom_point(aes(fill = Pval, shape = Pval.b), size = pt.size) 
			}
			obj1 <- obj1 +
					geom_abline(intercept=0, slope=1, color = 'gray') + 
					scale_shape_manual(values = c(21, 24)) +
                    xlab(xlab) +
					ylab(ylab) +
					labs(fill= ifelse(type == 'logP', 'log10(P-val)', 'Correlation'), shape = NULL) 
			
			if (scale == 'sqrt') {
				# obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
				obj1 <- obj1 + scale_y_sqrt(
								breaks = trans_breaks("sqrt", function(x) x^2),
								labels = trans_format("sqrt", math_format(.x^2))) +
						scale_x_sqrt(
								breaks = trans_breaks("sqrt", function(x) x^2),
								labels = trans_format("sqrt", math_format(.x^2)))
			}
			
			# To be revised
			if (scale == 'log') {
				obj1 <- obj1 + scale_y_log10(
								breaks = trans_breaks("log10", function(x) 10^x),
								labels = trans_format("log10", math_format(10^.x))) + 
						scale_x_log10(
								breaks = trans_breaks("log10", function(x) 10^x),
								labels = trans_format("log10", math_format(10^.x)))
			}
			if (scale == 'boxcox') {
				obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3)) +
						scale_x_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
			}
			
		} 
		
		plot.data.list[[LOI]] <- df2
		pdf(paste("Taxa_Correlation_", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)		
		print(obj1)
		dev.off()		
	}
	return(invisible(plot.data.list))
}

# Rev: 2018_03_08 add bootstrap standard error (normal approximation)
# Rev: 2019_01_03 handling zeros for log scales
taxa_barplot_aggregate <- function (prop, df, grp.name, strata=NULL, scale='sqrt',  xsize=10, ylab='Proportion', error='ci') {
	
	if (scale == 'log') {
		# Half of the non-zero proportion
		prop <- t(apply(prop, 1, function (x) {
					temp <- min(x[x != 0])
					x[x == 0] <- temp / 2
					x
				}))
	}
	
    grp <- factor(df[, grp.name])
	
	if (is.null(strata)) {
		df2 <- data.frame(Group=grp, t(prop))
		df2 <- melt(df2)
		colnames(df2) <- c('Group', 'Taxa', 'Value')
		
		# Could be revised
		temp1 <- aggregate(Value ~ Group + Taxa, df2, mean)

		if (error == 'se') {
			temp2 <- aggregate(Value ~ Group + Taxa, df2, function(x) {
						#sd(x) / sqrt(length(x))
						sd(sapply(1:50, function(i) mean(sample(x, repl=TRUE))))
					})
		}
		if (error == 'ci') {
			temp2 <- aggregate(Value ~ Group + Taxa, df2, function(x) {
						#sd(x) / sqrt(length(x))
						2 * sd(sapply(1:50, function (i) mean(sample(x, repl=TRUE))))
					})
		}
		
		
		df2 <- cbind(temp1, temp2[, 3])
		colnames(df2) <- c('Group', 'Taxa', 'Mean', 'SE')
		

		limits <- aes(ymax = Mean + SE, ymin = ifelse(Mean - SE > 0, Mean - SE, 0))
		dodge <- position_dodge(width=0.90)
		
		obj1 <- ggplot(df2, aes(x=Taxa, y=Mean, fill=Group)) + 
				geom_bar(position=dodge, stat="identity", alpha=0.75) + 
				geom_bar(position=dodge, stat="identity", alpha=0.75, colour="black", show.legend=FALSE, size=0.25) +
				geom_errorbar(limits, position=dodge, size=0.25, width=0.25) +
				labs(y=paste(ylab), x='') +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
				theme(legend.position="top", legend.title=element_blank())
		
		if (scale == 'sqrt') {
			# obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
			obj1 <- obj1 + scale_y_sqrt(
					breaks = trans_breaks("sqrt", function(x) x^2),
					labels = trans_format("sqrt", math_format(.x^2)))
		}
		
		# To be revised
		if (scale == 'log') {
			obj1 <- obj1 + scale_y_log10(
					breaks = trans_breaks("log10", function(x) 10^x),
					labels = trans_format("log10", math_format(10^.x)))
		}
		if (scale == 'boxcox') {
			obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
		}
		
	} else {
		grp2 <- factor(paste(strata, df[, strata]), levels=paste(strata, levels(df[, strata])))
		df2 <- data.frame(Group=grp, Strata=grp2, t(prop))
		df2 <- melt(df2)
		colnames(df2) <- c('Group', 'Strata', 'Taxa', 'Value')
		
		# Could be revised
		temp1 <- aggregate(Value ~ Group + Strata + Taxa, df2, mean)

		if (error == 'se') {
			temp2 <- aggregate(Value ~ Group + Strata +Taxa, df2, function(x) {
						#sd(x) / sqrt(length(x))
						sd(sapply(1:50, function(i) mean(sample(x, repl=TRUE))))
					})
		}
		if (error == 'ci') {
			temp2 <- aggregate(Value ~ Group + Strata + Taxa, df2, function(x) {
						#sd(x) / sqrt(length(x))
						2 * sd(sapply(1:50, function (i) mean(sample(x, repl=TRUE))))
					})
		}
		
		
		df2 <- cbind(temp1, temp2[, 4])
		colnames(df2) <- c('Group', 'Strata', 'Taxa', 'Mean', 'SE')
		
		limits <- aes(ymax = Mean + SE, ymin = ifelse(Mean - SE > 0, Mean - SE, 0))
		dodge <- position_dodge(width=0.90)
		
		obj1 <- ggplot(df2, aes(x=Taxa, y=Mean, fill=Group)) + 
				geom_bar(position=dodge, stat="identity", alpha=0.75) + 
				geom_bar(position=dodge, stat="identity", alpha=0.75, colour="black", show.legend=FALSE, size=0.25) +
				geom_errorbar(limits, position=dodge, size=0.25, width=0.25) +
				facet_wrap(~Strata, ncol=1) +
				labs(y=paste(ylab), x='') +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
				theme(legend.position="top", legend.title=element_blank())
		if (scale == 'sqrt') {
			# obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
			obj1 <- obj1 + scale_y_sqrt(
					breaks = trans_breaks("sqrt", function(x) x^2),
					labels = trans_format("sqrt", math_format(.x^2)))
		}
		if (scale == 'log') {
			obj1 <- obj1 + scale_y_log10(
					breaks = trans_breaks("log10", function(x) 10^x),
					labels = trans_format("log10", math_format(10^.x)))
		}
		if (scale == 'boxcox') {
			obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
		}
		
	}
	return(obj1)
	
}

# Rev: 2018_03_08 Add geom_point and jsize
# Rev: 2019_01_03 handling zeros for log scales
taxa_boxplot_aggregate <- function (prop, df, grp.name, strata=NULL, scale='none', xsize=10, jsize=NULL, ylab='Proportion') {
	
	
	grp <- factor(df[, grp.name])
	
	if (scale == 'log') {
		# Half of the non-zero proportion
		prop <- t(apply(prop, 1, function (x) {
							temp <- min(x[x != 0])
							x[x == 0] <- temp / 2
							x
						}))
	}
	
	if (is.null(jsize)) {
		nT <- nrow(prop) 
		jsize <- 2 / (nT %/% 20 + 1)
	}

	if (is.null(strata)) {
		df2 <- data.frame(Group=grp, t(prop))
		df2 <- melt(df2)
		colnames(df2) <- c('Group', 'Taxa', 'Value')
		
		dodge <- position_dodge(width=0.88)
		
		obj1 <- ggplot(df2, aes(x=Taxa, y=Value, fill=Group)) + 
				geom_boxplot(position=dodge, alpha = ifelse(jsize == 0, 0.75, 0.75), outlier.alpha=ifelse(jsize == 0, 0.5, 0),  lwd=0.25, fatten=1)
		if (jsize != 0) {
			obj1 <- obj1 + geom_point(position=position_jitterdodge(dodge.width=0.88), size=jsize, alpha=0.3)
		}
		obj1 <- obj1 +
				labs(y=paste(ylab), x='') +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
				theme(legend.position="top", legend.title=element_blank())
		if (scale == 'sqrt') {
			# obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
			obj1 <- obj1 + scale_y_sqrt(
					breaks = trans_breaks("sqrt", function(x) x^2),
					labels = trans_format("sqrt", math_format(.x^2)))
		}
		if (scale == 'log') {
			obj1 <- obj1 + scale_y_log10(
					breaks = trans_breaks("log10", function(x) 10^x),
					labels = trans_format("log10", math_format(10^.x)))
		}
		if (scale == 'boxcox') {
			obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
		}
		
	} else {
		grp2 <- factor(paste(strata, df[, strata]), levels=paste(strata, levels(df[, strata])))
		df2 <- data.frame(Group=grp, Strata=grp2, t(prop))
		df2 <- melt(df2)
		colnames(df2) <- c('Group', 'Strata', 'Taxa', 'Value')
		
		dodge <- position_dodge(width=0.88)
		
		obj1 <- ggplot(df2, aes(x=Taxa, y=Value, fill=Group)) + 
				geom_boxplot(position=dodge, alpha = ifelse(jsize == 0, 0.75, 0.75),  outlier.alpha=ifelse(jsize == 0, 0.5, 0),  lwd=0.25, fatten=1) 
		if (jsize != 0) {
			obj1 <- obj1 + geom_point(position=position_jitterdodge(dodge.width=0.88), size=jsize, alpha=0.2)
		}
		obj1 <- obj1 +
			#	geom_point(position=position_jitterdodge(dodge.width=0.95), size=jsize, alpha=0.3) +
				facet_wrap(~Strata, ncol=1) +
				labs(y=paste(ylab), x='') +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
				theme(legend.position="top", legend.title=element_blank())
		if (scale == 'sqrt') {
			# obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
			obj1 <- obj1 + scale_y_sqrt(
					breaks = trans_breaks("sqrt", function(x) x^2),
					labels = trans_format("sqrt", math_format(.x^2)))
		}
		# To be revised
		if (scale == 'log') {
			obj1 <- obj1 + scale_y_log10(
					breaks = trans_breaks("log10", function(x) 10^x),
					labels = trans_format("log10", math_format(10^.x)))
		}
		if (scale == 'boxcox') {
			obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
		}
		
	}
	return(obj1)
	
}

generate_taxa_biplot <- function (data.obj, taxa, trans='sqrt', grp.name, ann='', ...) {
	
    grp <- data.obj$meta.dat[, grp.name]

	prop <- NULL
	for (LOI2 in names(data.obj$abund.list)) {
		ct <- data.obj$abund.list[[LOI2]]
		ct <- ct + 1
		prop0 <- t(t(ct) / colSums(ct))
		prop <- rbind(prop, prop0[intersect(rownames(prop0), taxa), , drop=FALSE])	
	}
	colnames(prop) <- colnames(prop0)
	if (nrow(prop) != length(taxa)) {
		warnings('Some taxa not found in abundance lists! Please check the names!\n')
	}
	
	if (trans == 'normal') 	prop <- t(apply(prop, 1, function(x) qqnorm(x, plot=F)$x))
	if (trans == 'log') prop <- log(prop)
	if (trans == 'sqrt') prop <- sqrt(prop)
	if (trans == 'rank') prop <- t(apply(prop, 1, rank))
	
	wine.pca <- prcomp(t(prop), scale. = TRUE)
	g <- ggbiplot(wine.pca, obs.scale = 1, var.scale = 1, 
			groups = grp, ellipse = TRUE, circle = FALSE, ...) 
	g <- g + scale_color_discrete(name = '') + theme_bw()
	g <- g + theme(legend.direction = 'horizontal', legend.position = 'top') 
	pdf(paste0("Taxa_Biplot_", ann, ".pdf"), height=6, width=6)
	print(g)
	dev.off()
}

# Rev: 2017_05_19: new formula for width
#      Add presence and absence bar
# Rev: 2018_07_16
generate_taxa_barplot_aggregate <- function (data.obj,  grp.name, strata=NULL, scale='sqrt', pseudo.ct = 0.5,
		taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
		taxa.name='All',  prev=0.1, minp=0.002, ann='All', wids=NULL, heis=NULL, xsize=c(14, 12, 10, 9, 8)) {
	# To be completed
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	for (i in 1:length(taxa.levels)) {
		LOI <- taxa.levels[i]
		cat(LOI, "\n")
		
		ct <- data.obj$abund.list[[LOI]]
		prop <- t(t(ct / colSums(ct)))
		ind <- matrixStats::rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop)
		
		if (scale == 'log') {
			ct <- ct + pseudo.ct
			prop <- t(t(ct) / colSums(ct))
		}
		
		if ('All' %in% taxa.name) {
			prop <- prop[ind, , drop=FALSE]
		} else {
			# Rev: 2018_07_16
			prop <- prop[rownames(prop) %in% taxa.name, , drop=FALSE]
		}
			
		# decreasing
		prop <- prop[rev(order(rowMeans(prop))), ]
		
		if (is.null(wids) | is.null(heis)) {
			wid <- 7 * ifelse(nrow(prop) / 30 < 1, 1, nrow(prop) / 30)
			wid <- sqrt(nlevels(grp) / 2) * wid
	        hei <- 7
		} else {
			wid <- wids[i]
			hei <- heis[i]
		}
			
		pdf(paste("Taxa_Barplot_Aggregate", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)		
		obj1 <- taxa_barplot_aggregate(prop, df, grp.name, strata, scale, xsize[i]) 
		print(obj1)
		dev.off()		
	}
}

# Rev: 2017_05_19: new formula for width
# Rev: 2018_03_08  jsize
# Rev: 2018_07_16 
generate_taxa_boxplot_aggregate <- function (data.obj,  grp.name, strata=NULL, scale='sqrt', pseudo.ct = 0.5,
		taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
		taxa.name='All',  prev=0.1, minp=0.002, ann='All', wids=NULL, heis=NULL, xsize=c(14, 12, 10, 9, 8, 6, 4), jsize = 0) {
	# To be completed
	df <- data.obj$meta.dat
	grp <- factor(df[, grp.name])
	
	for (i in 1:length(taxa.levels)) {
		LOI <- taxa.levels[i]
		cat(LOI, "\n")
		
		ct <- data.obj$abund.list[[LOI]]
		prop <- t(t(ct) / colSums(ct))
		ind <- matrixStats::rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop)
		
		if (scale == 'log') {
			ct <- ct + pseudo.ct
			prop <- t(t(ct / colSums(ct)))
		}
		
		if ('All' %in% taxa.name) {
			prop <- prop[ind, , drop=FALSE]
		} else {
			# Rev: 2018_07_16
			prop <- prop[rownames(prop) %in% taxa.name, , drop=FALSE]
		}
		
		# decreasing
		prop <- prop[rev(order(rowMeans(prop))), ]
		
		if (is.null(wids) | is.null(heis)) {
			wid <- 7 * ifelse(nrow(prop) / 30 < 1, 1, nrow(prop) / 30)
			wid <- sqrt(nlevels(grp) / 2) * wid
		    
			hei <- 7
		} else {
			wid <- wids[i]
			hei <- heis[i]
		}
		
		pdf(paste("Taxa_Boxplot_Aggregate", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)		
		obj1 <- taxa_boxplot_aggregate (prop, df, grp.name, strata, scale, xsize[i], jsize=jsize) 
		print(obj1)
		dev.off()		
	}
}

# Add combined barplot with error bar and presence/absence bar (currently presence/absence bar is in generate_taxa_boxplot
# Rev: 2018_01_15 Add qt.outlier=0.97
# Rev: 2018_03_06 coord_cartesian
# Rev: 2018_05_31 Add x axis label control
# Rev: 2018_07_16
# Rev: 2019_01_03 add sqrtP and selecting taxa before adding pseudo-count
generate_taxa_barplot <- function (data.obj,  grp.name, strata=NULL, scale='P', pseudo.ct = 0.5,
		taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
		taxa.name='All', rm.outlier=T, qt.outlier=0.97, prev=0.1, minp=0.002, ann='All', x.axis.lab=FALSE, x.axis.size = 6) {
	# To be completed
	df <- data.obj$meta.dat
	grp <- factor(df[, grp.name])
	
	for (LOI in taxa.levels) {
		if (LOI == 'All') {
			if ('All' %in% taxa.name) stop('Taxa names are not allowed to be all when taxa level is also all!\n')
			headnames <- NULL
			prop <- NULL
			for (LOI2 in names(data.obj$abund.list)) {
				prop2 <- data.obj$abund.list[[LOI2]]
				taxa.name2 <- taxa.name[taxa.name %in% rownames(prop2)]
				if (length(taxa.name2) != 0) {
					if (scale == 'logP') {
						prop2 <- prop2 + pseudo.ct
					} 
					prop2 <- t(t(prop2) / colSums(prop2))			
					headnames2 <- sapply(strsplit(rownames(prop2), ";"), paste, collapse="\n")
					names(headnames2) <- rownames(prop2)
					prop <- rbind(prop, prop2[taxa.name2, , drop=FALSE])
					headnames <- c(headnames, headnames2)
				}
			}
			
		} else {
			cat(LOI, "\n")
			ct <- data.obj$abund.list[[LOI]]
			prop <- t(t(ct) / colSums(ct))
			ind <- matrixStats::rowMaxs(prop) > minp & rowSums(prop != 0) > prev * ncol(prop)
			
			if (scale == 'logP') {
				ct <- ct + pseudo.ct
				prop <- t(t(ct) / colSums(ct))
			} 
			
			headnames <- sapply(strsplit(rownames(prop), ";"), paste, collapse="\n")
			names(headnames) <- rownames(prop)
			
			if ('All' %in% taxa.name) {
				prop <- prop[ind, , drop=FALSE]
			} else {
				# Rev: 2018_07_16
				prop <- prop[rownames(prop) %in% taxa.name, , drop=FALSE]
			}
			
		}

		if (scale == 'logP') {
			prop <- log10(prop)
		}
		
		if (scale == 'sqrtP') {
			prop <- sqrt(prop)
		}
		
		
		hei <- 5
		if (is.null(strata)) {
			wid <- 5
		} else {
			wid <- 4 * nlevels(df[, strata])
		}
		
		if (scale == 'P') {
			ylab <- 'Proportion'
		} 
		
		if (scale == 'none') {
			ylab <- 'Value'
		} 
		
		if (scale == 'logP') {
			ylab <- 'log10(Proportion)'
		}
		
		if (scale == 'sqrtP') {
			ylab <- 'sqrt(Proportion)'
		}
		
		pdf(paste("Taxa_Barplot", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)
		if (is.null(strata)) {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				taxon.abund2 <- taxon.abund
				if (scale == 'P') {
					if (rm.outlier == T) {
						ylims <- c(0, quantile(taxon.abund, qt.outlier) * 1.25)
						taxon.abund[taxon.abund > quantile(taxon.abund, qt.outlier) * 1.25] <- quantile(taxon.abund, qt.outlier) * 1.25
					} else {
						ylims <- c(0, max(taxon.abund))
					}
					
				} else {
					ylims <- range(taxon.abund)
				}
				# df2 <- data.frame(x=factor(1:length(taxon.abund)), Value=taxon.abund, Group=grp)
		        df2 <- data.frame(x=colnames(prop), Value=taxon.abund, Group=grp)
				rownames(df2) <- colnames(prop)
				
				dodge <- position_dodge(width=0.99)
				obj <- ggplot(df2, aes(x=x, y=Value, fill=Group)) +
						geom_bar(position=dodge, stat='identity', alpha=0.75) + 
						facet_grid(. ~ Group, scales='free_x', space="free") +
						labs(y=ylab, title=headnames[taxon])
				
				if (rm.outlier) {
					obj <- obj + coord_cartesian(ylim = ylims)
				}
				
				if (!x.axis.lab) {
					obj <- obj + theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) 
				} else {
					obj <- obj + theme(axis.text.x = element_text(size = x.axis.size, angle = 90, vjust = 0.5, hjust=1))
				}
				
				obj <- obj  + 
						xlab('') +

						theme(legend.position="none")
				mean_df <- data.frame(Group=levels(grp), yint = tapply(taxon.abund2, grp, mean))
				obj <- obj + geom_hline(data = mean_df, aes(yintercept = yint))		
				mean_df <- data.frame(Group=levels(grp), yint = tapply(taxon.abund2, grp, median))
				obj <- obj + geom_hline(data = mean_df, aes(yintercept = yint), linetype=2)
				
				print(obj)
			}	
		} else {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				taxon.abund2 <- taxon.abund
				if (scale == 'P') {
					if (rm.outlier == T) {
						ylims <- c(0, quantile(taxon.abund, qt.outlier) * 1.25)
						taxon.abund[taxon.abund > quantile(taxon.abund, qt.outlier) * 1.25] <- quantile(taxon.abund, qt.outlier) * 1.25
					} else {
						ylims <- c(0, max(taxon.abund))
					}
					
				} else {
					ylims <- range(taxon.abund)
				}
				
				grp2 <- df[, strata]
				obj.list <- list()
				for (level in levels(grp2)) {
					ind <- grp2 %in% level
					# df2 <- data.frame(x=factor(1:length(taxon.abund[ind])), Value=taxon.abund[ind], Group=grp[ind])	
			        df2 <- data.frame(x=colnames(prop)[ind], Value=taxon.abund[ind], Group=grp[ind])	
					rownames(df2) <- colnames(prop)[ind]
					
					dodge <- position_dodge(width=0.99)
					obj0 <- ggplot(df2, aes(x=x, y=Value, fill=Group)) +
							geom_bar(position=dodge, stat='identity', alpha=0.75) + 
							facet_grid(. ~ Group, scales='free_x', space="free") +
							labs(y=ylab, title=headnames[taxon])
					
					if (rm.outlier) {
						obj0 <- obj0 + coord_cartesian(ylim = ylims)
					}
					
					if (!x.axis.lab) {
						obj0 <- obj0 + theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) 
					} else {
						obj0 <- obj0 + theme(axis.text.x = element_text(size = x.axis.size, angle = 90, vjust = 0.5, hjust=1))
					}
					
					
					obj0 <- obj0  + 
							xlab(paste(strata, level)) +
							theme(legend.position="none")
					mean_df <- data.frame(Group=levels(grp[ind]), yint = tapply(taxon.abund2[ind], grp[ind], mean))
					obj0 <- obj0 + geom_hline(data = mean_df, aes(yintercept = yint))		
					mean_df <- data.frame(Group=levels(grp[ind]), yint = tapply(taxon.abund2[ind], grp[ind], median))
					obj0 <- obj0 + geom_hline(data = mean_df, aes(yintercept = yint), linetype=2)
					obj.list[[level]] <- obj0		
				}
				multiplot(plotlist=obj.list, cols=nlevels(grp2))
			}	
			
		}
		dev.off()
	}		
}


twopart.test <- function(x1, x2, zero.p=0.2) {
	
	n1 <- length(x1)
	n2 <- length(x2)
	p1 <- mean(x1 != 0)
	p2 <- mean(x2 != 0)
	m1 <- sum(x1 != 0)
	m2 <- sum(x2 != 0)
	p12 <- (m1 + m2) / (n1 + n2)
	q12 <- 1 - p12
	
	if (q12 >= zero.p) {
		Z <- (abs(p1 - p2) - (1/(2*n1) + 1/(2*n2))) / sqrt(p12 * q12 * (1/n1 + 1/n2))
		x1 <- x1[x1!=0]
		x2 <- x2[x2!=0]
		R1 <- sum(rank(c(x1, x2))[1:length(x1)])
		ti <- as.vector(table(c(x1, x2)))
		W <- (abs(R1 - m1*(m1+m2+1)/2) - 1/2) / sqrt((m1*m2/12)*(m1+m2+1-sum(ti*(ti^2-1))/(m1+m2)/(m1+m2-1)))
		X2 <- Z^2 + W^2
		res <- list()
		res$stat <- X2
		res$p.value <- 1 - pchisq(X2, 2)
		res$Z <- Z
		res$W <- W
		res$test <- 'TwoPart'
	} else {
		res <- wilcox.test(x1, x2)
		res$test <- 'Wilcox'
	}
	res
}

getPermuteMatrix <- function (perm, N, strata = NULL) 
{
	if (length(perm) == 1) {
		perm <- how(nperm = perm)
	}
	if (!missing(strata) && !is.null(strata)) {
		if (inherits(perm, "how") && is.null(getBlocks(perm))) 
			setBlocks(perm) <- strata
	}
	if (inherits(perm, "how")) 
		perm <- shuffleSet(N, control = perm)
	if (is.null(attr(perm, "control"))) 
		attr(perm, "control") <- structure(list(within = list(type = "supplied matrix"), 
						nperm = nrow(perm)), class = "how")
	perm
}

perm_fdr_adj <- function (F0, Fp) {
	ord <- order(F0, decreasing = T)
	F0 <- F0[ord]
	perm.no <- ncol(Fp)
	Fp <- as.vector(Fp)
	Fp <- Fp[!is.na(Fp)]
	Fp <- sort(c(Fp, F0), decreasing = F)
	n <- length(Fp)
	m <- length(F0)
	FPN <- (n + 1) - match(F0, Fp) - 1:m
	p.adj.fdr <- FPN / perm.no / (1:m)
#		p.adj.fdr <- sapply(F0, function(x) sum(Fp >= 
#									x, na.rm=TRUE) / perm.no)/(1:length(F0))
	p.adj.fdr <- pmin(1, rev(cummin(rev(p.adj.fdr))))[order(ord)]
}

perm_fwer_adj <- function (F0, Fp) {
	ord <- order(F0, decreasing = T)
	F0 <- F0[ord]
	col.max <- colMaxs(Fp, na.rm=TRUE)
	p.adj.fwer <- sapply(F0, function(x) mean(col.max >= x))[order(ord)]
}

# New: 2018_06_20
# Westfall young
perm_fwer_adj2 <- function (F0, Fp) {
	ord <- order(F0, decreasing = T)
	m <- length(F0)
	F0 <- F0[ord]
	Fp <- Fp[ord, , drop=FALSE]
	col.max <- Fp[m, ]
	p.adj.fwer <- sapply(m:1, function(i) {
				x <- F0[i]
				y <- Fp[i, ]
				col.max <<- ifelse(y > col.max, y, col.max)
				mean(col.max >= x)
			})
	p.adj.fwer <- rev(p.adj.fwer)
	p.adj.fwer <- pmin(1, rev(cummin(rev(p.adj.fwer))))[order(ord)]
}

# Need to be further comprehensively tested
# a. Speed up, finished!
# b. Permutation method （response, covariate or residual permutation), residual will be default!
# c. Revise - add permutation stratified by subject, finished!
# Rev: 2017_02_02 permutation-based FDR control
# Rev: 2017_02_13 Add LMM-based permutation for block.perm=TRUE (type I error is controled, but power study hasn't been comprehensively studied)
# Rev: 2017_02_24 allow 'adj.name' to contain multiple covariates
# Still need to address NA's, currently simply remove NA's
# Rev: 2018_03_18 add winsorization support for permutation test - Only one side winsorization
# Rev: 2018_03_18 add effect size (partial R2) and Direction
# Rev: 2018_06_20 add Westfall young FWER adjustment
permute_differential_analysis <- function (meta.dat, prop, grp.name, adj.name=NULL, strata=NULL, 
		block.perm=FALSE, sqrt.trans=TRUE, resid.perm=TRUE, perm.no=999, winsor=FALSE, winsor.qt = 0.97, trace = FALSE) {
	# Square root transformation
	# User should take care of the normalization, transformation and addressing outliers
	opt.old <- options('contrasts')
	options(contrasts = c("contr.treatment", "contr.poly"))
	meta.dat <- droplevels(meta.dat)
	if (sqrt.trans) {
		Y <- sqrt(prop)
	} else {
		Y <- prop
	}
	row.names <- rownames(Y)
	
	if (winsor == TRUE) {
		# Only one side winsorization - be careful for different type of transformation
		Y <- t(apply(Y, 1, function (x) {
					cutoff <- quantile(x, winsor.qt)
					x[x >= cutoff] <- cutoff
					x
				}))
	}
	
	if (!is.null(strata)) {
		strata <- factor(strata)
	}
	
	# Prepare model matrix
	n <- ncol(prop)
	I <- diag(n)
	if (is.null(adj.name)) {
		M0 <- model.matrix(~ 1, meta.dat)
	} else {
		df0 <- meta.dat[, c(adj.name), drop=F]
		M0 <- model.matrix( ~., df0)
	}
	
#	P0 <- M0 %*% solve(t(M0) %*% M0) %*% t(M0)
	
	df1 <- meta.dat[, c(adj.name, grp.name), drop=F]
	M1 <- model.matrix( ~., df1)
	
	# Make sure it is contr.tretament
	M2 <- model.matrix (~.,  meta.dat[, c(grp.name), drop=F])
	
	# Remove the effect of confounder
	M2 <- as.matrix(resid(lm(M2[, -1] ~ M0 - 1)))
	
	# QR decompostion
	qrX0 <- qr(M0, tol = 1e-07)
	Q0 <- qr.Q(qrX0)
	Q0 <- Q0[, 1:qrX0$rank, drop=FALSE]
	
	qrX1 <- qr(M1, tol = 1e-07)
	Q1 <- qr.Q(qrX1)
	Q1 <- Q1[, 1:qrX1$rank, drop=FALSE]
	
	# Got residual
	if (resid.perm) {
		if (!block.perm) {
			# Permute the residual
			if (is.null(adj.name)) {
				Y <- t(resid(lm(as.formula(paste('t(Y) ~ 1')), meta.dat)))
			} else {
				Y <- t(resid(lm(as.formula(paste('t(Y) ~ ', paste(adj.name, collapse='+'))), meta.dat)))
			}
			
		} else {
			
			if (is.null(strata)) {
				stop('Block permutation requires strata!\n')
			} else {
				Y.r <- matrix(NA, nrow(Y), nlevels(strata))
				Y.e <- Y
				if(trace) cat('Fitting linear mixed effects model ...\n')
				for (j in 1:nrow(Y)) {
					# Linear mixed effects model
					yy <- Y[j, ]
					meta.dat$yy <- yy
					meta.dat$strata <- strata
					if (is.null(adj.name)) {
						# obj <- lme(as.formula(paste('yy ~ 1')), random =~ 1 |  strata, data=meta.dat, method='ML')
					     obj <- lmer(as.formula(paste('yy ~ 1 + (1 | strata)')),  data=meta.dat, REML=FALSE)
						# The order is the same as the levels
						# Y.r[j, ] <- random.effects(obj)[, 1]
					    Y.r[j, ] <- random.effects(obj)[[1]][, 1]
						Y.e[j, ] <- resid(obj)
					} else {
						# obj <- lme(as.formula(paste('yy ~ ', paste(adj.name, collapse='+'))), random =~ 1 |  strata, data=meta.dat, method='ML')
					    obj <- lmer(as.formula(paste('yy ~ ', paste(adj.name, collapse='+'), ' + (1 | strata)')),data=meta.dat, REML=FALSE)
						# Y.r[j, ] <- random.effects(obj)[, 1]
						Y.r[j, ] <- random.effects(obj)[[1]][, 1]
						Y.e[j, ] <- resid(obj)
					}
				}
				Y <- Y.r[, as.numeric(strata)] + Y.e
	#			Y <- Y - rowMeans(Y)
	#		    Y <- Y.e
			}
		}
	}
	
	TSS <- rowSums(Y^2)
	MSS1 <- rowSums((Y %*% Q1)^2)
	MSS0 <- rowSums((Y %*% Q0)^2)  # Not necessary, it's zero
	F0 <- (MSS1 - MSS0) /  (TSS - MSS1) 
	
	R2 <- (MSS1 - MSS0) / (TSS - MSS0)
	coefficients <- t(solve(t(M2) %*% M2) %*% t(M2) %*% t(Y))
	
#	P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
#	F0 <- diag(Y %*% (P1 - P0) %*% t(Y)) / diag(Y %*% (I - P1) %*% t(Y))
#	df3 <- df1
	
	if (block.perm == FALSE) {
		perm.ind <- vegan:::getPermuteMatrix(perm.no, n, strata = strata)
		perm.no <- nrow(perm.ind)
	}
	
	if(trace) cat('Permutation test ....\n')
	Fp <- sapply(1:perm.no, function(i) {
				if (i %% 100 == 0)  if(trace) cat('.')
				if (block.perm == FALSE) {
					Yp <- Y[, perm.ind[i, ]]
				} else {
					# Double permutation
					strata.p <- factor(strata, levels=sample(levels(strata)))
				    Yp <- Y.r[, as.numeric(strata.p)] + Y.e[, sample(ncol(Y))]
#				    Yp <- Y.e[, sample(ncol(Y))]
#					Yp <- Yp - rowMeans(Yp)
				}
				
#				df3[, grp.name] <- sample(df1[, grp.name])
#				M1 <- model.matrix( ~., df3)
#				P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
				MSS1p <- rowSums((Yp %*% Q1)^2)
				MSS0p <- rowSums((Yp %*% Q0)^2)    
			    if (block.perm == FALSE) {
					TSSp <- TSS
				} else {
					TSSp <- rowSums(Yp^2)
				}
				(MSS1p - MSS0p) /  (TSSp - MSS1p) 
			})
	
	
	if (mean(is.na(F0)) >= 0.1) {
		warning('More than 10% observed F stats have NA! Please check! \n')
	}
	
	if (mean(is.na(Fp)) >= 0.1) {
		warning('More than 10% permuted F stats have NA! Please check! \n')
	}
	
	na.ind <- is.na(F0)
	F0 <- F0[!na.ind]
	Fp <- Fp[!na.ind, ]
	
	p.raw <- cbind(Fp >= F0, 1)
	p.raw <- rowMeans(p.raw)
#	p.raw[is.na(p.raw)] <- 1   
	
	p.adj.fdr <- perm_fdr_adj(F0, Fp)
	p.adj.fwer <- perm_fwer_adj2(F0, Fp)
	
	# Pad back the NA values
	pad <- function (vec, ind) {
		vec0 <- numeric(length(ind))
		vec0[!ind] <- vec
		vec0[ind] <- NA
		vec0
	}
	
	options(contrasts=opt.old$contrast)
	
	F0 <- pad(F0, na.ind)
	p.raw <- pad(p.raw, na.ind)
	p.adj.fdr <- pad(p.adj.fdr, na.ind)
	p.adj.fwer <- pad(p.adj.fwer, na.ind)
	
	names(F0) <- names(p.raw) <- names(p.adj.fdr) <- names(p.adj.fwer) <- row.names
	return(list(F.stat=F0, R2=R2, coefficients=coefficients, p.raw=p.raw, p.adj.fdr=p.adj.fdr, p.adj.fwer=p.adj.fwer))
#	return(p.raw)
}
#


# New: 2017_08_17 Add a new variant of PERMANOVA with matrix decomposition
# Partially validated
permanova2 <- function (meta.dat, D, grp.name, adj.name=NULL, strata=NULL, 
		block.perm=FALSE, resid.perm=TRUE, perm.no=999, eig = c('All', 'Positive', 'Top'), var.exp = 0.90) {
	
	eig <- match.arg(eig)
	# Square root transformation
	# User should take care of the normalization, transformation and addressing outliers

	D <- -as.matrix(D)^2 / 2
	n <- nrow(D)
	G <- mean(D) + D - rowMeans(D) - matrix(rep(1, n), ncol=1) %*% colMeans(D)
	eig.obj <- eigen(G)
	if (eig == 'All') {
		Y <- eig.obj$vectors
		lambda <- eig.obj$values
		ind <- abs(lambda) > 1e-6
		lambda <- lambda[ind]
		Y <- Y[, ind, drop = FALSE]
	}
	
	if (eig == 'Positive') {
		lambda <- eig.obj$values
		ind <- lambda > 1e-6
		lambda <- lambda[ind]
		Y <- Y[, ind, drop = FALSE]
	}

	if (eig == 'Top') {
		Y <- eig.obj$vectors
		lambda <- eig.obj$values
		ind <- lambda > 1e-6
		lambda <- lambda[ind]
		Y <- Y[, ind, drop = FALSE]
	
		ind <- 1:(which(cumsum(lambda / sum(lambda)) > var.exp)[1])
		lambda <- lambda[ind]
		Y <- Y[, ind, drop = FALSE]
	}
	
	
    Y <- t(Y)
	
	if (!is.null(strata)) {
		strata <- factor(strata)
	}
	
	# Prepare model matrix
	n <- ncol(Y)
	I <- diag(n)
	if (is.null(adj.name)) {
		M0 <- model.matrix(~ 1, meta.dat)
	} else {
		df0 <- meta.dat[, c(adj.name), drop=F]
		M0 <- model.matrix( ~., df0)
	}
	
#	P0 <- M0 %*% solve(t(M0) %*% M0) %*% t(M0)
	
	df1 <- meta.dat[, c(adj.name, grp.name), drop=F]
	M1 <- model.matrix( ~., df1)
	
	# QR decompostion
	qrX0 <- qr(M0, tol = 1e-07)
	Q0 <- qr.Q(qrX0)
	Q0 <- Q0[, 1:qrX0$rank, drop=FALSE]
	
	qrX1 <- qr(M1, tol = 1e-07)
	Q1 <- qr.Q(qrX1)
	Q1 <- Q1[, 1:qrX1$rank, drop=FALSE]
	
	# Got residual
	if (resid.perm) {
		if (!block.perm) {
			# Permute the residual
			if (is.null(adj.name)) {
				Y <- t(resid(lm(as.formula(paste('t(Y) ~ 1')), meta.dat)))
			} else {
				Y <- t(resid(lm(as.formula(paste('t(Y) ~ ', paste(adj.name, collapse='+'))), meta.dat)))
			}
			
		} else {
			
			if (is.null(strata)) {
				stop('Block permutation requires strata!\n')
			} else {
				Y.r <- matrix(NA, nrow(Y), nlevels(strata))
				Y.e <- Y
#				cat('Fitting linear mixed effects model ...\n')
				for (j in 1:nrow(Y)) {
					# Linear mixed effects model
					yy <- Y[j, ]
					meta.dat$yy <- yy
					meta.dat$strata <- strata
					if (is.null(adj.name)) {
						
						obj <- lmer(as.formula(paste('yy ~ 1 + (1 | strata)')), data=meta.dat)
						# The order is the same as the levels
						Y.r[j, ] <- random.effects(obj)[[1]][, 1]
						Y.e[j, ] <- resid(obj)
					} else {
						obj <- lmer(as.formula(paste('yy ~ ', paste(adj.name, collapse='+'), '+ (1 | strata)')), data=meta.dat)
						Y.r[j, ] <- random.effects(obj)[[1]][, 1]
						Y.e[j, ] <- resid(obj)
					}
				}
				Y <- Y.r[, as.numeric(strata)] + Y.e
				#			Y <- Y - rowMeans(Y)
				#		    Y <- Y.e
			}
		}
	}
	
	TSS <- rowSums(Y^2)
	MSS1 <- rowSums((Y %*% Q1)^2)
	MSS0 <- rowSums((Y %*% Q0)^2)  # Not necessary, it's zero
	F0 <- sum(lambda * (MSS1 - MSS0)) /  sum(lambda * (TSS - MSS1))
	
#	P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
#	F0 <- diag(Y %*% (P1 - P0) %*% t(Y)) / diag(Y %*% (I - P1) %*% t(Y))
#	df3 <- df1
	
	if (block.perm == FALSE) {
		perm.ind <- vegan:::getPermuteMatrix(perm.no, n, strata = strata)
		perm.no <- nrow(perm.ind)
	}
	
#	cat('Permutation test ....\n')
	Fp <- sapply(1:perm.no, function(i) {
#				if (i %% 100 == 0) cat('.')
				if (block.perm == FALSE) {
					Yp <- Y[, perm.ind[i, ]]
				} else {
					# Double permutation
					strata.p <- factor(strata, levels=sample(levels(strata)))
					Yp <- Y.r[, as.numeric(strata.p)] + Y.e[, sample(ncol(Y))]
#				    Yp <- Y.e[, sample(ncol(Y))]
#					Yp <- Yp - rowMeans(Yp)
				}
				
#				df3[, grp.name] <- sample(df1[, grp.name])
#				M1 <- model.matrix( ~., df3)
#				P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
				MSS1p <- rowSums((Yp %*% Q1)^2)
				MSS0p <- rowSums((Yp %*% Q0)^2)    
				if (block.perm == FALSE) {
					TSSp <- TSS
				} else {
					TSSp <- rowSums(Yp^2)
				}
				sum(lambda * (MSS1p - MSS0p)) / sum(lambda * (TSSp - MSS1p)) 
			})
	
	p.value <- mean(c(Fp >= F0, 1))
  
	return(list(f0 = F0, f.perms = Fp, p.value = p.value))
#	return(p.raw)
}




# This function for nonparametric/permutaiton method
# Rev: 2017_02_16  Add normalization method; Add transformation; Remove rarefaction (only output warnings);
# Rev: 2017_10_30  Support filtering based on coefficient of variation
# Rev: 2018_01_30 Add ct.min and handle NA size factor
# Rev: 2018_03_15 Add winsorization support for permutatin test - default is FALSE
# Rev: 2018_10_03 TSS size factor calculation is based on specific taxa level
perform_differential_analysis <- function (data.obj, grp.name, adj.name=NULL, subject=NULL, 
		taxa.levels=c('Phylum', 'Order', 'Class', 'Family', 'Genus', 'Species'),
		method='perm', block.perm=FALSE, perm.no=999,
		norm='GMPR', norm.level='Species', intersect.no=4, ct.min = 2,
		transform='sqrt',
		winsor = FALSE, winsor.qt = 0.97,  
		prev=0.1, minp=0.002, medianp=NULL, cv=NULL,
		mt.method='fdr', cutoff=0.15, 
		ann='', seed=123, ...) {
	# To be completed
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	ind <- !is.na(grp)
	data.obj <- subset_data(data.obj, ind)
	grp <- grp[ind]
	df <- df[ind, ]
	
	if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
		data.obj$abund.list[['Species']] <- data.obj$otu.tab
		rownames(data.obj$abund.list[['Species']]) <- paste0("OTU", rownames(data.obj$otu.tab), ":", 
				data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
	}

	# Test for sequence-depth confounding
	dep <- colSums(data.obj$otu.tab)
	
	diff.seq.p <- summary(aov(dep ~ grp))[[1]][1, 'Pr(>F)']
	if (!is.na(diff.seq.p) & diff.seq.p <= 0.05) {
		warning(paste0(
		 '\nSignificant sequencing depth confounding with the variable of interest!\n',
		 'For nonparametric test/permutaiton test, there may be potentially many false postives (for those less prevalent taxa)!\n',
		 'Consider performing rarefaction first! (However, rarefaction will not completly solve the problem.)\n',
		 'May also try count-based models, which might have better false postive controls!\n'))
#		cat("For parametric test with sequence depth adjustment (DESeq2), please be cautious about the results!\n There may be potential residual sequence depth confounding!\n")
		# potential revising
	}
  
	# Calculate size.factor
	if (norm == 'Precalculated') {
		size.factor <- data.obj$size.factor
	}
	if (norm == 'GMPR') {
		if (norm.level %in% c('OTU', 'Species', 'ASV')) {
			tab <- data.obj$otu.tab
		} else {
			tab <- data.obj$abund.list[[norm.level]]
		}
		size.factor <- GMPR(tab, intersect.no, ct.min)
		
		# Rev: 2018_01_30
		ind <- !is.na(size.factor)
		data.obj <- subset_data(data.obj, ind)
		grp <- grp[ind]
		df <- df[ind, ]
		size.factor <- size.factor[ind]
	}	

	
	# Method-dependent processing
	if (is.null(method)) {
		if (nlevels(grp) == 2) {
			method <- 'wilcox'
		} else {
			method <- 'kruskal'
		}
	}
	if (method == 'wilcox.pair') {
		if (nlevels(grp) != 2) stop("Wilcox test requires two groups!\n")
		if (is.null(subject)) stop("Paired wilcox needs subject information!\n")
		subject <- factor(df[, subject])
		ind1 <- ind2 <- NULL
		for(sub in levels(subject)) {
			temp1 <- which(as.numeric(grp) == 1 & subject == sub)
			temp2 <- which(as.numeric(grp) == 2 & subject == sub)
			if (length(temp1) != 0 & length(temp2) != 0) {
				ind1 <- c(ind1, temp1[1])
				ind2 <- c(ind2, temp2[1])
			}
			
		}
		#	if (length(ind1) != length(ind2))  warning('Some subjects are not paired!\n')
	}
	
	if (method == 'perm.pair') {
		if (is.null(subject)) stop("Paired permutation test needs subject information!\n")
		subject <- factor(df[, subject])
	}
	
	if (method == 'perm') {
		if (!is.null(subject)) {
			subject <- factor(df[, subject])
		} 
	}
	
	ann <- paste0(ann, '.', method)
	
	if (is.factor(grp)) {
		
		pv.list <- qv.list <- qv.perm.list <- R2.list <- coef.list <- fc.list <- pc.list <- m.list <- nzm.list  <- prv.list <- list()
		res.final <- NULL
		for (LOI in taxa.levels) {
			cat(LOI, "\n")
			ct <- data.obj$abund.list[[LOI]]
			
			# Rev: 2018_10_03 Move here
			if (norm == 'TSS') {
				size.factor <- colSums(data.obj$abund.list[[LOI]])
			}
			
			# Filtering
			prop0 <- t(t(ct) / colSums(ct))
			if (!is.null(prev)) {
				prop0 <- prop0[rowSums(prop0!=0) > prev * ncol(prop0), , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			if (!is.null(minp)) {
				prop0 <- prop0[matrixStats::rowMaxs(prop0) > minp, , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			
			if (!is.null(medianp)) {
				nz.mean <- apply(prop0, 1, function(x) median(x[x!=0]))
				prop0 <- prop0[nz.mean > medianp, , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			
			if (!is.null(cv)) {
				prop0 <- prop0[rowSds(prop0) / rowMeans(prop0) > cv, , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			
			# Normalization
	        prop <- t(t(ct) / size.factor)
			
			# Transformation - Others are possible/may be explored in the future
	        if (transform == 'sqrt') {
				prop <- sqrt(prop)
			}
					 		
			if (method == 'perm') {
				set.seed(seed)
				perm.obj <- permute_differential_analysis(df, prop, grp.name, adj.name, strata=subject, 
						block.perm=block.perm, sqrt.trans=FALSE, perm.no=perm.no, winsor = winsor, winsor.qt = winsor.qt)
				pv.de2 <- perm.obj$p.raw
				names(pv.de2) <- rownames(prop)
			}
			
			# For legacy use
			if (method == 'perm.pair') {
				set.seed(seed)
				perm.obj <- permute_differential_analysis(df, prop, grp.name, adj.name, strata=subject,
						block.perm=FALSE, sqrt.trans=FALSE, perm.no=perm.no, winsor = winsor, winsor.qt = winsor.qt)
				pv.de2 <- perm.obj$p.raw
				names(pv.de2) <- rownames(prop)
			}
			
			pv.vec <- m.vec <- nzm.vec <- prv.vec <-  fc.vec <- pc.vec <-  NULL
			
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				taxon.abund0 <- prop0[taxon, ]
				
				pv <- fc <- pc <- m <- nzm <- prv <- NULL
				if (method == 'wilcox') {
					if (nlevels(grp) != 2) stop("Wilcox test requires two groups!\n")
					pv <- wilcox.test(taxon.abund ~ grp)$p.value
					
				}
				if (method == 'wilcox.pair') {
					pv <- wilcox.test(taxon.abund[ind1], taxon.abund[ind2], paired=T)$p.value
				}
				if (method == 'twopart') {
					if (nlevels(grp) != 2) stop("Two part test requires two groups!\n")
					grp1 <- taxon.abund[as.numeric(grp)==1]
					grp2 <- taxon.abund[as.numeric(grp)==2]
					pv <- twopart.test(grp1, grp2)$p.value
				}	
				if (method == 'kruskal') {
					if (nlevels(grp) <= 2) warning("Kruskal-wallis test requires three or more groups!\n")
					pv <- kruskal.test(taxon.abund ~ grp)$p.value
				}
		
				if (method == 'perm') {
					pv <- pv.de2[taxon]
				}
				
				if (method == 'perm.pair') {
					pv <- pv.de2[taxon]
				}
				m <- tapply(taxon.abund0, grp, function(x) mean(x))
				nzm <- tapply(taxon.abund0, grp, function(x) mean(x[x != 0]))
				prv <- tapply(taxon.abund0, grp, function(x) sum(x != 0))				
				
				# Rev: 2017_02_21 fc change baseline grp
				if (nlevels(grp) == 2) {
					grp.no <- table(grp)
					fc <- log2(m[2] / m[1])
					pc <- prv[2] / grp.no[2] / prv[1] * grp.no[1]
				} else {
					pc <- fc <- NA
				}
				
				pv.vec <- rbind(pv.vec, pv)
				fc.vec <- rbind(fc.vec, fc)
				m.vec <- rbind(m.vec, m)
				nzm.vec <- rbind(nzm.vec, nzm)
				pc.vec <- rbind(pc.vec, pc)
				prv.vec <- rbind(prv.vec, prv / table(grp))	
			}
			
			
			temp <- p.adjust(pv.vec[, 1], 'fdr')
			
			qv.vec <- matrix(temp, ncol=1)
			
			rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(pc.vec) <- rownames(m.vec) <- rownames(nzm.vec) <- rownames(prv.vec) <- rownames(prop)
			colnames(pv.vec) <- 'Pvalue'
			colnames(qv.vec) <- 'Qvalue'
			colnames(fc.vec) <- 'log2(AbundRatio)'
			colnames(pc.vec) <- 'PrevRatio'
			colnames(m.vec) <- paste(levels(grp), 'Mean')
			colnames(nzm.vec) <- paste(levels(grp), 'nzMean')
			colnames(prv.vec) <- paste(levels(grp), 'preval')
			
			pv.list[[LOI]] <- pv.vec
			qv.list[[LOI]] <- qv.vec
			fc.list[[LOI]] <- fc.vec
			pc.list[[LOI]] <- pc.vec
			m.list[[LOI]] <- m.vec
			nzm.list[[LOI]] <- nzm.vec
			prv.list[[LOI]] <- prv.vec
			
			if (method %in% c('perm', 'perm.pair')) {
				qv.perm.list[[LOI]] <- t(t(perm.obj$p.adj.fdr))
				R2.list[[LOI]] <- t(t(perm.obj$R2))
				coef.list[[LOI]] <- perm.obj$coefficients
			}
			
			if (method %in% c('perm', 'perm.pair')) {
				coef.mat <- perm.obj$coefficients
				colnames(coef.mat) <- paste0('coef_', colnames(coef.mat))
				res <- cbind(pv.vec, qv.vec, 'Qvalue.perm'=perm.obj$p.adj.fdr, R2=perm.obj$R2, coef.mat, 
						m.vec, nzm.vec, fc.vec, prv.vec, pc.vec)
				rownames(res) <- rownames(prop)
			} else {
				res <- cbind(pv.vec, qv.vec, m.vec, nzm.vec, fc.vec, prv.vec, pc.vec)
				rownames(res) <- rownames(prop)
			}

			write.csv(res, paste0("Taxa_DifferentialAbundanceAnalysis_", LOI, "_", ann, ".csv"))
			
			if (mt.method == 'fdr') {
				res.final <- rbind(res.final, res[res[, 'Qvalue'] <= cutoff, , drop=F])
			}
			if (mt.method == 'raw') {
				res.final <- rbind(res.final, res[res[, 'Pvalue'] <= cutoff, , drop=F])
			}
			if (mt.method == 'fdr.perm') {
				res.final <- rbind(res.final, res[res[, 'Qvalue.perm'] <= cutoff, , drop=F])
			}
			
		}
		
		if (!is.null(res.final)) {
			colnames(res.final) <- colnames(res)
			write.csv(res.final, paste0("Taxa_DifferentialAbundanceAnalysis_AllLevels_", mt.method, '_', cutoff, "_", ann, ".csv"))
		}
		
		return(list(pv.list=pv.list, fc.list=fc.list, pc.list=pc.list, qv.list=qv.list, qv.perm.list=qv.perm.list, R2.list=R2.list,
						coef.list=coef.list, m.list=m.list))
	} else {
		if (is.null(method)) {
			method <- 'Spearman'
		}
		# Continuous case - currently only has DESeq2 
		pv.list <- qv.list <-  fc.list <-  m.list <- R2.list <- coef.list <- qv.perm.list <- list()
		res.final <- NULL
		for (LOI in taxa.levels) {
			cat(LOI, "\n")
			ct <- data.obj$abund.list[[LOI]]
			
			# Rev: 2018_12_25
			if (norm == 'TSS') {
				size.factor <- colSums(data.obj$abund.list[[LOI]])
			}
			
			# Filtering
			prop0 <- t(t(ct) / colSums(ct))
			if (!is.null(prev)) {
				prop0 <- prop0[rowSums(prop0!=0) > prev * ncol(prop0), , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			if (!is.null(minp)) {
				prop0 <- prop0[matrixStats::rowMaxs(prop0) > minp, , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			
			if (!is.null(medianp)) {
				nz.mean <- apply(prop0, 1, function(x) median(x[x!=0]))
				prop0 <- prop0[nz.mean > medianp, , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			
			if (!is.null(cv)) {
				prop0 <- prop0[rowSds(prop0) / rowMeans(prop0) > cv, , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			
			# Normalization
			prop <- t(t(ct) / size.factor)
			
			# Transformation - Others are possible/may be explored in the future
			if (transform == 'sqrt') {
				prop <- sqrt(prop)
			}
			
			if (method == 'perm') {
				set.seed(seed)
				perm.obj <- permute_differential_analysis(df, prop, grp.name, adj.name, strata=subject, 
						block.perm=block.perm, sqrt.trans=FALSE, perm.no=perm.no, winsor = winsor, winsor.qt = winsor.qt)
				pv.de2 <- perm.obj$p.raw
				names(pv.de2) <- rownames(prop)
				# Place holder
				fc.de2 <- sapply(1:nrow(prop), function(i) cor(prop[i, ], df[, grp.name], method='spearman'))
			}
			
			if (method == 'Spearman') {
				if (!is.null(adj.name)) {
					stop("Spearman test can't adjust covariates!")
				}
				pv.de2 <- apply(prop, 1, function(x) {
							cor.test(x, grp, method='spearman')$p.value
						})
				names(pv.de2) <- rownames(prop)
				# Place holder
				fc.de2 <- sapply(1:nrow(prop), function(i) cor(prop[i, ], df[, grp.name], method='spearman'))
			}
			pv.vec <- matrix(pv.de2, ncol=1)	
			qv.vec <- matrix(p.adjust(pv.vec[, 1], 'fdr'), ncol=1)
			fc.vec <- matrix(fc.de2, ncol=1)
			m.vec <- matrix(rowMeans(prop0), ncol=1)
			
			rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(m.vec)  <- rownames(prop)
			colnames(pv.vec) <- 'Pvalue'
			colnames(qv.vec) <- 'Qvalue'
			colnames(fc.vec) <- 'SpearmanCorr'
			colnames(m.vec) <- 'Mean'
			
			
			pv.list[[LOI]] <- pv.vec
			qv.list[[LOI]] <- qv.vec
			fc.list[[LOI]] <- fc.vec
			m.list[[LOI]] <- m.vec
			
			if (method %in% c('perm', 'perm.pair')) {
				qv.perm.list[[LOI]] <- t(t(perm.obj$p.adj.fdr))
				R2.list[[LOI]] <- t(t(perm.obj$R2))
				coef.list[[LOI]] <- perm.obj$coefficients
			}
			
			if (method %in% c('perm', 'perm.pair')) {
				coef.mat <- perm.obj$coefficients
				colnames(coef.mat) <- paste0('coef_', colnames(coef.mat))
				res <- cbind(pv.vec, qv.vec, 'Qvalue.perm'=perm.obj$p.adj.fdr, R2=perm.obj$R2, coef.mat, fc.vec)
				rownames(res) <- rownames(prop)
			} else {
				res <- cbind(m.vec, fc.vec, pv.vec, qv.vec)
				rownames(res) <- rownames(prop)
			}
			
			
			write.csv(res, paste0("Taxa_DifferentialAbundanceAnalysis_", LOI, "_", ann, ".csv"))
			
			if (mt.method == 'fdr') {
				res.final <- rbind(res.final, res[res[, 'Qvalue'] <= cutoff, , drop=F])
			}
			if (mt.method == 'raw') {
				res.final <- rbind(res.final, res[res[, 'Pvalue'] <= cutoff, , drop=F])
			}
			if (mt.method == 'fdr.perm') {
				res.final <- rbind(res.final, res[res[, 'Qvalue.perm'] <= cutoff, , drop=F])
			}
		}
		
		if (!is.null(res.final)) {
			colnames(res.final) <- colnames(res)
			write.csv(res.final, paste0("Taxa_DifferentialAbundanceAnalysis_AllLevels_", mt.method, '_', cutoff, "_", ann, ".csv"))
		}
		return(list(pv.list=pv.list, fc.list=fc.list, qv.list=qv.list, m.list=m.list, qv.perm.list=qv.perm.list, coef.list=coef.list, R2.list=R2.list))
	}
	
}

bbmix.fit.MM <- function (ct, dep,  nIter = 10, winsor.qt = 1.0) {
  
  if (mean(ct == 0) >= 0.95 | sum(ct != 0) < 5) {
    stop('The number of nonzero is too small to fit the model! Consider removing it from testing!\n')
  }
  
  # Initialization
  prop0 <- ct / dep
  qt <- quantile(prop0, winsor.qt)
  prop0[prop0 >= qt] <- qt
  
  
  var1 <- var(prop0) / 2
  mean1 <- mean(prop0) / 2
  
  var2 <- var(prop0) / 2
  mean2 <- 3 * mean(prop0) / 2
  
  #	var1 <-  var(prop0[prop0 <= median(prop0)])
  #	mean1 <- mean(prop0) / 2
  #	
  #	var2 <- var(prop0[prop0 > median(prop0)])
  #	mean2 <- 3 * mean(prop0) / 2
  
  pi <- 0.5
  
  for (i in 1 : nIter) {
    
    
    shape1.1 <- ((1 - mean1) / var1 - 1 / mean1) * mean1 ^ 2
    shape1.2 <- shape1.1 * (1 / mean1 - 1)
    
    shape2.1 <- ((1 - mean2) / var2 - 1 / mean2) * mean2 ^ 2
    shape2.2 <- shape2.1 * (1 / mean2 - 1)
    
    m1 <- shape1.1 / (shape1.1 + shape1.2)
    s1 <- shape1.1 + shape1.2
    
    m2 <- shape2.1 / (shape2.1 + shape2.2)
    s2 <- shape2.1 + shape2.2
    
    f1 <- dbetabinom(ct, dep, m1, s1)
    f2 <- dbetabinom(ct, dep, m2, s2)
    
    q1 <-  pi * f1 /  (pi * f1 + (1 - pi) * f2)
    q2 <-  1 - q1
    
    pi <- mean(q1)
    
    # Rough estimation
    mean1 <- sum(prop0 * q1) / sum(q1)
    var1 <- sum((prop0 - mean1)^2 * q1) / sum(q1)
    
    mean2 <- sum(prop0 * q2) / sum(q2)
    var2 <- sum((prop0 - mean2)^2 * q2) / sum(q2)
  }
  
  return(list(shape1.1 = shape1.1, shape1.2 = shape1.2, shape2.1 = shape2.1, shape2.2 = shape2.2, pi = pi, q1 = q1))
  
}


perform_differential_analysis_zicoseq <- function (data.obj, grp.name, adj.name=NULL, subject=NULL, 
                                           taxa.levels=c('Phylum', 'Order', 'Class', 'Family', 'Genus'),
                                           prev.filter = 0.1, max.abund.filter = 0.002,
                                           is.winsor = TRUE, outlier.pct = 0.03,
                                           is.post.sample = FALSE, perm.no = 999,
                                           link.func = list(function (x) x^0.5),
                                           mt.method = c('fdr', 'raw'), cutoff = 0.10, 
                                           ann = '', width = 8, height = 7, seed = 123, ...) {
  mt.method <- match.arg(mt.method)
  set.seed(seed)
  
  meta.dat <- data.obj$meta.dat
  grp <- meta.dat[, grp.name]
  if (is.character(grp)) {
    grp <- factor(grp)
  }
  

  ind <- apply(meta.dat[, c(grp.name, adj.name, subject), drop = FALSE], 1, function (x) sum(is.na(x))) == 0


  data.obj <- subset_data(data.obj, ind)
  grp <- grp[ind]
  meta.dat <- meta.dat[ind, ]
  
  if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
    data.obj$abund.list[['Species']] <- data.obj$otu.tab
    rownames(data.obj$abund.list[['Species']]) <- paste0("S", substr(rownames(data.obj$otu.tab), 1, 10), "(", 
                                              substr(data.obj$otu.name[, 'Genus'], 1, 10), ")")
  }
  
  # Test for sequence-depth confounding
  dep <- colSums(data.obj$otu.tab)
  diff.seq.p <- summary(aov(dep ~ grp))[[1]][1, 'Pr(>F)']
  if (!is.na(diff.seq.p) & diff.seq.p <= 0.05 & is.post.sample == FALSE) {
    warning(paste0(
      '\nSignificant sequencing depth confounding with the variable of interest!\n',
      'Posterior sampling will be set to TRUE to reduce false positives!\n'))
      is.post.sample <- TRUE
  }
  
  ann <- paste0(ann, '.', grp.name, '.adj.', paste(adj.name, collapse = '.'))
  
  
  pv.list <- qv.list <- R2.list <- coef.list <- fc.list <- pc.list <- m.list <- nzm.list  <- prv.list <- list()
  res.final <- NULL
  for (LOI in taxa.levels) {
    cat(LOI, "\n")
    feature.dat <- data.obj$abund.list[[LOI]]
    # fit <- MicrobiomeStat::linda(feature.dat = feature.dat, meta.dat = meta.dat, formula = '~icd10_first_3_name + BMI + Sex + Batch')
    # fit1 <- MicrobiomeStat::linda(feature.dat = feature.dat, meta.dat = meta.dat, formula = '~icd10_first_3_name + BMI + Sex')
    # sum(fit$output[[1]]$padj<0.05)
    # sum(fit1$output[[1]]$padj<0.05)
    # adj.name <- c('BMI','Sex')
    # perm.obj <- GUniFrac::ZicoSeq(meta.dat = meta.dat, feature.dat = feature.dat, feature.dat.type = 'count',
    #                               is.winsor = is.winsor, outlier.pct = outlier.pct, 
    #                               prev.filter = prev.filter, max.abund.filter = max.abund.filter,
    #                               grp.name = grp.name, adj.name = adj.name, strata = subject, 
    #                               perm.no = perm.no, is.post.sample = is.post.sample, 
    #                               link.func = link.func,
    #                               return.feature.dat = TRUE)
    
    perm.obj <- ZicoSeq2(meta.dat = meta.dat, feature.dat = feature.dat, feature.dat.type = 'count',
                        is.winsor = is.winsor, outlier.pct = outlier.pct, 
                        prev.filter = prev.filter, max.abund.filter = max.abund.filter,
                        grp.name = grp.name, adj.name = adj.name, strata = subject, 
                        perm.no = perm.no, is.post.sample = is.post.sample, 
                        link.func = link.func,
                        return.feature.dat = TRUE, ...)
    
    cat('ZicoSeq finished!\n')
    if (mt.method == 'fdr') pvalue.type <- 'p.adj.fdr'
    if (mt.method == 'raw') pvalue.type <- 'p.raw'
    plot.obj <- ZicoSeq.plot2(ZicoSeq.obj = perm.obj, grp.name = grp.name, meta.dat = meta.dat, pvalue.type = pvalue.type, cutoff = cutoff)
    cat('Plot is able to be generated!\n')
    pdf(paste0("Taxa_ZicoSeq_", LOI, "_", ann, ".volcano.pdf"), width = width, height = height)
    print(plot.obj)
    dev.off()
    prop <- perm.obj$feature.dat
    prop <- t(t(prop) / colSums(prop))
    
    m.vec <- nzm.vec <- prv.vec <-  fc.vec <- pc.vec <-  NULL
    if (is.factor(grp)) {
      for (taxon in rownames(prop)) {
        taxon.abund0 <- prop[taxon, ]
        
        fc <- pc <- m <- nzm <- prv <- NULL
        m <- tapply(taxon.abund0, grp, function(x) mean(x))
        nzm <- tapply(taxon.abund0, grp, function(x) mean(x[x != 0]))
        prv <- tapply(taxon.abund0, grp, function(x) sum(x != 0))				
        
        # Rev: 2017_02_21 fc change baseline grp
        if (nlevels(grp) == 2) {
          grp.no <- table(grp)
          fc <- log2(m[2] / m[1])
          pc <- prv[2] / grp.no[2] / prv[1] * grp.no[1]
        } else {
          pc <- fc <- NA
        }
        
        fc.vec <- rbind(fc.vec, fc)
        m.vec <- rbind(m.vec, m)
        nzm.vec <- rbind(nzm.vec, nzm)
        pc.vec <- rbind(pc.vec, pc)
        prv.vec <- rbind(prv.vec, prv / table(grp))	
      }
    } else {
      m.vec <- matrix(rowMeans(prop), ncol = 1)
      prv.vec <- matrix(rowMeans(prop != 0), ncol = 1)
      nzm.vec <- matrix(apply(prop, 1, function (x) mean(x[x != 0])), ncol = 1)
      fc.vec <- matrix(sapply(1:nrow(prop), function(i) cor(prop[i, ], grp, method='spearman')), ncol = 1)
      pc.vec <- matrix(sapply(1:nrow(prop), function(i) cor(prop[i, ] != 0, grp, method='spearman')), ncol = 1)               
    }

    pv.vec <- matrix(perm.obj$p.raw, ncol = 1)
    qv.vec <- matrix(perm.obj$p.adj.fdr, ncol = 1)
    
    rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(pc.vec) <-
      rownames(m.vec) <- rownames(nzm.vec) <- rownames(prv.vec) <- rownames(prop)
    colnames(pv.vec) <- 'Pvalue'
    colnames(qv.vec) <- 'Qvalue'
    if (is.factor(grp)) {
      colnames(fc.vec) <- 'log2(AbundRatio)'
      colnames(pc.vec) <- 'PrevRatio'
      colnames(m.vec) <- paste(levels(grp), 'Mean')
      colnames(nzm.vec) <- paste(levels(grp), 'nzMean')
      colnames(prv.vec) <- paste(levels(grp), 'preval')
    } else {
      colnames(fc.vec) <- 'AbundSpearmanCor'
      colnames(pc.vec) <- 'PrevSpearmanCor'
      colnames(m.vec) <- 'Mean'
      colnames(nzm.vec) <- 'nzMean'
      colnames(prv.vec) <- 'Prev'
    }

    pv.list[[LOI]] <- pv.vec
    qv.list[[LOI]] <- qv.vec
    fc.list[[LOI]] <- fc.vec
    pc.list[[LOI]] <- pc.vec
    m.list[[LOI]] <- m.vec
    nzm.list[[LOI]] <- nzm.vec
    prv.list[[LOI]] <- prv.vec
    
    R2.list[[LOI]] <- perm.obj$R2[, 1, drop = FALSE]
    coef.list[[LOI]] <- t(perm.obj$coef.list[[1]][-1, , drop = FALSE])
    
    
    coef.mat <- t(perm.obj$coef.list[[1]][-1, , drop = FALSE])
    colnames(coef.mat) <- paste0('coef_', colnames(coef.mat))
    res <- cbind(pv.vec, qv.vec,  R2=perm.obj$R2[, 1], coef.mat, 
                 m.vec, nzm.vec, fc.vec, prv.vec, pc.vec)
    rownames(res) <- rownames(prop)
    
    write.csv(res, paste0("Taxa_ZicoSeq_", LOI, "_", ann, ".csv"))
    
    if (mt.method == 'fdr') {
      res.final <- rbind(res.final, res[round(res[, 'Qvalue'], 2) <= cutoff, , drop=F])
    }
    if (mt.method == 'raw') {
      res.final <- rbind(res.final, res[round(res[, 'Pvalue'], 2) <= cutoff, , drop=F])
    }
    
  }
  
  if (!is.null(res.final)) {
    colnames(res.final) <- colnames(res)
    write.csv(res.final, paste0("Taxa_ZicoSeq_AllLevels_", mt.method, '_', cutoff, "_", ann, ".csv"))
  }
  
  return(list(pv.list = pv.list, fc.list = fc.list, pc.list = pc.list, qv.list = qv.list,  R2.list = R2.list,
              coef.list = coef.list, m.list = m.list))
} 

## add 09/22/2023 for feature.dat.type support both count and proportion
perform_differential_analysis_zicoseq2 <- function (data.obj, grp.name, adj.name=NULL, subject=NULL, feature.dat.type='count',
                                                   taxa.levels=c('Phylum', 'Order', 'Class', 'Family', 'Genus'),
                                                   prev.filter = 0.1, max.abund.filter = 0.002,
                                                   is.winsor = TRUE, outlier.pct = 0.03,
                                                   is.post.sample = FALSE, perm.no = 999,
                                                   link.func = list(function (x) sign(x) * (abs(x))^0.5),
                                                   mt.method = c('fdr', 'raw'), cutoff = 0.10, 
                                                   ann = '', width = 8, height = 7, seed = 123, ...) {
  mt.method <- match.arg(mt.method)
  set.seed(seed)
  
  meta.dat <- data.obj$meta.dat
  grp <- meta.dat[, grp.name]
  if (is.character(grp)) {grp <- factor(grp)}
  
  
  ind <- apply(meta.dat[, c(grp.name, adj.name, subject), drop = FALSE], 1, function (x) sum(is.na(x))) == 0
  
  
  data.obj <- subset_data(data.obj, ind)
  grp <- grp[ind]
  meta.dat <- meta.dat[ind, ]
  
  if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
    data.obj$abund.list[['Species']] <- data.obj$otu.tab
    rownames(data.obj$abund.list[['Species']]) <- paste0("S", substr(rownames(data.obj$otu.tab), 1, 10), "(", 
                                                         substr(data.obj$otu.name[, 'Genus'], 1, 10), ")")
  }
  
  # Test for sequence-depth confounding
  dep <- colSums(data.obj$otu.tab)
  diff.seq.p <- summary(aov(dep ~ grp))[[1]][1, 'Pr(>F)']
  if (!is.na(diff.seq.p) & diff.seq.p <= 0.05 & is.post.sample == FALSE) {
    warning(paste0(
      '\nSignificant sequencing depth confounding with the variable of interest!\n',
      'Posterior sampling will be set to TRUE to reduce false positives!\n'))
    is.post.sample <- TRUE
  }
  
  ann <- paste0(ann, '.', grp.name, '.adj.', paste(adj.name, collapse = '.'))
  
  
  pv.list <- qv.list <- R2.list <- coef.list <- fc.list <- pc.list <- m.list <- nzm.list  <- prv.list <- list()
  res.final <- NULL
  for (LOI in taxa.levels) {
    cat(LOI, "\n")
    feature.dat <- as.matrix(data.obj$abund.list[[LOI]])
    
    perm.obj <- ZicoSeq2(meta.dat = meta.dat, feature.dat = feature.dat, feature.dat.type = feature.dat.type,
                        is.winsor = is.winsor, outlier.pct = outlier.pct, 
                        prev.filter = prev.filter, max.abund.filter = max.abund.filter,
                        grp.name = grp.name, adj.name = adj.name, strata = subject, 
                        perm.no = perm.no, is.post.sample = is.post.sample, 
                        link.func = link.func,
                        return.feature.dat = TRUE, ...)
    cat('ZicoSeq calculation finished!\n')
    if (mt.method == 'fdr') pvalue.type <- 'p.adj.fdr'
    if (mt.method == 'raw') pvalue.type <- 'p.raw'
    perm.obj$call$feature.dat.type = feature.dat.type
    plot.obj <- ZicoSeq.plot2(ZicoSeq.obj = perm.obj, grp.name = grp.name, meta.dat = meta.dat, pvalue.type = pvalue.type, cutoff = cutoff)
    cat('Generate ZicoSeq plot!\n')
    pdf(paste0("Taxa_ZicoSeq_", LOI, "_", ann, ".volcano.pdf"), width = width, height = height)
    print(plot.obj)
    dev.off()
    cat('ZicoSeq plot finished!\n')
    
    prop <- perm.obj$feature.dat
    prop <- t(t(prop) / colSums(prop))
    
    m.vec <- nzm.vec <- prv.vec <-  fc.vec <- pc.vec <-  NULL
    if (is.factor(grp)) {
      for (taxon in rownames(prop)) {
        taxon.abund0 <- prop[taxon, ]
        
        fc <- pc <- m <- nzm <- prv <- NULL
        m <- tapply(taxon.abund0, grp, function(x) mean(x))
        nzm <- tapply(taxon.abund0, grp, function(x) mean(x[x != 0]))
        prv <- tapply(taxon.abund0, grp, function(x) sum(x != 0))				
        
        # Rev: 2017_02_21 fc change baseline grp
        if (nlevels(grp) == 2) {
          grp.no <- table(grp)
          fc <- log2(m[2] / m[1])
          pc <- prv[2] / grp.no[2] / prv[1] * grp.no[1]
        } else {
          pc <- fc <- NA
        }
        
        fc.vec <- rbind(fc.vec, fc)
        m.vec <- rbind(m.vec, m)
        nzm.vec <- rbind(nzm.vec, nzm)
        pc.vec <- rbind(pc.vec, pc)
        prv.vec <- rbind(prv.vec, prv / table(grp))	
      }
    } else {
      m.vec <- matrix(rowMeans(prop), ncol = 1)
      prv.vec <- matrix(rowMeans(prop != 0), ncol = 1)
      nzm.vec <- matrix(apply(prop, 1, function (x) mean(x[x != 0])), ncol = 1)
      fc.vec <- matrix(sapply(1:nrow(prop), function(i) cor(prop[i, ], grp, method='spearman')), ncol = 1)
      pc.vec <- matrix(sapply(1:nrow(prop), function(i) cor(prop[i, ] != 0, grp, method='spearman')), ncol = 1)               
    }
    
    pv.vec <- matrix(perm.obj$p.raw, ncol = 1)
    qv.vec <- matrix(perm.obj$p.adj.fdr, ncol = 1)
    
    rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(pc.vec) <-
      rownames(m.vec) <- rownames(nzm.vec) <- rownames(prv.vec) <- rownames(prop)
    colnames(pv.vec) <- 'Pvalue'
    colnames(qv.vec) <- 'Qvalue'
    if (is.factor(grp)) {
      colnames(fc.vec) <- 'log2(AbundRatio)'
      colnames(pc.vec) <- 'PrevRatio'
      colnames(m.vec) <- paste(levels(grp), 'Mean')
      colnames(nzm.vec) <- paste(levels(grp), 'nzMean')
      colnames(prv.vec) <- paste(levels(grp), 'preval')
    } else {
      colnames(fc.vec) <- 'AbundSpearmanCor'
      colnames(pc.vec) <- 'PrevSpearmanCor'
      colnames(m.vec) <- 'Mean'
      colnames(nzm.vec) <- 'nzMean'
      colnames(prv.vec) <- 'Prev'
    }
    
    pv.list[[LOI]] <- pv.vec
    qv.list[[LOI]] <- qv.vec
    fc.list[[LOI]] <- fc.vec
    pc.list[[LOI]] <- pc.vec
    m.list[[LOI]] <- m.vec
    nzm.list[[LOI]] <- nzm.vec
    prv.list[[LOI]] <- prv.vec
    
    R2.list[[LOI]] <- perm.obj$R2[, 1, drop = FALSE]
    coef.list[[LOI]] <- t(perm.obj$coef.list[[1]][-1, , drop = FALSE])
    
    
    coef.mat <- t(perm.obj$coef.list[[1]][-1, , drop = FALSE])
    colnames(coef.mat) <- paste0('coef_', colnames(coef.mat))
    res <- cbind(pv.vec, qv.vec,  R2=perm.obj$R2[, 1], coef.mat, 
                 m.vec, nzm.vec, fc.vec, prv.vec, pc.vec)
    rownames(res) <- rownames(prop)
    
    write.csv(res, paste0("Taxa_ZicoSeq_", LOI, "_", ann, ".csv"))
    
    if (mt.method == 'fdr') {
      res.final <- rbind(res.final, res[round(res[, 'Qvalue'], 2) <= cutoff, , drop=F])
    }
    if (mt.method == 'raw') {
      res.final <- rbind(res.final, res[round(res[, 'Pvalue'], 2) <= cutoff, , drop=F])
    }
    
  }
  
  if (!is.null(res.final)) {
    colnames(res.final) <- colnames(res)
    write.csv(res.final, paste0("Taxa_ZicoSeq_AllLevels_", mt.method, '_', cutoff, "_", ann, ".csv"))
  }
  
  return(list(pv.list = pv.list, fc.list = fc.list, pc.list = pc.list, qv.list = qv.list,  R2.list = R2.list,
              coef.list = coef.list, m.list = m.list))
} 



#Owls <- transform(Owls, 
#		Nest=reorder(Nest,NegPerChick), 
#		logBroodSize=log(BroodSize), 
#		NCalls=SiblingNegotiation)
#
#fit_zinbinom <- glmmadmb(NCalls~(FoodTreatment+ArrivalTime)*SexParent+ 
#				offset(logBroodSize)+(1|Nest), 
#		data=Owls, 
#		zeroInflation=TRUE, 
#		family="nbinom")
confint.lme <- function (object, parm, level = 0.95, ...) {
	cf <- fixed.effects(object)
	pnames <- names(cf)
	if (missing(parm)) 
		parm <- pnames
	else if (is.numeric(parm)) 
		parm <- pnames[parm]
	a <- (1 - level)/2
	a <- c(a, 1 - a)
	pct <- stats:::format.perc(a, 3)
	fac <- qnorm(a)
	ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
					pct))
	ses <- sqrt(diag(vcov(object)))[parm]
	ci[] <- cf[parm] + ses %o% fac
	ci
}

# 2020_02_20 adj.name allow multple variables
# Rev: 2016_09_13 add glmmPQL-based overdipersed Poisson and binomial model
# glmmPQL default has overdispersion parameter. Regular Poisson and Binomial does not work.
perform_differential_analysis_para_single_RE <- function (taxon, ldep, grp.name, adj.name = NULL, subject, df, method='NB', LRT=FALSE) {
	# ldep: log depth (size factor)
	if (!is.null(adj.name) ) {
		if (sum(grepl(grp.name, c(adj.name)))) {
			stop('grp.name could not be part of adj.name, or there will be problem!\n')
		}
	}
	
	if (is.null(subject)) {
		stop('Random effects model require the specification of the subject parameter!\n')
	}	
	
	df$ldep <- ldep
	df$taxon <- taxon
	if (LRT) warning('Currently, only Wald test is implemented!\n')
	
	
	if (method %in% c('NB', 'ZINB1', 'ZINB0', 'B0')) {
		if (is.null(adj.name)) {
			grp.name.adj.name.subject <- grp.name.adj.name.subject <- paste(grp.name,  '+ (1|', subject, ')')
		} else {
			adj.name <- paste(adj.name, collapse = ' + ')
			grp.name.adj.name.subject <- paste(grp.name, '+', adj.name, '+ (1|', subject, ')')
		}
		
		if (method == 'NB') {
			m1.nb <- glmmadmb(as.formula(paste('taxon ~', grp.name.adj.name.subject, '+ offset(ldep)')), data = df, 
					zeroInflation=FALSE, family='nbinom')
		}
		
		# 'B0' uses the glmmadmb, which will be compared to 'B' method of glmmPQL
		if (method == 'B0') {
			taxon2 <- as.numeric(taxon != 0)
			df$taxon2 <- taxon2
			m1.nb <- glmmadmb(as.formula(paste('taxon2 ~', grp.name.adj.name.subject, '+ ldep')), data = df, 
					zeroInflation=FALSE, family='binomial')
		}
		
		if (method == 'ZINB1') {
			m1.nb <- glmmadmb(as.formula(paste('taxon ~', grp.name.adj.name.subject, '+ offset(ldep)')), data = df, 
					zeroInflation=TRUE, family='nbinom')
		}
		
		if (method == 'ZINB0') {
			m1.nb <- glmmadmb(as.formula(paste('taxon ~', grp.name.adj.name.subject, '+ offset(ldep)')), data = df, 
					zeroInflation=TRUE, family='truncnbinom')
		}
		
		code <- list(m1.conv=m1.nb$convmsg)	
		pv.nb <- wald.test(b = coef(m1.nb), Sigma = vcov(m1.nb), Terms = grep(grp.name, names(coef(m1.nb))))$result$chi2['P']
		method <- paste(method, 'Wald')
		
		
		aic.nb <- -2 * m1.nb$loglik + 2 * m1.nb$npar	
		coef.nb <- coef(m1.nb)
		ci.nb <- confint.default(m1.nb)
	} 
	
	if (method %in% c('OP', 'B', 'QB')) {
		if (is.null(adj.name)) {
			grp.name.adj.name <- grp.name
		} else {
			adj.name <- paste(adj.name, collapse = ' + ')
			grp.name.adj.name <- paste(grp.name, '+', adj.name)
		}
		
		subject <- paste('~ 1 |', subject)
		if (method == 'OP') {
			m1.nb <- glmmPQL(as.formula(paste('taxon ~', grp.name.adj.name, '+ offset(ldep)')), random = as.formula(subject), data = df, 
					verbose=F, family=quasipoisson())
		}
		
		# 'B0' uses the glmmadmb, which will be compared to 'B' method of glmmPQL
		if (method == 'B') {
			taxon2 <- as.numeric(taxon != 0)
			df$taxon2 <- taxon2
			m1.nb <- glmmPQL(as.formula(paste('taxon2 ~', grp.name.adj.name, '+ ldep')),  random = as.formula(subject), data = df, 
					verbose=F, family=binomial)
		}
		
		if (method == 'QB') {
			taxon2 <- as.numeric(taxon != 0)
			df$taxon2 <- taxon2
			m1.nb <- glmmPQL(as.formula(paste('taxon2 ~', grp.name.adj.name, '+ ldep')),  random = as.formula(subject), data = df, 
					verbose=F, family=quasibinomial())
		}
		
		pv.nb <- wald.test(b = fixed.effects(m1.nb), Sigma = vcov(m1.nb), Terms = grep(grp.name, names(coef(m1.nb))))$result$chi2['P']
		method <- paste(method, 'Wald')
		code <- NA   # placeholder
		aic.nb <- NA # NA placeholder
		coef.nb <- fixed.effects(m1.nb)
		ci.nb <- confint.lme(m1.nb)
		
	}
	
	fc.nb <- coef.nb[grep(grp.name, names(coef.nb))]
	obj <- ci.nb[grep(grp.name, rownames(ci.nb)), ]
	
	if (is.vector(obj)) {
		fc.lc.nb <- obj[1]
		fc.uc.nb <- obj[2]
	} else {
		fc.lc.nb <- obj[, 1]
		fc.uc.nb <- obj[, 2]
		names(fc.lc.nb) <- paste(names(fc.lc.nb), '2.5%')
		names(fc.uc.nb) <- paste(names(fc.uc.nb), '97.5%')
	}
	
	return(list(method=method, pv=pv.nb, lfc=fc.nb, lfc.lci=fc.lc.nb, lfc.uci=fc.uc.nb, aic=aic.nb, code=code))
}


# 2020_02_20 adj.name allow multple variables
perform_differential_analysis_para_single_FE <- function (taxon.abund, ldep, grp.name, adj.name=NULL, subject=NULL, df, method='NB', LRT=FALSE) {
	# ldep: log depth (size factor)
	if (!is.null(adj.name)) {
		if (sum(grepl(grp.name, c(adj.name)))) {
			stop('grp.name could not be part of adj.name or subject, or there will be problem!\n')
		}
	}
	
	if (!is.null(subject)) {
		warnings('Fixed effects model will ignore the subject variable! Please use randome effects model!\n')
	}
	if (LRT & method == 'OP') warning('Overdispersed Poisson does not support LRT! Wald test used!\n')
	
	if (is.null(adj.name)) {
		grp.name.adj.name <- grp.name
	} else {
		adj.name <- paste(adj.name, collapse = ' + ')
		grp.name.adj.name <- paste(grp.name, '+', adj.name)
	}
	if (method == 'NB') {
		m1.nb <- glm.nb(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep)')), data = df)
		if (LRT) {
			m0.nb <- update(m1.nb, as.formula(paste('. ~ . -',  grp.name)))
			code <- list(m1.conv=m1.nb$converged, m1.bound=m1.nb$boundary, m0.conv=m0.nb$converged, m0.bound=m0.nb$boundary)
			pv.nb <- anova(m1.nb, m0.nb)['Pr(Chi)'][2, ]
			method <- paste(method, 'LRT')
		} else {
			code <- list(m1.conv=m1.nb$converged, m1.bound=m1.nb$boundary)	
			pv.nb <- wald.test(b = coef(m1.nb), Sigma = vcov(m1.nb), Terms = grep(grp.name, names(coef(m1.nb))))$result$chi2['P']
			method <- paste(method, 'Wald')
		}

		aic.nb <- summary(m1.nb)$aic
		
		coef.nb <- coef(m1.nb)
		fc.nb <- coef.nb[grep(grp.name, names(coef.nb))]
		
		ci.nb <- confint.default(m1.nb)
		obj <- ci.nb[grep(grp.name, rownames(ci.nb)), ]
		
		if (is.vector(obj)) {
			fc.lc.nb <- obj[1]
			fc.uc.nb <- obj[2]
		} else {
			fc.lc.nb <- obj[, 1]
			fc.uc.nb <- obj[, 2]
			names(fc.lc.nb) <- paste(names(fc.lc.nb), '2.5%')
			names(fc.uc.nb) <- paste(names(fc.uc.nb), '97.5%')
		}
		return(list(method=method, pv=pv.nb, lfc=fc.nb, lfc.lci=fc.lc.nb, lfc.uci=fc.uc.nb, aic=aic.nb, code=code))
	}
	if (method == 'B') {
		taxon.abund2 <- as.numeric(taxon.abund != 0)
		m1.b <- glm(as.formula(paste('taxon.abund2 ~', grp.name.adj.name, '+ ldep')), data = df, family=binomial)
		if (LRT) {
			m0.b <- update(m1.b, as.formula(paste('. ~ . -',  grp.name)))
			code <- list(m1.conv=m1.b$converged, m1.bound=m1.b$boundary, m0.conv=m0.b$converged, m0.bound=m0.b$boundary)
			pv.b <- pchisq(2 * (logLik(m1.b) - logLik(m0.b)), df = df.residual(m0.b) - df.residual(m1.b), lower.tail=FALSE)
			method <- paste(method, 'LRT')
		} else {
			code <- list(m1.conv=m1.b$converged, m1.bound=m1.b$boundary)	
			pv.b <- wald.test(b = coef(m1.b), Sigma = vcov(m1.b), Terms = grep(grp.name, names(coef(m1.b))))$result$chi2['P']
			method <- paste(method, 'Wald')
		}

		aic.b <- summary(m1.b)$aic
		coef.b <- coef(m1.b)
		fc.b <- coef.b[grep(grp.name, names(coef.b))]
		
		ci.b <- confint.default(m1.b)
		obj <- ci.b[grep(grp.name, rownames(ci.b)), ]
		
		if (is.vector(obj)) {
			fc.lc.b <- obj[1]
			fc.uc.b <- obj[2]
		} else {
			fc.lc.b <- obj[, 1]
			fc.uc.b <- obj[, 2]
			names(fc.lc.b) <- paste(names(fc.lc.b), '2.5%')
			names(fc.uc.b) <- paste(names(fc.uc.b), '97.5%')
		}
		return(list(method=method, pv=pv.b, lfc=fc.b, lfc.lci=fc.lc.b, lfc.uci=fc.uc.b, aic=aic.b, code=code))
	}
	
	# Rev: 2016_09_13 add 'QB', No likelihood ratio test
	if (method == 'QB') {
		taxon.abund2 <- as.numeric(taxon.abund != 0)
		m1.b <- glm(as.formula(paste('taxon.abund2 ~', grp.name.adj.name, '+ ldep')), data = df, family=quasibinomial)
		
		code <- list(m1.conv=m1.b$converged, m1.bound=m1.b$boundary)	
		pv.b <- wald.test(b = coef(m1.b), Sigma = vcov(m1.b), Terms = grep(grp.name, names(coef(m1.b))))$result$chi2['P']
		method <- paste(method, 'Wald')
		
		
		aic.b <- summary(m1.b)$aic
		coef.b <- coef(m1.b)
		fc.b <- coef.b[grep(grp.name, names(coef.b))]
		
		ci.b <- confint.default(m1.b)
		obj <- ci.b[grep(grp.name, rownames(ci.b)), ]
		
		if (is.vector(obj)) {
			fc.lc.b <- obj[1]
			fc.uc.b <- obj[2]
		} else {
			fc.lc.b <- obj[, 1]
			fc.uc.b <- obj[, 2]
			names(fc.lc.b) <- paste(names(fc.lc.b), '2.5%')
			names(fc.uc.b) <- paste(names(fc.uc.b), '97.5%')
		}
		return(list(method=method, pv=pv.b, lfc=fc.b, lfc.lci=fc.lc.b, lfc.uci=fc.uc.b, aic=aic.b, code=code))
	}
	
	if (method == 'OP') {
		# No LRT
		m1.op <- glm(as.formula(paste('taxon.abund ~', grp.name.adj.name)), offset=ldep, data = df, family=quasipoisson)
		code <- list(m1.conv=m1.op$converged, m1.bound=m1.op$boundary)

		# pv.op <- pchisq(2 * (logLik(m1.op) - logLik(m0.op)), df = df.residual(m0.op) - df.residual(m1.op), lower.tail=FALSE) # LRT not applicable
		coef.op <- coef(m1.op)		
		pv.op <- wald.test(b = coef.op, Sigma = vcov(m1.op), Terms = grep(grp.name, names(coef.op)))$result$chi2['P']
		method <- paste(method, 'Wald')
		fc.op <- coef.op[grep(grp.name, names(coef.op))]
		
		ci.op <- confint.default(m1.op)
		obj <- ci.op[grep(grp.name, rownames(ci.op)), ]
		
		if (is.vector(obj)) {
			fc.lc.op <- obj[1]
			fc.uc.op <- obj[2]
		} else {
			fc.lc.op <- obj[, 1]
			fc.uc.op <- obj[, 2]
			names(fc.lc.op) <- paste(names(fc.lc.op), '2.5%')
			names(fc.uc.op) <- paste(names(fc.uc.op), '97.5%')
		}
		return(list(method=method, pv=pv.op, lfc=fc.op, lfc.lci=fc.lc.op, lfc.uci=fc.uc.op, aic=NULL, code=code))
	}
	
	if (method == 'ZINB0') {
		m1.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep)')),
				data = df, dist = "negbin", EM = TRUE)
		if (LRT) {
			if (is.null(adj.name)) {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~ offset(ldep)')),
						data = df, dist = "negbin", EM = TRUE)
			} else {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', adj.name, '+ offset(ldep)')),
						data = df, dist = "negbin", EM = TRUE)
			}
			code <- list(m1.conv=m1.zinb$converged, m0.conv=m0.zinb$converged)
			# LRT
			pv.zinb <-  pchisq(2 * (logLik(m1.zinb) - logLik(m0.zinb)), df = df.residual(m0.zinb) - df.residual(m1.zinb), lower.tail=FALSE)
			method <- paste(method, 'LRT')
		} else {
			code <- list(m1.conv=m1.zinb$converged)	
			pv.zinb <- wald.test(b = coef(m1.zinb), Sigma = vcov(m1.zinb), Terms = grep(grp.name, names(coef(m1.zinb))))$result$chi2['P']
			method <- paste(method, 'Wald')
		}
		
		aic.zinb <- -2 * logLik(m1.zinb) + 2 * (m1.zinb$n - m1.zinb$df.residual)
		
		coef.zinb <- coef(m1.zinb)
		fc.zinb <- coef.zinb[grep(grp.name, names(coef.zinb))]
		
		ci.zinb <- confint.default(m1.zinb)
		obj <- ci.zinb[grep(grp.name, rownames(ci.zinb)), ]
		
		if (is.vector(obj)) {
			fc.lc.zinb <- obj[1]
			fc.uc.zinb <- obj[2]
		} else {
			fc.lc.zinb <- obj[, 1]
			fc.uc.zinb <- obj[, 2]
			names(fc.lc.zinb) <- paste(names(fc.lc.zinb), '2.5%')
			names(fc.uc.zinb) <- paste(names(fc.uc.zinb), '97.5%')
		}
		return(list(method=method, pv=pv.zinb, lfc=fc.zinb, lfc.lci=fc.lc.zinb, lfc.uci=fc.uc.zinb, aic=aic.zinb, code=code))
	}
	
	
	if (method == 'ZINB1') {
		m1.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep) | ldep')),
				data = df, dist = "negbin", EM = TRUE)
		if (LRT) {
			if (is.null(adj.name)) {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~ offset(ldep) | ldep')),
						data = df, dist = "negbin", EM = TRUE)
			} else {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', adj.name, '+ offset(ldep) | ldep')),
						data = df, dist = "negbin", EM = TRUE)
			}
			code <- list(m1.conv=m1.zinb$converged, m0.conv=m0.zinb$converged)
			# LRT
			pv.zinb <-  pchisq(2 * (logLik(m1.zinb) - logLik(m0.zinb)), df = df.residual(m0.zinb) - df.residual(m1.zinb), lower.tail=FALSE)
			method <- paste(method, 'LRT')
		} else {
			code <- list(m1.conv=m1.zinb$converged)	
			pv.zinb <- wald.test(b = coef(m1.zinb), Sigma = vcov(m1.zinb), Terms = grep(grp.name, names(coef(m1.zinb))))$result$chi2['P']
			method <- paste(method, 'Wald')
		}
		
		aic.zinb <- -2 * logLik(m1.zinb) + 2 * (m1.zinb$n - m1.zinb$df.residual)
		
		coef.zinb <- coef(m1.zinb)
		fc.zinb <- coef.zinb[grep(grp.name, names(coef.zinb))]
		
		ci.zinb <- confint.default(m1.zinb)
		obj <- ci.zinb[grep(grp.name, rownames(ci.zinb)), ]
		
		if (is.vector(obj)) {
			fc.lc.zinb <- obj[1]
			fc.uc.zinb <- obj[2]
		} else {
			fc.lc.zinb <- obj[, 1]
			fc.uc.zinb <- obj[, 2]
			names(fc.lc.zinb) <- paste(names(fc.lc.zinb), '2.5%')
			names(fc.uc.zinb) <- paste(names(fc.uc.zinb), '97.5%')
		}
		return(list(method=method, pv=pv.zinb, lfc=fc.zinb, lfc.lci=fc.lc.zinb, lfc.uci=fc.uc.zinb, aic=aic.zinb, code=code))
	}
	
	
	if (method == 'ZINB2') {
		m2.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep) |', grp.name.adj.name, '+ ldep')),
				data = df, dist = "negbin", EM = TRUE)
		if (LRT) {
			if (is.null(adj.name)) {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~ offset(ldep) | ldep')),
						data = df, dist = "negbin", EM = TRUE)
			} else {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', adj.name, '+ offset(ldep) |', adj.name, ' + ldep')),
						data = df, dist = "negbin", EM = TRUE)
			}
			
			code <- list(m1.conv=m2.zinb$converged, m0.conv=m0.zinb$converged)
			# LRT
			pv2.zinb <-  pchisq(2 * (logLik(m2.zinb) - logLik(m0.zinb)), df = df.residual(m0.zinb) - df.residual(m2.zinb), lower.tail=FALSE)
			method <- paste(method, 'LRT')
		} else {
			code <- list(m2.conv=m2.zinb$converged)	
			pv2.zinb <- wald.test(b = coef(m2.zinb), Sigma = vcov(m2.zinb), Terms = grep(grp.name, names(coef(m2.zinb))))$result$chi2['P']
			method <- paste(method, 'Wald')
		}

		aic2.zinb <- -2 * logLik(m2.zinb) + 2 * (m2.zinb$n - m2.zinb$df.residual)
		
		coef.zinb <- coef(m2.zinb)
		fc2.zinb <- coef.zinb[grep(grp.name, names(coef.zinb))]
		
		ci.zinb <- confint.default(m2.zinb)
		obj <- ci.zinb[grep(grp.name, rownames(ci.zinb)), ]
		
		if (is.vector(obj)) {
			fc2.lc.zinb <- obj[1]
			fc2.uc.zinb <- obj[2]
		} else {
			fc2.lc.zinb <- obj[, 1]
			fc2.uc.zinb <- obj[, 2]
			names(fc2.lc.zinb) <- paste(names(fc2.lc.zinb), '2.5%')
			names(fc2.uc.zinb) <- paste(names(fc2.uc.zinb), '97.5%')
		}
		return(list(method=method, pv=pv2.zinb, lfc=fc2.zinb, lfc.lci=fc2.lc.zinb, lfc.uci=fc2.uc.zinb, aic=aic2.zinb, code=code))
	}

}

# Rev: 2016_09_14 Add Adaptive 0: switch OP and QB depending on the number of zeros. 
# Rev: 2017_10_30  Support filtering based on coefficient of variation
# Rev: 2018_01_30 Add ct.min and handle NA size factor
perform_differential_analysis_para <- function (data.obj,  grp.name, adj.name=NULL, subject=NULL, RE=FALSE, method='Adaptive0', zerop.cutoff=0.25, ZINB='ZINB1', LRT=FALSE, 
		taxa.levels=c('Phylum', 'Order', 'Class', 'Family', 'Genus', 'Species'), winsor=TRUE, winsor.qt=0.97, norm='GMPR', norm.level='Species', intersect.no=4, ct.min = 2,
		prev=0.1, minp=0.002, medianp=NULL, cv=NULL, mt.method='fdr', cutoff=0.15, ann='', ...) {
	# To be completed
	# subject holds the random effects formula
	if (!RE) {
		if (!(method %in% c('ZINB', 'B', 'QB', 'NB', 'OP', 'Adaptive0', 'Adaptive1', 'Adaptive2'))) stop('The speficied model is not supported!\n')
		perform_differential_analysis_para_single <- perform_differential_analysis_para_single_FE
		if (!is.null(subject)) warning('subject will not be used. Are you sure you want to run fixed effects model? ')
	} else {
		if (!(method %in% c('ZINB', 'B', 'B0', 'QB', 'NB', 'OP', 'Adaptive0', 'Adaptive1', 'Adaptive2'))) stop('The speficied model does not have random effects implementation!\n')
		if (ZINB != 'ZINB1') stop('Currently only ZINB1 is supported!\n')
		if (is.null(subject)) warning('subject is not supplied. Fixed effects model will be used instead!\n')
		perform_differential_analysis_para_single <- perform_differential_analysis_para_single_RE
	}

	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	ind <- !is.na(grp)
	data.obj <- subset_data(data.obj, ind)
	grp <- grp[ind]
	df <- df[ind, ]
	
	if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
		data.obj$abund.list[['Species']] <- data.obj$otu.tab
		rownames(data.obj$abund.list[['Species']]) <- paste0("OTU", rownames(data.obj$otu.tab), ":", 
				data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
	}
	
	dep <- colSums(data.obj$otu.tab)
	diff.seq.p <- summary(aov(dep ~ grp))[[1]][1, 'Pr(>F)']
	
	if (!is.na(diff.seq.p) & diff.seq.p <= 0.05) {
		cat("Signficant sequencing depth confounding!\n")
		cat("For parametric test with sequence depth adjustment, please be cautious about the results!\n")
		cat("There may be potential residual sequence depth confounding!\n")
	}
	
	pv.list <- qv.list <-  fc.list <- fc.lc.list <- fc.uc.list <- met.list <- list()
	res.final <- NULL
	
	if (norm == 'Precalculated') {
		dep <- data.obj$size.factor
	}
	
	if (norm == 'GMPR') {
		if (norm.level %in% c('OTU', 'Species', 'ASV')) {
			tab <- data.obj$otu.tab
		} else {
			tab <- data.obj$abund.list[[norm.level]]
		}
		dep <- GMPR(tab, intersect.no, ct.min)
		
		# Rev: 2018_01_30
		ind <- !is.na(dep)
		data.obj <- subset_data(data.obj, ind)
		grp <- grp[ind]
		df <- df[ind, ]
		dep <- dep[ind]
		
	}
	
	if (norm == 'TSS') {
		dep <- colSums(data.obj$otu.tab)
	}
	
	ldep <- log(dep)
	
	for (LOI in taxa.levels) {
		cat(LOI, "\n")

		taxon.ct <- data.obj$abund.list[[LOI]]
		

		if (winsor == TRUE) {
			# Addressing the outlier (97% percent) or at least one outlier
			
			taxon.ct.p <- t(t(taxon.ct) / dep)
			taxon.ct.p <- apply(taxon.ct.p, 1, function(x) {
						cutoff <- quantile(x, winsor.qt)
						x[x >= cutoff] <- cutoff
						x
					}
			)
			# column/row switch
			taxon.ct <- t(round(taxon.ct.p * dep))
			
		}
		
		prop <- t(t(taxon.ct) / colSums(taxon.ct))
#		if (!is.null(minp)) {
#			prop <- prop[matrixStats::rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
#			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
#		}
#		
#		if (!is.null(medianp)) {
#			nz.mean <- apply(prop, 1, function(x) median(x[x!=0]))
#			prop <- prop[nz.mean > medianp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
#			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
#		}
		
		if (!is.null(prev)) {
			prop <- prop[rowSums(prop!=0) > prev * ncol(prop), , drop=FALSE]	
			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
		}
		if (!is.null(minp)) {
			prop <- prop[matrixStats::rowMaxs(prop) > minp, , drop=FALSE]	
			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
		}
		
		if (!is.null(medianp)) {
			nz.mean <- apply(prop, 1, function(x) median(x[x!=0]))
			prop <- prop[nz.mean > medianp, , drop=FALSE]	
			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
		}
		
		if (!is.null(cv)) {
			prop <- prop[rowSds(prop) / rowMeans(prop) > cv, , drop=FALSE]	
			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
		}
			
		pv.vec <-  fc.vec <- fc.lc.vec <- fc.uc.vec <- met.vec <- conv.vec <- NULL
		obj <- NULL
		for (taxon in rownames(taxon.ct)) {
			cat('.')
			taxon.abund <- taxon.ct[taxon, ]
			
			######## Logistic regression ###############
			if (method == 'B0') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='B0', LRT, ...)) 
			if (method == 'B') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='B', LRT, ...)) 
			if (method == 'QB') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='QB', LRT, ...)) 
			######## Overdispersed Poisson regression #########	
			if (method == 'OP') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='OP', LRT, ...)) 
			######## Negative binomial regression #########
			if (method == 'NB') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='NB', LRT, ...)) 
			######## Zeroinflated negbinomial regression 1 ########
			if (method == 'ZINB') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method=ZINB, LRT, ...)) 
			
			# Adpative 0 selects OP and QB based on the zero proportion (Not optimal)
			if (method == 'Adaptive0') {
				temp <- mean(as.numeric(taxon.abund != 0))
				
				if (temp > zerop.cutoff) {
					error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='QB', LRT, ...)) 
				} else {
					error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='OP', LRT, ...)) 
				}
			}
			
			# Adpative 1 selects NB and ZIB based on AIC
			if (method == 'Adaptive1') {
				error1 <- try(obj1 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='NB', LRT, ...)) 
				error2 <- try(obj2 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method=ZINB, LRT, ...)) 
				if (class(error1) != 'try-error' & class(error2) != 'try-error') {
					if (obj1$aic < obj2$aic) {
						obj <- obj1
					} else {
						obj <- obj2
					}
					error <- error1
				} else {
					# pv == 0 indicates some problems in fitting
					if (class(error1) != 'try-error' & obj1$pv != 0) {
						obj <- obj1
						error <- error1
					} else {
						if (class(error2) != 'try-error' & obj1$pv != 0) {
							obj <- obj2
							error <- error2
						} else {
							error <- error2
						}
					}
				}	
			}
			
			# Adaptive 2 starts with NB model, if it fails, it switches ZINB
			if (method == 'Adaptive2') {
				error1 <- try(obj1 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='NB', LRT, ...)) 
				
				if (class(error1) == 'try-error' | obj1$pv == 0) {
					error2 <- try(obj2 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method=ZINB, LRT, ...)) 
					if (class(error2) != 'try-error') {
						obj <- obj2
						error <- error2
					} else {
						error <- error2
					}
				} else {
					error <- error1
					obj <- obj1
				}
				
			}
			
			# Random Effects model 
			# ZINB, B, NB, Adpative1 is implemented based on glmmADMB
			
            # Set P value NA for those not makes sense
			if (class(error) == "try-error" | abs(obj$lfc) > 100) {
				obj$pv <- obj$lfc <- obj$lfc.lci <- obj$lfc.uci <- obj$method <- NA
			}
			
			
			pv.vec <- rbind(pv.vec, obj$pv)
			fc.vec <- rbind(fc.vec, obj$lfc)
			fc.lc.vec <- rbind(fc.lc.vec, obj$lfc.lci)
			fc.uc.vec <- rbind(fc.uc.vec, obj$lfc.uci)
			met.vec <- rbind(met.vec, obj$method)
			
		}
		cat('\n')

		qv.vec <- matrix(p.adjust(pv.vec[, 1], 'fdr'), ncol=1)
		
		rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(fc.uc.vec) <- rownames(fc.lc.vec) <- rownames(met.vec) <- rownames(prop)
		colnames(pv.vec) <- 'Pvalue'
		colnames(qv.vec) <- 'Qvalue'
        colnames(met.vec) <- 'Method'

		
		pv.list[[LOI]] <- pv.vec
		qv.list[[LOI]] <- qv.vec
		fc.list[[LOI]] <- fc.vec
		fc.lc.list[[LOI]] <- fc.lc.vec
		fc.uc.list[[LOI]] <- fc.uc.vec
		met.list[[LOI]] <- met.vec
 
		res <- cbind(pv.vec, qv.vec, fc.vec, fc.lc.vec, fc.uc.vec, met.vec)
		rownames(res) <- rownames(prop)
		write.csv(res, paste0("Taxa_DifferentialAbundanceAnalysis_", LOI, "_", ann, ".csv"))
		
		if (mt.method == 'fdr') {
			res.final <- rbind(res.final, res[as.numeric(res[, 'Qvalue']) <= cutoff, , drop=F])
		}
		if (mt.method == 'raw') {
			res.final <- rbind(res.final, res[ as.numeric(res[, 'Pvalue']) <= cutoff, , drop=F])
		}
	}
	
	if (!is.null(res.final)) {
		colnames(res.final) <- colnames(res)
		res.final <- res.final[rowSums(is.na(res.final)) == 0, , drop=F]
		write.csv(res.final, paste0("Taxa_DifferentialAbundanceAnalysis_AllLevels_", mt.method, '_', cutoff, "_", ann, ".csv"))
	}
	return(list(pv.list=pv.list, qv.list=qv.list, fc.list=fc.list, fc.uc.list=fc.uc.list, fc.lc.list=fc.lc.list, met.list=met.list))
}


plot_effect_size <- function (month,  value, pos.lab, neg.lab, ylab, hjust1=1.3, hjust2=-0.3, lab.size=3, xsize=10) {
	month <- factor(month, levels=month)
	dtm <- data.frame(month=month, value=value)
	dtm$colour <- factor(ifelse(dtm$value < 0, neg.lab, pos.lab), levels=c(pos.lab, neg.lab))
	dtm$hjust <- ifelse(dtm$value > 0, hjust1, hjust2)
	obj <- ggplot(dtm, aes(month, value, label = month, hjust = hjust)) + 
			geom_text(aes(y = 0, colour = colour), size=lab.size) + 
			geom_bar(stat = "identity", aes(fill = colour)) +
			theme(axis.text.x = element_text(size=xsize)) +
            ylim(c(-max(abs(value))*1.1, max(abs(value))*1.1)) + 
			coord_flip() +
			scale_x_discrete(breaks = NULL) +
			labs(x = "", y = ylab) +
			theme(legend.position="top", legend.title=element_blank())
}


plot_effect_size2 <- function (fold.dat.plot1, ylabel='log(Fold change)', is.ln = TRUE, ord = TRUE, yint = 0) {
	
	if (is.ln) {
		fold.dat.plot1[, c('Estimate', 'LCI', 'UCI')] <- fold.dat.plot1[, c('Estimate', 'LCI', 'UCI')] * log2(exp(1))
	}
	Alphas <- seq(1, 99, 2) / 100
#Multiplier <- qnorm(1 - Alphas / 2)
	Multiplier <- seq(1, 0.01, len=50)
	zzTransparency <<- 1/(length(Multiplier)/4)
	
#	fold.dat.plot1 <- data.frame(IV=rownames(fold.dat.plot), Estimate=fold.dat.plot[, 1], LCI=fold.dat.plot[, 2], UCI=fold.dat.plot[, 3])
	fold.dat.plot1 <- data.frame(cbind(fold.dat.plot1, Scalar=rep(Multiplier, each = nrow(fold.dat.plot1))))
#	fold.dat.plot1[,  -1] <- apply(fold.dat.plot1[, -1], 2, function(x){as.numeric(as.character(x))})
	fold.dat.plot1$Emphasis <- by(1 - seq(0, 1, length = length(Multiplier) + 1)[-(length(Multiplier) + 1)],
			as.character(round(Multiplier, 5)), mean)[as.character(round(fold.dat.plot1$Scalar, 5))]
	
	fold.dat.plot1$IV <- factor(fold.dat.plot1$IV, unique(fold.dat.plot1$IV))
	
	if (ord) {
		fold.dat.plot1 <- fold.dat.plot1[order(as.character(fold.dat.plot1$IV)), ]
	}
	
	OutputPlot <- ggplot2::qplot(data = fold.dat.plot1, x = IV, y = Estimate,
			ymin = Estimate - (Estimate -LCI)*Scalar, ymax = Estimate + (UCI - Estimate)*Scalar,
			ylab = NULL, xlab = NULL, alpha = I(zzTransparency), colour = I(gray(0)), geom = "blank")
	
	OutputPlot <- OutputPlot + geom_hline(yintercept = yint, lwd = I(7/12), colour = I(hsv(0/12, 7/12, 7/12)), alpha = I(5/12))
	OutputPlot <- OutputPlot + geom_linerange(data = fold.dat.plot1, aes(size = 1/Emphasis), alpha = I(zzTransparency), colour = I(gray(0)))
	OutputPlot <- OutputPlot + scale_size_continuous() + guides(size=FALSE)
#OutputPlot <- OutputPlot + facet_grid(~ ModelName)
	OutputPlot <- OutputPlot + coord_flip() + geom_point(aes(x = IV, y = Estimate), colour = I(gray(0))) + theme_bw() + ylab(ylabel)
}

# Rev: 2016_11_25 Uniqufy
# Rev: 2016_04_18 strata


# New: 2020_02_20  for survival
plot_effect_size3 <- function (fold.dat.plot1, ylabel = 'HR',  ord.type = c('effect', 'taxa'), yint = 1) {
	
	ord.type <- match.arg(ord.type)
	
	Alphas <- seq(1, 99, 2) / 100
#Multiplier <- qnorm(1 - Alphas / 2)
	Multiplier <- seq(1, 0.01, len=50)
	zzTransparency <<- 1/(length(Multiplier)/4)
	
#	fold.dat.plot1 <- data.frame(IV=rownames(fold.dat.plot), Estimate=fold.dat.plot[, 1], LCI=fold.dat.plot[, 2], UCI=fold.dat.plot[, 3])
	fold.dat.plot1 <- data.frame(cbind(fold.dat.plot1, Scalar=rep(Multiplier, each = nrow(fold.dat.plot1))))
#	fold.dat.plot1[,  -1] <- apply(fold.dat.plot1[, -1], 2, function(x){as.numeric(as.character(x))})
	fold.dat.plot1$Emphasis <- by(1 - seq(0, 1, length = length(Multiplier) + 1)[-(length(Multiplier) + 1)],
			as.character(round(Multiplier, 5)), mean)[as.character(round(fold.dat.plot1$Scalar, 5))]
	
	fold.dat.plot1$IV <- factor(fold.dat.plot1$IV, unique(fold.dat.plot1$IV))
	
	if (ord.type == 'taxa') {
		fold.dat.plot1 <- fold.dat.plot1[order(as.character(fold.dat.plot1$IV)), ]
	}
	if (ord.type == 'effect') {
		fold.dat.plot1 <- fold.dat.plot1[order(as.character(fold.dat.plot1$Estimate)), ]
	}
	
	OutputPlot <- ggplot2::qplot(data = fold.dat.plot1, x = IV, y = Estimate,
			ymin = Estimate - (Estimate -LCI) * Scalar, ymax = Estimate + (UCI - Estimate)*Scalar,
			ylab = NULL, xlab = NULL, alpha = I(zzTransparency), colour = I(gray(0)), geom = "blank")
	
	OutputPlot <- OutputPlot + geom_hline(yintercept = yint, lwd = I(7/12), colour = I(hsv(0/12, 7/12, 7/12)), alpha = I(5/12))
	OutputPlot <- OutputPlot + geom_linerange(data = fold.dat.plot1, aes(size = 1/Emphasis), alpha = I(zzTransparency), colour = I(gray(0)))
	OutputPlot <- OutputPlot + scale_size_continuous() + guides(size=FALSE)
#OutputPlot <- OutputPlot + facet_grid(~ ModelName)
	OutputPlot <- OutputPlot + coord_flip() + geom_point(aes(x = IV, y = Estimate), colour = I(gray(0))) + theme_bw() + ylab(ylabel)
}


# Rev: 2020_03_08 Add support for numeric grp.name.n 
visualize_differential_analysis <- function (data.obj, diff.obj,  grp.name=NULL, grp.name.n = NULL,  strata=NULL, test='Nonpara', mt.method='fdr', scale='sqrt', cutoff=0.15,
                                             taxa.levels=c('Phylum', 'Family', 'Genus'), ord=TRUE, eff.type='logP', indivplot=TRUE, colFnsC=NULL, colFnsF=NULL, subject=NULL,
                                             xsize=10, jsize.agg = 0,  ann='', hei1=NULL, wid1=NULL, hei2=NULL, wid2=NULL) {
  
  # uniquefy names
  # For backward compatibility. Newer version will not need this and below. The old version has 'unclassified' which leads to duplicate names.
  # Newer version has 'Unclassified'. Case difference.
  
  # Check whether there is name duplication
  check.names <- NULL
  
  obj0 <- diff.obj[[1]]
  for (level in names(obj0)) {
    obj <- obj0[[level]]
    # rownames(obj) <- gsub('unclassified', paste0('Unclassified',substr(level, 1, 1)), rownames(obj))
    check.names <- c(check.names, rownames(obj))
  }
  
  if (sum(table(check.names) >= 2)) {
    data.obj <- uniquefy_taxa_names(data.obj)
    
    for (name1 in names(diff.obj)) {
      obj0 <- diff.obj[[name1]]
      for (level in names(obj0)) {
        obj <- obj0[[level]]
        # rownames(obj) <- gsub('unclassified', paste0('Unclassified',substr(level, 1, 1)), rownames(obj))
        rownames(obj) <- paste0(rownames(obj), substr(level, 1, 1))
        obj0[[level]] <- obj
      }
      diff.obj[[name1]] <- obj0
    }
    
  }
  
  #pv.list <- diff.obj$pv.list
  fc.list <- diff.obj$fc.list
  
  if (mt.method == 'fdr.perm') {
    qv.list <- diff.obj$qv.perm.list
  } else {
    qv.list <- diff.obj$qv.list
  }
  
  pv.list <- diff.obj$pv.list
  if (test == 'Para') {
    fc.lc.list <- diff.obj$fc.lc.list
    fc.uc.list <- diff.obj$fc.uc.list
  }
  
  R2.list <- diff.obj$R2.list
  
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  
  
  if (!is.null(grp.name.n)) {
    grp.n <- df[, grp.name.n]
  }
  
  ind <- !is.na(grp)
  data.obj <- subset_data(data.obj, ind)
  grp <- grp[ind]
  df <- df[ind, ]
  
  prop <- NULL
  eff <- eff.lc <- eff.uc <- R2 <- NULL
  taxa.names <- NULL
  if (is.null(taxa.levels)) {
    LOIs <- names(qv.list)
  } else {
    LOIs <- taxa.levels
    if (sum(!(taxa.levels %in% names(qv.list)))) {
      stop('Taxa levels are not contained in differential abundance analysis results!\n')
    }
  }
  for (LOI in LOIs) {
    pv.vec <- pv.list[[LOI]]
    fc.vec <- fc.list[[LOI]]
    #qv.vec <- qvalue(pv.vec[, 1])$qvalues
    qv.vec <- qv.list[[LOI]]
    
    if (!is.null(R2.list)) R2.vec <- R2.list[[LOI]]
    
    if (test == 'Para') {
      fc.lc.vec <- fc.lc.list[[LOI]]
      fc.uc.vec <- fc.uc.list[[LOI]]
    }
    
    if (mt.method %in% c('fdr', 'fdr.perm')) {
      taxa.name <- rownames(qv.vec)[qv.vec <= cutoff]
      taxa.name <- taxa.name[!is.na(taxa.name)]
    }
    
    if (mt.method == 'raw') {
      taxa.name <- rownames(pv.vec)[pv.vec <= cutoff]
      taxa.name <- taxa.name[!is.na(taxa.name)]
    }
    
    if (length(taxa.name) != 0) {
      prop0 <- data.obj$abund.list[[LOI]]
      prop0 <- t(t(prop0) / colSums(prop0))
      prop0 <-  prop0[taxa.name, , drop=F]
      if (ord == TRUE) {
        prop0 <- prop0[rev(order(rowMeans(prop0))), , drop=F]
      }
      prop <- rbind(prop, prop0)
      # currently using fold change
      if (test == 'Para') {
        eff <- rbind(eff, fc.vec[taxa.name, , drop=F])
        eff.lc <- rbind(eff.lc, fc.lc.vec[taxa.name, , drop=F])
        eff.uc <- rbind(eff.uc, fc.uc.vec[taxa.name, , drop=F])
      } else {
        if (eff.type %in% c('LFC', 'Spearman')) {
          eff <- c(eff, fc.vec[taxa.name, ])
        }
        if (eff.type == 'logP') {
          eff <- c(eff, sign(fc.vec[taxa.name, ]) * (-log10(pv.vec[taxa.name, ])))
        }
        if (!is.null(R2.list)) R2 <- c(R2, R2.vec[taxa.name, ])
        
      }
      taxa.names <- c(taxa.names, taxa.name)
    }
  }
  
  if (length(taxa.names) == 0) {
    cat('No differential taxa! \n')
  } else {
    if (length(taxa.names) >= 2) {
      if (is.null(wid1) | is.null(hei1)) {
        wid1 <- 7 * ifelse(nrow(prop) / 30 < 1, 1, nrow(prop) / 30)
        hei1 <- 7
      }
      
      if (!is.null(grp.name)) {
        obj <- taxa_barplot_aggregate(prop, df, grp.name, strata, scale, xsize)
        pdf(paste("Taxa_DifferentialAbundance_AbundanceBarplot", scale, mt.method, cutoff, ann, ".pdf", sep="_"), height=hei1, width=wid1)
        print(obj)
        dev.off()
        
        obj1 <- taxa_boxplot_aggregate(prop, df, grp.name, strata, scale, xsize, jsize = jsize.agg) 
        pdf(paste("Taxa_DifferentialAbundance_AbundanceBoxplot", scale, mt.method, cutoff, ann, ".pdf", sep="_"), height=hei1, width=wid1)
        print(obj1)
        dev.off()
        
      }
      
      # currently fold change
      if (test == 'Para') {
        rownames(eff) <- rownames(eff.lc) <- rownames(eff.uc) <- taxa.names
        for (k in 1:ncol(eff)) {
          fold.dat.plot1 <- data.frame(Estimate=eff[, k], LCI=eff.lc[, k], UCI=eff.uc[, k], IV=taxa.names)
          obj <- plot_effect_size2(fold.dat.plot1)
          pdf(paste("Taxa_DifferentialAbundance_EffectBarplot", colnames(eff)[k], mt.method, cutoff, ann, ".pdf", sep="_"), height=hei2, width=wid2)
          print(obj)
          dev.off()
        }
      } else {
        if (!is.null(eff)) {
          names(eff) <- taxa.names
          eff <- eff[!is.na(eff) & is.finite(eff)]
          eff <- sort(eff)
          taxa.names2 <- names(eff)
          if (is.null(wid2) | is.null(hei2)) {
            hei2 <- 4 + length(taxa.names2) / 20 * 3 
            wid2 <- 6
          }
          if (eff.type == 'LFC') {
            obj <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='Log2 fold change')
            pdf(paste("Taxa_DifferentialAbundance_LFCBarplot", mt.method, cutoff, ann, ".pdf", sep="_"), height=hei2, width=wid2)
            print(obj)
            dev.off()	
          }
          
          if (eff.type == 'Spearman') {
            obj <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='Spearman correlation')
            pdf(paste("Taxa_DifferentialAbundance_SpearmanBarplot", mt.method, cutoff, ann, ".pdf", sep="_"), height=hei2, width=wid2)
            print(obj)
            dev.off()	
          }
          
          if (eff.type == 'logP') {
            obj <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='-log10(P)')
            pdf(paste("Taxa_DifferentialAbundance_logPBarplot", mt.method, cutoff, ann, ".pdf", sep="_"), height=hei2, width=wid2)
            print(obj)
            dev.off()	
          }	
        } else {
          if (!is.null(R2.list))   {
            R2 <- sort(R2, decr=TRUE)
            taxa <- names(R2)
            taxa <- factor(taxa, levels=taxa)
            df <- data.frame(taxa = taxa, R2 = R2)
            obj <- ggplot(df, aes(taxa, R2)) + geom_col(fill='steelblue', col='gray', size=0.25) + ylab(expression(R^2)) + xlab('') +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
            pdf(paste("Taxa_DifferentialAbundance_R2_Barplot", mt.method, cutoff, ann, ".pdf", sep="_"), height=wid2, width=hei2)
            print(obj)
            dev.off()
            
          }
        }
        
        
      }
    }
    
    # create heatmp
    #			taxa.names2 <- taxa.names[!grepl('unclassified', taxa.names, ignore.case=T)]
    if (!is.null(grp.name)) {
      
      generate_taxa_heatmap(data.obj, taxa.levels='All', sam.ord=order(grp), taxa=taxa.names, meta.info=c(grp.name, strata), Colv=F, dendrogram='row', 
                            ann=paste0(mt.method, '_', cutoff, '_', ann, '_Unclustered'), colFnsC=colFnsC, colFnsF=colFnsF)
      generate_taxa_heatmap(data.obj, taxa.levels='All', taxa=taxa.names, meta.info=c(grp.name, strata),  ann=paste0(mt.method, '_', cutoff, '_', ann), colFnsC=colFnsC, colFnsF=colFnsF)
      generate_taxa_heatmap(data.obj, taxa.levels='All', taxa=taxa.names, meta.info=c(grp.name, strata), data.type='R', ann=paste0(mt.method, '_', cutoff, '_Rank_', ann), colFnsC=colFnsC, colFnsF=colFnsF)
      try(
        generate_taxa_biplot(data.obj, taxa=taxa.names, trans='sqrt', grp.name, ann=paste0(mt.method, '_', cutoff, '_', ann), varname.size = 1.5)	
      )	
    }
  }
  if (!is.null(grp.name)) {
    # Individual plots
    if (indivplot == TRUE) {
      generate_taxa_boxplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, ann=paste0(mt.method, '_', cutoff, '_', ann))
      if (!is.null(subject)) {
        generate_taxa_boxplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, subject=subject, ann=paste0(mt.method, '_', cutoff, '_', ann, '_Paired'))
      }
      generate_taxa_boxplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, scale='binary', ann=paste0(mt.method, '_', cutoff, '_', ann))
      generate_taxa_barplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, ann=paste0(mt.method, '_', cutoff, '_', ann))
      
      if (!is.null(grp.name.n)) {
        
        if (!is.null(strata)) {
          generate_taxa_scatterplot(data.obj, grp.name=grp.name.n, taxa.levels='All', strata=strata, combine.strata=FALSE, taxa.name=taxa.names, scale='sqrtP', 
                                    subject=subject, subject.pt=TRUE, ann=paste0(mt.method, '_', cutoff, '_s1_', ann))
          generate_taxa_scatterplot(data.obj, grp.name=grp.name.n, taxa.levels='All', strata=strata, combine.strata=TRUE,  taxa.name=taxa.names, scale='sqrtP',
                                    subject=subject,  subject.pt=TRUE, ann=paste0(mt.method, '_', cutoff, '_s2_', ann))
        } else {
          generate_taxa_scatterplot(data.obj, grp.name=grp.name.n, taxa.levels='All', taxa.name=taxa.names, scale='sqrtP', 
                                    subject=subject, subject.pt=TRUE, ann=paste0(mt.method, '_', cutoff, '_', ann))
        }
        
        
        
      }
    }
    
  }
}


# TODO: Add comment
# 
# Author: jchen1981
###############################################################################
###############################################################################
# New: 2020_02_20
perform_survival_analysis_para <- function (data.obj, obstime.name, delta.name, adj.name = NULL, subject = NULL, 
		RE = FALSE, method = 'cph', 
		taxa.levels = c('Phylum', 'Order', 'Class', 'Family', 'Genus', 'Species'), transform = c('sqrt', 'binary'),
		winsor = TRUE, winsor.qt = 0.97, 
		norm = 'GMPR', norm.level = 'Species', intersect.no=4, ct.min = 2,
		prev = 0.1, minp = 0.002, 
		mt.method = 'fdr', cutoff = 0.1049, ann = '', ...) {
	
	transform <- match.arg(transform)
	# To be completed
	# subject holds the random effects formula
	if (!RE) {
		if (!(method %in% c('cph'))) stop('The speficied method is not supported!\n')
		if (!is.null(subject)) warning('subject will not be used. Are you sure you want to run fixed effects model? ')
	} else {
		if (!(method %in% c('cph'))) stop('The speficied method does not have random effects implementation!\n')
		if (is.null(subject)) stop('subject is not supplied!\n')
	}
	
	# Regression model will handle NA automatically
	df <- data.obj$meta.dat
	obstime <- df[, obstime.name]
	delta <- df[, delta.name]
	
	
	if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
		data.obj$abund.list[['Species']] <- data.obj$otu.tab
		rownames(data.obj$abund.list[['Species']]) <- paste0("OTU", rownames(data.obj$otu.tab), ":", 
				data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
	}
	
	dep <- colSums(data.obj$otu.tab)
	diff.seq.p <- coef(summary(coxph(Surv(obstime, delta) ~  dep)))[1, 'Pr(>|z|)']
	
	if (!is.na(diff.seq.p) & diff.seq.p <= 0.05) {
		warning("Sequencing depth is signficantly correlated with surival times!\n")
	}
	
	if (norm == 'Precalculated') {
		dep <- data.obj$size.factor
	}
	
	if (norm == 'GMPR') {
		if (norm.level %in% c('OTU', 'Species', 'ASV')) {
			tab <- data.obj$otu.tab
		} else {
			tab <- data.obj$abund.list[[norm.level]]
		}
		
		dep <- GMPR(tab, intersect.no, ct.min)
		
		ind <- !is.na(dep)
		data.obj <- subset_data(data.obj, ind)
		obstime <- obstime[ind]
		delta <- delta[ind]
		df <- df[ind, ]
		dep <- dep[ind]	
	}
	
	if (norm == 'TSS') {
		dep <- colSums(data.obj$otu.tab)
	}
	
	ldep <- log(dep)
	
	pv.list <- qv.list <-  hr.list <- hr.lc.list <- hr.uc.list <- met.list <- list()
	res.final <- NULL
	
	for (LOI in taxa.levels) {
		cat(LOI, "\n")
		taxon.ct <- data.obj$abund.list[[LOI]]
		
		if (winsor == TRUE) {
			# Addressing the outlier (97% percent) or at least one outlier
			taxon.ct.p <- t(t(taxon.ct) / dep)
			taxon.ct.p <- apply(taxon.ct.p, 1, function(x) {
						cutoff <- quantile(x, winsor.qt)
						x[x >= cutoff] <- cutoff
						x
					}
			)
			# column/row switch
			taxon.ct <- t(round(taxon.ct.p * dep))
		}
		
		# original relative abundance
		prop <- t(t(taxon.ct) / colSums(taxon.ct))
		
		if (!is.null(prev)) {
			prop <- prop[rowSums(prop != 0) > prev * ncol(prop), , drop = FALSE]	
			taxon.ct <- taxon.ct[rownames(prop), , drop = FALSE]
		}
		
		if (!is.null(minp)) {
			prop <- prop[matrixStats::rowMaxs(prop) > minp, , drop = FALSE]	
			taxon.ct <- taxon.ct[rownames(prop), , drop = FALSE]
		}
		
#		if (!is.null(medianp)) {
#			nz.mean <- apply(prop, 1, function(x) median(x[x!=0]))
#			prop <- prop[nz.mean > medianp, , drop=FALSE]	
#			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
#		}
#		
#		if (!is.null(cv)) {
#			prop <- prop[rowSds(prop) / rowMeans(prop) > cv, , drop=FALSE]	
#			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
#		}
#		
		pv.vec <-  hr.vec <- hr.lc.vec <- hr.uc.vec <- met.vec <- conv.vec <- NULL
		obj <- NULL
		
		for (taxon in rownames(taxon.ct)) {
			cat('.')
			
			taxon.abund <- taxon.ct[taxon, ]
			# cox proportional hazard model
			if (!RE) {
				if (method == 'cph') {
					error <- try({
								obj <- survival_cph_fe(taxon.abund, dep, transform, obstime.name, delta.name, adj.name,  df, ...)
							}
					)
				}		
			}
			# Set P value NA for those not makes sense
			if (class(error) == "try-error") {
				obj$pv <- obj$hr <- obj$hr.lci <- obj$hr.uci <- obj$method <- NA
			}
			
			pv.vec <- rbind(pv.vec, obj$pv)
			hr.vec <- rbind(hr.vec, obj$hr)
			hr.lc.vec <- rbind(hr.lc.vec, obj$hr.lci)
			hr.uc.vec <- rbind(hr.uc.vec, obj$hr.uci)
			met.vec <- rbind(met.vec, obj$method)
			
		}
		cat('\n')
		qv.vec <- matrix(p.adjust(pv.vec[, 1], 'fdr'), ncol=1)
		
		rownames(pv.vec) <- rownames(qv.vec) <- rownames(hr.vec) <- rownames(hr.uc.vec) <- rownames(hr.lc.vec) <- rownames(met.vec) <- rownames(prop)
		colnames(pv.vec) <- 'Pvalue'
		colnames(qv.vec) <- 'Qvalue'
		colnames(met.vec) <- 'Method'
		
		pv.list[[LOI]] <- pv.vec
		qv.list[[LOI]] <- qv.vec
		hr.list[[LOI]] <- hr.vec
		hr.lc.list[[LOI]] <- hr.lc.vec
		hr.uc.list[[LOI]] <- hr.uc.vec
		met.list[[LOI]] <- met.vec
		
		res <- cbind(pv.vec, qv.vec, hr.vec, hr.lc.vec, hr.uc.vec, met.vec)
		
		rownames(res) <- rownames(prop)
		
		write.csv(res, paste0("Taxa_SurvivalAnalysis_", LOI, "_", ann, ".csv"))
		
		if (mt.method == 'fdr') {
			res.final <- rbind(res.final, res[as.numeric(res[, 'Qvalue']) <= cutoff, , drop=F])
		}
		
		if (mt.method == 'raw') {
			res.final <- rbind(res.final, res[as.numeric(res[, 'Pvalue']) <= cutoff, , drop=F])
		}
	}
	
	if (!is.null(res.final)) {
		colnames(res.final) <- colnames(res)
		res.final <- res.final[rowSums(is.na(res.final)) == 0, , drop=F]
		write.csv(res.final, paste0("Taxa_SurvivalAnalysis_AllLevels_", mt.method, '_', cutoff, "_", ann, ".csv"))
	}
	return(list(pv.list = pv.list, qv.list = qv.list, hr.list = hr.list, 
					hr.uc.list = hr.uc.list, hr.lc.list = hr.lc.list, met.list = met.list))
}

# New: 2020_02_20
# CPH for fixed effects
survival_cph_fe <- function(taxon.abund, depth, transform = c('sqrt', 'binary'), obstime.name, delta.name, adj.name,  df, ...) {
	
	if (transform == 'sqrt') abundance <- scale(sqrt(taxon.abund / depth))
	if (transform == 'binary') abundance <- as.numeric(taxon.abund == 0)
	
	if (length(grep('abundance', adj.name)) != 0) {
		stop('The adj.name should not contain "abundance"!\n')
	}
	
	obstime <- df[, obstime.name]
	delta <- df[, delta.name]
	if (is.null(adj.name)) {
		obj <- coxph(Surv(obstime, delta) ~ abundance) 
	} else {
		X <- scale(model.matrix(~. , df[, adj.name])[, -1])
		obj <- coxph(Surv(obstime, delta) ~ X + abundance) 
	}
	
	tab <- coef(summary(obj))
	tab <- tab[grep('abundance', rownames(tab)), ]
	
	res <- NULL
	res$pv <- tab["Pr(>|z|)"]
	res$hr <- tab["exp(coef)"]
	res$hr.lci <- exp(tab["coef"] - 1.96 * tab["se(coef)"])
	res$hr.uci <- exp(tab["coef"] + 1.96 * tab["se(coef)"])
	names(res$hr.lci) <- paste('2.5%')
	names(res$hr.uci) <- paste('97.5%')
	res$method <- 'cph(fe)'
	
	return(res)
}

visualize_survival_analysis <- function (data.obj, surv.obj,  mt.method = 'fdr', cutoff = 0.1049,
		taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'), 
		indivplot = FALSE, wid2 = NULL, hei2 = NULL, ann = '') {
	
	# uniquefy names
	# For backward compatibility. Newer version will not need this and below. The old version has 'unclassified' which leads to duplicate names.
	# Newer version has 'Unclassified'. Case difference.
	
	# Check whether there is name duplication
	check.names <- NULL
	
	obj0 <- surv.obj[[1]]
	for (level in names(obj0)) {
		obj <- obj0[[level]]
		# rownames(obj) <- gsub('unclassified', paste0('Unclassified',substr(level, 1, 1)), rownames(obj))
		check.names <- c(check.names, rownames(obj))
	}
	
	if (sum(table(check.names) >= 2)) {
		data.obj <- uniquefy_taxa_names(data.obj)
		
		for (name1 in names(surv.obj)) {
			obj0 <- surv.obj[[name1]]
			for (level in names(obj0)) {
				obj <- obj0[[level]]
				# rownames(obj) <- gsub('unclassified', paste0('Unclassified',substr(level, 1, 1)), rownames(obj))
				rownames(obj) <- paste0(rownames(obj), substr(level, 1, 1))
				obj0[[level]] <- obj
			}
			surv.obj[[name1]] <- obj0
		}
		
	}
	
	qv.list <- surv.obj$qv.list
	hr.list <- surv.obj$hr.list
	pv.list <- surv.obj$pv.list
	
	
	hr.lc.list <- surv.obj$hr.lc.list
	hr.uc.list <- surv.obj$hr.uc.list
	
	df <- data.obj$meta.dat
	
	eff <- eff.lc <- eff.uc <- NULL
	taxa.names <- NULL
	if (is.null(taxa.levels)) {
		LOIs <- names(qv.list)
	} else {
		LOIs <- taxa.levels
		if (sum(!(taxa.levels %in% names(qv.list)))) {
			stop('Taxa levels are not contained in survival analysis results!\n')
		}
	}
	for (LOI in LOIs) {
		pv.vec <- pv.list[[LOI]]
		hr.vec <- hr.list[[LOI]]
		qv.vec <- qv.list[[LOI]]
		
		
		hr.lc.vec <- hr.lc.list[[LOI]]
		hr.uc.vec <- hr.uc.list[[LOI]]
		
		if (mt.method %in% c('fdr')) {
			taxa.name <- rownames(qv.vec)[qv.vec <= cutoff]
			taxa.name <- taxa.name[!is.na(taxa.name)]
		}
		
		if (mt.method == 'raw') {
			taxa.name <- rownames(pv.vec)[pv.vec <= cutoff]
			taxa.name <- taxa.name[!is.na(taxa.name)]
		}
		
		if (length(taxa.name) != 0) {
			
			eff <- rbind(eff, hr.vec[taxa.name, , drop=F])
			eff.lc <- rbind(eff.lc, hr.lc.vec[taxa.name, , drop=F])
			eff.uc <- rbind(eff.uc, hr.uc.vec[taxa.name, , drop=F])
			
			taxa.names <- c(taxa.names, taxa.name)
		}
	}
	
	if (length(taxa.names) == 0) {
		cat('No taxa associated with survival times! \n')
	} else {
		if (length(taxa.names) >= 2) {
			if (is.null(hei2)) {
				hei2 <- 4 + length(taxa.names) / 20 * 3 
				wid2 <- 6
			}
			rownames(eff) <- rownames(eff.lc) <- rownames(eff.uc) <- taxa.names
			for (k in 1:ncol(eff)) {
				fold.dat.plot1 <- data.frame(Estimate=eff[, k], LCI = eff.lc[, k], UCI = eff.uc[, k], IV = taxa.names)
				obj <- plot_effect_size3(fold.dat.plot1)
				pdf(paste("Taxa_Survival_HRBarplot", colnames(eff)[k], mt.method, cutoff, ann, ".pdf", sep="_"), height = hei2, width = wid2)
				print(obj)
				dev.off()
			}
			
		} 
	}
	cat('Visualization completed!\n')
}

# Rev: 2016_12_07 when the names of some differential taxa are not consistent
# Rev: 2017_02_21  Color by log fold change m.list ->> fc.list
# To-do-list:  size by effects, accomodate OTU level
create_lefse_format <- function(data.obj, diff.obj, grp.name, mt.method='fdr', cutoff=0.15, prev=0.1, minp=0.002, 
		lefse.dir="/data2/microbiome/jeff/tools/nsegata-lefse/", ann="") {
	
	# To be improved - currently no effect size
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	levels(grp) <- paste0(1:nlevels(grp), levels(grp))

	qv.list <- diff.obj$qv.list
	pv.list <- diff.obj$pv.list
	fc.list <- diff.obj$fc.list
	
	otu.name.12 <- data.obj$otu.name
	otu.tab.12 <- data.obj$otu.tab
	tax.family.a <- NULL
	for (i in 1:6) {
		tax.family <- apply(otu.name.12[, 1:i, drop=F], 1, paste, collapse=".")
		phlan.tab <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
		rownames(phlan.tab) <- phlan.tab[, 1]
		phlan.tab <- as.matrix(phlan.tab [, -1, drop=F])
		phlan.tab <- t(t(phlan.tab) / colSums(phlan.tab))	
		phlan.tab <- phlan.tab[rowSums(phlan.tab!=0) > prev*ncol(phlan.tab) & matrixStats::rowMaxs(phlan.tab) > minp, , drop=F]
		tax.family.a <- c(tax.family.a, rownames(phlan.tab))
		#tax.family.a <- c(tax.family.a, unique(tax.family))
	}
	alias.a <- sapply(strsplit(tax.family.a, "\\."), function(x) {
				if (length(x) >= 2) {
					if (length(x) == 2) {
						return(x[2])
					} else {
						return(paste0(x[2], ";", x[length(x)]))
					}
				} else {
					return(x)
				}
			})
	
	taxa.names <- NULL
	abundant.grp.names <- NULL
	
	for (LOI in setdiff(names(qv.list), 'Species')) {
		qv.vec <- qv.list[[LOI]]
		pv.vec <- pv.list[[LOI]]
		#qv.vec <- qvalue(pv.vec[, 1])$qvalues
		fc.vec <- fc.list[[LOI]]
		if (mt.method == 'fdr') {
			taxa.name <- rownames(qv.vec)[qv.vec <= cutoff]
		#	abundant.grp.name <- apply(m.vec, 1, function(x) levels(grp)[which.max(x)])[qv.vec <= cutoff]
			abundant.grp.name <- apply(fc.vec, 1, function(x) levels(grp)[as.numeric(x >= 0) + 1])[qv.vec <= cutoff]
		}

		if (mt.method == 'raw') {
			taxa.name <- rownames(pv.vec)[pv.vec <= cutoff]
		#	abundant.grp.name <- apply(m.vec, 1, function(x) levels(grp)[which.max(x)])[pv.vec <= cutoff]
		    abundant.grp.name <- apply(fc.vec, 1, function(x) levels(grp)[as.numeric(x >= 0) + 1])[pv.vec <= cutoff]
		}
		
		if (length(taxa.name) != 0) {
			taxa.names <- c(taxa.names, taxa.name)
			abundant.grp.names <- c(abundant.grp.names, abundant.grp.name)
		}
	}

	# remove 'unclassified'
    abundant.grp.names <- abundant.grp.names[!grepl('unclassified', taxa.names, ignore.case=T)]
	taxa.names <- taxa.names[!grepl('unclassified', taxa.names, ignore.case=T)]
	
	# Rev: 2016_12_07
#	taxa.names <- intersect(taxa.names, alias.a)

	ind <- match(taxa.names, alias.a)
	na.ind <- !is.na(ind)

	phlan <- cbind(tax.family.a, '1.5', '\t', '\t', '-')

	phlan[ind[na.ind], 2:5] <- cbind('3.5',  abundant.grp.names[na.ind], '3.0', '-')
	write.table(phlan,  paste0('Lefse.LDA.', ann, '.txt'), row.names=F, col.names=F, quote=F, sep='\t')
#	cmd1 <- paste0('python ', lefse.dir,  'plot_cladogram.py Lefse.LDA.txt Taxa_DifferentialAbundance_', mt.method, cutoff, '_Cladogram_', ann, '.pdf', ' --format pdf')
#	cmd2 <- paste0("python ", lefse.dir, 'plot_res.py Lefse.LDA.txt Taxa_DifferentialAbundance_', mt.method, cutoff, '_LDA_', ann, '.pdf --format pdf')
#	system(cmd1)
#	system(cmd2)
}

perform_lefse_analysis <- function (data.obj,  grp.name, sub.grp.name=NULL, prev=0.1, minp=0.002,
		lefse.dir="/data2/microbiome/jeff/tools/nsegata-lefse/", ann="") {
	otu.name.12 <- data.obj$otu.name
	otu.tab.12 <- data.obj$otu.tab
	meta.dat <- data.obj$meta.dat
	
	levels(meta.dat[, grp.name]) <- paste0(1:nlevels(meta.dat[, grp.name]), levels(meta.dat[, grp.name]))
	
	tax.family <- apply(otu.name.12, 1, paste, collapse="|")
	lefse.tab <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
	rownames(lefse.tab) <- lefse.tab[, 1]
	lefse.tab <- as.matrix(lefse.tab [, -1])

	lefse.tab <- t(t(lefse.tab) / colSums(lefse.tab))
	
	lefse.tab <- lefse.tab[matrixStats::rowMaxs(lefse.tab) > minp & rowSums(lefse.tab!=0) > prev*ncol(lefse.tab), , drop=FALSE]

	if (is.null(sub.grp.name)) {
		header <- rbind(class=as.character(meta.dat[colnames(lefse.tab), grp.name]), 
				id=colnames(lefse.tab))
	} else {
		header <- rbind(class=as.character(meta.dat[colnames(lefse.tab), grp.name]), 
				subclass=as.character(meta.dat[colnames(lefse.tab), sub.grp.name]), 
				id=colnames(lefse.tab))
	}

	
	lefse.tab <- rbind(header, lefse.tab)
	write.table(lefse.tab, "lefse.txt", sep="\t", col.names=F, quote=F)
	
#	system(paste0('mkdir LefSe_', ann))
#	output <- paste0('LefSe_', ann, '/')
#	if (is.null(sub.grp.name)) {
#		cmd1 <- paste0("python ", lefse.dir, "format_input.py lefse.txt ", output, "temp.in ", "-c 1 -u 2 -o 1000000")
#	} else {
#		cmd1 <- paste0("python ", lefse.dir, "format_input.py lefse.txt ", output, "temp.in ", "-c 1 -s 2 -u 3 -o 1000000")
#	}
#
#	cmd2 <- paste0("python ", lefse.dir, "run_lefse.py ", output, "temp.in ",  output, "lda.res")
#	cmd3 <- paste0("python ", lefse.dir, "plot_res.py ",  output, "lda.res ", output, "lda.pdf --format pdf")
#	cmd4 <- paste0("python ", lefse.dir, "plot_cladogram.py ",  output, "lda.res ", output, "cladogram.pdf --format pdf --labeled_stop_lev 6 --abrv_stop_lev  6")
#	
#	system(cmd1)
#	system(cmd2)
#	system(cmd3)
#	system(cmd4)
}

plot.Boruta2 <- function (x, colCode = c("green", "yellow", "red", "blue"), sort = TRUE,
		whichShadow = c(TRUE, TRUE, TRUE), col = NULL, xlab = "Attributes", ids=NULL,
		ylab = "Importance", ...)
{
	if (class(x) != "Boruta")
		stop("This function needs Boruta object as an argument.")
	lz <- lapply(1:ncol(x$ImpHistory), function(i) x$ImpHistory[is.finite(x$ImpHistory[, i]), i])
	names(lz) <- colnames(x$ImpHistory)
	numShadow <- sum(whichShadow)
	lz <- lz[c(rep(TRUE, length(x$finalDecision)), whichShadow)]
	col <- Boruta:::generateCol(x, colCode, col, numShadow)
	if (sort) {
		ii <- order(sapply(lz, median))
		lz <- lz[ii]
		col <- col[ii]
	}
	names(lz) <- gsub("^X",  "", names(lz))
	if (is.null(ids)) {
		len <- sum(x$finalDecision %in% c('Confirmed', 'Tentative'))
		ind <- (length(lz) - len) : length(lz)
	} else {
		ind <- match(ids, names(lz))
	}
	
	boxplot(lz[ind], xlab = xlab, ylab = ylab, col = col[ind], ...)
	invisible(x)
	names(lz[ind])
}



createROC <- function (pv.list, lab.list, pos.lab='1', file.name='ROC.pdf', width = 6, height = 6) {
	require(ROCR)
	n <- length(pv.list)
	aucs <- numeric(n)
	names(aucs) <- names(pv.list)

#	cols <- scales::hue_pal()(n)
	cols <- rep(c('red', 'blue', 'orange', 'cyan', 'purple'), ceiling(n/5))[1:n]
	ltys <- rep(c(1, 2), ceiling(n/2))[1:n]
	pdf(file.name, height=6, width=6)
	for (i in 1:n) {
		
  		    cat("*")
			pv.mat <- pv.list[[i]]
			lab.mat <- lab.list[[i]]
			
			pred <- prediction(pv.mat, lab.mat==pos.lab)
			perf <- performance(pred, "tpr", "fpr")
			aucs[i] <- mean(unlist(performance(pred, 'auc')@y.values))
			plot(perf, avg="threshold", spread.estimate = 'stddev', spread.scale = 2,  col=cols[i], lty=ltys[i], lwd=1.5,  add=ifelse(i==1, FALSE, TRUE),  main='ROC curve')
				
		}
	abline(0, 1, col='black')
	legend("right", legend=paste0(names(pv.list), "(AUC:", round(aucs, 3), ")"), col=cols, lty=ltys, lwd=2,  bty="n")
	dev.off()
}

#invlogit <- function(x) {
#	rbinom(length(x), 1, exp(x) / (1 + exp(x)))
#}
#pv.list <- list(x1=rnorm(100), x2=rnorm(100))
#lab.list <- list(x1=invlogit(pv.list[['x1']]), x2=invlogit(pv.list[['x2']]))
#createROC(pv.list, lab.list)

# Rev: 2017_02_17 supply aug.var to decision tree
# Rev: 2017_11_08 Change formula
predictionRF <- function (data.obj,  resp.name, formula=NULL, taxa.level='Species', binary=FALSE, prev=0.1, minp=0.002, B=50, seed=123, 
		boruta.level=c('Confirmed', 'Tentative'), ann='',...) {

	#sink(paste0("Taxa_RandomForest_", taxa.level, ".txt"))
	date()
	response <- data.obj$meta.dat[, resp.name]
	
	if (!is.null(formula)) {
		# adj.var <- unlist(strsplit(unlist(strsplit(formula, '\\s*~\\s*'))[2], '\\s*\\+\\s*'))
		# adj <- data.obj$meta.dat[, adj.var, drop=F]
		adj <- model.matrix(as.formula(formula), data.obj$meta.dat)
		adj.var <- colnames(adj)
	} 

	if (taxa.level == 'Species') {
		if (taxa.level %in% names(data.obj$abund.list)) {
			ct <- data.obj$abund.list[[taxa.level]]
		} else {
			# Accomodate different version
			ct <- data.obj$otu.tab
			rownames(ct) <- paste0("OTU", rownames(ct), ":", data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
			data.obj$abund.list[['Species']] <- ct
		}
	} else {
		ct <- data.obj$abund.list[[taxa.level]]
	}

	prop <- t(t(ct) / colSums(ct))
	prop <- prop[matrixStats::rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
	prop <- t(prop)
	
	if (binary == TRUE) {
		prop <- (prop != 0)
	}
	
    original.names <- colnames(prop)
	set.seed(seed)
	
	if (is.null(formula)) {
		performance <- matrix(0, nrow=B, ncol=2)
		colnames(performance) <- c("RF_M","Guess")
		roc.list <- list(RF_M=NULL)
		lab.list <- roc.list
	} else {
		performance <- matrix(0, nrow=B, ncol=4)
		colnames(performance) <- c("RF_M", "RF_CF", "RF_M+CF", "Guess")
		roc.list <- list("RF_M"=NULL, "RF_CF"=NULL, "RF_M+CF"=NULL)
		lab.list <- roc.list
	}

	
	colnames(prop) <- gsub(";", "_", colnames(prop))
	colnames(prop) <- gsub(":", "_", colnames(prop))
	colnames(prop) <- gsub("-", "_", colnames(prop))
	colnames(prop) <- gsub("\\.", "_", colnames(prop))
	names(original.names) <- colnames(prop)
	if (!is.null(formula)) {
		names(adj.var) <- adj.var
		original.names <- c(original.names, adj.var)
	}
	
	
	if (!is.null(formula)) {
		padj <- cbind(prop, adj)
	} else {
		padj <- prop
	}
	
	I <- nrow(prop)
	cat('Begin to bootstrap ...\n')
	for(b in 1:B){
		if (b %% 50 == 0) cat (".")
		err <- try({
			bsample <- sample(1:I, I, replace=T)
			if (is.factor(response)) {
				rf1 <- randomForest(x=prop[bsample, ], y=response[bsample], xtest=prop[-bsample, ], ytest=response[-bsample], sampsize=table(response[bsample]), ...)
				performance[b, "RF_M"] <- mean(response[-bsample]!= rf1$test$predicted)
				performance[b, "Guess"]<- mean(response[-bsample] != levels(response)[which.max(tabulate(response[bsample]))])
				roc.list[['RF_M']] <- c(roc.list[['RF_M']], rf1$test$vote[, 1])
				lab.list[['RF_M']] <- c(lab.list[['RF_M']], response[-bsample])
				if (!is.null(formula)) {
					rf2 <- randomForest(x=adj[bsample, ], y=response[bsample], xtest=adj[-bsample, ], ytest=response[-bsample], sampsize=table(response[bsample]), ...)
					rf3 <- randomForest(x=padj[bsample, ], y=response[bsample], xtest=padj[-bsample, ], ytest=response[-bsample], sampsize=table(response[bsample]), ...)
					performance[b, "RF_CF"] <- mean(response[-bsample]!= rf2$test$predicted)
					performance[b, "RF_M+CF"] <- mean(response[-bsample]!= rf3$test$predicted)
					roc.list[['RF_CF']] <- c(roc.list[['RF_CF']], rf2$test$vote[, 1])
					roc.list[['RF_M+CF']] <- c(roc.list[['RF_M+CF']], rf3$test$vote[, 1])
					lab.list[['RF_CF']] <- c(lab.list[['RF_CF']], response[-bsample])
					lab.list[['RF_M+CF']] <- c(lab.list[['RF_M+CF']], response[-bsample])
				}
				
			} else {
				rf1 <- randomForest(x=prop[bsample, ], y=response[bsample], xtest=prop[-bsample, ], ytest=response[-bsample], ...)
				performance[b, "RF_M"] <- mean((response[-bsample] - rf1$test$predicted)^2)
				performance[b, "Guess"]<- mean((response[-bsample] - mean(response[bsample]))^2)
				
				if (!is.null(formula)) {
					
					rf2 <- randomForest(x=adj[bsample, ], y=response[bsample], xtest=adj[-bsample, ], ytest=response[-bsample],  ...)
					rf3 <- randomForest(x=padj[bsample, ], y=response[bsample], xtest=padj[-bsample, ], ytest=response[-bsample],  ...)
					performance[b, "RF_CF"] <- mean((response[-bsample] - rf2$test$predicted)^2)
					performance[b, "RF_M+CF"] <- mean((response[-bsample]- rf3$test$predicted)^2)
				}
			}
		})
		if (inherits(err, "try-error")) {
			next
		}

	}
	cat("\n")
	
	sink(paste0("Taxa_Random_forest_FriedmanTest_", taxa.level, '_', ann, ".txt"))
	if (is.null(formula)) {
		cat("Fridman.test p value: ", friedman.test(performance)$p.value, "\n")
	} else {
		cat("Fridman.test p value (M+CF vs CF): ", friedman.test(performance[, c('RF_M+CF', 'RF_CF')])$p.value, "\n")
	}
	sink()
	
	pdf(paste0("Taxa_Random_forest_misclassification_barplot_", taxa.level, '_', ann, ".pdf"), height=5, width=4)
	if (!is.null(formula)) {
		performance2 <- performance[, c('Guess', 'RF_CF', 'RF_M', 'RF_M+CF')]
	} else {
		performance2 <- performance
	}
	if (is.factor(response)) {
		boxplot(performance2, col="#4DAF4A", ylab='Classification error', las=2)
	} else {
		boxplot(performance2, col="#4DAF4A", ylab='PMSE', las=2)
	}	
	dev.off()
	
	# ROC curve - only compare to the first level of the factor 
    if (is.factor(response)) {
		lab.list <- lapply(lab.list, function(x) {
					y <- factor(x)
					levels(y) <- levels(response)
					y
				}
		)
		createROC(roc.list, lab.list, pos.lab=levels(response)[1], 
				file.name=paste0("Taxa_Random_forest_ROC_", taxa.level, '_', ann, ".pdf"))
	}
	
	if (is.factor(response)) {
		rf <- randomForest(x=padj, y=response, sampsize=table(response), importance=T, ...)
	} else {
		rf <- randomForest(x=padj, y=response,  importance=T, ...)
	}
	
	importance <- as.data.frame(rf$importance)
	
	if (is.factor(response)) {
		write.csv(importance[rev(order(importance$MeanDecreaseAccuracy)), ],
				paste0("Taxa_RandomForest_Ranking_MeanDecreaseAccuracy_", taxa.level,'_', ann, ".csv"))
		write.csv(importance[rev(order(importance$MeanDecreaseGini)), ], 
				paste0("Taxa_RandomForest_Ranking_MeanDecreaseGini_", taxa.level, '_', ann, ".csv"))
	} else {
		write.csv(importance[rev(order(importance[, '%IncMSE'])), ],
				paste0("Taxa_RandomForest_Ranking_IncMSE_", taxa.level, '_', ann, ".csv"))
		write.csv(importance[rev(order(importance$IncNodePurity)), ], 
				paste0("Taxa_RandomForest_Ranking_IncNodePurity_", taxa.level, '_', ann, ".csv"))
	}

	
# Boruta Feature Selection - All subset selection (can't solve the confounding problem)
	obj.Boruta <- Boruta(padj, response, doTrace = 2)	
	write.csv(obj.Boruta$finalDecision, paste0("Taxa_Random_forest_Boruta_Feature_Selection_", taxa.level, '_', ann, ".csv"))
	
	pdf(paste0("Taxa_Random_forest_Boruta_Feature_Selection_", taxa.level, '_', ann,  ".pdf"), height=6, width=10)
	par(mar=par('mar') + c(3, 0, 0, 0))
	plot(obj.Boruta, main = "Feature selection by Boruta", ylab="Importance z-score", lwd = 0.5, las = 3, xlab = "",
			cex=1 / (ncol(prop)/50), cex.axis=0.25*200/ncol(prop), yaxt='n')
	axis(2, cex.axis=1)
	dev.off()
	#sink()
	
	pdf(paste0("Taxa_Random_forest_Boruta_Feature_Selection_Significant_Only_", taxa.level, '_', ann, ".pdf"), height=6, width=6)
	par(mar=par('mar') + c(3, 0, 0, 0))
	otu.ids <- plot.Boruta2(obj.Boruta, main = "Feature selection by Boruta", lwd = 0.5, las = 3, 
			ylab="Importance z-score", xlab = "", cex.axis=1, yaxt='n')
	axis(2, cex.axis=1)
	dev.off()
	
	if ('Confirmed' %in% boruta.level) {
		taxa.names <- original.names[names(obj.Boruta$finalDecision)[obj.Boruta$finalDecision %in% c('Confirmed')]] 
		if (!is.null(formula)) {
			taxa.names <- setdiff(taxa.names, adj.var)
			aug.var <- adj.var
		} else {
			aug.var <- NULL
		}

		if (length(taxa.names) > 0) {
			if (length(taxa.names) > 1) {

				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12,
						margins=c(5, 15), ann=paste0('BorutaFeatures_P_Confirmed_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, data.type='R',
						margins=c(5, 15), ann=paste0('BorutaFeatures_R_Comfirmed_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, data.type='B',
						margins=c(5, 15), ann=paste0('BorutaFeatures_B_Comfirmed_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, Colv=F, dendrogram='row', sam.ord=order(response),
						margins=c(5, 15), ann=paste0('BorutaFeatures_P_Confirmed_Unclusterded_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, data.type='R', Colv=F, dendrogram='row', sam.ord=order(response),
						margins=c(5, 15), ann=paste0('BorutaFeatures_R_Comfirmed_Unclusterded_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, data.type='B', Colv=F, dendrogram='row', sam.ord=order(response),
						margins=c(5, 15), ann=paste0('BorutaFeatures_B_Comfirmed_Unclusterded_', ann))
			}

			if (is.factor(response)) {
				generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
				generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Confirmed_', ann))
				generate_taxa_barplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
				if (length(taxa.names) > 1) {
					generate_taxa_barplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
					generate_taxa_boxplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
					build.decision.tree(data.obj,  resp.name=resp.name, aug.var=aug.var, taxa.level=taxa.level, binary=binary, taxa=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann)) 
				}
				gen <- data.obj$abund.list[[taxa.level]]
				gen <- t(t(gen) / colSums(gen))
				dat <- t(gen[taxa.names, ,drop=FALSE])
				colnames(dat) <- paste0('V', 1:ncol(dat))
				response2 <- factor(-as.numeric(response)+2)
				
				mylr <- function(formula, train, test){
					model <- randomForest(formula, data=train)
					x <- predict(model, newdata=test, type='prob')[, 2]
				}
				
				ACC <- Daim(response2 ~., model=mylr, data=as.data.frame(dat), labpos="1",
						control=Daim.control(method="boot", number=100), cutoff="0.632+")
				pdf(paste0('BorutaFeatures_Confirmed_ROC_', taxa.level, '_0.632+', ann, '.pdf'), height=6, width=6)
				plot(ACC, main='Boruta taxa', method='0.632+', lwd=1.5)	
				abline(0, 1, col='black')
				x <- ACC
				legend("bottomright", legend = paste("AUC:", 
								formatC(Daim::auc(x$roc$sens632p, x$roc$spec632p), 
										digits = max(3, getOption("digits") - 3))),
						inset = 0.01)
				dev.off()	
			} else {
				data.obj$meta.dat[, paste0(resp.name, '_b')] <- factor(data.obj$meta.dat[, resp.name] < median(data.obj$meta.dat[, resp.name]))
				levels(data.obj$meta.dat[, paste0(resp.name, '_b')]) <- c('High', 'Low')
				generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
				generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Confirmed_', ann))
				generate_taxa_barplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
				if (length(taxa.names) > 1) {
					generate_taxa_barplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
					generate_taxa_boxplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
				}
			}
		}
	}

	if ('Tentative' %in% boruta.level) {
		taxa.names <- original.names[names(obj.Boruta$finalDecision)[obj.Boruta$finalDecision %in% c('Confirmed',  'Tentative')]] 
		if (!is.null(formula)) {
			taxa.names <- setdiff(taxa.names, adj.var)
			aug.var <- adj.var
		} else {
			aug.var <- NULL
		}
		
		if (length(taxa.names) > 0) {
			
			if (length(taxa.names) > 1)  {
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12,
						margins=c(5, 15), ann=paste0('BorutaFeatures_P_Tentative_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, data.type='R',
						margins=c(5, 15), ann=paste0('BorutaFeatures_R_Tentative_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, data.type='B',
						margins=c(5, 15), ann=paste0('BorutaFeatures_B_Tentative_', ann))
				
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, Colv=F, dendrogram='row', sam.ord=order(response),
						margins=c(5, 15), ann=paste0('BorutaFeatures_P_Tentative_Unclustered_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, data.type='R', Colv=F, dendrogram='row', sam.ord=order(response),
						margins=c(5, 15), ann=paste0('BorutaFeatures_R_Tentative_Unclustered_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, data.type='B', Colv=F, dendrogram='row', sam.ord=order(response),
						margins=c(5, 15), ann=paste0('BorutaFeatures_B_Tentative_Unclustered_', ann))
			}

			if (is.factor(response)) {		
				generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
				generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Tentative_', ann))
				generate_taxa_barplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
				if (length(taxa.names) > 1) {
					generate_taxa_barplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
					generate_taxa_boxplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
					build.decision.tree(data.obj,  resp.name=resp.name, aug.var=aug.var, taxa.level=taxa.level, binary=binary, taxa=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann)) 
				}
				
				gen <- data.obj$abund.list[[taxa.level]]
				gen <- t(t(gen) / colSums(gen))
				dat <- t(gen[taxa.names, ,drop=FALSE])
  			    colnames(dat) <- paste0('V', 1:ncol(dat))
				response2 <- factor(-as.numeric(response)+2)
				
				mylr <- function(formula, train, test){
					model <- randomForest(formula, data=train)
					x <- predict(model, newdata=test, type='prob')[, 2]
				}
				
				ACC <- Daim(response2 ~., model=mylr, data=as.data.frame(dat), labpos="1",
						control=Daim.control(method="boot", number=100), cutoff="0.632+")
				pdf(paste0('BorutaFeatures_Tentative_ROC_', taxa.level, '_0.632+', ann, '.pdf'), height=6, width=6)
				plot(ACC, main='Boruta taxa', method='0.632+', lwd=2)	
				abline(0, 1, col='black')
				x <- ACC
				legend("bottomright", legend = paste("AUC:", 
								formatC(Daim::auc(x$roc$sens632p, x$roc$spec632p), 
										digits = max(3, getOption("digits") - 3))),
						inset = 0.01)
				dev.off()	
				
			} else {
				data.obj$meta.dat[, paste0(resp.name, '_b')] <- factor(data.obj$meta.dat[, resp.name] < median(data.obj$meta.dat[, resp.name]))
				levels(data.obj$meta.dat[, paste0(resp.name, '_b')]) <- c('High', 'Low')
				generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
				generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Tentative_', ann))
				generate_taxa_barplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
				if (length(taxa.names) > 1) {
					generate_taxa_barplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
					generate_taxa_boxplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
				}
			}
		}
	}	
}

col.func <- colorRampPalette(c("blue", "cyan", "red", "yellow"))
#col.fun2 <- colorpanel(9, 'yellow', 'red')

# Rev: 2017_02_19 key alignment modified, rowsep, colsep bug
# Rev: 2017_12_19 row colors fixed
generate_taxa_heatmap <- function (data.obj, taxa.levels='Genus', taxa='All', meta.info, sam.ord=NULL, data.type='P',  prev=0.1, minp=0.002, 
		row.col.dat='Phyla', phy.no=4, sepwidth=0.01, colsep=NULL, rowsep=NULL, sepcolor='black',
		white='white', colFnsC=NULL, colFnsF=NULL, Rowv=T, Colv=T, dendrogram='both', margins=c(5, 15), in.grid=F,  is.labCol=T, cexCol=1, cexRow=NULL,
		omas=c(1, 1, 1, 8), width=12, height=6, ann='All', return.obj=FALSE, ...) {
	
	colsep0 <- colsep
	rowsep0 <- rowsep
	df <- data.obj$meta.dat

	# Determine the col/rowside color
	if (is.null(colFnsC)) {
		# colFnsC <- c(colorRampPalette(c('black', 'yellow', 'red'), colorRampPalette(c('black', 'green')), colorRampPalette(c('black', 'blue'))))
		# https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
		colFnsC <- c(colorRampPalette(colors = brewer.pal(9,"RdBu")), 
				colorRampPalette(colors = brewer.pal(9,"PRGn")), colorRampPalette(colors = brewer.pal(9,"RdYlBu")),  
		        colorRampPalette(colors = brewer.pal(9,"RdGy")), colorRampPalette(colors = brewer.pal(9,"PuOr")))
	} 
	if (is.null(colFnsF)) {
#		rainbow3 <- function(x) {
#			rainbow(x + 2)[1:x]
#		}
		jet3 <- function(x) {
			jet(x + 2)[1:x]
		}
#		colFnsF <- c(rainbow3, jet3)
#		colFnsF <- c(colorRampPalette(colors = brewer.pal(8, "Set1")), colorRampPalette(colors = brewer.pal(7, "Set2")), 
#				colorRampPalette(colors = brewer.pal(7,"Dark2")), colorRampPalette(colors = jet(8)))
		colFnsF <- function (x) {
			if (x <= 6) {
				return(brewer.pal(7, "Set2")[1:x])
			} else {
				return(colorRampPalette(colors = jet(8))(x))
			}
		}
		colFnsF <- c(colFnsF)
	}
	
	
	if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
		data.obj$abund.list[['Species']] <- data.obj$otu.tab
		rownames(data.obj$abund.list[['Species']]) <- paste0("A", substr(rownames(data.obj$otu.tab), 1, 6), ":", 
				data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
	}
	for (LOI in taxa.levels) {
		cat(LOI, "\n")
		
		if (LOI == 'All') {
			if ('All' %in% taxa) {
				stop("Please specify the taxa names that will be included in the heatmap!\n")
			} 
			prop <- NULL
			for (LOI2 in names(data.obj$abund.list)) {
				ct <- data.obj$abund.list[[LOI2]]
				prop0 <- t(t(ct) / colSums(ct))
				prop <- rbind(prop, prop0[intersect(rownames(prop0), taxa), , drop=FALSE])	
	
			}
			colnames(prop) <- colnames(prop0)
			if (nrow(prop) != length(taxa)) {
				warnings('Some taxa not found in abundance lists! Please check the names!\n')
			}
			
  		} else {
			ct <- data.obj$abund.list[[LOI]]
			prop <- t(t(ct) / colSums(ct))
			
			if ('All' %in% taxa) {
				prop <- prop[taxa, , drop=FALSE]
			} else {
				prop <- prop[matrixStats::rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
			}			
		}

		# Sort row and column
		if (is.null(sam.ord)) {
			prop <- prop[order(rownames(prop)), , drop=FALSE]
		} else {
			prop <- prop[order(rownames(prop)), sam.ord, drop=FALSE]
		}
		
		if (is.labCol) {
			labCol <- colnames(prop)
		}  else {
			labCol <- ''
		}
		# Deal with zeros
		if (data.type == 'B') {
			col.scheme <- c("lightyellow", "red")
			prop[, ] <- as.numeric(prop != 0)
			breaks <- c(-0.01, 0.01, 1.1)
		} 
		if (data.type == 'P' | data.type == 'logP'){
			#col.scheme <- c(white, colFunc(50))		
			col.scheme = c(white, brewer.pal(11, "Spectral"))
	        # col.scheme = c(white, colorRampPalette(c("blue", "yellow", "red"))(11))
			ind.temp <- prop != 0
			minp <- min(prop[prop!=0])/1.1
			prop[prop==0] <- minp
			prop <- log10(prop)
		#	breaks <- c(log10(minp)-0.01, seq(log10(minp)+0.01, 0, len=51))	
	        breaks <- c(log10(minp)-0.01, seq(log10(minp)+0.01, 0, len=12))
		}
		
		if (data.type == 'sqrtP'){
			#col.scheme <- c(white, colFunc(50))		
			col.scheme <-  c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026")
			# col.scheme = c(white, colorRampPalette(c("blue", "yellow", "red"))(11))
			prop <- sqrt(prop)
			breaks <- c(-0.01, seq(min(prop[prop != 0]) * 0.99, max(prop), len=9))
		}
		if (data.type == 'R'){
			col.scheme <- c(white, colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(ncol(prop)-1))
			#col.scheme <- c('white', colorRampPalette(c("green", "black", "red"))(ncol(prop)-1))
	
			prop <- t(apply(prop, 1, function(x) {
								temp <- rank(x[x!=0])
								s <- (ncol(prop) - 1) / (max(temp) - min(temp))
								temp <- 1 + (temp - min(temp)) * s
								x[x!=0] <- temp
								x
								}))
			breaks <- seq(0, ncol(prop), len=ncol(prop)+1)
		}
		
		if (nrow(prop) > 150) {
			labRow <- ''
			cexRow <- 1
		} else {
			labRow <- rownames(prop)
			if (is.null(cexRow)) {
				cexRow <- ifelse(0.5 * 60 / nrow(prop) > 1, 1, 0.5 * 60 / nrow(prop))
			}

		}
		
		if (is.null(colsep0)) {
			if (in.grid == T) {
				colsep <- c(0:ncol(prop))  
			} else {
				colsep <- c(0, ncol(prop))  
			}
		} else {
			colsep <- colsep0
		}
		
		if (is.null(rowsep0)) {
			if (in.grid == T) {
				rowsep <- c(0:nrow(prop))
			} else {
				rowsep <- c(0, nrow(prop))
			}
		} else {
			rowsep <- rowsep0
		}
		
		key.list <- list()
		colsidecol <- NULL
		i <- 0
		j <- 0
		for (keyID in meta.info) {
			if (is.null(sam.ord)) {
				x <- df[, keyID]
			} else {
				x <- df[sam.ord, keyID]
			}		

			if (is.factor(x) | is.character(x)) {
				x <- factor(x)
				key.list[[keyID]] <- list(breaks=levels(x), colors=(colFnsF[i+1][[1]])(nlevels(x)), base=NA, col.na=NA, right=F, include.lowest=F)	
				i <- (i + 1) %% length(colFnsF)
				colsidecol <- cbind(colsidecol, key.list[[keyID]]$colors[x])
			} else {
				key.list[[keyID]] <- makecmap(x, n=5, colFn=colFnsC[j+1][[1]])
				j <- (j + 1) %% length(colFnsC)
				colsidecol <- cbind(colsidecol, cmap(x, key.list[[keyID]]))
			}
		}
		colnames(colsidecol) <- meta.info
		
		# Rev: 2016_09_13
		colsidecol[is.na(colsidecol)] <- 'white'
		
		# add aunbdance key 
		if (data.type == 'B') {
			prop.cmap <- list(breaks=c("Absence", "Presence"), colors=c("lightyellow", "red"), base=NA, col.na=NA, right=F, include.lowest=F) 
			KeyName <- 'Abundance'
		} 
#		if (data.type == 'P'){
#			prop.cmap <- makecmap(prop[ind.temp], n = 5, colFn = colFunc)
#		}
#		if (data.type == 'R'){
#			prop.cmap <- makecmap(as.vector(prop), n = 5, colFn = colFunc)
#		}
		
		if (data.type %in% c('P', 'logP')){
			KeyName <- 'log(Proportion)'	
		}
		if (data.type %in% c('sqrtP')){
			KeyName <- 'sqrt(Proportion)'	
		}
		if (data.type == 'R'){
			KeyName <- 'Rank'
		}
			
		# add row col key
        # Strong assumption of the structure of 'row' names
        if (row.col.dat == 'Phyla') {
			if (LOI %in% c('Class', 'Order', 'Family', 'Genus', 'Species')) {
				if (LOI == 'Species') {
					phy <- sapply(strsplit(rownames(prop), ":"), function(x) x[2])
					phy <- sapply(strsplit(phy, ";"), function(x) x[1])
				} else {
					phy <- sapply(strsplit(rownames(prop), ";"), function(x) x[1])
				}
				
				if (sum(is.na(phy)) == 0) {
					temp <- sort(table(phy), decr=T)
					if (length(temp) > phy.no) {
						rare.phy <- names(temp)[-(1:phy.no)]
						phy[phy %in% rare.phy] <- 'Other'
						phy <- factor(phy, levels=c(names(temp)[(1:phy.no)], 'Other'))
					} else {
						phy <- factor(phy)
					}
					rowsidecol<- rainbow(nlevels(phy))[phy]
					rowsidecol <- rbind(rowsidecol, rowsidecol)
					rownames(rowsidecol) <- c('', '')
					phy.cmap <- list(breaks=levels(phy), colors=rainbow(nlevels(phy)), base=NA, col.na=NA, right=F, include.lowest=F) 
				} else {
					rowsidecol <- NULL
				}

			} else {
				rowsidecol <- NULL
			}
		} else {
			rowsidecol <- NULL
		}
		
		pdf(paste0('Taxa_Heatmap_', LOI, '_', ann, '.pdf'), width=width, height=height)
		par(oma = omas)

#		if (data.type == 'R' | data.type == 'B') {
#			dist2 <- dist
#		}
#		if (data.type == 'P') {
#			# Better clustering of taxa
##			dist2 <- function(x) {    
##				as.dist((1-cor(t(x)))/2)
##			}
#	       dist2 <- dist
#		}
#
#		# Rev: 2016_09_13 handle zero sd cases
#        if (sum(rowSds(prop) == 0) != 0 | sum(colSds(prop) == 0) != 0)  {
#			dist2 <- dist
#			warning('Zero sd produced! Euclidean distance is used instead!\n')
#		}
        hclust2 <- function (d) {
			return(hclust(d, method = 'ward.D'))
		}
		# Pearson correlation distance
		obj <- heatmap.3(prop, 
				Rowv=Rowv, 
				Colv=Colv, 
#				distfun = dist2,
		        hclustfun = hclust2,
				dendrogram=dendrogram,
				scale='none',
				col=col.scheme, 
				breaks=breaks, 
				symbreaks=F,
				trace='none',
				margins= margins, 
				colsep = colsep,  
				rowsep = rowsep,  
				sepcolor= sepcolor, 
				sepwidth=c(sepwidth, sepwidth),
				ColSideColors=colsidecol,
				RowSideColors=rowsidecol,
				cexRow=cexRow,
				labRow=labRow,
				labCol=labCol,
				cexCol=cexCol,
				key=(data.type != 'B'), density.info='none', symkey=F, KeyValueName=KeyName,	
				NumColSideColors= 0.5 *length(meta.info),
				NumRowSideColors= 0.5,
				...
		)

		par(cex=0.75)
		par(oma=c(0, 0, 1, 0)) 
		
		if (!is.null(rowsidecol) & row.col.dat == 'Phyla') {
			y.cord <- (1/(length(meta.info)+2)) * (0:((length(meta.info)+2) - 1))
			vkey2(phy.cmap, 'Phylum', x=0, y=-0.2, stretch=1.2)

		}

		y.cord <- (1/(length(meta.info))) * (0:(length(meta.info) - 1))
		k <- 1
		xs <- c(0.95, 1)
		for (keyID in meta.info) {
			x <- df[, keyID]
			if (is.factor(x) | is.character(x)) {
				x <- factor(x)
				vkey2(key.list[[keyID]], keyID, x=xs[(k-1) %% 2 + 1], y=y.cord[k], stretch=1.2)
			} else {
				vkey(key.list[[keyID]], keyID, x=xs[(k-1) %% 2 + 1], y=y.cord[k], stretch=1.2)
			}
			k <- k + 1
		}	
		
		if (data.type == 'B') {
			vkey2(prop.cmap, KeyName,  x=xs[(k-1) %% 2 + 1], y=1, stretch=1.2)		
		} 
#		if (data.type == 'P'){
#			vkey(prop.cmap, 'log10(Proportion)', y=1, stretch=1.2)
#		}	
#		if (data.type == 'R'){
#			vkey(prop.cmap, 'Rank', y=1, stretch=1.2)
#		}	
		dev.off()	
	}
	if (return.obj == TRUE) {
		return(obj)
	}
	
}

# New: 2017_03_02  General heatmap without taxa
# data.obj$data, data.obj$meta.dat
# Rev: 2017_11_20, Add 'LogC', taxa.as.row
generate_heatmap <- function (data.obj,  meta.info, sam.ord=NULL, data.type='P',
		row.col.dat=NULL,  sepwidth=0.01, colsep=NULL, rowsep=NULL, taxa.as.row = TRUE,
		colFnsC=NULL, colFnsF=NULL, Rowv=T, Colv=T, dendrogram='both', margins=c(5, 15),
		in.grid=F, sepcolor='black', is.labCol=T, cexCol=1, cexRow=NULL,
		omas=c(1, 1, 1, 8), width=12, height=6, ann='All', return.obj=FALSE, pdf=TRUE, ...) {
	
	colsep0 <- colsep
	rowsep0 <- rowsep
	df <- data.obj$meta.dat
	
	# Determine the col/rowside color
	# Determine the col/rowside color
	if (is.null(colFnsC)) {
		# colFnsC <- c(colorRampPalette(c('black', 'yellow', 'red'), colorRampPalette(c('black', 'green')), colorRampPalette(c('black', 'blue'))))
		# https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
		colFunsC <- c(colorRampPalette(colors = brewer.pal(9,"RdBu")), 
				colorRampPalette(colors = brewer.pal(9,"PRGn")), colorRampPalette(colors = brewer.pal(9,"RdYlBu")),  
				colorRampPalette(colors = brewer.pal(9,"RdGy")), colorRampPalette(colors = brewer.pal(9,"PuOr")))
	} 
	if (is.null(colFnsF)) {

		colFnsF <- function (x) {
			if (x <= 6) {
				return(brewer.pal(7, "Set2")[1:x])
			} else {
				return(colorRampPalette(colors = jet(8))(x))
			}
		}
		colFnsF <- c(colFnsF)
	}
	prop <- data.obj$data
	
	# Sort row and column
	if (!is.null(sam.ord)) {
		prop <- prop[, sam.ord, drop=FALSE]
	}
	
	if (is.labCol) {
		labCol <- colnames(prop)
	}  else {
		labCol <- ''
	}
	if (data.type == 'LogC') {
		col.scheme <- c('white', colorpanel(9, 'aliceblue', 'lightblue', 'blue'))
		breaks <- c(-0.1, 0.1, 2, 4, 6, 8, 10, 12, 14, 16, 18)
	}
	# Deal with zeros
	if (data.type == 'B') {
		col.scheme <- c("lightyellow", "red")
		breaks <- c(-0.01, 0.01, 1.1)
	} 
	if (data.type == 'P'){
		#
		#col.scheme = brewer.pal(11, "Spectral")
	     col.scheme <- bluered(11)
		breaks <- c(seq(min(prop), max(prop), len=12))
	}
	if (data.type == 'R'){
		col.scheme <- c(white, colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(ncol(prop)-1))
		breaks <- seq(0, ncol(prop), len=ncol(prop)+1)
	}
	
	if (nrow(prop) > 1500) {
		labRow <- ''
		cexRow <- 1
	} else {
		labRow <- rownames(prop)
		if (is.null(cexRow)) {
			cexRow <- ifelse(0.5 * 60 / nrow(prop) > 1, 1, 0.5 * 60 / nrow(prop))
		}
		
	}
	
	if (is.null(colsep0)) {
		if (in.grid == T) {
			colsep <- c(0:ncol(prop))  
		} else {
			colsep <- c(0, ncol(prop))  
		}
	} else {
		colsep <- colsep0
	}
	
	if (is.null(rowsep0)) {
		if (in.grid == T) {
			rowsep <- c(0:nrow(prop))
		} else {
			rowsep <- c(0, nrow(prop))
		}
	} else {
		rowsep <- rowsep0
	}
	
	key.list <- list()
	colsidecol <- NULL

	i <- 0
	j <- 0
	
	for (keyID in meta.info) {
		if (is.null(sam.ord)) {
			x <- df[, keyID]
		} else {
			x <- df[sam.ord, keyID]
		}		
		if (is.factor(x)) {
			key.list[[keyID]] <- list(breaks=levels(x), colors=(colFnsF[i+1][[1]])(nlevels(x)), base=NA, col.na=NA, right=F, include.lowest=F)	
			i <- (i + 1) %% length(colFnsF)
			colsidecol <- cbind(colsidecol, key.list[[keyID]]$colors[x])
		} else {
			key.list[[keyID]] <- makecmap(x, n=5, colFn=colFnsC[j+1][[1]])
			j <- (j + 1) %% length(colFnsC)
			colsidecol <- cbind(colsidecol, cmap(x, key.list[[keyID]]))
		}
	}
	colnames(colsidecol) <- meta.info
	
	# Rev: 2016_09_13
	colsidecol[is.na(colsidecol)] <- 'white'
	
	# add aunbdance key 
	if (data.type == 'B') {
		prop.cmap <- list(breaks=c("0", "1"), colors=c("lightyellow", "red"), base=NA, col.na=NA, right=F, include.lowest=F) 
		KeyName <- 'Value'
	} 
	
	if (data.type == 'P'){
		KeyName <- 'Value'	
	}
	if (data.type == 'R'){
		KeyName <- 'Rank'
	}
	
	if (data.type == 'LogC'){
		KeyName <- 'Log2(count+1)'
	}
	
	# add row col key
	if (!is.null(row.col.dat)) {
		
		phy <- df[, row.col.dat]
		phy <- factor(phy)
		
		rowsidecol<- rainbow(nlevels(phy))[phy]
		rowsidecol <- rbind(rowsidecol, rowsidecol)
		rownames(rowsidecol) <- c('', '')
		phy.cmap <- list(breaks=levels(phy), colors=rainbow(nlevels(phy)), base=NA, col.na=NA, right=F, include.lowest=F) 
	} else {
		rowsidecol <- NULL
	}
	if(pdf == TRUE) {
		pdf(paste0('Heatmap_', ann, '.pdf'), width=width, height=height)
	}
	
	par(oma = omas)
	
	
	# Pearson correlation distance
    if (taxa.as.row) {
		prop <- prop
	} else {
		prop <- t(prop)
	}
	obj <- heatmap.3(prop, 
			Rowv=Rowv, 
			Colv=Colv, 
#				distfun = dist2,
			dendrogram=dendrogram,
			scale='none',
			col=col.scheme, 
			breaks=breaks, 
			symbreaks=F,
			trace='none',
			margins= margins, 
			colsep = colsep,  
			rowsep = rowsep,  
			sepcolor= sepcolor, 
			sepwidth=c(sepwidth, sepwidth),
			ColSideColors=colsidecol,
			RowSideColors=rowsidecol,
			cexRow=cexRow,
			labRow=labRow,
			labCol=labCol,
			cexCol=cexCol,
			key=(data.type != 'B'), density.info='none', symkey=F, KeyValueName=KeyName,	
			NumColSideColors= 0.5 *length(meta.info),
			NumRowSideColors= 0.5,
			...
	)
	
	par(cex=0.75)
	par(oma=c(0, 0, 1, 0)) 
	
	if (!is.null(rowsidecol)) {
		y.cord <- (1/(length(meta.info)+2)) * (0:((length(meta.info)+2) - 1))
		vkey2(phy.cmap, rowsidecol, x=0, y=-0.2, stretch=1.2)
		
	}
	
	y.cord <- (1/(length(meta.info))) * (0:(length(meta.info) - 1))
	k <- 1
	xs <- c(0.95, 1)
	for (keyID in meta.info) {
		x <- df[, keyID]
		if (is.factor(x)) {
			vkey2(key.list[[keyID]], keyID, x=xs[(k-1) %% 2 + 1], y=y.cord[k], stretch=1.2)
		} else {
			vkey(key.list[[keyID]], keyID, x=xs[(k-1) %% 2 + 1], y=y.cord[k], stretch=1.2)
		}
		k <- k + 1
	}	
	
	if (data.type == 'B') {
		vkey2(prop.cmap, KeyName, x=xs[(k-1) %% 2 + 1], y=1, stretch=1.2)		
	} 
	
	if (pdf == TRUE) {
		dev.off()	
	}
	
	
	if (return.obj == TRUE) {
		return(obj)
	}
}

lighter <- function(color, factor = 0.2) {
	## converts color to hsv, multiplies v by factor, returns colors as hexcode
	x = rgb2hsv(col2rgb(color))
	v = pmax(pmin(x[3, ] + factor, 1), 0)
	hsv(h = x[1, ], s = x[2, ], v = v)
}

# Rev: 2016_06_13
# Still need to be revised: Multiple groups, the auto order should look similar for differnt groups;
# The legends do not work
# Rev: 2017_04_28 add col.map so the color scheme can be the same for different data sets
# Rev: 2017_11_09 add label for individual bar plot and change ordering 
# Rev: 2018_03_15 add par(lwd..)
# Rev: 2018_05_31 add return col.map and aggregate the taxa not in col.map into 'Other' Class
# Rev: 2018_06_07 many colors
generate_stacked_barplot <- function(data.obj, grp.name=NULL, taxa.levels=c('Phylum', 'Family', 'Genus'), agg.cutoff=0.005, 
		border=TRUE, order.auto=FALSE, order.method='abund', order.taxa=FALSE, cex.names=0.5, cex.top.lab = 0.75, separate=FALSE, cex.names2=1, indiv=TRUE, aggre=TRUE,
		hei1=6, wid1=9, hei2=6, wid2=9, margin=10, ann='', pdf=TRUE, col.map = NULL, many.colors=FALSE, ...) {
	if (is.null(grp.name)) {
		grp <- rep(1, nrow(data.obj$meta.dat))
	} else {
		grp <- factor(data.obj$meta.dat[, grp.name])
	}
	
	name.list <- list()
	col.list <- list()
	

	lwd.o <- par('lwd')
	# rearrange the order of taxa
	if (order.taxa) {
		for (taxa.level in taxa.levels) {
			abund0 <- data.obj$abund.list[[taxa.level]]
			abund0 <- t(t(abund0) / colSums(abund0))
			data.obj$abund.list[[taxa.level]] <- data.obj$abund.list[[taxa.level]][rev(order(rowMeans(abund0))), ]
		}
	}
	
	# create color mapping
	prop.list <- list()
	for (taxa.level in taxa.levels) {
		abund0 <- data.obj$abund.list[[taxa.level]]
		abund0 <- t(t(abund0) / colSums(abund0))
		
		abund1 <- abund0[rowMeans(abund0) >= agg.cutoff, , drop=F]
		abund2 <- abund0[rowMeans(abund0) < agg.cutoff, , drop=F]
		
		prop <- rbind(abund1, Other=colSums(abund2))
		colnames(prop) <- colnames(abund0)
#		name.list[[taxa.level]] <- rownames(abund1)
		prop.list[[taxa.level]] <- prop
		
		#rand.col <- c(sample(rainbow(nrow(prop)*2), nrow(prop)-1), 'gray')
		# Rev: 2017_04_28
		if (is.null(col.map[[taxa.level]])) {
			if (many.colors) {
				basic <- brewer.pal(12, "Paired")
				rand.col <- c(rep_len(c(lighter(basic, -0.2), lighter(basic, -0.1), basic, lighter(basic, 0.1),  lighter(basic, 0.2)), nrow(prop) - 1), 'gray')
			} else {
				rand.col <- c(rep_len(brewer.pal(12, "Paired"), nrow(prop) - 1), 'gray')
			}
			
			# rand.col <- c(brewer.pal(nrow(prop) - 1, "Paired"), 'gray')
			names(rand.col) <- rownames(prop)
			col.list[[taxa.level]] <- rand.col
		} else {
			# Rev: 2018_05_30
			if (sum(!(rownames(prop) %in% names(col.map[[taxa.level]]) ))) {
				warning('Color mapping does not contain all the taxa! Unmapped taxa will be aggreated into "Other" category!\n')
				umapped.taxa <- setdiff(rownames(prop), names(col.map[[taxa.level]]))
				Other2 <- colSums(prop[umapped.taxa, , drop=FALSE])
				prop <- prop[setdiff(rownames(prop), umapped.taxa), , drop=FALSE]
				if ('Other' %in% rownames(prop)) {
					prop['Other', ] <- prop['Other', ] + Other2
				} else {
					prop <- rbind(prop, Other = Other2)
					colnames(prop) <- colnames(data.obj$abund.list[[taxa.level]])
				}
				#		name.list[[taxa.level]] <- rownames(prop)
				prop.list[[taxa.level]] <- prop
			} 
			col.list[[taxa.level]] <- col.map[[taxa.level]]
		}
		
	}	
	
	
	if (indiv == TRUE) {
		if (pdf == TRUE)
			pdf(paste0("Taxa_Stacked_Barplot_Overall_Compo_", ann, ".pdf"), height=hei1, width=wid1)
		
		if (border == FALSE) {
			lty.o <- par("lty")
			par(lty = 0)
		}
		
		mar.o <- par(mar=par('mar') + c(0, margin, 0, 0))
		
		for (taxa.level in taxa.levels) {
#			abund0 <- data.obj$abund.list[[taxa.level]]
#			abund0 <- t(t(abund0) / colSums(abund0))
#			
			##			abund1 <- abund0[rowMeans(abund0) >= agg.cutoff, , drop=F]
			##			abund2 <- abund0[rowMeans(abund0) < agg.cutoff, , drop=F]
#			
#			abund1 <- abund0[name.list[[taxa.level]], , drop=FALSE]
#			abund2 <- 1 - colSums(abund1)
#			
#			prop <- rbind(abund1, Other=abund2)
#			colnames(prop) <- colnames(abund0)
			
			prop <- prop.list[[taxa.level]]
			
			#rand.col <- c(sample(rainbow(nrow(prop)*2), nrow(prop)-1), 'gray')
			#col.list[[taxa.level]] <- rand.col
			
			cex.legend = ifelse (nrow(prop) > 35, 35/nrow(prop)*0.75, 0.75)
			
			# Rev: 2016_09_13, better ordering within group based on single linkage
			# Rev: 2017_11_10, ordering rule the same across group
			if (order.auto) {
				
				temp <- rev(order(rowMeans(prop)))
				
				prop <- do.call(cbind, tapply(1:length(grp), factor(grp), function(i){
									if (length(i) >= 2) {
										prop.sub <- prop[, i, drop=FALSE]
										if (order.method == 'single') {
											ord <- hclust(vegdist(t(prop.sub)), method='single')$order
										} else {
											prop.dummy <- prop.sub
											prop.dummy[prop.dummy < 0.05] <- 0
											if (length(temp) >=3) {
												ord <- (order(-prop.dummy[temp[1], ], -prop.dummy[temp[2], ], -prop.dummy[temp[3], ]))
											} else {
												ord <- (order(-prop.dummy[temp[1], ]))
											}
										}
										# temp <- rev(order(rowMeans(prop.sub)))
										prop.sub <- prop.sub[, ord, drop=FALSE]
										return(prop.sub)
									} else {
										return(prop[, i, drop=FALSE])
									}
								}))
			} else {
				prop <- prop[, order(grp), drop=FALSE]
			}
			
			grp2 <- sort(grp)
			if (border == FALSE) {
				if (is.null(grp.name)) {
					par(lwd = 0.25)
					barplot(prop, col=col.list[[taxa.level]][rownames(prop)], ylab='Proportion', las=2, legend.text=rownames(prop), cex.names=cex.names, space=0,
							args.legend=list(x='left', bty='n',  cex=cex.legend, inset=c(-0.5, 0)), main=taxa.level, ...)
				}  else {
					
					par(lwd = 0.25)
					coord <- barplot(prop, ylim = c(0, 1.05), col=col.list[[taxa.level]][rownames(prop)], ylab='Proportion', las=2, legend.text=rownames(prop), cex.names=cex.names, space=0,
							args.legend=list(x='left', bty='n',  cex=cex.legend, inset=c(-0.5, 0)), main=taxa.level, ...)
					
					coord.text <- tapply(coord, grp2, mean)
					coord.st <- tapply(coord, grp2, max)[-(length(unique(grp2)))]
					coord.end <- tapply(coord, grp2, min)[-1]
					coord.line <- (coord.st +  coord.end) / 2
					
					text.names <- names(coord.text)
					
					
					for (i in 1 : length(coord.line)) {
						abline(v = coord.line[i], lty = 1, lwd = 0.5)	
					}
					for (i in 1 : length(coord.text)) {
						if (i %% 2 == 0) {
							pos <- 3
						} else {
							pos <- 3
						}
						text(coord.text[i], 1, text.names[i], pos = pos, cex = cex.top.lab)
					}
					
				}
			} else {
				if (is.null(grp.name)) {
					par(lwd = 0.25)
					barplot(prop, col=col.list[[taxa.level]][rownames(prop)], ylab='Proportion', las=2, legend.text=rownames(prop), cex.names=cex.names,
							args.legend=list(x='left', bty='n',  cex=cex.legend, inset=c(-0.5, 0)), main=taxa.level, ...)
				} else {
					par(lwd = 0.25)
					coord <- barplot(prop, ylim = c(0, 1.05), col=col.list[[taxa.level]][rownames(prop)], ylab='Proportion', las=2, legend.text=rownames(prop), cex.names=cex.names,
							args.legend=list(x='left', bty='n',  cex=cex.legend, inset=c(-0.5, 0)), main=taxa.level, ...)
					
					coord.text <- tapply(coord, grp2, mean)
					coord.st <- tapply(coord, grp2, max)[-(length(unique(grp2)))]
					coord.end <- tapply(coord, grp2, min)[-1]
					coord.line <- (coord.st +  coord.end) / 2
					
					text.names <- names(coord.text)
					for (i in 1 : length(coord.line)) {
						abline(v = coord.line[i], lty = 1, lwd = 0.5)	
					}
					for (i in 1 : length(coord.text)) {
						if (i %% 2 == 0) {
							pos <- 3
						} else {
							pos <- 3
						}
						text(coord.text[i], 1, text.names[i], pos = pos, cex = cex.top.lab)
					}
				}
				
				
			}
			
		}
		if (border == FALSE) {
			par(lty = lty.o)
		}
		
		par(mar = mar.o)
		if (pdf == TRUE)
			dev.off()
		
	}
	
	# Generate averaged over stack barplot
	
	if (!is.null(grp.name) & aggre == TRUE) {
		if (separate != TRUE) {
			if (pdf == TRUE)
				pdf(paste0("Taxa_Stacked_Barplot_Grouped_Compo_", ann, "_combine.pdf"), height=hei2, width=wid2)
			mar.o <- par(mar=par('mar') + c(0, margin/length(taxa.levels), 0, 0))
			mfrow.o <- par(mfrow=c(1, length(taxa.levels)))
			cex.legend <- 10
			for (taxa.level in taxa.levels) {
#				abund0 <- data.obj$abund.list[[taxa.level]]
#				abund0 <- t(t(abund0) / colSums(abund0))
#			
#				
#			    temp <- aggregate(t(abund0), by=list(grp), mean)
#				abund0 <- as.matrix(temp[, -1])
#				rownames(abund0) <- temp[, 1]
#				abund0 <- t(abund0)
#				
#				abund1 <- abund0[name.list[[taxa.level]], , drop=FALSE]
#				abund2 <- 1 - colSums(abund1)
#				
#				prop <- rbind(abund1, Other=abund2)
#				colnames(prop) <- colnames(abund0)
				
				temp <- aggregate(t(prop.list[[taxa.level]]), by=list(grp), mean)
				prop <- as.matrix(temp[, -1])
				rownames(prop) <- temp[, 1]
				prop <- t(prop)
				
				newsize <- ifelse (nrow(prop) > 35, 35/nrow(prop)*0.75, 0.75)
				cex.legend <- ifelse(cex.legend < newsize, cex.legend, newsize) 
				#prop <- prop[, order(grp)]
				par(lwd = 0.25)
				barplot(prop,  col=col.list[[taxa.level]][rownames(prop)], ylab='Proportion', las=2, cex.names=cex.names2, main=taxa.level, ...)
				#		legend.text=rownames(prop), args.legend=list(x='left', bty='n',  cex=cex.legend, inset=c(-2.2, 0)))
				
			}
			
			par(mar=c(0, 0, 0, 0))
			oma.o <- par(oma=c(0, 0, 0, 0))
			for (taxa.level in taxa.levels) {
#				abund0 <- data.obj$abund.list[[taxa.level]]
#				abund0 <- t(t(abund0) / colSums(abund0))
#				
#				abund0 <- t(apply(abund0, 1, function(x) {
#									tapply(x, grp, mean)
#								}))
#				
#				abund1 <- abund0[name.list[[taxa.level]], , drop=FALSE]
#				abund2 <- 1 - colSums(abund1)
#				
#				prop <- rbind(abund1, Other=abund2)
#				colnames(prop) <- colnames(abund0)
				
				
				#	cex.legend = ifelse (nrow(prop) > 35, 35/nrow(prop)*0.75, 0.75)
				plot(1, type="n", axes=FALSE, xlab="", ylab="")
				# Rev: 2016_09_10
				legend('left', legend=rev(rownames(prop.list[[taxa.level]])), bty='n',
						fill=rev(col.list[[taxa.level]][rownames(prop.list[[taxa.level]])]), cex=cex.legend)
				
			}
			par(mar=mar.o)
			par(oma=oma.o)
			par(mfrow=mfrow.o)
			if (pdf == TRUE)
				dev.off()
		} else {
			if (pdf == TRUE)
				pdf(paste0("Taxa_Stacked_Barplot_Grouped_Compo_", ann, "_Separate.pdf"), height=hei2, width=wid2)
			
			mar.o <- par(mar=par('mar') + c(0, margin, 0, 0))
			for (taxa.level in taxa.levels) {
#				abund0 <- data.obj$abund.list[[taxa.level]]
#				abund0 <- t(t(abund0) / colSums(abund0))
#				
#				abund0 <- t(apply(abund0, 1, function(x) {
#									tapply(x, grp, mean)
#								}))
#				
#				abund1 <- abund0[name.list[[taxa.level]], , drop=FALSE]
#				abund2 <- 1 - colSums(abund1)
#				
#				prop <- rbind(abund1, Other=abund2)
#				colnames(prop) <- colnames(abund0)
				temp <- aggregate(t(prop.list[[taxa.level]]), by=list(grp), mean)
				prop <- as.matrix(temp[, -1])
				rownames(prop) <- temp[, 1]
				prop <- t(prop)
				
				cex.legend <- ifelse (nrow(prop) > 35, 35/nrow(prop)*0.75, 0.75)
				par(lwd = 0.25)
				barplot(prop, col=col.list[[taxa.level]][rownames(prop)], ylab='Proportion', las=2, cex.names=cex.names2, main=taxa.level,
						legend.text=rownames(prop), args.legend=list(x='left', bty='n',  cex=cex.legend, inset=c(-0.5, 0)), ...)
				
			}
			par(mar=mar.o)
			if (pdf == TRUE)
				dev.off()
			
		}
		
	}
	
	par(lwd = lwd.o)
	cat('Completed!\n')
	return(invisible(list(col.map=col.list)))
}

build.decision.tree <- function(data.obj,  resp.name, taxa.level='Species', binary=FALSE, taxa, aug.var=NULL, ann='All') {
	ann <- paste(taxa.level, ann, sep="_")
	response <- data.obj$meta.dat[, resp.name]
	
	ct <- data.obj$abund.list[[taxa.level]]
	prop <- t(t(ct) / colSums(ct))
	prop <- prop[taxa, , drop=F] 	
	if (binary == TRUE) {
		prop <- (prop != 0)
	}
	
	dat <- as.data.frame(t(prop))
	# Rev: 2017_02_17 Add additional variables from meta dat
	if (!is.null(aug.var)) {
		dat <- cbind(dat, data.obj$meta.dat[, aug.var])
	}
	
	dat <- data.frame(dat, response)
	try(
		if (is.factor(response)) {
			fit <- rpart(response ~ ., method="class", data=dat)
			post(fit, file = paste0("Taxa_Unpruned_Classification_tree_", ann, ".ps"), title = "Unpruned Classification Tree")
			pfit<- prune(fit, cp= fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
			post(pfit, file = paste0("Taxa_Pruned_Classification_", ann, ".ps"), title = "Pruned Classification Tree")
			
		} else {
			fit <- rpart(response ~ ., method="anova", data=dat)
			post(fit, file = paste0("Taxa_Unpruned_Regression_tree_", ann, ".ps"), title = "Unpruned Regression Tree")
			pfit<- prune(fit, cp= fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
			post(pfit, file = paste0("Taxa_Pruned_Regression_tree_", ann, ".ps"), title = "Pruned RegressionTree")
		}
	)	
}

bootstrap.pwr <- function(pcs, formula, dat, ns=NULL, perm.no=199, iter.no=100) {
	if (is.null(ns)) {
		ns <- nrow(dat)*seq(1, 4, len=8)
	}
	pvs <- numeric(length(ns))
	names(pvs) <- paste(ns)
	for (n in ns) {
		cat('.')
		temp <- sapply(1:iter.no, function(i) {
			bt.ind <- sample(1:nrow(dat), n, repl=T)
			dat.bt <- dat[bt.ind, ]
			dist.bt <- dist(pcs[bt.ind, ])
			aov.tab <- adonis(as.formula(paste('dist.bt', formula)), dat=dat.bt, permutations=perm.no)$aov.tab
	        pv <- aov.tab[nrow(aov.tab)-2, ncol(aov.tab)]
		})
        pvs[paste(n)] <- mean(temp <= 0.05)
	}

	return(pvs)
}


perform_power_analysis <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), formula=NULL,  
		grp.name=NULL, adj.name=NULL, ann='', ...) {
	if (is.null(formula)) {
		if (is.null(adj.name)) {
			formula <- paste('~', grp.name)
		} else {
			formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
		}
	}	 
	df <- data.obj$meta.dat
    pvm <- NULL
	for (dist.name in dist.names) {
		cat('*')
		dist.mat <- dist.obj[[dist.name]]
		pcs <- cmdscale(dist.mat, k=nrow(dist.mat)-1)
		pvs <- bootstrap.pwr(pcs, formula, df, ...)
		pvm <- rbind(pvm, pvs)
	}
	colnames(pvm) <- names(pvs)
	rownames(pvm) <- dist.names
	
	pvdf <- melt(pvm)
	colnames(pvdf) <- c('Distance_type', 'Sample_size', 'Value')
	pdf(paste0("Power_curve_bootstrap_", ann, '.pdf')) 
	g.obj <- ggplot(pvdf, aes(x=Sample_size, y=Value)) +
				geom_point() +
				geom_line() +
				ylab('Power') +
				facet_wrap(~ Distance_type, ncol=2) +
				theme_bw()
	print(g.obj)
	dev.off()
	
	return(pvdf)
}

# More robust identification
k.nonoverlap.se <- function (tab, SE.factor=1) {
	max.v <- which.max(tab[, 'gap']) 
	
	for (i in (max.v-1):1) {
		if (tab[i, 'gap'] + SE.factor * tab[i, 'SE.sim'] < tab[i+1, 'gap'] - SE.factor * tab[i+1, 'SE.sim']) {
			break
		} 
	}
	
	return(i+1)
}


perform_cluster_analysis <- function (data.obj, dist.obj, dist.name=c('UniFrac'), k.best=NULL, method='pam', stat='gap', 
		grp.name=NULL, adj.name=NULL, subject=NULL, ann='', seed=1234) {
	
	df <- data.obj$meta.dat
	if (!is.null(grp.name)) {
		grp <- df[, grp.name]
	}
	if (!is.null(adj.name)) {
		adj <- df[, adj.name]
	}
	cat(dist.name, ' distance ...\n')
	mat <- dist.obj[[dist.name]]
	
	if (is.null(k.best)) {
		pc.obj <- cmdscale(as.dist(mat), k=ncol(mat) - 1, eig=T)
		eig <- pc.obj$eig
		eig <- eig[eig > 0]
		pvar <- eig / sum(eig)
		pc <- pc.obj$points
		
		cat('Assess cluster number by gap statistics ...\n')
		
		set.seed(seed)
		gs <-  gapstat_ord(pc, axes=1:ncol(pc), verbose=FALSE)
		plot_clusgap(gs)
		ggsave(paste0('cluster_assess_gap_statistic_', dist.name, ann, '.pdf') , width=6, height=6)
		
		
#		print(gs, method="Tibs2001SEmax")
#		print(gs, method="firstSEmax")
		
		tab <- gs$Tab
		write.csv(tab, paste0('gap_stat_', dist.name, '.csv'))
		
		k.best.gap <- k.nonoverlap.se(tab)
		
		cat('Gap statistic finds ', k.best.gap, ' clusters.\n')
		
		
		cat('Assess cluster number by average silhouette width ...\n')
		asw <- numeric(20)
		for (k in 2:20){
			asw[k] <- pam(as.dist(mat), k)$silinfo$avg.width
		}
		
		k.best.asw <- which.max(asw)
		cat("silhouette-optimal number of clusters:", k.best.asw, "\n")
		pdf(paste0('cluster_assess_asw_statistic_', dist.name, ann, '.pdf'), width=6, height=6)
		plot(1:20, asw, type= "h", main = "pam() clustering assessment",
				xlab= "k (# clusters)", ylab = "average silhouette width")
		axis(1, k.best.asw, paste("best",k.best.asw, sep="\n"), col = "red", col.axis = "red")
		dev.off()
		
		cat('ASW statistic finds ', k.best.asw, ' clusters.\n')
		
		if (stat == 'gap') {
			k.best <- k.best.gap
		} else {
			if (stat == 'asw') {
				k.best <- k.best.asw
			}
		}
		
	}
	
	if (k.best == 1) {
		cat('No robust cluster structure found!\n')
	} else {
		pam.obj <- pam(as.dist(mat), k=k.best)
		pam.class <- factor(pam.obj$clustering)
		write.csv(pam.class, paste0('Cluster.membership', dist.name, '.csv'))
		
		df$Cluster <- pam.class
		data.obj$meta.dat <- df
		if (is.null(grp.name)) {
			generate_ordination(data.obj, dist.obj, dist.name, grp.name='Cluster', ann=paste0('Cluster', dist.name, ann))
		} else {
			generate_ordination(data.obj, dist.obj, dist.name, grp.name='Cluster', strata=grp.name, ann=paste0('Cluster', dist.name, ann))
		}
		
		# Define cluster characteristics
	   diff.obj <- perform_differential_analysis(data.obj, grp.name='Cluster', 
			   taxa.levels=c('Genus'),
			   method='kruskal', mt.method='fdr', 
			   cutoff=0.01, prev=0.1, minp=0.002, ann='Cluster')
	   
	   visualize_differential_analysis(data.obj, diff.obj, grp.name='Cluster', taxa.levels=c('Genus'), indivplot=FALSE,
			   mt.method='fdr', cutoff=0.01, ann='Cluster')
		
		try(
				if (!is.null(grp.name)) {
					if (is.null(adj.name)) {
						form <- as.formula(paste('yy ~', grp.name))
					} else {
						form <- as.formula(paste('yy ~', paste(adj.name, collapse='+'), '+', grp.name))		
					}
					# Enrichment analysis - logistic regression - 1 vs other
					cat('Testing for association with the clusters - 1 vs other ...\n')
					sink(paste0('Cluster_association_test', dist.name, ann, '.txt'))
					if (is.null(subject)) {
						cat('Generalized linear model (logistic regression) is performed.\n')
						
						for (clus in levels(df$Cluster)[1:(nlevels(df$Cluster))]) {	
							cat('Test for enrichment in cluster', clus, '\n')
							y <- rep(0, length(df$Cluster)) 
							y[df$Cluster == clus] <- 1
							df$yy <- y
							prmatrix(summary(glm(form, data=df, family=binomial))$coefficients)	
						}
					} else {
						cat('Generalized linear mixed effects model (logistic regression) is performed.\n')
						for (clus in levels(df$Cluster)[1:(nlevels(df$Cluster))]) {
							cat('Test for enrichment in cluster', clus, '\n')
							y <- rep(0, length(df$Cluster)) 
							y[df$Cluster == clus] <- 1
							df$yy <- y
							prmatrix(summary(glmmPQL(form, data=df, random = as.formula(paste0('~ 1|', subject)), family=binomial, verbose=F))$tTable)
							
						}
					}
					sink()
				}
		)
		
	}
}


##########
B2M <- function(x) {
	x[x==0] <- min(x[x!=0])
	x[x==1] <- max(x[x!=1])
	log(x / (1-x))
}


fastDist <- function(X) {
	temp <- colSums(X^2)  
	D <- outer(temp, temp, "+") - 2 * t(X) %*% X
	diag(D) <- 0
	sqrt(D)
}

fastLM <- function(Y, M) {
	Y <- as.matrix(Y)
	XXI <- solve(t(M) %*% M)
	dof <- ncol(Y) - ncol(M)
	est <- XXI %*% t(M) %*% t(Y)
	resid <- t(Y) - M %*% est
	sigma <- sqrt(colSums(resid^2)/dof)
	Pvec <- 2*pt(-abs(t(est/(sqrt(diag(XXI))))/sigma), dof)
	return(Pvec)
}

matrix.paste0 <- function (m.list) {
	p <- length(m.list)
	
	m <- max(unlist(sapply(m.list, function (x) nrow(x))))
	n <- max(unlist(sapply(m.list, function (x) ncol(x))))
	res <- NULL
	for (i in 1:p) {
		mat <- m.list[[i]]
		if (!is.matrix(mat)) {
			mat <- matrix(mat, m, n)
		} else {
			if (!is.null(rownames(mat))) {
				row.names <- rownames(mat)
			}
			if (!is.null(colnames(mat))) {
				col.names <- colnames(mat)
			}
			
		}
		res <- paste0(res, mat)
	}
	return(matrix(res, m, n, dimnames=list(row.names, col.names)))
}
#
#normalize.func <- function (x) {
#	names <- paste(x)
#	tab <- table(x)
#	tab.name <- names(tab)
#	tab.x <- qqnorm(as.numeric(tab.name), plot.it = F)$x
#	names(tab.x) <- tab.name
#	as.numeric(tab.x[names])
#}

# Typer I error 1
calculateK <- function (G) {
	n <- nrow(G)
	m <- mean(diag(G))
	v <- var(G[lower.tri(G)])
	sigma <- v / m
	k <- m / sigma
	se <- 1.96 * sqrt(2 / n)
	list(K=k, K.lower = (sqrt(k) - se)^2, K.upper=(sqrt(k) + se)^2, sigma=sigma)
}

fPERMANOVA <- function (formula, data=NULL, contr.unordered = "contr.sum", contr.ordered = "contr.poly") {
	
	Terms <- terms(formula, data = data)
	lhs <- formula[[2]]
	lhs <- eval(lhs, data, parent.frame())
	
	formula[[2]] <- NULL
	rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
	op.c <- options()$contrasts
	options(contrasts = c(contr.unordered, contr.ordered))
	rhs <- model.matrix(formula, rhs.frame)
	options(contrasts = op.c)
	grps <- attr(rhs, "assign")
	
	# The first group includes the intercept, nterms count the intercept
	u.grps <- unique(grps)
	nterms <- length(u.grps)
	
	Z <- rhs[, grps %in% u.grps[1:(nterms-1)], drop=F]
	XZ <- rhs[, grps %in% u.grps[1:(nterms)], drop=F]
	
	n <- nrow(XZ)
	HZ <- Z %*% solve(t(Z) %*% (Z)) %*% t(Z)
	HXZ <- XZ %*% solve(t(XZ) %*% (XZ)) %*% t(XZ)
	HX <- HXZ - HZ
	HIX <- diag(n) - HXZ
	
	if (class(lhs) == 'dist') {
		D <- -as.matrix(lhs)^2 / 2
		
	} else {
		D <- -lhs^2 / 2
	}
	
	G <- mean(D) + D - rowMeans(D) - matrix(rep(1, n), ncol=1) %*% colMeans(D)
	G <- G / n
	obj <- calculateK(G)
	
	K <- obj$K
	
	sigma <- obj$sigma
	
	df1 <- ncol(XZ) - ncol(Z)
	df2 <- n - ncol(XZ)
	
	MSS <- sum(G * HX)
	RSS <- sum(G * HIX)
	TSS <- sum(diag(G))
#	RSS <- TSS - MSS
	
	f.stat <- (MSS / df1) / (RSS / df2)
	
#	p.value1 <- pnorm(f.stat, mean=1, sd=sqrt(2 / (K * df1)), lower.tail=F)
	p.value2 <- pchisq(f.stat * K * df1, df = K * df1, lower.tail=F)
	
	SumsOfSqs <- c(MSS, RSS, TSS)
	
	tab <- data.frame(Df = c(df1, df2, n - 1), SumsOfSqs = n * SumsOfSqs, 
			MeanSqs = c(n * MSS/df1, n * RSS/df2, NA), F.Model = c(f.stat, 
					NA, NA), R2 = c(MSS/TSS, NA, NA), 
			"Pr(>F)" = c(p.value2, NA, NA))
	
	rownames(tab) <- c("Model(Adjusted)", "Residuals", "Total")
	colnames(tab)[ncol(tab)] <- c("Pr(>F)")
	class(tab) <- c("anova", class(tab))
	attr(tab, "heading") <- c("F stat and P value of the last term is adjusted by previous terms!\n")
	
	out <- list(aov.tab = tab, K=K, sigma=sigma, call = match.call())
	out
}

# Determine the optimal depth to filter the samples
optim_dep <- function (data.obj,  ID.ctrl,  ID.other = NULL, taxa.level = 'Species', ann = '') {
	ID <- rownames(data.obj$meta.dat)
	if (is.null(ID.other)) {
		ID.other <- setdiff(ID, ID.ctrl)
	}
	
	dep <- colSums(data.obj$otu.tab)
	gen <- data.obj$abund.list[[taxa.level]]
	gen <- t(gen) / colSums(gen)
	
	dist.mat <- as.matrix(vegdist(gen))
	
	dist.mat <- dist.mat[ID.other, ID.ctrl]
	dep <- dep[ID.other]
	ord <- order(dep)
	dep <- dep[ord]
	dist.mat <- dist.mat[ord, ]
#	data <- data.frame(distance=as.vector(t(dist.mat)), depth=rep(log10(dep), each = ncol(dist.mat)))
	data <- data.frame(distance=rowMeans(dist.mat), depth=log10(dep))
	
	pdf(paste0('Distance2CtrlvsDepth.', taxa.level, '.', ann, '.pdf'), width = 6, height = 5)
	obj <- ggplot(data, aes(x = depth, y = distance)) +
			#	geom_boxplot() +
			geom_smooth() +
			geom_jitter() +
			geom_vline(xintercept=log10(2000), col='red') +
			xlab('Log10(depth)') +
			theme_bw()
	print(obj)
	dev.off()
}


# Determine the influential taxa - Need to test
contaminant_detection <- function (data.obj, ID.ctrl, ID.other = NULL, taxa.level = 'Species', ann = '') {
	ID <- rownames(data.obj$meta.dat)
	if (is.null(ID.other)) {
		ID.other <- setdiff(ID, ID.ctrl)
	}

	dep <- colSums(data.obj$otu.tab)
	gen <- data.obj$abund.list[[taxa.level]]
	gen <- t(gen) / colSums(gen)
	
	func <- function (gen, ID.other, ID.ctrl) {
		n.ctrl <- length(ID.ctrl)
		n.other <- length(ID.other)
		dist.mat <- as.matrix(vegdist(gen))
		dist1 <- rowMeans(dist.mat[ID.other, ID.ctrl])
		#dist2 <-  rowMeans(dist.mat[ID.other, ID.other]) * n.other / (n.other - 1)
		
		# mean(dist1) / mean(dist2)
		mean(dist1)
	}
	
	r0 <- func(gen, ID.other, ID.ctrl)
	r.vec <- numeric(ncol(gen))
	names(r.vec) <- colnames(gen)
	for (i in 1 : ncol(gen)) {
		if (i %% 10 == 0) cat('.')
#		gen2 <- gen[, -i]
#		gen2 <- gen2 / rowSums(gen2)
		gen2 <- gen
		gen2[, i] <- sample(gen2[, i])
		r.vec[i] <- func(gen2, ID.other, ID.ctrl)
	}
	return(list(r0=r0, r.vec=r.vec))
}

# New: 2018_07_13
generate_count_heatmap <- function (data.obj, ann='CountHeatmaps', width = 6, height = 20, 
		taxa.level = 'OTU', cutoff = 10, occurence = 2,
		breaks = c(-0.1, 0.1, 1, 2, 3, 4, 5, 6, 7,  8, 9, 10, 12, 14, 18),
		col = c('white', colorpanel(13, 'aliceblue', 'lightblue', 'blue')), 
		margins = c(10, 5),
		cexRow = 0.2, 
		cexCol = 0.4,
		...) {
	if (taxa.level %in% c('Species', 'OTU')) {
		otu.tab <- data.obj$otu.tab 
	} else {
		otu.tab <- data.obj$abund.list[[taxa.level]]
	}
	
	otu.tab <- otu.tab[rowSums(otu.tab) > cutoff, ]
	otu.tab <- otu.tab[rowSums(otu.tab != 0) > occurence, ]
	
	pdf(paste0(ann, '.', taxa.level, '.pdf'), width = width, height=height)
	heatmap.2(log2(otu.tab+1), scale = 'none', col = col,
			trace = 'none', breaks = breaks, key.xlab= 'log2(count + 1)', 
			density.info='none', margins = margins, cexRow = cexRow, cexCol =cexCol, ...)
	dev.off()
}


ZicoSeq.plot2 <- function (ZicoSeq.obj, meta.dat, grp.name = grp.name, 
                           pvalue.type = c("p.adj.fdr", "p.raw", "p.adj.fwer"), 
                           cutoff = 0.1, text.size = 10, out.dir = NULL, 
                           file.name = "ZicoSeq.plot.pdf", width = 10, height = 6) 
{
  # grp.name <- ZicoSeq.obj$call$grp.name # it does not give the correct grp.name, thus disable this code and add "grp.name" para
  if (ZicoSeq.obj$call$feature.dat.type == "proportion") {
    abundance <- ZicoSeq.obj$feature.dat
    prevalence <- apply(ZicoSeq.obj$feature.dat, 1, function(x) mean(x > 
                                                                       0))
  }
  if (ZicoSeq.obj$call$feature.dat.type == "count") {
    abundance <- t(t(ZicoSeq.obj$feature.dat)/colSums(ZicoSeq.obj$feature.dat))
    prevalence <- apply(ZicoSeq.obj$feature.dat, 1, function(x) mean(x > 
                                                                       0))
  }
  if (ZicoSeq.obj$call$feature.dat.type == "other") {
    abundance <- ZicoSeq.obj$feature.dat
    prevalence <- apply(abundance, 1, function(x) sd(x))
  }
  
  abundance <- apply(abundance, 1, function(x) mean(x))
  if (!(length(unique(meta.dat[, grp.name])) > 2 & !is.numeric(meta.dat[, grp.name]))) {
    coefs <- sapply(ZicoSeq.obj$coef.list, function(x) x[grep(grp.name, rownames(x)), ])
    colnames(coefs) <- paste0("coef_Func", 1:length(ZicoSeq.obj$coef.list))
  }
  if (!is.numeric(meta.dat[, grp.name]) & length(unique(meta.dat[, 
                                                                 grp.name])) == 2) {
    level2 <- rownames(ZicoSeq.obj$coef.list[[1]])[grep(paste0("^", 
                                                               grp.name), rownames(ZicoSeq.obj$coef.list[[1]]))]
    level2 <- gsub(grp.name, "", level2)
    base <- setdiff(unique(meta.dat[, grp.name]), level2)
  }
  ZicoSeq.obj$R2 <- ZicoSeq.obj$R2 + runif(length(ZicoSeq.obj$R2), 
                                           0, 1e-10)
  R2 <- matrixStats::rowMaxs(ZicoSeq.obj$R2)
  if (length(unique(meta.dat[, grp.name])) > 2 & !is.numeric(meta.dat[, 
                                                                      grp.name])) {
    plot.data <- data.frame(pvals = ZicoSeq.obj[[pvalue.type]], 
                            prevalence = prevalence, abundance = abundance, R2 = R2, 
                            taxa = rownames(ZicoSeq.obj$R2))
    plot.data[plot.data$pvals > cutoff, "taxa"] <- ""
  }else {
    signs <- t(sign(coefs))[t(ZicoSeq.obj$R2 == R2)]
    signs[signs == 0] <- 1
    plot.data <- data.frame(pvals = ZicoSeq.obj[[pvalue.type]], 
                            prevalence = prevalence, abundance = abundance, R2 = R2 * 
                              signs, taxa = rownames(ZicoSeq.obj$R2))
    plot.data[plot.data$pvals > cutoff, "taxa"] <- ""
  }
  if (is.numeric(meta.dat[, grp.name])) {
    title.name = paste0("Association between ", grp.name, 
                        " and  taxa abundance")
  }
  if (!is.numeric(meta.dat[, grp.name]) & length(unique(meta.dat[, 
                                                                 grp.name])) > 2) {
    title.name = paste0("Association between ", grp.name, 
                        " and  taxa abundance")
  }
  if (!is.numeric(meta.dat[, grp.name]) & length(unique(meta.dat[, 
                                                                 grp.name])) == 2) {
    title.name = paste0("Differential abundance between ", 
                        paste0(base, " (reference) and ", level2))
  }
  pvals <- R2 <- taxa <- NULL
  plot.obj <- ggplot(plot.data, aes(x = R2, y = -log10(pvals))) + 
    geom_point(aes(size = abundance, color = prevalence)) + 
    geom_vline(aes(xintercept = 0), color = "gray", linetype = "dashed") + 
    geom_hline(aes(yintercept = -log10(cutoff)), color = "gray", 
               linetype = "dashed") + scale_colour_gradient2(low = "white", 
                                                             high = "#006D2C") + scale_y_continuous(limits = c(0, 
                                                                                                               max(-log10(plot.data$pvals)) * 1.3)) + ggrepel::geom_text_repel(aes(label = taxa), 
                                                                                                                                                                               max.overlaps = Inf, color = "black") + labs(x = bquote(R^2), 
                                                                                                                                                                                                                           y = paste0("-log10(", pvalue.type, ")"), color = ifelse(ZicoSeq.obj$call$feature.dat.type == 
                                                                                                                                                                                                                                                                                     "other", "Standard deviation", "Prevalence"), size = ifelse(ZicoSeq.obj$call$feature.dat.type == 
                                                                                                                                                                                                                                                                                                                                                   "other", "Mean value", "Mean abundance")) + theme_bw() + 
    theme(axis.text = element_text(color = "black", size = text.size), 
          axis.title = element_text(color = "black", size = text.size), 
          legend.text = element_text(color = "black", size = text.size), 
          legend.title = element_text(color = "black", size = text.size)) + 
    ggtitle(title.name)
  if (!is.null(out.dir)) {
    print(plot.obj)
    ggsave(paste0(out.dir, file.name), width = width, height = height)
  }
  return(plot.obj)
}





# y is a distance matrix
partial.r.adj.dmanova <- function (x, z = NULL, y) {
  x <- as.data.frame(x)
  if (!is.null(z)) {
    z <- as.data.frame(z)
  }
  
  n <- nrow(x)
  x <- model.matrix( ~ . - 1, x)
  if (!is.null(z)) {
    z <- model.matrix( ~ . - 1, z)
  }
  df.x <- ncol(x)
  if (!is.null(z)) {
    df.z <- ncol(z)
  } else {
    df.z <- 0
  }
  
  if (!is.null(z)) {
    lm.obj1 <- dmanova(y ~ z)
    lm.obj2 <- dmanova(y ~ x + z)
    tss <- lm.obj1$aov.tab[2, 2]
    rss <- lm.obj2$aov.tab[2, 2]
  } else {
    lm.obj2 <- dmanova(y ~ x)
    tss <- lm.obj2$aov.tab[3, 2]
    rss <- lm.obj2$aov.tab[2, 2]
  }
  
  1 - (n - 1 - df.z) * rss / (n - 1 - df.x - df.z) / tss
}


# linear regression
# x variable of interest, could be vector or data.frame
# z variables to be adjusted, could be vector or data.frame
# y variable, vector
partial.r.adj <- function (x, z = NULL, y) {
  x <- as.data.frame(x)
  if (!is.null(z)) {
    z <- as.data.frame(z)
  }
  
  n <- length(y)
  x <- model.matrix( ~ . - 1, x)
  if (!is.null(z)) {
    z <- model.matrix( ~ . - 1, z)
  }
  df.x <- ncol(x)
  if (!is.null(z)) {
    df.z <- ncol(z)
  } else {
    df.z <- 0
  }
  
  if (!is.null(z)) {
    lm.obj1 <- summary(lm(y ~ z))
    lm.obj2 <- summary(lm(y ~ x + z))
  } else {
    lm.obj1 <- summary(lm(y ~ 1))
    lm.obj2 <- summary(lm(y ~ x))
  }
  tss <- sum(lm.obj1$residuals ^ 2)
  rss <- sum(lm.obj2$residuals ^ 2)
  1 - (n - 1 - df.z) * rss / (n - 1 - df.x - df.z) / tss
  
}



winzor <- function(x, outlier.pct=0.01, winsor.end = c("top","bottom","both")){
  if (winsor.end == "top") {
    qt <- quantile(x, 1 - outlier.pct)
    x[x >= qt] <- qt
  }
  if (winsor.end == "bottom") {
    qt <- quantile(x, outlier.pct)
    x[x <= qt] <- qt
  }
  if (winsor.end == "both") {
    qt1 <- quantile(x, 1 - outlier.pct/2)
    qt2 <- quantile(x, outlier.pct/2)
    x[x >= qt1] <- qt1
    x[x <= qt2] <- qt2
  }
  return(x)
}




perform_RF_prediction <- function(data, formula, grp.name, imbalanced=F,
                                  pkg = c('randomForest','ranger'), 
                                  positive=NULL, # the group you want to you predict
                                  Boruta=F, train = F){
  data[,grp.name] <- as.factor(data[,grp.name])
  if(is.null(positive)){positive=levels(data[,grp.name])[1]}
  if(pkg=='randomForest'){
    cat('-Using randomForest package-\n')
    if(train ==F){
      mtry <- round(sqrt(ncol(data)))
      rf <- randomForest(as.formula(formula),data = data, 
                         importance = TRUE, na.action = na.omit, 
                         ntree = 500, mtry = mtry, seed = 12)
      pred <- rf$votes
      OOB <- mean(rf$err.rate[,'OOB'])
    }else{
      idx <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
      train.idx <- idx==1
      test.idx <- idx==2
      mtry <- round(sqrt(ncol(data)))
      rf <- randomForest(as.formula(formula),data = data, 
                         importance = TRUE, na.action = na.omit, 
                         ntree = 500, mtry = mtry, seed = 12)
      pred <- rf$votes
      OOB <- mean(rf$err.rate[,'OOB'])
    }
  }
  
  if(pkg=='ranger'){
    cat('-Using ranger package-\n')
    if(imbalanced==T){
      ## add 02/07/2024 for imbalanced data issue
      class_wt <- 1 / prop.table(table(data[,grp.name]))
      frc <- as.vector(table(data[,grp.name])/nrow(data))
      rf <- ranger(as.formula(formula), data = data, probability = TRUE, importance = 'impurity', 
                   seed = 12, class.weights = class_wt, 
                   sample.fraction = frc, # stratified sampling
                   replace = T, keep.inbag = TRUE)
      
    }else{
      rf <- ranger(as.formula(formula), data = data, probability = TRUE, importance = 'impurity', seed = 12)
    }
    pred <- rf$predictions
    rownames(pred) <- rownames(data)
    OOB <- rf$prediction.error
  }
  
  cat('--Prediction finished--\n')
  pred.res <- merge(pred, data[,grp.name, drop = F], all = T, by = 0)
  row.names(pred.res) <- pred.res$Row.names
  pred.res$Row.names <- NULL
  
  roc.obj <- pROC::roc(pred.res[,grp.name], rf$predictions[,positive], ci=TRUE, plot=F, direction='>')
  ci.obj <- pROC::ci.se(roc.obj)
  dat.ci <- data.frame(x = as.numeric(rownames(ci.obj)),lower = ci.obj[, 1],upper = ci.obj[, 3])
  
  pred.v <- as.factor(ifelse(rf$predictions[,levels(data[,grp.name])[1]] > 0.5,levels(data[,grp.name])[1],levels(data[,grp.name])[2]))
  levels(pred.v) <- levels(data[,grp.name])
  acc <- calculate_performance(data[,grp.name], pred.v, positive = positive)
  
  plt.clinc <- pROC::ggroc(roc.obj, color = brewer.pal(8,'Dark2')[2]) + 
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") + 
    theme_bw() + 
    theme(axis.title = element_text(color = 'black', size = 16),
          panel.grid = element_blank(),
          axis.text = element_text(color = 'black', size = 16)) + 
    geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), fill = "orange", alpha= 0.2)+
    annotate(geom="text", x=0.25, y=0.1, label=gsub('\\(DeLong\\)','',capture.output(roc.obj$ci)),color="black") +
    ggtitle(paste0(grp.name,' [AUC: ',round(auc(roc.obj),2),', OOB:',round(OOB,2),']')) 
  
  if(Boruta==T){
    boruta_output <- Boruta(as.formula(formula), data = data, doTrace = 2)
    boruta.select <- boruta_output$finalDecision[boruta_output$finalDecision %in% c('Confirmed','Tentative')] %>% droplevels()
    return(list(roc.plot = plt.clinc, boruta_output=boruta_output, rf=rf, boruta.select=boruta.select, roc.obj = roc.obj, ci.obj=ci.obj,OOB=OOB, acc=acc, pred.res=pred.res))
  }else{
    return(list(plt.clinc=plt.clinc, rf=rf, roc.obj=roc.obj, ci.obj=ci.obj,OOB=OOB, acc=acc,pred.res=pred.res))
  }
}


generate_ROC <- function(rf.obj, ann=''){
  ## generate ROC by the length of the rf.obj model
  coord.all <- bind_rows(
    lapply(rf.obj, function(x){
      coord <- pROC::coords(x$roc.obj, "all", transpose = FALSE)
      return(coord[rev(seq(nrow(coord))),])
    }), .id = "source")
  
  dat.ci.all <- bind_rows(
    lapply(rf.obj, function(x){
      dat.ci <- data.frame(x = as.numeric(rownames(x$ci.obj)),lower = x$ci.obj[, 2],upper = x$ci.obj[, 3])
      return(dat.ci)
    }), .id = "source")
  
  CI <- sapply(rf.obj, function(x) gsub('\\ \\(DeLong\\)','',capture.output(x$roc.obj$ci)))
  CI_vector <- paste(names(CI), CI, sep = ": ")
  AUC <- round(sapply(rf.obj, function(x) x$roc.obj$auc),2)
  OOB <- round(sapply(rf.obj, function(x) x$OOB),2)
  AUC_OOB <- paste0('OOB/AUC [',paste(names(AUC),paste(OOB, AUC, sep = "/"), sep = ":", collapse = '; '),']')
  col <- c(brewer.pal(8,'Dark2'), brewer.pal(12,'Paired'),brewer.pal(8,'Accent')[1:4])[1:length(rf.obj)]
  names(col) <- names(rf.obj)
  plt <- ggplot() +
    geom_line(data = coord.all , aes(x = specificity, y = sensitivity, color = source), lwd = 0.5) +
    ggplot2::scale_x_reverse(lim=c(1, 0)) +
    geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0), color="grey", linetype="dashed") + 
    geom_ribbon(data = dat.ci.all, aes(x = x, ymin = lower, ymax = upper, fill = source), alpha= 0.1) +
    scale_fill_manual(values = col)+
    scale_color_manual(values = col)+
    theme_bw() + 
    theme(axis.title = element_text(color = 'black', size = 16),
          panel.grid = element_blank(),
          plot.title = element_text(size = 10),
          legend.text =  element_text(color = 'black', size = 14),
          axis.text = element_text(color = 'black', size = 16)) +
    annotate(geom="text", x=rep(0.3, length(rf.obj)), y=seq(0.2, 1, by = 0.05)[1:length(rf.obj)], label=stringr::str_wrap(CI_vector, width = 50),color="black") +
    ggtitle(paste(ann,AUC_OOB))+
    labs(color = '', fill ='')
  return(plt)
}


otu_filter <- function(feature.dat, prev.filter = 0.1, max.abund.filter = 0.002, mean.abund.filter = 0, verbose=T){
  temp <- t(t(feature.dat)/colSums(feature.dat))
  filter.ind <- rowMeans(temp != 0) >= prev.filter & rowMeans(temp) >= 
    mean.abund.filter & matrixStats::rowMaxs(temp) >= max.abund.filter
  names(filter.ind) <- rownames(feature.dat)
  if (verbose) 
    cat(sum(!filter.ind), " features are filtered!\n ",sum(filter.ind), " features are remained!\n ")
  filter.features <- rownames(feature.dat)[!filter.ind]
  feature.dat <- feature.dat[filter.ind, ]
  return(feature.dat)
}



test_match_order <- function(x,y) {
  
  if (all(x==y)) print('Perfect match in same order')
  
  if (!all(x==y) && all(sort(x)==sort(y))) print('Perfect match in wrong order')
  
  if (!all(x==y) && !all(sort(x)==sort(y))) print('No match')
}



calculate_performance <- function(predicted.labels, actual.labels, positive=NULL){
  if(is.null(positive)){positive = levels(actual.labels)[1]}
  actual.labels <- as.factor(actual.labels)
  predicted.labels <- factor(predicted.labels, levels = levels(actual.labels))
  tb <- (table(predicted.labels, actual.labels))
  negative <- colnames(tb)[!(colnames(tb) %in% positive)]
  FP <- tb[positive,negative]
  FN <- tb[negative,positive]
  TP <- tb[positive,positive]
  TN <- tb[negative,negative]
  recall <- sensitivity <- ifelse(TP==0, 0, TP/(TP+FN))
  precision <- ifelse(TP==0, 0,TP/(TP+FP))
  specificity <- ifelse(TN==0, 0,TN / (TN + FP))
  F1 <- ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))
  balanced_accuracy <- (sensitivity + specificity) / 2
  res <- c(F1, precision, recall, sensitivity, specificity,balanced_accuracy)
  names(res) <- c('F1', 'precision', 'recall', 'sensitivity', 'specificity','balanced_accuracy')
  return(res)
}





calculate_AUC_OOB <- function(rf.obj, data, grp.name){
  conf <- rf.obj[['confusion']][,-ncol(rf.obj[['confusion']])]
  oob <- 1 - (sum(diag(conf))/sum(conf))
  pred <- prediction(rf.obj[['votes']][,2], data[,variable], label.ordering = NULL)
  perf <- performance(pred, "tpr", "fpr")
  auc <- mean(unlist(performance(pred, 'auc')@y.values))
  return(c(AUC =auc, OOB=oob))
}



perform_RandomForest <- function(data.obj, grp.name, adj.name, taxa.level, 
                                 positive = NULL,
                                 prev=0.1, minp=0.002,meanp=0, 
                                 imbalanced=F, weighted.mt = c("VotingRule","SamplingStratification","ClassWeight"),
                                 Boruta=F,
                                 seed = 123){
  weighted.mt = match.arg(weighted.mt)
  
  
  set.seed(seed)
  tt <- otu_filter(data.obj$abund.list[[taxa.level]],prev.filter = prev, max.abund.filter = minp, mean.abund.filter = meanp)
  longname <- rownames(tt)
  names(longname) <- paste0('s',1:length(longname))
  rownames(tt) <- NULL
  rownames(tt) <- names(longname)
  data <- as.data.frame(merge(t(tt), data.obj$meta.dat[,c(adj.name, variable), drop =F], by = 0)) %>% column_to_rownames('Row.names') %>% droplevels()
  adj.name2 <- colnames(data)[!(colnames(data) %in% variable)]
  formula <- paste0(variable,' ~ ',paste0(adj.name2, collapse = '+'))
  
  if(imbalanced==F){
    rf.obj <- randomForest(as.formula(formula),data=data)
  }else{
    prop.tb <- prop.table(table(data[,grp.name]))
    rare.class.prevalence = unname(prop.tb[which.min(prop.tb)])
    
    if(weighted.mt =="VotingRule"){
      #1 - Balancing by voting rule, AUC of ROC will be unchanged...
      rf.obj = randomForest(as.formula(formula),data=data,cutoff=c(1-rare.class.prevalence,rare.class.prevalence))
    }
    
    if(weighted.mt =="SamplingStratification"){
      #2 - Balancing by sampling stratification
      nRareSamples = round(nrow(data) * rare.class.prevalence)
      rf.obj = randomForest(as.formula(formula),data=data,strata=data[,variable],sampsize=c(nRareSamples,nRareSamples))
    }
    
    if(weighted.mt =="ClassWeight"){
      #3 - Balancing by class-weight.
      w <- 1 / table(data[,grp.name]);w <- w/sum(w)
      rf.obj = randomForest(as.formula(formula),data=data,classwt=w)
    }
  }
  
  if(is.null(positive)){positive = levels(data[,grp.name])[1]}
  acc <- calculate_AUC_OOB(rf.obj, data = data, grp.name = grp.name)
  perf <- c(calculate_performance(actual.labels=rf.obj$y,predicted.labels=rf.obj$predicted, positive = positive), acc)
  
  roc.obj <- pROC::roc(rf.obj$y, rf.obj$votes[,positive], ci=TRUE, plot=T, direction='>')
  ci.obj <- pROC::ci.se(roc.obj)
  dat.ci <- data.frame(x = as.numeric(rownames(ci.obj)),lower = ci.obj[, 1],upper = ci.obj[, 3])
  OOB <- calculate_AUC_OOB(rf.obj, data = data, grp.name = grp.name)
  plt <- pROC::ggroc(roc.obj, color = brewer.pal(8,'Dark2')[2]) + 
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") + 
    theme_bw() + 
    theme(axis.title = element_text(color = 'black', size = 16),
          panel.grid = element_blank(),
          axis.text = element_text(color = 'black', size = 16)) + 
    geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), fill = "orange", alpha= 0.2)+
    annotate(geom="text", x=0.25, y=0.1, label=gsub('\\(DeLong\\)','',capture.output(roc.obj$ci)),color="black") +
    ggtitle(paste0(grp.name,' [AUC: ',round(OOB['AUC'],2),', OOB:',round(OOB['OOB'],2),']')) 
  
  
  if(Boruta ==T){
    boruta_output <- Boruta(as.formula(formula), data = data, doTrace = 2)
    boruta.select <- names(boruta_output$finalDecision[boruta_output$finalDecision %in% c('Confirmed','Tentative')] %>% droplevels())
    boruta.select.s <- boruta.select[!(boruta.select %in% adj.name)]
    boruta.select.sn <- (longname[names(longname) %in% boruta.select.s])[(boruta.select.s)]
    boruta.select <- c(boruta.select[(boruta.select %in% adj.name)],boruta.select.sn)
    return(list(perf=perf, rf.obj= rf.obj, plt = plt, boruta.select=boruta.select))
  }else{
    return(list(perf=perf, rf.obj= rf.obj, plt = plt))
  }
}

perform_xgboost <- function(data.obj, train.ratio = 0.8, grp.name, adj.name, positive=NULL, 
                            prev=0.1, minp=0.002, meanp=0, seed=123, nthread=2){
  library(xgboost)
  set.seed(seed)
  tt <- otu_filter(data.obj$abund.list[[taxa.level]],prev.filter = prev, max.abund.filter = minp, mean.abund.filter = meanp)
  longname <- rownames(tt)
  names(longname) <- paste0('s',1:length(longname))
  rownames(tt) <- NULL
  rownames(tt) <- names(longname)
  data <- as.data.frame(merge(t(tt), data.obj$meta.dat[,c(adj.name, grp.name), drop =F], by = 0)) %>% column_to_rownames('Row.names') %>% droplevels()
  
  data_case <- data[data[,grp.name]==unique(data[,grp.name])[1],]
  data_control <- data[data[,grp.name]!=unique(data[,grp.name])[1],]
  
  nRareSamples = round(0.8 * min(nrow(data_case), nrow(data_control)))
  train_case <- data_case[sample(nrow(data_case), size = nRareSamples), ]
  train_control <- data_control[sample(nrow(data_control), size = nRareSamples), ]
  
  data.train <- rbind(train_case, train_control)
  data.test <- data[!rownames(data) %in% rownames(data.train), ]
  
  train.matrix <- as.matrix(data.train[, -which(names(data.train) == grp.name)])
  test.matrix <- as.matrix(data.test[, -which(names(data.test) == grp.name)])
  
  if(is.null(positive)){positive <- levels(data.train[,grp.name])[1]}
  
  label <- ifelse(data.train[,grp.name]==positive,1,0)
  bst1 <- xgboost(data = train.matrix, label = label, max_depth = 2, eta = 1, nthread = nthread, nrounds = 2,objective = "binary:logistic")
  pred <- predict(bst1, test.matrix)
  pred_binary <- ifelse(pred > 0.5, 1, 0)
  true_labels <- ifelse(data.test[,grp.name]==positive,1,0)
  roc_curve <- pROC::roc(true_labels, pred, direction = '>')
  perf <- c(calculate_performance(actual.labels=true_labels,predicted.labels=pred_binary, positive = '1'), AUC=pROC::auc(roc_curve))
  return(perf)
}


# library(caret)
# library(xgboost)
# library(mlr)
# set.seed(seed)
# tt <- otu_filter(data.obj$abund.list[[taxa.level]],prev.filter = prev, max.abund.filter = minp, mean.abund.filter = meanp)
# longname <- rownames(tt)
# names(longname) <- paste0('s',1:length(longname))
# rownames(tt) <- NULL
# rownames(tt) <- names(longname)
# data <- as.data.frame(merge(t(tt), data.obj$meta.dat[,c(adj.name, grp.name), drop =F], by = 0)) %>% column_to_rownames('Row.names') %>% droplevels()
# 
# data_case <- data[data[,grp.name]==unique(data[,grp.name])[1],]
# data_control <- data[data[,grp.name]!=unique(data[,grp.name])[1],]
# 
# nRareSamples = round(0.8 * min(nrow(data_case), nrow(data_control)))
# train_case <- data_case[sample(nrow(data_case), size = nRareSamples), ]
# train_control <- data_control[sample(nrow(data_control), size = nRareSamples), ]
# 
# data.train <- rbind(train_case, train_control)
# data.test <- data[!rownames(data) %in% rownames(data.train), ]
# 
# #using one hot encoding
# labels <- data.train[,grp.name]
# ts_label <- data.test[,grp.name]
# new_tr <- model.matrix(~.+0,data = data.train %>% dplyr::select(-grp.name))
# new_ts <- model.matrix(~.+0,data = data.test %>% dplyr::select(-grp.name))
# 
# #convert factor to numeric
# labels <- as.numeric(labels)-1
# ts_label <- as.numeric(ts_label)-1
# 
# #preparing matrix
# dtrain <- xgb.DMatrix(data = new_tr,label = labels)
# dtest <- xgb.DMatrix(data = new_ts,label=ts_label)
# 
# 
# params <- list(booster = "gbtree", objective = "binary:logistic", eta=0.3, gamma=0, max_depth=6, min_child_weight=1, subsample=1, colsample_bytree=1)
# xgbcv <- xgb.cv( params = params, data = dtrain, nrounds = 100, nfold = 5, showsd = T, stratified = T, print.every.n = 10, early.stop.round = 20, maximize = F)
# xgb1 <- xgb.train (params = params, data = dtrain, nrounds = xgbcv$best_iteration, watchlist = list(val=dtest,train=dtrain), print.every.n = 10, early.stop.round = 10, maximize = F , eval_metric = "error")
# xgbpred <- predict(xgb1,dtest)
# xgbpred <- ifelse (xgbpred > 0.5,1,0)
# 
# acc <- confusionMatrix(as.factor(xgbpred), as.factor(ts_label), positive = "1")
# Precision <- unname(acc$byClass['Precision'])
# Recall <- unname(acc$byClass['Recall'])
# 2*(Precision*Recall/(Precision+Recall))

dmanova2 <- function (formula, data = NULL, positify = FALSE, contr.unordered = "contr.sum", 
                      contr.ordered = "contr.poly", returnG = FALSE) 
{
  Terms <- terms(formula, data = data)
  lhs <- formula[[2]]
  lhs <- eval(lhs, data, parent.frame())
  formula[[2]] <- NULL
  rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
  op.c <- options()$contrasts
  options(contrasts = c(contr.unordered, contr.ordered))
  rhs <- model.matrix(formula, rhs.frame)
  options(contrasts = op.c)
  grps <- attr(rhs, "assign")
  u.grps <- unique(grps)
  nterms <- length(u.grps)
  Z <- rhs[, grps %in% u.grps[1:(nterms - 1)], drop = F]
  XZ <- rhs[, grps %in% u.grps[1:(nterms)], drop = F]
  n <- nrow(XZ)
  if (inherits(lhs, "dist")) {
    D <- as.matrix(lhs)
  }
  else {
    D <- lhs
  }
  D <- -D^2/2
  G <- mean(D) + D - rowMeans(D) - matrix(rep(1, n), ncol = 1) %*% 
    colMeans(D)
  if (positify) {
    G <- as.matrix(nearPD(G)$mat)
  }
  XZi <- solve(t(XZ) %*% XZ)
  Zi <- solve(t(Z) %*% (Z))
  HZ <- Z %*% Zi %*% t(Z)
  HXZ <- XZ %*% XZi %*% t(XZ)
  HX <- HXZ - HZ
  HIXZ <- diag(n) - HXZ
  HIX <- diag(n) - HX
  df1 <- ncol(XZ) - ncol(Z)
  df2 <- n - ncol(XZ)
  MSS <- sum(G * HX) ##MSS <- sum(diag(G * HX))
  MSS_all <- sum(G * HXZ) ## MSS_all <- sum(diag(G * HXZ))
  RSS <- sum(G * HIXZ)
  TSS <- sum(diag(G))
  f.stat <- (MSS/df1)/(RSS/df2)
  GXZ <- G %*% XZ
  XZXZi <- XZ %*% XZi
  GXZtXZXZi <- GXZ %*% t(XZXZi)
  G.tilde <- G + XZXZi %*% (t(GXZ) %*% XZ) %*% t(XZXZi) - GXZtXZXZi - 
    t(GXZtXZXZi)
  obj <- calculateK(G.tilde, HIXZ, n - df2)
  K <- obj$K
  p.value <- pchisq(f.stat * K * df1, df = K * df1, lower.tail = F)
  SumsOfSqs <- c(MSS, RSS, TSS)
  tab <- data.frame(Df = c(df1, df2, n - 1), SumsOfSqs = SumsOfSqs, 
                    MeanSqs = c(MSS/df1, RSS/df2, NA), F.Model = c(f.stat, NA, NA), 
                    R2.full = c(MSS_all / TSS, NA, NA),
                    R2 = c(MSS / TSS,  NA, NA),
                    `Pr(>F)` = c(p.value, NA, NA))
  rownames(tab) <- c("Model(Adjusted)", "Residuals", "Total")
  colnames(tab)[ncol(tab)] <- c("Pr(>F)")
  class(tab) <- c("anova", class(tab))
  attr(tab, "heading") <- c("F stat and P value of the last term is adjusted by preceding terms!\n")
  if (returnG == TRUE) {
    out <- list(aov.tab = tab, df = K, G = G, call = match.call())
  }
  else {
    out <- list(aov.tab = tab, df = K, call = match.call())
  }
  out
}

calculateK <- function (G, H, df2) {
  n <- nrow(G)
  dG <- diag(G)
  
  mu1 <- sum(dG) / (n - df2)
  mu2 <- (sum(G^2) - sum(dG^2)) / ((n - df2)^2 + sum(H^4) - 2 * sum(diag(H^2)))
  
  k <- mu1^2 / mu2
  list(K = k, mu1 = mu1, mu2 = mu2)
}
