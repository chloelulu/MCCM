################################################################################
# File    :   fastcclasso.R
# Aim     :   estimating correlation inference for compositional data
# Referece:   Zhang, S., Fang, H., and Hu, T. (2024), “fastCCLasso: Fast 
#             and efficient algorithm for estimating correlation matrix 
#             from compositional data,” Submitted.
#  Rcode  :   https://github.com/ShenZhang-Statistics/fastCCLasso
#-------------------------------------------------------------------------------
# Author : Zhang Shen (Capital Normal University)
# Email  : zhangshen@cnu.edu.cn
# Date   : 2024-01-11
#-------------------------------------------------------------------------------
# Main function: fastCCLasso(xx, isCnt = FALSE, pseudo = 0.5, k_cv = 3, 
#            	 lam_min_ratio = 1E-4, k_max = 20, n_boot=100) 
#
#  Input:
#            xx ------ n x p data matrix (row/column is sample/variable)
#                      n samples & p compositional variables
#         isCnt ------ Is the compositional data matrix a count matrix? 
#                      Default: FALSE
#        pseudo ------ pseudo count if isCnt = TRUE
#                      Default: 0.5
#          k_cv ------ folds of cross validation
#                      Default: 3     
# lam_min_ratio ------ Smallest tuning parameter value provided as a fraction of maximum
#                      Default: 1e-4
#         k_max ------ maximum iterations for golden section method
#                      Default: 20
#        n_boot ------ Bootstrap times
#                      Default: 100
#  Output: 
#      A list structure contains:
#           rho ------ correlation estimation
#      cov_diag ------ variance estimation      
#   lambda_best ------ the optimal turning parameter
#       info_cv ------ information for cross validation
#        p_vals ------ p-values for elements of rho equal 0 or not
#        aa, bb ------ the diagonal elements of the weight matrics A and B in fastCCLasso, respectively.
#-------------------------------------------------------------------------------
fastCCLasso <- function(xx, isCnt = FALSE, pseudo = 0.5, k_cv = 3, 
	lam_min_ratio = 1E-4, k_max = 20, n_boot=100, aa=NULL, bb=NULL) {
	n = nrow(xx);
	p = ncol(xx);
	if(isCnt) {
		xx = xx + pseudo;
		xx = xx / rowSums(xx);
	};
	xx2 = log(xx) - rowMeans(log(xx));
	vxx2 = stats::var(xx2);
	
	if(is.null(aa)){
           aa = rep(1, p);
        }
        if(is.null(bb)){
           bb = 1 / diag(vxx2);
        }		
	
	#-Golden section method for the selection of lambda (log10 scale)
	xx = vxx2 * (aa * rep(bb, each = p) + bb * rep(aa, each = p))/2;
	diag(xx) = 0;
	lam_max = max(abs(xx));
	lam_int2 = log10(lam_max * c(lam_min_ratio, 1));
	a1 = lam_int2[1]; 
	b1 = lam_int2[2];
	
	#-Store lambda and corresponding cross validation's loss
	lams = NULL; 
	fvals = NULL;
	#-Two trial points in first 
	a2 = a1 + 0.382 * (b1 - a1); 
	b2 = a1 + 0.618 * (b1 - a1);
	fb2 = cvfastCCLasso(lambda = 10^b2, k_cv = k_cv, xx2 = xx2, 
		aa = aa, bb = bb);
	lams = c(lams, b2); 
	fvals = c(fvals, fb2);
	fa2 = cvfastCCLasso(lambda = 10^a2, k_cv = k_cv, xx2 = xx2, 
		aa = aa, bb = bb);
	lams = c(lams, a2); 
	fvals = c(fvals, fa2);
	#Error tolerance for convergence
	err_lam2 = 1e-1 * max(1, lam_int2);
	err_fval = 1e-4;
	err = b1 - a1;
	k = 0;
	
	while(err > err_lam2 && k < k_max) {
		fval_max = max(fa2, fb2);
		if(fa2 > fb2) {
			a1 = a2;
			a2 = b2;
			fa2 = fb2;
			b2 = a1 + 0.618 * (b1 - a1);
			fb2 = cvfastCCLasso(lambda = 10^b2, k_cv = k_cv, xx2 = xx2, 
				aa = aa, bb = bb);
			lams = c(lams, b2); 
			fvals = c(fvals, fb2);
		} else {
			b1 = b2;
			b2 = a2;
			fb2 = fa2;
			a2 = a1 + 0.382 * (b1 - a1);
			fa2 = cvfastCCLasso(lambda = 10^a2, k_cv = k_cv, xx2 = xx2, 
				aa = aa, bb = bb);
			lams = c(lams, a2);
			fvals = c(fvals, fa2);
		};
		fval_min = min(fa2, fb2);
		k = k + 1;
		err = b1 - a1;
		if(abs(fval_max - fval_min) / (1 + fval_min) <= err_fval) {
			break;
		};
	};
	info_cv = list(lams = lams, fvals = fvals, k = k + 2, 
		lam_int = 10^c(a1, b1));
	#if(a1 == lam_int2[1] || b1 == lam_int2[2]) {
	#	cat("WARNING:\n", "\tOptimal lambda is near boundary! ([", 10^a1, ",", 
	#		10^b1, "])\n", sep = "");
	#};
	lambda = 10^((a2 + b2)/2);
	fit_res = fastcclasso_sub(lambda = lambda, SS2 = vxx2, aa = aa, bb = bb);

	sigma_mod <- boot_fastCCLasso(xx2=xx2,sigma_hat= fit_res$sigma,
					 lambda=lambda,aa=aa,bb=bb,
					 n_boot = n_boot, max_iter=200,
					 stop_eps=1e-6)

	return(list(rho=sigma_mod$cor_w,
				cov_diag=sigma_mod$var_w,
				lambda_best=lambda,
				info_cv=info_cv,
				p_vals=sigma_mod$p_vals))
};
#---------------------------------------
#cross validation's loss of fastcclasso for single lambda
cvfastCCLasso <- function(lambda, k_cv, xx2, aa, bb) {
	n = nrow(xx2);
	p = ncol(xx2);
	n_b = floor(n / k_cv);
	cv.loss = 0;
	for(k in 1:k_cv) {
		ite = (n_b * (k - 1) + 1):(n_b * k);
		vxx2te = stats::var(xx2[ite, ]);
		vxx2tr = stats::var(xx2[-ite, ]);
		out = fastcclasso_sub(lambda = lambda, SS2 = vxx2tr, aa = aa, bb = bb);
		mm = out$sigma - out$ww - rep(out$ww, each = p) - vxx2te;
		cv.loss = cv.loss + mean(mm^2 * aa * rep(bb, each = p));
  };
  return(cv.loss);
};
#---------------------------------------
#fastcclasso for single lambda
fastcclasso_sub <- function(lambda, SS2, aa, bb, k_max = 200, x_tol = 1E-4) {
	p = ncol(SS2);
	cc = 1 / (aa * sum(bb) + bb * sum(aa));
	aa2 = aa * cc;
	bb2 = bb * cc;
	cab1 = 1 + sum(aa * bb2);
	caa = sum(aa * aa2);
	cbb = sum(bb * bb2);
	aabb = aa * rep(bb, each = p) + bb * rep(aa, each = p);
	lambda2 = 2 * lambda / aabb;
	ss2 = rowSums(SS2 * aabb);
	sigma = SS2;
	ww = colMeans(sigma) - mean(sigma)/2;
	k = 0;
	err = 1;
	while(err > x_tol && k < k_max) {
		# Update ww
		xx = rowSums(sigma * aabb) - ss2;
		ax1 = sum(aa2 * xx);
		bx1 = sum(bb2 * xx);
		ww2 = xx * cc + (aa2 * (cbb * ax1 - cab1 * bx1) + bb2 * (caa * bx1 -
			cab1 * ax1)) / (cab1^2 - caa * cbb);
		# Update sigma
		sigma2 = SS2 + ww2 + rep(ww2, each = p);
		oo = diag(sigma2);
		sigma2 = (sigma2 > lambda2) * (sigma2 - lambda2) + 
			(sigma2 < - lambda2) * (sigma2 + lambda2);
		diag(sigma2) = oo;
		# Check convergence
		err = max(abs(sigma2 - sigma)/(abs(sigma) + 1)); 
		k = k + 1;
		sigma = sigma2;
	};
	#if(k >= k_max) {
	#	cat("WARNING of fastcclasso_sub:\n", "\tMaximum Iteration:", k_max, 
	#		"&& Relative error:", err, "!\n");
	#};
	return(list(sigma = sigma, ww = ww2));
};
#---------------------------------------
boot_fastCCLasso <- function(xx2,sigma_hat,lambda,aa,bb, 
                             n_boot = 100,
                             max_iter=200, 
                             stop_eps=1e-6) {
  n <- nrow(xx2);
  p <- ncol(xx2);
  
  # Store the result of bootstrap
  cors_boot <- matrix(0, nrow = p * (p - 1)/2, ncol = n_boot + 1);
  vars_boot <- matrix(0, nrow = p, ncol = n_boot + 1);
  cors_mat <- matrix(0, p, p);
  ind_low <- lower.tri(cors_mat);
  
  # Bootstrap procedure
  sam_boot <- matrix(sample(1:n, size = n * n_boot, replace = T), 
                     ncol = n_boot);
  for(k in 1:n_boot) {
    ind_samp <- sam_boot[, k];
    S_samp <- stats::var(xx2[ind_samp,])
    cov_est<- fastcclasso_sub(lambda, SS2=S_samp, 
                               aa=aa, bb=bb,
                               k_max = 200, x_tol = stop_eps);
    vars_boot[, k] <- diag(cov_est$sigma);
    Is <- 1 / sqrt(vars_boot[, k]);
    cor_est <- Is * cov_est$sigma * rep(Is, each = p);
    cors_boot[, k] <- cor_est[ind_low];
  }
  
  vars_boot[, n_boot + 1] <- diag(sigma_hat);
  Is <- 1 / sqrt(vars_boot[, n_boot + 1]);
  cor_est <- Is * sigma_hat * rep(Is, each = p);  
  cors_boot[, n_boot + 1] <- cor_est[ind_low]; 
  
  # Variance estimation via bootstrap
  vars2 <- rowMeans(vars_boot);
  cors2mod <- rowMeans(cors_boot);
  cors2_mat <- diag(p);
  cors2_mat[ind_low] <- cors2mod;
  cors2_mat <- t(cors2_mat);
  cors2_mat[ind_low] <- cors2mod;
  # P-values with null distribution of correlation estimations of absolute data
  p_vals <- pt(cors2mod * sqrt((n - 2) / (1 - cors2mod^2)), df = n - 2);
  p_vals <- ifelse(p_vals <= 0.5, p_vals, 1 - p_vals);
  pval_mat <- diag(p);
  pval_mat[ind_low] <- p_vals;
  pval_mat <- t(pval_mat);
  pval_mat[ind_low] <- p_vals;
  
  return(list(var_w = vars2, cor_w = cors2_mat, p_vals = pval_mat));    
}
#-------------------------------------------------------------------------------
# call all methods
Callallmethods <- function(method,xMat,cv_k, lambda_min_ratio=1e-4,
                           Edge_eps=1e-4){
  # xMat: compositional data
  p <- dim(xMat)[2];
  S <- var(log(xMat) - rowMeans(log(xMat)));
  lambda_max <- max(max(S - diag(p)), -min(S - diag(p)));
  lambda_min <- lambda_min_ratio * lambda_max;
  lambda_int <- c(lambda_min,lambda_max);
  
  # method  
  if(method=="fastCCLasso"){
    begin_time <- proc.time();
    result <- fastCCLasso(xx=xMat,lam_min_ratio = lambda_min_ratio, 
                          k_cv = cv_k, k_max = 20);
    end_time <- proc.time();
    result_cor <- result$rho;
  }else if(method=="SparCC"){
    begin_time <- proc.time();
    result <- compute_corr_mod(fracs=xMat, iter=10, th=0.1);
    end_time <- proc.time();
    result_cor <- result$Cor.mat;
  }else if(method=="CCLasso"){    
    begin_time <- proc.time();
    result <- cclasso(xMat, counts = FALSE, pseudo = 0.5, k_cv = cv_k, 
                      lam_int = lambda_int, k_max = 20);
    end_time <- proc.time();	
    result_cor <- result$cor_w;
  }else if(method=="COAT"){
    begin_time <- proc.time();
    result <- coat(xMat, nFoler = cv_k, soft = 1);
    end_time <- proc.time();
    result_cor <- result$corr;
  }else{
    return(message("This method is not exist, please check it! "));
  } 
  result_cor[abs(result_cor)< Edge_eps] <- 0 ;
  
  return(list( runtime=as.numeric((end_time-begin_time)[3]),
               est_lower=result_cor[lower.tri(result_cor)] ,
               cor_est=result_cor));
}
#-------------------------------------------------------------------------------
# Aim     : Correlation inference for compositional data through lasso
# Referece: Fang, H., Huang, C., Zhao, H., and Deng, M. (2015),“CCLasso: 
#           correlation inference for compositional data through Lasso,” 
#           Bioinformatics, 31, 3172–3180.
# Rcode   : https://github.com/huayingfang/CCLasso/tree/master/R
#-------------------------------------------------------------------------------
# Author : Fang Huaying (Peking University)
# Email  : hyfang@pku.edu.cn
# Date   : 2016-01-08
# Version: 2.0
#-------------------------------------------------------------------------------
# Main function: cclasso(x, counts = FALSE, pseudo = 0.5, k_cv = 3, 
#                        lam_int = c(1e-4, 1), k_max = 20, n_boot = 20) 
#
#  Input:
#           x ------ n x p data matrix (row/column is sample/variable)
#                    n samples & p compositional variables
#      counts ------ Is the compositional data matrix a count matrix? 
#                    Default: FALSE
#      pseudo ------ pseudo count if counts = TRUE
#                    Default: 0.5
#        k_cv ------ folds of cross validation
#                    Default: 3     
#     lam_int ------ tuning parameter interval
#                    Default: [1e-4, 1]
#       k_max ------ maximum iterations for golden section method
#                    Default: 20
#      n_boot ------ Bootstrap times
#                    Default: 20
#  Output: 
#      A list structure contains:
#       var_w ------ variance estimation
#       cor_w ------ correlation estimation
#      p_vals ------ p-values for elements of cor_w equal 0 or not
#      lambda ------ final tuning parameter
#     info_cv ------ information for cross validation
#-------------------------------------------------------------------------------
cclasso <- function(x, counts = FALSE, pseudo = 0.5, k_cv = 3, 
                    lam_int = c(1e-4, 1), k_max = 20, n_boot = 20) {
  n <- nrow(x);
  p <- ncol(x);

  if(counts) {
    x <- x + pseudo;
    x <- x / rowSums(x);
  }
  x <- log(x);
  vx2 <- stats::var(x);

  # Diagonal weight for loss function
  rmean_vx2 <- rowMeans(vx2);
  wd <- 1/diag(vx2 - rmean_vx2 - rep(rmean_vx2, each = p) + mean(rmean_vx2));
  wd2 <- sqrt(wd);
  
  # Some global parameters for optimization with single lambda
  rho <- 1;
  u_f <- eigen(diag(p) - 1/p)$vectors;
  wd_u <- (t(u_f) %*% (wd * u_f))[-p, -p];
  wd_u_eig <- eigen(wd_u);
  d0_wd <- 1 / ( (rep(wd_u_eig$values, each = p-1) + wd_u_eig$values) / 
    (2 * rho) + 1 );
  u0_wd <- wd_u_eig$vectors;

  # Golden section method for the selection of lambda (log10 scale)
  sigma <- vx2;
  lam_int2 <- log10(range(lam_int));
  a1 <- lam_int2[1]; 
  b1 <- lam_int2[2];
  # Store lambda and corresponding cross validation's loss
  lams <- NULL; 
  fvals <- NULL;
  # Two trial points in first 
  a2 <- a1 + 0.382 * (b1 - a1); 
  b2 <- a1 + 0.618 * (b1 - a1);
  fb2 <- cv_loss_cclasso(lambda2 = 10^b2 / rho, x = x, k_cv = k_cv, 
    sigma = sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
    wd2 = wd2);
  lams <- c(lams, b2); 
  fvals <- c(fvals, fb2$cv_loss);
  fa2 <- cv_loss_cclasso(lambda2 = 10^a2 / rho, x = x, k_cv = k_cv, 
    sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
    wd2 = wd2);
  lams <- c(lams, a2); 
  fvals <- c(fvals, fa2$cv_loss);
  # Error tolerance for convergence
  err_lam2 <- 1e-1 * max(1, lam_int2);
  err_fval <- 1e-4;
    
  err <- b1 - a1;
  k <- 0;
  while(err > err_lam2 && k < k_max) {
    fval_max <- max(fa2$cv_loss, fb2$cv_loss);

    if(fa2$cv_loss > fb2$cv_loss) {
      a1 <- a2;      
      a2 <- b2;
      fa2 <- fb2;
      b2 <- a1 + 0.618 * (b1 - a1);
      fb2 <- cv_loss_cclasso(lambda2 = 10^b2 / rho, x = x, k_cv = k_cv, 
        sigma = fa2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
        wd2 = wd2);

      lams <- c(lams, b2); 
      fvals <- c(fvals, fb2$cv_loss);
    } else {
      b1 <- b2;
      b2 <- a2;
      fb2 <- fa2;
      a2 <- a1 + 0.382 * (b1 - a1);
      fa2 <- cv_loss_cclasso(lambda2 = 10^a2 / rho, x = x, k_cv = k_cv, 
        sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
        wd2 = wd2);

      lams <- c(lams, a2); 
      fvals <- c(fvals, fa2$cv_loss);
    }
    fval_min <- min(fa2$cv_loss, fb2$cv_loss);      

    k <- k + 1;
    err <- b1 - a1;
    if(abs(fval_max - fval_min) / (1 + fval_min) <= err_fval) {
      break;
    }
  }
  info_cv <- list(lams = lams, fvals = fvals, k = k + 2, 
    lam_int = 10^c(a1, b1)); 
  #if(a1 == lam_int2[1] || b1 == lam_int2[2]) {
  #  cat("WARNING:\n", "\tOptimal lambda is near boundary! ([", 10^a1, ",", 
   #   10^b1, "])\n", sep = "");
  #}

  lambda <- 10^((a2 + b2)/2);
  # Bootstrap for cclasso
  lambda2 <- lambda / rho;
  info_boot <- boot_cclasso(x = x, sigma = fb2$sigma, lambda2 = lambda2, 
    n_boot = n_boot, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);

  return(list(var_w = info_boot$var_w, cor_w = info_boot$cor_w, 
    p_vals = info_boot$p_vals, lambda = lambda, info_cv = info_cv));
}
#---------------------------------------
# Bootstrap for cclasso
boot_cclasso <- function(x, sigma, lambda2, n_boot = 20, 
                         wd, u_f, u0_wd, d0_wd) {
  n <- nrow(x);
  p <- ncol(x);

  # Store the result of bootstrap
  cors_boot <- matrix(0, nrow = p * (p - 1)/2, ncol = n_boot + 1);
  vars_boot <- matrix(0, nrow = p, ncol = n_boot + 1);
  cors_mat <- matrix(0, p, p);
  ind_low <- lower.tri(cors_mat);

  # Bootstrap procedure
  sam_boot <- matrix(sample(1:n, size = n * n_boot, replace = T), 
    ncol = n_boot);
  for(k in 1:n_boot) {
    ind_samp <- sam_boot[, k];
    sigma2 <- cclasso_sub(sigma = sigma, vx = var(x[ind_samp, ]), 
      lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
        
    vars_boot[, k] <- diag(sigma2);
    Is <- 1 / sqrt(vars_boot[, k]);
    cors_mat <- Is * sigma2 * rep(Is, each = p);
    cors_boot[, k] <- cors_mat[ind_low];
  }
  Is <- 1 / sqrt(diag(sigma));
  cors_mat <- sigma * Is * rep(Is, each = p);
  cors_boot[, n_boot + 1] <- cors_mat[ind_low];
  vars_boot[, n_boot + 1] <- diag(sigma);
  
  # Variance estimation via bootstrap
  vars2 <- rowMeans(vars_boot);
  # Correlations' relationship for artificial null sample
  tol_cor <- 1e-3;
  sam_art0 <- matrix(rnorm(n * p), nrow = n) * rep(sqrt(vars2), each = n);
  cors_art0 <- cor(sam_art0)[ind_low];
  sam_art <- sam_art0 - log(rowSums(exp(sam_art0)));
  sigma_art <- cclasso_sub(sigma = sigma, vx = var(sam_art), 
    lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
  Is <- 1 / sqrt(diag(sigma_art));
  cors_mat <- Is * sigma_art * rep(Is, each = p);
  cors_art2 <- cors_mat[ind_low];
  # Bias of estimation between absolute data and cclasso of compositional data
  cors0m2 <- log( ((1 + cors_art0) * (1 - cors_art2)) / ((1 + cors_art2) * 
    (1 - cors_art0)) );
  tmp <- abs(cors_art2) >= tol_cor;
  bias02 <- ifelse(sum(tmp), median(abs(cors0m2)[tmp]), 0);
  # Modification of estimation for cclasso
  cors2 <- log( (1 + cors_boot) / (1 - cors_boot) );    
  cors2mod <- (cors_boot >= tol_cor) * (cors2 + bias02) + 
    (cors_boot <= -tol_cor) * (cors2 - bias02);
  cors2mod <- 1 - rowMeans(2 / (exp(cors2mod) + 1));
  cors2_mat <- diag(p);
  cors2_mat[ind_low] <- cors2mod;
  cors2_mat <- t(cors2_mat);
  cors2_mat[ind_low] <- cors2mod;
  # P-values with null distribution of correlation estimations of absolute data
  p_vals <- pt(cors2mod * sqrt((n - 2) / (1 - cors2mod^2)), df = n - 2);
  p_vals <- ifelse(p_vals <= 0.5, p_vals, 1 - p_vals);
  pval_mat <- diag(p);
  pval_mat[ind_low] <- p_vals;
  pval_mat <- t(pval_mat);
  pval_mat[ind_low] <- p_vals;

  return(list(var_w = vars2, cor_w = cors2_mat, p_vals = pval_mat));    
}
#---------------------------------------
# cross validation's loss of cclasso for single lambda
cv_loss_cclasso <- function(lambda2, x, k_cv, sigma, 
                            wd, u_f, u0_wd, d0_wd, wd2) {
  n <- nrow(x);
  p <- ncol(x);

  n_b <- floor(n / k_cv);
  cv_loss <- 0;
  for(k in 1:k_cv) {
    itest <- (n_b * (k - 1) + 1):(n_b * k);
    vxk <- stats::var(x[itest, ]);
    vx2k <- stats::var(x[-itest, ]);

    sigma <- cclasso_sub(sigma = sigma, vx = vx2k, lambda2 = lambda2, 
      wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
    
    dsig <- sigma - vxk;
    tmp <- rowMeans(dsig);
    dsig <- dsig - tmp - rep(tmp, each = p) + mean(tmp);
    cv_loss <- cv_loss + base::norm(wd2 * dsig, "F")^2;
  }

  return(list(cv_loss = cv_loss, sigma = sigma));
}
#---------------------------------------
# cclasso for single lambda
cclasso_sub <- function(sigma, vx, lambda2, 
                        wd, u_f, u0_wd, d0_wd, 
                        k_max = 200, x_tol = 1e-4) {
  p <- ncol(sigma);
  sigma2 <- sigma;
  LAMBDA <- matrix(0, p, p);
  lambda2 <- matrix(lambda2, p, p);
  diag(lambda2) <- 0;

  k <- 0;
  err <- 1;
  while(err > x_tol && k < k_max) {
    # Update sigma
    x_sigma <- t(u_f) %*% ((sigma2 - vx) - LAMBDA) %*% u_f;
    x_sigma[-p,-p] <- u0_wd %*% ((t(u0_wd) %*% x_sigma[-p, -p] %*% u0_wd) * 
      d0_wd) %*% t(u0_wd);
    sigma_new <- vx + u_f %*% x_sigma %*% t(u_f);
    # Update sigma2
    A <- LAMBDA + sigma_new;
    sigma2_new <- (A > lambda2) * (A - lambda2) + (A < -lambda2) * 
      (A + lambda2);
    # Update Lambda
    LAMBDA <- LAMBDA + (sigma_new - sigma2_new);

    err <- max( abs(sigma_new - sigma)/(abs(sigma) + 1), 
      abs(sigma2_new - sigma2)/(abs(sigma2) + 1) ); 
    k <- k + 1;
    sigma <- sigma_new;
    sigma2 <- sigma2_new;
  }

  #if(k >= k_max) {
  #  cat("WARNING of cclasso_sub:\n", "\tMaximum Iteration:", k_max, 
  #    "&& Relative error:", err, "!\n");
  #}
  
  return(sigma);
}
#-------------------------------------------------------------------------------
# Aim :     Correlation inference for compositional data via composition-adjusted thresholding
# Referece: Cao, Y., Lin, W. and Li, H. (2019). Large covariance estimation for 
#           compositional data via composition-adjusted thresholding, Journal of
#           the American Statistical Association 114(526): 759–772.
# Rcode:    https://github.com/yuanpeicao/COAT
#-------------------------------------------------------------------------------
#  COAT estimate
#  Input:
#           x ------ n x p composition data matrix (row/column is sample/variable)
#     nFolder ------ number of the foler in cross validation
#        soft ------ soft = 1: soft thresholding; soft = 0: hard thrsholding
#  Output:
#       sigma ------ covariance estimation based on adaptive thresholding 
#        corr ------ correlation estimation based on adaptive thresholding
#        time ------ execution time
#-------------------------------------------------------------------------------
coat <- function(x, nFoler = 5, soft = 1){
  startTime <- proc.time()
  p <- ncol(x)
  clrX <- log(x) - rowSums(log(x)) %*%matrix(1,1,p) / p
  coatPred <- adaptThresoldCov(clrX, soft = soft)
  sigma <- coatPred$sigma
  corr <- coatPred$corr
  exeTimeClass <- proc.time() - startTime
  exeTime <- as.numeric(exeTimeClass[3])
  return(list(sigma = sigma, corr = corr, time = exeTime))
}

#---------------------------------------
#  Adaptive thresholding estimation of cov(x)
#  Input:
#           x ------ n x p data matrix (row/column is sample/variable)
#     nFolder ------ number of the foler in cross validation
#        soft ------ soft = 1: soft thresholding; soft = 0: hard thrsholding
#  Output:
#       sigma ------ covariance estimation based on adaptive thresholding 
#        corr ------ correlation estimation based on adaptive thresholding 
#---------------------------------------
adaptThresoldCov <- function(x, nFolder = 5, soft = 1){
  n <- nrow(x)
  p <- ncol(x)
  # Set the grid for the choice of tuning parameter
  nGrid <- 100
  gridInfo <- adaptThresholdRange(x)
  grid <- gridInfo$lower + (gridInfo$upper - gridInfo$lower)*rep(1:nGrid)/nGrid
  # Multi-folder cross validation
  part <- 1 + sample(c(1:n))%%nFolder
  error <- matrix(0, nFolder, nGrid)
  for (i in 1:nFolder){
    xTest <- x[which(part == i),]
    xTrain <- x[which(part != i),]
    gridInfoTrain <- adaptThresholdRange(xTrain)
    covTest <- cov(xTest)*(n-1)/n
    for (j in 1:nGrid){
      sigmaTrain <- adaptThreshold(gridInfoTrain$cov,gridInfoTrain$theta,grid[j],soft)
      error[i,j] <- (norm(sigmaTrain-covTest, "F"))
    }
  }
  errorSum <- colSums(error)
  lambda <- grid[which(errorSum == min(errorSum))][1]
  sigma <- adaptThreshold(gridInfo$cov,gridInfo$theta,lambda,soft)
  corr <- diag(diag(sigma)^(-0.5))%*%sigma%*%diag(diag(sigma)^(-0.5))
  return(list(sigma = sigma, corr = corr))
}
#---------------------------------------
#  Range of the tuning parameter
#  Input:
#           x ------ n x p data matrix (row/column is sample/variable)
#  Output:
#      A list structure contains:
#       upper ------ upper bound of tuning parameter
#       lower ------ lower bound of tuning parameter
#         cov ------ sample covariance of x
#       theta ------ sample variance of covariance
#---------------------------------------
adaptThresholdRange <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  cov <- cov(x)*(n-1)/n
  centered.x <- scale(x, scale = FALSE)
  theta <- (t(centered.x)^2)%*%(centered.x^2)/n - cov^2
  delta <- cov/(theta^0.5)
  delta <- abs(delta - diag(diag(delta)))
  upper <- max(delta)
  lower <- min(delta[which(delta != 0)])
  return(list(upper = upper, lower = lower, theta = theta, cov = cov))
}
#---------------------------------------
#  Apply adaptive thresholding to the sample covariance
#  Input:
#           cov ------ p x p covariance matrix
#         theta ------ p x p variance of covariance matrix
#        lambda ------ tuning parameter
#          soft ------ soft = 1: soft thresholding; soft = 0: hard thrsholding
#  Output:
#         sigma ------ p x p matrix, adaptive thresholding result
#---------------------------------------
adaptThreshold <- function(cov,theta,lambda,soft){
  covOffDiag <- cov - diag(diag(cov))
  thetaOffDiag <- theta - diag(diag(theta))
  sigmaTmp <- abs(covOffDiag) - lambda*thetaOffDiag^0.5
  sigmaTmp[which(sigmaTmp < 0)] <- 0
  if (soft == 1){
    sigma <- diag(diag(cov)) + sigmaTmp*sign(covOffDiag)
  }else{
    sigma <- cov
    sigma[which(sigmaTmp < 1e-10)] <- 0
    sigma <- sigma + diag(diag(cov))
  }
  return(sigma)
}
#-------------------------------------------------------------------------------
# Aim     : Correlation inference for compositional data
# Referece: Friedman, J. and Alm, E. J. (2012). Inferring correlation networks from 
#           genomic survey data, PLoS Computational Biology 8(9): e1002687.
# Rcode   : https://github.com/MPBA/r-sparcc/tree/master/R
#-------------------------------------------------------------------------------
# compute_corr_mod: modify the compute_corr function for fraction data
#          compute_corr_mod(fracs, iter=10, th=0.1) 
#
#  Input:
#       fracs ------ n x p data matrix (row/column is sample/variable)
#                    n samples & p compositional variables
#        iter ------ maximum iterations 
#                    Default: 10
#          th ------ thresholding
#                    Default: 0.1
#  Output: 
#      A list structure contains:
#       Vbase ------ variance estimation
#     Cor.mat ------ correlation estimation
#     Cov.mat ------ covariance estimation
#-------------------------------------------------------------------------------
## count matrix x should be samples on the rows and OTUs on the colums,
## assuming dim(x) -> samples by OTUs
sparcc <- function(x, max_iter=20, th=0.1, exiter=10){
  xdim <- dim(x)
  Vlist <- matrix(NA,nrow=xdim[2],ncol=max_iter)
  Corlist <- array(,dim=c(max_iter, xdim[2], xdim[2]))
  Covlist <- array(,dim=c(max_iter, xdim[2], xdim[2]))
  
  ## Cycle max_iter times for variability in variance estimation
  for (i in 1:max_iter){
    #cat("Iteration: %d\n",i)
    tmpres <- compute_corr(x, iter=exiter, th=th)
    Vlist[,i] <- tmpres[["Vbase"]]
    Corlist[i,,] <- tmpres[["Cor.mat"]]
    Covlist[i,,] <- tmpres[["Cov.mat"]]
  }
  
  ## Compute variance basis and correlation
  vdef <- apply(Vlist,1,median)
  cor_def <- apply(Corlist,2:3,median)
  
  ## Square root variances
  vdefsq <- vdef**0.5
  
  ## Compute covariance
  ttmp <- cor_def * vdefsq
  cov_def <- t(ttmp) * vdefsq
  
  ## Uncomment following lines for an alternative method
  ## x <- matrix(vdefsq,ncol=50,nrow=50, byrow=TRUE)
  ## y <- t(x)
  ## cov_def <- cor_def * x * y
  
  return(list(CORR=cor_def, COV=cov_def, VBASIS=vdef))
}

compute_corr <- function(x, iter=10, th=0.1){

  ## Compute relative fraction from dirichlet distribution
  ## NB think on different normalization for improvements
  fracs <- counts2frac(x)

  ## Compute the variation matrix
  V <- variation_mat(fracs)
  
  ## Compute the Sparcc correlation
  ## Initialize matrices
  ll1 <- basis_var(fracs, V)
  ll2 <- cor_from_basis(V, ll1[["Vbase"]])
  excluded <- NULL
  
  for (i in 1:iter){
    ## Search for excluded pairs
    ## ll2[[1]] -> Cor.mat
    ll3 <- exclude_pairs(ll2[["Cor.mat"]], ll1[["M"]], th = th, excluded = excluded)
    excluded <- ll3[["excluded"]]
    if (!ll3[["flag"]]){
      ll1 <- basis_var(fracs, V, M=ll3[["M"]], excluded=excluded)
      ll2 <- cor_from_basis(V, ll1[["Vbase"]])
    }
  }
  
  return(list(Vbase=ll1[["Vbase"]], Cor.mat=ll2[["Cor.mat"]], Cov.mat=ll2[["Cov.mat"]]))
}
#----------------------------------------------------------------------------------------
# modity compute_corr for fraction data
compute_corr_mod <- function(fracs, iter=10, th=0.1){ 
  ## Compute relative fraction from dirichlet distribution
  ## NB think on different normalization for improvements
  #fracs <- counts2frac(x)
  
  ## Compute the variation matrix
  V <- variation_mat(fracs)
  
  ## Compute the Sparcc correlation  
  ## Initialize matrices
  ll1 <- basis_var(fracs, V)
  ll2 <- cor_from_basis(V, ll1[["Vbase"]])
  excluded <- NULL
  
  for (i in 1:iter){
    ## Search for excluded pairs
    ## ll2[[1]] -> Cor.mat
    ll3 <- exclude_pairs(ll2[["Cor.mat"]], ll1[["M"]], th = th, excluded = excluded)
    excluded <- ll3[["excluded"]]
    if (!ll3[["flag"]]){
      ll1 <- basis_var(fracs, V, M=ll3[["M"]], excluded=excluded)
      ll2 <- cor_from_basis(V, ll1[["Vbase"]])
    }
  }
  
  return(list(Vbase=ll1[["Vbase"]], Cor.mat=ll2[["Cor.mat"]], Cov.mat=ll2[["Cov.mat"]]))
}

counts2frac <- function(x, method="dirichlet"){
  xsize <- dim(x)
  fracs <- matrix(1/xsize[2], nrow=xsize[1], ncol=xsize[2])
  if (method=="dirichlet"){
    fracs_t <- apply(x,1,function(y){rdirichlet(1,y + 1)})
    fracs <- t(fracs_t)
  }
  return(fracs)
}

variation_mat <- function(fracs){
  ## Initialize variation matrix
  V <- matrix(NA, ncol=dim(fracs)[2], nrow=dim(fracs)[2])
  ## Compute log for each OTU
  tmplog <- apply(fracs,2,log)
  idx <- combn(1:dim(fracs)[2],2)
  
  ## create matrix Ti,j
  ttmp <- tmplog[,idx[1,]] -   tmplog[,idx[2,]]
  
  ## Compute Variance
  vartmp <- apply(ttmp,2, var)
  
  ## Fill the variance matrix
  for (i in 1:length(vartmp)){
    V[idx[1,i],idx[2,i]] <- V[idx[2,i],idx[1,i]] <- vartmp[i]
  }
  diag(V) <- 1
  return(V)
}

basis_var <- function(fracs, V, Vmin=1e-4, excluded=NULL, Covmat=NULL, M=NULL){
  
  Vsize <- dim(V)
  Vvec <- apply(V,1,sum)
  
  ## Initialize Covmat matrix
  if (is.null(Covmat))
    Covmat <- matrix(0, nrow=Vsize[1], ncol=Vsize[2])
  
  Covvec <- apply(Covmat - diag(Covmat),1,sum)
  ## Initialize M matrix 
  if (is.null(M)){
    M <- matrix(1, nrow=Vsize[1], ncol=Vsize[2])
    diag(M) <- Vsize[1] - 1
  }
  Minv <- solve(M)
  Vbase <- Minv %*% (Vvec + 2*Covvec)
  Vbase[Vbase<0] <- Vmin
  return(list(Vbase=Vbase, M=M))
}

cor_from_basis <- function(V, Vbase){
  ## Compute the correlation from variation matrix and basis variations
  
  p <- dim(Vbase)[1]
  Cor.mat <- diag(rep(1,p))
  Cov.mat <- diag(Vbase[,1])
  
  idx <- combn(p,2)
  
  for (i in 1:(p-1)){
    idxslice <- idx[1,]==i
    cov.tmp <- .5 * (Vbase[i] + Vbase[idx[2,idxslice]] - V[i,idx[2,idxslice]])
    denom <- sqrt(Vbase[i]) * sqrt(Vbase[idx[2,idxslice]])
    cor.tmp <- cov.tmp / denom
    abscor <- abs(cor.tmp)
    if (any(abscor > 1)){
      idxthr <- abscor > 1
      
      ## Set the max correlation to -1,1
      cor.tmp[idxthr] <- sign(cor.tmp[idxthr])
      
      ## Compute the covariance basis
      cov.tmp[idxthr] <- cor.tmp[idxthr] * denom[idxthr] 
    }
    
    ## Fill the cor and cov matrix
    Cor.mat[i,idx[2,idxslice]] <- Cor.mat[idx[2,idxslice],i] <- cor.tmp
    Cov.mat[i,idx[2,idxslice]] <- Cov.mat[idx[2,idxslice],i] <- cov.tmp
  }
  return(list(Cor.mat=Cor.mat, Cov.mat=Cov.mat))
}


exclude_pairs <- function(Cor.mat, M, th=0.1, excluded=NULL){
  
  flag <- FALSE
  
  ## Remove autocorrelation
  cor.tmp <- abs(Cor.mat)
  diag(cor.tmp) <- diag(cor.tmp) - diag(Cor.mat)
  
  if (!is.null(excluded))
    cor.tmp[excluded,] <- 0
  
  ## Search highly correlated pairs
  mm <- max(cor.tmp)
  idxtorm <- which(cor.tmp==mm, arr.ind=TRUE)
  
  if (mm > th){
    ## Subtract 1 in in the M matrix where found highly correlated pairs
    for (i in 1:dim(idxtorm)[1]){
      M[idxtorm[i,1],idxtorm[i,2]] <- M[idxtorm[i,1],idxtorm[i,2]] - 1
    }
    
    ## Subtract one to the diagonal
    dd <- diag(M)[unique(c(idxtorm))]
    diag(M)[unique(c(idxtorm))] <- dd - 1
    excluded <- rbind(excluded, idxtorm)
  } else {
    excluded <- excluded
    flag <- TRUE
  }
  return(list(M=M, excluded=excluded, flag=flag))
}

normalize_matrix <- function(x){
  # Normalize by total count
  x <- apply(x,2,function(y){
    y/sum(y)})
  
  return(x)
}
#-------------------------------------------------------------------------------
