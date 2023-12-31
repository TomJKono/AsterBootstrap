# Functions related to running the parameteric bootstrap for variance
# estimation.

#' Estimate bootstrpped additive genetic variance for fitness with parallel
#' processing
bootstrap_VaW_par <- function(mod_obj, long_dat, model_spec, reff, blocking,
    n_iter, VaW_fun, n_cores=1, forumula=NULL, quiet=TRUE) {
    # Check that the "parallel" and "doParallel" packages are available
    if(!(requireNamespace("parallel", quietly=TRUE)
       || requireNamespace("doParallel", quietly=TRUE))) {
    stop('"parallel" and "doParallel" are required for parallel bootstrap.',
        call.=TRUE)
    }
    # Print a warning if someone is asking for more cores than there are
    # available on the system.
    tot_cores <- parallel::detectCores()
    if(n_cores > tot_cores) {
        warning(
            paste(
                "Requested more cores than R thinks you have on your ",
                "workstation! R thinks you have ",
                tot_cores,
                " cores available. This is not an error.",
                sep=""),
            immediate.=TRUE)
    }
    # heck that the number of iterations requested is greater than 0.
    if(!n_iter > 0) {
        stop("'n_iter' must be an integer greater than 0.", call.=TRUE)
    }
    # Cast the number of iterations to integer
    n_iter <- as.integer(n_iter)
    # Check that the number of cores requested is at least 1
    if(!n_cores >=1) {
        stop("'n_cores' must be an integer greater than or equal to 1.",
            call=TRUE)
    }
    n_cores <- as.integer(n_cores)
    # Make a variable that records which rows in the long data correspond to
    # the observations for the fitness variable
    if(!model_spec$fit %in% unique(long_dat$varb)) {
        stop(paste(
            "Fitness variable ",
            model_spec$fit,
            " does not exist in the long data! The long data has ",
            paste(unique(long_dat$varb), collapse=", "),
            sep=""),
        call.=TRUE)
    }
    # Try to get the 'alpha' and 'sigma' list elements from the model object.
    # If these do not exist, then we do not have a random-effects Aster model.
    elems <- c("alpha", "sigma", "random")
    if(!all(elems %in% names(mod_obj))) {
        stop(paste(
            "Random-effects elements not found in mod_obj!",
            "This indicates that mod_obj is not a random-effects Aster model.",
            "Please supply a random-effects model.",
            sep=" "),
        call.=TRUE)
    }
    # Now, let's run the bootstrap. We will use the foreach and doParallel
    # packages to distribute this across multiple local CPU cores. This
    # approach is not as efficient as mclapply(), but it is more portable for
    # those on Windows workstations.
    clust <- parallel::makeCluster(n_cores, outfile="")
    doParallel::registerDoParallel(clust)
    # Stop the cluster on exit, regardless of success or error
    on.exit(parallel::stopCluster(clust))
    # define a local %dopar%
    `%dopar%` <- foreach::`%dopar%`
    # Define a function that will obtain a single replicate estimate of the
    # additive genetic variance for fitness. In essence, this will fit a new
    # model with resampled values for the random effects, then use that
    # model to estimate additive genetic variance.
    single_boot <- function(boot_mod_obj, boot_dat, boot_model_spec, boot_reff,
        boot_blocking, vaw_function) {
        # Extract the original estimates of the 'alpha' and 'sigma' parameters
        # from the model object
        alpha_est <- boot_mod_obj$alpha
        sigma_est <- boot_mod_obj$sigma
        # Make a model matrix from the fixed and random effects
        obj_f_eff <- boot_mod_obj$fixed
        obj_r_eff <- boot_mod_obj$random
        modmat <- cbind(obj_f_eff, Reduce(cbind, obj_r_eff))
        # Extract the number of observations for the random effects - this will
        # be used in fitting a new model for a bootstrap replicate
        n_rand <- sapply(obj_r_eff, ncol)
        # Make up new values for the random effects by multiplying them by
        # a random [~N(0, 1)] value
        a_hat <- rep(sigma_est, times=n_rand)
        c_star <- rnorm(sum(n_rand), mean=0, sd=1)
        b_star <- a_hat * c_star
        eff_star <- c(alpha_est, b_star)
        # Transform the new effects to canonical parameter space so that they
        # can be used to generate new "data" from the Aster model distributions
        phi_star <- as.numeric(
            as.vector(boot_mod_obj$obj$origin) + modmat %*% eff_star)
        theta_star <- aster::astertransform(
            phi_star,
            boot_mod_obj$obj,
            to.cond="conditional",
            to.mean="canonical")
        # Use raster() to generate new data from an Aster distribution
        y_star <- aster::raster(
            theta_star,
            boot_model_spec$pred,
            boot_model_spec$fam,
            boot_mod_obj$obj$root)
        y_star <- as.vector(y_star)
        # Use the blocking factors supplied to make a new model with the
        # proper blocking variables
        if(is.null(boot_blocking)) {
            mod <- y_star ~ varb
        } else {
            block_string <- paste(boot_blocking, collapse="+")
            mod_string <- paste("y_star~varb+fit:(", block_string, ")", sep="")
            mod <- as.formula(mod_string)
        }
        # Fit a new model
        rout_boot <- random_effects_aster(
            boot_dat,
            boot_model_spec,
            reff=boot_reff,
            effects=c(alpha_est, c_star),
            sigma=sigma_est,
            formula=mod)
        # ADD HERE: call to function for VaW estimation
        vaw <- vaw_function(rout_boot, "Pat", "Father", boot_model_spec,
            1, 1)
        return(vaw)
    }
    # Apply this bootstrapping function across multiple cores now
    boot_est <- foreach::foreach(
        i=1:n_iter,
        combine=rbind,
        .multicombine=TRUE,
        .inorder=FALSE) %dopar% {
            devtools::load_all("~/Dropbox/GitHub/RIS/AsterBootstrap/")
            single_boot(
                boot_mod_obj=mod_obj,
                boot_dat=long_dat,
                boot_model_spec=model_spec,
                boot_reff=reff,
                boot_blocking=blocking,
                vaw_function=VaW_fun)
        }
    return(boot_est)
}
