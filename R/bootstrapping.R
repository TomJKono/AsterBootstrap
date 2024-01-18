# Functions related to running the parameteric bootstrap for variance
# estimation.

#' Estimate bootstrpped additive genetic variance for fitness with parallel
#' processing
single_bootstrap_VaW_par <- function(mod_obj, long_dat, model_spec, r_eff,
    blocking, n_iter, VaW_fun, boot_effect, seed, typical_ind, typical_fe_level,
    pkg_dir, n_cores=1, forumula=NULL, quiet=TRUE) {
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
    # Tell the cluster to load the AsterBootstrap package
    parallel::clusterCall(
        clust,
        function(pkg_dir) {
            devtools::load_all(pkg_dir, quiet=TRUE)
        },
        pkg_dir=pkg_dir
    )
    # Stop the cluster on exit, regardless of success or error
    on.exit(parallel::stopCluster(clust))
    # define a local %dopar%
    `%dopar%` <- foreach::`%dopar%`
    # Define a function that will obtain a single replicate estimate of the
    # additive genetic variance for fitness. In essence, this will fit a new
    # model with resampled values for the random effects, then use that
    # model to estimate additive genetic variance.
    single_boot <- function(boot_mod_obj, boot_dat, boot_model_spec, boot_reff,
        boot_blocking, vaw_function, boot_effect, boot_typical_ind,
        boot_typical_fe_level) {
        # Extract the original estimates of the 'alpha' and 'sigma' parameters
        # from the model object. 'alpha' are the MLE for the fixed effects, and
        # 'sigma' are the square roots of the MLE of the variance components.
        alpha_est <- boot_mod_obj$alpha
        sigma_est <- boot_mod_obj$sigma
        # Extract the penalized likelihood estimates for the random effects
        b_hat <- boot_mod_obj$b
        # Identify which random effects are the ones to replace for
        # bootstrapping.
        #   E.g., if we are bootstrapping the paternal random effect, first
        #   find effects that have 'Pat' in the name. Then, "invert" the
        #   matches because we will keep non-Pat effects as their original
        #   estimated values and replace Pat effects with simulated ones
        beff_idx <- !grepl(boot_effect, attributes(b_hat)$names)
        # Make a model matrix from the fixed and random effects
        obj_f_eff <- boot_mod_obj$fixed
        obj_r_eff <- boot_mod_obj$random
        n_rand <- sapply(obj_r_eff, ncol)
        modmat <- cbind(obj_f_eff, Reduce(cbind, obj_r_eff))
        # Check the random effects again - we need to drop any corresponding
        # entries that are missing in the original data, lest we get a mismatch
        # of dimension.
        missing_res <- sapply(names(boot_reff), function(re) {
            miss_rows <- which(is.na(boot_dat[[boot_reff[[re]]]]))
            return(as.numeric(miss_rows))
        })
        # Keep track of which rows have missing random effects. Remove them
        # the bootstrapped data if necessary
        missing_res <- unique(Reduce(c, missing_res))
        if(length(missing_res) > 0) {
            boot_dat <- boot_dat[-missing_res,]
        }
        # Use the sigma values from the original random effects Aster model as
        # a starting point for simulating new values
        reff_start <- rep(sigma_est, times=n_rand)
        # Make up new values for the random effects by multiplying by ~N(0, 1)
        sim_factor <- rnorm(sum(n_rand), mean=0, sd=1)
        b_star <- reff_start * sim_factor
        # Then, overwrite the simulated random effects with the original
        # estimates, but only for the random effects that are not being
        # bootstrapped.
        b_star[beff_idx] <- b_hat[beff_idx]
        # Put the fixed and random effects together to apply to the bootstrap
        # data
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
        # Use raster() to generate new data from an Aster distribution with the
        # new parameters and the existing graphical model
        y_star <- aster::raster(
            theta_star,
            boot_model_spec$pred,
            boot_model_spec$fam,
            boot_mod_obj$obj$root)
        y_star <- as.vector(y_star)
        # Use the blocking factors supplied to make a new model with the
        # proper blocking variables
        if(is.null(boot_blocking)) {
            mod <- resp ~ varb
        } else {
            block_string <- paste(boot_blocking, collapse="+")
            mod_string <- paste("resp~varb+fit:(", block_string, ")", sep="")
            mod <- as.formula(mod_string)
        }
        # Fit a new model. Suppress warnings - the user would have seen any
        # relevant warnings for the single-iteration model fitting
        boot_dat$resp <- y_star
        rout_boot <- random_effects_aster(
                boot_dat,
                boot_model_spec,
                r_eff=boot_reff,
                formula=mod,
                quiet=quiet)
        # call to function for VaW estimation
        vaw <- suppressWarnings(
            vaw_function(
                reaster_obj=rout_boot,
                sire_label="Pat",
                effect_label="Father",
                model_spec=boot_model_spec,
                fixed_eff_idx=boot_typical_fe_level,
                typical_ind_idx=boot_typical_ind)
            )
        return(vaw)
    }
    # Apply this bootstrapping function across multiple cores now
    boot_est <- foreach::foreach(
        i=1:n_iter,
        combine=rbind,
        .multicombine=TRUE,
        .inorder=FALSE) %dopar% {
            set.seed(seed + i)
            sb_est <- NULL
            while(is.null(sb_est)) {
                try(
                    sb_est <- single_boot(
                        boot_mod_obj=mod_obj,
                        boot_dat=long_dat,
                        boot_model_spec=model_spec,
                        boot_reff=r_eff,
                        boot_blocking=blocking,
                        vaw_function=VaW_fun,
                        boot_effect=boot_effect,
                        boot_typical_ind=typical_ind,
                        boot_typical_fe_level=typical_fe_level)
                    )
            }
            return(sb_est)
        }
    return(boot_est)
}
