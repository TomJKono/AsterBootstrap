# Functions related to fitting Aster models and generating point predictions
# from Aster models.

#' Fit a fixed-effects Aster model to a dataset
fixed_effects_aster <- function(dat, model_spec, blocking=NULL, formula=NULL,
    quiet=TRUE) {
    # If there are blocking factors specified, then coerce them to factor type
    if(!is.null(blocking)) {
        for(bf in blocking) {
            if(!bf %in% colnames(dat)) {
                stop(paste(
                    "Blocking factor '",
                    bf,
                    "' not found in the dataset!",
                    sep=""),
                call.=TRUE)
            }
            dat[[bf]] <- as.factor(dat[[bf]])
        }
    }

    # Build the formula. If the user has specified a custom formula, then use
    # that instead.
    if(is.null(formula)) {
        # The default formula will be a super basic one. It will use just the
        # response (resp), the graphical model variables (varb), the fitness
        # surrogate (fit), and whichever blocking factors were specified.
        if(is.null(blocking)) {
            mod <- resp ~ varb + fit
        } else {
            # We will only support additive "main effects" of blocking
            # factors. Interactions or more complicated linear combinations of
            # blocking factors will need a custom formula.
            block_string <- paste(blocking, collapse="+")
            mod_string <- paste("resp~varb+fit:(", block_string, ")", sep="")
            mod <- as.formula(mod_string)
        }
    } else {
        mod <- formula
    }

    # If the user asks for non-quiet output, we will print the formula
    if(!quiet) {
        cat("Formula:\n", file=stderr())
        cat(capture.output(print(mod)), file=stderr())
        cat("\n", file=stderr())
    }

    # Fit the model and return the model object
    aout <- aster::aster(
        formula=mod,
        pred=model_spec$pred,
        fam=model_spec$fam,
        varvar=varb,
        idvar=id,
        root=dat[[model_spec$initial]],
        data=dat)
    return(aout)
}


#' Fit a random-effects Aster model to a dataset
random_effects_aster <- function(dat, model_spec, r_eff=NULL, blocking=NULL,
    formula=NULL, effects=NULL, sigma=NULL, quiet=TRUE) {
    # If there are blocking factors specified, then coerce them to factor type
    if(!is.null(blocking)) {
        for(bf in blocking) {
            if(!bf %in% colnames(dat)) {
                stop(paste(
                    "Blocking factor '",
                    bf,
                    "' not found in the dataset!",
                    sep=""),
                call.=TRUE)
            }
            dat[[bf]] <- as.factor(dat[[bf]])
        }
    }
    # Check that both 'effects' and 'sigma' are either both defined (not NULL)
    # or are both NULL
    if(is.null(effects) & is.null(sigma)) {
        sigma_defined <- FALSE
    } else {
        sigma_defined <- TRUE
    }
    if((!is.null(effects) & is.null(sigma)) ||
        (is.null(effects) & !is.null(sigma))) {
        stop(
            "Both 'effects' and 'sigma' must be specified.",
            call.=TRUE)
    }
    # Check the random effects specification. If it is NULL (none specified),
    # then we will throw an error. Otherwise, it should be a list of random
    # effects to add.
    if(is.null(r_eff) | !is.list(r_eff)) {
        stop("Random effects must be specified as a list! See help page.",
             call.=TRUE)
    }
    # This is a bit ugly, but hopefully the list of random effects is not so
    # large that this for() loop becomes a problem
    rand_effs <- list()
    for(re in names(r_eff)) {
        # Check that the random effect variable is a single character string.
        if(!is.character(r_eff[[re]]) & length(r_eff[[re]]) != 1) {
            stop(
                paste("Random effect '",
                      re,
                      "' is not a single character string.",
                      sep=""),
                call.=TRUE)
        }
        # Check that the random effects are part of the data to which the
        # model will be fit
        if(!r_eff[[re]] %in% colnames(dat)) {
            stop(paste(
                "Random effect '",
                re,
                "' not found in the dataset!",
                sep=""),
            call.=TRUE)
        }
        # Remove data rows where the random effects columns have missing data
        if(any(is.na(dat[[r_eff[[re]]]]))) {
            miss_rows <- which(is.na(dat[[r_eff[[re]]]]))
            warning(paste(
                "Input data has missing values for random effect '",
                re,
                "' on rows ",
                paste(miss_rows, collapse=", "),
                ". Removing these rows.",
                sep=""),
            immediate.=TRUE,
            call.=TRUE
            )
            dat <- dat[-miss_rows,]
        }
        # Finally if we get here, then make the random effects formulae
        revar <- r_eff[[re]]
        re_f_string <- paste("~0+fit:", r_eff[[re]], sep="")
        rand_effs[[re]] <- as.formula(re_f_string)
    }

    # Build the formula. If the user has specified a custom formula, then use
    # that instead.
    if(is.null(formula)) {
        # The default formula will be a super basic one. It will use just the
        # response (resp), the graphical model variables (varb), the fitness
        # surrogate (fit), and whichever blocking factors were specified.
        if(is.null(blocking)) {
            mod <- resp ~ varb + fit
        } else {
            # We will only support additive "main effects" of blocking
            # factors. Interactions or more complicated linear combinations of
            # blocking factors will need a custom formula.
            block_string <- paste(blocking, collapse="+")
            mod_string <- paste("resp~varb+fit:(", block_string, ")", sep="")
            mod <- as.formula(mod_string)
        }
    } else {
        mod <- formula
    }

    # If the user asked for non-quiet output, print what we will be fitting
    # to the terminal.
    if(!quiet) {
        cat("Random effects:\n", file=stderr())
        cat(capture.output(print(rand_effs)), file=stderr())
        cat("\n\n", file=stderr())
        cat("Formula:\n", file=stderr())
        cat(capture.output(print(mod)), file=stderr())
        cat("\n", file=stderr())
    }

    # Then, fit the model and return the model object
    #   This is a bit of an ugly hack, but we have to write the input data into
    #   the global environment because the reaster() function tries to access
    #   the data from the global environment.
    .GlobalEnv$reaster_input_data__ <- dat
    if(!sigma_defined) {
        rout <- aster::reaster(
            fixed=mod,
            random=rand_effs,
            pred=model_spec$pred,
            fam=model_spec$fam,
            varvar=varb,
            idvar=id,
            root=reaster_input_data__[[model_spec$initial]],
            data=reaster_input_data__)
    } else {
        rout <- aster::reaster(
            fixed=mod,
            random=rand_effs,
            pred=model_spec$pred,
            fam=model_spec$fam,
            varvar=varb,
            idvar=id,
            root=reaster_input_data__[[model_spec$initial]],
            data=reaster_input_data__,
            effects=effects,
            sigma=sigma)
    }
    return(rout)
}

#' Calculate VaW from a random-effects Aster model using a paternal half-sibling
#' design.
VaW_paternal_halfsib <- function(reaster_obj, sire_label, effect_label, 
    model_spec, fixed_eff_idx, typical_ind_idx) {
    # Extract information about the graphical model from the model specification
    # object
    n_nodes <- length(model_spec$vars)
    fit_node <- which(model_spec$vars == model_spec$fit)
    bhat<- reaster_obj$b # random effects
    bhat.sire<- bhat[grep(sire_label, names(bhat))] # specifies Sire effects
    hoom.star <- predict(reaster_obj$obj,
                       newcoef=reaster_obj$alpha)
    hoom.star<- matrix(hoom.star, ncol=n_nodes)
    hoom.star<- hoom.star[,fit_node]
    # mapping function
    map <- function(b) {
        stopifnot(length(b) == 1)
        stopifnot(is.finite(b))
        alpha <- reaster_obj$alpha
        alpha[fixed_eff_idx] <- alpha[fixed_eff_idx] + b 
        # adding random effect to fixed effect
        hoom.star <- predict(reaster_obj$obj, newcoef=alpha)
        hoom.star<- matrix(hoom.star, ncol=n_nodes)
        return(hoom.star[typical_ind_idx, fit_node]) # return value of final node for typical individual
    }
    map.vector <- Vectorize(map)
    bhat.sire.mu<- map.vector(bhat.sire)
    hoom.star2<- predict(
        reaster_obj$obj,
        newcoef=reaster_obj$alpha,
        se.fit=TRUE,
        info.tol=1e-13)
    goom.star <- hoom.star2$gradient
    moom.star<- goom.star[,fixed_eff_idx]
    moom.star<- matrix(moom.star, ncol=n_nodes)
    # calcualtion for Va(w)
    boot_Va<- 4 * moom.star[typical_ind_idx, fit_node]^2 * reaster_obj$nu[effect_label] # final calcuation of VaW
    soutstar <- summary(reaster_obj)
    boot_SE<- 4 * moom.star[typical_ind_idx, fit_node]^2 * soutstar$nu[effect_label, "Std. Error"]
    vaandse <- c(boot_Va,boot_SE)
    return(vaandse)
    }

# TODO:
#   See TR658 for Pearson residuals. We will likely want to include some
#   way to plot these.
