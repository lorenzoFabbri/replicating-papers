keil2014 <- function(dat, idxs, mc_sample_size = 137000, type_intervention) {
  # Step 1: Model conditional probabilities in observed data
  ## Model A: platelet
  modA <- glm(
    formula = platnorm ~ all + cmv + male + age + agecurs1 +
      agecurs2 + gvhdm1 + daysgvhd + daysnorelapse + wait,
    data = subset(dat, platnormm1 == 0),
    family = "binomial"
  )
  ## Model B: relapse
  modB <- glm(
    formula = relapse ~ all + cmv + male + age + gvhdm1 +
      daysgvhd + platnormm1 + daysnoplatnorm + agecurs1 + agecurs2 +
      day + daysq + wait,
    data = subset(dat, relapsem1 == 0),
    family = "binomial"
  )
  ## Model D: censoring
  modD <- glm(
    formula = censlost ~ all + cmv + male + age + daysgvhd +
      daysnoplatnorm + daysnorelapse + agesq + day +
      daycurs1 + daycurs2,
    data = dat,
    family = "binomial"
  )
  ## Model E: death
  modE <- glm(
    formula = d ~ all + cmv + male + age + gvhd + platnorm +
      daysnoplatnorm + relapse + daysnorelapse + agesq +
      day + daycurs1 + daycurs2 + wait + day*gvhd + daycurs1*gvhd +
      daycurs2*gvhd,
    data = dat,
    family = "binomial"
  )
  ## Model H: exposure
  modH <- glm(
    formula = gvhd  ~ all + cmv + male + age + platnormm1 +
      daysnoplatnorm + relapsem1 + daysnorelapse +
      agecurs1 + agecurs2,
    data = subset(dat, gvhdm1 == 0),
    family = "binomial"
  )
  mods <- list(
    modA = modA,
    modB = modB,
    modD = modD,
    modE = modE,
    modH = modH
  )
  
  # Step 2: Generate time-varying exposures, covariates,
  #         and outcomes in Monte Carlo sample
  dat_baseline <- dat |>
    dplyr::filter(day == 1)
  dat_mc <- dat_baseline[sample(
    nrow(dat_baseline),
    size = mc_sample_size,
    replace = TRUE
  ), ]
  covariates <- c(
    "relapse", "gvhd", "platnorm",
    "relapsem1", "gvhdm1", "platnormm1",
    "daysnorelapse", "daysnogvhd", "daysnoplatnorm",
    "daysrelapse", "daysgvhd", "daysplatnorm"
  )
  
  # Loop over pseudo-subjects
  ret <- lapply(1:nrow(dat_baseline), function(idx) {
    # Simulate the "future" of one subject
    subject <- dat_mc[idx, ]
    simulate(subject, covariates, mods, type_intervention)
  }) # End loop over pseudo-subjects
  
  # Step 3: Estimate effect measure
  mod_survival <- survival::coxph(
    formula = survival::Surv(td * d) ~ gvhd,
    data = res
  )
  
  return(mod_survival)
}
################################################################################

simulate <- function(info, covariates, mods, type_intervention) {
  done <- 0
  
  # Simulate future from baseline covariates
  for (t in 1:1825) {
    # Splines for time t
    info[t, "daysq"] <- t**2
    info[t, "daycu"] <- t**3
    info[t, "daycurs1"] <- splines_t1(t)
    info[t, "daycurs2"] <- splines_t2(t)
    
    # Set baseline and time-varying covariates
    if (t == 1) {
      info <- info |>
        dplyr::mutate(dplyr::across(dplyr::all_of(covariates), ~ 0L))
      ##########################################################################
    } else {
      # Make duplicate of previous time point so that we have time-fixed covariates
      info[t, ] <- info[t-1, ]
      
      if (info[t-1, "relapse"] == 0) {
        info[t, "daysnorelapse"] <- as.integer(info[t-1, "daysnorelapse"] + 1)
      } else {
        info[t, "daysrelapse"] <- as.integer(info[t-1, "daysrelapse"] + 1)
      }
      if (info[t-1, "platnorm"] == 0) {
        info[t, "daysnoplatnorm"] <- as.integer(info[t-1, "daysnoplatnorm"] + 1)
      } else {
        info[t, "daysplatnorm"] <- as.integer(info[t-1, "daysplatnorm"] + 1)
      }
      if (info[t-1, "gvhd"] == 0) {
        info[t, "daysnogvhd"] <- as.integer(info[t-1, "daysnogvhd"] + 1)
      } else {
        info[t, "daysgvhd"] <- as.integer(info[t-1, "daysgvhd"] + 1)
      }
      
      if (info[t-1, "platnorm"] == 1) {
        info[t, "platnorm"] <- 1L
      } else {
        info[t, "platnorm"] <- as.integer(rbinom(
          n = 1,
          size = 1,
          prob = predict(mods$modA, newdata = info[t-1, ], type = "response")
        ))
      }
      
      if (info[t-1, "relapse"] == 1) {
        info[t, "relapse"] <- 1L
      } else {
        info_relapse_0 <- info[t-1, ]
        info_relapse_0$platnorm <- as.integer(info[t, "platnorm"])
        info[t, "relapse"] <- as.integer(rbinom(
          n = 1,
          size = 1,
          prob = predict(mods$modB, newdata = info_relapse_0, type = "response")
        ))
      }
      
      if (info[t-1, "gvhd"] == 1) {
        info[t, "gvhd"] <- 1L
      } else {
        info_gvhd_0 <- info[t-1, ]
        info_gvhd_0$platnorm <- as.integer(info[t, "platnorm"])
        info_gvhd_0$relapse <- as.integer(info[t, "relapse"])
        if (type_intervention == "natural_course") {
          info[t, "gvhd"] <- as.integer(rbinom(
            n = 1,
            size = 1,
            prob = predict(mods$modH, newdata = info_gvhd_0, type = "response")
          ))
        } else if (type_intervention == "always") {
          info[t, "gvhd"] <- 1L
        } else if (type_intervention == "never") {
          info[t, "gvhd"] <- 0L
        }
      }
      
      if (done == 0) {
        info_censlost_0 <- info[t-1, ]
        info_censlost_0$platnorm <- as.integer(info[t, "platnorm"])
        info_censlost_0$relapse <- as.integer(info[t, "relapse"])
        info_censlost_0$gvhd <- as.integer(info[t, "gvhd"])
        info[t, "censlost"] <- as.integer(rbinom(
          n = 1,
          size = 1,
          prob = predict(mods$modD, newdata = info_censlost_0, type = "response")
        ))
        if (type_intervention != "natural_course") {
          info[t, "censlost"] <- 0L
        }
        done <- as.integer(info[t, "censlost"])
        
        if (done == 0) {
          info_d_0 <- info[t-1, ]
          info_d_0$platnorm <- as.integer(info[t, "platnorm"])
          info_d_0$relapse <- as.integer(info[t, "relapse"])
          info_d_0$gvhd <- as.integer(info[t, "gvhd"])
          info_d_0$censlost <- as.integer(info[t, "censlost"])
          info[t, "d"] <- as.integer(rbinom(
            n = 1,
            size = 1,
            prob = predict(mods$modE, newdata = info_d_0, type = "response")
          ))
          done <- as.integer(info[t, "d"])
        } # End done death
        
        if (t >= 1825) { done <- 1L }
        if (info[t, "gvhd"] == 1 & info[t-1, "gvhd"] == 0) { info[t, "tg"] <- t }
        if (info[t, "relapse"] == 1 & info[t-1, "relapse"] == 0) { info[t, "tr"] <- t }
        if (info[t, "platnorm"] == 1 & info[t-1, "platnorm"] == 0) { info[t, "tp"] <- t }
        
        if (done == 1) {
          info[t, "td"] <- t
          if (info[t, "gvhd"] == 0) { info[t, "tg"] <- t+1 }
          if (info[t, "relapse"] == 0) { info[t, "tr"] <- t+1 }
          if (info[t, "platnorm"] == 0) { info[t, "tp"] <- t+1 }
          
          return(info)
        }
      } # End done censoring
    } # End set covariates
    
    info[t, "relapsem1"] <- as.integer(info[t, "relapse"])
    info[t, "gvhdm1"] <- as.integer(info[t, "gvhd"])
    info[t, "platnormm1"] <- as.integer(info[t, "platnorm"])
  } # End loop over time
  
  return(info)
}
################################################################################

splines_t1 <- function(t) {
  ((t > 83.6) * ((t - 83.6) / 83.6)**3) +
    ((t > 1862.2) * ((t - 1862.2) / 83.6)**3) * (947.0 - 83.6) -
    ((t > 947.0) * ((t - 947.0) / 83.6)**3) * (1862.2 - 83.6) / (1862.2 - 947.0)
}
################################################################################

splines_t2 <- function(t) {
  ((t > 401.4) * ((t - 401.4) / 83.6)**3) +
    ((t > 1862.2) * ((t - 1862.2) / 83.6)**3) * (947.0 - 401.4) -
    ((t > 947.0) * ((t - 947.0) / 83.6)**3) * (1862.2 - 401.4) / (1862.2 - 947.0)
}
################################################################################
