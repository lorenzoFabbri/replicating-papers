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
  ret <- lapply(seq_along(nrow(dat_baseline)), function(idx) {
    subject <- dat_mc |>
      dplyr::filter(id == idx)
    simulate(subject, covariates, mods, type_intervention)
  }) # End loop over pseudo-subjects
  
  # Step 3: Estimate effect measure
}
################################################################################

simulate <- function(info, covariates, mods, type_intervention) {
  for (t in 1:nrow(info)) {
    # Cubic splines for time t
    t_sq <- t**2
    t_c <- t**3
    daycurs1 <- splines_t1(t)
    daycurs2 <- splines_t2(t)
    
    # Set baseline and time-varying covariates
    if (t == 1) {
      info_0 <- info[t, ] |>
        dplyr::mutate(dplyr::across(dplyr::all_of(covariates), ~ 0))
      
      info[t, "platnorm"] <- predict(
        mods$modA, newdata = info_0, type = "response"
      )
      info[t, "relapse"] <- predict(
        mods$modB, newdata = info_0, type = "response"
      )
      info[t, "censlost"] <- predict(
        mods$modD, newdata = info_0, type = "response"
      )
      info[t, "d"] <- predict(
        mods$modE, newdata = info_0, type = "response"
      )
      
      if (type_intervention == "natural_course") {
        info[t, "gvhd"] <- predict(
          mods$modH, newdata = info_0, type = "response"
        )
      } else if (type_intervention == "always") {
        # TODO
      } else if (type_intervention == "never") {
        # TODO
      }
      
    } else {
      info_platnorm_0 <- info[t-1, ]
      info_relapse_0 <- info[t-1, ]
      info_relapse_0$platnorm <- info[t, ]$platnorm
      info_censlost_0 <- info[t, ]
      info_censlost_0$d <- info[t-1, ]$d
      info_censlost_0$censlost <- info[t-1, ]$censlost
      info_d_0 <- info[t, ]
      info_d_0$d <- info[t-1, ]$d
      info_gvhd_0 <- info[t-1, ]
      info_gvhd_0$platnorm <- info[t, ]$platnorm
      info_gvhd_0$relapse <- info[t, ]$relapse
    }
    
    relapsem1
    gvhdm1
    platnormm1
    # End set covariates
  } # End loop over time
}
################################################################################

splines_t1 <- function(t) {
  ((day > 83.6) * ((day - 83.6) / 83.6)**3) +
    ((day > 1862.2) * ((day - 1862.2) / 83.6)**3) * (947.0 - 83.6) -
    ((day > 947.0) * ((day - 947.0) / 83.6)**3) * (1862.2 - 83.6) / (1862.2 - 947.0)
}
################################################################################

splines_t2 <- function(t) {
  ((day > 401.4) * ((day - 401.4) / 83.6)**3) +
    ((day > 1862.2) * ((day - 1862.2) / 83.6)**3) * (947.0 - 401.4) -
    ((day > 947.0) * ((day - 947.0) / 83.6)**3) * (1862.2 - 401.4) / (1862.2 - 947.0)
}
################################################################################
