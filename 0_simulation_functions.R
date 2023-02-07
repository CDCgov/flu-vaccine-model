## File description and notes --------------------------------------------------
# Sinead E Morris
# Edit of original SEAIRVTModel.R from flumodels to include 2 vaccine classes
# (standard vaccines and higher-dose or adjuvanted vaccines)


## Main Function --------------------------------------------------------------------

SEAIRV_2vaccines <- function(population, populationFractions, contactMatrix, R0,
                          latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                          fractionLatentThatIsInfectious, relativeInfectivityAsymptomatic,
                          useCommunityMitigation, communityMitigationStartDay,
                          communityMitigationDuration, communityMitigationMultiplier,
                          fractionSymptomatic,  fractionSeekCare, 
                          fractionDiagnosedAndPrescribedOutpatient,
                          fractionAdhere, fractionAdmitted, fractionDiagnosedAndPrescribedInpatient, 
                          AVEi, AVEp,
                          vaccineAdministrationRatePerDay, 
                          vaccineAvailabilityByDay,
                          hd_vaccineAvailabilityByDay,
                          vaccineUptakeMultiplier, 
                          hd_vaccineUptakeMultiplier, 
                          VEs, VEi, VEp, 
                          hd_VEs, hd_VEi, hd_VEp,
                          vaccineEfficacyDelay,
                          boostDelay = 0, simulationLength, 
                          seedStartDay, tolerance, method,
                          useSecondCommunityMitigation=F, secondCommunityMitigationDuration, 
                          secondCommunityMitigationMultiplier, 
                          hospDuration, caseHospitalizationRatio, 
                          DelayOnsetToHosp, hospFatalityRatio, 
                          Seasonality=F,offset=-3.8*7,R0_amp=0, dayOfR0=75, 
                          vaccineUptakeMultipliers=list(1), 
                          vaccineUptakeMultiplierStartDay=c(0),
                          hd_vaccineUptakeMultipliers=list(1), 
                          hd_vaccineUptakeMultiplierStartDay=c(0)
                          ) {
  # Check inputs ----
  specifiedArguments <- names(match.call())[-1]
  
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  if(!"contactMatrix" %in% names(argumentList)){stop("contactMatrix must be supplied")}
  
  if(!"vaccineUptakeMultiplier" %in% names(argumentList) )
  {argumentList$vaccineUptakeMultiplier = vaccineUptakeMultipliers[[1]]}
  if(!"vaccineUptakeMultiplierStartDay" %in% names(argumentList) )
  {argumentList$vaccineUptakeMultiplierStartDay = vaccineUptakeMultiplierStartDay}
  
  if(!"hd_vaccineUptakeMultiplier" %in% names(argumentList) )
  {argumentList$hd_vaccineUptakeMultiplier = hd_vaccineUptakeMultipliers[[1]]}
  if(!"hd_vaccineUptakeMultiplierStartDay" %in% names(argumentList) )
  {argumentList$hd_vaccineUptakeMultiplierStartDay = hd_vaccineUptakeMultiplierStartDay}
  
  if(!"boostDelay" %in% names(argumentList) )
  {argumentList$boostDelay=boostDelay}

  parameters <- do.call("checkInputs.SEAIRV.2vaccines", argumentList) #Get parameters from checked inputs
  parameters$useSecondCommunityMitigation = useSecondCommunityMitigation
  parameters$Seasonality = Seasonality
  
  if(! "boostDelay" %in% names(parameters)){stop("boostDelay not in names(parameters)")}
  if(Seasonality)
  {
    if(missing(offset)|missing(dayOfR0)|missing(R0_amp))
    {stop("offset, R0_amp, and dayOfR0 must be specified when transmission seasonality is used")}
    else
    {
      parameters$offset = offset
      parameters$dayOfR0 = dayOfR0
      parameters$R0_amp = abs(R0_amp)
    }  }
  if (useSecondCommunityMitigation) {
    if (missing(secondCommunityMitigationDuration)) {
      stop("secondCommunityMitigationDuration must be specified when using second community mitigation.", 
           call. = FALSE)
    }
    else
    {parameters$secondCommunityMitigationDuration = secondCommunityMitigationDuration  }
    if (missing(secondCommunityMitigationMultiplier)) {
      stop("secondCommunityMitigationMultiplier must be specified when using second community mitigation.", 
           call. = FALSE)
    }
    else
    {checkNonNegative(secondCommunityMitigationMultiplier)
      parameters$secondCommunityMitigationMultiplier = secondCommunityMitigationMultiplier}
    if (!all(dim(secondCommunityMitigationMultiplier) == dim(contactMatrix))) {
      stop("Dimensions of secondCommunityMitigationMultiplier do not match those of contactMatrix", 
           call. = FALSE)
    }
    
  }
  
  initialState <- with(parameters, {
    c(S  = (1 - priorImmunity) * populationFractions,
      E  = 0 * populationFractions,
      A  = 0 * populationFractions,
      I  = 0 * populationFractions,
      R  = priorImmunity * populationFractions,
      Sv = 0 * populationFractions,
      Ev = 0 * populationFractions,
      Av = 0 * populationFractions,
      Iv = 0 * populationFractions,
      Rv = 0 * populationFractions,
      Svb = 0 * populationFractions,
      Evb = 0 * populationFractions,
      Avb = 0 * populationFractions,
      Ivb = 0 * populationFractions,
      Rvb = 0 * populationFractions,
      V  = 0 * populationFractions,
      Vb = 0 * populationFractions,
      vaccinatingPrime = rep(1, length(populationFractions)),
      vaccinatingBoost = rep(1, length(populationFractions)),
      Symp = 0 * populationFractions,
      Hosp = 0 * populationFractions,
      Discharged = 0 * populationFractions,
      Died = 0 * populationFractions)
  })
  
  rootFunction <- function(t, state, parameters) {
    stateList <- reconstructState.SEAIRV(state)
    with(append(stateList, parameters), {
      return(c(ifelse(vaccinatingPrime > 0, populationFractions - V - Vb - tolerance, 1),
               ifelse(vaccinatingBoost > 0, populationFractions - V - Vb - tolerance, 1),
               ifelse(populationFractions - V - Vb - tolerance > tolerance & (vaccinatingBoost + vaccinatingPrime) <= 0, 0, 1)))
    })
  }
  
  eventFunction <- function(t, state, parameters) {
    stateList <- reconstructState.SEAIRV(state)
    with(append(stateList, parameters), {
      state[flumodels:::getLabels("vaccinatingPrime", length(populationFractions))] <-
        ifelse(populationFractions - V - Vb > tolerance, 1, 0)
      state[flumodels:::getLabels("vaccinatingBoost", length(populationFractions))] <-
        ifelse(populationFractions - V - Vb > tolerance, 1, 0)
      return(state)
    })
  }
  
  rawOutput <- flumodels:::integrateModel(initialState = initialState,
                                          parameters = parameters,
                                          derivativeFunction = getDerivative.SEAIRTV,
                                          seedFunction = doSeed.SEAIRV,
                                          rootFunction = rootFunction,
                                          eventFunction = eventFunction)
  
  # Build the SEAIRVModel object to return ----
  model <- list(parameters = parameters,
                rawOutput = rawOutput)
  class(model) <- c("SEAIRTVModel")
  return(model)
}



## Functions to check inputs -------------------------------------------------

checkInputs.SEIR <- flumodels:::checkInputs.SEIR

checkInputs.Antiviral <- flumodels:::checkInputs.Antiviral

checkInputs.SEAIRV.2vaccines <- function(population, populationFractions, contactMatrix, R0,
                                latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                                fractionLatentThatIsInfectious, relativeInfectivityAsymptomatic,
                                useCommunityMitigation, communityMitigationStartDay,
                                communityMitigationDuration, communityMitigationMultiplier,
                                fractionSymptomatic,  fractionSeekCare, fractionDiagnosedAndPrescribedOutpatient,
                                fractionAdhere, fractionAdmitted, fractionDiagnosedAndPrescribedInpatient, AVEi, AVEp,
                                vaccineAdministrationRatePerDay, 
                                vaccineAvailabilityByDay,
                                hd_vaccineAvailabilityByDay,
                                vaccineUptakeMultiplier, 
                                hd_vaccineUptakeMultiplier, 
                                VEs, VEi, VEp, 
                                hd_VEs, hd_VEi, hd_VEp,
                                vaccineEfficacyDelay, boostDelay = 0,
                                simulationLength, seedStartDay, tolerance, method, 
                                useSecondCommunityMitigation, secondCommunityMitigationDuration, 
                                secondCommunityMitigationMultiplier, 
                                hospDuration, caseHospitalizationRatio, DelayOnsetToHosp, 
                                hospFatalityRatio,
                                Seasonality=F,offset=-3.8*7,R0_amp=0,dayOfR0=75,
                                vaccineUptakeMultipliers=list(1), 
                                vaccineUptakeMultiplierStartDay=c(0),
                                hd_vaccineUptakeMultipliers=list(1), 
                                hd_vaccineUptakeMultiplierStartDay=c(0)
                                ) 
{
  
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments

  if(missing(contactMatrix)){stop("contactMatrix must be supplied")}
  
  if(!missing(vaccineAdministrationRatePerDay) & missing(vaccineUptakeMultiplier) & (missing(vaccineUptakeMultipliers) | missing(vaccineUptakeMultiplierStartDay)))
  {stop("if vaccinating in model, and vaccineUptakeMultiplier isn't provided, both vaccineUptakeMultipliers and vaccineUptakeMultiplierStartDay must be")}
  if(missing(vaccineUptakeMultiplier) & !missing(vaccineAdministrationRatePerDay)){argumentList$vaccineUptakeMultiplier = vaccineUptakeMultipliers[[1]]}
  SEIRParameters <- do.call("checkInputs.SEIR", argumentList)
  
  
  if (missing(fractionLatentThatIsInfectious)) {
    stop("fractionLatentThatIsInfectious must be specified.", call. = FALSE)
  }
  flumodels:::checkBetween0and1(fractionLatentThatIsInfectious)
  flumodels:::checkPositive(relativeInfectivityAsymptomatic)
  
  
  
  for(k in 1:length(vaccineUptakeMultipliers))
  {
    flumodels:::checkNonNegative(vaccineUptakeMultipliers[[k]])
    flumodels:::checkDimensionsMatch(vaccineUptakeMultipliers[[k]], populationFractions)
  }
  SEIRParameters$vaccineUptakeMultipliers = vaccineUptakeMultipliers
  SEIRParameters$vaccineUptakeMultiplierStartDay = vaccineUptakeMultiplierStartDay
  
  for(k in 1:length(hd_vaccineUptakeMultipliers))
  {
    flumodels:::checkNonNegative(hd_vaccineUptakeMultipliers[[k]])
    flumodels:::checkDimensionsMatch(hd_vaccineUptakeMultipliers[[k]], populationFractions)
  }
  SEIRParameters$hd_vaccineUptakeMultipliers = hd_vaccineUptakeMultipliers
  SEIRParameters$hd_vaccineUptakeMultiplierStartDay = hd_vaccineUptakeMultiplierStartDay
  
  antiviralParameters <- do.call("checkInputs.Antiviral", argumentList)
  
  #Update arguments passed to checkInputs.Vaccine using SEIRParameters
  argumentList$population <- SEIRParameters$population
  argumentList$populationFractions <- SEIRParameters$populationFractions
  argumentList$seedStartDay <- SEIRParameters$seedStartDay
  argumentList$simulationLength <- SEIRParameters$simulationLength
  vaccineParameters <- do.call("checkInputs.2vaccines", argumentList)
  
  # Update SEIRParameters for beta, lambda1, and lambda2
  SEIRParameters$lambda1 = 1 / ((1-fractionLatentThatIsInfectious) * latentPeriod)
  SEIRParameters$lambda2 = 1 / (fractionLatentThatIsInfectious * latentPeriod)
  SEIRParameters$fractionLatentThatIsInfectious = fractionLatentThatIsInfectious
  SEIRParameters$relativeInfectivityAsymptomatic = relativeInfectivityAsymptomatic
  SEIRParameters$hospDuration = hospDuration
  SEIRParameters$caseHospitalizationRatio = caseHospitalizationRatio
  SEIRParameters$hospFatalityRatio = hospFatalityRatio
  SEIRParameters$DelayOnsetToHosp = DelayOnsetToHosp
  SEIRParameters$boostDelay = boostDelay
  # If there are no community mitigations, eg contact matrix among presymptomatic class is the same as among symptomatic class
  SEIRParameters$beta = R0 / 
    max(Mod(eigen(
      (infectiousPeriod + 1 / SEIRParameters$lambda2) * contactMatrix,
      symmetric = FALSE, only.values = TRUE
    )$values))
  
  SEIRParameters$beta_amp = R0_amp / 
    max(Mod(eigen(
      (infectiousPeriod + 1 / SEIRParameters$lambda2) * contactMatrix,
      symmetric = FALSE, only.values = TRUE
    )$values))
  #Return the parameters
  
  if(! "boostDelay" %in% names(SEIRParameters)){stop("boostDelay not in names(vaccineParameters)")}
  return(c(SEIRParameters, antiviralParameters, vaccineParameters))
}


checkInputs.2vaccines <- function(population, populationFractions, 
                                  seedStartDay, simulationLength,
                                  vaccineAdministrationRatePerDay = 0, 
                                  vaccineAvailabilityByDay = 0,
                                  hd_vaccineAvailabilityByDay = 0,
                                  vaccineUptakeMultiplier = 1,
                                  hd_vaccineUptakeMultiplier = 1,
                                  VEs = 0, VEi = 0, VEp = 0, 
                                  hd_VEs = 0, hd_VEi = 0, hd_VEp = 0, 
                                  boostDelay = 0, vaccineEfficacyDelay = 14, ...) {
  #Validate vaccine parameters
  flumodels:::checkNonNegativeNumber(vaccineAdministrationRatePerDay)
  flumodels:::checkNonNegative(vaccineAvailabilityByDay)
  flumodels:::checkNonNegative(hd_vaccineAvailabilityByDay)
  
  flumodels:::checkNonNegative(vaccineUptakeMultiplier)
  flumodels:::checkDimensionsMatch(vaccineUptakeMultiplier, populationFractions)
  flumodels:::checkNonNegative(hd_vaccineUptakeMultiplier)
  flumodels:::checkDimensionsMatch(hd_vaccineUptakeMultiplier, populationFractions)
  
  flumodels:::checkBetween0and1(VEs)
  flumodels:::checkDimensionsMatch(VEs, populationFractions)
  
  flumodels:::checkBetween0and1(VEi)
  flumodels:::checkDimensionsMatch(VEi, populationFractions)
  
  flumodels:::checkBetween0and1(VEp)
  flumodels:::checkDimensionsMatch(VEp, populationFractions)
  
  flumodels:::checkBetween0and1(hd_VEs)
  flumodels:::checkDimensionsMatch(hd_VEs, populationFractions)
  
  flumodels:::checkBetween0and1(hd_VEi)
  flumodels:::checkDimensionsMatch(hd_VEi, populationFractions)
  
  flumodels:::checkBetween0and1(hd_VEp)
  flumodels:::checkDimensionsMatch(hd_VEp, populationFractions)
  
  flumodels:::checkNonNegativeNumber(vaccineEfficacyDelay)
  flumodels:::checkNonNegativeNumber(boostDelay)
  
  #Compute the daily vaccination rates for priming and boost doses
  totalSimulationLength <- seedStartDay + simulationLength
  vaccinationRateByDay <- rep(0, totalSimulationLength)
  hd_vaccinationRateByDay <- rep(0, totalSimulationLength)
  
  currentVaccineAvailability <- 0
  hd_currentVaccineAvailability <- 0
  
  for (i in 1:totalSimulationLength) {
    if (i <= length(vaccineAvailabilityByDay)){
      currentVaccineAvailability <- currentVaccineAvailability + vaccineAvailabilityByDay[i] 
    }
    if (i <= length(hd_vaccineAvailabilityByDay)){
      hd_currentVaccineAvailability <- hd_currentVaccineAvailability + hd_vaccineAvailabilityByDay[i] 
    }
    
    # prioritize people receiving HD vaccine
    hd_vaccinationRateByDay[i] <- max( min(vaccineAdministrationRatePerDay, hd_currentVaccineAvailability), 
                                       0)
    
    vaccinationRateByDay[i] <- max( min(vaccineAdministrationRatePerDay - hd_vaccinationRateByDay[i], 
                                        currentVaccineAvailability), 
                                    0)
    
    currentVaccineAvailability <- currentVaccineAvailability - vaccinationRateByDay[i]
    hd_currentVaccineAvailability <- hd_currentVaccineAvailability - hd_vaccinationRateByDay[i]
  }
  vaccinationRateByDay <- vaccinationRateByDay / population 
  hd_vaccinationRateByDay <- hd_vaccinationRateByDay / population
  
  #Define vaccination rate functions
  vaccinationRate <- function(t) {
    if ((t < vaccineEfficacyDelay) || (t >= totalSimulationLength + vaccineEfficacyDelay)) {
      return(0)
    } else {
      return(vaccinationRateByDay[floor(t - vaccineEfficacyDelay + 1)])
    }
  }
  
  hd_vaccinationRate <- function(t) {
    if ((t < vaccineEfficacyDelay) || (t >= totalSimulationLength + vaccineEfficacyDelay)) {
      return(0)
    } else {
      return(hd_vaccinationRateByDay[floor(t - vaccineEfficacyDelay + 1)])
    }
  }
  
  #Compute the vaccination rate age multiplier
  vaccinationRateAgeMultiplier <- vaccineUptakeMultiplier * populationFractions
  totalMultiplier <- sum(vaccinationRateAgeMultiplier)
  if (totalMultiplier > 0) {
    vaccinationRateAgeMultiplier <- vaccinationRateAgeMultiplier / totalMultiplier
  } else {
    warning("vaccineUptakeMultiplier prevents vaccination from occurring.", call. = FALSE)
  }
  
  #Compute the vaccination rate age multiplier
  hd_vaccinationRateAgeMultiplier <- hd_vaccineUptakeMultiplier * populationFractions
  totalMultiplier <- sum(hd_vaccinationRateAgeMultiplier)
  if (totalMultiplier > 0) {
    hd_vaccinationRateAgeMultiplier <- hd_vaccinationRateAgeMultiplier / totalMultiplier
  } else {
    warning("hd_vaccineUptakeMultiplier prevents vaccination from occurring.", call. = FALSE)
  }
  
  #Return the parameters
  return(list(vaccinationRate = vaccinationRate,
              hd_vaccinationRate = hd_vaccinationRate,
              vaccinationRateAgeMultiplier = vaccinationRateAgeMultiplier,
              hd_vaccinationRateAgeMultiplier = hd_vaccinationRateAgeMultiplier,
              VEs = VEs, VEi = VEi, VEp = VEp, 
              hd_VEs = hd_VEs, hd_VEi = hd_VEi, hd_VEp = hd_VEp,
              vaccineEfficacyDelay = vaccineEfficacyDelay))
}


## Function with equations -----------------------------------------------------

# parameters should define populationFractions, contactMatrix, beta, lambda, gamma,
# AVEi.eff, VEs, VEi, and the function vaccinationRate(t)

# Note: the total population is normalized to be 1

getDerivative.SEAIRTV <- function(t, state, parameters) {
  
  stateList <- reconstructState.SEAIRV(state)
  with(append(stateList, parameters), {
    
    if(Seasonality)
    {
      day = t + 1
      if(length(day)!=1){stop(paste0("length day != 1, length day = ",length(day)))}
      if(length(beta)!=1){stop(paste0("length beta != 1, length beta = ",length(beta)))}
      if(length(R0_amp)!=1){stop(paste0("length R0_amp != 1, length R0_amp = ",length(R0_amp)))}
      if(length(offset)!=1){stop(paste0("length offset != 1, length offset = ",length(offset)))}
      
      beta_Max =  beta + beta_amp*0.5 - 0.5*beta_amp*cos((2*pi/(52*7))*(dayOfR0-(offset)))
      beta_time = 0.5*beta_amp*cos((2*pi/(52*7))*(day-(offset))) + (beta_amp/2)+beta_Max-beta_amp
      if(length(beta_time)!=1){stop(paste0("length beta_time != 1, length beta_time = ",length(beta_time)))}
      beta = beta_time
    }
    
    
    if (useCommunityMitigation) {
      if ((t >= communityMitigationStartDay) && (t < communityMitigationEndDay)) {
        contactMatrix <- communityMitigationMultiplier * contactMatrix
      } 
    }
    if(useSecondCommunityMitigation)
    {
      if((t>communityMitigationEndDay) & (t <= communityMitigationEndDay+secondCommunityMitigationDuration))
      {
        contactMatrix <- secondCommunityMitigationMultiplier * contactMatrix
      }
    }
    
    # replacing negative counts with 0
    S = (abs(S)+S)/2 
    Sv = (abs(Sv)+Sv)/2 
    Svb = (abs(Svb)+Svb)/2 
    
    E = (abs(E)+E)/2 
    Ev = (abs(Ev)+Ev)/2 
    Evb = (abs(Evb)+Evb)/2 
    
    A = (abs(A)+A)/2 
    Av = (abs(Av)+Av)/2 
    Avb = (abs(Avb)+Avb)/2 
    
    I = (abs(I)+I)/2 
    Iv = (abs(Iv)+Iv)/2 
    Ivb = (abs(Ivb)+Ivb)/2 
    
    R = (abs(R)+R)/2 
    Rv = (abs(Rv)+Rv)/2 
    Rvb = (abs(Rvb)+Rvb)/2 
    
    V = (abs(V)+V)/2 
    Vb = (abs(Vb)+Vb)/2 
    

    EligibleForHighDose = ifelse(populationFractions - V - Vb < 0, 0, populationFractions - V - Vb)
    EligibleForStDose = ifelse(populationFractions - V - Vb < 0, 0, populationFractions - V - Vb)
    
    # Standard vaccine -----
    vaxMultiplier = 
      vaccineUptakeMultipliers[[max( which(vaccineUptakeMultiplierStartDay <= 
                                           max(t- vaccineEfficacyDelay + 1,0)) )]] * populationFractions
    
    if (sum(vaxMultiplier)>0) { vaxMultiplier = vaxMultiplier/sum(vaxMultiplier) }
    else{ warning("vaccineUptakeMultipliers prevents vaccination from occurring.", call. = FALSE) }
    
    if (any(vaxMultiplier<0)) { stop(paste("any(vaxMultiplier<0), t =",t)) }
    
    isVaccinatingByAge <- (vaccinatingPrime > 0) & (populationFractions - V - Vb > 0)
    effectiveVaccinationMultiplier <- sum(ifelse(isVaccinatingByAge, 1, 0) * vaxMultiplier)
    
    if (effectiveVaccinationMultiplier > 0) {
      vaccinationRateByAge <- vaccinationRate(t) * vaxMultiplier / effectiveVaccinationMultiplier
    } else { vaccinationRateByAge <- 0 }
  
    # High dose vaccine -----
    hd_vaxMultiplier = 
      hd_vaccineUptakeMultipliers[[max( which(hd_vaccineUptakeMultiplierStartDay <= 
                                              max(t- vaccineEfficacyDelay + 1,0)) )]] * populationFractions
    
    if (sum(hd_vaxMultiplier)>0) { hd_vaxMultiplier = hd_vaxMultiplier/sum(hd_vaxMultiplier) }
    else{ warning("hd_vaccineUptakeMultipliers prevents vaccination from occurring.", call. = FALSE) }
    
    if (any(hd_vaxMultiplier<0)) { stop(paste("any(hd_vaxMultiplier<0), t =",t)) }
    
    hd_isVaccinatingByAge <- (vaccinatingBoost > 0) & (populationFractions - V - Vb > 0)
    hd_effectiveVaccinationMultiplier <- sum(ifelse(hd_isVaccinatingByAge, 1, 0) * hd_vaxMultiplier )
    
    if (hd_effectiveVaccinationMultiplier > 0) {
      hd_vaccinationRateByAge <- hd_vaccinationRate(t) * hd_vaxMultiplier / hd_effectiveVaccinationMultiplier
    } else {hd_vaccinationRateByAge <- 0}

    
    # Flows -----
    # Adjusted to account for VEp, which reduces the impact of AVEi since it reduces fractionSymptomatic
    forceOfInfection <- beta / populationFractions / 
      (fractionSymptomatic + relativeInfectivityAsymptomatic - fractionSymptomatic * relativeInfectivityAsymptomatic) * 
      as.vector((contactMatrix %*% (
        (A + I) * (1 - fractionSymptomatic) * relativeInfectivityAsymptomatic + 
          fractionSymptomatic * (A + (Av * (1 - VEi) * (1 - VEp)) + (Avb * (1 - hd_VEi) * (1 - hd_VEp)) ) + 
          (1 - AVEi.eff) * fractionSymptomatic * (I + (Iv * (1 - VEi) * (1 - VEp)) + (Ivb * (1 - hd_VEi) * (1 - hd_VEp)) ) + 
          (Av + Iv) * relativeInfectivityAsymptomatic * (1 - VEi) * (1 - fractionSymptomatic + fractionSymptomatic * VEp) +
          (Avb + Ivb) * relativeInfectivityAsymptomatic * (1 - hd_VEi) * (1 - fractionSymptomatic + fractionSymptomatic * hd_VEp)
      )))
    
    
    if(any(forceOfInfection<0)){stop("forceOfInfection < 0 , t = ", t)}
    forceOfInfection=ifelse(forceOfInfection<0,0,forceOfInfection)
    
    S_to_E <- S * forceOfInfection
    E_to_A <- lambda1 * E
    A_to_I <- lambda2 * A
    I_to_R <- gamma * I
    
    S_to_Sv <-  ifelse(isVaccinatingByAge, vaccinationRateByAge * S / EligibleForStDose, 0)
    Sv_to_Ev <- Sv * (1 - VEs) * forceOfInfection
    Ev_to_Av <- lambda1 * Ev
    Av_to_Iv <- lambda2 * Av
    Iv_to_Rv <- gamma * Iv
    
    S_to_Svb <- ifelse(hd_isVaccinatingByAge, hd_vaccinationRateByAge *  S / EligibleForHighDose, 0)
    Svb_to_Evb <- Svb * (1 - hd_VEs) * forceOfInfection
    Evb_to_Avb <- lambda1 * Evb
    Avb_to_Ivb <- lambda2 * Avb
    Ivb_to_Rvb <- gamma * Ivb

    
    if(any(S_to_Sv<0)){stop("any(S_to_Sv<0)")}
    if(any(Sv_to_Ev<0)){stop("any(Sv_to_Ev<0)")}
    if(any(Ev_to_Av<0)){stop("any(Ev_to_Av<0)")}
    if(any(Av_to_Iv<0)){stop("any(Av_to_Iv<0)")}
    if(any(Iv_to_Rv<0)){stop("any(Iv_to_Rv<0)")}
    
    if(any(S_to_Svb<0)){stop("any(S_to_Svb<0)")}
    if(any(Svb_to_Evb<0)){stop("any(Svb_to_Evb<0)")}
    if(any(Evb_to_Avb<0)){stop("any(Evb_to_Avb<0)")}
    if(any(Avb_to_Ivb<0)){stop("any(Avb_to_Ivb<0)")}
    if(any(Ivb_to_Rvb<0)){stop("any(Ivb_to_Rvb<0)")}
    
    
    # Derivatives -----
    #Non-vaccinated compartments
    dS <- -S_to_E - S_to_Sv - S_to_Svb
    dE <- S_to_E  - E_to_A
    dA <- E_to_A - A_to_I
    dI <- A_to_I - I_to_R
    dR <- I_to_R
    #Vaccinated compartments
    dSv <- S_to_Sv - Sv_to_Ev 
    dEv <- Sv_to_Ev - Ev_to_Av   
    dAv <- Ev_to_Av - Av_to_Iv
    dIv <- Av_to_Iv - Iv_to_Rv
    dRv <- Iv_to_Rv
    
    dSvb <- S_to_Svb - Svb_to_Evb
    dEvb <- Svb_to_Evb - Evb_to_Avb 
    dAvb <- Evb_to_Avb - Avb_to_Ivb
    dIvb <- Avb_to_Ivb - Ivb_to_Rvb
    dRvb <- Ivb_to_Rvb
    #Auxiliary vaccinated compartment
    dV <- ifelse(isVaccinatingByAge, vaccinationRateByAge, 0)
    #Auxiliary vaccinated compartment - people with a boosting dose
    dVb <- ifelse(hd_isVaccinatingByAge, hd_vaccinationRateByAge, 0)
    
    if (any(S==0 & dV>0) ) { stop("any(S==0 & dV>0)")  }
    if (any(S==0 & dVb>0)) { stop("any(S==0 & dVb>0)") }
    
    # Compartments tracking overall burden
    A_and_Av_to_Symp = lambda2 * (A + Av * (1 - VEp) + Avb * (1 - hd_VEp))  * fractionSymptomatic * caseHospitalizationRatio
    
    Symp_to_Hosp = (1/DelayOnsetToHosp) * Symp 
    Hosp_Discharge = Hosp/hospDuration
    
    dSymp = A_and_Av_to_Symp - Symp_to_Hosp
    dHosp = Symp_to_Hosp - Hosp_Discharge
    dDischarged = Hosp_Discharge
    dDied = Hosp_Discharge * hospFatalityRatio
    zeroVec <- 0 * populationFractions
    
    # Return derivative -------
    return(list(c(dS, dE, dA, dI, dR, dSv, dEv, dAv, dIv, dRv, dSvb, dEvb, dAvb, dIvb, dRvb, 
                  dV, dVb, zeroVec, zeroVec, dSymp, dHosp, dDischarged, dDied )))
  })
}

## Utility function that reconstructs the model state as a list so that equations can refer to compartments by name -----
reconstructState.SEAIRV <- function(state) {
  
  numberOfClasses <- length(state) / 23 #Each of the 11 classes are vectors of the same length
  S  <- state[                       1 :     numberOfClasses ]
  E  <- state[    (numberOfClasses + 1):(2 * numberOfClasses)]
  A  <- state[(2 * numberOfClasses + 1):(3 * numberOfClasses)]
  I  <- state[(3 * numberOfClasses + 1):(4 * numberOfClasses)]
  R  <- state[(4 * numberOfClasses + 1):(5 * numberOfClasses)]
  Sv <- state[(5 * numberOfClasses + 1):(6 * numberOfClasses)]
  Ev <- state[(6 * numberOfClasses + 1):(7 * numberOfClasses)]
  Av <- state[(7 * numberOfClasses + 1):(8 * numberOfClasses)]
  Iv <- state[(8 * numberOfClasses + 1):(9 * numberOfClasses)]
  Rv <- state[(9 * numberOfClasses + 1):(10 * numberOfClasses)]
  Svb <- state[(10 * numberOfClasses + 1):(11 * numberOfClasses)]
  Evb <- state[(11 * numberOfClasses + 1):(12 * numberOfClasses)]
  Avb <- state[(12 * numberOfClasses + 1):(13 * numberOfClasses)]
  Ivb <- state[(13 * numberOfClasses + 1):(14 * numberOfClasses)]
  Rvb <- state[(14 * numberOfClasses + 1):(15 * numberOfClasses)]
  V  <- state[(15 * numberOfClasses + 1):(16 * numberOfClasses)]
  Vb <- state[(16* numberOfClasses + 1):(17 * numberOfClasses)]
  vaccinatingPrime <- state[(17 * numberOfClasses + 1):(18 * numberOfClasses)]
  vaccinatingBoost <- state[(18 * numberOfClasses + 1):(19 * numberOfClasses)]
  Symp = state[(19 * numberOfClasses + 1):(20 * numberOfClasses)]
  Hosp <- state[(20 * numberOfClasses + 1):(21 * numberOfClasses)]
  Discharged = state[(21 * numberOfClasses + 1):(22 * numberOfClasses)]
  Died = state[(22 * numberOfClasses + 1):(23 * numberOfClasses)]
  return(as.list(environment()))
}

## Function that seeds infections in the SEAIR+V model -------------------------
   # parameters should define seedInfections, lambda1, lambda2, and gamma

doSeed.SEAIRV <- function(state, parameters) {
  
  stateList <- reconstructState.SEAIRV(state)
  with(append(stateList, parameters), {
    seedInfectionsFractions <- seedInfections / population
    S <- S - seedInfectionsFractions
    # E, A, I initialized such that A and I have derivatives of 0 at seed
    E <- E + seedInfectionsFractions * (gamma * lambda2) / 
      (lambda1 * lambda2 + gamma * (lambda1 + lambda2))
    A <- A + 0 * (gamma * lambda1) /
      (lambda1 * lambda2 + gamma * (lambda1 + lambda2))
    I <- I + 0 * (lambda1 * lambda2) /
      (lambda1 * lambda2 + gamma * (lambda1 + lambda2))
    
    #Return derivative
    return(c(S, E, A, I, R, Sv, Ev, Av, Iv, Rv,  Svb, Evb, Avb, Ivb, Rvb, V, Vb, vaccinatingPrime, vaccinatingBoost, Symp, Hosp, Discharged, Died))
  })
}


