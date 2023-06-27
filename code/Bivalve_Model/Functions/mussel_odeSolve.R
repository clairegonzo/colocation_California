mussel_odeSolve <- function(cellno, celldata, days) {
  # Interpolate to get daily values for Temp, Vs, F, mld (so datasets match 'Time')
  year_Oct2Mar <- seq(as.Date("2014-10-01"), as.Date("2016-03-31"), by = "day")
  growPeriod <- length(days) - 1
  dayno <- which(days[1] == year_Oct2Mar)[1]
  
  # # Make environmental parameters global for use in the derivative function
  # globalenv$mlds <- NULL
  # globalenv$Temps <- NULL
  # globalenv$Fs <- NULL
  # globalenv$Vs <- NULL
  # globalenv$Time <- NULL
  
  ## Mixed layer depth in meters(m)
  mld_monthly <- celldata[4:21]
  
  # Create an empty matrix to store the interpolated values
  mld_daily <- matrix(nrow = length(mld_monthly[[1]]), ncol = 547)
  # Convert each column to a numeric vector and apply pchip(), and fill the matrix with the interpolated values
  for (i in 1:nrow(mld_monthly)) {
    row_vector <- as.numeric(mld_monthly[i,])
    interp_vector <- pchip(1:18, row_vector, seq(1,18, by = 0.0311))
    mld_daily[i ,] <- interp_vector
  }
  mlds <- mld_daily[dayno:(dayno + growPeriod)]
  
  ## Temperature in Kelvin
  Temp_monthly <- celldata[22:39] + 273
  #holding matrix
  Temp_daily <- matrix(nrow = length(Temp_monthly[[1]]), ncol = 547)
  # Convert each column to a numeric vector and apply pchip(), and fill the matrix with the interpolated values
  for (i in 1:nrow(Temp_monthly)) {
    row_vector <- as.numeric(Temp_monthly[i,])
    interp_vector <- pchip(1:18, row_vector, seq(1,18, by = 0.0311))
    Temp_daily[i ,] <- interp_vector
  }
  Temps <- Temp_daily[dayno:(dayno + growPeriod)]
  
  ## POC in molC/cm3
  # Calculation: X (mgC/m3) * (1molC/12000mgC) * (1m3/1000000cm3)
  F_monthly <- celldata[58:75] / 12e9
  # holding matrix
  F_daily <- matrix(nrow = length(F_monthly[[1]]), ncol = 547)
  # Convert each column to a numeric vector and apply pchip(), and fill the matrix with the interpolated values
  for (i in 1:nrow(F_monthly)) {
    row_vector <- as.numeric(F_monthly[i,])
    interp_vector <- pchip(1:18, row_vector, seq(1,18, by = 0.0311))
    F_daily[i ,] <- interp_vector
  }
  Fs <- F_daily[dayno:(dayno + growPeriod)]
  
  ## Current speed in cm/d
  # Calculation: Input data = m/s; X (100cm/m) * (86400s/d)
  V_monthly <- celldata[40:57] * 8640000
  # holding matrix
  V_daily <- matrix(nrow = length(V_monthly[[1]]), ncol = 547)
  # Convert each column to a numeric vector and apply pchip(), and fill the matrix with the interpolated values  
  for (i in 1:nrow(V_monthly)) {
    row_vector <- as.numeric(V_monthly[i,])
    interp_vector <- pchip(1:18, row_vector, seq(1,18, by = 0.0311))
    V_daily[i ,] <- interp_vector
  }
  Vs <- V_daily[dayno:(dayno + growPeriod)]
  
  
  # Set time range/length of harvest cycle
  Time <- as.numeric(days)
  
  # Set initial conditions
  # Addmypet parameters from http://www.bio.vu.nl/thb/deb/deblab/add_my_pet/entries_web/i_results_Mytilus_galloprovincialis.html
  
  # physical length 0.3 mm = 0.03cm (approx from Bernard)
  # at birth: 0.00885056 (addmypet)
  L_w_init <- 0.03
  
  # mol reserve C/mol structural C
  # initial reserve mass at growth ceasing at birth = 1.49992e-010 molC (addymypet)
  # initial structural mass at birth (molC) = 8.86389e-011 (addmypet)
  mE_init <- (4.24e-10 / 3.22e-10)
  
  # mol gonadal (C mol)
  MR_init <- 0
  
  # POC (molC/cm3)
  F_init <- Fs[1]
  
  # initial number of mussels in farm (# indiv.)
  # Assumptions:
  # 100 lines per farm
  # 13,000 feet of fuzzy rope per line
  # 100 mussels per foot of fuzzy rope
  n_init <- 100 * 13000 * 100
  
  Init <- c(mE_init, L_w_init, MR_init, F_init, n_init)
  
  # ODE Solver
  options <- c(atol = 1e-7, rtol = 1e-7)
  sol <- pracma::ode23s(f = mussel_derivative, t0 = Time[1], tfinal = Time[length(Time)], y0 = Init, atol = options$atol, rtol = options$rtol)
  
  # Derived Values (weights of product)
  deltam <- 0.1989  # aspect ratio; see notes on size measures
  Mvdensity <- 0.0041841  # density of structure molC/cm3
  Mvdensity_w <- Mvdensity * deltam^3
  LW <- sol$y[, 2]  # apical length of an individual mussel (cm)
  MV <- 12 * Mvdensity_w * LW^3  # structural biomass of an individual mussel (gC)
  ME <- 12 * Mvdensity_w * sol$y[, 1] * LW^3  # reserve biomass of an individual mussel (gC)
  MR <- 12 * sol$y[, 3]  # reproductive biomass of an individual mussel (gC)
  M <- MV + ME + MR  # total biomass of an individual mussel (gC)
  Indiv <- sol$y[, 5]
  C_content <- 0.034  # Carbon weight/wet weight; Haamer, J. 1996.
  # Improving water quality in a eutrophied fjord system with mussel
  # farming. Ambio. Vol. 25. pp. 356-362.
  TM <- Indiv * M / (1000 * C_content)  # total biomass of all mussels on the farm (kg)
  Mwet <- M / C_content  # wet weight of individual mussel (g)
  
  # Available Food
  Fdens <- sol$y[, 4] * 12e9  # in mgC/m3
  # Calculation: X (molC/cm3) / [(1molC/12000mgC) * (1m3/1000000cm3)]
  
  # Functional Response
  Fh <- 1.21e-8
  f <- sol$y[, 4] / (Fh + sol$y[, 4])  # scaled functional response
  
  # Physical Conditions
  Temperature <- approx(x = Time, y = Temps, xout = sol$x, method = "linear")$y - 273  # Temperature in degrees C
  MixedLayerDepth <- approx(x = Time, y = mlds, xout = sol$x, method = "linear")$y  # Mixed layer depth in meters
  F_Ambient <- approx(x = Time, y = Fs, xout = sol$x, method = "linear")$y * 12e9  # Ambient POC concentration in mgC/m3
  CurrentSpeed <- approx(x = Time, y = Vs / 86400, xout = sol$x, method = "linear")$y  # Current speed in cm/s
  
  return(list(T = sol$x, TM = TM, Indiv = Indiv, MV = MV, ME = ME, MR = MR, M = M, Mwet = Mwet, LW = LW, Fdens = Fdens, F_Ambient = F_Ambient, f = f, Temperature = Temperature, MixedLayerDepth = MixedLayerDepth, CurrentSpeed = CurrentSpeed))
}
