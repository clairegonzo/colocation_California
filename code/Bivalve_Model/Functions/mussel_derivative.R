
mussel_derivative <- function(T, N) {
  # GLOBAL VARIABLES
  # globalenv()$Temps <- Temps
  # globalenv()$mlds <- mlds
  # globalenv()$Vs <- Vs
  # globalenv()$Fs <- Fs
  # globalenv()$Time <- Time
  
  # STATE VARIABLES
  mE <- N[1] # Reserve biomass (molC)
  L_w <- N[2] # Physical length (molC)
  MR <- N[3] # Reproductive biomass (molC)
  F <- N[4] # POC concentration (molC/cm3)
  n <- N[5] # Time-dependent abundance of mussels in the farm modified by mortality rate (dn) in units of (# indiv.)
  
  # ENVIRONMENTAL VARIABLES
  # Interpolate using cubic spline
  #~#~#~##~############ ERROR HERE #################~#~##~#~#
  # why do we even need to reintorpolate these values???? #
  Temp <- pracma::interp1(Time, Temps, T, method = "cubic") #consider using approx() or pchip() instead of interp1()
  V <- pracma::interp1(Time, Vs, T, method = "cubic")
  mld <- pracma::interp1(Time, mlds, T, method = "cubic") 
  F_in <- pracma::interp1(Time, Fs, T, method = "cubic")
  
  # TEMPERATURE DEPENDENT PARAMETERS
  # Conversion function to calculate rate parameters for actual temperature
  Arrh <- exp((Ta / Tref) - (Ta / Temp))
  
  v <- Arrh * v_ref # temperature adjusted energy conductance rate (cm/d)
  v_w <- v / deltam # temperature adjusted physical length conversion rate (cm/d)
  Mvdensity_w <- Mvdensity * deltam^3
  Lp_w <- Lp / deltam # shell length at puberty (approx 1.2 cm)
  kM <- kM_ref * Arrh # temperature adjusted maintenance rate coeff (1/d)
  jEM <- 1 * kM / yVE # 1.4672; kM/yVE [maintenance rate coeff/yield of structure from reserves] (1/d)
  kJ <- maint_ratio * jEM * yVE # = maint_ratio * kM [maintenance ratio * maintenance coefficient]
  
  # FUNCTIONAL RESPONSE
  # Farm dimensions
  farm_area <- 1000 * mld # Cross-sectional area of farm in m2
  
  # Food competition
  r <- V * farm_area * 10000 # V(cm/d) * farm_area(m2) * 10000(cm2/m2) = (cm3/d)
  supply <- F_in * r # F_in(molC/cm3) * r(cm3/d) = (molC/d)
  Fh <- 1.21e-8 # half saturation constant (molC/cm3)
  f <- F / (Fh + F) # scaled functional response
  Jx <- Jxmax * L_w^2 * f # scaled individual rate of food consumption (molC/d)
  consumption <- n * Jx # total food consumption by whole mussel population (molC/d)
  dF <- supply - consumption - r * F
  
  # GROWTH EQUATIONS
  ME <- Mvdensity_w * mE * L_w^3
  if ((kappa * v_w * mE) > (jEM * L_w)) {
    # Then there is growth
    dmE <- yEX * Jx / (Mvdensity_w * L_w^3) - v_w * mE / L_w # reserve dynamics (mol reserve C/mol structural C time)
    dL_w <- (kappa * v_w * mE - jEM * L_w) / (3 * kappa * mE + 3 / yVE) # shell growth
    if (L_w >= Lp_w) {
      Jcr <- kappar * ((1 - kappa) * mE * Mvdensity_w * (v_w * L_w^2 + yVE * jEM * L_w^3) / (kappa * mE * yVE + 1) - kJ * Ehp / mu_E) # gonad production rate (mol C/time)
    } else {
      Jcr <- 0
    }
  } else {
    if (v_w * mE > jEM * L_w) {
      dmE <- yEX * Jx / (Mvdensity_w * L_w^3) - v_w * mE / L_w
      dL_w <- 0
      if (L_w >= Lp_w) {
        Jcr <- kappar * (Mvdensity_w * L_w^3 * (mE * v_w / L_w - jEM) - kJ * Ehp / mu_E)
        if (Jcr < 0) {
          cat("Animal is starving; Gonads are resorbed for maturity maintenance\n")
          cat("Discuss with Roger or Erik if this happens frequently\n")
          if (MR < 0) {
            cat("WARNING: maturity maintenance is not fully paid\n")
            cat("Results not reliable after this point\n")
            Jcr <- 0
          }
        }
      } else {
        Jcr <- 0
        if (Mvdensity_w * L_w^3 * (mE * v_w / L_w - jEM) < kJ * Ehp / mu_E) {
          cat("WARNING: Juvenile may not mature properly\n")
          cat("Results may not be reliable\n")
          cat("Code needs to be changed to include maturation dynamics\n")
          cat("See Roger or Erik\n")
        }
      }
    } else {
      dL_w <- 0
      Jcr <- 0
      if (L_w >= Lp_w) {
        if (MR > 0) {
          dmE <- yEX * Jx / (Mvdensity_w * L_w^3) - v_w * mE / L_w
          Jcr <- Mvdensity_w * L_w^3 * (v_w * mE / L_w - jEM) - kJ * Ehp / mu_E
          cat("Animal is starving; gonad resorption for maintenance\n")
          cat("See Roger or Erik\n")
        } else {
          dmE <- yEX * Jx / (Mvdensity_w * L_w^3) - jEM
          cat("Animal is starving; extra reserve mobilization to pay somatic maintenance\n")
          cat("This is a debatable choice; Standard DEB prescribes death\n")
          cat("See Erik or Roger\n")
          if (ME < 0) {
            dmE <- 0
            cat("Animal starved to death\n")
          }
        }
      } else {
        dmE <- yEX * Jx / (Mvdensity_w * L_w^3) - jEM
        cat("Animal is starving; extra reserve mobilization to pay somatic maintenance\n")
        cat("This is a debatable choice; Standard DEB prescribes death\n")
        cat("See Erik or Roger\n")
        if (ME < 0) {
          dmE <- 0
          cat("Animal starved to death\n")
        }
      }
    }
  }
  dn <- 0 # currently no mortality; Alternative: dn <- -0.0005 * n * (1 - (L_w / 12))
  # calculation: 15% annual mortality rate (FAO); .15/300 = (probability of death/day);
  # logistically size-dependent.
  
  output <- c(dmE, dL_w, Jcr, dF, dn)
  return(output)
}

