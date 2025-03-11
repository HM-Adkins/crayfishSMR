#-------------------------------

# crayfishSMR

#-------------------------------

# Author: Hannah Adkins

#-----------------------------------------------------------------------------

# References
# Chabot, D., McKenzie, D. J., & Craig, J. F. (2016). Metabolic rate in fishes: Definitions, methods and significance for conservation physiology: editorial. Journal of Fish Biology, 88(1), 1–9. https://doi.org/10.1111/jfb.12873
# Chabot, D., Steffensen, J. F., & Farrell, A. P. (2016). The determination of standard metabolic rate in fishes: Measuring smr in fishes. Journal of Fish Biology, 88(1), 81–121. https://doi.org/10.1111/jfb.12845
# https://github.com/denis-chabot/fishMO2.git

#-----------------------------------------------------------------------------

# Load libraries
library(forcats)
library(dplyr)
library(devtools)
install_github("denis-chabot/fishMO2")
library(fishMO2)

# Citation - fishMO2
bibentry(
  bibtype = "Manual",
  title = "{fishMO2}: Calculate and plot the standard metabolic rate (SMR), the critical oxygen level (O2crit) and the specific dynamic action (SDA) and related variables in fishes and crustaceans",
  author = c(person("Denis", "Chabot")),
  year = 2020,
  url = "https://www.researchgate.net/project/fishMO2-a-R-package-to-calculate-and-plot-SMR-O2crit-and-SDA"
)


#-----------------------------------------------------------------------------

# SMR calculation
# Summarized definitions of output values from fishMO2 package
# mlnd - The Mean of the Lowest Normal Distribution, calculated by ⁠Mclust⁠.
# quant	- A vector of quantile values, with length = number of values in parameter ⁠q⁠.
# cl - A vector of same length as Z containing the number of the distribution assigned to each oxygen uptake value by ⁠Mclust.
# CVmlnd - The coefficient of variation (in % of MLND) of the oxygen uptake values assigned to the lowest normal distribution.
# Chabot et al. (2016) use CVmlnd to assess whether the MLND method should be used to estimate SMR or not.
# When it cannot (CV > 5.4%), they recommend using the quantile method with p = 0.2.

#----

# Function to calculate SMR (by multiple methods) for each ID
# n - number of individuals
# x - empty list to store results of function
# y - list of data frames where each data frame has a column labeled "MO2" containing the MO2 values for a single individual to be used in calculation of SMR
# Names of the data frames in 'y' should be the ID

calcSMR_byID <- function(n, x, y) {
  for (i in 1:n) {
    x[[i]] <- calcSMR(y[[i]]$MO2, q = c(0.2, 0.25))
  }
  names(x) <- names(y)
  x
}

# To save your results in 'x'
x <- calcSMR_byID()

#----

# Function to store SMR results in data frame and calculate the below values
# Sample size of the LND
# Sample size/mean/SD/CV of MO2
# q10/q30
# Difference between the MLND/q20 and q20/q25
# Percent difference between the MLND/q20 and q20/q25
# n - same as above
# x - same as above
# y - same as above
# z - data frame with an 'ID' column
# a - time over which SMR was calculated (i.e., 24 hr, day or night)

calcSMR_byID_table <- function(n, z, x, y, a) {
  for (i in 1:n) {
    z$time[i] <- a
    z$N_MO2[i] <- length(y[[i]]$MO2)
    z$mean_MO2[i] <- mean(y[[i]]$MO2)
    z$sd_MO2[i] <- sd(y[[i]]$MO2)
    z$N_dist[i] <- max(x[[i]]$cl)
    z$theSMRdistr[i] <- x[[i]]$theSMRdistr
    z$N_theSMRdistr[i] <- length(which(x[[i]]$cl == x[[i]]$theSMRdistr))
    z$MLND[i] <- x[[i]]$mlnd
    z$CVmlnd[i] <- x[[i]]$CVmlnd
    z$q10[i] <- quantile(y[[i]]$MO2, 0.1)
    z$q20[i] <- x[[i]]$quant[1]
    z$q25[i] <- x[[i]]$quant[2]
    z$q30[i] <- quantile(y[[i]]$MO2, 0.3)
  }
  z$CVmo2 <- (z$sd / z$mean) * 100
  z$diff_MLND_q20 <- z$MLND - z$q20
  z$perc_diff_MLND_q20 <- ((z$MLND - z$q20) / ((z$MLND + z$q20) / 2)) * 100
  z$diff_quartile <- z$q25 - z$q20
  z$perc_diff_quartile <- ((z$q25 - z$q20) / ((z$q25 + z$q20) / 2)) * 100
  z <- relocate(z, CVmo2, .before = N_dist)
  z
}

# To save your results in 'z'
z <- calcSMR_byID_table()

#----

# Function to create plots of frequency distribution with fitted normal distributions and lines at estimated MLND/q20/q25 for each ID
# n - same as above
# x - same as above
# y - same as above

plotMO2fdis_byID <- function(n, x, y) {
  par(mfrow = c(3, 2))
  for (i in 1:n) {
    plotMO2fdis(y[[i]]$MO2)
    abline(
      v = x[[i]]$mlnd,
      lty = 1,
      col = "black",
      lwd = 2
    )
    abline(
      v = x[[i]]$quant[1],
      lty = 2,
      col = "red",
      lwd = 2
    )
    abline(
      v = x[[i]]$quant[2],
      lty = 3,
      col = "green3",
      lwd = 2
    )
    legend(
      "topright",
      c("MLND", "q0.2", "q0.25"),
      title = IDs[i],
      lty = c(1, 2, 3),
      col = c("black", "red", "green3"),
      lwd = 2,
      bty = "n",
      y.intersp = 1.2
    )
  }
  par(def.par)
}

#----

# Function to create plots of MO2 over time with lines at estimated MLND/q20/q25 for each ID
# n - same as above
# x - same as above
# y - same as above

plotMO2_byID <- function(n, x, y) {
  par(mfrow = c(3, 2))
  for (i in 1:n) {
    plot(
      y[[i]]$MO2 ~ y[[i]]$Time_hr,
      main = IDs[i],
      ylim = c(0, 250),
      bty = "l",
      pch = x[[i]]$cl,
      ylab = yl,
      xlab = "Time (hr)"
    )
    abline(
      h = x[[i]]$mlnd,
      lty = 1,
      col = "black",
      lwd = 2
    )
    abline(
      h = x[[i]]$quant[1],
      lty = 2,
      col = "red",
      lwd = 2
    )
    abline(
      h = x[[i]]$quant[2],
      lty = 3,
      col = "green3",
      lwd = 2
    )
    legend(
      "topright",
      c("1", "2", "3"),
      col = c("black"),
      pch = c(1, 2, 3),
      bty = "n"
    )
    legend(
      "topright",
      c("MLND", "q0.2", "q0.25"),
      lty = c(1, 2, 3),
      col = c("black", "red", "green3"),
      lwd = 2,
      bty = "n",
      y.intersp = 1.2,
      inset = c(0.075, 0, 0, 0)
    )
  }
  par(def.par)
}

#-----------------------------------------------------------------------------

# Add SMR and SMR method columns to table 'z'
# Use SMR value estimated by MLND if CVmlnd is below the threshold value
# A different method (i.e., q20 or q25) is needed to estimate SMR when the CVmlnd is above the threshold value
# When the LND contains MO2 values likely associated with mild activity/stress (i.e., values slightly higher than SMR) the algorithm has trouble isolating values associated with SMR
# This results in the appearance of a distribution to the left of the MLND in some cases (Chabot et al. 2016)
# If only 1 distribution was fitted and CVMO2 > 15%, the SMR cannot be estimated

CVmlnd_thrshld <- 5.4

# Function to do the above
# n - same as above
# z - same as above

selSMR_by_CVmlnd <- function(n, z) {
  for (i in 1:n) {
    z$SMR[i] <- ifelse(z$CVmlnd[i] > CVmlnd_thrshld, "", z$MLND[i])
    z$SMR_method[i] <- ifelse(z$CVmlnd[i] > CVmlnd_thrshld, "", "MLND")
  }
  z
}

# To save your results in 'z'
z <- selSMR_by_CVmlnd()

#----

# Determine which method (i.e., q20 or q25) is needed to estimate SMR when the CVmlnd is above the threshold value
# If the difference between q20 and q25 exceeds 2% use the q20 value
# If the difference does not exceed 2% use the q25 value
# When the difference between q20 and q25 is low, the q20 is more influenced by outliers

# Function to create a data frame of IDs with a CVmlnd > threshold
# z - same as above

eval_subset <- function(z) {
  subset(z, z$SMR_method == "")
}

# To save your results
b <- eval_subset()

# Function to select the q20 or q25 as the SMR value and record the method used

selSMR_by_q <- function(n, z) {
  for (i in 1:n) {
    if (z$perc_diff_quartile[i] > 2) {
      z$SMR_method[i] <- "q20"
    }
    if (z$SMR_method[i] == "") {
      z$SMR_method[i] <- "q25"
    }
    if (z$perc_diff_quartile[i] > 2) {
      z$SMR[i] <- z$q20[i]
    }
    if (z$SMR[i] == "") {
      z$SMR[i] <- z$q25[i]
    }
  }
  z
}

# To save your results in 'z'
z <- selSMR_by_q()

#----

# Function to create plots for IDs with a CVmlnd > threshold
# b - data frame of IDs with a CVmlnd > threshold
# x - same as above
# y - same as above

eval_d <- eval_SMR <- list()

plot_eval <- function(b, y, x) {
  eval_IDs <- levels(fct_inorder(as.factor(b$ID)))
  eval_CVmlnd <- round(b$CVmlnd, 2)
  for (j in 1:length(eval_IDs)) {
    eval_d[[j]] <- y[[eval_IDs[j]]]
    eval_SMR[[j]] <- x[[eval_IDs[j]]]
  }
  par(
    mfrow = c(2, 3),
    mar = c(4, 4, 4, 2) + 0.1,
    oma = c(1, 1, 1, 1)
  )
  for (i in 1:length(eval_d)) {
    plot(
      eval_d[[i]]$MO2 ~ eval_d[[i]]$Time_hr,
      bty = "l",
      main = eval_IDs[i],
      pch = eval_SMR[[i]]$cl,
      ylab = yl,
      xlab = "Time (hr)"
    )
    abline(
      h = b$q20[i],
      lty = 2,
      col = "red",
      lwd = 2
    )
    abline(
      h = b$q25[i],
      lty = 3,
      col = "green3",
      lwd = 2
    )
    abline(
      h = b$MLND[i],
      lty = 2,
      col = "black",
      lwd = 2
    )

    plot(
      eval_d[[i]]$MO2 ~ eval_d[[i]]$Time_hr,
      ylim = c((eval_SMR[[i]]$mlnd * 0.7), (eval_SMR[[i]]$mlnd * 1.7)),
      bty = "l",
      pch = eval_SMR[[i]]$cl,
      ylab = "",
      xlab = "Time (hr)"
    )
    abline(
      h = eval_SMR[[i]]$quant[1],
      lty = 1,
      col = "red",
      lwd = 2
    )
    abline(
      h = eval_SMR[[i]]$quant[2],
      lty = 3,
      col = "green3",
      lwd = 2
    )
    abline(
      h = eval_SMR[[i]]$mlnd,
      lty = 2,
      col = "black",
      lwd = 2
    )
    rect(
      xleft = 10,
      ybottom = b$q20[i] - (b$q20[i] * 0.1),
      xright = 36,
      ytop = b$q20[i] + (b$q20[i] * 0.1),
      col = rgb(0.8, 0.1, 0.1, 0.2),
      border = NA
    )
    rect(
      xleft = 10,
      ybottom = b$q25[i] - (b$q25[i] * 0.1),
      xright = 36,
      ytop = b$q25[i] + (b$q25[i] * 0.1),
      col = rgb(0, 0.8, 0.1, 0.2),
      border = NA
    )

    plotMO2fdis(eval_d[[i]]$MO2,
      xlim = c(min(eval_d[[i]]$MO2), max(eval_d[[i]]$MO2) + (max(eval_d[[i]]$MO2) *
        0.2)),
      bty = "l"
    )
    abline(
      v = eval_SMR[[i]]$mlnd,
      lty = 1,
      col = "black",
      lwd = 2
    )
    abline(
      v = eval_SMR[[i]]$quant[1],
      lty = 2,
      col = "red",
      lwd = 2
    )
    abline(
      v = eval_SMR[[i]]$quant[2],
      lty = 3,
      col = "green3",
      lwd = 2
    )
    legend(
      "topright",
      c(
        paste("MLND", round(b$MLND[i], 2), sep = "= "),
        paste("q0.2", round(b$q20[i], 2), sep = "= "),
        paste("q0.25", round(b$q25[i], 2), sep = "= "),
        paste("CV of MLND", round(b$CVmlnd[i], 2), sep = "= ")
      ),
      lty = c(1, 2, 3),
      col = c("black", "red", "green3", gray(0, 0)),
      lwd = 2,
      bty = "n",
      y.intersp = 1.2,
      seg.len = 2
    )
  }
  par(def.par)
}
