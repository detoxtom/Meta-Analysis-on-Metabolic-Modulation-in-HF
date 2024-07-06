library(readxl)
library(metafor)
library(tidyverse)
library(ggplot2)

# 4-Level Meta Analyse
V <- vcalc(vi = vi, cluster=study_id, obs=esid, data=data, rho=c)

res <- rma.mv(yi=yi, V, random = ~ 1 | realstudy_id/intrastudy_id/esid, data=data, method = "REML",test = "knha")
res

regression <- rma.mv(yi=yi, V, mods = ~ Disease.model-1, random = ~ 1 | realstudy_id/intrastudy_id/esid, data=data, digits=3, method = "REML",test = "knha")
regression

#R**2
(res$sigma2[3]+res$sigma2[2]+res$sigma2[1]) - (regression$sigma2[3]+regression$sigma2[2]+regression$sigma2[1]) / (res$sigma2[3]+res$sigma2[2]+res$sigma2[1])
#Subsets
FAOi_mv <- rma.mv(yi=yi, V, subset = (Metabolic_Change =="FAO↑"), random = ~ 1 | realstudy_id/intrastudy_id/esid, data=data,method = "REML",test = "knha")
FAOd_mv <- rma.mv(yi=yi, V, subset = (Metabolic_Change =="FAO↓"), random = ~ 1 | realstudy_id/intrastudy_id/esid, data=data, method = "REML",test = "knha")
GOi_mv <- rma.mv(yi=yi, V, subset = (Metabolic_Change =="GO↑"), random = ~ 1 | realstudy_id/intrastudy_id/esid, data=data, method = "REML",test = "knha")
GOiFAOd_mv <- rma.mv(yi=yi, V, subset = (Metabolic_Change =="GO↑FAO↓"),random = ~ 1 | realstudy_id/intrastudy_id/esid, data=data, method = "REML",test = "knha")
GOdFAOi_mv <- rma.mv(yi=yi, V, subset = (Metabolic_Change =="GO↓FAO↑"), random = ~ 1 | realstudy_id/intrastudy_id/esid,  data=data, method = "REML",test = "knha")

FAOi_mv
FAOd_mv
GOi_mv
GOiFAOd_mv
GOdFAOi_mv
res

aggFAOi
aggFAOd
aggGOi
aggGOiFAOd
aggGOdFAOi
aggres


#Calculate 4-Level Variance
mlm.variance.distribution <- function(x) {
  
  m <- x
  
  # Check class
  if (!(class(m)[1] %in% c("rma.mv", "rma"))) {
    stop("x must be of class 'rma.mv'.")
  }
  
  # Check for four-level model
  if (m$sigma2s != 3) {
    stop("The model you provided does not seem to be a four-level model. This function can only be used for four-level models.")
  }
  

  
  # Get variance diagonal and calculate total variance
  n <- m$k.eff
  vector.inv.var <- 1/(diag(m$V))
  sum.inv.var <- sum(vector.inv.var)
  sum.sq.inv.var <- (sum.inv.var)^2
  vector.inv.var.sq <- 1/(diag(m$V)^2)
  sum.inv.var.sq <- sum(vector.inv.var.sq)
  num <- (n-1)*sum.inv.var
  den <- sum.sq.inv.var - sum.inv.var.sq
  est.samp.var <- num/den
  
  # Calculate variance proportions
  level1 <- ((est.samp.var) / (sum(m$sigma2) + est.samp.var) * 100)
  level2 <- ((m$sigma2[3]) / (sum(m$sigma2) + est.samp.var) * 100)
  level3 <- ((m$sigma2[2]) / (sum(m$sigma2) + est.samp.var) * 100)
  level4 <- ((m$sigma2[1]) / (sum(m$sigma2) + est.samp.var) * 100)
  
  # Prepare df for return
  Level <- c("Level 1", "Level 2", "Level 3", "Level 4")
  Variance <- c(level1, level2, level3, level4)
  df.res <- data.frame(Variance)
  colnames(df.res) <- c("% of total variance")
  rownames(df.res) <- Level
  I2 <- c("---", round(Variance[2:4], 2))
  df.res <- as.data.frame(cbind(df.res, I2))
  
  totalI2 <- sum(Variance[2:4])
  
  returnlist <- list(results = df.res,
                     totalI2 = totalI2)
  class(returnlist) <- c("mlm.variance.distribution", "list")
  
  invisible(returnlist)
  
  returnlist
}

#1: Sampling Error
#2: Between effect sizes
#3: Between individual interventions
#4: Between Paper

mlm.variance.distribution(res) #80.96 1.14
mlm.variance.distribution(FAOi_mv) #78,25 1.21
mlm.variance.distribution(FAOd_mv) #82.92 1.42
mlm.variance.distribution(GOi_mv) # 63.78 0.41
mlm.variance.distribution(GOiFAOd_mv) #65.00 0.49 
mlm.variance.distribution(GOdFAOi_mv) #78.74 0.89

GOdFAOi_mv$sigma2[3]+GOdFAOi_mv$sigma2[2]+GOdFAOi_mv$sigma2[1]


table(aggFAOi_data$Disease.model)
table(aggGOiFAOd_data$Disease.model)
table(aggGOi_data$Disease.model)
table(aggGOdFAOi_data$Disease.model)
table(aggFAOd_data$Disease.model)


#Aggregierte Subsets (gleicher Effekt wie rma.mv)

dFAOi_mv <- subset(data, Metabolic_Change== "FAO↑")
aggFAOi_data <- aggregate(dFAOi_mv, cluster=study_id, V=vcov(FAOi_mv, type="obs"), addk=TRUE)
aggFAOi <- rma(yi, vi, method="EE", data=aggFAOi_data, digits=3,test = "knha" )

dFAOd_mv <- subset(data, Metabolic_Change== "FAO↓")
aggFAOd_data <- aggregate(dFAOd_mv, cluster=study_id, V=vcov(FAOd_mv, type="obs"), addk=TRUE)
aggFAOd <- rma(yi, vi, method="EE", data=aggFAOd_data, digits=3,test = "knha")

dGOi_mv <- subset(data, Metabolic_Change== "GO↑")
aggGOi_data <- aggregate(dGOi_mv, cluster=study_id, V=vcov(GOi_mv, type="obs"), addk=TRUE)
aggGOi <- rma(yi, vi, method="EE", data=aggGOi_data, digits=3,test = "knha")

dGOiFAOd_mv <- subset(data, Metabolic_Change== "GO↑FAO↓")
aggGOiFAOd_data <- aggregate(dGOiFAOd_mv, cluster=study_id, V=vcov(GOiFAOd_mv, type="obs"), addk=TRUE)
aggGOiFAOd <- rma(yi, vi, method="EE", data=aggGOiFAOd_data, digits=3,test = "knha")

dGOdFAOi_mv <- subset(data, Metabolic_Change== "GO↓FAO↑")
aggGOdFAOi_data <- aggregate(dGOdFAOi_mv, cluster=study_id, V=vcov(GOdFAOi_mv, type="obs"), addk=TRUE)
aggGOdFAOi <- rma(yi, vi, method="EE", data=aggGOdFAOi_data, digits=3,test = "knha")

aggres_data <- aggregate(data, cluster=study_id, V=vcov(res, type="obs"), addk=TRUE)
aggres <- rma(yi, vi, method="EE", data=aggres_data, digits=3,test = "knha")
aggres


#Funnelplot Asymmetry (Regression test with inverse of sample size)
rma.mv(yi=yi, V, mod = 1/N2, subset = (Metabolic_Change =="FAO↑"), random = ~ 1 | realstudy_id/intrastudy_id/esid, data=data,  method = "REML",test = "knha",digits=3)
rma.mv(yi=yi, V, mod = 1/N2,subset = (Metabolic_Change =="FAO↓"), random = ~ 1 | realstudy_id/intrastudy_id/esid, data=data,  method = "REML",test = "knha",digits=3)
rma.mv(yi=yi, V, mod = 1/N2,subset = (Metabolic_Change =="GO↑"), random = ~ 1 | realstudy_id/intrastudy_id/esid, data=data,  method = "REML",test = "knha",digits=3)
rma.mv(yi=yi, V, mod = 1/N2,subset = (Metabolic_Change =="GO↑FAO↓"), random = ~ 1 | realstudy_id/intrastudy_id/esid, data=data, method = "REML",test = "knha",digits=3)
rma.mv(yi=yi, V, mod = 1/N2,subset = (Metabolic_Change =="GO↓FAO↑"), random = ~ 1 | realstudy_id/intrastudy_id/esid,  data=data,  method = "REML",test = "knha",digits=3)




regtest(FAOi, predictor= "ninv", model = "rma")
rma(yi, vi, mod = 1/N2, data = dFAOi, ni = dFAOi$N2, digits = 3, test = "knha")


#Funnelplots-----------------------------------------------------------------------------------------------------------
tiff("/Users/ai/Desktop/Promotion/plots/supplfunnel.tiff", height = 15, width = 20, units = "cm", compression = "lzw", res = 300)
par(mar = c(4,4,1,2))
funnel(GOdFAOi_mv, 
       yaxis="seinv",
       pch = 19, col="black",      #symbols
       xlim=c(-4, 5), ylim=c(.00001, 5), 
       xaxs="i", yaxs="i", las=1, 
       level=c(.10, .05, .01), 
       lty = c(0,1),                #lines for (ci, combined effct)
       steps= 5,
       # label= TRUE,
       refline = 0,
       shade=c("white", "gray55", "gray75"), back="grey90",hlines="white",,#"white"
       legend=TRUE, 
       refline2 = c(GOdFAOi_mv$beta),
       xlab = "Hedges' g ", ylab="Precision (1/Standard Error)")
dev.off()


tiff("/Users/ai/Desktop/Promotion/plots/GERMAN.tiff", height = 20, width = 20, units = "cm", compression = "lzw", res = 500)
par(mar=c(5,5,3,2)) # left, top, lef, rig
par(family ="Times New Roman")
x <- funnel(GOiFAOd_mv, 
       yaxis="seinv",
       pch = 19, col="black",      #symbols
       xlim=c(-4, 5), ylim=c(.00001, 4), 
       xaxs="i", yaxs="i", las=1, 
       level=c(.10, .05, .01), 
       lty = c(0,1),                #lines for (ci, combined effct)
       steps= 5,
       # label= TRUE,
       refline = 0,
       shade=c("white", "gray55", "gray75"), back="grey90",hlines="white",
       #legend=(TRUE), 
       # legend=list(show="pvals", bty="n"),
       refline2 = c(GOiFAOd_mv$beta),
       xlab = "Hedges' g ", ylab="Präzision (1/Standardfehler)")

legend("topright",
       c("0.10 < p ≤ 1.00","0.05 < p ≤ 0.10","0.01 < p ≤ 0.05","0.00 < p ≤ 0.01","Studien","95% Pseudo CIs"),
       pch = c(22, 22, 22, 22, 19, NA),  # set the symbol for each point
       pt.cex = c(2,2,2,2,1,1),
       inset = c(0.005,0.01),
       bty = "o",
       lty = c(NA, NA, NA, NA,NA, "dotted") ,
       lwd = c(1,1,1,1,1,1.5),
       pt.bg = c("white","gray55","gray75","gray90"),
       #bg = "#e6f5ff",
       #lty = c(NA, NA, NA, NA, NA, "dotted"),
       col = c("black","black","black","black")) 
dev.off()


#FOREST PLOT------------------------------------------------------------------------------------------------------------
mlabfun <- function(text, res, text2, text3) {
  text2_formatted <- formatC(text2, digits = 4)
  text3_formatted <- formatC(text3, digits = 3)
  list(bquote(bold(paste(.(text),
                         " (Q = ", .(formatC(res$QE, digits=2, format="f")),
                         ", p ", .(metafor:::.pval(res$QEp, digits=2, showeq=TRUE, sep=" 0")), "; ",
                         I^2, " =", .(text2_formatted),  "%, ",
                         σ^2, "= ", .(text3_formatted),")" ))))}

#--FAOI--
tiff("/Users/ai/Desktop/Promotion/plots/FAOi_mv.tiff", height = 15, width = 20, units = "cm", compression = "lzw", res = 300)
par(mar = c(4,4,1,2))
forest(aggFAOi,
       slab = paste(dFAOi$Authors, dFAOi$Year, sep = ", "), digits = c(2,1),
       cex = 0.8, fonts = "Times New Roman",    
       ylim = c(-1.5,18), xlim = c(-12,12),
       ilab = c(dFAOi$Disease.model), ilab.xpos = -6,
       mlab=mlabfun("Heterogeneity", FAOi_mv,78.25, 1.21),
       showweights = TRUE,
       #width = 0,
       efac = c(1,2),
       shade = c(-1.8,-1.5, -0.5),
       colshade = "grey87",
       addfit = TRUE,
       refline = c(FAOi_mv$b), col = "black",
       order = dFAOi$Year,
       at = c(-4,-2, 0, 2, 4, 6),
       xlab = "                              Favours metabolic treatment →")

par(xpd=NA)
rect(0, 16, 0, -3,border = "black")
text(-7.5,16.5, "Disease model", pos = 4, cex = 0.8, font = 2)
text(-12,16.5, "Author and Year", pos = 4, cex = 0.8, font = 2)
text(12,16.5, "Weight and Hedges' g [95% CI]", pos = 2, cex = 0.8, font = 2)
text(-3.5,17.5, "Increased Fatty Acid Oxidation", pos = 4, cex = 1, font = 2)
text(-12, -1.8, pos=4, cex=0.8, bquote(bold(paste("Test for Overall Effect: ",
                                                  t, " = ", .(fmtx(FAOi_mv$zval, digits=2)),
                                                  ", df = ", .(formatC(FAOi_mv$k - FAOi_mv$p),format="f"),
                                                  ", p ", .(metafor:::.pval(FAOi_mv$pval, digits=2, showeq=TRUE, sep=" 0")),))))
rect(-12, -2.5, -9, -4.7,border="grey50", col="grey87", lty="dashed")
text(-12, -3.2, pos=4, cex=0.8,  expression(paste("I" ^ 2,""[(2)]," = 30.33% ",sep = "")))
text(-12, -3.8, pos=4, cex=0.8,  expression(paste("I" ^ 2,""[(3)]," = 23.96% ", sep = "")))
text(-12, -4.4, pos=4, cex=0.8,  expression(paste("I" ^ 2,""[(4)]," = 23.96% ", sep = "")))
rect(7, -0.1, 12, -2,border = "grey87",col = "grey87")
text(7, -1.1,"100%   1.17 [ 0.58,  1.76]", pos=4, cex=0.8, font = 2)
dev.off()


tiff("/Users/ai/Desktop/Promotion/plots/FAOd_mv.tiff", height = 12, width = 20, units = "cm", compression = "lzw", res = 300)
par(mar = c(4,4,1,2))
forest(aggFAOd,
       slab = paste(dFAOd$Authors, dFAOd$Year, sep = ", "), digits = c(2,1),
       cex = 0.8, fonts = "Times New Roman",    
       ylim = c(-1.5 ,16), xlim = c(-15,13),
       ilab = c(dFAOd$Disease.model), ilab.xpos = -8,
       mlab=mlabfun("Heterogeneity", FAOd_mv,82.92,1.42),
       shade = c(-1.5, -0.5,-1.7),
       colshade = "grey87",
       showweights = TRUE,
       efac = c(1,2),
       at = c(-6,-4,-2, 0, 2, 4, 6),
       refline = c(FAOd_mv$b),
       order = dFAOd$Year,
       xlab = "                                                      Favours metabolic treatment →")

par(xpd=NA)
rect(0, 14, 0, -5,border = "black")
text(-8,14.5, "Disease model", cex = 0.8, font = 2)
text(-15,14.5, "Author and Year", pos = 4, cex = 0.8, font = 2)
text(13,14.5, "Weight and Hedges' g [95% CI]", pos = 2, cex = 0.8, font = 2)
text(-4.5,15.5, "Decreased Fatty Acid Oxidation", pos = 4, cex = 1, font = 2)
#text(0.1,-1.8, "Favours metabolic treatment →", pos = 4, cex = 0.8, font = 1)
text(-15, -1.8, pos=4, cex=0.8, bquote(bold(paste("Test for Overall Effect: ",
                                                  t, " = ", .(fmtx(FAOd_mv$zval, digits=2)),
                                                  ", df = ", .(formatC(FAOd_mv$k - FAOd_mv$p),format="f"),
                                                  ", p ", .(metafor:::.pval(FAOd_mv$pval, digits=2, showeq=TRUE, sep=" ")),))))
rect(5, -0.1, 13, -2,border = "grey87",col = "grey87")
text(6.95, -1.1,"100%    0.24 [-0.57,  1.05]", pos=4, cex=0.8, font = 2)

rect(-15, -2.4, -11.5, -4.6,border="grey50", col="grey87", lty="dashed")
text(-15, -3.1, pos=4, cex=0.8,  expression(paste("I" ^ 2,""[(2)]," =   9.65% ",sep = "")))
text(-15, -3.7, pos=4, cex=0.8,  expression(paste("I" ^ 2,""[(3)]," =   0.75% ", sep = "")))
text(-15, -4.3, pos=4, cex=0.8,  expression(paste("I" ^ 2,""[(4)]," = 72.51% ", sep = "")))
dev.off()


tiff("/Users/ai/Desktop/Promotion/plots/GOi_mv.tiff", height = 25, width = 20, units = "cm", compression = "lzw", res = 300)
par(mar = c(4,4,1,2))
forest(aggGOi,
       slab = paste(dGOi$Authors, dGOi$Year, sep = ", "), digits = c(2,1),
       cex = 0.8, fonts = "Times New Roman",
       ylim = c(-0.5,52), xlim = c(-15,12),
       ilab = c(dGOi$Disease.model), ilab.xpos = -8,
       mlab=mlabfun("Heterogeneity", GOi_mv,63.58,0.41),
       shade = c(-1.5, -0.5,-2.5),
       colshade = "grey87",
       refline = c(GOi_mv$b), col = "black",
       order = dGOi$Year,
       showweights = TRUE,
       at = c(-4, -2, 0, 2, 4,6),
       xlab = "                                      Favours metabolic treatment →")

rect(0, 50, 0, -5,border = "black")
text(-8,50.5, "Disease model", cex = 0.8, font = 2)
text(-15,50.5, "Author and Year", pos = 4, cex = 0.8, font = 2)
text(12,50.5, "Weight and Hedges' g [95% CI]", pos = 2, cex = 0.8, font = 2)
text(-4,52.5, "Increased Glucose Oxidation", pos = 4, cex = 1, font = 2)
#text(-0.1,-2, "Favours metabolic treatment →", pos = 4, cex = 0.8, font = 1)
text(-15, -2, pos=4, cex=0.8, bquote(bold(paste("Test for Overall Effect: ",
                                                t, " = ", .(fmtx(GOi_mv$zval, digits=2)),
                                                ", df = ", .(formatC(GOi_mv$k - GOi_mv$p),format="f"),
                                                ", p ", .(metafor:::.pval(GOi_mv$pval, digits=2, showeq=TRUE, sep=" 0")),))))
rect(5, -0.1, 13, -2,border = "grey87",col = "grey87")
text(6.2, -1.1,"100%    1.03 [ 0.79,  1.26]", pos=4, cex=0.8, font = 2)

par(xpd=NA)
rect(-15, -2.8, -11.5, -5.8,border="grey50", col="grey87", lty="dashed")
text(-15, -3.6, pos=4, cex=0.8,  expression(paste("I" ^ 2,""[(2)]," = 18.19% ",sep = "")))
text(-15, -4.5, pos=4, cex=0.8,  expression(paste("I" ^ 2,""[(3)]," = 36.06% ", sep = "")))
text(-15, -5.4, pos=4, cex=0.8,  expression(paste("I" ^ 2,""[(4)]," =   9.33% ", sep = "")))
dev.off()


#--GOiFAOd--


mlabfun00 <- function(text, res, text2, text3) {
  text2_formatted <- formatC(text2, digits = 2)
  text3_formatted <- formatC(text3, digits = 3)
  list(bquote(bold(paste(.(text),
                         " (Q = ", .(formatC(res$QE, digits=2, format="f")),
                         ", p ", .(metafor:::.pval(res$QEp, digits=2, showeq=TRUE, sep=" 0")), "; ",
                         I^2, " =", .(text2_formatted),  ".00%, ",
                         σ^2, "= ", .(text3_formatted),")" ))))}

tiff("/Users/ai/Desktop/Promotion/plots/GOiFAOd_mv.tiff", height = 15, width = 20, units = "cm", compression = "lzw", res = 300)
par(mar = c(4,4,1,2))
forest(aggGOiFAOd,
       slab = paste(dGOiFAOd$Authors, dGOiFAOd$Year, sep = ", "), digits = c(2,1),
       cex = 0.8, fonts = "Times New Roman",    
       ylim = c(-1.3,19), xlim = c(-15,12),
       ilab = c(dGOiFAOd$Disease.model), ilab.xpos = -8,
       mlab=mlabfun00("Heterogeneity", GOiFAOd_mv,65.00,0.49),
       showweights = TRUE,
       shade = c(-1.5, -0.5,-2.5),
       colshade = "grey87",
       at = c(-4,-2,0,2,4,6),
       refline = c(GOiFAOd_mv$b), col = "black",
       order = dGOiFAOd$Year,
       xlab = "                                       Favours metabolic treatment →")

rect(0, 17, 0, -5,border = "black")
text(-8,17.5, "Disease model", cex = 0.8, font = 2)
text(-15,17.5, "Author and Year", pos = 4, cex = 0.8, font = 2)
text(12,17.5, "Weight and Hedges' g [95% CI]", pos = 2, cex = 0.8, font = 2)
text(-5,18.5, "Decreased FAO and Increased GO", pos = 4, cex = 1, font = 2)
#text(-0.1,-1.8, "Favours metabolic treatment →", pos = 4, cex = 0.8, font = 1)
text(-15, -1.8, pos=4, cex=0.8, bquote(bold(paste("Test for Overall Effect: ",
                                                  t, " = ", .(fmtx(GOiFAOd_mv$zval, digits=2)),
                                                  ", df = ", .(formatC(GOiFAOd_mv$k - GOiFAOd_mv$p),format="f"),
                                                  ", p ", .(metafor:::.pval(GOiFAOd_mv$pval, digits=3, showeq=TRUE, sep=" ")),))))
rect(5, -0.1, 13, -2,border = "grey87",col = "grey87")
text(6.15, -1.1,"100%   0.44 [-0.003, 0.89]", pos=4, cex=0.8, font = 2)

par(xpd=NA)
rect(-15, -2.4, -11.5, -4.5,border="grey50", col="grey87", lty="dashed")
text(-15, -3, pos=4, cex=0.8,  expression(paste("I" ^ 2,""[(2)]," = 17.98% ",sep = "")))
text(-15, -3.6, pos=4, cex=0.8,  expression(paste("I" ^ 2,""[(3)]," = 23.51% ", sep = "")))
text(-15, -4.2, pos=4, cex=0.8,  expression(paste("I" ^ 2,""[(4)]," = 23.51% ", sep = "")))
dev.off()


tiff("/Users/ai/Desktop/Promotion/plots/GOdFAOi_mv.tiff", height = 10, width = 20, units = "cm", compression = "lzw", res = 300)
par(mar = c(4,4,1,2))
forest(aggGOdFAOi,
       slab = paste(dGOdFAOi$Authors, dGOdFAOi$Year, sep = ", "), digits = c(2,1),
       cex = 0.8, fonts = "Times New Roman",    
       ylim = c(-1.5,15), xlim = c(-15,13),
       ilab = c(dGOdFAOi$Disease.model), ilab.xpos = -8,
       mlab=mlabfun("Heterogeneity", GOdFAOi_mv,78.73,0.89),
       shade = c(-1.5, -0.5, -2.5),
       colshade = "grey87",
       showweights = TRUE,
       efac = c(1,2),
       refline = c(GOdFAOi_mv$b), col = "black",
       order = dGOdFAOi$Year,
       at = c(-4,-2, 0, 2, 4, 6),
       xlab = "                    Favours metabolic treatment →")

rect(0, 13, 0, -5,border = "black")
text(-8,13.5, "Disease model", cex = 0.8, font = 2)
text(-15,13.5, "Author and Year", pos = 4, cex = 0.8, font = 2)
text(13,13.5, "Weight and Hedges' g [95% CI]", pos = 2, cex = 0.8, font = 2)
text(-5,15, "Increased FAO and Decreased GO", pos = 4, cex = 1, font = 2)
#text(0.05,-2.2, "Favours metabolic treatment →", pos = 4, cex = 0.8, font = 1)
text(-15, -1.8, pos=4, cex=0.8, bquote(bold(paste("Test for Overall Effect: ",
                                                t, " = ", .(fmtx(GOdFAOi_mv$zval, digits=2)),
                                                ", df = ", .(formatC(GOdFAOi_mv$k - GOdFAOi_mv$p),format="f"),
                                                ", p ", .(metafor:::.pval(GOdFAOi_mv$pval, digits=2, showeq=TRUE, sep=" ")),))))
rect(5, -0.1, 13, -2,border = "grey87",col = "grey87")
text(6.85, -1.1,"100%   -0.03 [-0.77,  0.71]", pos=4, cex=0.8, font = 2)

par(xpd=NA)
rect(-15, -2.4, -11.5, -4.8,border="grey50", col="grey87", lty="dashed")
text(-15, -3.1, pos=4, cex=0.8,  expression(paste("I" ^ 2,""[(2)]," = 22.04% ",sep = "")))
text(-15, -3.8, pos=4, cex=0.8,  expression(paste("I" ^ 2,""[(3)]," =   0.00% ", sep = "")))
text(-15, -4.5, pos=4, cex=0.8,  expression(paste("I" ^ 2,""[(4)]," = 56.69% ", sep = "")))
dev.off()



#LEAVE ONE OUT----------------------------------------------------------------------------------------------------------

# FAOi
l1o <- leave1out(aggFAOi)
l1o
l1oforest <- aggFAOi_data
tiff("/Users/ai/Desktop/Promotion/plots/L1O_FAOi.tiff", height = 15, width = 20, units = "cm", compression = "lzw", res = 300)
par(mar = c(5,4,1,2))
x <- forest(l1o$estimate, 
            slab = paste(l1oforest$Authors, l1oforest$Year, sep = ", "), digits = c(2,1,2),
            sei=l1o$se,
            fonts = "Times New Roman",
            cex = 1.2,
            order = l1o$estimate, 
            xlab="Leave One Out Estimate", 
           ilab =  rep("< 0.01",times=15), #ilab = round(l1o$pval,digits = 3), 
            refline=coef(FAOi_mv))

text(0.36,17.8, "Increased Fatty Acid Oxidation", pos = 4, cex = 1.2, font = 2)
text(0,16.5, "p-value", cex = 1, font = 2)
text(-1.4,16.5, "Author and Year", pos = 4, cex = 1, font = 2)
text(2.9,16.5, "Hedges' g [95% CI]", pos = 2, cex = 1, font = 2)
dev.off()


# GOi
l1o <- leave1out(aggGOi)
l1o
l1oforest <- aggGOi_data
tiff("/Users/ai/Desktop/Promotion/plots/L1O_GOi.tiff", height = 25, width = 20, units = "cm", compression = "lzw", res = 300)
par(mar = c(5,4,1,2))
x <-forest(l1o$estimate, 
           slab = paste(l1oforest$Authors, l1oforest$Year, sep = ", "), digits = c(2,1,2),
           sei=l1o$se,
           fonts = "Times New Roman",
           cex = 1.2,
           order = l1o$estimate, 
           xlab="Leave One Out Estimate", 
           ilab = rep("< 0.001",times=49), #round(l1o$pval,digits = 5), # 
           refline=coef(GOi_mv))

text(0.63,51, "p-value", cex = 1, font = 2)
text(0.1,51, "Author and Year", pos = 4, cex = 1, font = 2)
text(1.7,51, "Hedges' g [95% CI]", pos = 2, cex = 1, font = 2)
text(0.68,52.5, "Increased Glucose Oxidation", pos = 4, cex = 1.2, font = 2)
dev.off()
x$xlim
x$ylim

# FAOd
l1o <- leave1out(aggFAOd)
l1oforest <- aggFAOd_data
tiff("/Users/ai/Desktop/Promotion/plots/L1O_FAOd.tiff", height = 15, width = 20, units = "cm", compression = "lzw", res = 300)
par(mar = c(5,4,1,2))
x <- forest(l1o$estimate, 
            slab = paste(l1oforest$Authors, l1oforest$Year, sep = ", "), digits = c(2,1,2),
            sei=l1o$se,
            #ylim = c(-3, 14),
            fonts = "Times New Roman",
            cex = 1.2,
            order = l1o$estimate, 
            xlab="Leave One Out Estimate", 
            ilab =  round(l1o$pval,digits = 3),# rep("p < 0.01",times=15), # 
            ilab.xpos = -1,
            refline=coef(FAOd_mv))

text(-1,15.5, "Decreased Fatty Acid Oxidation", pos = 4, cex = 1.2, font = 2)
text(-1,14.5, "p-value", cex = 1, font = 2)
text(-3,14.5, "Author and Year", pos = 4, cex = 1, font = 2)
text(2.6,14.5, "Hedges' g [95% CI]", pos = 2, cex = 1, font = 2)
dev.off()


# GOiFAOd
l1o <- leave1out(aggGOiFAOd)
l1oforest <- aggGOiFAOd_data
tiff("/Users/ai/Desktop/Promotion/plots/L1O_GOiFAOd.tiff", height = 15, width = 20, units = "cm", compression = "lzw", res = 300)
par(mar = c(5,4,1,2))
x <- forest(l1o$estimate, 
            slab = paste(l1oforest$Authors, l1oforest$Year, sep = ", "), digits = c(2,1,2),
            sei=l1o$se,
            fonts = "Times New Roman",
            cex = 1.2,
            order = l1o$estimate, 
            xlab="Leave One Out Estimate", 
            ilab =  round(l1o$pval,digits = 3), #rep("p < 0.01",times=15), 
            refline=coef(GOiFAOd_mv))

text(-0.3,18.7, "Decreased FAO and Increased GO", pos = 4, cex = 1.2, font = 2)
text(-0.43,17.5, "p-value", cex = 1, font = 2)
text(-1.5,17.5, "Author and Year", pos = 4, cex = 1, font = 2)
text(1.8,17.5, "Hedges' g [95% CI]", pos = 2, cex = 1, font = 2)
dev.off()
x$xlim
x$ylim

# GOdFAOi
l1o <- leave1out(aggGOdFAOi)
l1oforest <- aggGOdFAOi_data
tiff("/Users/ai/Desktop/Promotion/plots/L1O_GOdFAOi.tiff", height = 15, width = 20, units = "cm", compression = "lzw", res = 300)
par(mar = c(5,4,1,2))
x <- forest(l1o$estimate, 
            slab = paste(l1oforest$Authors, l1oforest$Year, sep = ", "), digits = c(2,1,2),
            sei=l1o$se,
            fonts = "Times New Roman",
            cex = 1.2,
            order = l1o$estimate, 
            xlab="Leave One Out Estimate", 
            ilab =  round(l1o$pval,digits = 3), #rep("p < 0.01",times=15), 
            refline=coef(GOdFAOi_mv))

text(-0.8,14.2, "Increased FAO and Decreased GO", pos = 4, cex = 1.2, font = 2)
text(-0.95,13.5, "p-value", cex = 1, font = 2)
text(-2.1,13.5, "Author and Year", pos = 4, cex = 1, font = 2)
text(1.5,13.5, "Hedges' g [95% CI]", pos = 2, cex = 1, font = 2)
dev.off()
x$xlim
x$ylim










#FULL FOREST PLOT-----------------------------------------------------------------------------------------------------------

GOiFAOi_mv <- rma.mv(yi=yi, V, subset = (Metabolic_Change =="GO↑FAO↑"),random = ~ 1 | realstudy_id/intrastudy_id/esid, data=data, method = "REML",test = "knha")
GOdFAOd_mv <- rma.mv(yi=yi, V, subset = (Metabolic_Change =="GO↓FAO↓"), random = ~ 1 | realstudy_id/intrastudy_id/esid,  data=data, method = "REML",test = "knha")
GOd_mv <- rma.mv(yi=yi, V, subset = (Metabolic_Change =="GO↓"), random = ~ 1 | realstudy_id/intrastudy_id/esid,  data=data, method = "REML",test = "knha")

dGOiFAOi_mv <- subset(data, Metabolic_Change== "GO↑FAO↑")
aggGOiFAOi_data <- aggregate(dGOiFAOi_mv, cluster=study_id, V=vcov(GOiFAOi_mv, type="obs"), addk=TRUE)
aggGOiFAOi <- rma(yi, vi, method="EE", data=aggGOiFAOi_data, digits=3,test = "knha" )
dGOdFAOd_mv <- subset(data, Metabolic_Change== "GO↓FAO↓")
aggGOdFAOd_data <- aggregate(dGOdFAOd_mv, cluster=study_id, V=vcov(GOdFAOd_mv, type="obs"), addk=TRUE)
aggGOdFAOd_data
aggGOdFAOd <- rma(yi, vi, method="EE", data=aggGOdFAOd_data, digits=3,test = "knha")
dGOd_mv <- subset(data, Metabolic_Change== "GO↓")
aggGOd_data <- aggregate(dGOd_mv, cluster=study_id, V=vcov(GOd_mv, type="obs"), addk=TRUE)
aggGOd <- rma(yi, vi, method="EE", data=aggGOd_data, digits=3,test = "knha")

mlm.variance.distribution(res) #,80.96, 1.14
mlm.variance.distribution(FAOi_mv) #,78,25, 1.21
mlm.variance.distribution(FAOd_mv) #,82.92, 1.42
mlm.variance.distribution(GOi_mv) # ,63.78, 0.41
mlm.variance.distribution(GOiFAOd_mv) #,65.00, 0.49
mlm.variance.distribution(GOdFAOi_mv) #,78.74, 0.89
mlm.variance.distribution(GOiFAOi_mv) # ,57.92, 0.37
mlm.variance.distribution(GOdFAOd_mv) #,60.50, 0.79
mlm.variance.distribution(GOd_mv) #,85.49, 1.93

GOiFAOi_mv$sigma2[3]+GOiFAOi_mv$sigma2[2]+GOiFAOi_mv$sigma2[1]


tiff("/Users/ai/Desktop/Promotion/plots/Forest Plot with subgroups24.tiff", height = 65, width = 40, units = "cm", compression = "lzw", res = 300)
par(mar = c(5,4,1,2))
forest(aggres,
       order = factor(aggres_data$Metabolic_Change, levels = c("FAO↑", "FAO↓", "GO↑", "GO↓", "GO↑FAO↓", "GO↓FAO↑","GO↑FAO↑","GO↓FAO↓")),
       slab = paste(aggres_data$Authors, aggres_data$Year, sep = ", "), digits = c(2,1),
       cex = 1, fonts = "Times New Roman",
       ylim = c(2,165), xlim = c(-15,12),
       ilab = c(aggres_data$Intervention), ilab.xpos = -8,
       xlab = "                    Favours metabolic treatment →",
       refline = c(res$b), col = "black",
       mlab=mlabfun("Heterogeneity", res, 80.96, 1.14),
       showweights = TRUE,
       shade = c(-1.8,-1.5, -0.5,-2,-3,-3.5,-4),
       colshade = "grey87",
       at = c(-6, -4, -2, 0, 2, 4, 6, 8),
       rows=c(158:144, 138:126, 120:72, 66:60, 54:39, 33:22, 16:11, 5:4)                                              #15 50 13 7 16 12 6 2
) 

par(cex = 1, font = 4)
rect(0, 163, 0, -4,border = "black")
text(-15, c(159.5, 139.5, 121.5, 67.5, 55.5, 34.5, 17.5, 6.5), pos = 4, c("FAO↑","FAO↓", "GO↑",  "GO↓", "FAO↓GO↑", "GO↓FAO↑","FAO↑GO↑","FAO↓GO↓"))
text(-8,164, "Intervention", cex = 1.2, font = 2)
text(-15,164, "Author and Year", pos = 4, cex = 1.2, font = 2)
text(12,164, "Weight Hedges' g [95% CI]", pos = 2, cex = 1.2, font = 2)

addpoly(FAOi_mv, efac  = (0.5),row = 142, cex = 1, mlab = mlabfun("Heterogeneity for Subgroup", FAOi_mv,78.25, 1.21), col="grey", border="gray")
addpoly(FAOd_mv, efac  = (0.5),row = 124, cex = 1, mlab = mlabfun("Heterogeneity for Subgroup", FAOd_mv,82.92, 1.42), col="gray", border="gray")
addpoly(GOi_mv, efac  = (0.5), row = 70, cex = 1, mlab = mlabfun("Heterogeneity for Subgroup", GOi_mv,63.78, 0.41), col="gray", border="gray")
addpoly(GOd_mv, efac  = (0.3),row = 58, cex = 1, mlab = mlabfun("Heterogeneity for Subgroup", GOd_mv,85.49, 1.93), col="gray", border="gray")
addpoly(GOiFAOd_mv, efac  = (0.5),row = 37, cex = 1, mlab = mlabfun("Heterogeneity for Subgroup", GOiFAOd_mv,65.01, 0.49), col="gray", border="gray")
addpoly(GOdFAOi_mv, efac  = (0.2),row = 20, cex = 1, mlab = mlabfun("Heterogeneity for Subgroup", GOdFAOi_mv,78.74, 0.89), col="gray", border="gray")
addpoly(GOiFAOi_mv, efac  = (0.5),row = 9, cex = 1, mlab = mlabfun("Heterogeneity for Subgroup", GOiFAOi_mv,57.92, 0.37), col="gray", border="gray")
addpoly(GOdFAOd_mv, efac  = (0.1),row = 2, cex = 1, mlab = mlabfun("Heterogeneity for Subgroup", GOdFAOd_mv,85.49, 1.93), col="gray", border="gray")

### test for subgroup differences
#res2 <- rma.mv(yi=yi, V, mods = ~Metabolic_Change, subset = (Metabolic_Change =="FAO↑"|Metabolic_Change =="FAO↓"|Metabolic_Change =="GO↑"|Metabolic_Change =="GO↑FAO↓"|Metabolic_Change =="GO↓FAO↑"),random = ~ 1 | realstudy_id/intrastudy_id/esid, data=data, method = "REML",test = "knha")
res2 <- rma.mv(yi=yi, V, mods = ~Metabolic_Change,random = ~ 1 | realstudy_id/intrastudy_id/esid, data=data, method = "REML",test = "knha")

text(-15, -3.9, pos=4, cex=1, bquote(bold(paste("Test for Subgroup Differences: ",
                                                Q[M], " = ", .(fmtx(res2$QM, digits=2)),
                                                ", df = ", .(res2$p - 1), ", ",
                                                .(fmtp(res2$QMp, digits=2, pname="p", add0=TRUE, sep=TRUE, equal=TRUE))))))

par(xpd=NA)
rect(-15, -4.8, -12.8, -8.8,border="grey50", col="grey87", lty="dashed")

text(-15, -5.8, pos=4, cex=1,  expression(paste("I" ^ 2,""[(2)]," = 30.33% ",sep = "")))
text(-15, -7, pos=4, cex=1,  expression(paste("I" ^ 2,""[(3)]," = 23.96% ", sep = "")))
text(-15, -8.2, pos=4, cex=1,  expression(paste("I" ^ 2,""[(4)]," = 23.96% ", sep = "")))

text(-15, -2.5, pos=4, cex=1, bquote(bold(paste("Test for Overall Effect: ",
                                                  t, " = ", .(fmtx(res$zval, digits=2)),
                                                  ", df = ", .(formatC(res$k - res$p),format="f"),
                                                  ", p ", .(metafor:::.pval(res$pval, digits=2, showeq=TRUE, sep=" 0")),))))


rect(7, -0.1, 12, -2,border = "grey87",col = "grey87")
text(8, -1.1,"      100%    0.69  [ 0.47,  0.92]", pos=4, cex=1, font = 2)
dev.off()








 var.comp = function(x){
  
  m = x
  
  # Check class
  if (!(class(m)[1] %in% c("rma.mv", "rma"))){
    stop("x must be of class 'rma.mv'.")
  }
  
  # Check for three level model
  if (m$sigma2s != 2){
    stop("The model you provided does not seem to be a three-level model. This function can only be used for three-level models.")
  }
  
  # Check for right specification (nested model)
  if (sum(grepl("/", as.character(m$random[[1]]))) < 1){
    stop("Model must contain nested random effects. Did you use the '~ 1 | cluster/effect-within-cluster' notation in 'random'? See ?metafor::rma.mv for more details.")
  }
  
  # Get variance diagonal and calculate total variance
  n = m$k.eff
  vector.inv.var = 1/(diag(m$V))
  sum.inv.var = sum(vector.inv.var)
  sum.sq.inv.var = (sum.inv.var)^2
  vector.inv.var.sq = 1/(diag(m$V)^2)
  sum.inv.var.sq = sum(vector.inv.var.sq)
  num = (n-1)*sum.inv.var
  den = sum.sq.inv.var - sum.inv.var.sq
  est.samp.var = num/den
  
  # Calculate variance proportions
  level1=((est.samp.var)/(m$sigma2[1]+m$sigma2[2]+est.samp.var)*100)
  level2=((m$sigma2[2])/(m$sigma2[1]+m$sigma2[2]+est.samp.var)*100)
  level3=((m$sigma2[1])/(m$sigma2[1]+m$sigma2[2]+est.samp.var)*100)
  
  # Prepare df for return
  Level=c("Level 1", "Level 2", "Level 3")
  Variance=c(level1, level2, level3)
  df.res=data.frame(Variance)
  colnames(df.res) = c("% of total variance")
  rownames(df.res) = Level
  I2 = c("---", round(Variance[2:3], 2))
  df.res = as.data.frame(cbind(df.res, I2))
  
  totalI2 = Variance[2] + Variance[3]
  
  
  # Generate plot
  df1 = data.frame("Level" = c("Sampling Error", "Total Heterogeneity"),
                   "Variance" = c(df.res[1,1], df.res[2,1]+df.res[3,1]),
                   "Type" = rep(1,2))
  
  df2 = data.frame("Level" = rownames(df.res),
                   "Variance" = df.res[,1],
                   "Type" = rep(2,3))
  
  df = as.data.frame(rbind(df1, df2))
  
  
  g = ggplot(df, aes(fill=Level, y=Variance, x=as.factor(Type))) +
    coord_cartesian(ylim = c(0,1), clip = "off") +
    geom_bar(stat="identity", position="fill", width = 1, color="black") +
    scale_y_continuous(labels = scales::percent)+
    theme(axis.title.x=element_blank(),
          axis.text.y = element_text(color="black"),
          axis.line.y = element_blank(),
          axis.title.y=element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_line(lineend = "round"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.background = element_rect(linetype="solid",
                                           colour ="black"),
          legend.title = element_blank(),
          legend.key.size = unit(0.75,"cm"),
          axis.ticks.length=unit(.25, "cm"),
          plot.margin = unit(c(1,3,1,1), "lines")) +
    scale_fill_manual(values = c("darkseagreen3", "deepskyblue3", "darkseagreen2",
                                 "deepskyblue1", "deepskyblue2")) +
    
    # Add Annotation
    
    # Total Variance
    annotate("text", x = 1.5, y = 1.05,
             label = paste("Total Variance:",
                           round(m$sigma2[1]+m$sigma2[2]+est.samp.var, 3))) +
    
    # Sampling Error
    annotate("text", x = 1, y = (df[1,2]/2+df[2,2])/100,
             label = paste("Sampling Error Variance: \n", round(est.samp.var, 3)), size = 3) +
    
    # Total I2
    annotate("text", x = 1, y = ((df[2,2])/100)/2-0.02,
             label = bquote("Total"~italic(I)^2*":"~.(round(df[2,2],2))*"%"), size = 3) +
    annotate("text", x = 1, y = ((df[2,2])/100)/2+0.05,
             label = paste("Variance not attributable \n to sampling error: \n", round(m$sigma2[1]+m$sigma2[2],3)), size = 3) +
    
    # Level 1
    annotate("text", x = 2, y = (df[1,2]/2+df[2,2])/100, label = paste("Level 1: \n",
                                                                       round(df$Variance[3],2), "%", sep=""), size = 3) +
    
    # Level 2
    annotate("text", x = 2, y = (df[5,2]+(df[4,2]/2))/100,
             label = bquote(italic(I)[Level2]^2*":"~.(round(df[4,2],2))*"%"), size = 3) +
    
    # Level 3
    annotate("text", x = 2, y = (df[5,2]/2)/100,
             label = bquote(italic(I)[Level3]^2*":"~.(round(df[5,2],2))*"%"), size = 3)
  
  returnlist = list(results = df.res,
                    totalI2 = totalI2,
                    plot = g)
  class(returnlist) = c("mlm.variance.distribution", "list")
  
  invisible(returnlist)
  
  returnlist
  
}
