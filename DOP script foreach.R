#####################################
# DOP test analysis - 2022
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: Oct 5 2022
#####################################
study <- "DOP"
###################### loading files
#loading files
data <- read.csv("phenotype.csv", skip = 3)
pipeline <- paste0(unique(data$studyYear), "_", study)

data$plantingDate <- stringr::str_replace_all(data$plantingDate, 
                      c( 
                     "February" =  "02",
                     "March" = "03",
                     "April" = "04",
                     "May" = "05"))
trials <- unique(data$plantingDate)  

# coerce col to factors or numeric
library(magrittr)
library(dplyr)
cols <- colnames(data)[1:30]
data %<>%
  mutate_each_(funs(factor(.)), cols)
cols <- colnames(data)[length(colnames(data))-1]
data %<>%
  mutate_each_(funs(as.numeric(.)), cols)

data$row <- as.numeric(data$rowNumber)
data$col <- as.numeric(data$colNumber)

#################################
# fitting a model for each trial
#################################
require(foreach)
require(doParallel)
require(doMC)
library(SpATS)
library(lme4)
library(car)
library(data.table)
library(ggplot2)

# setting the number of cores that will be used
detectCores()
registerDoParallel(cores = detectCores()) # type the number of cores you want to use
getDoParWorkers()

#running all the single-trials in parallel 
system.time(
results.st <- foreach(i = 1:length((trials)), 
                          .packages = c("SpATS", "car", "data.table"), 
                          .combine = "rbind", 
                          .multicombine = TRUE, 
                          .errorhandling = "remove",
                          .verbose = TRUE
 ) %dopar% {

      # subset the data  
    sample <- droplevels.data.frame(data[data$plantingDate == trials[i],])

      outlier <- boxplot.stats(sample$Yield.LSU_01.0000138)$out
      log.outliers <- data.frame(sample[sample$Yield.LSU_01.0000138 %in% outlier,])
      sample$Yield.LSU_01.0000138[sample$Yield.LSU_01.0000138 %in% outlier] <- NA
      
      nrow <- max(sample$row)
      ncol <- max(sample$col)
      nseg.row <- nrow
      nseg.col <- ncol
      
      fitF <- SpATS(response = "Yield.LSU_01.0000138", 
                                fixed = ~ 1, 
                                random = ~ replicate + rowNumber + colNumber, 
                                spatial = ~ PSANOVA(col, row, nseg = c(nseg.col, nseg.row)), 
                                genotype = "germplasmName", 
                                genotype.as.random = FALSE, 
                                data = sample)
      
      # obtaining the spatial trends - from raw data to BLUES and BLUPS
      #plot.SpATS(fitF) 
    
      # Estimate BLUEs
      blues <- predict.SpATS(fitF, which = "germplasmName")
      blues <- blues[,c(1,7,8)]
      colnames(blues)[2:3] <- c("BLUE", "sep_BLUE")
      
      # Now, run as random
      fitR <- SpATS(response = "Yield.LSU_01.0000138", 
                    fixed = ~ 1, 
                    random = ~ replicate + rowNumber + colNumber, 
                    spatial = ~ PSANOVA(col, row, nseg = c(nseg.col, nseg.row)), 
                    genotype = "germplasmName", 
                    genotype.as.random = TRUE, 
                    data = sample)

      # to obtain the heritability via the package function we can use 
      h2g <- getHeritability(fitR)
      # Broad-sense heritability based on Cullis method
      Vg <- fitR$var.comp["germplasmName"]
      ng <- length(unique(sample$germplasmName))
      C11_g <- fitR$vcov$C11_inv
      trC11_g <-sum(diag(C11_g))
      av2 <- 2/ng * (trC11_g - (sum(C11_g) - trC11_g) / ng-1) # mean var of a difference between genotypic BLUPS
      H2.Cullis <- 1 - av2 / (2 * Vg)
      # Estimate BLUPs for Grain Yield
      blups <- predict.SpATS(fitR, which = "germplasmName")
      blups <- blups[,c(7,8)]
      colnames(blups)[1:2] <- c("BLUP", "sep_BLUP")
      # Reliability
      rel <- mean(1 - blups$sep_BLUP^2 / fitR$var.comp["germplasmName"])
      # weights for ID's - adjust residual for further analysis
      vcov.mme <- fitR$vcov$C11_inv
      w <- diag(vcov.mme)
      
      output <- data.frame(blues,
                           w = w,
                           blups,
                           DOP = as.character(unique(sample$plantingDate)), 
                           h2g = h2g,
                           r = rel,
                           H.cullis = H2.Cullis,
                           trait = fitR$model$response[1]
                          )
        }
)

# saving the output file for single-trials analysis
output <- results.st
write.csv(output, paste0(pipeline, "_", "output_single-trials.csv"))

# estimate how many trials were eliminated by singularities
length(trials) - length(unique(output$DOP))
# eliminate trials with low Cullis heritability
unique(output$H.cullis)
output <- output[output$H.cullis > 0.50,]
# the final dimensions after QC
length(unique(output$DOP))

##################################################
# second step - joint analysis
##################################################

# Fitting genotype by environment model - joint analysis
library(sommer)
fitMET <- mmer(BLUE ~ 1,
               random= ~ DOP + germplasmName + germplasmName:DOP,
               weights = w,
               rcov= ~ units,
               data=output, 
               verbose = FALSE)

# Broad-sense heritability
nloc <- length(unique(output$DOP))
h2g.MET <- vpredict(fitMET, h2 ~ V2 / ( V2 + V3/nloc)) # MET level

# predicting the BLUP - main effect
BLUPs.main <- predict.mmer(object = fitMET, classify = "germplasmName")
# reliability
rel <- mean(1 - BLUPs.main$pvals$standard.error^2 / summary(fitMET)$varcomp[2,1])

# predicting the BLUP per enviroment
BLUPs.MET <- predict.mmer(object = fitMET, classify = c("germplasmName","DOP"))

rice.handbook.table <- round(reshape2::acast(BLUPs.MET$pvals, germplasmName ~  DOP, value.var = "predicted.value"))
write.csv(rice.handbook.table, paste0(pipeline, "_rice.handbook.table.csv"))

# barplot graph with confidence interval using main
data.plot <- BLUPs.main$pvals

limits <- aes(ymax = data.plot$predicted.value + data.plot$standard.error*1.96,
              ymin = data.plot$predicted.value - data.plot$standard.error*1.96)

p <- ggplot(data = data.plot, aes(x = reorder(germplasmName, -predicted.value), y = predicted.value, color = reorder(germplasmName, -predicted.value))) + 
  geom_bar(position = position_dodge(), 
         stat="identity") +
  scale_fill_brewer(palette="Set1") +
  geom_errorbar(limits, position = position_dodge(),
                width = 0.5) +
  labs(x = "Variety", y = "Yield") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  theme(legend.position = "bottom", legend.box = "horizontal") +
  annotate("text", x=length(unique(data.plot$germplasmName))/2, 
           y=max(data.plot$predicted.value), 
           label= paste("Reliability = ", round(sqrt(h2g.MET[1]), 2)*100, "%")) +
  labs(x = "", fill = "DOP") +
  labs(color = "") + 
  coord_cartesian(ylim=c(min(data.plot$predicted.value - data.plot$standard.error*1.96),
                         max(data.plot$predicted.value + data.plot$standard.error*1.96)))


ggsave(filename = paste0("./", pipeline, '_overall_performances.tiff'),
       plot = p,
       device = 'tiff',
       width = 280,
       height = 140,
       units = 'mm',
       dpi = 300)

data.plot$standard.error <- data.plot$standard.error*1.96
colnames(data.plot) <- c("Estimate", "Variety", "Yield", "Confidence_interval")
# saving the output file for single-trials analysis
write.csv(data.plot, paste0(pipeline, "_", "overall_performances.csv"))

#######################################
# third step - GGE MET analysis
#######################################
library(statgenGxE)

data.MET <- BLUPs.MET$pvals
colnames(data.MET)[c(2,4)] <- c("Variety", "Yield") 

## Create a TD object
dropsTD <- statgenSTA::createTD(data = data.MET, genotype = "Variety", trial = "DOP")
## Fit a model where trials are nested within scenarios.
dropsVarComp <- gxeVarComp(TD = dropsTD, trait = "Yield")
# Finlay-Wilkinson Analysis
dropsFW <- gxeFw(TD = dropsTD, trait = "Yield")

# reorganizing the data
test <- matrix(dropsFW$TD[[1]]$beta) %*% matrix(sort(dropsFW$envEffs[,2]), nrow = 1) + dropsFW$estimates[,2]
locations <- c(trials) 
locations <- as.character(locations[match(sort(dropsFW$envEffs[,2]), dropsFW$envEffs[,2])])
colnames(test) <- sort(dropsFW$envEffs[,2])
rownames(test) <- dropsFW$estimates[, 1]
test <- data.frame(Variety = rownames(test), test)
colnames(test)[2:ncol(test)] <- sort(dropsFW$envEffs[,2])
test <- melt(test)
test$variable <- as.numeric(as.character(test$variable))


## Create line plot for Finlay Wilkinson analysis.
q <- ggplot(data = test, aes(x = variable, 
                               y = value, 
                               group = Variety, 
                               colour = Variety)) + 
  geom_line() + geom_point() + 
  labs(x = "DOP", y = "Yield") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_vline(xintercept = meanEnV <- mean(dropsFW$envEffs[,2]), colour = "grey") +
  scale_x_continuous(breaks = round(dropsFW$envEffs[,2]), sec.axis = dup_axis(labels = locations))

ggsave(filename = paste0("./", pipeline, '_stability_adaptability_across_locations.tiff'),
       plot = q,
       device = 'tiff',
       width = 300,
       height = 400,
       units = 'mm',
       dpi = 300)

# GGE Biplot
library(metan)
model <- gge(data.MET, DOP, Variety, Yield, svp = "symmetrical")
a <- plot(model, type = 1)
b <- plot(model, type = 2)
c <- plot(model, type = 3)
d <- arrange_ggplot(a, b, c, tag_levels = "a")

ggsave(filename = paste0("./", pipeline, '_GGE_Biplots.tiff'),
       plot = d,
       device = 'tiff',
       width = 400,
       height = 250,
       units = 'mm',
       dpi = 300)

## barplot grouped by DOP
q <- ggplot(data = data.MET, aes(y = Yield, 
                                 x = DOP, 
                                 fill = Variety)) +
  geom_bar(position="dodge", stat="identity") +
  labs(x = "DOP", y = "Yield") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  coord_cartesian(ylim=c(min(data.MET$Yield - data.MET$standard.error*1.96),
                         max(data.MET$Yield + data.MET$standard.error*1.96)))


ggsave(filename = paste0("./", pipeline, '_performances _by_DOP.tiff'),
       plot = q,
       device = 'tiff',
       width = 400,
       height = 250,
       units = 'mm',
       dpi = 300)

###### the end ###########