###################################################################################
#COMBINED SITES SCALE 0 to 3
#################################################################################

scabdata <- read.csv("combinedimputed03.csv", header = TRUE)
View(scabdata)
library(agricolae)
summary(scabdata)
dim(scabdata)
scabscore <- scabdata[,5:7]
scabscore
dim(scabscore)
#############################################
#AREA UNDER THE CURVE DISEASE PROGRESSION
#############################################
library(epifitter)
dates <- c(33, 39, 49)
evaluation <- c(scaDISEASE.SCORE.3rd.JUNE = 33, DISEASE_SCORE.9TH.JUNE = 66, DISEASE_SCORE19THJUNE = 99)
AUDPC(evaluation,dates, y_proportion = TRUE, type = "absolute")
plot(dates,evaluation,type="h",ylim=c(0,100),col="red",axes=FALSE)
title(cex.main=0.8,main="Absolute or Relative AUDPC\nTotal area = 2640*(58-44)=26400")
2640*(49-39)
lines(dates,evaluation,col="red")
text(dates,evaluation+5,evaluation)
#text(18,20,"A = (44-30)*(2+1)/2")
(30-44)*(60+30)/2
#text(25,60,"B = (58-44)*(3+2)/2")
(58-44)*(90+60)/2
#text(25,40,"audpc = A+B = 1015")
#text(24.5,33,"relative = audpc/area = 0.725")
abline(h=0)
axis(1,dates)
axis(2,seq(0,100,30),las=2)
lines(rbind(c(14,40),c(14,100)),lty=8,col="green")
lines(rbind(c(14,100),c(28,100)),lty=8,col="green")
lines(rbind(c(28,90),c(28,100)),lty=8,col="green")

h1<-graph.freq(AUDPC,border="red",density=4,col="blue")
table.freq(h1)
############################################
####END AREA UNDER THE CURVE
###########################################
scabgeomeanscombined <- exp(rowMeans(log(scabscore)))
as.numeric(scabgeomeanscombined)
trtc <- scabdata[,1]
as.factor(trtc)
repsc <- scabdata[,4]
as.factor(repsc)
sitec <- scabdata[,3]
as.factor(sitec)
blockc <- scabdata[,8]
as.factor(blockc)
stand_count <- scabdata[,2]
datac <- cbind(trtc, sitec, repsc, blockc, scabgeomeanscombined, stand_count)
as.data.frame(datac)
str(datac)
analysisc <- lm(scabgeomeanscombined ~ trtc+repsc+sitec, scabdata)
anova(analysisc)
out3 <- LSD.test(analysisc, "trtc", p.adj="bonferroni")
out3
out3$groups$scabgeomeans
library(writexl)
combinedphenotype <- out3$groups
shapiro.test(combinedphenotype$scabgeomeanscombined)
library("dplyr")
library("ggpubr")
library("forecast")
ggqqplot(combinedphenotype$scabgeomeanscombined)
as.data.frame(combinedphenotype)
Resistants <- combinedphenotype[c(162,166:171,173:179),]
Susceptible <- combinedphenotype[c(1:161,163,164,165),]
Resistants
Susceptible
t.test(Resistants[,1], Susceptible[,1])
boxplot(Resistants[,1], Susceptible[,1], names = c("Resistant", "Susceptible"))

write.table(datac, file = "datac.txt" , append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
write.csv(datac, file = "datac.csv")
write.csv(combinedphenotype, file = "PHcombinedAnalysed.csv")          
write.xlsx(combinedphenotype, file = "combinedAnalysed000003.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE, showNA = TRUE,password = NULL)
hist(combinedphenotype$scabgeomeanscombined, xlab = "scab scores", 
     main = "SCAB MEAN DISTRIBUTION", col = "red"
     )
boxplot(combinedphenotype$scabgeomeanscombined, ylab = "mean scores", 
        main = "SCAB MEAN DISTRIBUTION"
        )

######################################################
##########HERITABILITY
######################################################
library(variability)
View(datac)
as.data.frame.array(datac)
gen.var(data = datac[5:6], datac[1], datac[3])





library(biganalytics)
library(bigmemory)
library(sommer)

#G <- read.big.matrix("SNPDATACOMPLETE.txt", type = "char," sep = "\t")
G <- read.csv("SNPDATACOMPLETE.csv", header = TRUE)
Y <- read.csv("datac.csv", header = TRUE)

ans1 <- mmer(scabgeomeanscombined~sitec,
             random= ~ trtc + sitec:trtc,
             rcov= ~ units,
             data=Y)
summary(ans1)



# <- as.data.frame(unclass(datac),
 #                  stringsAsFactors = TRUE)

Y
str(Y)
as.factor(Y$trtc)
as.factor(Y$sitec)
as.factor(Y$repsc)
str(Y)
#########################################


#######################################
library(inti)
library(agridat)

dt <- john.alpha
str(dt)

hr<-H2cal(data = dt
          , trait = "yield"
          , gen.name = "gen"
          , rep.n = 3
          , ran.model = "1 + rep + (1|rep:block) + (1|gen)"
          , fix.model = "0 + rep + (1|rep:block) + gen"
          , emmeans = TRUE
          , plot_diag = TRUE
          , outliers.rm = TRUE
          )

hr###################################


dt <- Y
dt
str(Y)
Y$scabgeomeanscombined <- as.numeric(as.character(Y$scabgeomeanscombined))
hr<-H2cal(data = dt
          , trait = "scabgeomeanscombined"
          , gen.name = "trtc"
          , rep.n = 3
          , ran.model = "1 + repsc + (1|repsc:blockc) + (1|trtc)"
          , fix.model = "0 + repsc + (1|repsc:blockc) + trtc"
          , emmeans = TRUE
          , plot_diag = TRUE
          , outliers.rm = TRUE
)
#hr$tabsmr %>% kable(caption = "Variance Component table")
hr$model %>% summary()
hr$tabsmr
install.packages("kableExtra")
library(kableExtra)
hr$tabsmr %>% kable(caption = "Variance Component table")
hr$blups

str(Y)
cols <- c("X", "trtc", "sitec", "repsc", "blockc")
Y[cols] <- lapply(Y[cols], factor)
str(Y)
str(datac)
#as.factor(Y$trtc)
as.factor(Y$sitec)
as.factor(Y$repsc)
as.numeric(Y$scabgeomeanscombined)

str(datac)
#library(DataCombine)
#FindReplace(Y, Y$sitec, KK.2021, from = "from", to = "to", exact = TRUE, vector = FALSE)
#replace(Y, Y$sitec==1, "KK.2021") 
#Y$sitec[Y$sitec == 1] <- "kk.2021" 
#Y$sitec[Y$sitec == 2] <- "bu.2021"
#Y$repsc[Y$repsc[1:178] ==1] <- "kk.2021.1"
#Y$repsc[Y$repsc[179:355] ==2] <- "kk.2021.2"
#Y$repsc[Y$repsc[356:505] ==3] <- "kk.2021.3"
#Y$repsc[Y$repsc[506:683] ==1] <- "bu.2021.1"
#Y$repsc[Y$repsc[684:857] ==2] <- "bu.2021.2"
#Y$repsc[Y$repsc[858:1031] ==3] <- "bu.2021.3"
#getOption(max.print(Y$repsc))
#Y$sitec[Y$repsc == 2] <- "bu.2021.2"
library(biganalytics)
library(bigmemory)
library(sommer)

str(Y)
Y$scabgeomeanscombined <- as.numeric(as.character(Y$scabgeomeanscombined))
ans1 <- mmer(scabgeomeanscombined~1, 
             random = trtc + sitec + trtc:sitec + sitec:repsc, 
             rcov= ~ units, 
             data = Y, 
             verbose = FALSE)
install.packages('sommer')
ans1

data("DT_example")
DT <- DT_example
DT
str(DT)
str(Y)
