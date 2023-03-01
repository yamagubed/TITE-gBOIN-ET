
DA1.b <- read.csv("../../data/simulation1_dose_allocation_gBOIN-ET.csv")
DA1.t <- read.csv("../../data/simulation1_dose_allocation_TITE-gBOIN-ET.csv")
OD1.b <- read.csv("../../data/simulation1_obd_gBOIN-ET.csv")
OD1.t <- read.csv("../../data/simulation1_obd_TITE-gBOIN-ET.csv")

DA2.b <- read.csv("../../data/simulation2_dose_allocation_gBOIN-ET.csv")
DA2.t <- read.csv("../../data/simulation2_dose_allocation_TITE-gBOIN-ET.csv")
OD2.b <- read.csv("../../data/simulation2_obd_gBOIN-ET.csv")
OD2.t <- read.csv("../../data/simulation2_obd_TITE-gBOIN-ET.csv")

# Table 1 #

colvar <- c("Z4","M","X1","X2","X3","X4","X5","X6")
WK1 <- OD1.b[(OD1.b$Z1==60)&(OD1.b$Z2==10)&(OD1.b$Z3=="Highest Efficacy"),colvar]
WK2 <- OD1.t[(OD1.t$Z1==60)&(OD1.t$Z2==10)&(OD1.t$Z3=="Highest Efficacy"),colvar]

Table <- NULL
for(s1 in 1:10){
 WK3  <- WK1[WK1$Z4==s1,]
 WK11 <- data.frame(
           Scenario = WK3$Z4,
           Design   = paste("   ",WK3$M),
           PCS      = paste("   ",format(round(WK3$X1,digits=1),nsmall=1)),
           PCA      = paste("   ",format(round(WK3$X2,digits=1),nsmall=1)),
           POS      = paste("   ",format(round(WK3$X3,digits=1),nsmall=1)),
           POA      = paste("   ",format(round(WK3$X4,digits=1),nsmall=1)),
           PET      = paste("   ",format(round(WK3$X5,digits=1),nsmall=1)),
           Dur      = paste("   ",format(round(WK3$X6,digits=0),nsmall=0)))

 WK3  <- WK2[WK2$Z4==s1,]
 WK21 <- data.frame(
           Scenario = WK3$Z4,
           Design   = paste("   ",WK3$M),
           PCS      = paste("   ",format(round(WK3$X1,digits=1),nsmall=1)),
           PCA      = paste("   ",format(round(WK3$X2,digits=1),nsmall=1)),
           POS      = paste("   ",format(round(WK3$X3,digits=1),nsmall=1)),
           POA      = paste("   ",format(round(WK3$X4,digits=1),nsmall=1)),
           PET      = paste("   ",format(round(WK3$X5,digits=1),nsmall=1)),
           Dur      = paste("   ",format(round(WK3$X6,digits=0),nsmall=0)))

 Table <- rbind(Table,WK11,WK21)
}

sink("../../output/Table 1.txt")
Table
sink()

# Table S2 #

colvar <- c("Z4","M","Z5","X1","X2","X3","X4","X5","X6","X7")
WK1 <- DA1.b[(DA1.b$Z1==60)&(DA1.b$Z2==10)&(DA1.b$Z3=="Highest Efficacy"),colvar]
WK2 <- DA1.t[(DA1.t$Z1==60)&(DA1.t$Z2==10)&(DA1.t$Z3=="Highest Efficacy"),colvar]

Table <- NULL
for(s1 in 1:10){
 WK3  <- WK1[WK1$Z4==s1,]
 WK11 <- data.frame(
           Scenario  = WK3$Z4,
           Design    = paste("   ",WK3$M),
           Parameter = paste("   ",WK3$Z5),
           ET        = c(paste("   ",format(round(WK3$X1[1],digits=0),nsmall=0)),paste("   "," ")),
           Dose1     = c(paste("   ",format(round(WK3$X2[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X2[2],digits=1),nsmall=1))),
           Dose2     = c(paste("   ",format(round(WK3$X3[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X3[2],digits=1),nsmall=1))),
           Dose3     = c(paste("   ",format(round(WK3$X4[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X4[2],digits=1),nsmall=1))),
           Dose4     = c(paste("   ",format(round(WK3$X5[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X5[2],digits=1),nsmall=1))),
           Dose5     = c(paste("   ",format(round(WK3$X6[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X6[2],digits=1),nsmall=1))),
           Dose6     = c(paste("   ",format(round(WK3$X7[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X7[2],digits=1),nsmall=1))))

 WK3  <- WK2[WK2$Z4==s1,]
 WK21 <- data.frame(
           Scenario  = WK3$Z4,
           Design    = paste("   ",WK3$M),
           Parameter = paste("   ",WK3$Z5),
           ET        = c(paste("   ",format(round(WK3$X1[1],digits=0),nsmall=0)),paste("   "," ")),
           Dose1     = c(paste("   ",format(round(WK3$X2[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X2[2],digits=1),nsmall=1))),
           Dose2     = c(paste("   ",format(round(WK3$X3[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X3[2],digits=1),nsmall=1))),
           Dose3     = c(paste("   ",format(round(WK3$X4[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X4[2],digits=1),nsmall=1))),
           Dose4     = c(paste("   ",format(round(WK3$X5[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X5[2],digits=1),nsmall=1))),
           Dose5     = c(paste("   ",format(round(WK3$X6[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X6[2],digits=1),nsmall=1))),
           Dose6     = c(paste("   ",format(round(WK3$X7[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X7[2],digits=1),nsmall=1))))

 Table <- rbind(Table,WK11,WK21)
}

sink("../../output/Table S2.txt")
Table
sink()

# Table S3 #

colvar <- c("Z4","M","X1","X2","X3","X4","X5","X6")
WK1 <- OD1.b[(OD1.b$Z1==60)&(OD1.b$Z2==5)&(OD1.b$Z3=="Highest Efficacy"),colvar]
WK2 <- OD1.t[(OD1.t$Z1==60)&(OD1.t$Z2==5)&(OD1.t$Z3=="Highest Efficacy"),colvar]

Table <- NULL
for(s1 in 1:10){
 WK3  <- WK1[WK1$Z4==s1,]
 WK11 <- data.frame(
           Scenario = WK3$Z4,
           Design   = paste("   ",WK3$M),
           PCS      = paste("   ",format(round(WK3$X1,digits=1),nsmall=1)),
           PCA      = paste("   ",format(round(WK3$X2,digits=1),nsmall=1)),
           POS      = paste("   ",format(round(WK3$X3,digits=1),nsmall=1)),
           POA      = paste("   ",format(round(WK3$X4,digits=1),nsmall=1)),
           PET      = paste("   ",format(round(WK3$X5,digits=1),nsmall=1)),
           Dur      = paste("   ",format(round(WK3$X6,digits=0),nsmall=0)))

 WK3  <- WK2[WK2$Z4==s1,]
 WK21 <- data.frame(
           Scenario = WK3$Z4,
           Design   = paste("   ",WK3$M),
           PCS      = paste("   ",format(round(WK3$X1,digits=1),nsmall=1)),
           PCA      = paste("   ",format(round(WK3$X2,digits=1),nsmall=1)),
           POS      = paste("   ",format(round(WK3$X3,digits=1),nsmall=1)),
           POA      = paste("   ",format(round(WK3$X4,digits=1),nsmall=1)),
           PET      = paste("   ",format(round(WK3$X5,digits=1),nsmall=1)),
           Dur      = paste("   ",format(round(WK3$X6,digits=0),nsmall=0)))

 Table <- rbind(Table,WK11,WK21)
}

sink("../../output/Table S3.txt")
Table
sink()

# Table S4 #

colvar <- c("Z4","M","Z5","X1","X2","X3","X4","X5","X6","X7")
WK1 <- DA1.b[(DA1.b$Z1==60)&(DA1.b$Z2==5)&(DA1.b$Z3=="Highest Efficacy"),colvar]
WK2 <- DA1.t[(DA1.t$Z1==60)&(DA1.t$Z2==5)&(DA1.t$Z3=="Highest Efficacy"),colvar]

Table <- NULL
for(s1 in 1:10){
 WK3  <- WK1[WK1$Z4==s1,]
 WK11 <- data.frame(
           Scenario  = WK3$Z4,
           Design    = paste("   ",WK3$M),
           Parameter = paste("   ",WK3$Z5),
           ET        = c(paste("   ",format(round(WK3$X1[1],digits=0),nsmall=0)),paste("   "," ")),
           Dose1     = c(paste("   ",format(round(WK3$X2[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X2[2],digits=1),nsmall=1))),
           Dose2     = c(paste("   ",format(round(WK3$X3[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X3[2],digits=1),nsmall=1))),
           Dose3     = c(paste("   ",format(round(WK3$X4[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X4[2],digits=1),nsmall=1))),
           Dose4     = c(paste("   ",format(round(WK3$X5[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X5[2],digits=1),nsmall=1))),
           Dose5     = c(paste("   ",format(round(WK3$X6[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X6[2],digits=1),nsmall=1))),
           Dose6     = c(paste("   ",format(round(WK3$X7[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X7[2],digits=1),nsmall=1))))

 WK3  <- WK2[WK2$Z4==s1,]
 WK21 <- data.frame(
           Scenario  = WK3$Z4,
           Design    = paste("   ",WK3$M),
           Parameter = paste("   ",WK3$Z5),
           ET        = c(paste("   ",format(round(WK3$X1[1],digits=0),nsmall=0)),paste("   "," ")),
           Dose1     = c(paste("   ",format(round(WK3$X2[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X2[2],digits=1),nsmall=1))),
           Dose2     = c(paste("   ",format(round(WK3$X3[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X3[2],digits=1),nsmall=1))),
           Dose3     = c(paste("   ",format(round(WK3$X4[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X4[2],digits=1),nsmall=1))),
           Dose4     = c(paste("   ",format(round(WK3$X5[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X5[2],digits=1),nsmall=1))),
           Dose5     = c(paste("   ",format(round(WK3$X6[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X6[2],digits=1),nsmall=1))),
           Dose6     = c(paste("   ",format(round(WK3$X7[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X7[2],digits=1),nsmall=1))))

 Table <- rbind(Table,WK11,WK21)
}

sink("../../output/Table S4.txt")
Table
sink()

# Table S5 #

colvar <- c("Z4","M","X1","X2","X3","X4","X5","X6")
WK1 <- OD1.b[(OD1.b$Z1==90)&(OD1.b$Z2==10)&(OD1.b$Z3=="Highest Efficacy"),colvar]
WK2 <- OD1.t[(OD1.t$Z1==90)&(OD1.t$Z2==10)&(OD1.t$Z3=="Highest Efficacy"),colvar]

Table <- NULL
for(s1 in 1:10){
 WK3  <- WK1[WK1$Z4==s1,]
 WK11 <- data.frame(
           Scenario = WK3$Z4,
           Design   = paste("   ",WK3$M),
           PCS      = paste("   ",format(round(WK3$X1,digits=1),nsmall=1)),
           PCA      = paste("   ",format(round(WK3$X2,digits=1),nsmall=1)),
           POS      = paste("   ",format(round(WK3$X3,digits=1),nsmall=1)),
           POA      = paste("   ",format(round(WK3$X4,digits=1),nsmall=1)),
           PET      = paste("   ",format(round(WK3$X5,digits=1),nsmall=1)),
           Dur      = paste("   ",format(round(WK3$X6,digits=0),nsmall=0)))

 WK3  <- WK2[WK2$Z4==s1,]
 WK21 <- data.frame(
           Scenario = WK3$Z4,
           Design   = paste("   ",WK3$M),
           PCS      = paste("   ",format(round(WK3$X1,digits=1),nsmall=1)),
           PCA      = paste("   ",format(round(WK3$X2,digits=1),nsmall=1)),
           POS      = paste("   ",format(round(WK3$X3,digits=1),nsmall=1)),
           POA      = paste("   ",format(round(WK3$X4,digits=1),nsmall=1)),
           PET      = paste("   ",format(round(WK3$X5,digits=1),nsmall=1)),
           Dur      = paste("   ",format(round(WK3$X6,digits=0),nsmall=0)))

 Table <- rbind(Table,WK11,WK21)
}

sink("../../output/Table S5.txt")
Table
sink()

# Table S6 #

colvar <- c("Z4","M","Z5","X1","X2","X3","X4","X5","X6","X7")
WK1 <- DA1.b[(DA1.b$Z1==90)&(DA1.b$Z2==10)&(DA1.b$Z3=="Highest Efficacy"),colvar]
WK2 <- DA1.t[(DA1.t$Z1==90)&(DA1.t$Z2==10)&(DA1.t$Z3=="Highest Efficacy"),colvar]

Table <- NULL
for(s1 in 1:10){
 WK3  <- WK1[WK1$Z4==s1,]
 WK11 <- data.frame(
           Scenario  = WK3$Z4,
           Design    = paste("   ",WK3$M),
           Parameter = paste("   ",WK3$Z5),
           ET        = c(paste("   ",format(round(WK3$X1[1],digits=0),nsmall=0)),paste("   "," ")),
           Dose1     = c(paste("   ",format(round(WK3$X2[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X2[2],digits=1),nsmall=1))),
           Dose2     = c(paste("   ",format(round(WK3$X3[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X3[2],digits=1),nsmall=1))),
           Dose3     = c(paste("   ",format(round(WK3$X4[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X4[2],digits=1),nsmall=1))),
           Dose4     = c(paste("   ",format(round(WK3$X5[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X5[2],digits=1),nsmall=1))),
           Dose5     = c(paste("   ",format(round(WK3$X6[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X6[2],digits=1),nsmall=1))),
           Dose6     = c(paste("   ",format(round(WK3$X7[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X7[2],digits=1),nsmall=1))))

 WK3  <- WK2[WK2$Z4==s1,]
 WK21 <- data.frame(
           Scenario  = WK3$Z4,
           Design    = paste("   ",WK3$M),
           Parameter = paste("   ",WK3$Z5),
           ET        = c(paste("   ",format(round(WK3$X1[1],digits=0),nsmall=0)),paste("   "," ")),
           Dose1     = c(paste("   ",format(round(WK3$X2[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X2[2],digits=1),nsmall=1))),
           Dose2     = c(paste("   ",format(round(WK3$X3[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X3[2],digits=1),nsmall=1))),
           Dose3     = c(paste("   ",format(round(WK3$X4[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X4[2],digits=1),nsmall=1))),
           Dose4     = c(paste("   ",format(round(WK3$X5[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X5[2],digits=1),nsmall=1))),
           Dose5     = c(paste("   ",format(round(WK3$X6[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X6[2],digits=1),nsmall=1))),
           Dose6     = c(paste("   ",format(round(WK3$X7[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X7[2],digits=1),nsmall=1))))

 Table <- rbind(Table,WK11,WK21)
}

sink("../../output/Table S6.txt")
Table
sink()

# Table 2 #

colvar <- c("Z4","M","X1","X2","X3","X4","X5","X6")
WK1 <- OD2.b[(OD2.b$Z1==60)&(OD2.b$Z2==10)&(OD2.b$Z3=="Highest Efficacy"),colvar]
WK2 <- OD2.t[(OD2.t$Z1==60)&(OD2.t$Z2==10)&(OD2.t$Z3=="Highest Efficacy"),colvar]

Table <- NULL
for(s1 in 1:10){
 WK3  <- WK1[WK1$Z4==s1,]
 WK11 <- data.frame(
           Scenario = WK3$Z4,
           Design   = paste("   ",WK3$M),
           PCS      = paste("   ",format(round(WK3$X1,digits=1),nsmall=1)),
           PCA      = paste("   ",format(round(WK3$X2,digits=1),nsmall=1)),
           POS      = paste("   ",format(round(WK3$X3,digits=1),nsmall=1)),
           POA      = paste("   ",format(round(WK3$X4,digits=1),nsmall=1)),
           PET      = paste("   ",format(round(WK3$X5,digits=1),nsmall=1)),
           Dur      = paste("   ",format(round(WK3$X6,digits=0),nsmall=0)))

 WK3  <- WK2[WK2$Z4==s1,]
 WK21 <- data.frame(
           Scenario = WK3$Z4,
           Design   = paste("   ",WK3$M),
           PCS      = paste("   ",format(round(WK3$X1,digits=1),nsmall=1)),
           PCA      = paste("   ",format(round(WK3$X2,digits=1),nsmall=1)),
           POS      = paste("   ",format(round(WK3$X3,digits=1),nsmall=1)),
           POA      = paste("   ",format(round(WK3$X4,digits=1),nsmall=1)),
           PET      = paste("   ",format(round(WK3$X5,digits=1),nsmall=1)),
           Dur      = paste("   ",format(round(WK3$X6,digits=0),nsmall=0)))

 Table <- rbind(Table,WK11,WK21)
}

sink("../../output/Table 2.txt")
Table
sink()

# Table S7 #

colvar <- c("Z4","M","Z5","X1","X2","X3","X4","X5","X6","X7")
WK1 <- DA2.b[(DA2.b$Z1==60)&(DA2.b$Z2==10)&(DA2.b$Z3=="Highest Efficacy"),colvar]
WK2 <- DA2.t[(DA2.t$Z1==60)&(DA2.t$Z2==10)&(DA2.t$Z3=="Highest Efficacy"),colvar]

Table <- NULL
for(s1 in 1:10){
 WK3  <- WK1[WK1$Z4==s1,]
 WK11 <- data.frame(
           Scenario  = WK3$Z4,
           Design    = paste("   ",WK3$M),
           Parameter = paste("   ",WK3$Z5),
           ET        = c(paste("   ",format(round(WK3$X1[1],digits=0),nsmall=0)),paste("   "," ")),
           Dose1     = c(paste("   ",format(round(WK3$X2[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X2[2],digits=1),nsmall=1))),
           Dose2     = c(paste("   ",format(round(WK3$X3[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X3[2],digits=1),nsmall=1))),
           Dose3     = c(paste("   ",format(round(WK3$X4[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X4[2],digits=1),nsmall=1))),
           Dose4     = c(paste("   ",format(round(WK3$X5[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X5[2],digits=1),nsmall=1))),
           Dose5     = c(paste("   ",format(round(WK3$X6[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X6[2],digits=1),nsmall=1))),
           Dose6     = c(paste("   ",format(round(WK3$X7[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X7[2],digits=1),nsmall=1))))

 WK3  <- WK2[WK2$Z4==s1,]
 WK21 <- data.frame(
           Scenario  = WK3$Z4,
           Design    = paste("   ",WK3$M),
           Parameter = paste("   ",WK3$Z5),
           ET        = c(paste("   ",format(round(WK3$X1[1],digits=0),nsmall=0)),paste("   "," ")),
           Dose1     = c(paste("   ",format(round(WK3$X2[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X2[2],digits=1),nsmall=1))),
           Dose2     = c(paste("   ",format(round(WK3$X3[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X3[2],digits=1),nsmall=1))),
           Dose3     = c(paste("   ",format(round(WK3$X4[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X4[2],digits=1),nsmall=1))),
           Dose4     = c(paste("   ",format(round(WK3$X5[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X5[2],digits=1),nsmall=1))),
           Dose5     = c(paste("   ",format(round(WK3$X6[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X6[2],digits=1),nsmall=1))),
           Dose6     = c(paste("   ",format(round(WK3$X7[1],digits=0),nsmall=0)),paste("   ",format(round(WK3$X7[2],digits=1),nsmall=1))))

 Table <- rbind(Table,WK11,WK21)
}

sink("../../output/Table S7.txt")
Table
sink()

