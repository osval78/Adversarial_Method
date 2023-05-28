rm(list=ls())
library(BGLR) #load BLR
library(dplyr)
library(caret)
library(plyr)
library(tidyr)
library(SKM)
library(STPGA)
library(reshape2)
#####Loading data set#############
load("Groundnut.RData")
ls()
Pheno=Pheno
Geno=Geno
head(Pheno)
Traits_to_evaluate=colnames(Pheno)[c(2:5)]
#Traits_to_evaluate
ToP_Selection=20
Threshold_Opt_trn=0.8
# Variables definitions ---------------------------------------------------------
#######Model to be implemented
model_name <- "BGBLUP"
Name_data_set="Groundnut"
Version=1

########Folds and iterations to use 
cv_folds_number <- 5
tuning_folds_number <-10
iterations_number <- 8000
burn_in <- 2500

# Data preparation -------------------------------------------------------------
Pheno <- Pheno %>% arrange(Line)

# Select the unique lines both in pheno and geno and sort them
final_geno_lines <- intersect(Pheno$Line, rownames(Geno)) %>% sort()
Geno <- Geno[final_geno_lines, final_geno_lines]
dim(Geno)


Pheno$Env="EYT_3"
Pheno_E=Pheno
Predictions_All_traits_Fold=data.frame()
Summary_AUC_Fold=data.frame()
Summary_Accuracy_Testing_Fold=data.frame()

for (trait in Traits_to_evaluate) {
  #trait=Traits_to_evaluate[1]
  SKM::echo("*** Trait: %s ***", trait)
  BLUEs_trait=Pheno_E[,trait]
  Q.20=quantile(BLUEs_trait,prob=(1-ToP_Selection/100))
  
  Predictions_All_traits=data.frame()
  Summary_AUC=data.frame()
  Summary_Accuracy_Testing=data.frame()
  
  for (r in 1:5){
    
    set.seed(r)
    All_pos=1:length(BLUEs_trait)
    size_TST=which(BLUEs_trait>Q.20)
    Tst_Set=sample(All_pos,length(size_TST))
    TRN_Set=All_pos[-Tst_Set]
    
    y_Bin=rep(0,nrow(Pheno_E))
    
    y_Bin[-TRN_Set]=1
    Bin_name=paste("Bin",trait, r,sep="_")
    Pheno_E=cbind(Pheno_E,Y_Bin=y_Bin)
    #head(Pheno_E)
    colnames(Pheno_E)=c(colnames(Pheno_E)[-ncol(Pheno_E)],Bin_name[1])
    #head(Pheno_E)
    # Cross validation -----------------------------------------------------------
    
    Predictions=data.frame()
    
    PhenoTuning <- Pheno_E
    y_tuning_Bin <- PhenoTuning %>% pull(Bin_name)
    tuning_lines <- as.character(PhenoTuning$Line)
    GenoTuning <- Geno[tuning_lines,tuning_lines]
    GenoTuning=data.matrix(GenoTuning)
    y_tuning_Bin_Factor=as.factor(y_tuning_Bin)
    ZL=model.matrix(~0+Line,data=PhenoTuning)
    GRM=ZL%*%GenoTuning%*%t(ZL)
    ETATuning1=list(Line=list(model='RKHS',K=GRM))
    
    inner_folds <- SKM::cv_kfold_strata(y_tuning_Bin_Factor, k = tuning_folds_number)
    
    Predictions = data.frame(Trait=trait,
                             Line = Pheno_E$Line,
                             Fold = NA,
                             Observed = y_tuning_Bin,
                             Predicted = NA,
                             prob=NA,
                             Opt_trn=NA)
    
    for (k in seq_along(inner_folds)) {
      #         k=1
      SKM::echo("\t\t\t*** InnerFold: %s / %s ***", k, length(inner_folds))
      inner_fold <- inner_folds[[k]]
      
      y_na <-as.factor(y_tuning_Bin)
      y_na[inner_fold$testing] <- NA
      
      model <- BGLR::BGLR(
        y = y_na,
        ETA = ETATuning1,
        response_type = "ordinal",
        nIter = iterations_number,
        burnIn = burn_in,
        verbose = FALSE
      )
      
      prob_val <- model$probs[inner_fold$testing,2]
      Observed <-y_tuning_Bin[inner_fold$testing]
      Numbers_One=length(which(y_tuning_Bin==1))
      Threshold=Numbers_One/length(y_tuning_Bin)
      Predicted=ifelse(prob_val>Threshold,1,0)
      Threshold_Opt_trn_adj=(Threshold_Opt_trn*Threshold)/0.5
      Opt_trn2=ifelse(prob_val>Threshold_Opt_trn_adj,1,0)
      Opt_trn=(1-Observed)*Opt_trn2
      
      Predictions$Fold[inner_fold$testing] = k
      Predictions$Predicted[inner_fold$testing] = Predicted
      Predictions$prob[inner_fold$testing] = prob_val
      Predictions$Opt_trn[inner_fold$testing] = Opt_trn      
    }
    
    ##########Predictions of testing with full data
    yy=BLUEs_trait
    yy_na=yy
    yy_na[-TRN_Set]=NA
    model1 <- BGLR::BGLR(
      y = yy_na,
      ETA = ETATuning1,
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE
    )
    
    Predicted_trait=model1$yHat
    Predictions$Predicted_trait=Predicted_trait
    Predictions$Observed_trait=yy
    MSE=mse(yy[-TRN_Set],Predicted_trait[-TRN_Set])
    COR=cor(yy[-TRN_Set],Predicted_trait[-TRN_Set])
    NRMSE=nrmse(yy[-TRN_Set],Predicted_trait[-TRN_Set])
    
    ##########Predictions of testing with sample of full data
    Opt_trn=which(Predictions$Opt_trn==1)
    
    All_sample=c(Opt_trn,Tst_Set)
    PhenoTuning2=Pheno_E[All_sample,]
    PhenoTuning2=droplevels(PhenoTuning2)
    ZL2=model.matrix(~0+Line,data=PhenoTuning2)
    GID_Samples=unique( PhenoTuning2$Line)
    GenoTuning2=Geno[ GID_Samples, GID_Samples]
    GRM2=ZL2%*%GenoTuning2%*%t(ZL2)
    ETATuning2=list(Line=list(model='RKHS',K=GRM2))
    yy2=PhenoTuning2[,trait]
    yy_na2=yy2
    TRN_Set2=1:length(Opt_trn)
    yy_na2[-TRN_Set2]=NA
    
    model2 <- BGLR::BGLR(
      y = yy_na2,
      ETA = ETATuning2,
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE
    )
    
    Predicted_trait2=model2$yHat
    MSE2=mse(yy2[-TRN_Set2],Predicted_trait2[-TRN_Set2])
    COR2=cor(yy2[-TRN_Set2],Predicted_trait2[-TRN_Set2])
    NRMSE2=nrmse(yy2[-TRN_Set2],Predicted_trait2[-TRN_Set2])
    
    ##########CDMEAN and PEVMEAN methods
    n=length(yy)
    Perc1=0.5
    K= GRM
    svd_G=svd(K,nu=round(n*Perc1,0),nv=round(n*Perc1,0))
    str(svd_G)
    PC5_Wheat=K%*%svd_G$v
    #plot(PC5_Wheat[,1],PC5_Wheat[,2])
    
    
    rownames(PC5_Wheat)=1:n
    Candidates=TRN_Set
    NTrn=n*0.6
    
    ListTrain1<-GenAlgForSubsetSelection(P=PC5_Wheat,Candidates=Candidates,Test=Tst_Set,ntoselect=NTrn, InitPop=NULL,
                                         npop=20, nelite=10, mutprob=.5, mutintensity = 1,
                                         niterations=5000,minitbefstop=10, tabu=F,tabumemsize = 0,plotiters=F,
                                         lambda=1e-5,errorstat="CDMEAN", mc.cores=3)
    
    ListTrain2<-GenAlgForSubsetSelection(PC5_Wheat,Candidates=Candidates,Test=Tst_Set,ntoselect=NTrn, InitPop=NULL,
                                         npop=20, nelite=10, mutprob=.5, mutintensity = 1,
                                         niterations=5000,minitbefstop=10, tabu=F,tabumemsize = 0,plotiters=F,
                                         lambda=1e-5,errorstat="PEVMEAN", mc.cores=3)
    
    
    Opt_trn_A=as.vector(ListTrain1[[1]]);
    
    All_sample3=c(Opt_trn_A,Tst_Set)
    PhenoTuning3=Pheno_E[All_sample3,]
    PhenoTuning3=droplevels(PhenoTuning3)
    ZL3=model.matrix(~0+Line ,data=PhenoTuning3)
    GID_Samples3=unique( PhenoTuning3$Line)
    GenoTuning3=Geno[ GID_Samples3, GID_Samples3]
    GRM3=ZL3%*%GenoTuning3%*%t(ZL3)
    
    ETATuning3=list(Line=list(model='RKHS',K=GRM3))
    
    yy3=PhenoTuning3[,trait]
    yy_na3=yy3
    TRN_Set3=1:length(Opt_trn_A)
    yy_na3[-TRN_Set3]=NA
    model3 <- BGLR::BGLR(
      y = yy_na3,
      ETA = ETATuning3,
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE
    )
    Predicted_trait3=model3$yHat
    #Predictions$Predicted_trait=Predicted_trait
    #Predictions$Observed_trait=yy
    MSE3=mse(yy3[-TRN_Set3],Predicted_trait3[-TRN_Set3])
    COR3=cor(yy3[-TRN_Set3],Predicted_trait3[-TRN_Set3])
    NRMSE3=nrmse(yy3[-TRN_Set3],Predicted_trait3[-TRN_Set3])
    
    Opt_trn4=as.vector(ListTrain2[[1]]);
    
    All_sample4=c(Opt_trn4,Tst_Set)
    PhenoTuning4=Pheno_E[All_sample4,]
    PhenoTuning4=droplevels(PhenoTuning4)
    ZL4=model.matrix(~0+Line ,data=PhenoTuning4)
    GID_Samples4=unique( PhenoTuning4$Line)
    GenoTuning4=Geno[ GID_Samples4, GID_Samples4]
    GRM4=ZL4%*%GenoTuning4%*%t(ZL4)
    
    ETATuning4=list(Line=list(model='RKHS',K=GRM4))
    
    yy4=PhenoTuning4[,trait]
    yy_na4=yy4
    TRN_Set4=1:length(Opt_trn4)
    yy_na4[-TRN_Set4]=NA
    model4 <- BGLR::BGLR(
      y = yy_na4,
      ETA = ETATuning4,
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE
    )
    Predicted_trait4=model4$yHat
    #Predictions$Predicted_trait=Predicted_trait
    #Predictions$Observed_trait=yy
    MSE4=mse(yy4[-TRN_Set4],Predicted_trait4[-TRN_Set4])
    COR4=cor(yy4[-TRN_Set4],Predicted_trait4[-TRN_Set4])
    NRMSE4=nrmse(yy4[-TRN_Set4],Predicted_trait4[-TRN_Set4])
    
    ########Average and Max distance
    K_trn_tst=(GRM[TRN_Set,-TRN_Set])
    K_trn_tst_Max=apply(K_trn_tst,1,max)
    K_trn_tst_Mean=apply(K_trn_tst,1,mean)
    Q_max=quantile(K_trn_tst_Max,prob=0.75)  
    Q_mean=quantile(K_trn_tst_Mean,prob=0.75)    
    Max_pos=which(K_trn_tst_Max>Q_max)
    Trn_Set_Max=TRN_Set[Max_pos]
    Mean_pos=which(K_trn_tst_Mean>Q_mean)
    Trn_Set_Mean=TRN_Set[Mean_pos]
    
    Opt_trn5=Trn_Set_Max
    
    All_sample5=c(Opt_trn5,Tst_Set)
    PhenoTuning5=Pheno_E[All_sample5,]
    PhenoTuning5=droplevels(PhenoTuning5)
    ZL5=model.matrix(~0+Line ,data=PhenoTuning5)
    GID_Samples5=unique( PhenoTuning5$Line)
    GenoTuning5=Geno[ GID_Samples5, GID_Samples5]
    GRM5=ZL5%*%GenoTuning5%*%t(ZL5)
    
    ETATuning5=list(Line=list(model='RKHS',K=GRM5))
    
    yy5=PhenoTuning5[,trait]
    yy_na5=yy5
    TRN_Set5=1:length(Opt_trn5)
    yy_na5[-TRN_Set5]=NA
    model5 <- BGLR::BGLR(
      y = yy_na5,
      ETA = ETATuning5,
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE
    )
    Predicted_trait5=model5$yHat
    #Predictions$Predicted_trait=Predicted_trait
    #Predictions$Observed_trait=yy
    MSE5=mse(yy5[-TRN_Set5],Predicted_trait5[-TRN_Set5])
    COR5=cor(yy5[-TRN_Set5],Predicted_trait5[-TRN_Set5])
    NRMSE5=nrmse(yy5[-TRN_Set5],Predicted_trait5[-TRN_Set5])
    
    Opt_trn6=Trn_Set_Mean
    
    All_sample6=c(Opt_trn6,Tst_Set)
    PhenoTuning6=Pheno_E[All_sample6,]
    PhenoTuning6=droplevels(PhenoTuning6)
    ZL6=model.matrix(~0+Line ,data=PhenoTuning6)
    GID_Samples6=unique( PhenoTuning6$Line)
    GenoTuning6=Geno[ GID_Samples6, GID_Samples6]
    GRM6=ZL6%*%GenoTuning6%*%t(ZL6)
    
    ETATuning6=list(Line=list(model='RKHS',K=GRM6))
    
    yy6=PhenoTuning6[,trait]
    yy_na6=yy6
    TRN_Set6=1:length(Opt_trn6)
    yy_na6[-TRN_Set6]=NA
    model6 <- BGLR::BGLR(
      y = yy_na6,
      ETA = ETATuning6,
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE
    )
    Predicted_trait6=model6$yHat
    #Predictions$Predicted_trait=Predicted_trait
    #Predictions$Observed_trait=yy
    MSE6=mse(yy6[-TRN_Set6],Predicted_trait6[-TRN_Set6])
    COR6=cor(yy6[-TRN_Set6],Predicted_trait6[-TRN_Set6])
    NRMSE6=nrmse(yy6[-TRN_Set6],Predicted_trait6[-TRN_Set6])
    
    
    Summary_Accuracy_Testing=rbind(Summary_Accuracy_Testing,data.frame(Trait=trait, MSE= MSE,COR=COR, NRMSE= NRMSE,MSE_AV= MSE2,COR_AV=COR2, NRMSE_AV= NRMSE2,MSE_CD= MSE3,COR_CD=COR3, NRMSE_CD= NRMSE3,MSE_PEV= MSE4,COR_PEV=COR4, NRMSE_PEV= NRMSE4,MSE_AVG= MSE5,COR_AVG=COR5, NRMSE_AVG= NRMSE5,MSE_Max= MSE6,COR_Max=COR6, NRMSE_Max= NRMSE6 ))
    
    Predictions_All_traits=rbind(Predictions_All_traits,Predictions)
    
    classes <- c("No", "Yes")  
    Predictions_Adapted=data.frame(Fold=Predictions$Fold,Observed=as.factor(ifelse(Predictions$Observed=="1", "Yes", "No")),Predicted=as.factor(ifelse(Predictions$Predicted=="1", "Yes", "No")), Yes=Predictions$prob,No=1-Predictions$prob)
    
    Summary <- Predictions_Adapted %>%
      group_by(Fold) %>%
      dplyr::summarise(
        ROC_AUC = roc_auc(
          Observed,
          select_at(across(), classes),
          positive_class = "Yes"
        ))
    Summary
    
    Global <- Summary %>%
      summarise_all(mean)
    Global
    
    Final <- vctrs::vec_rbind(Summary, Global) %>%
      mutate(Fold = as.character(Fold)) 
    #%>%
    #  round(digits = 4)
    
    Final$Fold[nrow(Final)] <- "Global"
    
    
    AUC=Global[-1]
    Summary_AUC=rbind(Summary_AUC,data.frame(Trait=trait,AUC=AUC))   
    
    
  }
  
  TT=apply(Summary_Accuracy_Testing[,-1],2,mean)
  TT1=data.frame(Trait=Summary_Accuracy_Testing[1,1], MSE=TT[1],COR=TT[2], NRMSE=TT[3],MSE_AV= TT[4],COR_AV=TT[5], NRMSE_AV=TT[6],MSE_CD=TT[7],COR_CD=TT[8], NRMSE_CD=TT[9],MSE_PEV=TT[10],COR_PEV=TT[11], NRMSE_PEV= TT[12],MSE_AVG=TT[13],COR_AVG=TT[14], NRMSE_AVG=TT[15],MSE_Max= TT[16],COR_Max=TT[17], NRMSE_Max=TT[18] )
  
  RR=mean( Summary_AUC[,-1])
  RR1=data.frame(Trait= Summary_AUC[1,1],  AUC=RR[1])
  
  
  Summary_AUC_Fold=rbind(Summary_AUC_Fold,RR1)
  Predictions_All_traits_Fold=rbind(Predictions_All_traits_Fold,Predictions_All_traits)
  Summary_Accuracy_Testing_Fold=rbind(Summary_Accuracy_Testing_Fold,TT1) 
}
write.csv(Summary_AUC_Fold,file="Summary_AUC_Groundnut_Random_V3.csv")
write.csv(Summary_Accuracy_Testing_Fold,file="Summary_Accuracy_Testing_Groundnut_Random_V3.csv")
write.csv(Predictions_All_traits_Fold,file="Predictions_All_traits_Groundnut_Random_V3.csv")






