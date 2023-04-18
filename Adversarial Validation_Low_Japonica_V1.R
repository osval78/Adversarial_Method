rm(list=ls())
library(BGLR) #load BLR
library(dplyr)
library(caret)
library(plyr)
library(tidyr)
library(SKM)

library(reshape2)
#####Loading data set#############
load("Japonica.RData")
ls()
Pheno=Pheno
dim(Pheno)
Geno=Geno
head(Pheno)
#Traits_to_evaluate
Traits_to_evaluate=colnames(Pheno)[c(2:5)]
ToP_Selection=20
Threshold_Opt_trn=0.55
# Variables definitions ---------------------------------------------------------
#######Model to be implemented
model_name <- "BGBLUP"
Name_data_set="Japonica"
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


Pheno$Env="Japonica"
Pheno_E=Pheno
Predictions_All_traits=data.frame()
Summary_all_traits=data.frame()
Summary_all_traits_Testing=data.frame()
#######Fold loop for each trait
for (trait in Traits_to_evaluate) {
# trait=Traits_to_evaluate[1]
  SKM::echo("*** Trait: %s ***", trait)
  BLUEs_trait=Pheno_E[,trait]
  Q.20=quantile(BLUEs_trait,prob=(ToP_Selection/100))
  ######Training testing partitions
   TRN_Set=which(BLUEs_trait> Q.20)
   Tst_Set=which(BLUEs_trait< Q.20)
   y_Bin=rep(0,nrow(Pheno_E))

   y_Bin[-TRN_Set]=1
   Bin_name=paste("Bin",trait, sep="_")
   Pheno_E=cbind(Pheno_E,Y_Bin=y_Bin)
   #head(Pheno_E)
   colnames(Pheno_E)=c(colnames(Pheno_E)[-ncol(Pheno_E)],Bin_name[1])

   ####Preprosesing steps#
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
      ######Inner ten fold cross validation
      inner_folds <- SKM::cv_kfold_strata(y_tuning_Bin_Factor, k = tuning_folds_number)
     
      Predictions=data.frame()
      
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
        
      Predictions <- rbind(
         Predictions,
         data.frame(
            Trait=trait,
            Line = Pheno_E$Line[inner_fold$testing],
            Fold = k,
            Observed = Observed,
            Predicted =Predicted,
            prob=prob_val,
            Opt_trn=Opt_trn

         )
      )
      
     }
      
      ##########Predictions of testing with full data
      yy=BLUEs_trait
      yy_na=yy
      yy_na[-TRN_Set]=NA
      model <- BGLR::BGLR(
        y = yy_na,
        ETA = ETATuning1,
        nIter = iterations_number,
        burnIn = burn_in,
        verbose = FALSE
      )
      Predicted_trait=model$yHat
      Predictions$Predicted_trait=Predicted_trait
      Predictions$Observed_trait=yy
      MSE=mse(yy[-TRN_Set],Predicted_trait[-TRN_Set])
      COR=cor(yy[-TRN_Set],Predicted_trait[-TRN_Set])
      NRMSE=nrmse(yy[-TRN_Set],Predicted_trait[-TRN_Set])
      
      ##########Predictions of testing with optimal training
      Opt_trn=which(Predictions$Opt_trn==1)
      
      All_sample=c(Opt_trn,Tst_Set)
      PhenoTuning2=Pheno_E[All_sample,]
      PhenoTuning2=droplevels(PhenoTuning2)
      ZL2=model.matrix(~0+Line,data=PhenoTuning2)
      GID_Samples=unique( PhenoTuning2$Line)
      GenoTuning2=GenoTuning[ GID_Samples, GID_Samples]
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

      Summary_all_traits_Testing=rbind(Summary_all_traits_Testing,data.frame(Trait=trait, MSE= MSE,COR=COR, NRMSE= NRMSE,MSE_Opt= MSE2,COR_Opt=COR2, NRMSE_Opt= NRMSE2 ))
      
      Predictions_All_traits=rbind(Predictions_All_traits,Predictions)
      #########Computing metrics for binary data
      classes <- c("No", "Yes")
      Predictions_Adapted=data.frame(Fold=Predictions$Fold,Observed=as.factor(ifelse(Predictions$Observed=="1", "Yes", "No")),Predicted=as.factor(ifelse(Predictions$Predicted=="1", "Yes", "No")), Yes=Predictions$prob,No=1-Predictions$prob)
      
      Summary <- Predictions_Adapted %>%
        group_by(Fold) %>%
        dplyr::summarise(
          ROC_AUC = roc_auc(
            Observed,
            select_at(across(), classes),
            positive_class = "Yes"
          )
        )
      Summary
      
      Global <- Summary %>%
        summarise_all(mean)
      
      Final <- vctrs::vec_rbind(Summary, Global) %>%
        mutate(Fold = as.character(Fold)) 

      Final$Fold[nrow(Final)] <- "Global"
   Summary_all_traits=rbind(Summary_all_traits,data.frame(Trait=trait,Global[-1] ))   
  

}
Summary_all_traits
Predictions_All_traits
Summary_all_traits_Testing

write.csv(Summary_all_traits,file="Summary_all_traits_Japonica_Low_V1.csv")
write.csv(Summary_all_traits_Testing,file="Summary_all_traits_Testing_Japonica_Low_V1.csv")
write.csv(Predictions_All_traits,file="Predictions_All_traits_Japonica_Low_V1.csv")








