library(tidyverse)
library(ComplexHeatmap)

data_path <- "https://raw.githubusercontent.com/pughlab/ORACLE_TCR/main/Data"

#Reading the clinical data -------------------------------------------------------------------------

ClinicalData_fname <- "ORACLE_ClinicalData_UpdatedFeb2024.csv"
ClinicalData <- readr::read_csv(file.path( data_path , ClinicalData_fname)) %>%
        dplyr::select(Patient_id , IRAE , Age ,
                      Sex , Ethnicity , Smoking , 
                      `Pack years` , Stage ,
                      Histology , Mutation ,
                      Chemo , RTResponse , DurvaResponse ,
                      PFS_durva)%>%
        mutate(DurvaResponseLongevity = case_when(
                DurvaResponse == "CR" ~ "CR" ,
                DurvaResponse == "PR" ~ "PR" ,
                DurvaResponse == "SD" & PFS_durva > 11 ~ "Long-term SD" ,
                DurvaResponse == "SD" & PFS_durva < 11 ~ "Short-term SD" ,
                DurvaResponse == "PD" ~ "PD"  ))

#Reading the sample inventory data -----------------------------------------------------------------

inventory_fname <- "ORACLE_capTCRseq_Inventory.csv"
inventory <- readr::read_csv(file.path( data_path , inventory_fname))

#Defing the patients order for the heatmap ---------------------------------------------------------

PatientOrder <- (ClinicalData %>%
                          left_join(inventory %>%
                                            group_by(Patient_id) %>%
                                            summarise(N = n()) %>%
                                            ungroup() ,
                                    by = "Patient_id") %>%
                          arrange( desc (N) ,
                                   match(Histology , 
                                         c("Adeno" , "Squamous" , "LargeCell" , "NOS")) ,
                                   match(Stage , c("3A" , "3B" , "3C")) ,
                                   match(Smoking , 
                                         c("Current" , "Light" , "Former" , "Never")) ,
                                   match(Sex , c("F" , "M")) ,
                                   Ethnicity ,
                                   Chemo ,
                                   RTResponse ,
                                   DurvaResponseLongevity))$Patient_id

#Converting the sample inventory data to 0:1 matrix for heatmap body design ------------------------

sample_inventory_matrix <- inventory  %>%
        dplyr::select(Patient_id ,
                      Cycle ,
                      Source) %>%
        arrange(match(Source , c("Peripheral Blood" , 
                                 "Plasma" ,
                                 "Tumour"))) %>%
        mutate(Cycle_Source = paste(Cycle , Source , sep = "_")) %>%
        dplyr::select(-Source) %>%
        pivot_wider(values_from = Cycle ,
                    names_from = Cycle_Source) %>%
        arrange(match(Patient_id , PatientOrder)) %>%
        column_to_rownames(var = "Patient_id") %>%
        as.matrix()

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# Desining the color pallettes and annotations for the legends in the heatmap ----------------------

color_pal <- data.frame(
        Cycle = c("T0" , "T2" , "T3" , "T4" , "T5"),
        Color = c("#DAC5A5" ,  
                  "#6E3A06" , 
                  "#B0C291" , 
                  "#335A30" , 
                  "#E90000") ,
        Annotation = c("Baseline, T0" , 
                       "Post-CRT, T2" , 
                       "Pre-IO, T3" , 
                       "On-IO, T4" , 
                       "Progression, T5") )

#---------------------------

IRAE_ColorPal <- data.frame(
        IRAE = c( "Gr1Pneumonitis" , "Gr2Pneumonitis" , "Gr3Pneumonitis" ,
                  
                  "Gr1Colitis" , "Gr2Colitis" , 
                  
                  "Gr1Hypothyroid" , "Gr2Hypothyroid" ,
                  
                  "Gr2Hepatitis" , "Gr3Hepatitis" , 
                  
                  "Gr1Dermatitis" ,
                  "Gr2AdrenalInsufficiency" ,
                  "Gr2Arthritis" ,
                  "Gr2SiccaSymptoms" ,
                  "Gr2BullousPemphigoid" ,
                  "Gr4Diabetes" ,
                  "InfusionReaction" ),
        Color = c( "#9ecae1" , "#2171b5" , "#08306b" ,
                   
                   "#FFFFAC" , "#F4DD15" , 
                   
                   "#CFCFC4" , "#696969" ,
                   
                   "#9F7E61" , "#4E2F1C" , 
                   
                   "#AA2600" ,
                   "#DC6602" ,
                   "#FCA505" ,
                   "#9E8310" ,
                   "#51561F" ,
                   "#832C38" ,
                   "#000000" ))

#---------------------------

Histology_ColorPal <- data.frame(
        Histology = c("Adeno" , "Squamous" , "LargeCell" , "NOS"),
        Color = c("#255E63" , "#D24F08" , "#629B46" , "#F3CB66") ,
        Annotation = c("Adenocarcinoma" , "Squamous cell carcinoma" , "Large cell carcinoma" , "Other") )

#---------------------------

Smoking_ColorPal <- data.frame(
        SmokingStatus = c("Current" , "Light" , "Former" , "Never"),
        Color = c("#910000" , "#ED8111" , "#FFE503" , "#B0BF1A")  )

#---------------------------

Chemo_ColorPal <- data.frame(
        Regimen = c("CARBO/ETOP" , "CARBO/TAXOL" , "CARBO/PEM" , "CIS/ETOP" , "CIS/PEM"),
        Color = c("#000000" , "#902B21" , "#EDD998" , "#CC9F38" , "#A2A67C")  )

#---------------------------

Response_ColorPal <- data.frame(
        Response = c("CR" ,
                     "PR" ,
                     "Long-term SD" ,
                     "Short-term SD" ,
                     "PD" ),
        Color = c ("#033483" ,
                   "#A7D2E9" ,
                   "#7DB290" ,
                   "#FEA500" ,
                   "#C5231B" ) )

#---------------------------

Ethnicity_ColorPal <- data.frame(
        Ethnicity = c("Caucasian" ,
                      "Asian" ,
                      "African" ,
                      "Arabic" ,
                      "Other" ,
                      "Unknown" ),
        Color = c ("#E7E6DF" ,
                   "#FFF0A3" ,
                   "#2F0B00" ,
                   "#84563B" ,
                   "#51661C" ,
                   "grey" ) )

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# Visualization ------------------------------------------------------------------------------------

ht_opt$ROW_ANNO_PADDING = unit(1, "cm")
ht_opt$legend_gap = unit(c(1 , 1), "cm")
ht_opt$HEATMAP_LEGEND_PADDING = unit(1.5, "cm")
ht_opt$legend_title_gp = gpar(fontsize = 12, fontface = "plain") 
ht_opt$legend_labels_gp = gpar(fontsize = 10, fontface = "plain") 



draw (
        Heatmap(sample_inventory_matrix ,
                name = "Timepoint" , # Legend title
                
                # Defining the colors for the sample timepoints-----------------------------------------------
                col = structure(names = color_pal$Cycle, color_pal$Color) ,
                na_col = "transparent" ,
                
                # Splitting the columns based on the sample types---------------------------------------------
                
                column_split = c(rep("Peripheral Blood" , 5) ,
                                 rep("Plasma" , 5) ,
                                 rep("Tumour" , 2) ) ,
                
                # Setting the column names for the splitted columns-------------------------------------------
                column_title = c("Peripheral Blood" , 
                                 "Plasma" ,
                                 "Tumour") ,
                column_title_side = "top" ,
                column_title_gp = gpar(fontsize = 12, fontfamily = "Helvetica" , fontface = "plain") ,
                
                #Setting the Gap between the splitted columns:
                column_gap = unit(8, "mm"),
                border = FALSE , # Do you want the splitted columns have borders around them?
                
                show_column_names = TRUE,
                column_labels = structure(c (rep( c( "T0" , "T2" , "T3" , "T4" , "T5") , 2),
                                              "T0" , "T3"), 
                                          names = colnames(sample_inventory_matrix)) ,
                column_names_rot = 0 ,
                column_names_gp = gpar(fontsize = 8) ,
                column_names_centered = TRUE ,
                
                # Defing Annotation columns for each row -----------------------------------------------------
                #Starts here: --------------------------------------------------------------------------------
                #---------------------------------------------------------------------------------------------
                #---------------------------------------------------------------------------------------------

                right_annotation = rowAnnotation(
                        #---------------------------
                        Histology = (ClinicalData %>% 
                                             dplyr::select(Patient_id , Histology) %>%
                                             arrange(match(Patient_id , rownames(sample_inventory_matrix))))$Histology ,
                        #---------------------------
                        Stage = (ClinicalData %>% 
                                         dplyr::select(Patient_id , Stage) %>%
                                         arrange(match(Patient_id , rownames(sample_inventory_matrix))))$Stage ,
                        #---------------------------
                        Sex = (ClinicalData %>% 
                                       dplyr::select(Patient_id , Sex) %>%
                                       arrange(match(Patient_id , rownames(sample_inventory_matrix))))$Sex ,
                        #---------------------------
                        `Smoking status` = (ClinicalData %>% 
                                                    dplyr::select(Patient_id , Smoking) %>%
                                                    arrange(match(Patient_id , rownames(sample_inventory_matrix))))$Smoking ,
                        #---------------------------
                        Ethnicity = (ClinicalData %>% 
                                             dplyr::select(Patient_id , Ethnicity) %>%
                                             arrange(match(Patient_id , rownames(sample_inventory_matrix))))$Ethnicity ,
                        #---------------------------
                        `Chemotherapy regimen` = (ClinicalData %>% 
                                                          dplyr::select(Patient_id , Chemo) %>%
                                                          arrange(match(Patient_id , rownames(sample_inventory_matrix))))$Chemo ,
                        #---------------------------
                        `Response to RT` = (ClinicalData %>% 
                                                    dplyr::select(Patient_id , RTResponse) %>%
                                                    arrange(match(Patient_id , rownames(sample_inventory_matrix))))$RTResponse ,
                        #---------------------------
                        `Response to ICB` = (ClinicalData %>% 
                                                     dplyr::select(Patient_id , DurvaResponseLongevity) %>%
                                                     arrange(match(Patient_id , rownames(sample_inventory_matrix))))$DurvaResponseLongevity ,
                        #---------------------------
                        IRAE = as.matrix ((ClinicalData %>%
                                                   dplyr::select(Patient_id , IRAE) %>%
                                                   separate_rows(IRAE , sep = "-and") %>%
                                                   
                                                   rowwise() %>%
                                                   mutate(IRAEType = case_when(
                                                           isTRUE(grepl (IRAE , pattern = "Gr")) ~ substr (IRAE , start = 4 , stop = 100) ,
                                                           TRUE ~ IRAE )) %>%
                                                   ungroup() %>%
                                                   arrange(match(IRAEType , c("Pneumonitis" , "Colitis" , "Hypothyroid" , "Hepatitis" , 
                                                                              "Dermatitis" , "AdrenalInsufficiency" , "Arthritis" , 
                                                                              "SiccaSymptoms" , "BullousPemphigoid" , "Diabetes" ,
                                                                              "InfusionReaction" , "None" , "Unknown"))) %>%
                                                   pivot_wider(
                                                           values_from = IRAE ,
                                                           names_from = IRAEType
                                                   ) %>%
                                                   arrange(match(Patient_id , rownames(sample_inventory_matrix)))) [, -c(1 , 13 , 14)]) ,
                        #---------------------------
                        Mutations = as.matrix((ClinicalData %>%
                                                       dplyr::select(Patient_id , Mutation) %>%
                                                       separate_rows(Mutation , sep = "-and") %>%
                                                       rowwise() %>%
                                                       mutate(MotherMutation = stringr::str_split(Mutation , pattern = "-")[[1]][1]) %>%
                                                       ungroup() %>%
                                                       arrange(match(MotherMutation ,
                                                                     c("TP53" , "KRAS" , "EGFR" , "ALK" , 
                                                                       "CHEK2" , "PIK3CA" , "RB1" ,
                                                                       "TMB4" , "MSS" , 
                                                                       "None" , "Unknown"))) %>%
                                                       mutate(Mutation = 1) %>%
                                                       pivot_wider(names_from = MotherMutation , 
                                                                   values_from = Mutation) %>%
                                                       mutate_all(~ replace_na ( . , 0 )) %>%
                                                       arrange(match(Patient_id , rownames(sample_inventory_matrix))) %>%
                                                       mutate_all(~ as.character(.))) [, -1]) ,
                        #---------------------------
                        # Defing Colors for each Annotation ----------------------------------------------------------
                        col = list(Histology = setNames(as.character(Histology_ColorPal$Color), 
                                                        Histology_ColorPal$Histology),
                                   
                                   `Smoking status` = setNames(as.character(Smoking_ColorPal$Color), 
                                                               Smoking_ColorPal$SmokingStatus) ,
                                   
                                   `Response to ICB` = setNames(as.character(Response_ColorPal$Color), 
                                                                Response_ColorPal$Response),
                                   
                                   `Chemotherapy regimen` = setNames(as.character(Chemo_ColorPal$Color), 
                                                                     Chemo_ColorPal$Regimen),
                                   
                                   IRAE = setNames(as.character(IRAE_ColorPal$Color), 
                                                   IRAE_ColorPal$IRAE),
                                   
                                   Ethnicity = setNames(as.character(Ethnicity_ColorPal$Color), 
                                                        Ethnicity_ColorPal$Ethnicity),
                                   
                                   `Response to RT` = c("CR" = "#033483" , 
                                                        "PR" = "#A7D2E9" , 
                                                        "SD" = "#7DB290" , 
                                                        "PD" = "#C5231B") ,
                                   
                                   Stage = c("3A" = "#E4D8E0" ,                                                 
                                             "3B" = "#8F7193" , 
                                             "3C" = "#2B222C" ) ,
                                   
                                   Sex = c("M" = "#000000" , "F" = "#DAC7A9") ,
                                   
                                   Mutations = c("1" = "#000000" , "0" = "transparent")
                        ),
                        border = c(rep(FALSE , 8) , TRUE , TRUE) ,
                        gap = unit(0.1, "cm") ,
                        na_col = "transparent" ,
                        annotation_name_rot = 30 ,
                        annotation_name_gp = gpar(fontsize = 8) ,
                        
                        
                        # Defining the order of the variables in annotation legends: ------
                        # Starts here: ----------------------------------------------------
                        annotation_legend_param = list(
                                
                                Histology = list(at = Histology_ColorPal$Histology , 
                                                 labels = Histology_ColorPal$Annotation ) ,
                                
                                Sex = list(at = c("F" , "M"), 
                                           labels = c("Female" , "Male")),
                                
                                Ethnicity = list(at = Ethnicity_ColorPal$Ethnicity),
                                
                                `Smoking status` = list(at = Smoking_ColorPal$SmokingStatus),
                                
                                `Response to RT` = list(at = c("CR" , "PR" , "SD" , "PD")) ,
                                
                                `Response to ICB` = list(at = Response_ColorPal$Response) ,
                                
                                IRAE = list(at = IRAE_ColorPal$IRAE ) 
                        ) 
                        # ------------------------------------------------------- Ends here
                ) ,
                
                #---------------------------------------------------------------------------------------------
                #---------------------------------------------------------------------------------------------
                #Ends here... --------------------------------------------------------------------------------
                
                
                # Heatmap aeshetics---------------------------------------------------------------------------
                row_names_side = "left" ,
                row_names_gp = gpar(fontsize = 8),
                
                rect_gp = gpar(col = "white", lwd = 4), #Setting lines around each rectangle
                heatmap_legend_param = list(
                        labels = color_pal$Annotation ,
                        at = color_pal$Cycle ,
                        ncol = 5 ,
                        by_row = TRUE ) ,
                
                
                
                # Size of the heatmap-------------------------------------------------------------------------
                
                width = unit(8, "cm"), 
                height = unit(21, "cm") ,
                use_raster = TRUE ,
                raster_quality = 5
        ) ,
        legend_grouping = "original" ,
        heatmap_legend_side = "bottom"
)


