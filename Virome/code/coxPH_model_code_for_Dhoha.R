


species.name <- colnames(PA.table.species.all.prev)[26:81]
new.species.name <- gsub("/|-","_",species.name)
colnames(PA.table.species.all.prev)[26:81] <- new.species.name
PA.table.species.all.prev <- PA.table.species.all.prev %>% mutate(across(starts_with("Viruses"), ~ factor(., levels = c(0, 1))))


PA.table.species.all.prev$Dx.Status <- as.numeric(ifelse(PA.table.species.all.prev$Dx.Status == "CONTROL",0,1))
table(PA.table.species.all.prev %>%  distinct(patientID,.keep_all = TRUE) %>% pull(Dx.Status))
#  0  1 
# 36 36 


species.stats.list <- list()

for (i in new.species.name) {
  result <- tryCatch({
    
    
    print(paste("Processing i =", i))
    # Example code that may cause the specific error
    cox.fit <- coxph(Surv(month, Dx.Status) ~ get(i) + Country + Sex + Delivery.Mode + HLA.Category + feeding_first_year + Age.at.Gluten.Introduction..months., data = PA.table.species.all.prev)
    
    
    fit.tidy <- tidy(cox.fit)
    
    write.csv(fit.tidy,paste0("~/Handley Lab Dropbox/16S/Celiac/plot/COX/grand_final/species/",i,"_cox.csv"))
    
    
    
    # Get the summary of the Cox model
    cox_summary <- summary(cox.fit)
    
    # Extract the coefficients, standard errors, z-values, and p-values
    model_df <- data.frame(
      term = rownames(cox_summary$coefficients),         # The names of the coefficients (terms)
      coef = cox_summary$coefficients[, 1],              # Coefficients (estimates)
      exp_coef = exp(cox_summary$coefficients[, 1]),    # Hazard ratios (exp(coef))
      se = cox_summary$coefficients[, 2],                # Standard errors
      z_value = cox_summary$coefficients[, 3],           # Z-values
      p_value = cox_summary$coefficients[, 4],           # P-values
      ci_lower = exp(cox_summary$conf.int[, 3]),        # Lower bound of 95% CI for HR
      ci_upper = exp(cox_summary$conf.int[, 4])         # Upper bound of 95% CI for HR
    )
    
    
    write.csv(model_df,paste0("~/Handley Lab Dropbox/16S/Celiac/plot/COX/grand_final/species/",i,"_cox_fill_stats.csv"))
    
    print("test1")
    
    
    df <- data.frame(summary(cox.fit)$coefficients) %>%
      filter(row.names(.) == rownames(data.frame(summary(cox.fit)$coefficients))[1]) %>%
      rownames_to_column("taxa") %>% 
      mutate(taxa = i)
    
    print(dim(df))
    
    species.stats.list[[i]] <- df
    
    print("test2")
    
    NULL  # Return a value if no error occurs
  }, error = function(e) {
    # Check if the error message matches the specific condition
    if (grepl("Downdated VtV is not positive definite", e$message)||grepl("Response is constant", e$message)) {
      message(paste("Caught specific error at iteration", i, "- skipping to next iteration"))
      return(NULL)  # Return NULL to indicate the error occurred and skip the current iteration
    }
    stop(e)  # Rethrow the error if it doesn't match the condition
  })
  
  # If result is NULL, skip to the next iteration
  if (is.null(result)) {
    next  # Skip to the next iteration if the specific error occurred
  }
}



species.cox.res <- bind_rows(species.stats.list) %>%
  mutate(`P-Value` = ifelse(`Pr...z..` < 0.05,"< 0.05","> 0.05")) %>%
  mutate(`P-adjust` = p.adjust(`Pr...z..`)) %>%
  mutate(`P-adjust-sig` = ifelse(`P-adjust` < 0.05,"< 0.05","> 0.05")) %>%
  mutate("taxa_species" = unlist(lapply(strsplit(taxa,";"),"[",7))) %>%
  mutate("taxa_genus" = unlist(lapply(strsplit(taxa,";"),"[",6))) %>%
  mutate("taxa_family" = unlist(lapply(strsplit(taxa,";"),"[",5))) %>%
  mutate("taxa_genus_species" = paste0(unlist(lapply(strsplit(taxa,";"),"[",6)),"_",unlist(lapply(strsplit(taxa,";"),"[",7)))) %>%
  mutate(taxa_species = gsub("_","-",taxa_species))



p.species.stats <- ggplot(species.cox.res,aes(x = `exp.coef.`,y=taxa_species,color = `P-Value`)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0,linetype="dotted") +
  scale_color_manual(values = c("< 0.05" = "red", "> 0.05" = "black")) +
  geom_text_repel(data = subset(species.cox.res, `Pr...z..` < 0.05),aes(label = format(round(`Pr...z..`,3),nsmall = 2)),show.legend = FALSE) + #, aes(label = `P-adjust`)
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12,face = "bold")) +
  labs(x = "exp(Coef)",y = "") +
  theme_bw()



ggsave(p.species.stats,file="~/Handley Lab Dropbox/16S/Celiac/plot/COX/grand_final/species/species.cox.res.pdf",width = 10, height = 15,dpi=600)


species_to_genus <- unlist(split(species.cox.res$taxa_genus,species.cox.res$taxa_species))

library(scales) 

custom_trans <- trans_new(
  name = "custom",
  transform = function(x) {
    # Piecewise transformation
    ifelse(x < 0.05, x * 20,  # Map 0.00-0.05 to 0-1
           ifelse(x < 0.25, (x - 0.05) * 5 + 1,  # Map 0.05-0.25 to 1-2
                  (x - 0.25) * 1 + 2))  # Map 0.25+ proportionally
  },
  inverse = function(x) {
    # Reverse the transformation
    ifelse(x < 1, x / 20,  # Map 0-1 back to 0.00-0.05
           ifelse(x < 2, (x - 1) / 5 + 0.05,  # Map 1-2 back to 0.05-0.25
                  (x - 2) / 1 + 0.25))  # Map 2+ back to 0.25+
  }
)


p.species.stats.Pvalue <- ggplot(species.cox.res %>% filter(`Pr...z..` < 0.25),aes(x = `Pr...z..`,y=reorder(taxa_species, -`Pr...z..`),color = taxa_family)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0.05,linetype="dotted") +
  scale_color_manual(breaks = rev(unique(species.cox.res %>% filter(`Pr...z..` < 0.25) %>% pull(taxa_family))),values = c("#FF5733","#33FF57","#3357FF","#FF33A1","#FFFF33","#FF33FF","#DFFF00","#ab2be2","#FF6347","#FF8C00"), #"#33FFFF",
                     labels = rev(unique(unique(species.cox.res %>% filter(`Pr...z..` < 0.25) %>% pull(taxa_family))))) +  #values = c("< 0.05" = "red", "> 0.05" = "black")
  #geom_text_repel(data = subset(species.cox.res, `Pr...z..` < 0.05),aes(label = format(round(`Pr...z..`,3),nsmall = 2)),show.legend = FALSE) + #, aes(label = `P-adjust`)
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12,face = "bold")) +
  labs(x = "",y = "",color = "Family") + 
  scale_y_discrete(labels = species_to_genus) +
  scale_x_continuous(
    trans = custom_trans,  # Apply the custom transformation
    breaks = c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25),  # Specify breaks
    labels = c("0.00", "0.05", "0.10", "0.15", "0.20", "0.25")  # Specify labels
  ) +
  theme_bw()


ggsave(p.species.stats.Pvalue,file="~/Handley Lab Dropbox/16S/Celiac/plot/COX/grand_final/species/species.cox.res.pValue.genusName.update.pdf",width = 8, height = 6,dpi=600)


sig.species.adjust <- species.cox.res %>% 
  filter(`Pr...z..` < 0.05) %>%
  pull(taxa)


PA.table.species.all.prev <- PA.table.species.all.prev %>% mutate(across(starts_with("Viruses"), ~ as.numeric(as.character(.))))


PA.table.sig.adjust <- PA.table.species.all.prev %>%
  dplyr::select(patientID,Country,month,Dx.Status,disease_onset_details,all_of(sig.species.adjust)) %>%
  pivot_longer(cols = starts_with("Viruses"),
               names_to = "Species", 
               values_to = "Count") %>%
  group_by(month,Dx.Status,Species) %>%
  mutate(total_count = sum(Count))



PA.table.sig.adjust$Species <- lapply(strsplit(PA.table.sig.adjust$Species,";"),function(x) paste0(x[5],";",x[6],";",x[7]))



PA.table.sig.adjust$Dx.Status <- ifelse(PA.table.sig.adjust$Dx.Status == 0,"CONTROL","CELIAC")
sig_Species_PA.heatmap.adjust <- ggplot(PA.table.sig.adjust, aes(x = as.character(Species), y = as.character(month), fill = total_count)) +
  geom_tile(color = "black", size = 0.5) +  # Adding a white border between tiles
  #scale_fill_manual(values = c("0" = "white", "1" = "blue")) +  # Custom colors for absence and presence
  scale_fill_gradient(low="white", high="black") +
  geom_text(aes(label = total_count), color = "white", size = 3) + 
  facet_wrap(~Dx.Status) +
  labs(title = "PA numbers",
       x = "Groups", y = "Timepoints", fill = "total_PA") +
  theme_minimal() +
  theme(axis.text.x = element_text(size=15,face="bold",angle = 90,hjust = 1),
        axis.text.y = element_text(size=15,face="bold"),
        strip.text.x = element_text(size = 15, face = "bold"))


ggsave(sig_Species_PA.heatmap.adjust,file ="~/Handley Lab Dropbox/16S/Celiac/plot/COX/species_prev/Species_PA.heatmap.pdf",width = 15,height = 15)



# check1 <- PA.table.species.all.prev %>% dplyr::select(patientID,Dx.Status,month,`Viruses;unclassified Viruses phylum;unclassified Viruses class;unclassified Viruses order;Anelloviridae;Gammatorquevirus;Torque teno midi virus 8`) %>%
#   dplyr::rename("Torque_teno_midi_virus8" = "Viruses;unclassified Viruses phylum;unclassified Viruses class;unclassified Viruses order;Anelloviridae;Gammatorquevirus;Torque teno midi virus 8")
#
#
# check.anello2.PA.number <- ggplot(check1,aes(x = patientID,y =as.character(month),fill = Torque_teno_midi_virus8)) +
#   geom_tile(color = "black", size = 0.5) +  # Adding a white border between tiles
#   #scale_fill_manual(values = c("0" = "white", "1" = "blue")) +  # Custom colors for absence and presence
#   scale_fill_gradient(low="white", high="black") +
#   geom_text(aes(label =Torque_teno_midi_virus8), color = "red", size = 3) +
#   facet_wrap(~Dx.Status) +
#   labs(title = "PA numbers",
#        x = "Groups", y = "Timepoints", fill = "total_PA") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(size=15,face="bold",angle = 90,hjust = 1),
#         axis.text.y = element_text(size=15,face="bold"),
#         strip.text.x = element_text(size = 15, face = "bold"))

# ggsave(check.anello2.PA.number,file = "~/Handley Lab Dropbox/16S/Celiac/plot/COX/species_prev/anello_Torque_teno_midi_virus8_PAcheck.pdf",width = 25,height = 8,dpi = 600)


