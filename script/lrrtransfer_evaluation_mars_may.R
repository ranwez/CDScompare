
#devtools::install_github("kassambara/ggpubr")
library(ggplot2)

# lecture du fichier de comparaison de la comparaison des résultats de Mars avec la référence expertisée
tab_mars = read.csv("/home/vetea/Documents/premiere_annee/second_semestre/stage/annot_CSC/results/exp_best_mars.csv", 
                    header=TRUE, 
                    sep=",", 
                    dec=".")
# modification des noms de colonnes du jeu de données pour éviter des bugs
colnames(tab_mars) = c("chromosome", "cluster", "reference_locus", "alternative_locus", "matches", "mismatches", "identity", "reference_start", "reference_end", "alternative_start", "alternative_end", "ref_mRNA", "alt_mRNA", "EI_mismatch_zones", "RF_mismatch_zones", "EI_mismatches", "RF_mismatches", "ref_mRNA_number", "alt_mRNA_number")
# histogramme des valeurs d'identité de la comparaison
hist(tab_mars$identity, 
     main="Histogramme de l'identité pour la comparaison EXP-Mars", 
     xlab="Identité", 
     ylab="effectif")

# lecture du fichier de comparaison de la comparaison des résultats de Mai avec la référence expertisée
tab_mai = read.csv("/home/vetea/Documents/premiere_annee/second_semestre/stage/annot_CSC/results/exp_best_mai.csv", 
                   header=TRUE, 
                   sep=",", 
                   dec=".")
# modification des noms de colonnes du jeu de données pour éviter des bugs
colnames(tab_mai) = c("chromosome", "cluster", "reference_locus", "alternative_locus", "matches", "mismatches", "identity", "reference_start", "reference_end", "alternative_start", "alternative_end", "ref_mRNA", "alt_mRNA", "EI_mismatch_zones", "RF_mismatch_zones", "EI_mismatches", "RF_mismatches", "ref_mRNA_number", "alt_mRNA_number")
# histogramme des valeurs d'identité de la comparaison
hist(tab_mai$identity, 
     main="Histogramme de l'identité pour la comparaison EXP-Mai", 
     xlab="Identité", 
     ylab="effectif")


# lecture du fichier combinant les comparaisons de Mars et Mai avec la référence
tab_fusion = read.csv("/home/vetea/Documents/premiere_annee/second_semestre/stage/annot_CSC/results/mars_may_combined.csv", 
                      header=TRUE, 
                      sep=",", 
                      dec=".")

# création d'un vecteur contenant la différence d'identité entre Mars et Mai
diff = tab_fusion$identity_may - tab_fusion$identity_mars
# histogramme de la différence d'identité entre Mars et Mai 
hist(diff, 
     main="Histogramme de l'identité pour l'écart Mars-Mai", 
     xlab="Identité", 
     ylab="effectif", 
     breaks = 10)

# test de normalité
shapiro.test(diff) # p-value significative = pas de normalité

# test de différence des valeurs appariées d'identité entre Mars et Mai
# n>30 pour les deux distributions, pas besoin de normalité apparemment...
t.test(tab_fusion$identity_may, 
       tab_fusion$identity_mars, 
       paired=TRUE) 

# graphe à points des valeurs d'identité pour chaque locus entre Mars (x) et Mai (y)
tab_fusion$is_canonical <- as.factor(tab_fusion$is_canonical)
ggplot(tab_fusion, aes(x=identity_mars, y=identity_may, color=is_canonical)) +
  geom_point()

# avec les graphes de densité marginaux

library(ggpubr)

# Scatter plot colored by groups ("Species")
sp <- ggscatter(tab_fusion, x = "identity_mars", y = "identity_may",
                color = "is_canonical", palette = "jco",
                size = 3, alpha = 0.6)+
  border()                                         
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(tab_fusion, "identity_mars", fill = "is_canonical",
                   palette = "jco")
yplot <- ggdensity(tab_fusion, "identity_may", fill = "is_canonical", 
                   palette = "jco")+
  rotate()
# Cleaning the plots
yplot <- yplot + clean_theme() 
xplot <- xplot + clean_theme()
# Arranging the plot
ggarrange(xplot, NULL, sp, yplot, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)

write.csv(tab_fusion, "Comparaison_transfert_Mars_Mai.csv")

