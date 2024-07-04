
library(ggplot2)

# lecture du fichier de comparaison des méthodes 'cdna' et 'cdna_exon'
tab_cdna = read.csv("/home/vetea/Documents/premiere_annee/second_semestre/stage/CDScompR/results/cdna_cdna_exon.csv", 
                    header=TRUE, 
                    sep=",", 
                    dec=".")
# modification des noms de colonnes du jeu de données pour éviter des bugs
colnames(tab_cdna) = c("chromosome", "cluster", "reference_locus", "alternative_locus", "matches", "mismatches", "identity", "reference_start", "reference_end", "alternative_start", "alternative_end", "ref_mRNA", "alt_mRNA", "EI_mismatch_zones", "RF_mismatch_zones", "EI_mismatches", "RF_mismatches", "ref_mRNA_number", "alt_mRNA_number")
# histogramme des valeurs d'identité de la comparaison
hist(tab_cdna$identity, 
     main="Histogramme de l'identité pour la comparaison cdna - cdna-exon", 
     xlab="Identité", 
     ylab="effectif")

# lecture du fichier de comparaison des méthodes 'prot' et 'prot_exon'
tab_prot = read.csv("/home/vetea/Documents/premiere_annee/second_semestre/stage/CDScompR/results/prot_prot_exon.csv", 
                    header=TRUE, 
                    sep=",", 
                    dec=".")
# modification des noms de colonnes du jeu de données pour éviter des bugs
colnames(tab_prot) = c("chromosome", "cluster", "reference_locus", "alternative_locus", "matches", "mismatches", "identity", "reference_start", "reference_end", "alternative_start", "alternative_end", "ref_mRNA", "alt_mRNA", "EI_mismatch_zones", "RF_mismatch_zones", "EI_mismatches", "RF_mismatches", "ref_mRNA_number", "alt_mRNA_number")
# histogramme des valeurs d'identité de la comparaison
hist(tab_prot$identity, 
     main="Histogramme de l'identité pour la comparaison prot - prot-exon", 
     xlab="Identité", 
     ylab="effectif")

# lecture du fichier de comparaison de toutes les méthodes (contre la référence expertisée)
tab_methods = read.csv("/home/vetea/Documents/premiere_annee/second_semestre/stage/CDScompR/results/may_methods_combined.csv", 
                    header=TRUE, 
                    sep=",", 
                    dec=".")

ggplot(tab_methods, aes(x=ident_cdna, y=ident_cdna_exon)) +
  geom_point()

diff_cdna <- tab_methods$ident_cdna -tab_methods$ident_cdna_exon
summary(diff_cdna)

# vérification de l'égalité du résultat de 'best' avec le maximum des résultats des méthodes
tab_methods$max <- pmax(tab_methods$ident_cdna, tab_methods$ident_cdna_exon, tab_methods$ident_prot, tab_methods$ident_prot_exon, tab_methods$ident_mapping)
tab_methods$best_is_max <- tab_methods$ident_best == tab_methods$max
tab_methods$diff_best_max <- tab_methods$max - tab_methods$ident_best
print(all(tab_methods$best_is_max)) # affiche 'TRUE' si 'best' correspond pour chaque locus au
                                    # maximum des résultats de méthodes

write.csv(tab_methods, "Comparaison_des_methodes_lrrtransfer.csv")

# tests appariés pour déterminer si l'une des méthodes approxime les résultats de 'best'

t.test(tab_methods$ident_best, 
       tab_methods$ident_cdna, 
       paired=TRUE)

t.test(tab_methods$ident_best, 
       tab_methods$ident_cdna_exon, 
       paired=TRUE)

t.test(tab_methods$ident_best, 
       tab_methods$ident_prot, 
       paired=TRUE)

t.test(tab_methods$ident_best, 
       tab_methods$ident_prot_exon, 
       paired=TRUE)

# mapping renvoie une erreur, mais de toute façon, les résultats sont mauvais

# t.test(tab_methods$ident_best, 
#        tab_methods$ident__mapping, 
#        paired=TRUE)

# graphe de visualisation de la densité des résultats de 'best' et de 'prot' pour
# vérifier la correspondance
plot(density(tab_methods$ident_best))
lines(density(tab_methods$ident_prot))

