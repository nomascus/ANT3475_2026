# Notes pour l'atelier

x <- 3 * (14/2)
x
X

num_vecteur <- c(1, 4, 3, -1)
num_vecteur

chr_vecteur <- c("bleu", "rouge")
chr_vecteur

x <- c(1,2,3,4,5)
x <- 1:5
y <- 6
x * y
x + y


# On commence par créer les vecteurs.
site_id <- c("A1.01", "A1.02", "B1.01", "B1.02") # l'identifiant du site
pH_sol <- c(5.6, 7.3, 4.1, 6.0) #la mesure du pH du sol à chaque site
num_sp <- c(17, 23, 15, 7) #le nombre d'espèces observées
traitement <- c("Fert", "Fert", "Non_fert", "Non_fert") #traitement appliqué

# On peut grouper tous ces vecteurs en un data frame avec la fonction data.frame().
mon_df <- data.frame(site_id, pH_sol, num_sp, traitement)

# On l’affiche à la console!
mon_df
View(mon_df)

impair <- c(1,3,5,7,9)

impair[2]
impair[c(3,5)]
impair[c(-1)]
impair[6]
impair[impair > 4]

chr_vecteur[chr_vecteur == "rouge"]

mon_df[mon_df$traitement == 'Fert',]

mon_df[1,]
mon_df[,3]
mon_df[2,4]

mon_df$traitement
mon_df[c("site_id","traitement")]

install.packages(c("ggplot2", "tidyr", "dplyr", "palmerpenguins", "lubridate"))

## Telechargez les fichiers de StudiUM
setwd("endroit et nom de repertoire que tu as choisi")
CO2 <- read.csv("co2_good.csv", header = TRUE)
CO2[,2]
CO2[3,4]
CO2[CO2$Treatment == "chilled",]
