# a script to calculate the probability of detecting a father by sub-sampling a clutch of offspring
# assumes only one mother per clutch
# assumes enough polymorphim / power to assign each offspring

# setwd("/Users/tbroquet/Documents/Recherche/Jaera_projet2/Data/Manip1_sexualselection/")

prob <- function(x,n){
1-((1-x)^n)                  # probability of detecting a father contributing a fraction x of a clutch by sampling n larvae
}

N <- function(x,n){          # expected number of offspring in a sample of size n
round(x*n)
}

n <- c(5,10,12,15,20)             # sample size
x <- seq(0,0.5,0.001)         # actual contribution of a father
y1 <- sapply(n,prob,x=x)       # calculate fonction "prob" defined above for all n and x values
y2 <- sapply(n,N,x=x)          # idem with function "N"
# 
# pdf("detection_power.pdf",width=6,height=6,useDingbats=T)
#   matplot(x,y1,type='l',lwd=2,lty=c(3,1,2,1,2),
#           col=c("black","black","black", "grey","grey"),
#           xlab="Actual father contribution",
#           ylab="Probability of detecting at least one offspring",
#           cex.lab=0.8,
#           las=1
#         )
#   legend(0.15,0.4,legend=n,lwd=2,lty=c(3,1,2,1,2),bty="n",title="sample size", col=c("black","black","black","grey","grey"), cex=0.8)
#   abline(v=0.1)
# 
# dev.off()


library(tidyverse)
p <- cbind(x, y1) %>%
  as_tibble() %>%
  rename(V5=V2,
         V10=V3,
         V12=V4,
         V15=V5,
         V20=V6) %>%
  pivot_longer(cols = contains("V"), names_to = "Sampled_offspring", values_to = "Proba") %>%
  mutate(Sampled_offspring = factor(str_remove_all(Sampled_offspring, "V"),
                                    levels = n)) %>%
  ggplot(aes(x = x, y = Proba)) +
  geom_hline(yintercept = 0.9, lwd = 1.1, colour = "blue", lty = 1) +
  geom_line(aes(lty = Sampled_offspring, colour = Sampled_offspring), lwd = 1.1) +
  # annotate("text", x = 0.01, y = 0.95, label = "0.9", size = 8, colour = "blue") +
  scale_colour_manual(name = "Sampled offspring",
                      values = c("black", "black", "black", "grey", "grey")) +
  scale_linetype_manual(name = "Sampled offspring",
                        values = c(3, 2, 1, 1, 2)) +
  labs(x = "Actual father contribution",
       y = "Probability of detecting at least one offspring") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 20))


p
ggsave("../../../Chapter_1_Sexual_selection/Figures/Suppfigs/Proba_detecting_contribs.png", p, scale = 3, width = 1600, height = 1000, units = "px")


