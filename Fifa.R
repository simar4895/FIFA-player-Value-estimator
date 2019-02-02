library(pacman)

p_load(tidyverse)

fifa <- read.csv("complete.csv")

summary(fifa)

fifa_new <- fifa %>% select(c(7,17:18,20,28))

summary(fifa_new)

write.csv(fifa_new,"fifa.csv",row.names = F)

hist(fifa_new$age)
hist(fifa_new$eur_value)
hist(fifa$eur_wage)
hist(fifa_new$overall)
hist(fifa_new$international_reputation)

plot(fifa_new$eur_wage,fifa_new$eur_value)
plot(fifa_new$age,fifa_new$eur_value)
plot(fifa_new$overall,fifa_new$eur_value)
plot(fifa_new$international_reputation,fifa_new$eur_value)

a <- lm(eur_value~eur_wage+international_reputation+overall+age+I(age^2),data = fifa_new)
summary(a)

