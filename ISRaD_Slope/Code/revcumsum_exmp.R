c <- c(2.17, 0.98, 1.17, 0.52, 0.32, 0.24, 0.18, 0.14, 0.19, 0.14, 0.27, 0.26, 0.3)
a <- c(16,-17, -75,-145, -125,  -162, -201, -221, -289, -308, -322, -349, -340)

df <- data.frame(c, a)

library(tidyverse)

df %>% 
  arrange(desc(a)) %>% 
  mutate(diff = rev(cumsum(rev(c*100/sum(c))))) %>% 
  ggplot(aes(x = a, y = diff)) +
  geom_path() +
  geom_vline(xintercept = 0, linetype = "dashed")
   

