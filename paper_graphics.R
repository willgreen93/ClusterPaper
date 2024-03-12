library(dplyr)
library(ggplot2)

vl_func <- function(a, b, vl_max, vi, t){
  t1 <- (vl_max-vi)/a
  vl <- ifelse(t<t1, vi+a*t, vl_max-b*(t-t1))
  return(pmax(vl, vi))
}

vl_func(a=1.5, b=1, vl_max=8.5, vi=0, t=seq(0,15,0.1))

vl_df <- data.frame(t=seq(0,15,0.1)) %>%
  mutate(vl=vl_func(a=2, b=1, vl_max=8.5, vi=0, t=t)) %>%
  mutate(thresh=ifelse(t==1.5, vl, NA)) %>%
  mutate(lag=dnorm(x=t, mean=4.5, sd=1)) %>%
  mutate(peak=ifelse(vl==max(vl), vl, NA))

point_s = 3

# tdist
ggplot(vl_df, aes(x=t)) +
  theme_bw() +
  geom_rect(aes(xmin=point_s, xmax=point_s+2.0, ymin=0, ymax=Inf), fill="red", alpha=0.002) +
  geom_rect(aes(xmin=point_s+0.4, xmax=point_s+1.6, ymin=0, ymax=Inf), fill="red", alpha=0.006) +
  geom_rect(aes(xmin=point_s+0.8, xmax=point_s+1.2, ymin=0, ymax=Inf), fill="red", alpha=0.01) + 
  geom_line(aes(y=vl), size=1.5) +
  geom_point(aes(y=thresh), size=5, color="blue") + 
  labs(y="Viral Load", x="time") +
  theme(panel.grid=element_blank()) + 
  geom_segment(
    aes(x = 1.5, xend = 4.0, y = 3, yend = 3),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    linetype = "dashed", color = "blue", size=0.5
  )

##thresh and thresh peak  
ggplot(vl_df, aes(x=t)) +
  theme_bw() +
  geom_rect(aes(xmin=0.5, xmax=3, ymin=2, ymax=5), fill="blue", alpha=0.002) +
  geom_rect(aes(xmin=0.5, xmax=3, ymin=2.5, ymax=4.5), fill="blue", alpha=0.006) +
  geom_rect(aes(xmin=0.5, xmax=3, ymin=3, ymax=4), fill="blue", alpha=0.01) + 
  geom_rect(aes(xmin=point_s, xmax=point_s+2.0, ymin=0, ymax=Inf), fill="red", alpha=0.002) +
  geom_rect(aes(xmin=point_s+0.4, xmax=point_s+1.6, ymin=0, ymax=Inf), fill="red", alpha=0.006) +
  geom_rect(aes(xmin=point_s+0.8, xmax=point_s+1.2, ymin=0, ymax=Inf), fill="red", alpha=0.01) + 
  geom_line(aes(y=vl), size=2) +
  #geom_point(aes(y=thresh), size=5, color="red") + 
  labs(y="Viral Load", x="time") +
  theme(panel.grid=element_blank()) 
  
  #geom_segment(
  #  aes(x = 1.75, xend = 5.0, y = 3.5, yend = 3.5),
  #  arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
  #  linetype = "dashed", color = "red", size=0.5
  #) #+
  #geom_segment(x = 5, xend=5, y=0, yend=7.7, color="red", linetype="dashed", size=1)

point_s <- 4
ggplot(vl_df, aes(x=t)) +
  theme_bw() +
  geom_rect(aes(xmin=point_s, xmax=point_s+2.0, ymin=0, ymax=Inf), fill="red", alpha=0.002) +
  geom_rect(aes(xmin=point_s+0.3, xmax=point_s+1.6, ymin=0, ymax=Inf), fill="red", alpha=0.006) +
  geom_rect(aes(xmin=point_s+0.6, xmax=point_s+1.1, ymin=0, ymax=Inf), fill="red", alpha=0.01) + 
  geom_line(aes(y=vl), size=1.5) +
  geom_point(aes(y=peak), size=5, color="blue") + 
  labs(y="Viral Load", x="time") +
  theme(panel.grid=element_blank()) + 
  geom_segment(
    aes(x = 4.3, xend = 5, y = 8.5, yend = 8.5),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    linetype = "dashed", color = "blue", size=0.5
  )
