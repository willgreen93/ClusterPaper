library(dplyr)
library(ggplot2)
library(scales)

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
symp1 <- ggplot(vl_df, aes(x=t)) +
  theme_bw() +
  geom_rect(aes(xmin=point_s, xmax=point_s+2.0, ymin=0, ymax=Inf), fill="red", alpha=0.002) +
  geom_rect(aes(xmin=point_s+0.4, xmax=point_s+1.6, ymin=0, ymax=Inf), fill="red", alpha=0.006) +
  geom_rect(aes(xmin=point_s+0.8, xmax=point_s+1.2, ymin=0, ymax=Inf), fill="red", alpha=0.01) + 
  geom_line(aes(y=vl), size=1.5) +
  geom_point(aes(y=thresh), size=5, color="blue") + 
  labs(y="", x="") +
  theme(panel.grid=element_blank()) + 
  geom_segment(
    aes(x = 1.5, xend = 4.0, y = 3, yend = 3),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    linetype = "dashed", color = "blue", size=0.5
  )

##thresh and thresh peak  
point_s <- 4

symp2 <- ggplot(vl_df, aes(x=t)) +
  theme_bw() +
  geom_rect(aes(xmin=0.5, xmax=3, ymin=2, ymax=5), fill="blue", alpha=0.002) +
  geom_rect(aes(xmin=0.5, xmax=3, ymin=2.5, ymax=4.5), fill="blue", alpha=0.006) +
  geom_rect(aes(xmin=0.5, xmax=3, ymin=3, ymax=4), fill="blue", alpha=0.01) + 
  geom_rect(aes(xmin=point_s, xmax=point_s+2.0, ymin=0, ymax=Inf), fill="red", alpha=0.002) +
  geom_rect(aes(xmin=point_s+0.4, xmax=point_s+1.6, ymin=0, ymax=Inf), fill="red", alpha=0.006) +
  geom_rect(aes(xmin=point_s+0.8, xmax=point_s+1.2, ymin=0, ymax=Inf), fill="red", alpha=0.01) + 
  geom_line(aes(y=vl), size=2) +
  labs(y="", x="Time since infection (days)") +
  theme(panel.grid=element_blank()) 
  
symp3 <- ggplot(vl_df, aes(x=t)) +
  theme_bw() +
  geom_rect(aes(xmin=point_s, xmax=point_s+2.0, ymin=0, ymax=Inf), fill="red", alpha=0.002) +
  geom_rect(aes(xmin=point_s+0.3, xmax=point_s+1.6, ymin=0, ymax=Inf), fill="red", alpha=0.006) +
  geom_rect(aes(xmin=point_s+0.6, xmax=point_s+1.1, ymin=0, ymax=Inf), fill="red", alpha=0.01) + 
  geom_line(aes(y=vl), size=1.5) +
  geom_point(aes(y=peak), size=5, color="blue") + 
  labs(y="Viral load", x="") +
  theme(panel.grid=element_blank()) + 
  geom_segment(
    aes(x = 4.3, xend = 5, y = 8.5, yend = 8.5),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    linetype = "dashed", color = "blue", size=0.5
  )

comb_symp_rel <- ggpubr::ggarrange(symp3, symp2, symp1, nrow=1, labels=c("A", "B", "C"))
ggsave("images/vk_symp_relationship.png", comb_symp_rel, width=10, height=3)

infecteds <- data.frame(t=seq(10,40,1)) %>%
  mutate(infecteds=ifelse(t<25, exp(t/2), exp((50-t)/2)))

p_inf <- ggplot(infecteds, aes(x=t, y=infecteds)) +
  geom_col(fill=c("brown1")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y="Infection incidence (M)", x=NULL) +
  lims(x=c(15,40)) +
  scale_y_continuous(labels = function(x) paste0(x/1e6, ""))

vls2 <- infecteds %>% mutate(infecteds=round(infecteds)) %>% tidyr::uncount(infecteds) %>%
  mutate(a=1+exp(rnorm(n=nrow(.), mean=-1, sd=0.3)),
         b=0.5+exp(rnorm(n=nrow(.), mean=-1, sd=0.2)),
         vlmax=rnorm(n=nrow(.), 3.6, 0.3),
         vt=truncnorm::rtruncnorm(n=nrow(.), a=0, mean=1.5, sd=0.4),
         symp_onset=round((vlmax-vt)/a),
         test_day=1,#floor(runif(nrow(.), 1, 5)),
         inf_test_delay=as.character(symp_onset+test_day),
         symp_calendar_day=t+symp_onset,
         test_calendar_day=symp_calendar_day+test_day) %>% 
  mutate(vl=round(vl_func(2, 1, 4, 0, symp_onset+test_day))) %>%
  mutate(vl=as.character(vl), symp_onset=as.character(symp_onset)) %>%
  dplyr::filter(symp_onset != 4)

p_symp <- ggplot(vls2, aes(x=symp_calendar_day)) +
  geom_bar(aes(fill=symp_onset)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values=c("cadetblue1", "cadetblue2", "cadetblue3", "cadetblue4")) +
  scale_y_continuous(labels = function(x) paste0(x/1e6, "")) +
  lims(x=c(15,40)) +
  labs(x=NULL, y="Symptom onset (M)", fill="Delay from infection to onset")# +
  #ggtitle("Number of people with symptom onset, by calendar day")

p_symp_test1 <- ggplot(vls2 %>% dplyr::filter(test_day==1), aes(x=test_calendar_day)) +
  geom_bar(aes(fill=inf_test_delay)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values=c("darkgoldenrod1", "darkgoldenrod2", "darkgoldenrod3", "darkgoldenrod4")) +
  scale_y_continuous(labels = function(x) paste0(x/1e6, "")) +
  lims(x=c(15,40)) +
  labs(x=NULL, y="People tested (M)", fill="Delay infection to test") #+
  #ggtitle("Number of people tested one day after symptom onset, by time taken for symptom onset from infection")
  #facet_wrap(~test_day)

vls3 <- vls2 %>% 
  group_by(test_calendar_day, test_day, inf_test_delay) %>%
  summarise(count = n()) %>%
  group_by(test_calendar_day, test_day) %>%
  mutate(proportion = count / sum(count))

p_symp2 <- ggplot(vls3 %>% dplyr::filter(test_calendar_day<40, test_calendar_day>15, test_day==1), aes(x=test_calendar_day)) +
  geom_col(aes(y=proportion, fill=inf_test_delay)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  geom_vline(aes(xintercept=25.5), linetype="dotted") +
  scale_fill_manual(values=c("darkgoldenrod1", "darkgoldenrod2", "darkgoldenrod3", "darkgoldenrod4")) +
  labs(x=NULL, y="Proportion") #+
  #ggtitle("Proportion tested by delay from infection to symptom onset")
  #facet_wrap(~test_day)

vls4 <- vls2 %>% 
  group_by(test_calendar_day, test_day, vl) %>%
  summarise(count = n()) %>%
  group_by(test_calendar_day, test_day) %>%
  mutate(proportion = count / sum(count))

p_vl <- ggplot(vls4 %>% dplyr::filter(test_calendar_day<40, test_calendar_day>15, test_day==1), aes(x=test_calendar_day)) +
  geom_col(aes(y=proportion, fill=vl)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_vline(aes(xintercept=25.5), linetype="dotted") +
  scale_fill_manual(values=c("darkolivegreen1", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4")) +
  labs(x="day", y="Proportion", fill="Viral load") #+
  #ggtitle("Proportion of viral loads recorded by calendar day")
  #facet_wrap(~test_day)

COMB <- ggpubr::ggarrange(p_inf, p_symp, p_symp_test1, p_symp2, p_vl, nrow=5, ncol=1, align="hv", labels="AUTO")
ggsave("images/epidemic_phase_bias_schematic.png", COMB, height=12, width=8)


## methods image
output_dem <- readRDS("alpha/alpha_N_thresh_peak_Cov_ignored_7,8_corr6_bounded_a_b_vg=6_340000_2024-03-12_mult.rds")

incidence <- data.frame(dates = as.Date(dimnames(output_dem$x$data_array)[[2]]), 
                        symptoms=infecteds_generator(parameters=output_dem$MCMC_output[340000,],
                                            knots=output_dem$x$knots, population=output_dem$x$population,
                                            max_day=output_dem$x$max_day, form=output_dem$x$form))

incidence_plot <- ggplot(incidence[30:nrow(incidence),], aes(x=dates, symptoms)) +
  geom_col(fill="lightpink2", color="transparent") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  labs(x="Date", y="Symptom onset incidence")

#ggsave("images/methods_incidence_plot.png", incidence_plot, width=5, height=5)

cond_test_prob <- c(output_dem$MCMC_output[340000,c("test0", "test1", "test2","test3", "test4", "test5")],test6=1)
daily_test_prob <- vector(length=length(cond_test_prob))
daily_test_prob[1] <- cond_test_prob[1]
for(i in 2:length(daily_test_prob)) daily_test_prob[i] <- prod(cond_test_prob[i], 1-(cond_test_prob[0:(i-1)]))

test <- data.frame(test_days=0:(dim(output_dem$x$data_array)[3]-1),
                   test_prob=daily_test_prob)

test_plot <- ggplot(test, aes(x=test_days, y=test_prob)) +
  geom_col(fill="cadetblue1") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  labs(x="Days after onset", y="Probability of testing")

#ggsave("images/methods_test_plot.png", test_plot, width=5, height=5)

output_dem$MCMC_output[315000,c("a_sigma", "b_sigma",  "inc2")] <- c(0.4, 0.5 ,1)

p_array <- p_array_func(parameters=output_dem$MCMC_output[315000,],
                        output_dem$x$knots, vk=1, vg=output_dem$x$vg, max_day=output_dem$x$max_day,
                        population=output_dem$x$population, output_dem$x$data_array, output_dem$x$test_pop,
                        output_dem$x$ncores, output_dem$x$form, output_dem$x$symp_delay_lim, stoch=0.5) 

q_i_k <- p_array[[2]] %>%
  as.data.frame() %>%
  dplyr::rename("1"="negative") %>%
  mutate(day=0:(nrow(.)-1)) %>%
  tidyr::gather(vl, value, 1:(ncol(.)-1)) %>%
  dplyr::filter(day < 7) %>%
  mutate(day=factor(day, levels=0:40),
         value=round(value*1000), 
         vl=as.numeric(vl)) %>%
  tidyr::uncount(value) 

summary_stats <- q_i_k %>%
  group_by(day) %>%
  summarize(
    median = median(vl),
    q1 = quantile(vl, 0.25),
    q3 = quantile(vl, 0.75)
  )

q_ik_plot <- ggplot(data=q_i_k, aes(x=day)) +
  geom_violin(aes(y=vl), fill="darkgoldenrod1") +
  geom_point(data = summary_stats, aes(y = median), color = "black", size = 3) +
  geom_errorbar(data = summary_stats, aes(ymin = q1, ymax = q3), width = 0.1, color = "black") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  labs(y="Probability density of viral load", x="Day after onset")

q_ik_plot

#ggsave("images/methods_qik_plot.png", q_ik_plot, width=5, height=5)

ggsave("images/methods_plot.png", ggarrange(incidence_plot, test_plot, q_ik_plot, nrow=1), width=10, height=3.5)

N_i_j_k <- p_array[[1]] %>% 
  as.data.frame.table() %>% magrittr::set_colnames(c("vl", "day", "test_day", "freq")) %>%
  dplyr::filter(day==100) %>%
  mutate(test_day=paste0("Test day = ", test_day))

N_ijk_plot <- ggplot(N_i_j_k %>% dplyr::filter(test_day!="Test day = 6"), aes(x=vl, y=freq)) +
  geom_col(aes(fill=test_day)) +
  theme_bw() +
  facet_wrap(~test_day) +
  theme(panel.grid=element_blank(),
        strip.text=element_text(size=12, colour="black"),
        legend.position = "none") +
  labs(x="Viral load", y="Number", fill="Day of test after onset") +
  scale_x_discrete(breaks=seq(2,10,2))

ggsave("images/methods_N_ijk_plot.png", N_ijk_plot, width=7, height=7)
  
#





#####
p_a <- ggplot(vls, aes(x=t)) +
  geom_bar(fill="lightblue") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  labs(x="Calendar time", y="Number of infections")

p_a

p_b <- ggplot(vls, aes(x=calendar_symp_day)) +
  geom_bar(aes(fill=symp_onset_t_rounded)) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  #scale_fill_manual(values=c("darkgoldenrod1", "darkgoldenrod2", "darkgoldenrod3", "darkgoldenrod4")) +
  labs(x="Calendar time", y="Number of people with symptom onset", fill="Days since infection") + 
  theme(
    legend.position = c(.99, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(fill="transparent")
  )

p_b

p_c <- ggplot(vls %>% filter(test_time==3), aes(x=calendar_test_day)) +
  geom_bar(aes(fill=vl_test)) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  scale_fill_manual(values=c("aquamarine1", "aquamarine2", "aquamarine3", "aquamarine4","aliceblue")) +
  labs(x="Calendar time", y="Number of people with symptom onset", fill="Viral load") +
  theme(
    legend.position = c(.99, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(fill="transparent")
  ) 

p_c

vls_prop <- vls %>% 
  group_by(test_time, calendar_symp_day, vl_test) %>%
  summarise(count = n()) %>%
  group_by(calendar_symp_day, test_time) %>%
  mutate(proportion = count / sum(count))

p_c2 <- ggplot(vls_prop %>% filter(calendar_symp_day>15, calendar_symp_day<3), aes(x=calendar_symp_day)) +
  geom_col(aes(y=proportion, fill=vl_test)) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  scale_fill_manual(values=c("aquamarine1", "aquamarine2", "aquamarine3", "aquamarine4","aliceblue")) +
  labs(x="Calendar time", y="Number of people with symptom onset", fill="Viral load") +
  theme(
    legend.position = c(.99, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.background = element_rect(fill="transparent")
  ) +
  facet_wrap(~test_time)

p_c2


p_d <- ggplot(vls %>% filter(calendar_test_day %in% c(24, 30)) %>% mutate(calendar_test_day=factor(calendar_test_day)), aes(x=calendar_test_day)) +
  geom_bar(aes(fill=vl_test)) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  scale_fill_manual(values=c("aquamarine1", "aquamarine2", "aquamarine3", "aquamarine4")) +
  labs(x="Calendar time", y="Number of people with symptom onset")

p_phase_shift <- ggpubr::ggarrange(p_a, p_b, p_c, nrow=1)

ggsave("images/phase_shift_demonstration.png", p_phase_shift, width=15, height=5)

vls <- vl_func(2, 1, 4, 0, t=seq(0,6))
plot(vls, x=seq(0,6))
p_symp <- c(0, 0.25, 0.25, 0.25, 0.25)


