library(dplyr)
library(ggplot2)

ct_data_raw <- readRDS("ct_cleaned.rds")

wt_delta_cutoff <- as.Date("2021-03-05")
alpha_omicron_cutoff <- as.Date("2021-06-23")

ct_data0 <- ct_data_raw %>%
  mutate(strain=ifelse(specimen_date<wt_delta_cutoff & sneg==0, "wt",
                       ifelse(specimen_date>=wt_delta_cutoff & sneg==0, "delta",
                              ifelse(specimen_date<alpha_omicron_cutoff & sneg==1, "alpha",
                                     ifelse(specimen_date>=alpha_omicron_cutoff & sneg==1, "omicron", "NA")))))

strain_start_raw <- setNames(ct_data0 %>% group_by(specimen_date, strain) %>% summarise(count=n()) %>% filter(count > 50) %>% ungroup() %>% group_by(strain) %>% filter(specimen_date==min(specimen_date)) %>% dplyr::select(specimen_date) %>% pull(), c("wt", "alpha", "delta", "omicron"))
strain_end <- setNames(ct_data0 %>% group_by(specimen_date, strain) %>% summarise(count=n()) %>% filter(count > 50) %>% ungroup() %>% group_by(strain) %>% filter(specimen_date==max(specimen_date)) %>% dplyr::select(specimen_date) %>% pull(), c("wt", "alpha", "delta", "omicron"))
strain_start <- strain_start_raw - (7-as.numeric(strain_end-strain_start_raw)%%7)

ct_data <- ct_data_raw %>% 
  mutate(vl=pmin(round((13.698-0.328*p2ch2cq)*2)/2,11)) %>%
  mutate(strain=ifelse(sneg==0 & specimen_date>=strain_start[1] & specimen_date<=strain_end[1], "wt",
                       ifelse(sneg==1 & specimen_date>=strain_start[2] & specimen_date<=strain_end[2], "alpha",
                              ifelse(sneg==0 & specimen_date>=strain_start[3] & specimen_date<=strain_end[3], "delta",
                                     ifelse(sneg==1 & specimen_date>=strain_start[4]+161 & specimen_date<=strain_end[4], "omicron", NA)))))

start_dates <- ct_data %>% group_by(strain) %>% summarise(date=min(onsetdate))

strain_start_final <- setNames(ct_data %>% group_by(strain) %>% 
   filter(specimen_date==min(specimen_date)) %>% 
   group_by(strain, specimen_date) %>% 
   summarise(count=n()) %>% dplyr::select(specimen_date) %>% pull(), c("alpha", "delta", "omicron", "wt"))
 
strain_end_final <- setNames(ct_data %>% group_by(strain) %>% 
   filter(specimen_date==max(specimen_date)) %>% 
   group_by(strain, specimen_date) %>% 
   summarise(count=n()) %>% dplyr::select(specimen_date) %>% pull(), c("alpha", "delta", "omicron", "wt"))

all_start_final <- c(all=min(ct_data$specimen_date, na.rm=T))
all_end_final <- c(all=max(ct_data$specimen_date, na.rm=T))

saveRDS(list(start=c(all_start_final, strain_start_final), end=c(all_end_final, strain_end_final)), "strain_dates.rds")

ct_data

ct_plot <- ct_data %>% group_by(strain, specimen_date) %>% summarise(count=n())
  
pillar2_cases <- ggplot(ct_plot %>% filter(strain!="NA"), aes(x=specimen_date, y=count)) + geom_line(aes(color=strain)) + theme_classic() +
  theme(legend.title=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(angle=0, vjust=0.5)) +
  labs(y="Count")

ggsave(filename=paste0("images/pillar2_cases.png"), plot = pillar2_cases, device = "png", width=20/3, height=10/3)

ct_data_boxplot <- ct_data %>% mutate(week=floor((onsetdate-min(onsetdate, na.rm=T))/7),week_strain=paste(strain,week))

ggplot(ct_data_boxplot %>% filter(is.na(strain)==F), aes(x=onsetdate, y=vl, group=week)) +
  geom_boxplot() +
  #geom_boxplot(aes(color=strain)) +
  #facet_wrap(~strain, scales="free") +
  scale_x_date(labels = scales::date_format("%b-%y")) +
  theme_classic()

thresh <- 1.2
cases_data <- read.csv("data/data_2023-Jun-15.csv") %>% mutate(date=as.Date(date)) %>% arrange(date) %>% dplyr::select(date, newCasesBySpecimenDate) %>%
  mutate(rolling_cases = rollmean(newCasesBySpecimenDate,7, fill=NA)) %>% mutate(weekly_change=lead(rolling_cases,3)/lag(rolling_cases, 3)) %>%
  mutate(growth_rate=ifelse(weekly_change > thresh, "Growing epidemic", ifelse(weekly_change > 1/thresh, "medium", "Shrinking epidemic"))) %>% dplyr::select(date, weekly_change, growth_rate) %>%
  rename(specimen_date=date)


ct_from_onset <- ct_data %>% dplyr::select(onsetdate, specimen_date, onset2specimen, strain, p2ch1cq) %>% filter(is.na(strain)==F) %>% 
  group_by(strain, onset2specimen, specimen_date) %>% mutate(ctraw=0.14298+1.07095*p2ch1cq, strain=factor(strain, levels=c("wt", "alpha", "delta", "omicron"))) %>% 
  filter(ctraw < 40, is.na(ctraw)==F, onset2specimen <= 10) %>% left_join(cases_data, by="specimen_date") %>% ungroup() 

ct_from_onset_strain <- ct_from_onset %>% group_by(strain, onset2specimen) %>% summarise(ct=mean(ctraw), se=sd(ctraw)/sqrt(n())) %>% mutate(se_low=ct-se, se_high=ct+se)
ct_from_onset_growth <- ct_from_onset %>% group_by(strain, onset2specimen, growth_rate) %>% summarise(ct=mean(ctraw), se=sd(ctraw)/sqrt(n())) %>% mutate(se_low=ct-se, se_high=ct+se)

ct_from_onset_strain_plot <- ggplot(ct_from_onset_strain, aes(x=onset2specimen)) +
  geom_line(aes(y=ct, color=strain)) +
  geom_ribbon(aes(ymin=se_low, ymax=se_high, fill=strain), alpha=0.5) +
  labs(x="Time from onset (days)", y="Ct\nvalue") +
  theme_classic() +
  theme(legend.title=element_blank(),
        axis.title.y = element_text(angle=0, vjust=0.5)) +
  scale_x_continuous(breaks=seq(0,10))

ct_from_onset_strain_plot

ggsave(filename=paste0("images/ct_from_onset_variant.png"), plot = ct_from_onset_strain_plot, device = "png", width=5, height=10/3)

ct_from_onset_growth_plot <- ggplot(ct_from_onset_growth %>% filter(growth_rate != "medium"), aes(x=onset2specimen, color=growth_rate, fill=growth_rate)) +
  geom_line(aes(y=ct)) +
  geom_ribbon(aes(ymin=se_low, ymax=se_high), alpha=0.2) +
  facet_wrap(~strain, nrow=1) +
  theme_classic() +
  labs(x="Time from onset (days)", y="Ct\nvalue") +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(angle=0, vjust=0.5)) +
  scale_x_continuous(breaks=seq(0,10))
  
ct_from_onset_growth_plot

ggsave(filename=paste0("images/ct_from_onset_growth.png"), plot = ct_from_onset_growth_plot, device = "png", width=10, height=10/3)

ggplot(ct_data %>% filter(onsetdate < as.Date("2020-12-01")), aes(x=p2ch1cq, y=p2ch2cq)) +
  geom_point()

## strain specific data array
array_func <- function(ct_df_raw, strain_interest){
  
  if(strain_interest!="all") ct_df <- ct_df_raw %>% filter(strain==strain_interest) 
  if(strain_interest=="all") ct_df <- ct_df_raw %>% filter(specimen_date<=max(specimen_date,na.rm=T)-3) 
  
  ct_df <- ct_df %>% mutate(day_number=as.numeric(specimen_date-min(ct_df$specimen_date, na.rm=T))+35)
  
  max_day <- max(ct_df$day_number, na.rm=T)
  
  ct_values <- sort(unique(ct_df_raw$vl))
  day_numbers <- 1:max_day
  dates <- seq(min(ct_df$specimen_date, na.rm=T)-34,max(ct_df$specimen_date, na.rm=T),1)
  onset2symp <- 0:13
  
  data_array <- array(NA, dim=c(length(ct_values), length(day_numbers), length(onset2symp)))
  
  skeleton <- data.frame(day_number=day_numbers)
  
  for(i in 1:dim(data_array)[1]){ 
    ct_data_filtered <- ct_df %>% filter(vl==i/2)
    for(k in 1:dim(data_array)[3]){
      filtered <- left_join(skeleton, ct_data_filtered %>% filter(onset2specimen==k-1) %>% group_by(day_number, vl, onset2specimen) %>% summarise(count=n()), by="day_number")
      data_array[i,,k] <- filtered$count
    }
  }
  
  data_array[is.na(data_array)] = 0
  
  dimnames(data_array) <- list(paste0("vl_",ct_values), paste0(dates),paste0("dtt_",onset2symp))
  
  saveRDS(data_array, paste0("data/data_array_vl_",strain_interest,"_N",".rds"))
  
  return(data_array)
}

array_func(ct_df_raw=ct_data, strain_interest="all")
array_func(ct_data, strain_interest="wt")
array_func(ct_data, strain_interest="alpha")
array_func(ct_data, strain_interest="delta")
array_func(ct_df=ct_data, strain_interest="omicron")

readRDS("data/data_array_vl_all_N.rds") %>% ncol() %% 7
readRDS("data/data_array_vl_wt_N.rds") %>% ncol() 
readRDS("data/data_array_vl_alpha_N.rds") %>% ncol() 
readRDS("data/data_array_vl_delta_N.rds") %>% ncol() 
readRDS("data/data_array_vl_omicron_N.rds") %>% ncol() 
