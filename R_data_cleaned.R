library(dplyr)

REACT_2_5 <- read.csv("data/REACT_2_5_aggregated.csv") %>% mutate(round = as.numeric(substr(round,7,nchar(round))), date = as.Date(date, format="%d/%m/%Y"))
REACT_6_10 <- read.csv("data/REACT_6_10_aggregated.csv")[,1:7] %>% mutate(round = as.numeric(substr(round,7,nchar(round))), date = as.Date(date, format="%d/%m/%Y"))
REACT_11_19 <- read.csv("data/REACT_11_19_aggregated.csv")[2:8] %>% magrittr::set_colnames(colnames(REACT_6_10)) %>% mutate(date=as.Date(date))

strain_start <- c(wt="2020-03-16", alpha="2020-12-31", delta="2021-05-25",omicron="2021-12-21") 
strain_end <- c(wt="2020-11-24", alpha="2021-04-30", delta="2021-12-08", omicron="2023-03-01")

REACT_df <- rbind(REACT_11_19, REACT_6_10, REACT_2_5) %>% 
  filter(is.na(date)==FALSE, is.na(result_ipsos)==FALSE) %>%
  mutate(n_gene_ct_value = ifelse(n_gene_ct_value==0 | n_gene_ct_value> 40, 40, n_gene_ct_value)) %>%
  filter(n_gene_ct_value > 6) %>% 
  mutate(n_gene_vl_value = round(2*pmin(13.574-0.3011*n_gene_ct_value,11))/2) %>%
  dplyr::arrange(date) %>% 
  group_by(date, n_gene_vl_value) %>% 
  summarise(count=as.numeric(n())) %>% rowwise() %>%
  mutate(strain=names(strain_start)[which(date>=strain_start & date <= strain_end)[1]]) %>%
  ungroup()

ggplot(REACT_df) +
  geom_bar(aes(x=n_gene_vl_value)) +
  theme_bw()

ggplot(REACT_df %>% group_by(date, strain) %>% filter(n_gene_vl_value!=1.5) %>% summarise(count=sum(count)), aes(x=date, y=count)) +
  geom_point(aes(color=strain)) +
  theme_bw()

data_prepare <- function(REACT_df=REACT_df, strain_interest){
  if(strain_interest != "all") REACT_df <- REACT_df %>% filter(strain==strain_interest)
  
  dates <- as.Date(min(REACT_df$date):(max(REACT_df$date)+30))
  
  REACT_raw_combined <- REACT_df %>% mutate(date_number=as.numeric(date-min(date)+31))
  
  days <- 1:max(REACT_raw_combined$date_number)
  
  vl_values <- sort(unique(REACT_df$n_gene_vl_value))
  
  REACT_data_matrix <- matrix(NA, nrow=length(vl_values), ncol=length(days))
  
  skeleton <- data.frame(n_gene_vl_value = vl_values)
  
  for(j in 1:length(days)){
    print(j)
    hold <- REACT_raw_combined %>% filter(date_number==days[j])
    hold2 <- left_join(skeleton, hold, "n_gene_vl_value")
    
    REACT_data_matrix[,j] <- hold2 %>% dplyr::select(count) %>% mutate(count=ifelse(is.na(count), 0, count)) %>% pull()
  } 
  
  colnames(REACT_data_matrix) <- paste(dates)
  rownames(REACT_data_matrix) <- paste0("vl_",vl_values)
  
  #saveRDS(REACT_data_matrix, paste0("REACT_data_matrix_",strain_interest,"2.rds"))
  saveRDS(REACT_data_matrix, paste0("data/REACT_data_matrix_",strain_interest,".rds"))
  return(REACT_data_matrix)
  
}

data_prepare(REACT_df=REACT_df, "all")
data_prepare(REACT_df=REACT_df, "alpha")
data_prepare(REACT_df=REACT_df, "delta")
data_prepare(REACT_df=REACT_df, "omicron")
data_prepare(REACT_df=REACT_df, "wt")

REACT_data_plot <- rbind(REACT_11_19, REACT_6_10, REACT_2_5) %>% 
  filter(is.na(date)==FALSE, is.na(result_ipsos)==FALSE) %>%
  dplyr::mutate(n_gene_ct_value = ifelse(n_gene_ct_value==0, 40, round(n_gene_ct_value,0))) %>%
  dplyr::mutate(n_gene_ct_value = pmin(40,round(as.numeric(n_gene_ct_value), 0))) %>% 
  arrange(date)  %>% 
  group_by(date, n_gene_ct_value) %>% 
  summarise(count=as.numeric(n())) %>% rowwise() %>%
  mutate(strain=names(strain_start)[which(date>=strain_start & date <= strain_end)[1]]) %>% 
  mutate(strain=ifelse(is.na(strain), "interim", strain),
         strain=factor(strain, levels=c("wild_type", "alpha", "delta", "omicron", "interim"))) %>%
  filter(is.na(n_gene_ct_value) == FALSE) %>% rowwise() %>%
  ungroup() 

REACT_ct_dist_plot <- ggplot(REACT_data_plot %>% filter(n_gene_ct_value < 40), aes(x=n_gene_ct_value, y=count)) +
  geom_col(aes(fill=strain)) +
  theme_classic() +
  labs(y="Count", x="Ct value") +
  scale_x_reverse() +
  theme(legend.title = element_blank())

ggsave(filename=paste0("images/REACT_ct_dist_plot.png"), plot = REACT_ct_dist_plot, device = "png", width=10, height=5.5)

REACT_data_plot$count %>% sum()
REACT_data_plot %>% filter(n_gene_ct_value<40) %>% select(count) %>% sum()

positives_data <- REACT_data_plot %>% mutate(positive=ifelse(n_gene_ct_value>=35, "negative", "positive")) %>%
  group_by(date, positive, strain) %>% summarise(count=sum(count)) %>% 
  tidyr::spread(positive, count) %>%
  mutate(positive=ifelse(is.na(positive), 0, positive),
         total=positive+negative,
         proportion=positive/(total), 
         strain=factor(strain, levels=c("wild_type", "alpha", "delta", "omicron", "interim")))

REACT_positives_plot <- ggplot(positives_data, aes(x=date, y=proportion)) +
  scale_x_date(labels = scales::date_format("%b-%y")) +
  geom_col(aes(fill=strain)) +
  theme_classic() +
  labs(y="Proportion\npositive") +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(angle=0, vjust=0.5),
        legend.title = element_blank())

ggsave(filename=paste0("images/REACT_positives_plot.png"), plot = REACT_positives_plot, device = "png", width=10, height=5)

positives_data$proportion %>% max()


REACT_totals_plot <- ggplot(positives_data, aes(x=date, y=total, fill=strain)) +
  geom_col() +
  theme_classic() +
  labs(y="Total\nsamples") +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(angle=0, vjust=0.5),
        legend.title = element_blank())

ggsave(filename=paste0("images/REACT_totals_plot.png"), plot = REACT_totals_plot, device = "png", width=10, height=5)

boxplot_data <- data.frame(date=as.Date(rep(REACT_data_plot$date, REACT_data_plot$count)),
                ct_value=rep(REACT_data_plot$n_gene_ct_value, REACT_data_plot$count),
                strain=rep(REACT_data_plot$strain, REACT_data_plot$count)) %>%
  mutate(week=as.numeric(floor((date-min(date))/7))) %>%
  filter(ct_value != 40) %>%
  mutate(month=format(date,"%b-%y"),
         week_strain=paste(week, strain)) %>%
  mutate(strain=factor(strain, levels=c("wild_type", "alpha", "delta", "omicron", "interim")))

boxplot_plot2 <- ggplot(boxplot_data, aes(x=date, y=ct_value, group=week_strain)) + 
  scale_x_date(labels = scales::date_format("%b-%y")) +
  geom_boxplot(aes(color=strain)) +
  theme_classic() + 
  labs(y="Ct value") +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(angle=0, vjust=0.5))

ggsave(filename=paste0("images/boxplot_plot2.png"), plot=boxplot_plot2, device="png", width=10, height=5.5)

growth_ct <- left_join(positives_data, boxplot_data %>% group_by(date, strain) %>% summarise(mean_ct_value=mean(ct_value)), by="date") %>% ungroup() %>%
  mutate(roll_mean = rollmean(proportion, 7, fill=NA),
         change = roll_mean/lag(roll_mean,7))

ggplot(growth_ct, aes(x=mean_ct_value, y=change)) + geom_point()
