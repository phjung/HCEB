############################
## 1. GEO DATA SETTING
############################

library("tidyverse")
load("data/RData/data_TB.RData")

############################
## 1.4 TEEN BIRTH RATES AND HYPERTENSION RATES
############################
############################
## 1.4.0 CV: Birth Count, Fertile Women, Total Women
############################

census.tract@data <- census.tract@data %>%
  mutate(CV_BC_2010=(MOE_BC_2010/1.645)/(0.1+BC_2010), # Birth Count
         CV_BC_2011=(MOE_BC_2010/1.645)/(0.1+BC_2011),
         CV_BC_2012=(MOE_BC_2010/1.645)/(0.1+BC_2012),
         CV_BC_2013=(MOE_BC_2010/1.645)/(0.1+BC_2013),
         CV_BC_2014=(MOE_BC_2010/1.645)/(0.1+BC_2014),
         
         CV_WF_2010=(MOE_WF_2010/1.645)/(0.1+WF_2010), # Fertile Women
         CV_WF_2011=(MOE_WF_2010/1.645)/(0.1+WF_2011),
         CV_WF_2012=(MOE_WF_2010/1.645)/(0.1+WF_2012),
         CV_WF_2013=(MOE_WF_2010/1.645)/(0.1+WF_2013),
         CV_WF_2014=(MOE_WF_2010/1.645)/(0.1+WF_2014),
         
         CV_WT_2010=(MOE_WT_2010/1.645)/(0.1+WT_2010), # Total Women
         CV_WT_2011=(MOE_WT_2010/1.645)/(0.1+WT_2011),
         CV_WT_2012=(MOE_WT_2010/1.645)/(0.1+WT_2012),
         CV_WT_2013=(MOE_WT_2010/1.645)/(0.1+WT_2013),
         CV_WT_2014=(MOE_WT_2010/1.645)/(0.1+WT_2014)
  )

############################
## 1.4.1 Teen Birth Rate by Total Brith
############################

### Raw rate

census.tract@data <- census.tract@data %>%
  mutate(rate_raw_TB_2010=TB_2010/(0.1+BC_2010), # Raw rate - TB
         rate_raw_TB_2011=TB_2011/(0.1+BC_2011),
         rate_raw_TB_2012=TB_2012/(0.1+BC_2012),
         rate_raw_TB_2013=TB_2013/(0.1+BC_2013),
         rate_raw_TB_2014=TB_2014/(0.1+BC_2014),
         
         rate_std_raw_TB_2010 = rate_raw_TB_2010 - sum(TB_2010)/(sum(0.1+BC_2010)), # Std Raw rate
         rate_std_raw_TB_2011 = rate_raw_TB_2011 - sum(TB_2011)/(sum(0.1+BC_2011)),
         rate_std_raw_TB_2012 = rate_raw_TB_2012 - sum(TB_2012)/(sum(0.1+BC_2012)),
         rate_std_raw_TB_2013 = rate_raw_TB_2013 - sum(TB_2013)/(sum(0.1+BC_2013)),
         rate_std_raw_TB_2014 = rate_raw_TB_2014 - sum(TB_2014)/(sum(0.1+BC_2014)))

## Weight 

weight_BC <- census.tract@data %>%
  mutate(n_TB_2010 = sum(0.1+BC_2010),
         n_TB_2011 = sum(0.1+BC_2011),
         n_TB_2012 = sum(0.1+BC_2012),
         n_TB_2013 = sum(0.1+BC_2013),
         n_TB_2014 = sum(0.1+BC_2014)) %>%
  
  mutate(m_TB_2010 = sum(TB_2010)/n_TB_2010,
         m_TB_2011 = sum(TB_2011)/n_TB_2011,
         m_TB_2012 = sum(TB_2012)/n_TB_2012,
         m_TB_2013 = sum(TB_2013)/n_TB_2013,
         m_TB_2014 = sum(TB_2014)/n_TB_2014) %>%
  
  mutate(s2_TB_2010 = sum((0.1+BC_2010)*(rate_raw_TB_2010-m_TB_2010)^2)/n_TB_2010,
         s2_TB_2011 = sum((0.1+BC_2011)*(rate_raw_TB_2011-m_TB_2011)^2)/n_TB_2011,
         s2_TB_2012 = sum((0.1+BC_2012)*(rate_raw_TB_2012-m_TB_2012)^2)/n_TB_2012,
         s2_TB_2013 = sum((0.1+BC_2013)*(rate_raw_TB_2013-m_TB_2013)^2)/n_TB_2013,
         s2_TB_2014 = sum((0.1+BC_2014)*(rate_raw_TB_2014-m_TB_2014)^2)/n_TB_2014) %>%
  
  mutate(A_TB_2010 = s2_TB_2010-m_TB_2010/(n_TB_2010/233),
         A_TB_2011 = s2_TB_2011-m_TB_2011/(n_TB_2011/233),
         A_TB_2012 = s2_TB_2012-m_TB_2012/(n_TB_2012/233),
         A_TB_2013 = s2_TB_2013-m_TB_2013/(n_TB_2013/233),
         A_TB_2014 = s2_TB_2014-m_TB_2014/(n_TB_2014/233),
         
         B_TB_2010 = m_TB_2010/(0.1+BC_2010),
         B_TB_2011 = m_TB_2011/(0.1+BC_2011),
         B_TB_2012 = m_TB_2012/(0.1+BC_2012),
         B_TB_2013 = m_TB_2013/(0.1+BC_2013),
         B_TB_2014 = m_TB_2014/(0.1+BC_2014),
         
         C_TB_2010 = rate_raw_TB_2010^2*CV_BC_2010^2,
         C_TB_2011 = rate_raw_TB_2011^2*CV_BC_2011^2,
         C_TB_2012 = rate_raw_TB_2012^2*CV_BC_2012^2,
         C_TB_2013 = rate_raw_TB_2013^2*CV_BC_2013^2,
         C_TB_2014 = rate_raw_TB_2014^2*CV_BC_2014^2) %>%
  
  mutate(w_EB_TB_2010 = A_TB_2010/(A_TB_2010+B_TB_2010),
         w_EB_TB_2011 = A_TB_2011/(A_TB_2011+B_TB_2011),
         w_EB_TB_2012 = A_TB_2012/(A_TB_2012+B_TB_2012),
         w_EB_TB_2013 = A_TB_2013/(A_TB_2013+B_TB_2013),
         w_EB_TB_2014 = A_TB_2014/(A_TB_2014+B_TB_2014),
         
         w_EB2_TB_2010 = A_TB_2010/(A_TB_2010+B_TB_2010+C_TB_2010),
         w_EB2_TB_2011 = A_TB_2011/(A_TB_2011+B_TB_2011+C_TB_2011),
         w_EB2_TB_2012 = A_TB_2012/(A_TB_2012+B_TB_2012+C_TB_2012),
         w_EB2_TB_2013 = A_TB_2013/(A_TB_2013+B_TB_2013+C_TB_2013),
         w_EB2_TB_2014 = A_TB_2014/(A_TB_2014+B_TB_2014+C_TB_2014)) %>%
  
  dplyr::select(m_TB_2010,
                m_TB_2011,
                m_TB_2012,
                m_TB_2013,
                m_TB_2014,
                w_EB_TB_2010,
                w_EB_TB_2011,
                w_EB_TB_2012,
                w_EB_TB_2013,
                w_EB_TB_2014,
                w_EB2_TB_2010,
                w_EB2_TB_2011,
                w_EB2_TB_2012,
                w_EB2_TB_2013,
                w_EB2_TB_2014)

## EB Rate: EB rate, EB2 rate

census.tract@data <- census.tract@data %>%
  mutate(rate_EB_TB_2010 = weight_BC$m_TB_2010+rate_std_raw_TB_2010*weight_BC$w_EB_TB_2010,
         rate_EB_TB_2011 = weight_BC$m_TB_2011+rate_std_raw_TB_2011*weight_BC$w_EB_TB_2011,
         rate_EB_TB_2012 = weight_BC$m_TB_2012+rate_std_raw_TB_2012*weight_BC$w_EB_TB_2012,
         rate_EB_TB_2013 = weight_BC$m_TB_2013+rate_std_raw_TB_2013*weight_BC$w_EB_TB_2013,
         rate_EB_TB_2014 = weight_BC$m_TB_2014+rate_std_raw_TB_2014*weight_BC$w_EB_TB_2014,
         
         rate_std_EB_TB_2010 = rate_std_raw_TB_2010*weight_BC$w_EB_TB_2010,
         rate_std_EB_TB_2011 = rate_std_raw_TB_2011*weight_BC$w_EB_TB_2011,
         rate_std_EB_TB_2012 = rate_std_raw_TB_2012*weight_BC$w_EB_TB_2012,
         rate_std_EB_TB_2013 = rate_std_raw_TB_2013*weight_BC$w_EB_TB_2013,
         rate_std_EB_TB_2014 = rate_std_raw_TB_2014*weight_BC$w_EB_TB_2014,
         
         rate_EB2_TB_2010 = weight_BC$m_TB_2010+rate_std_raw_TB_2010*weight_BC$w_EB2_TB_2010,
         rate_EB2_TB_2011 = weight_BC$m_TB_2011+rate_std_raw_TB_2011*weight_BC$w_EB2_TB_2011,
         rate_EB2_TB_2012 = weight_BC$m_TB_2012+rate_std_raw_TB_2012*weight_BC$w_EB2_TB_2012,
         rate_EB2_TB_2013 = weight_BC$m_TB_2013+rate_std_raw_TB_2013*weight_BC$w_EB2_TB_2013,
         rate_EB2_TB_2014 = weight_BC$m_TB_2014+rate_std_raw_TB_2014*weight_BC$w_EB2_TB_2014,
         
         rate_std_EB2_TB_2010 = rate_std_raw_TB_2010*weight_BC$w_EB2_TB_2010,
         rate_std_EB2_TB_2011 = rate_std_raw_TB_2011*weight_BC$w_EB2_TB_2011,
         rate_std_EB2_TB_2012 = rate_std_raw_TB_2012*weight_BC$w_EB2_TB_2012,
         rate_std_EB2_TB_2013 = rate_std_raw_TB_2013*weight_BC$w_EB2_TB_2013,
         rate_std_EB2_TB_2014 = rate_std_raw_TB_2014*weight_BC$w_EB2_TB_2014)

############################
## 1.4.2 Teen Birth Rate by Fertile Women
############################

### Raw rate

census.tract@data <- census.tract@data %>%
  mutate(rate_raw_WF_2010=TB_2010/(0.1+WF_2010), # Raw rate - TB
         rate_raw_WF_2011=TB_2011/(0.1+WF_2011),
         rate_raw_WF_2012=TB_2012/(0.1+WF_2012),
         rate_raw_WF_2013=TB_2013/(0.1+WF_2013),
         rate_raw_WF_2014=TB_2014/(0.1+WF_2014),
         
         rate_std_raw_WF_2010 = rate_raw_WF_2010 - sum(TB_2010)/(sum(0.1+WF_2010)), # Std Raw rate
         rate_std_raw_WF_2011 = rate_raw_WF_2011 - sum(TB_2011)/(sum(0.1+WF_2011)),
         rate_std_raw_WF_2012 = rate_raw_WF_2012 - sum(TB_2012)/(sum(0.1+WF_2012)),
         rate_std_raw_WF_2013 = rate_raw_WF_2013 - sum(TB_2013)/(sum(0.1+WF_2013)),
         rate_std_raw_WF_2014 = rate_raw_WF_2014 - sum(TB_2014)/(sum(0.1+WF_2014)))

## Weight 

weight_WF <- census.tract@data %>%
  mutate(n_WF_2010 = sum(0.1+WF_2010),
         n_WF_2011 = sum(0.1+WF_2011),
         n_WF_2012 = sum(0.1+WF_2012),
         n_WF_2013 = sum(0.1+WF_2013),
         n_WF_2014 = sum(0.1+WF_2014)) %>%
  
  mutate(m_WF_2010 = sum(TB_2010)/n_WF_2010,
         m_WF_2011 = sum(TB_2011)/n_WF_2011,
         m_WF_2012 = sum(TB_2012)/n_WF_2012,
         m_WF_2013 = sum(TB_2013)/n_WF_2013,
         m_WF_2014 = sum(TB_2014)/n_WF_2014) %>%
  
  mutate(s2_WF_2010 = sum((0.1+WF_2010)*(rate_raw_WF_2010-m_WF_2010)^2)/n_WF_2010,
         s2_WF_2011 = sum((0.1+WF_2011)*(rate_raw_WF_2011-m_WF_2011)^2)/n_WF_2011,
         s2_WF_2012 = sum((0.1+WF_2012)*(rate_raw_WF_2012-m_WF_2012)^2)/n_WF_2012,
         s2_WF_2013 = sum((0.1+WF_2013)*(rate_raw_WF_2013-m_WF_2013)^2)/n_WF_2013,
         s2_WF_2014 = sum((0.1+WF_2014)*(rate_raw_WF_2014-m_WF_2014)^2)/n_WF_2014) %>%
  
  mutate(A_WF_2010 = s2_WF_2010-m_WF_2010/(n_WF_2010/233),
         A_WF_2011 = s2_WF_2011-m_WF_2011/(n_WF_2011/233),
         A_WF_2012 = s2_WF_2012-m_WF_2012/(n_WF_2012/233),
         A_WF_2013 = s2_WF_2013-m_WF_2013/(n_WF_2013/233),
         A_WF_2014 = s2_WF_2014-m_WF_2014/(n_WF_2014/233),
         
         B_WF_2010 = m_WF_2010/(0.1+WF_2010),
         B_WF_2011 = m_WF_2011/(0.1+WF_2011),
         B_WF_2012 = m_WF_2012/(0.1+WF_2012),
         B_WF_2013 = m_WF_2013/(0.1+WF_2013),
         B_WF_2014 = m_WF_2014/(0.1+WF_2014),
         
         C_WF_2010 = rate_raw_WF_2010^2*CV_WF_2010^2,
         C_WF_2011 = rate_raw_WF_2011^2*CV_WF_2011^2,
         C_WF_2012 = rate_raw_WF_2012^2*CV_WF_2012^2,
         C_WF_2013 = rate_raw_WF_2013^2*CV_WF_2013^2,
         C_WF_2014 = rate_raw_WF_2014^2*CV_WF_2014^2) %>%
  
  mutate(w_EB_WF_2010 = A_WF_2010/(A_WF_2010+B_WF_2010),
         w_EB_WF_2011 = A_WF_2011/(A_WF_2011+B_WF_2011),
         w_EB_WF_2012 = A_WF_2012/(A_WF_2012+B_WF_2012),
         w_EB_WF_2013 = A_WF_2013/(A_WF_2013+B_WF_2013),
         w_EB_WF_2014 = A_WF_2014/(A_WF_2014+B_WF_2014),
         
         w_EB2_WF_2010 = A_WF_2010/(A_WF_2010+B_WF_2010+C_WF_2010),
         w_EB2_WF_2011 = A_WF_2011/(A_WF_2011+B_WF_2011+C_WF_2011),
         w_EB2_WF_2012 = A_WF_2012/(A_WF_2012+B_WF_2012+C_WF_2012),
         w_EB2_WF_2013 = A_WF_2013/(A_WF_2013+B_WF_2013+C_WF_2013),
         w_EB2_WF_2014 = A_WF_2014/(A_WF_2014+B_WF_2014+C_WF_2014)) %>%
  
  dplyr::select(m_WF_2010,
                m_WF_2011,
                m_WF_2012,
                m_WF_2013,
                m_WF_2014,
                w_EB_WF_2010,
                w_EB_WF_2011,
                w_EB_WF_2012,
                w_EB_WF_2013,
                w_EB_WF_2014,
                w_EB2_WF_2010,
                w_EB2_WF_2011,
                w_EB2_WF_2012,
                w_EB2_WF_2013,
                w_EB2_WF_2014)

## EB Rate: EB rate, EB2 rate

census.tract@data <- census.tract@data %>%
  mutate(rate_EB_WF_2010 = weight_WF$m_WF_2010+rate_std_raw_WF_2010*weight_WF$w_EB_WF_2010,
         rate_EB_WF_2011 = weight_WF$m_WF_2011+rate_std_raw_WF_2011*weight_WF$w_EB_WF_2011,
         rate_EB_WF_2012 = weight_WF$m_WF_2012+rate_std_raw_WF_2012*weight_WF$w_EB_WF_2012,
         rate_EB_WF_2013 = weight_WF$m_WF_2013+rate_std_raw_WF_2013*weight_WF$w_EB_WF_2013,
         rate_EB_WF_2014 = weight_WF$m_WF_2014+rate_std_raw_WF_2014*weight_WF$w_EB_WF_2014,
         
         rate_std_EB_WF_2010 = rate_std_raw_WF_2010*weight_WF$w_EB_WF_2010,
         rate_std_EB_WF_2011 = rate_std_raw_WF_2011*weight_WF$w_EB_WF_2011,
         rate_std_EB_WF_2012 = rate_std_raw_WF_2012*weight_WF$w_EB_WF_2012,
         rate_std_EB_WF_2013 = rate_std_raw_WF_2013*weight_WF$w_EB_WF_2013,
         rate_std_EB_WF_2014 = rate_std_raw_WF_2014*weight_WF$w_EB_WF_2014,
         
         rate_EB2_WF_2010 = weight_WF$m_WF_2010+rate_std_raw_WF_2010*weight_WF$w_EB2_WF_2010,
         rate_EB2_WF_2011 = weight_WF$m_WF_2011+rate_std_raw_WF_2011*weight_WF$w_EB2_WF_2011,
         rate_EB2_WF_2012 = weight_WF$m_WF_2012+rate_std_raw_WF_2012*weight_WF$w_EB2_WF_2012,
         rate_EB2_WF_2013 = weight_WF$m_WF_2013+rate_std_raw_WF_2013*weight_WF$w_EB2_WF_2013,
         rate_EB2_WF_2014 = weight_WF$m_WF_2014+rate_std_raw_WF_2014*weight_WF$w_EB2_WF_2014,
         
         rate_std_EB2_WF_2010 = rate_std_raw_WF_2010*weight_WF$w_EB2_WF_2010,
         rate_std_EB2_WF_2011 = rate_std_raw_WF_2011*weight_WF$w_EB2_WF_2011,
         rate_std_EB2_WF_2012 = rate_std_raw_WF_2012*weight_WF$w_EB2_WF_2012,
         rate_std_EB2_WF_2013 = rate_std_raw_WF_2013*weight_WF$w_EB2_WF_2013,
         rate_std_EB2_WF_2014 = rate_std_raw_WF_2014*weight_WF$w_EB2_WF_2014)

############################
## 1.4.3 Teen Birth Rate by Teen Women (ACS)
############################

### Raw rate

census.tract@data <- census.tract@data %>%
  mutate(rate_raw_WT_2010=TB_2010/(0.1+WT_2010), # Raw rate - TB
         rate_raw_WT_2011=TB_2011/(0.1+WT_2011),
         rate_raw_WT_2012=TB_2012/(0.1+WT_2012),
         rate_raw_WT_2013=TB_2013/(0.1+WT_2013),
         rate_raw_WT_2014=TB_2014/(0.1+WT_2014),
         
         rate_std_raw_WT_2010 = rate_raw_WT_2010 - sum(TB_2010)/(sum(0.1+WT_2010)), # Std Raw rate
         rate_std_raw_WT_2011 = rate_raw_WT_2011 - sum(TB_2011)/(sum(0.1+WT_2011)),
         rate_std_raw_WT_2012 = rate_raw_WT_2012 - sum(TB_2012)/(sum(0.1+WT_2012)),
         rate_std_raw_WT_2013 = rate_raw_WT_2013 - sum(TB_2013)/(sum(0.1+WT_2013)),
         rate_std_raw_WT_2014 = rate_raw_WT_2014 - sum(TB_2014)/(sum(0.1+WT_2014)))

## Weight 

weight_WT <- census.tract@data %>%
  mutate(n_WT_2010 = sum(0.1+WT_2010),
         n_WT_2011 = sum(0.1+WT_2011),
         n_WT_2012 = sum(0.1+WT_2012),
         n_WT_2013 = sum(0.1+WT_2013),
         n_WT_2014 = sum(0.1+WT_2014)) %>%
  
  mutate(m_WT_2010 = sum(TB_2010)/n_WT_2010,
         m_WT_2011 = sum(TB_2011)/n_WT_2011,
         m_WT_2012 = sum(TB_2012)/n_WT_2012,
         m_WT_2013 = sum(TB_2013)/n_WT_2013,
         m_WT_2014 = sum(TB_2014)/n_WT_2014) %>%
  
  mutate(s2_WT_2010 = sum((0.1+WT_2010)*(rate_raw_WT_2010-m_WT_2010)^2)/n_WT_2010,
         s2_WT_2011 = sum((0.1+WT_2011)*(rate_raw_WT_2011-m_WT_2011)^2)/n_WT_2011,
         s2_WT_2012 = sum((0.1+WT_2012)*(rate_raw_WT_2012-m_WT_2012)^2)/n_WT_2012,
         s2_WT_2013 = sum((0.1+WT_2013)*(rate_raw_WT_2013-m_WT_2013)^2)/n_WT_2013,
         s2_WT_2014 = sum((0.1+WT_2014)*(rate_raw_WT_2014-m_WT_2014)^2)/n_WT_2014) %>%
  
  mutate(A_WT_2010 = s2_WT_2010-m_WT_2010/(n_WT_2010/233),
         A_WT_2011 = s2_WT_2011-m_WT_2011/(n_WT_2011/233),
         A_WT_2012 = s2_WT_2012-m_WT_2012/(n_WT_2012/233),
         A_WT_2013 = s2_WT_2013-m_WT_2013/(n_WT_2013/233),
         A_WT_2014 = s2_WT_2014-m_WT_2014/(n_WT_2014/233),
         
         B_WT_2010 = m_WT_2010/(0.1+WT_2010),
         B_WT_2011 = m_WT_2011/(0.1+WT_2011),
         B_WT_2012 = m_WT_2012/(0.1+WT_2012),
         B_WT_2013 = m_WT_2013/(0.1+WT_2013),
         B_WT_2014 = m_WT_2014/(0.1+WT_2014),
         
         C_WT_2010 = rate_raw_WT_2010^2*CV_WT_2010^2,
         C_WT_2011 = rate_raw_WT_2011^2*CV_WT_2011^2,
         C_WT_2012 = rate_raw_WT_2012^2*CV_WT_2012^2,
         C_WT_2013 = rate_raw_WT_2013^2*CV_WT_2013^2,
         C_WT_2014 = rate_raw_WT_2014^2*CV_WT_2014^2
         
         #C_WT_2010 = (A_WT_2010+m_WT_2010)*CV_WT_2010^2,
         #C_WT_2011 = (A_WT_2011+m_WT_2011)*CV_WT_2011^2,
         #C_WT_2012 = (A_WT_2012+m_WT_2012)*CV_WT_2012^2,
         #C_WT_2013 = (A_WT_2013+m_WT_2013)*CV_WT_2013^2,
         #C_WT_2014 = (A_WT_2014+m_WT_2014)*CV_WT_2014^2
         ) %>%
         
  
  
  mutate(w_EB_WT_2010 = A_WT_2010/(A_WT_2010+B_WT_2010),
         w_EB_WT_2011 = A_WT_2011/(A_WT_2011+B_WT_2011),
         w_EB_WT_2012 = A_WT_2012/(A_WT_2012+B_WT_2012),
         w_EB_WT_2013 = A_WT_2013/(A_WT_2013+B_WT_2013),
         w_EB_WT_2014 = A_WT_2014/(A_WT_2014+B_WT_2014),
         
         w_EB2_WT_2010 = A_WT_2010/(A_WT_2010+B_WT_2010+C_WT_2010),
         w_EB2_WT_2011 = A_WT_2011/(A_WT_2011+B_WT_2011+C_WT_2011),
         w_EB2_WT_2012 = A_WT_2012/(A_WT_2012+B_WT_2012+C_WT_2012),
         w_EB2_WT_2013 = A_WT_2013/(A_WT_2013+B_WT_2013+C_WT_2013),
         w_EB2_WT_2014 = A_WT_2014/(A_WT_2014+B_WT_2014+C_WT_2014)) %>%
  
  dplyr::select(m_WT_2010,
                m_WT_2011,
                m_WT_2012,
                m_WT_2013,
                m_WT_2014,
                w_EB_WT_2010,
                w_EB_WT_2011,
                w_EB_WT_2012,
                w_EB_WT_2013,
                w_EB_WT_2014,
                w_EB2_WT_2010,
                w_EB2_WT_2011,
                w_EB2_WT_2012,
                w_EB2_WT_2013,
                w_EB2_WT_2014)

## EB Rate: EB rate, EB2 rate

census.tract@data <- census.tract@data %>%
  mutate(rate_EB_WT_2010 = weight_WT$m_WT_2010+rate_std_raw_WT_2010*weight_WT$w_EB_WT_2010,
         rate_EB_WT_2011 = weight_WT$m_WT_2011+rate_std_raw_WT_2011*weight_WT$w_EB_WT_2011,
         rate_EB_WT_2012 = weight_WT$m_WT_2012+rate_std_raw_WT_2012*weight_WT$w_EB_WT_2012,
         rate_EB_WT_2013 = weight_WT$m_WT_2013+rate_std_raw_WT_2013*weight_WT$w_EB_WT_2013,
         rate_EB_WT_2014 = weight_WT$m_WT_2014+rate_std_raw_WT_2014*weight_WT$w_EB_WT_2014,
         
         rate_std_EB_WT_2010 = rate_std_raw_WT_2010*weight_WT$w_EB_WT_2010,
         rate_std_EB_WT_2011 = rate_std_raw_WT_2011*weight_WT$w_EB_WT_2011,
         rate_std_EB_WT_2012 = rate_std_raw_WT_2012*weight_WT$w_EB_WT_2012,
         rate_std_EB_WT_2013 = rate_std_raw_WT_2013*weight_WT$w_EB_WT_2013,
         rate_std_EB_WT_2014 = rate_std_raw_WT_2014*weight_WT$w_EB_WT_2014,
         
         rate_EB2_WT_2010 = weight_WT$m_WT_2010+rate_std_raw_WT_2010*weight_WT$w_EB2_WT_2010,
         rate_EB2_WT_2011 = weight_WT$m_WT_2011+rate_std_raw_WT_2011*weight_WT$w_EB2_WT_2011,
         rate_EB2_WT_2012 = weight_WT$m_WT_2012+rate_std_raw_WT_2012*weight_WT$w_EB2_WT_2012,
         rate_EB2_WT_2013 = weight_WT$m_WT_2013+rate_std_raw_WT_2013*weight_WT$w_EB2_WT_2013,
         rate_EB2_WT_2014 = weight_WT$m_WT_2014+rate_std_raw_WT_2014*weight_WT$w_EB2_WT_2014,
         
         rate_std_EB2_WT_2010 = rate_std_raw_WT_2010*weight_WT$w_EB2_WT_2010,
         rate_std_EB2_WT_2011 = rate_std_raw_WT_2011*weight_WT$w_EB2_WT_2011,
         rate_std_EB2_WT_2012 = rate_std_raw_WT_2012*weight_WT$w_EB2_WT_2012,
         rate_std_EB2_WT_2013 = rate_std_raw_WT_2013*weight_WT$w_EB2_WT_2013,
         rate_std_EB2_WT_2014 = rate_std_raw_WT_2014*weight_WT$w_EB2_WT_2014)

############################
## 1.4.4 Teen Birth Rate by Teen Women (Census)
############################

### Raw rate

census.tract@data <- census.tract@data %>%
  mutate(rate_raw_WT_DC_2010=TB_2010/(0.1+WT_2010_DC), # Raw rate - TB
         rate_std_raw_WT_DC_2010 = rate_raw_WT_DC_2010 - sum(TB_2010)/(sum(0.1+WT_2010_DC))) # Std Raw rate

## Weight 

weight_WT_DC <- census.tract@data %>%
  mutate(n_WT_DC_2010 = sum(0.1+WT_2010_DC)) %>%
  mutate(m_WT_DC_2010 = sum(TB_2010)/n_WT_DC_2010) %>%
  mutate(s2_WT_DC_2010 = sum((0.1+WT_2010_DC)*(rate_raw_WT_DC_2010-m_WT_DC_2010)^2)/n_WT_DC_2010) %>%
  mutate(A_WT_DC_2010 = s2_WT_DC_2010-m_WT_DC_2010/(n_WT_DC_2010/233),
         B_WT_DC_2010 = m_WT_DC_2010/(0.1+WT_2010_DC)) %>%
  mutate(w_EB_WT_DC_2010 = A_WT_DC_2010/(A_WT_DC_2010+B_WT_DC_2010)) %>%
  dplyr::select(m_WT_DC_2010,
                w_EB_WT_DC_2010)

## EB Rate: EB rate

census.tract@data <- census.tract@data %>%
  mutate(rate_EB_WT_DC_2010 = weight_WT_DC$m_WT_DC_2010+rate_std_raw_WT_DC_2010*weight_WT_DC$w_EB_WT_DC_2010,
         rate_std_EB_WT_DC_2010 = rate_std_raw_WT_DC_2010*weight_WT_DC$w_EB_WT_DC_2010)

############################
save.image("data/RData/data_TB.RData")

