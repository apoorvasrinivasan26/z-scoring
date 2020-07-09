library(tidyverse)
library(sqldf)
library(robustHD)
library(nlme)
library(splines)
library(ggplot2)
require(rms)
library(rms)
library(MASS)
library(faraway)
library(devtools)
library(SnakeCharmR)
library(foreign)


##reading in new data - 10/28/19

aging_eprime_final = read.spss("AFSP ALL SITES EPrime for Zscoring wDemog.sav",to.data.frame=T,
                     use.value.lables=F) %>%
  filter(use == 1) 

dropsubnew = c(202200, 220051, 220116, 20769, 208806)

aging_eprime_final = aging_eprime_final[!((aging_eprime_final$subject) %in% dropsubnew),]


dim(aging_eprime_final) 
# 304 170

length(unique(aging_eprime_final$subject[which(aging_eprime_final$grp == 0)]))
#controls: 93


aging_eprime_final = aging_eprime_final %>%
  janitor::clean_names() %>%
mutate(sex = as.factor(sex),
         grp = as.factor(grp)) %>%
  rename(agesqr = agesq) %>%
         rename(toteduc = educa) 



aging_eprime_final = sqldf("SELECT *,
                           CASE 
                           WHEN age >60 THEN '61-80'
                           WHEN age >= 46 THEN '46-60'
                           WHEN age >= 31 THEN '31-45'
                           WHEN age >= 16 THEN '16-30'
                           END AS age_group
                           FROM aging_eprime_final;
                           ")
dim(aging_eprime_final)

aging_control = aging_eprime_final[which(aging_eprime_final$grp == 0),] 
dim(aging_control)
#93 171
aging_cases = aging_eprime_final %>%
  filter(!(grp == 0))
dim(aging_cases)
#211 171

###chrt_tc1


m_chrt_tc1 = gls(log(chrt_tc1) ~ age + agesqr,  data = aging_control[which(aging_control$chrt_tc1<1200),], weights = varExp(form=~age), na.action = na.omit)
summary(m_chrt_tc1)

aging_agg = aging_control %>%
  filter(!is.na(chrt_tc1)) %>%
  filter(chrt_tc1 <=1200) %>%
  group_by(age_group) %>%
  summarise(sd_res_chrt_tc1 = sd(resid(m_chrt_tc1)),
            mean_chrt_tc1 = mean(log(chrt_tc1)))


aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")  
dim(aging_control)


aging_control$p_chrt_tc1 = predict(m_chrt_tc1, aging_control)

aging_control = aging_control %>%
  mutate(z_chrt_tc1 = ((log(chrt_tc1) - p_chrt_tc1) / sd_res_chrt_tc1)) 

aging_cases$p_chrt_tc1 = predict(m_chrt_tc1, aging_cases)

aging_cases = aging_cases %>%
  mutate(z_chrt_tc1 = ((log(chrt_tc1) - p_chrt_tc1) / sd_res_chrt_tc1)) 


aging_exp = rbind(aging_cases, aging_control)

##chrt_c1


m_chrt_c1 = gls(chrt_c1 ~ age + agesqr + toteduc + sex, data = aging_control, weights = varExp(form=~age), na.action = na.omit) 
summary(m_chrt_c1)

aging_agg = aging_control %>%
  filter(!is.na(chrt_c1)) %>%
  filter(chrt_c1 >=50) %>%
  group_by(age_group) %>%
  summarise(sd_res_chrt_c1 = sd(resid(m_chrt_c1)),
            mean_chrt_c1 = mean(chrt_c1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 


aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 


aging_control$p_chrt_c1 = predict(m_chrt_c1, aging_control)

aging_control = aging_control %>%
  mutate(z_chrt_c1 = ((chrt_c1 - p_chrt_c1)/ sd_res_chrt_c1)) 

aging_cases$p_chrt_c1 = predict(m_chrt_c1, aging_cases)

aging_cases = aging_cases %>%
  mutate(z_chrt_c1 = ((chrt_c1 - p_chrt_c1)/ sd_res_chrt_c1)) 


aging_exp = rbind(aging_cases, aging_control)

##rtall_wor1


#m = gls(log(rtall_wor1) ~  toteduc + age + agesqr + sex ,  data = aging_control, weights = varExp(form=~age), na.action = na.omit)
#summary(m)

aging_agg = aging_control %>%
  filter(!is.na(rtall_wor1)) %>%
  group_by(age_group) %>%
  summarise(sd_rtall_wor1 = sd(log(rtall_wor1)),
            mean_rtall_wor1 = mean(log(rtall_wor1)))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)


aging_control = aging_control %>%
  mutate(z_rtall_wor1 = ((log(rtall_wor1) - mean_rtall_wor1 )/ sd_rtall_wor1)) 



aging_cases = aging_cases %>%
  mutate(z_rtall_wor1 = ((log(rtall_wor1) - mean_rtall_wor1 )/ sd_rtall_wor1)) 


aging_exp = rbind(aging_cases, aging_control)


##rtall_col1


aging_agg = aging_control %>%
  filter(!is.na(rtall_col1)) %>%
  group_by(age_group) %>%
  summarise(sd_rtall_col1 = sd(log(rtall_col1)),
            mean_rtall_col1 = mean(log(rtall_col1)))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)


aging_control = aging_control %>%
  mutate(z_rtall_col1 = ((log(rtall_col1) - mean_rtall_col1 )/ sd_rtall_col1)) 



aging_cases = aging_cases %>%
  mutate(z_rtall_col1 = ((log(rtall_col1) - mean_rtall_col1 )/ sd_rtall_col1)) 


aging_exp = rbind(aging_cases, aging_control)

##rtall_cwd1


aging_agg = aging_control %>%
  filter(!is.na(rtall_cwd1)) %>%
  group_by(age_group) %>%
  summarise(sd_rtall_cwd1 = sd(log(rtall_cwd1)),
            mean_rtall_cwd1 = mean(log(rtall_cwd1)))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)


aging_control = aging_control %>%
  mutate(z_rtall_cwd1 = ((log(rtall_cwd1) - mean_rtall_cwd1 )/ sd_rtall_cwd1)) 



aging_cases = aging_cases %>%
  mutate(z_rtall_cwd1 = ((log(rtall_cwd1) - mean_rtall_cwd1 )/ sd_rtall_cwd1)) 


aging_exp = rbind(aging_cases, aging_control)

##cstravg1

m_cstravg1 = gls(cstravg1 ~ sex + age, data = aging_control[which(aging_control$cstravg1 <= 1.5),], weights = varExp(form=~age), na.action = na.omit) 
summary(m_cstravg1)

aging_agg = aging_control %>%
  filter(!is.na(cstravg1)) %>%
  filter(cstravg1 <=1.5) %>%
  group_by(age_group) %>%
  summarise(sd_res_cstravg1 = sd(resid(m_cstravg1)),
            mean_cstravg1 = mean(cstravg1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 


aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 



aging_control$p_cstravg1 = predict(m_cstravg1, aging_control)

aging_control = aging_control %>%
  mutate(z_cstravg1 = ((cstravg1 - p_cstravg1)/ sd_res_cstravg1)) 



aging_cases$p_cstravg1 = predict(m_cstravg1, aging_cases)

aging_cases = aging_cases %>%
  mutate(z_cstravg1 = ((cstravg1 - p_cstravg1)/ sd_res_cstravg1)) 


aging_exp = rbind(aging_cases, aging_control)

##5 missing sex variables. making missing values its own category

dim(aging_exp)

stnd_resid_cstravg1 = residuals(m_cstravg1, type = "pearson") 



##ct_wor1



m_ct_wor1 = gls(ct_wor1 ~ age + agesqr, data = aging_control, weights = varExp(form=~age), na.action = na.omit) 
summary(m_ct_wor1)

aging_agg = aging_control %>%
  filter(!is.na(ct_wor1)) %>%
  group_by(age_group) %>%
  summarise(sd_res_ct_wor1 = sd(resid(m_ct_wor1)),
            mean_ct_wor1 = mean(ct_wor1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 
dim(aging_cases)                        


aging_control$p_ct_wor1 = predict(m_ct_wor1, aging_control)

aging_control = aging_control %>%
  mutate(z_ct_wor1 = (ct_wor1 - p_ct_wor1 / sd_res_ct_wor1))


aging_cases$p_ct_wor1 = predict(m_ct_wor1, aging_cases)

aging_cases = aging_cases %>%
  mutate(z_ct_wor1 = (ct_wor1 - p_ct_wor1 / sd_res_ct_wor1)) 



aging_exp = rbind(aging_cases, aging_control)


dim(aging_exp)
stnd_resid_ict_wor1 = residuals(m_ict_wor1, type = "pearson") 

length(which(stnd_resid_ict_wor1 >= 4))


## 1 outlier


###ict_col1


#aging_control$ict_col1 = (41- aging_control$ct_col1)
#aging_cases$ict_col1 = (41- aging_cases$ct_col1)

m_ct_col1 = gls(ct_col1 ~ age , data = aging_control, weights = varExp(form=~age), na.action = na.omit) 
summary(m_ct_col1)

aging_agg = aging_control %>%
  filter(!is.na(ct_col1)) %>%
  group_by(age_group) %>%
  summarise(sd_res_ct_col1 = sd(resid(m_ct_col1)),
            mean_ct_col1 = mean(ct_col1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 


aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 




aging_control$p_ct_col1 = predict(m_ct_col1, aging_control)

aging_control = aging_control %>%
  mutate(z_ct_col1 = ((ct_col1 - p_ct_col1 )/ sd_res_ct_col1)) 



aging_cases$p_ct_col1 = predict(m_ct_col1, aging_cases)

aging_cases = aging_cases %>%
  mutate(z_ct_col1 = ((ct_col1 - p_ct_col1 )/ sd_res_ct_col1)) 


aging_exp = rbind(aging_cases, aging_control)



##ct_cwd1

#aging_control$ict_cwd1 = (89- aging_control$ct_cwd1)
#aging_cases$ict_cwd1 = (89- aging_cases$ct_cwd1)

m_ct_cwd1 = gls(ct_cwd1 ~ age + toteduc, data = aging_control, na.action = na.omit, weights = varExp(form=~age)) 
summary(m_ct_cwd1)


aging_agg = aging_control %>%
  filter(!is.na(ct_cwd1)) %>%
  group_by(age_group) %>%
  summarise(sd_res_ct_cwd1 = sd(resid(m_ct_cwd1)),
            mean_ct_cwd1 = mean(ct_cwd1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 


aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 



aging_control$p_ct_cwd1 = predict(m_ct_cwd1, aging_control)

aging_control = aging_control %>%
  mutate(z_ct_cwd1 = ((ct_cwd1 - p_ct_cwd1 )/ sd_res_ct_cwd1)) 




aging_cases$p_ct_cwd1 = predict(m_ct_cwd1, aging_cases)

aging_cases = aging_cases %>%
  mutate(z_ct_cwd1 = ((ct_cwd1 - p_ct_cwd1 )/ sd_res_ct_cwd1)) 


aging_exp = rbind(aging_cases, aging_control)


##rtall_congruent1




aging_agg = aging_control %>%
  filter(!is.na(rtall_congruent1)) %>%
  group_by(age_group) %>%
  summarise(sd_rtall_congruent1 = sd(log(rtall_congruent1)),
            mean_rtall_congruent1 = mean(log(rtall_congruent1)))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)


aging_control = aging_control %>%
  mutate(z_rtall_congruent1 = ((log(rtall_congruent1) - mean_rtall_congruent1 )/ sd_rtall_congruent1)) 



aging_cases = aging_cases %>%
  mutate(z_rtall_congruent1 = ((log(rtall_congruent1) - mean_rtall_congruent1 )/ sd_rtall_congruent1)) 


aging_exp = rbind(aging_cases, aging_control)


##rtall_incongruent1



aging_agg = aging_control %>%
  filter(!is.na(rtall_incongruent1)) %>%
  group_by(age_group) %>%
  summarise(sd_rtall_incongruent1 = sd(log(rtall_incongruent1)),
            mean_rtall_incongruent1 = mean(log(rtall_incongruent1)))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)


aging_control = aging_control %>%
  mutate(z_rtall_incongruent1 = ((log(rtall_incongruent1) - mean_rtall_incongruent1 )/ sd_rtall_incongruent1)) 



aging_cases = aging_cases %>%
  mutate(z_rtall_incongruent1 = ((log(rtall_incongruent1) - mean_rtall_incongruent1 )/ sd_rtall_incongruent1)) 


aging_exp = rbind(aging_cases, aging_control)

##ccistravg1


aging_agg = aging_control %>%
  filter(!is.na(ccistravg1)) %>%
  group_by(age_group) %>%
  summarise(sd_res_ccistravg1 = sd(ccistravg1),
            mean_ccistravg1 = mean(ccistravg1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 


aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 




aging_control = aging_control %>%
  mutate(z_ccistravg1 = ((ccistravg1 - mean_ccistravg1) / sd_res_ccistravg1)) 


aging_cases = aging_cases %>%
  mutate(z_ccistravg1 = ((ccistravg1 - mean_ccistravg1) / sd_res_ccistravg1)) 



aging_exp = rbind(aging_cases, aging_control)


dim(aging_exp)

stnd_resid_ccistravg1 = residuals(m_ccistravg1, type = "pearson") 


##ctall_congruent1



m_ct_congruent1 = gls(ct_congruent1 ~age + agesqr,  data = aging_control, na.action = na.omit, weights = varExp(form=~age)) 
summary(m_ct_congruent1)


aging_agg = aging_control %>%
  filter(!is.na(ct_congruent1)) %>%
  group_by(age_group) %>%
  summarise(sd_res_ct_congruent1 = sd(resid(m_ct_congruent1)),
            mean_ct_congruent1 = mean(ct_congruent1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 

dim(aging_control)
aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 


dim(aging_cases)

aging_control$p_ct_congruent1 = predict(m_ct_congruent1, aging_control)

aging_control = aging_control %>%
  mutate(z_ct_congruent1 = ((ct_congruent1 - p_ct_congruent1 )/ sd_res_ct_congruent1)) 


aging_cases$p_ct_congruent1 = predict(m_ct_congruent1, aging_cases)

aging_cases = aging_cases %>%
  mutate(z_ct_congruent1 = ((ct_congruent1 - p_ct_congruent1 )/ sd_res_ct_congruent1)) 



aging_exp = rbind(aging_cases, aging_control)



###ct_incongruent

#aging_control$ict_incongruent1 = (52- aging_control$ct_incongruent1)
#aging_cases$ict_incongruent1 = (52- aging_cases$ct_incongruent1)

m_ct_incongruent1 = gls(ct_incongruent1 ~ age  , data = aging_control, na.action = na.omit, weights = varExp(form=~age)) 
summary(m_ct_incongruent1)      


aging_agg = aging_control %>%
  filter(!is.na(ct_incongruent1)) %>%
  group_by(age_group) %>%
  summarise(sd_res_ct_incongruent1 = sd(resid(m_ct_incongruent1)),
            mean_ct_incongruent1 = mean(ct_incongruent1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)


aging_control$p_ct_incongruent1 = predict(m_ct_incongruent1, aging_control)

aging_control = aging_control %>%
  mutate(z_ct_incongruent1 = ((ct_incongruent1 - p_ct_incongruent1) / sd_res_ct_incongruent1)) 


aging_cases$p_ct_incongruent1 = predict(m_ct_incongruent1, aging_cases)

aging_cases = aging_cases %>%
  mutate(z_ct_incongruent1 = ((ct_incongruent1 - p_ct_incongruent1 )/ sd_res_ct_incongruent1)) 



aging_exp = rbind(aging_cases, aging_control)

###rtall_ll1


aging_agg = aging_control %>%
  filter(!is.na(rtall_ll1)) %>%
  group_by(age_group) %>%
  summarise(sd_rtall_ll1 = sd(log(rtall_ll1)),
            mean_rtall_ll1 = mean(log(rtall_ll1)))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)

aging_control = aging_control %>%
  mutate(z_rtall_ll1 = ((log(rtall_ll1) - mean_rtall_ll1 )/ sd_rtall_ll1)) 

aging_cases = aging_cases %>%
  mutate(z_rtall_ll1 = ((log(rtall_ll1) - mean_rtall_ll1 )/ sd_rtall_ll1)) 

aging_exp = rbind(aging_cases, aging_control)


##rtall_lh1



aging_agg = aging_control %>%
  filter(!is.na(rtall_lh1)) %>%
  group_by(age_group) %>%
  summarise(sd_rtall_lh1 = sd(log(rtall_lh1)),
            mean_rtall_lh1 = mean(log(rtall_lh1)))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)


aging_control = aging_control %>%
  mutate(z_rtall_lh1 = ((log(rtall_lh1) - mean_rtall_lh1 )/ sd_rtall_lh1)) 



aging_cases = aging_cases %>%
  mutate(z_rtall_lh1 = ((log(rtall_lh1) - mean_rtall_lh1 )/ sd_rtall_lh1)) 




aging_exp = rbind(aging_cases, aging_control)

##rtall_hl1

aging_agg = aging_control %>%
  filter(!is.na(rtall_hl1)) %>%
  group_by(age_group) %>%
  summarise(sd_rtall_hl1 = sd(log(rtall_hl1)),
            mean_rtall_hl1 = mean(log(rtall_hl1)))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)


aging_control = aging_control %>%
  mutate(z_rtall_hl1 = ((log(rtall_hl1) - mean_rtall_hl1 )/ sd_rtall_hl1)) 



aging_cases = aging_cases %>%
  mutate(z_rtall_hl1 = ((log(rtall_hl1) - mean_rtall_hl1 )/ sd_rtall_hl1)) 


aging_exp = rbind(aging_cases, aging_control)


##rtall_hh1

aging_agg = aging_control %>%
  filter(!is.na(rtall_hh1)) %>%
  group_by(age_group) %>%
  summarise(sd_rtall_hh1 = sd(log(rtall_hh1)),
            mean_rtall_hh1 = mean(log(rtall_hh1)))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)


aging_control = aging_control %>%
  mutate(z_rtall_hh1 = ((log(rtall_hh1) - mean_rtall_hh1 )/ sd_rtall_hh1)) 



aging_cases = aging_cases %>%
  mutate(z_rtall_hh1 = ((log(rtall_hh1) - mean_rtall_hh1 )/ sd_rtall_hh1)) 


aging_exp = rbind(aging_cases, aging_control)


##rtall_neu1


aging_agg = aging_control %>%
  filter(!is.na(rtall_neut1)) %>%
  group_by(age_group) %>%
  summarise(sd_rtall_neut1 = sd(rtall_neut1),
            mean_rtall_neut1 = mean(rtall_neut1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)



aging_control = aging_control %>%
  mutate(z_rtall_neut1 = ((rtall_neut1 - mean_rtall_neut1 )/ sd_rtall_neut1)) 



aging_cases = aging_cases %>%
  mutate(z_rtall_neut1 = ((rtall_neut1 - mean_rtall_neut1 )/ sd_rtall_neut1)) 

#aging_control = aging_control %>%
# mutate(p_rtall_neut1 = "NA")

aging_exp = rbind(aging_cases, aging_control)


##rtall_pos1




aging_agg = aging_control %>%
  filter(!is.na(rtall_pos1)) %>%
  group_by(age_group) %>%
  summarise(sd_rtall_pos1 = sd(log(rtall_pos1)),
            mean_rtall_pos1 = mean(log(rtall_pos1)))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)



aging_control = aging_control %>%
  mutate(z_rtall_pos1 = ((log(rtall_pos1) - mean_rtall_pos1 )/ sd_rtall_pos1)) 

aging_cases = aging_cases %>%
  mutate(z_rtall_pos1 = ((log(rtall_pos1) - mean_rtall_pos1 )/ sd_rtall_pos1)) 


aging_exp = rbind(aging_cases, aging_control)

##rtall_neg1



aging_agg = aging_control %>%
  filter(!is.na(rtall_neg1)) %>%
  group_by(age_group) %>%
  summarise(sd_rtall_neg1 = sd(log(rtall_neg1)),
            mean_rtall_neg1 = mean(log(rtall_neg1)))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)



aging_control = aging_control %>%
  mutate(z_rtall_neg1 = ((log(rtall_neg1) - mean_rtall_neg1 )/ sd_rtall_neg1)) 



aging_cases = aging_cases %>%
  mutate(z_rtall_neg1 = ((log(rtall_neg1) - mean_rtall_neg1 )/ sd_rtall_neg1)) 


aging_exp = rbind(aging_cases, aging_control)


##rtall_sui1





aging_agg = aging_control %>%
  filter(!is.na(rtall_sui1)) %>%
  group_by(age_group) %>%
  summarise(sd_rtall_sui1 = sd(log(rtall_sui1)),
            mean_rtall_sui1 = mean(log(rtall_sui1)))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)



aging_control = aging_control %>%
  mutate(z_rtall_sui1 = ((log(rtall_sui1) - mean_rtall_sui1 )/ sd_rtall_sui1)) 



aging_cases = aging_cases %>%
  mutate(z_rtall_sui1 = ((log(rtall_sui1) - mean_rtall_sui1 )/ sd_rtall_sui1)) 


aging_exp = rbind(aging_cases, aging_control)


##rtall_nlex1



aging_agg = aging_control %>%
  filter(!is.na(rtall_nlex1)) %>%
  group_by(age_group) %>%
  summarise(sd_rtall_nlex1 = sd(log(rtall_nlex1)),
            mean_rtall_nlex1 = mean(log(rtall_nlex1)))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)


aging_control = aging_control %>%
  mutate(z_rtall_nlex1 = ((log(rtall_nlex1) - mean_rtall_nlex1 )/ sd_rtall_nlex1)) 



aging_cases = aging_cases %>%
  mutate(z_rtall_nlex1 = ((log(rtall_nlex1) - mean_rtall_nlex1 )/ sd_rtall_nlex1)) 


aging_exp = rbind(aging_cases, aging_control)


##ct_neut1



aging_agg = aging_control %>%
  filter(!is.na(ct_neut1)) %>%
  group_by(age_group) %>%
  summarise(sd_ct_neut1 = sd(ct_neut1),
            mean_ct_neut1 = mean(ct_neut1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)



aging_control = aging_control %>%
  mutate(z_ct_neut1 = ((ct_neut1 - mean_ct_neut1 )/ sd_ct_neut1)) 



aging_cases = aging_cases %>%
  mutate(z_ct_neut1 = ((ct_neut1 - mean_ct_neut1 )/ sd_ct_neut1)) 



aging_exp = rbind(aging_cases, aging_control)


##ct_pos1


aging_agg = aging_control %>%
  filter(!is.na(ct_pos1)) %>%
  group_by(age_group) %>%
  summarise(sd_ct_pos1 = sd(ct_pos1),
            mean_ct_pos1 = mean(ct_pos1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)

aging_control = aging_control %>%
  mutate(z_ct_pos1 = ((ct_pos1 - mean_ct_pos1 )/ sd_ct_pos1)) 

aging_cases = aging_cases %>%
  mutate(z_ct_pos1 = ((ct_pos1 - mean_ct_pos1 )/ sd_ct_pos1)) 

aging_exp = rbind(aging_cases, aging_control)


##ct_neg1

m =  gls(log(ct_neg1) ~ age, data = aging_control, weights = varExp(form=~age), na.action = na.omit) 
summary(m)

aging_agg = aging_control %>%
  filter(!is.na(ct_neg1)) %>%
  group_by(age_group) %>%
  summarise(sd_ct_neg1 = sd(ct_neg1),
            mean_ct_neg1 = mean(ct_neg1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)


aging_control = aging_control %>%
  mutate(z_ct_neg1 = ((ct_neg1 - mean_ct_neg1 )/ sd_ct_neg1)) 



aging_cases = aging_cases %>%
  mutate(z_ct_neg1 = ((ct_neg1 - mean_ct_neg1 )/ sd_ct_neg1)) 



aging_exp = rbind(aging_cases, aging_control)



##ct_sui1


m_ct_sui1 = gls(ct_sui1 ~ age + agesqr, data = aging_control, weights = varExp(form=~age), na.action = na.omit) 
summary(m_ct_sui1)

aging_agg = aging_control %>%
  filter(!is.na(ct_sui1)) %>%
  filter(ct_sui1 >=15) %>%
  group_by(age_group) %>%
  summarise(sd_res_ct_sui1 = sd(resid(m_ct_sui1)),
            mean_ct_sui1 = mean(ct_sui1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")    

aging_control$p_ct_sui1 = predict(m_ct_sui1, aging_control)

aging_control = aging_control %>%
  mutate(z_ct_sui1 = ((ct_sui1 - p_ct_sui1) / sd_res_ct_sui1)) 


aging_cases$p_ct_sui1 = predict(m_ct_sui1, aging_cases)

aging_cases = aging_cases %>%
  mutate(z_ct_sui1 = ((ct_sui1 - p_ct_sui1 )/ sd_res_ct_sui1)) 



aging_exp = rbind(aging_cases, aging_control)


##ct_nlex1

m_ct_nlex1 = gls(ct_nlex1~   age ,  data = aging_control, weights = varExp(form=~age), na.action = na.omit)
summary(m_ct_nlex1)


aging_agg = aging_control %>%
  filter(!is.na(ct_nlex1)) %>%
  group_by(age_group) %>%
  summarise(sd_res_ct_nlex1 = sd(resid(m_ct_nlex1)),
            mean_ct_nlex1 = mean(ct_nlex1))


aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")  
dim(aging_cases)


aging_control$p_ct_nlex1 = predict(m_ct_nlex1, aging_control)

aging_control = aging_control %>%
  mutate(z_ct_nlex1 = ((ct_nlex1 - p_ct_nlex1) / sd_res_ct_nlex1)) 


aging_cases$p_ct_nlex1 = predict(m_ct_nlex1, aging_cases)

aging_cases = aging_cases %>%
  mutate(z_ct_nlex1 = ((ct_nlex1 - p_ct_nlex1) / sd_res_ct_nlex1)) 



aging_exp = rbind(aging_cases, aging_control)


##anb_rt1


m_anb_rt1 = gls(log(anb_rt1) ~ age ,  data = aging_control[which(aging_control$anb_rt1 >7),], weights = varExp(form=~age), na.action = na.omit)
summary(m_anb_rt1)


aging_agg = aging_control %>%
  filter(!is.na(anb_rt1)) %>%
  filter(anb_rt1 > 7) %>%
  group_by(age_group) %>%
  summarise(sd_res_anb_rt1 = sd(resid(m_anb_rt1)),
            mean_anb_rt1 = mean(log(anb_rt1)))


aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")  
dim(aging_cases)


aging_control$p_anb_rt1 = predict(m_anb_rt1, aging_control)

aging_control = aging_control %>%
  mutate(z_anb_rt1 = ((log(anb_rt1) - p_anb_rt1) / sd_res_anb_rt1)) 


aging_cases$p_anb_rt1 = predict(m_anb_rt1, aging_cases)

aging_cases = aging_cases %>%
  mutate(z_anb_rt1 = ((log(anb_rt1) - p_anb_rt1) / sd_res_anb_rt1)) 



aging_exp = rbind(aging_cases, aging_control)


##anb_ct1

m_anb_ct1 = gls(anb_ct1 ~ toteduc ,  data = aging_control, weights = varExp(form=~age), na.action = na.omit)
summary(m_anb_ct1)

aging_agg = aging_control %>%
  filter(!is.na(anb_ct1)) %>%
  group_by(age_group) %>%
  summarise(sd_res_anb_ct1 = sd(resid(m_anb_ct1)),
            mean_anb_ct1 = mean(anb_ct1))


aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")  
dim(aging_cases)


aging_control$p_anb_ct1 = predict(m_anb_ct1, aging_control)

aging_control = aging_control %>%
  mutate(z_anb_ct1 = ((anb_ct1 - p_anb_ct1) / sd_res_anb_ct1)) 


aging_cases$p_anb_ct1 = predict(m_anb_ct1, aging_cases)

aging_cases = aging_cases %>%
  mutate(z_anb_ct1 = ((anb_ct1 - p_anb_ct1) / sd_res_anb_ct1)) 



aging_exp = rbind(aging_cases, aging_control)


##anb_srt1

aging_agg = aging_control %>%
  filter(!is.na(anb_srt1)) %>%
  group_by(age_group) %>%
  summarise(sd_anb_srt1 = sd(anb_srt1),
            mean_anb_srt1 = mean(anb_srt1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)

aging_control = aging_control %>%
  mutate(z_anb_srt1 = ((anb_srt1 - mean_anb_srt1 )/ sd_anb_srt1)) 



aging_cases = aging_cases %>%
  mutate(z_anb_srt1 = ((anb_srt1 - mean_anb_srt1 )/ sd_anb_srt1)) 



aging_exp = rbind(aging_cases, aging_control)


##anb_sct1



aging_agg = aging_control %>%
  filter(!is.na(anb_sct1)) %>%
  group_by(age_group) %>%
  summarise(sd_anb_sct1 = sd(anb_sct1),
            mean_anb_sct1 = mean(anb_sct1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)

aging_control = aging_control %>%
  mutate(z_anb_sct1 = ((anb_sct1 - mean_anb_sct1 )/ sd_anb_sct1)) 



aging_cases = aging_cases %>%
  mutate(z_anb_sct1 = ((anb_sct1 - mean_anb_sct1 )/ sd_anb_sct1)) 

#aging_control = aging_control %>%
# mutate(p_anb_sct1 = "NA")

aging_exp = rbind(aging_cases, aging_control)



##anb_nrt1

m_anb_nrt1 = gls(log(anb_nrt1) ~ age ,  data = aging_control, weights = varExp(form=~age), na.action = na.omit)
summary(m_anb_nrt1)


aging_agg = aging_control %>%
  filter(!is.na(anb_nrt1)) %>%
  group_by(age_group) %>%
  summarise(sd_res_anb_nrt1 = sd(resid(m_anb_nrt1)),
            mean_anb_nrt1 = mean(log(anb_nrt1)))


aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")  
dim(aging_cases)


aging_control$p_anb_nrt1 = predict(m_anb_nrt1, aging_control)

aging_control = aging_control %>%
  mutate(z_anb_nrt1 = ((log(anb_nrt1) - p_anb_nrt1) / sd_res_anb_nrt1)) 


aging_cases$p_anb_nrt1 = predict(m_anb_nrt1, aging_cases)

aging_cases = aging_cases %>%
  mutate(z_anb_nrt1 = ((log(anb_nrt1) - p_anb_nrt1) / sd_res_anb_nrt1)) 



aging_exp = rbind(aging_cases, aging_control)




##anb_nct1

m_anb_nct1 = gls(anb_nct1 ~  toteduc ,  data = aging_control, weights = varExp(form=~age), na.action = na.omit)
summary(m_anb_nct1)


aging_agg = aging_control %>%
  filter(!is.na(anb_nct1)) %>%
  group_by(age_group) %>%
  summarise(sd_res_anb_nct1 = sd(resid(m_anb_nct1)),
            mean_anb_nct1 = mean(anb_nct1))


aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")  
dim(aging_cases)


aging_control$p_anb_nct1 = predict(m_anb_nct1, aging_control)

aging_control = aging_control %>%
  mutate(z_anb_nct1 = ((anb_nct1 - p_anb_nct1) / sd_res_anb_nct1)) 


aging_cases$p_anb_nct1 = predict(m_anb_nct1, aging_cases)

aging_cases = aging_cases %>%
  mutate(z_anb_nct1 = ((anb_nct1 - p_anb_nct1) / sd_res_anb_nct1)) 



aging_exp = rbind(aging_cases, aging_control)


##rl_true_errors


aging_agg = aging_control %>%
  filter(!is.na(rl_true_errors1)) %>%
  group_by(age_group) %>%
  summarise(sd_rl_true_errors1 = sd(rl_true_errors1),
            mean_rl_true_errors1 = mean(rl_true_errors1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)

aging_control = aging_control %>%
  mutate(z_rl_true_errors1 = ((rl_true_errors1 - mean_rl_true_errors1 )/ sd_rl_true_errors1)) 

aging_cases = aging_cases %>%
  mutate(z_rl_true_errors1 = ((rl_true_errors1 - mean_rl_true_errors1 )/ sd_rl_true_errors1)) 

aging_exp = rbind(aging_cases, aging_control)

##rl_true_pphase_errors



aging_agg = aging_control %>%
  filter(!is.na(rl_true_pphase_errors1)) %>%
  group_by(age_group) %>%
  summarise(sd_rl_true_pphase_errors1 = sd(rl_true_pphase_errors1),
            mean_rl_true_pphase_errors1 = mean(rl_true_pphase_errors1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)

aging_control = aging_control %>%
  mutate(z_rl_true_pphase_errors1 = ((rl_true_pphase_errors1 - mean_rl_true_pphase_errors1 )/ sd_rl_true_pphase_errors1)) 

aging_cases = aging_cases %>%
  mutate(z_rl_true_pphase_errors1 = ((rl_true_pphase_errors1 - mean_rl_true_pphase_errors1 )/ sd_rl_true_pphase_errors1)) 

aging_exp = rbind(aging_cases, aging_control)


##rl_true_persev_errors1



aging_agg = aging_control %>%
  filter(!is.na(rl_true_persev_errors1)) %>%
  group_by(age_group) %>%
  summarise(sd_rl_true_persev_errors1 = sd(rl_true_persev_errors1),
            mean_rl_true_persev_errors1 = mean(rl_true_persev_errors1))



aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")                 

dim(aging_cases)

aging_control = aging_control %>%
  mutate(z_rl_true_persev_errors1 = ((rl_true_persev_errors1 - mean_rl_true_persev_errors1 )/ sd_rl_true_persev_errors1)) 

aging_cases = aging_cases %>%
  mutate(z_rl_true_persev_errors1 = ((rl_true_persev_errors1 - mean_rl_true_persev_errors1 )/ sd_rl_true_persev_errors1)) 

aging_exp = rbind(aging_cases, aging_control)



##rl_rt_correct1

m_rl_rt_correct1 = gls(rl_rt_correct1 ~  age, data = aging_control, na.action = na.omit, weights = varExp(form=~age))
summary(m_rl_rt_correct1)

aging_agg = aging_control %>%
  filter(!is.na(rl_rt_correct1)) %>%
  group_by(age_group) %>%
  summarise(sd_res_rl_rt_correct1 = sd(resid(m_rl_rt_correct1)),
            mean_rl_rt_correct1 = mean(log(rl_rt_correct1)))


aging_control =   sqldf("SELECT *
                        FROM aging_control as ac
                        LEFT JOIN aging_agg as agg
                        ON ac.age_group = agg.age_group")                 
dim(aging_control)

aging_cases = sqldf("SELECT *
                    FROM aging_cases as aca
                    LEFT JOIN aging_agg as agg
                    ON aca.age_group = agg.age_group")  
dim(aging_cases)


aging_control$p_rl_rt_correct1 = predict(m_rl_rt_correct1, aging_control)

aging_control = aging_control %>%
  mutate(z_rl_rt_correct1 = ((log(rl_rt_correct1) - p_rl_rt_correct1) / sd_res_rl_rt_correct1)) 


aging_cases$p_rl_rt_correct1 = predict(m_rl_rt_correct1, aging_cases)

aging_cases = aging_cases %>%
  mutate(z_rl_rt_correct1 = ((log(rl_rt_correct1) - p_rl_rt_correct1) / sd_res_rl_rt_correct1)) 



aging_exp = rbind(aging_cases, aging_control)

##finaldat

finaldat = aging_exp %>%
  dplyr::select(subject, age_group, age, sex, toteduc,  starts_with("z"))


head(finaldat)
dim(finaldat)

write.csv(finaldat, "new_z_scores_final.csv")


attach(aging_eprime_final)

plot(age, rt_ll1, main="Scatterplot Example",
     xlab="Age ", ylab="rt_LL1", pch=19)



attach(finaldat)
plot(age, z_rtall_ll1, main="Scatterplot Example",
     xlab="Age ", ylab="z_rtall_ll1", pch=19)



