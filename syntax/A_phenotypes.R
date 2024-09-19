################################################################################
# To use:
# 1. Put this script into the folder with the WLS .dta file
# 2. Change the file name on line 19
# 3. Change familypub on line 20 to family identifier in private data
# 4. If needed, uncomment and use line 16 to set working directory
# 5. Edit variable names on line 26 to line up with restricted data
# 6. All of the ADHD code only runs with restricted data, not public data
# 7. Check sample output .csv and let me know if it looks good
################################################################################

library(haven)
library(dplyr)
library(tidyr)

pwd()
setwd('..')

# idpub + rytpe or familypub + personid uniquely identify
d <- read_dta('./data/wls/wls_bl_14.03.stata/wls_bl_14_03.dta') |>
  mutate(id = paste0(as.character(familypub), '-', as.character(personid))) |>
  mutate(female = case_when(z_sexrsp == 1 ~ 0,
                            z_sexrsp == 2 ~ 1,
                            TRUE ~ NA_real_),
         birthdate = if_else(z_brdxdy >= 0, as.numeric(z_brdxdy), NA_real_))

id <- select(d, id, female, idpub, rtype, familypub, personid)
resid_controls <- select(d, id, female, birthdate)

# Earlier versions residualized phenotypes relative to age and sex, now they are
# centered within wave. to change this behavior, set this flag to TRUE
resid_phen <- FALSE


#################
# anthopometric #
#################

# bmi

bmi <- d |>
  mutate(bmi_92 = if_else(z_mx011rec >= 0, as.numeric(z_mx011rec), NA_real_),
         bmi_03 = if_else(z_ix011rec >= 0, as.numeric(z_ix011rec), NA_real_),
         bmi_11 = if_else(z_jx011rec >= 0, as.numeric(z_jx011rec), NA_real_)) |>
  select(id, birthdate, female, bmi_92, bmi_03, bmi_11)

bmi_92 <- bmi |>
  select(id, birthdate, female, bmi_92) |>
  na.omit()

bmi_03 <- bmi |>
  select(id, birthdate, female, bmi_03) |>
  na.omit()

bmi_11 <- bmi |>
  select(id, birthdate, female, bmi_11) |>
  na.omit()

if(resid_phen){
  m_bmi_92 <- lm(bmi_92 ~ I(birthdate^2)*female + birthdate*female, bmi_92)
  m_bmi_03 <- lm(bmi_03 ~ I(birthdate^2)*female + birthdate*female, bmi_03)
  m_bmi_11 <- lm(bmi_11 ~ I(birthdate^2)*female + birthdate*female, bmi_11)
} else {
  m_bmi_92 <- lm(bmi_92 ~ 1, bmi_92)
  m_bmi_03 <- lm(bmi_03 ~ 1, bmi_03)
  m_bmi_11 <- lm(bmi_11 ~ 1, bmi_11)
}

bmi_92$bmi_92_resid <- resid(m_bmi_92)
bmi_03$bmi_03_resid <- resid(m_bmi_03)
bmi_11$bmi_11_resid <- resid(m_bmi_11)

bmi <- bmi_92 |>
  full_join(bmi_03, by=c('id', 'female', 'birthdate')) |>
  full_join(bmi_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(bmi = mean(value, na.rm=TRUE))

# height in meters

height <- d |>
  mutate(height = 0.0254*pmax(as.numeric(z_ix010rec), as.numeric(z_jx010rec),
                              na.rm=TRUE)) |>
  mutate(height = if_else(height <= 0, NA_real_, height)) |>
  select(id, height) |>
  left_join(resid_controls, by='id') |>
  na.omit()

if(resid_phen){
  m_height <- lm(height ~ I(birthdate^2)*female + birthdate*female, height)
} else {
  m_height <- lm(height ~ 1, height)
}

height$height <- resid(m_height)
height <- select(height, id, height)

###########################
# cognition and education #
###########################

# cognitive performance (1957-1964 graduates)

cognitive_performance <- d |>
  mutate(iq = as.numeric(z_gwiiq_bm)) |>
  mutate(cognitive_performance = if_else(iq >=0, iq, NA_real_)) |>
  select(id, cognitive_performance) |>
  left_join(resid_controls, by='id') |>
  na.omit()

if(resid_phen){
  m_cognitive_performance <- lm(cognitive_performance ~ I(birthdate^2)*female + birthdate*female, cognitive_performance)
} else {
  m_cognitive_performance <- lm(cognitive_performance ~ 1, cognitive_performance)
}

cognitive_performance$cognitive_performance <- resid(m_cognitive_performance)
cognitive_performance <- select(cognitive_performance, id, cognitive_performance)

# educational attainment

ed_attainment <- d |>
  mutate(ed_75 = if_else(z_edeqyr >= 0, as.numeric(z_edeqyr), NA_real_),
         ed_92 = if_else(z_rb003red >= 0, as.numeric(z_rb003red), NA_real_),
         ed_03 = if_else(z_gb103red >= 0, as.numeric(z_gb103red), NA_real_),
         ed_11 = if_else(z_hb103red >= 0, as.numeric(z_hb103red), NA_real_)) |>
  select(id, starts_with('ed_')) |>
  pivot_longer(-id) |>
  na.omit() |>
  group_by(id) |>
  summarize(ed_attainment = max(value)) |>
  left_join(resid_controls, by='id') |>
  na.omit()

if(resid_phen){
  m_ed_attainment <- lm(ed_attainment ~ I(birthdate^2)*female + birthdate*female, ed_attainment)
} else {
  m_ed_attainment <- lm(ed_attainment ~ 1, ed_attainment)
}

ed_attainment$ed_attainment <- resid(m_ed_attainment)
ed_attainment <- select(ed_attainment, id, ed_attainment)


####################################
# fertility and sexual development #
####################################

# age first birth - in months

age_first_birth <- d |>
  mutate(afb_75 = if_else(z_agrkd1 >= 0, as.numeric(z_agrkd1), NA_real_),
         dfb_11 = if_else(z_hd01401 >= 0, as.numeric(z_hd01401), NA_real_),
         dob = 12*birthdate) |>
  mutate(afb_11 = dfb_11 - dob) |>
  mutate(age_first_birth = case_when(!is.na(afb_11) ~ afb_11,
                                     !is.na(afb_75) ~ afb_75,
                                     TRUE ~ NA_real_)) |>
  select(id, age_first_birth) |>
  left_join(resid_controls, by='id') |>
  na.omit()

if(resid_phen){
  m_age_first_birth <- lm(age_first_birth ~ I(birthdate^2)*female + birthdate*female, age_first_birth)
} else {
  m_age_first_birth <- lm(age_first_birth ~ 1, age_first_birth)
}

age_first_birth$age_first_birth <- resid(m_age_first_birth)
age_first_birth <- select(age_first_birth, id, age_first_birth)

# age first menses

age_first_menses <- d |>
  filter(female ==1) |>
  mutate(age_first_menses =
           if_else(z_in190rer >= 0, as.numeric(z_in190rer), NA_real_)) |>
  select(id, age_first_menses) |>
  left_join(resid_controls, by='id') |>
  na.omit()

if(resid_phen){
  m_age_first_menses <- lm(age_first_menses ~ I(birthdate^2)*female + birthdate*female, age_first_menses)
} else {
  m_age_first_menses <- lm(age_first_menses ~ 1, age_first_menses)
}

age_first_menses$age_first_menses <- resid(m_age_first_menses)
age_first_menses <- select(age_first_menses, id, age_first_menses)


# number ever born

number_ever_born <- d |>
  mutate(number_ever_born =
           if_else(z_hd301kd >= 0, as.numeric(z_hd301kd), NA_real_)) |>
  select(id, number_ever_born) |>
  left_join(resid_controls, by='id') |>
  na.omit()

if(resid_phen){
  m_number_ever_born <- lm(number_ever_born ~ I(birthdate^2)*female + birthdate*female, number_ever_born)
} else {
  m_number_ever_born <- lm(number_ever_born ~ 1, number_ever_born)
}

number_ever_born$number_ever_born <- resid(m_number_ever_born)
number_ever_born <- select(number_ever_born, id, number_ever_born)

###############################
# health and health behaviors #
###############################

# alcohol misuse

alcohol_misuse <- d |>
  mutate(am_1_92 = case_when(z_ru025re == 1 ~ 1,
                             z_ru025re == 2 ~ 0,
                             TRUE ~ NA_real_),
         am_1_03 = case_when(z_gu025re == 1 ~ 1,
                             z_gu025re == 2 ~ 0,
                             TRUE ~ NA_real_),
         am_1_11 = case_when(z_hu025re == 1 ~ 1,
                             z_hu025re == 2 ~ 0,
                             TRUE ~ NA_real_),

         am_2_92 = case_when(z_ru026re > 4 ~ 1,
                             z_ru026re >= 0 ~ 0,
                             TRUE ~ NA_real_),
         am_2_03 = case_when(z_gu026re > 4 ~ 1,
                             z_gu026re >= 0 ~ 0,
                             TRUE ~ NA_real_),
         am_2_11 = case_when(z_hu026re > 4 ~ 1,
                             z_hu026re >= 0 ~ 0,
                             TRUE ~ NA_real_),

         am_3_92 = case_when(z_ru027re > 4 ~ 1,
                             z_ru027re >= 0 ~ 0,
                             TRUE ~ NA_real_),
         am_3_03 = case_when(z_gu027re > 4 ~ 1,
                             z_gu027re >= 0 ~ 0,
                             TRUE ~ NA_real_),
         am_3_11 = case_when(z_hu027re > 4 ~ 1,
                             z_hu027re >= 0 ~ 0,
                             TRUE ~ NA_real_),

         am_4_92 = case_when(z_ru029re > 1 ~ 1,
                             z_ru029re >= 0 ~ 0,
                             TRUE ~ NA_real_),
         am_4_03 = case_when(z_gu029re > 1 ~ 1,
                             z_gu029re >= 0 ~ 0,
                             TRUE ~ NA_real_),
         am_4_11 = case_when(z_hu029re > 1 ~ 1,
                             z_hu029re >= 0 ~ 0,
                             TRUE ~ NA_real_),

         am_5_92 = case_when(z_ru030re == 1 ~ 1,
                             z_ru030re == 2 ~ 0,
                             TRUE ~ NA_real_),
         am_5_03 = case_when(z_gu030re == 1 ~ 1,
                             z_gu030re == 2 ~ 0,
                             TRUE ~ NA_real_),
         am_5_11 = case_when(z_hu030re == 1 ~ 1,
                             z_hu030re == 2 ~ 0,
                             TRUE ~ NA_real_),

         am_6_92 = case_when((z_ru035re == 1) & (z_ru036re == 1) ~ 1,
                             (z_ru035re == 1) & (z_ru036re == 3) ~ 1,
                             z_ru035re > 0 ~ 0,
                             TRUE ~ NA_real_),
         am_6_03 = case_when((z_gu035re == 1) & (z_gu036re == 1) ~ 1,
                             (z_gu035re == 1) & (z_gu036re == 3) ~ 1,
                             z_gu035re > 0 ~ 0,
                             TRUE ~ NA_real_),
         am_6_11 = case_when((z_hu035re == 1) & (z_hu036re == 1) ~ 1,
                             (z_hu035re == 1) & (z_hu036re == 3) ~ 1,
                             z_hu035re > 0 ~ 0,
                             TRUE ~ NA_real_),

          am_8_92 = case_when(z_ru032re == 1 ~ 1,
                              z_ru032re == 2 ~ 0,
                              TRUE ~ NA_real_),
          am_8_03 = case_when(z_gu032re == 1 ~ 1,
                              z_gu032re == 2 ~ 0,
                              TRUE ~ NA_real_),
          am_8_11 = case_when(z_hu032re == 1 ~ 1,
                              z_hu032re == 2 ~ 0,
                              TRUE ~ NA_real_),


          am_9_92 = case_when(z_ru033re == 1 ~ 1,
                              z_ru033re == 2 ~ 0,
                              TRUE ~ NA_real_),
          am_9_03 = case_when(z_gu033re == 1 ~ 1,
                              z_gu033re == 2 ~ 0,
                              TRUE ~ NA_real_),
          am_9_11 = case_when(z_hu033re == 1 ~ 1,
                              z_hu033re == 2 ~ 0,
                              TRUE ~ NA_real_)) |>
  select(id, starts_with('am_'))


am_missing <- alcohol_misuse |>
  pivot_longer(-id, names_to=c('item', 'year'), values_to='resp', names_sep=5) |>
  group_by(id, year) |>
  summarize(am_missing = if_else(sum(is.na(resp)) > 4, 1, 0),
            .groups='drop')

am_no <- alcohol_misuse |>
  pivot_longer(-id, names_to=c('item', 'year'), values_to='resp', names_sep=5) |>
  filter(item == 'am_1_') |>
  mutate(am_no = if_else(resp == 0, 1, 0)) |>
  select(id, year, am_no)



am <- alcohol_misuse |>
  pivot_longer(-id, names_to=c('item', 'year'), values_to='resp', names_sep=5) |>
  filter(item != 'am_1_') |>
  group_by(id, year) |>
  summarize(mean = mean(resp, na.rm=TRUE),
            .groups='drop') |>
  full_join(am_missing, by=c('id', 'year')) |>
  full_join(am_no, by=c('id', 'year')) |>
  mutate(am = case_when(am_missing == 1 ~ NA_real_,
                        am_no == 1 ~ 0,
                        TRUE ~ mean)) |>
  select(id, year, am) |>
  left_join(resid_controls, by='id')


am_92 <- am |>
  filter(year == '92') |>
  na.omit()

am_03 <- am |>
  filter(year == '03') |>
  na.omit()

am_11 <- am |>
  filter(year == '11') |>
  na.omit()

if(resid_phen){
  m_am_92 <- lm(am ~ I(birthdate^2)*female + birthdate*female, am_92)
  m_am_03 <- lm(am ~ I(birthdate^2)*female + birthdate*female, am_03)
  m_am_11 <- lm(am ~ I(birthdate^2)*female + birthdate*female, am_11)
} else {
  m_am_92 <- lm(am ~ 1, am_92)
  m_am_03 <- lm(am ~ 1, am_03)
  m_am_11 <- lm(am ~ 1, am_11)
}

am_92$am_92_resid <- resid(m_am_92)
am_03$am_03_resid <- resid(m_am_03)
am_11$am_11_resid <- resid(m_am_11)

alcohol_misuse <- am_92 |>
  full_join(am_03, by=c('id', 'female', 'birthdate')) |>
  full_join(am_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(alcohol_misuse = mean(value, na.rm=TRUE))

# allergies - cat, dust, pollen, hayfever

allergies <- d |>
  mutate(allergy_cat = case_when(is.na(z_jx109rer) ~ NA_real_,
                                 z_jx109rer < 0 ~ NA_real_,
                                 z_jx410rer == -2 ~ 0,
                                 z_jx410rer < 0 ~ NA_real_,
                                 z_jx410rer == 22 ~ NA_real_,
                                 z_jx410rer == 20 ~ 1,
                                 TRUE ~ 0),
         allergy_dust = case_when(is.na(z_jx109rer) ~ NA_real_,
                                  z_jx109rer < 0 ~ NA_real_,
                                  z_jx410rer == -2 ~ 0,
                                  z_jx410rer < 0 ~ NA_real_,
                                  z_jx410rer == 22 ~ NA_real_,
                                  z_jx410rer == 6 ~ 1,
                                  TRUE ~ 0),
         allergy_pollen = case_when(is.na(z_jx109rer) ~ NA_real_,
                                    z_jx109rer < 0 ~ NA_real_,
                                    z_jx410rer == -2 ~ 0,
                                    z_jx410rer < 0 ~ NA_real_,
                                    z_jx410rer == 22 ~ NA_real_,
                                    z_jx410rer == 3 ~ 1,
                                    TRUE ~ 0),
         hayfever = case_when(is.na(z_jx109rer) ~ NA_real_,
                              z_jx109rer < 0 ~ NA_real_,
                              z_jx410rer == -2 ~ 0,
                              z_jx410rer < 0 ~ NA_real_,
                              z_jx410rer == 22 ~ NA_real_,
                              z_jx410rer == 14 ~ 1,
                              TRUE ~ 0)) |>
  select(id, starts_with('allergy'), hayfever) |>
  left_join(resid_controls, by='id')

## cat

cat <- allergies |>
  select(id, female, birthdate, allergy_cat) |>
  na.omit()

if(resid_phen){
  m_cat <- lm(allergy_cat ~ I(birthdate^2)*female + birthdate*female, cat)
} else {
  m_cat <- lm(allergy_cat ~ 1, cat)
}

cat$allergy_cat <- resid(m_cat)
cat <- select(cat, id, allergy_cat)

## dust

dust <- allergies |>
  select(id, female, birthdate, allergy_dust) |>
  na.omit()

if(resid_phen){
  m_dust <- lm(allergy_dust ~ I(birthdate^2)*female + birthdate*female, dust)
} else {
  m_dust <- lm(allergy_dust ~ 1, dust)
}

dust$allergy_dust <- resid(m_dust)
dust <- select(dust, id, allergy_dust)

## pollen

pollen <- allergies |>
  select(id, female, birthdate, allergy_pollen) |>
  na.omit()

if(resid_phen){
  m_pollen <- lm(allergy_pollen ~ I(birthdate^2)*female + birthdate*female, pollen)
} else {
  m_pollen <- lm(allergy_pollen ~ 1, pollen)
}

pollen$allergy_pollen <- resid(m_pollen)
pollen <- select(pollen, id, allergy_pollen)

## hayfever

hayfever <- allergies |>
  select(id, female, birthdate, hayfever) |>
  na.omit()

if(resid_phen){
  m_hayfever <- lm(hayfever ~ I(birthdate^2)*female + birthdate*female, hayfever)
} else {
  m_hayfever <- lm(hayfever ~ 1, hayfever)
}

hayfever$hayfever <- resid(m_hayfever)
hayfever <- select(hayfever, id, hayfever)

## combine

allergies <- cat |>
  full_join(dust, by='id') |>
  full_join(pollen, by='id') |>
  full_join(hayfever, by='id')

# asthma

asthma <- d |>
  mutate(asthma = case_when(z_mx085rer == 1 ~ 1,
                            z_ix085rer == 1 ~ 1,
                            z_jx085rer == 1 ~ 1,
                            z_mx085rer == 2 ~ 0,
                            z_ix085rer == 2 ~ 0,
                            z_jx085rer == 2 ~ 0,
                            TRUE ~ NA_real_)) |>
  mutate(asthma_eczema_rhinitis = case_when(asthma == 1 ~ 1,
                                            z_jx410rer == 14 ~ 1,
                                            asthma == 0 ~ 0,
                                            z_jx109rer >= 0 ~ 0,
                                            TRUE ~ NA_real_)) |>
  select(id, female, birthdate, asthma, asthma_eczema_rhinitis)

ast <- asthma |>
  select(id, female, birthdate, asthma) |>
  na.omit()

aer <- asthma |>
  select(id, female, birthdate, asthma_eczema_rhinitis) |>
  na.omit()

if(resid_phen){
  m_ast <- lm(asthma ~ I(birthdate^2)*female + birthdate*female, ast)
  m_aer <- lm(asthma_eczema_rhinitis ~ I(birthdate^2)*female + birthdate*female, aer)
} else {
  m_ast <- lm(asthma ~ 1, ast)
  m_aer <- lm(asthma_eczema_rhinitis ~ 1, aer)
}

ast$asthma <- resid(m_ast)
aer$asthma_eczema_rhinitis <- resid(m_aer)

ast <- select(ast, id, asthma)
aer <- select(aer, id, asthma_eczema_rhinitis)

asthma <- full_join(ast, aer, by='id')

# adhd

# adhd <- d |>
#   mutate(adhd = case_when(z_hh125rek == 1 ~ 1,
#                           z_hh125rek == 2 ~ 0,
#                           z_hh125rek == -2 ~ 0,
#                           TRUE ~ NA_real_)) |>
#   select(id, adhd) |>
#   left_join(resid_controls, by='id') |>
#   na.omit()
#
# if(resid_phen){
#   m_adhd <- lm(adhd ~ I(birthdate^2)*female + birthdate*female, adhd)
# } else {
#   m_adhd <- lm(adhd ~ 1, adhd)
# }
#
# adhd$adhd <- resid(m_adhd)
# adhd <- select(adhd, id, adhd)

# cigarettes per day -

cpd <- d |>
  mutate(cpd_92 = case_when(mx015rer < 0 ~ NA_real_,
                            mx015rer > 4 ~ NA_real_,
                            TRUE ~ as.numeric(mx015rer)),
         cpd_03 = case_when(z_ix015rer < 0 ~ NA_real_,
                            z_ix015rer > 10 ~ NA_real_,
                            TRUE ~ as.numeric(z_ix015rer)),
         cpd_11 = case_when(z_jx015rer < 0 ~ NA_real_,
                            z_jx015rer > 10 ~ NA_real_,
                            TRUE ~ as.numeric(z_jx015rer)))


cpd_92 <- cpd |>
  select(id, birthdate, female, cpd_92) |>
  na.omit()

cpd_03 <- cpd |>
  select(id, birthdate, female, cpd_03) |>
  na.omit()

cpd_11 <- cpd |>
  select(id, birthdate, female, cpd_11) |>
  na.omit()

if(resid_phen){
  m_cpd_92 <- lm(cpd_92 ~ I(birthdate^2)*female + birthdate*female, cpd_92)
  m_cpd_03 <- lm(cpd_03 ~ I(birthdate^2)*female + birthdate*female, cpd_03)
  m_cpd_11 <- lm(cpd_11 ~ I(birthdate^2)*female + birthdate*female, cpd_11)
} else {
  m_cpd_92 <- lm(cpd_92 ~ 1, cpd_92)
  m_cpd_03 <- lm(cpd_03 ~ 1, cpd_03)
  m_cpd_11 <- lm(cpd_11 ~ 1, cpd_11)
}

cpd_92$cpd_92_resid <- resid(m_cpd_92)
cpd_03$cpd_03_resid <- resid(m_cpd_03)
cpd_11$cpd_11_resid <- resid(m_cpd_11)

cigarettes_per_day <- cpd_92 |>
  full_join(cpd_03, by=c('id', 'female', 'birthdate')) |>
  full_join(cpd_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(cigarettes_per_day = mean(value, na.rm=TRUE))

# copd

copd <- d |>
  mutate(copd = case_when(z_mx089rer == 1 ~ 1,
                          z_ix089rer == 1 ~ 1,
                          z_jx089rer == 1 ~ 1,
                          z_mx089rer == 2 ~ 0,
                          z_ix089rer == 2 ~ 0,
                          z_jx089rer == 2 ~ 0,
                          TRUE ~ NA_real_)) |>
  select(id, copd) |>
  left_join(resid_controls, by='id') |>
  na.omit()

if(resid_phen){
  m_copd <- lm(copd ~ I(birthdate^2)*female + birthdate*female, copd)
} else {
  m_copd <- lm(copd ~ 1, copd)
}

copd$copd <- resid(m_copd)
copd <- select(copd, id, copd)

# depressive symptoms

ds <- d |>
  mutate(ds_92 = if_else(z_mu001rec >= 0, as.numeric(z_mu001rec), NA_real_),
         ds_03 = if_else(z_iu001rec >= 0, as.numeric(z_mu001rec), NA_real_),
         ds_11 = if_else(z_ju001rec >= 0, as.numeric(z_mu001rec), NA_real_))

ds_92 <- ds |>
  select(id, birthdate, female, ds_92) |>
  na.omit()

ds_03 <- ds |>
  select(id, birthdate, female, ds_03) |>
  na.omit()

ds_11 <- ds |>
  select(id, birthdate, female, ds_11) |>
  na.omit()

if(resid_phen){
  m_ds_92 <- lm(ds_92 ~ I(birthdate^2)*female + birthdate*female, ds_92)
  m_ds_03 <- lm(ds_03 ~ I(birthdate^2)*female + birthdate*female, ds_03)
  m_ds_11 <- lm(ds_11 ~ I(birthdate^2)*female + birthdate*female, ds_11)
} else {
  m_ds_92 <- lm(ds_92 ~ 1, ds_92)
  m_ds_03 <- lm(ds_03 ~ 1, ds_03)
  m_ds_11 <- lm(ds_11 ~ 1, ds_11)
}

ds_92$ds_92_resid <- resid(m_ds_92)
ds_03$ds_03_resid <- resid(m_ds_03)
ds_11$ds_11_resid <- resid(m_ds_11)

depressive_symptoms <- ds_92 |>
  full_join(ds_03, by=c('id', 'female', 'birthdate')) |>
  full_join(ds_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(depressive_symptoms = mean(value, na.rm=TRUE))


# depressive symptoms (alternate)

reverse <- c("z_mu006rer", "z_mu009rer", "z_mu016rer", "z_mu018rer",
             "z_iu006rer", "z_iu009rer", "z_iu016rer", "z_iu018rer",
             "z_ju006rer", "z_ju009rer", "z_ju016rer", "z_ju018rer")

ds_92 <- d |>
  select(id, z_mu003rer, z_mu004rer, z_mu005rer, z_mu006rer, z_mu007rer,
         z_mu008rer, z_mu009rer, z_mu010rer, z_mu011rer, z_mu012rer, z_mu013rer,
         z_mu014rer, z_mu015rer, z_mu016rer, z_mu017rer, z_mu018rer, z_mu019rer,
         z_mu020rer, z_mu021rer, z_mu022rer) |>
  pivot_longer(-id) |>
  mutate(recoded = case_when(value < 0 ~ NA_real_,
                             is.na(value) ~ NA_real_,
                             name %in% reverse & value < 1 ~ 3,
                             name %in% reverse & value <= 2 ~ 2,
                             name %in% reverse & value <= 4 ~ 1,
                             name %in% reverse ~ 0,
                             value < 1 ~ 0,
                             value <= 2 ~ 1,
                             value <= 4 ~ 2,
                             TRUE ~ 3),
         year = '92') |>
  na.omit() |>
  group_by(id, year) |>
  summarize(ds = sum(recoded, na.rm=TRUE), .groups='drop') |>
  left_join(resid_controls, by='id') |>
  na.omit()


ds_03 <- d |>
  select(id, z_iu003rer, z_iu004rer, z_iu005rer, z_iu006rer, z_iu007rer,
         z_iu008rer, z_iu009rer, z_iu010rer, z_iu011rer, z_iu012rer, z_iu013rer,
         z_iu014rer, z_iu015rer, z_iu016rer, z_iu017rer, z_iu018rer, z_iu019rer,
         z_iu020rer, z_iu021rer, z_iu022rer) |>
  pivot_longer(-id) |>
  mutate(recoded = case_when(value < 0 ~ NA_real_,
                             is.na(value) ~ NA_real_,
                             name %in% reverse & value < 1 ~ 3,
                             name %in% reverse & value <= 2 ~ 2,
                             name %in% reverse & value <= 4 ~ 1,
                             name %in% reverse ~ 0,
                             value < 1 ~ 0,
                             value <= 2 ~ 1,
                             value <= 4 ~ 2,
                             TRUE ~ 3),
         year = '03') |>
  na.omit() |>
  group_by(id, year) |>
  summarize(ds = sum(recoded, na.rm=TRUE), .groups='drop') |>
  left_join(resid_controls, by='id') |>
  na.omit()


ds_11 <- d |>
  select(id, z_ju003rer, z_ju004rer, z_ju005rer, z_ju006rer, z_ju007rer,
         z_ju008rer, z_ju009rer, z_ju010rer, z_ju011rer, z_ju012rer, z_ju013rer,
         z_ju014rer, z_ju015rer, z_ju016rer, z_ju017rer, z_ju018rer, z_ju019rer,
         z_ju020rer, z_ju021rer, z_ju022rer) |>
  pivot_longer(-id) |>
  mutate(recoded = case_when(value < 0 ~ NA_real_,
                             is.na(value) ~ NA_real_,
                             name %in% reverse & value < 1 ~ 3,
                             name %in% reverse & value <= 2 ~ 2,
                             name %in% reverse & value <= 4 ~ 1,
                             name %in% reverse ~ 0,
                             value < 1 ~ 0,
                             value <= 2 ~ 1,
                             value <= 4 ~ 2,
                             TRUE ~ 3),
         year = '11') |>
  na.omit() |>
  group_by(id, year) |>
  summarize(ds = sum(recoded, na.rm=TRUE), .groups='drop') |>
  left_join(resid_controls, by='id') |>
  na.omit()

if(resid_phen){
  m_ds_92 <- lm(ds ~ I(birthdate^2)*female + birthdate*female, ds_92)
  m_ds_03 <- lm(ds ~ I(birthdate^2)*female + birthdate*female, ds_03)
  m_ds_11 <- lm(ds ~ I(birthdate^2)*female + birthdate*female, ds_11)
} else {
  m_ds_92 <- lm(ds ~ 1, ds_92)
  m_ds_03 <- lm(ds ~ 1, ds_03)
  m_ds_11 <- lm(ds ~ 1, ds_11)
}

ds_92$ds_92_resid <- resid(m_ds_92)
ds_03$ds_03_resid <- resid(m_ds_03)
ds_11$ds_11_resid <- resid(m_ds_11)

depressive_symptoms_alt <- ds_92 |>
  full_join(ds_03, by=c('id', 'female', 'birthdate')) |>
  full_join(ds_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(depressive_symptoms_alt = mean(value, na.rm=TRUE))



# drinks per week

dpw <- d |>
  mutate(days_92 = if_else(z_ru026re >= 0, as.numeric(z_ru026re)/4, NA_real_),
         drinks_92 =if_else(z_ru027re >= 0, as.numeric(z_ru027re), NA_real_),
         days_03 = if_else(z_gu026re >= 0, as.numeric(z_gu026re)/4, NA_real_),
         drinks_03 =if_else(z_gu027re >= 0, as.numeric(z_gu027re), NA_real_),
         days_11 = if_else(z_hu026re >= 0, as.numeric(z_hu026re)/4, NA_real_),
         drinks_11 =if_else(z_hu027re >= 0, as.numeric(z_hu027re), NA_real_)) |>
  mutate(dpw_92 = days_92*drinks_92,
         dpw_03 = days_03*drinks_03,
         dpw_11 = days_11*drinks_11)

dpw_92 <- dpw |>
  select(id, birthdate, female, dpw_92) |>
  na.omit()

dpw_03 <- dpw |>
  select(id, birthdate, female, dpw_03) |>
  na.omit()

dpw_11 <- dpw |>
  select(id, birthdate, female, dpw_11) |>
  na.omit()

if(resid_phen){
  m_dpw_92 <- lm(dpw_92 ~ I(birthdate^2)*female + birthdate*female, dpw_92)
  m_dpw_03 <- lm(dpw_03 ~ I(birthdate^2)*female + birthdate*female, dpw_03)
  m_dpw_11 <- lm(dpw_11 ~ I(birthdate^2)*female + birthdate*female, dpw_11)
} else {
  m_dpw_92 <- lm(dpw_92 ~ 1, dpw_92)
  m_dpw_03 <- lm(dpw_03 ~ 1, dpw_03)
  m_dpw_11 <- lm(dpw_11 ~ 1, dpw_11)
}

dpw_92$dpw_92_resid <- resid(m_dpw_92)
dpw_03$dpw_03_resid <- resid(m_dpw_03)
dpw_11$dpw_11_resid <- resid(m_dpw_11)

drinks_per_week <- dpw_92 |>
  full_join(dpw_03, by=c('id', 'female', 'birthdate')) |>
  full_join(dpw_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(drinks_per_week = mean(value, na.rm=TRUE))


# ever smoker

ever_smoker <- d |>
  mutate(ever_smoker = case_when(z_mx012rer == 1 ~ 1,
                                 z_ix012rer == 1 ~ 1,
                                 z_jx012rer == 1 ~ 1,
                                 z_mx012rer == 2 ~ 0,
                                 z_ix012rer == 2 ~ 0,
                                 z_jx012rer == 2 ~ 0,
                                 TRUE ~ NA_real_)) |>
  select(id, ever_smoker) |>
  left_join(resid_controls, by='id') |>
  na.omit()

if(resid_phen){
  m_ever_smoker <- lm(ever_smoker ~ I(birthdate^2)*female + birthdate*female, ever_smoker)
} else {
  m_ever_smoker <- lm(ever_smoker ~ 1, ever_smoker)
}

ever_smoker$ever_smoker <- resid(m_ever_smoker)
ever_smoker <- select(ever_smoker, id, ever_smoker)

# physical activity

pa <- d |>
  mutate(light_92 = case_when(z_mx005rer == 1 ~ 15,
                              z_mx005rer == 2 ~ 6,
                              z_mx005rer == 3 ~ 2,
                              z_mx005rer == 4 ~ 0.5,
                              TRUE ~ NA_real_),
         heavy_92 = case_when(z_mx006rer == 1 ~ 15,
                              z_mx006rer == 2 ~ 6,
                              z_mx006rer == 3 ~ 2,
                              z_mx006rer == 4 ~ 0.5,
                              TRUE ~ NA_real_),
         light1_03 = if_else(z_iz165rer >= 0, as.numeric(z_iz165rer), NA_real_),
         light2_03 = if_else(z_iz168rer >= 0, as.numeric(z_iz168rer), NA_real_),
         heavy1_03 = if_else(z_iz171rer >= 0, as.numeric(z_iz171rer), NA_real_),
         heavy2_03 = if_else(z_iz174rer >= 0, as.numeric(z_iz174rer), NA_real_),
         light1_11 = if_else(z_jz165rer >= 0, as.numeric(z_jz165rer), NA_real_),
         light2_11 = if_else(z_jz168rer >= 0, as.numeric(z_jz168rer), NA_real_),
         heavy1_11 = if_else(z_jz171rer >= 0, as.numeric(z_jz171rer), NA_real_),
         heavy2_11 = if_else(z_jz174rer >= 0, as.numeric(z_jz174rer), NA_real_)) |>
  mutate(pa_92 = 2*light_92 + 8*heavy_92,
         pa_03 = 2*(light1_03+light2_03) + 8*(heavy1_03+heavy2_03),
         pa_11 = 2*(light1_11+light2_11) + 8*(heavy1_11+heavy2_11))


pa_92 <- pa |>
  select(id, birthdate, female, pa_92) |>
  na.omit()

pa_03 <- pa |>
  select(id, birthdate, female, pa_03) |>
  na.omit()

pa_11 <- pa |>
  select(id, birthdate, female, pa_11) |>
  na.omit()

if(resid_phen){
  m_pa_92 <- lm(pa_92 ~ I(birthdate^2)*female + birthdate*female, pa_92)
  m_pa_03 <- lm(pa_03 ~ I(birthdate^2)*female + birthdate*female, pa_03)
  m_pa_11 <- lm(pa_11 ~ I(birthdate^2)*female + birthdate*female, pa_11)
} else {
  m_pa_92 <- lm(pa_92 ~ 1, pa_92)
  m_pa_03 <- lm(pa_03 ~ 1, pa_03)
  m_pa_11 <- lm(pa_11 ~ 1, pa_11)
}

pa_92$pa_92_resid <- resid(m_pa_92)
pa_03$pa_03_resid <- resid(m_pa_03)
pa_11$pa_11_resid <- resid(m_pa_11)

physical_activity <- pa_92 |>
  full_join(pa_03, by=c('id', 'female', 'birthdate')) |>
  full_join(pa_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(physical_activity = mean(value, na.rm=TRUE))

# self rated health

srh <- d |>
  mutate(srh_92 = if_else(z_mx001rer > 0, as.numeric(z_mx001rer), NA_real_),
         srh_03 = if_else(z_ix001rer > 0, as.numeric(z_ix001rer), NA_real_),
         srh_11 = if_else(z_jx001rer > 0, as.numeric(z_jx001rer), NA_real_))


srh_92 <- srh |>
  select(id, birthdate, female, srh_92) |>
  na.omit()

srh_03 <- srh |>
  select(id, birthdate, female, srh_03) |>
  na.omit()

srh_11 <- srh |>
  select(id, birthdate, female, srh_11) |>
  na.omit()

if(resid_phen){
  m_srh_92 <- lm(srh_92 ~ I(birthdate^2)*female + birthdate*female, srh_92)
  m_srh_03 <- lm(srh_03 ~ I(birthdate^2)*female + birthdate*female, srh_03)
  m_srh_11 <- lm(srh_11 ~ I(birthdate^2)*female + birthdate*female, srh_11)
} else {
  m_srh_92 <- lm(srh_92 ~ 1, srh_92)
  m_srh_03 <- lm(srh_03 ~ 1, srh_03)
  m_srh_11 <- lm(srh_11 ~ 1, srh_11)
}

srh_92$srh_92_resid <- resid(m_srh_92)
srh_03$srh_03_resid <- resid(m_srh_03)
srh_11$srh_11_resid <- resid(m_srh_11)

self_rated_health <- srh_92 |>
  full_join(srh_03, by=c('id', 'female', 'birthdate')) |>
  full_join(srh_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(self_rated_health = mean(value, na.rm=TRUE))


##############################
# personality and well-being #
##############################

# big 5: mean impute missing components

b5 <- d |>
  mutate(extra_92 = if_else(z_mh001rei >= 0, as.numeric(z_mh001rei), NA_real_),
         extra_03 = if_else(z_ih001rei >= 0, as.numeric(z_ih001rei), NA_real_),
         extra_11 = if_else(z_jh001rei >= 0, as.numeric(z_jh001rei), NA_real_),
         neuro_92 = if_else(z_mh025rei >= 0, as.numeric(z_mh025rei), NA_real_),
         neuro_03 = if_else(z_ih025rei >= 0, as.numeric(z_ih025rei), NA_real_),
         neuro_11 = if_else(z_jh025rei >= 0, as.numeric(z_jh025rei), NA_real_),
         open_92 = if_else(z_mh032rei >= 0, as.numeric(z_mh032rei), NA_real_),
         open_03 = if_else(z_ih032rei >= 0, as.numeric(z_ih032rei), NA_real_),
         open_11 = if_else(z_jh032rei >= 0, as.numeric(z_jh032rei), NA_real_))

## extraversion

extra_92 <- b5 |>
  select(id, birthdate, female, extra_92) |>
  na.omit()

extra_03 <- b5 |>
  select(id, birthdate, female, extra_03) |>
  na.omit()

extra_11 <- b5 |>
  select(id, birthdate, female, extra_11) |>
  na.omit()

if(resid_phen){
  m_extra_92 <- lm(extra_92 ~ I(birthdate^2)*female + birthdate*female, extra_92)
  m_extra_03 <- lm(extra_03 ~ I(birthdate^2)*female + birthdate*female, extra_03)
  m_extra_11 <- lm(extra_11 ~ I(birthdate^2)*female + birthdate*female, extra_11)
} else {
  m_extra_92 <- lm(extra_92 ~ 1, extra_92)
  m_extra_03 <- lm(extra_03 ~ 1, extra_03)
  m_extra_11 <- lm(extra_11 ~ 1, extra_11)
}

extra_92$extra_92_resid <- resid(m_extra_92)
extra_03$extra_03_resid <- resid(m_extra_03)
extra_11$extra_11_resid <- resid(m_extra_11)

extra <- extra_92 |>
  full_join(extra_03, by=c('id', 'female', 'birthdate')) |>
  full_join(extra_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(extraversion = mean(value, na.rm=TRUE))

## neuroticism

neuro_92 <- b5 |>
  select(id, birthdate, female, neuro_92) |>
  na.omit()

neuro_03 <- b5 |>
  select(id, birthdate, female, neuro_03) |>
  na.omit()

neuro_11 <- b5 |>
  select(id, birthdate, female, neuro_11) |>
  na.omit()


if(resid_phen){
  m_neuro_92 <- lm(neuro_92 ~ I(birthdate^2)*female + birthdate*female, neuro_92)
  m_neuro_03 <- lm(neuro_03 ~ I(birthdate^2)*female + birthdate*female, neuro_03)
  m_neuro_11 <- lm(neuro_11 ~ I(birthdate^2)*female + birthdate*female, neuro_11)
} else {
  m_neuro_92 <- lm(neuro_92 ~ 1, neuro_92)
  m_neuro_03 <- lm(neuro_03 ~ 1, neuro_03)
  m_neuro_11 <- lm(neuro_11 ~ 1, neuro_11)
}

neuro_92$neuro_92_resid <- resid(m_neuro_92)
neuro_03$neuro_03_resid <- resid(m_neuro_03)
neuro_11$neuro_11_resid <- resid(m_neuro_11)

neuro <- neuro_92 |>
  full_join(neuro_03, by=c('id', 'female', 'birthdate')) |>
  full_join(neuro_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(neuroticism = mean(value, na.rm=TRUE))

## openness

open_92 <- b5 |>
  select(id, birthdate, female, open_92) |>
  na.omit()

open_03 <- b5 |>
  select(id, birthdate, female, open_03) |>
  na.omit()

open_11 <- b5 |>
  select(id, birthdate, female, open_11) |>
  na.omit()

if(resid_phen){
  m_open_92 <- lm(open_92 ~ I(birthdate^2)*female + birthdate*female, open_92)
  m_open_03 <- lm(open_03 ~ I(birthdate^2)*female + birthdate*female, open_03)
  m_open_11 <- lm(open_11 ~ I(birthdate^2)*female + birthdate*female, open_11)
} else {
  m_open_92 <- lm(open_92 ~ 1, open_92)
  m_open_03 <- lm(open_03 ~ 1, open_03)
  m_open_11 <- lm(open_11 ~ 1, open_11)
}

open_92$open_92_resid <- resid(m_open_92)
open_03$open_03_resid <- resid(m_open_03)
open_11$open_11_resid <- resid(m_open_11)

open <- open_92 |>
  full_join(open_03, by=c('id', 'female', 'birthdate')) |>
  full_join(open_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(openness = mean(value, na.rm=TRUE))

## combine big 5

big_5 <- extra |>
  full_join(neuro, by='id') |>
  full_join(open, by='id')


# life satisfaction: family, finance, work

ls <- d |>
  mutate(life_satisfaction_family = if_else(z_gb040re > 0, as.numeric(z_gb040re), NA_real_),
         ls_finance_03 = if_else(z_gp226re > 0, as.numeric(z_gp226re), NA_real_),
         ls_finance_11 = if_else(z_hp226re > 0, as.numeric(z_hp226re), NA_real_),
         ls_job_92 = if_else(z_rg044jjc > 0, as.numeric(z_rg044jjc), NA_real_),
         ls_job_03 = if_else(z_gg044jjc > 0, as.numeric(z_gg044jjc), NA_real_),
         ls_job_11 = if_else(z_hg044jjc > 0, as.numeric(z_hg044jjc), NA_real_))

## family

ls_family <- ls |>
  select(id, birthdate, female, life_satisfaction_family) |>
  na.omit()

if(resid_phen){
  m_ls_family <- lm(life_satisfaction_family ~ I(birthdate^2)*female + birthdate*female, ls_family)
} else {
  m_ls_family <- lm(life_satisfaction_family ~ 1, ls_family)
}

ls_family$life_satisfaction_family <- resid(m_ls_family)
ls_family <- select(ls_family, id, life_satisfaction_family)

## work

ls_job_92 <- ls |>
  select(id, birthdate, female, ls_job_92) |>
  na.omit()

ls_job_03 <- ls |>
  select(id, birthdate, female, ls_job_03) |>
  na.omit()

ls_job_11 <- ls |>
  select(id, birthdate, female, ls_job_11) |>
  na.omit()

if(resid_phen){
  m_ls_job_92 <- lm(ls_job_92 ~ I(birthdate^2)*female + birthdate*female, ls_job_92)
  m_ls_job_03 <- lm(ls_job_03 ~ I(birthdate^2)*female + birthdate*female, ls_job_03)
  m_ls_job_11 <- lm(ls_job_11 ~ I(birthdate^2)*female + birthdate*female, ls_job_11)
} else {
  m_ls_job_92 <- lm(ls_job_92 ~ 1, ls_job_92)
  m_ls_job_03 <- lm(ls_job_03 ~ 1, ls_job_03)
  m_ls_job_11 <- lm(ls_job_11 ~ 1, ls_job_11)
}

ls_job_92$ls_job_92_resid <- resid(m_ls_job_92)
ls_job_03$ls_job_03_resid <- resid(m_ls_job_03)
ls_job_11$ls_job_11_resid <- resid(m_ls_job_11)

ls_job <- ls_job_92 |>
  full_join(ls_job_03, by=c('id', 'female', 'birthdate')) |>
  full_join(ls_job_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(life_satisfaction_job = mean(value, na.rm=TRUE))

## finances

ls_finance_03 <- ls |>
  select(id, birthdate, female, ls_finance_03) |>
  na.omit()

ls_finance_11 <- ls |>
  select(id, birthdate, female, ls_finance_11) |>
  na.omit()

if(resid_phen){
  m_ls_finance_03 <- lm(ls_finance_03 ~ I(birthdate^2)*female + birthdate*female, ls_finance_03)
  m_ls_finance_11 <- lm(ls_finance_11 ~ I(birthdate^2)*female + birthdate*female, ls_finance_11)
} else {
  m_ls_finance_03 <- lm(ls_finance_03 ~ 1, ls_finance_03)
  m_ls_finance_11 <- lm(ls_finance_11 ~ 1, ls_finance_11)
}

ls_finance_03$ls_finance_03_resid <- resid(m_ls_finance_03)
ls_finance_11$ls_finance_11_resid <- resid(m_ls_finance_11)

ls_finance <- ls_finance_03 |>
  full_join(ls_finance_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(life_satisfaction_finance = mean(value, na.rm=TRUE))

## combine

life_satisfaction <- ls_family |>
  full_join(ls_finance, by='id') |>
  full_join(ls_job, by='id') |>
  select(id, starts_with('life_satisfaction'))


# loneliness

lone <- d |>
  mutate(lone_92 = if_else(z_mu008rer >= 0, as.numeric(z_mu008rer), NA_real_),
         lone_03 = if_else(z_iu008rer >= 0, as.numeric(z_iu008rer), NA_real_),
         lone_11 = if_else(z_ju008rer >= 0, as.numeric(z_ju008rer), NA_real_))

lone_92 <- lone |>
  select(id, birthdate, female, lone_92) |>
  na.omit()

lone_03 <- lone |>
  select(id, birthdate, female, lone_03) |>
  na.omit()

lone_11 <- lone |>
  select(id, birthdate, female, lone_11) |>
  na.omit()

if(resid_phen){
  m_lone_92 <- lm(lone_92 ~ I(birthdate^2)*female + birthdate*female, lone_92)
  m_lone_03 <- lm(lone_03 ~ I(birthdate^2)*female + birthdate*female, lone_03)
  m_lone_11 <- lm(lone_11 ~ I(birthdate^2)*female + birthdate*female, lone_11)
} else {
  m_lone_92 <- lm(lone_92 ~ 1, lone_92)
  m_lone_03 <- lm(lone_03 ~ 1, lone_03)
  m_lone_11 <- lm(lone_11 ~ 1, lone_11)
}

lone_92$lone_92_resid <- resid(m_lone_92)
lone_03$lone_03_resid <- resid(m_lone_03)
lone_11$lone_11_resid <- resid(m_lone_11)

loneliness <- lone_92 |>
  full_join(lone_03, by=c('id', 'female', 'birthdate')) |>
  full_join(lone_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(loneliness = mean(value, na.rm=TRUE))

# religious attendance

ra <- d |>
  mutate(ra_75 = case_when(z_bkxrl3 < 0 ~ NA_real_,
                           z_bkxrl3 == 6 ~ 5,
                           TRUE ~ as.numeric(z_bkxrl3)),
         ra_11 = case_when(z_jl004rec < 0 ~ NA_real_,
                           z_jl004rec <= 1 ~ 5,
                           z_jl004rec <= 3 ~ 4,
                           z_jl004rec == 4 ~ 3,
                           z_jl004rec <= 6 ~ 2,
                           TRUE ~ 1))

ra_75 <- ra |>
  select(id, birthdate, female, ra_75) |>
  na.omit()

ra_11 <- ra |>
  select(id, birthdate, female, ra_11) |>
  na.omit()

if(resid_phen){
  m_ra_75 <- lm(ra_75 ~ I(birthdate^2)*female + birthdate*female, ra_75)
  m_ra_11 <- lm(ra_11 ~ I(birthdate^2)*female + birthdate*female, ra_11)
} else {
  m_ra_75 <- lm(ra_75 ~ 1, ra_75)
  m_ra_11 <- lm(ra_11 ~ 1, ra_11)
}

ra_75$ra_75_resid <- resid(m_ra_75)
ra_11$ra_11_resid <- resid(m_ra_11)

religious_attendance <- ra_75 |>
  full_join(ra_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(religious_attendance = mean(value, na.rm=TRUE))


# risk tolerance

risk_tolerance <- d |>
  select(id,
         z_jstk02re, z_jstk03re, z_jstk04re, z_jstk05re, z_jstk06re, z_jstk07re,
         z_jstk08re, z_jstk09re, z_jstk10re, z_jstk11re, z_jstk12re, z_jstk13re,
         z_jstk14re, z_jstk15re, z_jstk16re, z_jstk17re, z_jstk18re, z_jstk19re,
         z_jstk20re, z_jstk21re, z_jstk22re) |>
  pivot_longer(-id) |>
  mutate(resp = case_when(value == 1 ~ 1,
                          value == 2 ~ 0,
                          TRUE ~ NA_real_)) |>
  group_by(id) |>
  summarize(risk = mean(resp, na.rm=TRUE)) |>
  full_join(resid_controls, by='id')

risk <- na.omit(risk_tolerance)

if(resid_phen){
  m_risk <- lm(risk ~ I(birthdate^2)*female + birthdate*female, risk)
} else {
  m_risk <- lm(risk ~ 1, risk)
}

risk$risk_tolerance <- resid(m_risk)
risk_tolerance <- select(risk, id, risk_tolerance)

# subjective well-being

swb <- d |>
  mutate(happy_92 = if_else(z_mu006rer >= 0, as.numeric(z_mu006rer), NA_real_),
         happy_03 = if_else(z_iu006rer >= 0, as.numeric(z_iu006rer), NA_real_),
         happy_11 = if_else(z_ju006rer >= 0, as.numeric(z_ju006rer), NA_real_),
         enjoy_92 = if_else(z_mu009rer >= 0, as.numeric(z_mu009rer), NA_real_),
         enjoy_03 = if_else(z_iu009rer >= 0, as.numeric(z_iu009rer), NA_real_),
         enjoy_11 = if_else(z_ju009rer >= 0, as.numeric(z_ju009rer), NA_real_)) |>
  select(id, starts_with('happy_'), starts_with('enjoy_')) |>
  pivot_longer(-id, names_to=c('item', 'year'), names_sep = 6) |>
  group_by(id, year) |>
  summarize(swb = mean(value, na.rm=TRUE), .groups='drop') |>
  full_join(resid_controls, by='id')

swb_92 <- swb |>
  filter(year == '92') |>
  na.omit()

swb_03 <- swb |>
  filter(year == '03') |>
  na.omit()

swb_11 <- swb |>
  filter(year == '11') |>
  na.omit()

if(resid_phen){
  m_swb_92 <- lm(swb ~ I(birthdate^2)*female + birthdate*female, swb_92)
  m_swb_03 <- lm(swb ~ I(birthdate^2)*female + birthdate*female, swb_03)
  m_swb_11 <- lm(swb ~ I(birthdate^2)*female + birthdate*female, swb_11)
} else {
  m_swb_92 <- lm(swb ~ 1, swb_92)
  m_swb_03 <- lm(swb ~ 1, swb_03)
  m_swb_11 <- lm(swb ~ 1, swb_11)
}

swb_92$swb_92_resid <- resid(m_swb_92)
swb_03$swb_03_resid <- resid(m_swb_03)
swb_11$swb_11_resid <- resid(m_swb_11)

subjective_well_being <- swb_92 |>
  full_join(swb_03, by=c('id', 'female', 'birthdate')) |>
  full_join(swb_11, by=c('id', 'female', 'birthdate')) |>
  select(id, ends_with('resid')) |>
  pivot_longer(-id) |>
  group_by(id) |>
  summarize(subjective_well_being = mean(value, na.rm=TRUE))



####################
# construct output #
####################

out <- id |>
  full_join(bmi, by='id') |>
  full_join(height, by='id') |>
  full_join(cognitive_performance, by='id') |>
  full_join(ed_attainment, by='id') |>
  full_join(age_first_birth, by='id') |>
  full_join(age_first_menses, by='id') |>
  full_join(number_ever_born, by='id') |>
  full_join(alcohol_misuse, by='id') |>
  full_join(allergies, by='id') |>
  full_join(asthma, by='id') |>
  #full_join(adhd, by='id') |>
  full_join(cigarettes_per_day, by='id') |>
  full_join(copd, by='id') |>
  full_join(depressive_symptoms, by='id') |>
  full_join(depressive_symptoms_alt, by='id') |>
  full_join(drinks_per_week, by='id') |>
  full_join(ever_smoker, by='id') |>
  full_join(physical_activity, by='id') |>
  full_join(self_rated_health, by='id') |>
  full_join(big_5, by='id') |>
  full_join(life_satisfaction, by='id') |>
  full_join(loneliness, by='id') |>
  full_join(religious_attendance, by='id') |>
  full_join(risk_tolerance, by='id') |>
  full_join(subjective_well_being, by='id')

write.csv(out, './data/clean/wls_phenotypes.csv')
