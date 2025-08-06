library(gtsummary)
library(gt)
meta <- read.csv("data/metadata/metadata_c.csv")
# Defining variables
meta1 <- meta %>%
  mutate(status = dplyr::recode(status,
                                sc ="Pre-HIV",
                                nc = "Non-HIV")) %>%
  mutate(Drinking = dplyr::recode(DRINKCAT,
    `3` = "Heavy drinking",
    `2` = "Heavy drinking",
    `1` = "Not heavy drinking",
    `0` = "Not heavy drinking"
  )) %>%
  mutate(Smoking = dplyr::recode(SMOKE,
    `1` = "Never",
    `2` = "Former",
    `3` = "Current"
  )) %>%
  mutate(Education = dplyr::recode(EDUCA,
    `1` = "No degree",
    `2` = "No degree",
    `3` = "No degree",
    `4` = "No degree",
    `5` = "Undergrad",
    `6` = "Postgrad",
    `7` = "Postgrad"
  )) %>%
  mutate(Race = dplyr::recode(RACE,
    `1` = "White",
    `2` = "White",
    `3` = "Black",
    `4` = "Black",
    .default = "Other"
  )) %>%
  mutate(Location = dplyr::recode(loc,
    `1` = "Baltimore",
    `2` = "Chicago",
    `3` = "Pittsburgh",
    `4` = "Los Angeles"
  )) %>%
  mutate(Substance_use = dplyr::recode(substance_use,
    substance_use = "yes",
    no_substance_use = "no"
  )) %>%
  mutate(sexual_activity = dplyr::recode(group1,
    g1 = "0",
    g2 = "1", g3 = "2-5", g4 = "6 or more"
  ))

# Summary table
tbl <- meta1 %>%
  tbl_summary(
    by = status,
    include = c("status", "age", "Location", "Race", "Education", "Drinking", "Smoking", "abx_use", "Substance_use", "sexual_activity"),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})"
    ),
    digits = list(
      all_categorical() ~ 1,
      all_continuous() ~ 1
    ),
    label = list(
      age ~ "Age",
      Substance_use ~ "Any substance use",
      abx_use ~ "Antibiotic use",
      sexual_activity ~ "Number of receptive anal intercourse partners"
    )
  ) %>%
  add_p(list(all_continuous() ~ "wilcox.test", all_categorical() ~ "fisher.test"))


tbl %>%
  as_gt() %>%
  gt::gtsave("Outputs/Summary_Table.docx")
