
compexpt <- function(df) {
  decline_list <- unique(setDT(bene_lvlBU)[!CRNT_NUM %in% cclf8$CRNT_NUM, CRNT_NUM])

  if (any(colnames(df) == "CLM_PMT_AMT")) {
    temp <- xref2[df, on = c(PRVS_NUM = "BENE_MBI_ID"), .(BENE_MBI_ID,
      CRNT_NUM = if_else(is.na(CRNT_NUM), BENE_MBI_ID, CRNT_NUM),
      CLM_FROM_DT, CLM_PMT_AMT,
      CLM_MNTH = paste(str_sub(month.name[as.integer(str_sub(CLM_FROM_DT, start = 6L, end = 7L))], 1L, 3L),
        str_sub(CLM_FROM_DT, start = 3L, end = 4L),
        sep = "_"
      )
    )][CRNT_NUM %in% bene_lvlBU$CRNT_NUM,
      .(tot_exp = sum(CLM_PMT_AMT)),
      by = .(CRNT_NUM, CLM_MNTH)
    ]
  } else {
    temp <- xref2[df, on = c(PRVS_NUM = "BENE_MBI_ID"), .(BENE_MBI_ID,
      CRNT_NUM = if_else(is.na(CRNT_NUM), BENE_MBI_ID, CRNT_NUM),
      CLM_FROM_DT, CLM_LINE_CVRD_PD_AMT,
      CLM_MNTH = paste(str_sub(month.name[as.integer(str_sub(CLM_FROM_DT, start = 6L, end = 7L))], 1L, 3L),
        str_sub(CLM_FROM_DT, start = 3L, end = 4L),
        sep = "_"
      ),
      CRNT_NUM = if_else(is.na(CRNT_NUM), BENE_MBI_ID, CRNT_NUM)
    )][CRNT_NUM %in% bene_lvlBU$CRNT_NUM,
      .(tot_exp = sum(CLM_LINE_CVRD_PD_AMT)),
      by = .(CRNT_NUM, CLM_MNTH)
    ]
  }


  cost_temp <- merge(temp, bene_lvlBU, by = c("CRNT_NUM", "CLM_MNTH"), all = TRUE)[
    ,

    `:=`(eligibility = case_when(
      attrib == 1 ~ "ESRD",
      attrib == 2 ~ "Disabled",
      attrib == 3 ~ "AgedDual",
      attrib == 4 ~ "AgedNonDual",
      attrib == 0 ~ "Not_Eligible"
    ))
  ][,
    .(
      tot_exp = sum(tot_exp, na.rm = TRUE),
      ESRD = sum(attrib == 1),
      Disabled = sum(attrib == 2),
      AgedDual = sum(attrib == 3),
      AgedNonDual = sum(attrib == 4),
      Not_Eligible = sum(attrib == 0)
    ),
    by = .(CRNT_NUM, eligibility, TIN, CHAPTER)
  ][,
    membermonths := Reduce(`+`, lapply(.SD, function(x) {
      replace(
        x,
        which(is.na(x)), 0
      )
    })),
    .SDcols = c("ESRD", "Disabled", "AgedDual", "AgedNonDual", "Not_Eligible")
  ][
    ,
    `:=`(val = paste(tot_exp, membermonths, sep = "_"))
  ]

  cost <- cost_temp[, dcast(cost_temp, CRNT_NUM + TIN + CHAPTER ~ eligibility, value.var = "val")] %>%
    separate(ESRD, into = c("ESRD_exp", "ESRD.p"), sep = "_") %>%
    separate(Disabled, into = c("Disabled_exp", "Disabled.p"), sep = "_") %>%
    separate(AgedDual, into = c("AgedDual_exp", "AgedDual.p"), sep = "_") %>%
    separate(AgedNonDual, into = c("AgedNonDual_exp", "AgedNonDual.p"), sep = "_") %>%
    select(CRNT_NUM, TIN, CHAPTER, ESRD.p, Disabled.p, AgedDual.p, AgedNonDual.p, ESRD_exp, Disabled_exp, AgedDual_exp, AgedNonDual_exp) %>%
    mutate_at(vars(ESRD.p:AgedNonDual_exp), list(~ replace_na(., 0))) %>%
    mutate_at(vars(ESRD.p:AgedNonDual_exp), list(~ as.numeric(.))) %>%
    mutate_at(vars(ESRD.p:AgedNonDual.p), list(~ round(. / 12, 10))) %>%

    # annualize

    mutate(
      annl.ESRD = ifelse(ESRD.p == 0, 0, ESRD_exp / ESRD.p),
      annl.Disabled = ifelse(Disabled.p == 0, 0, Disabled_exp / Disabled.p),
      annl.AgedDual = ifelse(AgedDual.p == 0, 0, AgedDual_exp / AgedDual.p),
      annl.AgedNonDual = ifelse(AgedNonDual.p == 0, 0, AgedNonDual_exp / AgedNonDual.p),


      # truncate and apply completion factor
      ESRD.cost = ifelse(abs(annl.ESRD) > ESRD.Tr, sign(annl.ESRD) * ESRD.Tr * cmplt, annl.ESRD * cmplt),
      Disabled.cost = ifelse(abs(annl.Disabled) > DSBLD.Tr, sign(annl.Disabled) * DSBLD.Tr * cmplt, annl.Disabled * cmplt),
      AgedDual.cost = ifelse(abs(annl.AgedDual) > AGD.Tr, sign(annl.AgedDual) * AGD.Tr * cmplt, annl.AgedDual * cmplt),
      AgedNonDual.cost = ifelse(abs(annl.AgedNonDual) > AGND.Tr, sign(annl.AgedNonDual) * AGND.Tr * cmplt, annl.AgedNonDual * cmplt),
      Eligibility = rowSums(select(., ESRD.p:AgedNonDual.p)),

      # truncate and apply completion factor
      
      # ESRD.cost = annl.ESRD * (1 - ESRD.Tr) * cmplt,
      # Disabled.cost = annl.Disabled * (1 - DSBLD.Tr) * cmplt,
      # AgedDual.cost = annl.AgedDual * (1 - AGD.Tr) * cmplt,
      # AgedNonDual.cost = annl.AgedNonDual * (1 - AGND.Tr) * cmplt,
      # Eligibility = rowSums(select(., ESRD.p:AgedNonDual.p)),
      
      
      ESRD.cost.wt = ESRD.cost * ESRD.p, # weighted by enrollment proportion
      Disabled.cost.wt = Disabled.cost * Disabled.p,
      AgedDual.cost.wt = AgedDual.cost * AgedDual.p,
      AgedNonDual.cost.wt = AgedNonDual.cost * AgedNonDual.p
    ) %>%
    filter(Eligibility != 0) %>%
    mutate(
      wtcost = rowSums(select(., ESRD.cost.wt:AgedNonDual.cost.wt), na.rm = TRUE) / Eligibility,
      decline = ifelse(CRNT_NUM %in% decline_list, 1, 0)
    ) %>% # annualized & truncated totals)

    group_by(CHAPTER) %>%
    summarise(
      decline_yrs_ESRD = sum(ESRD.p[decline == 1], na.rm = TRUE),
      decline_yrs_DIS = sum(Disabled.p[decline == 1], na.rm = TRUE),
      decline_yrs_AGD = sum(AgedDual.p[decline == 1], na.rm = TRUE),
      decline_yrs_AGND = sum(AgedNonDual.p[decline == 1], na.rm = TRUE),

      ESRD_yrs = sum(ESRD.p, na.rm = TRUE) - decline_yrs_ESRD,
      Disabled_yrs = sum(Disabled.p, na.rm = TRUE) - decline_yrs_DIS,
      AgedDual_yrs = sum(AgedDual.p, na.rm = TRUE) - decline_yrs_AGD,
      AgedNonDual_yrs = sum(AgedNonDual.p, na.rm = TRUE) - decline_yrs_AGND,

      years = ESRD_yrs + Disabled_yrs + AgedDual_yrs + AgedNonDual_yrs,

      p.ESRD = ESRD_yrs / years, # proportions
      p.Disabled = Disabled_yrs / years,
      p.AgedDual = AgedDual_yrs / years,
      p.AgedNonDual = AgedNonDual_yrs / years,


      tot.ESRD = ifelse(ESRD_yrs != 0, sum(ESRD.cost.wt, na.rm = TRUE), 0),
      adj_esrd = (tot.ESRD + ((tot.ESRD / ESRD_yrs) * decline_yrs_ESRD * exp_ratio)) / (ESRD_yrs + decline_yrs_ESRD), # adjust for missing

      tot.Disabled = ifelse(Disabled_yrs != 0, sum(Disabled.cost.wt, na.rm = TRUE), 0),
      adj_dis = (tot.Disabled + ((tot.Disabled / Disabled_yrs) * decline_yrs_DIS * exp_ratio)) / (Disabled_yrs + decline_yrs_DIS),

      tot.AgedDual = ifelse(AgedDual_yrs != 0, sum(AgedDual.cost.wt, na.rm = TRUE), 0),
      adj_agd = (tot.AgedDual + ((tot.AgedDual / AgedDual_yrs) * decline_yrs_AGD * exp_ratio)) / (AgedDual_yrs + decline_yrs_AGD),

      tot.AgedNonDual = ifelse(AgedNonDual_yrs != 0, sum(AgedNonDual.cost.wt, na.rm = TRUE), 0),
      adj_agnd = (tot.AgedNonDual + ((tot.AgedNonDual / AgedNonDual_yrs) * decline_yrs_AGND * exp_ratio)) / (AgedNonDual_yrs + decline_yrs_AGND),

      comp_exp = 1.00 * (adj_esrd * p.ESRD +
        adj_dis * p.Disabled +
        adj_agd * p.AgedDual +
        adj_agnd * p.AgedNonDual)
    ) %>%
    select(CHAPTER, comp_exp) %>%
    spread(key = CHAPTER, value = comp_exp)

  return(cost)
}
