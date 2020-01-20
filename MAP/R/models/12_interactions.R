##interactions


load(file.path(mpath,'scantwo.scans.elr.short.rsave'))
elr_gg <- gg
elr_gg_step2 <- gg_step2
elr_sc2_normal_imp <- sc2_normal_imp
elr_sc2_normal_imp_perms <- sc2_normal_imp_perms

elr_int <- summary(elr_sc2_normal_imp, thresholds=c(0, Inf, 4, Inf, Inf), what="int")


load(file.path(mpath,'scantwo.scans.nbh.short.rsave'))
nbh_gg <- gg
nbh_gg_step2 <- gg_step2
nbh_sc2_normal_imp <- sc2_normal_imp
nbh_sc2_normal_imp_perms <- sc2_normal_imp_perms

nbh_int <- summary(nbh_sc2_normal_imp, thresholds=c(0, Inf, 4, Inf, Inf), what="int")
