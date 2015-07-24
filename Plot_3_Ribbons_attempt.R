load("Latest_workspace.RData")
source("Output_ribbon_data_for_3_datasets.R")
require(ggplot2)
source("Data_NZoo_week13.R")
NZoo.cases.report1 = cases.by.age.group.5.group

source("Data_NZoo.R")
NZoo.cases.report2 = cases.by.age.group.5.group


NZoo.ribbons = ribbon.plot.NZoo.3.sims(NZoo.sims1, NZoo.cases.report1, 13,
                                       NZoo.sims2, NZoo.cases.report2[1:10, ], 13,
                                       NZoo.sims3, NZoo.cases.report2, 17,
                                       4)
NZoo.ribbons


pdf("3ribbons.pdf", height = 5, width = 6)
NZoo.ribbons
dev.off()



Foum.ribbons <- ribbon.plot.Non.NZoo.3.sims(Foum.sims1.bounded, Foum.ests1.bounded$Cases_by_age, 13,
                                            Foum.sims2.bounded, Foum.ests2.bounded$Cases_by_age, 13,
                                            Foum.sims3.bounded, Foum.ests3.bounded$Cases_by_age, 17,
                                            10, "Foumbadou")

Foum.ribbons
Gueasso.ribbons <- ribbon.plot.Non.NZoo.3.sims(Gueasso.sims1.bounded, Gueasso.ests1.bounded$Cases_by_age, 13,
                                               Gueasso.sims2.bounded, Gueasso.ests2.bounded$Cases_by_age, 13,
                                               Gueasso.sims3.bounded, Gueasso.ests3.bounded$Cases_by_age, 17,
                                               11, "Gueasso")


Gama.ribbons <- ribbon.plot.Non.NZoo.3.sims(Gama.sims1.bounded, Gama.ests1.bounded$Cases_by_age, 13,
                                            Gama.sims2.bounded, Gama.ests2.bounded$Cases_by_age, 13,
                                            Gama.sims3.bounded, Gama.ests3.bounded$Cases_by_age, 17,
                                            11, "Gama Berema")


Kokota.ribbons <- ribbon.plot.Non.NZoo.3.sims(Kokota.sims1.bounded, Kokota.ests1.bounded$Cases_by_age, 13,
                                              Kokota.sims2.bounded, Kokota.ests2.bounded$Cases_by_age, 13,
                                              Kokota.sims3.bounded, Kokota.ests3.bounded$Cases_by_age, 17,
                                              9,"Kokota")


Laine.ribbons <- ribbon.plot.Non.NZoo.3.sims(Laine.sims1.bounded, Laine.ests1.bounded$Cases_by_age, 13,
                                             Laine.sims2.bounded, Laine.ests2.bounded$Cases_by_age, 13,
                                             Laine.sims3.bounded, Laine.ests3.bounded$Cases_by_age, 17,
                                             12, "Laine")

Cu.Lola.ribbons <- ribbon.plot.Non.NZoo.3.sims(Cu.Lola.sims1.bounded, Cu.Lola.ests1.bounded$Cases_by_age, 13,
                                               Cu.Lola.sims2.bounded, Cu.Lola.ests2.bounded$Cases_by_age, 13,
                                               Cu.Lola.sims3.bounded, Cu.Lola.ests3.bounded$Cases_by_age, 17,
                                               12, "Cu Lola")

# Cu.Lola.ribbons.no.obs.bound <- ribbon.plot.Non.NZoo.3.sims(Cu.Lola.sims1.bounded, Cu.Lola.ests1.bounded$Cases_by_age, 13,
#                                                Cu.Lola.sims2.bounded, Cu.Lola.ests2.bounded$Cases_by_age, 13,
#                                                Cu.Lola.sims3.bounded, Cu.Lola.est.bounded.case.obs3$Cases_by_age, 17,
#                                                12, "Cu Lola")


pdf("Compare_predictions_ribbons_3_datasets_July_23.pdf", width = 13, height = 15)
multiplot(NZoo.ribbons, Gueasso.ribbons, Kokota.ribbons, Cu.Lola.ribbons,
          Foum.ribbons, Gama.ribbons, Laine.ribbons, cols=2)
dev.off()

require(grid)
