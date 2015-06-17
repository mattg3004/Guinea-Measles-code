load("June_16.RData")
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



pdf("3ribbons.pdf", height = 7, width = 7)
NZoo.ribbons
dev.off()



Foum.ribbons <- ribbon.plot.Non.NZoo.3.sims(Foum.sims1, Foum.ests1$Cases_by_age, 13,
                                            Foum.sims2, Foum.ests2$Cases_by_age, 13,
                                            Foum.sims3, Foum.ests3$Cases_by_age, 17,
                                            10)


Gueasso.ribbons <- ribbon.plot.Non.NZoo.3.sims(Gueasso.sims1, Gueasso.ests1$Cases_by_age, 13,
                                               Gueasso.sims2, Gueasso.ests2$Cases_by_age, 13,
                                               Gueasso.sims3, Gueasso.ests3$Cases_by_age, 17,
                                               11)




Gama.ribbons <- ribbon.plot.Non.NZoo.3.sims(Gama.sims1, Gama.ests1$Cases_by_age, 13,
                                            Gama.sims2, Gama.ests2$Cases_by_age, 13,
                                            Gama.sims3, Gama.ests3$Cases_by_age, 17,
                                            11)


Kokota.ribbons <- ribbon.plot.Non.NZoo.3.sims(Kokota.sims1, Kokota.ests1$Cases_by_age, 13,
                                              Kokota.sims2, Kokota.ests2$Cases_by_age, 13,
                                              Kokota.sims3, Kokota.ests3$Cases_by_age, 17,
                                              9)


Laine.ribbons <- ribbon.plot.Non.NZoo.3.sims(Laine.sims1, Laine.ests1$Cases_by_age, 13,
                                             Laine.sims2, Laine.ests2$Cases_by_age, 13,
                                             Laine.sims3, Laine.ests3$Cases_by_age, 17,
                                             12)

Cu.Lola.ribbons <- ribbon.plot.Non.NZoo.3.sims(Cu.Lola.sims1, Cu.Lola.ests1$Cases_by_age, 13,
                                               Cu.Lola.sims2, Cu.Lola.ests2$Cases_by_age, 13,
                                               Cu.Lola.sims3, Cu.Lola.ests3$Cases_by_age, 17,
                                               12)
Cu.Lola.ribbons +  geom_text(data = NULL, x = 5, y = 30, label = "plot mpg vs. wt")



pdf("Compare_predictions_ribbons_3_datasets.pdf", width = 13, height = 15)
multiplot(NZoo.ribbons, Gueasso.ribbons, Kokota.ribbons, Cu.Lola.ribbons,
          Foum.ribbons, Gama.ribbons, Laine.ribbons, cols=2)
dev.off()
