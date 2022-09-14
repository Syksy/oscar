###
#
# Raw R script used to precompute the examples in the vignette, so vignette compilation doesn't take too long
#
###

library(oscar)
data(ex)

# Cox model fit for example data
tmp <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family="cox")

png("tmp_visu.png", width=700, height=700)
visu(tmp)
dev.off()

# Cross-validation

cv.tmp <- cv.casso(tmp, fold=10, seed=0)

png("tmp_cv_visu.png", width=700, height=700)
cv.visu(cv.tmp)
dev.off()

# Bootstrapping

bs.tmp <- bs.oscar(tmp, bootstrap=10, seed=0)

png("tmp_bs_visu.png", width=700, height=700)
bs.visu(bs.tmp)
dev.off()

# kmax

tmp.kmax10<- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, kmax=10, family="cox")

png("tmp_kmax10_visu.png")
visu(tmp.kmax10)
dev.off()

# Save it as a convenient .RData to be loaded during vignette construction
save.image("vignette_workspace.RData")

