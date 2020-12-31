# This is the script to do batch correction for Chris' lungNet data Author: Yuan
# Tian

library("sva")

champ.SlimCombat <- function(beta, batchname="Slide", variablename="Sample_Name", pd, noMod = FALSE) {
    beta <- logit2(beta)
    batch <- pd[, batchname]
    
    formdf <- as.formula(paste(" ~", variablename))
    print(formdf)
    
    mod <- model.matrix(formdf, data = pd)
    message("Generate mod success. Started to run ComBat, which is quite slow...")
    
    if (noMod) {
        combat <- ComBat(dat = beta, batch = batch, mod = NULL, par.prior = TRUE)
    } else {
        combat <- ComBat(dat = beta, batch = batch, mod = mod, par.prior = TRUE)
    }
    
    combat = ilogit2(combat)
    message("Combat Success")
    return(combat)
}
