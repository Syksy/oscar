.onAttach <- function(libname, pkgname){
	packageStartupMessage(paste("If you use oscar-methodology, please consider citing the following paper:",
	"Halkola AS, Joki K, Mirtti T, M\U{E4}kel\U{E4} MM, Aittokallio T, Laajala TD (2023).", 
	"OSCAR: Optimal subset cardinality regression using the L0-pseudonorm with applications to prognostic modelling of prostate cancer.", 
	"PLoS Comput Biol 19(3): e1010333.", 
	"https://doi.org/10.1371/journal.pcbi.1010333", sep="\n")
	)
}
