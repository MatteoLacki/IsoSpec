# add this file to the tests

library(enviPat)
library(stringr)
library(IsoSpecR)
library(tidyverse)
library(microbenchmark)
library(jsonlite)

data(isotopicData)
# store args in function env for multiple calls without the passing of args
envipat_call_factory = function(...) function() isopattern(...)
is.error = function(x) inherits(x, "try-error")

# test the call for results and their timing
test_envipat = function(..., times = 100, timing = T){
  envipat_call = envipat_call_factory(...)
  out = list()
  call_result = try(envipat_call())
  if(is.error(call_result)){
    return(call_result)
  } else{
    call_result = call_result[[1]]
    if(timing) out$timing = microbenchmark(envipat_call(), times = times, unit = 'us')
    call_result = tbl_df(call_result)
    call_result = call_result[,1:2]
    colnames(call_result) = c('mass', 'prob')
    call_result = arrange(tbl_df(call_result), desc(prob))
    out$call_result = call_result
    return(out)
  }
}


run_test = function(molecule, threshold, out_path, timing=T){
  res = test_envipat(isotopes = isotopicData$IsoSpec, 
                     chemforms = molecule, 
                     threshold = threshold, 
                     verbose = F, 
                     rel_to = 2,
                     timing = timing)
  if(is.error(res)){
    return(res)
  } else {
    file_name = file.path(out_path, paste0(molecule, "_", threshold, ".tsv"))
    readr::write_tsv(res$call_result, file_name, col_names = F)
    return(res$timing)
  }
}

out_path  = "/Users/matteo/Projects/isospec/IsoSpec/tests/envipat_results"
thresholds = seq(10^{-5}, .05, by=.01)

data(chemforms)
atoms = c("C", "H", "N", "O", "S")
envipat_mols = chemforms[!str_detect(chemforms, "\\[")]
envipat_path= file.path(out_path, "envipat_mols")
mist_mols = c("P1", "P2", "H1", "H2", "O1", "O2", "H2O1", "C0", "P0", "C100O0P100", "C100", "P100", "C1", "H10C10O10N10S5", "Se1", "Se10", "Sn1", "Sn4", "Sn4C1", "C2H6O1", "C1000", "C1H1O2N2Se1Sn1P1", "P1C1Sn1", "Se5", "Sn5", "Se2Sn2C2O2N2S2B2He2U2Na2Cl2")
mist_path = file.path(out_path, "mist_mols")
human = fromJSON("/Users/matteo/Projects/furious_fastas/py_data/human_fastas_small.json")
human_path= file.path(out_path, "human_uniprot_mols")
human = human[,atoms]
human = sapply(
  split(human, 1:nrow(human)),
  function(x){
    a = atoms[x > 0]
    x = x[x > 0]
    paste(a, x, sep='', collapse='')
  },
  USE.NAMES = F)

timeout = R.utils::withTimeout
run_no_timing = function(..., secs) try(timeout(run_test(...),
                                                timeout = secs),
                                        silent = T)

run_tests = function(molecules, out_path, thresholds, secs = 10, mc.cores = 5)
  parallel::mclapply(thresholds, 
                     function(thr)  sapply(molecules, run_no_timing,
                                           threshold = thr,
                                           out_path  = out_path,
                                           secs      = secs),
                     mc.cores = 5)

mist_timings = run_tests(mist_mols, mist_path, thresholds)
envipat_timings = run_tests(envipat_mols, envipat_path, thresholds)
human_timings = run_tests(human[1:2000], human_path, thresholds[1])

save(mist_timings,
     envipat_timings,
     human_timings,
     file = paste0(out_path,"timings_envipat.Rda"))
