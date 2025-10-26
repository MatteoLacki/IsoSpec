#' An example of how to add your own elements.
#' 
#' This can be used, for instance, with isotopically labelled molecules.
#'
#' @export
custom_isotopes_example = function(){
  print(isotopicData)

  print('Here are some isotopes')
  i = isotopicData$IsoSpecShortZero
  print(i)
  
  # convention: D for deuterium, M for N15
  foney_elements = data.frame(
  	element = c('D', 'M'),
  	isotope = c('D', 'M'),
  	mass    = c(i[i$isotope == 'H2', 'mass'], i[i$isotope == 'N15', 'mass']),
  	abundance = c(1, 1),
  	ratioC  = c(NA, NA) # this is not important actually, I just keep it for consistency
  )	
  
  isotopes = rbind(i, foney_elements)
  print(isotopes)
  
  # Your atom counts:
  # C  C13  H   N   N15 O
  # 37  6   71  9   4   13
  # translate into:
  your_molecule = c(C=37, H=71, N=9, O=13, D=6, M=4)
  
  res = IsoSpecify(your_molecule, .999, isotopes, showCounts=T)
  print(res)
  
  return(invisible())
}