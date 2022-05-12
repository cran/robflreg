MMcontrol <-
function(bdp=0.50, eff=0.95, shapeEff=FALSE, convTol.MM=1e-7, maxIt.MM=50, fastScontrols=Scontrol(...), ...) {
  
  return(list(eff=eff, bdp=bdp, shapeEff=shapeEff, convTol.MM=convTol.MM, maxIt.MM=maxIt.MM, fastScontrols=fastScontrols))
}
