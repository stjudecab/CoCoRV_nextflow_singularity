estimateAFJoint = function(controlID, controlAC, controlAN,
                           caseID, caseAC, caseAN) {
  # AN AC for missing variants
  controlANForNoVariant = round(mean(controlAN))
  caseANForNoVariant = round(mean(caseAN))
  
  controlAF = numeric(length(controlID))
  caseAF = numeric(length(caseID))
  
  # variants both in case and control
  commonVariant = intersect(controlID, caseID)
  if (length(commonVariant) > 0) {
    commonCaseIndex = match(commonVariant, caseID)
    commonControlIndex = match(commonVariant, controlID)
    commonAF = (caseAC[commonCaseIndex] + controlAC[commonControlIndex]) / 
      (caseAN[commonCaseIndex] + controlAN[commonControlIndex])
    controlAF[commonControlIndex] = commonAF
    caseAF[commonCaseIndex] = commonAF
  }
  
  # variants case only
  caseOnlyVariant = setdiff(caseID, controlID)
  if (length(caseOnlyVariant) > 0) {
    caseOnlyIndex = match(caseOnlyVariant, caseID)
    caseOnlyAF = (caseAC[caseOnlyIndex]) / 
      (caseAN[caseOnlyIndex] + controlANForNoVariant)
    caseAF[caseOnlyIndex] = caseOnlyAF
  }
  
  # variants control only
  controlOnlyVariant = setdiff(controlID, caseID)
  if (length(controlOnlyVariant) > 0) {
    controlOnlyIndex = match(controlOnlyVariant, controlID)
    controlOnlyAF = (controlAC[controlOnlyIndex]) / 
      (controlAN[controlOnlyIndex] + caseANForNoVariant)
    controlAF[controlOnlyIndex] = controlOnlyAF
  }
  
  return(list(controlAF = controlAF, caseAF = caseAF))
}