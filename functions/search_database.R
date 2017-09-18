search.datebase = function(sample.peaklist, reference.spectra, tolerance, 
                           method = c('eu', 'ieu', 'cosine')) {
  
  scores = do.call(rbind, lapply(reference.spectra, function(peaklist) {
    try(
      score(sample.peaklist, peaklist, tolerance, method = method),
      c(0, 0, 0)
    )
  }))
  
  result = lapply(method, function(m) {
    ref.order = order(scores[, m], decreasing = TRUE)
    s = cbind(
      name = names(reference.spectra)[ref.order],
      score = scores[, m][ref.order]
    )
    rownames(s) = NULL
    s
  })
  names(result) = method
  
  result
}