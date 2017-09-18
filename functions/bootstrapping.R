search.datebase.bootstrap = function(
  sample.peaklist, reference.spectra, tolerance, 
  method = c('eu', 'ieu', 'cosine'),
  bootstrap.times = 100
) {
  
  peak.numbers = nrow(sample.peaklist)
  bootstrap.peak.indexs = lapply(1:bootstrap.times, function(i) {
    indexs = sample.int(peak.numbers, peak.numbers, replace = TRUE)
  })
  
  scores = lapply(reference.spectra, function(peaklist) {
    peaks1 = normalize.intensity(sample.peaklist) 
    peaks2 = normalize.intensity(peaklist)
    common.peaks = find.common.peaks(peaks1, peaks2, tolerance)
    
    common.peak.indexs = sapply(1:peak.numbers, function(peak.idx) {
      common.idx = which((common.peaks[[1]][, 1] == common.peaks[[3]][peak.idx, 1]) & 
                           (common.peaks[[1]][, 2] == common.peaks[[3]][peak.idx, 2]))
      if(length(common.idx) == 0)
        common.idx = 0
      common.idx[1]
    })
    
    scores.origin = sapply(method, function(m) {
      do.call(m, args = list(common.peaks = common.peaks))
    })
    names(scores.origin) = method
    
    scores.bootstrap = lapply(bootstrap.peak.indexs, function(bootidx) {
      cp = local({
        p3 = common.peaks[[3]][bootidx, ]
        p4 = common.peaks[[4]]
        cp.idx = common.peak.indexs[bootidx]
        cp.idx = cp.idx[cp.idx > 0]
        if(length(cp.idx) == 1) {
          cp1 = matrix(common.peaks[[1]][cp.idx, ], nrow = 1)
          cp2 = matrix(common.peaks[[2]][cp.idx, ], nrow = 1)
        }
        else {
          cp1 = common.peaks[[1]][cp.idx, ]
          cp2 = common.peaks[[2]][cp.idx, ]
        }
        list(cp1, cp2, p3, p4, tolerance = tolerance)
      })
      scores.bootstrap = sapply(method, function(m) {
        do.call(m, args = list(common.peaks = cp))
      })
      names(scores.bootstrap) = method
      scores.bootstrap
    })
    
    scores = lapply(method, function(m) {
      list(
        score = scores.origin[m],
        scores.bootstrap = sapply(scores.bootstrap, function(s) s[m])
      )
    })
    names(scores) = method
    scores
  })
  
  order.scores = function(sc) {
    ref.order = order(sc, decreasing = TRUE)
    s = cbind(
      name = names(reference.spectra)[ref.order],
      score = sc[ref.order]
    )
    rownames(s) = NULL
    s
  }
  
  result = lapply(method, function(m) {
    scores.origin = sapply(scores, function(s) s[[m]]$score)
    sample = order.scores(scores.origin)

    bootstrap = lapply(1:bootstrap.times, function(i) {
      scores.bootstrap = sapply(scores, function(s) s[[m]]$scores.bootstrap[i])
      order.scores(scores.bootstrap)
    })
    
    list(
      sample = sample,
      bootstrap = bootstrap
    )
  })
  names(result) = method
  
  result
}

get.bootstrap.score = function(bootstrap.result, classFunc = function(x) x) {
  bootstrap.names = sapply(bootstrap.result, function(x) x[1, 'name'])
  bootstrap.classes = classFunc(bootstrap.names)
  scores = table(bootstrap.classes) / length(bootstrap.classes)
  scores[order(scores, decreasing = T)]
}

get.bootstrap.confidence.score = function(search.database.result, classFunc = function(x) x) {
  bootstrap.names = sapply(search.database.result$bootstrap, function(x) x[1, 'name'])
  bootstrap.classes = classFunc(bootstrap.names)
  sum(bootstrap.classes == classFunc(search.database.result$sample[1, 'name'])) / length(bootstrap.classes)
}
