peaklists = local({
  filenames = list.files(pattern = '.txt')
  peaklists = lapply(filenames, read.table)
  names(peaklists) = sub('.txt', '', filenames)
  lapply(peaklists, normalize.intensity)
})

peaklist.names = local({
  split = strsplit(names(peaklists), ' ')
  genus = sapply(split, function(s) s[1])
  species = sapply(split, function(s) s[2])
  strain = sapply(split, function(s) do.call(paste, as.list(s[-1:-2])))
  cbind(genus = genus, species = species, strain = strain)
})

local({
  genus.species = paste(peaklist.names[, 1], peaklist.names[, 2])
  count = table(unlist(genus.species))
  species = names(count)[count > 1]
  indexs = genus.species %in% species
  peaklists <<- peaklists[indexs]
  peaklist.names <<- peaklist.names[indexs, ]
})


scores = lapply(1:length(peaklists), function(i) {
  result = search.datebase.bootstrap(
    peaklists[[i]],
    peaklists[-i],
    tolerance = 2000
  )
  
  res = lapply(result, function(res) {
    list(
      res$sample,
      genus = get.bootstrap.score(
        res$bootstrap, 
        function(x) sapply(strsplit(x, ' '), function(x) paste(x[1]))
      ),
      species = get.bootstrap.score(
        res$bootstrap, 
        function(x) sapply(strsplit(x, ' '), function(x) paste(x[1], x[2]))
      )
      # genus = get.bootstrap.confidence.score(
      #   res, 
      #   function(x) sapply(strsplit(x, ' '), function(x) paste(x[1]))
      # ),
      # species = get.bootstrap.confidence.score(
      #   res, 
      #   function(x) sapply(strsplit(x, ' '), function(x) paste(x[1], x[2]))
      # )
    )
  })
  
  print(i)
  res
})

results = lapply(c('cosine', 'eu', 'ieu'), function(method) {
  sample = names(peaklists)
  
  original = local({
    match = sapply(scores, function(x) x[[method]][[1]][1, 1])
    score = sapply(scores, function(x) x[[method]][[1]][1, 2])
  
    genus = sapply(strsplit(sample, ' '), function(x) x[1]) == 
      sapply(strsplit(match, ' '), function(x) x[1])
    species = sapply(strsplit(sample, ' '), function(x) paste(x[1], x[2])) == 
      sapply(strsplit(match, ' '), function(x) paste(x[1], x[2]))
    
    cbind(
      sample = sample,
      match = match,
      score = score,
      genus = genus,
      species = species
    )
  })

  bootstrap.species = local({
    match = sapply(scores, function(x) names(x[[method]]$species)[1])
    score = sapply(scores, function(x) x[[method]]$species[1])
  
    species = sapply(strsplit(sample, ' '), function(x) paste(x[1], x[2])) == 
      sapply(strsplit(match, ' '), function(x) paste(x[1], x[2]))
  
    cbind(
      sample = sample,
      match = match,
      score = score,
      species = species
    )
  })

  bootstrap.genus = local({
    match = sapply(scores, function(x) names(x[[method]]$genus)[1])
    score = sapply(scores, function(x) x[[method]]$genus[1])
    
    genus = sapply(strsplit(sample, ' '), function(x) x[1]) == 
      sapply(strsplit(match, ' '), function(x) x[1])
    
    cbind(
      sample = sample,
      match = match,
      score = score,
      genus = genus
    )
  })
  
  list(original = original,
       bootstrap.species = bootstrap.species,
       bootstrap.genus = bootstrap.genus)
})


write.csv(results[[1]]$original, 'result_cosine.csv', row.names = FALSE)
write.csv(results[[1]]$bootstrap.species, 'result_cosine_bootstrap_species.csv', row.names = FALSE)
write.csv(results[[1]]$bootstrap.genus, 'result_cosine_bootstrap_genus.csv', row.names = FALSE)

write.csv(results[[2]]$original, 'result_eu.csv', row.names = FALSE)
write.csv(results[[2]]$bootstrap.species, 'result_eu_bootstrap_species.csv', row.names = FALSE)
write.csv(results[[2]]$bootstrap.genus, 'result_eu_bootstrap_genus.csv', row.names = FALSE)

write.csv(results[[3]]$original, 'result_ieu.csv', row.names = FALSE)
write.csv(results[[3]]$bootstrap.species, 'result_ieu_bootstrap_species.csv', row.names = FALSE)
write.csv(results[[3]]$bootstrap.genus, 'result_ieu_bootstrap_genus.csv', row.names = FALSE)


pROC::plot.roc(
  results[[1]]$original[, 'species'], 
  predictor = as.numeric(results[[1]]$original[, 'score']), 
  print.auc = T, print.auc.y = 0.4, col = 'tomato3'
)
pROC::plot.roc(
  results[[2]]$original[, 'species'], 
  predictor = as.numeric(results[[2]]$original[, 'score']), 
  print.auc = T, print.auc.y = 0.5, col = 'navy', add = T
)
pROC::plot.roc(
  results[[3]]$original[, 'species'], 
  predictor = as.numeric(results[[3]]$original[, 'score']), 
  print.auc = T, print.auc.y = 0.6, col = 'springgreen3', add = T
)

pROC::plot.roc(
  results[[1]]$original[, 'genus'], 
  predictor = as.numeric(results[[1]]$original[, 'score']), 
  print.auc = T, print.auc.y = 0.4, col = 'tomato3'
)
pROC::plot.roc(
  results[[2]]$original[, 'genus'], 
  predictor = as.numeric(results[[2]]$original[, 'score']), 
  print.auc = T, print.auc.y = 0.5, col = 'navy', add = T
)
pROC::plot.roc(
  results[[3]]$original[, 'genus'], 
  predictor = as.numeric(results[[3]]$original[, 'score']), 
  print.auc = T, print.auc.y = 0.6, col = 'springgreen3', add = T
)

pROC::plot.roc(
  results[[1]]$bootstrap.species[, 'species'], 
  predictor = as.numeric(results[[1]]$bootstrap.species[, 'score']), 
  print.auc = T, print.auc.y = 0.4, col = 'tomato3'
)
pROC::plot.roc(
  results[[2]]$bootstrap.species[, 'species'], 
  predictor = as.numeric(results[[2]]$bootstrap.species[, 'score']), 
  print.auc = T, print.auc.y = 0.5, col = 'navy', add = T
)
pROC::plot.roc(
  results[[3]]$bootstrap.species[, 'species'], 
  predictor = as.numeric(results[[3]]$bootstrap.species[, 'score']), 
  print.auc = T, print.auc.y = 0.6, col = 'springgreen3', add = T
)

pROC::plot.roc(
  results[[1]]$bootstrap.genus[, 'genus'], 
  predictor = as.numeric(results[[1]]$bootstrap.genus[, 'score']), 
  print.auc = T, print.auc.y = 0.4, col = 'tomato3'
)
pROC::plot.roc(
  results[[2]]$bootstrap.genus[, 'genus'], 
  predictor = as.numeric(results[[2]]$bootstrap.genus[, 'score']), 
  print.auc = T, print.auc.y = 0.5, col = 'navy', add = T
)
pROC::plot.roc(
  results[[3]]$bootstrap.genus[, 'genus'], 
  predictor = as.numeric(results[[3]]$bootstrap.genus[, 'score']), 
  print.auc = T, print.auc.y = 0.6, col = 'springgreen3', add = T
)


plot.sensitivity.errorrate = function(data, upper = T, ...) {
  x = if (upper)
    do.call(rbind, lapply(seq(0, 1, 0.01), function(t) {
      acc = sum(data[, 2] >= t & data[, 1]) / sum(data[, 2] >= t)
      s = sum(data[, 2] >= t & data[, 1]) / sum(data[, 1])
      c(t, acc, s)
    }))
  else
    do.call(rbind, lapply(seq(0, 1, 0.01), function(t) {
      acc = sum(data[, 2] <= t & data[, 1]) / sum(data[, 2] <= t)
      s = sum(data[, 2] <= t & data[, 1]) / sum(data[, 1])
      c(t, acc, s)
    }))
  plot(c(0, 0), xlim = c(0, 1), ylim = c(0, 1), ...)
  lines(x = x[, 1], y = x[, 3], lwd = 2, col = 'springgreen4')
  lines(x = x[, 1], y = 1 - x[, 2], lwd = 2, col = 'tomato3')
  cbind(threshold = x[, 1],
        sensitivity = x[, 3],
        error.rate = 1 - x[, 2])
}

plot.sensitivity.errorrate(cbind(
  as.logical(results[[1]]$original[, 'species']), 
  predictor = as.numeric(results[[1]]$original[, 'score'])
), xlab = 'Score', ylab = 'Value')
plot.sensitivity.errorrate(cbind(
  as.logical(results[[1]]$bootstrap.species[, 'species']), 
  predictor = as.numeric(results[[1]]$bootstrap.species[, 'score'])
), xlab = 'Score', ylab = 'Value')

plot.sensitivity.errorrate(cbind(
  as.logical(results[[2]]$original[, 'species']), 
  predictor = as.numeric(results[[2]]$original[, 'score'])
), xlab = 'Score', ylab = 'Value')
plot.sensitivity.errorrate(cbind(
  as.logical(results[[2]]$bootstrap.species[, 'species']), 
  predictor = as.numeric(results[[2]]$bootstrap.species[, 'score'])
), xlab = 'Score', ylab = 'Value')

plot.sensitivity.errorrate(cbind(
  as.logical(results[[3]]$original[, 'species']), 
  predictor = as.numeric(results[[3]]$original[, 'score'])
), xlab = 'Score', ylab = 'Value')
plot.sensitivity.errorrate(cbind(
  as.logical(results[[3]]$bootstrap.species[, 'species']), 
  predictor = as.numeric(results[[3]]$bootstrap.species[, 'score'])
), xlab = 'Score', ylab = 'Value')


plot.sensitivity.errorrate(cbind(
  as.logical(results[[1]]$original[, 'genus']), 
  predictor = as.numeric(results[[1]]$original[, 'score'])
), xlab = 'Score', ylab = 'Value')
plot.sensitivity.errorrate(cbind(
  as.logical(results[[1]]$bootstrap.genus[, 'genus']), 
  predictor = as.numeric(results[[1]]$bootstrap.genus[, 'score'])
), xlab = 'Score', ylab = 'Value')

plot.sensitivity.errorrate(cbind(
  as.logical(results[[2]]$original[, 'genus']), 
  predictor = as.numeric(results[[2]]$original[, 'score'])
), xlab = 'Score', ylab = 'Value')
plot.sensitivity.errorrate(cbind(
  as.logical(results[[2]]$bootstrap.genus[, 'genus']), 
  predictor = as.numeric(results[[2]]$bootstrap.genus[, 'score'])
), xlab = 'Score', ylab = 'Value')

plot.sensitivity.errorrate(cbind(
  as.logical(results[[3]]$original[, 'genus']), 
  predictor = as.numeric(results[[3]]$original[, 'score'])
), xlab = 'Score', ylab = 'Value')
plot.sensitivity.errorrate(cbind(
  as.logical(results[[3]]$bootstrap.genus[, 'genus']), 
  predictor = as.numeric(results[[3]]$bootstrap.genus[, 'score'])
), xlab = 'Score', ylab = 'Value')

