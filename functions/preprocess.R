normalize.intensity = function(data, max.intensity = 100) {
  normalized.intensity = data[, 2] / max(data[, 2]) * max.intensity
  normalized.intensity[normalized.intensity < 0] = 0
  new.data = data
  new.data[, 2] = normalized.intensity
  new.data
}

crop.mz = function(data, mz.lower = 4000, mz.upper = 12000) {
  data[data[, 1] >= mz.lower & data[, 1] <= mz.upper, ]
}