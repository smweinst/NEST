#' function to calculate enrichment score given a map of brain-phenotype associations & network labels
#'@param stat.map - vector (length V) of location (v)-specific associations
#'@param net.map - vector (also length V) consisting of all 0's, 1's. It is assumed that 1's denote locations inside the network of interest, and 0's denote all other locations
#'@param one.sided - default is TRUE, so that the running sum statistic increases in increments proportional to the magnitude of T(v) at a given location in the network; setting one.sided = FALSE will change the increment for locations in the network to increase by 1/(# locations in network)
#'@param save.detail - default is FALSE, so that the calculated enrichment score is the only information returned. it can be set to FALSE to return additional information used to calculate the enrichment score
#'@export

enrichScore = function(stat.map, net.map, one.sided = TRUE, save.detail = FALSE){

  V = length(stat.map) # number of locations (e.g., vertices)

  # sort stat.map (strong positive to strong negative)
  stat.map.order = order(stat.map, decreasing = TRUE)
  stat.map.sorted = stat.map[stat.map.order]

  net.map.sorted = net.map[stat.map.order]

  # p.hit = running sum increment when there is a 'hit' (location in network)
  if (one.sided == TRUE){ # one-sided test is default -- test whether stat values inside network are 'more extreme' than those outside the network
    p.hit.numerator = cumsum(abs(stat.map.sorted)*net.map.sorted)
  } else{ # two-sided setting -- test whether statistic differs (in distribution) inside vs. outside the network; no distinction between enrichment vs. under-enrichment
    p.hit.numerator = cumsum(net.map.sorted) #cumsum(rep(1, length(stat.map.sorted)))
  }
  p.hit = p.hit.numerator/p.hit.numerator[V]

  # p.miss = running sum increment when there is a 'miss' (location not in network)
  p.miss.numerator = cumsum(1-net.map.sorted)
  p.miss = p.miss.numerator/p.miss.numerator[V]

  running.sum = p.hit - p.miss

  # ES = enrichment score (test statistic for NEST)
  ES = max(abs(running.sum))

  if (save.detail == FALSE){
    return(ES)
  } else{
    return(list(ES = ES,
                p.hit = p.hit,
                p.miss = p.miss,
                running.sum = running.sum))
  }

}
