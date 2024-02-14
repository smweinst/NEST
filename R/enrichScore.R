# function to compute enrichment score
# input arguments:
## stat.map
## net.map - vector (with same dimension as stat.map) with 1s indicating locations inside network and 0s indicating locations outside network\
## one.sided
## save.detail

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

# OLD VERSION:
# enrichScore = function(L, L.network.labels, p, networks, save_vertex_level=FALSE){
#   # L = unsorted statistics at each vertex/location
#   # L.network.labels = labels of networks corresponding to each location
#   # p = exponent in enrichment score calculation (p=1 if score should be weighted by statistic; p = 0 if they get equal weight)
#   # netEnrich = network labels that we want enrichment score for
#   
#   # order locations by statistic
#   L.order = order(L, decreasing = TRUE) # added decreasing = TRUE 02/22/2023
#   L.sorted = L[L.order]
#   L.sorted.p = abs(L.sorted)^p # absolute value and exponentiate (potential contribution of each vertex to enrichment score if vertex belongs to network of interest)
#   L.network.labels.sorted = L.network.labels[L.order]
#   #L.with.labels.sorted = cbind(L = L.sorted.p,net = L.network.labels[L.order]) 
#   V = length(L) # number of vertices/locations
#   
#   ES = vector(mode="numeric", length = length(networks))
#   P.hit = P.miss = list()
#   out = list()
#   par(mfrow=c(2,1))
#   for (curr.net in networks){ # loop over the networks
#     
#     curr.net.ind = ifelse(L.network.labels.sorted==curr.net, TRUE, FALSE)
#     curr.net.ind[which(is.na(curr.net.ind))] = FALSE # add this in in case some network labels are NA (01/17/2024)
#     P.hit.numerator = cumsum(L.sorted.p*curr.net.ind)#*curr.net.ind
#     P.hit[[paste0("net",curr.net)]] = P.hit.numerator/P.hit.numerator[V] # denominator = final cumulative sum
#     
#     P.miss.numerator = cumsum(1-curr.net.ind)#*(1-curr.net.ind)
#     P.miss[[paste0("net",curr.net)]] = P.miss.numerator/P.miss.numerator[V] # denominator = final cumulative sum (=number of misses)
#     
#     running_sum = P.hit[[paste0("net",curr.net)]] - P.miss[[paste0("net",curr.net)]]
#     
#     #ES.v = (P.hit[[paste0("net",curr.net)]]-P.miss[[paste0("net",curr.net)]])
#     ES[curr.net] = max(abs(running_sum))#*sign(running_sum[which.max(abs(running_sum))])
#     
#     if (save_vertex_level==TRUE){
#       out[[paste0("net",curr.net)]] = list(P.hit = P.hit[[paste0("net",curr.net)]],
#                                            P.miss = P.miss[[paste0("net",curr.net)]],
#                                            net.ind = curr.net.ind,
#                                            ES = ES[curr.net],
#                                            running_sum = running_sum)
#     } else{
#       out[[paste0("net",curr.net)]] = ES[curr.net]
#     }
#     
#     # ES.v = (P.hit[[paste0("net",curr.net)]]-P.miss[[paste0("net",curr.net)]])
#     # ES[curr.net] = max(abs(ES.v))
#     # ES.ind = which.max(abs(ES.v))
#     # plot(ES.v, main = paste0("network ", curr.net), ylim=c(-0.45,0.45), type = 'l', lwd=3,xlab = "rank", ylab = "P.hit-P.miss")
#     # abline(h=0, lty=3)
#     # abline(v=ES.ind, col = "cornflowerblue",lwd=2)
#     # 
#   }
#   
#   return(out)
#   
# }
