
# pso = P(stay | optimal)
# pss = P(stay | suboptimal)
pso <- .895
pss <- .195




di <- function(pso, pss){
  c <- pso - pss
  d <- (pso+pss-1)/(1-pso+pss)
  c <- ifelse(c<0,0,c)
  d <- ifelse(d<0,0,d)
  cbind(c=c,  d=d)
}
di(pso, pss)


oo <-.8 # seq(0,1,.02)
ss <- seq(0,1,.01)
g <- expand.grid(oo,ss)
plot(NA, xlim=0:1,ylim=0:1, xlab="c", ylab="d",asp=1)
for(i in 1:length(cc))
  lines(di(oo[i],ss),type="l", col = 1) #paste0("gray", i))

