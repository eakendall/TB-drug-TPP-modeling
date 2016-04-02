q <- c(0.025,0.25,0.5,0.75,0.975)
quantile(alldrDST$tbdeaths, q)
quantile(alldrDST$inc, q)

quantile(alldrout$tbdeaths, q)
quantile(alldrout$inc, q)
quantile(alldrout$rrdeaths, q)
quantile(alldrout$rrinc, q)

quantile(1- novelwide1$tbdeaths10alloptimal/alldrDST$tbdeaths10, q)
quantile(1- novelwide1$inc10alloptimal/alldrDST$inc10, q)

quantile(1- novelwide3$rrdeaths10alloptimal/alldrout$rrdeaths10, q)
quantile(1- novelwide3$tbdeaths10alloptimal/alldrout$tbdeaths10, q)
quantile(1- novelwide3$rrinc10alloptimal/alldrout$rrinc10, q)

quantile(1- novelwide1$tbdeaths10allintermediate/alldrDST$tbdeaths10, q)

quantile(1- novelwide3$rrdeaths10allintermediate/alldrout$rrdeaths10, q)

dsmin["efficacy",c(1,3,5)]
dsmax["efficacy",c(1,3,5)]
