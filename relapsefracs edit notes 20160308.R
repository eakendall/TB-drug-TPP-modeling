Tperiod.lengths <- c(1,1,1,1,2,3,3,8)/12
regimens <- c("s","r","n")
durations <- c(6, 20, 4)/12; names(durations) <- regimens 
Rnames <- c("R0", "Rc", "Rn", "Rcn", "Rr", "Rrc", "Rrn", "Rrcn")



a <- array(999, dim=c(6,5)); rownames(a)<- paste0("duration",1:6); colnames(a) <- c(paste0("finished",0:3,".5"),"finished5")
for (p in 1:5) for (d in 1:6) 
{
  durations[3] <- d/12
  a[d,p] <- relapsefracs(p)[1,3]
}
a  

Old version: 
  finished0.5 finished1.5 finished2.5 finished3.5 finished4.5 finished5.5
duration2           1  0.06090294  0.04060196  0.04060196  0.04060196  0.04060196
duration3           1  0.15225735  0.05413595  0.04060196  0.04060196  0.04060196
duration4           1  0.17509595  0.12941875  0.05075245  0.04060196  0.04060196
duration5           1  1.00000000  0.15225735  0.06496314  0.04060196  0.04060196
duration6           1  1.00000000  0.16748308  0.13703161  0.05413595  0.04060196

New version: 
  finished0.5 finished1.5 finished2.5 finished3.5 finished4.5 finished5.5
duration2   0.3653348  0.05125992  0.02050397  0.02050397  0.02050397  0.02050397
duration3   0.5768899  0.10764584  0.04100794  0.02050397  0.02050397  0.02050397
duration4   0.6826674  0.14224628  0.07304539  0.03588195  0.02050397  0.02050397
duration5   0.7461339  0.23840179  0.10764584  0.05741111  0.02050397  0.02050397
duration6   0.7884449  0.36533483  0.13071280  0.08457887  0.04100794  0.02050397

more accurately:
  finished0.5 finished1.5 finished2.5 finished3.5 finished5
duration1  0.09806163  0.01903855  0.01903855  0.01903855 0.01903855
duration2  0.34746851  0.05438356  0.01903855  0.01903855 0.01903855
duration3  0.56497901  0.09806163  0.04260189  0.01903855 0.01903855
duration4  0.67373426  0.12198392  0.07413933  0.03671106 0.01903855
duration5  0.73898741  0.21696222  0.09806163  0.06145257 0.01903855
duration6  0.78248950  0.34746851  0.11400982  0.08211343 0.04260189


2 vs 4 month regimen: relapses are reduced by ltfu*(0.08-0.03 + 0.05-0.02 + 0.023-0.019) = 0.08*ltfu ~ 0.24% of total patients
4 vs 6: ltfu*(1-0.08 + 0.075-0.052 + 0.06-0.02 + 2*(0.024-0.02)) ~ 1*ltfu ~ 3% of total patients, so theres actually a lot more relapse among drug suscetibles as a result of the longer duration for 6 vs 4 than for 4 vs 2 under these approximations and assumptions. And no one with any resistance is getting treated with the novel regimen, so their outcomes aren't changing with novel reigmen durations. Therefore, we see a big impact from reduced ltfu and from better outcomes in especially tfu's during month 2. 


(For comparison, the relapse that occurs regardless of ltfu is 2% here, and the relapse due to ltfu for even a 2mo regimen is ltfu*(1+0.03) ~ 3%. )

Effects that can make a longer novel regimen worse: 
1) relapses that do occur, occur sooner.
2) cured individuals are susceptible to reinfection sooner.

But a 6 mo regimen should have ~1.4% more relapses than a 4mo regimen (3% * (0.1+0.22+0.06+0.05+0.02+0.02))
compared to a (3% + (0.32+0.09+0.05+0.02)) ~ also ~1.4% difference for a 4 vs 2 mo regimen,

               
               
