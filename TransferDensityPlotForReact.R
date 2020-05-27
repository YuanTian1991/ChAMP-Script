#load("./BMIQ_Result.rda")

getPoints <- function(i)
{
tmp <- density(myNorm[,i])

#points <- list()
#for(i in 1:length(tmp$x))
#{
#    points[[i]] <- list(x=tmp$x[i],y=tmp$y[i])
#}

oneSample <- list(label=colnames(myNorm)[i],
                  fill=FALSE,
                        backgroundColor='rgba(75,192,192,0.4)',
      borderColor='rgba(75,192,192,1)',
      borderCapStyle='butt',
      borderDashOffset=0.0,
      borderJoinStyle='miter',
      pointBorderColor='rgba(75,192,192,1)',
      pointBackgroundColor='#fff',
      pointBorderWidth=1,
      pointHoverRadius=1,
      pointHoverBackgroundColor='rgba(75,192,192,1)',
      pointHoverBorderColor='rgba(220,220,220,1)',
      pointHoverBorderWidth=2,
      pointRadius=1,
      pointHitRadius=2,
      data=tmp$y[seq(1,512,by=8)]
                  )

}

datasets <- list()
for(i in 1:ncol(myNorm))
{
    datasets[[i]] <- getPoints(i)
}

#apply(myNorm,2,function(x) getPoints(x))

library(jsonlite)
#toJSON(points,auto_unbox = FALSE)

sink("myjson.json")
toJSON(datasets,auto_unbox=TRUE)
sink()
