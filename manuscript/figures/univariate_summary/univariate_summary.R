#!/usr/bin/env Rscript


library(RSQLite)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(sensitivity)
library(scales)
library(xtable)
library(gtable)
library(plyr)

resultsDir = '../../analysis/univariate_analyses/'
plotDir = './'
resultsDb = paste(resultsDir,'results.sqlite',sep='')

textSize = 6
pointSize = 0.5
lineSize = .5
plot_themes  = 	theme_classic() +
				theme(axis.line = element_line(size=.3)) +
				theme(axis.ticks = element_line(size=.3)) +
				theme(axis.ticks.length = unit(-.07,'cm')) +
				theme(axis.ticks.margin = unit(.11,'cm')) +
				theme(axis.title.x=element_text(size=textSize)) + 
				theme(axis.text.x=element_text(size= textSize)) + 
				theme(axis.title.y=element_text(size= textSize)) +
				theme(axis.text.y=element_text(size= textSize)) +
				theme(plot.title=element_text(size=textSize)) +
				theme(plot.margin=unit(c(3,3,3,3),'mm')) +
				theme(legend.title=element_text(size=textSize)) +
				theme(legend.text=element_text(size=textSize)) +
				theme(legend.position ='bottom') +
				theme(legend.direction='horizontal') +
				theme(legend.margin = unit(-.5,'cm')) +
				theme(panel.border = element_rect(colour = "black", fill=NA, size=.5)) +
				theme(axis.line = element_blank())

parNames = c('relativeR0','seasonalAmplitude','relativeN','relativeTurnover','tropicFractionI0')
defaults = c(1.0,1.0,1.0,1.0,0.0)
parDefaults = data.frame('relativeN'=1.0,'tropicFractionI0'=1.0,'relativeR0'=1.0,'relativeTurnover'=1.0,'seasonalAmplitude'=0.0)

comboDb = dbConnect(SQLite(), dbname = paste(resultsDir,resultsDb,sep='/'))
initExtension(comboDb)

antigenDF = dbGetQuery(comboDb, paste('SELECT',paste(parNames,collapse=','),', sweepID,tropicsTrunkPro, tropicsAgLag, extinct, excessDiversity, fluLike FROM pooled_results',sep=' '))
headerNames = c('\"Relative R\"[0]','\"Seasonal amplitude\"','\"Relative population size\"','\"Relative turnover\"','\"Fraction of initial infecteds in tropics\"')
headers = data.frame(parNames,headerNames)
fit = lm(tropicsAgLag~relativeN+ tropicFractionI0 +relativeR0+relativeTurnover+seasonalAmplitude,data=antigenDF)

plots=c()

for(i in 1:length(parNames)){
	parName = parNames[i]
	print(parNames[i])
	#exclude non-sweep parameters
	plotDF = antigenDF[which(antigenDF$sweepID==parName),]

	trunkProFit = lm(tropicsTrunkPro~plotDF[,i],data=plotDF)
	pearsonR = cor.test(plotDF[,i],plotDF$tropicsTrunkPro)$estimate
	CI = cor.test(plotDF[,i],plotDF$tropicsTrunkPro)$conf.int
	
	trunkProPlot = ggplot(plotDF, aes_string(x=parName,y="tropicsTrunkPro")) + 
		geom_point(size = pointSize) + 
		xlab(NULL) +
		ylab(NULL) +
		ylim(c(0,1)) + 
		plot_themes +	
		theme(axis.text.x=element_blank()) +
		theme(axis.text.y=element_blank()) +
		stat_smooth(method="lm",se=FALSE,size=lineSize) +
		geom_hline(aes(yintercept=1/3),linetype=2) +
		ggtitle(paste('r = ',round(pearsonR,2),'\n95% CI: (',round(CI[1],2),',',round(CI[2],2),')',sep=''))	
		
	if(i==1){
		trunkProPlot = trunkProPlot +
			ylab("Tropics trunk proportion") +
			plot_themes +
			theme(axis.text.x = element_blank()) 
		trunkProPlotRow = ggplotGrob(trunkProPlot)
	}
	else{
		trunkProPlotRow = cbind(trunkProPlotRow,ggplotGrob(trunkProPlot),size='first')
	}
}	

for(i in 1:length(parNames)){
	parName = parNames[i]
	print(parNames[i])
	plotDF = antigenDF[which(antigenDF$sweepID==parName),]

	pearsonR = cor.test(plotDF[,i],plotDF$tropicsAgLag)$estimate
 	CI = cor.test(plotDF[,i],plotDF$tropicsAgLag)$conf.int

	agLagPlot = ggplot(plotDF, aes_string(x=parName,y="tropicsAgLag")) + 
		geom_point(size = pointSize) + 
		xlab(parse(text=paste(headers$headerNames[headers$parNames==parName]))) + 
		ylab(NULL) + 
		ylim(c(-0.4,0.4)) + 
		plot_themes +	
		theme(axis.text.y=element_blank()) +
		stat_smooth(method="lm",se=FALSE,size=lineSize) +
		geom_hline(aes(yintercept=0.0),linetype=2) +
		ggtitle(paste('r = ',round(pearsonR,2),'\n95% CI: (',round(CI[1],2),',',round(CI[2],2),')',sep=''))	

	if(i==1){
		agLagPlot = agLagPlot +
			ylab("Tropics antigenic lead") +
			plot_themes 
		agLagPlotRow = ggplotGrob(agLagPlot)
	}
	else{
		agLagPlotRow = cbind(agLagPlotRow, ggplotGrob(agLagPlot), size='first')
	}

}

for(parName in parNames){
	plotDF = antigenDF[which(antigenDF$sweepID==parName),]
	plotDF = ddply(plotDF, .(relativeR0, seasonalAmplitude, relativeN, relativeTurnover, tropicFractionI0), numcolwise(sum))

	dfm = melt(plotDF, id.vars=c(parName),measure.vars=c("extinct","excessDiversity","fluLike"), variable.name="behavior",value.name="Count")
	fluCounts = ggplot(dfm,aes_string(x=parName,y="Count")) +
					geom_bar(aes_string(fill="behavior"),stat="identity") + 
					plot_themes + 
					theme(legend.position='none') +
					ylab('') +
					xlab('') +
					theme(axis.text.y=element_blank()) 

	if(parName == parNames[1]){
		fluCounts = fluCounts + 
				ylab('Counts') + 				
				theme(axis.title.y = element_text(size = textSize)) +
				theme(axis.text.y = element_text(size = textSize,hjust = -.4)) 
		statusPlotRow = ggplotGrob(fluCounts)
	}
	else{
		statusPlotRow = cbind(statusPlotRow,ggplotGrob(fluCounts),size='first')
	}
}

figureS1 = rbind(trunkProPlotRow,agLagPlotRow,statusPlotRow,size='first')
pdf(paste(plotDir,'univariate_summary.pdf',sep='/'),7,5)
grid.draw(figureS1)

dev.off()
dbDisconnect(comboDb)
