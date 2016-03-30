#!/usr/bin/env Rscript
library(RSQLite)
library(ggplot2)
library(reshape2)
library(plyr)
library(gridExtra)
library(scales)
library(grid)

resultsDir = '../../analysis/sensitivity_analysis/'
plotDir = './'
resultsDb = paste(resultsDir,'results.sqlite',sep='')

parNames = c('meanStep','sdStep','muPhenotype')

reps = 10
comboDb = dbConnect(SQLite(), dbname = resultsDb)
initExtension(comboDb)
textSize = 6
plot_themes  = 	theme_classic() +
				theme(axis.line = element_line(size=.3)) +
				theme(axis.ticks = element_line(size=.3)) +
				theme(axis.ticks.length = unit(-.07,'cm')) +
				theme(axis.ticks.margin = unit(.11,'cm')) +
				theme(axis.title.x=element_text(size=textSize+1)) + 
				theme(axis.text.x=element_text(size= textSize)) + 
				theme(axis.title.y=element_text(size= textSize+1)) +
				theme(axis.text.y=element_text(size= textSize)) +
				theme(plot.title=element_text(size=textSize+2)) +
				theme(plot.margin=unit(c(1,1,1,1),'mm')) +
				theme(legend.title=element_text(size=textSize)) +
				theme(legend.text=element_text(size=textSize)) +
				theme(legend.position ='none') +
				theme(panel.border = element_rect(colour = "black", fill=NA, size=.5)) +
				theme(axis.line = element_blank()) 
								
blank_axes = theme(axis.title.x=element_blank()) + 
			theme(axis.text.x=element_blank()) + 
			theme(axis.title.y=element_blank()) + 
			theme(axis.text.y=element_blank()) 
			
xlabels = 	theme(axis.title.x=element_text(size=textSize+1, margin=margin(t = -2, b=10))) + 
			theme(axis.text.x=element_text(size= textSize,margin=margin(t = 2, b=12)))
antigenDF = dbGetQuery(comboDb, paste('SELECT',paste(parNames,collapse=','),', runId, trunkProportions, antigenicLagAG1, meanFluxRate, meanAnnualIncidence, extinct, excessDiversity, fluLike FROM pooled_results',sep=' '))

summaryDF = ddply(antigenDF, .(meanStep, sdStep, muPhenotype), summarise, tropicsAgLag = mean(antigenicLagAg1,na.rm=TRUE), tropicsTrunkPro = mean(trunkProportions, na.rm=TRUE), meanAnnualIncidence = mean(meanAnnualIncidence,na.rm=TRUE), meanFluxRate = mean(meanFluxRate,na.rm=TRUE), extinct=sum(extinct),excessDiversity= sum(excessDiversity),fluLike=sum(fluLike),total=sum(extinct,excessDiversity,fluLike))

summaryDF$meanFluxRate[is.na(summaryDF$meanFluxRate)] = -1
summaryDF$meanAnnualIncidence[is.na(summaryDF$meanAnnualIncidence)] = -1


parName1 = 'meanStep'
parName2 = 'sdStep'

MUPHENOTYPES = c(5e-5, 1e-4, 2e-4)
FLUPLOTS = c()
plots = c()

for(plotMu in MUPHENOTYPES){

	plotDF = subset(summaryDF, muPhenotype == plotMu)
	
	fluStatPlot = ggplot(plotDF, aes_string(x=parName1,y=parName2)) + 
		geom_tile(aes(fill=fluStat)) + 
		scale_fill_gradient2("Flu like",limits=c(0,reps),low='#2A4A7F',high="#800000",midpoint=reps/2) +
		xlab(expression(beta['mean'])) + 
		ylab(parName2) + 
		plot_themes +
		ggtitle(paste('mu = ',plotMu)) + #,'\nFlu-like: R0 = ',plotR0)) +
			scale_x_continuous(expand = c(0,0)) +
			scale_y_continuous(expand = c(0,0)) 
	
	fluLikePlot = ggplot(plotDF, aes_string(x=parName1,y=parName2)) + 
		geom_tile(aes(fill=fluLike)) + 
		scale_fill_gradient2("Flu like counts",limits=c(0,reps),low='#2A4A7F',high="#800000",midpoint=reps/2) +
		ylab(substitute(delta[sd])) + 
		plot_themes +
		theme(axis.title.x=element_blank()) + 
		theme(axis.text.x=element_blank()) + 
		ggtitle(bquote(mu ==.(plotMu))) + theme(plot.title=element_text(hjust=-.3)) +#,'\nFlu-like',sep=''))) +
			scale_x_continuous(expand = c(0,0)) +
			scale_y_continuous(expand = c(0,0)) 
			
	extinctPlot = ggplot(plotDF, aes_string(x=parName1,y=parName2)) + 
		geom_tile(aes(fill=extinct)) + 
		scale_fill_gradient2("Extinction counts",limits=c(0,reps),low='#2A4A7F',high="#800000",midpoint=reps/2) +
		#xlab(substitute(delta[mean])) + 
		#ylab(substitute(delta[sd])) + 
 			plot_themes +
 			blank_axes +
		scale_x_continuous(expand = c(0,0)) +
		scale_y_continuous(expand = c(0,0)) 

	diversityPlot = ggplot(plotDF, aes_string(x=parName1,y=parName2)) + 
		geom_tile(aes(fill=excessDiversity)) + 
		scale_fill_gradient2("Excess diversity counts",limits=c(0,reps),low='#2A4A7F',high="#800000",midpoint=reps/2) +
		plot_themes +
		blank_axes +
		scale_x_continuous(expand = c(0,0)) +
		scale_y_continuous(expand = c(0,0)) 

	meanFluxPlot = ggplot(plotDF, aes_string(x=parName1,y=parName2)) + 
		geom_tile(aes(fill=as.numeric(meanFluxRate))) + 
		scale_fill_gradientn('Mean antigenic flux',
							colours=c('#2A4A7F','white',"#800000"),
							values = rescale(c(0,1.0,3.0)),
							guide = 'colorbar', limits = c(0,3.0)) +	
		plot_themes +
		blank_axes +
		scale_x_continuous(expand = c(0,0)) +
		scale_y_continuous(expand = c(0,0)) 
		
	incidencePlot = ggplot(plotDF, aes_string(x=parName1,y=parName2)) + 
		geom_tile(aes(fill=as.numeric(meanAnnualIncidence))) + 
		scale_fill_gradientn('Incidence',
							colours=c('#2A4A7F','white',"#800000"),
							values = rescale(c(0,.12,.35)),
							guide = 'colorbar', limits = c(0,.35)) +	
		plot_themes +
		blank_axes +
		scale_x_continuous(expand = c(0,0)) +
		scale_y_continuous(expand = c(0,0)) 
		 
	if(plotMu == MUPHENOTYPES[length(MUPHENOTYPES)]){
		legend_themes = theme(legend.title=element_text(size=textSize)) +
			theme(legend.text=element_text(size=textSize)) +
			theme(legend.position ='bottom') +
			theme(legend.direction='horizontal') +
			theme(legend.key.height=unit(.2,'cm')) +
			theme(plot.margin=unit(c(1,1,1,1),'mm')) +
			theme(legend.margin = unit(0,'mm')) 
		fluLikePlot = fluLikePlot + legend_themes +	guides(fill=guide_colorbar(barwidth=4,barheight=0.2,title.position='bottom')) + xlabels + xlab(substitute(delta[mean]))
		extinctPlot = extinctPlot + legend_themes +guides(fill=guide_colorbar(barwidth=4,barheight=0.2,title.position='bottom')) + xlabels + xlab(substitute(delta[mean]))
		diversityPlot = diversityPlot + legend_themes +guides(fill=guide_colorbar(barwidth=4,barheight=0.2,title.position='bottom')) + xlabels + xlab(substitute(delta[mean]))
		meanFluxPlot = meanFluxPlot+ legend_themes +guides(fill=guide_colorbar(barwidth=4,barheight=0.2,title.position='bottom')) + xlabels + xlab(substitute(delta[mean]))
		incidencePlot = incidencePlot + legend_themes +guides(fill=guide_colorbar(barwidth=4,barheight=0.2,title.position='bottom')) + xlabels + xlab(substitute(delta[mean]))

	}
	thisRow = cbind(ggplotGrob(fluLikePlot), ggplotGrob(extinctPlot), ggplotGrob(diversityPlot), ggplotGrob(meanFluxPlot), ggplotGrob(incidencePlot),size = 'first')
	print(plotMu)
	if(plotMu == MUPHENOTYPES[1]){
		plots = thisRow
	}
	else{
		plots =  rbind(plots,thisRow,size='first')
	}

}


pdf(paste(plotDir,'sensitivity_analysis.pdf',sep='/'),7,5)
grid.draw(plots)
dev.off()

dbDisconnect(comboDb)