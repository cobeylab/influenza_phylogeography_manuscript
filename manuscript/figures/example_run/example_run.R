library(RSQLite)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(png)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
######################

figureWidth = 5 #inches
aspectRatio = 0.8

resultsDir = '../../analysis/default_params/'
plotDir = './'
exampleDir = paste(resultsDir,'6',sep='')		
tips = read.csv(paste(exampleDir,'/out.tips',sep=''))
branches = read.csv(paste(exampleDir,'/out.branches',sep=''),sep=',',header=FALSE)
timeseries = read.csv(paste(exampleDir,'/out.timeseries',sep=''),sep='\t')
prevalence = cbind(timeseries$date,data.frame(timeseries$northI/timeseries$northN, timeseries$tropicsI/timeseries$tropicsN, timeseries$southI/timeseries$southN)*10^5)
names(prevalence) = c('year','north','tropics','south')

textSize = 9
plot_themes  = 	theme_classic() +
				theme(axis.line = element_line(size=.2)) +
				theme(axis.ticks = element_line(size=.2)) +
				theme(axis.ticks.length = unit(.05,'cm')) +
				theme(axis.title.x=element_text(size=textSize)) + 
				theme(axis.text.x=element_text(size= textSize)) + 
				theme(axis.title.y=element_text(size= textSize)) +
				theme(axis.text.y=element_text(size= textSize)) +
				theme(plot.title=element_text(size=textSize)) +
				theme(plot.margin=unit(c(.1,0,-1.0,0),'cm')) +
				theme(legend.title=element_text(size=textSize)) +
				theme(legend.text=element_text(size=textSize)) +
				theme(legend.position ='bottom') +
				theme(legend.direction='horizontal') +
				theme(legend.margin = unit(0,'cm')) +
				theme(text = element_text(family='serif')) 

zoom_themes = theme(legend.position="none", axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),axis.title.y=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_rect(color='black', fill="white"),
            plot.margin = unit(c(0,0,0,0),"mm"))

branches = within(branches, V9<-data.frame(do.call('rbind',strsplit(as.character(branches$V9),'}\t{',fixed=TRUE))))
branches = within(branches, V17<-data.frame(do.call('rbind',strsplit(as.character(branches$V17),'}\t',fixed=TRUE))))
virus = data.frame(branches[,1],branches[,2],branches[,3],branches[,4],branches[,5],branches[,6],branches[,7],branches[,8],branches[,9][1])
parents = data.frame(branches[,9][2],branches[,10],branches[,11],branches[,12],branches[,13],branches[,14],branches[,15],branches[,16],branches[,17][1])
names(virus) = names(tips)
names(parents) = names(tips)
branches =  data.frame(virus,parents)
names(branches)=c(names(tips),paste('parent',names(tips),sep=''))

tips$distance = (tips$ag1^2+tips$ag2^2)^.5

northTips = tips[which(tips$loc==0),]
tropicsTips = tips[which(tips$loc==1),]
southTips = tips[which(tips$loc==2),]

fit = loess(distance~year,data=tips,span=0.2,degree=1)

j = order(tips$year)
regionColors = c(rgb(0.765, 0.728, 0.274), rgb(0.324, 0.609, 0.708), rgb(0.857, 0.131, 0.132))
plotDF = data.frame(tips$year,tips$distance,tips$loc)
names(plotDF) = c('year','distance','loc')
ydXlim = c(20,23); ydYlim = c(18,22)

yeardistancePlotFull = ggplot(data = plotDF, aes(x=year, y=distance)) +
				geom_point(colour = regionColors[plotDF$loc +1], size = 0.1) +
				geom_smooth(method='loess',colour = 'black',size=0.2,se=FALSE,span=0.25) +
				ylab('Antigenic distance from founder') +
				xlab('year') +
				plot_themes +
				ggtitle(expression(paste("(",italic(b),")"))) +  
				theme(plot.title=element_text(hjust=-.25)) +
				geom_rect(aes(xmin = ydXlim[1], xmax= ydXlim[2],ymin= ydYlim[1], ymax= ydYlim[2]),color='black',fill='transparent', size = .2) +
				theme(plot.margin=unit(c(.1,0,-1.0,.5),'cm')) 


yeardistancePlotZoom = ggplot(data = plotDF, aes(x=year, y=distance)) +
				geom_point(colour = regionColors[plotDF$loc +1], size = 0.1) +
				geom_smooth(method='loess',colour = 'black',size=0.5,se=FALSE,span=0.25) +
				coord_cartesian(xlim= ydXlim, ylim= ydYlim) + zoom_themes
				
yeardistancePlot = yeardistancePlotFull + annotation_custom(grob=ggplotGrob(yeardistancePlotZoom),xmin=0,xmax=16,ymin=18.6,ymax=40)					


treeXlim = c(20,23); treeYlim = c(1600,1750)
treePlotFull = ggplot(data = branches) + 
			geom_segment(data = branches, mapping = aes(x = year, y = layout, xend = parentyear,yend = layout),colour = regionColors[branches$loc+1], size = .1) +
			geom_segment(data = branches, mapping = aes(x = parentyear, y = layout, xend = parentyear,yend = parentlayout), colour = regionColors[branches$loc+1], size= .1) +
			plot_themes +
			ylab(NULL) +
			xlab('year') +
			theme(axis.line.y = element_blank()) +
			theme(axis.ticks.y = element_blank()) +
			theme(axis.text.y = element_blank()) +
			ggtitle(expression(paste("(",italic(a),")"))) +  
			theme(plot.title=element_text(hjust=0)) +
			geom_rect(aes(xmin = treeXlim[1], xmax=treeXlim[2],ymin=treeYlim[1], ymax=treeYlim[2]),color='black',fill='transparent', size = .3)

trunkBranches = branches[branches$trunk==1,]
sideBranches = branches[branches$trunk==0,]
treePlotZoom = ggplot(data = trunkBranches) + 
			geom_segment(data = trunkBranches, mapping = aes(x = year, y = layout, xend = parentyear,yend = layout),colour = regionColors[trunkBranches$loc+1], size = .6,lineend='round') +
			geom_segment(data = trunkBranches, mapping = aes(x = parentyear, y = layout, xend = parentyear,yend = parentlayout), colour = regionColors[trunkBranches$loc+1], size= .6,lineend='round') +
			geom_segment(data = sideBranches, mapping = aes(x = year, y = layout, xend = parentyear,yend = layout),colour = regionColors[sideBranches $loc+1], size = .1, alpha=0.5) +
			geom_segment(data = sideBranches, mapping = aes(x = parentyear, y = layout, xend = parentyear,yend = parentlayout), colour = regionColors[sideBranches$loc+1], size= .1,alpha=0.5) +
			coord_cartesian(xlim=treeXlim, ylim=treeYlim) + zoom_themes
			
treePlot = treePlotFull + annotation_custom(grob=ggplotGrob(treePlotZoom),xmin=0,xmax=20.5,ymin=2175,ymax=3200)					
	
prevalencePlot = ggplot(data = prevalence) +
					geom_line(mapping = aes(x=year, y=north), colour = regionColors[1], size = 0.3) + 
					geom_line(mapping = aes(x=year, y=tropics), colour = regionColors[2], size =0.3) + 
					geom_line(mapping = aes(x=year, y=south), colour = regionColors[3], size =0.3) +
					plot_themes +
					xlab('year') +
					ylab(expression('prevalence per 10'^5)) +
					theme(aspect.ratio=0.25) +
					theme(plot.margin=unit(c(0,2,-1,0),'cm')) +
					ggtitle(expression(paste("(",italic(c),")"))) +  
					theme(plot.title=element_text(hjust=-.15,margin=margin(t = -1, b=12))) 

pdf(paste(plotDir,'example_run.pdf',sep=''),width = figureWidth,height = figureWidth*aspectRatio)
multiplot(treePlot, yeardistancePlot, prevalencePlot, layout = matrix(c(1,2,3,3), nrow=2, byrow=TRUE))
dev.off()