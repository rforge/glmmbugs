\name{mergeBugsData-methods}
\docType{methods}
\alias{mergeBugsData}
\alias{mergeBugsData-methods}
\alias{mergeBugsData,data.frame-method}
\alias{mergeBugsData,SpatialPolygonsDataFrame-method}

\title{Merge results from BUGS into a data.frame or SPDF}
\description{
 merge the result from bugs function
}
\usage{
\S4method{mergeBugsData}{data.frame}(
x, bugsSummary, by.x = NULL, newcol = "mean", ...
)
\S4method{mergeBugsData}{SpatialPolygonsDataFrame}(
x, bugsSummary, by.x = NULL, newcol = "mean", ...
)
}

\arguments{
  \item{x}{spatial polygon object i.e population data set (popdata)}
  \item{bugsSummary}{ posterior distribution result from summaryChain function }
  \item{by.x}{ the common term from the spatial polygon object and the bugs function result  }
  \item{newcol}{the summary statistic that to be merged back to the data frame  }
  \item{\dots}{ additional arguments }
}

\author{ Patrick Brown }


\examples{
if(require('diseasemapping')){
	data('popdata')
	newdata = c("3560102"=2, "3560104"=3)
	popdatatry = mergeBugsData(popdata, newdata, by.x="CSDUID")
}

}
