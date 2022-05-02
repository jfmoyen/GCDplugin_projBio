\encoding{utf8}
\name{projbioplot}
\alias{projbioplot}
\alias{projbiocoords}
\title{
    Calculates ternary coordinates projected from biotite (and plots the ternary diagram)
}

\description{
    The function \code{projbiocoords} calculates the coordinates used to define the projection, 
    whereas \code{projbioplot} calls \code{projbiocoords} and does the other operations needed 
    to define a template that will be used by \code{\link{plotDiagram}}.
}

\usage{
    projbioplot(mins=FALSE,addWR=FALSE,ticks=FALSE,xmin=-2,xmax=1,ymin=-0.5,ymax=0.5)
    projbiocoords(where=WR,add=FALSE)
}
\arguments{
  \item{mins}{Boolean. If \code{TRUE}, the composition of ideal minerals will be plotted on 
  the diagram, namely "q","fsp","an","an50","cz","Ep","sill","opx","cpx","olv","grs-Gt","Gt","NaCrd","Crd","bio","ms","MgHbl","Edn" and "Pgs". }

  \item{addWR}{Boolean. If \code{TRUE}, the results of the calculations will be added to WR as new columns with colnames = \code{c("ms1","fsp","CaAl","bio")} }
  
    \item{ticks}{Boolean. If \code{TRUE}, tick marks will be added to the side of the diagram (similar to the \code{ticks} option in \code{\link{ternary}} ) }
	
	 \item{xmin, xmax, ymin, ymax}{Numeric. Bounds of the plotting area, in rectangular coordinates. Note that the vertical axis (left side) is at x=0; (ms1 = 3 al + 2 (Na + K) is at x=0 and y=0.5, and CaAl = Ca + Al is at x=0 and y=-0.5); the fsp (=Al + (Na + K) ) is at \eqn{x=\sqrt{3}/2} and y=0. Plagioclase an50 is at \eqn{x=-\sqrt{3}/2} and y=0, and most points should fall to the right of it. The defaults are therefore sensible.}
	 
	   \item{where}{A data matrix containing whole rock analyses, to be projected. By default \code{WR}}

  \item{add}{Boolean. If \code{TRUE}, the results of the calculations will be added to WR as new columns with colnames = \code{c("ms1","fsp","CaAl","bio")} } 
}

\value{
    For \code{projbioplot}, nothing. This function is meant only to be called via \code{\link{plotDiagram}}.

    \code{projbiocoords} returns a matrix of 4 columns containing the new coordinates,  \code{(ms1,fsp,CaAl,bio)}.
}

\details{
    The "projection from biotite" of \emph{Moyen et al. (2016)} is based on four coordinates based on molar proportions:

    \deqn{ms1 = Al - Ca - NK}	
    \deqn{fsp = -2 Al -2 Ca + 3 NK -1/3 FM}
	\deqn{CaAl = Ca}
	\deqn{bio = 1/3 FM}
    
    where
    Al = molar Al

    Ca = molar Ca

    FM = molar Fe + molar Mg

    NK = molar Na + molar K

    molar proportions are calculated by \code{\link{millications}}.

    They are primarily used for plotting the relevant diagram (ms1-fsp-CaAl, note that bio 
    is dropped during plotting as this is projected from bio).

    \code{projbiocoords} calculates the coordinates. It also has the side effect of (globally) 
    affecting its results to \code{results}, where they are available to \code{\link{addResults}}, 
    \code{\link{r2clipboard}(results)}, \code{\link{HTMLTableResults}}, etc.

    \code{projbioplot} is used purely for defining a Figaro template. The values of the four 
    coordinates are calculated (and (globally) assigned to \code{results}). 
    If called with \code{\link{plotDiagram}}, the diagram is plotted, which is the real 
    use of the function.

    If used as a plugin, a GUI function \code{.projbioGUI} is also supplied 
    (and linked to the menu item), in charge of gathering the missing arguments and 
    calling \code{\link{plotDiagram}}.

    \if{html}{\figure{projbio.png}}
    \if{latex}{\figure{projbio.png}{options: width=5.5in}}
}

\seealso{
    \link{plotDiagram} \link{millications}
}

\references{
    Moyen JF, Laurent O, Chelle-Michou C, Couzinie S, Vanderhaeghe O, Zeh A, 
    Villaros A, Gardien V, 2016. Collision vs. subduction-related magmatism: 
    two contrasting ways of granite formation and implications for crustal growth. 
    Lithos in press.
    \href{http://dx.doi.org/10.1016/j.lithos.2016.09.018}{doi:10.1016/j.lithos.2016.09.018} 

    and in particular the supplementary item SE4 "Multivariate statistics and 
    projection for granitic rocks".
}

\author{
    Jean-François Moyen \email{jfmoyen@gmail.com}
}

\examples{
    data(sazava)
    accessVar("sazava")
    projbiocoords()
    plotDiagram("projbioplot")

    plotDiagram("projbioplot",mins=T,xmin=-1)
}

\keyword{hplot}
\concept{Calculation}
\concept{Plot}
\concept{Projection}
\concept{Granites}
\concept{Cations}
\concept{Menu: Plugins}

