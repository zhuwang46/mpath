\name{be.zeroinfl}
\alias{be.zeroinfl}

\title{conduct backward stepwise variable elimination for zero inflated count regression}

\description{
 conduct backward stepwise variable elimination for zero inflated count regression from zeroinfl function}

\usage{
be.zeroinfl(object, data, dist=c("poisson", "negbin", "geometric"), alpha=0.05, 
            trace=FALSE)
}

\arguments{
  \item{object}{an object from function zeroinfl}
  \item{data}{argument controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{dist}{one of the distributions in zeroinfl function; defaults to the distribution of \code{object}.}
  \item{alpha}{significance level of term elimination}
  \item{trace}{logical value, if TRUE, print detailed calculation results}
}
\details{
 Conduct backward stepwise term elimination for zero-inflated count regression.
 Terms in the count and zero components are compared using Wald tests. For a
 term with multiple coefficients, such as a factor, a joint Wald test is used
 and the whole term is retained or removed together. Main effects required by
 interaction terms are retained. Elimination stops when no eligible term has
 a test p-value greater than \code{alpha}, or when a reduced model cannot be
 fitted.}
\value{An object of class \code{zeroinfl} for the selected model. Terms whose
 Wald tests cannot be computed are retained.}
\references{
Zhu Wang, Shuangge Ma, Ching-Yun Wang, Michael Zappitelli, Prasad Devarajan and Chirag R. Parikh (2014)
\emph{EM for Regularized Zero Inflated Regression Models with Applications to Postoperative Morbidity after Cardiac Surgery in Children}, \emph{Statistics in Medicine}. 33(29):5192-208.

Zhu Wang, Shuangge Ma and Ching-Yun Wang (2015) \emph{Variable selection for zero-inflated and overdispersed data with application to health care demand in Germany}, \emph{Biometrical Journal}. 57(5):867-84.
}

\author{Zhu Wang <zwang145@uthsc.edu>}

\keyword{models}
\keyword{regression}

 
