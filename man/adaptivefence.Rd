\name{adaptivefence}
\alias{adaptivefence}
\title{Adaptive Fence model selection}
\usage{
adaptivefence(mf, f, ms, d, lf, bootset = list(), bs = NULL, grid = 51,
  pickf = md, sizef = size_default)
}
\arguments{
  \item{mf}{function for fitting the model}

  \item{f}{formular of full model}

  \item{ms}{list of formular of candidates models}

  \item{d}{data}

  \item{lf}{function for lack of fit}

  \item{bootset}{bootstrap setting B: number of bootstrap
  sample mode: parametric("p") or nonparametric("np")}

  \item{bs}{bootstrap sample, default setting is
  automatically nonparametric bootstrap (B = 100)}

  \item{grid}{grid for c}

  \item{pickf}{function for choosing the model for fixed c}

  \item{sizef}{function for determinant model size}
}
\value{
list with whatever
}
\description{
Adaptive Fence model selection
}

