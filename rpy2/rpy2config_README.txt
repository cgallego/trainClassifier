rpy2config_README for running rpy2 from Winpython:

1) Set up enviroment variables as:
R_HOME=C:\Program Files\R\R-2.15.3
R_LIBS_USER=C:\Users\windows\R\win-library\2.15


2) within spyder python interpreter:

from rpy2.robjects.packages import importr
utils = importr("utils")
caret = importr("caret")
graphics = importr('graphics')
grdevices = importr('grDevices')
base = importr('base')
stats = importr('stats')

import array

x = array.array('i', range(10))
y = stats.rnorm(10)

grdevices.X11()

graphics.par(mfrow = array.array('i', [2,2]))
graphics.plot(x, y, ylab = "foo/bar", col = "red")

kwargs = {'ylab':"foo/bar", 'type':"b", 'col':"blue", 'log':"x"}
graphics.plot(x, y, **kwargs)


m = base.matrix(stats.rnorm(100), ncol=5)
pca = stats.princomp(m)
graphics.plot(pca, main="Eigen values")
stats.biplot(pca, main="biplot")


3) to port code from R:

Finding where an R symbol is coming from
====================
import rpy2.robjects.packages as rpacks
env = rpacks.wherefrom('expand.grid')
env.do_slot('name')[0]

import rpy2.robjects as robjects
x = robjects.r.seq(1, 5)
tuple(x)

x.names = robjects.StrVector('abcde')
print(x)

x[0]

==========
x = robjects.ListVector({'a': 1, 'b': 2, 'c': 3})
x.rx2[{'i': x.names.index('b')}] = 9

4) to import code from R and call it from python:
====================
http://rpy.sourceforge.net/rpy2/doc-dev/html/robjects_rpackages.html


====================
library(devtools)
source_url('https://raw.github.com/hadley/stringr/master/R/c.r')
source('Z:/Cristina/MassNonmass/Section1 - ExperimentsUpToDate/RMarkdown/classifierComp/example.r')
 
or downloading content
====================
import urllib2
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

bioc_url = urllib2.urlopen('https://raw.github.com/hadley/stringr/master/R/c.r')
rcode = ''.join(bioc_url.readlines())

rfun = SignatureTranslatedAnonymousPackage(rcode, "rfun")


or for simpler functions:
====================
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

string = """
square <- function(x) {
    return(x^2)
}

cube <- function(x) {
    return(x^3)
}
"""

powerpack = SignatureTranslatedAnonymousPackage(string, "powerpack")
