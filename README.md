MethodGames
===========

A short discussion of how one might evaluate a new method for discovering patterns. 

This is the text and code for Bowers, Jake. forthcoming. "Method Games,"
*Sociological Methodology*. 

The results of the machine method game are in ```socmeth.Rout```.

To rerun this game, or to modify it for your own use, you will need to install
the [R statistical computing environment](http://www.r-project.org).

If you are using OS X or some other version of unix, you can rerun the
game in ```socmeth.R``` to recreate ```socmeth.Rout``` by typing the
following at the command line:

```bash
R CMD BATCH socmeth.R
```

Be prepared to wait a while. The code is designed to return NA and to keep
running if the problem becomes intractable for a given method. This means one
should not worry about errors --- the proportion of NAs will be reported below
the report of the proportion of successes. In this version, we see, for
example, that 5/800=.00625 of the versions of the easy data problem (5
variables, all involved in the truth) were intractable by the ```SIS()```
function.

If you are using Windows, you will need to edit the code to avoid the
```parallel``` library and to substitute ```lapply()``` or ```replicate()```
because, as far as I know, the shared memory multicore parallel capabilities
of R have not yet been ported to Windows. You can also setup your own local or remote cluster
using distributed memory parallelism using the ```snow``` package and the ```parLapply()``` function.
