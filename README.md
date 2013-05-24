MethodGames
===========

A short discussion of how one might evaluate a new method for discovering patterns.

The results of the simulation are in ```socmeth.Rout```.

If you are using OS X or some other version of unix, you can rerun the
simulation in ```socmeth.R``` to recreate ```socmeth.Rout``` by typing the following at the command line:

```
R CMD BATCH socmeth.R
```

Be prepared to wait a while.

If you are using Windows, you will need to edit the code to not use the ```parallel``` library and to substitute ```lapply()``` or ```replicate()``` because, as far as I know, the shared memory multicore parallel capabilities of R have not yet been ported to Windows. You can also setup your own cluster using distributed memory parallelism using the ```snow``` package and ```parLapply()```.
