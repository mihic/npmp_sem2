Build instructions
------------------

cabal sandbox init
cabal install --dependencies-only
cabal configure
cabal build

Running
-------

cabal run -- +RTS -N4 -H2096m

