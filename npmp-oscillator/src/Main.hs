module Main where

import Genetski
import Osebek

import System.Random

import Numeric.LinearAlgebra

import Graphics.Rendering.Chart.Easy
import Graphics.Rendering.Chart.Backend.Cairo


--problem = [(LinMod 1.0 1,Lin 1.0),(Gen 100.0 [] [(0,1.0)],Lin 1.0)]
--
--pf :: Double -> [Double] -> [Double]
--pf t [p0,p1] =
--  [ p1 - p0
--  , (100*heaviside (1-p0)) - 2*p1
--  ]

linearnaDegradacija = [ (Gen 0 [] [], Lin 0.1) ]

encimskaModifikacija = [ (Gen 0 [] [], Lin 0.0)
                       , (EnzMod 2.0 30.0 0, Lin 0.0)
                       ]

represilator3 = [ (Gen 100.0 [] [(2,1.0)], Lin 1.0)
                , (Gen 100.0 [] [(0,1.0)], Lin 1.0)
                , (Gen 100.0 [] [(1,1.0)], Lin 1.0)
                ]

sinus amplituda perioda x = amplituda /2 * (sin ((x*2/perioda*pi)+pi/2) + 1)
hvsinus amplituda perioda x =
    let y = amplituda * (sin ((x*2/perioda*pi)+pi/2))
    in  if y > 0
            then y
            else 0


main = do
    g <- getStdGen
    --let g = mkStdGen 249
    let amplituda = 100
        perioda = 13.5
        ts = linspace 10000 (0::Double,3*perioda)
        f = sinus amplituda perioda
--    let (populacija,g') = ustvariPopulacijo g 50
        populacija = replicate 50 (replicate 3 (Gen 1.0 [] [],Lin 1.0))
        --populacija = replicate 50 represilator3
        iteracija = iterate (\(pop,g) -> izboljsajPopulacijo g (oceni1 ts f) pop) $ (populacija,g)
        bests = dropWhile (\x -> (oceni1 ts f x)>1.0) $ map (head . fst) iteracija
        populacija2 = replicate 50 $ head bests
        iteracija2 = iterate (\(pop,g) -> izboljsajPopulacijo g (oceni_slabo ts f) pop) $ (populacija2,g)
        bests2 = map (head . fst) iteracija2

--        ocene = map (oceni_slabo ts f) bests
    --mapM_ (putStrLn . show) ocene

    let sol = simuliraj ts $ bests2 !! 400
    --let sol = simuliraj ts $ bests !! 0

    toFile def{_fo_format=PS} "najboljsi.ps" $ do
        layout_title .= "najboljsi"
        setColors [opaque blue, opaque red, opaque green, opaque black,
                   opaque orange, opaque cyan, opaque blueviolet,
                   opaque lime ,opaque magenta, opaque chocolate]
        mapM_ plot $
              zipWith (\i l -> line ("P" ++ show i) [l])
                      [0..] $
                      map (zip $ toList ts) $ map toList $ toColumns sol
        plot $ line "kriterij" [zip (toList ts) (toList $ cmap f ts)]
    return ()
