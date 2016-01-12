module Genetski where

import Debug.Trace

import Osebek
import SpremeniKoef
import ToODE

import Numeric.GSL
import Numeric.LinearAlgebra

import Data.List
import Data.Ord
import System.Random

import Control.Parallel
import Control.Parallel.Strategies

simuliraj :: Vector Double -> Osebek -> Matrix Double
simuliraj ts osebek = ("\n"++(show osebek)++"\n") `trace` odeSolve (toODE osebek) (100 : replicate (length osebek - 1) 0.5) ts

simuliraj_slabo :: Double -> Int -> Osebek -> [Double]
simuliraj_slabo delta num osebek = 
    let ode = toODE osebek 0 --cas v toODE je brezveze
        zacetni_pogoji = (100 : replicate (length osebek - 1) (0.5::Double))
        rezultati = scanl (\pog _ -> zipWith (\a b -> a + b*delta) pog $ ode pog
                            ) zacetni_pogoji [0..num-1]
    in map head rezultati

data Smer = Gor Int | Dol Int | Naravnost
stej_ekstreme :: Double -> [Double] -> Int
stej_ekstreme delta xs= res
  where
    (res,_,_) = foldl' stej (0,Dol 0,head xs) xs

    stej :: (Int,Smer,Double) -> Double -> (Int,Smer,Double)
    stej (cnt,p_smer,prev) x =
        (n_cnt,smer,x)
      where
        
        smer =
            case ((x-prev)/delta,p_smer) of
                (y, Dol 100)  | y>1 -> Gor 0
                (y, Gor n)    | y>1 -> Gor (min 100 (n+1))
                (y, Dol n)    | y>1 -> Dol (n+1)

                (y,Gor (-100)) | y<1 -> Dol 0
                (y,Dol n)      | y<1 -> Dol (max (-100) (n-1))
                (y,Gor n)      | y<1 -> Dol (n-1)


                _       -> p_smer
        n_cnt =
            case (p_smer,smer) of
                (Gor _,Dol _) -> cnt + 1
                (Dol _,Gor _) -> cnt + 1
                _         -> cnt
oceni1 :: Vector Double -> (Double -> Double) -> Osebek -> Double -- (Double -> Double) je odvec
oceni1 ts _ osebek =
    let lts = toList ts
        delta = lts !! 1 - lts !! 0 --to je slabo, lahko bi bil argument
        num = length lts
        p0s = simuliraj_slabo delta num osebek
    in ((fromIntegral $ stej_ekstreme delta p0s) - 5.0)**2



oceni :: Vector Double -> (Double -> Double) -> Osebek -> Double
oceni ts f osebek =
    let ys0 = cmap f ts
        ys = head $ toColumns $ simuliraj ts osebek
        razlike = ys - ys0
    in  razlike <.> razlike

oceni_slabo :: Vector Double -> (Double -> Double) -> Osebek -> Double
oceni_slabo ts f osebek = 
    let lts = toList ts 
        delta = lts !! 1 - lts !! 0 --to je slabo, lahko bi bil argument
        num = length lts
        ys0 = map f $ lts
        p0s = simuliraj_slabo delta num osebek
        razlike = zipWith (-) ys0 p0s
        kvadrati = zipWith (*) razlike razlike
--    in  trace ((show osebek)++"\n") $ sum kvadrati
    in  sum kvadrati

ustvariPopulacijo :: RandomGen g => g -> Int -> ([Osebek],g)
ustvariPopulacijo g n =
    iterate (\(osebki,g) -> let (osebek,g') = random g in (osebek:osebki,g')) ([],g) !! n

-- mutacija v prvi fazi
mutirajPopulacijo :: RandomGen g => g -> [Osebek] -> ([Osebek],g)
mutirajPopulacijo g [] = ([],g)
mutirajPopulacijo g (osebek:osebki) =
    let (osebki',g') = mutirajPopulacijo g osebki
        (c,g'') = randomR (0,1::Double) g'
        (osebek',g''') =
            case c of
                x | x < -1.09 -> randDodajProtein g'' osebek
                x | x < -1.19 -> randOdstraniProtein g'' osebek
                x | x < 0.3 -> randSpremeniProtein g'' osebek
                _            -> randSpremeniKoef g'' osebek
    in  (osebek':osebek:osebki',g''')


-- izboljsava v prvi fazi
izboljsajPopulacijo :: (RandomGen g) => g -> (Osebek -> Double) -> [Osebek] -> ([Osebek],g)
izboljsajPopulacijo g oceni populacija =
    let (mutPop,g') = mutirajPopulacijo g populacija
        graded = zip (parMap rdeepseq oceni mutPop) mutPop
        sortedtuples = sortBy (comparing fst) graded
        sortPop = map snd sortedtuples
        --sortPop = sortOn oceni mutPop
        newPop = take (length populacija) sortPop
    in  (show $ oceni $ head newPop) `trace` (newPop,g')

{-
mutirajPopulacijo :: RandomGen g => g -> [Osebek] -> ([Osebek],g)
mutirajPopulacijo g [] = ([],g)
mutirajPopulacijo g (osebek:osebki) =
    let (osebki',g') = mutirajPopulacijo g osebki
        (c,g'') = randomR (0,1::Double) g'
        (osebek',g''') =
            case c of
                x | x < 0.05 -> randDodajProtein g'' osebek
                x | x < 0.10 -> randOdstraniProtein g'' osebek
                x | x < 0.40 -> randSpremeniProtein g'' osebek
                _            -> randSpremeniKoef g'' osebek
    in  (osebek:osebek':osebki',g''')


-- izboljsava v drugi fazi
izboljsajPopulacijo :: (RandomGen g) => g -> (Osebek -> Double) -> [Osebek] -> ([Osebek],g)
izboljsajPopulacijo g oceni populacija =
    let (mutPop,g') = mutirajPopulacijo g populacija
        graded = zip (parMap rdeepseq oceni mutPop) mutPop
        sortedtuples = sortBy (comparing fst) graded
        sortPop = map snd sortedtuples
        --sortPop = sortOn oceni mutPop
        newPop = take (length populacija) sortPop
    in  (newPop,g')
-}
