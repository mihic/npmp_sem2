module SpremeniKoef (randSpremeniKoef) where

import Osebek

import System.Random


countParamsDeg i = case i of
                   Lin _ -> 1
                   Enz _ _ -> 2
                   Act _ _ -> 1

countParamsIzr d = case d of
                   LinMod _ _ -> 1
                   EnzMod _ _ _  -> 2
                   ProtProt _ _ _ _ _ _ -> 2
                   Gen _ akts reps -> 1 + length akts + length reps

countParamsSingle (i,d) = countParamsIzr i + countParamsDeg d
countParams ps = foldl (\acc x -> acc + countParamsSingle x) 0 ps


spremeniKoef :: Int -> Double -> Osebek -> Osebek
spremeniKoef idx c (p:ps) = 
    let 
      cnt = idx - countParamsSingle p
    in
      if cnt < 0 then
        (changeParameterSingle idx c p):ps
      else p:(spremeniKoef cnt c ps)

randSpremeniKoef :: RandomGen g => g -> Osebek -> (Osebek,g)
randSpremeniKoef g osebek =
    let n = countParams osebek
        (idx,g') = randomR (0,n-1) g
        (c,g'') = randomR (0.0,2.0) g'
    in  (spremeniKoef idx c osebek,g'')


changeParameterSingle idx c (i,d) =
    let lenI = countParamsIzr i
    in
        if lenI > idx then
            (changeParameterIzr idx c i,d)
        else
            (i,changeParameterDeg (idx-lenI) c d)

changeParameterIzr idx c i = 
    case i of
        LinMod a b-> LinMod (a*c) b
        EnzMod a b n -> if idx==0 then EnzMod (a*c) b n else EnzMod a (b*c) n
        ProtProt x1 x2 y1 y2 z1 z2 ->
            if idx==0 then ProtProt x1 x2 (y1*c) y2 z1 z2
            else ProtProt x1 x2 y1 (y2*c) z1 z2
        Gen a akts reps -> 
            if idx==0 then Gen (a*c) akts reps
            else
                if (idx-1)<(length akts) then 
                    let tmp = if (idx-1)==0 then 
                            let (x,y) = head akts
                            in (x,c*y) : (tail akts)
                        else 
                            let (x,y) = akts !! 1
                            in [head akts,(x,c*y)]
                    in
                        Gen a tmp reps
                else 
                    let tmp = if (idx-1-length akts)==0 then
                            let (x,y) = head reps
                            in (x,c*y) : (tail reps)
                        else 
                            let (x,y) = reps !! 1
                            in [head reps,(x,c*y)]
                    in Gen a akts tmp


changeParameterDeg idx c d =
    case d of
        Lin a -> Lin $ a*c
        Enz a b -> if idx==0 then Enz (a*c) b else Enz a (b*c)
        Act a b -> Act (a*c) b
