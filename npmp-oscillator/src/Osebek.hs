{-# LANGUAGE TypeSynonymInstances,FlexibleInstances #-}
module Osebek where

import System.Random


data Izr = Gen    Double [(Int,Double)] [(Int,Double)] -- Gen alpha activators repressors
         | LinMod Double Int                           -- LinMod beta p1
         | EnzMod Double Double Int                    -- EnzMod beta km p1
         | ProtProt Int Int Double Double Int Int      -- ProtProt m n kon koff p1 p2
    deriving Show

instance Random Izr where
    randomR (LinMod _ lo,LinMod _ hi) g =
        let (tip,g') = randomR (0,1::Double) g
        in  case tip of
                x|x<0.7 -> (Gen 1.0 [] [],g')
                x|x<0.8 ->
                    let (idx1,g'') = randomR (lo,hi) g'
                    in  (LinMod 1.0 idx1,g'')
                x|x<0.9 ->
                    let (idx1,g'') = randomR (lo,hi) g'
                    in  (EnzMod 1.0 1.0 idx1,g'')
                _ ->
                    let (idx1,g'') = randomR (lo,hi) g'
                        (idx2,g''') = randomR (lo,hi) g''
                    in  (ProtProt 1 1 1.0 1.0 idx1 idx2, g''')
    random = randomR (LinMod undefined 0,LinMod undefined 9)
        
data Deg = Lin Double        -- Lin delta
         | Enz Double Double -- Enz delta km
         | Act Double Int    -- Act delta p2
    deriving Show

instance Random Deg where
    randomR (Act _ lo,Act _ hi) g =
        let (tip,g') = randomR (0,1::Double) g
        in  case tip of
                x|x<0.8 -> (Lin 1.0,g')
                x|x<0.9 -> (Enz 1.0 1.0,g')
                _ ->
                    let (idx2,g'') = randomR (lo,hi) g'
                    in  (Act 1.0 idx2,g'')
    random = randomR (Act undefined 0,Act undefined 9)

type Protein = (Izr,Deg)

instance Random Protein where
    randomR ((LinMod _ lo,_),(LinMod _ hi,_)) g =
        let (izr,g') = randomR (LinMod undefined lo,LinMod undefined hi) g
            (deg,g'') = randomR (Act undefined lo,Act undefined hi) g
        in  ((izr,deg),g'')
    random = randomR ((LinMod undefined 0,undefined),(LinMod undefined 9,undefined))

type Osebek = [Protein]

instance Random Osebek where
    randomR (lo,hi) g =
        iterate f ([],g') !! n
      where
        (n,g') = randomR (length lo,length hi) g
        loidx = (LinMod undefined 0, undefined :: Deg) 
        hiidx = (LinMod undefined (n-1), undefined :: Deg)
        f (osebek,g) =
            let (protein,g') = randomR (loidx,hiidx) g
            in  ((protein:osebek),g')
    random = randomR ([undefined],[undefined | _ <- [0..9]])
        

--newtype Sprememba = Sprememba {unsprememba :: RandomGen g => Int -> Protein -> (Protein,g)}
{-
instance Random Sprememba where
    randomR (_,_) = random
    random g =
        let (c,g') = randomR (0,5::Int) g
        in  case c of
                0 ->
                    let f n p =
                        let (p',g') = randSpremeniIzr g p
                        in 
            
  -}      


----------------------------------------
-- Funkcije za spreminjanje proteinov --
----------------------------------------

spremeniIzr :: Izr -> Protein -> Protein
spremeniIzr izr (_,deg) = (izr,deg)

randSpremeniIzr :: RandomGen g => g -> Int -> Int -> Protein -> (Protein , g)
randSpremeniIzr g n idx p =
    let (izr,g') = randomR (LinMod 0 0,LinMod 0 (n-1)) g
    in (spremeniIzr izr p,g')

spremeniDeg :: Deg -> Protein -> Protein
spremeniDeg deg (izr,_) = (izr,deg)

randSpremeniDeg :: RandomGen g => g -> Int -> Int -> Protein -> (Protein , g)
randSpremeniDeg g n idx p =
    let (deg,g') = randomR (Act 0 0,Act 0 (n-1)) g
    in  (spremeniDeg deg p,g')

dodajAktivator :: (Int,Double) -> Protein -> Protein
dodajAktivator a (izr,deg) =
    case izr of
        Gen alpha as rs | length as <= 2 ->
            (Gen alpha (a:as) rs,deg)
        _  ->
            (izr,deg)

randDodajAktivator :: RandomGen g => g -> Int -> Int -> Protein -> (Protein , g)
randDodajAktivator g n idx p =
    let (a,g') = randomR (0,n-1) g
        --a' = if a >= idx then a+1 else a
    in  (dodajAktivator (a,1.0) p,g')

dodajRepresor :: (Int,Double) -> Protein -> Protein
dodajRepresor r (izr,deg) =
    case izr of
        Gen alpha as rs | length rs <= 2 ->
            (Gen alpha as (r:rs),deg)
        _  ->
            (izr,deg)

randDodajRepresor :: RandomGen g => g -> Int -> Int -> Protein -> (Protein , g)
randDodajRepresor g n idx p =
    let (r,g') = randomR (0,n-1) g
        --r' = if r >= idx then r+1 else r
    in  (dodajRepresor (r,1.0) p,g')

odstraniAktivator :: Int -> Protein -> Protein
odstraniAktivator idxa (izr,deg) =
    case izr of
        Gen alpha as rs ->
            (Gen alpha (delete idxa as) rs,deg)
        _ ->
            (izr,deg)

randOdstraniAktivator :: RandomGen g => g -> Protein -> (Protein , g)
randOdstraniAktivator g p =
    let (idxa,g') = randomR (0,1) g
    in  (odstraniAktivator idxa p,g')

odstraniRepresor :: Int -> Protein -> Protein
odstraniRepresor idxr (izr,deg) =
    case izr of
        Gen alpha as rs ->
            (Gen alpha as (delete idxr rs),deg)
        _ ->
            (izr,deg)

randOdstraniRepresor :: RandomGen g => g -> Protein -> (Protein , g)
randOdstraniRepresor g p =
    let (idxr,g') = randomR (0,1) g
    in  (odstraniRepresor idxr p,g')


--------------------------------------
-- Funkcije za spreminjanje osebkov --
--------------------------------------

dodajProtein :: Protein -> Osebek -> Osebek
dodajProtein protein osebek =
    osebek ++ [protein]

randDodajProtein :: RandomGen g => g -> Osebek -> (Osebek,g)
randDodajProtein g osebek =
    let n = length osebek
        (izr,g') = randomR (LinMod undefined 0,LinMod undefined (n-1)) g
        (deg,g'') = randomR (Act undefined 0, Act undefined (n-1)) g'
    in  if n < 10
            then (dodajProtein (izr,deg) osebek,g'')
            else (osebek,g)

spremeniProtein :: (Protein -> Protein) -> Int -> Osebek -> Osebek
spremeniProtein f idx osebek =
    let (l1,p:l2) = splitAt idx osebek
    in  l1 ++ f p : l2

-- TODO: tole lahko spremenim v tip: RandomGen g => g -> (Osebek -> Osebek , g),
--       kar bi bilo mogoce bols
--       to bi lahko naredil tako, da za idx izberem randomR (0,maxBound) g in potem
--       kot rezultat vrnem funkcijo, ki na tem pravilno naredi `mod`.
--       Isto velja tudi za randSpremeniKoef in randOdstraniProtein
randSpremeniProtein :: RandomGen g => g -> Osebek -> (Osebek,g)
randSpremeniProtein g osebek =
    let n = length osebek
        (idx,g') = randomR (0,n-1) g
        (l1,p:l2) = splitAt idx osebek
        (c,g'') = randomR (0,1::Double) g 
        (p',g''') =
            case c of
                x|x< -0.10 -> randSpremeniIzr g'' n idx p
                x|x< -0.20 -> randSpremeniDeg g'' n idx p
                x|x< 0.25 -> randDodajAktivator g'' n idx p
                x|x< 0.50 -> randDodajRepresor g'' n idx p
                x|x< 0.75 -> randOdstraniAktivator g'' p
                _         -> randOdstraniRepresor g'' p
    in  --(spremeniProtein f idx osebek,g''')
        (l1++p':l2,g''')

odstraniProtein :: Int -> Osebek -> Osebek
odstraniProtein idx osebek = 
    map (\x -> popravi x) $ delete idx osebek
  where
    premakni brisaniIdx idx =
        if idx > brisaniIdx
            then idx - 1
            else idx
    popravi (i,d) =
        (popraviIzr i, popraviDeg d)
    popraviIzr i = case i of
        Gen alpha as rs ->
            Gen alpha
                (map (\(a,kd) -> (premakni idx a,kd)) $ filter ((idx /=) . fst) as)
                (map (\(r,kd) -> (premakni idx r,kd)) $ filter ((idx /=) . fst) rs)
        LinMod beta p1 ->
            if p1 /= idx
                then LinMod beta $ premakni idx p1
                else Gen 1.0 [] []
        EnzMod beta km p1 ->
            if p1 /= idx
                then EnzMod beta km $ premakni idx p1
                else Gen 1.0 [] []
        ProtProt m n kon koff p1 p2 ->
            if p1 /= idx && p2 /= idx
                then ProtProt m n kon koff (premakni idx p1) (premakni idx p2)
                else Gen 1.0 [] []
    popraviDeg d = case d of
        Act delta p2 ->
            if p2 /= idx
                then Act delta $ premakni idx p2
                else Lin 1.0
        _ -> d

randOdstraniProtein :: RandomGen g => g -> Osebek -> (Osebek,g)
randOdstraniProtein g osebek =
    let n = length osebek
        (idx,g') = randomR (0,n-1) g
    in  if n == 1
            then (osebek,g)
            else (odstraniProtein idx osebek, g')
    

----------------------
-- Pomozne funkcije --
----------------------

delete n xs =
    case splitAt n xs of
        (ys,z:zs) -> ys ++ zs
        _         -> xs

