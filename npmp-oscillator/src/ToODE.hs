module ToODE (toODE) where

import Osebek


toVecIzr vec idx izr t ps0 =
  let p = ps0 !! idx
  in case izr of
      Gen alpha as rs ->
        vec `update` ( idx
                     ,   alpha
                       * (product $ map (\(idxa,kd) -> heaviside ((ps0!!idxa)-kd)) as)
                       * (product $ map (\(idxr,kd) -> heaviside (kd-(ps0!!idxr))) rs)
                     )
      LinMod beta idx1 ->
        let p1 = ps0 !! idx1
        in  vec
            `update` (idx ,  beta*p1)
            `update` (idx1, -beta*p1)
      EnzMod beta km idx1 ->
        let p1 = ps0 !! idx1
        in  vec
            `update` (idx ,  beta * p1 / (km + p1))
            `update` (idx1, -beta * p1 / (km + p1))
      ProtProt m n kon koff idx1 idx2 ->
        let p1 = ps0 !! idx1
            p2 = ps0 !! idx2
        in vec
           `update` (idx1, -kon
                          * (fromIntegral m) * p1^m
                          * (fromIntegral n) * p2^n
                          + koff * (fromIntegral m) * p)
           `update` (idx2, -kon
                          * (fromIntegral m) * p1^m
                          * (fromIntegral n) * p2^n
                          + koff * (fromIntegral n) * p)
           `update` (idx ,  kon
                          * (fromIntegral m) * p1^m
                          * (fromIntegral n) * p2^n
                          - koff * p)
toVecDeg vec idx deg t ps0 =
  let p = ps0!!idx
  in case deg of
      Lin delta ->
        vec `update` (idx, -delta * p)
      Enz delta km ->
        vec `update` (idx, -delta * p / (km + p))
      Act delta idx2 ->
        vec `update` (idx, -delta * p * (ps0!!idx2))

toODE ps t ps0 =
  foldl prispevek
        (replicate (length ps) 0)
        (zip [0..] ps)
  where
    prispevek acc (idx,(c,d)) =
      toVecDeg (toVecIzr acc idx c t ps0) idx d t ps0


----------------------
-- Pomozne funkcije --
----------------------

update [] _         = []
update (x:xs) (0,v) = x+v : xs
update (x:xs) (i,v) = x : update xs ((i-1),v)
infixl 5 `update`

heaviside x 
    | x < 0     = 0
    | x == 0    = 0.5
    | otherwise = 1
