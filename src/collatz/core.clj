(ns collatz.core
  (:gen-class)
  (:use seesaw.core
        seesaw.chooser
        ;;mikera.image.core
        ;;clojure.math.numeric-tower
        )
  (:require [clojure.math.numeric-tower :as math]
            )
  (:import [org.apache.commons.math3.distribution NormalDistribution]
           [org.apache.commons.math3.analysis.function Log]
           )
  )

(native!)

(defn plus3 "add 3 bits" [a b c]
  (let [sum (+ a b c)
        newbit (mod sum 2)
        carry (if (> sum 1) 1 0)]
    (list carry newbit)))

;; (defn nextstepaux "iterate sequence" [s c acc]
;;   (if (= (count s) 1)
;;     (let [[a] s
;;           [carry newbit] (plus3 a 0 c)]
;;       (concat (list carry newbit) acc))
;;     (let [[a b] s
;;           [carry newbit] (plus3 a b c)
;;           news (drop 1 s)
;;           newacc (conj acc newbit)]
;;       (nextstepaux news carry newacc))))

;; (defn nextstep "wrapper" [s]
;;   (nextstepaux (reverse s) 1 '()))

(defn nextstep "iterate through sequence" [seq]
  (loop [s (reverse seq)
         c 1
         acc '()]
    (if (= (count s) 1)
      (let [[a] s
            [carry newbit] (plus3 a 0 c)]
        (concat (list carry newbit) acc))
      (let [[a b] s
            [carry newbit] (plus3 a b c)
            news (drop 1 s)
            newacc (conj acc newbit)]
        (recur news carry newacc)))))

(defn trimL "remove given char from left of list" [r [h & tail :as rem]]
  (loop [r r
         [h & tail :as rem] rem]
    (if (= h r)
      (recur r tail)
      rem)))

(defn trimR [r l]
  (reverse (trimL r (reverse l))))

(defn trim0 [l]
  (trimL 0 (trimR 0 l)))

;; (defn lifeaux "calc life using accumulator" [s t]
;;   (if (< (count s) 67)
;;     t
;;     (lifeaux (trim0 (nextstep s)) (+ t 1))))
  
(defn life "calculate the life of a seq till length is < 67" [seq]
  (loop [s seq
         t 0]
    (if (< (count s) 67) t
        (recur (trim0 (nextstep s)) (+ t 1)))))

(defn rand-seq "generate random seq of 0/1 of length len" [len]
  (map rand-int (repeat len 2)))

(defn mutatebit "change bit with given prob" [p b]
  (if (< (rand) p ) (- 1 b) b))
(defn mutate "change bits of seq with given probability" [p seq]
  (map (fn [b] (mutatebit p b)) seq))

(defn gene->seq "derive seq from gene" [g]
  (concat (trimL 0 g) '(1)))

(defn normalcdf [m sd x]
  (-> (new NormalDistribution m sd) (.cumulativeProbability x)))

(defn log "natural log" [x]
  (-> (new Log) (.value x)))

(def leftVel (- (/ (log 3.0) (log 2.0)) 1))

(defn probSum "prob that sum of given number of GeoRVs is less than given"
  [sum num]
  (let [m (* num (- 1 leftVel))
        sd (math/sqrt (* 2 num))]
    (normalcdf m sd sum)))

(defn meansdLife "returns mean and sd of time to shrink by given amount"
  [len]
  (let [l (map inc (range (* 100 len)))
        m1list (map #(probSum len %) l)
        m1 (inc (reduce + m1list))
        m2list (map #(* (+ 1 (* 2 %)) (nth m1list (dec %))) l)
        m2 (inc (reduce + m2list))
        sd (math/sqrt (- m2 (* m1 m1)))]
    (list m1 sd)))

(def meansdLife-mem (memoize meansdLife))

(defn rarity "# sd's above mean of lifetime" [seq]
  (let [len (count seq)]
    (if (= (count seq) 1) 0
        (let [n (life seq)
              [m sd] (meansdLife-mem (- len 66))]
          (if (= sd 0) (println seq))
          (/ (- n m) sd)))))

(defn fitness [g] 
  (rarity (gene->seq g)))

(defn rand-pool-1 "make one entry in a pool, given length" [len]
  (let [g (rand-seq len)
        f (fitness g)]
    (list f g)))

(defn rand-pool "make random pool of given size and length"
  ([len size] (rand-pool len size []))
  ([len size res] (if (= (count res) size)
                    res
                    (rand-pool len size (conj res (rand-pool-1 len))))))

(defn get-best-in-pool [pool]
  (let [mx (atom -10)
        bg (atom nil)]
    (doseq [[f g :as cur] pool]
      (when (>= f @mx)
        (reset! mx f)
        (reset! bg cur)))
    @bg))

;; (defn detect-aux [a b]
;;   (* b (inc a)))
  

;; (defn detect-ones [l]
;;   (let [ones (reductions detect-aux l)
;;         mx (apply max ones)
;;         ind (

(defn find-overlap [s1 s2]
  (let [bi (atom -1)
        prev (atom 0)
        best (atom 0)
        tab (apply mapv vector (list s1 s2 (range (count s1))))]
    (doseq [[a b i] tab]
      (if (= a b)
        (do
          (swap! prev inc)
          (when (> @prev @best)
            (reset! best @prev)
            (reset! bi i)))
        (reset! prev 0)))
      @bi))
            
(defn competition "determine whether new gene enters pool" [pool best gene len popl]
  (let [newscore (fitness gene)
        newone (list newscore gene)
        randi (rand-int popl)
        oldscore (first (nth pool randi))
        bestscore (first best)
        newbest (if (< newscore bestscore) best newone)
        test (< (rand) (/ oldscore (+ oldscore newscore)))
        newpool (if test pool
                    (assoc (vec pool) randi newone))]
    (when (>= newscore bestscore)
      (println (str newbest)))
    (list newbest newpool)))

(defn met-evolve [[best pool :as pair]]
  (let [len (count (second (first pool)))
        popl (count pool)
        g1 (vec (second (nth pool (rand-int popl))))
        g2 (vec (second (nth pool (rand-int popl))))
        cut (find-overlap g1 g2)
        splice (concat (subvec g1 0 cut) (subvec g2 cut))
        mrate (/ 3 len)
        newgene (mutate mrate splice)
        old? (some #(= newgene (second %)) pool)]
    (if old?
      pair
      (competition pool best newgene len popl))))

(defn generations "iterate met-evolve for n generations applied to st"
  [st n]
  (loop [st st
         n n]
    (if (= n 0) st
        (recur (met-evolve st) (dec n)))))
              
(defn -main
  "I don't do a whole lot ... yet."
  [& args]
  (println "Hello, World!"))
