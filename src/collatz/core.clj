(ns collatz.core
  (:gen-class)
  (:use seesaw.core
        seesaw.chooser
        )
  (:require [clojure.math.numeric-tower :as math]
             [clojure.java.io :as io]
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

(defn rand-state [len pop]
  (let [p (rand-pool len pop)
        b (get-best-in-pool p)]
    (list b p)))

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

(def gene-evo)
(def extract-worst)

(defn competition
  "determine whether new gene enters pool"
  [[best pool] gene len popl]
  (let [newscore (fitness gene)
        newp (math/expt 2 newscore)
        newone (list newscore gene)
        ;;randi (rand-int popl)
        [worst oks] (extract-worst pool)
        ;;oldscore (first (nth pool randi))
        oldscore (first worst)
        oldp (math/expt 2 oldscore)
        bestscore (first best)
        newbest (if (< newscore bestscore) best newone)
        rat (/ oldp (+ oldp newp))
        check (< (rand) rat)
        newpool (if check pool
                    ;;(assoc (vec pool) randi newone)
                    (conj oks newone)
                    )]
    (when (> newscore bestscore)
      (println (str newbest))
      (future (gene-evo gene)))
    (list newbest newpool)))

(defn mate [[best pool :as pair] len popl mrate g1 g2]
  (let [cut (find-overlap g1 g2)
        splice (concat (subvec g1 0 cut) (subvec g2 cut))
        newgene (mutate mrate splice)
        old? (some #(= newgene (second %)) pool)]
    (if old?
      pair
      (competition pair newgene len popl))))

(def seq-diff)

(defn met-evolve [[best pool :as pair]]
  (let [len (count (second (first pool)))
        popl (count pool)
        mrate (/ 3 len)
        g1 (vec (second (nth pool (rand-int popl))))
        g2 (vec (second (nth pool (rand-int popl))))
        close? (< (seq-diff g1 g2) (* 0.1 len))]
    (if close? pair
        (mate pair len popl mrate g1 g2))))

(defn generations "iterate met-evolve for n generations applied to st"
  [st n]
  (loop [st st
         n n]
    (if (= n 0) st
        (recur (met-evolve st) (dec n)))))

(defn alts "count alternations on right, looking for 1 at start" [seq]
  (loop [[a & tail] (reverse seq)
         l 0
         b 1]
    (if (= a b)
      (recur tail (inc l) (- 1 b))
      l)))

(defn padline "pad line of bits to given length with r 0s on right"
  [len [seq r]]
  (let [l (- len r (count seq))]
    (concat (repeat l 0) seq (repeat r 0))))

(defn evo "make tab of evolution of seq (from gene)" [gene]
  (let [rec (loop [seq (gene->seq gene)
                   r 0
                   tab '()]
              (if (< (count seq) 67) (conj tab (list seq r))
                  (let [z (dec (alts seq))]
                    (recur (trim0 (nextstep seq)) (+ r z) (conj tab (list seq r))))))
        len (apply max (map #(+ (count (first %)) (second %)) rec))]
        (map #(padline len %) rec)))

(defn mk-gp-cmd "make string for gnuplot to display file fn" [fn h w]
  (list "gnuplot" "-p" "-e" (str "``set " "xrange " "[1:" (str w) "]; " "set " "yrange " "[1:" (str h) "]; " "plot '" fn "' matrix with image``")))

(defn tab->str [tab]
  (let [lines (map #(clojure.string/join "\t" %) tab)]
    (clojure.string/join "\n" lines)))

(set! *print-length* nil)
(defn save-tab
  [tab fnm]
  (let [fw (io/file fnm)
        clobber (atom false)]
    (if (.exists fw)
      (let [dial (dialog :content (str "File " fw " exists, do you want to overwrite it?")
                         :success-fn (fn [e] (reset! clobber true))
                         :type :warning
                         :option-type :yes-no)]
        (pack! dial)
        (show! dial))
      (reset! clobber true))
    (if (and fw @clobber)
      (spit fw (tab->str tab)))))

(defn viewtab
  ([tab] (viewtab tab (str (java.time.LocalDateTime/now))))
  ([tab fnm]
   (save-tab tab fnm)
   (let [h (count tab)
         w (count (first tab))
         cmd (mk-gp-cmd fnm h w)]
     (apply clojure.java.shell/sh cmd))))
   
(defn gene-evo [g]
  (viewtab (evo g)))

(defn get-stats
  ([b p]
   (let [fs (map first p)
         [mn m mx] (get-stats fs)]
     (list mn m mx (first b))))
  ([lst]
   (let [mn (apply min lst)
         mx (apply max lst)
         m (/ (apply + lst) (count lst))]
     (list mn m mx))))

(def ave-dist-pool)

(defn gogen [st n]
  (swap! st generations n)
  (println (apply get-stats @st))
  (println (ave-dist-pool (second @st))))

(defn join-states [[[f1 b1 :as bg1] p1 :as s1] & ss]
  (if (empty? ss)
    s1
    (let [[[[f2 b2 :as bg2] p2] & sss] ss
          bnew (if (> f1 f2) bg1 bg2)
          pnew (concat p1 p2)]
      (apply join-states (list bnew pnew) sss))))

(defn gotarget [st targ]
  (let [len (count (second (first @st)))]
    (while (< (first (first @st)) targ)
      (swap! st generations 1000)
      ;; (let [d (ave-dist-pool (second @st))]
      ;;   (when (< (/ d len) 0.45)
      ;;     (let [extra (rand-state len 2)]
      ;;       (swap! st join-states extra)
      ;;       (println "d: " d ", population: " (count (second @st))))))
      )
    (println (apply get-stats @st))
    ))

(defn extract-worst "returns worst entry of pool and the pool without that entry"
  [pool]
  (loop [res '()
         worst (first pool)
         rem (rest pool)]
    (if (= (count rem) 0)
      (list worst res)
      (let [[h & t] rem]
        (if (< (first worst) (first h))
          (recur (conj res h) worst t)
          (recur (conj res worst) h t))))))

(defn seq-diff
  "count number of places where seq's differ; must have same length"
  [s1 s2]
  (count (filter #(= % false) (map = s1 s2))))

(defn mean [l] (/ (apply + l) (count l)))

(defn total-dist [p pool]
  (apply + (map #(seq-diff p %) pool)))

(defn ave-dist-pool
  "calculate the average distance across the given pool"
  ([pool]
   (let [len (count pool)
        div (/ (* len (dec len)) 2.0)
        genes (map second pool)
        sum (loop [s 0
                   [h & t] genes]
              ;;(println (count t))
              (if (empty? t)
                s
                (recur (+ s (total-dist h t)) t)))]
     (/ sum div)))
  ([p1 p2]
   (let [g1 (map second p1)
         g2 (map second p2)
         pop1 (count g1)
         pop2 (count g2)
         div (* pop1 pop2)
         dlist (map #(total-dist % g2) g1)]
     (/ (apply + dlist) div))))

(defn ave-dist
  ([state] (ave-dist-pool (second state)))
  ([s1 s2] (ave-dist-pool (second s1) (second s2))))
  
    
(defn save-state [st]
  (let [ststr (prn-str st)
        fw (choose-file :type :save)
        clobber (atom false)]
    (if (.exists fw)
      (let [dial (dialog :content (str "File " fw " exists, do you want to overwrite it?")
                         :success-fn (fn [e] (reset! clobber true))
                         :type :warning
                         :option-type :yes-no)]
        (pack! dial)
        (show! dial))
      (reset! clobber true))
    (if (and fw @clobber)
      (spit fw ststr))))

(defn read-state []
  (let [stfile (choose-file :type :open
                             :selection-mode :files-only)]
    (when stfile
      (read-string (slurp stfile)))))

(defn gene-sensitivity [gene]
  (let [len (count gene)
        gv (vec gene)
        vars (map #(assoc gv % (- 1 (get gv %))) (range len))
        fs (map #(list (fitness %1) %2) vars (range len))]
    fs))

(defn new6 [len]
  (let [pool (rand-pool len len)
        state (atom (list (get-best-in-pool pool) pool))]
    (gotarget state 6)
    @state))

(defn -main
  "I don't do a whole lot ... yet."
  [& args]
  (println "Hello, World!"))
