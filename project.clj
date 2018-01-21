(defproject collatz "0.1.0-SNAPSHOT"
  :description "Collatz conjecture explorer"
  :url ""
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.9.0"]
                 [seesaw "1.4.5"]
                 [org.apache.commons/commons-math3 "3.6"]
                 [org.clojure/math.numeric-tower "0.0.4"]
                 [net.mikera/imagez "0.12.0"]
                 ;;[ch.cern/cernlib "2005-1.4"]
                 ]
  :main ^:skip-aot collatz.core
  :target-path "target/%s"
  :profiles {:uberjar {:aot :all}}
  :bin {:name "collatz"
        :bin-path "~/bin"
        :bootclasspath false}
  )
