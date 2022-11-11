(ns user
  (:require [hiccup.page :as hiccup]
            [nextjournal.clerk.config :as config]
            [nextjournal.clerk :as clerk]
            [nextjournal.clerk.view]))

;; To get everything running, first follow the README instructions:
;;
;;```
;; shadow-cljs watch sicm-browser
;;```
;;
;; Then jack in, come here and run the commands in the comment.
;;
;; Better rendering for slides. NOTE that for now, until
;; https://github.com/babashka/sci/issues/832 is fixed, this has to go in
;; `demo.clj` or some non-user namespace. It has to be here so that SCI gets the
;; change before setting its dynamic var.

#_[sicmutils.expression.render :as xr]
#_(alter-var-root
   #'xr/*TeX-vertical-down-tuples*
   (constantly true))

;; Same with my `[sicmutils.env :refer :all]` to get the REPL working.

(defn rebind [^clojure.lang.Var v f]
  (let [old (.getRawRoot v)]
    (.bindRoot v (f old))))

(defonce _ignore
  (rebind
   #'nextjournal.clerk.view/include-viewer-css
   (fn [old]
     (fn []
       (concat
        (list
         (hiccup/include-css
          "https://unpkg.com/mathlive@0.83.0/dist/mathlive-static.css")
         (hiccup/include-css
          "https://unpkg.com/mathlive@0.83.0/dist/mathlive-fonts.css"))
        (old))))))

(def notebooks
  ["src/demo.clj"])

(defn start! []
  (swap! config/!resource->url
         merge
         {"/js/viewer.js" "http://localhost:9000/js/main.js"})
  (clerk/serve!
   {:browse? true
    :watch-paths ["src"]
    :port 7777})
  (Thread/sleep 500)
  (clerk/show! "src/demo.clj"))

(defn github-pages! [_]
  (swap! config/!resource->url merge {"/js/viewer.js" "/sicmutils-clerk/js/main.js"})
  (clerk/build!
   {:paths notebooks
    :bundle? false
    :browse? false
    :out-path "public"}))

(defn publish-local!
  ([] (publish-local! nil))
  ([_]
   (swap! config/!resource->url merge {"/js/viewer.js" "/js/main.js"})
   (clerk/build!
    {:paths notebooks
     :bundle? false
     :browse? false
     :out-path "public"})))

(comment
  (start!)
  (clerk/serve! {:browse? true})
  (publish-local!))

(comment
  ;; If the watcher's not running, call clerk/show on files to be rendered:

  ;; intro:
  (clerk/show! "src/demo.clj")

  ;; Mathbox basics:
  (clerk/show! "src/cube_controls.clj")

  ;; functions:
  (clerk/show! "src/functions.clj")

  ;; symbolic physics:
  (clerk/show! "src/einstein.clj")

  ;; vega, symbolic, double-pendulum
  (clerk/show! "src/pendulum.clj")

  ;; mathbox physics:
  (clerk/show! "src/oscillator.clj")
  (clerk/show! "src/ellipsoid.clj")
  (clerk/show! "src/double_ellipsoid.clj")

  ;; browser/client comms:
  (clerk/show! "src/live_oscillator.clj"))
