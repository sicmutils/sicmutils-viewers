(ns emmy.portal.tex
  "TODO should I name this katex?"
  (:require [emmy.portal.css :as css]
            ["katex" :as k]
            [portal.ui.api :as p]
            [portal.ui.rpc :as rpc]
            ["react" :as react]))

(css/inject!
 "https://cdn.jsdelivr.net/npm/katex@0.16.7/dist/katex.min.css")

(defn tex [text]
  [:span
   {:dangerouslySetInnerHTML
    {:__html
     (.renderToString katex text
                      #js {:throwOnError false
                           :displayMode true})}}])

(defn expand [v]
  (let [[form set-form!] (react/useState nil)]
    (react/useEffect
     (fn []
       (-> (rpc/call 'emmy.expression.render/->TeX v)
           (.then set-form!))
       js/undefined)
     #js [])
    (when form
      [tex form])))

(defn show-expr [x]
  (if (string? x)
    [tex x]
    [:f> expand x]))

(p/register-viewer!
 {:name :emmy.portal/tex
  :component show-expr
  ;; TODO modify predicate?
  :predicate (constantly true)})
