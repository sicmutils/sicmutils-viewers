(ns emmy.mafs.core
  "Server-side rendering functions for the components declared in the
  [`mafs.core`](https://cljdoc.org/d/org.mentat/mafs.cljs/CURRENT/api/mafs.core)
  namespace of the [`Mafs.cljs` project](https://mafs.mentat.org)."
  (:refer-clojure :exclude [vector])
  (:require [emmy.viewer :as ev]))

(defn ^:no-doc default-viewer [child]
  (ev/render
   ['mafs.core/Mafs
    ['mafs.coordinates/Cartesian]
    child]))

(defn ^:no-doc tagged
  ([v] (ev/tagged v default-viewer))
  ([v viewer]
   (ev/tagged v viewer)))

(def Theme
  {:red "var(--mafs-red)"
   :orange "var(--mafs-orange)"
   :green "var(--mafs-green)"
   :blue "var(--mafs-blue)"
   :indigo "var(--mafs-indigo)"
   :violet "var(--mafs-violet)"
   :pink "var(--mafs-pink)"
   :yellow "var(--mafs-yellow)"})

(defn mafs
  "
  - `:width`
  - `:height`
  - `:pan`
  - `:viewbox`
  - `:preserve-aspect-ratio`
  - `:ssr`
  "
  [& children]
  (let [[opts children] (ev/split-opts children)]
    (tagged
     (into ['mafs.core/Mafs opts] children)
     ev/render)))

(defn point
  "
  - `:x`
  - `:y`
  - `:color`
  - `:opacity`
  - `:svg-circle-opts`"
  [opts]
  (tagged ['mafs.core/Point opts]))

(defn polygon
  "
  - `:points`
  - `:svg-polygon-opts`
  - `:color`
  - `:weight`
  - `:fill-opacity`
  - `:stroke-opacity`
  - `:stroke-style`
  "
  [opts]
  (tagged ['mafs.core/Polygon opts]))

(defn ellipse
  "
  - `:center`
  - `:radius`
  - `:svg-ellipse-opts`
  - `:color`
  - `:weight`
  - `:fill-opacity`
  - `:stroke-opacity`
  - `:stroke-style`
  "
  [opts]
  (tagged ['mafs.core/Ellipse opts]))

(defn circle
  "
  - `:center`
  - `:radius`
  - `:svg-ellipse-opts`
  - `:color`
  - `:weight`
  - `:fill-opacity`
  - `:stroke-opacity`
  - `:stroke-style`
  "
  [opts]
  (tagged ['mafs.core/Circle opts]))

(defn text
  "
  - `:x`
  - `:y`
  - `:attach`
  - `:attach-distance`
  - `:size`
  - `:color`
  - `:svg-text-opts`
  "
  [& children]
  (let [[opts children] (ev/split-opts children)]
    (tagged
     (into ['mafs.core/Text opts] children))))

(defn vector
  "
  - `:tail`
  - `:tip`
  - `:svg-line-opts`
  - `:color`
  - `:opacity`
  - `:weight`
  - `:style`
  "
  [opts]
  (tagged ['mafs.core/Vector opts]))

(defn transform
  "
  - `:matrix`
  - `:translate`
  - `:scale`
  - `:rotate`
  - `:shear`
  "
  [& children]
  (let [[opts children] (ev/split-opts children)]
    (tagged
     (into ['mafs.core/Transform opts] children))))

(defn movable-point
  "This version takes an atom and, optionally, a path into the atom.
  - `:atom`
  - `:path` optional
  - `:constrain`
  - `:color`

  also discuss

  - `:point`
  - `:on-move`"
  [opts]
  (tagged ['mafs.core/MovablePoint opts]))
