(ns demo.mathbox
  (:require [mathbox]
            #_[demo.bench]
            [mathbox.primitives :as box]
            [nextjournal.clerk.sci-viewer :as sv]
            [reagent.core :as r :include-macros true]
            [sicmutils.env :as e]
            [sicmutils.expression.compile :as xc]
            [sicmutils.mechanics.lagrange :as l]
            [sicmutils.numerical.ode :as ode]
            [sicmutils.structure :as struct]
            ["three/src/math/Color.js" :refer [Color]]))

;; ## Components

(defn format-number [x]
  (-> x (.toFixed 2) (.replace #"\.0+$" "")))

(defn Cartesian [opts & children]
  (let [opts (update opts :range
                     (fn [{:keys [x y z]}]
                       [(or x [0 1])
                        (or y [0 1])
                        (or z [0 1])]))]
    (into [box/Cartesian opts] children)))

(defn Volume
  "Like the original, but allows for a 'dimensions' key."
  [{:keys [dimensions] :as opts}]
  (let [[width height depth] dimensions
        opts (-> (dissoc opts :dimensions)
                 (assoc :width width
                        :height height
                        :depth depth))]
    [box/Volume opts]))

;; ## Demo-Specific Components
;;
;; These are higher-level components that we'll include for the demos.

(defn ColorCube
  [{:keys [dimensions size opacity]}]
  [:<>
   [Volume
    {:dimensions dimensions
     :items 1
     :channels 4
     :live false
     :expr (fn [emit x y z]
             (emit x y z opacity))}]
   [box/Point
    {:points "<"
     :colors "<"
     :color 0xffffff
     :size size}]])

(defn Function1 [{:keys [samples f] :or {samples 256}}]
  (let [f' (xc/sci-eval f)]
    [:<>
     [box/Interval
      {:width samples
       :channels 2
       :expr (fn [emit x _ time]
               (emit x (f' x time)))}]
     [box/Line {:color 0x3090ff :width 4}]
     [box/Point {:color 0x3090ff :size 8}]]))

(defn Lagrangian-updater
  "hardcoded at first for this use case."
  [state-derivative initial-state & {:keys [compile?]
                                     :or {compile? false}}]
  (let [{:keys [integrator equations]}
        (ode/integration-opts (constantly state-derivative)
                              []
                              initial-state
                              {:epsilon 1e-6
                               :compile? compile?})]
    (fn [[t :as state] t2]
      (let [s (into-array (flatten state))
            output (.solve integrator equations t s t2 nil)]
        (struct/unflatten (.-y ^js output) state)))))

(defn state-deriv* [g m l]
  (l/Lagrangian->state-derivative
   (l/L-pendulum g m l)))

(defn Lagrangian-collector
  "hardcoded at first for this use case."
  [state-derivative initial-state {:keys [compile? parameters]
                                   :or {compile? false}}]
  (let [{:keys [integrator equations]}
        (ode/integration-opts state-derivative
                              (if parameters
                                @parameters
                                [])
                              initial-state
                              {:epsilon 1e-6
                               :compile? compile?})]
    (set! (.-denseOutput integrator) true)
    (fn [state t2 step-size emit]
      (.solve integrator equations 0 state t2
              (.grid integrator step-size
                     (fn [_ y]
                       (emit (aget y 1) (aget y 2))))))))

(defn Tail [{:keys [length] :as opts}]
  [:<>
   [box/Spread {:height [0 0 -0.02] :alignHeight -1}]
   ;; Ah, this is the color channel, and fades out the tail as you go.
   [box/Array
    {:width length
     :channels 4
     :expr (fn [emit i]
             (emit 1 1 1 (- 1 (/ i 16))))}]
   [box/Transpose {:order "yx"}]
   [box/Point
    (-> (dissoc opts :length)
        (assoc :points "<<<"
               :colors "<"))]])

(defn Comet
  "Path is a function of i, t
  dimensions is how many you want to emit
  history is tail length,
  rest of options go to the final point

  Note that i think we have to emit with xzy?? weird..."
  [{:keys [dimensions path length] :as opts}]
  [:<>
   [box/Array
    {:history length
     :channels dimensions
     :expr (fn [emit i t]
             (apply emit (path i t)))}]
   [Tail
    (dissoc opts :dimensions :path)]])

;; NOTE: it comes across the wire ALREADY COMPILED and simplified down
;; to source.
(defn Mass [{:keys [state->xyz L initial-state var-name]}]
  (let [render-fn   (xc/sci-eval state->xyz)
        state-deriv (xc/sci-eval L)
        my-updater  (Lagrangian-updater state-deriv initial-state :compile? false)]
    (r/with-let [!state (r/atom initial-state)]
      [Comet
       {:dimensions 3
        :length 16
        :color 0x3090ff
        :size 10
        :opacity 0.99
        :path
        (fn [_ t]
          (swap! !state #(my-updater % t))
          (when var-name
            (sv/clerk-eval
             (list 'clojure.core/reset!
                   var-name
                   (mapv (fn rec [x]
                           (if (sequential? x)
                             (mapv rec x)
                             x))
                         @!state))))
          (let [[x y z] (render-fn @!state)]
            [x z y]))}])))

(defn DoubleMass
  "Obviously these should be merged!"
  [{:keys [state->xyz L initial-state]}]
  (let [render-fn   (xc/sci-eval state->xyz)
        state-deriv (xc/sci-eval L)
        my-updater  (Lagrangian-updater state-deriv initial-state)]
    (r/with-let [!state (r/atom initial-state)]
      [:<>
       [box/Array
        {:items 1
         ;; because we have two items to emit.
         :width 2
         :channels 3
         :history 16
         :expr
         (fn [emit _ t]
           (swap! !state #(my-updater % t))
           (let [[x1 y1 z1 x2 y2 z2] (render-fn @!state)]
             (emit x1 z1 y1)
             (emit x2 z2 y2)))}]
       [Tail
        {:length 16
         :color 0x3090ff
         :size 10
         :zIndex 1}]])))

(def ^:private two-pi (* 2 Math/PI))

(defn Ellipse [{:keys [a b c]}]
  [:<>
   [box/Area
    {:width 64
     :height 64
     :rangeX [0 two-pi]
     :rangeY [0 two-pi]
     :axes [1 3]
     :expr (fn [emit theta phi _i _j _time]
             (let [sin-theta (Math/sin theta)
                   cos-theta (Math/cos theta)]
               ;; xzy
               (emit (* a sin-theta (Math/cos phi))
                     (* c cos-theta)
                     (* b sin-theta (Math/sin phi)))))
     :items 1
     :channels 3}]
   [box/Surface
    {:shaded true
     :opacity 0.2
     :lineX true
     :lineY true
     :points "<"
     :color 0xffffff
     :width 1}]])

;; ## Hamiltonian Example

(defn sq [x] (* x x))
(def normalize (e/principal-value Math/PI))

(def m 1)
(def g 9.8)
(def l 1)

(defn V [theta]
  (* -1 m g l (Math/cos theta)))

(def potential V)

(defn DoublePendulum
  "For later, here's how to extend this."
  []
  [:<>
   [box/Array
    {:width 2
     :channels 2
     :items 3
     :expr (fn [emit _i now]
             (emit -1 0)
             (emit 0 0)
             (emit (Math/sin now) (- (Math/cos now))))}]
   ;; attach a bob between the two.
   [box/Vector {:color 0xffffff :width 2}]
   ;; JUST attach a point to the second of the two items, ie [1, 2)
   [box/Slice {:items [0 1]}]
   [box/Point {:color 0x909090 :size 4}]

   [box/Slice {:items [1 3]}]
   [box/Point {:color 0xffffff :size 10}]])

(def L
  (l/L-pendulum 9.8 1 1))

(def state-deriv
  (l/Lagrangian->state-derivative L))

(defn Pendulum [!state]
  (r/with-let [updater (Lagrangian-updater state-deriv @!state :compile? true)]
    [:<>
     [box/Array
      {:channels 2
       :items 2
       :expr (fn [emit _i now]
               (swap! !state #(updater % now))
               (let [[_ theta] @!state
                     theta     (normalize theta)]
                 (emit 0 0)
                 (emit (Math/sin theta)
                       (- (Math/cos theta)))))}]
     ;; attach a bob between the two.
     [box/Vector {:color 0xffffff :width 2}]

     [box/Slice {:items [0 1]}]
     [box/Point {:color 0x909090 :size 4}]

     [box/Slice {:items [1 2]}]
     [box/Point {:color 0xffffff :size 10}]]))



(defn WellAxes []
  [:<>
   [box/Axis
    {:axis "x"
     :color 0xffffff}]
   [box/Scale
    {:axis "x"
     :divide 5
     :unit Math/PI
     :base 2
     :start true
     :end true}]
   [box/Format
    {:expr (fn [x] (str (format-number (/ x Math/PI)) "π"))
     :font ["Helvetica"]}]
   [box/Label
    {:color 0xffffff
     :background 0x000000
     :depth 0.5
     :zIndex 1
     :zOrder 5
     :size 10}]
   [box/Axis
    {:axis "y" :color 0xffffff}]
   [box/Scale
    {:axis "y"
     :divide 5
     :unit 1
     :base 10
     :start true
     :end true
     :zero false}]
   [box/Format
    {:expr format-number
     :font ["Helvetica"]}]
   [box/Label
    {:color 0xffffff
     :background 0x000000
     :depth 0.5
     :zIndex 1
     :zOrder 5
     :size 10
     :offset [20 0]}]
   ])

(defn PotentialLine [V]
  [:<>
   ;; This is the potential well. Gotta redo this to make more sense.
   [box/Interval
    {:width 128
     :channels 2
     :live false
     :expr (fn [emit theta]
             (emit theta (V theta)))}]
   [box/Line {:color 0x3090ff}]])

(defn Well [!state]
  [:<>
   [box/Grid
    {:color 0x808080
     :unitX Math/PI
     :baseX 2}]
   [WellAxes]
   [PotentialLine V]

   ;; this is the bead traveling with history along the potential.
   [Comet
    {:dimensions 2
     :length 16
     :color 0xa0d0ff
     :size 5
     :opacity 0.99
     :path
     (fn [_ _]
       (let [[_ theta] @!state]
         [theta (V theta)]))}]])

(defn PhaseAxes []
  [:<>
   [box/Axis
    {:axis "x"
     :color 0xffffff}]
   [box/Scale
    {:axis "x"
     :divide 5
     :unit 1
     :base 10
     :start true
     :end true}]
   [box/Format
    {:expr format-number
     :font ["Helvetica"]}]
   [box/Label
    {:color 0xffffff
     :background 0x000000
     :depth 0.5
     :zIndex 1
     :zOrder 5
     :size 10}]
   [box/Axis
    {:axis "y"
     :color 0xffffff}]
   [box/Scale
    {:axis "y"
     :divide 5
     :unit 1
     :base 10
     :start true
     :end true
     :zero false}]
   [box/Format
    {:expr format-number
     :font ["Helvetica"]}]
   [box/Label
    {:color 0xffffff
     :background 0x000000
     :depth 0.5
     :zIndex 1
     :zOrder 5
     :size 10
     :offset [20 0]}]])

(defn PhaseVectors [!params items]
  (let [dt      3e-2
        updater (Lagrangian-collector state-deriv* [0 0 0]
                                      {:compile? true
                                       :parameters !params})]
    [:<>
     [box/Area
      {:width 16
       :height 16
       :channels 2
       :items items
       :centeredX true
       :centeredY true
       :live true
       :expr
       (fn [emit q p _ _ t]
         (updater #js [0 (normalize (+ q t)) p] (* dt (+ 0.01 (dec items))) dt emit))}]
     [box/Vector
      {:color 0x3090ff
       :size 5
       :end true}]]))

(defn Phase [!state !params items]
  [:<>
   [box/Grid {:color 0x808080}]
   [PhaseAxes]
   [PhaseVectors !params items]
   [Comet
    {:dimensions 2
     :length 16
     :color 0xa0d0ff
     :size 10
     :opacity 0.99
     :path
     (fn [_ _]
       (let [[_ q p] @!state]
         [q p]))}]])

;; SO what have I learned??
;;
;;
(defn Hamilton []
  (r/with-let [!state  (r/atom [0 3 0])
               !params (r/atom [9.8 1 1])
               !items  (r/atom 51)]
    [:<>
     [:input
      {:type :range
       :min 1
       :max 10
       :step 1
       :value (nth @!params 0)
       :on-change
       (fn [target]
         (let [v (.. target -target -value)]
           (swap! !params assoc 0 (js/parseFloat v))))}]
     [:input
      {:type :range
       :min 1
       :max 10
       :step 1
       :value (nth @!params 1)
       :on-change
       (fn [target]
         (let [v (.. target -target -value)]
           (swap! !params assoc 1 (js/parseFloat v))))}]
     [:input
      {:type :range
       :min 1
       :max 10
       :step 1
       :value (nth @!params 2)
       :on-change
       (fn [target]
         (let [v (.. target -target -value)]
           (swap! !params assoc 2 (js/parseFloat v))))}]
     [:input
      {:type :range
       :min 1
       :max 50
       :step 1
       :value @!items
       :on-change
       (fn [target]
         (let [v (.. target -target -value)]
           (reset! !items (js/parseInt v))))}]
     [mathbox/Mathbox
      {:style {:height "600px" :width "100%"}
       :options {:plugins ["core" "controls" "cursor" "stats"]}
       :init  (fn [mb]
                (let [three    (.-three mb)
                      renderer (.-renderer three)]
                  (.setClearColor renderer (Color. 0x000000) 1.0)
                  (.camera mb #js {:proxy true
                                   :position #js [0 0 20]
                                   :fov 90})))}
      [box/Layer
       [box/Unit {:scale 720 :focus 1}
        #_[box/Cartesian
           {:id "pendulum"
            :range [[-1 1] [-1 1]]
            :scale [0.25 0.25]
            :position [-0.5 0.35 0]}
           [Pendulum !state]]

        #_[box/Cartesian
           {:id "well"
            ;; TODO fix our `normalize` so we don't map pi back to negative pi.
            :range [[(- Math/PI) (- Math/PI 0.00001)]
                    [-10 10]]
            :scale [0.48 0.25]
            :position [-0.5 -0.25 0]}
           [Well !state]]

        [box/Cartesian
         {:id "phase"
          :range [[-4 4] [-8 8]]
          :scale [0.6 0.6]
          #_#_:position [0.6 0]}
         [Phase !state !params @!items]
         ]]]]]))
