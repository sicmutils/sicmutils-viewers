;; ## Stationary Integrals of a Simple Bézier Curve

(ns examples.minimal-bezier
  {:nextjournal.clerk/toc true}
  (:refer-clojure
   :exclude [+ - * / = zero? compare
             numerator denominator ref partial
             infinite? abs])
  (:require [emmy.clerk :as ec]
            [emmy.env :refer :all :exclude [point]]
            [emmy.mafs :as mafs]
            [nextjournal.clerk :as clerk]))

{::clerk/width :wide}

^{::clerk/visibility {:code :hide :result :hide}}
(ec/install!)

(def show (comp clerk/tex ->TeX simplify))

;; GJS has a demonstration in SICM using an interpolating
;; polynomial to represent an unknown path between two points.
;; We think we can do a little better with a Bézier curve.
;; Let's set up some of the math for such curves. For the
;; present, we will think of a curve as a function of
;; $t\in[0,1]$. For starters, a fixed point is just a
;; constant function.

(def point constantly)

;; Given two curves, we can "lerp", or linearly interpolate,
;; between them, by making a new function of $t$ out of two
;; existing ones, in which we smoothly blend from one to another.
;; This is also called the *convex combination* of the objects.
;;
;; We'll use Emmy's generic arithmetic to do this. It knows how
;; to work on vectors with no fuss.

(defn lerp [v w]
  (fn [t]
    (+ (* (- 1 t) (v t))
       (* t (w t)))))

;; We should expect that to create a smooth line between two
;; fixed points:

(def l1
  (let [p (point [0 0])
        q (point [1 1])]
    (lerp p q)))


(mapv l1 (range 0 1 0.1))

;; looks good.

;; To get this into the Mafs world, we want a JS function. How can
;; we get one?
;;

(show ((lerp (point [3 4]) (point [5 6])) 't))

;; seems to be working. OK. Now let's add the control points.

(l1 't)

;; making a bezier is basically a tableau of lerps.

(defn bez [p1 c1 c2 p2]
  (let [lerp1 (lerp p1 c1)
        lerp2 (lerp c1 c2)
        lerp3 (lerp c2 p2)
        lerp12 (lerp lerp1 lerp2)
        lerp23 (lerp lerp2 lerp3)]
    (lerp lerp12 lerp23)))

;; it looks something like this:

;; ```
;; p1 l1 l12 l123.
;; c1 l2 l23.
;; c2 l3.
;; p2.
;; ```

;; so we might consider making that explicit,
;; redoing it all as an iterated lazy sequence of the
;; inputs, and taking the head of the 3rd column (counting
;; from zero).

;; Returning to the experiment: Let's make a Bézier that
;; goes from (0, 1) to (1, 0), with control points at
;; (0, 0) and (1, 1).

(def bez1
  (bez (point [0 1])
       (point [0 0])
       (point [1 1])
       (point [1 0])))

(show (bez1 't))

;; We need to visualize this to find out how it works. Let's
;; start with something simple and non-interactive, just to see
;; if the curve we have computed makes any sense.

(defn bezier-viewer
  [{:keys [view-box p1 c1 c2 c2-color p2 f]
    :or {c2-color :blue}}]
  (mafs/mafs
   {:view-box view-box}
   (mafs/cartesian {:subdivisions 4})
   (mafs/parametric {:xy f :t [0 1]})
   (mafs/point p1 {:color :red})
   (mafs/point c1 {:color :blue})
   (mafs/point c2 {:color c2-color})
   (mafs/point p2 {:color :red})))

;; We can apply it here:

(bezier-viewer
 {:view-box {:x [-2 2] :y [-2 2]}
  :p1 [0 1] :p2 [1 0]
  :c1 [0 0] :c2 [1 1]
  :f bez1})

;; That looks right. The endpoints are in red, and the control points in blue
;; show how the Bézier curve has forced the tangent lines to be vertical at
;; the approaches to the endpoints. Now let's consider how to get to the
;; interesting part of the investigation: finding the control points that
;; extremize integrals of the bezier curve.  For example, where should we put
;; the control points if we want the bezier curve to have minimum length?

(defn bez2 [[c_1x c_1y c_2x c_2y]]
  (bez (point [0 1])
       (point [c_1x c_1y])
       (point [c_2x c_2y])
       (point [1 0])))

(defn objective [state]
  (let [f (bez2 state)
        Df (D f)]
    (definite-integral
      #(sqrt (square (Df %))) 0 1
      {:method :bulirsch-stoer-closed})))

(objective [0 0 1 1])
(objective [0.51 0.4 0.5 0.5])

;; That says the length of the curve with the control points as given is
;; about 1.68. We know the right answer is $\sqrt 2$. How to get there?

(def result1
  (multidimensional-minimize objective [0 0 1 1] {:info? true}))

;; Let's plug that data back into the viewer and see what we find.

(let [[c_1x c_1y c_2x c_2y] (:result result1)]
  (bezier-viewer
   {:view-box {:x [-0.5 5] :y [-0.5 2]}
    :p1 [0 1] :p2 [1 0]
    :c1 [c_1x c_1y]
    :c2 [c_2x c_2y]
    :f (bez2 (:result result1))}))

;; That does look like a straight line, doesn't it.
;; Now, by changing the objective function, we should be able
;; to find other interesting mathematical results. Let's try
;; the [brachistochrone problem](https://en.wikipedia.org/wiki/Brachistochrone_curve).
;;
;; To do this we want to find the $f$ minimizing
;; $$\int_0^1 ds\over{v}$$
;; where $ds^2 = dx^2 + dy^2$ and v is the speed of
;; descent, which is proportional to $y$. Traditional
;; approaches to this problem have the y axis pointing
;; downward, so that velocity increases as the object
;; falls.

;; There's a singularity at y=0 that we want to avoid, so
;; for this one we shift everything up by one unit.
;; Lame?

(defn bez3 [[c_1x c_1y c_2x c_2y]]
  (bez (point [0 2])
       (point [c_1x c_1y])
       (point [c_2x c_2y])
       (point [4 1])))

(def result2
  (multidimensional-minimize
   (fn [state]
     (let [f (bez3 state)
           Df (D f)]
       (definite-integral
         (fn [t]
           (let [[_x y] (f t)
                 df (Df t)]
             (sqrt (/ (square df) y))
             ))
         0 1 {:method :bulirsch-stoer-closed})))
   [0.5 0.5 0.5 0.5]
   {:info? true
    :nonzero-delta 0.5
    :zero-delta 0.1}))

;; Let's see what happens when we plug that in:

(let [[c_1x c_1y c_2x c_2y] (:result result2)]
  (bezier-viewer
   {:view-box {:x [-0.5 5] :y [-0.5 2]}
    :p1 [0 2] :p2 [4 1]
    :c1 [c_1x c_1y]
    :c2 [c_2x c_2y]
    :c2-color :green
    :f (bez3 (:result result2))}))

;; The physics folks tell us the answer is a cycloid, and the picture
;; in Wikipedia shows that this 1:4 shape ought to dip below the zero
;; point (if we look at this graph upside down.) But it's not entirely
;; convincing. We might need more control points for this, and/or a
;; background curve that represents the analytical solution.

;; Let's grasp the nettle and make a function that will allow the
;; addition of another control point.

;; Alternately, we now have what we need to wire up movable points.
