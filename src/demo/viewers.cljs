(ns demo.viewers
  (:require [demo.jsxgraph]
            [demo.mathbox]
            [demo.mathlive]
            [jsxgraph.core]
            [mathbox]
            [mathbox.primitives]
            [mathlive.core]
            [nextjournal.clerk.sci-env]
            [nextjournal.clerk.static-app]
            [nextjournal.clerk.trim-image]
            [sci.core :as sci]
            [sci.ctx-store]
            [sicmutils.env.sci :as sicm.sci]
            [sicmutils.value]))

(sci.ctx-store/swap-ctx!
 sci/merge-opts
 {:classes {'Math js/Math}
  :namespaces
  (assoc
   (:namespaces sicm.sci/context-opts)
   ;; These are the three provided by this library.
   'demo.mathbox  (sci/copy-ns demo.mathbox (sci/create-ns 'demo.mathbox))
   'demo.jsxgraph (sci/copy-ns demo.jsxgraph (sci/create-ns 'demo.jsxgraph))
   'demo.mathlive (sci/copy-ns demo.mathlive (sci/create-ns 'demo.mathlive))

   ;; These come from external deps, now.
   'mathbox (sci/copy-ns mathbox (sci/create-ns 'mathbox))
   'jsxgraph.core (sci/copy-ns jsxgraph.core (sci/create-ns 'jsxgraph.core))
   'mathlive.core (sci/copy-ns mathlive.core (sci/create-ns 'mathlive.core))
   'mathbox.primitives
   (sci/copy-ns mathbox.primitives (sci/create-ns 'mathbox.primitives)))
  :aliases {'e 'sicmutils.env
            'mb 'demo.mathbox
            'box 'mathbox.primitives
            'jsx 'jsxgraph.core
            'ml  'mathlive.core}})
