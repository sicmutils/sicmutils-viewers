goog.provide('sicmutils.numerical.unimin.bracket');
sicmutils.numerical.unimin.bracket.epsilon = 1.0E-21;
/**
 * Returns the points ordered as f(a) < f(b)
 */
sicmutils.numerical.unimin.bracket.ascending_by = (function sicmutils$numerical$unimin$bracket$ascending_by(f,a,b){
var fa = (f.cljs$core$IFn$_invoke$arity$1 ? f.cljs$core$IFn$_invoke$arity$1(a) : f.call(null,a));
var fb = (f.cljs$core$IFn$_invoke$arity$1 ? f.cljs$core$IFn$_invoke$arity$1(b) : f.call(null,b));
if((fa < fb)){
return new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [a,fa], null),new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [b,fb], null)], null);
} else {
return new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [b,fb], null),new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [a,fa], null)], null);
}
});
/**
 * Accepts three pairs of `[x, (f x)]`, fits a quadratic function to all three
 *   points and returns the step from `xb` (the coordinate of the second argument)
 *   to the minimum of the fitted quadratic.
 * 
 *   Returns the numerator and denominator `p` and `q` of the required step. If `q`
 *   is 0, then the supplied points were colinear.
 * 
 *   `q` is guaranteed to be `>= 0`, while `p` might be negative.
 * 
 *   See these notes for the derivation of this method:
 *   http://fourier.eng.hmc.edu/e176/lectures/NM/node25.html
 */
sicmutils.numerical.unimin.bracket.parabolic_pieces = (function sicmutils$numerical$unimin$bracket$parabolic_pieces(p__102143,p__102144,p__102145){
var vec__102146 = p__102143;
var xa = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102146,(0),null);
var fa = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102146,(1),null);
var vec__102149 = p__102144;
var xb = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102149,(0),null);
var fb = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102149,(1),null);
var vec__102152 = p__102145;
var xc = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102152,(0),null);
var fc = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102152,(1),null);
var _PERCENT_ = (function (){var tmp1 = ((xb - xa) * (fb - fc));
var tmp2 = ((xb - xc) * (fb - fa));
var v = (tmp2 - tmp1);
var p = (((xb - xc) * tmp2) - ((xb - xa) * tmp1));
var q = (2.0 * v);
if((q > (0))){
return new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [sicmutils.generic.negate.cljs$core$IFn$_invoke$arity$1(p),q], null);
} else {
return new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [p,sicmutils.generic.abs.cljs$core$IFn$_invoke$arity$1(q)], null);
}
})();
if(cljs.core.truth_((function (p1__102138_SHARP_){
return (cljs.core.second(p1__102138_SHARP_) >= (0));
}))){
} else {
throw (new Error("Assert failed: (fn* [p1__102138#] (>= (second p1__102138#) 0))"));
}

return _PERCENT_;
});
/**
 * Fits a parabola through all three points, and returns the coordinate of the
 *   minimum of the parabola.
 * 
 *   If the supplied points are colinear, returns a point that takes a large jump
 *   in the direction of the downward slope of the line.
 */
sicmutils.numerical.unimin.bracket.parabolic_step = (function sicmutils$numerical$unimin$bracket$parabolic_step(a,p__102161,c){
var vec__102163 = p__102161;
var xb = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102163,(0),null);
var b = vec__102163;
var two_eps = (2.0 * sicmutils.numerical.unimin.bracket.epsilon);
var vec__102166 = sicmutils.numerical.unimin.bracket.parabolic_pieces(a,b,c);
var p = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102166,(0),null);
var q = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102166,(1),null);
var q__$1 = (((q < two_eps))?two_eps:q);
return (xb + (p / q__$1));
});
/**
 * Returns a function that performs steps of bracket extension.
 * 
 *   :grow-limit is the maximum factor that the parabolic interpolation can jump
 *   the function.
 */
sicmutils.numerical.unimin.bracket.bracket_step_fn = (function sicmutils$numerical$unimin$bracket$bracket_step_fn(f,p__102183){
var map__102184 = p__102183;
var map__102184__$1 = cljs.core.__destructure_map(map__102184);
var grow_limit = cljs.core.get.cljs$core$IFn$_invoke$arity$3(map__102184__$1,new cljs.core.Keyword(null,"grow-limit","grow-limit",2130208943),110.0);
return (function (a,p__102192,p__102193){
var vec__102194 = p__102192;
var xb = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102194,(0),null);
var fb = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102194,(1),null);
var b = vec__102194;
var vec__102197 = p__102193;
var xc = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102197,(0),null);
var fc = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102197,(1),null);
var c = vec__102197;
var wlim = (xb + (grow_limit * (xc - xb)));
var w = sicmutils.numerical.unimin.bracket.parabolic_step(a,b,c);
if((((xb <= w)) && ((w <= xc)))){
var fw = (f.cljs$core$IFn$_invoke$arity$1 ? f.cljs$core$IFn$_invoke$arity$1(w) : f.call(null,w));
if((fw < fc)){
return new cljs.core.PersistentVector(null, 3, 5, cljs.core.PersistentVector.EMPTY_NODE, [b,new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [w,fw], null),c], null);
} else {
if((fw > fb)){
return new cljs.core.PersistentVector(null, 3, 5, cljs.core.PersistentVector.EMPTY_NODE, [a,b,new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [w,fw], null)], null);
} else {
var new_c = sicmutils.numerical.unimin.golden.extend_pt(xc,xb);
return new cljs.core.PersistentVector(null, 3, 5, cljs.core.PersistentVector.EMPTY_NODE, [b,c,new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [new_c,(f.cljs$core$IFn$_invoke$arity$1 ? f.cljs$core$IFn$_invoke$arity$1(new_c) : f.call(null,new_c))], null)], null);

}
}
} else {
if((((xc <= w)) && ((w <= wlim)))){
var fw = (f.cljs$core$IFn$_invoke$arity$1 ? f.cljs$core$IFn$_invoke$arity$1(w) : f.call(null,w));
if((fw < fc)){
var new_c = sicmutils.numerical.unimin.golden.extend_pt(w,xc);
return new cljs.core.PersistentVector(null, 3, 5, cljs.core.PersistentVector.EMPTY_NODE, [c,new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [w,fw], null),new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [new_c,(f.cljs$core$IFn$_invoke$arity$1 ? f.cljs$core$IFn$_invoke$arity$1(new_c) : f.call(null,new_c))], null)], null);
} else {
return new cljs.core.PersistentVector(null, 3, 5, cljs.core.PersistentVector.EMPTY_NODE, [b,c,new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [w,fw], null)], null);
}
} else {
if((((xc <= wlim)) && ((wlim <= w)))){
return new cljs.core.PersistentVector(null, 3, 5, cljs.core.PersistentVector.EMPTY_NODE, [b,c,new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [wlim,(f.cljs$core$IFn$_invoke$arity$1 ? f.cljs$core$IFn$_invoke$arity$1(wlim) : f.call(null,wlim))], null)], null);
} else {
var new_c = sicmutils.numerical.unimin.golden.extend_pt(xc,xb);
return new cljs.core.PersistentVector(null, 3, 5, cljs.core.PersistentVector.EMPTY_NODE, [b,c,new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [new_c,(f.cljs$core$IFn$_invoke$arity$1 ? f.cljs$core$IFn$_invoke$arity$1(new_c) : f.call(null,new_c))], null)], null);

}
}
}
});
});
/**
 * Generates an interval `[lo, hi]` that is guaranteed to contain a minimum of the
 *   function `f`, along with a candidate point `[mid, (f mid)]` that the user can
 *   use to start a minimum search.
 * 
 *   Returns a dictionary of the form:
 * 
 *   {:lo `lower end of the bracket`
 * :mid `candidate point`
 * :hi `upper end of the bracket`
 * :fncalls `# of fn evaluations so far`
 * :iterations `total iterations`}
 * 
 *   `:lo`, `:mid` and `:hi` are each pairs of the form `[x, (f x)]`.
 * 
 *   The implementation works by growing the bounds using either:
 * 
 *   - a step outside the bounds that places one bound at the golden-ratio cut
 *   point between the new bounds, or
 *   - a parabola with a minimum interpolated outside the current bounds, bounded b
 *   a max.
 * 
 *   This implementation was ported from `scipy.optimize.optimize.bracket`:
 *   https://github.com/scipy/scipy/blob/v1.5.2/scipy/optimize/optimize.py#L2450
 * 
 *   `bracket-min` supports the following optional keyword arguments:
 * 
 *   `:xa` the initial guess for the lower end of the bracket. Defaults to 0.0.
 * 
 *   `:xb` the initial guess for the upper end of the bracket. Defaults to 1.0. (If
 *   these points aren't supplied in sorted order they'll be switched.)
 * 
 *   `:grow-limit` The maximum factor that the parabolic interpolation can jump the
 *   function. Defaults to 110.0.
 * 
 *   `:maxiter` Maximum number of iterations allowed for the minimizer. Defaults to
 *   1000.
 * 
 *   `:maxfun` Maximum number of times the function can be evaluated before exiting.
 *   Defaults to 1000.
 *   
 */
sicmutils.numerical.unimin.bracket.bracket_min = (function sicmutils$numerical$unimin$bracket$bracket_min(var_args){
var G__102209 = arguments.length;
switch (G__102209) {
case 1:
return sicmutils.numerical.unimin.bracket.bracket_min.cljs$core$IFn$_invoke$arity$1((arguments[(0)]));

break;
case 2:
return sicmutils.numerical.unimin.bracket.bracket_min.cljs$core$IFn$_invoke$arity$2((arguments[(0)]),(arguments[(1)]));

break;
default:
throw (new Error(["Invalid arity: ",cljs.core.str.cljs$core$IFn$_invoke$arity$1(arguments.length)].join('')));

}
});

(sicmutils.numerical.unimin.bracket.bracket_min.cljs$core$IFn$_invoke$arity$1 = (function (f){
return sicmutils.numerical.unimin.bracket.bracket_min.cljs$core$IFn$_invoke$arity$2(f,cljs.core.PersistentArrayMap.EMPTY);
}));

(sicmutils.numerical.unimin.bracket.bracket_min.cljs$core$IFn$_invoke$arity$2 = (function (f,p__102210){
var map__102211 = p__102210;
var map__102211__$1 = cljs.core.__destructure_map(map__102211);
var opts = map__102211__$1;
var xa = cljs.core.get.cljs$core$IFn$_invoke$arity$3(map__102211__$1,new cljs.core.Keyword(null,"xa","xa",-251918777),0.0);
var xb = cljs.core.get.cljs$core$IFn$_invoke$arity$3(map__102211__$1,new cljs.core.Keyword(null,"xb","xb",-1641297359),1.0);
var maxiter = cljs.core.get.cljs$core$IFn$_invoke$arity$3(map__102211__$1,new cljs.core.Keyword(null,"maxiter","maxiter",-1867223754),(1000));
var maxfun = cljs.core.get.cljs$core$IFn$_invoke$arity$3(map__102211__$1,new cljs.core.Keyword(null,"maxfun","maxfun",-50301814),(1000));
var vec__102216 = sicmutils.util.counted.cljs$core$IFn$_invoke$arity$1(f);
var f_counter = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102216,(0),null);
var f__$1 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102216,(1),null);
var step = sicmutils.numerical.unimin.bracket.bracket_step_fn(f__$1,opts);
var stop_fn = (function (_,p__102228,p__102229,iteration){
var vec__102230 = p__102228;
var ___$1 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102230,(0),null);
var fb = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102230,(1),null);
var vec__102233 = p__102229;
var ___$2 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102233,(0),null);
var fc = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102233,(1),null);
return (((iteration > maxiter)) || ((((cljs.core.deref(f_counter) > maxfun)) || ((fb <= fc)))));
});
var complete = (function (p__102238,b,p__102239,iterations){
var vec__102240 = p__102238;
var xa__$1 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102240,(0),null);
var a = vec__102240;
var vec__102243 = p__102239;
var xc = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102243,(0),null);
var c = vec__102243;
var m = new cljs.core.PersistentArrayMap(null, 5, [new cljs.core.Keyword(null,"lo","lo",-931799889),a,new cljs.core.Keyword(null,"mid","mid",-2123385246),b,new cljs.core.Keyword(null,"hi","hi",-1821422114),c,new cljs.core.Keyword(null,"fncalls","fncalls",-71376182),cljs.core.deref(f_counter),new cljs.core.Keyword(null,"iterations","iterations",-1402710890),iterations], null);
if((xc < xa__$1)){
return cljs.core.assoc.cljs$core$IFn$_invoke$arity$variadic(m,new cljs.core.Keyword(null,"lo","lo",-931799889),c,cljs.core.prim_seq.cljs$core$IFn$_invoke$arity$2([new cljs.core.Keyword(null,"hi","hi",-1821422114),a], 0));
} else {
return m;
}
});
var vec__102219 = sicmutils.numerical.unimin.bracket.ascending_by(f__$1,xa,xb);
var vec__102222 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102219,(0),null);
var xb__$1 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102222,(0),null);
var b = vec__102222;
var vec__102225 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102219,(1),null);
var xa__$1 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102225,(0),null);
var a = vec__102225;
var xc = sicmutils.numerical.unimin.golden.extend_pt(xb__$1,xa__$1);
var fc = (f__$1.cljs$core$IFn$_invoke$arity$1 ? f__$1.cljs$core$IFn$_invoke$arity$1(xc) : f__$1.call(null,xc));
var G__102250 = new cljs.core.PersistentVector(null, 3, 5, cljs.core.PersistentVector.EMPTY_NODE, [a,b,new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [xc,fc], null)], null);
var vec__102251 = G__102250;
var a__$1 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102251,(0),null);
var b__$1 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102251,(1),null);
var c = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102251,(2),null);
var iteration = (0);
var G__102250__$1 = G__102250;
var iteration__$1 = iteration;
while(true){
var vec__102257 = G__102250__$1;
var a__$2 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102257,(0),null);
var b__$2 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102257,(1),null);
var c__$1 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102257,(2),null);
var iteration__$2 = iteration__$1;
if(stop_fn(a__$2,b__$2,c__$1,iteration__$2)){
return complete(a__$2,b__$2,c__$1,iteration__$2);
} else {
var G__102328 = step(a__$2,b__$2,c__$1);
var G__102329 = (iteration__$2 + (1));
G__102250__$1 = G__102328;
iteration__$1 = G__102329;
continue;
}
break;
}
}));

(sicmutils.numerical.unimin.bracket.bracket_min.cljs$lang$maxFixedArity = 2);

/**
 * Identical to bracket-min, except brackets a maximum of the supplied fn.
 */
sicmutils.numerical.unimin.bracket.bracket_max = (function sicmutils$numerical$unimin$bracket$bracket_max(var_args){
var G__102261 = arguments.length;
switch (G__102261) {
case 1:
return sicmutils.numerical.unimin.bracket.bracket_max.cljs$core$IFn$_invoke$arity$1((arguments[(0)]));

break;
case 2:
return sicmutils.numerical.unimin.bracket.bracket_max.cljs$core$IFn$_invoke$arity$2((arguments[(0)]),(arguments[(1)]));

break;
default:
throw (new Error(["Invalid arity: ",cljs.core.str.cljs$core$IFn$_invoke$arity$1(arguments.length)].join('')));

}
});

(sicmutils.numerical.unimin.bracket.bracket_max.cljs$core$IFn$_invoke$arity$1 = (function (f){
return sicmutils.numerical.unimin.bracket.bracket_max.cljs$core$IFn$_invoke$arity$2(f,cljs.core.PersistentArrayMap.EMPTY);
}));

(sicmutils.numerical.unimin.bracket.bracket_max.cljs$core$IFn$_invoke$arity$2 = (function (f,opts){
var _f = cljs.core.comp.cljs$core$IFn$_invoke$arity$2(sicmutils.generic.negate,f);
return sicmutils.numerical.unimin.bracket.bracket_min.cljs$core$IFn$_invoke$arity$2(_f,opts);
}));

(sicmutils.numerical.unimin.bracket.bracket_max.cljs$lang$maxFixedArity = 2);

/**
 *  Given a function f, a starting point and a step size, try to bracket a local
 *   extremum for f.
 * 
 *   Return a list (retcode a b c fa fb fc iter-count) where a < b < c, and fa, fb,
 *   fc are the function values at these points. In the case of a minimum, fb
 *   <= (min fa fc); the opposite inequality holds in the case of a maximum.
 * 
 *   iter-count is the number of function evaluations required. retcode is 'okay if
 *   the search succeeded, or 'maxcount if it was abandoned.
 *   
 */
sicmutils.numerical.unimin.bracket.bracket_min_scmutils = (function sicmutils$numerical$unimin$bracket$bracket_min_scmutils(var_args){
var G__102266 = arguments.length;
switch (G__102266) {
case 1:
return sicmutils.numerical.unimin.bracket.bracket_min_scmutils.cljs$core$IFn$_invoke$arity$1((arguments[(0)]));

break;
case 2:
return sicmutils.numerical.unimin.bracket.bracket_min_scmutils.cljs$core$IFn$_invoke$arity$2((arguments[(0)]),(arguments[(1)]));

break;
default:
throw (new Error(["Invalid arity: ",cljs.core.str.cljs$core$IFn$_invoke$arity$1(arguments.length)].join('')));

}
});

(sicmutils.numerical.unimin.bracket.bracket_min_scmutils.cljs$core$IFn$_invoke$arity$1 = (function (f){
return sicmutils.numerical.unimin.bracket.bracket_min_scmutils.cljs$core$IFn$_invoke$arity$2(f,cljs.core.PersistentArrayMap.EMPTY);
}));

(sicmutils.numerical.unimin.bracket.bracket_min_scmutils.cljs$core$IFn$_invoke$arity$2 = (function (f,p__102267){
var map__102268 = p__102267;
var map__102268__$1 = cljs.core.__destructure_map(map__102268);
var start = cljs.core.get.cljs$core$IFn$_invoke$arity$3(map__102268__$1,new cljs.core.Keyword(null,"start","start",-355208981),(0));
var step = cljs.core.get.cljs$core$IFn$_invoke$arity$3(map__102268__$1,new cljs.core.Keyword(null,"step","step",1288888124),(10));
var maxiter = cljs.core.get.cljs$core$IFn$_invoke$arity$3(map__102268__$1,new cljs.core.Keyword(null,"maxiter","maxiter",-1867223754),(1000));
var vec__102269 = sicmutils.util.counted.cljs$core$IFn$_invoke$arity$1(f);
var f_counter = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102269,(0),null);
var f__$1 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102269,(1),null);
var stop_fn = (function (p__102281,p__102282,p__102283,iteration){
var vec__102284 = p__102281;
var _ = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102284,(0),null);
var fa = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102284,(1),null);
var vec__102287 = p__102282;
var ___$1 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102287,(0),null);
var fb = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102287,(1),null);
var vec__102290 = p__102283;
var ___$2 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102290,(0),null);
var fc = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102290,(1),null);
return (((iteration > maxiter)) || ((fb <= (function (){var x__4339__auto__ = fa;
var y__4340__auto__ = fc;
return ((x__4339__auto__ < y__4340__auto__) ? x__4339__auto__ : y__4340__auto__);
})())));
});
var complete = (function (p__102294,b,p__102295,iterations){
var vec__102297 = p__102294;
var xa = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102297,(0),null);
var a = vec__102297;
var vec__102300 = p__102295;
var xc = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102300,(0),null);
var c = vec__102300;
var m = new cljs.core.PersistentArrayMap(null, 6, [new cljs.core.Keyword(null,"lo","lo",-931799889),a,new cljs.core.Keyword(null,"mid","mid",-2123385246),b,new cljs.core.Keyword(null,"hi","hi",-1821422114),c,new cljs.core.Keyword(null,"fncalls","fncalls",-71376182),cljs.core.deref(f_counter),new cljs.core.Keyword(null,"converged?","converged?",1779059976),(iterations <= maxiter),new cljs.core.Keyword(null,"iterations","iterations",-1402710890),iterations], null);
if((xc < xa)){
return cljs.core.assoc.cljs$core$IFn$_invoke$arity$variadic(m,new cljs.core.Keyword(null,"lo","lo",-931799889),c,cljs.core.prim_seq.cljs$core$IFn$_invoke$arity$2([new cljs.core.Keyword(null,"hi","hi",-1821422114),a], 0));
} else {
return m;
}
});
var run = (function (p__102307,b,p__102308,iter){
while(true){
var vec__102309 = p__102307;
var xa = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102309,(0),null);
var a = vec__102309;
var vec__102312 = p__102308;
var xc = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102312,(0),null);
var c = vec__102312;
if(stop_fn(a,b,c,iter)){
return complete(a,b,c,iter);
} else {
var xd = (xc + (xc - xa));
var G__102337 = b;
var G__102338 = c;
var G__102339 = new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [xd,(f__$1.cljs$core$IFn$_invoke$arity$1 ? f__$1.cljs$core$IFn$_invoke$arity$1(xd) : f__$1.call(null,xd))], null);
var G__102340 = (iter + (1));
p__102307 = G__102337;
b = G__102338;
p__102308 = G__102339;
iter = G__102340;
continue;
}
break;
}
});
var vec__102272 = sicmutils.numerical.unimin.bracket.ascending_by(f__$1,start,(start + step));
var vec__102275 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102272,(0),null);
var xb = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102275,(0),null);
var b = vec__102275;
var vec__102278 = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102272,(1),null);
var xa = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__102278,(0),null);
var a = vec__102278;
var xc = (xb + (xb - xa));
return run(a,b,new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [xc,(f__$1.cljs$core$IFn$_invoke$arity$1 ? f__$1.cljs$core$IFn$_invoke$arity$1(xc) : f__$1.call(null,xc))], null),(0));
}));

(sicmutils.numerical.unimin.bracket.bracket_min_scmutils.cljs$lang$maxFixedArity = 2);

/**
 * Identical to bracket-min-scmutils, except brackets a maximum of the supplied
 *   fn.
 */
sicmutils.numerical.unimin.bracket.bracket_max_scmutils = (function sicmutils$numerical$unimin$bracket$bracket_max_scmutils(var_args){
var G__102317 = arguments.length;
switch (G__102317) {
case 1:
return sicmutils.numerical.unimin.bracket.bracket_max_scmutils.cljs$core$IFn$_invoke$arity$1((arguments[(0)]));

break;
case 2:
return sicmutils.numerical.unimin.bracket.bracket_max_scmutils.cljs$core$IFn$_invoke$arity$2((arguments[(0)]),(arguments[(1)]));

break;
default:
throw (new Error(["Invalid arity: ",cljs.core.str.cljs$core$IFn$_invoke$arity$1(arguments.length)].join('')));

}
});

(sicmutils.numerical.unimin.bracket.bracket_max_scmutils.cljs$core$IFn$_invoke$arity$1 = (function (f){
return sicmutils.numerical.unimin.bracket.bracket_max_scmutils.cljs$core$IFn$_invoke$arity$2(f,cljs.core.PersistentArrayMap.EMPTY);
}));

(sicmutils.numerical.unimin.bracket.bracket_max_scmutils.cljs$core$IFn$_invoke$arity$2 = (function (f,opts){
var _f = cljs.core.comp.cljs$core$IFn$_invoke$arity$2(sicmutils.generic.negate,f);
return sicmutils.numerical.unimin.bracket.bracket_min.cljs$core$IFn$_invoke$arity$2(_f,opts);
}));

(sicmutils.numerical.unimin.bracket.bracket_max_scmutils.cljs$lang$maxFixedArity = 2);


//# sourceMappingURL=sicmutils.numerical.unimin.bracket.js.map
