goog.provide('sicmutils.polynomial.factor');
/**
 * Given a [[Polynomial]] `p`, returns a sequence of factors of in order of
 *   increasing power.
 * 
 *   The first element is a constant factor, the next is a factor with power 1, and
 *   so on.
 */
sicmutils.polynomial.factor.split_polynomial = (function sicmutils$polynomial$factor$split_polynomial(p){
var answer = (function sicmutils$polynomial$factor$split_polynomial_$_answer(tracker,const$){
var final$ = cljs.core.peek(tracker);
if(sicmutils.value.number_QMARK_(final$)){
return cljs.core.into.cljs$core$IFn$_invoke$arity$2(new cljs.core.PersistentVector(null, 1, 5, cljs.core.PersistentVector.EMPTY_NODE, [final$], null),cljs.core.subvec.cljs$core$IFn$_invoke$arity$2(cljs.core.conj.cljs$core$IFn$_invoke$arity$2(cljs.core.pop(tracker),(1)),(1)));
} else {
return cljs.core.into.cljs$core$IFn$_invoke$arity$2(new cljs.core.PersistentVector(null, 1, 5, cljs.core.PersistentVector.EMPTY_NODE, [const$], null),cljs.core.subvec.cljs$core$IFn$_invoke$arity$2(tracker,(1)));
}
});
var m = (0);
var h = p;
var tracker = cljs.core.PersistentVector.EMPTY;
var old_s = p;
var old_m = (1);
while(true){
if(sicmutils.value.one_QMARK_(m)){
return answer(tracker,h);
} else {
var gg = sicmutils.polynomial.gcd.gcd_Dp(h);
var new_s = sicmutils.polynomial.evenly_divide(h,(sicmutils.polynomial.gcd.gcd.cljs$core$IFn$_invoke$arity$2 ? sicmutils.polynomial.gcd.gcd.cljs$core$IFn$_invoke$arity$2(h,gg) : sicmutils.polynomial.gcd.gcd.call(null,h,gg)));
var new_m = (sicmutils.polynomial.gcd.gcd.cljs$core$IFn$_invoke$arity$2 ? sicmutils.polynomial.gcd.gcd.cljs$core$IFn$_invoke$arity$2(gg,new_s) : sicmutils.polynomial.gcd.gcd.call(null,gg,new_s));
var facts = sicmutils.polynomial.evenly_divide(old_s,new_s);
var doublefacts = (sicmutils.polynomial.gcd.gcd.cljs$core$IFn$_invoke$arity$2 ? sicmutils.polynomial.gcd.gcd.cljs$core$IFn$_invoke$arity$2(facts,old_m) : sicmutils.polynomial.gcd.gcd.call(null,facts,old_m));
var singlefacts = sicmutils.polynomial.evenly_divide(new_s,new_m);
var G__100853 = new_m;
var G__100854 = sicmutils.polynomial.evenly_divide(h,sicmutils.polynomial.mul(new_m,new_s));
var G__100855 = cljs.core.conj.cljs$core$IFn$_invoke$arity$variadic(tracker,doublefacts,cljs.core.prim_seq.cljs$core$IFn$_invoke$arity$2([singlefacts], 0));
var G__100856 = new_s;
var G__100857 = new_m;
m = G__100853;
h = G__100854;
tracker = G__100855;
old_s = G__100856;
old_m = G__100857;
continue;
}
break;
}
});
/**
 * Given some sequence of polynomial factors ordered by increasing power,
 *   symbolically evaluates each power and generates a symbolic expression
 *   representing the product of all factors.
 * 
 *   For example:
 * 
 *   ```clojure
 *   (factors->expression ['c 'x 'y 1 'z])
 *   ;;=> (* c x (expt y 2) (expt z 4))
 *   ```
 */
sicmutils.polynomial.factor.factors__GT_expression = (function sicmutils$polynomial$factor$factors__GT_expression(factors){
var expt = sicmutils.numsymb.symbolic_operator(new cljs.core.Symbol(null,"expt","expt",2123850345,null));
return cljs.core.cons(new cljs.core.Symbol(null,"*","*",345799209,null),cljs.core.map_indexed.cljs$core$IFn$_invoke$arity$2((function (i,f){
if((i === (0))){
return f;
} else {
return (expt.cljs$core$IFn$_invoke$arity$2 ? expt.cljs$core$IFn$_invoke$arity$2(f,i) : expt.call(null,f,i));
}
}),factors));
});
/**
 * Simplifier that flattens nested products, converts singleton calls
 *   like `(* x) => x`, and squashes no-argument products like `(*)` into a
 *   constant `1`.
 */
sicmutils.polynomial.factor.simplify_product = pattern.rule.rule_simplifier.cljs$core$IFn$_invoke$arity$variadic(cljs.core.prim_seq.cljs$core$IFn$_invoke$arity$2([sicmutils.simplify.rules.associative.cljs$core$IFn$_invoke$arity$variadic(cljs.core.prim_seq.cljs$core$IFn$_invoke$arity$2([new cljs.core.Symbol(null,"*","*",345799209,null)], 0)),sicmutils.simplify.rules.unary_elimination.cljs$core$IFn$_invoke$arity$variadic(cljs.core.prim_seq.cljs$core$IFn$_invoke$arity$2([new cljs.core.Symbol(null,"*","*",345799209,null)], 0)),sicmutils.simplify.rules.constant_elimination(new cljs.core.Symbol(null,"*","*",345799209,null),(1)),pattern.rule.rule_STAR_(pattern.rule.pattern_STAR_.cljs$core$IFn$_invoke$arity$2((new cljs.core.List(null,new cljs.core.Symbol(null,"*","*",345799209,null),null,(1),null)),pattern.rule._EQ__GT_),(function (G__100766){
var r__20948__auto__ = (1);
var or__4253__auto__ = r__20948__auto__;
if(cljs.core.truth_(or__4253__auto__)){
return or__4253__auto__;
} else {
return pattern.consequence.succeed(r__20948__auto__);
}
}))], 0));
/**
 * Given a polynomial `p`, and a sequence of variables `vars` (one for each
 *   indeterminate in `p`), returns a symbolic expression representing the product
 *   of all factors of `p`.
 * 
 *   Optionally accepts a `simplify` function that will be called on each factor of
 *   exponent 0, 1, 2 etc. Defaults to `identity`.
 */
sicmutils.polynomial.factor.poly__GT_factored_expression = (function sicmutils$polynomial$factor$poly__GT_factored_expression(var_args){
var G__100778 = arguments.length;
switch (G__100778) {
case 2:
return sicmutils.polynomial.factor.poly__GT_factored_expression.cljs$core$IFn$_invoke$arity$2((arguments[(0)]),(arguments[(1)]));

break;
case 3:
return sicmutils.polynomial.factor.poly__GT_factored_expression.cljs$core$IFn$_invoke$arity$3((arguments[(0)]),(arguments[(1)]),(arguments[(2)]));

break;
default:
throw (new Error(["Invalid arity: ",cljs.core.str.cljs$core$IFn$_invoke$arity$1(arguments.length)].join('')));

}
});

(sicmutils.polynomial.factor.poly__GT_factored_expression.cljs$core$IFn$_invoke$arity$2 = (function (p,vars){
return sicmutils.polynomial.factor.poly__GT_factored_expression.cljs$core$IFn$_invoke$arity$3(p,vars,cljs.core.identity);
}));

(sicmutils.polynomial.factor.poly__GT_factored_expression.cljs$core$IFn$_invoke$arity$3 = (function (p,vars,simplify){
try{var factors = cljs.core.map.cljs$core$IFn$_invoke$arity$2((function (factor){
var G__100794 = sicmutils.polynomial.__GT_expression(factor,vars);
return (simplify.cljs$core$IFn$_invoke$arity$1 ? simplify.cljs$core$IFn$_invoke$arity$1(G__100794) : simplify.call(null,G__100794));
}),sicmutils.polynomial.factor.split_polynomial(p));
var G__100798 = sicmutils.polynomial.factor.factors__GT_expression(factors);
return (sicmutils.polynomial.factor.simplify_product.cljs$core$IFn$_invoke$arity$1 ? sicmutils.polynomial.factor.simplify_product.cljs$core$IFn$_invoke$arity$1(G__100798) : sicmutils.polynomial.factor.simplify_product.call(null,G__100798));
}catch (e100785){if((e100785 instanceof Error)){
var _ = e100785;
taoensso.timbre._log_BANG_.cljs$core$IFn$_invoke$arity$11(taoensso.timbre._STAR_config_STAR_,new cljs.core.Keyword(null,"warn","warn",-436710552),"sicmutils.polynomial.factor",null,120,new cljs.core.Keyword(null,"p","p",151049309),new cljs.core.Keyword(null,"auto","auto",-566279492),(new cljs.core.Delay((function (){
return new cljs.core.PersistentVector(null, 1, 5, cljs.core.PersistentVector.EMPTY_NODE, ["Factorization choked! Simplifying the unfactored polynomial."], null);
}),null)),null,-1593612079,null);

var G__100787 = sicmutils.polynomial.__GT_expression(p,vars);
return (sicmutils.polynomial.factor.simplify_product.cljs$core$IFn$_invoke$arity$1 ? sicmutils.polynomial.factor.simplify_product.cljs$core$IFn$_invoke$arity$1(G__100787) : sicmutils.polynomial.factor.simplify_product.call(null,G__100787));
} else {
throw e100785;

}
}}));

(sicmutils.polynomial.factor.poly__GT_factored_expression.cljs$lang$maxFixedArity = 3);

/**
 * Given some symbolic expression containing only polynomial operations, returns a
 *   factored version of the expression with basic simplifications applied.
 * 
 *   Optionally accepts a `simplify` function that will be called on each factor of
 *   exponent 0, 1, 2 etc. Defaults to `identity`.
 * 
 *   NOTE prefer [[factor]], as [[factor]] can handle expressions with
 *   non-polynomial operations. The trigonometric functions, for example.
 */
sicmutils.polynomial.factor.factor_expression = (function sicmutils$polynomial$factor$factor_expression(var_args){
var G__100803 = arguments.length;
switch (G__100803) {
case 1:
return sicmutils.polynomial.factor.factor_expression.cljs$core$IFn$_invoke$arity$1((arguments[(0)]));

break;
case 2:
return sicmutils.polynomial.factor.factor_expression.cljs$core$IFn$_invoke$arity$2((arguments[(0)]),(arguments[(1)]));

break;
default:
throw (new Error(["Invalid arity: ",cljs.core.str.cljs$core$IFn$_invoke$arity$1(arguments.length)].join('')));

}
});

(sicmutils.polynomial.factor.factor_expression.cljs$core$IFn$_invoke$arity$1 = (function (expr){
return sicmutils.polynomial.factor.factor_expression.cljs$core$IFn$_invoke$arity$2(expr,cljs.core.identity);
}));

(sicmutils.polynomial.factor.factor_expression.cljs$core$IFn$_invoke$arity$2 = (function (expr,simplify){
var unwrapped = sicmutils.expression.expression_of(expr);
var cont = (function (p1__100799_SHARP_,p2__100800_SHARP_){
return sicmutils.polynomial.factor.poly__GT_factored_expression.cljs$core$IFn$_invoke$arity$3(p1__100799_SHARP_,p2__100800_SHARP_,simplify);
});
return sicmutils.polynomial.expression__GT_.cljs$core$IFn$_invoke$arity$2(unwrapped,cont);
}));

(sicmutils.polynomial.factor.factor_expression.cljs$lang$maxFixedArity = 2);

/**
 * Expression analyzer, identical to [[polynomial/analyzer]] except
 *   the symbolic expressions returned are in factored form.
 */
sicmutils.polynomial.factor.analyzer = (function (){var symgen = sicmutils.expression.analyze.monotonic_symbol_generator("-f-");
return sicmutils.expression.analyze.make_analyzer.cljs$core$IFn$_invoke$arity$2((function (){
if((typeof sicmutils !== 'undefined') && (typeof sicmutils.polynomial !== 'undefined') && (typeof sicmutils.polynomial.factor !== 'undefined') && (typeof sicmutils.polynomial.factor.t_sicmutils$polynomial$factor100805 !== 'undefined')){
} else {

/**
* @constructor
 * @implements {sicmutils.expression.analyze.ICanonicalize}
 * @implements {cljs.core.IMeta}
 * @implements {cljs.core.IWithMeta}
*/
sicmutils.polynomial.factor.t_sicmutils$polynomial$factor100805 = (function (symgen,meta100806){
this.symgen = symgen;
this.meta100806 = meta100806;
this.cljs$lang$protocol_mask$partition0$ = 393216;
this.cljs$lang$protocol_mask$partition1$ = 0;
});
(sicmutils.polynomial.factor.t_sicmutils$polynomial$factor100805.prototype.cljs$core$IWithMeta$_with_meta$arity$2 = (function (_100807,meta100806__$1){
var self__ = this;
var _100807__$1 = this;
return (new sicmutils.polynomial.factor.t_sicmutils$polynomial$factor100805(self__.symgen,meta100806__$1));
}));

(sicmutils.polynomial.factor.t_sicmutils$polynomial$factor100805.prototype.cljs$core$IMeta$_meta$arity$1 = (function (_100807){
var self__ = this;
var _100807__$1 = this;
return self__.meta100806;
}));

(sicmutils.polynomial.factor.t_sicmutils$polynomial$factor100805.prototype.sicmutils$expression$analyze$ICanonicalize$ = cljs.core.PROTOCOL_SENTINEL);

(sicmutils.polynomial.factor.t_sicmutils$polynomial$factor100805.prototype.sicmutils$expression$analyze$ICanonicalize$expression__GT_$arity$4 = (function (_,expr,cont,v_compare){
var self__ = this;
var ___$1 = this;
return sicmutils.polynomial.expression__GT_.cljs$core$IFn$_invoke$arity$3(expr,cont,v_compare);
}));

(sicmutils.polynomial.factor.t_sicmutils$polynomial$factor100805.prototype.sicmutils$expression$analyze$ICanonicalize$__GT_expression$arity$3 = (function (_,p,vars){
var self__ = this;
var ___$1 = this;
return sicmutils.polynomial.factor.poly__GT_factored_expression.cljs$core$IFn$_invoke$arity$2(p,vars);
}));

(sicmutils.polynomial.factor.t_sicmutils$polynomial$factor100805.prototype.sicmutils$expression$analyze$ICanonicalize$known_operation_QMARK_$arity$2 = (function (_,o){
var self__ = this;
var ___$1 = this;
return sicmutils.polynomial.analyzer.sicmutils$expression$analyze$ICanonicalize$known_operation_QMARK_$arity$2(null,o);
}));

(sicmutils.polynomial.factor.t_sicmutils$polynomial$factor100805.getBasis = (function (){
return new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [new cljs.core.Symbol(null,"symgen","symgen",-880919559,null),new cljs.core.Symbol(null,"meta100806","meta100806",2051147628,null)], null);
}));

(sicmutils.polynomial.factor.t_sicmutils$polynomial$factor100805.cljs$lang$type = true);

(sicmutils.polynomial.factor.t_sicmutils$polynomial$factor100805.cljs$lang$ctorStr = "sicmutils.polynomial.factor/t_sicmutils$polynomial$factor100805");

(sicmutils.polynomial.factor.t_sicmutils$polynomial$factor100805.cljs$lang$ctorPrWriter = (function (this__4491__auto__,writer__4492__auto__,opt__4493__auto__){
return cljs.core._write(writer__4492__auto__,"sicmutils.polynomial.factor/t_sicmutils$polynomial$factor100805");
}));

/**
 * Positional factory function for sicmutils.polynomial.factor/t_sicmutils$polynomial$factor100805.
 */
sicmutils.polynomial.factor.__GT_t_sicmutils$polynomial$factor100805 = (function sicmutils$polynomial$factor$__GT_t_sicmutils$polynomial$factor100805(symgen__$1,meta100806){
return (new sicmutils.polynomial.factor.t_sicmutils$polynomial$factor100805(symgen__$1,meta100806));
});

}

return (new sicmutils.polynomial.factor.t_sicmutils$polynomial$factor100805(symgen,cljs.core.PersistentArrayMap.EMPTY));
})()
,symgen);
})();
/**
 * Accepts a single symbolic expression and returns a factored version
 *  of that expression.
 * 
 *  Differs from [[factor-expression]] in that it can handle any expression, not
 *  just expressions limited to polynomial operations.
 */
sicmutils.polynomial.factor.factor = sicmutils.expression.analyze.default_simplifier(sicmutils.polynomial.factor.analyzer);
/**
 * Given an unwrapped symbolic expression of the form `(sqrt x)`, returns a new,
 *   unsimplified symbolic expression with any even power removed from underneath
 *   the square root.
 * 
 *   For example:
 * 
 *   ```clojure
 *   (process-sqrt
 *  '(sqrt (* x (expt y 2) (expt z 4))))
 *   ;;=> (* (sqrt x) y (expt z 2))
 *   ```
 */
sicmutils.polynomial.factor.process_sqrt = (function sicmutils$polynomial$factor$process_sqrt(expr){
var fact_exp = (function (){var G__100814 = cljs.core.first((sicmutils.numsymb.operands.cljs$core$IFn$_invoke$arity$1 ? sicmutils.numsymb.operands.cljs$core$IFn$_invoke$arity$1(expr) : sicmutils.numsymb.operands.call(null,expr)));
return (sicmutils.polynomial.factor.factor.cljs$core$IFn$_invoke$arity$1 ? sicmutils.polynomial.factor.factor.cljs$core$IFn$_invoke$arity$1(G__100814) : sicmutils.polynomial.factor.factor.call(null,G__100814));
})();
var expt = sicmutils.numsymb.symbolic_operator(new cljs.core.Symbol(null,"expt","expt",2123850345,null));
var _STAR_ = sicmutils.numsymb.symbolic_operator(new cljs.core.Symbol(null,"*","*",345799209,null));
var sqrt = sicmutils.numsymb.symbolic_operator(new cljs.core.Symbol(null,"sqrt","sqrt",370479598,null));
var even_QMARK_ = (function (n){
return ((sicmutils.value.native_integral_QMARK_(n)) && (cljs.core.even_QMARK_(n)));
});
var factors = (cljs.core.truth_(sicmutils.numsymb.product_QMARK_(fact_exp))?(sicmutils.numsymb.operands.cljs$core$IFn$_invoke$arity$1 ? sicmutils.numsymb.operands.cljs$core$IFn$_invoke$arity$1(fact_exp) : sicmutils.numsymb.operands.call(null,fact_exp)):new cljs.core.PersistentVector(null, 1, 5, cljs.core.PersistentVector.EMPTY_NODE, [fact_exp], null));
var odds = (1);
var evens = (1);
while(true){
if(cljs.core.empty_QMARK_(factors)){
if(cljs.core._EQ_.cljs$core$IFn$_invoke$arity$2(evens,(1))){
} else {
sicmutils.util.logic.assume_BANG_.cljs$core$IFn$_invoke$arity$2(cljs.core.sequence.cljs$core$IFn$_invoke$arity$1(cljs.core.seq(cljs.core.concat.cljs$core$IFn$_invoke$arity$2((new cljs.core.List(null,new cljs.core.Symbol(null,"non-negative?","non-negative?",1248641804,null),null,(1),null)),(new cljs.core.List(null,evens,null,(1),null))))),new cljs.core.Symbol(null,"root-out-squares","root-out-squares",-278750312,null));
}

var G__100841 = (sqrt.cljs$core$IFn$_invoke$arity$1 ? sqrt.cljs$core$IFn$_invoke$arity$1(odds) : sqrt.call(null,odds));
var G__100842 = evens;
return (_STAR_.cljs$core$IFn$_invoke$arity$2 ? _STAR_.cljs$core$IFn$_invoke$arity$2(G__100841,G__100842) : _STAR_.call(null,G__100841,G__100842));
} else {
var vec__100843 = factors;
var seq__100844 = cljs.core.seq(vec__100843);
var first__100845 = cljs.core.first(seq__100844);
var seq__100844__$1 = cljs.core.next(seq__100844);
var f = first__100845;
var more = seq__100844__$1;
if(cljs.core.truth_(sicmutils.numsymb.expt_QMARK_(f))){
var vec__100846 = (sicmutils.numsymb.operands.cljs$core$IFn$_invoke$arity$1 ? sicmutils.numsymb.operands.cljs$core$IFn$_invoke$arity$1(f) : sicmutils.numsymb.operands.call(null,f));
var b = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__100846,(0),null);
var e = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__100846,(1),null);
if((!(even_QMARK_(e)))){
var G__100873 = more;
var G__100874 = (_STAR_.cljs$core$IFn$_invoke$arity$2 ? _STAR_.cljs$core$IFn$_invoke$arity$2(f,odds) : _STAR_.call(null,f,odds));
var G__100875 = evens;
factors = G__100873;
odds = G__100874;
evens = G__100875;
continue;
} else {
var G__100876 = more;
var G__100877 = odds;
var G__100878 = (function (){var power = cljs.core.quot(e,(2));
if((power > (1))){
var G__100849 = evens;
var G__100850 = (expt.cljs$core$IFn$_invoke$arity$2 ? expt.cljs$core$IFn$_invoke$arity$2(b,power) : expt.call(null,b,power));
return (_STAR_.cljs$core$IFn$_invoke$arity$2 ? _STAR_.cljs$core$IFn$_invoke$arity$2(G__100849,G__100850) : _STAR_.call(null,G__100849,G__100850));
} else {
if(cljs.core._EQ_.cljs$core$IFn$_invoke$arity$2(power,(1))){
return (_STAR_.cljs$core$IFn$_invoke$arity$2 ? _STAR_.cljs$core$IFn$_invoke$arity$2(evens,b) : _STAR_.call(null,evens,b));
} else {
return evens;

}
}
})();
factors = G__100876;
odds = G__100877;
evens = G__100878;
continue;
}
} else {
var G__100884 = more;
var G__100885 = (_STAR_.cljs$core$IFn$_invoke$arity$2 ? _STAR_.cljs$core$IFn$_invoke$arity$2(f,odds) : _STAR_.call(null,f,odds));
var G__100886 = evens;
factors = G__100884;
odds = G__100885;
evens = G__100886;
continue;
}
}
break;
}
});
/**
 * Given an unwrapped symbolic expression, returns a new symbolic expression with
 *   any perfect square (exponent with an even power) removed from underneath any
 *   `sqrt` that appears in the expression.
 * 
 *   To use [[root-out-squares]] with a wrapped symbolic expression,
 *   use [[sicmutils.expression/fmap]].
 */
sicmutils.polynomial.factor.root_out_squares = (function sicmutils$polynomial$factor$root_out_squares(expr){
return clojure.walk.prewalk((function (t){
if(cljs.core.truth_(sicmutils.numsymb.sqrt_QMARK_(t))){
return sicmutils.polynomial.factor.process_sqrt(t);
} else {
return t;
}
}),expr);
});

//# sourceMappingURL=sicmutils.polynomial.factor.js.map
