goog.provide('cljs.tools.reader.impl.commons');
/**
 * Checks whether the reader is at the start of a number literal
 */
cljs.tools.reader.impl.commons.number_literal_QMARK_ = (function cljs$tools$reader$impl$commons$number_literal_QMARK_(reader,initch){
return ((cljs.tools.reader.impl.utils.numeric_QMARK_(initch)) || (((((("+" === initch)) || (("-" === initch)))) && (cljs.tools.reader.impl.utils.numeric_QMARK_(reader.cljs$tools$reader$reader_types$Reader$peek_char$arity$1(null))))));
});
/**
 * Read until first character that doesn't match pred, returning
 * char.
 */
cljs.tools.reader.impl.commons.read_past = (function cljs$tools$reader$impl$commons$read_past(pred,rdr){
var ch = rdr.cljs$tools$reader$reader_types$Reader$read_char$arity$1(null);
while(true){
if((pred.cljs$core$IFn$_invoke$arity$1 ? pred.cljs$core$IFn$_invoke$arity$1(ch) : pred.call(null,ch))){
var G__74440 = rdr.cljs$tools$reader$reader_types$Reader$read_char$arity$1(null);
ch = G__74440;
continue;
} else {
return ch;
}
break;
}
});
/**
 * Advances the reader to the end of a line. Returns the reader
 */
cljs.tools.reader.impl.commons.skip_line = (function cljs$tools$reader$impl$commons$skip_line(reader){
while(true){
if(cljs.tools.reader.impl.utils.newline_QMARK_(reader.cljs$tools$reader$reader_types$Reader$read_char$arity$1(null))){
} else {
continue;
}
break;
}

return reader;
});
cljs.tools.reader.impl.commons.int_pattern = /^([-+]?)(?:(0)|([1-9][0-9]*)|0[xX]([0-9A-Fa-f]+)|0([0-7]+)|([1-9][0-9]?)[rR]([0-9A-Za-z]+)|0[0-9]+)(N)?$/;
cljs.tools.reader.impl.commons.ratio_pattern = /([-+]?[0-9]+)\/([0-9]+)/;
cljs.tools.reader.impl.commons.float_pattern = /([-+]?[0-9]+(\.[0-9]*)?([eE][-+]?[0-9]+)?)(M)?/;
cljs.tools.reader.impl.commons.match_int = (function cljs$tools$reader$impl$commons$match_int(s){
var m = cljs.core.vec(cljs.core.re_find(cljs.tools.reader.impl.commons.int_pattern,s));
if((!(((m.cljs$core$IFn$_invoke$arity$1 ? m.cljs$core$IFn$_invoke$arity$1((2)) : m.call(null,(2))) == null)))){
return (0);
} else {
var negate_QMARK_ = ("-" === (m.cljs$core$IFn$_invoke$arity$1 ? m.cljs$core$IFn$_invoke$arity$1((1)) : m.call(null,(1))));
var a = (((!(((m.cljs$core$IFn$_invoke$arity$1 ? m.cljs$core$IFn$_invoke$arity$1((3)) : m.call(null,(3))) == null))))?new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [(m.cljs$core$IFn$_invoke$arity$1 ? m.cljs$core$IFn$_invoke$arity$1((3)) : m.call(null,(3))),(10)], null):(((!(((m.cljs$core$IFn$_invoke$arity$1 ? m.cljs$core$IFn$_invoke$arity$1((4)) : m.call(null,(4))) == null))))?new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [(m.cljs$core$IFn$_invoke$arity$1 ? m.cljs$core$IFn$_invoke$arity$1((4)) : m.call(null,(4))),(16)], null):(((!(((m.cljs$core$IFn$_invoke$arity$1 ? m.cljs$core$IFn$_invoke$arity$1((5)) : m.call(null,(5))) == null))))?new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [(m.cljs$core$IFn$_invoke$arity$1 ? m.cljs$core$IFn$_invoke$arity$1((5)) : m.call(null,(5))),(8)], null):(((!(((m.cljs$core$IFn$_invoke$arity$1 ? m.cljs$core$IFn$_invoke$arity$1((7)) : m.call(null,(7))) == null))))?new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [(m.cljs$core$IFn$_invoke$arity$1 ? m.cljs$core$IFn$_invoke$arity$1((7)) : m.call(null,(7))),parseInt((m.cljs$core$IFn$_invoke$arity$1 ? m.cljs$core$IFn$_invoke$arity$1((6)) : m.call(null,(6))))], null):new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [null,null], null)
))));
var n = (a.cljs$core$IFn$_invoke$arity$1 ? a.cljs$core$IFn$_invoke$arity$1((0)) : a.call(null,(0)));
if((n == null)){
return null;
} else {
var bn = parseInt(n,(a.cljs$core$IFn$_invoke$arity$1 ? a.cljs$core$IFn$_invoke$arity$1((1)) : a.call(null,(1))));
var bn__$1 = ((negate_QMARK_)?((-1) * bn):bn);
if(cljs.core.truth_(isNaN(bn__$1))){
return null;
} else {
return bn__$1;
}
}
}
});
cljs.tools.reader.impl.commons.match_ratio = (function cljs$tools$reader$impl$commons$match_ratio(s){
var m = cljs.core.vec(cljs.core.re_find(cljs.tools.reader.impl.commons.ratio_pattern,s));
var numerator = (m.cljs$core$IFn$_invoke$arity$1 ? m.cljs$core$IFn$_invoke$arity$1((1)) : m.call(null,(1)));
var denominator = (m.cljs$core$IFn$_invoke$arity$1 ? m.cljs$core$IFn$_invoke$arity$1((2)) : m.call(null,(2)));
var numerator__$1 = (cljs.core.truth_(cljs.core.re_find(/^\+/,numerator))?cljs.core.subs.cljs$core$IFn$_invoke$arity$2(numerator,(1)):numerator);
return (parseInt(numerator__$1) / parseInt(denominator));
});
cljs.tools.reader.impl.commons.match_float = (function cljs$tools$reader$impl$commons$match_float(s){
var m = cljs.core.vec(cljs.core.re_find(cljs.tools.reader.impl.commons.float_pattern,s));
if((!(((m.cljs$core$IFn$_invoke$arity$1 ? m.cljs$core$IFn$_invoke$arity$1((4)) : m.call(null,(4))) == null)))){
return parseFloat((m.cljs$core$IFn$_invoke$arity$1 ? m.cljs$core$IFn$_invoke$arity$1((1)) : m.call(null,(1))));
} else {
return parseFloat(s);
}
});
cljs.tools.reader.impl.commons.matches_QMARK_ = (function cljs$tools$reader$impl$commons$matches_QMARK_(pattern__$1,s){
var vec__74435 = cljs.core.re_find(pattern__$1,s);
var match = cljs.core.nth.cljs$core$IFn$_invoke$arity$3(vec__74435,(0),null);
return (match === s);
});
cljs.tools.reader.impl.commons.match_number = (function cljs$tools$reader$impl$commons$match_number(s){
if(cljs.tools.reader.impl.commons.matches_QMARK_(cljs.tools.reader.impl.commons.int_pattern,s)){
return cljs.tools.reader.impl.commons.match_int(s);
} else {
if(cljs.tools.reader.impl.commons.matches_QMARK_(cljs.tools.reader.impl.commons.float_pattern,s)){
return cljs.tools.reader.impl.commons.match_float(s);
} else {
if(cljs.tools.reader.impl.commons.matches_QMARK_(cljs.tools.reader.impl.commons.ratio_pattern,s)){
return cljs.tools.reader.impl.commons.match_ratio(s);
} else {
return null;
}
}
}
});
/**
 * Parses a string into a vector of the namespace and symbol
 */
cljs.tools.reader.impl.commons.parse_symbol = (function cljs$tools$reader$impl$commons$parse_symbol(token){
if(((("" === token)) || (((/:$/.test(token) === true) || (/^::/.test(token) === true))))){
return null;
} else {
var ns_idx = token.indexOf("/");
var ns = (((ns_idx > (0)))?cljs.core.subs.cljs$core$IFn$_invoke$arity$3(token,(0),ns_idx):null);
if((!((ns == null)))){
var ns_idx__$1 = (ns_idx + (1));
if((ns_idx__$1 === cljs.core.count(token))){
return null;
} else {
var sym = cljs.core.subs.cljs$core$IFn$_invoke$arity$2(token,ns_idx__$1);
if((((!(cljs.tools.reader.impl.utils.numeric_QMARK_(cljs.core.nth.cljs$core$IFn$_invoke$arity$2(sym,(0)))))) && ((((!(("" === sym)))) && (((/:$/.test(ns) === false) && ((((sym === "/")) || (((-1) === sym.indexOf("/"))))))))))){
return new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [ns,sym], null);
} else {
return null;
}
}
} else {
if((((token === "/")) || (((-1) === token.indexOf("/"))))){
return new cljs.core.PersistentVector(null, 2, 5, cljs.core.PersistentVector.EMPTY_NODE, [null,token], null);
} else {
return null;
}
}
}
});
cljs.tools.reader.impl.commons.read_comment = (function cljs$tools$reader$impl$commons$read_comment(var_args){
var args__4870__auto__ = [];
var len__4864__auto___74443 = arguments.length;
var i__4865__auto___74444 = (0);
while(true){
if((i__4865__auto___74444 < len__4864__auto___74443)){
args__4870__auto__.push((arguments[i__4865__auto___74444]));

var G__74445 = (i__4865__auto___74444 + (1));
i__4865__auto___74444 = G__74445;
continue;
} else {
}
break;
}

var argseq__4871__auto__ = ((((1) < args__4870__auto__.length))?(new cljs.core.IndexedSeq(args__4870__auto__.slice((1)),(0),null)):null);
return cljs.tools.reader.impl.commons.read_comment.cljs$core$IFn$_invoke$arity$variadic((arguments[(0)]),argseq__4871__auto__);
});

(cljs.tools.reader.impl.commons.read_comment.cljs$core$IFn$_invoke$arity$variadic = (function (rdr,_){
return cljs.tools.reader.impl.commons.skip_line(rdr);
}));

(cljs.tools.reader.impl.commons.read_comment.cljs$lang$maxFixedArity = (1));

/** @this {Function} */
(cljs.tools.reader.impl.commons.read_comment.cljs$lang$applyTo = (function (seq74438){
var G__74439 = cljs.core.first(seq74438);
var seq74438__$1 = cljs.core.next(seq74438);
var self__4851__auto__ = this;
return self__4851__auto__.cljs$core$IFn$_invoke$arity$variadic(G__74439,seq74438__$1);
}));

cljs.tools.reader.impl.commons.throwing_reader = (function cljs$tools$reader$impl$commons$throwing_reader(msg){
return (function() { 
var G__74447__delegate = function (rdr,_){
return cljs.tools.reader.impl.errors.reader_error.cljs$core$IFn$_invoke$arity$variadic(rdr,cljs.core.prim_seq.cljs$core$IFn$_invoke$arity$2([msg], 0));
};
var G__74447 = function (rdr,var_args){
var _ = null;
if (arguments.length > 1) {
var G__74448__i = 0, G__74448__a = new Array(arguments.length -  1);
while (G__74448__i < G__74448__a.length) {G__74448__a[G__74448__i] = arguments[G__74448__i + 1]; ++G__74448__i;}
  _ = new cljs.core.IndexedSeq(G__74448__a,0,null);
} 
return G__74447__delegate.call(this,rdr,_);};
G__74447.cljs$lang$maxFixedArity = 1;
G__74447.cljs$lang$applyTo = (function (arglist__74449){
var rdr = cljs.core.first(arglist__74449);
var _ = cljs.core.rest(arglist__74449);
return G__74447__delegate(rdr,_);
});
G__74447.cljs$core$IFn$_invoke$arity$variadic = G__74447__delegate;
return G__74447;
})()
;
});

//# sourceMappingURL=cljs.tools.reader.impl.commons.js.map
