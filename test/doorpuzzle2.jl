using Spot

f = ltl"G room & G!w1 & G!w2 & G!w3 & G!w4 & (!d5 U k5) & (!d4 U k4)"
#f = ltl"G! w"
a = translate(LTLTranslator(), f)
