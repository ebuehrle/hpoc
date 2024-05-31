using Spot

#f = ltl"G room & G!w1 & G!w2 & G!w3 & G!w4 & (!d5 U k5) & (!d4 U k4)"
#f = ltl"G((a1 | a2) -> F(b1 | b2))"
#f = ltl"G(r1a) & G(r2a) & G((r1o1 | r1o2 | r1o3) -> F(r1t1 | r1t2)) & G((r2o1 | r2o2 | r2o3) -> F(r2t1 | r2t2)) & F(r1o1 | r2o1)"# & F(r1o2 | r2o2) & F(r1o3 | r2o3)"# & F(r1o4 | r2o4) & F(r1t1 | r2t1) & F(r1t2 | r2t2)#
#f = ltl"(G(a -> F(b))) & (G(c -> F(d)))"
#f = ltl"G((a1 | a2) -> F(b)) & G(c -> F(d))"
#f = ltl"G(r1a) & G(r2a) & G((r1o1 | r1o2 | r1o3) -> F(r1t1 | r1t2)) & G((r2o1 | r2o2 | r2o3) -> F(r2t1 | r2t2)) & F(r1o1 | r2o1) & F(r1o2 | r2o2) & F(r1o3 | r2o3)"# & F(r1o4 | r2o4) & F(r1t1 | r2t1) & F(r1t2 | r2t2)
#f = ltl"G((a1|a2) -> F(b1|b2)) & F(a1|a2)"
#f = ltl"G! w"
f = ltl"G(a -> F(b))"
a = translate(LTLTranslator(deterministic=true, state_based_acceptance=true, buchi=true), f)
get_rabin_acceptance(a)
