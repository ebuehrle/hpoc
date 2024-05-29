using Spot

k1 = "(k11 & k12 & k13 & k14)"
k2 = "(k21 & k22 & k23 & k24)"
k3 = "(k31 & k32 & k33 & k34)"
k4 = "(k41 & k42 & k43 & k44)"
k5 = "(k51 & k52 & k53 & k54)"

d1 = "(d11 & d12 & d13 & d14)"
d2 = "(d21 & d22 & d23 & d24)"
d3 = "(d31 & d32 & d33 & d34)"
d4 = "(d41 & d42 & d43 & d44)"
d5 = "(d51 & d52 & d53 & d54)"

w1 = "(w11 & w12 & w13 & w14)"
w2 = "(w21 & w22 & w23 & w24)"
w3 = "(w31 & w32 & w33 & w34)"
w4 = "(w41 & w42 & w43 & w44)"
w5 = "(w51 & w52 & w53 & w54)"

room = "(r1 & r2 & r3 & r4)"

f = SpotFormula("G$(room) & G!$(w1) & G!$(w2) & G!$(w3) & G!$(w4) & (!$(d5) U $(k5)) & (!$(d4) U $(k4))")
#f = SpotFormula("G$(room) & G!$(w1) & G!$(w2) & G!$(w3) & G!$(w4) & G!$(w5) & (!$(d5) U $(k5)) & (!$(d4) U $(k4)) & (!$(d3) U $(k3)) & (!$(d2) U $(k2)) & (!$(d1) U $(k1))")
a = translate(LTLTranslator(), f)
