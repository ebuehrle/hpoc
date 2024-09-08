julia --project stlcg-1.jl > img/stlcg-1.log
mv gmp.sdpa img/stlcg-1.sdpa
julia --project stlcg-2.jl > img/stlcg-2.log
mv gmp.dat-s img/stlcg-2.dat-s
julia --project doorpuzzle-1.jl > img/doorpuzzle-1.log
mv gmp.dat-s img/doorpuzzle-1.dat-s
julia --project rover-2.jl > img/rover-2.log
mv gmp.dat-s img/rover-2.dat-s
julia --project stlcg-1-miqp.jl > img/stlcg-1-miqp.log
julia --project stlcg-2-miqp.jl > img/stlcg-2-miqp.log
