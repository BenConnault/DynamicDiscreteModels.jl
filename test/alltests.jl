using Base.Test

tests = [
		"dynamicdiscretemodel", 
        "toymodel",
        "utils",
         ]


for t in tests
    tfile = string(t, ".jl")
    println(" * $(tfile) ...")
    include(tfile)
end