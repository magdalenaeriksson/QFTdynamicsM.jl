# perform calculation in structure constructor
struct teststruct{d,dpone}
    a::Array{Float64,d}
    b::Array{Float64,dpone}

    function teststruct(d)
        c = d + 5
        return new{d,d+1}(zeros(c),zeros(c,c))
    end
end