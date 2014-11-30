module Utils

export Interval

type Interval{T<:Real}
    lower::T
    upper::T
end

Base.in{T<:Real}(x::T, interval::Interval{T}) = interval.lower <= x <= interval.upper
==(a::Interval, b::Interval) = a.lower == b.lower && a.upper == b.upper
⊆(a::Interval, b::Interval) = b.lower <= a.lower <= a.upper <= b.upper
⊈(a::Interval, b::Interval) = !(a ⊆ b)
⊊(a::Interval, b::Interval) = a ⊆ b && (b.lower < a.lower || a.upper < b.upper)

end  # module Utils
