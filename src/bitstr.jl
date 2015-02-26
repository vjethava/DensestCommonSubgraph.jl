export find_set_bits,  gray_to_binary, binary_to_gray, find_set_bits!, binary_to_gray!, int_to_binary_bool!

function find_set_bits!{T<:Number}(S::Array{T, 1}, code::Array{Bool, 1})
  ns = 0
  n = length(code)
  @inbounds begin
    for i=1:n
      if code[i] == true
        ns = ns + 1
        S[ns] = i
        #      println("$ns -> $S")
      end
    end
  end
  #  println("\tfind_set_bits!\n\t\tcode: $code\n\t\tS: $S\n\t\tns: $ns")
  return ns
end

function find_set_bits(code::String, set_char::Char ='1')
  n = length(code)
  S = zeros(Integer, n, 1)
  c = 0
  @inbounds begin
    for i=[1:length(code)]
      if code[i] == set_char
        c = c + 1
        S[c] = i
      end
    end
  end
  V = S[1:c]
  return V
end

function binary_to_gray!(gray_code::Array{Bool, 1}, binary_code::Array{Bool, 1})
  n = length(binary_code)
  @inbounds begin
    for i=(n-1):-1:1
      if (binary_code[i] != binary_code[i+1])
        gray_code[i] = true
      else
        gray_code[i] = false
      end
    end
  end
  gray_code[n] = binary_code[n]
  #  println("\tbinary: $binary_code\n\tgray: $gray_code")
end

function int_to_binary_bool!(binary_code::Array{Bool, 1}, d::Integer, n::Integer)
  i = d::Integer # local copy
  for j=[1:n;]
    if i > 0
      binary_code[j] = i % 2
      i = i >> 1
    else
      binary_code[j] = false
    end
  end
  #  println("int_to_binary_bool! i: $d n: $n code: $binary_code")
end

function binary_to_gray(binary_code::String)
  n= length(binary_code)
  gray_code  = ""
  for i=[n:-1:2;]
    # @inbounds c_char = binary_code[(i-1):i]
    # println(n, "-> " , i, " -> ", c_char, " : ", gray_code)
    if binary_code[i-1] != binary_code[i] 
      gray_code = "1" * gray_code
    else
      gray_code =  "0" * gray_code
    end
  end
  gray_code =  string(binary_code[1]) * gray_code
  return gray_code
end

function gray_to_binary(gray_code::String)
  n = length(gray_code)
  binary_code = string(gray_code[1])
  for i=2:n
    if gray_code[i] == binary_code[i-1]
      binary_code = binary_code * "0"
    else
      binary_code = binary_code * "1"
    end
  end
  return binary_code
end
