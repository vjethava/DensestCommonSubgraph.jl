
export dcs_parallel
using DensestCommonSubgraph
using Devectorize

function reduce_dcs(sA::(Array{Int64, 1}, Float64), sB::(Array{Int64, 1}, Float64))
  if sA[2] > sB[2]
    return sA
  else
    return sB
  end
end

function dcs_parallel{T}(G::Array{SparseMatrixCSC{T, Int64}, 1}, n::Number)
  @parallel reduce_dcs for j=[1:n;]
    dcs_exhaustive_worker(G, n, j)
  end
end


function dcs_exhaustive_worker{T}(G::Array{SparseMatrixCSC{T, Int64}, 1}, num_parts::Number, curr_part::Number)
  m = length(G)
  n = size(G[1], 1)

  prev_state = zeros(Bool, n)
  prev_nodes = zeros(Int64, n)
  prev_num_nodes = 0

  binary_state = zeros(Bool, n)

  curr_state = zeros(Bool, n)
  curr_nodes = zeros(Int64, n)
  curr_num_nodes = 0

  opt_idx = 0
  opt_deg = 0.0

  curr_ne = zeros(m, 1)
  curr_deg = zeros(m, 1)

  # Which part of state space to look at
  if VERSION.minor == 3
    j_begin = ifloor(2^(n)/num_parts*(curr_part - 1)) + 1
    j_end = min(iceil(2^(n)/num_parts*(curr_part)), 2^n - 1)
  else # avoids deprecated warning!
    j_begin = floor(Integer, 2^(n)/num_parts*(curr_part - 1)) + 1
    j_end = min(ceil(Integer, 2^(n)/num_parts*(curr_part)), 2^n - 1)
  end
  # update previous state
  if j_begin > 1
    int_to_binary_bool!(binary_state, j_begin - 1, n)
    binary_to_gray!(prev_state, binary_state)
    prev_num_nodes = find_set_bits!(prev_nodes, prev_state)
    for k=1:m
      for i=1:prev_num_nodes
        for j=(i+1):prev_num_nodes
          ni = prev_nodes[i]
          nj = prev_nodes[j]
          curr_ne[k] = curr_ne[k] + G[k][ni, nj] + G[k][nj, ni]
        end
      end
      curr_deg[k] = curr_ne[k]/prev_num_nodes
    end
    opt_deg = minimum(curr_deg)
    opt_idx = j_begin - 1
    # println("ne: $(curr_ne') deg: $(curr_deg') nodes: $(prev_nodes[1:prev_num_nodes])")
  end

  # start the main loop
  for j=j_begin:j_end
    int_to_binary_bool!(binary_state, j, n)
    binary_to_gray!(curr_state, binary_state)
    curr_num_nodes = find_set_bits!(curr_nodes, curr_state)

    # Find which node has been added or deleted
    changed_idx = 1
    while (curr_state[changed_idx] == prev_state[changed_idx])
      changed_idx = changed_idx + 1
    end

    # Has the node been added or deleted?
    if curr_state[changed_idx]
      m_sign = 1
    else
      m_sign = -1
    end
    @inbounds begin
      for k=1:m
        delta_deg = 0
        if prev_num_nodes > 0
          for ll=1:prev_num_nodes
            other_node = prev_nodes[ll]
            delta_deg = delta_deg + G[k][other_node, changed_idx] + G[k][changed_idx, other_node]
          end
        end
        curr_ne[k] = curr_ne[k] + m_sign  * delta_deg
        curr_deg[k] = curr_ne[k] / curr_num_nodes
        #      println("G[$k] j: $j: $(curr_ne[k]), $(curr_deg[k])")
      end
    end
    @devec curr_min_deg = minimum(curr_deg)
    if curr_min_deg > opt_deg
      opt_deg = curr_min_deg
      opt_idx = j
    end

    for k=[1:n;]
      prev_state[k] = curr_state[k]
      prev_nodes[k] = curr_nodes[k]
    end
    prev_num_nodes = curr_num_nodes
  end
  int_to_binary_bool!(binary_state, opt_idx, n)
  binary_to_gray!(curr_state, binary_state)
  curr_num_nodes = find_set_bits!(curr_nodes, curr_state)
  opt_S = curr_nodes[1:curr_num_nodes]
  # println("$(lpad(curr_part, 3))/$(rpad(num_parts, 3)) min_deg: $(rpad(opt_deg, 25)) S: $opt_S")


  return opt_S , opt_deg
end

