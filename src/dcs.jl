export dcs_exhaustive
using Devectorize


function dcs_exhaustive{T}(G::Array{SparseMatrixCSC{T, Int64}, 1})
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

  for j=[1:(2^n-1);]
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
      for k=[1:m;]
        delta_deg = 0
        if prev_num_nodes > 0
          for ll=[1:prev_num_nodes;]
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
  println("S: $opt_S min_deg: $opt_deg")
  return opt_S , opt_deg
end



