"""
    function cleanLabels!(df::DataFrame)

  # Description
Strip braces and underscores from symbolic labels.  Can be used in formatting output 
or for plotting.
- `df` a DataFrame with the symbolic CFs stored in a column called CF.

  # Returns
- `df` DataFrame with the new labels
"""
function cleanLabels(df::DataFrame)
  subs = Dict("{" => "", "}" => "","_"=>"")
  subs2 = Dict("r1"=>"g1","r2"=>"g2")
  
  for (k,v) in subs
    df.CF .= replace.(df.CF,k => v)
  end
  for (k,v) in subs2
    df.CF .= replace.(df.CF,k => v)
  end
  return(df)
end

