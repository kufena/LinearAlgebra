namespace LinearAlgebra2

// A module representing real vectors.
module public Vector =

  // A type of vectors - a size and a list.
  // Of course, there's no link here between the size of the vector and
  // the length of the list.  I guess the size is separate so we don't have
  // an O(1) way of determing its size rather than using List.length.
  //
  // In matrix, we use row or column vectors.  Maybe we should be putting that
  // information in the vector, not the matrix?
  type public VectorT = VectorT of (int * ((float) list))

  // dot product of two vectors to get a scalar.
  let public dotProd x y =
    let mf (a,b) = a * b
    match x with
    | (VectorT (n1,c1)) ->
                      match y with 
                      | (VectorT (n2,c2)) -> if (n1 = n2) then (Some (List.zip c1 c2 |> (List.map mf) |> List.sum ) )
                                                          else None

  // scalar multiplication of a vector.
  let public scalarMult x y =
    let xmult a = x * a
    match y with
    | (VectorT (n,c)) -> VectorT (n, List.map xmult c)

  let toVec l = VectorT (List.length l, l)