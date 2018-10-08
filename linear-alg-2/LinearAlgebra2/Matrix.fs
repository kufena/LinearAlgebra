namespace LinearAlgebra2

open LinearAlgebra2.Vector

// A module for matrices, as a kind of DSL for linear algebra.
// A type of matrix can be in row or column vector form.
// There's a cross product for vectors, and other functions on matrices.
// Also some bits to go from row to column form.
//
// I guess the idea here is not to target performance but clarity.
// Although getting some GPU stuff in there might be interesting in the future.
module public Matrix =
  
  type public MatrixStyle = RowVectors | ColVectors

  type public MatrixT = MatrixT of ( MatrixStyle *
                                   int *
                                   int *
                                   (VectorT) list
                                 )

  let StyleSwitch x : MatrixStyle = match x with 
                                    |   RowVectors -> ColVectors
                                    |   ColVectors -> RowVectors

  let public crossProd v1 v2 =
    let subprod l1 x = List.map (fun v -> x * v) l1
    let mkrow n l = VectorT (n,l)
    match v1 with
      | VectorT (n1, c1) ->
        match v2 with
          | VectorT (n2, c2) ->
            if (n1 = n2) 
            then Some (MatrixT (RowVectors, n1, n2, List.map (subprod c2) c1 |> (List.map (mkrow n1)) ))
            else None

  // Take a list of lists to a zip-list type scenario where
  // the head of each list make a list, the second element of
  // each list makes a list, and so on.
  // This is a partial function - probably ought to return an option.
  let rec zippo li =
    match li with
    | [x] -> List.map (fun c -> [c]) x
    | (a :: x) -> List.zip a (zippo x) |> List.map (fun (h,t) -> h::t)
  
  // Take a list of vectors to a list of lists.
  // This is a partial function - probably ought to return an option.
  let rec toListofList li =
    match li with
    | [(VectorT (n, vcomps))] -> [vcomps]
    | ((VectorT (n, vcomps)) :: xs) -> vcomps :: (toListofList xs)
  
  // Make n-sized vectors from a list of lists.
  let revectorize n lili = List.map (fun li -> VectorT (n,li)) lili

  // Takes a matrix and returns an equivalent in row vector style.
  let public toRow mat =
    match mat with
      | (MatrixT (style, m, n, vecs)) ->
        match style with
          | RowVectors -> mat
          | ColVectors -> toListofList vecs |> zippo |> (revectorize m) |> (fun nvecs -> MatrixT (RowVectors, n, m, nvecs))

  // Takes a matrix and returns an equivalent in column matrix style.
  let public toCol mat =
    match mat with
      | (MatrixT (style, m, n, vecs)) ->
        match style with
          | ColVectors -> mat
          | RowVectors -> toListofList vecs |> zippo |> (revectorize m) |> (fun nvecs -> MatrixT (ColVectors, n, m, nvecs))
