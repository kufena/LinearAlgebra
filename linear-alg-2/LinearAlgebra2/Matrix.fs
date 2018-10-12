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
  
  // a wrap list function which takes a list of numbers to an index and
  // a wrapped rest of list.  So [1;2;3] would go to (1,[2;3]) and (2,[3;1])
  // and (3,[1;2])
  // This is predominantly for use when working out sub-matrices when doing
  // determinants of matrices bigger than 2x2.
  let rec wraps l acc =
    match l with
    | [] -> []
    | (a :: x) -> (a, (-1.0) ** ((float) ((List.length x))), List.append x acc) :: (wraps x (List.append acc [a]))
  
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

  let public index m (p,q) =

    let (style,rows,cols,vecs) = 
      match m with
      | MatrixT (st, rw, co, vc) -> (st, rw, co, vc)

    if (p >= rows || q >= cols) then
      None
    else
      match style with
      | ColVectors -> match (List.item q vecs) with
                      | VectorT (rs, ls) -> Some (List.item p ls)
      | RowVectors -> match (List.item p vecs) with
                      | VectorT (cs, ls) -> Some (List.item q ls)

  // 2x2 determinants
  // Does the normal (ab - cd) determinant from anywhere inside a matrix.
  // Ideally, we would check to make sure we are on the bottom two rows, but
  // it's fun to see the inner determinants.

  let simple2x2determinant vecs c1 c2 row =
    let (VectorT (n1,v1)) = List.item c1 vecs
    let (VectorT (n2,v2)) = List.item c2 vecs

    let a = List.item row v1
    let b = List.item row v2
    let c = List.item (row+1) v1
    let d = List.item (row+1) v2

    (a * d) - (b * c)

  // This uses a list of wraps to find determinants of sub-matrices and multiplies
  // them by a row/col position value.  If we've got to 2x2 already (from sz) then
  // we just return the appropriate 2x2 from the precomputed list.  Otherwise, we
  // use the wraps to recursively produce determinants for sub matrices, then sum
  // them all at the end.
  let rec computedeterminant row sz vecs mywraps =

    printf "%d %d" row sz

    // find a specific 2x2 det using our column position in the wrap.
    let wrapTo2x2 q =
      match q with
      | (r, sg, l) -> simple2x2determinant vecs r (List.head l) row

    // for one wrap, compute the determinant by creating the wraps for our column
    // list, and the multiply the result using the value from the vectors.
    let handleWrap wr =
      match wr with
      | (c, sg, clist) -> sg * computedeterminant (row + 1) (sz - 1) vecs (wraps clist [])

    if (sz = 2)
    // our size is 2 - just get the 2x2 - there should only be one.
    then List.head mywraps |> wrapTo2x2
    // otherwise get a value for each wrap and combine by addition.
    else List.map handleWrap mywraps |> List.fold (+) 0.0
  

  // A more general determinant calculator.
  let public determinant m p q sz =
    let colMatrix = toCol m

    let (style,rows,cols,vecs) = 
      match colMatrix with
      | MatrixT (st, rw, co, vc) -> (st, rw, co, vc)

    if (rows <> cols) // not sqaure
    then None
    else if (cols = 2) // simple - just do it
         then Some (simple2x2determinant vecs 0 1 0)
         else Some (computedeterminant 0 cols vecs (wraps [0 .. (cols-1)] []))

  