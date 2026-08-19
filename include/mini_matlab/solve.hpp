#pragma once

#include "mini_matlab/matrix.hpp"

namespace mini_matlab {

// Thrown when elimination reaches a pivot that is zero, or small enough
// relative to the scale of the original matrix to be indistinguishable from
// zero in double precision.
//
// This detects a *structurally* singular system. It does NOT detect
// ill-conditioning: a matrix can have every pivot equal to 1 and still be
// numerically singular (unit diagonal, -1 above, kappa ~ 2^(n-1)). Deciding
// that requires a condition estimator, which this engine does not have.
class SingularMatrix : public MatrixError {
public:
    using MatrixError::MatrixError;
};

// Solves a * x = b by Gaussian elimination with partial pivoting, returning x.
// MATLAB spells this `a \ b`.
//
// `a` must be square and n-by-n; `b` must have n rows and may have any number
// of columns (each column is an independent right-hand side).
//
// Both parameters are taken by value: the algorithm destroys them, and the
// caller's matrices -- which in the REPL are entries in the symbol table --
// must survive unchanged. Passing a temporary moves rather than copies.
//
// Throws DimensionMismatch if `a` is not square or the shapes disagree,
// SingularMatrix if no usable pivot exists.
Matrix solve(Matrix a, Matrix b);

} // namespace mini_matlab
