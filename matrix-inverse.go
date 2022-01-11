package kriging

import (
	"gonum.org/v1/gonum/mat"
)

func matrixInverse(x []float64, n int) ([]float64, bool) {
	a := mat.NewDense(n, n, x)
	var ia mat.Dense

	err := ia.Inverse(a)
	if err != nil {
		return ia.RawMatrix().Data, false
	}

	return ia.RawMatrix().Data, true
}
