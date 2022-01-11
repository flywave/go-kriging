package kriging

import (
	"math"

	vec3d "github.com/flywave/go3d/float64/vec3"
)

func minFloat64(t []vec3d.T, j int) float64 {
	min := float64(0)
	for i := 0; i < len(t); i++ {
		if min == 0 || min > t[i][j] {
			min = t[i][j]
		}
	}

	return min
}

func maxFloat64(t []vec3d.T, j int) float64 {
	max := float64(0)
	for i := 0; i < len(t); i++ {
		if max < t[i][j] {
			max = t[i][j]
		}
	}

	return max
}

func pipFloat64(t []Point, x, y float64) bool {
	c := false
	for i, j := 0, len(t)-1; i < len(t); j, i = i, i+1 {
		if ((t[i][1] > y) != (t[j][1] > y)) && (x < (t[j][0]-t[i][0])*(y-t[i][1])/(t[j][1]-t[i][1])+t[i][0]) {
			c = !c
		}
	}

	return c
}

func exp(x float64) float64 {
	if x == 0 {
		return 1
	}
	return math.Exp(x)
}

func pow2(x float64) float64 {
	return x * x
}

func pow3(x float64) float64 {
	return x * x * x
}
