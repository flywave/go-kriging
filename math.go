package kriging

import (
	"math"

	vec3d "github.com/flywave/go3d/float64/vec3"
)

func degToRad(angle float64) float64 {
	return angle * math.Pi / 180
}

func minFloat64(t []vec3d.T, k int) float64 {
	if len(t) == 0 {
		return 0
	}
	min := t[0][k]
	for i := 1; i < len(t); i++ {
		if t[i][k] < min {
			min = t[i][k]
		}
	}

	return min
}

func maxFloat64(t []vec3d.T, k int) float64 {
	if len(t) == 0 {
		return 0
	}
	max := t[0][k]
	for i := 1; i < len(t); i++ {
		if t[i][k] > max {
			max = t[i][k]
		}
	}

	return max
}
