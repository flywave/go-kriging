package kriging

import (
	"math"

	vec3d "github.com/flywave/go3d/float64/vec3"
)

func degToRad(angle float64) float64 {
	return angle * math.Pi / 180
}

func minFloat64(t []vec3d.T, k int) float64 {
	min := float64(0)
	for i := 0; i < len(t); i++ {
		if min == 0 || min > t[i][k] {
			min = t[i][k]
		}
	}

	return min
}

func maxFloat64(t []vec3d.T, k int) float64 {
	max := float64(0)
	for i := 0; i < len(t); i++ {
		if max < t[i][k] {
			max = t[i][k]
		}
	}

	return max
}
