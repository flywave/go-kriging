package kriging

import (
	"math"
)

func degToRad(angle float64) float64 {
	return angle * math.Pi / 180
}

func radToDeg(angle float64) float64 {
	return angle * 180 / math.Pi
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
