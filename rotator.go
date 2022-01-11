package kriging

import (
	"math"

	mat2d "github.com/flywave/go3d/float64/mat2"
	vec2d "github.com/flywave/go3d/float64/vec2"
)

func DegToRad(angle float64) float64 {
	return angle * math.Pi / 180
}

func RadToDeg(angle float64) float64 {
	return angle * 180 / math.Pi
}

type Rotator struct {
	Degrees float64
}

func ZERO() Rotator {
	return Rotator{0}
}

func (r *Rotator) Add(degrees float64) {
	r.Degrees += degrees
}

func (r *Rotator) AddScaled(degrees, scale float64) {
	r.Degrees += degrees * scale
}

func (r Rotator) RotateVector(v vec2d.T) vec2d.T {
	v2 := v
	mat := r.RotationMatrix()
	mat.TransformVec2(&v2)
	return v2
}

func (r Rotator) RotationMatrix() (m mat2d.T) {
	rad := DegToRad(r.Degrees)

	c := math.Cos(rad)
	s := math.Sin(rad)

	m[0][0] = c
	m[0][1] = -s
	m[1][0] = s
	m[1][1] = c

	return m
}
