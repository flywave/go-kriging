package kriging

import (
	"math"

	vec2d "github.com/flywave/go3d/float64/vec2"
	vec3d "github.com/flywave/go3d/float64/vec3"
)

type Convex struct {
	vertices []vec3d.T
	hull     []vec2d.T
	edges    []Edge
}

type Edge struct {
	Start  vec2d.T
	End    vec2d.T
	Normal vec2d.T
}

func NewConvex(vertices []vec3d.T) *Convex {
	c := Convex{vertices, nil, nil}
	return &c
}

func (c *Convex) Rect() vec2d.Rect {
	r := vec2d.Rect{Min: vec2d.MaxVal, Max: vec2d.MinVal}
	for i := range c.hull {
		r.Extend(&c.hull[i])
	}
	return r
}

func (c *Convex) Hull() []vec2d.T {
	if c.hull == nil {
		minX, maxX := c.getExtremePoints()
		c.hull = append(c.quickHull(c.vertices, maxX, minX), c.quickHull(c.vertices, minX, maxX)...)
	}

	return c.hull
}

func (c *Convex) Edges() []Edge {
	if c.edges == nil {
		hull := c.Hull()
		for i, start := range hull {
			nextIndex := i + 1
			if len(hull) <= nextIndex {
				nextIndex = 0
			}
			end := hull[nextIndex]
			r := Rotator{90}
			normal := r.RotateVector(vec2d.Sub(&start, &end))
			normal.Normalize()
			c.edges = append(c.edges, Edge{
				start,
				end,
				normal})
		}
	}
	return c.edges
}

func (c *Convex) Support(dir vec2d.T, rot Rotator) (bestVertex vec2d.T) {
	bestProjection := -math.MaxFloat64

	for _, vertex := range c.Hull() {
		v := rot.RotateVector(vertex)
		v2 := vec2d.T{dir[0], dir[1]}
		projection := vec2d.Dot(&v, &v2)

		if bestProjection < projection {
			bestVertex = rot.RotateVector(vec2d.T{vertex[0], vertex[1]})
			bestProjection = projection
		}
	}

	return bestVertex
}

func (c *Convex) quickHull(points []vec3d.T, start, end vec2d.T) []vec2d.T {
	pointDistanceIndicators := c.getLhsPointDistanceIndicatorMap(points, start, end)
	if len(pointDistanceIndicators) == 0 {
		return []vec2d.T{end}
	}

	farthestPoint := c.getFarthestPoint(pointDistanceIndicators)

	newPoints := []vec3d.T{}
	for point := range pointDistanceIndicators {
		newPoints = append(newPoints, point)
	}

	return append(
		c.quickHull(newPoints, farthestPoint, end),
		c.quickHull(newPoints, start, farthestPoint)...)
}

func Subtract(lhs vec3d.T, rhs vec2d.T) vec2d.T {
	return vec2d.T{lhs[0] - rhs[0], lhs[1] - rhs[1]}
}

func Subtract2(lhs vec2d.T, rhs vec2d.T) vec2d.T {
	return vec2d.T{lhs[0] - rhs[0], lhs[1] - rhs[1]}
}

func Add(lhs vec3d.T, rhs vec2d.T) vec2d.T {
	return vec2d.T{lhs[0] + rhs[0], lhs[1] + rhs[1]}
}

func OnTheRight(v vec2d.T, o vec2d.T) bool {
	return Cross(v, o) < 0
}

func (c *Convex) InHull(position vec3d.T, rotation Rotator, point vec2d.T) bool {
	for _, edge := range c.Edges() {
		if !OnTheRight(Subtract2(point, Add(position, rotation.RotateVector(edge.Start))), Subtract2(Add(position, rotation.RotateVector(edge.End)), Add(position, rotation.RotateVector(edge.Start)))) {
			return false
		}
	}

	return true
}

func (c *Convex) getExtremePoints() (minX, maxX vec2d.T) {
	minX = vec2d.T{math.MaxFloat64, 0}
	maxX = vec2d.T{-math.MaxFloat64, 0}

	for _, p := range c.vertices {
		if p[0] < minX[0] {
			minX = vec2d.T{p[0], p[1]}
		}

		if maxX[0] < p[0] {
			maxX = vec2d.T{p[0], p[1]}
		}
	}

	return minX, maxX
}

func (c *Convex) getLhsPointDistanceIndicatorMap(points []vec3d.T, start, end vec2d.T) map[vec3d.T]float64 {
	pointDistanceIndicatorMap := make(map[vec3d.T]float64)

	for _, point := range points {
		distanceIndicator := c.getDistanceIndicator(point, start, end)
		if distanceIndicator > 0 {
			pointDistanceIndicatorMap[point] = distanceIndicator
		}
	}

	return pointDistanceIndicatorMap
}

func Cross(lhs, rhs vec2d.T) float64 {
	return (lhs[0] * rhs[1]) - (lhs[1] * rhs[0])
}

func (c *Convex) getDistanceIndicator(point vec3d.T, start, end vec2d.T) float64 {
	point2d := vec2d.T{point[0], point[1]}
	vLine := vec2d.Sub(&end, &start)

	vPoint := vec2d.Sub(&point2d, &start)

	return Cross(vLine, vPoint)
}

func (c *Convex) getFarthestPoint(pointDistanceIndicatorMap map[vec3d.T]float64) (farthestPoint vec2d.T) {
	maxDistanceIndicator := -math.MaxFloat64
	for point, distanceIndicator := range pointDistanceIndicatorMap {
		if maxDistanceIndicator < distanceIndicator {
			maxDistanceIndicator = distanceIndicator
			farthestPoint = vec2d.T{point[0], point[1]}
		}
	}

	return farthestPoint
}
