package kriging

import (
	"testing"

	vec2d "github.com/flywave/go3d/float64/vec2"
	vec3d "github.com/flywave/go3d/float64/vec3"

	"github.com/stretchr/testify/assert"
)

func TestNewConvex(t *testing.T) {
	a := assert.New(t)

	vertices := []vec3d.T{{0, 0, 0}, {100, 0, 0}, {100, -10, 0}, {150, 100, 0}, {100, 200, 0}, {0, 210, 0}, {-50, 100, 0}, {30, 30, 0}, {75, 30, 0}}
	hull := []vec2d.T{{-50, 100}, {0, 0}, {100, -10}, {150, 100}, {100, 200}, {0, 210}}

	c := NewConvex(vertices)

	a.Equal(hull, c.Hull())
}

func TestEdge(t *testing.T) {
	a := assert.New(t)

	vertices := []vec3d.T{
		{0, 0, 0},
		{100, 0, 0},
		{0, 100, 0},
		{100, 100, 0}}

	c := NewConvex(vertices)

	edges := c.Edges()
	for i, edge := range edges {
		nextIndex := i + 1
		if len(edges) <= nextIndex {
			nextIndex = 0
		}

		nextEdge := edges[nextIndex]
		a.True(OnTheRight(Subtract2(nextEdge.End, nextEdge.Start), Subtract2(edge.End, edge.Start)))
	}
}

func TestInHull(t *testing.T) {
	a := assert.New(t)

	vertices := []vec3d.T{
		{0, 0, 0},
		{100, 0, 0},
		{0, 100, 0},
		{100, 100, 0}}

	c := NewConvex(vertices)

	a.True(c.InHull(vec3d.Zero, zRotator(), vec2d.T{50, 50}))
	a.False(c.InHull(vec3d.Zero, zRotator(), vec2d.T{50, -50}))
}

func TestSupport(t *testing.T) {
	a := assert.New(t)

	c := NewConvex(
		[]vec3d.T{
			{0, 0, 0},
			{100, 0, 0},
			{0, 100, 0},
			{100, 100, 0}})

	a.Equal(c.Support(vec2d.T{1, 1}, zRotator()), vec2d.T{100, 100})
}
