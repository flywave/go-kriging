package kriging

import (
	vec2d "github.com/flywave/go3d/float64/vec2"
	vec3d "github.com/flywave/go3d/float64/vec3"

	"github.com/flywave/go-cog"
	"github.com/flywave/go-geo"
	"github.com/flywave/go-geoid"
	"github.com/flywave/go-geom"
)

const (
	default_no_data     = float64(-9999)
	default_no_data_str = "-9999"
)

func (c FlatPoints) contains(x, y float64) bool {
	contains := false
	k := c.Len()
	j := k - 2
	for i := 0; i < k-1; i++ {
		xi, yi := c.Take(i)
		xj, yj := c.Take(j)

		if ((yi > y) != (yj > y)) && (x < (xj-xi)*(y-yi)/(yj-yi)+xi) {
			contains = !contains
		}
		j = i
	}
	return contains
}

var epsg4326 geo.Proj

func init() {
	epsg4326 = geo.NewProj(4326)
}

type process struct {
	heightModel  geoid.VerticalDatum
	heightOffset float64
	pixelSize    [2]float64
	inputProj    geo.Proj
	input        *geom.FeatureCollection
	inputPos     []vec3d.T
	model        KrigingModel
	nodata       string
	convexHull   FlatPoints
	bounds       vec2d.Rect
	src          cog.TileSource
	output       string
}
