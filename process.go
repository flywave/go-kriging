package kriging

import (
	vec2d "github.com/flywave/go3d/float64/vec2"
	vec3d "github.com/flywave/go3d/float64/vec3"

	"github.com/flywave/go-geo"
	"github.com/flywave/go-geoid"
	"github.com/flywave/go-geom"
)

const (
	default_no_data     = float64(-9999)
	default_no_data_str = "-9999"
)

var epsg4326 geo.Proj

func init() {
	epsg4326 = geo.NewProj(4326)
}

type Interpolation struct {
	heightModel  geoid.VerticalDatum
	heightOffset float64
	pixelSize    [2]float64
	inputProj    geo.Proj
	input        *geom.FeatureCollection
	inputPos     []vec3d.T
	model        ModelType
	nodata       string
	convexHull   *Convex
	kriging      *Kriging
	bounds       vec2d.Rect
	output       string
	background   *string
}

type Options struct {
	HeightModel  geoid.VerticalDatum
	HeightOffset float64
	PixelSize    [2]float64
	InputSrs     *string
	Input        *geom.FeatureCollection
	Output       string
	Background   *string
	Model        *ModelType
}

func NewInterpolation(opts Options) *Interpolation {
	inter := &Interpolation{
		input:        opts.Input,
		heightModel:  opts.HeightModel,
		heightOffset: opts.HeightOffset,
		pixelSize:    opts.PixelSize,
		output:       opts.Output,
		background:   opts.Background,
		nodata:       default_no_data_str,
	}

	if opts.InputSrs != nil {
		inter.inputProj = geo.NewProj(opts.InputSrs)
	}

	if opts.Model != nil {
		inter.model = *opts.Model
	} else {
		inter.model = Gaussian
	}

	return inter
}

func (p *Interpolation) extractPosion() []vec3d.T {
	ret := make([]vec3d.T, 0, 1000)
	for _, feas := range p.input.Features {
		g_ := feas.Geometry
		switch g := g_.(type) {
		case geom.Point3:
			ret = append(ret, vec3d.T{g.X(), g.Y(), g.Z()})
		case geom.MultiPoint3:
			for _, pos := range g.Points() {
				ret = append(ret, vec3d.T{pos.X(), pos.Y(), pos.Z()})
			}
		case geom.LineString3:
			for _, pos := range g.Subpoints() {
				ret = append(ret, vec3d.T{pos.X(), pos.Y(), pos.Z()})
			}
		case geom.MultiLine3:
			for _, li := range g.Lines() {
				for _, pos := range li.Subpoints() {
					ret = append(ret, vec3d.T{pos.X(), pos.Y(), pos.Z()})
				}
			}
		case geom.Polygon3:
			for _, sli := range g.Sublines() {
				for _, pos := range sli.Subpoints() {
					ret = append(ret, vec3d.T{pos.X(), pos.Y(), pos.Z()})
				}
			}
		case geom.MultiPolygon3:
			for _, poly := range g.Polygons() {
				for _, sli := range poly.Sublines() {
					for _, pos := range sli.Subpoints() {
						ret = append(ret, vec3d.T{pos.X(), pos.Y(), pos.Z()})
					}
				}
			}
		}
	}
	return ret
}

func (p *Interpolation) Process() error {
	p.inputPos = p.extractPosion()
	p.computeConvexHull()
	return nil
}

func (p *Interpolation) computeConvexHull() []vec2d.T {
	p.convexHull = NewConvex(p.inputPos)
	return p.convexHull.Hull()
}

func (p *Interpolation) computeKriging() error {
	p.kriging = New(p.inputPos)
	_, err := p.kriging.Train(p.model, 0, 100)
	return err
}

func (p *Interpolation) cacleGrid() *Grid {
	return nil
}

func (p *Interpolation) writeTiff(grid *Grid) error {
	return nil
}

func (p *Interpolation) convertHeight() {
}
