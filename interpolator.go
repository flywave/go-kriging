package kriging

import (
	"errors"
	"image"
	"math"

	"github.com/flywave/go-cog"
	"github.com/flywave/go-geom/general"

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

const (
	BILINEAR   = "bilinear"
	HYPERBOLIC = "hyperbolic"
)

var epsg4326 geo.Proj

func init() {
	epsg4326 = geo.NewProj(4326)
}

type Interpolator interface {
	Interpolate(southWestHeight, southEastHeight, northWestHeight, northEastHeight, x, y float64) float64
}

type BilinearInterpolator struct {
	Interpolator
}

func Lerp(value1, value2, amount float64) float64 { return value1 + (value2-value1)*amount }

func (i *BilinearInterpolator) Interpolate(southWestHeight, southEastHeight, northWestHeight, northEastHeight, x, y float64) float64 {
	sw := southWestHeight
	se := southEastHeight
	nw := northWestHeight
	ne := northEastHeight

	hi_linear := Lerp(Lerp(nw, sw, y), Lerp(ne, se, y), x)

	return hi_linear
}

type HyperbolicInterpolator struct {
	Interpolator
}

func (i *HyperbolicInterpolator) Interpolate(southWestHeight, southEastHeight, northWestHeight, northEastHeight, x, y float64) float64 {
	h1 := southWestHeight
	h2 := southEastHeight
	h3 := northWestHeight
	h4 := northEastHeight
	a00 := h1
	a10 := h2 - h1
	a01 := h3 - h1
	a11 := h1 - h2 - h3 + h4
	hi_hyperbolic := a00 + a10*x + a01*y + a11*x*y
	return hi_hyperbolic
}

type KrigingInterpolator struct {
	heightModel  geoid.VerticalDatum
	heightOffset float64
	pixelSize    *[2]float64
	filterSize   [3]uint32
	inputProj    geo.Proj
	input        *geom.FeatureCollection
	inputPos     []vec3d.T
	model        ModelType
	nodata       string
	convexHull   *Convex
	kriging      *Kriging
	bounds       vec2d.Rect
	output       string
	background   *cog.Reader
	interpolator string
}

type Options struct {
	HeightModel  geoid.VerticalDatum
	HeightOffset float64
	PixelSize    *[2]float64
	InputSrs     *string
	Input        *geom.FeatureCollection
	Output       string
	Background   *string
	Model        *ModelType
	Interpolator *string
	FilterSize   *[3]uint32
}

func NewKrigingInterpolator(opts Options) *KrigingInterpolator {
	inter := &KrigingInterpolator{
		input:        opts.Input,
		heightModel:  opts.HeightModel,
		heightOffset: opts.HeightOffset,
		pixelSize:    opts.PixelSize,
		output:       opts.Output,
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

	if opts.Interpolator == nil {
		inter.interpolator = BILINEAR
	} else {
		inter.interpolator = *opts.Interpolator
	}

	if opts.Background != nil {
		inter.background = cog.Read(*opts.Background)
	}

	if opts.FilterSize != nil {
		inter.filterSize = *opts.FilterSize
	} else {
		inter.filterSize = [3]uint32{1024, 1024, 512}
	}

	return inter
}

func (p *KrigingInterpolator) extractPosion() []vec3d.T {
	ret := make([]vec3d.T, 0, 1000)

	for _, feas := range p.input.Features {
		g_ := feas.Geometry
		switch g := g_.(type) {
		case *general.Point:
			if p.inputProj != nil && !p.inputProj.Eq(epsg4326) {
				pos2 := []vec2d.T{{g.X(), g.Y()}}
				pos2 = p.inputProj.TransformTo(epsg4326, pos2)
				ret = append(ret, vec3d.T{pos2[0][0], pos2[0][1], g.Data()[2]})
			} else {
				ret = append(ret, vec3d.T{g.X(), g.Y(), g.Data()[2]})
			}
		case *general.MultiPoint:
			for _, pos := range g.Points() {
				if p.inputProj != nil && !p.inputProj.Eq(epsg4326) {
					pos2 := []vec2d.T{{pos.X(), pos.Y()}}
					pos2 = p.inputProj.TransformTo(epsg4326, pos2)
					ret = append(ret, vec3d.T{pos2[0][0], pos2[0][1], pos.Data()[2]})
				} else {
					ret = append(ret, vec3d.T{pos.X(), pos.Y(), pos.Data()[2]})
				}
			}
		case *general.LineString:
			for _, pos := range g.Subpoints() {
				if p.inputProj != nil && !p.inputProj.Eq(epsg4326) {
					pos2 := []vec2d.T{{pos.X(), pos.Y()}}
					pos2 = p.inputProj.TransformTo(epsg4326, pos2)
					ret = append(ret, vec3d.T{pos2[0][0], pos2[0][1], pos.Data()[2]})
				} else {
					ret = append(ret, vec3d.T{pos.X(), pos.Y(), pos.Data()[2]})
				}
			}
		case *general.MultiLine:
			for _, li := range g.Lines() {
				for _, pos := range li.Subpoints() {
					if p.inputProj != nil && !p.inputProj.Eq(epsg4326) {
						pos2 := []vec2d.T{{pos.X(), pos.Y()}}
						pos2 = p.inputProj.TransformTo(epsg4326, pos2)
						ret = append(ret, vec3d.T{pos2[0][0], pos2[0][1], pos.Data()[2]})
					} else {
						ret = append(ret, vec3d.T{pos.X(), pos.Y(), pos.Data()[2]})
					}
				}
			}
		case *general.Polygon:
			for _, sli := range g.Sublines() {
				for _, pos := range sli.Subpoints() {
					if p.inputProj != nil && !p.inputProj.Eq(epsg4326) {
						pos2 := []vec2d.T{{pos.X(), pos.Y()}}
						pos2 = p.inputProj.TransformTo(epsg4326, pos2)
						ret = append(ret, vec3d.T{pos2[0][0], pos2[0][1], pos.Data()[2]})
					} else {
						ret = append(ret, vec3d.T{pos.X(), pos.Y(), pos.Data()[2]})
					}
				}
			}
		case *general.MultiPolygon:
			for _, poly := range g.Polygons() {
				for _, sli := range poly.Sublines() {
					for _, pos := range sli.Subpoints() {
						if p.inputProj != nil && !p.inputProj.Eq(epsg4326) {
							pos2 := []vec2d.T{{pos.X(), pos.Y()}}
							pos2 = p.inputProj.TransformTo(epsg4326, pos2)
							ret = append(ret, vec3d.T{pos2[0][0], pos2[0][1], pos.Data()[2]})
						} else {
							ret = append(ret, vec3d.T{pos.X(), pos.Y(), pos.Data()[2]})
						}
					}
				}
			}
		}
	}
	return ret
}

func (p *KrigingInterpolator) filter(inputPos []vec3d.T) ([]vec3d.T, error) {
	min, max, _ := MinMaxVec3(inputPos)

	vg := NewVoxelGrid(vec3d.T{(max[0] - min[0]) / float64(p.filterSize[0]), (max[1] - min[1]) / float64(p.filterSize[1]), (max[2] - min[2]) / float64(p.filterSize[2])})

	res, err := vg.Filter(inputPos)

	if err != nil {
		return nil, err
	}
	return res, nil
}

func (p *KrigingInterpolator) Process() error {
	pos := p.extractPosion()

	pos, err := p.filter(pos)

	if err != nil {
		return err
	}

	p.inputPos = pos

	p.convertHeight()
	p.computeConvexHull()
	p.computeKriging()

	grid := p.cacleGrid()

	if grid == nil {
		return errors.New("gen grid error")
	}

	p.resample(grid)

	tiledata, si, bbox, srs := grid.GetDate()

	rect := image.Rect(0, 0, int(si[0]), int(si[1]))

	src := cog.NewSource(tiledata, &rect, cog.CTLZW)

	return cog.WriteTile(p.output, src, bbox, srs, si, nil)
}

func (p *KrigingInterpolator) computeConvexHull() []vec2d.T {
	p.convexHull = NewConvex(p.inputPos)
	return p.convexHull.Hull()
}

func (p *KrigingInterpolator) computeKriging() error {
	p.kriging = New(p.inputPos)
	_, err := p.kriging.Train(p.model, 0, 100)
	return err
}

func (p *KrigingInterpolator) cacleGrid() *Grid {
	if p.convexHull == nil {
		return nil
	}
	var width, height int
	if p.background != nil {
		ps := p.background.GetPixelSize(0)
		p.pixelSize = &ps
		p.bounds = p.background.GetBounds(0)
		si := p.background.GetSize(0)
		width, height = int(si[0]), int(si[1])
		epsgcode, err := p.background.GetEPSGCode(0)
		if err != nil {
			return nil
		}
		if epsgcode != 4326 {
			proj := geo.NewProj(epsgcode)
			p.bounds = proj.TransformRectTo(epsg4326, p.bounds, 16)
		}
	} else {
		if p.pixelSize == nil {
			conf := geo.DefaultTileGridOptions()
			conf[geo.TILEGRID_SRS] = epsg4326
			conf[geo.TILEGRID_RES_FACTOR] = 2.0
			conf[geo.TILEGRID_TILE_SIZE] = []uint32{512, 512}
			conf[geo.TILEGRID_ORIGIN] = geo.ORIGIN_UL

			grid := geo.NewTileGrid(conf)

			ps := [2]float64{grid.Resolutions[13], grid.Resolutions[16]}
			p.pixelSize = &ps
		}
		p.bounds = p.convexHull.Rect()

		width, height = int((p.bounds.Max[0]-p.bounds.Min[0])/p.pixelSize[0]), int((p.bounds.Max[1]-p.bounds.Min[1])/p.pixelSize[1])
	}
	grid := CaclulateGrid(width, height, geo.NewGeoReference(p.bounds, epsg4326))
	return grid
}

func (p *KrigingInterpolator) resample(grid *Grid) error {
	if p.background == nil {
		for i := range grid.Coordinates {
			if p.convexHull.InHull(vec3d.Zero, ZERO(), vec2d.T{grid.Coordinates[i][0], grid.Coordinates[i][1]}) {
				grid.Coordinates[i][2] = p.kriging.Predict(grid.Coordinates[i][0], grid.Coordinates[i][1])
			} else {
				grid.Coordinates[i][2] = default_no_data
			}
		}
	} else {
		var interpolator Interpolator

		if p.interpolator == HYPERBOLIC {
			interpolator = &HyperbolicInterpolator{}
		} else {
			interpolator = &BilinearInterpolator{}
		}

		georef := geo.NewGeoReference(p.bounds, epsg4326)

		for i := range grid.Coordinates {
			if p.convexHull.InHull(vec3d.Zero, ZERO(), vec2d.T{grid.Coordinates[i][0], grid.Coordinates[i][1]}) {
				grid.Coordinates[i][2] = p.kriging.Predict(grid.Coordinates[i][0], grid.Coordinates[i][1])
			} else {
				grid.Coordinates[i][2] = p.GetElevation(grid.Coordinates[i][0], grid.Coordinates[i][1], georef, interpolator)
			}
		}
	}
	return nil
}

func (p *KrigingInterpolator) convertHeight() {
	if (p.heightModel == geoid.HAE && p.heightOffset == 0) || p.heightModel == geoid.UNKNOWN {
		return
	}
	for i := range p.inputPos {
		if p.heightModel == geoid.HAE {
			p.inputPos[i][2] += p.heightOffset
		} else {
			gid := geoid.NewGeoid(p.heightModel, false)
			p.inputPos[i][2] = gid.ConvertHeight(p.inputPos[i][0], p.inputPos[i][1], p.inputPos[i][2], geoid.GEOIDTOELLIPSOID)
		}
	}
}

func getAverageExceptForNoDataValue(noData, valueIfAllBad float64, values ...float64) float64 {
	withValues := []float64{}

	epsilon := math.Nextafter(1, 2) - 1

	for _, v := range values {
		if math.Abs(v-noData) > epsilon {
			withValues = append(withValues, v)
		}
	}
	if len(withValues) > 0 {
		sum := 0.0
		for _, v := range withValues {
			sum += v
		}

		return sum / float64(len(withValues))
	} else {
		return valueIfAllBad
	}
}

const (
	NO_DATA_OUT = 0
)

func (s *KrigingInterpolator) getBackgroundElevation(x, y int) float64 {
	data := s.background.Data[0].([]float64)
	si := s.background.GetSize(0)
	if x >= int(si[0]) {
		x = int(si[0] - 1)
	}
	if x < 0 {
		x = 0
	}
	if y >= int(si[1]) {
		y = int(si[1] - 1)
	}
	if y < 0 {
		y = 0
	}
	return data[y*int(si[0])+x]
}

func (s *KrigingInterpolator) GetElevation(lon, lat float64, georef *geo.GeoReference, interpolator Interpolator) float64 {
	heightValue := 0.0

	si := s.background.GetSize(0)

	var yPixel, xPixel, xInterpolationAmount, yInterpolationAmount float64

	dataEndLat := georef.GetOrigin()[1] + float64(s.pixelSize[1])*float64(si[1])

	if float64(s.pixelSize[1]) > 0 {
		yPixel = ((dataEndLat - lat) / float64(s.pixelSize[1]))
	} else {
		yPixel = (lat - dataEndLat) / float64(s.pixelSize[1])
	}
	xPixel = (lon - georef.GetOrigin()[0]) / float64(s.pixelSize[0])

	epsilon := math.Max(float64(s.pixelSize[0])/10, float64(s.pixelSize[1])/10)

	_, xInterpolationAmount = math.Modf(float64(xPixel))
	_, yInterpolationAmount = math.Modf(float64(yPixel))

	xOnDataPoint := math.Abs(xInterpolationAmount) < epsilon
	yOnDataPoint := math.Abs(yInterpolationAmount) < epsilon

	if xOnDataPoint && yOnDataPoint {
		x := int(math.Floor(xPixel))
		y := int(math.Floor(yPixel))
		heightValue = s.getBackgroundElevation(x, y)
	} else {
		xCeiling := int(math.Ceil(xPixel))
		xFloor := int(math.Floor(xPixel))
		yCeiling := int(math.Ceil(yPixel))
		yFloor := int(math.Floor(yPixel))

		northWest := s.getBackgroundElevation(xFloor, yFloor)
		northEast := s.getBackgroundElevation(xCeiling, yFloor)
		southWest := s.getBackgroundElevation(xFloor, yCeiling)
		southEast := s.getBackgroundElevation(xCeiling, yCeiling)

		avgHeight := getAverageExceptForNoDataValue(default_no_data, NO_DATA_OUT, southWest, southEast, northWest, northEast)

		if northWest == default_no_data {
			northWest = avgHeight
		}
		if northEast == default_no_data {
			northEast = avgHeight
		}
		if southWest == default_no_data {
			southWest = avgHeight
		}
		if southEast == default_no_data {
			southEast = avgHeight
		}

		heightValue = interpolator.Interpolate(southWest, southEast, northWest, northEast, xInterpolationAmount, yInterpolationAmount)
	}

	return heightValue
}
