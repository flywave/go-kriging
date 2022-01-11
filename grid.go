package kriging

import (
	"sort"

	vec2d "github.com/flywave/go3d/float64/vec2"
	vec3d "github.com/flywave/go3d/float64/vec3"

	"github.com/flywave/go-geo"
)

type Coordinates []vec3d.T

func (s Coordinates) Len() int {
	return len(s)
}

func (s Coordinates) Less(i, j int) bool {
	if s[i][1] == s[j][1] {
		return s[i][0] < s[j][0]
	} else {
		return s[i][1] < s[j][1]
	}
}

func (s Coordinates) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}

type Grid struct {
	Width       int
	Height      int
	Coordinates Coordinates
	Count       int
	Minimum     float64
	Maximum     float64
	box         *vec3d.Box
	srs         geo.Proj
}

func NewGrid(width, height int) *Grid {
	return &Grid{Width: width, Height: height, Count: width * height, Minimum: 15000, Maximum: -15000}
}

func caclulatePixelSize(width, height int, bbox vec2d.Rect) []float64 {
	pixelSize := []float64{0, 0}
	pixelSize[0] = (bbox.Max[0] - bbox.Min[0]) / float64(width)
	pixelSize[1] = (bbox.Max[1] - bbox.Min[1]) / float64(height)
	return pixelSize
}

func CaclulateGrid(width, height int, georef *geo.GeoReference) *Grid {
	grid := NewGrid(width, height)

	grid.Count = grid.Width * grid.Height
	grid.srs = georef.GetSrs()

	coords := make(Coordinates, 0, grid.Count)

	pixelSize := caclulatePixelSize(grid.Width, grid.Height, georef.GetBBox())

	for y := grid.Height - 1; y >= 0; y-- {
		latitude := georef.GetOrigin()[1] + (float64(pixelSize[1]) * float64(y))
		for x := 0; x < grid.Width; x++ {
			longitude := georef.GetOrigin()[0] + (float64(pixelSize[0]) * float64(x))
			coords = append(coords, vec3d.T{longitude, latitude, 0})
		}
	}

	grid.Coordinates = coords
	return grid
}

func (h *Grid) GetRect() vec2d.Rect {
	bbox := h.GetBBox()
	return vec2d.Rect{Min: vec2d.T{bbox.Min[0], bbox.Min[1]}, Max: vec2d.T{bbox.Max[0], bbox.Max[1]}}
}

func (h *Grid) GetBBox() vec3d.Box {
	if h.box == nil {
		r := vec3d.Box{Min: vec3d.MaxVal, Max: vec3d.MinVal}
		for i := range h.Coordinates {
			r.Extend(&h.Coordinates[i])
		}
		return r
	}
	return *h.box
}

func (h *Grid) GetRange() float64 {
	return h.Maximum - h.Minimum
}

func (h *Grid) Sort() {
	sort.Sort(h.Coordinates)
}

func (h *Grid) Value(row, column int) float64 {
	return h.Coordinates[row*h.Width+column][2]
}

func (h *Grid) GetTileDate() ([]float64, vec3d.Box, geo.Proj) {
	tiledata := make([]float64, h.Width*h.Height)

	row, col := h.Height, h.Width

	for x := 0; x < col; x++ {
		for y := 0; y < row; y++ {
			tiledata[y*h.Width+x] = h.Value(y, x)
		}
	}

	return tiledata, h.GetBBox(), h.srs
}
