package kriging

import "image/color"

type ModelType string

const (
	Gaussian    ModelType = "gaussian"
	Exponential ModelType = "exponential"
	Spherical   ModelType = "spherical"
)

type DistanceList [][2]float64

func (t DistanceList) Len() int {
	return len(t)
}

func (t DistanceList) Less(i, j int) bool {
	return t[i][0] < t[j][0]
}

func (t DistanceList) Swap(i, j int) {
	tmp := t[i]
	t[i] = t[j]
	t[j] = tmp
}

type GridMatrices struct {
	Data        [][]float64 `json:"data"`
	Width       float64     `json:"width"`
	Xlim        [2]float64  `json:"xLim"`
	Ylim        [2]float64  `json:"yLim"`
	Zlim        [2]float64  `json:"zLim"`
	NodataValue float64     `json:"nodataValue"`
}

type ContourRectangle struct {
	Contour     []float64  `json:"contour"`
	XWidth      int        `json:"xWidth"`
	YWidth      int        `json:"yWidth"`
	Xlim        [2]float64 `json:"xLim"`
	Ylim        [2]float64 `json:"yLim"`
	Zlim        [2]float64 `json:"zLim"`
	XResolution float64    `json:"xResolution"`
	YResolution float64    `json:"yResolution"`
}

type Point [2]float64 // example [103.614373, 27.00541]

type Ring []Point

type PolygonCoordinates []Ring

type PolygonGeometry struct {
	Type        string `json:"type" default:"Polygon"` // Polygon
	Coordinates []Ring `json:"coordinates,omitempty"`  // coordinates
}

type GridLevelColor struct {
	Value [2]float64 `json:"value"` // range [0, 5]
	Color color.RGBA `json:"color"` // RGBA {255, 255, 255, 255}
}

type PredictDate struct {
	X     int
	Y     int
	Value float64
}
