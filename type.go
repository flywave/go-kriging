package kriging

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
