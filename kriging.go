package kriging

import (
	"errors"
	"math"
	"sort"
	"sync"

	vec3d "github.com/flywave/go3d/float64/vec3"
)

type Kriging struct {
	pos []vec3d.T

	Nugget float64 `json:"nugget"`
	Range  float64 `json:"range"`
	Sill   float64 `json:"sill"`
	A      float64 `json:"A"`
	N      int     `json:"n"`

	K []float64 `json:"K"`
	M []float64 `json:"M"`

	model KrigingModel
}

func New(pos []vec3d.T) *Kriging {
	return &Kriging{pos: pos}
}

type KrigingModel func(float64, float64, float64, float64, float64) float64

func krigingKrigingGaussian(h, nugget, range_, sill, A float64) float64 {
	x := -(1.0 / A) * ((h / range_) * (h / range_))
	return nugget + ((sill-nugget)/range_)*
		(1.0-exp(x))
}

func krigingKrigingExponential(h, nugget, range_, sill, A float64) float64 {
	x := -(1.0 / A) * (h / range_)
	return nugget + ((sill-nugget)/range_)*
		(1.0-exp(x))
}

func krigingKrigingSpherical(h, nugget, range_, sill, A float64) float64 {
	if h > range_ {
		return nugget + (sill-nugget)/range_
	} else {
		x := h / range_
		return nugget + ((sill-nugget)/range_)*
			(1.5*(x)-0.5*(pow3(x)))
	}
}

func (kri *Kriging) Train(model ModelType, sigma2 float64, alpha float64) (*Kriging, error) {
	kri.Nugget = 0.0
	kri.Range = 0.0
	kri.Sill = 0.0
	kri.A = float64(1) / float64(3)
	kri.N = 0.0

	switch model {
	case Gaussian:
		kri.model = krigingKrigingGaussian
		break
	case Exponential:
		kri.model = krigingKrigingExponential
		break
	case Spherical:
		kri.model = krigingKrigingSpherical
		break
	}

	var i, j, k, l, n int
	n = len(kri.pos)

	var distance DistanceList
	distance = make([][2]float64, (n*n-n)/2)

	i = 0
	k = 0
	for ; i < n; i++ {
		for j = 0; j < i; {
			distance[k] = [2]float64{}
			distance[k][0] = math.Sqrt(pow2(kri.pos[i][0]-kri.pos[j][0]) + pow2(kri.pos[i][1]-kri.pos[j][1]))
			distance[k][1] = math.Abs(kri.pos[i][2] - kri.pos[j][2])
			j++
			k++
		}
	}
	sort.Sort(distance)
	kri.Range = distance[(n*n-n)/2-1][0]

	var lags int
	if ((n*n - n) / 2) > 30 {
		lags = 30
	} else {
		lags = (n*n - n) / 2
	}

	tolerance := kri.Range / float64(lags)

	lag := make([]float64, lags)
	semi := make([]float64, lags)
	if lags < 30 {
		for l = 0; l < lags; l++ {
			lag[l] = distance[l][0]
			semi[l] = distance[l][1]
		}
	} else {
		i = 0
		j = 0
		k = 0
		l = 0
		for i < lags && j < ((n*n-n)/2) {
			for {
				if distance[j][0] > (float64(i+1) * tolerance) {
					break
				}
				lag[l] += distance[j][0]
				semi[l] += distance[j][1]
				j++
				k++
				if j >= ((n*n - n) / 2) {
					break
				}
			}

			if k > 0 {
				lag[l] = lag[l] / float64(k)
				semi[l] = semi[l] / float64(k)
				l++
			}
			i++
			k = 0
		}
		if l < 2 {
			return nil, errors.New("not enough points")
		}
	}

	n = l
	kri.Range = lag[n-1] - lag[0]
	X := make([]float64, 2*n)
	for i := 0; i < len(X); i++ {
		X[i] = 1
	}
	Y := make([]float64, n)
	var A = kri.A
	for i = 0; i < n; i++ {
		switch model {
		case Gaussian:
			X[i*2+1] = 1.0 - exp(-(1.0/A)*pow2(lag[i]/kri.Range))
			break
		case Exponential:
			X[i*2+1] = 1.0 - exp(-(1.0/A)*lag[i]/kri.Range)
			break
		case Spherical:
			X[i*2+1] = 1.5*(lag[i]/kri.Range) - 0.5*pow3(lag[i]/kri.Range)
			break
		}
		Y[i] = semi[i]
	}

	var Xt = matrixTranspose(X, n, 2)
	var Z = matrixMultiply(Xt, X, 2, n, 2)
	Z = matrixAdd(Z, matrixDiag(float64(1)/alpha, 2), 2, 2)
	var cloneZ = make([]float64, len(Z))
	copy(cloneZ, Z)
	if matrixChol(Z, 2) {
		matrixChol2inv(Z, 2)
	} else {
		Z, _ = matrixInverse(cloneZ, 2)
	}

	var W = matrixMultiply(matrixMultiply(Z, Xt, 2, 2, n), Y, 2, n, 1)

	kri.Nugget = W[0]
	kri.Sill = W[1]*kri.Range + kri.Nugget
	kri.N = len(kri.pos)

	n = len(kri.pos)
	K := make([]float64, n*n)
	for i = 0; i < n; i++ {
		for j = 0; j < i; j++ {
			K[i*n+j] = kri.model(
				math.Sqrt(pow2(kri.pos[i][0]-kri.pos[j][0])+pow2(kri.pos[i][1]-kri.pos[j][1])),
				kri.Nugget,
				kri.Range,
				kri.Sill,
				kri.A)
			K[j*n+i] = K[i*n+j]
		}
		K[i*n+i] = kri.model(0, kri.Nugget,
			kri.Range,
			kri.Sill,
			kri.A)
	}

	var C = matrixAdd(K, matrixDiag(sigma2, n), n, n)
	var cloneC = make([]float64, len(C))
	copy(cloneC, C)
	if matrixChol(C, n) {
		matrixChol2inv(C, n)
	} else {
		// TODO false
		C, _ = matrixInverse(cloneC, n)
	}

	copy(K, C)
	t := make([]float64, n)

	for i := range kri.pos {
		t[i] = kri.pos[i][2]
	}

	var M = matrixMultiply(C, t, n, n, 1)
	kri.K = K
	kri.M = M

	return kri, nil
}

func (kri *Kriging) Predict(x, y float64) float64 {
	k := make([]float64, kri.N)
	for i := 0; i < kri.N; i++ {
		x_ := x - kri.pos[i][0]
		y_ := y - kri.pos[i][1]
		h := math.Sqrt(pow2(x_) + pow2(y_))
		k[i] = kri.model(
			h,
			kri.Nugget, kri.Range,
			kri.Sill, kri.A,
		)
	}

	return matrixMultiply(k, kri.M, 1, kri.N, 1)[0]
}

func (kri *Kriging) Grid(polygon PolygonCoordinates, width float64) *GridMatrices {
	n := len(polygon)
	if n == 0 {
		return &GridMatrices{}
	}

	var nodataValue float64 = -9999

	xlim := [2]float64{polygon[0][0][0], polygon[0][0][0]}
	ylim := [2]float64{polygon[0][0][1], polygon[0][0][1]}

	for i := 0; i < n; i++ {
		for j := 0; j < len(polygon[i]); j++ {
			if polygon[i][j][0] < xlim[0] {
				xlim[0] = polygon[i][j][0]
			}
			if polygon[i][j][0] > xlim[1] {
				xlim[1] = polygon[i][j][0]
			}
			if polygon[i][j][1] < ylim[0] {
				ylim[0] = polygon[i][j][1]
			}
			if polygon[i][j][1] > ylim[1] {
				ylim[1] = polygon[i][j][1]
			}
		}
	}

	x := int(math.Ceil((xlim[1] - xlim[0]) / width))
	y := int(math.Ceil((ylim[1] - ylim[0]) / width))

	A := make([][]float64, x+1)
	for i := 0; i <= x; i++ {
		A[i] = make([]float64, y+1)
	}

	for i := 0; i < n; i++ {
		currentPolygon := polygon[i]
		var lxlim [2]float64
		var lylim [2]float64
		lxlim[0] = currentPolygon[0][0]
		lxlim[1] = lxlim[0]
		lylim[0] = currentPolygon[0][1]
		lylim[1] = lylim[0]
		for j := 1; j < len(currentPolygon); j++ {
			if currentPolygon[j][0] < lxlim[0] {
				lxlim[0] = currentPolygon[j][0]
			}
			if currentPolygon[j][0] > lxlim[1] {
				lxlim[1] = currentPolygon[j][0]
			}
			if currentPolygon[j][1] < lylim[0] {
				lylim[0] = currentPolygon[j][1]
			}
			if currentPolygon[j][1] > lylim[1] {
				lylim[1] = currentPolygon[j][1]
			}
		}

		var a, b [2]int

		a[0] = int(math.Floor(((lxlim[0] - math.Mod(lxlim[0]-xlim[0], width)) - xlim[0]) / width))
		a[1] = int(math.Ceil(((lxlim[1] - math.Mod(lxlim[1]-xlim[1], width)) - xlim[0]) / width))
		b[0] = int(math.Floor(((lylim[0] - math.Mod(lylim[0]-ylim[0], width)) - ylim[0]) / width))
		b[1] = int(math.Ceil(((lylim[1] - math.Mod(lylim[1]-ylim[1], width)) - ylim[0]) / width))

		var wg sync.WaitGroup
		predictCh := make(chan *PredictDate, (b[1]-b[0])*(a[1]-a[0]))
		var parallelPredict = func(j, k int, polygon []Point, xTarget, yTarget float64) {
			predictDate := &PredictDate{X: j, Y: k}
			predictDate.Value = kri.Predict(xTarget,
				yTarget,
			)
			predictCh <- predictDate
			defer wg.Done()
		}

		var xTarget, yTarget float64
		for j := a[0]; j <= a[1]; j++ {
			xTarget = xlim[0] + float64(j)*width
			for k := b[0]; k <= b[1]; k++ {
				yTarget = ylim[0] + float64(k)*width

				if pipFloat64(currentPolygon, xTarget, yTarget) {
					wg.Add(1)
					go parallelPredict(j, k, currentPolygon, xTarget, yTarget)
				} else {
					A[j][k] = nodataValue
				}
			}
		}

		go func() {
			wg.Wait()
			close(predictCh)
		}()

		for predictDate := range predictCh {
			if predictDate.Value != 0 {
				j := predictDate.X
				k := predictDate.Y
				A[j][k] = predictDate.Value
			}

		}
	}

	gridMatrices := &GridMatrices{
		Xlim:        xlim,
		Ylim:        ylim,
		Zlim:        [2]float64{minFloat64(kri.pos, 2), maxFloat64(kri.pos, 2)},
		Width:       width,
		Data:        A,
		NodataValue: nodataValue,
	}
	return gridMatrices
}

func (kri *Kriging) Contour(xWidth, yWidth int) *ContourRectangle {
	xlim := [2]float64{minFloat64(kri.pos, 0), maxFloat64(kri.pos, 0)}
	ylim := [2]float64{minFloat64(kri.pos, 1), maxFloat64(kri.pos, 1)}
	zlim := [2]float64{minFloat64(kri.pos, 2), maxFloat64(kri.pos, 2)}
	xl := xlim[1] - xlim[0]
	yl := ylim[1] - ylim[0]
	gridW := xl / float64(xWidth)
	gridH := yl / float64(yWidth)
	var contour []float64

	var xTarget, yTarget float64

	for j := 0; j < yWidth; j++ {
		yTarget = ylim[0] + float64(j)*gridW
		for k := 0; k < xWidth; k++ {
			xTarget = xlim[0] + float64(k)*gridH
			contour = append(contour, kri.Predict(xTarget, yTarget))
		}
	}

	contourRectangle := &ContourRectangle{
		Contour:     contour,
		XWidth:      xWidth,
		YWidth:      yWidth,
		Xlim:        xlim,
		Ylim:        ylim,
		Zlim:        zlim,
		XResolution: 1,
		YResolution: 1,
	}

	return contourRectangle
}

func (kri *Kriging) ContourWithBBox(bbox [4]float64, width float64) *ContourRectangle {
	xlim := [2]float64{bbox[0], bbox[2]}
	ylim := [2]float64{bbox[1], bbox[3]}
	zlim := [2]float64{minFloat64(kri.pos, 2), maxFloat64(kri.pos, 2)}

	geoXWidth := xlim[1] - xlim[0]
	geoYWidth := ylim[1] - ylim[0]

	xWidth := int(math.Ceil(width))
	yWidth := int(math.Ceil(float64(xWidth) / (geoXWidth / geoYWidth)))

	var xResolution = geoXWidth * 1.0 / float64(xWidth)
	var yResolution = geoYWidth * 1.0 / float64(yWidth)

	var xTarget, yTarget float64
	var contour []float64

	for j := 0; j < yWidth; j++ {
		yTarget = bbox[1] + float64(j)*yResolution
		for k := 0; k < xWidth; k++ {
			xTarget = bbox[0] + float64(k)*xResolution
			contour = append(contour, kri.Predict(xTarget, yTarget))
		}
	}
	contourRectangle := &ContourRectangle{
		Contour:     contour,
		XWidth:      xWidth,
		YWidth:      yWidth,
		Xlim:        xlim,
		Ylim:        ylim,
		Zlim:        zlim,
		XResolution: xResolution,
		YResolution: yResolution,
	}

	return contourRectangle
}
