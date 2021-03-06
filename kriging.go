package kriging

import (
	"errors"
	"math"
	"sort"

	vec3d "github.com/flywave/go3d/float64/vec3"
)

type Kriging struct {
	pos []vec3d.T

	nugget float64
	rangex float64
	sill   float64
	A      float64
	n      int

	K []float64
	M []float64

	model KrigingModel
}

func New(pos []vec3d.T) *Kriging {
	return &Kriging{pos: pos}
}

type KrigingModel func(float64, float64, float64, float64, float64) float64

func krigingKrigingGaussian(h, nugget, range_, sill, A float64) float64 {
	x := -(1.0 / A) * ((h / range_) * (h / range_))
	return nugget + ((sill-nugget)/range_)*(1.0-math.Exp(x))
}

func krigingKrigingExponential(h, nugget, range_, sill, A float64) float64 {
	x := -(1.0 / A) * (h / range_)
	return nugget + ((sill-nugget)/range_)*(1.0-math.Exp(x))
}

func krigingKrigingSpherical(h, nugget, range_, sill, A float64) float64 {
	if h > range_ {
		return nugget + (sill-nugget)/range_
	} else {
		x := h / range_
		return nugget + ((sill-nugget)/range_)*(1.5*(x)-0.5*(math.Pow(x, 3)))
	}
}

func (kri *Kriging) Train(model ModelType, sigma2 float64, alpha float64) (*Kriging, error) {
	kri.nugget = 0.0
	kri.rangex = 0.0
	kri.sill = 0.0
	kri.A = float64(1) / float64(3)
	kri.n = 0.0

	switch model {
	case Gaussian:
		kri.model = krigingKrigingGaussian
	case Exponential:
		kri.model = krigingKrigingExponential
	case Spherical:
		kri.model = krigingKrigingSpherical
	}

	var i, j, k, l, n int
	n = len(kri.pos)

	distance := make([][2]float64, (n*n-n)/2)

	i = 0
	k = 0
	for ; i < n; i++ {
		for j = 0; j < i; {
			distance[k] = [2]float64{}
			distance[k][0] = math.Pow(
				math.Pow(kri.pos[i][0]-kri.pos[j][0], 2)+
					math.Pow(kri.pos[i][1]-kri.pos[j][1], 2), 0.5)
			distance[k][1] = math.Abs(kri.pos[i][2] - kri.pos[j][2])
			j++
			k++
		}
	}
	sort.Sort(DistanceList(distance))
	kri.rangex = distance[(n*n-n)/2-1][0]

	var lags int
	if ((n*n - n) / 2) > 30 {
		lags = 30
	} else {
		lags = (n*n - n) / 2
	}

	tolerance := kri.rangex / float64(lags)

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
	kri.rangex = lag[n-1] - lag[0]
	X := make([]float64, 2*n)
	for i := 0; i < len(X); i++ {
		X[i] = 1
	}
	Y := make([]float64, n)
	var A = kri.A
	for i = 0; i < n; i++ {
		switch model {
		case Gaussian:
			X[i*2+1] = 1.0 - math.Exp(-(1.0/A)*math.Pow(lag[i]/kri.rangex, 2))
		case Exponential:
			X[i*2+1] = 1.0 - math.Exp(-(1.0/A)*lag[i]/kri.rangex)
		case Spherical:
			X[i*2+1] = 1.5*(lag[i]/kri.rangex) - 0.5*math.Pow(lag[i]/kri.rangex, 3)
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

	kri.nugget = W[0]
	kri.sill = W[1]*kri.rangex + kri.nugget
	kri.n = len(kri.pos)

	n = len(kri.pos)
	K := make([]float64, n*n)
	for i = 0; i < n; i++ {
		for j = 0; j < i; j++ {
			K[i*n+j] = kri.model(
				math.Pow(math.Pow(kri.pos[i][0]-kri.pos[j][0], 2)+
					math.Pow(kri.pos[i][1]-kri.pos[j][1], 2), 0.5),
				kri.nugget,
				kri.rangex,
				kri.sill,
				kri.A)
			K[j*n+i] = K[i*n+j]
		}
		K[i*n+i] = kri.model(0, kri.nugget,
			kri.rangex,
			kri.sill,
			kri.A)
	}

	var C = matrixAdd(K, matrixDiag(sigma2, n), n, n)
	var cloneC = make([]float64, len(C))
	copy(cloneC, C)
	if matrixChol(C, n) {
		matrixChol2inv(C, n)
	} else {
		matrixSolve(cloneC, n)
		C = cloneC
	}

	K = C
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
	k := make([]float64, kri.n)
	for i := 0; i < kri.n; i++ {
		x_ := x - kri.pos[i][0]
		y_ := y - kri.pos[i][1]
		h := math.Pow(math.Pow(x_, 2)+math.Pow(y_, 2), 0.5)
		k[i] = kri.model(
			h,
			kri.nugget, kri.rangex,
			kri.sill, kri.A,
		)
	}

	return matrixMultiply(k, kri.M, 1, kri.n, 1)[0]
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
