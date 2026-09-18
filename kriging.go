package kriging

import (
	"errors"
	"fmt"
	"math"
	"runtime"
	"sort"
	"sync"

	vec3d "github.com/flywave/go3d/float64/vec3"
	"gonum.org/v1/gonum/mat"
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

	model     KrigingModel
	modelType ModelType
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

	kri.modelType = model
	switch model {
	case Gaussian:
		kri.model = krigingKrigingGaussian
	case Exponential:
		kri.model = krigingKrigingExponential
	case Spherical:
		kri.model = krigingKrigingSpherical
	default:
		// 没有 default 时 kri.model 会保持 nil，下面并行构造 K 矩阵的
		// worker goroutine 里调用 nil 函数会直接崩掉整个进程（不可恢复）
		return nil, fmt.Errorf("unknown kriging model type: %q", model)
	}

	var i, j, k, l, n int
	n = len(kri.pos)
	if n < 3 {
		return nil, fmt.Errorf("kriging needs at least 3 points, got %d", n)
	}

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
	for ii := 0; ii < len(X); ii++ {
		X[ii] = 1
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

	numCPU := runtime.GOMAXPROCS(0)
	var wg sync.WaitGroup
	rowsPerWorker := (n + numCPU - 1) / numCPU
	for wi := 0; wi < numCPU; wi++ {
		start := wi * rowsPerWorker
		end := start + rowsPerWorker
		if end > n {
			end = n
		}
		wg.Add(1)
		go func(start, end int) {
			defer wg.Done()
			nugget := kri.nugget
			rangex := kri.rangex
			sill := kri.sill
			a := kri.A
			model := kri.model
			diag := model(0, nugget, rangex, sill, a)
			for i := start; i < end; i++ {
				xi, yi := kri.pos[i][0], kri.pos[i][1]
				K[i*n+i] = diag
				for j := 0; j < i; j++ {
					dx := xi - kri.pos[j][0]
					dy := yi - kri.pos[j][1]
					h := math.Sqrt(dx*dx + dy*dy)
					val := model(h, nugget, rangex, sill, a)
					K[i*n+j] = val
					K[j*n+i] = val
				}
			}
		}(start, end)
	}
	wg.Wait()

	C := matrixAdd(K, matrixDiag(sigma2, n), n, n)

	t := make([]float64, n)
	for ii := range kri.pos {
		t[ii] = kri.pos[ii][2]
	}

	cMat := mat.NewSymDense(n, C)
	var chol mat.Cholesky
	if chol.Factorize(cMat) {
		zVec := mat.NewVecDense(n, t)
		var mVec mat.VecDense
		if err := chol.SolveVecTo(&mVec, zVec); err == nil {
			kri.M = mVec.RawVector().Data
		}
	} else {
		matrixSolve(C, n)
		kri.M = matrixMultiply(C, t, n, n, 1)
	}
	kri.K = nil

	return kri, nil
}

func (kri *Kriging) Predict(x, y float64) float64 {
	nugget := kri.nugget
	rangex := kri.rangex
	sill := kri.sill
	a := kri.A
	pos := kri.pos
	M := kri.M
	invRange := 1.0 / rangex
	sillMinusNugget := (sill - nugget) * invRange

	var result float64

	switch kri.modelType {
	case Gaussian:
		for i := 0; i < kri.n; i++ {
			dx := x - pos[i][0]
			dy := y - pos[i][1]
			h2 := (dx*dx + dy*dy) * invRange * invRange
			result += (nugget + sillMinusNugget*(1.0-math.Exp(-(1.0/a)*h2))) * M[i]
		}
	case Exponential:
		oneOverA := 1.0 / a
		for i := 0; i < kri.n; i++ {
			dx := x - pos[i][0]
			dy := y - pos[i][1]
			h := math.Sqrt(dx*dx+dy*dy) * invRange
			result += (nugget + sillMinusNugget*(1.0-math.Exp(-oneOverA*h))) * M[i]
		}
	case Spherical:
		for i := 0; i < kri.n; i++ {
			dx := x - pos[i][0]
			dy := y - pos[i][1]
			h := math.Sqrt(dx*dx + dy*dy)
			if h > rangex {
				result += (nugget + sillMinusNugget) * M[i]
			} else {
				xr := h * invRange
				result += (nugget + sillMinusNugget*(1.5*xr-0.5*xr*xr*xr)) * M[i]
			}
		}
	}

	return result
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
		yTarget = ylim[0] + float64(j)*gridH
		for k := 0; k < xWidth; k++ {
			xTarget = xlim[0] + float64(k)*gridW
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
