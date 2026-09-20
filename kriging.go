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

// dedupByTolerance 合并平面距离小于 eps 的点（z 取均值），保持首次出现顺序。
// 用格网哈希把复杂度从 O(n²) 降到近似 O(n)；eps <= 0 时只合并坐标完全相同的点。
func dedupByTolerance(pos []vec3d.T, eps float64) []vec3d.T {
	if eps <= 0 || math.IsNaN(eps) || math.IsInf(eps, 0) {
		// 无容差：只合并坐标完全相同的点
		type acc struct {
			z   float64
			num int
		}
		seen := make(map[[2]float64]int, len(pos))
		out := make([]vec3d.T, 0, len(pos))
		var accs []acc
		for _, p := range pos {
			key := [2]float64{p[0], p[1]}
			if idx, ok := seen[key]; ok {
				accs[idx].z += p[2]
				accs[idx].num++
			} else {
				seen[key] = len(out)
				out = append(out, p)
				accs = append(accs, acc{z: p[2], num: 1})
			}
		}
		for i := range out {
			if accs[i].num > 1 {
				out[i][2] = accs[i].z / float64(accs[i].num)
			}
		}
		return out
	}
	type acc struct {
		z   float64
		num int
	}
	cells := make(map[[2]int][]int) // 格网 -> 保留点下标
	out := make([]vec3d.T, 0, len(pos))
	accs := make([]acc, 0, len(pos))
	for _, p := range pos {
		cx, cy := int(math.Floor(p[0]/eps)), int(math.Floor(p[1]/eps))
		merged := false
		for dx := -1; dx <= 1 && !merged; dx++ {
			for dy := -1; dy <= 1 && !merged; dy++ {
				for _, idx := range cells[[2]int{cx + dx, cy + dy}] {
					q := out[idx]
					d := math.Hypot(p[0]-q[0], p[1]-q[1])
					if d <= eps {
						accs[idx].z += p[2]
						accs[idx].num++
						merged = true
						break
					}
				}
			}
		}
		if !merged {
			cells[[2]int{cx, cy}] = append(cells[[2]int{cx, cy}], len(out))
			out = append(out, p)
			accs = append(accs, acc{z: p[2], num: 1})
		}
	}
	for i := range out {
		if accs[i].num > 1 {
			out[i][2] = accs[i].z / float64(accs[i].num)
		}
	}
	return out
}

func bboxDiagonal(pos []vec3d.T) float64 {
	minX, minY := pos[0][0], pos[0][1]
	maxX, maxY := minX, minY
	for _, p := range pos[1:] {
		if p[0] < minX {
			minX = p[0]
		}
		if p[0] > maxX {
			maxX = p[0]
		}
		if p[1] < minY {
			minY = p[1]
		}
		if p[1] > maxY {
			maxY = p[1]
		}
	}
	return math.Hypot(maxX-minX, maxY-minY)
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

	// 平面坐标近乎重合的点（如同一位置多次测量、过密采样）会让 K 矩阵的
	// 行/列到求解器精度下完全相同，系统奇异。按容差合并，z 取均值。
	// 容差取包围盒对角线的 1e-3：远小于数据范围，又足以吸收测量噪声级的重合。
	deduped := dedupByTolerance(kri.pos, bboxDiagonal(kri.pos)*1e-3)
	if len(deduped) != len(kri.pos) {
		if len(deduped) < 3 {
			return nil, fmt.Errorf("kriging needs at least 3 distinct points, got %d", len(deduped))
		}
		kri.pos = deduped
		n = len(deduped)
	}

	distance := make([][2]float64, (n*n-n)/2)

	i = 0
	k = 0
	for ; i < n; i++ {
		for j = 0; j < i; {
			dx := kri.pos[i][0] - kri.pos[j][0]
			dy := kri.pos[i][1] - kri.pos[j][1]
			dz := kri.pos[i][2] - kri.pos[j][2]
			distance[k] = [2]float64{math.Sqrt(dx*dx + dy*dy), math.Abs(dz)}
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
	// 退化保护：所有点对平面距离几乎相同（如等边三角形三点）时 range 为 0，
	// 后续 lag/rangex 与 1/rangex 会产生 NaN/Inf。退回用最大点对距离。
	if !(kri.rangex > 0) || math.IsNaN(kri.rangex) || math.IsInf(kri.rangex, 0) {
		kri.rangex = distance[(len(kri.pos)*len(kri.pos)-len(kri.pos))/2-1][0]
		if !(kri.rangex > 0) {
			return nil, errors.New("kriging: all sample points are coincident")
		}
	}
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
	// 拟合是无约束线性回归，可能给出物理上非法的参数（负块金、负基台）。
	// 负块金会让 K 矩阵失去半正定性，求解必然失败；钳制到有效域。
	if kri.nugget < 0 {
		kri.nugget = 0
	}
	if kri.sill < kri.nugget {
		kri.sill = kri.nugget
	}
	kri.n = len(kri.pos)

	// 斜率（sill−nugget）退化为 0：数据没有可用的空间结构（如 z 全等的
	// 平面点集），γ≡0，K 矩阵全零不可解。降级为常数（均值）预测器：
	// 取 γ≡1，M_i = mean(z)/n，则 Predict = Σ γ·M_i ≡ mean(z)。
	if kri.sill-kri.nugget <= 0 {
		mean := 0.0
		for i := range kri.pos {
			mean += kri.pos[i][2]
		}
		mean /= float64(len(kri.pos))
		kri.nugget = 1
		kri.sill = 1
		kri.M = make([]float64, kri.n)
		for i := range kri.M {
			kri.M[i] = mean / float64(kri.n)
		}
		kri.K = nil
		return kri, nil
	}

	n = len(kri.pos)

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

	// 解 K·M = t。优先 Cholesky（对称正定，快），失败退回通用 LU；
	// 解出后校验残差，病态系统（如 gaussian 小 range 大点集、σ²=0 无正则）
	// 会解出 1e30 量级的垃圾权重，必须显式报错而不是静默输出。
	zVec := mat.NewVecDense(n, t)
	var mVec mat.VecDense
	solved := false
	cMat := mat.NewSymDense(n, C)
	var chol mat.Cholesky
	if chol.Factorize(cMat) {
		if err := chol.SolveVecTo(&mVec, zVec); err == nil {
			solved = true
		}
	}
	if !solved {
		if err := mVec.SolveVec(mat.NewDense(n, n, C), zVec); err != nil {
			return nil, fmt.Errorf("kriging: solving kriging system failed: %w", err)
		}
	}
	M := mVec.RawVector().Data
	if relResidual(C, M, t, n) > 1e-6 {
		return nil, errors.New("kriging: ill-conditioned kriging system, " +
			"predicted weights are unreliable; increase sigma2, change the variogram model, or decimate the sample points")
	}
	kri.M = M
	kri.K = nil

	return kri, nil
}

// relResidual 计算 ||C·M − t||∞ / max(||t||∞, 1)。
// M 含 NaN/Inf 或残差过大都返回 +Inf，调用方据此判定系统病态。
func relResidual(C, M, t []float64, n int) float64 {
	tNorm := 0.0
	for i := 0; i < n; i++ {
		if a := math.Abs(t[i]); a > tNorm {
			tNorm = a
		}
	}
	if tNorm < 1 {
		tNorm = 1
	}
	rNorm := 0.0
	for i := 0; i < n; i++ {
		s := 0.0
		for j := 0; j < n; j++ {
			s += C[i*n+j] * M[j]
		}
		r := s - t[i]
		if math.IsNaN(r) || math.IsInf(r, 0) {
			return math.Inf(1)
		}
		if a := math.Abs(r); a > rNorm {
			rNorm = a
		}
	}
	return rNorm / tNorm
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
