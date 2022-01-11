package kriging

import (
	"math"
	"sort"
	"sync"

	geo "github.com/paulmach/go.geo"

	"github.com/furstenheim/SimpleRTree"
	"github.com/paulmach/go.geo/reducers"
)

const DEFAULT_SEGLENGTH = 0.001

type concaver struct {
	rtree            *SimpleRTree.SimpleRTree
	seglength        float64
	options          *Options
	closestPointsMem []closestPoint
	searchItemsMem   []searchItem
	flatPointBuffer  []float64
	rtreePool        *sync.Pool
}
type Options struct {
	Seglength                   float64
	EstimatedRatioConcaveConvex int
	ConcaveHullPool             *sync.Pool
}

type concaveHullPoolElement struct {
	fpbMem           []float64
	closestPointsMem []closestPoint
	searchItemsMem   []searchItem
	rtreePool        *sync.Pool
	convexHullPool   *sync.Pool
	pointsCopy       FlatPoints
}

func Compute(points FlatPoints) (concaveHull FlatPoints) {
	return ComputeWithOptions(points, nil)
}
func ComputeWithOptions(points FlatPoints, o *Options) (concaveHull FlatPoints) {
	sort.Sort(lexSorter(points))
	return ComputeFromSortedWithOptions(points, o)
}
func ComputeFromSorted(points FlatPoints) (concaveHull FlatPoints) {
	return ComputeFromSortedWithOptions(points, nil)
}

func ComputeFromSortedWithOptions(points FlatPoints, o *Options) (concaveHull FlatPoints) {
	var pointsCopy FlatPoints
	var rtreeOptions SimpleRTree.Options
	var isConcaveHullPoolElementsSet bool
	var poolEl *concaveHullPoolElement
	var rtreePool *sync.Pool
	var convexHullPool *sync.Pool
	if o != nil && o.ConcaveHullPool != nil {
		concaveHullPoolElementsCandidate := o.ConcaveHullPool.Get()
		if concaveHullPoolElementsCandidate != nil {
			isConcaveHullPoolElementsSet = true
			poolEl = concaveHullPoolElementsCandidate.(*concaveHullPoolElement)
			rtreePool = poolEl.rtreePool
			convexHullPool = poolEl.convexHullPool
		} else {
			rtreePool = &sync.Pool{}
			convexHullPool = &sync.Pool{}
		}

	}
	if isConcaveHullPoolElementsSet && cap(poolEl.pointsCopy) >= len(points) {
		pointsCopy = poolEl.pointsCopy[0:0]
	} else {
		pointsCopy = make(FlatPoints, 0, len(points))
	}
	pointsCopy = append(pointsCopy, points...)
	rtreeOptions.RTreePool = rtreePool
	rtreeOptions.UnsafeConcurrencyMode = true
	rtree := SimpleRTree.NewWithOptions(rtreeOptions)
	var wg sync.WaitGroup
	wg.Add(2)
	go func() {
		points = NewFromSortedArrayWithOptions(points, cOptions{Pool: convexHullPool}).(FlatPoints)
		wg.Done()
	}()

	func() {
		rtree.LoadSortedArray(SimpleRTree.FlatPoints(pointsCopy))
		wg.Done()
	}()
	wg.Wait()
	var c concaver
	c.seglength = DEFAULT_SEGLENGTH
	if o != nil && o.Seglength != 0 {
		c.seglength = o.Seglength
	}
	c.rtree = rtree
	if isConcaveHullPoolElementsSet {
		c.closestPointsMem = poolEl.closestPointsMem
		c.searchItemsMem = poolEl.searchItemsMem
		c.flatPointBuffer = poolEl.fpbMem[0:0]

	} else {
		c.closestPointsMem = make([]closestPoint, 0, 2)
		c.searchItemsMem = make([]searchItem, 0, 2)
		estimatedProportionConcave2Convex := 4
		if c.options != nil && c.options.EstimatedRatioConcaveConvex != 0 {
			estimatedProportionConcave2Convex = c.options.EstimatedRatioConcaveConvex
		}
		c.flatPointBuffer = make([]float64, 0, (2 * points.Len() * estimatedProportionConcave2Convex))
	}

	result := c.computeFromSorted(points)
	rtree.Destroy()
	if o != nil && o.ConcaveHullPool != nil {
		o.ConcaveHullPool.Put(
			&concaveHullPoolElement{
				rtreePool:        rtreePool,
				convexHullPool:   convexHullPool,
				searchItemsMem:   c.searchItemsMem,
				closestPointsMem: c.closestPointsMem,
				fpbMem:           c.flatPointBuffer,
				pointsCopy:       pointsCopy,
			},
		)
	}
	return result
}

func (c *concaver) computeFromSorted(convexHull FlatPoints) (concaveHull FlatPoints) {
	if convexHull.Len() < 3 {
		return convexHull
	}

	x0, y0 := convexHull.Take(0)
	concaveHullBuffer := c.flatPointBuffer
	concaveHullBuffer = append(concaveHullBuffer, x0, y0)
	for i := 0; i < convexHull.Len(); i++ {
		x1, y1 := convexHull.Take(i)
		var x2, y2 float64
		if i == convexHull.Len()-1 {
			x2, y2 = convexHull.Take(0)
		} else {
			x2, y2 = convexHull.Take(i + 1)
		}
		sideSplit := c.segmentize(x1, y1, x2, y2)
		for _, p := range sideSplit {
			concaveHullBuffer = append(concaveHullBuffer, p.x, p.y)
		}
	}
	concaveHull = make([]float64, 0, len(concaveHullBuffer))
	concaveHull = append(concaveHull, concaveHullBuffer...)
	path := reducers.DouglasPeucker(geo.NewPathFromFlatXYData(concaveHull), c.seglength)

	concaveHull = concaveHull[0:0]
	reducedPoints := path.Points()

	for _, p := range reducedPoints {
		concaveHull = append(concaveHull, p.Lng(), p.Lat())
	}
	return concaveHull
}

func (c *concaver) segmentize(x1, y1, x2, y2 float64) (points []closestPoint) {
	dist := math.Sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))
	nSegments := math.Ceil(dist / c.seglength)
	factor := 1 / nSegments
	vX := factor * (x2 - x1)
	vY := factor * (y2 - y1)

	closestPoints := c.closestPointsMem[0:0]
	closestPoints = append(closestPoints, closestPoint{index: 0, x: x1, y: y1})
	closestPoints = append(closestPoints, closestPoint{index: int(nSegments), x: x2, y: y2})

	if nSegments < 2 {
		return closestPoints[1:]
	}

	stack := c.searchItemsMem[0:0]
	stack = append(stack, searchItem{left: 0, right: int(nSegments), lastLeftIndex: 0, lastRightIndex: 1})
	for len(stack) > 0 {
		var item searchItem
		item, stack = stack[len(stack)-1], stack[:len(stack)-1]
		if item.right-item.left <= 1 {
			continue
		}
		index := (item.left + item.right) / 2
		fIndex := float64(index)
		currentX := x1 + vX*fIndex
		currentY := y1 + vY*fIndex
		lx := closestPoints[item.lastLeftIndex].x
		ly := closestPoints[item.lastLeftIndex].y
		rx := closestPoints[item.lastRightIndex].x
		ry := closestPoints[item.lastRightIndex].y

		d1 := (currentX-lx)*(currentX-lx) + (currentY-ly)*(currentY-ly)
		d2 := (currentX-rx)*(currentX-rx) + (currentY-ry)*(currentY-ry)
		x, y, _, found := c.rtree.FindNearestPointWithin(currentX, currentY, math.Min(d1, d2))
		if !found {
			continue
		}
		isNewLeft := x != lx || y != ly
		isNewRight := x != rx || y != ry

		if isNewLeft && isNewRight {
			newResultIndex := len(closestPoints)
			closestPoints = append(closestPoints, closestPoint{index: index, x: x, y: y})
			stack = append(stack, searchItem{left: item.left, right: index, lastLeftIndex: item.lastLeftIndex, lastRightIndex: newResultIndex})
			stack = append(stack, searchItem{left: index, right: item.right, lastLeftIndex: newResultIndex, lastRightIndex: item.lastRightIndex})
		} else if isNewLeft {
			stack = append(stack, searchItem{left: item.left, right: index, lastLeftIndex: item.lastLeftIndex, lastRightIndex: item.lastRightIndex})
		} else {
			stack = append(stack, searchItem{left: index, right: item.right, lastLeftIndex: item.lastLeftIndex, lastRightIndex: item.lastRightIndex})
		}
	}
	closestPointSorter(closestPoints).cpSort()
	c.searchItemsMem = stack
	c.closestPointsMem = closestPoints
	return closestPoints[1:]
}

type closestPoint struct {
	index int
	x, y  float64
}

type searchItem struct {
	left, right, lastLeftIndex, lastRightIndex int
}

type FlatPoints []float64

func (fp FlatPoints) Len() int {
	return len(fp) / 2
}

func (fp FlatPoints) Slice(i, j int) Interface {
	return fp[2*i : 2*j]
}

func (fp FlatPoints) Swap(i, j int) {
	fp[2*i], fp[2*i+1], fp[2*j], fp[2*j+1] = fp[2*j], fp[2*j+1], fp[2*i], fp[2*i+1]
}

func (fp FlatPoints) Take(i int) (x1, y1 float64) {
	return fp[2*i], fp[2*i+1]
}

type Interface interface {
	Take(i int) (x, y float64)
	Len() int
	Swap(i, j int)
	Slice(i, j int) Interface
}

type cOptions struct {
	Pool *sync.Pool
}
type poolElStruct struct {
	lowerIndexes, upperIndexes []int
}

func cNew(points Interface) Interface {
	sort.Sort(pointSorter{i: points})
	return NewFromSortedArray(points)
}

func NewFromSortedArray(points Interface) Interface {
	o := cOptions{}
	return NewFromSortedArrayWithOptions(points, o)
}
func NewWithOptions(points Interface, o cOptions) Interface {
	sort.Sort(pointSorter{i: points})
	return NewFromSortedArrayWithOptions(points, o)
}

func NewFromSortedArrayWithOptions(points Interface, o cOptions) Interface {
	n := points.Len()
	if n < 3 {
		return points
	}
	var w sync.WaitGroup
	var lowerIndexes []int
	var upperIndexes []int
	var isPooledMemReceived bool
	if o.Pool != nil {
		poolElCandidate := o.Pool.Get()
		if poolElCandidate != nil {
			isPooledMemReceived = true
			poolEl := poolElCandidate.(*poolElStruct)
			lowerIndexes = poolEl.lowerIndexes[0:0]
			upperIndexes = poolEl.upperIndexes[0:0]
		}
	}
	if !isPooledMemReceived {
		lowerIndexes = make([]int, 0, 5)
		upperIndexes = make([]int, 0, 5)
	}
	lowerIndexes = append(lowerIndexes, 0, 1)
	upperIndexes = append(upperIndexes, n-1, n-2)
	w.Add(2)
	func() {
		for i := 2; i < n; i++ {
			x, y := points.Take(i)
			m := len(lowerIndexes)
			for m > 1 {
				x2, y2 := points.Take(lowerIndexes[m-2])
				x3, y3 := points.Take(lowerIndexes[m-1])
				if isOrientationPositive(x2, y2, x3, y3, x, y) {
					break
				}
				lowerIndexes = lowerIndexes[:m-1]
				m -= 1
			}
			lowerIndexes = append(lowerIndexes, i)
		}

		w.Done()
	}()
	func() {
		for i := n - 3; i >= 0; i-- {
			x, y := points.Take(i)
			m := len(upperIndexes)
			for m > 1 {
				x2, y2 := points.Take(upperIndexes[m-2])
				x3, y3 := points.Take(upperIndexes[m-1])
				if isOrientationPositive(x2, y2, x3, y3, x, y) {
					break
				}
				upperIndexes = upperIndexes[:m-1]
				m -= 1
			}
			upperIndexes = append(upperIndexes, i)

		}
		w.Done()
	}()
	w.Wait()

	upperIndexes = upperIndexes[:len(upperIndexes)-1]
	lowerIndexes = lowerIndexes[:len(lowerIndexes)-1]
	allIndexes := append(lowerIndexes, upperIndexes...)

	result := sortByIndexes(points, allIndexes)
	if o.Pool != nil {
		o.Pool.Put(&poolElStruct{lowerIndexes: allIndexes, upperIndexes: upperIndexes})
	}
	return result
}

func isOrientationPositive(x1, y1, x2, y2, x3, y3 float64) (isPositive bool) {
	// compute determinant to obtain the orientation
	// |x1 - x3 x2 - x3 |
	// |y1 - y3 y2 - y3 |
	return (x1-x3)*(y2-y3)-(y1-y3)*(x2-x3) > 0
}

type pointSorter struct {
	i Interface
}

func (s pointSorter) Less(i, j int) bool {
	x1, y1 := s.i.Take(i)
	x2, y2 := s.i.Take(j)
	return x1 < x2 || (x1 == x2 && y1 < y2)
}

func (s pointSorter) Swap(i, j int) {
	s.i.Swap(i, j)
}

func (s pointSorter) Len() int {
	return s.i.Len()
}

func sortByIndexes(points Interface, indices []int) Interface {
	s := indexSorter{indices: indices, points: points}
	sort.Sort(s)
	for i, index := range indices {
		points.Swap(i, index)
	}
	return points.Slice(0, len(indices))
}

type indexSorter struct {
	indices []int
	points  Interface
}

func (s indexSorter) Less(i, j int) bool {
	return s.indices[i] < s.indices[j]
}

func (s indexSorter) Swap(i, j int) {
	s.indices[i], s.indices[j] = s.indices[j], s.indices[i]
	s.points.Swap(s.indices[i], s.indices[j])
}

func (s indexSorter) Len() int {
	return len(s.indices)
}
