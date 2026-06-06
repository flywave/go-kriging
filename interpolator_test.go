package kriging

import (
	"io"
	"os"
	"testing"

	"github.com/flywave/go-geom/general"
)

func TestInterpolator1(t *testing.T) {

	f, _ := os.Open("./test.json")

	json, _ := io.ReadAll(f)

	fcs, _ := general.UnmarshalFeatureCollection(json)
	m := ModelType("spherical")

	bg := "./test.tif"

	opts := Options{
		Input:      fcs,
		Output:     "./out1.tif",
		Model:      &m,
		Background: &bg,
	}

	ker := NewKrigingInterpolator(opts)

	_, _, err := ker.Process()

	if err != nil {
		t.FailNow()
	}

}
