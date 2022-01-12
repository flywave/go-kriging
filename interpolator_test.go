package kriging

import (
	"io/ioutil"
	"os"
	"testing"

	"github.com/flywave/go-geom/general"
)

func TestInterpolator(t *testing.T) {

	f, _ := os.Open("./test.json")

	json, _ := ioutil.ReadAll(f)

	fcs, _ := general.UnmarshalFeatureCollection(json)
	m := ModelType("spherical")

	bg := "./test.tif"

	opts := Options{
		Input:      fcs,
		Output:     "./out.tif",
		Model:      &m,
		Background: &bg,
	}

	ker := NewKrigingInterpolator(opts)

	err := ker.Process()

	if err != nil {
		t.FailNow()
	}

}
