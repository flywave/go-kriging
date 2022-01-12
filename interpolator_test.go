package kriging

import (
	"io/ioutil"
	"os"
	"testing"

	"github.com/flywave/go-geom/general"
)

func TestInterpolator1(t *testing.T) {

	f, _ := os.Open("./test.json")

	json, _ := ioutil.ReadAll(f)

	fcs, _ := general.UnmarshalFeatureCollection(json)
	m := ModelType("spherical")

	opts := Options{
		Input:  fcs,
		Output: "./out1.tif",
		Model:  &m,
	}

	ker := NewKrigingInterpolator(opts)

	err := ker.Process()

	if err != nil {
		t.FailNow()
	}

}

func TestInterpolator(t *testing.T) {

	f, _ := os.Open("./test.json")

	json, _ := ioutil.ReadAll(f)

	fcs, _ := general.UnmarshalFeatureCollection(json)
	m := ModelType("spherical")

	bg := "./test.tif"

	opts := Options{
		Input:      fcs,
		Output:     "./out2.tif",
		Model:      &m,
		Background: &bg,
	}

	ker := NewKrigingInterpolator(opts)

	err := ker.Process()

	if err != nil {
		t.FailNow()
	}

}
