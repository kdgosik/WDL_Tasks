package main

import (
    "fmt"
    "encoding/json"
    "io/ioutil"
)

type WDL struct {
	// id used in firecloud's api
	Id       string `json:"id"`
	// version of the wdl
	Version  string `json:"version"`
	// cid created from below
	Hash     string `json:"hash"`
	// // actual wdl content (may not need / just use api)
	// Output string `json:"output"`
}

func main() {

	content, err := ioutil.ReadFile("wdl_hash.json")
	if err != nil {
		fmt.Println(err)
	}

    var out []WDL
    err = json.Unmarshal(content, &out)
    if err != nil {
        panic(err)
    }
    fmt.Println(out)
}
