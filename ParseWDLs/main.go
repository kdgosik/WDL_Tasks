package main

import (
	"fmt"
	cid "github.com/ipfs/go-cid"
	mh "github.com/multiformats/go-multihash"
	//  ipld "github.com/ipfs/go-ipld-format"
	blocks "github.com/ipfs/go-block-format"
	"regexp"
)

type WDL struct {
	// id used in firecloud's api
	Id   string `json:"id"`
	// cid created from below
	Hash string `json:"hash"`
	// // actual wdl content (may not need / just use api)
	// Output string `json:"output"`
}

type Link struct {
	// parent WDL id or hash
	Parent string `json:"parent"`
	// child WDL id or hash
	Child  string `json:"child"`
	// Combined hash (concatenate then hash again)
	Hash   string `json:"hash"`
}

// inputs WDL and returns a hash of the content
func HashWDL(id string, wdl string) WDL {

	// using cid to create it from scratch
	pref := cid.Prefix{
		Version:  1,
		Codec:    cid.Raw,
		MhType:   mh.SHA2_256,
		MhLength: -1, // default length
	}

	// hashes the wdl
	c, err := pref.Sum([]byte(wdl))
	if err != nil {
		fmt.Println(err)
	}

	// cretes a block given the hash above
	hash, _ := blocks.NewBlockWithCid([]byte(wdl), c)

	// outputs as a WDL struct
	out := WDL{
		Id:   id,
		Hash: hash.Cid().String(),
	}

	return out
}


// inputs a parent and a child wdl and returns a hash
func HashLink(parent string, child string) Link {
	// using cid to create it from scratch
	pref := cid.Prefix{
		Version:  1,
		Codec:    cid.Raw,
		MhType:   mh.SHA2_256,
		MhLength: -1, // default length
	}

	// creates a concatenated hash
	c, err := pref.Sum([]byte(parent + child))
	if err != nil {
		fmt.Println(err)
	}

	// creates a link struct of concatenated hash
	hash, _ := blocks.NewBlockWithCid([]byte(parent+child), c)

	out := Link{
		Parent: parent,
		Child:  child,
		Hash:   hash.Cid().String(),
	}

	return out
}


func main() {

// from api_get_workflows.go  gets all tool urls
	tool_urls := getToolsURL()

// loops through urls and returns a slice of WDL hashes
	var wdl_hashes []WDL
	for _, url := range tool_urls[:5] {
		reg, _ := regexp.Compile(`https://api.firecloud.org/ga4gh/v1/tools/(.+?)/versions/`)
		reg_result := reg.FindAllStringSubmatch(url, -1)
		fmt.Println(reg_result)
		id := reg_result[0][1]
		wf := getWorkflow(url)

		wdl_hashes = append(wdl_hashes, HashWDL(id, wf))

	}
	fmt.Println(wdl_hashes)


// Example to be DELETED
/*

	// using cid to create it from scratch
	pref := cid.Prefix{
		Version:  1,
		Codec:    cid.Raw,
		MhType:   mh.SHA2_256,
		MhLength: -1, // default length
	}

	// And then feed it some data
	c, err := pref.Sum([]byte("Hello World!"))
	if err != nil {
		fmt.Println(err)
	}

	fmt.Println("Created CID: ", c)

	// using go-block-format to create block
	out, _ := blocks.NewBlockWithCid([]byte("Hello World!"), c)

	// now need to extend to a Node interface in ipld-format
	fmt.Println(out)
	fmt.Println(out.RawData())
	fmt.Println(string(out.RawData()))
	fmt.Println(out.Cid())
	fmt.Println(out.String())
	fmt.Println(out.Loggable())
*/
}
