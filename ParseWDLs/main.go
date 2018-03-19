package main

import (
	"fmt"
	cid "github.com/ipfs/go-cid"
	mh "github.com/multiformats/go-multihash"
	//  ipld "github.com/ipfs/go-ipld-format"
	blocks "github.com/ipfs/go-block-format"
)

type WDL struct {
  // id used in firecloud's api
	Id     string `json:"id"`
  // cid created from below
	Hash   string `json:"hash"`
  // actual wdl content (may not need / just use api)
	Output string `json:"output"`
}

type Link struct {
  // parent WDL id or hash
	Parent    string `json:"parent"`
  // child WDL id or hash
	Child     string `json:"child"`
  // Combined hash (concatenate then hash again)
	Hash      string `json:"hash"`
}


// rough code structure
func  HashWDL(id string) {

  // using cid to create it from scratch
  pref := cid.Prefix{
    Version:  1,
    Codec:    cid.Raw,
    MhType:   mh.SHA2_256,
    MhLength: -1, // default length
  }

  wdl := get("https://firecloud.api/" + id)

  c, err := pref.Sum([]byte(wdl))
  	if err != nil {
  		fmt.Println(err)
  	}

  hash, _ := blocks.NewBlockWithCid([]byte(wdl), c)

  out = WDL{
    Id: id,
    Hash: hash.Cid()
  }

  return out
}


func HashLink(parent string, child string) {
  // using cid to create it from scratch
  pref := cid.Prefix{
    Version:  1,
    Codec:    cid.Raw,
    MhType:   mh.SHA2_256,
    MhLength: -1, // default length
  }

  c, err := pref.Sum([]byte(parent + child))
    if err != nil {
      fmt.Println(err)
    }

  hash, _ := blocks.NewBlockWithCid([]byte(parent + child), c)

  out = Link{
    Parent: parent
    Child: child,
    Hash: hash.Cid()
  }

  return out
}

func main() {

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

}
