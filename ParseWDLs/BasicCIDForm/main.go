package main

import (
	"github.com/ipfs/go-cid"
  mh "github.com/multiformats/go-multihash"
	"html/template"
	"log"
	"net/http"
	"net/url"
  "strings"
)


func displayHandler(w http.ResponseWriter, r *http.Request) {
	err := r.ParseForm()
	if err != nil {
		log.Fatalln(err)
	}

// input data from the Form
	data := struct {
		Method        string
		URL           *url.URL
		Submissions   map[string][]string
		Header        http.Header
		Host          string
		ContentLength int64
		Output        string
	}{
		r.Method,
		r.URL,
		r.Form,
		r.Header,
		r.Host,
		r.ContentLength,
		"",
	}

// check if there is any submission data
	if len(data.Submissions) > 0 {
// put submission text into a slice of strings
    var out []string
    for _, v := range data.Submissions {
      out = append(out, v[0])
    }

// Create a cid manually by specifying the 'prefix' parameters
		pref := cid.Prefix{
			Version:  1,
			Codec:    cid.Raw,
			MhType:   mh.SHA2_256,
			MhLength: -1, // default length
		}

		// And then feed it some data
    c, err := pref.Sum([]byte(strings.Join(out, "")))
		if err != nil {
			log.Println(err)
		}

// assign the Output the CID as a string
		data.Output = c.String()
  }

// execute template to update with Form data
	tpl.ExecuteTemplate(w, "index.gohtml", data)
}

/*
func saveHandler(w http.ResponseWriter, r *http.Request, title string) {
	body := r.FormValue("body")
	p := &Page{Title: title, Body: []byte(body)}
	err := p.save()
	if err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		return
	}
	http.Redirect(w, r, "/view/"+title, http.StatusFound)
}
*/

var tpl *template.Template

// initialize template
func init() {
	tpl = template.Must(template.ParseFiles("index.gohtml"))
}

func main() {
  http.HandleFunc("/", displayHandler) // setting router rule
	http.ListenAndServe(":8080", nil)
}
