package main

import (
	"encoding/json"
	"fmt"
	"github.com/ipfs/go-cid"
	mh "github.com/multiformats/go-multihash"
	"html/template"
	"io/ioutil"
	"log"
	"net/http"
	"net/url"
	"regexp"
	"strings"
)

// tool class struct of json object
type ToolClass struct {
	Id          string `json:"id"`
	Name        string `json:"name"`
	Description string `json:"description"`
}

// version struct
type Version struct {
	Id              string   `json:"id"`
	Name            string   `json:"name"`
	Dockerfile      bool     `json:"dockerfile"`
	Url             string   `json:"url"`
	Image           string   `json:"image"`
	VerifiedSource  string   `json:"verified-source"`
	MetaVersion     string   `json:"meta-version"`
	DescriptionType []string `json:"description-type"`
	Verified        bool     `json:"verified"`
}

// tool struct, main body of API call
type Tool struct {
	Organization   string    `json:"organization"`
	Author         string    `json:author`
	Toolname       string    `json:toolname`
	Url            string    `json:url`
	Description    string    `json:description`
	Versions       []Version `json:versions`
	VerifiedSource string    `json:"verified-source"`
	Signed         bool      `json:signed`
	MetaVersion    string    `json:"meta-version"`
	Id             string    `json:id`
	Toolclass      ToolClass `json:toolclass`
	Contains       []string  `json:contains`
	Verified       bool
}


// API call to firecloud to get all info on tools in firecloud
func getToolsURL() []string {

	url := "https://api.firecloud.org/ga4gh/v1/tools"
	// Build the request
	req, err := http.NewRequest("GET", url, nil)
	if err != nil {
		fmt.Println("NewRequest: ", err)
	}
	// For control over HTTP client headers,
	// redirect policy, and other settings,
	// create a Client
	// A Client is an HTTP client
	client := &http.Client{}

	// Send the request via a client
	// Do sends an HTTP request and
	// returns an HTTP response
	resp, err := client.Do(req)
	if err != nil {
		fmt.Println("Do: ", err)
	}
	// Callers should close resp.Body
	// when done reading from it
	// Defer the closing of the body
	defer resp.Body.Close()

	// Fill the record with the data from the JSON
	var records []Tool

	// Use json.Decode for reading streams of JSON data
	if err := json.NewDecoder(resp.Body).Decode(&records); err != nil {
		fmt.Println(err)
	}

	// https://api.firecloud.org/ga4gh/v1/tools/erictdawson:wham/versions/1/plain-WDL/descriptor
	var url_out []string
	for _, val := range records {
		url_out = append(url_out, "https://api.firecloud.org/ga4gh/v1/tools/"+val.Id+"/versions/"+val.MetaVersion+"/plain-WDL/descriptor")
	}

	return url_out

}


// just get id from url
func getID(tool_urls []string) []string {
	var out []string
	for _, v := range tool_urls {
		reg, _ := regexp.Compile(`https://api.firecloud.org/ga4gh/v1/tools/(.+?)/versions/`)
		reg_result := reg.FindAllStringSubmatch(v, -1)
		out = append(out, reg_result[0][1])
	}

	return out
}



// takes a url and returns a workflow as a string
func getWorkflow(url string) string {

	res, err := http.Get(url)
	if err != nil {
		fmt.Println(err)
	}
	body, err := ioutil.ReadAll(res.Body)
	res.Body.Close()
	if err != nil {
		fmt.Println(err)
	}
	return string(body)
}


// handler for data in go template
func displayHandler(w http.ResponseWriter, r *http.Request) {
	err := r.ParseForm()
	if err != nil {
		log.Fatalln(err)
	}

	tool_urls := getToolsURL()

	// input data from the Form
	data := struct {
		Method        string
		URL           *url.URL
		Submissions   map[string][]string
		Header        http.Header
		Host          string
		ContentLength int64
		Choice        []string
		WDL           []string
		Hash          string
	}{
		r.Method,
		r.URL,
		r.Form,
		r.Header,
		r.Host,
		r.ContentLength,
		getID(tool_urls),
		[]string{""},
		"",
	}




	// check if there is any submission data
	if len(data.Submissions) > 0 {

		// needs verified it is working properly
			var out []string
			for _, tool := range tool_urls {
				idx := strings.Index(tool, data.Submissions["wdl"][0])

				if idx > 0 {
					out = append(out, getWorkflow(tool))
				}

			}
			
		// put submission text into a slice of strings
		//var out []string
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

		// assign wdl to display?
		data.WDL = out

		// assign the Output the CID as a string
		data.Hash = c.String()
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
