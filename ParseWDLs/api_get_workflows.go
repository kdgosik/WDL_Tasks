package main

import (
	"encoding/json"
	"fmt"
	"net/http"
	"io/ioutil"
	"strings"
)


// tool class struct of json object
type ToolClass struct {
	Id          		string	 	`json:"id"`
	Name        		string 		`json:"name"`
	Description 		string 		`json:"description"`
}

// version struct
type Version struct {
	Id 							string 		`json:"id"`
	Name 						string 		`json:"name"`
	Dockerfile 			bool 			`json:"dockerfile"`
	Url 						string 		`json:"url"`
	Image 					string 		`json:"image"`
	VerifiedSource 	string 		`json:"verified-source"`
	MetaVersion 		string 		`json:"meta-version"`
	DescriptionType []string 	`json:"description-type"`
	Verified 				bool 			`json:"verified"`
}

// tool struct, main body of API call
type Tool struct {
	Organization   string    	`json:"organization"`
	Author         string    	`json:author`
	Toolname       string    	`json:toolname`
	Url            string    	`json:url`
	Description    string    	`json:description`
	Versions       []Version  `json:versions`
	VerifiedSource string    	`json:"verified-source"`
	Signed         bool	    	`json:signed`
	MetaVersion    string    	`json:"meta-version"`
	Id             string    	`json:id`
	Toolclass      ToolClass 	`json:toolclass`
	Contains       []string  	`json:contains`
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
	var out []string
	for _, val := range records {
		out = append(out, "https://api.firecloud.org/ga4gh/v1/tools/"+val.Id+"/versions/"+val.MetaVersion+"/plain-WDL/descriptor")
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


// converts a workflow into a slice of strings
func workflowToSlice(workflow string) []string {

	w := strings.FieldsFunc(workflow, func(r rune) bool {
		switch r {
		case '\n', '\t', ' ':
			return true
		}
		return false
	})
	//fmt.Printf("%q\n", w)
	//fmt.Println(w)
	return w
}
