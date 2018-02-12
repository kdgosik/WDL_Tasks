package main

import (
	"bufio"
	"io/ioutil"
	"fmt"
  "os"
	"net/http"
	"strings"
)

func main() {

	file, err := os.Open("Methods_APIs.txt")
	if err != nil {
		fmt.Println(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		var s string

		res, err := http.Get(scanner.Text())
		if err != nil {
			fmt.Println(err)
		}

		//https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_cga:WXS_hg19_MutationCalling_QC_v1-3_BETA/versions/1/plain-WDL/descriptor
		out_path := strings.Replace(scanner.Text(), "https://api.firecloud.org/ga4gh/v1/tools/", "", 1)
		out_path = strings.Replace(out_path, "/versions/1/plain-WDL/descriptor", "", 1)

		// read
		bs, err := ioutil.ReadAll(res.Body)
		if err != nil {
			fmt.Println(err)
		}
		s = string(bs)

		OutFile, err := os.Create("workflows/" + out_path + ".txt")
		if err != nil {
			fmt.Println(err)
		}
		defer OutFile.Close()

		OutFile.WriteString(s)

	}

}
