package main

import (
	"fmt"
	"regexp"
	"os"
)

func main() {

	tool_urls := getToolsURL()

	for _, url := range tool_urls {
		wf := getWorkflow(url)

		reg, _ := regexp.Compile(`\:(.+?)\/`)
		reg_result := reg.FindAllStringSubmatch(url, -1)
		path := reg_result[1][1]

		file, err := os.Create("RawWDLs/" + path + ".wdl")
		if err != nil {
			fmt.Println(err)
		}

		file.WriteString(wf)
		file.Close()
	}

}
