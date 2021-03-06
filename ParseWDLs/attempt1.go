package main

import (
	"fmt"
//	"strings"
	"regexp"
	"os"
	//	"encoding/json"
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
		defer file.Close()

		file.WriteString(wf)

	}

	/*
	   	wf := getWorkflow(tool_urls[0])
	   	w := workflowToSlice(wf)
	   	s := strings.Join(w, " ")


	   	// using regex to find task components
	   	task := findWDLComponent(s, `task (\w+)`)
	   	input := findWDLComponent(s, `{(.*?) command `)
	   	command := findWDLComponent(s, `command {(.*?)} (output|runtime|meta)`)
	   	output := findWDLComponent(s, `output {(.*?)}`)
	   	meta := findWDLComponent(s, `meta {(.*?)}`)
	   	runtime := findWDLComponent(s, `runtime {(.*?)}`)

	   // needs work
	   //	inputs := parseInputs(input)


	   	fmt.Println("task: ", task)
	   	fmt.Println("inputs: ", input)
	   	fmt.Println("command: ", command)
	   	fmt.Println("output: ", output)
	   	fmt.Println("meta: ", meta)
	   	fmt.Println("runtime: ", runtime)
	*/

	/*
		inputs_out, err := json.Marshal(inputs)
		if err != nil {
			fmt.Println(err)
		}
		fmt.Println(string(inputs_out))

		fmt.Println("--------------------------------")
		fmt.Println("--------------------------------")
		fmt.Println("--------------------------------")
		// still need to incorpate a loop for multiple tasks in a workflow
		out := Task{}
		out.Id = rand.Int()
		out.Name = task[0]
		out.Inputs = input
		out.Command = command[0]
		out.Output = output[0]

		out1, err := json.Marshal(out)
		if err != nil {
			fmt.Println(err)
		}
		fmt.Println(string(out1))

	*/
}
