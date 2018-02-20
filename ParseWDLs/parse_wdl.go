package main

import (
	"regexp"
	"strings"
)

type Input struct {
	Id   int    `json:"id"`
	Name string `json:"name"`
	Type string `json:"type"`
}

type Task struct {
	Id          int      `json:"id"`
	Name        string   `json:"name"`
	Inputs      []string `json:"inputs"`
	Command     string   `json:"command"`
	Output      string   `json:"output"`
	Runtime     string   `json:"runtime"`
	Meta        string   `json:"meta"`
	Connections []string `json:"connections"`
}

// locate WDL Component
func findWDLComponent(s string, expr string) []string {
	// locating tasks
	reg, _ := regexp.Compile(expr)
	reg_result := reg.FindAllStringSubmatch(s, -1)

	var out []string
	for _, v1 := range reg_result {
		for k2, v2 := range v1 {
			if k2 == 1 {
				tmp := strings.TrimPrefix(v2, " ")
				tmp = strings.TrimSuffix(tmp, " ")
				out = append(out, tmp)
			}
		}
	}
	return out
}


func parseInputs(input []string) [][]Input {
	var inputs [][]Input
	for _, v := range input {
		var input_slice []Input
		tmp := strings.Split(v, " ")
		for k2, _ := range tmp {
			input_obj := Input{}
			if k2%2 == 0 {
				input_obj.Id = rand.Int()
				input_obj.Name = tmp[k2+1]
				input_obj.Type = tmp[k2]
				input_slice = append(input_slice, input_obj)
			}
		}
		inputs = append(inputs, input_slice)
	}

	return inputs
}

/*
func main() {

		//f1, err := ioutil.ReadFile("../RetrieveWDLsAPI/workflows/alex_methods_new4:redact_me_brad.txt")
		//f1, err := ioutil.ReadFile("../RetrieveWDLsAPI/workflows/broadgdac-aggregate_clinical.txt")
		//f1, err := ioutil.ReadFile("../RetrieveWDLsAPI/workflows/broad_cga-mutect_pipeline_v1.txt") // command regex is not working properly
		//f1, err := ioutil.ReadFile("../RetrieveWDLsAPI/workflows/broad_cga_pipelines-mutect_pipeline_v1.txt") // gives error becuase of command <<<
		f1, err := ioutil.ReadFile("../RetrieveWDLsAPI/workflows/broadinstitute_gtex-samtools_view_v1-0_BETA.txt")
		if err != nil {
			fmt.Println(err)
		}

		xs := string(f1)
		w := strings.FieldsFunc(xs, func(r rune) bool {
			switch r {
			case '\n', '\t', ' ':
				return true
			}
			return false
		})

	tool_urls := getToolsURL()
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

	var inputs [][]Input
	for _, v := range input {
		var input_slice []Input
		tmp := strings.Split(v, " ")
		for k2, _ := range tmp {
			input_obj := Input{}
			if k2%2 == 0 {
				input_obj.Id = rand.Int()
				input_obj.Name = tmp[k2+1]
				input_obj.Type = tmp[k2]
				input_slice = append(input_slice, input_obj)
			}
		}
		inputs = append(inputs, input_slice)
	}

	fmt.Println("task: ", task)
	fmt.Println("inputs: ", input)
	fmt.Println("command: ", command)
	fmt.Println("output: ", output)
	fmt.Println("meta: ", meta)
	fmt.Println("runtime: ", runtime)
	fmt.Println("--------------------------------")
	fmt.Println("--------------------------------")
	fmt.Println("--------------------------------")

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

}
*/
