package main

import (
	"fmt"
	"io/ioutil"
	"os"
	"os/exec"
	"regexp"
	"strings"
	 "encoding/csv"
	//"encoding/json"
)


func check(e error) {
	if e != nil {
			panic(e)
	}
}

// node class struct of json object
type Node struct {
	Name        		string 		`json:"name"`
	Type            string    `json:type`
	Description 		string 		`json:"description"`
	Hash            string    `json:"hash"`
	PreviousHash    string   	`json:previousHash`
}

func main() {
	fout, err := os.Create("wdl_input_fields.csv")
	check(err)
	defer fout.Close()
	w := csv.NewWriter(fout)
	defer w.Flush()


	files, err := ioutil.ReadDir("ValidatedWDLs")
	check(err)

	var row [][]string
	for _, file := range files {
fmt.Println(file.Name())
		inputs, err := exec.Command("java", "-jar", "jars/wdltool-0.14.jar", "inputs", "ValidatedWDLs/"+file.Name()).Output()
		if err != nil {
			fmt.Printf("Command finished with error: %v", err)
		} else {
			trimmed := strings.Replace(string(inputs), "{", "", 1)
			trimmed = strings.Replace(trimmed, "}", "", 1)
			trimmed = strings.Replace(trimmed, " ", "", -1)

			tmp := strings.Split(trimmed, ",")

			for _, v := range tmp {
				var out []string
				out = append(out, file.Name())
				out = append(out, strings.TrimSpace(v))

				reg, _ := regexp.Compile(`\.([^.]*?)":`)
				reg_result := reg.FindAllStringSubmatch(v, -1)
				fmt.Println(reg_result)
				out = append(out, reg_result[0][1])

				reg2, _ := regexp.Compile(`\:"(.*?)"`)
				reg2_result := reg2.FindAllStringSubmatch(v, -1)
				out = append(out, reg2_result[0][1])
				fmt.Println(reg2_result)
				row = append(row, out)
			}

		}

	}

	for _, v := range row {
		err := w.Write(v)
		check(err)
	}

}
