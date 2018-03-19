package main

import (
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"os/exec"
)

func main() {
	files, err := ioutil.ReadDir("RawWDLs")
	if err != nil {
		fmt.Println(err)
	}

	for _, file := range files {

		valid, err := exec.Command("java", "-jar", "jars/wdltool-0.14.jar", "validate", "RawWDLs/"+file.Name()).Output()
		if err != nil {
			fmt.Printf("Command finished with error: %v", err)
		} else {
			inputs, _ := exec.Command("java", "-jar", "jars/wdltool-0.14.jar", "inputs", "RawWDLs/"+file.Name()).Output()
			if len(inputs) == 0 {
				fmt.Printf("There are no inputs to the WDL")
			} else {
				fmt.Println(string(valid) + "It Worked!")

				in, err := os.Open("RawWDLs/" + file.Name())
				if err != nil {
					fmt.Println(err)
				}

				out, err := os.Create("ValidatedWDLs/" + file.Name())
				if err != nil {
					fmt.Println(err)
				}

				_, err = io.Copy(out, in)
				if err != nil {
					fmt.Println(err)
				}

				in.Close()
				out.Close()
			}
		}
	}
}
