package main

import (
	"fmt"
	"io/ioutil"
	"strings"
)

func main() {

	fmt.Println("-----------------------------------")
	fmt.Println("Start of Program")

	f1, err := ioutil.ReadFile("../RetrieveWDLsAPI/workflows/alex_methods_new4:redact_me_brad.txt")
	if err != nil {
		fmt.Println(err)
	}

	s := string(f1)
	w := strings.FieldsFunc(s, func(r rune) bool {
		switch r {
		case '\n', '\t', ' ':
			return true
		}
		return false
	})
	fmt.Printf("%q\n", w)
	fmt.Println(w)
}
