package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	filepath2 "path/filepath"
	"strings"
)

func readLipidMaps(filePath string, moleculesDir string) {

	_ = os.MkdirAll(moleculesDir, 0755)

	file, err := os.Open(filePath)
	if err != nil {
		fmt.Println("Failed to open file: " + filePath)
		log.Fatal(err)
	}
	// Initialize scanner
	thisScanner := bufio.NewScanner(file)

	var molID string

	for thisScanner.Scan() {
		// get next line
		line := thisScanner.Text()
		//
		if strings.Contains(line, "<LM_ID>") {
			thisScanner.Scan()
			molID = thisScanner.Text()
		} else if strings.Contains(line, "<SMILES>") {
			thisScanner.Scan()
			molSMILES := thisScanner.Text()
			outPath := filepath2.Join(moleculesDir,molID + ".smi")
			outFile, err := os.Create(outPath)
			if err != nil {
				fmt.Println("Failed to create file: " + outPath)
				log.Fatal(err)
			}
			_, err = outFile.WriteString(molSMILES + "\n")

		}
	}
}


