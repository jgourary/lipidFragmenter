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


	file, err := os.Open(filePath)
	if err != nil {
		fmt.Println("Failed to open file: " + filePath)
		log.Fatal(err)
	}
	// Initialize scanner
	thisScanner := bufio.NewScanner(file)

	// initialize storage
	thisScanner.Scan()
	molName := thisScanner.Text()
	outDir := filepath2.Join(moleculesDir,molName)
	os.MkdirAll(outDir,0755)
	outPath := filepath2.Join(outDir,molName + ".sdf")
	outFile, err := os.Create(outPath)
	if err != nil {
		fmt.Println("Failed to open file: " + outPath)
		log.Fatal(err)
	}
	_, err = outFile.WriteString(molName + "\n")


	for thisScanner.Scan() {
		// get next line
		line := thisScanner.Text()
		//
		if strings.Contains(line, "$$$$") {
			thisScanner.Scan()
			molName := thisScanner.Text()
			outDir := filepath2.Join(moleculesDir,molName)
			os.MkdirAll(outDir,0755)
			outPath := filepath2.Join(outDir,molName + ".sdf")
			outFile, err = os.Create(outPath)
			if err != nil {
				fmt.Println("Failed to open file: " + outPath)
				log.Fatal(err)
			}
			_, err = outFile.WriteString(molName + "\n")

		} else {
			_, err = outFile.WriteString(line + "\n")
		}
	}
}


